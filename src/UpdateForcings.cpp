/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2020 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Model.h"
#include "Radiation.h"
#include "Properties.h"
#include "Forcings.h"
#include <limits.h>

double EstimateRelativeHumidity(const relhum_method method,
                                const force_struct &F);
double EstimateAirPressure     (const airpress_method method,
                                const force_struct &F,
                                const double       &elev);

//////////////////////////////////////////////////////////////////
/// \brief Updates HRU forcing functions
/// \details Interpolate meteorological information from Gauge stations and
///  assign to each HRU
///  Additionally estimates missing forcings from avaialable data and
///  corrects for orographic effects
///  \remark Called prior to each computational timestep and presumes constant
///  forcing functions over global time step
///
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
//
void CModel::UpdateHRUForcingFunctions(const optStruct &Options,
                                       const time_struct &tt)
{

  force_struct        F;
  static force_struct *Fg=NULL; 
  double              elev,slope;
  int                 mo,yr;
  int                 k,g,nn;
  double              mid_day,model_day, time_shift;
  double              wt;

  //Reserve static memory (only gets called once in course of simulation)
  if (Fg==NULL){
    Fg=new force_struct [_nGauges];
  }

  double t  = tt.model_time;
  mo        = tt.month;
  yr        = tt.year;
  nn        = (int)((tt.model_time+TIME_CORRECTION)/Options.timestep);//current timestep index.
 
  time_shift= Options.julian_start_day-floor(Options.julian_start_day);
  model_day = floor(tt.model_time+time_shift+TIME_CORRECTION); //model time of 00:00 of current day
  mid_day   = floor(tt.julian_day+TIME_CORRECTION)+0.5;//mid day

  CForcingGrid *pGrid_pre        = NULL;            // forcing grids
  CForcingGrid *pGrid_rain       = NULL;
  CForcingGrid *pGrid_snow       = NULL;
  CForcingGrid *pGrid_pet        = NULL;
  CForcingGrid *pGrid_tave       = NULL;
  CForcingGrid *pGrid_daily_tmin = NULL;
  CForcingGrid *pGrid_daily_tmax = NULL;
  CForcingGrid *pGrid_daily_tave = NULL;
  CForcingGrid *pGrid_recharge   = NULL;
  CForcingGrid *pGrid_precip_temp= NULL;

  // see if gridded forcing is read from a NetCDF
  bool pre_gridded            = ForcingGridIsInput(F_PRECIP)         ;
  bool rain_gridded           = ForcingGridIsInput(F_RAINFALL)       ;
  bool snow_gridded           = ForcingGridIsInput(F_SNOWFALL)       ;
  bool pet_gridded            = ForcingGridIsInput(F_PET)            ;
  bool temp_ave_gridded       = ForcingGridIsInput(F_TEMP_AVE)       ;
  bool temp_daily_min_gridded = ForcingGridIsInput(F_TEMP_DAILY_MIN) ;
  bool temp_daily_max_gridded = ForcingGridIsInput(F_TEMP_DAILY_MAX) ;
  bool temp_daily_ave_gridded = ForcingGridIsInput(F_TEMP_DAILY_AVE) ;
  bool recharge_gridded       = ForcingGridIsInput(F_RECHARGE)       ;
  bool precip_temp_gridded    = ForcingGridIsInput(F_PRECIP_TEMP)    ;

  //Extract data from gauge time series
  for (g=0;g<_nGauges;g++)
  {
    ZeroOutForcings(Fg[g]);

    if ( !(pre_gridded || snow_gridded || rain_gridded) )
    {
      if(_pGauges[g]->TimeSeriesExists(F_PRECIP)){//if precip exists, others exist
        Fg[g].precip          =_pGauges[g]->GetForcingValue(F_PRECIP,nn);     //mm/d
        Fg[g].precip_daily_ave=_pGauges[g]->GetForcingValue(F_PRECIP,model_day,1);
        Fg[g].precip_5day     =_pGauges[g]->GetForcingValue(F_PRECIP,t-5.0,5.0)*5.0;
        Fg[g].snow_frac       =_pGauges[g]->GetAverageSnowFrac(nn);
      }
    }
    if ( !(temp_daily_min_gridded && temp_daily_max_gridded) )
    {
      if(_pGauges[g]->TimeSeriesExists(F_TEMP_AVE)){//if temp_ave exists, others exist
        Fg[g].temp_ave        =_pGauges[g]->GetForcingValue(F_TEMP_AVE,nn);
        Fg[g].temp_daily_ave  =_pGauges[g]->GetForcingValue(F_TEMP_DAILY_AVE,nn);
        Fg[g].temp_daily_min  =_pGauges[g]->GetForcingValue(F_TEMP_DAILY_MIN,nn);
        Fg[g].temp_daily_max  =_pGauges[g]->GetForcingValue(F_TEMP_DAILY_MAX,nn);
        Fg[g].temp_ave_unc    =Fg[g].temp_daily_ave;
        Fg[g].temp_min_unc    =Fg[g].temp_daily_min;
        Fg[g].temp_max_unc    =Fg[g].temp_daily_max;
        if(Fg[g].temp_daily_max < Fg[g].temp_daily_min)
        {
          WriteWarning("UpdateHRUForcingFunctions: max_temp<min_temp at gauge: "+_pGauges[g]->GetName() + " on " + tt.date_string,Options.noisy);
        }
      }
    }
    Fg[g].temp_month_max  =_pGauges[g]->GetMonthlyMaxTemp  (mo);
    Fg[g].temp_month_min  =_pGauges[g]->GetMonthlyMinTemp  (mo);
    Fg[g].temp_month_ave  =_pGauges[g]->GetMonthlyAveTemp  (mo);
    Fg[g].PET_month_ave   =_pGauges[g]->GetMonthlyAvePET   (mo);

    //if (Options.uses_full_data)
    //{ // \todo [optimize] should add Options.uses_full_data to be calculated if LW,SW,etc. methods=USE_DATA or otherwise need data streams
    Fg[g].LW_radia_net    =_pGauges[g]->GetForcingValue    (F_LW_RADIA_NET,nn);
    Fg[g].SW_radia        =_pGauges[g]->GetForcingValue    (F_SW_RADIA,nn);
    Fg[g].SW_radia_net    =_pGauges[g]->GetForcingValue    (F_SW_RADIA_NET,nn);
    Fg[g].LW_incoming     =_pGauges[g]->GetForcingValue    (F_LW_INCOMING,nn);
    Fg[g].ET_radia        =_pGauges[g]->GetForcingValue    (F_ET_RADIA,nn);
    Fg[g].SW_radia_unc    =Fg[g].SW_radia;

    Fg[g].PET             =_pGauges[g]->GetForcingValue    (F_PET,nn);
    Fg[g].OW_PET          =_pGauges[g]->GetForcingValue    (F_OW_PET,nn);
    Fg[g].potential_melt  =_pGauges[g]->GetForcingValue    (F_POTENTIAL_MELT,nn);

    Fg[g].air_pres        =_pGauges[g]->GetForcingValue    (F_AIR_PRES,nn);
    Fg[g].air_dens        =_pGauges[g]->GetForcingValue    (F_AIR_DENS,nn);
    Fg[g].rel_humidity    =_pGauges[g]->GetForcingValue    (F_REL_HUMIDITY,nn);
    Fg[g].cloud_cover     =_pGauges[g]->GetForcingValue    (F_CLOUD_COVER,nn);
    Fg[g].wind_vel        =_pGauges[g]->GetForcingValue    (F_WIND_VEL,nn);
    //}

    if(!(recharge_gridded)){
      Fg[g].recharge        =_pGauges[g]->GetForcingValue    (F_RECHARGE,nn);
    }
    if(!(precip_temp_gridded)) {
      Fg[g].precip_temp     =_pGauges[g]->GetForcingValue    (F_TEMP_AVE,nn);
    }
  }
  if (_nGauges > 0) {g_debug_vars[4]=_pGauges[0]->GetElevation(); }//RFS Emulation cheat

  //Generate HRU-specific forcings from gauge data
  //---------------------------------------------------------------------
  double ref_elev_temp;
  double ref_elev_precip;
  double ref_measurement_ht; //m above land surface
  for (k = 0; k < _nHydroUnits; k++)
  {
    if(_pHydroUnits[k]->IsEnabled())
    {
      elev  = _pHydroUnits[k]->GetElevation();
      slope = _pHydroUnits[k]->GetSlope();

      ZeroOutForcings(F);
      ref_elev_temp=ref_elev_precip=0.0;
      ref_measurement_ht=0.0;

      //not gauge-based
      if(tt.day_changed)
      {
        F.day_angle  = CRadiation::DayAngle(mid_day,yr,Options.calendar);
        F.day_length = CRadiation::DayLength(_pHydroUnits[k]->GetLatRad(),CRadiation::SolarDeclination(F.day_angle));
      }

      //interpolate forcing values from gauges
      //-------------------------------------------------------------------
      for(g = 0; g < _nGauges; g++)
      {
        wt=_aGaugeWtPrecip[k][g];
        if(wt != 0.0){
          if(!(pre_gridded || snow_gridded || rain_gridded))
          {
            F.precip           += wt * Fg[g].precip;
            F.precip_daily_ave += wt * Fg[g].precip_daily_ave;
            F.precip_5day      += wt * Fg[g].precip_5day;
            F.snow_frac        += wt * Fg[g].snow_frac;
            ref_elev_precip    += wt * _pGauges[g]->GetElevation();
          }
        }
        wt=_aGaugeWtTemp[k][g];
        if(wt != 0.0) {
          if(!(temp_ave_gridded || (temp_daily_min_gridded && temp_daily_max_gridded) || temp_daily_ave_gridded))
          {
            F.temp_ave         += wt * Fg[g].temp_ave;
            F.temp_daily_ave   += wt * Fg[g].temp_daily_ave;
            F.temp_daily_min   += wt * Fg[g].temp_daily_min;
            F.temp_daily_max   += wt * Fg[g].temp_daily_max;
            F.temp_month_min   += wt * Fg[g].temp_month_min;
            F.temp_month_max   += wt * Fg[g].temp_month_max;
            F.temp_month_ave   += wt * Fg[g].temp_month_ave;
            ref_elev_temp      += wt * _pGauges[g]->GetElevation();
          }
        }
        wt=_aGaugeWeights[k][g];
        if(wt != 0.0){
          F.rel_humidity   += wt * Fg[g].rel_humidity;
          F.air_pres       += wt * Fg[g].air_pres;
          F.air_dens       += wt * Fg[g].air_dens;
          F.wind_vel       += wt * Fg[g].wind_vel;
          F.cloud_cover    += wt * Fg[g].cloud_cover;
          F.ET_radia       += wt * Fg[g].ET_radia;
          F.LW_incoming    += wt * Fg[g].LW_incoming;
          F.LW_radia_net   += wt * Fg[g].LW_radia_net;
          F.SW_radia       += wt * Fg[g].SW_radia;
          F.SW_radia_net   += wt * Fg[g].SW_radia_net;
          F.PET_month_ave  += wt * Fg[g].PET_month_ave;
          F.PET            += wt * Fg[g].PET;
          F.OW_PET         += wt * Fg[g].OW_PET;
          F.potential_melt += wt * Fg[g].potential_melt;
          if(!(recharge_gridded)){
            F.recharge       += wt * Fg[g].recharge;
          }
          if(!(precip_temp_gridded)) {
            F.precip_temp    += wt * Fg[g].precip_temp;
          }
          ref_measurement_ht+=wt*_pGauges[g]->GetMeasurementHt();
        }
      }

      //-------------------------------------------------------------------
      //  Gridded data support
      //-------------------------------------------------------------------
      bool   new_chunk1;                                // true if new chunk was read, otherwise false
      bool   new_chunk2;                                // true if new chunk was read, otherwise false
      bool   new_chunk3;                                // true if new chunk was read, otherwise false
      bool   new_chunk4;                                // true if new chunk was read, otherwise false

      //Override forcing functions with gridded data, if present
      // see if gridded forcing is available (either from NetCDF or derived)
      pre_gridded            = ForcingGridIsInput(F_PRECIP); 
      rain_gridded           = ForcingGridIsInput(F_RAINFALL);
      snow_gridded           = ForcingGridIsInput(F_SNOWFALL);
      pet_gridded            = ForcingGridIsInput(F_PET);
      temp_ave_gridded       = ForcingGridIsInput(F_TEMP_AVE);
      temp_daily_min_gridded = ForcingGridIsInput(F_TEMP_DAILY_MIN);
      temp_daily_max_gridded = ForcingGridIsInput(F_TEMP_DAILY_MAX);
      temp_daily_ave_gridded = ForcingGridIsInput(F_TEMP_DAILY_AVE);
      recharge_gridded       = ForcingGridIsInput(F_RECHARGE);
      precip_temp_gridded    = ForcingGridIsInput(F_PRECIP_TEMP);

      // find the correct grid
      if(pre_gridded)             { pGrid_pre         = GetForcingGrid(F_PRECIP); }
      if(rain_gridded)            { pGrid_rain        = GetForcingGrid(F_RAINFALL); }
      if(snow_gridded)            { pGrid_snow        = GetForcingGrid(F_SNOWFALL); }
      if(pet_gridded)             { pGrid_pet         = GetForcingGrid(F_PET); }
      if(temp_ave_gridded)        { pGrid_tave        = GetForcingGrid(F_TEMP_AVE); }
      if(temp_daily_min_gridded)  { pGrid_daily_tmin  = GetForcingGrid(F_TEMP_DAILY_MIN); }
      if(temp_daily_max_gridded)  { pGrid_daily_tmax  = GetForcingGrid(F_TEMP_DAILY_MAX); }
      if(temp_daily_ave_gridded)  { pGrid_daily_tave  = GetForcingGrid(F_TEMP_DAILY_AVE); }
      if(recharge_gridded)        { pGrid_recharge    = GetForcingGrid(F_RECHARGE); }
      if(precip_temp_gridded)     { pGrid_precip_temp = GetForcingGrid(F_PRECIP_TEMP); }
      // ---------------------
      // (1A) read gridded precip/snowfall/rainfall and populate additional time series
      // ---------------------
      if(pre_gridded || snow_gridded || rain_gridded)
      {
        // read data (actually new chunk is only read if timestep is not covered by old chunk anymore)
        new_chunk1 = false; new_chunk2 = false; new_chunk3 = false;
        if(pre_gridded)  { new_chunk1 = pGrid_pre-> ReadData(Options,tt.model_time); }
        if(snow_gridded) { new_chunk2 = pGrid_snow->ReadData(Options,tt.model_time); }
        if(rain_gridded) { new_chunk3 = pGrid_rain->ReadData(Options,tt.model_time); }

        // populate derived data from the ones just read (e.g., precip -> (rain,snow), or (rain,snow)->(precip)
        if(new_chunk1 || new_chunk2 || new_chunk3) { 
          GenerateGriddedPrecipVars(Options);//only call if new data chunk found
        } 

        pGrid_pre  = GetForcingGrid((F_PRECIP)); 
        pGrid_rain = GetForcingGrid((F_RAINFALL)); 
        pGrid_snow = GetForcingGrid((F_SNOWFALL)); 

        F.precip           = pGrid_pre->GetWeightedValue           (k,tt.model_time,Options.timestep);
        F.precip_daily_ave = pGrid_pre->GetDailyWeightedValue      (k,tt.model_time,Options.timestep,Options);
        F.snow_frac        = pGrid_snow->GetWeightedAverageSnowFrac(k,tt.model_time,Options.timestep,pGrid_rain);
        F.precip_5day      = NETCDF_BLANK_VALUE;

        ref_elev_precip    = pGrid_pre->GetRefElevation(k);
        if (ref_elev_precip==RAV_BLANK_DATA){ ref_elev_precip = _pHydroUnits[k]->GetElevation();} //disabling orographic effects if no elevation given (warning in Forcing grid init)
      }

      // ---------------------
      // (2a) read gridded temperature (average or min/max) and populate additional time series
      // ---------------------
      if(temp_ave_gridded || (temp_daily_min_gridded && temp_daily_max_gridded) || temp_daily_ave_gridded)
      {
        // read data (actually new chunk is only read if timestep is not covered by old chunk anymore)
        new_chunk1 = new_chunk2 = new_chunk3 = new_chunk4 = false;
        if(temp_ave_gridded)       { new_chunk1 =       pGrid_tave->ReadData(Options,tt.model_time); }
        if(temp_daily_ave_gridded) { new_chunk2 = pGrid_daily_tave->ReadData(Options,tt.model_time); }
        if(temp_daily_min_gridded) { new_chunk3 = pGrid_daily_tmin->ReadData(Options,tt.model_time); }
        if(temp_daily_max_gridded) { new_chunk4 = pGrid_daily_tmax->ReadData(Options,tt.model_time); }

        // populate derived data from the ones just read
        if(new_chunk1 || new_chunk2 || new_chunk3 || new_chunk4) { 
          GenerateGriddedTempVars(Options);//only generate if new data chunk found
        }
        if (new_chunk3!=new_chunk4){
          ExitGracefully("CModel::UpdateHRUForcingFunctions: Gridded Min and max temperature have to have same time discretization.",BAD_DATA);
        }

        pGrid_tave       = GetForcingGrid((F_TEMP_AVE)); 
        pGrid_daily_tave = GetForcingGrid((F_TEMP_DAILY_AVE)); 
        pGrid_daily_tmin = GetForcingGrid((F_TEMP_DAILY_MIN)); 
        pGrid_daily_tmax = GetForcingGrid((F_TEMP_DAILY_MAX)); 

        F.temp_ave         = pGrid_tave      ->GetWeightedValue(k,tt.model_time,Options.timestep);
        F.temp_daily_ave   = pGrid_daily_tave->GetWeightedValue(k,tt.model_time,Options.timestep);
        F.temp_daily_min   = pGrid_daily_tmin->GetWeightedValue(k,tt.model_time,Options.timestep);
        F.temp_daily_max   = pGrid_daily_tmax->GetWeightedValue(k,tt.model_time,Options.timestep);
        
				F.temp_month_ave   = NOT_SPECIFIED;
				F.temp_month_min   = NOT_SPECIFIED;
				F.temp_month_max   = NOT_SPECIFIED;

        ref_elev_temp      = pGrid_pre->GetRefElevation(k);
        if(ref_elev_temp==RAV_BLANK_DATA) { ref_elev_temp = _pHydroUnits[k]->GetElevation(); } //disabling orographic effects if no elevation given (warning in Forcing grid init)
      }

      // ---------------------
      // (3) read gridded recharge
      // ---------------------
      if(recharge_gridded)
      {
        pGrid_recharge   = GetForcingGrid(F_RECHARGE); 
        // read data (actually new chunk is only read if timestep is not covered by old chunk anymore)
        pGrid_recharge-> ReadData(Options,tt.model_time);

        F.recharge   = pGrid_recharge->GetWeightedValue(k,tt.model_time,Options.timestep);
      }
      // ---------------------
      // (3) read gridded precipitation temperature
      // ---------------------
      if(precip_temp_gridded)
      {
        pGrid_precip_temp   = GetForcingGrid(F_PRECIP_TEMP);
        // read data (actually new chunk is only read if timestep is not covered by old chunk anymore)
        pGrid_precip_temp-> ReadData(Options,tt.model_time);

        F.precip_temp   = pGrid_precip_temp->GetWeightedValue(k,tt.model_time,Options.timestep);
      }


      //-------------------------------------------------------------------
      //  Temperature Corrections
      //-------------------------------------------------------------------
      F.temp_ave_unc = F.temp_daily_ave;
      F.temp_min_unc = F.temp_daily_min;
      F.temp_max_unc = F.temp_daily_max;
      
      CorrectTemp(Options,F,elev,ref_elev_temp,tt);

      //-------------------------------------------------------------------
      //  Copy Daily values from current day, earlier time steps
      //    (done after temperature corrections so that the uncorrected values are used in calculating the lapse rate for temp_ave)
      //-------------------------------------------------------------------
      if(!tt.day_changed){
        _pHydroUnits[k]->CopyDailyForcings(F);
      }

      //-------------------------------------------------------------------
      //  Subdaily Corrections
      //-------------------------------------------------------------------
      // always calculated, never from gauge
      F.subdaily_corr = CalculateSubDailyCorrection(F,Options,elev,ref_elev_temp,tt,k);

      //-------------------------------------------------------------------
      //  Air Pressure, Density, relative humidity
      //-------------------------------------------------------------------
      F.air_pres = EstimateAirPressure(Options.air_pressure,F,elev);

      F.air_dens = GetAirDensity(F.temp_ave,F.air_pres);

      F.rel_humidity = EstimateRelativeHumidity(Options.rel_humidity,F);

      //-------------------------------------------------------------------
      // Snow fraction Calculations
      //-------------------------------------------------------------------
      F.snow_frac =EstimateSnowFraction(Options.rainsnow,&F,Options);
      
      //-------------------------------------------------------------------
      //  Precip Corrections
      //-------------------------------------------------------------------

      //--Gauge Corrections------------------------------------------------
      if(!(pre_gridded || snow_gridded || rain_gridded)) //Gauge Data
      {        
				double gauge_corr; double rc,sc;
				F.precip=F.precip_5day=F.precip_daily_ave=0.0;
        int p=_pHydroUnits[k]->GetSubBasinIndex();
        rc=_pSubBasins[p]->GetRainCorrection();
        sc=_pSubBasins[p]->GetSnowCorrection();
        for(int g=0; g<_nGauges; g++)
        {
          gauge_corr= F.snow_frac*sc*_pGauges[g]->GetSnowfallCorr() + (1.0-F.snow_frac)*rc*_pGauges[g]->GetRainfallCorr();
          wt=_aGaugeWtPrecip[k][g];

          F.precip         += wt*gauge_corr*Fg[g].precip;
          F.precip_daily_ave+=wt*gauge_corr*Fg[g].precip_daily_ave;
          F.precip_5day    += wt*gauge_corr*Fg[g].precip_5day;
        }
      }
      else //Gridded Data
      {
				double grid_corr;
        double rain_corr=pGrid_pre->GetRainfallCorr();
        double snow_corr=pGrid_pre->GetSnowfallCorr();
        grid_corr= F.snow_frac*snow_corr + (1.0-F.snow_frac)*rain_corr;
				
        F.precip          *=grid_corr;
        F.precip_daily_ave*=grid_corr;
        F.precip_5day      =NETCDF_BLANK_VALUE;
      }

      //--Orographic corrections-------------------------------------------
      CorrectPrecip(Options,F,elev,ref_elev_precip,k,tt);

      //-------------------------------------------------------------------
      //  Wind Velocity
      //-------------------------------------------------------------------

      F.wind_vel = EstimateWindVelocity(Options,F,ref_measurement_ht,k);
      
      //-------------------------------------------------------------------
      //  Cloud Cover
      //-------------------------------------------------------------------

      F.cloud_cover = EstimateCloudCover(Options,F,k);

      //-------------------------------------------------------------------
      //  Radiation Calculations
      //-------------------------------------------------------------------

      F.SW_radia = CRadiation::EstimateShortwaveRadiation(Options,&F,_pHydroUnits[k],tt,F.ET_radia);
      F.SW_radia_unc = F.SW_radia;
      F.SW_radia *= CRadiation::SWCloudCoverCorrection(Options,&F,elev);
      F.SW_radia *= CRadiation::SWCanopyCorrection(Options,_pHydroUnits[k]);

      if (Options.SW_radia_net == NETSWRAD_CALC) 
      {
        F.SW_radia_net = F.SW_radia*(1 - _pHydroUnits[k]->GetTotalAlbedo());
      }//otherwise, uses data

      F.LW_radia_net=CRadiation::EstimateLongwaveRadiation(GetStateVarIndex(SNOW),Options,&F,_pHydroUnits[k],F.LW_incoming);

      //-------------------------------------------------------------------
      //  Potential Melt Rate
      //-------------------------------------------------------------------

      F.potential_melt=EstimatePotentialMelt(&F,Options,_pHydroUnits[k],tt);

      //-------------------------------------------------------------------
      //  PET Calculations
      //-------------------------------------------------------------------
      // last but not least - needs all of the forcing params calculated above

      F.PET   =EstimatePET(F,_pHydroUnits[k],ref_measurement_ht,ref_elev_temp,Options.evaporation   ,Options,tt,false);
      F.OW_PET=EstimatePET(F,_pHydroUnits[k],ref_measurement_ht,ref_elev_temp,Options.ow_evaporation,Options,tt,true);

      CorrectPET(Options,F,_pHydroUnits[k],elev,ref_elev_temp,k);

      //-------------------------------------------------------------------
      // Direct evaporation of rainfall
      //-------------------------------------------------------------------
      if(Options.direct_evap) {
        double reduce=min(F.PET,F.precip*(1.0-F.snow_frac));
        if(F.precip-reduce>0.0) { F.snow_frac=1.0-(F.precip*(1.0-F.snow_frac)-reduce)/(F.precip-reduce); }
        F.precip-=reduce;
        F.PET   -=reduce;
      }

      //-------------------------------------------------------------------
      // Update
      //-------------------------------------------------------------------
      _pHydroUnits[k]->UpdateForcingFunctions(F);

    }//end if (!_pHydroUnits[k]->IsDisabled())
  }//end for k=0; k<nHRUs...

   //delete static arrays (only called once)=========================
  if(t>=Options.duration-Options.timestep)
  {
    if(DESTRUCTOR_DEBUG) { cout<<"DELETING STATIC ARRAY IN UPDATEHRUFORCINGFUNCTIONS"<<endl; }
    delete [] Fg; Fg=NULL;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Estimates air pressure given elevation [kPa]
/// \param method [in] Method of calculating air pressure
/// \param &F [in] Forcing functons structure
/// \param &elev [in] Elevation at which air pressure is to be estimated
/// \return Estimated air pressure at elev [kPa]
//
double EstimateAirPressure     (const airpress_method method,
                                const force_struct &F,
                                const double       &elev)
{
  if (method==AIRPRESS_DATA){
    return F.air_pres;
  }
  else if (method==AIRPRESS_BASIC)
  {
    return AMBIENT_AIR_PRESSURE*pow(1.0-0.0065*elev/(ZERO_CELSIUS+F.temp_ave),5.26);
  }
  else if (method==AIRPRESS_UBC)
  {
    return KPA_PER_ATM*(1.0-0.0001*elev);//UBC_WM (c) Michael Quick
  }
  else if (method==AIRPRESS_CONST){
    return AMBIENT_AIR_PRESSURE;
  }
  else{
    return AMBIENT_AIR_PRESSURE;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Estimates relative humidity = ea/esat
/// \param method [in] Method of calculating relative humidity
/// \param &F [in] Reference to forcing function structure
/// \return Estimated relative humidity
//
double EstimateRelativeHumidity(const relhum_method method,
                                const force_struct &F)
{
  if (method==RELHUM_CONSTANT)
  {
    return 0.5;
  }
  else if (method==RELHUM_MINDEWPT)
  { //uses minimum daily temperature as proxy for dew point temperature
    double dew_point_temp=F.temp_daily_min;
    return min(GetSaturatedVaporPressure(dew_point_temp)/GetSaturatedVaporPressure(F.temp_ave),1.0);
  }
  else if (method==RELHUM_DATA)
  {
    return F.rel_humidity;
  }
  return 0.5;
}

//////////////////////////////////////////////////////////////////
/// \brief Estimates wind velocity [m/s]
/// \details UBC method adapted from UBC Watershed model source code,
/// \copyright (c) Michael Quick
///
/// \param &Options [in] Global model options information
/// \param &F [in] Reference to forcing functions for HRU k
/// \param k [in] index of HRU
/// \return Estimated wind velocity in HRU k [m/s]
double CModel::EstimateWindVelocity(const optStruct    &Options,
                                    const force_struct &F,
                                    const double       &wind_measurement_ht,
                                    const int           k)
{
  //---------------------------------------------------------------------
  if (Options.wind_velocity==WINDVEL_CONSTANT)
  {
    return 2.0;//m/s (global average)
  }
  //---------------------------------------------------------------------
  else if (Options.wind_velocity==WINDVEL_DATA)
  {
    return F.wind_vel;
  }
  //---------------------------------------------------------------------
  else if (Options.wind_velocity==WINDVEL_UBCWM)
  {
    double TED, A1term;
    TED=max(F.temp_daily_max-F.temp_daily_min,0.0);

    const double rf_elev=2000;
    double elev=_pHydroUnits[k]->GetElevation();
    double Fc  =_pHydroUnits[k]->GetSurfaceProps()->forest_coverage;
    double P0TEDL=CGlobalParams::GetParams()->UBC_lapse_params.P0TEDL;
    double P0TEDU=CGlobalParams::GetParams()->UBC_lapse_params.P0TEDU;
    double A0term=CGlobalParams::GetParams()->UBC_lapse_params.max_range_temp;
    if (elev>=rf_elev)
    {
      A1term=25.0-P0TEDL*0.001*rf_elev-P0TEDU*0.001*(elev-rf_elev);
    }
    else{
      A1term=25.0-P0TEDL*0.001*elev;
    }
    lowerswap(A1term,A0term);

    const double max_wind_speed=8.0;//P0VBMX [km/h]
    lowerswap(TED,A1term);

    double wt=min(TED/25.0,1.0);

    double wind_vel = (1-wt)*max_wind_speed +(wt)*1;

    upperswap(wind_vel,1.0);
    lowerswap(wind_vel,max_wind_speed - 1.0); //why the minus one? who knows...

    //elevation correction
    wind_vel*=max(sqrt(elev/1000.0),1.0);

    //forest correction
    const double F0WIND=0.7; //ratio of wind in forest vs wind in open
    wind_vel*=((Fc)*F0WIND+(1.0-Fc)*1.0);

    wind_vel*=M_PER_KM/SEC_PER_HR; //KPH_TO_M_PER_S;
    return wind_vel;
  }
  else
  {
    return 2.0;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns estimate of cloud cover
/// \param &Options [in] Global model options information
/// \param &F [in] Reference to forcing functions for HRU
/// \param k [in] index of HRU
/// \return Estimated cloud cover [0..1] in HRU k
//
double CModel::EstimateCloudCover  (const optStruct &Options,
                                    const force_struct &F,
                                    const int k)
{
  double cover=0.0;
  switch(Options.cloud_cover)
  {
  case(CLOUDCOV_NONE):
  {
    cover=0.0;
    break;
  }
  case(CLOUDCOV_DATA):
  {
    cover=F.cloud_cover;
    break;
  }
  case (CLOUDCOV_UBCWM):
  {
    double range=(F.temp_max_unc-F.temp_min_unc); //uses uncorrected station temperature
    double cloud_min_range(0.0),cloud_max_range(0.0);

    for (int g=0;g<_nGauges;g++){
      cloud_min_range+=_aGaugeWtTemp[k][g]*_pGauges[g]->GetCloudMinRange();//[C] A0FOGY in UBC_WM
      cloud_max_range+=_aGaugeWtTemp[k][g]*_pGauges[g]->GetCloudMaxRange();//[C] A0SUNY in UBC_WM
    }
    cover=1.0-(range-cloud_min_range)/(cloud_max_range-cloud_min_range);
    lowerswap(cover,1.0);
    upperswap(cover,0.0);
    break;
  }
  }
  return cover;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns sub-daily correction for daily snowmelt or PET calculations
/// \param &Options [in] Global model options information
/// \param &F [in] Reference to forcing functions for HRU
/// \param &elev [in] elevation of HRU
/// \param &ref_elev_temp [in] reference elevation of temperature gauge
/// \param &tt [in] current time structure
/// \param k [in] index of HRU
/// \return sub-daily correction for daily snowmelt or PET calculations
//
double CModel::CalculateSubDailyCorrection(const force_struct &F,
                                           const optStruct    &Options,
                                           const double       &elev,
                                           const double       &ref_elev_temp,
                                           const time_struct  &tt,
                                           const int          k)
{
  if (Options.timestep>=1.0){return 1.0;}

  //-----------------------------------------------------
  if (Options.subdaily==SUBDAILY_NONE)
  {
    return 1.0;
  }
  //-----------------------------------------------------
  else if (Options.subdaily==SUBDAILY_SIMPLE)
  { //tested and working - simple and elegant
    double DL=F.day_length;
    double dawn=0.5-0.5*DL;
    double dusk=0.5+0.5*DL;
    double t=tt.julian_day-floor(tt.julian_day); //time of day [d]
    double dt=Options.timestep;

    if      ((t>dawn) && (t+dt<=dusk)){
      return -0.5*(cos(PI*(t+dt-dawn)/DL)-cos(PI*(t-dawn)/DL))/Options.timestep;
    }
    else if ((t<dawn) && (t+dt>=dawn)){
      return -0.5*(cos(PI*(t+dt-dawn)/DL)-1                  )/Options.timestep;
    }
    else if ((t<dusk) && (t+dt>=dusk)){
      return -0.5*(-1                   -cos(PI*(t-dawn)/DL))/Options.timestep;
    }
    return 0.0;
  }
  //-----------------------------------------------------
  else if (Options.subdaily==SUBDAILY_UBC)
  {
    if ( ForcingGridIsAvailable(F_PRECIP)         ||
         ForcingGridIsAvailable(F_TEMP_AVE)       ||
         ForcingGridIsAvailable(F_TEMP_DAILY_MIN) ||
         ForcingGridIsAvailable(F_TEMP_DAILY_MAX)  ) {
      ExitGracefully( "CModel::CalculateSubDailyCorrection: Option SUBDAILY_UBC is not implemented when gridded inputs are given.", BAD_DATA);
    }
    //this is not pretty (and somewhat expensive), due to the need to correct all daily temperatures for every timestep, but it works
    double time_shift=Options.julian_start_day-floor(Options.julian_start_day)+TIME_CORRECTION;
    int nn_start=(int)(floor(tt.model_time+time_shift    )/Options.timestep);//index of first timestp in day
    int nn_end  =(int)(floor(tt.model_time+time_shift+1.0)/Options.timestep);//index of first timestp in following day
    force_struct Ftmp;
    time_struct  tt_tmp;
    double start_of_day;
    tt_tmp=tt;
    tt_tmp.day_changed=true;

    double sum=0.0;
    for (int nnn=nn_start;nnn<nn_end;nnn++)
    {
      tt_tmp.model_time=nnn*Options.timestep;

      start_of_day=floor(tt_tmp.model_time+time_shift);
      ZeroOutForcings(Ftmp);
      for (int g=0;g<_nGauges;g++)
      {
        Ftmp.precip_daily_ave+=_aGaugeWtPrecip[k][g]*_pGauges[g]->GetForcingValue(F_PRECIP,start_of_day,1);

        Ftmp.temp_ave        +=_aGaugeWtTemp[k][g]*_pGauges[g]->GetForcingValue(F_TEMP_AVE,nnn);
        Ftmp.temp_daily_max  +=_aGaugeWtTemp[k][g]*_pGauges[g]->GetForcingValue(F_TEMP_DAILY_MAX,nnn);
        Ftmp.temp_daily_min  +=_aGaugeWtTemp[k][g]*_pGauges[g]->GetForcingValue(F_TEMP_DAILY_MIN,nnn);
      }
      CorrectTemp(Options,Ftmp,elev,ref_elev_temp,tt_tmp);
      sum+=max(Ftmp.temp_ave,0.0);
    }

    if (sum==0){return 0.0;}
    else       {return max(F.temp_ave,0.0)/sum/Options.timestep;}
  }
  return 1.0;
}

//////////////////////////////////////////////////////////////////
/// \brief Checks whether a forcing grid (by name) is read actually
///        from a NetCDF or is a derived grid.
///        In case grid does not exist false is returned.
///
/// \param ftype [in]  type of forcing grid, e.g. F_TEMP_DAILY_MIN
//
bool CModel::ForcingGridIsInput(const forcing_type &ftype) const
{
  bool is_input = false; // class variable determining if grid is read from NetCDF (false) or is derived (true)

  int f=GetForcingGridIndexFromType(ftype);
  if(f != DOESNT_EXIST) {
    is_input = !(_pForcingGrids[f]->GetIsDerived());
  }

  return is_input;
}

//////////////////////////////////////////////////////////////////
/// \brief Checks whether a forcing grid (by type) is available (read from NetCDF or derived)
///        or not (means it is read from a non-gridded time series).
///
/// \param ftype [in]  type of forcing grid, e.g. "F_TEMP_DAILY_MIN"
//
bool CModel::ForcingGridIsAvailable(const forcing_type &ftype) const
{
  return (GetForcingGridIndexFromType(ftype) != DOESNT_EXIST);
}
