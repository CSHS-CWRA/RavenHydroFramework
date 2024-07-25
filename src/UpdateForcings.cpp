/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2024 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Model.h"
#include "Radiation.h"
#include "Properties.h"
#include "Forcings.h"
#include <limits.h>

double EstimateRelativeHumidity(const relhum_method method,
                                const force_struct &F,
                                const CHydroUnit   *pHRU);
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
  double              elev;
  int                 mo,yr;
  int                 k,g,nn;
  double              mid_day,model_day, time_shift;
  double              wt;
  bool                rvt_file_provided = (strcmp(Options.rvt_filename.c_str(), "") != 0);

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
  CForcingGrid *pGrid_owpet      = NULL;
  CForcingGrid *pGrid_windspeed  = NULL;
  CForcingGrid *pGrid_relhum     = NULL;
  CForcingGrid *pGrid_SW_net     = NULL;
  CForcingGrid *pGrid_LW_inc     = NULL;
  CForcingGrid *pGrid_SW         = NULL;
  CForcingGrid *pGrid_tave       = NULL;
  CForcingGrid *pGrid_daily_tmin = NULL;
  CForcingGrid *pGrid_daily_tmax = NULL;
  CForcingGrid *pGrid_daily_tave = NULL;
  CForcingGrid *pGrid_recharge   = NULL;
  CForcingGrid *pGrid_precip_temp= NULL;
  CForcingGrid *pGrid_precip_conc= NULL;

  // see if gridded forcing is read from a NetCDF
  bool pre_gridded            = ForcingGridIsInput(F_PRECIP)         ;
  bool rain_gridded           = ForcingGridIsInput(F_RAINFALL)       ;
  bool snow_gridded           = ForcingGridIsInput(F_SNOWFALL)       ;
  bool temp_ave_gridded       = ForcingGridIsInput(F_TEMP_AVE)       ;
  bool temp_daily_min_gridded = ForcingGridIsInput(F_TEMP_DAILY_MIN) ;
  bool temp_daily_max_gridded = ForcingGridIsInput(F_TEMP_DAILY_MAX) ;
  bool temp_daily_ave_gridded = ForcingGridIsInput(F_TEMP_DAILY_AVE) ;
  bool precip_temp_gridded    = ForcingGridIsInput(F_PRECIP_TEMP)    ;
  bool precip_conc_gridded    = ForcingGridIsInput(F_PRECIP_CONC)    ;

  bool pet_gridded            = ForcingGridIsInput(F_PET)            && (Options.evaporation   ==PET_DATA);
  bool owpet_gridded          = ForcingGridIsInput(F_OW_PET)         && (Options.ow_evaporation==PET_DATA);
  bool windspeed_gridded      = ForcingGridIsInput(F_WIND_VEL)       && (Options.wind_velocity ==WINDVEL_DATA);
  bool relhum_gridded         = ForcingGridIsInput(F_REL_HUMIDITY)   && (Options.rel_humidity  ==RELHUM_DATA);
  bool SWnet_gridded          = ForcingGridIsInput(F_SW_RADIA_NET)   && (Options.SW_radia_net  ==NETSWRAD_DATA);
  bool LWinc_gridded          = ForcingGridIsInput(F_LW_INCOMING)    && (Options.LW_incoming   ==LW_INC_DATA);
  bool SW_gridded             = ForcingGridIsInput(F_SW_RADIA)       && (Options.SW_radiation  ==SW_RAD_DATA);
  bool recharge_gridded       = ForcingGridIsInput(F_RECHARGE)       && (Options.recharge      ==RECHARGE_DATA);

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

    Fg[g].LW_radia_net    =_pGauges[g]->GetForcingValue    (F_LW_RADIA_NET,nn);
    Fg[g].SW_radia        =_pGauges[g]->GetForcingValue    (F_SW_RADIA,nn);
    Fg[g].SW_radia_net    =_pGauges[g]->GetForcingValue    (F_SW_RADIA_NET,nn);
    Fg[g].SW_radia_subcan =_pGauges[g]->GetForcingValue    (F_SW_RADIA_SUBCAN,nn);
    Fg[g].SW_subcan_net   =_pGauges[g]->GetForcingValue    (F_SW_SUBCAN_NET,nn);
    Fg[g].LW_incoming     =_pGauges[g]->GetForcingValue    (F_LW_INCOMING,nn);
    Fg[g].ET_radia        =_pGauges[g]->GetForcingValue    (F_ET_RADIA,nn);
    Fg[g].SW_radia_unc    =Fg[g].SW_radia;
    Fg[g].ET_radia_flat   =Fg[g].ET_radia;

    Fg[g].PET             =_pGauges[g]->GetForcingValue    (F_PET,nn);
    Fg[g].OW_PET          =_pGauges[g]->GetForcingValue    (F_OW_PET,nn);
    Fg[g].potential_melt  =_pGauges[g]->GetForcingValue    (F_POTENTIAL_MELT,nn);

    Fg[g].air_pres        =_pGauges[g]->GetForcingValue    (F_AIR_PRES,nn);
    Fg[g].air_dens        =_pGauges[g]->GetForcingValue    (F_AIR_DENS,nn);
    Fg[g].rel_humidity    =_pGauges[g]->GetForcingValue    (F_REL_HUMIDITY,nn);
    Fg[g].cloud_cover     =_pGauges[g]->GetForcingValue    (F_CLOUD_COVER,nn);
    Fg[g].wind_vel        =_pGauges[g]->GetForcingValue    (F_WIND_VEL,nn);
    Fg[g].recharge        =_pGauges[g]->GetForcingValue    (F_RECHARGE,nn);
    Fg[g].precip_temp     =_pGauges[g]->GetForcingValue    (F_TEMP_AVE,nn);
    Fg[g].precip_conc     =_pGauges[g]->GetForcingValue    (F_PRECIP_CONC,nn);
  }
  if (_nGauges > 0) {g_debug_vars[4]=_pGauges[0]->GetElevation(); }//UBCWM RFS Emulation cheat

  //Generate HRU-specific forcings from gauge data
  //---------------------------------------------------------------------
  double ref_elev_temp;
  double ref_elev_precip;
  double ref_measurement_ht; //m above land surface
  for (k = 0; k < _nHydroUnits; k++)
  {
    elev  = _pHydroUnits[k]->GetElevation();

    ZeroOutForcings(F);
    ref_elev_temp=ref_elev_precip=0.0;
    ref_measurement_ht=0.0;

    //not gauge-based
    if(tt.day_changed)
    {
      F.day_angle  = CRadiation::DayAngle(mid_day,yr,Options.calendar);
      F.day_length = CRadiation::DayLength(_pHydroUnits[k]->GetLatRad(),CRadiation::SolarDeclination(F.day_angle));
    }

    if(_pHydroUnits[k]->IsEnabled())
    {
      ApplyLocalParamOverrrides(k,false);

      //interpolate forcing values from gauges
      //-------------------------------------------------------------------
      for(g = 0; g < _nGauges; g++)
      {
        wt=_aGaugeWtPrecip[k][g];
        if(wt != 0.0) {
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
        if(wt != 0.0) {
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
          F.SW_radia_subcan+= wt * Fg[g].SW_radia_subcan;
          F.SW_subcan_net  += wt * Fg[g].SW_subcan_net;
          F.PET_month_ave  += wt * Fg[g].PET_month_ave;
          F.potential_melt += wt * Fg[g].potential_melt;
          F.PET            += wt * Fg[g].PET;
          F.OW_PET         += wt * Fg[g].OW_PET;
          F.recharge       += wt * Fg[g].recharge;
          F.precip_temp    += wt * Fg[g].precip_temp;
          F.precip_conc    += wt * Fg[g].precip_conc;
          ref_measurement_ht+=wt*_pGauges[g]->GetMeasurementHt();
        }
      }

      // if in BMI without RVT file, precip and temp values are expected to have been given before this point
      if (Options.in_bmi_mode && !rvt_file_provided) {
        F.precip           = _pHydroUnits[k]->GetForcingFunctions()->precip + 0.0;
        F.precip_daily_ave = _pHydroUnits[k]->GetForcingFunctions()->precip + 0.0;  // TODO: check
        F.precip_5day      = F.precip * 5;                                          // TODO: check
        F.temp_ave         = _pHydroUnits[k]->GetForcingFunctions()->temp_ave + 0.0;
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
      owpet_gridded          = ForcingGridIsInput(F_OW_PET);
      windspeed_gridded      = ForcingGridIsInput(F_WIND_VEL);
      relhum_gridded         = ForcingGridIsInput(F_REL_HUMIDITY);
      SWnet_gridded          = ForcingGridIsInput(F_SW_RADIA_NET);
      LWinc_gridded          = ForcingGridIsInput(F_LW_INCOMING);
      SW_gridded             = ForcingGridIsInput(F_SW_RADIA);
      temp_ave_gridded       = ForcingGridIsInput(F_TEMP_AVE);
      temp_daily_min_gridded = ForcingGridIsInput(F_TEMP_DAILY_MIN);
      temp_daily_max_gridded = ForcingGridIsInput(F_TEMP_DAILY_MAX);
      temp_daily_ave_gridded = ForcingGridIsInput(F_TEMP_DAILY_AVE);
      recharge_gridded       = ForcingGridIsInput(F_RECHARGE);
      precip_temp_gridded    = ForcingGridIsInput(F_PRECIP_TEMP);
      precip_conc_gridded    = ForcingGridIsInput(F_PRECIP_CONC);

      // find the correct grid
      if(pre_gridded)             { pGrid_pre         = GetForcingGrid(F_PRECIP); }
      if(rain_gridded)            { pGrid_rain        = GetForcingGrid(F_RAINFALL); }
      if(snow_gridded)            { pGrid_snow        = GetForcingGrid(F_SNOWFALL); }
      if(pet_gridded)             { pGrid_pet         = GetForcingGrid(F_PET); }
      if(owpet_gridded)           { pGrid_owpet       = GetForcingGrid(F_OW_PET); }
      if(windspeed_gridded)       { pGrid_windspeed   = GetForcingGrid(F_WIND_VEL); }
      if(relhum_gridded)          { pGrid_relhum      = GetForcingGrid(F_REL_HUMIDITY); }
      if(SWnet_gridded)           { pGrid_SW_net      = GetForcingGrid(F_SW_RADIA_NET); }
      if(LWinc_gridded)           { pGrid_LW_inc      = GetForcingGrid(F_LW_INCOMING); }
      if(SW_gridded)              { pGrid_SW          = GetForcingGrid(F_SW_RADIA); }
      if(temp_ave_gridded)        { pGrid_tave        = GetForcingGrid(F_TEMP_AVE); }
      if(temp_daily_min_gridded)  { pGrid_daily_tmin  = GetForcingGrid(F_TEMP_DAILY_MIN); }
      if(temp_daily_max_gridded)  { pGrid_daily_tmax  = GetForcingGrid(F_TEMP_DAILY_MAX); }
      if(temp_daily_ave_gridded)  { pGrid_daily_tave  = GetForcingGrid(F_TEMP_DAILY_AVE); }
      if(recharge_gridded)        { pGrid_recharge    = GetForcingGrid(F_RECHARGE); }
      if(precip_temp_gridded)     { pGrid_precip_temp = GetForcingGrid(F_PRECIP_TEMP); }
      if(precip_conc_gridded)     { pGrid_precip_conc = GetForcingGrid(F_PRECIP_CONC); }

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

        F.precip           = pGrid_pre->GetWeightedValue(k,tt.model_time,Options.timestep);
        F.precip_daily_ave = pGrid_pre->GetDailyWeightedValue(k,tt.model_time,Options.timestep,Options);
        F.snow_frac        = pGrid_snow->GetWeightedAverageSnowFrac(k,tt.model_time,Options.timestep,pGrid_rain);
        F.precip_5day      = NETCDF_BLANK_VALUE;

        ref_elev_precip    = pGrid_pre->GetRefElevation(k);
        if(ref_elev_precip==RAV_BLANK_DATA) { ref_elev_precip = _pHydroUnits[k]->GetElevation(); } //disabling orographic effects if no elevation given (warning in Forcing grid init)
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
        if(new_chunk3!=new_chunk4) {
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

        ref_elev_temp      = pGrid_tave->GetRefElevation(k);
        if(ref_elev_temp==RAV_BLANK_DATA) { ref_elev_temp = _pHydroUnits[k]->GetElevation(); } //disabling orographic effects if no elevation given (warning in Forcing grid init)
      }

      // ---------------------
      // (3) read gridded recharge
      // ---------------------
      if(recharge_gridded)
      {
        pGrid_recharge   = GetForcingGrid(F_RECHARGE);
        pGrid_recharge-> ReadData(Options,tt.model_time);
        F.recharge   = pGrid_recharge->GetWeightedValue(k,tt.model_time,Options.timestep);
      }
      // ---------------------
      // (4) read gridded precipitation temperature
      // ---------------------
      if(precip_temp_gridded)
      {
        pGrid_precip_temp   = GetForcingGrid(F_PRECIP_TEMP);
        pGrid_precip_temp-> ReadData(Options,tt.model_time);
        F.precip_temp   = pGrid_precip_temp->GetWeightedValue(k,tt.model_time,Options.timestep);
      }
      if(precip_conc_gridded)
      {
        pGrid_precip_conc   = GetForcingGrid(F_PRECIP_CONC);
        pGrid_precip_conc-> ReadData(Options,tt.model_time);
        F.precip_conc   = pGrid_precip_conc->GetWeightedValue(k,tt.model_time,Options.timestep);
      }
      // ---------------------
      // (5) read gridded PET, Others
      // ---------------------
      if(pet_gridded)
      {
        pGrid_pet   = GetForcingGrid(F_PET);
        pGrid_pet-> ReadData(Options,tt.model_time);
        F.PET   = pGrid_pet->GetWeightedValue(k,tt.model_time,Options.timestep);
      }
      if(owpet_gridded) {
        pGrid_owpet   = GetForcingGrid(F_OW_PET);
        pGrid_owpet-> ReadData(Options,tt.model_time);
        F.OW_PET   = pGrid_owpet->GetWeightedValue(k,tt.model_time,Options.timestep);
      }
      if (windspeed_gridded) {
        pGrid_windspeed = GetForcingGrid(F_WIND_VEL);
        pGrid_windspeed->ReadData(Options, tt.model_time);
        F.wind_vel = pGrid_windspeed->GetWeightedValue(k, tt.model_time, Options.timestep);
      }
      if (relhum_gridded) {
        pGrid_relhum = GetForcingGrid(F_REL_HUMIDITY);
        pGrid_relhum->ReadData(Options, tt.model_time);
        F.rel_humidity = pGrid_relhum->GetWeightedValue(k, tt.model_time, Options.timestep);
      }
      if (SWnet_gridded) {
        pGrid_SW_net = GetForcingGrid(F_SW_RADIA_NET);
        pGrid_SW_net->ReadData(Options, tt.model_time);
        F.SW_radia_net = pGrid_SW_net->GetWeightedValue(k, tt.model_time, Options.timestep);
      }
      if (LWinc_gridded) {
        pGrid_LW_inc = GetForcingGrid(F_LW_INCOMING);
        pGrid_LW_inc->ReadData(Options, tt.model_time);
        F.LW_incoming = pGrid_LW_inc->GetWeightedValue(k, tt.model_time, Options.timestep);
      }
      if (SW_gridded) {
        pGrid_SW = GetForcingGrid(F_SW_RADIA);
        pGrid_SW->ReadData(Options, tt.model_time);
        F.SW_radia = pGrid_SW->GetWeightedValue(k, tt.model_time, Options.timestep);
      }

      //-------------------------------------------------------------------
      //  Temperature Corrections
      //-------------------------------------------------------------------
      double tc;
      int p = _pHydroUnits[k]->GetSubBasinIndex();
      tc = _pSubBasins[p]->GetTemperatureCorrection();

      //--Gauge Corrections------------------------------------------------
      if (Options.in_bmi_mode && !rvt_file_provided)  // temperature was given by the BMI and no gauge corrections are to be applied
      {
        F.temp_daily_ave = F.temp_daily_max = F.temp_daily_min = F.temp_ave;  // TODO: check if this is acceptable
      }
      else if (!(temp_ave_gridded || (temp_daily_min_gridded && temp_daily_max_gridded) || temp_daily_ave_gridded)) //Gauge Data
      {
          double gauge_corr;
          F.temp_ave = F.temp_daily_ave = F.temp_daily_max = F.temp_daily_min = 0.0; // leave out monthly for now
          for (g = 0; g < _nGauges; g++)
          {
              gauge_corr = tc + _pGauges[g]->GetTemperatureCorr();
              wt = _aGaugeWtTemp[k][g];

              F.temp_ave       += wt * (gauge_corr + Fg[g].temp_ave);
              F.temp_daily_ave += wt * (gauge_corr + Fg[g].temp_daily_ave);
              F.temp_daily_max += wt * (gauge_corr + Fg[g].temp_daily_max);
              F.temp_daily_min += wt * (gauge_corr + Fg[g].temp_daily_min);

          }
      }
      else //Gridded Data
      {
          double grid_corr;
          grid_corr = pGrid_pre->GetTemperatureCorr();
          if ((tc+grid_corr)!=0.0){
            F.temp_ave       += tc + grid_corr;
            F.temp_daily_ave += tc + grid_corr;
            F.temp_daily_max += tc + grid_corr;
            F.temp_daily_min += tc + grid_corr;
          }
      }

      F.temp_ave_unc = F.temp_daily_ave;
      F.temp_min_unc = F.temp_daily_min;
      F.temp_max_unc = F.temp_daily_max;

      CorrectTemp(Options,F,elev,ref_elev_temp,tt);

      ApplyForcingPerturbation(F_TEMP_AVE, F, k, Options, tt);

      //-------------------------------------------------------------------
      //  Copy Daily values from current day, earlier time steps
      //    (done after temperature corrections so that the uncorrected values are used in calculating the lapse rate for temp_ave)
      //-------------------------------------------------------------------
      if(!tt.day_changed) {
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

      F.rel_humidity = EstimateRelativeHumidity(Options.rel_humidity,F,_pHydroUnits[k]);

      //-------------------------------------------------------------------
      // Snow fraction Calculations
      //-------------------------------------------------------------------
      F.snow_frac =EstimateSnowFraction(Options.rainsnow,_pHydroUnits[k],&F,Options);

      //-------------------------------------------------------------------
      //  Precip Corrections
      //-------------------------------------------------------------------
      double rc,sc;
      // int p=_pHydroUnits[k]->GetSubBasinIndex(); // defined for temp correction already
      rc=_pSubBasins[p]->GetRainCorrection();
      sc=_pSubBasins[p]->GetSnowCorrection();

      //--Gauge Corrections------------------------------------------------
      if(!(pre_gridded || snow_gridded || rain_gridded)) //Gauge or BMI-injected Data
      {
        double gauge_corr;
        F.precip=F.precip_5day=F.precip_daily_ave=0.0;
        if ((!Options.in_bmi_mode) || rvt_file_provided) {
          // Gauge-based precip and snowfall correction
          for(g=0; g<_nGauges; g++)
          {
            gauge_corr= F.snow_frac*sc*_pGauges[g]->GetSnowfallCorr() + (1.0-F.snow_frac)*rc*_pGauges[g]->GetRainfallCorr();
            wt=_aGaugeWtPrecip[k][g];
            F.precip         += wt*gauge_corr*Fg[g].precip;
            F.precip_daily_ave+=wt*gauge_corr*Fg[g].precip_daily_ave;
            F.precip_5day    += wt*gauge_corr*Fg[g].precip_5day;
          }
        }
		else {
          // Gauge-less precip and snowfall correction
          gauge_corr         = (F.snow_frac * sc) + ((1.0-F.snow_frac)*rc);
          F.precip           = gauge_corr * _pHydroUnits[k]->GetForcingFunctions()->precip;
          F.precip_daily_ave = gauge_corr * _pHydroUnits[k]->GetForcingFunctions()->precip;  // TODO: check if this is acceptable
          F.precip_5day      = gauge_corr * _pHydroUnits[k]->GetForcingFunctions()->precip * 5;   // TODO: check if this is acceptable
        }
      }
      else //Gridded Data
      {
        double grid_corr;
        double rain_corr=pGrid_pre->GetRainfallCorr();
        double snow_corr=pGrid_pre->GetSnowfallCorr();
        grid_corr= F.snow_frac*sc*snow_corr + (1.0-F.snow_frac)*rc*rain_corr;

        F.precip          *=grid_corr;
        F.precip_daily_ave*=grid_corr;
        F.precip_5day      =NETCDF_BLANK_VALUE;
      }

      //--Orographic corrections-------------------------------------------
      CorrectPrecip(Options,F,elev,ref_elev_precip,k,tt);

      ApplyForcingPerturbation(F_PRECIP  , F, k, Options, tt);
      ApplyForcingPerturbation(F_RAINFALL, F, k, Options, tt);
      ApplyForcingPerturbation(F_SNOWFALL, F, k, Options, tt);

      //-------------------------------------------------------------------
      //  Wind Velocity
      //-------------------------------------------------------------------

      F.wind_vel = EstimateWindVelocity(Options,_pHydroUnits[k],F,ref_measurement_ht,k);

      //ApplyForcingPerturbation(F_WIND_VEL, F, k, Options, tt);

      //-------------------------------------------------------------------
      //  Cloud Cover
      //-------------------------------------------------------------------

      F.cloud_cover = EstimateCloudCover(Options,F,k);

      //-------------------------------------------------------------------
      //  Radiation Calculations
      //-------------------------------------------------------------------

      F.SW_radia = CRadiation::EstimateShortwaveRadiation(this, &F, _pHydroUnits[k], tt, F.ET_radia, F.ET_radia_flat);
      F.SW_radia_unc = F.SW_radia;
      F.SW_radia *= CRadiation::SWCloudCoverCorrection(this, &F, elev);

      F.SW_radia_subcan = F.SW_radia * CRadiation::SWCanopyCorrection(this, _pHydroUnits[k]);

      if(Options.SW_radia_net == NETSWRAD_CALC) //(default)
      {
        F.SW_radia_net  = F.SW_radia       *(1-_pHydroUnits[k]->GetTotalAlbedo(false));
        F.SW_subcan_net = F.SW_radia_subcan*(1-_pHydroUnits[k]->GetTotalAlbedo(true ));
      }//otherwise, uses data

      F.LW_radia_net = CRadiation::EstimateLongwaveRadiation(GetStateVarIndex(SNOW), this, &F, _pHydroUnits[k], F.LW_incoming);

      //-------------------------------------------------------------------
      //  Potential Melt Rate
      //-------------------------------------------------------------------

      F.potential_melt=EstimatePotentialMelt(&F,Options.pot_melt,Options,_pHydroUnits[k],tt);

      //-------------------------------------------------------------------
      //  PET Calculations
      //-------------------------------------------------------------------
      // last but not least - needs all of the forcing params calculated above
      if(!pet_gridded) //Gauge Data
      {
        F.PET   =EstimatePET(F,_pHydroUnits[k],ref_measurement_ht,ref_elev_temp,Options.evaporation,Options,tt,false);
      }
      if (!owpet_gridded)
      {
        F.OW_PET=EstimatePET(F,_pHydroUnits[k],ref_measurement_ht,ref_elev_temp,Options.ow_evaporation,Options,tt,true);
      }
      CorrectPET(Options,F,_pHydroUnits[k],elev,ref_elev_temp,k);

      //-------------------------------------------------------------------
      // Irrigation
      //-------------------------------------------------------------------
      F.irrigation=0.0; //FOR NOW
      //F.irrigation=GetDemandOptimizer()->EstimateIrrigation(F,_pHydroUnits[k],Options,tt);

      //-------------------------------------------------------------------
      // Direct evaporation of rainfall / irrigation
      //-------------------------------------------------------------------
      if(Options.direct_evap)
      {
        double rainfall=F.precip*(1.0-F.snow_frac);
        double reduce  =min(F.PET,rainfall+F.irrigation);
        double rainfrac=rainfall/(rainfall+F.irrigation);
        if ((rainfall+F.irrigation)==0){rainfrac=1.0;}
        if(F.precip+F.irrigation-reduce>0.0) { F.snow_frac=1.0-(rainfall-reduce*rainfrac)/(F.precip-reduce*rainfrac); }
        F.precip    -=reduce*rainfrac;
        F.PET       -=reduce*rainfrac;
        F.irrigation-=reduce*(1.0-rainfrac);
      }

      ApplyLocalParamOverrrides(k,true);
    }//end if (!_pHydroUnits[k]->IsDisabled())

    //-------------------------------------------------------------------
    // Update
    //-------------------------------------------------------------------
    _pHydroUnits[k]->UpdateForcingFunctions(F);

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
                                const force_struct &F,
                                const CHydroUnit   *pHRU)
{
  double relhum=0.5;
  if (method==RELHUM_CONSTANT)
  {
    relhum=0.5;
  }
  else if (method==RELHUM_MINDEWPT)
  { //uses minimum daily temperature as proxy for dew point temperature
    double dew_point_temp=F.temp_daily_min;
    relhum=min(GetSaturatedVaporPressure(dew_point_temp)/GetSaturatedVaporPressure(F.temp_ave),1.0);
  }
  else if (method==RELHUM_DATA)
  {
    relhum=F.rel_humidity;
  }

  return min(relhum * pHRU->GetSurfaceProps()->relhum_corr,1.0);

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
//
double CModel::EstimateWindVelocity(const optStruct    &Options,
                                    const CHydroUnit   *pHRU,
                                    const force_struct &F,
                                    const double       &wind_measurement_ht,
                                    const int           k)
{
  double wind_vel=2.0; //[m/s]
  //---------------------------------------------------------------------
  if (Options.wind_velocity==WINDVEL_CONSTANT)
  {
    wind_vel=2.0;//[m/s] (global average)
  }
  //---------------------------------------------------------------------
  else if (Options.wind_velocity==WINDVEL_DATA)
  {
    wind_vel=F.wind_vel;
  }
  //---------------------------------------------------------------------
  else if(Options.wind_velocity==WINDVEL_UBC_MOD)
  { //simplifed version of the UBCWM algorithm
    double v_min=pHRU->GetSurfaceProps()->min_wind_speed;
    double v_max=pHRU->GetSurfaceProps()->max_wind_speed;
    double Tmin =F.temp_daily_min;
    double Tmax =F.temp_daily_max;

    double wt=min(0.04*(Tmax-Tmin),1.0);

    wind_vel=max(0.0,(1-wt)*v_min+(wt)*v_max);
  }
  //---------------------------------------------------------------------
  else if(Options.wind_velocity==WINDVEL_SQRT)
  {
    double b = this->_pGlobalParams->GetParams()->windvel_icept;
    double m = this->_pGlobalParams->GetParams()->windvel_scale;
    double Tmin =F.temp_daily_min;
    double Tmax =F.temp_daily_max;

    wind_vel=max(0.0,m*(Tmax-Tmin)+b);
  }
  //---------------------------------------------------------------------
  else if(Options.wind_velocity==WINDVEL_LOG)
  {
    double b = this->_pGlobalParams->GetParams()->windvel_icept;
    double m = this->_pGlobalParams->GetParams()->windvel_scale;
    double Tmin = F.temp_daily_min;
    double Tmax = F.temp_daily_max;

    wind_vel=max(0.0,exp(m*(Tmax-Tmin)+b)-1.0);
  }
  //---------------------------------------------------------------------
  else if (Options.wind_velocity==WINDVEL_UBCWM)
  {
    double TED, A1term;
    TED=max(F.temp_daily_max-F.temp_daily_min,0.0);

    const double rf_elev=2000;
    double elev = _pHydroUnits[k]->GetElevation();
    double Fc   = _pHydroUnits[k]->GetSurfaceProps()->forest_coverage;
    double P0TEDL = this->_pGlobalParams->GetParams()->UBC_lapse_params.P0TEDL;
    double P0TEDU = this->_pGlobalParams->GetParams()->UBC_lapse_params.P0TEDU;
    double A0term = this->_pGlobalParams->GetParams()->UBC_lapse_params.max_range_temp;
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

    wind_vel = (1-wt)*max_wind_speed +(wt)*1;

    upperswap(wind_vel,1.0);
    lowerswap(wind_vel,max_wind_speed - 1.0); //why the minus one? who knows...

    //elevation correction
    wind_vel*=max(sqrt(elev/1000.0),1.0);

    //forest correction
    const double F0WIND=0.7; //ratio of wind in forest vs wind in open
    wind_vel*=((Fc)*F0WIND+(1.0-Fc)*1.0);

    wind_vel*=M_PER_KM/SEC_PER_HR; //KPH_TO_M_PER_S;
  }

  return wind_vel * pHRU->GetSurfaceProps()->wind_vel_corr;
}

double CModel::WindspeedAtHeight (const double       &z,
                                  const optStruct    &Options,
                                  const CHydroUnit   *pHRU,
                                  const force_struct &F,
                                  const double       &z_ref)
{
  //---------------------------------------------------------------------
  if      (Options.wind_profile==WINDPROF_UNIFORM)
  {
    return F.wind_vel;//m/s
  }
  //---------------------------------------------------------------------
  else if (Options.wind_profile==WINDPROF_LOGARITHMIC)
  {
    double z_0=pHRU->GetVegVarProps()->roughness;     //[m]
    double zpd=pHRU->GetVegVarProps()->zero_pln_disp; //[m]

    return F.wind_vel*log((z-zpd)/z_0)/log((z_ref-zpd)/z_0);
  }
  //---------------------------------------------------------------------
  else if(Options.wind_profile==WINDPROF_MAHAT)
  {
    double LAI       =pHRU->GetVegVarProps()->LAI;
    double h_raw     =pHRU->GetVegVarProps()->height;
    double snow_depth=pHRU->GetSnowDepth();
    double z0ground  =pHRU->GetSurfaceProps()->roughness;
    double nd=1.0;
    double y=1;
    //y is an integer indicating one of the three basic forest profiles [e.g., Massman, 1982; Meyers et al., 1998]:
    //y = 1 for young pine, y=2 for leafed deciduous tree, and y=3 for old pine with long stems and clumping at the top.

    double z0snow = 0.01;         // roughness length of snow[m]
    double z0s = z0ground;        //surface roughness

    double h,zz,d,z0c,u_star,Ua,Uc,Ub,Uz;
    if(snow_depth>z0ground) {z0s=z0snow;}
    // Generate coordinates; correct heights relative to surface(wind profile,arbitrarily,goes to double veg height)
    h = h_raw - snow_depth;                                 // height of canopy above snow depth[m]
    zz = z - snow_depth;                                     // height above snow[m]
    d   = h*(0.05+0.5*pow(LAI,0.20)+0.050*(y-1));           // zero-plane displacement height[m] FLAG!the exponent is 0.02 in paper,0.2 in source code
    //d      = h*0.64                                       // zero-plane displacement height[m] simple
    z0c = h*(0.23-0.1*pow(LAI,0.25)-0.015*(y-1));           // roughness length for top of canopy
    // z0c+=z0snow;                                    // add nominal depth so it's never 0

    // if canopy taller than reference height (and not buried in snow), calculate canopy corrections
    if((h_raw > z_ref) && (h > snow_depth))
    {
      // Wind at different heights
      u_star = VON_KARMAN * F.wind_vel/log((h+zz-d)/z0c);                  // compute the friction velocity(m/s)
      Ua=u_star/VON_KARMAN*log((zz-d)/z0c);                                // above canopy(Logarithmic) // z >= h
      Uc=u_star/VON_KARMAN*log((h -d)/z0c)*exp(-nd*LAI*((h-zz    )/h));     // within canopy(exponential)
      Ub=u_star/VON_KARMAN*log((h -d)/z0c)*exp(-nd*LAI*((h-z_ref)/h))*log((zz-z0s)/z0s)/log(z_ref/z0s);    // assuming zero-plane displacement of ground is 0...

      // Assemble them(logarithmic when z > h; exponential when z <= h & z > z_ms' logarithmic as z_ref -> 0)
      if       (zz>=h)               {Uz=Ua;}
      else if ((zz< h) && (zz>z_ref)){Uz=Uc;}
      else                           {Uz=Ub;}
    }
    else {
      // Otherwise,single logarithmic profile
      Uz = F.wind_vel*log((zz-z0s)/z0s)/log((10.0+h)/z0s);    // 10 m measurement height assumign zero-plane displacement of ground is 0...
    }
    return max(0.0,Uz);
  }
  return 0.0;
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


//////////////////////////////////////////////////////////////////
/// \brief Estimates fraction of precipitation that is snow based temperature over timestep
/// \ref Brook90 approximation is really only valid for daily time step \cite Dingman1994
///
/// \param method [in] Method of converting total precipitation to rain/snow
/// \param pHRU [in] HRU unit
/// \param *F [in & out] Forcing functions for a specific HRU
/// \param &Options [in] Global model options information
/// \return Double value for fraction of precipitation which is snow
//
double CModel::EstimateSnowFraction(const rainsnow_method method,
                                    const CHydroUnit* pHRU,
                                    const force_struct* F,
                                    const optStruct& Options)
{
	//-----------------------------------------------------------
	if (method == RAINSNOW_DATA)
	{
	    return F->snow_frac;
	}
	//-----------------------------------------------------------
	else if (method == RAINSNOW_DINGMAN)
	{ //from Brook90 model, Dingman pg 109
	    double temp = this->_pGlobalParams->GetParams()->rainsnow_temp;
	    if (F->temp_daily_max <= temp) { return 1.0; }
	    if (F->temp_daily_min >= temp) { return 0.0; }
	    return (temp - F->temp_daily_min) / (F->temp_daily_max - F->temp_daily_min);
	}
	//-----------------------------------------------------------
	else if (method == RAINSNOW_THRESHOLD)
	{ //abrupt threshold temperature (e.g., HYMOD)
	  double temp = this->_pGlobalParams->GetParams()->rainsnow_temp;
	  if (F->temp_ave <= temp) { return 1.0; }
	  else                     { return 0.0; }
	}
  //-----------------------------------------------------------
  else if ((method == RAINSNOW_HBV) || (method == RAINSNOW_UBCWM))
  {//linear variation based upon daily average temperature
      double frac;
      double delta = this->_pGlobalParams->GetParams()->rainsnow_delta;
      double temp  = this->_pGlobalParams->GetParams()->rainsnow_temp;

      if      (F->temp_daily_ave <= (temp - 0.5 * delta)) { frac = 1.0; }
      else if (F->temp_daily_ave >= (temp + 0.5 * delta)) { frac = 0.0; }//assumes only daily avg. temp is included
      else {
          frac = 0.5 + (temp - F->temp_daily_ave) / delta;
      }
      if    (method == RAINSNOW_UBCWM) { return frac; }
      //HBV-EC implementation - correction only applied to rain portion of snow (assumes snow data provided)
      else if (method == RAINSNOW_HBV) { return frac * (1.0 - F->snow_frac) + F->snow_frac; }
  }
  //-----------------------------------------------------------
  else if (method == RAINSNOW_HSPF) // Also, from HydroComp (1969)
  {
      double temp = this->_pGlobalParams->GetParams()->rainsnow_temp;
      double snowtemp;
      double dewpt = GetDewPointTemp(F->temp_ave, F->rel_humidity);

      snowtemp = temp + (F->temp_ave - dewpt) * (0.5808 + 0.0144 * F->temp_ave);

      if (F->temp_ave < snowtemp) { return 1.0; }
      else                        { return 0.0; }
  }
  //-----------------------------------------------------------
  else if (method == RAINSNOW_HARDER) // From Harder & Pomeroy 2013, Ported from CRHM (Pomeroy, 2007)
  {
    // \todo[funct] replace with T_icebulb=GetIcebulbTemp(T,rel_hum);
    double Tk, D, lambda, pta, L, snowfrac, T_icebulb;
    double rel_hum;
    rel_hum = F->rel_humidity;
    Tk = F->temp_ave + ZERO_CELSIUS;

    D = 0.0000206 * pow(Tk / ZERO_CELSIUS, 1.75);
    lambda = 0.000063 * Tk + 0.00673;
    pta = 18.01528 * ((rel_hum) * 0.611 * exp((17.3 * F->temp_ave) / (237.3 + F->temp_ave))) / (0.00831441 * Tk) / 1000.0;

    if (F->temp_ave > FREEZING_TEMP) { L = 1000.0 * (2501.0 - 2.361 * F->temp_ave); }
    else { L = 1000.0 * (2834.1 - 0.29 * F->temp_ave - 0.004 * F->temp_ave * F->temp_ave); }

    double crit = ALMOST_INF;
    double T_old(250.0), T_guess(250.0);
    double T1, T2;
    double dT = 0.002; //deg C
    double a, b, c;
    while (crit > 0.0001)
    { //ice bulb temp found using the newton-raphson method
        dT = 0.002 * T_old;
        T1 = T_old + dT / 2.0;
        T2 = T_old - dT / 2.0;
        a = (-T_old + Tk + (L * D / lambda) * (pta - (18.01528 * (0.611 * exp((17.3 * (T_old - ZERO_CELSIUS)) / (237.3 + (T_old - ZERO_CELSIUS)))) / (0.00831441 * T_old) / 1000)));
        b = (-T1    + Tk + (L * D / lambda) * (pta - (18.01528 * (0.611 * exp((17.3 * (T1    - ZERO_CELSIUS)) / (237.3 + (T1    - ZERO_CELSIUS)))) / (0.00831441 * T1   ) / 1000)));
        c = (-T2    + Tk + (L * D / lambda) * (pta - (18.01528 * (0.611 * exp((17.3 * (T2    - ZERO_CELSIUS)) / (237.3 + (T2    - ZERO_CELSIUS)))) / (0.00831441 * T2   ) / 1000)));

        T_guess = T_old - (a / ((b - c) / dT));
        crit = fabs(T_old - T_guess);
        T_old = T_guess;
    } // end while

    T_icebulb = T_guess - ZERO_CELSIUS;

    snowfrac = 1.0;
    if (T_icebulb > -10.0) {
        snowfrac = 1.0 - 1.0 / (1.0 + 2.50286 * pow(0.125, T_icebulb));
    }
    return snowfrac;
  }
  else if (method == RAINSNOW_WANG)
  {
    //from Wang et al., 2019. A wet-bulb temperature-based rain-snow partitioning scheme improves snowpack prediction over the drier western United States. Geophysical Research Letters, 46(23), pp.13825-13835.

    double Ta     =F->temp_ave;
    double rel_hum=F->rel_humidity;
    double P      =F->air_pres;

    double Twet=GetWetBulbTemperature(P,rel_hum,Ta);
    return  1.0/(1.0+6.99e-5*exp(2.0*(Twet+3.97)));
  }
  else if (method == RAINSNOW_SNTHERM89)
  {
    //from Jordan, R.: A one-dimensional temperature model for a snow cover: Technical documentation for SNTHERM.89., 1991. (as used in Noah-MP-3.6)
    double Ta=F->temp_ave;

    if      (Ta>2.5){return 0.0;}
    else if (Ta>2.0){return 0.6;}
    else if (Ta>0.5){return 1.0-(0.4/1.5)*(Ta-0.5); }
    else            {return 1.0;}
  }
  return 0.0;
  //Should check/add:
  //Stefan W. Kienzle,A new temperature based method to separate rain and snow,
  ///< Hydrological Processes 22(26),p5067-5085,2008,http://dx.doi.org/10.1002/hyp.7131 \cite kienzle2008HP
}
