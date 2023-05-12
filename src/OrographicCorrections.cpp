/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2020 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Model.h"

double UBCPreciptiationByElev(const double temp,
                              const double gauge_elev,
                              const double band_elev,
                              const UBC_lapse PP);

//////////////////////////////////////////////////////////////////
/// \brief Corrects gauge temperature for elevation effects
/// \remark UBCWM Orographic Corrections adapted from UBC Watershed model
///  code,
/// \copyright (c) Michael Quick
///
/// \param &Options [in] Global model options information
/// \param &F [out] Forcing functions for HRU
/// \param elev [in] elevation of this HRU
/// \param ref_elev [in] Reference temperature elevation (usually met station elevation)
/// \param &tt [in] current time strucure
//
void   CModel::CorrectTemp(const optStruct   &Options,
                           force_struct      &F,
                           const double       elev,
                           const double       ref_elev,
                           const time_struct &tt)
{
  
  //---------------------------------------------------------------------------
  if ((Options.orocorr_temp==OROCORR_SIMPLELAPSE) ||
      (Options.orocorr_temp==OROCORR_HBV        ))
  {
    double lapse=CGlobalParams::GetParams()->adiabatic_lapse;//[C/km]
    lapse/=1000.0;//convert to C/m
    F.temp_ave-=lapse*(elev-ref_elev);

    if(tt.day_changed)
    {
      F.temp_daily_ave-=lapse*(elev-ref_elev);
      F.temp_daily_min-=lapse*(elev-ref_elev);
      F.temp_daily_max-=lapse*(elev-ref_elev);
      F.temp_month_max-=lapse*(elev-ref_elev);
      F.temp_month_min-=lapse*(elev-ref_elev);
      if (Options.orocorr_temp!=OROCORR_HBV){
        F.temp_month_ave-=lapse*(elev-ref_elev); //not corrected in HBV
      }
    }
  }
  //---------------------------------------------------------------------------
  else if (Options.orocorr_temp==OROCORR_UBCWM)
  {

    double V=0.0, wt;
    double adiabatic;
    double lo_lapse_max,hi_lapse_max,lo_lapse_min,hi_lapse_min,lo_lapse_ave,hi_lapse_ave,hi_lapse_sub,lo_lapse_sub;

    //calculate temperature lapse rates
    //--------------------------------------------------------------------
    const   global_struct *globals=CGlobalParams::GetParams();
    UBC_lapse lapse_params=globals->UBC_lapse_params;

    if (lapse_params.A0PPTP > 0){
      V = min(F.precip_daily_ave/ lapse_params.A0PPTP,1.0);
    }
    adiabatic=(V)*globals->wet_adiabatic_lapse+(1.0-V)* globals->adiabatic_lapse;
    wt=((F.temp_daily_max-F.temp_daily_min) / lapse_params.max_range_temp);//A0TERM

    //when temperature range is small (-->0), lapse approaches adiabatic/wet adiabatic
    lo_lapse_max =(wt)*lapse_params.A0TLXM+(1-wt)*adiabatic;
    hi_lapse_max =(wt)*lapse_params.A0TLXH+(1-wt)*adiabatic;
    lo_lapse_min =(wt)*lapse_params.A0TLNM+(1-wt)*adiabatic;
    hi_lapse_min =(wt)*lapse_params.A0TLNH+(1-wt)*adiabatic;

    lo_lapse_ave = (lo_lapse_max + lo_lapse_min)/2.0;
    hi_lapse_ave = (hi_lapse_max + hi_lapse_min)/2.0;
    hi_lapse_sub = hi_lapse_min+(F.temp_ave-F.temp_daily_min)*(hi_lapse_max-hi_lapse_min)/(F.temp_daily_max-F.temp_daily_min);
    lo_lapse_sub = lo_lapse_min+(F.temp_ave-F.temp_daily_min)*(lo_lapse_max-lo_lapse_min)/(F.temp_daily_max-F.temp_daily_min);

    if(tt.day_changed)
    {
      //calculate elevation-adjusted temperatures
      //--------------------------------------------------------------------
      double max_temp2K, min_temp2K,ave_temp2K;
      //lapse (station to 2K)
      if (ref_elev >= 2000) {
        max_temp2K=F.temp_daily_max-hi_lapse_max * (2000 - ref_elev) / 1000;
        min_temp2K=F.temp_daily_min-hi_lapse_min * (2000 - ref_elev) / 1000;
        ave_temp2K=F.temp_daily_ave-hi_lapse_ave * (2000 - ref_elev) / 1000;
      }
      else{
        max_temp2K=F.temp_daily_max-lo_lapse_max * (2000 - ref_elev) / 1000;
        min_temp2K=F.temp_daily_min-lo_lapse_min * (2000 - ref_elev) / 1000;
        ave_temp2K=F.temp_daily_ave-lo_lapse_ave * (2000 - ref_elev) / 1000;
      }
      //lapse (2K to HRU elevation)
      if (elev >= 2000 ){
        F.temp_daily_max = max_temp2K - hi_lapse_max * (elev - 2000) / 1000;
        F.temp_daily_min = min_temp2K - hi_lapse_min * (elev - 2000) / 1000;
        F.temp_daily_ave = ave_temp2K - hi_lapse_ave * (elev - 2000) / 1000;
      }
      else{
        F.temp_daily_max = max_temp2K - lo_lapse_max * (elev - 2000) / 1000;
        F.temp_daily_min = min_temp2K - lo_lapse_min * (elev - 2000) / 1000;
        F.temp_daily_ave = ave_temp2K - lo_lapse_ave * (elev - 2000) / 1000;
      }
      upperswap(F.temp_daily_max, F.temp_daily_min);
    }
    if (Options.timestep>=1.0){
      F.temp_daily_ave = 0.5*(F.temp_daily_max + F.temp_daily_min);
      F.temp_ave       =F.temp_daily_ave;
    }
    else{//sub-daily timestep
      double temp2K;
      if (ref_elev >= 2000) {
        temp2K=F.temp_ave-hi_lapse_sub * (2000 - ref_elev) / 1000;
      }
      else{
        temp2K=F.temp_ave-lo_lapse_sub * (2000 - ref_elev) / 1000;
      }
      if (elev >= 2000 ){
        F.temp_ave=temp2K-hi_lapse_sub * (elev - 2000) / 1000;
      }
      else{
        F.temp_ave=temp2K-lo_lapse_sub * (elev - 2000) / 1000;
      }
    }

  }
  //---------------------------------------------------------------------------
  else if (Options.orocorr_temp==OROCORR_UBCWM2)
  {
    double elev1=_pGauges[0]->GetElevation();
    double elev2=_pGauges[1]->GetElevation();
    double t=tt.model_time;

    if(tt.day_changed)
    {
      double XLAPSXA,XLAPSXB;
      double XLAPSNA,XLAPSNB;
      double gaugemax1=_pGauges[0]->GetForcingValue(F_TEMP_DAILY_MAX,t,Options.timestep);
      double gaugemax2=_pGauges[1]->GetForcingValue(F_TEMP_DAILY_MAX,t,Options.timestep);
      double gaugemin1=_pGauges[0]->GetForcingValue(F_TEMP_DAILY_MIN,t,Options.timestep);
      double gaugemin2=_pGauges[1]->GetForcingValue(F_TEMP_DAILY_MIN,t,Options.timestep);
      XLAPSXA = (gaugemax2 - gaugemax1) / (elev2 - elev1);
      XLAPSXB = gaugemax1-XLAPSXA*elev1;
      XLAPSNA = (gaugemin2 - gaugemin1) / (elev2 - elev1);
      XLAPSNB = gaugemin1-XLAPSNA*elev1;

      //Lapse rate calculated from first two stations
      F.temp_daily_max = XLAPSXB + XLAPSXA * elev;
      F.temp_daily_min = XLAPSNB + XLAPSNA * elev;

      upperswap(F.temp_daily_max,F.temp_daily_min);
      F.temp_daily_ave=0.5*(F.temp_daily_max+F.temp_daily_min);
    }

    if (Options.timestep>=1.0){
      F.temp_ave=F.temp_daily_ave;
    }
    else{//sub-daily timestep
      double XLAPSTA,XLAPSTB;
      double gaugetemp1=_pGauges[0]->GetForcingValue(F_TEMP_AVE,t,Options.timestep);
      double gaugetemp2=_pGauges[1]->GetForcingValue(F_TEMP_AVE,t,Options.timestep);
      XLAPSTA = (gaugetemp2 - gaugetemp1) / (elev2 - elev1);
      XLAPSTB = gaugetemp1-XLAPSTA*elev1;
      F.temp_ave = XLAPSTB + XLAPSTA * elev;
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects gauge precipitation for elevation effects
/// \remark UBCWM Orographic Corrections adapted from UBC Watershed model
///  code,
/// \copyright (c) Michael Quick
///
/// \param &Options [in] Global model options information
/// \param &F [out] Forcing functions in HRU k
/// \param elev [in] HRU elevation
/// \param ref_elev [in] Reference elevation (usually met station elevation)
/// \param k [in] index of HRU in model
/// \param &tt [in] current time structure
//
void CModel::CorrectPrecip(const optStruct     &Options,
                           force_struct  &F,
                           const double         elev,
                           const double         ref_elev ,
                           const int            k,
                           const time_struct    &tt) const
{


  //---------------------------------------------------------------------------
  if (Options.orocorr_precip==OROCORR_SIMPLELAPSE)
  {
    double lapse=CGlobalParams::GetParams()->precip_lapse;
    lapse/=1000; //[mm/d/km]->[mm/d/m]
    if (F.precip > REAL_SMALL){
      F.precip           = max(F.precip           + lapse*(elev - ref_elev), 0.0);
      F.precip_5day      = max(F.precip_5day      + lapse*(elev - ref_elev), 0.0);
      F.precip_daily_ave = max(F.precip_daily_ave + lapse*(elev - ref_elev), 0.0);
    }
  }
  //---------------------------------------------------------------------------
  else if (Options.orocorr_precip==OROCORR_HBV)
  {
    double corr_upper=CGlobalParams::GetParams()->HBVEC_lapse_upper/M_PER_KM;
    double corr      =CGlobalParams::GetParams()->HBVEC_lapse_rate/M_PER_KM;
    double lapse_elev=CGlobalParams::GetParams()->HBVEC_lapse_elev;
    double add=0.0;
    if (elev>lapse_elev){
      add=(corr_upper-corr)*(elev-lapse_elev);
    }
    F.precip          *=max(1.0+corr*(elev-ref_elev)+add,0.0);
    F.precip_5day     *=max(1.0+corr*(elev-ref_elev)+add,0.0);
    F.precip_daily_ave*=max(1.0+corr*(elev-ref_elev)+add,0.0);
  }
  //---------------------------------------------------------------------------
  else if ((Options.orocorr_precip==OROCORR_UBCWM) || (Options.orocorr_precip==OROCORR_UBCWM2))
  {
    if(tt.day_changed)
    {
      //Current cheat on temp (assumes first HRU is in band 1)
      double band1_temp=_pHydroUnits[0]->GetForcing(F_TEMP_DAILY_AVE);
      //if(!Options.keepUBCWMbugs) {
      //  double band1_temp=F.temp_daily_ave; //More desirable?
      // }
      /// \todo [bug] temp should not be the temp at band 1!! (this is really unacceptable, but part of UBCWM)

      //orographic corrections
      _pHydroUnits[k]->SetPrecipMultiplier(UBCPreciptiationByElev(band1_temp,ref_elev,elev,CGlobalParams::GetParams()->UBC_lapse_params));
    }
    F.precip           *= _pHydroUnits[k]->GetPrecipMultiplier();
    F.precip_daily_ave *= _pHydroUnits[k]->GetPrecipMultiplier();

    //gauge and annual corrections
    //F.precip*=_pHydroUnits[k]->GetSurfaceParams()->precip_adj;//(1.0+P0PADJ)

    // \todo [re-org] should move to interception -needed because UBCWM doesn't explicitly model canopy
    double P0PINX=_pHydroUnits[k]->GetVegetationProps()->max_intercept_rate;
    double P0PINT=_pHydroUnits[k]->GetVegetationProps()->rain_icept_pct;
    double icept_factor=_pHydroUnits[k]->GetSurfaceProps()->UBC_icept_factor;

    F.precip          -=P0PINT*min(F.precip,P0PINX)*icept_factor;
    F.precip_daily_ave-=P0PINT*min(F.precip,P0PINX)*icept_factor;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates precipitation at given elevation with additional temperature correction
/// \remark UBCWM Orographic Corrections adapted from UBC Watershed model
///  code,
/// \note Utility method called by UBCPrecipitationByElev
/// \copyright (c) Michael Quick
///
/// \param PP [in] UBC lapse parameters
/// \param grad [in] gradient of lapse
/// \param delta_elev [in] Change in elevation [m]
/// \param ttband [in] gauge temperature
/// \return Precipitation at a given electation with additional temperature correction
//
double UBCOroAdjust(const UBC_lapse PP,const double grad, const double delta_elev, const double ttband)
{
  double F, adj_factor;
  double temp=max(ttband,0.0);
  F=max(1.0-PP.A0STAB*temp*temp,0.0);
  adj_factor=pow((1.0+0.01*grad*F),0.01*delta_elev);
  return adj_factor;
}

//////////////////////////////////////////////////////////////////
/// \brief Return a correction factor for precipitation at a specific elevation
///
/// \param temp [in] Temperature at gauge
/// \param gauge_elev [in] Gauge elevation
/// \param band_elev [in] Band elevation
/// \param PP [in] UBC lapse parameters
/// \return Correction factor for precipitation at specified elevation
//
double UBCPreciptiationByElev(const double temp,
                              const double gauge_elev,
                              const double band_elev,
                              const UBC_lapse PP)
{
  double p_hi,p_mid,p_low;

  if     (gauge_elev<=PP.E0LMID){
    p_low=                      UBCOroAdjust(PP,PP.P0GRADL,PP.E0LLOW-gauge_elev ,temp);
    p_mid=p_low*UBCOroAdjust(PP,PP.P0GRADL,PP.E0LMID-PP.E0LLOW  ,temp);
    p_hi =p_mid*UBCOroAdjust(PP,PP.P0GRADM,PP.E0LHI -PP.E0LMID  ,temp);
  }
  else if(gauge_elev<=PP.E0LHI){
    p_hi =                      UBCOroAdjust(PP,PP.P0GRADM,PP.E0LHI -gauge_elev ,temp);
    p_mid=                      UBCOroAdjust(PP,PP.P0GRADM,PP.E0LMID-gauge_elev ,temp);
    p_low=p_mid*UBCOroAdjust(PP,PP.P0GRADL,PP.E0LLOW-PP.E0LMID  ,temp);
  }
  else{
    p_hi =                      UBCOroAdjust(PP,PP.P0GRADU,PP.E0LHI -gauge_elev ,temp);
    p_mid=p_hi *UBCOroAdjust(PP,PP.P0GRADM,PP.E0LMID-PP.E0LHI   ,temp);
    p_low=p_mid*UBCOroAdjust(PP,PP.P0GRADL,PP.E0LLOW-PP.E0LMID  ,temp);
  }
  if     (band_elev<=PP.E0LMID){
    return p_low*UBCOroAdjust(PP,PP.P0GRADL,band_elev-PP.E0LLOW ,temp);}
  else if(band_elev<=PP.E0LHI){
    return p_mid*UBCOroAdjust(PP,PP.P0GRADM,band_elev-PP.E0LMID ,temp);}
  else{
    return p_hi *UBCOroAdjust(PP,PP.P0GRADU,band_elev-PP.E0LHI  ,temp);}
}

////////////////////////////////////////////////////////////////
/// \brief Corrects gauge PET for elevation effects
/// \remark UBCWM Orographic Corrections adapted from UBC Watershed model code,
/// \copyright (c) Michael Quick
///
/// \param &Options [in] Global model options information
/// \param &F [out] Forcing functions in HRU k
/// \param elev [in] HRU elevation
/// \param ref_elev [in] Reference elevation (usually met station elevation)
/// \param k [in] index of HRU
//
void CModel::CorrectPET(const optStruct &Options,
                        force_struct &F,
                        const CHydroUnit *pHRU,
                        const double elev,
                        const double ref_elev_temp,
                        const int k)
{
  double factor;
  //---------------------------------------------------------------------------
  if (Options.orocorr_PET==OROCORR_HBV)
  {
    factor=GLOBAL_PET_CORR*max(1.0-HBV_PET_ELEV_CORR*(elev-ref_elev_temp),0.0);
    F.PET   *=factor;
    F.OW_PET*=factor;
  }
  //---------------------------------------------------------------------------
  else if (Options.orocorr_PET==OROCORR_PRMS)
  {
    double sat_vap_max,sat_vap_min,c1,ch;
    double max_month_temp(0.0),min_month_temp(0.0);
    for (int g=0;g<_nGauges;g++)
    {
      double max=-ALMOST_INF;
      double min=ALMOST_INF;
      double tmp;
      for(int i=0;i<12;i++) {//get min/max annual temp
        tmp=_pGauges[g]->GetMonthlyAveTemp(i);
        upperswap(max,tmp);
        lowerswap(min,tmp);
      }
      max_month_temp+=_aGaugeWtTemp[k][g]*max;
      min_month_temp+=_aGaugeWtTemp[k][g]*min;
    }

    sat_vap_max=GetSaturatedVaporPressure(max_month_temp);
    sat_vap_min=GetSaturatedVaporPressure(min_month_temp);
    c1=68.0-(3.6*(FEET_PER_METER*(elev-ref_elev_temp))/1000);
    ch=50/(sat_vap_max-sat_vap_min)*MB_PER_KPA;

    factor=1.0/(c1+(13.0*ch));

    F.PET   *=factor;
    F.OW_PET*=factor;
  }
  //---------------------------------------------------------------------------
  else if (Options.orocorr_PET==OROCORR_UBCWM)
  {
    //NOTE - should not use with orocorr_precip==OROCORR_UBCWM - duplicates correction
    //RFS Correction
    /*const double A0PELA=0.9;
      double PET_corr=pHRU->GetSurfaceProps()->forest_PET_corr;
      double Fc      =pHRU->GetSurfaceProps()->forest_coverage;
      double forest_corr=((PET_corr)*Fc + 1.0*(1.0-Fc));
      double xevap2=max(A0PELA*0.001*(ref_elev-ref_elev_temp),0.0);
      //F.PET+=xevap2*forest_corr; //RFS correction - handled in UBCWM estimate PET command instead
      */
  }

  //sub-daily correction
  F.PET*=F.subdaily_corr;

  //Snow cover
  if(Options.snow_suppressPET) 
  {
    double SWE     =0.0;
    double snow_cov=1.0;
    int iSnow=GetStateVarIndex(SNOW);
    int iSC  =GetStateVarIndex(SNOW_COVER);
    if(iSnow!=DOESNT_EXIST) { SWE     =pHRU->GetStateVarValue(iSnow); }
    if(iSC  !=DOESNT_EXIST) { snow_cov=pHRU->GetStateVarValue(iSC); }
      
    if(SWE>0.1) { F.PET*=(1.0-snow_cov); }
  }
}
