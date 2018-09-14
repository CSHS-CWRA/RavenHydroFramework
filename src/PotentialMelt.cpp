/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Model.h"
double UBC_DailyPotentialMelt( const optStruct &Options,
                               const force_struct &F,
                               const CHydroUnit *pHRU,
                               const time_struct &tt);
//////////////////////////////////////////////////////////////////
/// \brief Estimates potential melt rate [mm/d]
/// \details can be positive or negative
///
/// \param *F [in] Forcing functions for a particular HRU over current time step
/// \param Options [in] global options structure
/// \param *pHRU [in] pointer to HRU for which radiation is calculated
/// \param tt [in] current time structure
/// \return double potential melt rate [mm/d]
//
double CModel::EstimatePotentialMelt(const force_struct *F,
                                     const optStruct &Options,
                                     const CHydroUnit *pHRU,
                                     const time_struct &tt)
{
  double albedo;
  albedo=pHRU->GetSnowAlbedo();

  //----------------------------------------------------------
  if (Options.pot_melt==POTMELT_DATA)
  {
    return F->potential_melt; //already read into gauge time series & interpolated
  }
  //----------------------------------------------------------
  else if (Options.pot_melt==POTMELT_DEGREE_DAY)
  {
    double Ma=pHRU->GetSurfaceProps()->melt_factor;
    return Ma*(F->temp_daily_ave-FREEZING_TEMP);
  }
  //----------------------------------------------------------
  else if (Options.pot_melt==POTMELT_HBV)
  {
    const surface_struct *surf=pHRU->GetSurfaceProps();

    double Ma    =surf->melt_factor;
    double Ma_min=surf->min_melt_factor;
    double AM    =surf->HBV_melt_asp_corr;
    double MRF   =surf->HBV_melt_for_corr;
    double Fc    =surf->forest_coverage;
    double slope_corr;

    //annual variation
    Ma=Ma_min+(Ma-Ma_min)*0.5*(1.0-cos(F->day_angle-WINTER_SOLSTICE_ANG));

    //Below code to replicate bug in HBV-EC
    //double OCT_1=(365.0-31-30-31-(136))/365.0*2.0*PI;//no idea what the 136 correction is, but it is needed
    //Ma=Ma_min+(Ma-Ma_min)*0.5*(1.0-cos(2.0*(F->day_angle-OCT_1)-WINTER_SOLSTICE_ANG));//TMP DEBUG: Northern hemisphere only (2 is temporary to match HBV!)

    //forest correction
    Ma*=((1.0-Fc)*1.0+(Fc)*MRF);

    //slope correction (forest only, for some reason)
    slope_corr=((1.0-Fc)*1.0+(Fc)*sin(pHRU->GetSlope()));

    //aspect corrections
    Ma*=max(1.0-AM*slope_corr*cos(pHRU->GetAspect()),0.0);//north facing slopes (aspect=0) have lower melt rates

    return Ma*(F->temp_daily_ave-FREEZING_TEMP);
  }
  //----------------------------------------------------------
  else if (Options.pot_melt==POTMELT_EB)
  {
    double rainfall  = F->precip*(1-F->snow_frac);
    double air_temp  = F->temp_ave;
    double air_pres  = F->air_pres;
    double rel_humid = F->rel_humidity;
    double wind_vel  = F->wind_vel;

    double ref_ht    = 2.0; // [m];
    double roughness = pHRU->GetVegVarProps()->roughness ;
    double surf_temp = pHRU->GetSnowTemperature();

    double K = (1-albedo)* F->SW_radia;   //net short wave radiation to snowpack [MJ/m2/d]
    double L = F->LW_radia;       //net long wave radiation [MJ/m2/d]
    double H = GetSensibleHeatSnow(air_temp,surf_temp,wind_vel,ref_ht,roughness); //[MJ/m2/d]
    double LE= GetLatentHeatSnow  (air_pres,air_temp,surf_temp,rel_humid,wind_vel,ref_ht,roughness); //[MJ/m2/d]
    double R = GetRainHeatInput   (surf_temp,rainfall ,rel_humid); //[MJ/m2/d]

    double S = K + L + H + LE + R;                //total heat input rate [MJ/m2/d]
    S*=(MM_PER_METER/DENSITY_WATER/LH_FUSION);    //[MJ/m2/d-->mm/d]
    return S;
  }
  //----------------------------------------------------------
  else if (Options.pot_melt==POTMELT_RESTRICTED)
  {
    double rad,convert;
    double Ma=pHRU->GetSurfaceProps()->melt_factor;
    double tmp_rate=threshPositive(Ma*(F->temp_daily_ave-FREEZING_TEMP));

    rad     = ((1-albedo)*F->SW_radia+F->LW_radia);//[MJ/m2/d]
    convert = MM_PER_METER/DENSITY_WATER/LH_FUSION;//for converting radiation to mm/d [mm/m]/[kg/m3]/[MJ/kg] = [mm-m2/MJ]

    return tmp_rate+ threshPositive(rad*convert);
  }
  //----------------------------------------------------------
  else if (Options.pot_melt==POTMELT_DD_RAIN)
  {
    //degree-day with rain-on-snow correction
    double Ma=pHRU->GetSurfaceProps()->melt_factor;
    double rainfall=F->precip*(1-F->snow_frac);
    double adv_melt= SPH_WATER/LH_FUSION*max(F->temp_ave,0.0)*rainfall; //[mm/d]
    return Ma*(F->temp_daily_ave-FREEZING_TEMP)+adv_melt;
  }
  //----------------------------------------------------------
  else if (Options.pot_melt==POTMELT_UBCWM)
  {
    return UBC_DailyPotentialMelt(Options,*F,pHRU,tt);
  }
  //----------------------------------------------------------
  else if (Options.pot_melt == POTMELT_USACE)
  {
    // Potential Melt
    double Ma;

    // parameters for snowmelt during rain
    double k1 = pHRU->GetSurfaceProps()->wind_exposure;
    double Pr = pHRU->GetForcingFunctions()->precip*(1 - pHRU->GetForcingFunctions()->snow_frac);
    double Ta = pHRU->GetForcingFunctions()->temp_daily_ave;
    double v  = pHRU->GetForcingFunctions()->wind_vel * SEC_PER_HR/M_PER_KM; // [m/s->km/hr]
    double N  = pHRU->GetForcingFunctions()->cloud_cover; // Estimated cloud cover
    double TdP = pHRU->GetForcingFunctions()->temp_daily_min; // difference between dew point temp measured at 3 m and snow surface

    // parameters for snowmelt without rain
    double kP = 1.1; // basin shortwave radiation melt factor
    double SF = pHRU->GetSurfaceProps()->forest_coverage;
    double I  = F->SW_radia / kP;
    double a  = pHRU->GetSnowAlbedo(); // Average snow surface albedo
    double TaP = Ta; // Difference between air temp measured at 3 m and snow surface
    double k2  = k1; // Convection-condensation melt factor

    // Cloud Base Temperature
    double TcP; // Difference between the cloud base temp and snow surface
    double cloud_base = (TaP - TdP) * 400 / FEET_PER_METER; // Estimation of cloud base height in M
    double lapse = CGlobalParams::GetParams()->adiabatic_lapse;//[C/km]
    lapse = lapse / 1000.0; //[C/m]
    TcP = TaP - cloud_base * lapse;

    // rain on snow
    if (Pr > 0) {
      //Heavily Forested
      if (SF >= 0.8) {Ma = (3.38 + 0.0126*Pr)*Ta + 1.3;}
      else           {Ma = (1.33 + 0.239*k1*v + 0.0126*Pr)*Ta + 2.3;}
    }
    // no rain
    else {

      if      (SF >= 0.8){ Ma = 3.38*(0.53*TaP + 0.47*TdP);}      // Heavily Forested (> 80 %)
      else if (SF >= 0.6){ Ma = k2*(0.239*v)*(0.22*TaP + 0.78*TdP) + SF*1.33*TaP;}// Forested (60 - 80 %)
      else if (SF >= 0.1){ Ma = kP*(1 - SF)*(2.43*I)*(1 - a) + k2*(0.239*v)*(0.22*TaP + 0.78*TdP) + SF*(1.33*TaP);}// Partly Forested (10 - 60 %)
      else               { Ma = kP*(1 - SF)*(3.08*I)*(1 - a) + (1 - N)*(0.969*TaP - 21.34)+N*(1.33*TcP) + k2*(0.239*v)*(0.22*TaP + 0.78*TdP);}// Open (< 10%)
    }
    return Ma;
  }
  //----------------------------------------------------------
  else if(Options.pot_melt == POTMELT_CRHM_EBSM)
  {
    double Qn_ebsm=(F->SW_radia_net+F->LW_radia);
    double Qh_ebsm=-0.92+0.076*F->wind_vel+0.19*F->temp_daily_max;
    double Qp_ebsm=(F->precip*(1.0-F->snow_frac))*F->temp_ave*4.2/1000.0;

    return  Qn_ebsm + Qh_ebsm + Qp_ebsm; 
  }
  return 0.0;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns potential melt rate [mm/d]
/// \brief Adapted from UBC Watershed model source code,
/// \copyright (c) Michael Quick
///
/// \param &Options [in] Global model options information
/// \param &F [in] Forcing functions
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &tt [in] Current model time
/// \return Potential melt rate [mm/d]
//
double UBC_DailyPotentialMelt(const optStruct &Options,
                              const force_struct &F,
                              const CHydroUnit *pHRU,
                              const time_struct &tt)
{
  double CONVM,CONDM,RAINMELT,potential_melt;
  double SHORM,ONGWMO;
  double pressure,vel,RI,RM;
  double Fc,orient;

  const double P0CONV= pHRU->GetSurfaceProps()->conv_melt_mult; //convection melt multiplier
  const double P0COND= pHRU->GetSurfaceProps()->cond_melt_mult;  //condensation melt multiplier
  const double P0RAIN= pHRU->GetSurfaceProps()->rain_melt_mult;  //rain melt multiplier

  double albedo   =pHRU->GetSnowAlbedo();
  double snowSWE  =pHRU->GetSnowSWE();

  if ((pHRU->GetHRUType()==HRU_GLACIER) && (snowSWE<=0))
  {
    double minalb=CGlobalParams::GetParams()->min_snow_albedo;
    albedo=1-(1.0-albedo)*(1-minalb);//UBCWM RFS implmentation
  }

  Fc    =pHRU->GetSurfaceProps()->forest_coverage;
  orient=1.0-fabs(pHRU->GetAspect()/PI-1.0);        //=0 for north, 1.0 for south

  //calculate shortwave energy component.
  SHORM=(1.0-albedo)*F.SW_radia;     //[MJ/m2/d]
  SHORM*=(MM_PER_METER/DENSITY_WATER/LH_FUSION);   //[MJ/m2/d-->mm/d]

  //longwave energy component
  ONGWMO=F.LW_radia;                                //[MJ/m2/d]
  ONGWMO*=(MM_PER_METER/DENSITY_WATER/LH_FUSION);  //[MJ/m2/d-->mm/d]

  //calculate reduction factor RM
  pressure=F.air_pres/KPA_PER_ATM;                  //[kPa-->atm]
  vel     =F.wind_vel/M_PER_KM*SEC_PER_HR;          //[m/s->km/hr]

  RI=(0.095*F.temp_daily_ave)/(vel*vel);//Richardson Number
  //above is a linearization of
  //RI=(2*GRAVITATIONAL_CONST*ZA/(F.temp_daily_ave+ZERO_CELSIUS)*F.temp_daily_ave/vel/vel;

  RM=min(max(1.0-7.7*RI,0.0),1.6);

  //convective melt
  CONVM=P0CONV*pressure*F.temp_daily_ave*RM*vel; //[mm/d]

  //condensation melt
  CONDM=P0COND*max(F.temp_daily_min,0.0)*RM*vel;  //[mm/d]
  //temp_daily_min should be dew point temp
  //to account for a (sketchy) UBC implementation detail

  CONDM*=((1.0-Fc)*pressure+(Fc)*1.0);  //[mm/d]

  //melt due to rain heat input
  if (Options.keepUBCWMbugs){//factor of 10 off
    RAINMELT=0.00125            *max(F.temp_daily_ave,0.0)*F.precip_daily_ave*(1.0-F.snow_frac); //[mm/d]
  }
  else{
    RAINMELT= P0RAIN*SPH_WATER/LH_FUSION*max(F.temp_daily_ave,0.0)*F.precip_daily_ave*(1.0-F.snow_frac); //[mm/d]
  }
  potential_melt= SHORM + ONGWMO + CONVM + CONDM + RAINMELT;             //[mm/d]

  return potential_melt;
}
