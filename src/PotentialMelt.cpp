/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2023 the Raven Development Team
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
                                     const potmelt_method method,
                                     const optStruct &Options,
                                     const CHydroUnit *pHRU,
                                     const time_struct &tt)
{

  //----------------------------------------------------------
  if (method==POTMELT_DATA)
  {
    return F->potential_melt; //already read into gauge time series & interpolated
  }
  //----------------------------------------------------------
  else if (method==POTMELT_DEGREE_DAY)
  {
    double Ma=pHRU->GetSurfaceProps()->melt_factor;
    double melt_temp=pHRU->GetSurfaceProps()->DD_melt_temp;
    return Ma*(F->temp_daily_ave-melt_temp)*F->subdaily_corr;
  }
  //----------------------------------------------------------
  else if ((method==POTMELT_HBV) || (method==POTMELT_HBV_ROS))
  {
    const surface_struct *surf=pHRU->GetSurfaceProps();

    double Ma    =surf->melt_factor;
    double Ma_min=surf->min_melt_factor;
    double AM    =surf->HBV_melt_asp_corr;
    double MRF   =surf->HBV_melt_for_corr;
    double Fc    =surf->forest_coverage;
    double melt_temp=surf->DD_melt_temp;
    double slope_corr;

    //annual variation
    double adj=0;
    if (pHRU->GetLatRad()<0.0){adj=PI;}//southern hemisphere phase shift
    Ma=Ma_min+(Ma-Ma_min)*0.5*(1.0-cos(F->day_angle-WINTER_SOLSTICE_ANG-adj));

    //forest correction
    Ma*=((1.0-Fc)*1.0+(Fc)*MRF);

    //slope correction (forest only, for some reason)
    slope_corr=((1.0-Fc)*1.0+(Fc)*sin(pHRU->GetSlope()));

    //aspect corrections
    Ma*=max(1.0-AM*slope_corr*cos(pHRU->GetAspect()-adj),0.0);//north facing slopes (aspect=0) have lower melt rates

    double rain_energy(0.0);
    if(method==POTMELT_HBV_ROS) {
      double ROS_mult= pHRU->GetSurfaceProps()->rain_melt_mult;  //rain melt multiplier
      rain_energy=ROS_mult*SPH_WATER/LH_FUSION*max(F->temp_ave,0.0)*(F->precip*(1.0-F->snow_frac)); //[mm/d]
    }

    return Ma*(F->temp_daily_ave-melt_temp)*F->subdaily_corr+rain_energy;
  }
  //-------------------------------------------------------------------------------------
  else if (method==POTMELT_BLENDED)
  {
    ExitGracefullyIf(_PotMeltBlends_N==0,
      "EstimatePotentialMelt: POTMELT_BLENDED specified, but no weights were provided using :BlendedPotMeltWeights command",BAD_DATA);
    double melt=0;
    for(int i=0; i<_PotMeltBlends_N;i++) {
      potmelt_method etyp=_PotMeltBlends_type[i];
      double           wt=_PotMeltBlends_wts[i];
      melt+=wt*EstimatePotentialMelt(F,etyp,Options,pHRU,tt);

      if(rvn_isnan(melt)) {
        ExitGracefully("EstimatePotentialMelt: NaN value produced in POTMELT_BLENDED calculation by one or more POTMELT routines",RUNTIME_ERR);return false;
      }
    }
    return melt;
  }
  //----------------------------------------------------------
  else if (method==POTMELT_EB)
  {
    double rainfall  = F->precip*(1-F->snow_frac);
    double air_temp  = F->temp_ave;
    double air_pres  = F->air_pres;
    double rel_humid = F->rel_humidity;
    double wind_vel  = F->wind_vel;

    double ref_ht    = 2.0; // [m];
    double roughness = pHRU->GetSurfaceProps()->roughness ;
    double surf_temp = pHRU->GetSnowTemperature();

    double K = F->SW_radia_net;   //net short wave radiation to snowpack [MJ/m2/d]
    double L = F->LW_radia_net;   //net long wave radiation [MJ/m2/d]
    double H = GetSensibleHeatSnow(air_temp,surf_temp,wind_vel,ref_ht,roughness); //[MJ/m2/d]
    double LE= GetLatentHeatSnow  (air_pres,air_temp,surf_temp,rel_humid,wind_vel,ref_ht,roughness); //[MJ/m2/d]
    double R = GetRainHeatInput   (surf_temp,air_temp,rainfall ,rel_humid); //[MJ/m2/d]

    double S = K + L + H + LE + R;                //total heat input rate [MJ/m2/d]
    S*=(MM_PER_METER/DENSITY_WATER/LH_FUSION);    //[MJ/m2/d-->mm/d]
    return S;
  }
  //----------------------------------------------------------
  else if (method==POTMELT_RESTRICTED)
  {
    double rad,convert;
    double melt_temp=pHRU->GetSurfaceProps()->DD_melt_temp;
    double Ma       =pHRU->GetSurfaceProps()->melt_factor;
    double tmp_rate=threshPositive(Ma*(F->temp_daily_ave-melt_temp));

    rad     = (F->SW_radia_net+F->LW_radia_net);//[MJ/m2/d]
    convert = MM_PER_METER/DENSITY_WATER/LH_FUSION;//for converting radiation to mm/d [mm/m]/[kg/m3]/[MJ/kg] = [mm-m2/MJ]

    return tmp_rate*F->subdaily_corr+ threshPositive(rad*convert);
  }
  //----------------------------------------------------------
  else if (method==POTMELT_DD_RAIN)
  {
    //degree-day with rain-on-snow correction
    double melt_temp=pHRU->GetSurfaceProps()->DD_melt_temp;
    double Ma       =pHRU->GetSurfaceProps()->melt_factor;
    double rainfall=F->precip*(1-F->snow_frac);
    double adv_melt= SPH_WATER/LH_FUSION*max(F->temp_ave,0.0)*rainfall; //[mm/d]
    return Ma*(F->temp_daily_ave-melt_temp)*F->subdaily_corr+adv_melt;
  }
  //----------------------------------------------------------
  else if (method==POTMELT_UBCWM)
  {
    return UBC_DailyPotentialMelt(Options,*F,pHRU,tt)*F->subdaily_corr;
  }
  //----------------------------------------------------------
  else if (method == POTMELT_USACE)
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
    return Ma*F->subdaily_corr;
  }
  //----------------------------------------------------------
  else if(method == POTMELT_CRHM_EBSM)
  {
    double Qn_ebsm=(F->SW_radia_net+F->LW_radia_net);
    double Qh_ebsm=(-0.92+0.076*F->wind_vel+0.19*F->temp_daily_max);
    double Qp_ebsm=(F->precip*(1.0-F->snow_frac))*max(F->temp_ave,0.0)*SPH_WATER*DENSITY_WATER/MM_PER_METER;

    return  (Qn_ebsm + Qh_ebsm + Qp_ebsm)/LH_FUSION/DENSITY_WATER*MM_PER_METER; //[MJ/m2/d]->[mm/d]
  }
  //----------------------------------------------------------
  else if(method == POTMELT_HMETS)
  {
    double Ma;
    double Ma_min   =pHRU->GetSurfaceProps()->min_melt_factor;
    double Ma_max   =pHRU->GetSurfaceProps()->max_melt_factor;
    double melt_temp=pHRU->GetSurfaceProps()->DD_melt_temp;
    double Kcum     =pHRU->GetSurfaceProps()->DD_aggradation;
    int iCumMelt=GetStateVarIndex(CUM_SNOWMELT);
    
    double cum_melt=0.0;
    if(iCumMelt!=DOESNT_EXIST){cum_melt= pHRU->GetStateVarValue(iCumMelt); }

    Ma=min(Ma_max,Ma_min*(1+Kcum*cum_melt));

    return Ma*(F->temp_daily_ave-melt_temp)*F->subdaily_corr;
  }
  //----------------------------------------------------------
  else if(method == POTMELT_RILEY)
  {
    //From Riley, J.P., E.K. Israelsen, and K.O. Eggleston, some approaches to snowmelt prediction, 
    //Actes du Colloque de Banff sur le role de la niege et de la glace en hydrologie, AISH Pub, Vol 2, no 107, p956-971, 1972
    //As documented in HYDROTEL 2.1 manual
    double Ma       =pHRU->GetSurfaceProps()->melt_factor;
    double melt_temp=pHRU->GetSurfaceProps()->DD_melt_temp;
    double rainfall =F->precip*(1.0-F->snow_frac);
    double slope_fact=F->ET_radia/F->ET_radia_flat;
    double snoalb  = pHRU->GetSnowAlbedo(); // Average snow surface albedo
    double adv_melt= SPH_WATER/LH_FUSION*max(F->temp_ave,0.0)*rainfall; //[mm/d]

    return Ma*slope_fact*(F->temp_ave-melt_temp)*(1.0-snoalb)+adv_melt;
  }
  //----------------------------------------------------------
  else if(method == POTMELT_NONE)
  {
    return 0.0;
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
  double Fc;

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

  //calculate shortwave energy component.
  SHORM=(1.0-albedo)*F.SW_radia_subcan;            //[MJ/m2/d] //Uses subcanopy correction
  SHORM*=(MM_PER_METER/DENSITY_WATER/LH_FUSION);   //[MJ/m2/d-->mm/d]

  //longwave energy component
  ONGWMO=F.LW_radia_net;                           //[MJ/m2/d]
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
