/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2018 the Raven Development Team
  ----------------------------------------------------------------
  Routines for calculating radiation
  ----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Radiation.h"
#include "GlobalParams.h"
#include "Model.h"  // required for definition of CModel

//////////////////////////////////////////////////////////////////
/// \brief Estimates incoming clear sky shortwave radiation SW_radia and extraterrestrial radiation ET_rad [MJ/m2/d]
/// \details Returns the clear sky total solar radiation (total incident shortwave radiation)
/// \remark Units are Kcs in Dingman text
///
/// \param Options [in] global options structure
/// \param *F [in] Forcing functions for a particular HRU over current time step
/// \param *pHRU [in] pointer to HRU for which radiation is calculated
/// \param tt [in] current time structure
/// \param ET_rad [out] ET radiation [MJ/m2/d]
/// \return Double value representing shortwave radiation [MJ/m2/d]
//
double CRadiation::EstimateShortwaveRadiation(CModel* pModel,
                                              const force_struct *F,
                                              const CHydroUnit   *pHRU,
                                              const time_struct  &tt,
                                              double             &ET_rad,
                                              double             &ET_rad_flat)
{
  const optStruct* Options = pModel->GetOptStruct();
  double latrad=pHRU->GetLatRad();
  double lateq =pHRU->GetLatEq();
  switch(Options->SW_radiation)
  {
  //--------------------------------------------------------
  case(SW_RAD_NONE):
  {
    return 0.0;
  }
  //--------------------------------------------------------
  case(SW_RAD_DATA):
  //--------------------------------------------------------
  case(SW_RAD_DEFAULT):
  {
    double dew_pt    =GetDewPointTemp(F->temp_ave,F->rel_humidity);
    double slope     =pHRU->GetSlope();
    double solar_noon=pHRU->GetSolarNoon();
    double aspect    =pHRU->GetAspect();

    double SWrad= ClearSkySolarRadiation(tt.julian_day, Options->timestep,
                                         latrad, lateq, slope, aspect,
                                         F->day_angle, F->day_length,
                                         solar_noon, dew_pt, ET_rad, ET_rad_flat,
                                         (Options->timestep >= 1.0));
    if (Options->SW_radiation==SW_RAD_DATA){return F->SW_radia;} //ensures ET_rad still calculated!
    else                                   {return SWrad;}

  }
  //--------------------------------------------------------
  case(SW_RAD_UBCWM):
  {
    //calculates DAILY MEAN shortwave
    double solar_rad;
    double elev=pHRU->GetElevation();
    solar_rad= UBC_SolarRadiation(latrad,elev,F->day_angle, F->day_length,ET_rad);
    ET_rad_flat=ET_rad;
    double orient=1.0-fabs(pHRU->GetAspect()/PI-1.0);        //=0 for north, 1.0 for south
    if(pHRU->GetLatRad()<0.0) { orient=1.0-orient; }//southern hemisphere phase shift
    double shortwave_corr_S = InterpolateMo(pModel->GetGlobalParams()->GetParams()->UBC_s_corr, tt,Options->month_interp, *Options);
    double shortwave_corr_N = InterpolateMo(pModel->GetGlobalParams()->GetParams()->UBC_n_corr, tt,Options->month_interp, *Options);
    double shortwave_corr=((orient)* shortwave_corr_S + (1.0-orient)*shortwave_corr_N);
    ET_rad*=shortwave_corr;
    return solar_rad * shortwave_corr;
  }
  //--------------------------------------------------------
  case(SW_RAD_VALIANTZAS) :
  {
    double declin = SolarDeclination(F->day_angle);
    double ws = acos(-tan(latrad)*tan(declin));
    double ecc = EccentricityCorr(F->day_angle);
    ET_rad = 37.59*ecc*(ws*sin(latrad)*sin(declin) + sin(ws)*cos(latrad)*cos(declin)); //Extraterrestrial radiation
    ET_rad_flat=ET_rad; //no slope correction in this method
    return  min(0.75+0.00002*pHRU->GetElevation(),1.0)*ET_rad;
  }

  }
  return 0;
}

//////////////////////////////////////////////////////////////////
/// \brief Estimates net longwave radiation [MJ/m2/d]
///
/// \param Options [in] global options structure
/// \param *F [in] Forcing functions for a particular HRU over current time step
/// \param *pHRU [in] pointer to HRU for which radiation is calculated
/// \return Double value representing longwave radiation [MJ/m2/d]
//
double CRadiation::EstimateLongwaveRadiation(const int iSnow,
                                             CModel* pModel,
                                             const force_struct  *F,
                                             const CHydroUnit    *pHRU,
                                             double        &LW_incoming)
{
  const optStruct *Options = pModel->GetOptStruct();
  switch(Options->LW_incoming)
  {
    case(LW_INC_DATA):
    {
      LW_incoming=F->LW_incoming;
      break;
    }
    case(LW_INC_DEFAULT):
    {
      LW_incoming=0.0; //calculated below with LW_DEFAULT, or not at all for algorithms which can calculate net LW only
      break;
    }
    //--------------------------------------------------------
    case(LW_INC_SICART):
    {
      //Calculates incoming long-wave terrain emmission radiation using sky view factor
      //Sicart et al. (2005), Incoming longwave radiation to melting snow: observations, sensitivity and estimation in
      //northern environments, Hydrological Processes, 20, 3697-3708
      //Ported over from CRHM (Pomeroy et al., 2007), added forest cover term

      double Fc        =pHRU->GetSurfaceProps()->forest_coverage;
      double svf       =pHRU->GetSurfaceProps()->sky_view_factor;
      double epsilon_s =0.98;//vegetation emissivity
      double tau       =0.8; //atmospheric transmittance
      if(F->ET_radia >= 0.001) { tau =  (F->SW_radia_unc/F->ET_radia); } //Based upon comments, not code

      double Tair      =F->temp_ave + ZERO_CELSIUS;
      double ea        =F->rel_humidity*GetSaturatedVaporPressure(F->temp_ave);//kPa

      double epsilon_air= 1.24 * pow(ea*HPA_PER_KPA/Tair,1.0/7.0); //Brutsaert 1975

      svf=((svf)*Fc + 1.0 * (1 - Fc)); //correct for partial forest cover

      double L_0,L_F; //incoming from atmos [MJ/m2/d], incoming from vegetation [MJ/m2/d]
      L_0 = epsilon_air*(1.0 + 0.44*F->rel_humidity - 0.18*tau)*WATT_TO_MJ_PER_D;//eqn 9 of Sicart et al (2005)
      L_F = epsilon_s*STEFAN_BOLTZ*pow(Tair,4.0); // eqn 6 of Sicart et al (2005)

      LW_incoming = svf*L_0 + (1.0-svf)*L_F;
      break;
    }
    //--------------------------------------------------------
    case (LW_INC_SKYVIEW):
    {
      //simple sky view factor with air emissivity from
      //Prata, A. J. A new long wave formula for estimating downward clear sky radiation at the surface.
      //Quartely Journal of the Royal Meteorological Society, v. 122, p. 1127 - 1151,1996.
      double Fc = pHRU->GetSurfaceProps()->forest_coverage;
      double svf =pHRU->GetSurfaceProps()->sky_view_factor;//0.7
      double epsilon_s=0.95;
      double epsilon_air=0.18;
      double Tair=F->temp_ave+ZERO_CELSIUS;

      double e_0  = F->rel_humidity*GetSaturatedVaporPressure(F->temp_ave); //[kPa]
      double eta  = 46.5*(e_0*HPA_PER_KPA/Tair); //[g/cm2] - Prata 1996 eqn 16
      epsilon_air = 1.0 - (1 + eta) * exp(-sqrt(1.2+3*eta)); //Prata 1996 eqn 17

      //epsilon_air = 1.24 * pow(e_0*HPA_PER_KPA/Tair,1.0/7.0); //Brutsaert 1975

      svf=((svf)*Fc + 1.0 * (1 - Fc)); //correct for partial forest cover

      LW_incoming =  epsilon_s   * (1 - svf) * STEFAN_BOLTZ * pow(Tair,4);
      LW_incoming += epsilon_air * (    svf) * STEFAN_BOLTZ * pow(Tair,4);
      break;
    }
    //--------------------------------------------------------
    case (LW_INC_DINGMAN):
    {
      //from Dingman eqns. 5-40
      double forest_cover=pHRU->GetSurfaceProps()->forest_coverage;

      double Tair =F->temp_ave+ZERO_CELSIUS;  //[K]
      double Tsurf=F->temp_ave+ZERO_CELSIUS;  //[K] //\todo better way to do this.
      //double Tsurf=pHRU->GetSurfaceTemperature()+ZERO_CELSIUS;//[K] (not yet functional - also surface depends upon whether this is used for ET or melt)

      double ea=F->rel_humidity*GetSaturatedVaporPressure(F->temp_ave);

      double emiss_eff=1.72*pow(ea/Tair,1/7.0);/// effective clear sky emmisivity - Brutsaert, 1975 in units of kPa/K

      //emiss_eff = 1.08 *(1.0- exp(-(10*ea)*Tair/2016.0));            ///< Satterlund (1979) via Brook90 documentation \cite Satterlund1979WRR
      //emiss_eff = 1.0- 0.261*exp(-0.000777*F->temp_ave*F->temp_ave); ///< Idso and Jackson (1969) via Brook90 documentation \cite Idso1969JoGR
      //emiss_eff = 0.0000092*Tair*Tair;                               ///< Swinbank (1963) via Brook90 documentation \cite Swinbank1963QJotRMS

      double eps_at=emiss_eff*(1.0+0.22*F->cloud_cover*F->cloud_cover); //cloud cover correction

      eps_at=(1-forest_cover)*eps_at+(forest_cover)*1.0; //treats forest as blackbody - neglects sky view factor

      LW_incoming=STEFAN_BOLTZ*eps_at*pow(Tair,4);
      break;
    }
  }

  switch(Options->LW_radiation)
  {
  //--------------------------------------------------------
  case(LW_RAD_DATA):
  {
    return F->LW_radia_net;
  }
  //--------------------------------------------------------
  case(LW_RAD_NONE):
  {
    return 0.0;
  }
  //--------------------------------------------------------
  case(LW_RAD_DEFAULT):
  {
    //from Dingman eqns. 5-40
    double emissivity=0.99;
    double forest_cover=pHRU->GetSurfaceProps()->forest_coverage;
    double Tair =F->temp_ave+ZERO_CELSIUS;  //[K]
    double Tsurf=F->temp_ave+ZERO_CELSIUS;  //[K] //\todo better way to do this.
    //double Tsurf=pHRU->GetSurfaceTemperature()+ZERO_CELSIUS;//[K] (not yet functional)

    double ea   =F->rel_humidity*GetSaturatedVaporPressure(F->temp_ave);

    double emiss_eff=1.72*pow(ea/Tair,1/7.0);/// effective clear sky emmisivity \ref Brutsaert 1975, via Dingman clear sky emmissivity

    double eps_at   =emiss_eff*(1.0+0.22*F->cloud_cover*F->cloud_cover);

    eps_at=(1-forest_cover)*eps_at+(forest_cover)*1.0; //treats forest as blackbody - neglects sky view factor

    if (Options->LW_incoming==LW_INC_DEFAULT) { //now, DEFAULT==DINGMAN
      LW_incoming=STEFAN_BOLTZ*eps_at*pow(Tair,4);
    }
    return LW_incoming-STEFAN_BOLTZ*emissivity*pow(Tsurf,4);
  }
  //--------------------------------------------------------
  case (LW_RAD_UBCWM):
  {
    double cloud=F->cloud_cover;
    double LW_open=(1-cloud)*(-20+0.94*F->temp_daily_ave)+(cloud)*(1.24*F->temp_daily_min);//[mm/d]

    double P0BLUE = pModel->GetGlobalParams()->GetParams()->UBC_LW_forest_fact;//actually P0BLUE * P0LWVF , [mm/d/K]
    double LW_forest=P0BLUE*F->temp_daily_ave; //[mm/d]

    double tmp=(LH_FUSION*DENSITY_WATER/MM_PER_METER); //[mm/d]-->[MJ/m2/d]

    double Fc=pHRU->GetSurfaceProps()->forest_coverage;

    return tmp*((1.0-Fc)*LW_open+(Fc)*LW_forest);
  }
  //--------------------------------------------------------
  case (LW_RAD_HSPF):
  {
    double shade=0.0; //Shade factor (land surface parameter?)
    double LW;

    if (F->temp_ave>0){LW=(shade)*0.468*(F->temp_ave)+(1.0-shade)*(0.360*F->temp_ave-6.6);}
    else              {LW=(shade)*0.360*(F->temp_ave)+(1.0-shade)*(0.306*F->temp_ave-6.6);}
    if (LW<0){LW*=(1.0-F->cloud_cover);}

    LW*=HR_PER_DAY*MJ_PER_M2_LANGLEY;//convert Langleys/hr->MH/m2/day
    return LW;
  }
  //--------------------------------------------------------
  case(LW_RAD_VALIANTZAS):
  {
    //from Valiantzas, 2006, Simplified versions for the Penman evaporation equation using routine weather data, Journal of Hydrology, 331
    double f,eps,ea;
    ea =F->rel_humidity*GetSaturatedVaporPressure(F->temp_ave);                    //[kPa]

    f=max(1.35*(F->SW_radia/F->SW_radia_unc)-0.35,0.0); //cloud cover adjustment Valiantzas (2006) eqn 40
    eps= 0.34 - 0.14*sqrt(ea);                          //net emissivity Valiantzas (2006) eqn 41

    return -f*eps*STEFAN_BOLTZ*pow(F->temp_ave+ZERO_CELSIUS,4.0); //net LW Valiantzas (2006) eqn 3
  }
  //--------------------------------------------------------
  /*case (LW_RAD_UNFAO):
  {
    ///< from Crop evapotranspiration - Guidelines for computing crop water requirements - FAO Irrigation and drainage paper 56 \cite Allen1998FR
    //http://www.fao.org/docrep/X0490E/X0490E00.htm
    double ea =F->rel_humidity*GetSaturatedVaporPressure(F->temp_ave); //[kPa]
    double tmp=(0.34-0.14*pow(ea,0.5))*(1.35*(F->SW_radia/F->SW_uncorr)-0.35);
    return STEFAN_BOLTZ*tmp*(pow(F->temp_daily_min,4)+pow(F->temp_daily_max,4));
  }*/
  //--------------------------------------------------------
  /*
  case(LW_RAD_BRUTSAERT):
  { //Brutsaert (1982) equation for effective clear sky emissivity
    double cloud_corr=0.2+0.8*min(max((F->SW_radia/F->SW_radia_unc - 0.25) / 0.5,1.0),0.0);
    return STEFAN_BOLTZ*emmissivity*cloud_corr*(emiss_eff*pow(Tair,4)-pow(Tair,4));
  }
  */
  //--------------------------------------------------------
   /*
   case (LW_RAD_BRUNT):
   { //from CRHM (Pomeroy, 2007) netall routine
     // based upon Brunt (1944), Physical and dynamical meteorology: Cambridge Univ. Press
     double LW_net = -0.85;
     double emiss=0.97;
     double ea =F->rel_humidity*GetSaturatedVaporPressure(F->temp_ave);                    //[kPa]
     if(clear_sky_SW > 0.0)
     {
        LW_net = -0.85 - emiss*STEFAN_BOLTZ*pow(F->temp_ave+ZERO_CELSIUS,4.0)*(0.39-0.093*sqrt(ea))*(0.26+0.81*(F->SW_radia/F->SW_radia_unc));
        LW_incoming = 0.85;
     }
     return LW_net;
   }
  */
  /*case(LW_RAD_SWAT):
    {
    // calculate net long-wave radiation
    double emiss,cloud_corr;
    double sat_vap,ea;
    sat_vap=GetSaturatedVaporPressure(F->temp_ave); //[kPa]
    ea =F->rel_humidity*sat_vap;                    //[kPa]
    emiss      = (0.34 - 0.139 * sqrt(ea));                  // net emissivity  equation 2.2.20 in SWAT manual
    cloud_corr = 0.1+ 0.9 * (F->SW_radia / F->SW_radia_unc);  // cloud cover factor equation 2.2.19 SWAT
    LW_incoming=(0.98-emiss*cloud_corr) * STEFAN_BOLTZ * pow(F->temp_ave+ZERO_CELSIUS,4);
    return -emiss * cloud_corr * STEFAN_BOLTZ * pow(F->temp_ave+ZERO_CELSIUS,4);   // net long-wave radiation equation 2.2.21 SWAT
    }*/

  }
  return 0.0;
}

/////////////////////////////////////////////////////////////////////
/// \brief Estimates the effect of cloud cover on short wave radiation.
/// \details Caluculates the reduction of shortwave radiation due to cloud cover
/// \remark The default is to apply the UBCWM correction. The UBCWM equation is identical
/// \to Dingman (2008) Eq. 5-31 but uses a different semantic: The UBCWM 'cloud penetration factor POCAST'
/// \is identical to  'cloud height' in Dingman (2008) Eq. 5-31.
/// \param Options [in] global options structure
/// \param F [in] forcings structure
/// \param elev [in] elevation of HRU [masl]
/// \return Double Cloud cover correction factor for shortwave radiation

double CRadiation::SWCloudCoverCorrection(CModel* pModel,
                                          const force_struct *F,
                                          const double &elev)
{
  const optStruct *Options = pModel->GetOptStruct();
  switch(Options->SW_cloudcovercorr)
  {
  case(SW_CLOUD_CORR_UBCWM):
  {
    double POCAST = pModel->GetGlobalParams()->GetParams()->UBC_cloud_penet;
    return 1.0 -(1.0 - POCAST) * F->cloud_cover;
  }
  case(SW_CLOUD_CORR_DINGMAN): // Dingman (2008) Eq. 5-30 (does not require parameters)
  {
    return 0.355 + 0.68 * (1 - F->cloud_cover);
  }
  case(SW_CLOUD_CORR_ANNANDALE):
  // Annandale et al (2002), Software for missing data error analysis of Penman-Monteith
  // reference evapotranspiration, Irrigation Science 21(2),pp57-67
  {
    double kRs=0.16;   // interior = ~0.16/coastal~=0.19
    if(F->temp_daily_max>F->temp_daily_min) {
      return kRs*(1.0 + 2.7E-5*elev)*sqrt(F->temp_daily_max - F->temp_daily_min);
    }
    else {return 0.0;}
  }
  default: // apply no cloud cover correction (default)
  {
    return 1.0;
  }
  }
}


/////////////////////////////////////////////////////////////////////
/// \brief Calculates the ratio of solar radiation under forest canopy relative to open.
/// \details Caluculates the transmittance of both direct and diffuse radiation
/// \remark The default is to apply no correction so that in the absence of a canopy no canopy
/// \       related paramters need to be supplied (other than the required ones)
///
/// \param Options [in] global options structure
/// \param *pHRU [in] pointer to HRU for which radiation is calculated
/// \return Double Canopy correction factor for shortwave radiation

double CRadiation::SWCanopyCorrection(CModel *pModel,
                                      const CHydroUnit *pHRU)
{
  const optStruct *Options = pModel->GetOptStruct();
  double Fc=pHRU->GetSurfaceProps()->forest_coverage;

  switch(Options->SW_canopycorr)
  {
  case(SW_CANOPY_CORR_STATIC):
  {
    double LAI = pHRU->GetVegVarProps()->LAI;
    double SAI = pHRU->GetVegVarProps()->SAI;
    double extinction = pHRU->GetVegetationProps()->svf_extinction;

    return (1.0-Fc)*1.0+Fc*exp(-extinction*(LAI+SAI)); // Dingman (2008) Eq. 5-33
  }
  case(SW_CANOPY_CORR_UBCWM):  // A simple factor that switches on when forest cover is greater than zero
  {
    double shortwave_corr = 1.0;
    double UBC_correction_fact = pModel->GetGlobalParams()->GetParams()->UBC_exposure_fact; //Sun exposure factor of forested areas
    shortwave_corr*=((Fc)*UBC_correction_fact+(1.0-Fc)*1.0);
    return shortwave_corr;
  }
  case(SW_CANOPY_CORR_DYNAMIC):  // stub for Jost and Moore (2010) equations
  {
    ExitGracefully("SWCanopyCorrection::SW_CANOPY_CORR_DYNAMIC",STUB);
    return 0.0;
  }
  case(SW_CANOPY_CORR_NONE):
  {
    return 1.0;
  }
  default: // apply no correction for canopy (default)
  {
    return 1.0;
  }
  }
  return 1.0;
}
//////////////////////////////////////////////////////////////////
/// \brief Calculates day angle [rad]
/// \param &julian_day [in] Julian identifier of a specfic day
/// \param year        [in] Year in which this day occurs
/// \param calendar    [in] enum int of calendar used
/// \return The day angle of this day [rad]
//
double CRadiation::DayAngle(const double&julian_day,
                            const int    year,
                            const int    calendar)
{
  double leap=0.0;
  if (IsLeapYear(year,calendar)){leap=1.0;}
  return 2.0*PI*(julian_day/(365.0+leap));
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates day length, given a latitude and a solar declination [d]
/// \param lat [in] Latitude [rad]
/// \param declin [in] Solar declination [rad]
/// \return Day length corresponding to these parameters [d]
//
double CRadiation::DayLength(const double lat,    //latitude (in radians)
                             const double declin) //solar declination (in radians)
{
  double hd = -tan(declin) * tan(lat);
  if      (hd>= 1.0) {return 0.0;}//sun stays below horizon
  else if (hd<=-1.0) {return 1.0;}//sun stays above horizon
  return acos(hd)/PI;
}

//////////////////////////////////////////////////////////////////
/// \brief  function thatdetermines solar noon, in days, on a surface with arbitrary slope/aspect
/// \ref from Dingman E-22 \cite Dingman1994
/// \param lat [in] Latitude [rad]
/// \param slope [in] slope [rad]
/// \param aspect [in] aspect [rad counterclock from north]
///
double CRadiation::CalculateSolarNoon(const double &latrad,
                                      const double &slope,
                                      const double &aspect)
{
  double solar_noon; //[d]

  double denom=cos(slope)*cos(latrad) - sin(slope)*sin(latrad)*cos(aspect);
  if (denom==0.0){denom = REAL_SMALL;}
  solar_noon = -atan(sin(slope)*sin(-aspect)/denom)/EARTH_ANG_VEL; // Dingman eqn E-22
  if (solar_noon> 0.5){solar_noon-=1.0;}
  if (solar_noon<-0.5){solar_noon+=1.0;}

  return solar_noon;
}
//////////////////////////////////////////////////////////////////
/// \brief  function thatdetermines solar noon, in days, on a surface with arbitrary slope/aspect
/// \ref from Dingman E-22 \cite Dingman1994
///
double CRadiation::CalculateEquivLatitude(const double &latrad,
                                          const double &slope,
                                          const double &aspect)
{
  return asin(min(cos(slope)*sin(latrad) + sin(slope)*cos(latrad)*cos(aspect),1.0)); //Dingman eqn E-23
}
//////////////////////////////////////////////////////////////////
/// \brief Returns solar declination [rad]
/// \ref Uses Dingman eqn E-3 \cite Dingman1994
/// \param day_angle Day angle, in radians
/// \return Solar declination [rad]
//
double CRadiation::SolarDeclination(const double day_angle)
{
  /*return      0.006918-
    0.399912*cos(1.0*day_angle)+
    0.070257*sin(1.0*day_angle)-
    0.006758*cos(2.0*day_angle)+
    0.000907*sin(2.0*day_angle)-
    0.002697*cos(3.0*day_angle)+
    0.001480*sin(3.0*day_angle);*///Dingman eqn E-3
  //return asin(0.39785* sin(4.868961+0.017203*julian_day+
  //                             0.033446*sin(6.224111+0.017202*julian_day)));//Brook90
  //return 0.409*sin(day_angle-1.39); //Valiantzas, 2006
  return 0.40928*sin(day_angle-1.384);//UBCWM
}


//////////////////////////////////////////////////////////////////
/// \brief Returns the eccentricity correction, E0 [-]
/// \param day_angle [in] Day angle [rad]
/// \return Eccentricity correction, E0 [-]
//
double CRadiation::EccentricityCorr(const double day_angle)
{
  /// \ref Dingman eqn E-2
  return  1.000110+
    0.034221*cos(1.0*day_angle)+
    0.001280*sin(1.0*day_angle)+
    0.000719*cos(2.0*day_angle)+
    0.000077*sin(2.0*day_angle);
  //return  pow(1.0 - 0.0167 * cos(0.0172 * (julian_day - 3)),-2);//Brook90
  //return  1+0.033*cos(day_angle);//WATFLOOD, Valiantzas (2006)
}
//////////////////////////////////////////////////////////////////
/// \brief LOCAL function that applies a correction to extraterrestrial radiation
/// \details Applies multiplication factor to solar constant and eccentricity factor for daily incoming radiation
/// \ref from Dingman E-25 or E-6/E-7 \cite Dingman1994
/// \return Correction factor for shortwave radiation
///
double CRadiation::RadCorr (const double declin,     //declination, radians
                            const double solar_noon, //time of solar noon, days
                            const double lat_eq,     //equivalent latitude, radians
                            const double t_sunrise,  //time of sunrise, days before solar noon
                            const double t_sunset,   //time of sunrise, days after solar noon,
                            const double t_sol,      //time of day w.r.t. solar noon, days (not adjusted for slope)
                            const bool   avg_daily)
{
  if(t_sunset<=t_sunrise){return 0.0;}
  if (avg_daily)
  {
    //averaged daily correction - just integrating the point value
    return sin(declin)*sin(lat_eq)*(t_sunset-t_sunrise) +
           cos(declin)*cos(lat_eq)*(sin(EARTH_ANG_VEL*(t_sunset -solar_noon))-
                                    sin(EARTH_ANG_VEL*(t_sunrise-solar_noon)))/EARTH_ANG_VEL;
  }
  else if ((t_sol<=t_sunrise) || (t_sol>=t_sunset))//nighttime
  {
    return 0.0;
  }
  else{//eqn E-6 Dingman (instantaneous)
    return sin(declin)*sin(lat_eq) +
           cos(declin)*cos(lat_eq)*cos(EARTH_ANG_VEL*(t_sol-solar_noon));
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates extraterrestrial radiation on an inclined plane
/// \details Returns total incoming solar radation, K_{ET}, in [MJ/m2/d], corrected for slope and aspect
/// \ref from Swift (1976), Lee (1964) via Dingman \cite Swift1976WRR \cite Lee1964IJoB
///
/// \param latrad [in] Latitude in radians
/// \param lateq [in] Equivalent latitude for slope (Dingman E-23)
/// \param declin [in] Solar declination in radians
/// \param ecc [in] Eccentricity correction
/// \param solar_noon, //effective solar noon corrected for slope [days]
/// \param day_length, //length of day (not corrected for slope) [days]
/// \param t_sol [in] Time of day with respecto solar noon [d]
/// \param avg_daily [in] True if average daily total incoming solar radiation is to be computed instead
/// \return Double total extraterrestrial radiation
//
double CRadiation::CalcETRadiation( const double latrad,     //latitude in radians
                                    const double lateq,      //equivalent latitude for slope/aspect in radians
                                    const double declin,     //solar declination in radians
                                    const double ecc,        //eccentricity correction [-]
                                    const double slope,      //in radians
                                    const double solar_noon, //effective solar noon corrected for slope [days]
                                    const double day_length, //length of day (not corrected for slope) [days]
                                    const double t_sol,      //time of day w.r.t. solar noon [days] [-0.5..0.5]
                                    const bool   avg_daily)  //true if this is to be averaged over the day
{
  double ExTerCorr;      //Extraterrestrial radiation correction [-]
  double half_day;       //half day length [days]
  double lat_eq;         //slope-corrected equivalent latitude [rad]
  double t_rise,t_set;   //time of sunrise, sunset (w.r.t. solar noon) [days] [-0.5..0.5]
  double t_rise1,t_set1;
  double t_rise2,t_set2;
  bool   TwoDawns;

  //---slope corrections-------------------------------------------
  lat_eq =latrad;
  if (slope!=0.0){
    lat_eq = lateq;     //equivalent latitude for slope/aspect (Dingman E-23)
  }

  //calculate sunrise & sunset--------------------------------------------
  half_day=0.5*day_length;
  t_rise=-half_day;     //uncorrected sunrise/set: Dingman E-5a,b
  t_set = half_day;

  //sunrise/set corrected for slope & aspect
  if(slope !=0){half_day=0.5*DayLength(lat_eq, declin);}

  t_rise1=max(t_rise,solar_noon-half_day);//Dingman E-24a,b
  t_set1 =min(t_set ,solar_noon+half_day);
  if (t_rise1>t_set1) {t_rise1=t_set1=0.0;}//Correct for those endless nights...

  ExTerCorr=RadCorr(declin, solar_noon, lat_eq, t_rise1, t_set1,t_sol,avg_daily);

  //Corrections for two sunrises--------------------------------------------
  TwoDawns = false;
  double fnoon=0; //fake noon
  t_rise2=solar_noon+1.0-half_day;
  t_set2 =solar_noon-1.0+half_day;
  if      (t_rise2<t_set ){t_set2 =t_set; fnoon=solar_noon+1.0; TwoDawns=true;}
  else if (t_set2 >t_rise){t_rise2=t_rise;fnoon=solar_noon-1.0; TwoDawns=true;}

  if (TwoDawns){
    //cout<<"TwoDawns: "<<t_set1-t_rise1<<" "<<t_set2-t_rise2<<" slope, lat"<< slope*RADIANS_TO_DEGREES<<" "<<latrad*RADIANS_TO_DEGREES<<endl;
    if (t_rise2>t_set2) {t_rise2=t_set2=0.0;}//Correct for those endless nights...
    ExTerCorr+=RadCorr(declin, fnoon, lat_eq, t_rise2,t_set2,t_sol,avg_daily);
  }
  //---end two sunrise/sunset correction-------------------------------------

  return SOLAR_CONSTANT*ecc*ExTerCorr;//[MJ/m2/d] E-25 dingman
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates extraterrestrial radiation on an inclined or horizontal plane
/// \details Returns total incoming solar radation, K_{ET}, in [MJ/m2/d], corrected for slope and aspect
/// \ref from Allen, Trezza & Tasumi (2006)
///
/// \param latrad [in] Latitude in radians
/// \param lateq [in] Equivalent latitude for slope (Dingman E-23)
/// \param declin [in] Solar declination in radians
/// \param ecc [in] Eccentricity correction

/// \param solar_noon, //effective solar noon corrected for slope [days]
/// \param day_length, //length of day (not corrected for slope) [days]
/// \param t_sol [in] Time of day with respecto solar noon [d]

/// \param avg_daily [in] True if average daily total incoming solar radiation is to be computed instead

/// \return Double total extraterrestrial radiation [MJ/m2/d]
//
double CRadiation::CalcETRadiation2(const double &latrad,     //latitude in radians
                                    const double &aspect,     //aspect in radians clockwise from south
                                    const double &declin,     //solar declination in radians
                                    const double &ecc,        //eccentricity correction [-]
                                    const double &slope,      //in radians
                                    const double &t1,         //starttime of timeestep w.r.t. solar noon [days] [-0.5..0.5]
                                    const double &t2,         //endtime of timeestep w.r.t. solar noon [days] [-0.5..0.5] (must be > t1)
                                    const bool   avg_daily)
{
  double a, b, c;
  double ws1,ws2;        //sunrise and sunset on horizontal surface [rad]
  double noon_cos_theta; //cos of azimuth at solar noon [-]

  double revasp=-(aspect+PI);  if(revasp>2*PI){revasp-=2*PI;}//revised aspect

  double cosslope =cos(slope),   sinslope=sin(slope);
  double cosdeclin=cos(declin), sindeclin=sin(declin);
  double coslatrad=cos(latrad), sinlatrad=sin(latrad);
  double cosrevasp=cos(revasp), sinrevasp=sin(revasp);

  a = sindeclin * (coslatrad*sinslope*cosrevasp - sinlatrad*cosslope);//Eqn 11a, constant
  b = cosdeclin * (coslatrad*cosslope + sinlatrad*sinslope*cosrevasp);//Eqn 11b, constant
  c = cosdeclin * sinslope * sinrevasp;                               //Eqn 11c, constant

  double quad;
  quad = b*b + c*c - a*a;        //Quadratic function for eqn 13
  quad = max(0.0001, quad);      //Limit quadratic function to greater than 0 as per step A.4.i.

  //Find cos_theta at solar noon (used instead of sin(slope) > solar angle)
  noon_cos_theta =  (sindeclin * sinlatrad + cosdeclin * coslatrad)* cosslope
                  + (cosdeclin * sinlatrad - sindeclin * coslatrad)* sinslope * cosrevasp;

  // Find sunrise, sunset on flat surface
  //=================================================================================
  ws2 = acos(-sindeclin/cosdeclin*sinlatrad/coslatrad); //Eqn 8, Sunset [rad]
  ws1 = -ws2;                                           //Sunrise [rad]

  if(latrad>0){// (N. hemisphere)
    if      (declin + latrad > PI / 2.0)       { ws1 = -PI;  ws2 = PI; }//sun never setting
    else if (fabs(declin - latrad) > PI / 2.0) { ws1 = 0.0;  ws2 = 0.0;}//sun never rising
  }
  else if(latrad<0){//(S. hemisphere)
    if      (declin + latrad < -PI / 2.0)      { ws1 = -PI;  ws2 = PI; }//sun never setting
    else if (fabs(declin - latrad) < -PI / 2.0){ ws1 = 0;    ws2 = 0;  }//sun never rising
  }

  //Step A.2. Find beginning integration limit
  //=================================================================================
  double sin_w1,w1,w1x,cos_w1,cos_w1x,cos_ws1,w1_24(0.0);

  sin_w1 = (a * c - b * sqrt(quad)) / (b*b + c*c);    //Eqn 13a
  sin_w1 = max(-1.0, sin_w1);                         //Limit sin_w1 to between 1 and -1 as per step A.4.iii.
  sin_w1 = min( 1.0, sin_w1);                         //Limit sin_w1 to between 1 and -1 as per step A.4.iii.

  w1 = asin(sin_w1);
  w1x = -PI - w1;
  cos_w1  = -a + b * cos(w1)  + c * sin(w1);         //Eqn 14
  cos_w1x = -a + b * cos(w1x) + c * sin(w1x);        //Eqn 14

  // Apply conditions in A.2.iv and A.2.v
  cos_ws1 = -a + b * cos(ws1) + c * sin(ws1) ;       //Eqn 14
  if      ((cos_ws1 <= cos_w1) && (cos_w1  < 0.001)) {w1_24 = w1;}
  else if (cos_w1x > 0.001)                          {w1_24 = ws1;}
  else if (w1x <= ws1)                               {w1_24 = ws1;}
  else if (w1x > ws1)                                {w1_24 = w1x;}

  w1_24 = max(w1_24, ws1);                           //Sunrise [rad], Step A.2.vi

  //Step A.3. Find end integration limit
  //=================================================================================
  double sin_w2,w2,w2x,cos_w2,cos_w2x,cos_ws2,w2_24(0.0);
  sin_w2 = (a * c + b * sqrt(quad)) / (b *b + c *c); //Eqn 13b
  sin_w2 = max(-1.0, sin_w2);                        //Limit sin_w2 to between 1 and -1 as per step A.4.iii.
  sin_w2 = min( 1.0, sin_w2);                        //Limit sin_w2 to between 1 and -1 as per step A.4.iii.

  w2 = asin(sin_w2);
  w2x = PI - w2;
  cos_w2  = -a + b * cos(w2 ) + c * sin(w2 );        //Eqn 14
  cos_w2x = -a + b * cos(w2x) + c * sin(w2x);        //Eqn 14

  //Apply conditions in A.3.iv and A.3.v
  cos_ws2 = -a + b * cos(ws2) + c * sin(ws2) ;         //Eqn 14
  if      ((cos_ws2 <= cos_w2) && (cos_w2 < 0.001)) {w2_24 = w2;}
  else if (cos_w2x > 0.001)                         {w2_24 = ws2;}
  else if (w2x >= ws2)                              {w2_24 = ws2;}
  else if (w2x < ws2)                               {w2_24 = w2x;}

  w2_24 = min(w2_24, ws2);                           //Sunset [rad], Step A.3.vi

  if(w2_24 < w1_24) {//If sunrise is before sunset, then set equal (slope is always shaded) as per A.4.ii.
    w1_24 = w2_24;
  }

  double cos_theta(0.0); //average cosine of azimuth during interval from w1_24 to w2_24
  bool doublesunset(false);

  //Check for double sunrise/double sunset
  //if cos theta is negative at noon potential that double sunrise/sunset occurs
  //=================================================================================
  double cos_theta1,cos_theta2;
  if(noon_cos_theta < 0)
  {
    double sinA,sinB,AA,BB,w2_24b,w1_24b,cos_w2b,cos_w1b;
    // Double sunrise/double sunset can occur
    sinA = (a * c + b * sqrt(quad)) / (b*b + c*c);   //Eqn 44a
    sinA = max(-1.0,sinA);                           //Limit to between -1 and 1
    sinA = min( 1.0,sinA);
    sinB = (a * c - b * sqrt(quad)) / (b*b + c*c);   //Eqn 44b
    sinB = max(-1.0,sinB);                             //Limit to between -1 and 1
    sinB = min( 1.0,sinB);
    AA = asin(sinA);
    BB = asin(sinB);
    w2_24b = min(AA,BB);                             //Eqn 45
    w1_24b = max(AA,BB);                             //Eqn 46
    cos_w2b = -a + b * cos(w2_24b) + c * sin(w2_24b);//Eqn 47a
    cos_w1b = -a + b * cos(w1_24b) + c * sin(w1_24b);//Eqn 47b

    if ((cos_w2b < -0.001) || (cos_w2b > 0.001)) {
      w2_24b = -PI - w2_24b;                         //Eqn 48a
    }

    if ((cos_w1b < -0.001) || (cos_w1b > 0.001)) {
      w1_24b = PI - w1_24b;                          //Eqn 48b
    }

    //Confirm if double sunrise/sunset occurs
    double X;
    //From equation 49a/49 b. Does not match what is siad in point g but makes more sense.
    if((w2_24b >= w1_24) && (w1_24b <= w2_24))
    {

      X = sin(declin) * sin(latrad) * cos(slope) * (w1_24b - w2_24b)
        - sin(declin) * cos(latrad) * sin(slope) * cos(revasp) * (w1_24b - w2_24b)
        + cos(declin) * cos(latrad) * cos(slope) * (sin(w1_24b) - sin(w2_24b))
        + cos(declin) * sin(latrad) * sin(slope) * cos(revasp) * (sin(w1_24b) - sin(w2_24b))
        - cos(declin) * sin(slope) * sin(revasp) * (cos(w1_24b) - cos(w2_24b));

      if(X < 0) {
        //Dbl sunrise/sunset occurs

        if(!avg_daily){
          w1_24 =max(w1_24,t1*2.0*PI);
          w2_24b=min(w2_24b,t2*2.0*PI);
          if(w2_24b<w1_24){w2_24b=w1_24;}  //nighttime during this interval
        }

        cos_theta1 = sin(declin) * sin(latrad) * cos(slope) * (w2_24b - w1_24)
        - sin(declin) * cos(latrad) * sin(slope) * cos(revasp) * (w2_24b - w1_24)
        + cos(declin) * cos(latrad) * cos(slope) * (sin(w2_24b) - sin(w1_24))
        + cos(declin) * sin(latrad) * sin(slope) * cos(revasp) * (sin(w2_24b) - sin(w1_24))
        - cos(declin) * sin(slope) * sin(revasp) * (cos(w2_24b) - cos(w1_24)); // Eqn 5 - integrate between w1_24 and w2_24b

        if(!avg_daily){
          w1_24b=max(w1_24b,t1*2.0*PI);
          w2_24 =min(w2_24,t2*2.0*PI);
          if(w2_24<w1_24b){w2_24=w1_24b;} //nighttime during this interval
        }
        cos_theta2 = sin(declin) * sin(latrad) * cos(slope) * (w2_24 - w1_24b)
        - sin(declin) * cos(latrad) * sin(slope) * cos(revasp) * (w2_24 - w1_24b)
        + cos(declin) * cos(latrad) * cos(slope) * (sin(w2_24) - sin(w1_24b))
        + cos(declin) * sin(latrad) * sin(slope) * cos(revasp) * (sin(w2_24) - sin(w1_24b))
        - cos(declin) * sin(slope) * sin(revasp) * (cos(w2_24) - cos(w1_24b)); // Eqn 5 - integrate between w1_24b and w2_24

        doublesunset=true;
        cos_theta = cos_theta1 + cos_theta2;  //Eqn 51
        //day_length = w2_24 - w1_24b + w2_24b - w1_24;
      }
    }
  }// end double sunset condition

  if(!doublesunset){
    //change integration limits to start & end of timestep, if needed
    if(!avg_daily){
      w1_24=max(w1_24,t1*2.0*PI);
      w2_24=min(w2_24,t2*2.0*PI);
      if(w2_24<=w1_24){ w2_24=w1_24; } //nighttime
    }

    cos_theta = (sindeclin * sinlatrad * cosslope - sindeclin * coslatrad * sinslope * cosrevasp) * (w2_24 - w1_24)
              +( cosdeclin * coslatrad * cosslope + cosdeclin * sinlatrad * sinslope * cosrevasp) * (sin(w2_24) - sin(w1_24))
              -( cosdeclin * sinslope  * sinrevasp                                              ) * (cos(w2_24) - cos(w1_24));  //Eqn 5

    //day_length=w2_24-w1_24;
  }

  return max(0.0,SOLAR_CONSTANT*cos_theta*ecc)/2.0/PI/(t2-t1);

}
//////////////////////////////////////////////////////////////////
/// \brief Calculates scattering transmissivity
/// \ref Bolsenga 1964 (from Dingman E-10-E-13, tau_sa) \cite Bolsenga1964
/// \param dew_pt [in] Dew point temperature [C]
/// \param opt_air_mass [in] Optical air mass [-]
/// \return Transmissivity from scattering
//
double CRadiation::CalcScatteringTransmissivity(const double dew_pt,//dew point temp, [C]
                                                const double opt_air_mass)//[-]
{
  double precip_WV;                     //precipitable water vaport content [cm]
  double a,b;
  precip_WV=1.12*exp(0.0614*dew_pt);    //Dingman E-10 (Bolsenga 1964)
  a=-0.124-0.0207*precip_WV;            /// \ref Dingman E-12
  b=-0.0682-0.0248*precip_WV;           /// \ref Dingman E-13
  return exp(a+b*opt_air_mass);         /// \ref Dingman E-11
}

/*****************************************************************
   Calculate Diffuse Scattering Transmissivity
------------------------------------------------------------------
/// returns the transmissivity due to scattering
/// Bolsenga 1964 (from Dingman E-16-E-18, tau_s) \cite Bolsenga1964
*****************************************************************/
double CRadiation::CalcDiffScatteringTransmissivity(const double dew_pt,//dew point temp, [C]
                                                    const double opt_air_mass)//[-]
{
  double precip_WV;                     //precipitable water vapor content [cm]
  double a,b;
  precip_WV=1.12*exp(0.0614*dew_pt);    //Dingman E-10 (Bolsenga 1964)
  a=-0.0363-0.0084*precip_WV;           //Dingman E-17
  b=-0.0572-0.0173*precip_WV;           //Dingman E-18
  return exp(a+b*opt_air_mass);         //Dingman E-16
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates daily optical air mass
/// \details returns the daily optical air mass integral [d]
/// \ref from X. Yin, Metorol. Atmos. Phys. (63) 1997 \cite Yin1997MaAP
/// \docminor some documentation needed here
///
/// \param latrad [in] Latitude [rad]
/// \param declin [in] Declination [rad]
/// \param cos_omt [in] ??
/// \param c1 [in] ??
/// \param c2 [in] ??
/// \param c3 [in] ??
/// \param t [in] Time of day w.r.t. solar noon
/// \return Daily optical air mass
//
double CRadiation::MassFunct(const double latrad, const double declin, const double cos_omt,
                             const double c1,const double c2, const double c3,const double t)
{
  /// \ref eqn 2.6b X. Yin
  double a,b,tmp;
  a=c1+sin(latrad)*sin(declin);
  b=   cos(latrad)*cos(declin);

  if (a>b)
  {
    tmp=(b+a*cos_omt)/(a+b*cos_omt);
    return 1.0/EARTH_ANG_VEL*c2/sqrt(a*a-b*b)*acos(tmp)+c3*t;
  }
  else if (a<b)
  {
    tmp = sqrt(b+a)*sqrt(1.0+cos_omt)+sqrt(b-a)*sqrt(1.0-cos_omt);
    tmp/=(sqrt(b+a)*sqrt(1.0+cos_omt)-sqrt(b-a)*sqrt(1.0-cos_omt));
    return 1.0/EARTH_ANG_VEL*c2/sqrt(b*b-a*a)*log(tmp)+c3*t;
  }
  else
  {
    return 1.0/EARTH_ANG_VEL*c2/a*tan(0.5*acos(cos_omt))+c3*t;
  }
}

//////////////////////////////////////////////////////////////////////
/// \brief Calculates optical air mass [-]
/// \param latrad [in] Latitude [rad]
/// \param declin [in] Declination [rad]
/// \param day_length [in] Day length [days]
/// \param t_sol [in] time of day w.r.t. solar noon [-0.5..0.5]
/// \param avg_daily [in] True if daily average is to be calculated
/// \return Optical air mass [-] \cite linacre1992
//
double CRadiation::OpticalAirMass(const double latrad,//in radians
                                  const double declin,//in radians
                                  const double day_length, //[days]
                                  const double t_sol, //time of day w.r.t. solar noon [days] [-0.5..0.5]
                                  const bool   avg_daily)//in radians
{
  double t1;                                                            //halfday length, days
  double Zn,Z0,t80;

  const double c1a= 0.00830700; //fitting constants
  const double c2a= 1.02129227;
  const double c3a=-0.01287829;
  const double c1b= 0.03716000;
  const double c2b= 1.53832220;
  const double c3b=-1.69726057;
  const double m90=39.7;        //air mass @90deg

  static const double deg80=80.0/180.0*PI;
  static const double cos80=cos(deg80);

  t1=0.5*day_length;
  if (t1==0.0){return m90;}

  Zn=acos(sin(latrad)*sin(declin)-cos(latrad)*cos(declin));
  if (fabs(latrad+declin)<=0.5*PI){Zn=0.5*PI;} //eqn 2.9

  if (avg_daily)
  {
    Z0=0.5*PI;
    if (fabs(latrad-declin)<=0.5*PI){Z0=latrad-declin;}//eqn 2.10

    t80=(1.0/EARTH_ANG_VEL)*acos((cos80-sin(latrad)*sin(declin))/(cos(latrad)*cos(declin)));//eqn. 2.8

    double cos_omt=cos(EARTH_ANG_VEL*t1);
    if (Zn<=deg80)
    {
      return MassFunct(latrad,declin,cos_omt,c1a,c2a,c3a,t1)/t1;
    }
    else if (Z0>=deg80)
    {
      return MassFunct(latrad,declin,cos_omt,c1b,c2b,c3b,t1)/t1;
    }
    else
    {
      double cos_omt80=cos(EARTH_ANG_VEL*t80);
      double sum=0;
      sum+=MassFunct(latrad,declin,cos_omt80,c1a,c2a,c3a,t80)/t1;
      sum+=MassFunct(latrad,declin,cos_omt  ,c1b,c2b,c3b, t1)/t1;
      sum-=MassFunct(latrad,declin,cos_omt80,c1b,c2b,c3b,t80)/t1;
      return sum;
    }
  }
  else
  {
    double cosZ=sin(latrad)*sin(declin)+cos(latrad)*cos(declin)*cos(EARTH_ANG_VEL*t_sol);
    if (cosZ<0){return m90;}
    if (cosZ>cos80){return c2a/(c1a+cosZ)+c3a;}
    else                                         {return c2b/(c1b+cosZ)+c3b;}
  }

  //elevation correction:
  //return m*(-elev/8000.0);//elev is masl (Linacre: Climate, Data, and Resources, 1992)
}

/////////////////////////////////////////////////////////////////
/// \brief Calculates total incident radiation, in [MK/m2/d] using Dingman (2002) over time interval tstep
/// \details Returns clear sky solar radiation
/// \remark must be corrected with albedo to obtain net radiation
/// \remark Kcs in Dingman text
/// \param julian_day [in] Julian representation of a day in a year
/// \param tstep [in] time step, in days
/// \param latrad [in] Latitude in radians
/// \param lateq [in] Equivalent latitude for slope
/// \param slope [in] Slope [rad]
/// \param aspect [in] aspect [rad from north]
/// \param day_angle [in] Day angle [rad]
/// \param day_length [in] Day length [days]
/// \param solar_noon [in] solar noon correction
/// \param dew_pt [in] Dew point temperature [C]
/// \param ET_radia [out] ET radiation without atmospheric corrections
/// \param avg_daily [in] True if average daily total incident radiation is to be computed instead
/// \return Double clear sky solar radiation
//
double CRadiation::ClearSkySolarRadiation(const double &julian_day,
                                          const double &tstep,
                                          const double &latrad,    //[rad]
                                          const double &lateq,     //[rad]
                                          const double &slope,     //[rad]
                                          const double &aspect,    //rad]
                                          const double &day_angle,
                                          const double &day_length,
                                          const double &solar_noon,//[days]
                                          const double &dew_pt,    //dew point temp, [C]
                                                double &ET_radia, //Extraterrestrial radiation on slope[MJ/m2/d]
                                                double &ET_radia_flat, //ET radiation on flat ground [MJ/m2/d]
                                          const bool    avg_daily) //true if average daily is to be computed
{
  double Ket,Ketp;          //daily solar radiation with (Ket) and without (Ketp) slope correction, [MJ/m2/d]
  double tau;               //total atmospheric transmissivity [-]
  double tau2;              //50% of solar beam attenuation from vapor and dust [-]
  double Mopt;              //optical air mass [-]
  double gamma_dust=0.025;  //attenuation due to dust
  double declin;            //solar declination
  double ecc;               //eccentricity correction [-]
  double t_sol,t_sol2;      //time of day w.r.t. solar noon (start & end of timestep) [d]

  declin= SolarDeclination(day_angle);
  ecc   = EccentricityCorr(day_angle);

  t_sol=julian_day-floor(julian_day)-0.5;
  t_sol2=t_sol+tstep;

  Mopt =OpticalAirMass(latrad,declin,day_length,t_sol,avg_daily);
  tau  =CalcScatteringTransmissivity(dew_pt,Mopt)-gamma_dust;               //Dingman E-9
  tau2 =0.5*(1.0-CalcDiffScatteringTransmissivity(dew_pt,Mopt)+gamma_dust); //Dingman E-15

  //Ketp =CalcETRadiation(latrad,lateq ,declin,ecc,slope,solar_noon,day_length,t_sol,avg_daily); //old dingman approach
  //Ket  =CalcETRadiation(latrad,latrad,declin,ecc,0.0  ,0.0,       day_length,t_sol,avg_daily);

  Ketp =CalcETRadiation2(latrad,aspect,declin,ecc,slope,t_sol,t_sol2,avg_daily);
  Ket  =CalcETRadiation2(latrad,aspect,declin,ecc,0.0  ,t_sol,t_sol2,avg_daily);

  ET_radia=Ketp;
  ET_radia_flat=Ket;

  //DingmanE-26 (from E-8,E-14,E-19,E-20)
  return tau *Ketp+                               //direct solar radiation on surface
         tau2*Ket+                                //diffuse radiation
         tau2*GLOBAL_ALBEDO*(tau2+tau)*Ketp;      //backscattered radiation
}


//////////////////////////////////////////////////////////////////
/// \brief Returns the clear sky total radiation [MJ/m2/d]
/// \details Returns the clear sky total solar radiation (total incident solar (shortwave)
///  radiation), on a horizontal plane in MJ/m2/d using the UBC Watershed Model approach
///  No eccentricity correction used, simple approximation of air mass effects
/// \remark Adapted from UBC Watershed model source code,(c) Michael Quick
///
/// \param latrad [in] Latitude [rad]
/// \param elev [in] Elevation [masl]
/// \param day_angle [in] Day angle [rad]
/// \param day_length [in] Day length [days]
/// \param ET_radia [out] ET radiation without atmospheric corrections
/// \return Clear sky total radiation [MJ/m2/d]
//
double CRadiation::UBC_SolarRadiation(const double latrad,   //[rad]
                                      const double elev,   //[masl]
                                      const double day_angle,
                                      const double day_length,//[d]
                                      double &ET_radia )
{
  double declin;        //declination in radians
  double daylength;     //daylength [d]
  double solar_rad;     //solar radiation [MJ/m2/d]
  double horizon_corr;  //horizon [deg]
  double air_mass,A1,Io;
  double cosZ;
  double turbidity;     //atmospheric turbidity
  int nIncrements=10;

  horizon_corr = 20.0;//P0BAHT;
  turbidity = 2.0;    //P0CLAR;

  declin   =SolarDeclination(day_angle);

  //mountain barrier correction
  daylength = day_length - (horizon_corr/360)/cos(latrad);

  //integrate solar radiation over day
  solar_rad=0.0;
  ET_radia =0.0;
  double dt=daylength/nIncrements; //[d]
  for (double t=(-0.5*daylength)+dt/2;t<(0.5*daylength);t+=dt)
  {
    cosZ=sin(latrad)*sin(declin)+cos(latrad)*cos(declin)*cos(EARTH_ANG_VEL*t);
    Io = SOLAR_CONSTANT*cosZ ;       //ET radiation outside earths atmosphere
    if (Io>0.0)
    {
      air_mass =1.0/cosZ;                //relative thickness of air mass at sea level
      air_mass*=(1.0 - 0.0001*elev);     //adjustment due to altitude
      A1 = 0.128 - 0.023452*log(air_mass); //molecular scattering coefficient
      if (A1 > 0) {
        solar_rad+=Io*exp(-A1*turbidity*air_mass)*dt;
      }
      ET_radia+=Io*dt;
    }
  }

  return solar_rad;
}
