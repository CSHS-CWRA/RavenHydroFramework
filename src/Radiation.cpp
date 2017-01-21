/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright © 2008-2014 the Raven Development Team

  Routines for calculating radiation
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Radiation.h"
#include "GlobalParams.h"

//////////////////////////////////////////////////////////////////
/// \brief Estimates net shortwave radiation [MJ/m2/d]
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
double CRadiation::EstimateShortwaveRadiation(const optStruct    &Options,
                                              const force_struct *F,
                                              const CHydroUnit   *pHRU,
                                              const time_struct  &tt,
                                                    double       &ET_rad)	
{
  double latrad=pHRU->GetLatRad();
  double lateq =pHRU->GetLatEq();
  switch(Options.SW_radiation)
  {
    //--------------------------------------------------------
    case(SW_RAD_DATA):
		{
      return F->SW_radia;
      break;
    }
    //--------------------------------------------------------
    case(SW_RAD_DEFAULT)://Dingman
    {
      double dew_pt;
      dew_pt=GetDewPointTemp(F->temp_ave,F->rel_humidity);
      
      double slope=pHRU->GetSlope();
      double solar_noon= pHRU->GetSolarNoon();
      return ClearSkySolarRadiation(tt.julian_day,latrad,lateq,slope,F->day_angle,F->day_length,solar_noon,dew_pt,ET_rad,(Options.timestep>=1));
      break;
    }
    //--------------------------------------------------------
    case(SW_RAD_UBCWM):
    {  
      if (tt.day_changed) //no need to do calculations every timestep
      {
				double solar_rad;
				double elev=pHRU->GetElevation();
				solar_rad= UBC_SolarRadiation(latrad,elev,F->day_angle, F->day_length,ET_rad);

				double orient=1.0-fabs(pHRU->GetAspect()/PI-1.0);        //=0 for north, 1.0 for south 
				double shortwave_corr_S=InterpolateMo(CGlobalParams::GetParams()->UBC_s_corr,tt,Options);
				double shortwave_corr_N=InterpolateMo(CGlobalParams::GetParams()->UBC_n_corr,tt,Options);
				double shortwave_corr=((orient)* shortwave_corr_S + (1.0-orient)*shortwave_corr_N);

				return solar_rad * shortwave_corr;
      }
      else
      {
        return F->SW_radia_unc;
      }

      break;
    }
    //--------------------------------------------------------
    case(SW_RAD_VALIANTZAS) :
    {
      double declin = SolarDeclination(F->day_angle);
      double ws = acos(-tan(latrad)*tan(declin));
      double ecc = EccentricityCorr(F->day_angle);
      ET_rad = 37.59*ecc*(ws*sin(latrad)*sin(declin) + sin(ws)*cos(latrad)*cos(declin)); //Extraterrestrial radiation
      return  min(0.75+0.00002*pHRU->GetElevation(),1.0)*ET_rad;
      break;
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
double CRadiation::EstimateLongwaveRadiation(const optStruct     &Options,
                                             const force_struct  *F,
                                             const CHydroUnit    *pHRU)
{
  switch(Options.LW_radiation)
  {
    //--------------------------------------------------------
    case(LW_RAD_DATA):
		{
      return F->LW_radia;
      break;
    }
    //--------------------------------------------------------
    case(LW_RAD_DEFAULT):
    {
      //from Dingman eqns. 5-40
      double emissivity=0.99; // \todo [funct]: remove hardcoded emissivity
      double forest_cover=pHRU->GetSurfaceProps()->forest_coverage;
      double Tair =F->temp_ave+ZERO_CELSIUS;  //[K]
      //double Tsurf=pHRU->GetSurfaceTemperature()+ZERO_CELSIUS;//[K] (not yet functional)
      double Tsurf=F->temp_ave+ZERO_CELSIUS;  //[K] //\todo better way to do this.

      double ea=GetSaturatedVaporPressure(F->temp_ave);

      double emiss_eff; //effective clear sky emmisivity
      emiss_eff=1.72*pow(ea/Tair,1/7.0);/// \ref Kustas, 1994, via Dingman clear sky emmissivity \cite Moran1994WRR

      //emiss_eff = 1.24*pow(10.0*ea/Tair,1/7.0);                ///< Brutsaert, 1982, via Brook90 documentation \cite Brutsaert1982
 	    //emiss_eff = 1.08 *(1.0- exp(-(10*ea)*Tair/2016.0));      ///< Satterlund (1979) via Brook90 documentation \cite Satterlund1979WRR
 	    //emiss_eff = 1.0- 0.261*exp(-0.000777*F->temp_ave*F->temp_ave); ///< Idso and Jackson (1969) via Brook90 documentation \cite Idso1969JoGR
 	    //emiss_eff = 0.0000092*Tair*Tair;                         ///< Swinbank (1963) via Brook90 documentation \cite Swinbank1963QJotRMS

      double eps_at=emiss_eff*(1.0+0.22*F->cloud_cover*F->cloud_cover);

      eps_at=(1-forest_cover)*eps_at+(forest_cover)*1.0; //treats forest as blackbody

      return STEFAN_BOLTZ*emissivity*(eps_at*pow(Tair,4)-pow(Tsurf,4));
      break;
    }
    //--------------------------------------------------------
    case (LW_RAD_UBCWM):
    {
      double cloud=F->cloud_cover;
      double LW_open=(1-cloud)*(-20+0.94*F->temp_daily_ave)+(cloud)*(1.24*F->temp_daily_min);

      double P0BLUE=CGlobalParams::GetParams()->UBC_LW_forest_fact;//actually P0BLUE * P0LWVF , [mm/d/K]
      double LW_forest=P0BLUE*F->temp_daily_ave;

      double tmp=(LH_FUSION*DENSITY_WATER/MM_PER_METER); //[mm/d]-->[MJ/m2/d]

      double Fc=pHRU->GetSurfaceProps()->forest_coverage;

      return tmp*((1.0-Fc)*LW_open+(Fc)*LW_forest);

      break;
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
      break;
    }
    //--------------------------------------------------------
    case(LW_RAD_VALIANTZAS):
    {
      double f;
      double eps;
      double sat_vap,ea;
      sat_vap=GetSaturatedVaporPressure(F->temp_ave);
      ea =F->rel_humidity*sat_vap;
      eps= 0.34 - 0.14*sqrt(ea);
      f=(1.35*(F->SW_radia/F->SW_radia_unc)-0.35);
      return f*eps*STEFAN_BOLTZ*pow(F->temp_ave+ZERO_CELSIUS,4.0);
      break;
    }
    /*case (LW_RAD_UNFAO):
    {
      ///< from Crop evapotranspiration - Guidelines for computing crop water requirements - FAO Irrigation and drainage paper 56 \cite Allen1998FR
      //http://www.fao.org/docrep/X0490E/X0490E00.htm
      double tmp,ea;
      ea=F->rel_humidity*GetSaturatedVaporPressure(F->temp_ave); //actual vapor pressure
      tmp=(0.34-0.14*pow(ea,0.5))*(1.35*(F->SW_radia/F->CS_radia)-0.35);
      return STEFAN_BOLTZ*tmp*(pow(F->temp_daily_min,4)+pow(F->temp_daily_max,4));
    }*/
  /*
    case(LW_RAD_BRUTSAERT):
    {//Brook90-Brutsaert 
      const double C1=0.25;
      const double C2=0.5;
      const double C3=0.2;
      //Brutsaert (1982) equation for effective clear sky emissivity  
      double cloud_corr; //cloud correction
      double ratio=solar_radiation/potential_solar_rad;//if solar_radiation is specified
      cloud_corr=C3+(1.0-C3)*min(max((ratio - C1) / C2,1.0),0.0);
      return STEFAN_BOLTZ*emmissivity*cloud_corr*(emiss_eff*pow(Tair,4)-pow(Tair,4)); 
    }
  */
    /*case(LW_RAD_SWAT):
    {
      // calculate net long-wave radiation
      double rbo,rout;
      double incoming =F->SW_radia;
      double rmx=??;
	    rbo = -(0.34 - 0.139 * sqrt(sat_vap*F->rel_humidity));			  // net emissivity  equation 2.2.20 in SWAT manual - SHOULD CHECK - is RH 0-100 in SWAT???
	    if (rmx < 1e-4){rto = 0.0;}
	    else           {rto = 0.9 * (F->SW_radia / rmx) + 0.1;}				// cloud cover factor equation 2.2.19 SWAT
	    rout= rbo * rto * 4.9e-9 * pow(F->temp_ave+ZERO_CELSIUS,4);	  // net long-wave radiation equation 2.2.21 SWAT
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
/// \param *pHRU [in] pointer to HRU for which radiation is calculated
/// \return Double Cloud cover correction factor for shortwave radiation

double CRadiation::SWCloudCoverCorrection(const optStruct    &Options,
							                            const force_struct *F)
{
	switch(Options.SW_cloudcovercorr)
	{
	  case(SW_CLOUD_CORR_UBCWM):
	  {
		  double POCAST = CGlobalParams::GetParams()->UBC_cloud_penet;
		  return 1.0 -(1.0 - POCAST) * F->cloud_cover;
		  break;
	  }
	  case(SW_CLOUD_CORR_DINGMAN): // Dingman (2008) Eq. 5-30 (does not require parameters)
		  {
			  return 0.355 + 0.68 * (1 - F->cloud_cover);
			  break;
		  }
	  default: // apply no cloud cover correction (default)
		{
			return 1.0;
			break;
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

double CRadiation::SWCanopyCorrection(const optStruct  &Options, 
									                    const CHydroUnit *pHRU)
{
  switch(Options.SW_canopycorr)
  {
    case(SW_CANOPY_CORR_STATIC):
	  {
		  double LAI = pHRU->GetVegVarProps()->LAI;
		  double SAI = pHRU->GetVegVarProps()->SAI;
		  double extinction = pHRU->GetVegetationProps()->svf_extinction;
		
		  return exp(-extinction*(LAI+SAI)); // Dingman (2008) Eq. 5-33
		  break;
	  }
	  case(SW_CANOPY_CORR_UBCWM):  // A simple factor that switches on when forest cover is greater than zero
	  {
		  double shortwave_corr = 1.0;
		  double Fc = pHRU->GetSurfaceProps()->forest_coverage; // forest cover
		  double UBC_correction_fact = CGlobalParams::GetParams()->UBC_exposure_fact;	//Sun exposure factor of forested areas
		  shortwave_corr*=((Fc)*UBC_correction_fact+(1.0-Fc)*1.0);
		  return shortwave_corr;
		  break;
	  }
	  case(SW_CANOPY_CORR_DYNAMIC):  // stub for Jost and Moore (2010) equations
	  {
		  ExitGracefully("SWCanopyCorrection::SW_CANOPY_CORR_DYNAMIC",STUB);
	  }
	  case(SW_CANOPY_CORR_NONE):  
	  {
		  return 1.0;
	  }
	  default: // apply no correction for canopy (default)
	  {
		  return 1.0;
		  break;
	  }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Calculates day angle [rad]
/// \param &julian_day [in] Julian identifier of a specfic day
/// \param year [in] Year in which this day occurs
/// \return The day angle of this day [rad]
//
double CRadiation::DayAngle(const double&julian_day,
														const int    year)
{
	double leap=0.0;
	if (IsLeapYear(year)){leap=1.0;}
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
	if			(hd>= 1.0) {return 0.0;}//sun stays below horizon
	else if (hd<=-1.0) {return 1.0;}//sun stays above horizon
	return acos(hd)/PI;
}

//////////////////////////////////////////////////////////////////
/// \brief LOCAL function that applies a correction to extraterrestrial radiation
/// \details Applies multiplication factor to solar constant and eccentricity factor for daily incoming radiation
/// \ref from Dingman E-25 or E-6/E-7 \cite Dingman1994
/// \return Correction factor for shortwave radiation
///
double CRadiation::RadCorr (const double declin,		 //declination, radians
								            const double solar_noon, //time of solar noon, days
								            const double lat_eq,     //equivalent latitude, radians
								            const double t_sunrise,  //time of sunrise, days before solar noon
								            const double t_sunset,   //time of sunrise, days after solar noon,
								            const double t_sol,			 //time of day w.r.t. solar noon, days (not adjusted for slope)
								            const bool   avg_daily)	 //in radians
{
	if (avg_daily)
	{
		//averaged daily correction
		return sin(declin)*sin(lat_eq)*(t_sunset-t_sunrise) + 
					 cos(declin)*cos(lat_eq)*(sin(EARTH_ANG_VEL*(t_sunset -solar_noon))-
					                          sin(EARTH_ANG_VEL*(t_sunrise-solar_noon)))/EARTH_ANG_VEL;
	}
	else if ((t_sol<=t_sunrise) || (t_sol>=t_sunset))//nighttime
	{
		return 0.0; 
	}
	else{//eqn E-6 Dingman
		return sin(declin)*sin(lat_eq) + 
			     cos(declin)*cos(lat_eq)*cos(EARTH_ANG_VEL*(t_sol-solar_noon));
	}
}

//////////////////////////////////////////////////////////////////
/// \brief Returns solar declination [rad]
/// \ref Uses Dingman eqn E-3 \cite Dingman1994
/// \param day_angle Day angle, in radians
/// \return Solar declination [rad]
//
double CRadiation::SolarDeclination(const double day_angle)
{
	/*return	0.006918-
					0.399912*cos(1.0*day_angle)+
					0.070257*sin(1.0*day_angle)-
					0.006758*cos(2.0*day_angle)+
					0.000907*sin(2.0*day_angle)-
					0.002697*cos(3.0*day_angle)+
					0.001480*sin(3.0*day_angle);*///Dingman eqn E-3
	//return asin(0.39785* sin(4.868961+0.017203*julian_day+
	//	                         0.033446*sin(6.224111+0.017202*julian_day)));//Brook90
	//return 0.4903*sin(day_angle-1.405); //WATFLOOD
  //return 0.409*sin(day_angle-1.39); //Valiantzas, 2006
  return 0.40928*sin(day_angle-1.384);//UBCWM 
}

//////////////////////////////////////////////////////////////////
/// \brief Returns solar zenith angle [rad]
/// \param lat_eq [in] Equivalent latitude [rad]
/// \param declin [in] Solar declination [rad]
/// \param t [in] Absolute time [d]
/// \return Solar zenith angle [rad]
//
double CRadiation::SolarZenithAngle(const double lat_eq,//equivalent latitude, radians 
                                    const double declin,//solar declination, radians
                                    const double t)     //absolute time, days
{
  double Z,H;
  H=PI*(t-0.5);//hour angle (=0 at noon, -PI..PI over day)
  Z = acos( sin(lat_eq)*sin(declin) + cos(lat_eq)*cos(declin)*cos(H));
  return Z;
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
/// \brief Calculates extraterrestrial radiation on an inclined plane
/// \details Returns total incoming solar radation, K_{ET}, in [MJ/m2/d], corrected for slope and aspect
/// \note Has been tested and is working for zero-aspect slope
/// \ref from Swift (1976), Lee (1964) via Dingman \cite Swift1976WRR \cite Lee1964IJoB
///
/// \param latrad [in] Latitude in radians
/// \param lateq [in] Equivalent latitude for slope (Dingman E-23)
/// \param declin [in] Solar declination in radians
/// \param ecc [in] Eccentricity correction
/// \param slope [in] Slope [rad]
/// \param aspect [in] Aspect [rad]
/// \param t_sol [in] Time of day with respecto solar noon [d]
/// \param avg_daily [in] True if average daily total incoming solar radiation is to be computed instead
/// \return Double clear sky total radiation
//
double CRadiation::CalcETRadiation( const double latrad,//latitude in radians
                                    const double lateq, //equivalent latitude for slope in radians
														        const double declin,//solar declination in radians
														        const double ecc,		//eccentricity correction [-]
														        const double slope, //in radians
                                    const double solar_noon, //effective solar noon correction for slope [days]
                                    const double day_length, // length of day (not corrected for slope) [days]
														        const double t_sol, //time of day w.r.t. solar noon [days] [-0.5..0.5]
														        const bool   avg_daily)//in radians
{
	double ExTerCorr;  	//Extraterrestrial radiation correction [-]	
	double hd;					//half day length [days]
	double t_rise,t_set;//time of sunrise, sunset (w.r.t. solar noon)[days] [-0.5..0.5]
	double lat_eq;			//slope-corrected equivalent latitude [rad]
	double t_rise1,t_set1;
	double t_rise2,t_set2;
	bool	 TwoDawns;

	//---slope corrections-------------------------------------------
	lat_eq =latrad;
	if (slope!=0.0){
		lat_eq = lateq;	//equivalent latitude for slope (Dingman E-23)
  }
		
	//calculate sunrise & sunset--------------------------------------------
	hd=0.5*day_length;
	t_rise=-hd;			//actual sunrise/set: Dingman E-5a,b
	t_set = hd;	
  
  //sunrise/set corrected for slope & aspect
	if(slope !=0){
    hd=0.5*DayLength(lat_eq, declin);}

	t_rise1=max(t_rise,solar_noon-hd);//Dingman E-24a,b
	t_set1 =min(t_set ,solar_noon+hd);
	if (t_rise1>t_set1) {t_rise1=t_set1=0.0;}//Correct for those endless nights...

	ExTerCorr=RadCorr(declin, solar_noon, lat_eq, t_rise1, t_set1,t_sol,avg_daily);

	//Corrections for two sunrises 
	TwoDawns = false;
	t_rise2=solar_noon+1.0-hd;
	t_set2 =solar_noon-1.0+hd;
	if      (t_rise2<t_set ){t_set2 =t_set; TwoDawns=true;}
	else if (t_set2 >t_rise){t_rise2=t_rise;TwoDawns=true;}
	if (TwoDawns){
		ExTerCorr+=RadCorr(declin, solar_noon, lat_eq, t_rise2,t_set2,t_sol,avg_daily);
	}
	//---end sunrise/sunset correction----------------------------------------

	return SOLAR_CONSTANT*ecc*ExTerCorr;//[MJ/m2/d] E-25 dingman

	//Watflood approach (slope/aspect insensitive):
	/*
	double hdrad=0.5*DayLength(latrad,declin)/EARTH_ANG_VEL; 
	double tmp  = 15.392*ecc*(hdrad*sin(latrad)*sin(declin)+cos(latrad)*cos(declin)*sin(hdrad));//[mm/d]
	//[MJ/m2/d] =[mm/d]*(J/kg)*(MJ/J)*(kg/m3)*(m/mm)
	return tmp*(LH_VAPOR*DENSITY_WATER/MM_PER_METER);//[MJ/m2/d]
	*/
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
	double precip_WV;			//precipitable water vaport content [cm]
	double a,b;
	precip_WV=1.12*exp(0.0614*dew_pt);				 //Dingman E-10 (Bolsenga 1964)
	a=-0.124-0.0207*precip_WV;								 /// \ref Dingman E-12
	b=-0.0682-0.0248*precip_WV;								 /// \ref Dingman E-13
	return exp(a+b*opt_air_mass);              /// \ref Dingman E-11
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
	double precip_WV;			//precipitable water vapor content [cm]
	double a,b;
	precip_WV=1.12*exp(0.0614*dew_pt);				 //Dingman E-10 (Bolsenga 1964)
	a=-0.0363-0.0084*precip_WV;								 //Dingman E-17
	b=-0.0572-0.0173*precip_WV;								 //Dingman E-18
	return exp(a+b*opt_air_mass);//Dingman E-16
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
  double t1;								//halfday length, days
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
		else					 {return c2b/(c1b+cosZ)+c3b;}
	}

  //elevation correction:
  //return m*(-elev/8000.0);//elev is masl (Linacre: Climate, Data, and Resources, 1992) 
}

/////////////////////////////////////////////////////////////////
/// \brief Calculates total incident radiation, in [MK/m2/d] using Dingman (2002)
/// \details Returns clear sky solar radiation 
/// \remark must be corrected with albedo to obtain net radiation
/// \remark Kcs in Dingman text
/// \param julian_day [in] Julian representation of a day in a year
/// \param latrad [in] Latitude in radians
/// \param lateq [in] Equivalent latitude for slope
/// \param slope [in] Slope [rad]
/// \param day_angle [in] Day angle [rad]
/// \param day_length [in] Day length [days]
/// \param solar_noon [in] solar noon correction
/// \param dew_pt [in] Dew point temperature [C]
/// \param ET_radia [out] ET radiation without atmospheric corrections
/// \param avg_daily [in] True if average daily total incident radiation is to be computed instead
/// \return Double clear sky solar radiation
//
double CRadiation::ClearSkySolarRadiation(const double julian_day,
											const double latrad,   //[rad]
											const double lateq, //[rad]
											const double slope, //[rad]
											const double day_angle,
											const double day_length,
											const double solar_noon,//[days]
											const double dew_pt,//dew point temp, [C]
														double &ET_radia, //ET radiation [MJ/m2/d]
													const bool   avg_daily)	//true if average daily is to be computed
	{
	double Ket,Ketp;				//daily solar radiation with (Ket) and without (Ketp) slope correction, [MJ/m2/d]
	double tau;							//total atmospheric transmissivity [-]
	double tau2;						//50% of solar beam attenuation from vapor and dust [-]
	double Mopt;						//optical air mass [-]
	double gamma_dust=0.025;//attenuation due to dust
	double declin;					//solar declination
	double ecc;							//eccentricity correction [-]
	double t_sol;					  //time of day w.r.t. solar noon [d]
	
	declin=	SolarDeclination(day_angle);
	ecc   = EccentricityCorr(day_angle);

	t_sol=julian_day-floor(julian_day)-0.5;
  Mopt =OpticalAirMass(latrad,declin,day_length,t_sol,avg_daily);
	tau  =CalcScatteringTransmissivity(dew_pt,Mopt)-gamma_dust;							  //Dingman E-9
	tau2 =0.5*(1.0-CalcDiffScatteringTransmissivity(dew_pt,Mopt)+gamma_dust); //Dingman E-15

	Ketp =CalcETRadiation(latrad,lateq,declin,ecc,slope,solar_noon,day_length,t_sol,avg_daily);
  Ket  =CalcETRadiation(latrad,lateq,declin,ecc,0.0  ,0.0, day_length,t_sol,avg_daily);
 
  ET_radia=Ketp;

  //DingmanE-26 (from E-8,E-14,E-19,E-20)
	return tau *Ketp+												 //direct solar radiation on surface
		     tau2*Ket+												 //diffuse radiation
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
  double declin; //declination in radians
  double daylength;     //daylength [d]
  double solar_rad;     //solar radiation [MJ/m2/d]
  double horizon_corr;  //horizon [deg] 
  double air_mass,A1,Io;
  double cosZ;
  double turbidity;    //atmospheric turbidity
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



