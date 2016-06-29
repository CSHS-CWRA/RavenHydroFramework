/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright © 2008-2014 the Raven Development Team
----------------------------------------------------------------*/
#include "Properties.h"
#include "GlobalParams.h"

//////////////////////////////////////////////////////////////////
/// \brief Returns thermal conductivity of snow [W/m/K]
/// \ref from Jordan (1991)/CLM 3.0 eqn 6.65 \cite Jordan1991
///
/// \param &sno_dens [in] Snow density [kg/m^3]
/// \return Thermal conductivity of snow [W/m/K]
//
double GetSnowThermCond(const double &sno_dens)
{
	const double A1=7.75e-15;
  const double A2=1.105e-6;
	return TC_AIR+(A1*sno_dens+A2*pow(sno_dens,2.0))*(TC_ICE-TC_AIR);
}

//////////////////////////////////////////////////////////////////
/// \brief Estimates fraction of precipitation that is snow based temperature over timestep
/// \ref Brook90 approximation is really only valid for daily time step \cite Dingman1994
///
/// \param method [in] Method of converting total precipitation to rain/snow
/// \param *F [in & out] Forcing functions for a specific HRU
/// \param &Options [in] Global model options information
/// \return Double value for fraction of precipitation which is snow
//
double EstimateSnowFraction( const rainsnow_method method,
                             const force_struct       *F,
                             const optStruct       &Options)
{
	
  
  //-----------------------------------------------------------
  if (method==RAINSNOW_DATA)
  {
    return F->snow_frac;
  }
  //-----------------------------------------------------------
  else if (method==RAINSNOW_DINGMAN)
  { //from Brook90 model, Dingman pg 109
    double temp =CGlobalParams::GetParams()->rainsnow_temp;
    if (F->temp_daily_max<=temp){return 1.0;}
    if (F->temp_daily_min>=temp){return 0.0;}
	  return (temp-F->temp_daily_min)/(F->temp_daily_max-F->temp_daily_min);
  }
  //-----------------------------------------------------------
  else if ((method==RAINSNOW_HBV) || (method==RAINSNOW_UBCWM))
  {//linear variation based upon daily average temperature
    double frac;
    double delta=CGlobalParams::GetParams()->rainsnow_delta;
    double temp =CGlobalParams::GetParams()->rainsnow_temp;
    if      (F->temp_daily_ave<=(temp-0.5*delta)){frac=1.0;}
    else if (F->temp_daily_ave>=(temp+0.5*delta)){frac=0.0;}//assumes only daily avg. temp is included
    else    {
      frac= 0.5+(temp-F->temp_daily_ave)/delta;
    }
    if      (method==RAINSNOW_UBCWM){return frac;}
    //HBV-EC implementation - correction only applied to rain portion of snow
    else if (method==RAINSNOW_HBV){return frac*(1.0-F->snow_frac)+F->snow_frac;}
  }
  //-----------------------------------------------------------
  else if (method==RAINSNOW_HSPF)
  {
    double temp =CGlobalParams::GetParams()->rainsnow_temp;
    double snowtemp;
    double dewpt=GetDewPointTemp(F->temp_ave,F->rel_humidity);

    snowtemp=temp+(F->temp_ave-dewpt)*(0.3072+0.0144*F->temp_ave);

    if (F->temp_ave<snowtemp){return 1.0;}
    else                     {return 0.0;}
  }
  return 0.0;
	//Should check/add:
	//Stefan W. Kienzle,A new temperature based method to separate rain and snow, 
	///< Hydrological Processes 22(26),p5067-5085,2008,http://dx.doi.org/10.1002/hyp.7131 \cite kienzle2008HP
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates fresh snow density [kg/m^3] based on air temperature
/// 
/// \param &air_temp [in] Air temperature
/// \return Double value for fresh snow desnsity [kg/m^3]
//
double CalcFreshSnowDensity(const double &air_temp)
{
  //return FRESH_SNOW_DENS;
  
	// \ref WATCLASS 
  if (air_temp<=FREEZING_TEMP) {
    return 67.92+51.25*exp((air_temp-FREEZING_TEMP)/2.59);
  }
	else{
    return min((119.17+20.0*(air_temp-FREEZING_TEMP)),200.0);
  } 
	///< or (Anderson, 1976 -CLM manual eqn 7.18) \cite anderson1976
	/*if (air_temp>FREEZING_TEMP+2){
		return 169.15;
	}
	else if (air_temp<FREEZING_TEMP-15){
		return 50;
	}
	else{
		return 50+1.7*pow(air_temp-FREEZING_TEMP+15,1.5)
	}*/
}

/////////////////////////////////////////////////////////////////////
/// \brief Calculates and returns snow density [kg/m^3]
/// \remark From Dingman eqn 5-12
///
/// \param &snowSWE [in] Snow water equivalent [mm]
/// \param &snow_depth [in] Snow depth [mm]
/// \return Snow density [kg/m^3]
//
double GetSnowDensity(const double &snowSWE,const double &snow_depth)
{
	return (snowSWE/snow_depth)*DENSITY_WATER; // [kg m^-3]
}

/////////////////////////////////////////////////////////////////////
/// \brief Calculates snow depth from SWE and density [mm]
/// \param &snowSWE [in] Snow water equivalent [mm]
/// \param &snow_density [in] Density of snow [kg/m^3]
/// \return Snow depth [mm]
//
double GetSnowDepth(const double &snowSWE,//[mm]
                    const double &snow_density)//[kg/m^3]
{
	return (snowSWE/snow_density)*DENSITY_WATER; //[mm]
}

//////////////////////////////////////////////////////////////////////
/// \brief Calculates snow liquid holding capacity
/// \param &SWE [in] Snow water equivalent [mm]
/// \param &snow_depth [in] Depth of snow [mm]
/// \param &Options [in] Global model options information
/// \return Snow liquid capacity [-]
//
double CalculateSnowLiquidCapacity(const double &SWE,const double &snow_depth, const optStruct &Options)
{
  double liq_cap;
  //if (snow_depth>0.0){
  //  liq_cap=CGlobalParams::GetParams()->snow_SWI*(1.0-SWE/snow_depth);
  //}
  
  liq_cap= CGlobalParams::GetParams()->snow_SWI*SWE; //HBV-EC, Brook90, UBCWM, GAWSER

  return liq_cap;
  
  //from Dingman eqn 5-15
  //JRC: returns a negative value if snow_density<275 kg/m3 and 
  //doesnt seem to make sense - as snow becomes denser, 
  //holding capacity should reduce!
//  double snow_density=GetSnowDensity(SWE,snow_depth);
//	return -0.0735*(snow_density/DENSITY_WATER) + 
//         2.67e-4*(pow(snow_density,2)/DENSITY_WATER); //unitless
}


//////////////////////////////////////////////////////////////////
/// \brief Calculates sensible heat exchange with atmosphere for snow
/// \ref from Dingman pg 197-8 \cite Dingman 1994
///
/// \param &air_temp [in] Air temperature [C]
/// \param &surf_temp [in] Surface temperature [C]
/// \param &V [in] Velocity [m/s]
/// \param ref_ht [in] Reference height [m]
/// \param &rough Roughness coefficient [m]
/// \return Sensible heat exchange with atmosphere from snow [MJ/m^2/d]
//
double GetSensibleHeatSnow(const double &air_temp, //[C]
                           const double &surf_temp,//[C]
                           const double &V, 
                           const double &ref_ht,   //[m]
                           const double &rough)    //[m]
{
	double temp_var = pow(VON_KARMAN/log(ref_ht/rough),2);

	return DENSITY_AIR*SPH_AIR*temp_var*V*(air_temp-surf_temp); //sensible heat [MJ/m2/d]
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates latent heat exchange with atmosphere for snow
/// \ref from Dingman pg 198 \cite Dingman1994
///
/// \param &P [in] Pressure [kPa]
/// \param &air_temp [in] Air temperature [C]
/// \param &surf_temp [in] Surface temperature [C]
/// \param &rel_humid [in] Relative humidity [0..1]
/// \param &V [in] Velocity [m/s]
/// \param ref_ht [in] Reference height [m]
/// \param &rough Roughness coefficient [m]
/// \return Latent heat exchange with atmosphere from snow [MJ/m^2/d]
//
double GetLatentHeatSnow(const double &P,
                         const double &air_temp,
                         const double &surf_temp,
                         const double &rel_humid,
                         const double &V,
                         const double &ref_ht,
                         const double &rough)
{	
	double numer,denom,temp_var,vap_pres,surf_pres,LE;
	
	vap_pres  = GetVaporPressure(air_temp,rel_humid);
	surf_pres = GetVaporPressure(surf_temp,rel_humid);

	numer = AIR_H20_MW_RAT * DENSITY_AIR * pow(VON_KARMAN,2);
	denom = P * pow(log(ref_ht/rough),2);

	temp_var = numer/denom;
	LE = LH_VAPOR*temp_var*V*(vap_pres-surf_pres); //latent heat [MJ/m2/d]

	return LE;//[MJ/m2/d]
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates heat input from rain
/// \details Calculates heat input from rain by first calculating vapor pressure
/// and dew point temperature, and then by selecting an equationa ccording to the surface temperature
/// \ref from Dingman pg 199-200 \cite Dingman1994
///
/// \param &surf_temp [in] Surface temperature [C]
/// \param &rain_rate [in] Rate at which rain is accumulating? [mm/d]
/// \param &rel_humid [in] Relative humidity [0..1]
/// \return Heat input from rain [MJ/m^2/d]
//
double GetRainHeatInput(const double &surf_temp, 
                        const double &rain_rate,
                        const double &rel_humid)
{
	double R=0,rain_temp,vap_pres;

	vap_pres = GetVaporPressure(surf_temp,rel_humid); //[kPa]
	rain_temp = GetDewPointTemp(vap_pres); //[C]

	if(surf_temp < FREEZING_TEMP)
	{
		R =	DENSITY_WATER*(rain_rate/MM_PER_METER)*(SPH_WATER*(rain_temp - FREEZING_TEMP) + LH_FUSION);
	}
	if(surf_temp == FREEZING_TEMP)
	{
		R =	DENSITY_WATER*(rain_rate/MM_PER_METER)*(SPH_WATER*(rain_temp - FREEZING_TEMP));
	}

	return R; // [MJ m^-2 d^-1]
}


