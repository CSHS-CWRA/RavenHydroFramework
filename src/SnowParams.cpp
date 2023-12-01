/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2023 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Properties.h"
#include "GlobalParams.h"
#include "HydroUnits.h"
#include "Model.h"

//////////////////////////////////////////////////////////////////
/// \brief Returns thermal conductivity of snow [MJ/d/m/K]
/// \ref from Jordan (1991)/CLM 3.0 eqn 6.65 \cite Jordan1991
/// \param &sno_dens [in] Snow density [kg/m^3]
/// \return Thermal conductivity of snow [MJ/d/m/K]
//
double GetSnowThermCond(const double &sno_dens)
{
  const double A1=7.75e-15;
  const double A2=1.105e-6;
  return WATT_TO_MJ_PER_D*(TC_AIR+(A1*sno_dens+A2*pow(sno_dens,2.0))*(TC_ICE-TC_AIR));

  //Sturm et al. 1997. The thermal conductivity of seasonal snow Journal of Glaciology, Vol. 43, No. 143
  //Snow conductivity used in CRHM
  /*
  double rho=sno_dens*GPCM3_PER_KGPM3;//density in g/cm3
  if(rho < 0.156){return WATT_TO_MJ_PER_D*(0.023 + 0.234*rho*rho);}
  else           {return WATT_TO_MJ_PER_D*(0.138-1.01*rho+3.233*rho*rho);}
  */
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

  // \ref WATCLASS , also USACE (1956) (from Hedstrom & Pomeroy)
  if (air_temp<=FREEZING_TEMP) {
    return 67.92+51.25*exp((air_temp-FREEZING_TEMP)/2.59); //range: 67.92-119.17
  }
  else{
    return min((119.17+20.0*(air_temp-FREEZING_TEMP)),200.0); //range: 119.18-200
  }

  //From Hydrotel 2.1 documentation, eqn 3.1
  //return max(50.0,min(151.0,151+10.63*air_temp+0.2767*air_temp*air_temp)); //range: 50 to 151

  ///< or (Anderson, 1976 -CLM manual eqn 7.18) \cite anderson1976
  /*if      (air_temp>+2.0 ){return 169.15;}
    else if (air_temp<-15.0){return 50.0;}
    else                    {return 50+1.7*pow(air_temp+15,1.5) //range: 50-169.15
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

//////////////////////////////////////////////////////////////////////
/// \brief Calculates snow liquid holding capacity
/// \param &SWE [in] Snow water equivalent [mm]
/// \param &snow_depth [in] Depth of snow [mm]
/// \param &Options [in] Global model options information
/// \return Snow liquid capacity [-]
//
double CalculateSnowLiquidCapacity(const double &SWE,const double &snow_depth, const CModelABC* pModel)
{
  double liq_cap;
  //if (snow_depth>0.0){
  //  liq_cap=CGlobalParams::GetParams()->snow_SWI*(1.0-SWE/snow_depth);
  //}

  liq_cap = pModel->GetGlobalParams()->GetParams()->snow_SWI*SWE; //HBV-EC, Brook90, UBCWM, GAWSER

  return liq_cap;

  //from Dingman eqn 5-15
  //JRC: returns a negative value if snow_density<275 kg/m3 and
  //doesnt seem to make sense - as snow becomes denser,
  //holding capacity should reduce!
//  double snow_density=GetSnowDensity(SWE,snow_depth);
//      return -0.0735*(snow_density/DENSITY_WATER) +
//         2.67e-4*(pow(snow_density,2)/DENSITY_WATER); //unitless
}


//////////////////////////////////////////////////////////////////
/// \brief Calculates sensible heat exchange with atmosphere for snow
/// \ref from Dingman pg 197-8 \cite Dingman 1994
///
/// \param &air_temp [in] Air temperature [C]
/// \param &surf_temp [in] Surface temperature [C]
/// \param &V [in] Velocity [m/s]
/// \param ref_ht [in] Measurement height [m]
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

  return DENSITY_AIR*SPH_AIR*temp_var*V*SEC_PER_DAY*(air_temp-surf_temp); //sensible heat [MJ/m2/d]

  //        [kg/m3]  [MJ/kg/K] [-] [m/s] [s/d]  [K] -> [MJ/m2/d]
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates latent heat exchange with atmosphere for snow
/// \ref from Dingman pg 198 \cite Dingman1994 (consistent with SUBLIM_BULK_AERO)
///
/// \param &P [in] Pressure [kPa]
/// \param &air_temp [in] Air temperature [C]
/// \param &surf_temp [in] Surface temperature [C]
/// \param &rel_humid [in] Relative humidity [0..1]
/// \param &V [in] wind speed [m/s]
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
  double temp_var,vap_pres,surf_pres,LH,CE;

  vap_pres  = GetSaturatedVaporPressure(air_temp )*rel_humid;
  surf_pres = GetSaturatedVaporPressure(surf_temp)*1.0; //assume saturated

  CE = pow(VON_KARMAN,2)/pow(log(ref_ht/rough),2);

  temp_var = AIR_H20_MW_RAT * DENSITY_AIR * CE / P;

  LH = LH_SUBLIM*temp_var*V*SEC_PER_DAY*(vap_pres-surf_pres); //latent heat [MJ/m2/d]

  return LH;//[MJ/m2/d]
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates heat input from rain
/// \details Calculates heat input from rain by first calculating vapor pressure
/// and dew point temperature, and then by selecting an equationa ccording to the surface temperature
/// \ref from Dingman pg 199-200 \cite Dingman1994
///
/// \param &surf_temp [in] Surface temperature [C]
/// \param &air_temp [in] air temperature [C]
/// \param &rainfall [in] rainfall rate [mm/d]
/// \param &rel_humid [in] Relative humidity [0..1]
/// \return Heat input from rain [MJ/m^2/d]
//
double GetRainHeatInput(const double &surf_temp,
                        const double &air_temp,
                        const double &rainfall,
                        const double &rel_humid)
{
  double R=0,rain_temp,vap_pres;

  vap_pres = GetSaturatedVaporPressure(air_temp)*rel_humid; //[kPa]
  rain_temp = GetDewPointTemp(vap_pres); //[C]

  if(surf_temp < FREEZING_TEMP)
  {
    R = DENSITY_WATER*(rainfall/MM_PER_METER)*(SPH_WATER*(rain_temp - FREEZING_TEMP) + LH_FUSION);
  }
  if(surf_temp >= FREEZING_TEMP)
  {
    R = DENSITY_WATER*(rainfall/MM_PER_METER)*(SPH_WATER*(rain_temp - FREEZING_TEMP));
  }
  return R; // [MJ/m2/d]
}
