/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ----------------------------------------------------------------
  Sublimation
  ----------------------------------------------------------------*/
#include "HydroProcessABC.h"
#include "SnowMovers.h"
#include "Model.h"

/*****************************************************************
   Sublimation Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the sublimation constructor
/// \param sub_type [in] Sublimation-modelling technique selected
//
CmvSublimation::CmvSublimation(sublimation_type sub_type, CModelABC *pModel)
  :CHydroProcessABC(SUBLIMATION, pModel)
{
  type =sub_type;

  int iSnow,iAtmos,iSnowDepth;
  iSnow     =pModel->GetStateVarIndex(SNOW);
  iAtmos    =pModel->GetStateVarIndex(ATMOSPHERE);
  iSnowDepth=pModel->GetStateVarIndex(SNOW_DEPTH);

  ExitGracefullyIf(iSnow==DOESNT_EXIST,
                   "CmvSublimation: snow state variable required",BAD_DATA);

  model_depth_change=(iSnowDepth!=DOESNT_EXIST);

  if (model_depth_change){
    CHydroProcessABC::DynamicSpecifyConnections(2);//nConnections=2
    iFrom[0]=iSnow;      iTo[0]=iAtmos;     //rates[0]: SNOW->ATMOSPHERE
    iFrom[1]=iSnowDepth; iTo[1]=iSnowDepth; //rates[1]: SNOW_DEPTH change
  }
  else{
    CHydroProcessABC::DynamicSpecifyConnections(1);//nConnections=1
    iFrom[0]=iSnow;      iTo[0]=iAtmos;     //rates[0]: SNOW->ATMOSPHERE
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvSublimation::~CmvSublimation(){}

//////////////////////////////////////////////////////////////////
/// \brief Initializes sublimation modelling object
//
void CmvSublimation::Initialize(){}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for sublimation algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by sublimation algorithm (size of aP[] and aPC[])
//
void CmvSublimation::GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const
{
  nP=1;
  aP[0]="SNOW_TEMPERATURE"; aPC[0]=CLASS_GLOBAL;
  if (type==SUBLIM_SVERDRUP)
  {
    nP+=1;
    aP[1]="SNOW_ROUGHNESS";   aPC[1]=CLASS_GLOBAL;
  }
  else if (type==SUBLIM_PBSM)
  {
    nP=0;
    //algorithm not complete
  }
  else
  {
    nP=0; //most have no params
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
///
/// \param sub_type [in] Sublimation modelling type selected
/// \param *aSV [out] Array of state variable types needed by sublimation algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by sublimation algorithm (size of aSV[] and aLev[] arrays)
//
void CmvSublimation::GetParticipatingStateVarList(sublimation_type sub_type,sv_type *aSV, int *aLev, int &nSV)
{
  nSV=2;
  aSV[0]=SNOW;       aLev[0]=DOESNT_EXIST;
  aSV[1]=ATMOSPHERE; aLev[1]=DOESNT_EXIST;
  //might also effect SNOW_DEPTH at some point
}

double  SublimationRate  (const double      *state_vars,
                          const CHydroUnit  *pHRU,
                          const CModelABC   *pModel,
                          const time_struct &tt,
                          const double      &wind_vel,
                          sublimation_type   type)
{
  const optStruct* Options = pModel->GetOptStruct();
  double Ta      =pHRU->GetForcingFunctions()->temp_ave;
  double rel_hum =pHRU->GetForcingFunctions()->rel_humidity;

  //-----------------------------------------------------------------
  if(type==SUBLIM_SVERDRUP)
  {
    double roughness   = pModel->GetGlobalParams()->GetParams()->snow_roughness;  // [m]
    double air_density = pHRU->GetForcingFunctions()->air_dens;       // [kg/m3]
    double air_pres    = pHRU->GetForcingFunctions()->air_pres;       // [kPa]
    double Tsnow       = pHRU->GetSnowTemperature();

    double es          = GetSaturatedVaporPressure(Tsnow);     // [kPa]
    double ea          = GetSaturatedVaporPressure(Ta)*rel_hum;// [kPa]

    double C_E=VON_KARMAN*VON_KARMAN/(log(8.0/roughness))/(log(8.0/roughness));//ref ht =8 m
    double sublim=air_density*C_E*wind_vel *  0.623 * (es - ea)/air_pres*SEC_PER_DAY/DENSITY_WATER*MM_PER_METER; //[mm/d]
    //[kg/m3]*[m/s]*[s/d]/[kg/m3]*[mm/m]=[mm/d]
    return sublim;
  }
  //-----------------------------------------------------------------
  else if(type==SUBLIM_KUZMIN)
  {
    double Tave    =pHRU->GetForcingFunctions()->temp_daily_ave;
    double Tsnow   =pHRU->GetSnowTemperature();
    double ea      =GetSaturatedVaporPressure(Tave)*rel_hum*MB_PER_KPA;
    double es      =GetSaturatedVaporPressure(Tsnow)*1.0*MB_PER_KPA;

    return ((0.18 + 0.098 *wind_vel)*(es-ea)); //[mm/day]
  }
  //-----------------------------------------------------------------
  else if(type==SUBLIM_KUCHMENT_GELFAN)
  {
    //KutchmentGelfan1996
    // Effectively identical to KUZMIN but clearer to understand
    //supports negative sublimation (freezing of condensate) -uses hourly temperature
    double Tave    =pHRU->GetForcingFunctions()->temp_ave;
    double Tsnow   =pHRU->GetSnowTemperature();
    double ea      =GetSaturatedVaporPressure(Tave)*rel_hum*MB_PER_KPA;
    double es      =GetSaturatedVaporPressure(Tsnow)*1.0*MB_PER_KPA;

    double Q_E=32.82 *(0.18+0.098*wind_vel)*(es-ea)*WATT_TO_MJ_PER_D; //MJ/m2/d
    return Q_E/LH_SUBLIM/DENSITY_WATER*(MM_PER_METER); //[mm/d]
  }
  //-----------------------------------------------------------------
  else if(type==SUBLIM_BULK_AERO)
  { ///bulk aerodynamix flux formulation described in Hock & Holmgren [2005]
    //Hock, R., & Holmgren, B. (2005). A distributed surface energy-balance model for complex topography and
    //its application to Storglaciaren, Sweden. Journal of Glaciology, 51(172), 25-36. doi:10.3189/172756505781829566
    double air_pres=pHRU->GetForcingFunctions()->air_pres; //air pressure [kPa]
    double air_dens=pHRU->GetForcingFunctions()->air_dens;
    double Tsnow   =pHRU->GetSnowTemperature();

    double ea   =GetSaturatedVaporPressure(Ta   )*rel_hum;  //Vapour pressure of air[kPa]
    double es   =GetSaturatedVaporPressure(Tsnow)*1.0;      //snow surface vapour pressure(assume saturated)[kPa]

    double z_ref=2.0;            //reference height [m]
    double z0   =0.00005;        //aerodynamic roughness length for snow [m]
    double z0e  =0.00005;        //roughness length for temperature[m]

    //Calculate transfer coefficient for latent heat--> Assumes neutral conditions; probably ok for melt.
    //Flag to chat about: could add Richardson Numbers to this; Monin-Ohbukhov doesn't work great high relief due to gravity winds
    double C_E = (VON_KARMAN*VON_KARMAN) / (log(z_ref / z0) * log(z_ref / z0e));

    //Calculate Latent Heat Flux,then convert units to get sublimation rate
    double Q_E = air_dens * LH_SUBLIM * C_E * wind_vel * ((AIR_H20_MW_RAT * (es - ea)) / air_pres)*SEC_PER_DAY; //MJ/m2/d
    return  Q_E / (LH_SUBLIM * DENSITY_WATER)*MM_PER_METER; //[mm/d]
  }
  //-----------------------------------------------------------------
  else if(type==SUBLIM_CENTRAL_SIERRA)
  {
    double T         =pHRU->GetForcingFunctions()->temp_daily_ave;
    double sat_vap   = GetSaturatedVaporPressure(pHRU->GetSnowTemperature())*MB_PER_KPA;   // [millibar]
    double vap_pres  = GetSaturatedVaporPressure(T)*rel_hum*MB_PER_KPA;                    //[millibar]
    double wind_ht   = 2.0*FEET_PER_METER;                                                 // height of wind measurement [feet]
    double vap_ht    = 2.0*FEET_PER_METER;                                                 // height of vapor pressure measurement [feet]
    double wind_v    = wind_vel*MPH_PER_MPS;                                               // wind speed [miles per hour]

    return ((0.0063*(pow((wind_ht*vap_ht),(-1/6)))*wind_v*(sat_vap-vap_pres))*MM_PER_INCH);      // [mm/day]
  }
  //-----------------------------------------------------------------
  else if(type==SUBLIM_PBSM)
  {
    ExitGracefully("CmvSublimation:PBSM",STUB);
    return 0.0;
  }
  return 0.0;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns rates of loss of water to sublimation [mm/d]
///
/// \param *state_vars [in] Array of current State variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Rate of snow loss due to sublimation [mm/day]
//
void CmvSublimation::GetRatesOfChange(const double      *state_vars,
                                      const CHydroUnit  *pHRU,
                                      const optStruct   &Options,
                                      const time_struct &tt,
                                      double      *rates) const
{

  double wind_vel=pHRU->GetForcingFunctions()->wind_vel;          // [m/s]

  rates[0] = SublimationRate(state_vars, pHRU, pModel, tt, wind_vel, type); //Uses wind vel @ 2m

  if(_nConnections==2) {// simulating snow depth
    int iSnowDepth=pModel->GetStateVarIndex(SNOW_DEPTH);
    double SD=pHRU->GetStateVarValue(iSnowDepth);
    double SWE=pHRU->GetStateVarValue(iFrom[0]);
    if(SWE>0) {
      rates[1]=rates[0]*SD/SWE; //only works for single layer snow models.
    }
  }
};
//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow cannot drain "from" compartment over timestep
///
/// \param *state_vars [in] Array of current State variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Rate of snow loss due to sublimation [mm/day]

//
void  CmvSublimation::ApplyConstraints( const double             *state_vars,
                                        const CHydroUnit *pHRU,
                                        const optStruct  &Options,
                                        const time_struct &tt,
                                        double     *rates) const
{
  if (state_vars[iFrom[0]]<=0){rates[0]=0.0;}//reality check

  //cant remove more than is there
  rates[0]=threshMin(rates[0],state_vars[iFrom[0]]/Options.timestep,0.0);
}
