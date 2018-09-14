/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2018 the Raven Development Team
  ----------------------------------------------------------------
  Sublimation
  ----------------------------------------------------------------*/
#include "HydroProcessABC.h"
#include "SnowMovers.h"

/*****************************************************************
   Sublimation Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the sublimation constructor
/// \param sub_type [in] Sublimation-modelling technique selected
//
CmvSublimation::CmvSublimation(sublimation_type sub_type):
  CHydroProcessABC(SUBLIMATION)
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
    CHydroProcessABC::DynamicSpecifyConnections(2);//nConnections=1
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
  if (type==SUBLIM_SVERDRUP)
  {
    nP=2;
    aP[0]="SNOW_ROUGHNESS";   aPC[0]=CLASS_GLOBAL;
    aP[1]="SNOW_TEMPERATURE"; aPC[1]=CLASS_GLOBAL;
  }
  else if (type==SUBLIM_KUZMIN)
  {
    nP=0;
  }
  else if (type==SUBLIM_CENTRAL_SIERRA)
  {
    nP=1;
    aP[0]="SNOW_TEMPERATURE"; aPC[0]=CLASS_GLOBAL;
  }
  else if (type==SUBLIM_PBSM)
  {
    nP=0;
    //algorithm not complete
  }
  else if (type==SUBLIM_WILLIAMS)
  {
    nP=1;
    aP[0]="SNOW_TEMPERATURE"; aPC[0]=CLASS_GLOBAL;
  }
  else if(type==SUBLIM_CRHM_MARKS)
  {
    nP=0;
  }
  else
  {
    ExitGracefully("CmvSublimation::GetParticipatingParamList: undefined sublimation algorithm",BAD_DATA);
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

//////////////////////////////////////////////////////////////////
/// \brief Returns rates of loss of water to sublimation [mm/d]
/// \todo [add funct] SUBLIM_SVERDRUP needs to be revised
///
/// \param *state_vars [in] Array of current State variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Rate of snow loss due to sublimation [mm/day]
//
void CmvSublimation::GetRatesOfChange(const double               *state_vars,
                                      const CHydroUnit *pHRU,
                                      const optStruct    &Options,
                                      const time_struct &tt,
                                      double     *rates) const
{
  double Ta    =pHRU->GetForcingFunctions()->temp_ave;

  if(type==SUBLIM_SVERDRUP)
  {
    double numer,denom,sat_vap,vap_pres,air_density,roughness,rel_humid,wind_vel,air_pres;

    roughness   = (CGlobalParams::GetParams()->snow_roughness)/MM_PER_CM;                          // [cm]
    air_density = (pHRU->GetForcingFunctions()->air_dens)*GRAMS_PER_KG/CM3_PER_METER3;             // [grams/cm3]
    rel_humid   = pHRU->GetForcingFunctions()->rel_humidity;                                       
    wind_vel    = (pHRU->GetForcingFunctions()->wind_vel)*CM_PER_METER;                            // [cm/s]
    air_pres    = (pHRU->GetForcingFunctions()->air_pres)*MB_PER_KPA;                              // [millibars]
    sat_vap     = (GetSaturatedVaporPressure(pHRU->GetSnowTemperature()))*MB_PER_KPA;              // [millibars]
    vap_pres    = (GetVaporPressure((pHRU->GetForcingFunctions()->temp_ave),rel_humid))*MB_PER_KPA;// [millibars]

    numer = 0.623 * air_density * (pow(VON_KARMAN,2)) * wind_vel * (sat_vap - vap_pres);
    denom = air_pres * (pow((log(800/roughness)),2));

    rates[0]    = ((numer/denom));  //[gm/cm2/sec]
    rates[0] = rates[0]*SEC_PER_DAY/DENSITY_WATER/GRAMS_PER_KG*MM_PER_METER/METER_PER_CM/METER_PER_CM; //[mm/day]

  }
  else if(type==SUBLIM_KUZMIN)                          // needs to be tested still
  {
    double sat_vap,vap_pres;
    sat_vap             = (GetSaturatedVaporPressure(Ta))*MB_PER_KPA;   // [millibar]
    vap_pres    = (GetVaporPressure(Ta,(pHRU->GetForcingFunctions()->rel_humidity)))*MB_PER_KPA;        //[millibar]

    rates[0]= ((0.18 + 0.098 * (pHRU->GetForcingFunctions()->wind_vel))*(sat_vap-vap_pres)); //[mm/day]
  }
  else if(type==SUBLIM_CENTRAL_SIERRA)                  // needs to be tested still
  {
    double sat_vap,vap_pres,wind_ht,vap_ht,wind_v;

    sat_vap   = (GetSaturatedVaporPressure(pHRU->GetSnowTemperature()))*MB_PER_KPA;   // [millibar]
    vap_pres  = (GetVaporPressure(Ta,(pHRU->GetForcingFunctions()->rel_humidity)))*MB_PER_KPA;        //[millibar]
    wind_ht   = 2.0*FEET_PER_METER;    // height of wind measurement [feet]
    vap_ht    = 2.0*FEET_PER_METER;    // height of vapor pressure measurement [feet]
    wind_v    = (pHRU->GetForcingFunctions()->wind_vel   )*MPH_PER_MPS;                       // wind speed [miles per hour]

    rates[0]    = ((0.0063*(pow((wind_ht*vap_ht),(-1/6)))*(sat_vap-vap_pres)*wind_v)*MM_PER_INCH);      // [mm/day]

  }
  else if(type==SUBLIM_PBSM)
  {
    ExitGracefully("CmvSublimation:PBSM",STUB);
  }
  else if(type==SUBLIM_WILLIAMS)
  {
    double sat_vap,vap_pres,air_density,rel_humid,wind_vel;

    air_density = (pHRU->GetForcingFunctions()->air_dens)*GRAMS_PER_KG/CM3_PER_METER3;  // [grams/cm3]
    rel_humid   = (pHRU->GetForcingFunctions()->rel_humidity);                                                                                  // []
    wind_vel    = (pHRU->GetForcingFunctions()->wind_vel)*CM_PER_METER;                                 // [cm/s]
    sat_vap     = (GetSaturatedVaporPressure(pHRU->GetSnowTemperature()))*MB_PER_KPA;           // [millibars]
    vap_pres    = (GetVaporPressure((pHRU->GetForcingFunctions()->temp_ave),rel_humid))*MB_PER_KPA;     // [millibars]

    ExitGracefully("CmvSublimation:WILLIAMS: not working correctly",STUB);//***** Resolve errors in output, way too high/low!

    rates[0]= ((0.00011 * air_density * wind_vel * (vap_pres - sat_vap)) * MIN_PER_DAY * (1/(DENSITY_ICE * GRAMS_PER_KG *(1/CM3_PER_METER3))) * MM_PER_CM) / Options.timestep;          // [mm/day]
    // 0.00011 sec/min/mb
  }
  else if(type==SUBLIM_CRHM_MARKS) 
  {
    double T       =pHRU->GetForcingFunctions()->temp_daily_ave;
    double rel_hum =pHRU->GetForcingFunctions()->rel_humidity;
    double wind_vel=pHRU->GetForcingFunctions()->wind_vel;
    double ea      =GetSaturatedVaporPressure(T)*rel_hum;

    rates[0] =0.08*(0.18+0.098*wind_vel)*(6.11-ea*10.0)/LH_SUBLIM;
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
  if (rates[0]<0)             {rates[0]=0.0;}//positivity constraint

  //cant remove more than is there
  rates[0]=threshMin(rates[0],state_vars[iFrom[0]]/Options.timestep,0.0);
}
