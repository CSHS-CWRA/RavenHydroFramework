/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------
  Snow temp evolution
  ----------------------------------------------------------------*/
#include "HydroProcessABC.h"
#include "SnowMovers.h"
#include "Model.h"

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the Snow temp evolution constructor
/// \param ste_type [in] Model of snow temp selected
//
CmvSnowTempEvolve::CmvSnowTempEvolve(snowtemp_evolve_type ste_type, CModelABC *pModel)
  :CHydroProcessABC(SNOWTEMP_EVOLVE, pModel)
{
  _type=ste_type;
  CHydroProcessABC::DynamicSpecifyConnections(1);
  iFrom[0]=pModel->GetStateVarIndex  (SNOW_TEMP);
  iTo  [0]=iFrom[0];
}
//////////////////////////////////////////////////////////////////
/// \brief Implementation fo the default destructor
//
CmvSnowTempEvolve::~CmvSnowTempEvolve(){}

//////////////////////////////////////////////////////////////////
/// \brief Initializes snow temp modelling object
//
void CmvSnowTempEvolve::Initialize(){}

//////////////////////////////////////////////////////////////////
/// \brief Returns rates of change of state variables
///
/// \param *state_vars [in]  Array of current state variables in HRU
/// \param *pHRU       [in]  Reference to pertinent HRU
/// \param &Options    [in]  Global model options information
/// \param &tt         [in]  Current model time
/// \param *rates      [out] Rate of loss from "from" compartment [mm/d]
//
void CmvSnowTempEvolve::GetRatesOfChange( const double      *state_vars,
                                          const CHydroUnit  *pHRU,
                                          const optStruct   &Options,
                                          const time_struct &tt,
                                          double            *rates) const
{
  double Tsnow=state_vars[iFrom[0]];
  double Tair=pHRU->GetForcingFunctions()->temp_ave;
  if (_type==SNOTEMP_NEWTONS)
  {
    //linear heat transfer coefficient
    double alpha = pModel->GetGlobalParams()->GetParams()->airsnow_coeff; // = (1-x6) as used in original cema neige =[1/d]
    rates[0] = alpha*(Tair-Tsnow);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Applies constraints to calculations in GetRatesOfChange
///
/// \param *state_vars [in]  Array of current state variables in HRU
/// \param *pHRU       [in]  Reference to pertinent HRU
/// \param &Options    [in]  Global model options information
/// \param &tt         [in]  Current model time
/// \param *rates      [out] Rate of loss from "from" compartment [mm/d]
//
void CmvSnowTempEvolve::ApplyConstraints( const double      *state_vars,
                                          const CHydroUnit  *pHRU,
                                          const optStruct   &Options,
                                          const time_struct &tt,
                                          double            *rates) const
{
  //can't exceed a temperature of freezing
  double Tsnow=state_vars[iFrom[0]];
  rates[0]=min(rates[0],(FREEZING_TEMP-Tsnow)/Options.timestep);
}

/////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for snow melt algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by algorithm (size of aP[] and aPC[])
//
void CmvSnowTempEvolve::GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const
{
  nP=0;
  if (_type==SNOTEMP_NEWTONS)
  {
    aP[nP]="AIRSNOW_COEFF"; aPC[nP]=CLASS_GLOBAL; nP++;
  }
  else
  {
    ExitGracefully("CmvSnowTempEvolve::GetParticipatingParamList: undefined snow melt algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
///
/// \param ste_type [in]  modelling type selected
/// \param *aSV     [out] Array of state variable types needed by algorithm
/// \param *aLev    [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV     [out] Number of state variables required by algorithm (size of aSV[] and aLev[] arrays)
//
void CmvSnowTempEvolve::GetParticipatingStateVarList(snowtemp_evolve_type ste_type,sv_type *aSV, int *aLev, int &nSV)
{
  nSV=1;
  aSV[0]=SNOW_TEMP; aLev[0]=DOESNT_EXIST;
}

/*double Qs(double air_pres,double T){
  return 1.0; //TMP DEBUG
}

void GetSnowTemperature(force_struct *F, CHydroUnit *pHRU, const double ref_elev){
  const double EMISS=0.97;//?
  const double  LS=??;
  const double SBC=??;
  const double CP=??; //heat cap?

  double air_temp, air_dens;
  double LWi;
  double u1,zref,z0snow,ra;

  air_temp = F->temp_ave + ZERO_CELSIUS;
  air_dens = F->air_pres*1000/(UNIV_GAS_CONST*air_temp);
  u1       = F->wind_vel; // Wind speed (m/s)
  LWi      = F->LW_incoming;

  zref   = pHRU->GetVegVarProps()->reference_height;
  z0snow =1.0; //JRC - MUST ATTEND TO
  ra     = (log(zref/z0snow)*log(ref_elev/z0snow))/sqrt(VON_KARMAN)/u1;

  double Qss=Qs(F->air_pres,air_temp);
  double delta = AIR_H20_MW_RAT*LS*Qss/(UNIV_GAS_CONST*sqrt(air_temp));

  double q = F->rel_humidity*Qs(F->air_pres,air_temp); // specific humidity (kg/kg)

  double snow_temp = air_temp + (EMISS*(LWi - SBC*pow(air_temp,4.0)) + LS*(q - Qss)*air_dens/ra)/(4*EMISS*SBC*pow(air_temp,3.0) + (CP + LS*delta)*air_dens/ra);

  return min(snow_temp-ZERO_CELSIUS,FREEZING_TEMP);
}*/
