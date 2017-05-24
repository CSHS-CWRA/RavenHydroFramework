/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------
  Snow temp evolution
  ----------------------------------------------------------------*/
#include "HydroProcessABC.h"
#include "SnowMovers.h"

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the Snow temp evolution constructor
/// \param ste_type [in] Model of snow temp selected
//
CmvSnowTempEvolve::CmvSnowTempEvolve(snowtemp_evolve_type  ste_type):
  CHydroProcessABC(SNOWTEMP_EVOLVE)
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
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Rate of loss from "from" compartment [mm/d]
//
void CmvSnowTempEvolve::GetRatesOfChange( const double            *state_vars,
                                          const CHydroUnit  *pHRU,
                                          const optStruct         &Options,
                                          const time_struct &tt,
                                          double      *rates) const
{
  double Tsnow=state_vars[iFrom[0]];
  double Tair=pHRU->GetForcingFunctions()->temp_ave;
  if (_type==SNOTEMP_NEWTONS)
  {
    //linear heat transfer coefficient
    double alpha=CGlobalParams::GetParams()->airsnow_coeff; // = (1-x6) as used in original cema niege =[1/d]
    rates[0]=alpha*(Tair-Tsnow);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Applies constraints to calculations in GetRatesOfChange
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Rate of loss from "from" compartment [mm/d]
//
void CmvSnowTempEvolve::ApplyConstraints( const double      *state_vars,
                                          const CHydroUnit  *pHRU,
                                          const optStruct         &Options,
                                          const time_struct &tt,
                                          double      *rates) const
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
/// \param ste_type [in] modelling type selected
/// \param *aSV [out] Array of state variable types needed by algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by algorithm (size of aSV[] and aLev[] arrays)
//
void CmvSnowTempEvolve::GetParticipatingStateVarList(snowtemp_evolve_type ste_type,sv_type *aSV, int *aLev, int &nSV)
{
  nSV=1;
  aSV[0]=SNOW_TEMP; aLev[0]=DOESNT_EXIST;
}


