/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team

  Open Water Evap
  Lake Evap
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "OpenWaterEvap.h"

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the Open water evaporation constructor
/// \param owtype [in] Model of open water evaporation used
//
CmvOWEvaporation::CmvOWEvaporation(owevap_type owtype, const int i_from)
  :CHydroProcessABC(OPEN_WATER_EVAPORATION)
{
  type =owtype;

  CHydroProcessABC::DynamicSpecifyConnections(1);//nConnections=1
  iFrom[0]=i_from;
  iTo  [0]=pModel->GetStateVarIndex(ATMOSPHERE);    //rates[0]: PONDED_WATER/DEPRESSION->ATMOSPHERE
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvOWEvaporation::~CmvOWEvaporation(){}

//////////////////////////////////////////////////////////////////
/// \brief Verify that iFrom[] - iTo[] connection is DEPRESSION-ATMOSPHERE
//
void CmvOWEvaporation::Initialize()
{
  sv_type typ=pModel->GetStateVarType(iFrom[0]);
  ExitGracefullyIf((typ!=DEPRESSION) && (typ!=PONDED_WATER),
                   "CmvOWEvaporation::Initialize:Open Water evaporation must come from depression unit",BAD_DATA);
  ExitGracefullyIf(pModel->GetStateVarType(iTo[0])!=ATMOSPHERE,
                   "CmvOWEvaporation::Initialize:Open Water evaporation must go to atmosphere",BAD_DATA);
}
//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for baseflow algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by baseflow algorithm (size of aP[] and aPC[])
//
void CmvOWEvaporation::GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const
{
  nP=1;
  aP[0]="OW_PET_CORR";   aPC[0]=CLASS_LANDUSE;
}
//////////////////////////////////////////////////////////////////
/// \brief Sets reference to state variable types needed by evaporation algorithm
///
/// \param owtype [in] Model of open water evaporation used
/// \param *aSV [out] Array of state variable types needed by ow evaporation algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by ow evaporation algorithm (size of aSV[] and aLev[] arrays)
//
void CmvOWEvaporation::GetParticipatingStateVarList(owevap_type owtype, sv_type *aSV, int *aLev, int &nSV)
{
  nSV=1;
  aSV[0]=ATMOSPHERE;  aLev[0]=DOESNT_EXIST;

  //other state var is specified in constructor
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rate of loss of water from open water to atmosphere [mm/d]
/// \details if type==OPEN_WATER_EVAP, evaporation is calculated using fraction of PET method
///
/// \param *state_vars [in] Current array of state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time strucutre
/// \param *rates [out] rates[0] is rate of loss from open water to atmosphere [mm/day]
//
void CmvOWEvaporation::GetRatesOfChange( const double                   *state_vars,
                                         const CHydroUnit        *pHRU,
                                         const optStruct   &Options,
                                         const time_struct &tt,
                                         double      *rates) const
{
  double OWPET;
  OWPET = pHRU->GetForcingFunctions()->OW_PET;            //open water PET rate [mm/d]
  
  if(pHRU->IsLinkedToReservoir()){return;}//reservoir-linked HRUs handle ET via reservoir MB

  if (type==OPEN_WATER_EVAP)//-------------------------------------
  {
    rates[0]=pHRU->GetSurfaceProps()->ow_PET_corr*OWPET;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate ET cannot drain depression storage over timestep
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &t [in] Current model time structure
/// \param *rates [out] rates[0] is rate of loss from open water to atmosphere [mm/day]
//
void  CmvOWEvaporation::ApplyConstraints( const double            *state_vars,
                                          const CHydroUnit  *pHRU,
                                          const optStruct   &Options,
                                          const time_struct &t,
                                          double      *rates) const
{

  if (rates[0]<0)             {rates[0]=0.0;}//positivity constraint

  //can't remove more than is there
  rates[0]=threshMin(rates[0],state_vars[iFrom[0]]/Options.timestep,0.0);
  if (state_vars[iFrom[0]]<=0){rates[0]=0.0;}//reality check

  g_debug_vars[1]=rates[0];//Used Evap

}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the lake evaporation constructor
/// \param lktype [in] Model of open water evaporation used
/// \param fromIndex [in] Index of lake from which water evaporates
//
CmvLakeEvaporation::CmvLakeEvaporation(lakeevap_type lktype, const int fromIndex)
  :CHydroProcessABC(LAKE_EVAPORATION)
{
  type =lktype;
  ExitGracefullyIf(fromIndex==DOESNT_EXIST,
                   "CmvLakeEvaporation Constructor: invalid 'from' compartment specified",BAD_DATA);

  CHydroProcessABC::DynamicSpecifyConnections(1);//nConnections=1
  iFrom[0]=fromIndex;
  iTo  [0]=pModel->GetStateVarIndex(ATMOSPHERE);;     //rates[0]: LAKE->ATMOSPHERE
}
//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvLakeEvaporation::~CmvLakeEvaporation(){}

//////////////////////////////////////////////////////////////////
/// \brief Verify that process moves water to atmosphere
//
void CmvLakeEvaporation::Initialize()
{
  ExitGracefullyIf(pModel->GetStateVarType(iTo[0])!=ATMOSPHERE,
                   "CmvLakeEvaporation::Initialize:Open Water evaporation must go to atmosphere",BAD_DATA);
}
//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for baseflow algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by baseflow algorithm (size of aP[] and aPC[])
//
void CmvLakeEvaporation::GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const
{
  if(type==LAKE_EVAP_BASIC)//-------------------------------------
  {
    nP=1;
    aP[0]="LAKE_PET_CORR"; aPC[0]=CLASS_LANDUSE;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Sets reference to state variable types needed by evaporation algorithm
/// \note "From" compartment specified by :LakeStorage command
///
/// \param lktype [in] Model of open water evaporation used
/// \param *aSV [out] Array of state variable types needed by ow evaporation algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by lake evaporation algorithm (size of aSV[] and aLev[] arrays)
//
void CmvLakeEvaporation::GetParticipatingStateVarList(lakeevap_type lktype, sv_type *aSV, int *aLev, int &nSV)
{
  nSV=1;
  aSV[0]=ATMOSPHERE;  aLev[0]=DOESNT_EXIST;

  //'from' compartment specified by :LakeStorage command
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rate of loss of water from lake to atmosphere [mm/d]
/// \details  if type==LAKE_EVAP_BASIC, evaporation is calculated using fraction of PET method
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time structure
/// \param *rates [out] rates[0] is rate of loss from lake to atmosphere [mm/d]
//
void CmvLakeEvaporation::GetRatesOfChange(const double                  *state_vars,
                                          const CHydroUnit      *pHRU,
                                          const optStruct         &Options,
                                          const time_struct &tt,
                                          double      *rates) const
{
  if ((!pHRU->IsLake()) && (pModel->GetLakeStorageIndex()!=iFrom[0])){return;}

  double OWPET;
  OWPET = pHRU->GetForcingFunctions()->OW_PET;          //calls PET rate [mm/d]

  if(pHRU->IsLinkedToReservoir()){return;}//reservoir-linked HRUs handle ET via reservoir MB

  if (type==LAKE_EVAP_BASIC)//-------------------------------------
  {
    rates[0]= pHRU->GetSurfaceProps()->lake_PET_corr*OWPET;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow ccannot drain lake over timestep
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time structure
/// \param *rates [out] rates[0] is rate of loss from lake to atmosphere [mm/d]
//
void  CmvLakeEvaporation::ApplyConstraints( const double                *state_vars,
                                            const CHydroUnit  *pHRU,
                                            const optStruct   &Options,
                                            const time_struct &t,
                                            double      *rates) const
{
  if (state_vars[iFrom[0]]<=0){rates[0]=0.0;}//reality check
  if (rates[0]<0)             {rates[0]=0.0;}//positivity constraint

  //can't remove more than is there
  rates[0]=threshMin(rates[0],state_vars[iFrom[0]]/Options.timestep,0.0);
}
