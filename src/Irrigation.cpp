/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2026 the Raven Development Team
----------------------------------------------------------------*/
#include "Irrigation.h"

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the irrigation constructor
/// \param itype    [in] Irrigation method
/// \param to_index [in] Index of state variable receiving irrigation water
/// \param *pModel  [in] Model reference
//
CmvIrrigation::CmvIrrigation(irrigation_type itype,
                             int             to_index,
                             CModelABC      *pModel)
  :CHydroProcessABC(IRRIGATION_INPUT, pModel)
{
  ExitGracefullyIf(to_index==DOESNT_EXIST,
                   "CmvIrrigation Constructor: invalid 'to' compartment specified",BAD_DATA);

  CHydroProcessABC::DynamicSpecifyConnections(2);
  iFrom[0]=pModel->GetStateVarIndex(ATMOS_PRECIP);   iTo[0]=pModel->GetStateVarIndex(IRRIGATION_SRC);
  iFrom[1]=pModel->GetStateVarIndex(IRRIGATION_SRC); iTo[1]=to_index;
  _type=itype;
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvIrrigation::~CmvIrrigation(){}

//////////////////////////////////////////////////////////////////
/// \brief Initializes irrigation object
//
void CmvIrrigation::Initialize(){}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
///
/// \param itype [in] Irrigation method
/// \param *aSV [out] Reference to array of state variables needed by irrigation algorithm
/// \param *aLev [out] Array of levels of multilevel state variables (or DOESNT_EXIST if single level)
/// \param &nSV [out] Number of participating state variables (length of aSV and aLev arrays)
//
void CmvIrrigation::GetParticipatingStateVarList(irrigation_type itype, sv_type *aSV, int *aLev, int &nSV)
{
  nSV=1;
  aSV [0]=IRRIGATION_SRC;  aLev[0]=DOESNT_EXIST;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rates of irrigation water supply [mm/d]
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Rates of water movement [mm/d]
//
void CmvIrrigation::GetRatesOfChange(const double      *state_vars,
                                     const CHydroUnit  *pHRU,
                                     const optStruct   &Options,
                                     const time_struct &tt,
                                     double            *rates) const
{
  if(_type==IRRIG_FROMFILE){
    rates[0]=pHRU->GetForcingFunctions()->irrigation; //[mm/d]
  }
  else{
    ExitGracefully("CmvIrrigation::GetRatesOfChange: undefined irrigation type",BAD_DATA);
  }
  rates[1]=rates[0];
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Rates of water movement [mm/d]
//
void CmvIrrigation::ApplyConstraints(const double      *state_vars,
                                     const CHydroUnit  *pHRU,
                                     const optStruct   &Options,
                                     const time_struct &tt,
                                     double            *rates) const
{
  rates[0]=max(rates[0],0.0);
  rates[1]=rates[0];
}
