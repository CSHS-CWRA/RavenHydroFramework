/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ------------------------------------------------------------------
  Capillary Rise
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "SoilWaterMovers.h"

/*****************************************************************
   Recharge Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the standard constructor
/// \param cr_type [in] Model of capillary rise selected
/// \param In_index [in] Soil storage unit index from which water is lost
/// \param Out_index [in] Soil storage unit index to which water rises
//
CmvRecharge::CmvRecharge(int to_index)
  :CHydroProcessABC(RECHARGE)
{
  ExitGracefullyIf(to_index==DOESNT_EXIST,
                   "CmvRecharge Constructor: invalid 'to' compartment specified",BAD_DATA);

  CHydroProcessABC::DynamicSpecifyConnections(1);
  iFrom[0]=pModel->GetStateVarIndex(ATMOS_PRECIP); 
  iTo[0]=to_index;
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of default destructor
//
CmvRecharge::~CmvRecharge(){}

//////////////////////////////////////////////////////////////////
/// \brief Verifies iFrom - iTo connectivity
/// \details Ensures that water rises from soil layer/ groundwater
/// to another soil layer or groundwater unit
//
void CmvRecharge::Initialize()
{
  ExitGracefullyIf((pModel->GetStateVarType(iTo[0])!=SOIL) &&
                   (pModel->GetStateVarType(iTo[0])!=GROUNDWATER),
                   "CmvRecharge::Initialize:Recharge must go to soil or groundwater unit",BAD_DATA);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for capillary rise algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by capillary rise  algorithm (size of aP[] and aPC[])
//
void CmvRecharge::GetParticipatingParamList(string  *aP , class_type *aPC , int &nP) const
{
  nP=0;
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
/// \details User specifies from and to compartments, levels not known before construction
///
/// \param *aSV [out] Reference to state variable types needed by cr algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by cr algorithm (size of aSV[] and aLev[] arrays)
//
void CmvRecharge::GetParticipatingStateVarList(sv_type *aSV, int *aLev, int &nSV)
{
  nSV=0;
  //soil and atmos_precip will already be in model by default
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rate of loss of water from soil or aquifer to another soil layer [mm/d]
/// \details if type==CRISE_HBV, flow is linearly and inversely proportional to saturation of recieving unit \cite Bergstroem1992
///
/// \param *storage [in] Reference to soil storage from which water rises
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of loss of water from soil to another soil layer [mm/day]
//
void   CmvRecharge::GetRatesOfChange( const double      *storage,
                                      const CHydroUnit  *pHRU,
                                      const optStruct   &Options,
                                      const time_struct &tt,
                                      double            *rates) const
{
  if (pHRU->GetHRUType()!=HRU_STANDARD){return;}//Lakes & glaciers

  rates[0]=pHRU->GetForcingFunctions()->recharge;
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow cannot drain "from" compartment over timestep.
/// Presumes overfilling of "to" compartment is handled using cascade
///
/// \param *storage [in] state variable array for this HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of loss from soil to other soil layer [mm/day]
//
void   CmvRecharge::ApplyConstraints( const double     *storage,
                                           const CHydroUnit *pHRU,
                                           const optStruct  &Options,
                                           const time_struct &tt,
                                           double     *rates) const
{
  if (pHRU->GetHRUType()!=HRU_STANDARD){return;}//Lakes & glaciers

  //exceedance of max "to" compartment
  //water flow simply slows (or stops) so that receptor will not overfill during tstep
  rates[0]=threshMin(rates[0],
                     (pHRU->GetStateVarMax(iTo[0],storage,Options)-storage[iTo[0]])/Options.timestep,0.0);
}
