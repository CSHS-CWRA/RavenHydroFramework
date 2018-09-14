/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ------------------------------------------------------------------
  Interflow
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "SoilWaterMovers.h"

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the interflow constructor
///
/// \param itype [in] Model of interflow selected
/// \param In_index [in] Index of storage unit from which interflow is released
//
CmvInterflow::CmvInterflow(interflow_type       itype,
                           int                                            In_index)
  :CHydroProcessABC(INTERFLOW)
{
  type =itype;
  ExitGracefullyIf(In_index==DOESNT_EXIST,
                   "CmvInterflow Constructor: invalid 'from' compartment specified",BAD_DATA);
  CHydroProcessABC::DynamicSpecifyConnections(1);
  iFrom[0]=In_index; iTo[0]=pModel->GetStateVarIndex(SURFACE_WATER);
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvInterflow::~CmvInterflow(){}

//////////////////////////////////////////////////////////////////
/// \brief Validate iFrom connection of the interflow; ensures 'from' compartment is soil
//
void CmvInterflow::Initialize()
{
  if      (pModel->GetStateVarType(iFrom[0])!=SOIL) {
    ExitGracefully("CmvInterflow::Invalid 'from' type",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for interflow algorithm
/// \param *aPC [out] Array of class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by interflow algorithm (size of aP[] and aPC[])
//
void CmvInterflow::GetParticipatingParamList(string  *aP , class_type *aPC , int &nP) const
{
  nP=0;
  if (type==INTERFLOW_PRMS)
  {
    nP=4;
    aP[0]="MAX_INTERFLOW_RATE"; aPC[0]=CLASS_SOIL;
    aP[1]="POROSITY";           aPC[1]=CLASS_SOIL;
    aP[2]="FIELD_CAPACITY";     aPC[2]=CLASS_SOIL;
    aP[3]="SAT_WILT";           aPC[3]=CLASS_SOIL;
  }
  else
  {
    ExitGracefully("CmvInterflow::GetParticipatingParamList: undefined interflow algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
/// \note Out connection is always to surface water, From connection must be a user specified soil unit
///
/// \param itype [in] Model of interflow used
/// \param *aSV [out] Array of state variable types needed by interflow algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by interflow evaporation algorithm (size of aSV[] and aLev[] arrays)
//
void CmvInterflow::GetParticipatingStateVarList(interflow_type  itype,
                                                sv_type *aSV, int *aLev, int &nSV)
{
  //user specified compartments, out is always to SW
  nSV=2;
  aSV [0]=SOIL;  aSV [1]=SURFACE_WATER;
  aLev[0]=0;     aLev[1]=DOESNT_EXIST;
  //aLev[0] is actually user specified, but 1 layer is needed if 5 are
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rate of loss of water from soil to surface water [mm/d]
/// \details  if type==INTERFLOW_PRMS, uses PRMS approach
///
/// \param *state_vars [in] Array of current state variables
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Rate of loss of water from soil compartment to surface water [mm/d]
//
void   CmvInterflow::GetRatesOfChange(const double                      *state_vars,
                                      const CHydroUnit      *pHRU,
                                      const optStruct               &Options,
                                      const time_struct &tt,
                                      double                      *rates) const
{
  if (pHRU->GetHRUType()!=HRU_STANDARD){return;}//Lakes & glaciers

  double stor = state_vars[iFrom[0]];
  int    m    = pModel->GetStateVarLayer(iFrom[0]);

  if (type==INTERFLOW_PRMS)
  {
    double tens_stor,K,max_stor;

    max_stor  = pHRU  ->GetSoilCapacity(m);                     //maximum storage of soil layer [mm]
    tens_stor = pHRU  ->GetSoilTensionStorageCapacity(m);       // [mm]
    K         = pHRU  ->GetSoilProps(m)->max_interflow_rate;    //interflow rate [mm/d]

    rates[0] = K * max(stor - tens_stor ,0.0)/(max_stor-tens_stor);
  }
  else{
    rates[0]=0.0;
    ExitGracefully("CmvInterflow::GetRatesOfChange: undefined interflow type",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow cannot drain "from" compartment over timestep.
///
/// \param *state_vars [in] Array of current state variables
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Rate of loss of water from soil compartment to surface water [mm/d]
//
void   CmvInterflow::ApplyConstraints( const double              *state_vars,
                                       const CHydroUnit *pHRU,
                                       const optStruct        &Options,
                                       const time_struct &tt,
                                       double     *rates) const
{
  if (pHRU->GetHRUType()!=HRU_STANDARD){return;}//Lakes & glaciers

  //cant remove more than is there
  //shouldn't be able to remove below wilt_cap, either...
  rates[0]=threshMin(rates[0],state_vars[iFrom[0]]/Options.timestep,0.0);
}
