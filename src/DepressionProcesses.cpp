/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ------------------------------------------------------------------
  DepressionOverflow
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "DepressionProcesses.h"

/*****************************************************************
   DepressionOverflow Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of abstraction constructor
/// \param absttype [in] Selected model of abstraction
//
CmvDepressionOverflow::CmvDepressionOverflow(depflow_type dtype)
  :CHydroProcessABC(DEPRESSION_OVERFLOW)
{
  type=dtype;

  CHydroProcessABC::DynamicSpecifyConnections(1);
  //abstraction (ponded-->depression)
  iFrom[0]=pModel->GetStateVarIndex(DEPRESSION);
  iTo  [0]=pModel->GetStateVarIndex(SURFACE_WATER);

}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvDepressionOverflow::~CmvDepressionOverflow(){}

//////////////////////////////////////////////////////////////////
/// \brief Initializes abstraction object
//
void   CmvDepressionOverflow::Initialize()
{
}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for abstraction algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by abstraction algorithm (size of aP[] and aPC[])
//
void CmvDepressionOverflow::GetParticipatingParamList(string  *aP , class_type *aPC , int &nP) const
{
  if (type==DFLOW_THRESHPOW)
  {
    nP=4;
    aP[0]="DEP_THRESHHOLD";                     aPC[0]=CLASS_LANDUSE;
    aP[1]="DEP_N";                              aPC[1]=CLASS_LANDUSE;
    aP[2]="DEP_MAX_FLOW";       aPC[2]=CLASS_LANDUSE;
    aP[3]="DEP_MAX";            aPC[3]=CLASS_LANDUSE;

  }
  else
  {
    ExitGracefully("CmvDepressionOverflow::GetParticipatingParamList: undefined depression overflow algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variable list
///
/// \param absttype [in] Selected abstraction model
/// \param *aSV [out] Reference to array of state variables needed by abstraction algorithm
/// \param *aLev [out] Array of levels of multilevel state variables (or DOESNT_EXIST of single level)
/// \param &nSV [out] Number of participating state variables (length of aSV and aLev arrays)
//
void CmvDepressionOverflow::GetParticipatingStateVarList(depflow_type dtype, sv_type *aSV, int *aLev, int &nSV)
{
  nSV=2;
  aSV[0]=DEPRESSION;     aLev[0]=DOESNT_EXIST;
  aSV[1]=SURFACE_WATER;  aLev[1]=DOESNT_EXIST;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rates of water loss to abstraction
//
/// \param *state_vars [in] Array of state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time at which abstraction is to be calculated
/// \param *rates [out] rates[0]= rate of abstraction [mm/d]
//
void   CmvDepressionOverflow::GetRatesOfChange( const double                    *state_vars,
                                                const CHydroUnit       *pHRU,
                                                const optStruct        &Options,
                                                const time_struct &tt,
                                                double     *rates) const
{

  double stor=state_vars[iFrom[0]];

  //----------------------------------------------------------------------------
  if (type==DFLOW_THRESHPOW)
  {
    double max_flow   =pHRU->GetSurfaceProps()->dep_max_flow;
    double n          =pHRU->GetSurfaceProps()->dep_n;
    double thresh_stor=pHRU->GetSurfaceProps()->dep_threshhold;
    double max_stor   =pHRU->GetSurfaceProps()->dep_max;

    if (max_stor<thresh_stor){
      rates[0] = max_flow;
    }
    else {
      rates[0] = max_flow*pow(max((stor - thresh_stor) / (max_stor - thresh_stor), 0.0), n);
    }
  }

}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow cannot drain "from" compartment over timestep.
///
/// \param *state_vars [in] Array of state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time at which abstraction is to be calculated
/// \param *rates [out] rates[0]= rate of abstraction [mm/d]
//
void   CmvDepressionOverflow::ApplyConstraints(const double              *state_vars,
                                               const CHydroUnit *pHRU,
                                               const optStruct      &Options,
                                               const time_struct &tt,
                                               double     *rates) const
{
  //cant remove more than is there
  rates[0]=threshMin(rates[0],max(state_vars[iFrom[0]]/Options.timestep,0.0),0.0);

}
