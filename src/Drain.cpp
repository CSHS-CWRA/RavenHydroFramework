/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
------------------------------------------------------------------
  GW Drain Condition
----------------------------------------------------------------*/

#include "GWSWProcesses.h"

/*****************************************************************
   Drain Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the drain constructor
/// \param pGWModel [in] pointer to groundwater model
//
CGWDrain::CGWDrain(CGroundwaterModel *pGWModel,
                   CModel            *pModel)
  :CGWSWProcessABC(pGWModel, DRAIN, pModel)
{
  _aElevations   = NULL;
  _aConductances = NULL;
  _aHeads    = NULL;
  int iSW;
  iSW=pModel->GetStateVarIndex(SURFACE_WATER);   // NEEDS TO BE AN ARGUMENT

  CHydroProcessABC::DynamicSpecifyConnections(1);

  iFrom[0]=pModel->GetStateVarIndex(GROUNDWATER);
  iTo[0]=iSW;
  setProcName("RAVEN DRAINS");
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
//
CGWDrain::~CGWDrain()
{
  delete[] _aElevations;
  delete[] _aConductances;
  delete[] _aHeads;
}

//////////////////////////////////////////////////////////////////
/// \brief Initializes arrays
//
void CGWDrain::InitializeDrainClass(int nDrains)
{
  _nnodes        = nDrains;
  _inodes        = new int   [nDrains];
  _aElevations   = new double[nDrains];
  _aConductances = new double[nDrains];
  _aHeads    = new double[nDrains];
}

//////////////////////////////////////////////////////////////////
/// \brief Adds Drain Node
//
void CGWDrain::addDrain(int drainID, int nodeID, double drain_elev, double drain_cond)
{
  _inodes       [drainID] = nodeID;
  _aElevations  [drainID] = drain_elev;
  _aConductances[drainID] = drain_cond;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for drain algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by drain algorithm (size of aP[] and aPC[])
//
void CGWDrain::GetParticipatingParamList(string  *aP , class_type *aPC , int &nP) const
{

}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
///
/// \param se_type [in] Model of drain used
/// \param *aSV [out] Array of state variable types needed by Drain algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by drain algorithm (size of aSV[] and aLev[] arrays)
//
void CGWDrain::GetParticipatingStateVarList(sv_type *aSV, int *aLev, int &nSV)
{
  nSV=1;
  aSV [0]=GROUNDWATER;  aLev[0]=DOESNT_EXIST;
  // To SV is set by user
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rates of loss from set of soil layers & aquifer to SW due to Drains [mm/day]
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of loss from "from" compartment [mm/day]
//
void CGWDrain::GetRatesOfChange (const double		   *state_vars,
                                 const CHydroUnit  *pHRU,
                                 const optStruct	 &Options,
                                 const time_struct &tt,
                                       double      *rates) const
{
  if (pHRU->GetHRUType() != HRU_STANDARD) { return; } // CHECK
#ifdef _MODFLOW_USG_
  set<int> HRU_nodes;
  int      index;
  double   head, area;
  double   weight;
  double   drain_flow, AMAT_change, RHS_change;

  //-- Drainage is any amount of head over drain elevation controlled by conductance layer
  long long int HRUid = pHRU->GetHRUID();
  rates[0] = 0.0;


  //-- Loop over drain nodes in HRU
  HRU_nodes = getHRUnodes((int)(HRUid));

  for (int n=0; n< HRU_nodes.size(); n++)
  {
    index = _mNodeIndex.at(n);  // location of node in arrays

    // Replicating MFUSG GWF2DRN7U1FM subroutine

    //-- See if node head is above elevation
    head = pGWModel->GetNodeHead(n);
    _aHeads[index] = head;
    if (head > _aElevations[index])
    {
      //-- Get HRU-node overlap weight
      weight = pGWModel->GetOverlapWeight(HRUid, n);

      //-- Calculate rate of water volume leaving through drain
      //-- (Convention is positive leaving aquifer)
      drain_flow = weight * _aConductances[index] * (head - _aElevations[index]);  //modflow drainage  [m3/d]

      //-- Add to Flow Equation
      AMAT_change = weight * _aConductances[index];
      RHS_change  = weight * _aConductances[index] * _aElevations[index];

      //-- Node Area
      area = pGWModel->GetNodeArea(n);

      //-- Rate to SW SV
      rates[0] += (drain_flow / area) * MM_PER_METER;

      //-- Rate to GW Modle
      pGWModel->AddToGWEquation(n, AMAT_change, RHS_change);
    }
  }
  // Important for not double-counting fluxes added to GWSV
  SendRateToFluxTracker(pHRU->GetGlobalIndex(), rates);
#endif
}

//////////////////////////////////////////////////////////////////
/// \brief Fills arrays[_nnodes] with drain rates to(+)/from(-) the GW model.
///
/// \param pModel [in] pointer to groundwater model
/// \param nodes [out] array of drain nodes
/// \param rates [out] rates of flux to GW model
//
void CGWDrain::calcGWBudget(const CModel *pModel, int *nodes, double *rates)
{
  // Replicates MFUSG Subroutine GWF2DRN7U1BD

  // Loop over all drains
  for (int i = 0; i < _nnodes; i++)
  {
    nodes[i] = _inodes[i];
    // HRU-cell overlap weights are not needed here, since nothing is being calculated at the HRU level

    // Update Heads
    _aHeads[i] = pGWModel->GetNodeHead(nodes[i]);

    //-- Check if head over drain elevation
    if (_aHeads[i] > _aElevations[i])
    {
      // outflow = C * (elev_water - elev_drain)
      rates[i] = _aConductances[i] * (_aElevations[i] - _aHeads[i]);
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow cannot drain "from" compartment over timestep
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of loss from "from" compartment [mm/day]
//
void   CGWDrain::ApplyConstraints( const double		 *state_vars,
                                      const CHydroUnit *pHRU,
                                      const optStruct	 &Options,
                                      const time_struct &tt,
                                            double     *rates) const
{
  if (pHRU->GetHRUType()!=HRU_STANDARD){return;}//Lake/Glacier case

  //cant remove more than is there  // GWUPDATE needs to be from MODEL
  //rates[0]=threshMin(rates[0],state_vars[iFrom[0]]/Options.timestep,0.0);
  //cout<<"dr: "<<rates[0]<<endl;
}

//////////////////////////////////////////////////////////////////
/// \brief Updates SV rates of change given new head. Placeholder method for child processes to override.
///
/// \param *state_vars [in] Array of state variables stored in current HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model option information
/// \param &tt [in] Current model time
/// \param *rates [out] Rates of water/energy moved between storage locations / rates of change of modified state variables (size: _nConnections+_nCascades)
//
void CGWDrain::UpdateRatesOfChange(const double      *state_vars,
                                   const CHydroUnit  *pHRU,
                                   const optStruct   &Options,
                                   const time_struct &tt,
                                         double      *rates) const
{
  //-- Would it be worth it to store the previous rate of change to avoid recalculation?
  set<int> HRU_nodes;
  int      index;
  double   old_head, area, weight, flow_change;
  long long int  HRUid = pHRU->GetHRUID();
  rates[0] = 0.0;

  //-- Loop over drain nodes in HRU
  HRU_nodes = getHRUnodes((int)(HRUid)); //todo[funct] - support long long int

  for (unsigned int n=0; n< HRU_nodes.size(); n++)
  {
    index = _mNodeIndex.at(n);  // location of node in arrays

    //-- Store old head
    old_head = _aHeads[index];

    //-- Update head, see if node head is above elevation
    _aHeads[index] = pGWModel->GetNodeHead(n);

    if (_aHeads[index] > _aElevations[index])
    {
      // But what if old head wasn't? Set to Drain Elevation
      if (old_head < _aElevations[index]) { old_head = _aElevations[index]; }

      //-- Get HRU-node overlap weight & node area
      weight = pGWModel->GetOverlapWeight((int)(HRUid), n);
      area = pGWModel->GetNodeArea(n);

      //-- Calc drain flow change
      flow_change = weight * _aConductances[index] * ((_aHeads[index] - _aElevations[index]) - (old_head - _aElevations[index]));

      //-- Flux Change
      rates[0] += (flow_change / area) * MM_PER_METER;
    }
  }
}
