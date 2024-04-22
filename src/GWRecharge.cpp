/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2021 the Raven Development Team
------------------------------------------------------------------
*/
#include "GWSWProcesses.h"

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the gw recharge constructor
/// \param pGWModel [in] pointer to groundwater model
/// \param type [in] recharge type (enum)
/// \param From_index [in] Index of state variable/ storage unit where recharge is coming from
//
CGWRecharge::CGWRecharge(CGroundwaterModel *pGWModel,
                         gwrecharge_type    type,
                         int                From_index,
                         CModel            *pModel)
  :CGWSWProcessABC(pGWModel, GWRECHARGE, pModel)
{
  int iSW;
  _aFluxes = NULL;
  iSW = pModel->GetStateVarIndex(SURFACE_WATER);

  _recharge_type = type;
  _totalRecharge = 0.0;

  // ...better to just use different constructor?
  CHydroProcessABC::DynamicSpecifyConnections(1);

  iFrom[0] = From_index;
    iTo[0] = pModel->GetStateVarIndex(GROUNDWATER);
  setProcName("RAVEN RECHARGE");
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
//
CGWRecharge::~CGWRecharge()
{
  delete[] _aFluxes;
}

//////////////////////////////////////////////////////////////////
/// \brief Initializes arrays
///
/// \param nNodes [in] number of nodes the recharge process is moving water to
//
void CGWRecharge::InitializeRechargeClass(int nNodes)
{
  _nnodes        = nNodes;
  _inodes        = new int[nNodes];
  _totalRecharge = 0.0;
  if (_recharge_type == RECHARGE_FLUX) {
    _aFluxes     = new double[nNodes];
  }
  else if (_recharge_type == RECHARGE_HRU_DATA) {
    // Assumes all first layer is recharge (consistent with CGroundwaterModel::InitMFUSG initialization)
    for (int i = 0; i < _nnodes; i++) {
      _inodes[i] = i + 1; //indexed from 1, matching MFUSG
    }
  }
  else {
    _aFluxes     = NULL;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Adds Node when type is RECHARGE_FLUX
///
/// \param index  [in] location in arrays of where to store other parameters (why isn't the process tracking this itself?! -LS)
/// \param nodeID [in] node number in GW model where recharge will be moved
/// \param flux   [in] flux to be moved [mm/day] to node
//
void CGWRecharge::addNodeFlux(int index, int nodeID, double flux)
{
  _inodes [index] = nodeID;
  _aFluxes[index] = flux;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for drain algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by drain algorithm (size of aP[] and aPC[])
//
void CGWRecharge::GetParticipatingParamList(string  *aP, class_type *aPC, int &nP) const
{
  nP=0;
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
///
/// \param *aSV [out] Array of state variable types needed by recharge algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by recharge algorithm (size of aSV[] and aLev[] arrays)
//
void CGWRecharge::GetParticipatingStateVarList(sv_type *aSV, int *aLev, int &nSV)
{
  nSV=1;
  aSV [0]=GROUNDWATER;  aLev[0]=DOESNT_EXIST;
  // From SV is set by user
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rate of loss to subsurface
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of loss from "from" compartment [mm/day]
//
void CGWRecharge::GetRatesOfChange(const double		   *state_vars,
                                   const CHydroUnit  *pHRU,
                                   const optStruct	 &Options,
                                   const time_struct &tt,
                                   double      *rates) const
{
#ifdef _MODFLOW_USG_
  if (pHRU->GetHRUType() == HRU_STANDARD)
  {
    set<int> HRU_nodes;
    double   area, weight;
    int      index;
    double   node_rech;
    double   AMAT_change=0.0, RHS_change=0.0;

    //-- Initialize
    long long int HRUid = pHRU->GetHRUID();
    rates[0]  = 0.0;

    // Loop over nodes with recharge in HRU
    HRU_nodes = getHRUnodes((int)(HRUid)); //todo[funct] - support long long int

    //-- Recharge as a forcing
    if (_recharge_type == RECHARGE_HRU_DATA) {
      rates[0] = pHRU->GetForcingFunctions()->recharge; // [mm/day]

      //-- Loop over recharge nodes adding to gw flow equation variables
      for (auto n : HRU_nodes) {
        // Correct to topmost node with water
        int active_node = n;  // To avoid changing what we're looping over
        pGWModel->GetTopActiveNode(active_node);
        weight = pGWModel->GetOverlapWeight(HRUid, active_node);
        area   = pGWModel->GetNodeArea(active_node);
        AMAT_change = 0.0;
        RHS_change  = (rates[0] * area * weight) / MM_PER_METER;  // [m3/d]
        pGWModel->AddToGWEquation(active_node, AMAT_change, RHS_change);
      }
    }
    //-- Recharge controlled by a conductance value from source
    else if (_recharge_type == RECHARGE_FLUX) {

      for (auto n : HRU_nodes) {
        // Correct to topmost node with water
        int active_node = n;  // To avoid changing what we're looping over
        pGWModel->GetTopActiveNode(active_node);
        index       = _mNodeIndex.at(active_node);  // location of node in array
        weight      = pGWModel->GetOverlapWeight(HRUid, active_node);
        area        = pGWModel->GetNodeArea(active_node);
        node_rech   = _aFluxes[index] * weight; // [m/d]
        // SW SV rate
        rates[0] += node_rech * MM_PER_METER; // [mm/d]
        // GW Rate
        AMAT_change = 0.0;
        RHS_change = node_rech * area; // [m^3/d]
        pGWModel->AddToGWEquation(active_node, AMAT_change, RHS_change);
      }
    }

    else { ExitGracefully("CmvRecharge::GetRatesOfChange: undefined recharge type", BAD_DATA); }

  }
  else {//Lake/Glacier case (?) //GWUPDATE
    return;
  }
  // Important for not double-counting fluxes added to GWSV
  SendRateToFluxTracker(pHRU->GetGlobalIndex(), rates);
  #endif
}

//////////////////////////////////////////////////////////////////
/// \brief Fills arrays[_nnodes] with process rates to(+)/from(-) the GW model.
///        Placeholder method for child processes to override.
///
/// \param pModel [in] pointer to groundwater model
/// \param nodes [out] array of nodes where process is active
/// \param rates [out] rates of flux to GW model
//
void CGWRecharge::calcGWBudget(const CModel *pModel, int *nodes, double *rates)
{
#ifdef _MODFLOW_USG_
  int n, topnode;
  double weight, area;
  vector<int> nodeHRUs;
  CHydroUnit  *pHRU;

  //-- Loop over all recharge nodes
  for (int i = 0; i < _nnodes; i++)
  {
    n = _inodes[i];
    topnode = pGWModel->GetTopActiveNode(n); // recharge goes to top active node
    nodes[i] = topnode;
    rates[i] = 0.0;
    area = pGWModel->GetNodeArea(n);  // m

    //-- Recharge via DATA is by HRU
    if (_recharge_type == RECHARGE_HRU_DATA) {
      nodeHRUs = pGWModel->GetHRUsByNode(n);

      for (auto HRUid : nodeHRUs) {
        pHRU = pModel->GetHRUByID(HRUid);
        weight = pGWModel->GetOverlapWeight(HRUid, n);
        rates[i] += (pHRU->GetForcingFunctions()->recharge/MM_PER_METER * area * weight); // mm/d to m3/d
      }
    }

    //-- Recharge flux is by cell
    else if (_recharge_type == RECHARGE_FLUX) {
      rates[i] += _aFluxes[i] * area; // m/d to m3/d
    }
  }
#endif
}

//////////////////////////////////////////////////////////////////
/// \brief
/// \details
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of loss from "from" compartment [mm/day]
//
void   CGWRecharge::ApplyConstraints(const double		 *state_vars,
  const CHydroUnit *pHRU,
  const optStruct	 &Options,
  const time_struct &tt,
  double     *rates) const
{
  if (pHRU->GetHRUType() != HRU_STANDARD){ return; }//Lake/Glacier case

  //cant remove more than is there  // GWUPDATE needs to be from MODEL
  //rates[0]=threshMin(rates[0],state_vars[iFrom[0]]/Options.timestep,0.0);
  //cout<<"dr: "<<rates[0]<<endl;
}

gwrecharge_type CGWRecharge::getRechargeType(){ return(_recharge_type); }
