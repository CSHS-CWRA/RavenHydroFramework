/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2023 the Raven Development Team, Leland Scantlebury
------------------------------------------------------------------*/
#include "RavenInclude.h"
#include "MFUSGpp.h"
#include "GroundwaterModel.h"
#include "GWSWProcesses.h"
#include "GWRiverConnection.h"

// Constructor
CGroundwaterModel::CGroundwaterModel(CModel *pModel)
{
  _pModel = pModel;
  _pRiver = new CGWRiverConnection(this);  // setup automatically

  _nProcesses = 0;
  _nNodes = 0;
  _nLayers = 0;
  _aNodesPerLayer = NULL;

  _nStressPeriods = 0;
  _aStressPeriodLength = NULL;
  _aNumSteps = NULL;

  _maxiter = 0;
  _niunit = 0;
  _Converged = false;
  _doSolve = 0;

  _iStressPeriod = 0;
  _iStep = -1;

  _aProcessesAMAT = NULL;
  _aProcessesRHS = NULL;
  _aGWSWFluxes = NULL;
}

// Destructor
CGroundwaterModel::~CGroundwaterModel()
{
  free(_aNodesPerLayer);
  free(_aNumSteps);
  free(_aStressPeriodLength);
  delete [] _aProcessesAMAT;
  delete [] _aProcessesRHS;
  delete [] _aGWSWFluxes;
  // Shut down Modflow-USG, if active (necessary check?)
#ifdef _MODFLOW_USG_
  if (_nNodes > 0) { MFUSG::mf_shutdown(); }
#endif
}

//////////////////////////////////////////////////////////////////
/// \brief Initializes GW Model.
/// \details Is run AFTER InitMFUSG (requires values from reading MF model)
///
/// \param &Options [in] Global model options information
//
void CGroundwaterModel::Initialize(const optStruct  &Options) {
  // Allocate Modflow-USG variables
  _iStressPeriod = 0;
  _iStep = -1; // To be zero during first increment. Is there a smarter way to do this?
  _doSolve = 0;

  //-- Turn off Modflow-USG packages that Raven emulates
  MaskMFUSGpackages();

  //-- Setup first stress period
  InitializeStressPeriod();

  //-- Initialize PBJ package, if active
  _pRiver->InitializePBJ();

  //-- Initialize process LHS/RHS arrays
  _aProcessesAMAT = new double[_nNodes + 1];
  _aProcessesRHS  = new double[_nNodes + 1];
  for (int i = 1; i < _nNodes + 1; i++) {
    _aProcessesAMAT[i] = 0;
    _aProcessesRHS [i] = 0;
  }

  //-- Initialize Flux Tracker Array
  _aGWSWFluxes = new double[_pModel->GetNumHRUs()];
  for (int k=0; k < _pModel->GetNumHRUs(); k++)
  {
    _aGWSWFluxes[k] = 0.0;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Runs Fortran routine that reads in Modflow-USG model.
/// \details Currently is set to run during GW file parsing when name file is
///          known. Empty variables are passed and filled within routine.
///
/// \param namefile [in] filepath to MODFLOW name file
//
void CGroundwaterModel::InitMFUSG(string namefile)
{
#ifdef _MODFLOW_USG_
  // Jump through hoops to pass string to Fortran (more on the Fortran side, too!)
  int fileLen = namefile.size() + 1;
  char *cnamefile;
  cnamefile = new char[fileLen];
  strcpy(cnamefile, namefile.c_str());

  // Run Modflow USG initialize routines
  MFUSG::init_mfusg(cnamefile, &fileLen, &_nNodes, &_nLayers, &_nStressPeriods, &_maxiter, &_niunit);
  delete [] cnamefile;

  // Set default Recharge setting values
  SetRCHOptCode(3);          // See function, default is use top active node

  // Need to be contiguous arrays for Fortran
  _aNodesPerLayer      = (   int*)malloc(_nLayers        *sizeof(int   )); //new int[_nLayers];
  _aNumSteps           = (   int*)malloc(_nStressPeriods *sizeof(int   )); //new int[_nStressPeriods];
  _aStressPeriodLength = (double*)malloc(_nStressPeriods *sizeof(double)); //new double[_nStressPeriods];

  //-- Get memory from MFUSG
  MFUSG::get_mf_memory(_aNodesPerLayer, _aNumSteps, _aStressPeriodLength);
  for (int i = _nLayers; i > 0; i--) {
    _aNodesPerLayer[i] =  _aNodesPerLayer[i]- _aNodesPerLayer[i - 1]; //deaccumulate array
  }


  // Initialize recharge GW-SW process if HRU-data
  // Needs to know how many nodes (All in layer 1)
  // (unsure if I like this being here -LS)
  CGWRecharge* pRecharge = dynamic_cast<CGWRecharge*>(GetProcess(GWRECHARGE));
  if ((pRecharge != NULL) && (pRecharge->getRechargeType()==RECHARGE_HRU_DATA)){
    pRecharge->Initialize();
  }
#else
  ExitGracefully("CGroundwaterModel::InitMFUSG: cannot initialize ModflowUSG unless using Raven w/USG support", RUNTIME_ERR);
#endif
}

//////////////////////////////////////////////////////////////////
/// \brief Used to clear out the Raven-side AMAT/RHS accumulators
///  of entires for each node (1 to _nNodes)
void CGroundwaterModel::ClearMatrix()
{
  //-- Zero AMAT/RHS
  for (int n = 1; n < _nNodes + 1; n++) // To align with MFUSG, nodes are indexed from 1
  {
    _aProcessesAMAT[n] = 0.0;
    _aProcessesRHS [n] = 0.0;
  }
  //-- Zero fluxes
  for (int k=0; k < _pModel->GetNumHRUs(); k++)
  {
    _aGWSWFluxes[k] = 0.0;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Turns off MFUSG packages that Raven emulates
/// \details Sets IUNIT to zero in Modflow-USG to prevent certain packages from being active regardless of
///          if they're in the model namefile. Numbers correspond to MFUSGLib Extern.f90 CUNIT declaration
///          Currently turns off Drain and Recharge (see notes below on smart additions)
//
void CGroundwaterModel::MaskMFUSGpackages()
{
#ifdef _MODFLOW_USG_
  // Unit numbers for each package can be determined from CUNIT array in Extern.f90
  int nunits = 2;
  int off_units[2] = { 3, 8 };   // 3-DRN, 8-RCH // Eventual? 4-RIV, 18-STR, 22-LAK, 40-DRT, 44-SFR, 46-GAGE, 39-ETS,
  for (int i = 0; i < nunits; i++){
    MFUSG::mask_iunit(&off_units[i]);
  }
#endif
}

//////////////////////////////////////////////////////////////////
/// \brief Manages time steps between Raven and MFUSG
/// \details Handles changes in stress periods, which resets the MF time step to 1
///          Important for keeping the models in sync.
///
/// \param t [in] model time
//
void CGroundwaterModel::UpdateTStep(const double &t)
{
  double elapsed_from_start;
  // Eventually this method will need to account for differences in model stepping

  // Increment timestep counter
  _iStep++;

  // Check if we have a new stress period
  if (_iStep + 1 > _aNumSteps[_iStressPeriod]) {
    // Next stress period
    _iStressPeriod++;
    _iStep = 0;
    ExitGracefullyIf(_iStressPeriod > _nStressPeriods,
      "Raven has exceeded the number of stress periods in the Modflow Model.", RUNTIME_ERR);
    // Setup the new stress period
    InitializeStressPeriod();
  }

  // Setup new time step
  InitTS();

  // A little check to ensure Raven step matches MF
  // (Is this necessary? Wouldn't this be programmer error?)
  elapsed_from_start = _iStep;
  for (int i = 0; i < _iStressPeriod; i++) {
    elapsed_from_start += _aNumSteps[i];
  }
  ExitGracefullyIf(elapsed_from_start != t,
    "Modflow and Raven timesteps out of sync.", RUNTIME_ERR);
}

//////////////////////////////////////////////////////////////////
/// \brief Runs Fortran routine that reads in/sets up MFUSG active
///        package variables for the next stress period
//
void CGroundwaterModel::InitializeStressPeriod()
{
  // Fortran periods start at 1, ours start at zero
  // Quick conversion
  int mf_period = _iStressPeriod + 1;
#ifdef _MODFLOW_USG_
  MFUSG::init_stress_period(&mf_period);
#endif
}

//////////////////////////////////////////////////////////////////
/// \brief Runs Fortran routine that sets up MFUSG time step
///        (certain transient packages have additional routines)
//
void CGroundwaterModel::InitTS()
{
  // Fortran time steps start at 1, ours start at zero
  // Quick conversion
  int mf_step = _iStep + 1;
#ifdef _MODFLOW_USG_
  MFUSG::init_time_step(&mf_step, &_doSolve);
#endif
}

//////////////////////////////////////////////////////////////////
/// \brief Called every iteration to move Raven-side AMAT/RHS matrix entries
///        into the MFUSG solution.
/// \details Water/conductances will not be added/subtracted from inactive
/// nodes (i.e., IBOUND = 0). This is handled by MFUSGLib.
void CGroundwaterModel::InsertProcessesAsPackages()
{
#ifdef _MODFLOW_USG_
  double nodeAMAT=0.0;
  double nodeRHS =0.0;
  double area;
  // Loop over nodes

  ofstream outfile;
  for (int n = 1; n < _nNodes + 1; n++) // To align with MFUSG, nodes are indexed from 1
  {
    // Get stored values in GroundwaterModel AMAT/RHS variables
    nodeAMAT = _aProcessesAMAT[n];
    nodeRHS  = _aProcessesRHS [n];
    area        = GetNodeArea(n);
    if ((nodeAMAT == 0.0) && (nodeRHS == 0.0)) { continue; } // Cycles to next node
    // Add to MFUSG AMAT/RHS
    MFUSG::add_to_flow_eq(&n, &nodeAMAT, &nodeRHS);
  }
#endif
}

//////////////////////////////////////////////////////////////////
/// \brief Solves the groundwater flow equation, iterating until
///        convergence criteria is met
//
void CGroundwaterModel::Solve(const double &t)
{
  // Iteratively solve for new heads
#ifdef _MODFLOW_USG_

  UpdateTStep(t);

  int iter = 0, kiter = 0;  // Backtracking necessitates two tracking variables
  int ibflag = 0;
  _Converged  = false;

  while (iter < _maxiter) // Mimicks MF ancient Fortran control structure
  {
    kiter = iter + 1;
    MFUSG::formulate_boundaries(&kiter); // Fortran is a 1-index language

    // Add in Raven "Packages" contribution to MFUSG LHS/RHS
    InsertProcessesAsPackages();

    MFUSG::solve_matrix(&ibflag, &_Converged);
    // Check for backtracking
    if (ibflag == 1) { continue; }
    // Check for convergence
    if (_Converged) { break; }
    // No convergence, no backtracking = next iteration
    iter++;
  }
#endif
}
//////////////////////////////////////////////////////////////////
/// \brief Runs MFUSG routines run after a solution is found, adds in Raven budget entries
///  \details Primarly computes the groundwater flow budget given the new flow solution. Also
///           where the model crashes if no solution is found!
///
/// \param timestep [in] model timestep
//
void CGroundwaterModel::PostSolve(const double &timestep)
{
#ifdef _MODFLOW_USG_
  MFUSG::post_solution_updates();
  MFUSG::calc_boundary_budgets();
  // Add in Raven "Packages" contribution to MFUSG Flow Budget
  UpdateProcessBudgets(timestep);
  // Write out results
  MFUSG::print_save_data();

  // Was convergence not reached during this time step?
  ExitGracefullyIf(!_Converged, "Modflow-USG Failed to Converge", RUNTIME_ERR);
#endif
}

//////////////////////////////////////////////////////////////////
/// \brief Sends water budget (inflow/outflows) of GWSW Processes to MFUSG
/// \details Allows Raven GWSW Processes to have accurate representation in the
///          MFUSG List file water budgets for each time step. Calls MFUSG
///          function that imitates how a MFUSG Package would normally
///          report its water budget items, summing cell-by-cell entries.
///          Note rates passed on MFUSG via add_to_flow_budget should be RATES [L3/T]
///          The MFUSG routine multiplies these by the time change (DELT) internally
///
/// \param timestep [in] model timestep
//
void CGroundwaterModel::UpdateProcessBudgets(const double &timestep)
{
#ifdef _MODFLOW_USG_
  // Initialize variables to be filled by processes in loop
  int     HRUID, ibound, proc_nnodes, topnode, n;
  int    *proc_nodes;
  double *rates;
  double  hru_overlap_area, flow = 0.0;
  char   *procname;
  char    fluxname[16] = "RAVEN FLUX";
  vector<int> HRUnodes;

  // First, the Raven Groundwater Compartment Flux (non-GWSW process) volume
  proc_nodes  = new int   [_nNodes];
  rates       = new double[_nNodes];
  proc_nnodes = _nNodes;
  for (int i = 0; i < _nNodes; i++) {
    proc_nodes[i] = i+1;  //MFUSG uses 1-index
    rates     [i] = 0.0;  //initialize
  }
  // Loop over HRUs, fill/update arrays
  for (int k = 0; k < _pModel->GetNumHRUs(); k++)
  {
    HRUID = _pModel->GetHydroUnit(k)->GetID();
    HRUnodes = GetNodesByHRU(HRUID);
    hru_overlap_area = 0;

    for (int j=0; j < HRUnodes.size(); j++) {
      n = HRUnodes[j];
      // Get top active node - Raven Flux (recharge) is delivered to water table
      topnode = GetTopActiveNode(n);

      hru_overlap_area = GetNodeArea(n) * GetOverlapWeight(HRUID, n); // Get overlap weight only defined for top layer
      rates[topnode-1] += hru_overlap_area * _aGWSWFluxes[k];          // Rates is one off - node 1 is at rates[0]
    }
  }
  // Send to MFUSG
  MFUSG::add_to_flow_budget(&proc_nnodes, proc_nodes, rates, fluxname);
  delete [] proc_nodes;
  delete [] rates;

  // Loop over processes, running their "calcGWBudget" method.
  // Method fills the above arrays, allowing us to then pass them to MFUSG
  // _mGWSWProcesses is map of (process_type, CGWSWProcessABC).
  // using C++11 style map iterator https://stackoverflow.com/questions/26281979/c-loop-through-map
  for (auto const &proc : _mGWSWProcesses)
  {
    proc_nnodes = proc.second->getNumNodes();
    proc_nodes  = new int   [proc_nnodes];
    rates       = new double[proc_nnodes];

    proc.second->calcGWBudget(_pModel, proc_nodes, rates);
    procname = proc.second->getProcName();
    MFUSG::add_to_flow_budget(&proc_nnodes, proc_nodes, rates, procname);

    // Cleanup
    delete [] proc_nodes;
    delete [] rates;
  }
#endif
}

//////////////////////////////////////////////////////////////////
/// \brief Adds a hydrological process to GWModel
/// \details Adds a hydrological process that moves water to or from the GW model from a state variable.
/// Crashes if passed a Hydro Process
///
/// \param ptype    [in] process type (enum)
/// \param *newProc [in] pointer to GWSWProcess to be added to GW Model process map.
///
//
void CGroundwaterModel::AddProcess(process_type ptype, CHydroProcessABC* newProc)
{
  //-- Cast HydroProcessABC to CGWSWProcessABC
  CGWSWProcessABC* pGWSWProc = dynamic_cast<CGWSWProcessABC*>(newProc);

  ExitGracefullyIf(pGWSWProc == NULL, "CGroundwaterModel AddProcess:: Non-GWSW process (NULL)", BAD_DATA);

  _mGWSWProcesses.insert(pair<process_type, CGWSWProcessABC*>(ptype, pGWSWProc));
  _nProcesses++;
}

//////////////////////////////////////////////////////////////////
/// \brief Adds to AMAT/RHS arrays to be passed to MFUSG
/// \details Adds head_coefficent and flow_rate to the LHS (AMAT) and RHS GW Flow Equation,
///  respectively, for the specified node
///  Does not support off-diagonal AMAT coefficients (e.g., connection between n and n+1)
///  For a complete explanation about formulating boundaries with these equations,
///  see the documentation for add_to_flow_eq() in Extern.f90 in the MFUSGLib
///
/// \param n                [in] Node ID (valid 1 to _nNodes)
/// \param head_coefficient [in] Addition to location n in the coefficient (AMAT) matrix
/// \param flow_rate        [in] Addition to location n in the RHS matrix
///
/// LS NOTE 2021 - Consider an error for if n=0. SHOULD BE INDEXED FROM 1!
/// That's how it's passed to MFUSG
void CGroundwaterModel::AddToGWEquation(int n, double head_coefficient, double flow_rate)
{
  _aProcessesAMAT[n] += head_coefficient;
  _aProcessesRHS [n] += flow_rate;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of active GW SW Processes in the GW Model
//
int CGroundwaterModel::GetNumProcesses  () const { return _nProcesses; }

//////////////////////////////////////////////////////////////////
/// \brief Returns specific hydrologic process denoted by type
///
/// \param ptype    [in] process type (enum)
/// \return pointer to GWSW Process requested
//
CGWSWProcessABC* CGroundwaterModel::GetProcess(process_type ptype)
{
  CGWSWProcessABC* pReturn = NULL;
  if (_mGWSWProcesses.find(ptype) != _mGWSWProcesses.end()) {
    pReturn = _mGWSWProcesses.find(ptype)->second;
  }
  return pReturn;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns hydrologic process map indexed by type
///
map<process_type, CGWSWProcessABC*>* CGroundwaterModel::GetProcessMap() { return &_mGWSWProcesses; }

//////////////////////////////////////////////////////////////////
/// \brief Returns GWRiverConnection for GW Model
///
/// \return pointer GW River Connection object
//
CGWRiverConnection* CGroundwaterModel::GetRiverConnection() {return _pRiver;}

//////////////////////////////////////////////////////////////////
/// \brief Returns drain process pointer
///
/// \return pointer to drain gwsw process (error if doesn't exit)
//
CGWDrain* CGroundwaterModel::GetDrainProcess()
{
  CGWDrain* pDrain = dynamic_cast<CGWDrain*>(GetProcess(DRAIN));
  ExitGracefullyIf(pDrain == NULL, "ParseGWFile:: No drain processes in model", BAD_DATA);
  return(pDrain);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns recharge process pointer
///
/// \return pointer to recharge gwsw process (error if doesn't exit)
//
CGWRecharge* CGroundwaterModel::GetRechProcess()
{
  CGWRecharge* pRecharge = dynamic_cast<CGWRecharge*>(GetProcess(GWRECHARGE));
  ExitGracefullyIf(pRecharge == NULL, "ParseGWFile:: No Recharge processes in model", BAD_DATA);
  return(pRecharge);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns total number nodes in the GW model
///
int CGroundwaterModel::GetNumNodes() { return _nNodes; }

//////////////////////////////////////////////////////////////////
/// \brief Returns number nodes in the given layer of the groundwater model (1-index)
//
int CGroundwaterModel::GetNumNodesByLay(int layer)
{
  return _aNodesPerLayer[layer-1];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns total number of time steps set up in the GW model
///
int CGroundwaterModel::GetTotalTSteps()
{
  int total_steps = 0;
  for (int i = 0; i < _nStressPeriods; i++) {
    total_steps += _aNumSteps[i];
  }
  return total_steps;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns MFUSG head for the given node
double CGroundwaterModel::GetNodeHead(int n)
{
#ifdef _MODFLOW_USG_
  return (MFUSG::get_node_head(&n));
#else
  return 0;
#endif
}

//////////////////////////////////////////////////////////////////
/// \brief Returns MFUSG area for the given node
double CGroundwaterModel::GetNodeArea(int n)
{
#ifdef _MODFLOW_USG_
  return(MFUSG::get_node_area(&n));
#else
  return 0.0;
#endif
}

//////////////////////////////////////////////////////////////////
/// \brief Sets the weight value for the OverlapWeight sparse matrix [HRUID][node]
///
/// \param HRUID    [in] HRU identifier
/// \param node     [in] node id (1 to _nNodes)
/// \param weight   [in] Percentage of node that is covered by the HRU corresponding to HRUID
/// \param &Options [in] Global model options information
///
void CGroundwaterModel::SetOverlapWeight(const int HRUID, const int node, const double weight, const optStruct& Options)
{
#ifdef _MODFLOW_USG_
  bool newConnection = true;
  //-- Check if already exists (Should this be moved to the parser?)
  pair<int,int> p(HRUID,node);
  if (_mOverlapWeights.find(p) != _mOverlapWeights.end())
  {
    newConnection = false;
    string warn = "GroundwaterModel: Overlap Weight for Cell " + to_string(node) + " with HRU " + to_string(HRUID) + "defined multiple times.";
    WriteWarning(warn, Options.noisy);
  }
  //-- Overlap Weight Matrix currently implemented using maps(!)
  _mOverlapWeights[p] = weight;

  //-- Also add to HRUsByNode and NodeByHRU maps (if new!)
  if (newConnection)
  {
    _mHRUsByNode[node].push_back(HRUID);
    _mNodesByHRU[HRUID].push_back(node);
  }
#endif
};

//////////////////////////////////////////////////////////////////
/// \brief Gets the weight value for the OverlapWeight sparse matrix.
/// \details Returns 0.0 if no weight was exists. Only contains nodes
/// in layer 1. The intent for nodes in the lower layers is to use
/// GetTopActiveNode() for the correct node to move water to, but use
/// the layer 1 for finding HRU connections.
///
/// \param HRUID    [in] HRU identifier
/// \param node     [in] node id (1 to _aNodesPerLayer[0])
///
double CGroundwaterModel::GetOverlapWeight(const int HRUID, const int node)
{
  //-- Error checking for nodes in deeper layers
  // (Partially to help LS find where loops go past lay 1)
  if (node > _aNodesPerLayer[0]) {
    ExitGracefully("CGroundwaterModel::GetOverlapWeight: node-HRU overlap requested outside of the top layer (layer 1)", RUNTIME_ERR);
  }

  //-- Overlap Weight Matrix currently implemented using maps(!)
  map<pair<int, int>, double>::iterator weight_iter;    // <HRUid, node>, weight
  double weight = 0.0;

  //-- Use map.find() to get iterator to the element
  pair<int,int> p(HRUID,node);
  weight_iter = _mOverlapWeights.find(p);

  //-- Get weight if exists, otherwise return default (0.0)
  if (weight_iter != _mOverlapWeights.end()) {
    weight = weight_iter->second;
  }

  //-- Return weight (double)
  return weight;
};

//////////////////////////////////////////////////////////////////
/// \brief Returns a pointer to the Overlap Weight "matrix"
/// \return OverlapWeights pointer
///
map<pair<int, int>, double> *CGroundwaterModel::GetOverlapWeights()
{
  return &_mOverlapWeights;
}

//////////////////////////////////////////////////////////////////
/// \brief Get the HRUs that overlaps with a specific cell
/// \return vector of HRU ids
///
/// \param node     [in] node id (1 to _nNodes)
///
vector<int> CGroundwaterModel::GetHRUsByNode(const int node)
{
  vector<int> HRUs;
  //-- Obtain vector of HRUIDs for given cell
  map<int, vector<int> >::iterator map_query;
  map_query = _mHRUsByNode.find(node);

  //-- If not empty set to the return value
  if (map_query != _mHRUsByNode.end()) {
    HRUs = map_query->second;
  }

  return HRUs;
}

//////////////////////////////////////////////////////////////////
/// \brief Get the cells that overlap with a specific HRU
/// \return vector of node ids
///
/// \param HRUID    [in] HRU identifier
///
vector<int> CGroundwaterModel::GetNodesByHRU(const int HRUID)
{
  vector<int> Cells;
  //-- Obtain vector of nodes given HRU
  map<int, vector<int> >::iterator map_query;
  map_query = _mNodesByHRU.find(HRUID);

  //-- If not empty set to the return value
  if (map_query != _mNodesByHRU.end()) {
    Cells = map_query->second;
  }

  return Cells;
}

//////////////////////////////////////////////////////////////////
/// \brief Sets MFUSG Recharge option code, determinding where recharge gets delivered
/// \details A value of 1 sends recharge always to layer 1, 3 sends it to the top active
/// node in the vertical column. Used every time GetTopActiveNode is called.
/// Could potentially be used for GWSW Processes. See MFUSG I/O manual for more details.
///
/// \param code     [in] Value to set MODFLOW variable NRCHOP
///
void CGroundwaterModel::SetRCHOptCode(int code) {
#ifdef _MODFLOW_USG_
  MFUSG::set_nrchop(&code);
#endif
}

//////////////////////////////////////////////////////////////////
/// \brief Runs MF-USG function to find highest (top most layer) node, given a node # in the top layer
/// \return top active node
/// \details Controls recharge - see SetRCHOptCode & MFUSG::set_NRCHOP
/// Could potentially be used for GWSW Processes
///
/// \param node     [in] node id (1 to _nNodes)
///
const int CGroundwaterModel::GetTopActiveNode(const int node)
{
  int topnode = node;
#ifdef _MODFLOW_USG_
  MFUSG::first_active_below(&topnode);
#endif
  return topnode;
}

//////////////////////////////////////////////////////////////////
/// \brief Adds a flux to the GWSWFlux counter
/// \details Called by GWSW Processes to keep track of how much water they add to the HRU
/// Groundwater storage compartment.
///
/// \param k     [in] HRU global index
/// \param moved [in] amount of water [mm] moved by GWSW process to the HRU GW compartment
///
/// LS Note: Maybe this should actually be a method that loops over GWSW Processes in _mGWSWProcesses
/// and calls a method belonging to GWSW Process to return their GW flux.
/// The goal of such a change would be to make GWSW Processes not require a pointer to the model
/// (independent like HydroProcesses)
void CGroundwaterModel::AddFlux(const int k, const double moved)
{
  _aGWSWFluxes[k] += moved;
}

//////////////////////////////////////////////////////////////////
/// \brief Adds the non-GWSW Process GW flux (e.g., HydroProcess flux) to the Raven GW equation accumulators
/// \details Takes the HRU groundwater compartment current value, subtracts the contributions
/// from GWSW Processes (tracked in _GWSWFluxes for each HRU, see CGroundwaterModel::AddFlux and
/// CGWSWProcessABC::SendRateToFluxTracker) and adds that to the LHS and RHS arrays Raven keeps
/// for adding during MF iterations.
/// The separation is important for tracking and/or updating fluxes to/from GWSW Process-connected compartments
///
/// \param *pHRU    [in] Pointer to HRU
/// \param GWVal    [in] Depth of water [mm] in HRU GW Compartment
///
void CGroundwaterModel::FluxToGWEquation(const CHydroUnit *pHRU, double GWVal)
{
  long long int HRUID;
  int k;
  double      flux, weight, area, node_rech;
  int         n,active_node;
  vector<int> HRUnodes;

  k               = pHRU->GetGlobalIndex();
  HRUID           = pHRU->GetHRUID();
  HRUnodes        = GetNodesByHRU((int)(HRUID)); //\todo[funct] - support long

  //-- Calc Non-GWSW Process contribution    //testing flux = 9.027605E-04;
  flux = (GWVal - _aGWSWFluxes[k]) / MM_PER_METER; // [mm/T] to [m/T]

  //-- This flux should never be below zero, right? Something to think about [checking for]
  //-- Distribute flow among HRU cells
  for (unsigned int j=0; j< HRUnodes.size(); j++)
  {
    // Correct to topmost active node
    n           = HRUnodes[j];
    active_node = GetTopActiveNode(n);
    weight      = GetOverlapWeight((int)(HRUID), n); //todo[funct] - support long long int
    area        = GetNodeArea(n);
    node_rech   = flux * area * weight; // [m3/d] * overlap weight-corrected
    // GW Rate
    AddToGWEquation(active_node, 0.0, node_rech);
  }

  //-- Update HRU Flux for later GW Budget update
  _aGWSWFluxes[k] = flux;
}

//////////////////////////////////////////////////////////////////
/// \brief Set Raven Cell-by-Cell Unit Number for MODFLOW-USG
void CGroundwaterModel::SetCBCUnitNo(int unit_no)
{
#ifdef _MODFLOW_USG_
  MFUSG::set_irvncb(&unit_no);
#endif
}
