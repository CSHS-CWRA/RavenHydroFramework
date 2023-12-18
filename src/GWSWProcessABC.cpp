/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2021 the Raven Development Team
----------------------------------------------------------------*/

#include "GWSWProcesses.h"  //Declared here
#include "GroundwaterModel.h"

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the constructor for an abstract GW-SW process
/// \param pGWModel [in] pointer to groundwater model
/// \param ptype    [in] process type (enum)
//
CGWSWProcessABC::CGWSWProcessABC(CGroundwaterModel *pGWM,
                                 const process_type ptype,
                                 CModelABC *pModel)
  :CHydroProcessABC(ptype, pModel)
{
  pGWModel = pGWM;
  _nnodes  = 0;
  _inodes  = NULL;
  _pProcName = new char[16];
  setProcName("RAVEN DEFAULT");
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor for an abstract GW-SW process
//
CGWSWProcessABC::~CGWSWProcessABC()
{
  delete [] _inodes;
  delete [] _pProcName;
}

//////////////////////////////////////////////////////////////////
/// \brief Initializes process, generates Node-HRU maps
/// \details _inodes, _nnodes must be populated already when this runs
//
void CGWSWProcessABC::Initialize()
{
  CHydroProcessABC::Initialize();
  int nodeid;
  vector<int> HRUs;
  //-- Need to populate NodesByHRU map
  for (int i = 0; i < _nnodes; i++)
  {
    //-- Get vector of HRUs connected to node
    nodeid = _inodes[i];
    HRUs = pGWModel->GetHRUsByNode(nodeid);
    for (vector<int>::iterator HRU = HRUs.begin(); HRU != HRUs.end(); ++HRU) {
      _mNodesByHRU[*HRU].insert(_inodes[i]);
      _mNodeIndex[_inodes[i]] = i;
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns the number of GW model nodes involved in the process
//
int CGWSWProcessABC::getNumNodes() { return(_nnodes); }

//////////////////////////////////////////////////////////////////
/// \brief Returns character array of GW SW Process name (for MFUSG)
//
char *CGWSWProcessABC::getProcName() { return(_pProcName); }

//////////////////////////////////////////////////////////////////
/// \brief Returns set of nodes, for the given HRU, that have the GW SW Process
/// \return Set of nodes to iterate over, empty if none in specified HRU
///
/// \param HRUid [in] HRU identifier
//
set<int>  CGWSWProcessABC::getHRUnodes(int HRUid) const
{
  if (_mNodesByHRU.find(HRUid) != _mNodesByHRU.end()) {
    return(_mNodesByHRU.find(HRUid)->second);
  }
  else {
    return(set<int>()); // Empty
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets GW SW Process name, performs string to char conversion for MFUSG
///
/// \param name [in] the name of the process [max 16 characters] to identify it in the MFUSG budget
//
void CGWSWProcessABC::setProcName(string name)
{
  ExitGracefullyIf(name.length() > 16, "GWSW Process Name Too Long", RUNTIME_ERR);
  strcpy(_pProcName, name.c_str());
}

//////////////////////////////////////////////////////////////////
/// \brief Fills arrays[_nnodes] with process rates to(+)/from(-) the GW model.
///        Placeholder method for child processes to override.
///
/// \param pModel [in] pointer to groundwater model
/// \param nodes [out] array of nodes where process is active
/// \param rates [out] rates of flux to GW model
//
void CGWSWProcessABC::calcGWBudget(const CModel *pModel, int *nodes, double *rates)
{
  // Fill with zeros, this is a placeholder method for child processes
  for (int i=1; i < _nnodes; i++) {
    nodes [i] =  _inodes[i];
    rates [i] = 0.0;
  }
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
void CGWSWProcessABC::UpdateRatesOfChange(const double      *state_vars,
                                          const CHydroUnit  *pHRU,
                                          const optStruct   &Options,
                                          const time_struct &tt,
                                                double      *rates) const
{
  for (int q=0;q<_nConnections;q++)
  {
    rates[q]=0.0;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Finds rate(s) of change being send to the GW SV and adds to the GW Model flux tracker
/// \details The Groundwater SV ends up containing both HydroProcess and GWSWProcess fluxes
///    However, this results in a double counting since GWSWProcesses also add to the AMAT and RHS
///    arrays of the GW Model. Therefore, the fluxes added by GWSWProcesses need to be tracked so
///    they can be removed from the GW SV prior to using it to calculate net recharge.
///    This method makes calling GroundwaterModel::AddFlux() easier by automatically figuring out
///    the appropriate rate and sign to use for the passed flux (+ = to aquifer)
///
/// \param k      [in]  HRU global index
/// \param *rates [out] Rates of water/energy moved between storage locations / rates of change of modified state variables
//
void CGWSWProcessABC::SendRateToFluxTracker(const int k, const double *rates) const
{
  int GWSVID = pModel->GetStateVarIndex(GROUNDWATER);
  for (int q=0;q<_nConnections;q++)
  {
    if (iFrom[0] == GWSVID) {
      pGWModel->AddFlux(k, -1*rates[0]);
    }
    else if (iTo[q] == GWSVID) {
      pGWModel->AddFlux(k, rates[0]);
    }
  }
}
