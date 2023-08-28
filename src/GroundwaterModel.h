/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team, Leland Scantlebury
------------------------------------------------------------------
  Class CGroundwater
----------------------------------------------------------------*/
#ifndef GROUNDWATER_H
#define GROUNDWATER_H

#include <vector>
#include <map>
#include "RavenInclude.h"
#include "Model.h"

class CModel;
class CGWSWProcessABC;
class CGWDrain;
class CGWRecharge;
class CGWRiverConnection;

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for groundwater model
///
/// A note about HRU-node connections:
/// It is assumed, primarily, Raven will be interacting with the top layer of
/// the groundwater model. For this reason, the Node-HRU maps only define the
/// the connections to nodes in the first layer (_mNodesByHRU, _mHRUsByNode)
/// This is also true for _mOverlapWeights. The assumption is that any process
/// connected to lower nodes shouldn't be needing weights (e.g., all the water
/// goes to a specific node). This could change in the future.
/// For consistency with MODFLOW, processes moving water to the groundwater
/// table can check if the node they are sending water to is the top
/// active node using GetTopActiveNode() (see usage in UpdateProcessBudgets()).
/// This is further mediated by the recharge option settings (NRCHOP,
/// SetRCHOptCode) which defaults to 3 - "top active"
/// Thus, when finding the correct node to apply recharge to use the original
/// Layer 1 node for queries into HRU overlap calculations.
/// It's an unstructured model - figuring out the node above/below another is
/// not always arbitrary! (see rvn_first_active_below in Extern.f)
///
class CGroundwaterModel
{
private:/*----------------------------------------------------*/

  CModel                             *_pModel;
  CGWRiverConnection                 *_pRiver;

  int                                 _nNodes;  ///< Nodes in the groundwater model
  int                                _nLayers;  ///< number of layers in the groundwater model
  int                       * _aNodesPerLayer;  ///< Array of number of nodes per layer [size:_nLayers]

  int                         _nStressPeriods;  ///< Total # of Stress periods
  int                          _iStressPeriod;  ///< Current stress period (kper in MFUSG)

  int                             *_aNumSteps;  ///< Number of steps in each stress period [size:_nStressPeriods]
  int                                  _iStep;  ///< Current time step (kstp in MFUSG)

  int                                _maxiter;  ///< Maximum outer interations for GW Solver
  int                                 _niunit;  ///< Number of items in the Modflow IUNIT array
  int                                _doSolve;  ///< Equivalent to MFUSG var "idoflow" - 0 means don't solve (Steady state or transport step).
  bool                             _Converged;  ///< Model Convergence success (from MFUSG)

  int                             _nProcesses;  ///< GW-SW Processes in Raven
  map<process_type,
    CGWSWProcessABC*>         _mGWSWProcesses;  ///< map of pointers to GWSW Processes with type as key (IMPORTANT: are also in Model _pProcesses).

  double                *_aStressPeriodLength;  ///< Length of each stress period [d] (perlen in MFUSG)
  double                     *_aProcessesAMAT;  ///< Array of conductance terms to be added to MFUSG LHS, by node [size _nNodes+1] (to index from 1)
  double                      *_aProcessesRHS;  ///< Array of volumes to added to MFUSG RHS, by node [size _nNodes+1] (to index from 1)
  double                        *_aGWSWFluxes;  ///< Sum of fluxes (one for each HRU) to/from the GW compartment.
                                                ///  Used to determine non-GWSW process flux to GW model

  map<pair<int,int>, double> _mOverlapWeights;  ///< 2D "matrix" of weights relating which HRUs are connected to which GW cells
                                                ///  =A_ki/A_i where A_ki is area of overlap between HRU k and node i and A_i is area of node i
                                                ///  Only contains weights for layer 1 (see HRU-node comment at top)
                                                ///  Dimensions : [_nHydroUnits][_aNodesPerLayer[0]] (variable)
  map<int,vector<int> >          _mNodesByHRU;  ///< Map of length _nNodes with vectors of the HRUs connected to that node
  map<int,vector<int> >          _mHRUsByNode;  ///< Map of length _nHRU with vectors of the nodes connected to that HRU

  void                   UpdateTStep         (const double &t);
  void                   UpdateProcessBudgets(const double& timestep);

public:/*-------------------------------------------------------*/

  //-- Constructors:
  CGroundwaterModel(CModel *pModel);
  ~CGroundwaterModel();
  void                   Initialize(const optStruct &Options);

  //-- General Methods
  void                   ClearMatrix();

  //-- MFUSG Running Methods
  void                   InitMFUSG                (string namefile);
  void                   MaskMFUSGpackages        ();

  void                   InitializeStressPeriod   ();
  void                   InitTS                   ();
  void                   InsertProcessesAsPackages();
  void                   Solve                    (const double &t);
  void                   PostSolve                (const double &timestep);  // Revise miserable undescriptive name
  const int              GetTopActiveNode         (const int node);

  //-- Setters
  void                   AddProcess       (process_type ptype, CHydroProcessABC *pGWSWProc);
  void                   AddToGWEquation  (int n, double head_coefficient, double flow_rate); // added to AMAT and RHS, respectively.
  void                   SetOverlapWeight (const int        HRUID,
                                           const int        node,
                                           const double     weight,
                                           const optStruct &Options);
  void                   AddFlux          (const int k, const double moved);
  void                   FluxToGWEquation (const CHydroUnit *pHRU, double GWVal);
  void                   SetCBCUnitNo     (int unit_no);
  void                   SetRCHOptCode    (int code);

  //-- Getters
  int                    GetNumNodes      ();
  int                    GetNumNodesByLay (int layer);
  double                 GetNodeHead      (int n);
  double                 GetNodeArea      (int n);
  int                    GetTotalTSteps   ();
  double                 GetOverlapWeight (int HRUID, int node);
  int                    GetNumProcesses  () const;
  CGWSWProcessABC*       GetProcess       (process_type ptype);
  map<process_type,
    CGWSWProcessABC*>*   GetProcessMap    ();
  CGWRiverConnection*    GetRiverConnection ();
  CGWDrain*              GetDrainProcess  ();
  CGWRecharge*           GetRechProcess   ();
  map<pair<int,int>,
    double>             *GetOverlapWeights();
  vector<int>            GetHRUsByNode    (const int cellid);
  vector<int>            GetNodesByHRU    (const int HRUID);
};

#endif
