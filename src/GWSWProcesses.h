/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2021 the Raven Development Team
----------------------------------------------------------------*/

#ifndef GWSWPROCESSES_H
#define GWSWPROCESSES_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"
#include "GroundwaterModel.h"

#include <set>
#include <vector>
#include <map>

class CGroundwaterModel;

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for groundwater-surface water processes (abstract base class)
/// \details Data abstraction for processes that move water between the groundwater
///   model and the surface water model. Inherits HydroProcessABC values and methods.
///   The important difference between the classes is that GWSWProcess has access to the GW
///   Model object and the HRU-cell overlap weight matrix
///
///   GWSWProcesses have some additional requirements within the GetRatesOfChange() method:
///    - Must calculate flux to individual cells, passed to the GW model through the GW model
///      method AddToGWEquation()
///    - Need to pass rate of change (+ = to aquifer) to GWModel AddFlux() using
///      SendRateToFluxTracker() method in this class
///
///   Other important differences are the calcGWBudget method, which is used
///   to inform the GW model of the total flux of the process, and the UpdateRatesOfChange
///   method, which can be used to correct head-dependent SV fluxes after the model has
///   run (experimental).
//
class CGWSWProcessABC : public CHydroProcessABC
{
protected:/*-------------------------------------------------------*/
	CGroundwaterModel *pGWModel;       ///< pointer to groundwater model

	int          _nnodes;                     ///< Number of gw model nodes with process
	int         *_inodes;                     ///< array of nodes with process
  char        *_pProcName;                  ///< Name of process to be displayed in MFUSG Budget

	map<int, set<int> >  _mNodesByHRU;         ///< Map of the set of nodes with process, keyed by HRU
	map<int, int>        _mNodeIndex;           ///< Hash to easily find index of a node within _inode [nodeid, index]

public: /*-------------------------------------------------------*/
	// Constructors/Destructor
  CGWSWProcessABC(CGroundwaterModel  *pGWModel,
                  const process_type ptype,
                  CModelABC          *pModel);
	virtual ~CGWSWProcessABC();

	virtual void Initialize();

  // Main Methods
  virtual void calcGWBudget(const CModel *pModel, int *nodes, double *rates);
  virtual void UpdateRatesOfChange(const double      *state_vars,
                                   const CHydroUnit  *pHRU,
                                   const optStruct   &Options,
                                   const time_struct &tt,
                                         double      *rates) const;
  void SendRateToFluxTracker(const int k, const double *rates) const;

  // Setters
  void setProcName(string name);

	// Getters
	int   getNumNodes ();
  char* getProcName ();
  set<int>  getHRUnodes(int HRUid) const;

};

/*****************************************************************
                             Drain
------------------------------------------------------------------
*****************************************************************/

////////////////////////////////////////////////////////////////////
/// \brief Data abstraction for groundwater rising above a given elevation and leaving the subsurface
//
class CGWDrain : public CGWSWProcessABC
{
private:/*------------------------------------------------------*/
	double    *_aElevations;
	double    *_aConductances;
  double    *_aHeads;

public:/*-------------------------------------------------------*/
	// Constructors/destructors:
  CGWDrain(CGroundwaterModel *pGWM,
           CModel            *pModel);
	~CGWDrain();

	// HydroProcess inherited functions
	void        InitializeDrainClass        (int nDrains);
	void        GetRatesOfChange            (const double		  *state_vars,
		                                       const CHydroUnit  *pHRU,
		                                       const optStruct	  &Options,
		                                       const time_struct &tt,
		                                       double            *rates) const;
	void        ApplyConstraints            (const double      *state_vars,
		                                       const CHydroUnit  *pHRU,
		                                       const optStruct	  &Options,
		                                       const time_struct &tt,
		                                       double            *rates) const;
	void        GetParticipatingParamList   (string  *aP , class_type *aPC , int &nP) const;
	static void GetParticipatingStateVarList(sv_type *aSV, int *aLev, int &nSV);

  // GWSWProcess Inherited Functions
  void calcGWBudget       (const CModel *pModel, int *nodes, double *rates);
  void UpdateRatesOfChange(const double      *state_vars,
                           const CHydroUnit  *pHRU,
                           const optStruct   &Options,
                           const time_struct &tt,
                                 double      *rates) const;

	// Drain Setters
	void addDrain(int drainID, int nodeID, double drain_elev, double drain_cond);

	// Drain Getters
	// getElev(int drainID);
	// getCond(int drainID);

};

/*****************************************************************
                            Recharge
------------------------------------------------------------------
*****************************************************************/

enum gwrecharge_type {
  RECHARGE_HRU_DATA,     // RECHARGE_DATA in file
  RECHARGE_FLUX
};

////////////////////////////////////////////////////////////////////
/// \brief Data abstration of loss of water to groundwater
//
class CGWRecharge : public CGWSWProcessABC
{
private:/*------------------------------------------------------*/
  gwrecharge_type  _recharge_type;
  double           _totalRecharge;
  double          *_aFluxes;

public:/*-------------------------------------------------------*/
  // Constructors/destructors:
  CGWRecharge(CGroundwaterModel *pGWM,
              gwrecharge_type    type,
              int                From_index,
              CModel            *pModel);          // Needs to set iFrom
  ~CGWRecharge();

  // HydroProcess inherited functions
  void InitializeRechargeClass(int nNodes);
  void GetRatesOfChange(const double		  *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct	  &Options,
                        const time_struct &tt,
                        double      *rates) const;
  void ApplyConstraints(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct	  &Options,
                        const time_struct &tt,
                        double      *rates) const;
  void        GetParticipatingParamList(string  *aP, class_type *aPC, int &nP) const;
  static void GetParticipatingStateVarList(sv_type *aSV, int *aLev, int &nSV);

  // GWSWProcess Inherited Functions
  void calcGWBudget(const CModel *pModel, int *nodes, double *rates);

  // Recharge Setters
  void     addNodeFlux      (int index, int nodeID, double flux);

  // Recharge Getters
  gwrecharge_type getRechargeType();
};

#endif
