/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2022 the Raven Development Team
----------------------------------------------------------------*/
#ifndef ENKF_H
#define ENKF_H

#include "ModelEnsemble.h"

struct force_perturb 
{
  forcing_type forcing;
  disttype     distribution;
  double       distpar[3];   ///< distribution parameters - conditional on distribution
  int          kk;           //< HRU group index (or DOESNT_EXIST if perturbation should apply everywhere)
  adjustment   adj_type;     //< additive or multiplicative adjustment
};
struct obs_perturb
{
  sv_type      state;
  disttype     distribution;
  double       distpar[3];   ///< distribution parameters - conditional on distribution
  int          kk;           //< HRU group index (or DOESNT_EXIST if perturbation should apply everywhere)
  adjustment   adj_type;     //< additive or multiplicative adjustment
};
////////////////////////////////////////////////////////////////////
/// \brief Data abstraction for EnKF model ensemble run
//
class CEnKFEnsemble : public CEnsemble
{
private:
  int            _nEnKFMembers;     ///< number of members in ensemble (=1/2 of _nMembers, because pre/post runs are separate)

  double       **_state_matrix;     ///< current array of state vectors, X [size: _nEnKFMembers x_nStateVars]
  int            _nStateVars;       ///< total number of assimilated state variables

  double       **_obs_matrix;       ///< matrix of perturbed observations, D [size: _nEnKFMembers x _nObsDatapoints]
  double       **_output_matrix;    ///< matrix of simulated values corresponding to observations, HA [size: _nEnKFMembers x _nObsDatapoints]
  double       **_noise_matrix;     ///< matrix of observational noise [size: _nEnKFMembers x _nObsDatapoints]
  int            _nObsDatapoints;   ///< number of valid datapoints available for assimilation

  force_perturb**_pPerturbations;   ///< array of pointers to perturbation data; defines which forcing functions to perturb and how [size: _nPerturbations]
  int            _nPerturbations;   ///< number of forcing functions to perturb 

  obs_perturb  **_pObsPerturbations;///< array of pointers to observation perturbation data [size _nObsPerturbations] 
  int            _nObsPerturbations;///< number of observational perturbations. If observation does not have perturbation, it is assumed "perfect" data

  sv_type       *_aAssimStates;     ///< assimilated state variable [size: _nAssimStates]
  int           *_aAssimLayers;     ///< layer of assimilated state variable [size: _nAssimStates] 
  int           *_aAssimGroupID;    ///< index (kk) of HRU group or SubBasin group  [size: _nAssimStates]
  int            _nAssimStates;     ///< number of assimilated state variable types (!= number of total state vars)

  double         _t_assim_start;    //< assimilation start time 

  int            _data_horizon;     //< number of data points back in time to be used in assimilation (=1 for classic EnKF, more for variational methods)
  int            _nTimeSteps;       //< number of timesteps in hindcast period (prior to t0)

  bool           _warm_ensemble;    //< true if initial conditions should be read from previous ensemble run (default: false)
  string         _warm_runname;     //< run name of warm solution files (e.g., run1_solution.rvc)

  int           *_aObsIndices;      //< indices of CModel::pObsTS array corresponding to assimilation time series [size:_nObs]
  int            _nObs;             //< number of time series to be assimilated

  ofstream      _ENKFOUT;           ///< output file stream

  void AssimilationCalcs();         //< determines the final state matrix after assimilation
  void UpdateStates    (CModel *pModel,optStruct& Options,const int e);
  void AddToStateMatrix(CModel* pModel,optStruct& Options,const int e);
public:
  CEnKFEnsemble(const int num_members,const optStruct &Options);
  ~CEnKFEnsemble();

  double GetStartTime(const int e) const;

  void SetToWarmEnsemble(string runname);
  void SetDataHorizon(const int nTimesteps);
  void AddPerturbation   (forcing_type type, disttype distrib, double *distpars, int group_index, adjustment adj);
  void AddObsPerturbation(sv_type      type, disttype distrib, double *distpars, adjustment adj);
  void AddAssimilationState(sv_type sv, int layer, int assim_groupID);

  void Initialize       (const CModel* pModel,const optStruct &Options);
  void UpdateModel      (CModel *pModel,optStruct &Options,const int e);
  void StartTimeStepOps (CModel* pModel,optStruct& Options,const time_struct &tt,const int e);
  void CloseTimeStepOps (CModel* pModel,optStruct& Options,const time_struct& tt,const int e); //called at end of each timestep
  void FinishEnsembleRun(CModel *pModel,optStruct &Options,const int e);
};
#endif

