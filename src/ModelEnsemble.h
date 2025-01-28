/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2023 the Raven Development Team
----------------------------------------------------------------*/
#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include "RavenInclude.h"
#include "Model.h"
#include "SoilAndLandClasses.h"


struct param_dist
{
  //this structure defines a parameter distribution
  string       param_name;   ///< parameter name
  class_type   param_class;  ///< parametre type e.g., CLASS_SOIL
  string       class_group;  ///< specific class (e.g., SILTY_SAND) this refers to
  double       default_val;  ///< base parameter value
  disttype     distribution; ///< statistical distribution
  double       distpar[3];   ///< distribution parameters
  //                         //   for DIST_UNIFORM, distpar[0]=min, distpar[1]=max
  //                         //   for DIST_NORMAL, distpar[0]=mean, distpar[1]=std_dev
  //                         //   for DIST_NORMAL, distpar[0]=mean of ln, distpar[1]=std_dev of ln
  //                         //   for DIST_GAMMA, distpar[0]=shape, distpar[1]=scale
  //transformation trans; e.g., log transform

};
double SampleFromDistribution(disttype distribution,double distpar[3]);

////////////////////////////////////////////////////////////////////
/// \brief Data abstraction for model ensemble run
/// \details Stores ensemble member details and functionality on how to modify each ensemble member
//
class CEnsemble
{
protected:/*------------------------------------------------------*/
  ensemble_type _type;           ///< ensemble type (e.g., ENSEMBLE_MONTECARLO)

  int           _nMembers;       ///< number of ensemble members

  string       *_aOutputDirs;    ///< array of output directory names [size: _nMembers]
  string       *_aRunNames;      ///< array of output runnames [size: _nMembers]
  string       *_aSolutionFiles; ///< array of input solution filenames [size: _nMembers]

  int           _rand_seed;      ///< random seed

  bool          _disable_output; ///< true if output from ensemble should be turned off (default: false)

public:/*-------------------------------------------------------*/
  CEnsemble(const int num_members, const optStruct &Options);
  ~CEnsemble();

  //Acessor Function
  int            GetNumMembers() const;
  ensemble_type  GetType() const;

  virtual double GetStartTime(const int e) const;

  bool           DontWriteOutput() const;

  //Manipulator Functions
  void SetRandomSeed     (const unsigned int seed);
  void SetOutputDirectory(const string OutDirString);
  void SetRunNames       (const string RunNames);
  void SetSolutionFiles  (const string SolFiles);

  virtual void Initialize       (const CModel* pModel,const optStruct &Options); //called prior to ALL ensemble runs
  virtual void UpdateModel      (CModel *pModel,optStruct &Options,const int e); //called prior to each ensemble run
  virtual void StartTimeStepOps (CModel* pModel,optStruct &Options,const time_struct &tt,const int e) {} //called at start of every timestep
  virtual void CloseTimeStepOps (CModel* pModel,optStruct &Options,const time_struct &tt,const int e) {} //called at end of each timestep
  virtual void FinishEnsembleRun(CModel *pModel,optStruct &Options,const time_struct &tt,const int e) {} //called after all ensembles run
};

////////////////////////////////////////////////////////////////////
/// \brief Data abstraction for Monte Carlo model ensemble run
//
class CMonteCarloEnsemble : public CEnsemble
{
private:

  int          _nParamDists; ///< number of parameter distributions for sampling
  param_dist **_pParamDists; ///< array of pointers to parameter distributions


public:
  CMonteCarloEnsemble(const int num_members,const optStruct &Options);
  ~CMonteCarloEnsemble();

  void AddParamDist(const param_dist *dist);

  void Initialize(const CModel* pModel,const optStruct &Options);
  void UpdateModel(CModel *pModel,optStruct &Options, const int e);
};

////////////////////////////////////////////////////////////////////
/// \brief Data abstraction for DDS model ensemble run
//
class CDDSEnsemble : public CEnsemble
{
private:
  double       _r_val;       ///< perturbation value
  double      *_BestParams;  ///< vector of best parameter values
  double      *_TestParams;  ///< vector of test parameter values
  double       _Fbest;       ///< best obj function val

  int          _nParamDists; ///< number of parameter distributions for sampling
  param_dist **_pParamDists; ///< array of pointers to parameter distributions

  long long    _calib_SBID;  ///< observation hydrograph subbasin ID
  diag_type    _calib_Obj;   ///< diagnostic used as objective function (e.g., DIAG_NASH_SUTCLIFFE)
  string       _calib_Period;///< name of calibration period (e.g., CALIB)

  ofstream     _DDSOUT;      ///< output file stream

  double PerturbParam(const double &x_best,
                      const double &upperbound,
                      const double &lowerbound);

public:
  CDDSEnsemble(const int num_members,const optStruct &Options);
  ~CDDSEnsemble();

  void SetPerturbationValue(const double &perturb);
  void SetCalibrationTarget(const long long SBID, const diag_type object_diag, const string period);
  void AddParamDist(const param_dist *dist);

  void Initialize(const CModel* pModel,const optStruct &Options);
  void UpdateModel(CModel *pModel,optStruct &Options,const int e);
  void FinishEnsembleRun(CModel *pModel,optStruct &Options,const time_struct &tt,const int e);
};
#endif
