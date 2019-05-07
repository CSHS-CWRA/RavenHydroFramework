/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2019 the Raven Development Team
----------------------------------------------------------------*/
#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include "RavenInclude.h"
#include "Model.h"
#include "SoilAndLandClasses.h"

enum disttype {
  DIST_UNIFORM,     ///< Uniform distribution
  DIST_NORMAL       ///< Gaussian (normal) distribution
};
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
  //transformation trans; e.g., log transform
        
};
double SampleFromDistribution(param_dist *dist);

////////////////////////////////////////////////////////////////////
/// \brief Data abstraction for model ensemble run
/// \details Stores ensemble member details and functionality on how to modify each ensemble member
//
class CEnsemble
{
protected:/*------------------------------------------------------*/
  ensemble_type _type;        ///< ensemble type (e.g., ENSEMBLE_MONTECARLO)
  
  int           _nMembers;    ///< number of ensemble members

  string       *_aOutputDirs; ///< array of output directory names [size: _nMembers] 
  string       *_aRunNames;   ///< array of output runnames [size: _nMembers]

  int           _rand_seed;   ///< random seed 

public:/*-------------------------------------------------------*/
  CEnsemble(const int num_members, const optStruct &Options);
  ~CEnsemble();

  //Acessor Function
  int           GetNumMembers();
  ensemble_type GetType();

  //Manipulator Functions
  void SetRandomSeed     (const int seed);
  void SetOutputDirectory(const string OutDirString);
  void SetRunNames       (const string RunNames);

  virtual void Initialize(const optStruct &Options);
  virtual void UpdateModel(CModel *pModel,optStruct &Options,const int e);
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

  void Initialize(const optStruct &Options);
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

  double PerturbParam(const double &x_best,
                      const double &upperbound,
                      const double &lowerbound);

public:
  CDDSEnsemble(const int num_members,const optStruct &Options);
  ~CDDSEnsemble();

  void SetPerturbationValue(const double &perturb);
  void AddParamDist(const param_dist *dist);

  void Initialize(const optStruct &Options);
  void UpdateModel(CModel *pModel,optStruct &Options,const int e);
  void FinishEnsembleRun(CModel *pModel,optStruct &Options,const int e);
};
#endif
