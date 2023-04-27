/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2023 the Raven Development Team
  ----------------------------------------------------------------
  class definitions:
  CProcessGroup
  ----------------------------------------------------------------*/

#ifndef PROCESSGROUP_H
#define PROCESSGROUP_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"
/*enum procgrp_type{
  GRP_SOLVER,
  GRP_WEIGHTEDAVG
};*/
///////////////////////////////////////////////////////////////////
/// \brief groups processes into single process
/// \remark all state variables are potential receptacles for precipitation
//
class CProcessGroup: public CHydroProcessABC
{
private:/*------------------------------------------------------*/

  CHydroProcessABC **_pSubProcesses;  ///< array of subprocesses
  int                _nSubProcesses;  ///< number of subprocesses
  string             _name;

  double            *_aWeights;       ///< array of weights
  //procgrp_type      _type;

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CProcessGroup(string name);
  ~CProcessGroup();

  //inherited functions
  void Initialize();
  void GetRatesOfChange(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double      *rates) const;
  void ApplyConstraints(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double      *rates) const;

  static void GetParticipatingStateVarList(sv_type *aSV, int *aLev, int &nSV);
  void        GetParticipatingParamList   (string *aP, class_type *aPC, int &nP) const;

  //accessor functions
  int GetGroupSize() const;

  //manipulator functions
  void AddProcess(CHydroProcessABC *pProc);
  void SetWeights(const double *aWts, const int nVal);
  void CalcWeights(const double *aVals,const int nVal);
};

#endif