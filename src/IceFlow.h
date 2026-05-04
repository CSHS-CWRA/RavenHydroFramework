/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2026 the Raven Development Team
  ----------------------------------------------------------------
  class definitions:
  CmvLatIceFlow
  ----------------------------------------------------------------*/

#ifndef ICEFLOW_H
#define ICEFLOW_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"
#include "LateralExchangeABC.h"

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for the redistribution of snow based on slope and snow SWE

class CmvLatIceFlow: public CLateralExchangeProcessABC
{
private:/*------------------------------------------------------*/
  double *_aInitGlacierHeight; //initial glacier height, in mm [size: nHRUs]

  CLatConnect **_pLatConnect;  //< array of pointers to HRU source/recipient pairs [size:_nLatConnect]
  int           _nLatConnect;  //< number of HRU source/recipient pairs

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvLatIceFlow(CModel *pModel);
  ~CmvLatIceFlow();

  void SetHRUWidth(const int k, const double &w);
  void AddConnectionArray(CLatConnect **pConnections, const int nConnections);

  //inherited functions
  void Initialize();
  void InitializePostRVC(const optStruct &Options);

  static void GetParticipatingStateVarList(sv_type *aSV, int *aLev, int &nSV);

  void GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const;

  void GetLateralExchange(const double * const *state_vars, //array of all SVs for all HRUs, [k][i]
                          const CHydroUnit * const *pHRUs,
                          const optStruct   &Options,
                          const time_struct &tt,
                                double      *exchange_rates) const;//purely virtual

};



#endif
