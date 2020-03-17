/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2020 the Raven Development Team
  ----------------------------------------------------------------
  class definitions:
  CmvAdvection
  ----------------------------------------------------------------*/

#ifndef ADVECTION_H
#define ADVECTION_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"
#include "Model.h"


////////////////////////////////////////////////////////////////////
/// \brief Calculates the advection of a substance
//
class CmvAdvection: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  const CTransportModel* pTransModel;

  int     _constit_ind;         ///< index of constituent being advected

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvAdvection(string constit_name, CTransportModel *pTransportModel);
  ~CmvAdvection();

  //inherited functions
  void Initialize();
  void GetRatesOfChange(const double              *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double      *rates) const;
  void ApplyConstraints(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &t,
                        double      *rates) const;

  void        GetParticipatingParamList   (string  *aP, class_type *aPC, int &nP) const;

  /*static void GetParticipatingStateVarList(sv_type *aSV,
    int     *aLev,
    int     &nSV);*/ //not used
};
#endif
