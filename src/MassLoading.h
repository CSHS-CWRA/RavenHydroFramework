/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ----------------------------------------------------------------
  class definitions:
  CmvMassLoading
  ----------------------------------------------------------------*/

#ifndef MASSLOADING_H
#define MASSLOADING_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"
#include "Model.h"


////////////////////////////////////////////////////////////////////
/// \brief Calculates the mass loading of a constituent (Neumann sources)
//
class CmvMassLoading: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  const CTransportModel* pTransModel;

  int     _constit_ind;         ///< index of constituent being loaded

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvMassLoading(string          constit_name,
                 CTransportModel *pTransportModel,
                 CModelABC       *pModel);
  ~CmvMassLoading();

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

  void GetParticipatingParamList(string     *aP,
                                 class_type *aPC,
                                 int        &nP) const;
};
#endif
