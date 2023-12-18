/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2023 the Raven Development Team
  ----------------------------------------------------------------
  class definitions:
  CmvFrozenLake
  ----------------------------------------------------------------*/

#ifndef FROZEN_LAKE_H
#define FROZEN_LAKE_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"
#include "Transport.h"

///////////////////////////////////////////////////////////////////
/// \brief Method of calculating freezing of the lake surface
//
enum lakefreeze_type{
  LFREEZE_BASIC,      ///< basic - based upon potential melt energy, moderated by SWE
  LFREEZE_THERMAL      ///< based upon enthaly balance model
};

///////////////////////////////////////////////////////////////////
/// \brief Calculates the depth of frozen soil or lake ice
//
class CmvFrozenLake: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  lakefreeze_type _type; ///< algorithm choice

  const CTransportModel *_pTransModel; ///< transport model pointer (if needed for enthalpy model)

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvFrozenLake(const lakefreeze_type type,
                const CTransportModel *_pTransModel,
                CModelABC             *pModel);
  ~CmvFrozenLake();

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
                        const time_struct &t,
                              double      *rates) const;

  void        GetParticipatingParamList   (string  *aP, class_type *aPC, int &nP) const;
  static void GetParticipatingStateVarList(lakefreeze_type  dtype,
                                           sv_type *aSV,
                                           int     *aLev,
                                           int     &nSV);
};

#endif
