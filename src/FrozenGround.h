/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2026 the Raven Development Team
  ----------------------------------------------------------------
  class definitions:
  CmvFrozenGround
  ----------------------------------------------------------------*/

#ifndef FROZEN_GROUND_H
#define FROZEN_GROUND_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"
#include "ModelABC.h"

///////////////////////////////////////////////////////////////////
/// \brief Method of calculating freezing of the ground surface
//
enum groundfreeze_type{
  FREEZE_STEFAN,  ///< stefan model 
  FREEZE_THERMAL, ///< based upon enthaly balance model
  FREEZE_RANKINEN ///< from Rankinen et al., HESS, 2004
};

///////////////////////////////////////////////////////////////////
/// \brief Calculates the depth of frozen soil or lake ice
//
class CmvFrozenGround: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  groundfreeze_type _type; ///< algorithm choice

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvFrozenGround(const groundfreeze_type dtype, CModelABC *pModel);
  ~CmvFrozenGround();

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
  static void GetParticipatingStateVarList(groundfreeze_type  dtype,
                                           sv_type *aSV,
                                           int     *aLev,
                                           int     &nSV);
};

#endif
