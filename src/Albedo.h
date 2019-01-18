/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------
  class definitions:
  CmvSnowAlbedo
  CmvSoilAlbedo
  ----------------------------------------------------------------*/

#ifndef ALBEDO_H
#define ALBEDO_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"

////////////////////////////////////////////////////////////////////
/// \brief enumerated type for snow albedo evolution algorithms
enum snowalb_type{
  SNOALB_UBCWM,         ///< UBCWM method
  SNOALB_CRHM_ESSERY,   ///< ported from CRHM (Pomeroy et al., 2007) ClassalbedoRichard
  SNOALB_BAKER          ///< Baker et al (1990) as ported from CRHM
};

////////////////////////////////////////////////////////////////////
/// \brief Calculates the evolution of snow albedo
//
class CmvSnowAlbedoEvolve: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  snowalb_type type; ///< Type of albedo algorithm selected

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvSnowAlbedoEvolve(const snowalb_type snalbtype);
  ~CmvSnowAlbedoEvolve();

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
  static void GetParticipatingStateVarList(snowalb_type  snalbtype,
                                           sv_type *aSV,
                                           int     *aLev,
                                           int     &nSV);
};
#endif
