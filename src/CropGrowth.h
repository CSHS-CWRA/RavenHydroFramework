/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------
  class definitions:
  CmvCropHeatUnitEvolve
  ----------------------------------------------------------------*/

#ifndef CROPGROWTH_H
#define CROPGROWTH_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"

///////////////////////////////////////////////////////////////////
/// \brief Algorithms for modelling the evolution of cumulative crop heat units
//
enum CHUevolve_type{
  CHU_ONTARIO ///< Ontario CHU evolution method
};

////////////////////////////////////////////////////////////////////
/// \brief Calculates the evolution of Crop Heat Units
//
class CmvCropHeatUnitEvolve: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  CHUevolve_type type; ///< Method of modeling CHUs selected

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvCropHeatUnitEvolve(const CHUevolve_type snalbtype);
  ~CmvCropHeatUnitEvolve();

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
  static void GetParticipatingStateVarList(CHUevolve_type  snalbtype,
                                           sv_type *aSV,
                                           int     *aLev,
                                           int     &nSV);
};
#endif
