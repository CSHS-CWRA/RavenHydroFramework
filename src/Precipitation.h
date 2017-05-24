/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------
  class definitions:
  CmvPrecipitation
  ----------------------------------------------------------------*/

#ifndef PRECIPITATION_H
#define PRECIPITATION_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"
#include "Infiltration.h"

///////////////////////////////////////////////////////////////////
/// \brief Calculates distribution of snow and rain precipitation to surface storage
/// \remark all state variables are potential receptacles for precipitation
//
class CmvPrecipitation: public CHydroProcessABC
{
private:/*------------------------------------------------------*/

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvPrecipitation();
  ~CmvPrecipitation();

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
                        const time_struct &tt,
                        double      *rates) const;

  static void GetParticipatingStateVarList(sv_type *aSV, int *aLev, int &nSV);
  void        GetParticipatingParamList   (string *aP, class_type *aPC, int &nP) const;

};
#endif
