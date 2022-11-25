/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2022 the Raven Development Team
  ----------------------------------------------------------------*/
#include "RavenInclude.h"

#ifndef FORCINGS_H
#define FORCINGS_H

//functions defined in Forcings.cpp

forcing_type GetForcingTypeFromString(const string &f_string);
string                ForcingToString(const forcing_type ftype);
double           GetForcingFromString(const string &forcing_string, const force_struct &f);
double             GetForcingFromType(const forcing_type &ftype, const force_struct &f);
string            GetForcingTypeUnits(      forcing_type ftype);
void                  ZeroOutForcings(force_struct &F);
void            CopyDailyForcingItems(force_struct& Ffrom, force_struct& Fto); 
#endif
