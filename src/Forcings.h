#include "RavenInclude.h"

#ifndef FORCINGS_H
#define FORCINGS_H



forcing_type GetForcingTypeFromString(const string &f_string);
string       ForcingToString         (const forcing_type ftype);
double GetForcingFromString          (const string &forcing_string, const force_struct &f);
string GetForcingTypeUnits           (      forcing_type ftype);
void   ZeroOutForcings    (force_struct &F);//in CommonFunctions.cpp
#endif
