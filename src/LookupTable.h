/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2024 the Raven Development Team
  ----------------------------------------------------------------
  LookupTable.h
  ----------------------------------------------------------------*/
#ifndef LOOKUPTABLE_H
#define LOOKUPTABLE_H

#include "RavenInclude.h"

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for general lookup table y(x)
//
class CLookupTable
{
 private:
  string  _name;
  double *_aX;
  double *_aY;
  int     _nItems;

 public:
  CLookupTable(string name, double *x, double *y, int NM);
  ~CLookupTable();

  string GetName() const;
  double GetValue(const double &x) const;
  double GetSlope(const double &x) const;
};
#endif
