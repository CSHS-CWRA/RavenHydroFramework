/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2024 the Raven Development Team
----------------------------------------------------------------*/
#include "LookupTable.h"

//////////////////////////////////////////////////////////////////
// Constructor/Destructor
//
CLookupTable::CLookupTable(string name, double* x, double* y, int N)
{
  _name=name;
  _aX=new double [N];
  _aY=new double [N];
  _nItems=N;
  for (int i = 0; i < _nItems; i++) {
    _aX[i]=x[i];
    _aY[i]=y[i];
  }
}
CLookupTable::~CLookupTable() {
  delete [] _aX;
  delete [] _aY;
}
//////////////////////////////////////////////////////////////////
// Accessors
//
string CLookupTable::GetName() const
{
  return _name;
}
double CLookupTable::GetValue(const double& x) const
{
  if (x==RAV_BLANK_DATA){return RAV_BLANK_DATA;}
  return InterpolateCurve(x,_aX,_aY,_nItems,false);
}
double CLookupTable::GetSlope(const double& x) const
{
  if (x==RAV_BLANK_DATA){return RAV_BLANK_DATA;}
  double dx=0.001*(_aX[_nItems-1]-_aX[0]);
  return (1.0/dx)*(InterpolateCurve(x+dx,_aX,_aY,_nItems,false)-InterpolateCurve(x,_aX,_aY,_nItems,false));
}
