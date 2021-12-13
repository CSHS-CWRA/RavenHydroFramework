/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ----------------------------------------------------------------*/

#include "RavenInclude.h"
#include "StageDischargeRelations.h"

// STAGE DISCHARGE TABLE CLASS
/////////////////////////////////////////////////////////////////
/// \brief Stage discharge table constructor / destructor
//
CStageDischargeTable::CStageDischargeTable(const string name, const double *h, const double *Q, const int N) 
  :CStageDischargeRelationABC(name) 
{
  _Np=N;
  ExitGracefullyIf(N<=1,"CStageDischargeTable constructor : table must have at least 2 entries",BAD_DATA);
  _aStage=NULL;
  _aQ=NULL;
  _aStage=new double [_Np];
  _aQ    =new double [_Np];
  ExitGracefullyIf(_aQ==NULL,"CStageDischargeTable constructor",OUT_OF_MEMORY);
  for (int i = 0; i < _Np; i++) {
    _aStage[i]=h[i];
    _aQ    [i]=Q[i];
  }
  for (int i = 1; i < _Np; i++) {
    if ((_aStage[i] - _aStage[i-1]) <= 0.0) {
      ExitGracefully("CStageDischargeTable constructor: (stage,discharge) pairs must be ordered by increasing stage",BAD_DATA);
    }
    if ((_aQ[i] - _aQ[i-1]) < 0.0) {
      ExitGracefully("CStageDischargeTable constructor: discharge must monotonically increase with stage",BAD_DATA);
    }
  }
}
CStageDischargeTable::~CStageDischargeTable() 
{
  delete [] _aStage;
  delete [] _aQ;
 }
/////////////////////////////////////////////////////////////////
/// \param h, stage in absolute or relative height units [m]
/// \returns discharge, in [m3/s]
//
double CStageDischargeTable::GetDischarge(const double& h) const 
{
  return InterpolateCurve(h,_aStage,_aQ,_Np,true);
}


// BASIC WEIR CLASS

/////////////////////////////////////////////////////////////////
/// \brief Basic rectangular weir constructor / destructor
//
CBasicWeir::CBasicWeir(const string name, const double elev, const double width, const double coeff) : 
             CStageDischargeRelationABC(name) 
{
  _crest_elev=elev;
  _crestwidth=width;
  _weir_coeff=coeff;
}
CBasicWeir::~CBasicWeir(){}

/////////////////////////////////////////////////////////////////
/// \param h, stage in absolute height units [m]
/// \returns discharge, in [m3/s]
//
double CBasicWeir::GetDischarge(const double& h) const 
{
 return 2.0/3.0*_weir_coeff*sqrt(2*GRAVITY)*_crestwidth*pow(max(h-_crest_elev,0.0),1.5);
}
