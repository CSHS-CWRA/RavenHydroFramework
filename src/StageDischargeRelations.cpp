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
double CStageDischargeTable::GetDischarge(const double &h, const double &hstart, const double &Qstart, const double &rivdepth, const double &drefelev) const
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
double CBasicWeir::GetDischarge(const double &h, const double &hstart, const double &Qstart, const double &rivdepth, const double &drefelev) const
{
 return 2.0/3.0*_weir_coeff*sqrt(2*GRAVITY)*_crestwidth*pow(max(h-_crest_elev,0.0),1.5);
}

// SLUICE GATE CLASS

/////////////////////////////////////////////////////////////////
/// \brief Sluice gate constructor / destructor
//
CSluiceGate::CSluiceGate(const string name, const double elev, const double width, const double height, const double coeff, int numgates) :
             CStageDischargeRelationABC(name)
{
  _bottom_elev=elev;
  _gatewidth=width;
  _gateheight=height;
  _gate_coeff=coeff;
  _num_gates=numgates;
}
CSluiceGate::~CSluiceGate(){}

/////////////////////////////////////////////////////////////////
/// \param h, stage in absolute height units [m]
/// \returns discharge, in [m3/s]
//
double CSluiceGate::GetDischarge(const double &h, const double &hstart, const double &Qstart, const double &rivdepth, const double &drefelev) const
{
  // free flowing condition (not restrained by backwater)
  // Q = numgates * CWB*sqrt(2gH)
 return _num_gates*_gate_coeff*_gatewidth*_gateheight*sqrt(2*GRAVITY*max(h-_bottom_elev,0.0));
}

// ORIFICE CLASS

/////////////////////////////////////////////////////////////////
/// \brief Orifice constructor / destructor
//
COrifice::COrifice(const string name, const double elev, const double diameter, const double coeff, int numopenings) :
             CStageDischargeRelationABC(name)
{
  _bottom_elev=elev;
  _diameter=diameter;
  _coeff=coeff;
  _num_openings=numopenings;
}
COrifice::~COrifice(){}

/////////////////////////////////////////////////////////////////
/// \param h, stage in absolute height units [m]
/// \returns discharge, in [m3/s]
//
double COrifice::GetDischarge(const double &h, const double &hstart, const double &Qstart, const double &rivdepth, const double &drefelev) const
{
 // free flowing condition (not restrained by backwater)
  // Q = submerged_correction * numopenings * CA*sqrt(2gH)
  // A = (1/4)*PI*d^2

 return _num_openings*_coeff*(PI*pow(_diameter,2)/4)*sqrt(2*GRAVITY*max(h-_bottom_elev,0.0));
}

// BASIC PUMP CLASS

/////////////////////////////////////////////////////////////////
/// \brief Basic pump constructor / destructor
//
CBasicPump::CBasicPump(const string name, const double flow, const double on_elev, const double off_elev) :
             CStageDischargeRelationABC(name)
{
  _flow=flow;
  _on_elev=on_elev;
  _off_elev=off_elev;
}
CBasicPump::~CBasicPump(){}

/////////////////////////////////////////////////////////////////
/// \param h [in] stage in absolute height units [m]
/// \param Q_start [in] outflow at start of timestep [m3/s]
/// \returns discharge, in [m3/s]
//
double CBasicPump::GetDischarge(const double &h, const double &hstart, const double &Qstart, const double &rivdepth, const double &drefelev) const
{
  /// xxx update to use h and hstart to detect decreasing flow?

  if (h > _on_elev) { // stage > on_elev, turn on pump
    return _flow;
  }
  else if (h <= _off_elev) { // stage < off_elev, turn off pump
    return 0;
  } else if (h <= _on_elev && Qstart != 0) { // stage < on_elev but stage > off_elev and pump already on (Qstart !=0), keep pump running
    return _flow;
  } else { // catch all but this case should not occur
    // warning??
    return 0;
  }
}

/// STRUCTURES TO ADD
/*
- generalized seepage (same impementation but generalized for multiple possible instances etc)
- radial gates - https://www.hec.usace.army.mil/confluence/hmsdocs/hmstrm/modeling-reservoirs/spillways
- Ogee spillway - coefficient corrected in simulation based on design head
- pumps - min run time etc
- culvert outlets - https://www.hec.usace.army.mil/confluence/hmsdocs/hmstrm/modeling-reservoirs/outlets
-- now that downstream depth is included, the backwater calculation will be possible (although still rough since Raven is not a hydrualic model)
- non-level dam top (specify dam top x-z points)

*/
