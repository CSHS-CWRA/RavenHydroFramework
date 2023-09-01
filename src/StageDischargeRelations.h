/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ----------------------------------------------------------------
  StageDischargeRelations.h
  ------------------------------------------------------------------
  defines abstract stage discharge relations
  CStageDichargeRelationABC
  CStageDischargeTable
  CBasicWeir
  CCircularCulvert
  ...
  ----------------------------------------------------------------*/
#ifndef STAGE_DISCHARGE_H

/*****************************************************************
    Class CStageDischargeRelation
------------------------------------------------------------------
    Data Abstraction for abstract stage discharge relation
  ******************************************************************/
class CStageDischargeRelationABC
{
protected:
  string _name;

  CStageDischargeRelationABC (const string name){_name=name;}
  ~CStageDischargeRelationABC(){}

public:
  string         GetName     ()                      {return _name;}
  // virtual double GetDischarge(const double &h) const {return 0.0;}
  virtual double GetDischarge(const double &h, const double &hstart, const double &Qstart, const double &rivdepth, const double &drefelev) const {return 0.0;};
};

/*****************************************************************
    Class CStageDischargeTable
------------------------------------------------------------------
    Data Abstraction for tabular stage discharge relation
  ******************************************************************/
class CStageDischargeTable : public CStageDischargeRelationABC
{
public:
  int     _Np;                  ///< number of points on rating curve
  double* _aStage;              ///< Base rating curve for stage elevation [m]
  double* _aQ;                  ///< Rating curve for overflow (e.g., weir) flow rates [m3/s]

public:
  CStageDischargeTable(const string name, const double *h, const double *Q, const int N);
  ~CStageDischargeTable();

  double GetDischarge(const double &h, const double &hstart, const double &Qstart, const double &rivdepth, const double &drefelev) const;
};

/*****************************************************************
    Class CBasicWeir
------------------------------------------------------------------
    Data Abstraction for basic rectangular weir
  ******************************************************************/
class CBasicWeir : public CStageDischargeRelationABC
{
public:
  double _crest_elev;         ///< elevation of crest [m]
  double _crestwidth;         ///< width of crest [m]
  double _weir_coeff;         ///< weir coefficient [-]

public:
  CBasicWeir(const string name, const double elev, const double width, const double coeff);
  ~CBasicWeir();

  double GetDischarge(const double &h, const double &hstart, const double &Qstart, const double &rivdepth, const double &drefelev) const;
};

/*****************************************************************
    Class CSluiceGate
------------------------------------------------------------------
    Data Abstraction for sluice gate
  ******************************************************************/
class CSluiceGate : public CStageDischargeRelationABC
{
public:
  double _bottom_elev;        ///< bottom elevation of crest [m]
  double _gatewidth;          ///< width of gate [m]
  double _gateheight;         ///< height of gate when raised [m]
  double _gate_coeff;         ///< gate coefficient [-]
  int    _num_gates;          ///< number of parallel (identical) gates [-]

public:
  CSluiceGate(const string name, const double elev, const double width, const double height, const double coeff, int numgates);
  ~CSluiceGate();

  double GetDischarge(const double &h, const double &hstart, const double &Qstart, const double &rivdepth, const double &drefelev) const;
};

/*****************************************************************
    Class COrifice
------------------------------------------------------------------
    Data Abstraction for orifice (submerged opening, assumed circular)
  ******************************************************************/
class COrifice : public CStageDischargeRelationABC
{
public:
  double _bottom_elev;        ///< bottom elevation of orifice [m]
  double _diameter;           ///< diameter of circular orifice [m]
  double _coeff;              ///< discharge coefficient [-]
  int    _num_openings;       ///< number of identical openings [-]

public:
  COrifice(const string name, const double elev, const double diameter, const double coeff, int numopenings);
  ~COrifice();

  double GetDischarge(const double &h, const double &hstart, const double &Qstart, const double &rivdepth, const double &drefelev) const;
};

/*****************************************************************
    Class CBasicPump
------------------------------------------------------------------
    Data Abstraction for basic pump
  ******************************************************************/
class CBasicPump : public CStageDischargeRelationABC
{
public:
  double _flow;           ///< fixed flow for pump [cms]
  double _on_elev;        ///< diameter of circular orifice [m]
  double _off_elev;       ///< discharge coefficient [-]

public:
  CBasicPump(const string name, const double flow, const double on_elev, const double off_elev);
  ~CBasicPump();

  double GetDischarge(const double &h, const double &hstart, const double &Qstart, const double &rivdepth, const double &drefelev) const;
};
#endif
