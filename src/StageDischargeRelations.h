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
  virtual double GetDischarge(const double &h) const {return 0.0;}
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

  double GetDischarge(const double &h) const;
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

  double GetDischarge(const double &h) const;

};
#endif