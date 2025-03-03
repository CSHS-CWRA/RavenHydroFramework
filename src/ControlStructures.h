/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ----------------------------------------------------------------
  ControlStructures.h
  ------------------------------------------------------------------
  defines abstract control structures, flow regimes, and regime conditions
  ----------------------------------------------------------------*/
#ifndef CONTROL_STRUCTURES_H
#define CONTROL_STRUCTURES_H

#include "RavenInclude.h"
#include "StageDischargeRelations.h"

/*****************************************************************
 RegimeCondition: defines single condition under which regime applies
*****************************************************************/
struct RegimeCondition
{
  string     variable;     ///< one of {DATE, DOY, STAGE, FLOW, ...}
  long long  basinID;      ///< defaults to -1 (this basin)
  comparison compare;      ///< e.g., greater than, between, etc.
  double     compare_val1; ///< primary comparison value
  double     compare_val2; ///< secondary comparison value (used for between)
  RegimeCondition() { variable = "DATE";  basinID=DOESNT_EXIST; compare=COMPARE_BETWEEN;compare_val1=compare_val2=0.0;}
};
//example :Condition DATE IS_BETWEEN 2010-07-01 2012-06-31 (converted to model time in compare_vals)
//example :Condition DAY_OF_YEAR IS_BETWEEN 213 320
//example :Condition STAGE IS_GREATER_THAN 340

/*****************************************************************
 RegimeConstraint: defines single constraint to be applied to outflow
*****************************************************************/
struct RegimeConstraint
{
  string     variable; // one of {FLOW, FLOW_DELTA,...}
  comparison compare_typ;
  double     compare_val1;
  double     compare_val2;
};
//example :Constraint FLOW LESS_THAN 4000
//example :Constraint FLOW_DELTA BETWEEN 30 40 [m3/s/d]

class CModel;
/*****************************************************************
    Class COutflowRegime
------------------------------------------------------------------
    Data Abstraction for flow regime that defines conditions for operating in this regime
  ******************************************************************/
class COutflowRegime
{
private:
  string             _name;

  const CModel      *_pModel;              ///< pointer to access model for querying flows

  CStageDischargeRelationABC *_pCurve;      ///< pointer to stage discharge curve (or other Q=f(h) relation)

  RegimeCondition  **_pConditions;         ///< conditions for operation *e.g., date range, etc., ALL conditions must be met
  int                _nConditions;         ///< number of conditions for operation under this regime

  RegimeConstraint **_pConstraints;        ///< constraints on operation, sorted by priority
  int                _nConstraints;        ///< number of constraints

  bool IsConditionMet(int i) const;

public:/*-------------------------------------------------------*/
  //Constructors:
  COutflowRegime(const CModel *pMod, const string name);
  ~COutflowRegime();

  string             GetName() const;
  void    CheckConditionsAndConstraints() const;
  const CStageDischargeRelationABC *GetCurve() const;

  void              SetCurve(CStageDischargeRelationABC *pCurv);
  void    AddRegimeCondition(RegimeCondition  *pCond);
  void   AddRegimeConstraint(RegimeConstraint *pCond);

  bool      AreConditionsMet(const time_struct &tt) const;
  double          GetOutflow(const double &h, const double &hstart, const double &Qstart, const long long &target_SBID, const double &drefelev) const;
};


/*****************************************************************
    Class CControlStructure
------------------------------------------------------------------
    Data Abstraction for reservoir outflow control structure
  ******************************************************************/
class CControlStructure
{
private:/*-------------------------------------------------------*/
  string           _name;                ///< structure name
  long long        _SBID;                ///< subbasin ID
  long long        _target_SBID;         ///< outflow directed to this subbasin ID
  double           _dRefElev;            ///< downstream reference elevation [m] (used in some structures for backwater, limiting flow, etc.)

  int              _nRegimes;            ///< number of flow regimes
  COutflowRegime **_aRegimes;            ///< array of pointers to flow regimes with unique Q=f(h,Q,...)

public:/*-------------------------------------------------------*/
  //Constructors:
  CControlStructure(const string name,const long long SBID, const long long downID);
  ~CControlStructure();

  void      AddRegime       (COutflowRegime *pRegime);
  void      SetTargetBasin  (const long long SBID);

  void      SetDownstreamRefElevation (const double dRefElev);

  string    GetName         () const;
  long long GetTargetBasinID() const;

  double    GetOutflow(const double &stage,const double &stage_start,const double &Q_start,const time_struct &tt) const;

  string    GetCurrentRegimeName(const time_struct &tt) const;
};
#endif
