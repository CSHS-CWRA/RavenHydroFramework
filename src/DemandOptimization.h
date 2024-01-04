/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2024 the Raven Development Team
  ----------------------------------------------------------------
  DemandOptimization.h
  ----------------------------------------------------------------*/
#ifndef DEMAND_OPTIMIZATION_H
#define DEMAND_OPTIMIZATION_H

#include "RavenInclude.h"
#include <stdio.h>
#include "Model.h"

#ifdef _LPSOLVE_
namespace lp_lib  {
#include "../lib/lp_solve/lp_lib.h"
}
#endif 

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

///////////////////////////////////////////////////////////////////
/// \brief different expression types  
//
enum exptype 
{
  EXP,     //< default expression
  EXP_OP,  //< operator
  EXP_INV  //< inverse operator 
};
///////////////////////////////////////////////////////////////////
/// \brief different expression term types  
//
enum termtype
{
  TERM_DV,        //< decision variable
  TERM_TS,        //< time series @ts(name,n)  
  TERM_LT,        //< lookup table @lookup(x)
  TERM_CONST,     //< constant 
  TERM_HISTORY,   //< bracketed - !Q123[-2] 
  TERM_MAX,       //< @max(x,y) 
  TERM_MIN,       //< @min(x,y)
  TERM_CONVERT,   //< @convert(x,units)
  TERM_UNKNOWN    //< unknown 
};
///////////////////////////////////////////////////////////////////
/// \brief decision variable types 
//
enum dv_type 
{
  DV_QOUT,    //< outflow from reach
  DV_QOUTRES, //< outflow from reservoir 
  DV_STAGE,   //< reservoir stage
  DV_DELIVERY,//< delivery of water demand 
  DV_SLACK,   //< slack variable for goal satisfaction
  DV_USER     //< user specified decision variable 
};
// data structures used by CDemandOptimizer class
// -------------------------------------------------------------------

//////////////////////////////////////////////////////////////////
/// expression term 
///    individual term in expression
// 
struct expressionTerm 
{
  termtype      type;            //< type of expression 
  double        mult;            //< multiplier of expression (+/- 1, depending upon operator and location in exp)
  bool          reciprocal;      //< true if term is in denominator 
  
  double        value;           //< constant value or conversion multiplier
  CTimeSeries  *pTS;             //< pointer to time series, if this is a named time series
  CLookupTable *pLT;             //< pointer to lookup table, if this is a named time series
  bool          is_nested;       //< true if nested within (i.e., an argument to) another term 
  int           timeshift;       //< for time series (+ or -) or lookback value (+) 
  int           DV_ind;          //< index of decision variable 
  int           nested_ind1;     //< index k of first argument (e.g., for lookup table with term entry)
  int           nested_ind2;     //< index k of second argument (e.g., for min/max functions)
  string        nested_exp1;     //< contents of first argument to function - can be expression
  string        nested_exp2;     //< contents of second argument to function
 
  string        origexp;         //< original string expression 

  expressionTerm(); 
};

//////////////////////////////////////////////////////////////////
/// expression structure
///   abstraction of (A*B*C)+(D*E)+(-F)+(G*H) <= 0  
///   parenthetical collections are groups of terms -example has 4 groups with [3,2,1,2] Terms per group
//
struct expressionStruct //full expression 
{ 
  expressionTerm  ***pTerms;      //< 2D irregular array of pointers to expression terms size:[nGroups][nExpPerGrp[j]]
  int                nGroups;     //< total number of terms groups in expression 
  int               *nTermsPerGrp;//< number of terms per group [size: nGroups]
  comparison         compare;     //< comparison operator (==, <, >)

  expressionStruct();
  ~expressionStruct();
};


//////////////////////////////////////////////////////////////////
// decision variable  
//
struct decision_var
{
  string  name;      //< decision variable names: Qxxxx or Dxxxx where xxxx is SBID 
  dv_type dvar_type; //< decision variable type: e.g., DV_QOUT or DV_DELIVERY
  int     p_index;   //< raw subbasin index p (not SBID) of decision variable (or DOESNT_EXIST if not linked to SB) 
  int     dem_index; //< demand index in subbasin p (or DOESNT_EXIST if type not DV_DELIVERY) 
  double  value;     //< solution value for decision variable
  double  min;       //< minimum bound (default=0) 
  double  max;       //< maximum bound (default unbounded)

  decision_var(string nam, int p, dv_type typ){name=nam; p_index=p; value=0.0;min=0.0;max=ALMOST_INF;dvar_type=typ;dem_index=DOESNT_EXIST;}
};
//////////////////////////////////////////////////////////////////
// decision variable constraint condition 
//
struct dv_condition
{
  string      dv_name; //< decision variable name (e.g., Q1023) or "MONTH"
  double      value;   //< conditional value 
  double      value2;  //< second conditional (if COMPARE_BETWEEN)
  comparison  compare; //> comparison operator, e.g., COMPARE_IS_EQUAL

  dv_condition(){
    dv_name="";
    value=0.0;
    value2=0.0;
    compare=COMPARE_IS_EQUAL;
  }
};
//////////////////////////////////////////////////////////////////
// decision variable constraint or goal (soft constraint)
//
struct dv_constraint
{
  string            name;          //< constraint name 
  expressionStruct *pExpression;   //< constraint expression 
  
  bool              is_goal;       //< true if constraint is soft (goal rather than constraint)
  int               priority;      //< priority (default==1, for goals only)
  double            penalty_under; //< penalty if under specified value (for goals only)
  double            penalty_over;  //< penalty if over value (for goals only)
  int               slack_ind1;    //< slack index (0.._nSlackVars) of under/over slack index for goal, or DOESNT_EXIST if constraint
  int               slack_ind2;    //< slack index (0.._nSlackVars) of over slack index for target goal, or DOESNT_EXIST if constraint

  double            penalty_value; //< (from solution) penalty incurred by not satisfying goal (or zero for constraint) 

  int               nConditions;   //< number of conditional statments 
  dv_condition    **pConditions;   //< array of pointers to conditional statements 
  bool              conditions_satisfied; 

  dv_constraint();
  ~dv_constraint();
  void AddCondition(dv_condition *pCondition);
};
///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for demand optimization 
//
class CDemandOptimizer
{
private: /*------------------------------------------------------*/

  CModel          *_pModel;            //< pointer to model

  int              _nDecisionVars;      //< total number of decision variables considered
  decision_var   **_pDecisionVars;      //< array of pointers to decision variable structures [size:_nDecisionVars]

  int              _nConstraints;       //< number of user-defined enforced constraints/goals in management model 
  dv_constraint  **_pConstraints;       //< array of pointers to user-defined enforced constraints/goals in management model

  int              _nEnabledSubBasins;  //< number of enabled subbasins in model 
  
  int              _nReservoirs;        //< local storage of number of simulated lakes/reservoirs 
  int             *_aResIndices;        //< storage of enabled reservoir indices (0:_nReservoirs or DOESNT_EXIST) [size:_nSubBasins] 

  int              _nDemands;           //< local storage of number of demand locations (:IrrigationDemand/:WaterDemand + :ReservoirExtraction time series)
  int             *_aDemandIDs;         //< local storage of demand IDs [size:_nDemands]
  string          *_aDemandAliases;     //< list of demand aliases [size: _nDemands] 
  double          *_aDemandPenalties;   //< array of priority weights [size: _nDemands]
  double          *_aDelivery;          //< array of delivered water for each demand [m3/s] [size: _nDemands]
  double          *_aCumDelivery;       //< array of cumulative deliveries of demand since _aCumDivDate of this year [m3] [size: _nDemands] 
  int             *_aCumDelDate;        //< julian date to calculate cumulative deliveries from {default: Jan 1)[size: _nDemands]
  
  //int            _nDemandGroups;      //< number of demand groups 
  //CDemandGroup **_pDemandGroups;      //< array of pointers to demand groups   
  
  double          *_aSlackValues;       //< array of slack variable values [size: _nSlackVars]
  int              _nSlackVars;         //< number of slack variables 

  int             _nUserDecisionVars;   //< number of user-specified decision variables 
  int             _nUserConstants;      //< number of user-specified named constants
  string         *_aUserConstNames;     //< array of names of user-specified constants
  double         *_aUserConstants;      //< array of values of user-specified constants
  int             _nUserTimeSeries;     //< number of user variable time series
  CTimeSeries   **_pUserTimeSeries;     //< array of pointers to user variable time series
  int             _nUserLookupTables;   //< number of user variable lookup tables 
  CLookupTable  **_pUserLookupTables;   //< array of pointers to user variable lookup tables 

  int             _nHistoryItems;      //< size of flow/demand/stage history that needs to be stored (in time steps) (from :LookbackDuration)
  double        **_aQhist;             //< history of subbasin discharge [size: nSBs * _nHistoryItems]
  double        **_aDhist;             //< history of actual diversions [size: nSBs * _nHistoryItems]
  double        **_ahhist;             //< history of actual reservoir stages  (or 0 for non-reservoir basins) [size: nSBs * _nHistoryItems]

  ofstream        _DEMANDOPT;          //< ofstream for DemandOptimization.csv
  ofstream        _GOALSAT;            //< ofstream for GoalSatisfaction.csv

  bool     ConvertToExpressionTerm(const string s, expressionTerm* term, const int lineno, const string filename)  const;
  int               GetDVColumnInd(const dv_type typ, const int counter) const;

#ifdef _LPSOLVE_
  void           AddConstraintToLP(const int i, lp_lib::lprec *pLinProg, const time_struct &tt,int *col_ind, double *row_val) const;
#endif 
  double              EvaluateTerm(expressionTerm **pTerms,const int k, const int nn) const;

  void             CheckConditions(const int ii, const time_struct &tt); 

public: /*------------------------------------------------------*/
  CDemandOptimizer(CModel *pMod);
  ~CDemandOptimizer();

  int    GetDemandIndexFromName(const string dname) const;
  double GetNamedConstant      (const string s) const; 
  double GetTotalWaterDemandDelivery(const int p) const;

  void   SetHistoryLength      (const int n);
  void   SetCumulativeDate     (const int julian_date, const string demandID);
  
  void   AddDecisionVar        (const decision_var *pDV);
  void   SetDecisionVarBounds  (const string name, const double &min, const double &max);
  dv_constraint *AddConstraint (const string name, expressionStruct *exp, const bool soft_constraint);
  void   AddUserConstant       (const string name, const double &val);
  void   AddUserTimeSeries     (const CTimeSeries *pTS);
  void   AddUserLookupTable    (const CLookupTable *pLUT);
   
  //void MultiplyDemand        (const string dname, const double &mult);
  //void MultiplyGroupDemand   (const string groupname, const double &mult);
  //void SetDemandPenalty      (const string dname, const double &pen);
  //void SetDemandPriority     (const string dname, const int &prior);

  expressionStruct *ParseExpression(const char **s, const int Len, const int lineno, const string filename) const;

  void   Initialize            (CModel *pModel, const optStruct &Options);
  void   InitializePostRVMRead (CModel *pModel, const optStruct &Options); 
  void   SolveDemandProblem    (CModel *pModel, const optStruct &Options, const double *aSBrunoff, const time_struct &tt);

  void   WriteOutputFileHeaders(const optStruct &Options);
  void   WriteMinorOutput      (const optStruct &Options,const time_struct &tt);
  void   CloseOutputStreams    ();
};
#endif