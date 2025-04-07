/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2025 the Raven Development Team
  ----------------------------------------------------------------
  DemandOptimization.h
  ----------------------------------------------------------------*/
#ifndef DEMAND_OPTIMIZATION_H
#define DEMAND_OPTIMIZATION_H

#include "RavenInclude.h"
#include <stdio.h>
#include "Model.h"
#include "LookupTable.h"
#include "Demands.h"
#include "Expression.h"

#ifdef _LPSOLVE_
namespace lp_lib  {
#ifdef _WIN32
#include "../lib/lp_solve/lp_lib.h"
#else
#include "../lib/lp_solve_unix/lp_lib.h"
#endif
}
#endif

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
/// \brief decision variable types
//
enum dv_type
{
  DV_QOUT,    //< outflow from reach
  DV_QOUTRES, //< outflow from reservoir
  DV_STAGE,   //< reservoir stage
  DV_DSTAGE,  //< change in reservoir stage over time step
  DV_BINRES,  //< binary integer value for above/beneath stage-discharge sill
  DV_DELIVERY,//< delivery of water demand
  DV_RETURN,  //< return flows to reach
  DV_USER,    //< user specified decision variable
  DV_SLACK    //< slack variable for goal satisfaction
};

//////////////////////////////////////////////////////////////////
// decision variable
//
struct decision_var
{
  string  name;      //< decision variable names: Qxxxx or Dxxxx where xxxx is SBID
  dv_type dvar_type; //< decision variable type: e.g., DV_QOUT or DV_DELIVERY
  int     p_index;   //< raw subbasin index p (not SBID) of decision variable (or DOESNT_EXIST if not linked to SB) (target basin for DV_RETURN, source basin for DV_DELIVERY)
  int     dem_index; //< demand index ii in subbasin p (or DOESNT_EXIST if type not DV_DELIVERY/DV_RETURN) (ii..nDemands in SB p_index, usually <=1)
  int     loc_index; //< local index (rescount or subbasincount or demand count)
  double  value;     //< solution value for decision variable
  double  min;       //< minimum bound (default=0)
  double  max;       //< maximum bound (default unbounded)

  decision_var(string nam, int p, dv_type typ, int loc_ind)
  {
    name=nam;
    p_index=p;
    loc_index=loc_ind;
    value=0.0;min=-ALMOST_INF;max=ALMOST_INF;dvar_type=typ;dem_index=DOESNT_EXIST;
  }
};
//////////////////////////////////////////////////////////////////
// non-linear variable
//
struct nonlin_var
{
  string name;      //< name of nonlinear variable !Qxxx or !UserVar
  string target;    //< decision variable name Qxxx or UserVar
  double guess_val; //< current guess of value
  int    DV_index;  //index of linear decision variable (0:_nDecisionVars-1) corresponding to nonlinear variable

  nonlin_var(string nam, string targ) {
    name=nam; target=targ; DV_index=DOESNT_EXIST; guess_val=0.0;
  }
};
//////////////////////////////////////////////////////////////////
// goal/constraint condition
//
struct exp_condition
{
  string      dv_name;      //< decision variable name (e.g., Q1023) or "MONTH" or "DATE" or "DAY_OF_YEAR" or "YEAR"
  double      value;        //< conditional value
  double      value2;       //< second conditional (if COMPARE_BETWEEN)
  string      date_string;  //< conditional value (if date)
  string      date_string2; //< second conditional (if DATE COMPARE_BETWEEN)
  comparison  compare;      //> comparison operator, e.g., COMPARE_IS_EQUAL
  int         p_index;      //> subbasin or demand index of LHS of condition expression (or DOESNT_EXIST)

  expressionStruct *pExp;   //< condition expression (or NULL if not used)

  exp_condition(){
    dv_name="";
    value=value2=0.0;
    date_string=date_string2="";
    compare=COMPARE_IS_EQUAL;
    p_index=DOESNT_EXIST;
    pExp=NULL;
  }
};
//////////////////////////////////////////////////////////////////
// operating regime
//
struct op_regime
{
  string            reg_name;             //< regime name

  expressionStruct *pExpression;          //< goal expression

  int               nConditions;          //< number of conditional statments
  exp_condition   **pConditions;          //< array of pointers to conditional statements

  double            penalty_under;        //< penalty associated with not satisfying goal (or RAV_BLANK_DATA, if default is to be used)
  double            penalty_over;         //< penalty associated with not satisfying goal (or RAV_BLANK_DATA, if default is to be used)

  op_regime(string name) {
    reg_name=name;
    pExpression=NULL;
    nConditions=0;
    pConditions=NULL;
    penalty_under=RAV_BLANK_DATA;
    penalty_over =RAV_BLANK_DATA;
  }
};
//////////////////////////////////////////////////////////////////
// management goal (or constraint) or DV definition (a special constraint defining a user specified DV)
//
struct managementGoal
{
  string            name;          //< goal or constraint name

  bool              is_goal;       //< true if constraint is soft (goal rather than constraint)
  bool              is_nonlinear;  //< current expression has non-linear terms
  int               priority;      //< priority (default==1, for goals only)
  double            penalty_under; //< DEFAULT penalty if under specified value (for goals only)
  double            penalty_over;  //< DEFAULT penalty if over value (for goals only)
  int               slack_ind1;    //< slack index (0.._nSlackVars) of under/over slack index for goal, or DOESNT_EXIST if constraint
  int               slack_ind2;    //< slack index (0.._nSlackVars) of over slack index for target goal, or DOESNT_EXIST if constraint

  op_regime       **pOperRegimes;  //< array of pointers to operating regimes, which are chosen from conditionals and determine active expression [size:nOperRegimes]
  int               nOperRegimes;  //< size of operating regime array

  int               active_regime;        //< currently active operating regime (or DOESNT_EXIST if none)
  bool              conditions_satisfied; //< true if any operating regime satisfied during current timestep
  bool              ever_satisfied;       //< true if any operating regime ever satisfied during simulation (for warning at end of sim)

  bool              use_stage_units;      //< true if goal penalty needs to be converted to units of s/m3
  double            units_correction;     //< penalty correction for units different than m3/s [default: 1]
  int               reservoir_index;      //< p index of reservoir for stage conversion

  managementGoal();
  ~managementGoal();
  expressionStruct *GetCurrentExpression() const{return pOperRegimes[nOperRegimes-1]->pExpression; }
  void AddOperatingRegime(op_regime *pOR, bool first);
  void AddOpCondition(exp_condition *pCondition); //adds to most recent operating regime
  void AddExpression (expressionStruct *pExp);   //adds to most recent operating regime
};

//////////////////////////////////////////////////////////////////
// workflow variable definition
//
struct workflowVar
{
  string            name;          //< workflow variable name

  double            current_val;   //< current value of workflow variable (evaluated at start of time step)

  op_regime       **pOperRegimes;  //< array of pointers to operating regimes, which are chosen from conditionals and determine active expression [size:nOperRegimes]
  int               nOperRegimes;  //< size of operating regime array

  workflowVar();
  ~workflowVar();
  expressionStruct *GetCurrentExpression() const{return pOperRegimes[nOperRegimes-1]->pExpression; }
  void AddOperatingRegime(op_regime *pOR, bool first);
  void AddOpCondition(exp_condition *pCondition); //adds to most recent operating regime
  void AddExpression (expressionStruct *pExp);   //adds to most recent operating regime
};

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for demand optimization
//
class CDemandGroup;
class CDemandOptimizer
{
private: /*------------------------------------------------------*/

  CModel          *_pModel;             //< pointer to model

  int              _nDecisionVars;      //< total number of decision variables considered
  decision_var   **_pDecisionVars;      //< array of pointers to decision variable structures [size:_nDecisionVars]

  int              _nWorkflowVars;      //< total number of workflow variables considered
  workflowVar    **_pWorkflowVars;      //< array of pointers to workflow variables [size: _nWorkflowVars]

  int              _nNonLinVars;        //> total number of non-linear variables (e.g., ?Q130)
  nonlin_var     **_pNonLinVars;        //> array of pointers to non-linear variable pairs
  int              _maxIterations;      //> maximum iterations in iterative scheme (default:5)
  double           _iterTolerance;      //> iterative solve tolerance (%)
  double           _relaxCoeff;         //> iterative solve relaxation coefficient (default:1.0=no relaxation)

  int              _nGoals;             //< number of user-defined enforced constraints/goals in management model
  managementGoal **_pGoals;             //< array of pointers to user-defined enforced constraints/goals in management model

  int              _nEnabledSubBasins;  //< number of enabled subbasins in model
  int             *_aSBIndices;         //< local index of enabled subbasins (0:_nEnabledSubBasins) [size: _nSubBasins]

  bool            *_aDisableSDCurve;    //< array of booleans to determine whether reservoir stage-discharge curve is overridden [size: _nSubBasins]

  int              _nReservoirs;        //< local storage of number of enabled lakes/reservoirs
  int             *_aResIndices;        //< storage of enabled reservoir indices (0:_nReservoirs or DOESNT_EXIST) [size:_nSubBasins]

  CDemand        **_pDemands;           //< array of pointers to water demand instances
  int              _nDemands;           //< local storage of number of demand locations
  int              _nReturns;           //< total number of demands with return flows (<=_nDemands)

  double          *_aDelivery;          //< array of delivered water for each water demand [m3/s] [size: _nDemands]
  double          *_aCumDelivery;       //< array of cumulative deliveries of demand since _aCumDivDate of this year [m3] [size: _nDemands]

  int            **_aUpstreamDemands;   //< demand indices (d) upstream (inclusive) of subbasin p [size: nSBs][size: _aUpCount] (only restricted demands)
  int             *_aUpCount;           //< number of demands upstream (inclusive) of subbasin p [size: nSBs]

  int              _nDemandGroups;      //< number of demand groups
  CDemandGroup   **_pDemandGroups;      //< array of pointers to demand groups

  int              _nEnviroFlowGoals;   //< number of *active* environmental flow goals
  double          *_aSlackValues;       //< array of slack variable values [size: _nSlackVars]
  int              _nSlackVars;         //< number of slack variables

  int             _nUserDecisionVars;   //< number of user-specified decision variables
  decision_var  **_pUserDecisionVars;   //< array of pointers to user-specified decision vars [size: _nUserDecisionVars] (overlaps content of _pDecisionVars)

  int             _nUserConstants;      //< number of user-specified named constants
  string         *_aUserConstNames;     //< array of names of user-specified constants
  double         *_aUserConstants;      //< array of values of user-specified constants

  int             _nUserTimeSeries;     //< number of user variable time series
  CTimeSeries   **_pUserTimeSeries;     //< array of pointers to user variable time series

  int             _nUserLookupTables;   //< number of user variable lookup tables
  CLookupTable  **_pUserLookupTables;   //< array of pointers to user variable lookup tables

  int             _nHistoryItems;      //< size of flow/demand/stage history that needs to be stored (in time steps) (from :LookbackDuration)
  double        **_aQhist;             //< history of subbasin discharge [size: _nEnabledSBs * _nHistoryItems]
  double        **_aDhist;             //< history of actual delivery [size: _nEnabledSBs * _nHistoryItems]
  double        **_aIhist;             //< history of reservoir inflows [size: _nEnabledSBs * _nHistoryItems]
  double        **_ahhist;             //< history of actual reservoir stages  (or 0 for non-reservoir basins) [size: _nEnabledSBs * _nHistoryItems]

  int             _nSolverResiduals;   //< number of residuals of solution (equal to lp_lib::nRows)
  double         *_aSolverResiduals;   //< vector of residuals of solution [size: _nSolverResiduals]
  string         *_aSolverRowNames;    //< vector of solver row names [size: _nSolverResiduals]

  ofstream        _MANOPT;          //< ofstream for ManagementOptimization.csv
  ofstream        _GOALSAT;            //< ofstream for GoalSatisfaction.csv
  ofstream        _RESID;              //< ofstream for OptimizationResidual.csv

  int             _do_debug_level;      //< =1 if debug info is to be printed to screen, =2 if LP matrix also printed (full debug), 0 for nothing

  //Called during simualtion
  void         UpdateHistoryArrays();
  void     UpdateWorkflowVariables(const time_struct &tt,const optStruct &Options);
  bool     ConvertToExpressionTerm(const string s, expressionTerm* term, const int lineno, const string filename)  const;
  int               GetDVColumnInd(const dv_type typ, const int counter) const;
  double              EvaluateTerm(expressionTerm **pTerms,const int k, const double &t) const;
  bool        EvaluateConditionExp(const expressionStruct* pE,const double &t) const;
  bool     CheckOpRegimeConditions(const op_regime *pOperRegime, const time_struct &tt, const optStruct &Options) const;


#ifdef _LPSOLVE_
  void            WriteLPSubMatrix(lp_lib::lprec *pLinProg, string filename, const optStruct &Options) const; //for debugging
  void           AddConstraintToLP(const int i, const int k, lp_lib::lprec *pLinProg, const time_struct &tt,int *col_ind, double *row_val, const bool update, const int lpgoalrow) const;
  void      IncrementAndSetRowName(lp_lib::lprec *pLinProg,int &rowcount,const string &name);
#endif


  //Called during initialization
  void              AddDecisionVar(const decision_var *pDV);
  bool        UserTimeSeriesExists(string TSname) const;
  void     AddReservoirConstraints(const optStruct &Options);
  void     IdentifyUpstreamDemands();
  bool          VariableNameExists(const string &name) const;
  int    GetNLIndexFromGuessString(const string &name) const;

public: /*------------------------------------------------------*/
  CDemandOptimizer(CModel *pMod);
  ~CDemandOptimizer();

  int           GetDemandIndexFromName(const string dname) const;
  double        GetNamedConstant      (const string s) const;
  double        GetUnitConversion     (const string s) const;
  int           GetUserDVIndex        (const string s) const;
  double        GetWorkflowVariable   (const string s, int &index) const;
  workflowVar  *GetWorkflowVarStruct  (int index);
  workflowVar  *GetWorkflowVarStruct  (const string s);
  int           GetNumUserDVs         () const;
  int           GetDebugLevel         () const;
  int           GetIndexFromDVString  (string s) const;

  CDemandGroup *GetDemandGroup        (const int ii);
  CDemandGroup *GetDemandGroupFromName(const string name);
  int           GetNumDemandGroups    () const;
  CDemand      *GetWaterDemand        (const int d);
  int           GetNumWaterDemands    () const;

  void   SetHistoryLength      (const int n);
  void   SetCumulativeDate     (const int julian_date, const string demandID);
  void   SetDebugLevel         (const int lev);
  void   SetDemandAsUnrestricted(const string dname);
  void   OverrideSDCurve       (const int p);
  void   SetMaxIterations      (const int Nmax);
  void   SetSolverTolerance    (const double tol);
  void   SetRelaxationCoeff    (const double relax);

  void   AddGoalOrConstraint   (const managementGoal *pGoal);

  void   AddUserDecisionVar    (const decision_var *pDV);
  void   SetUserDecisionVarBounds(const string name, const double &min, const double &max);
  void   AddUserConstant       (const string name, const double &val);
  void   AddWorkflowVariable   (const workflowVar *pCV);
  void   AddUserTimeSeries     (const CTimeSeries *pTS);
  void   AddUserLookupTable    (const CLookupTable *pLUT);
  void   AddNonLinVariable     (const string name, const string targetDV);

  void   AddWaterDemand        (CDemand *pDem);
  void   AddDemandGroup        (const string groupname);

  //void MultiplyDemand        (const string dname, const double &mult);
  //void MultiplyGroupDemand   (const string groupname, const double &mult);
  void SetDemandPenalty        (const string dname, const double &pen);
  //void SetDemandPriority     (const string dname, const int &prior);

  double     EvaluateExpression(const expressionStruct* pE,const double &t,bool RHS_only) const;

  expressionStruct *ParseExpression(const char **s, const int Len, const int lineno, const string filename) const;
  exp_condition    *ParseCondition (const char **s, const int Len, const int lineno, const string filename) const;

  void   Initialize            (CModel *pModel, const optStruct &Options);
  void   InitializePostRVMRead (CModel *pModel, const optStruct &Options);
  void   InitializeDemands     (CModel *pModel, const optStruct &Options);
  void   PrepDemandProblem     (CModel *pModel, const optStruct &Options, const time_struct &tt);
  void   SolveManagementProblem(CModel *pModel, const optStruct &Options, const double *aSBrunoff, const time_struct &tt);

  void   WriteOutputFileHeaders(const optStruct &Options);
  void   WriteMinorOutput      (const optStruct &Options,const time_struct &tt);
  void   CloseOutputStreams    ();

  void   Closure               (const optStruct &Options);
};
#endif
