/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2025 the Raven Development Team
  ----------------------------------------------------------------
  enum termtype
  expressionTerm struct
  expressionStruct structure
  ----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "TimeSeries.h"
#include "LookupTable.h"

#pragma once
///////////////////////////////////////////////////////////////////
/// \brief different expression term types
//
enum termtype
{
  TERM_DV,        //< decision variable !Q123 or named
  TERM_TS,        //< time series @ts(name,n)
  TERM_LT,        //< lookup table @lookup(x)
  TERM_DLT,       //< lookup table derivative @dlookup(x)
  TERM_HRU,       //< state variable @HRU_var(SNOW,2345)
  TERM_SB,        //< state variable @SB_var(SNOW,234)
  TERM_CONST,     //< constant
  TERM_WORKFLOW,  //< workflow variable (constant updated each time step)
  TERM_HISTORY,   //< bracketed - !Q123[-2]
  TERM_ITER,      //< decision variable guess ?Q123 or ?Q.RalphRiver
  TERM_MAX,       //< @max(x,y)
  TERM_MIN,       //< @min(x,y)
  TERM_CUMUL_TS,  //< @cumul(ts_name,duration) //MAY WANT @cumul(ts_name,duration,n) to handle time shift, e.g., 3 days to 10 days ago?
  TERM_POW,       //< @pow(x,n)
  TERM_CUMUL,     //< cumulative delivery !C123
  TERM_UNKNOWN    //< unknown
};

// -------------------------------------------------------------------
// data structures used by CDemandOptimizer class:
//   -expressionTerm
//   -expressionStruct (built from expressionTerms)
//   -decision_var
//   -exp_condition (may be defined using expressionStruct)
//   -op_regime (has conditions)
//   -workflowVar (defined using expressionStruct)
//   -managementGoal (built using multiple operating regimes)
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
  int           p_index;         //< subbasin index p (for history variables, e.g, !Q324[-2] or @SB_var() command )
  int           HRU_index;       //< HRU index k (for @HRU_var command)
  int           SV_index;        //< state variable index i (for @SB_var or @HRU_var command)

  expressionTerm(); //defined in DemandExpressionHandling.cpp
};

//////////////////////////////////////////////////////////////////
/// expression structure
///   abstraction of (A*B*C)+(D*E)-(F)+(G*H) <= 0
///   parenthetical collections are groups of terms -example has 4 groups with [3,2,1,2] terms per group
//
struct expressionStruct //full expression
{
  expressionTerm  ***pTerms;      //< 2D irregular array of pointers to expression terms size:[nGroups][nTermsPerGrp[j]]
  int                nGroups;     //< total number of terms groups in expression
  int               *nTermsPerGrp;//< number of terms per group [size: nGroups]
  comparison         compare;     //< comparison operator (==, <, >)

  bool               has_nonlin;

  string             origexp;     //< original string expression

  expressionStruct();
  ~expressionStruct();
};
