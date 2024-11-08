/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2024 the Raven Development Team
----------------------------------------------------------------*/
#include "DemandOptimization.h"

const int MAX_EXP_GROUPS     =50;
const int MAX_TERMS_PER_GROUP=10;

//////////////////////////////////////////////////////////////////
/// constructor and destructor of expression term structure
//
expressionTerm::expressionTerm()
{
  type =TERM_UNKNOWN;
  mult =1.0;
  reciprocal=false;
  origexp="";

  DV_ind=DOESNT_EXIST;
  value =0.0;
  pTS   =NULL;
  pLT   =NULL;
  timeshift=0;
  p_index  =DOESNT_EXIST;
  HRU_index=DOESNT_EXIST;
  SV_index =DOESNT_EXIST;

  is_nested  =false;
  nested_ind1=DOESNT_EXIST;
  nested_ind2=DOESNT_EXIST;
  nested_exp1 = "";
  nested_exp2 = "";
}

//////////////////////////////////////////////////////////////////
/// constructor and destructor of expression structure
//
expressionStruct::expressionStruct()
{
  pTerms=NULL;
  nTermsPerGrp=NULL;
  nGroups=0;
  compare=COMPARE_IS_EQUAL;
}
expressionStruct::~expressionStruct()
{
  for (int i = 0; i < nGroups; i++) {
    for (int j = 0; j < nTermsPerGrp[i]; j++) {
      delete pTerms[i][j]; pTerms[i][j]=NULL;
    }
    delete [] pTerms[i]; pTerms[i]=NULL;
  }
  delete [] pTerms; pTerms=NULL;
}

//////////////////////////////////////////////////////////////////
/// constructor, destructor, and member functions of managementGoal structure
//
managementGoal::managementGoal()
{
  name="";

  is_goal=false;
  penalty_under=1.0;
  penalty_over =1.0;
  slack_ind1=DOESNT_EXIST;
  slack_ind2=DOESNT_EXIST;

  nOperRegimes=1;
  pOperRegimes=new op_regime* [1];
  pOperRegimes[0] = new op_regime("[DEFAULT]");

  active_regime=DOESNT_EXIST;
  ever_satisfied=false;
  conditions_satisfied=false;

  use_stage_units=false;
  units_correction=1.0;
  reservoir_index=DOESNT_EXIST;
  overrides_SDcurve=false;
}
managementGoal::~managementGoal()
{
  delete [] pOperRegimes; nOperRegimes=0;
}

void managementGoal::AddOperatingRegime(op_regime* pOR, bool first)
{
  if ((nOperRegimes==1) && (first)){
    pOperRegimes[0]->reg_name=pOR->reg_name;
  }
  else{
    if(!DynArrayAppend((void**&)(pOperRegimes),(void*)(pOR),nOperRegimes)) {
      ExitGracefully("management_constraint::AddOperatingRegime: adding NULL operating regime",BAD_DATA_WARN);
    }
  }
}
void managementGoal::AddOpCondition(exp_condition* pCond)
{
  if(!DynArrayAppend((void**&)(pOperRegimes[nOperRegimes-1]->pConditions),(void*)(pCond),pOperRegimes[nOperRegimes-1]->nConditions)) {
    ExitGracefully("management_constraint::AddOpCondition: adding NULL condition",BAD_DATA_WARN);
  }
}
void managementGoal::AddExpression(expressionStruct* pExp)
{
  pOperRegimes[nOperRegimes-1]->pExpression=pExp;
  ExitGracefullyIf(pExp==NULL,"managementGoal::AddExpression: NULL Expression",RUNTIME_ERR);
}

//////////////////////////////////////////////////////////////////
/// \brief retrieves value of named constant from list of user constants
/// \params s [in] - string
/// \returns value of named constant, or BLANK if not found
//
double CDemandOptimizer::GetNamedConstant(const string s) const
{
  for (int i = 0; i < _nUserConstants; i++) {
    if (s==_aUserConstNames[i]){return _aUserConstants[i]; }
  }
  return RAV_BLANK_DATA;
}

//////////////////////////////////////////////////////////////////
/// \brief retrieves value of named constant from list of user constants
/// \params s [in] - string
/// \returns value of unit conversion multiplier, or BLANK if not found
//
double CDemandOptimizer::GetUnitConversion(const string s) const
{
  if (s == "ACREFTD_TO_CMS") {return 1.0/ACREFTD_PER_CMS;}
  if (s == "CMS_TO_ACREFTD") {return ACREFTD_PER_CMS;    }
  if (s == "INCHES_TO_MM"  ) {return MM_PER_INCH;        }
  if (s == "MM_TO_INCHES"  ) {return 1.0/MM_PER_INCH;    }
  if (s == "FEET_TO_METER" ) {return 1.0/FEET_PER_METER; }
  if (s == "METER_TO_FEET" ) {return FEET_PER_METER;     }
  if (s == "CMS_TO_CFS"    ) {return 1.0/CFS_PER_CMS;    }
  if (s == "CFS_TO_CMS"    ) {return CFS_PER_CMS;        }
  return RAV_BLANK_DATA;
}
//////////////////////////////////////////////////////////////////
/// \brief retrieves value of control variable from list of control variables
/// \params s [in] - string
/// \param index [out] - index of found control variable, or DOESNT_EXIST if not found
/// \returns value of control var, or BLANK if not found
//
double CDemandOptimizer::GetControlVariable(const string s, int &index) const
{
  index=DOESNT_EXIST;
  for (int i = 0; i < _nControlVars; i++) {
    if (s==_pControlVars[i]->name){index=i; return _pControlVars[i]->current_val; }
  }
  return RAV_BLANK_DATA;
}

//////////////////////////////////////////////////////////////////
/// \brief retrieves value of decision variable index
/// \params s [in] - string
/// \returns index of named user-specified decision var, or DOESNT_EXIST if not found
//
int CDemandOptimizer::GetUserDVIndex(const string s) const
{
  int ct=0;
  for (int i = 0; i < _nDecisionVars; i++) {
    if (_pDecisionVars[i]->dvar_type==DV_USER){
      if (_pDecisionVars[i]->name==s){return ct;}
      ct++;
    }
  }
  return DOESNT_EXIST;
}

//////////////////////////////////////////////////////////////////
/// \brief retrieves index of native decision variable starting with !
/// \params s [in] - string
/// \returns index of named decision variable, or DOESNT_EXIST if not found
/// index is subbasin index p for subbasin-based DVs, or demand index ii for demand-linked DVs
/// supports !Qxxx, !Q.name, !hxxx, !Ixxx, !Dxxx, !Cxxx, $Bxxx, $Exxx, !dxxx, !Rxxx
//
int CDemandOptimizer::GetIndexFromDVString(string s) const //String in format !Qxxx or !Q.name but not !Qxx[n]
{
  if ((s[1] == 'Q') || (s[1] == 'h') || (s[1]=='I') || (s[1]=='B') || (s[1]=='E')) //Subbasin-indexed \todo[funct] - support q = dQ/dt
  {
    if (s[2] == '.') {
      string name=s.substr(3);
      for (int p = 0; p < _pModel->GetNumSubBasins(); p++) {
        if (name == _pModel->GetSubBasin(p)->GetName()) {return p;}
      }
      return DOESNT_EXIST;
    }
    else{
      long   SBID=s_to_l(s.substr(2).c_str());
      return _pModel->GetSubBasinIndex(SBID);
    }
  }
  else if ((s[1]=='D') || (s[1]=='C') || (s[1]=='R') || (s[1] == 'd'))
  {
    string demandID;
    if (s[2] == '.') {demandID=s.substr(3);} //!D.FarmerBobsDemand
    else             {demandID=s.substr(2);} //!D1234
    return GetDemandIndexFromName(demandID);
  }
  return DOESNT_EXIST;
}

//////////////////////////////////////////////////////////////////
/// \brief returns number of user-specified decision variables
//
int CDemandOptimizer::GetNumUserDVs() const
{
  return _nUserDecisionVars;
}
//////////////////////////////////////////////////////////////////
/// \brief enum to string routines
//
string TermTypeToString(termtype t)
{
  if (t==TERM_DV      ){return "TERM_DV"; }
  if (t==TERM_TS      ){return "TERM_TS"; }
  if (t==TERM_LT      ){return "TERM_LT"; }
  if (t==TERM_HRU     ){return "TERM_HRU";}
  if (t==TERM_SB      ){return "TERM_SB"; }
  if (t==TERM_CONST   ){return "TERM_CONST"; }
  if (t==TERM_CONTROL ){return "TERM_CONTROL"; }
  if (t==TERM_HISTORY ){return "TERM_HISTORY"; }
  if (t==TERM_MAX     ){return "TERM_MAX"; }
  if (t==TERM_MIN     ){return "TERM_MIN"; }
  if (t==TERM_CONVERT ){return "TERM_CONVERT"; }
  if (t==TERM_CUMUL   ){return "TERM_CUMUL";}
  if (t==TERM_CUMUL_TS){return "TERM_CUMUL_TS";}
  if (t==TERM_UNKNOWN ){return "TERM_UNKNOWN"; }
  return "TERM_UNKNOWN";
}
string DVTypeToString(dv_type t)
{
  if (t==DV_QOUT    ){return "DV_QOUT";     }
  if (t==DV_QOUTRES ){return "DV_QOUTRES";  }
  if (t==DV_STAGE   ){return "DV_STAGE";    }
  if (t==DV_DSTAGE  ){return "DV_DSTAGE";   }
  if (t==DV_BINRES  ){return "DV_BINRES";   }
  if (t==DV_DELIVERY){return "DV_DELIVERY"; }
  if (t==DV_RETURN  ){return "DV_RETURN";   }
  if (t==DV_SLACK   ){return "DV_SLACK";    }
  if (t==DV_USER    ){return "DV_USER";     }
  return "DV_UNKNOWN";
}
//////////////////////////////////////////////////////////////////
/// \brief writes expType as string
/// called by SummarizeExpression
//
string expTypeToString(termtype &typ){
  switch (typ)
  {
    case(TERM_DV):      return "TERM_DV"; break;
    case(TERM_TS):      return "TERM_TS"; break;
    case(TERM_LT):      return "TERM_LT"; break;
    case(TERM_CONST):   return "TERM_CONST"; break;
    case(TERM_CONTROL): return "TERM_CONTROL"; break;
    case(TERM_HISTORY): return "TERM_HISTORY"; break;
    case(TERM_MAX):     return "TERM_MAX"; break;
    case(TERM_MIN):     return "TERM_MIN"; break;
    case(TERM_CONVERT): return "TERM_CONVERT"; break;
    case(TERM_CUMUL):   return "TERM_CUMUL"; break;
    case(TERM_UNKNOWN): return "TERM_UNKNOWN"; break;
  }
  return "?";
}
//////////////////////////////////////////////////////////////////
/// \brief parses individual expression string and converts to expression term structure
/// \param s        [in] - string
/// \param term    [out] - expression structure
/// \param lineno   [in] - line number of original expression in input file filename, referenced in errors
/// \param filename [in] - name of input file, referenced in errors
/// \returns expression term, false if entire expression should be ignored
/// only called during parse of expression term by ParseExpression() - no need to optimize
//
bool CDemandOptimizer::ConvertToExpressionTerm(const string s, expressionTerm* term, const int lineno, const string filename)  const
{
  size_t NPOS=std::string::npos;
  bool   has_brackets=false;
  int    n=0;
  string warn;
  string warnstring = " at line #" + to_string(lineno) + " in file "+filename;

  int    index;
  string tmp(s);
  term->origexp=tmp;
  term->p_index=DOESNT_EXIST;

  size_t i1=s.find("[");
  size_t i2=s.find("]");
  size_t i3=s.find("(");

  if ((i1 != NPOS) && (i2 != NPOS) && (i3==NPOS)) { //array NOT nested within function
    string contents=s.substr(i1+1,i2-i1-1);
    n=s_to_i(contents.c_str());
    //cout<<" BRACKET CONTENTS: ["<<contents<<"]"<<" n: "<<n<<endl;
    has_brackets=true;
  }

  //----------------------------------------------------------------------
  if ( ((s[0] == '!') || (s[0] == '$')) && (has_brackets)) //historical lookup (e.g., !Q102[-2])
  {
    ExitGracefullyIf(n>=0,              "ConvertToExpressionTerm:: term in square brackets must be <0 for historical variables.",BAD_DATA_WARN);
    ExitGracefullyIf(n<-_nHistoryItems, "ConvertToExpressionTerm:: term in square brackets exceeds the magnitude of the :LookbackDuration", BAD_DATA_WARN);

    if ((s[1] == 'Q') || (s[1] == 'D') || (s[1]=='I') || (s[1] == 'h') || (s[1] == 'B') || (s[1] == 'E'))
    {
      if (n<0)
      {
        int p=GetIndexFromDVString(s.substr(0,i1));
        if (p == DOESNT_EXIST) {
          warn="ConvertToExpressionTerm: invalid subbasin ID ("+s.substr(0,i1)+") in expression"+warnstring;
          ExitGracefully(warn.c_str(), BAD_DATA_WARN);
          return false;
        }
        if (!_pModel->GetSubBasin(p)->IsEnabled()) {
          //if (!_suppress_disabled_warnings){
          warn="ConvertToExpressionTerm: history expression "+warnstring+" includes reference to disabled subbasin. Expression will be ignored";
          WriteWarning(warn,true);
          return false;
        }
        term->type=TERM_HISTORY;
        term->timeshift=-n;
        term->p_index=p;
        return true;

      }
    }
    else {
      warn="ConvertToExpression:: Unparseable term in history expression starting with !/$ - only !Q, !I, !D, !h, $B, or $E currently supported. "+warnstring;
      ExitGracefully(warn.c_str(), BAD_DATA_WARN);
      return false;
    }
  }
  //----------------------------------------------------------------------
  else if      (s[0]=='!')  //decision variable e.g., !Q234, !I32, or !D.Matts_Brewery
  {
    if ((s[1] == 'Q') || (s[1] == 'h') || (s[1]=='I')) //Subbasin-indexed 
    {
      int p=GetIndexFromDVString(s);
      if (p == DOESNT_EXIST) {
        warn="ConvertToExpressionTerm: invalid subbasin ID ("+s+") in expression"+warnstring;
        ExitGracefully(warn.c_str(),BAD_DATA_WARN);
        return false;
      }
      if (!_pModel->GetSubBasin(p)->IsEnabled()) {
        //if (!_suppress_disabled_warnings){
        warn="ConvertToExpressionTerm: reference to disabled subbasin "+warnstring+". Constraint/goal will be disabled";
        WriteWarning(warn.c_str(),true);
        return false;
      }
      if (_aSBIndices[p] == DOESNT_EXIST) {
        //if (!_suppress_disabled_warnings){
        warn="ConvertToExpressionTerm: disabled subbasin ID in expression "+warnstring+ ". Expression will be ignored";
        WriteWarning(warn.c_str(),true);
        return false;
      }
      if ((s[1]=='I') && (_pModel->GetSubBasin(p)->GetReservoir()==NULL)) {
        warn="ConvertToExpressionTerm: reservoir inflow in expression for basin without reservoir "+warnstring+ ". Expression will be ignored";
        WriteWarning(warn.c_str(),true);
        return false;
      }
      if ((s[1]=='h') && (_pModel->GetSubBasin(p)->GetReservoir()==NULL)) {
        warn="ConvertToExpressionTerm: reservoir stage in expression for basin without reservoir "+warnstring+ ". Expression will be ignored";
        WriteWarning(warn.c_str(),true);
        return false;
      }
      if      ((s[1]=='Q') && (_aResIndices[p] != DOESNT_EXIST))
      {
        term->DV_ind=GetDVColumnInd(DV_QOUTRES, _aResIndices[p]);
      }
      else if ((s[1]=='Q') && (_aSBIndices[p]  != DOESNT_EXIST))
      {
        term->DV_ind=GetDVColumnInd(DV_QOUT,    _aSBIndices[p]);
      }
      else if (s[1]=='I')
      {
        term->DV_ind=GetDVColumnInd(DV_QOUT,   _aSBIndices[p]);
      }
      else if (s[1]=='h')
      {
        term->DV_ind=GetDVColumnInd(DV_STAGE,  _aResIndices[p]);
      }
      term->p_index=p;//not needed?
    }
    else if ((s[1]=='D') || (s[1]=='C') || (s[1]=='R') || (s[1]=='d'))   //demand indexed
    {
      int d=GetIndexFromDVString(s);

      if (d == DOESNT_EXIST) {
        warn="ConvertToExpressionTerm: invalid demand ID or demand from disabled subbasin used in expression " +warnstring +": goal/constraint will be ignored";
        WriteWarning(warn.c_str(),true);
        return false;
      }
      if      (s[1]=='D') //Delivery
      {
        term->DV_ind=GetDVColumnInd(DV_DELIVERY,d);
      }
      else if (s[1]=='R')
      {
        term->DV_ind=GetDVColumnInd(DV_RETURN,d);
      }
      else if (s[1]=='C')
      {
        term->type=TERM_CUMUL;
        term->p_index=d;
        return true;
      }
      else if (s[1] == 'd') //Demand
      {
        term->type=TERM_CONST;
        term->value=_pDemands[d]->GetDemand();
        return true;
      }
    }
    else{
      warn="ConvertToExpression:: Unparseable term in expression starting with ! - only !Q, !I, !D, !d, !R, !C, or !h currently supported. "+warnstring;
      ExitGracefully(warn.c_str(), BAD_DATA_WARN);
    }

    term->type=TERM_DV;
    return true;
  }
  //----------------------------------------------------------------------
  else if (s[0] == '$')
  {
    if ((s[1] == 'B') || (s[1] == 'E')) //Subbasin-indexed linked to time series
    {
      int p=GetIndexFromDVString(s);
      if (p == DOESNT_EXIST) {
        warn="ConvertToExpressionTerm: invalid subbasin ID ("+s+") in expression"+warnstring;
        ExitGracefully(warn.c_str(),BAD_DATA_WARN);
        return false;
      }
      if (!_pModel->GetSubBasin(p)->IsEnabled()) {
        warn="ConvertToExpressionTerm: reference to disabled subbasin "+warnstring+". Constraint/goal will be disabled";
        WriteWarning(warn.c_str(),true);
        return false;
      }
      if (_aSBIndices[p] == DOESNT_EXIST) {
        warn="ConvertToExpressionTerm: disabled subbasin ID in expression "+warnstring+ ". Expression will be ignored";
        WriteWarning(warn.c_str(),true);
        return false;
      }
      term->type=TERM_HISTORY;
      term->timeshift=0;
      term->p_index=p;
      return true;
    }
  }
  //----------------------------------------------------------------------
  else if (s.substr(0, 4) == "@ts(")//time series (e.g., @ts(my_time_series,n)
  {
    string name;
    int    index;
    size_t is = s.find("@ts(");
    size_t ie = s.find(",",is);
    size_t ip = s.find_last_of(")");
    if (ie == NPOS) {
      warn="ConvertToExpressionTerm: missing comma in @ts expression"+warnstring;
      ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    }
    if (ip == NPOS) {
      warn="ConvertToExpressionTerm: missing end parentheses in @ts expression"+warnstring;
      ExitGracefully(warn.c_str(), BAD_DATA_WARN);
    }
    if ((is != NPOS) && (ie != NPOS))
    {
      bool found=false;
      name  = s.substr(is+4,ie-(is+4));
      index = s_to_i(s.substr(ie+1,ip-(ie+1)).c_str());
      term->pTS=NULL;
      for (int i = 0; i < _nUserTimeSeries; i++) {
        if (StringToUppercase(_pUserTimeSeries[i]->GetName()) == StringToUppercase(name)) {
          term->pTS = _pUserTimeSeries[i];
          found=true;
        }
      }
      term->timeshift=index;
      term->type     =TERM_TS;
      if (!found) {
        warn="ConvertToExpression: unrecognized time series name in @ts command"+warnstring;
        ExitGracefully(warn.c_str(), BAD_DATA_WARN);
      }
    }
  }
  //----------------------------------------------------------------------
  else if (s.substr(0, 7) == "@cumul(") //cumulative time series (e.g., @cumul(my_time_series,duration))
  {
    string name;
    int    index;
    size_t is = s.find("@cumul(");
    size_t ie = s.find(",",is);
    size_t ip = s.find_last_of(")");
    if (ie == NPOS) {
      warn="ConvertToExpressionTerm: missing comma in @cumul expression"+warnstring;
      ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    }
    if (ip == NPOS) {
      warn="ConvertToExpressionTerm: missing end parentheses in @cumul expression"+warnstring;
      ExitGracefully(warn.c_str(), BAD_DATA_WARN);
    }
    if ((is != NPOS) && (ie != NPOS))
    {
      bool found=false;
      name  = s.substr(is+7,ie-(is+7));
      index = s_to_i(s.substr(ie+1,ip-(ie+1)).c_str());
      term->pTS=NULL;
      for (int i = 0; i < _nUserTimeSeries; i++) {
        if (StringToUppercase(_pUserTimeSeries[i]->GetName()) == StringToUppercase(name)) {
          term->pTS = _pUserTimeSeries[i];
          found=true;
        }
      }
      term->timeshift=index;
      term->type     =TERM_CUMUL_TS;
      if (!found) {
        warn="ConvertToExpression: unrecognized time series name in @cumul command"+warnstring;
        ExitGracefully(warn.c_str(), BAD_DATA_WARN);
      }
    }
  }
  //----------------------------------------------------------------------
  else if (s.substr(0, 8) == "@lookup(")//lookup table (e.g., @lookup(my_table,EXPRESSION) -NOW SUPPORTS NESTED FUNCTIONS
  {
    string name;
    string x_in;
    size_t is = s.find("@lookup(");
    size_t ie = s.find(",",is);
    size_t ip = s.find_last_of(")");
    if (ie == NPOS) {
      warn="ConvertToExpressionTerm: missing comma in @lookup expression"+warnstring;
      ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    }
    if (ip == NPOS) {
      warn="ConvertToExpressionTerm: missing end parentheses in @lookup expression"+warnstring;
      ExitGracefully(warn.c_str(), BAD_DATA_WARN);
    }
    if ((is != NPOS) && (ie != NPOS))
    {
      bool found=false;
      name = s.substr(is+8,ie-(is+8));
      x_in = s.substr(ie+1,ip-(ie+1));
      term->pLT=NULL;
      for (int i = 0; i < _nUserLookupTables; i++) {
        if (StringToUppercase(_pUserLookupTables[i]->GetName()) == StringToUppercase(name)) {
          term->pLT = _pUserLookupTables[i];
          found=true;
        }
      }
      term->nested_exp1 =x_in;
      term->type        =TERM_LT;
      if (!found) {
        warn="ConvertToExpression: unrecognized lookup table name in @lookup command"+warnstring;
        ExitGracefully(warn.c_str(), BAD_DATA_WARN);
        return false;
      }
    }
  }
 //----------------------------------------------------------------------
 else if (s.substr(0, 9) == "@HRU_var(") // HRU state var (e.g., @HRU_var(SNOW,[id])
  {
    string sv_name;
    long long int HRUID;
    size_t is = s.find("@HRU_var(");
    size_t ie = s.find(",",is);
    size_t ip = s.find_last_of(")");
    if (ie == NPOS) {
      warn="ConvertToExpressionTerm: missing comma in @HRU_var expression"+warnstring;
      ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    }
    if (ip == NPOS) {
      warn="ConvertToExpressionTerm: missing end parentheses in @HRU_var expression"+warnstring;
      ExitGracefully(warn.c_str(), BAD_DATA_WARN);
    }
    if ((is != NPOS) && (ie != NPOS))
    {
      bool found=false;
      sv_name = s.substr(is+9,ie-(is+9));
      HRUID   = s_to_i(s.substr(ie+1,ip-(ie+1)).c_str());

      int lay=DOESNT_EXIST;
      sv_type sv=_pModel->GetStateVarInfo()->StringToSVType(sv_name,lay,false);
      if (sv== UNRECOGNIZED_SVTYPE)  {
        warn="ConvertToExpression: requested state variable in in @HRU_var command"+warnstring;
        ExitGracefully(warn.c_str(), BAD_DATA_WARN);
        return false;
      }
      if (_pModel->GetHRUByID(HRUID) == NULL) {
        warn="ConvertToExpression: requested state variable in in @HRU_var command"+warnstring;
        ExitGracefully(warn.c_str(), BAD_DATA_WARN);
        return false;
      }

      term->SV_index=_pModel->GetStateVarIndex(sv,lay);
      term->HRU_index=_pModel->GetHRUByID(HRUID)->GetGlobalIndex();
      term->type=TERM_HRU;
    }
  }
  //----------------------------------------------------------------------
  else if (s.substr(0, 8) == "@SB_var(") // SubBasin state var (e.g., @SB_var(SNOW,[id])
  {
    string sv_name;
    long   SBID;
    size_t is = s.find("@SB_var(");
    size_t ie = s.find(",",is);
    size_t ip = s.find_last_of(")");
    if (ie == NPOS) {
      warn="ConvertToExpressionTerm: missing comma in @SB_var expression"+warnstring;
      ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    }
    if (ip == NPOS) {
      warn="ConvertToExpressionTerm: missing end paretheses in @SB_var expression"+warnstring;
      ExitGracefully(warn.c_str(), BAD_DATA_WARN);
    }
    if ((is != NPOS) && (ie != NPOS))
    {
      bool found=false;
      sv_name = s.substr(is+8,ie-(is+8));
      SBID    = s_to_l(s.substr(ie+1,ip-(ie+1)).c_str());

      int lay=DOESNT_EXIST;
      sv_type sv=_pModel->GetStateVarInfo()->StringToSVType(sv_name,lay,false);

      if (sv== UNRECOGNIZED_SVTYPE) {
        warn="ConvertToExpression: requested state variable in in @SB_var command "+warnstring;
        ExitGracefully(warn.c_str(), BAD_DATA_WARN);
        return false;
      }
      if (_pModel->GetSubBasinByID(SBID) == NULL) {
        warn="ConvertToExpression: requested state variable in in @SB_var command "+warnstring;
        ExitGracefully(warn.c_str(), BAD_DATA_WARN);
        return false;
      }
      term->SV_index=_pModel->GetStateVarIndex(sv,lay);
      term->p_index=_pModel->GetSubBasinByID(SBID)->GetGlobalIndex();
      term->type=TERM_SB;
    }
  }
  //----------------------------------------------------------------------
  else if (s.substr(0, 5) == "@max(") //max function (e.g., @max(exp1,exp2)
  {
    string name;
    string x_in,y_in;
    size_t is = s.find("@max(");
    size_t ie = s.find(",",is);
    size_t ip = s.find_last_of(")");
    if (ie == NPOS) {
      warn="ConvertToExpressionTerm: missing comma in @max expression"+warnstring;
      ExitGracefully(warn.c_str(),BAD_DATA_WARN);
      return false;
    }
    if (ip == NPOS) {
      warn="ConvertToExpressionTerm: missing end parentheses in @max expression"+warnstring;
      ExitGracefully(warn.c_str(), BAD_DATA_WARN);
      return false;
    }
    if ((is != NPOS) && (ie != NPOS))
    {
      x_in = s.substr(is+5,ie-(is+5));
      y_in = s.substr(ie+1,ip-(ie+1));
      term->nested_exp1 =x_in;
      term->nested_exp2 =y_in;
      term->type        =TERM_MAX;
    }
  }
  //----------------------------------------------------------------------
  else if (s.substr(0, 5) == "@min(") //max function (e.g., @min(exp1,exp2)
  {
    string name;
    string x_in,y_in;
    size_t is = s.find("@min(");
    size_t ie = s.find(",",is); //\todo[funct] - handle nested functions as first argument (this only works with second argument)
    size_t ip = s.find_last_of(")");
    if (ie == NPOS) {
      warn="ConvertToExpressionTerm: missing comma in @min expression"+warnstring;
      ExitGracefully(warn.c_str(),BAD_DATA_WARN);
      return false;
    }
    if (ip == NPOS) {
      warn="ConvertToExpressionTerm: missing end parentheses in @min expression"+warnstring;
      ExitGracefully(warn.c_str(), BAD_DATA_WARN);
      return false;
    }
    if ((is != NPOS) && (ie != NPOS))
    {
      x_in = s.substr(is+5,ie-(is+5));
      y_in = s.substr(ie+1,ip-(ie+1));
      term->nested_exp1 =x_in;
      term->nested_exp2 =y_in;
      term->type        =TERM_MIN;
    }
  }
  //----------------------------------------------------------------------
  else if (GetUnitConversion(s)!=RAV_BLANK_DATA) // named unit conversion
  {
    term->type=TERM_CONST;
    term->value=GetUnitConversion(s);
  }
  //----------------------------------------------------------------------
  else if (GetNamedConstant(s)!=RAV_BLANK_DATA) // named constant
  {
    term->type=TERM_CONST;
    term->value=GetNamedConstant(s);
  }
  //----------------------------------------------------------------------
  else if (GetControlVariable(s,index) != RAV_BLANK_DATA) // control variable
  {
    
    term->type=TERM_CONTROL;
    term->value=GetControlVariable(s,index);// initially zero
    term->DV_ind=index;
  }
  //----------------------------------------------------------------------
  else if (GetUserDVIndex(s)!=DOESNT_EXIST) // user-defined decision variable
  {
    term->DV_ind=GetDVColumnInd(DV_USER, GetUserDVIndex(s));
    term->type=TERM_DV;
  }
  //----------------------------------------------------------------------
  else if (is_numeric(s.c_str()))  // numeric value (e.g., 132.4)
  {
    term->type=TERM_CONST;
    term->value=s_to_d(s.c_str());
  }
  //----------------------------------------------------------------------
  else {
    warn="ConvertToExpression:: Unparseable term "+ s + " in expression" + warnstring;
    ExitGracefully(warn.c_str(), BAD_DATA_WARN);
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////
/// \brief Parses demand constraint expression of form (e.g.,) A * B(x) * C + D * E - F = G * H(t) [NO PARENTHESES!]
/// \param s        [in] - array of strings of [size: Len]
/// \param Len      [in] - length of string array
/// \param lineno   [in] - line number of original expression in input file filename, referenced in errors
/// \param filename [in] - name of input file, referenced in errors
/// \returns expressionStruct: a 2D array of pointers to grouped terms (e.g., [0]:[A,B,C], [1]:[D,E], [2]:[F], [3]:[G,H] for above example)
//
expressionStruct *CDemandOptimizer::ParseExpression(const char **s,
                                                    const int    Len,
                                                    const int    lineno,
                                                    const string filename) const
{
  int rhs_ind=-1;
  exptype          type       [MAX_EXP_GROUPS];
  expressionTerm  *terms      [MAX_EXP_GROUPS][MAX_TERMS_PER_GROUP];
  int              termspergrp[MAX_EXP_GROUPS];
  comparison       comp=COMPARE_BETWEEN;
  int              nComps=0;
  size_t           strlen;

  //identify all strings as either operators or terms
  //determine nature of comparison =,<,>
  for (int i = 1; i < Len; i++) //starts at 1 - first term is :Expression or :DefineDecisionVariable
  {
    strlen=to_string(s[i]).length();
    type[i]=EXP;
    if ((s[i][0]=='+') || ((strlen==1) && (s[i][0]=='-')) || (s[i][0]=='*') || (s[i][0]=='/') || (s[i][0]=='=') || (s[i][0]=='<') || (s[i][0]=='>')){ 
      type[i] = EXP_OP;
      if ((i > 1) && (type[i - 1] == EXP_OP)) {
        ExitGracefully("ParseExpression: cannot have consecutive math operators in an expression.",BAD_DATA_WARN);
      }
    }
    if (s[i-1][0]=='/'){type[i]=EXP_INV; }

    if ((s[i][0]=='=') || (s[i][0]=='~')){
      comp=COMPARE_IS_EQUAL;    rhs_ind=i; nComps++;
    }
    else if (s[i][0]=='<'){//works for <=, too
      comp=COMPARE_LESSTHAN;    rhs_ind=i; nComps++;
    }
    else if (s[i][0]=='>'){ //works for >=, too
      comp=COMPARE_GREATERTHAN; rhs_ind=i; nComps++;
    }
  }
  if (nComps > 1) {
    ExitGracefully("ParseExpression: more than one comparison operator (e.g., <, >, =) in expression",BAD_DATA_WARN);
  }
  if (nComps==0) {
    ExitGracefully("ParseExpression: missing comparison operator (e.g., <, >, =) in expression",BAD_DATA_WARN);
  }

  //TMP DEBUG
  /*cout << "parse expression (2) comparison type = " << comp << " rhs ind= " << rhs_ind << " types=";
  for (int i = 1; i < Len; i++) //starts at 1 - first term is :Expression or :DefineDecisionVariable
  {
    if (type[i]==EXP_OP ){ cout<<" OP"; }
    if (type[i]==EXP    ){ cout<<" EXP"; }
    if (type[i]==EXP_INV){ cout<<" INV"; }
  }
  cout<<endl;*/
  //END TMP DEBUG

  int  j=0;
  int  k=0;
  int  DVcount=0;
  bool valid;
  for (int i = 1; i < Len; i++)
  {
    if (type[i] != EXP_OP)
    {
      terms[j][k] = new expressionTerm();
      terms[j][k]->mult=1.0;
      valid=ConvertToExpressionTerm(to_string(s[i]),terms[j][k],lineno,filename);
      if (!valid){return NULL; }

      //TMP DEBUG
      //cout<<"   TERM["<<k<<" in "<< j << "]: " << s[i] << " : type = " << TermTypeToString(terms[j][k]->type) << " index : " << terms[j][k]->DV_ind << endl;

      if ((i>rhs_ind)      && (k==0)){terms[j][k]->mult*=-1.0; } //only multiply leading term by -1 to move to LHS
      if ((s[i-1][0]=='-') && (k==0)){terms[j][k]->mult*=-1.0; } //and reverse if preceded with minus sign

      if (terms[j][k]->type == TERM_DV) {
        DVcount++;
        if (DVcount > 1) {
          ExitGracefully("ParseExpression: more than one decision variable found in product term within constraint expression. Cannot handle non-linear expressions.",BAD_DATA_WARN);
          return NULL;
        }
      }
      terms[j][k]->reciprocal=(type[i] == EXP_INV);

      //Handle nested expressions
      if (terms[j][k]->nested_exp1 != "")
      {
        terms[j][k  ]->nested_ind1=k+1;

        terms[j][k+1] = new expressionTerm();
        terms[j][k+1]->mult=1.0;
        terms[j][k+1]->is_nested=true;
        valid=ConvertToExpressionTerm(terms[j][k]->nested_exp1,terms[j][k+1],lineno, filename);
        if (!valid){return NULL;}

        if (terms[j][k + 1]->type == TERM_DV) {
          ExitGracefully("ParseExpression: decision variable found nested in functional expression (lookup, min, or max). Cannot handle non-linear expressions.",BAD_DATA_WARN);
          return NULL;
        }

        if (terms[j][k]->nested_exp2 != "")
        {
          terms[j][k  ]->nested_ind2=k+2;

          terms[j][k+2] = new expressionTerm();
          terms[j][k+2]->mult=1.0;
          terms[j][k+2]->is_nested=true;
          valid=ConvertToExpressionTerm(terms[j][k]->nested_exp2, terms[j][k + 2],lineno, filename);
          if (!valid){return NULL;}

          if (terms[j][k + 2]->type == TERM_DV) {
            ExitGracefully("ParseExpression: decision variable found nested in functional expression (min or max). Cannot handle non-linear expressions.",BAD_DATA_WARN);
            return NULL;
          }
          k++;
        }
        k++;
      }//end nested expressions
      k++;
    }
    else if ((i != 1) && (s[i][0]!='*') && (s[i][0]!='/')) { // +, -, or = operator - expression product term finished and expression doesn't start with + or -
      termspergrp[j]=k;
      j++; k=0; DVcount=0;
    }
  }
  termspergrp[j]=k; //last group
  j++;

  //copy to new expressionStruct structure
  //---------------------------------------------------------------
  int nGroups=j;
  expressionStruct *tmp;
  tmp=new expressionStruct();
  tmp->nGroups=nGroups;
  tmp->compare=comp;

  string origexp="";
  for (int i = 1; i < Len; i++) {
    origexp+=to_string(s[i])+" ";
  }
  tmp->origexp=origexp;

  tmp->pTerms=new expressionTerm **[nGroups];
  tmp->nTermsPerGrp=new int[nGroups];
  for (int j = 0; j < nGroups; j++) {
    tmp->pTerms[j] = new expressionTerm * [termspergrp[j]];
    tmp->nTermsPerGrp[j]=termspergrp[j];
    for (int k = 0; k < termspergrp[j]; k++) {
      tmp->pTerms[j][k]=terms[j][k];
    }
  }

  return tmp;
}


//////////////////////////////////////////////////////////////////
/// \brief writes out explicit details about expressionStructure contents
/// called in ParseManagementFile under Noisy mode
//
void SummarizeExpression(const char **s, const int Len, expressionStruct* exp)
{
  if (exp==NULL){return;}
  cout<<"*"<<endl;
  cout<<"EXPRESSION: "<<endl;
  for (int i = 0; i < Len; i++) {cout<<s[i]<<" ";}cout<<endl;
  cout<<"--# groups="<<exp->nGroups<<"-->"<<endl;
  for (int i = 0; i < exp->nGroups; i++) {
    cout << "  GROUP "<<i+1<<"(" <<exp->nTermsPerGrp[i]<<" terms): " << endl;
    for (int j = 0; j < exp->nTermsPerGrp[i]; j++) {
      cout<<"   TERM "<<j+1<<": " << exp->pTerms[i][j]->origexp<<" ("<< expTypeToString(exp->pTerms[i][j]->type) <<") mult: "<<exp->pTerms[i][j]->mult<<endl;
    }
  }
  cout<<"*"<<endl<<endl;
}
//////////////////////////////////////////////////////////////////
/// \brief checks if conditions of operating regime k or goal ii are satisfied
/// returns true if *all* conditions of operating regime k of goal ii are satisfied
///
bool CDemandOptimizer::CheckGoalConditions(const int ii, const int k, const time_struct &tt, const optStruct &Options) const
{
  double dv_value;
  managementGoal *pC=_pGoals[ii];

  //Check if conditionals are satisfied
  for (int j = 0; j < pC->pOperRegimes[k]->nConditions; j++)
  {
    exp_condition *pCond=pC->pOperRegimes[k]->pConditions[j];
    if (pCond->pExp != NULL) {
      if(!EvaluateConditionExp(pCond->pExp,tt.model_time))
      {
        return false;
      }
    }
    else
    {
      dv_value=0.0;
      if      (pCond->dv_name == "DAY_OF_YEAR")
      {
        dv_value=tt.julian_day;
      }
      else if (pCond->dv_name == "MONTH")
      {
        dv_value=(double)(tt.month);
      }
      else if (pCond->dv_name == "DATE")
      {
        dv_value=0;
        //handled below
      }
      else if (pCond->dv_name[0] == '!') //decision variable
      {
        char   tmp =pCond->dv_name[1];

        long   ind=pCond->p_index;

        if      (tmp == 'Q') {
          dv_value=_pModel->GetSubBasinByID(ind)->GetOutflowRate();
        }
        else if (tmp == 'h') {
          dv_value=_pModel->GetSubBasinByID(ind)->GetReservoir()->GetResStage();
        }
        else if (tmp == 'I') {
          dv_value=_pModel->GetSubBasinByID(ind)->GetOutflowArray()[_pModel->GetSubBasinByID(ind)->GetNumSegments()-1];
        }
        else if (tmp == 'C') {
          dv_value=_aCumDelivery[ind];
        }
        else if (tmp == 'D') {
          dv_value=_pModel->GetSubBasinByID(_pDemands[ind]->GetSubBasinID())->GetDemandDelivery(_pDemands[ind]->GetLocalIndex());
        }
        else if (tmp == 'R') {
          dv_value=_pModel->GetSubBasinByID(_pDemands[ind]->GetSubBasinID())->GetReturnFlow(_pDemands[ind]->GetLocalIndex());
        }
        else if (tmp == 'd') {
          dv_value=_pDemands[ind]->GetDemand();
        }
        else {
          ExitGracefully("Invalid decision variable in condition statement (letter after ! not supported)",BAD_DATA);
        }
        //todo: support !q, !B, !E, !F, !T
      }
      else {//handle user specified DVs and control variables
        int i=GetUserDVIndex(pCond->dv_name);
        if (i == DOESNT_EXIST) {
          bool found=false;
          for (int j = 0; j < _nControlVars; j++) {
            if (_pControlVars[j]->name == pCond->dv_name) {
              dv_value =_pControlVars[i]->current_val;
              found=true;
            }
          }

          if (!found){
            ExitGracefully("CheckGoalConditions: Unrecognized varible on left hand side of :Condition statement ",BAD_DATA_WARN);
            return false;
          }
        }
        else {
          dv_value =_pDecisionVars[i]->value;
        }
      }

      comparison comp=pCond->compare;
      double        v=pCond->value;
      double       v2=pCond->value2;

      if (pCond->dv_name == "DATE") //dates require special handling; remainder of values are just doubles
      {
        dv_value=0;

        time_struct time1=DateStringToTimeStruct(pCond->date_string ,"00:00:00",Options.calendar);
        v=-TimeDifference(time1.julian_day,time1.year,tt.julian_day,tt.year,Options.calendar);//v negative if day is after day1

        if (comp == COMPARE_BETWEEN) {
          time_struct time2=DateStringToTimeStruct(pCond->date_string2,"00:00:00",Options.calendar);
          v2=-TimeDifference(time2.julian_day,time2.year,tt.julian_day,tt.year,Options.calendar);//v2 negative if day is after day2
        }
      }

      if (comp == COMPARE_BETWEEN)
      {
        if ((pCond->dv_name == "DAY_OF_YEAR") || (pCond->dv_name =="MONTH")) { //handles wraparound
          if ( v2 < v ){
            if ((dv_value > v ) && (dv_value < v2)){return false;}
          }
          else {
            if ((dv_value > v2) || (dv_value < v )){return false;}
          }
        }
        else {
          if ((dv_value > v2) || (dv_value< v)){return false;}//don't apply condition if ANY conditionals unmet
        }
      }
      else if (comp == COMPARE_GREATERTHAN) {
        if (dv_value < v){return false;}//don't apply condition
      }
      else if (comp == COMPARE_LESSTHAN) {
        if (dv_value > v){return false;}//don't apply condition
      }
      else if (comp == COMPARE_IS_EQUAL) {
        if (fabs(dv_value - v) > PRETTY_SMALL){return false;}//usually integer value (e.g., month)
      }
      else if (comp == COMPARE_NOT_EQUAL) {
        if (fabs(dv_value - v) < PRETTY_SMALL){return false;}
      }
    }
  }
  return true; //all conditionals satisfied
}

//////////////////////////////////////////////////////////////////
/// adds constraint ii to LP solve problem statement
/// \params ii [in] - index of constraint in _pGoals array
/// \param k - index of operating regime in _pGoals[ii]
/// \param pLinProg [out] - pointer to valid lpsolve structure to be modified
/// \param tt [in] - time structure
/// \param *col_ind [in] - empty array (with memory reserved) for storing column indices
/// \param *row_val [in] - empty array (with memory reserved) for storing row values
//
#ifdef _LPSOLVE_
void CDemandOptimizer::AddConstraintToLP(const int ii, const int kk, lp_lib::lprec* pLinProg, const time_struct &tt, int *col_ind, double *row_val) const
{
  double coeff;
  int    i=0;
  int    retval;
  double RHS,term;
  bool   group_has_dv;
  int    DV_ind;
  bool   constraint_valid=true;

  managementGoal   *pC=_pGoals[ii];
  expressionStruct *pE;

  if (kk!=DOESNT_EXIST)
  {
    pE= pC->pOperRegimes[kk]->pExpression; //active expression

    RHS=0.0;
    for (int j = 0; j < pE->nGroups; j++)
    {
      coeff=1.0;
      group_has_dv=false;
      DV_ind=DOESNT_EXIST;

      for (int k = 0; k < pE->nTermsPerGrp[j]; k++)
      {
        if (pE->pTerms[j][k]->type == TERM_DV)
        {
          DV_ind=pE->pTerms[j][k]->DV_ind;
          group_has_dv=true;
          coeff *= (pE->pTerms[j][k]->mult);
        }
        else if (!(pE->pTerms[j][k]->is_nested))
        {
          term=EvaluateTerm(pE->pTerms[j], k, tt.model_time);

          if (term==RAV_BLANK_DATA){constraint_valid=false;}
          if (pE->pTerms[j][k]->reciprocal == true)
          {
            coeff /= (pE->pTerms[j][k]->mult) * term;
            ExitGracefullyIf(term==0.0, "AddConstraintToLP: Divide by zero error in evaluating :Expression with division term",BAD_DATA);
          }
          else {
            coeff *= (pE->pTerms[j][k]->mult) * term;
          }
        }
      }

      if (!group_has_dv) {
        RHS-=coeff; //term group goes on right hand side
      }
      else {        //term group multiplies decision variable
        col_ind[i]=DV_ind;
        row_val[i]=coeff;
        i++;
      }
    }
  }
  else {
    pE= pC->pOperRegimes[0]->pExpression; //inactive expression
  }

  if ((kk==DOESNT_EXIST) || (!constraint_valid)) // INACTIVE/INVALID GOAL/CONSTRAINT
  {
    if (!pC->is_goal) {
      col_ind[i]=1;
      row_val[i]=0.0;
      i++;
    }
    RHS=0; //leads to inert equation s1-s2=0, s1>=0, or s1<=0 OR 0*dv1=0 if constraint
  }

  int constr_type=ROWTYPE_EQ;
  int nSlack=1;
  coeff=1.0;
  if      (pE->compare==COMPARE_IS_EQUAL)   {constr_type=ROWTYPE_EQ; coeff=-1.0; nSlack=2;}
  else if (pE->compare==COMPARE_LESSTHAN)   {constr_type=ROWTYPE_LE; coeff=-1.0;}
  else if (pE->compare==COMPARE_GREATERTHAN){constr_type=ROWTYPE_GE; coeff=+1.0;}

  //Add slack variables to goals
  if (pC->is_goal) {
    col_ind[i]=GetDVColumnInd(DV_SLACK, pC->slack_ind1);
    row_val[i]=coeff;
    i++;
    if (nSlack == 2) {
      col_ind[i]=GetDVColumnInd(DV_SLACK, pC->slack_ind2);
      row_val[i]=-coeff;
      i++;
    }
  }

  retval = lp_lib::add_constraintex(pLinProg,i,row_val,col_ind,constr_type,RHS);
  ExitGracefullyIf(retval==0,"AddConstraintToLP::Error adding user-specified constraint/goal",RUNTIME_ERR);
}
#endif

//////////////////////////////////////////////////////////////////
/// evaluates difference between right hand side and left hand side of expression, or value of RHS if RHS_only==true and statement only has one term left of (in)equality
/// does not support active decision variables in expression (only to be used for conditionals and demand/return expressions)
/// \params pE [in] - conditional expression
/// \param t [in] - current model time
/// \param RHS_only [in] - true if only Right hand side of expression is to be evaluated, else RHS-LHS is returned
/// repeatedly calls EvaluateTerm()
/// \returns RHS if RHS_only or RHS-LHS if !RHS_only; returns BLANK if any expression is blank (usually time series with blank value)
//
double CDemandOptimizer::EvaluateExpression(const expressionStruct* pE,const double &t,bool RHS_only) const
{
  double coeff;
  double term;
  double RHS=0.0;
  bool constraint_valid=true;
  for (int j = 0; j < pE->nGroups; j++)
  {
    coeff=1.0;
    if (RHS_only && j==0) {
      //skip first term (assumes expression is of form  A = B + C - D /E, i.e., only one term on left
    }
    else{
      for (int k = 0; k < pE->nTermsPerGrp[j]; k++)
      {
        if (pE->pTerms[j][k]->type == TERM_DV)
        {
          string warn="EvaluateConditionalExp: conditional expressions or demand/return expressions cannot contain decision variables. Problematic command: "+pE->origexp;
          ExitGracefully(warn.c_str(), BAD_DATA);
        }
        else if (!(pE->pTerms[j][k]->is_nested))
        {
          term=EvaluateTerm(pE->pTerms[j], k, t);

          if (term==RAV_BLANK_DATA){return RAV_BLANK_DATA;} //returns blank if missing data in expression
          if (pE->pTerms[j][k]->reciprocal == true)
          {
            coeff /= (pE->pTerms[j][k]->mult) * term;
            ExitGracefullyIf(term==0.0, "EvaluateConditionExp: Divide by zero error in evaluating :Condition or :DemandExpression expression with division term",BAD_DATA);
          }
          else {
            coeff *= (pE->pTerms[j][k]->mult) * term;
          }
        }
      }

      RHS-=coeff; //term group goes on right hand side
    }
  }
  return RHS;
}

//////////////////////////////////////////////////////////////////
/// evaluates conditional expression
/// \params pE [in] - conditional expression
/// \param t [in] - current model time
/// \returns true if conditional expression is satisfied
//
bool CDemandOptimizer::EvaluateConditionExp(const expressionStruct* pE,const double &t) const
{

  double RHSminusLHS=EvaluateExpression(pE,t,false);

  if      (fabs(RHSminusLHS-RAV_BLANK_DATA)<REAL_SMALL) {return true;}//condition assumed satisfied if no data in conditional
  if      (pE->compare==COMPARE_IS_EQUAL)               {return (fabs(RHSminusLHS)<REAL_SMALL);}
  else if (pE->compare==COMPARE_LESSTHAN)               {return (RHSminusLHS>0);}
  else if (pE->compare==COMPARE_GREATERTHAN)            {return (RHSminusLHS<0);}
  return false;
}

//////////////////////////////////////////////////////////////////
/// evaluates term in constraint/goal expression to numerical value (unless term is DV)
/// \params pTerms [in] - pointer to array of terms in expression group (e.g., (A*B*C(D,E) stored as [A,B,C,D,E])
/// \param k [in] index of term to be evaluated (e.g., 2 would correspond to C)
/// \param t [in] current model time
///  this is called recursively for nested terms in expression group
/// returns RAV_BLANK_DATA if entire term should be invalidated / not applied
//
double CDemandOptimizer::EvaluateTerm(expressionTerm **pTerms,const int k, const double &t) const
{
  int nshift;
  double duration;
  double x,y;
  int kk;
  int p;
  const expressionTerm *pT=pTerms[k];

  if      (pT->type == TERM_DV)
  {
    //do nothing - can't evaluate
    ExitGracefully("CDemandOptimizer::EvaluateTerm: trying to expand out active decision variable",RUNTIME_ERR);
    return -1;
  }
  else if (pT->type == TERM_TS)
  {
    nshift=pT->timeshift;
    return pT->pTS->GetValue(t +(double)(nshift));
  }
  else if (pT->type == TERM_CUMUL_TS)
  {
    duration=(double)(pT->timeshift);
    return pT->pTS->GetAvgValue(t-duration,duration)*(duration); //assumes daily timestep
  }
  else if (pT->type == TERM_LT)
  {
    kk=pT->nested_ind1;
    x=EvaluateTerm(pTerms,kk,t);
    return pT->pLT->GetValue(x); //decision variable can't be in lookup table
  }
  else if (pT->type == TERM_HRU)
  {
    int i=pT->SV_index;
    int k=pT->HRU_index;
    return _pModel->GetHydroUnit(k)->GetStateVarValue(i); //This will be start of timestep value
  }
  else if (pT->type == TERM_SB)
  {
    int i=pT->SV_index;
    int p=pT->p_index;
    return _pModel->GetSubBasin(p)->GetAvgStateVar(i);
  }
  else if (pT->type == TERM_CONST) 
  {
    return pT->value;
  }
  else if (pT->type == TERM_CONTROL)
  {
    int i=pT->DV_ind;
    return _pControlVars[i]->current_val;
  }
  else if (pT->type == TERM_CUMUL)  //!C123
  {
    int d=pT->p_index; //todo - should probably use something else here - misleading code
    return _aCumDelivery[d];
  }
  else if (pT->type == TERM_HISTORY) //e.g., !Q100[-3]
  {
    nshift=pT->timeshift-1;
    char   tmp =pT->origexp[1];//e.g., Q

    p=pT->p_index;

    if      (tmp=='Q'){return _aQhist[_aSBIndices[p]][nshift]; }
    else if (tmp=='h'){return _ahhist[_aSBIndices[p]][nshift]; }
    else if (tmp=='D'){return _aDhist[_aSBIndices[p]][nshift]; }
    else if (tmp=='I'){return _aIhist[_aSBIndices[p]][nshift]; }
    else if (tmp=='B'){return _pModel->GetSubBasin(p)->GetSpecifiedInflow(t+pT->timeshift*1.0); } //ASSUMES DAILY TIMESTEP!
    else if (tmp=='E'){return _pModel->GetSubBasin(p)->GetEnviroMinFlow  (t+pT->timeshift*1.0); } //ASSUMES DAILY TIMESTEP!
    else {
      ExitGracefully("CDemandOptimizer::EvaluateTerm: Invalid history variable ",BAD_DATA);
    }
  }
  else if (pT->type == TERM_MAX) {
    x=EvaluateTerm(pTerms,pT->nested_ind1,t);
    y=EvaluateTerm(pTerms,pT->nested_ind2,t);
    return max(x,y);
  }
  else if (pT->type == TERM_MIN) {
    x=EvaluateTerm(pTerms,pT->nested_ind1,t);
    y=EvaluateTerm(pTerms,pT->nested_ind2,t);
    return min(x,y);
  }
  else if (pT->type == TERM_CONVERT) {
    x=EvaluateTerm(pTerms,pT->nested_ind1,t);
    return x*pT->value;
  }
  return 0;
}
