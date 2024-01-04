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
/// constructor, destructor, and member functions of dv_constraint structure
//
dv_constraint::dv_constraint()
{
  name="";
  pExpression=NULL;

  is_goal=false;
  penalty_under=1.0;
  penalty_over=1.0;
  slack_ind1=DOESNT_EXIST;
  slack_ind2=DOESNT_EXIST;
  
  penalty_value=0;

  nConditions=0;
  pConditions=NULL;
}
dv_constraint::~dv_constraint()
{
  delete pExpression; 
  delete [] pConditions; nConditions=0;
}
void dv_constraint::AddCondition(dv_condition* pCond)
{
  if(!DynArrayAppend((void**&)(pConditions),(void*)(pCond),nConditions)) {
    ExitGracefully("dv_constraint::AddCondition: adding NULL condition",BAD_DATA_WARN);
  }
}
//////////////////////////////////////////////////////////////////
//Move to LookupTable.h/.cpp
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
string CLookupTable::GetName() const
{
  return _name;
}
double CLookupTable::GetValue(const double& x) const
{
  return InterpolateCurve(x,_aX,_aY,_nItems,false);
}
double CLookupTable::GetSlope(const double& x) const
{
  double dx=0.001*(_aX[_nItems-1]-_aX[0]);
  return (1.0/dx)*(InterpolateCurve(x+dx,_aX,_aY,_nItems,false)-InterpolateCurve(x,_aX,_aY,_nItems,false));
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
/// \brief parses individual expression string and converts to expression term structure
/// \param s [in] - string
/// \param term [out] - expression structure
/// \param lineno [in] - line number of original expression in input file filename, referenced in errors
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

  string tmp(s);
  term->origexp=tmp;

  size_t i1=s.find("[");
  size_t i2=s.find("]");
  size_t i3=s.find("(");
  
  if ((i1 != NPOS) && (i2 != NPOS) && (i3==NPOS)) { //array NOT nested within function
    string contents=s.substr(i1+1,i2-i1); cout<<" BRACKET CONTENTS: ["<<contents<<"]";
    n=s_to_i(contents.c_str());
    has_brackets=true;
  }

  //----------------------------------------------------------------------
  if      (s[0]=='!')  //decision variable e.g., !Q234 or !D42a 
  {
    int p=DOESNT_EXIST;
    int d=DOESNT_EXIST;
    
    if ((s[1] == 'Q') || (s[1] == 'h')) 
    {
      string sbid=s.substr(2,s.size() - 2); cout << "EXP BASIN ID: " << sbid<<endl;
      long SBID=s_to_l(sbid.c_str());
      p=_pModel->GetSubBasinIndex(SBID);
      if (p == DOESNT_EXIST) {
        warn="ConvertToExpressionTerm: invalid subbasin ID in expression"+warnstring;
        ExitGracefully(warn.c_str(),BAD_DATA_WARN);
      }
      if (!_pModel->GetSubBasin(p)->IsEnabled()) {
        warn="ConvertToExpressionTerm: reference to disabled subbasin "+warnstring+". Constraint/goal will be disabled";
        WriteWarning(warn.c_str(),true);
        return false;
      }
    } 
    else {
      string demandID;
      if (s[1] == '.') {demandID=s.substr(2);} //!D1234
      else             {demandID=s.substr(1);} //!D.FarmerBobsDemand
      if (d == DOESNT_EXIST) {
        warn = "ConvertToExpressionTerm: invalid demand ID or demand from disabled subbasin used in expression " +
               warnstring +" goal/constraint will be ignored";
        WriteWarning(warn.c_str(),true);
        return false; 
      }
    }
    
    if      (s[1]=='Q')
    {
      if (_aResIndices[p] != DOESNT_EXIST) {
        term->DV_ind = GetDVColumnInd(DV_QOUTRES, _aResIndices[p]);  
      } 
      else {
        term->DV_ind=GetDVColumnInd(DV_QOUT,p);  
      }
      //TODO - need to handle reservoir outflow 
    } 
    else if (s[1]=='D')
    {
      term->DV_ind=GetDVColumnInd(DV_DELIVERY,d);
    } 
    else if (s[1]=='h'){
      term->DV_ind=GetDVColumnInd(DV_STAGE,_aResIndices[p]);  
    } 
    term->type=TERM_DV;
  } 
  //----------------------------------------------------------------------
  else if (s.find("@ts(") != NPOS) //time series (e.g., @ts(my_time_series,n)
  {
    string name;
    int    index;
    size_t is = s.find("@ts(");
    size_t ie = s.find(",",is);
    size_t ip = s.find(")",ie);
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
      name  = s.substr(is+4,ie-is);
      index = s_to_i(s.substr(ie+1,ip-ie).c_str());
      term->pTS=NULL;
      for (int i = 0; i < _nUserTimeSeries; i++) {
        if (_pUserTimeSeries[i]->GetName() == name) {
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
  else if (s.find("@lookup(") != NPOS) //lookup table (e.g., @lookup(my_table,EXPRESSION)
  {
    string name;
    string x_in;
    size_t is = s.find("@lookup(");
    size_t ie = s.find(",",is);
    size_t ip = s.find(")",ie);
    if (ie == NPOS) {
      warn="ConvertToExpressionTerm: missing comma in @lookup expression"+warnstring;
      ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    }
    if (ip == NPOS) {
      warn="ConvertToExpressionTerm: missing end paretheses in @lookup expression"+warnstring;
      ExitGracefully(warn.c_str(), BAD_DATA_WARN);
    }
    if ((is != NPOS) && (ie != NPOS)) 
    {
      bool found=false;
      name = s.substr(is+8,ie-is);
      x_in = s.substr(ie+1,ip-ie);
      term->pLT=NULL;
      for (int i = 0; i < _nUserLookupTables; i++) {
        if (_pUserLookupTables[i]->GetName() == name) {
          term->pLT = _pUserLookupTables[i];
          found=true;
        }
      }  
      term->nested_exp1 =x_in;
      term->type        =TERM_LT;
      if (!found) {
        warn="ConvertToExpression: unrecognized lookup table name in @lookup command"+warnstring;
        ExitGracefully(warn.c_str(), BAD_DATA_WARN);
      }
    }
  }
  //----------------------------------------------------------------------
  else if (s.find("@max(") != NPOS) //max function (e.g., @max(exp1,exp2)
  {
    string name;
    string x_in,y_in;
    size_t is = s.find("@max(");
    size_t ie = s.find(",",is);
    size_t ip = s.find(")",ie);
    if (ie == NPOS) {
      warn="ConvertToExpressionTerm: missing comma in @max expression"+warnstring;
      ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    }
    if (ip == NPOS) {
      warn="ConvertToExpressionTerm: missing end paretheses in @max expression"+warnstring;
      ExitGracefully(warn.c_str(), BAD_DATA_WARN);
    }
    if ((is != NPOS) && (ie != NPOS)) 
    {
      x_in = s.substr(is+5,ie-is);
      y_in = s.substr(ie+1,ip-ie);
      term->nested_exp1 =x_in;
      term->nested_exp2 =y_in;
      term->type        =TERM_MAX;
    }
  }
  //----------------------------------------------------------------------
  else if (s.find("@min(") != NPOS) //max function (e.g., @min(exp1,exp2)
  {
    string name;
    string x_in,y_in;
    size_t is = s.find("@min(");
    size_t ie = s.find(",",is);
    size_t ip = s.find(")",ie);
    if (ie == NPOS) {
      warn="ConvertToExpressionTerm: missing comma in @min expression"+warnstring;
      ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    }
    if (ip == NPOS) {
      warn="ConvertToExpressionTerm: missing end paretheses in @min expression"+warnstring;
      ExitGracefully(warn.c_str(), BAD_DATA_WARN);
    }
    if ((is != NPOS) && (ie != NPOS)) 
    {
      x_in = s.substr(is+5,ie-is);
      y_in = s.substr(ie+1,ip-ie);
      term->nested_exp1 =x_in;
      term->nested_exp2 =y_in;
      term->type        =TERM_MIN;
    }
  }
  //----------------------------------------------------------------------
  else if (s.find("@convert(") != NPOS) //conversion (e.g., @convert(x,ACREFTPERDAY_TO_CMS)
  {
    string x_in;
    string units;
    size_t is = s.find("@convert(");
    size_t ie = s.find(",",is);
    size_t ip = s.find(")",ie);
    if (ie == NPOS) {
      warn="ConvertToExpressionTerm: missing comma in @convert expression"+warnstring;
      ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    }
    if (ip == NPOS) {
      warn="ConvertToExpressionTerm: missing end parentheses in @convert expression"+warnstring;
      ExitGracefully(warn.c_str(), BAD_DATA_WARN);
    }
    if ((is != NPOS) && (ie != NPOS)) 
    {
      bool found=false;
      x_in  = s.substr(is+8,ie-is);
      units = s.substr(ie+1,ip-ie);
 
      if (units == "ACREFTD_TO_CMS") {term->value=1.0/ACREFTD_PER_CMS;     found=true;}
      if (units == "CMS_TO_ACREFTD") {term->value=ACREFTD_PER_CMS;     found=true;}
      if (units == "INCHES_TO_MM"  ) {term->value=MM_PER_INCH;           found=true;}
      if (units == "MM_TO_INCHES"  ) {term->value=1.0/MM_PER_INCH;       found=true;}
      term->nested_exp1 =x_in;
      term->type     =TERM_CONVERT;
      if (!found) {
        warn="ConvertToExpression: unrecognized time series name in @convert command"+warnstring;
        ExitGracefully(warn.c_str(), BAD_DATA_WARN);
      }
    }
  }
  //----------------------------------------------------------------------
  else if (has_brackets) // historical lookup (e.g., Q102[-2])
  {
    ExitGracefullyIf(n>=0,"ConvertToExpressionTerm:: term in square brackets must be <=0 for historical variables.",BAD_DATA_WARN);
    if (n<0)
    {
      term->type=TERM_HISTORY;
      term->timeshift=-n;

      string sbid=s.substr(1,i1-1); cout<<"EXP[] BASIN ID: "<<sbid<<endl;
      long   SBID=s_to_l(sbid.c_str());
      int       p=_pModel->GetSubBasinIndex(SBID);
      if (p == DOESNT_EXIST) {
        ExitGracefully("ConvertToExpressionTerm: invalid subbasin ID in expression",BAD_DATA_WARN);
      }
      if (!_pModel->GetSubBasin(p)->IsEnabled()) {
        WriteWarning("ConvertToExpressionTerm: expression includes reference to disabled subbasin. Constraint/goal will be disabled",true);
        //return false;
      }
      //TODO check if subbasin is enabled
    }
  }
  //----------------------------------------------------------------------
  else if (GetNamedConstant(s)!=RAV_BLANK_DATA) // named constant 
  {
    term->type=TERM_CONST;
    term->value=GetNamedConstant(s);
  }
  //----------------------------------------------------------------------
  else if (is_numeric(s.c_str()))  // numeric value (e.g., 132.4)
  {
    term->type=TERM_CONST;
    term->value=s_to_d(s.c_str());
  }
  //----------------------------------------------------------------------
  else {
    warn="ConvertToExpression:: Unparseable term in expression"+warnstring;
    ExitGracefully(warn.c_str(), BAD_DATA_WARN);
  }
  return true;
}

//////////////////////////////////////////////////////////////////
/// \brief Parses demand constraint expression of form (e.g.,) A * B(x) * C + D * E - F = G * H(t) [NO PARENTHESES!]
/// \params s [in] - array of strings of [size: Len]
/// \param Len [in] - length of string array
/// \param lineno [in] - line number of original expression in input file filename, referenced in errors 
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

  //identify all strings as either operators or terms
  //determine nature of comparison =,<,>
  for (int i = 1; i < Len; i++) //starts at 1 - first term is :Expression or :DefineDecisionVariable
  {
    type[i]=EXP;
    if ((s[i][0]=='+') || (s[i][0]=='*') || (s[i][0]=='=') || (s[i][0]=='/') || (s[i][0]=='<') || (s[i][0]=='>')){
      type[i] = EXP_OP; 
    } 
    if (s[i-1][0]=='/'){type[i]=EXP_INV; }

    if ((s[i][0]=='=') || (s[i][0]=='~')){
      comp=COMPARE_IS_EQUAL;    rhs_ind=i; 
    }
    if (s[i][0]=='<'){//works for <=, too 
      comp=COMPARE_LESSTHAN;    rhs_ind=i; 
    }
    if (s[i][0]=='>'){ //works for >=, too
      comp=COMPARE_GREATERTHAN; rhs_ind=i; 
    }
  }
  if (comp == COMPARE_BETWEEN) {
    ExitGracefully("ParseExpression: missing comparison operator (e.g., <, >, =, in expression",BAD_DATA_WARN);
  }

  int  j=0;
  int  k=0;
  int  DVcount=0;
  bool valid;
  for (int i = 1; i < Len; i++) 
  {
    if (type[i] != EXP_OP) 
    {
      terms[j][k] = new expressionTerm();
      valid=ConvertToExpressionTerm(to_string(s[i]),terms[j][k],lineno,filename);
      if (!valid){return NULL;}

      if (i>rhs_ind     ){terms[j][k]->mult*=-1; }
      if (s[i-1][0]=='-'){terms[j][k]->mult*=-1; }
      if (terms[j][k]->type == TERM_DV) {
        DVcount++;
        if (DVcount > 1) {
          ExitGracefully("ParseExpression: more than one decision variable found in product term within constraint expression. Cannot handle non-linear expressions.",BAD_DATA_WARN);
        }
      }
      if (type[i] == EXP_INV) {
        terms[j][k]->reciprocal=true;
      }
      //Handle nested expressions 
      if (terms[j][k]->nested_exp1 != "") 
      {
        terms[j][k] = new expressionTerm();
        valid=ConvertToExpressionTerm(terms[j][k]->nested_exp1, terms[j][k + 1],lineno, filename);
        if (!valid){return NULL;}
        terms[j][k]->nested_ind1=k+1;
        terms[j][k+1]->is_nested=true;
        //TODO - THROW WARNING IF DV IS IN NESTED EXP
        k++;
      }
      if (terms[j][k]->nested_exp2 != "") 
      {
        terms[j][k] = new expressionTerm();
        valid=ConvertToExpressionTerm(terms[j][k]->nested_exp2, terms[j][k + 1],lineno, filename);
        if (!valid){return NULL;}
        terms[j][k]->nested_ind2=k+1;
        terms[j][k+1]->is_nested=true;
        //TODO - THROW WARNING IF DV IS IN NESTED EXP
        k++;
      }
      k++;
    }
    else if (i!=1){ //operator - expression product term finished and expression doesn't start with + or - 
      termspergrp[j]=k;
      j++; k=0; DVcount=0;
    }
  }

  //copy to new structure 
  int nGroups=j;
  expressionStruct *tmp;
  tmp=new expressionStruct();
  tmp->nGroups=nGroups;
  tmp->compare=comp;

  tmp->pTerms=new expressionTerm **[nGroups];
  tmp->nTermsPerGrp=new int[nGroups];
  for (int j = 0; j < nGroups; j++) {
    tmp->pTerms[j] = new expressionTerm * [nGroups];
    tmp->nTermsPerGrp[j]=termspergrp[k];
    for (int k = 0; k < termspergrp[j]; k++) {
      tmp->pTerms[j][k]=terms[j][k];
    }
  }
  return tmp;
}

void CDemandOptimizer::CheckConditions(const int ii, const time_struct &tt)
{ 
  double dv_value;
  int p;
  dv_constraint *pC=_pConstraints[ii];
  _pConstraints[ii]->conditions_satisfied=false;
  //Check if conditionals are satisfied 
  for (int j = 0; j < pC->nConditions; j++) 
  {
    dv_value=0.0;
    if (pC->pConditions[j]->dv_name == "DAY_OF_YEAR") {
      dv_value=tt.julian_day;
    }
    else if (pC->pConditions[j]->dv_name == "MONTH") {
      dv_value=(double)(tt.month);
    }
    else if (pC->pConditions[j]->dv_name[0] == '!') { //decision variable 
      char   tmp =pC->pConditions[j]->dv_name[1];
      string tmp2=pC->pConditions[j]->dv_name.substr(2);
      long ind=s_to_l(tmp2.c_str());

      if      (tmp == 'Q') {
        dv_value=_pModel->GetSubBasinByID(ind)->GetOutflowRate();
      }
      else if (tmp == 'h') {
        dv_value=_pModel->GetSubBasinByID(ind)->GetReservoir()->GetResStage();
      }
      else if (tmp == 'C') {
        p=_pModel->GetSubBasinIndex(ind);
        dv_value=_aCumDelivery[p];
      }
      else if (tmp == 'D') {
        //int d=GetDemandIndex(tmp2);
        //int p=_aDemandBasin[d];
        //dv_value = _pModel->GetSubBasinByID(p)->GetDemandDelivery();
        //TODO 
      }
      //TODO : handle user specified DVs
    }

    comparison comp=pC->pConditions[j]->compare;
    if (comp == COMPARE_BETWEEN) {
      if ((dv_value > pC->pConditions[j]->value2) || (dv_value< pC->pConditions[j]->value)){return;}//don't apply condition if ANY conditionals unmet
    }
    else if (comp == COMPARE_GREATERTHAN) {
      if (dv_value< pC->pConditions[j]->value){return;}//don't apply condition
    }
    else if (comp == COMPARE_LESSTHAN) {
      if (dv_value> pC->pConditions[j]->value){return;}//don't apply condition
    }
  }
  _pConstraints[ii]->conditions_satisfied=true;
}
    //////////////////////////////////////////////////////////////////
/// adds constraint ii to LP solve problem statement 
/// \params ii [in] - index of constraint in _pConstraints array
/// \param pLinProg [out] - pointer to valid lpsolve structure to be modified
/// \param tt [in] - time structure  
/// \param *col_ind [in] - empty array (with memory reserved) for storing column indices 
/// \param *row_val [in] - empty array (with memory reserved) for storing row values
//
#ifdef _LPSOLVE_
void CDemandOptimizer::AddConstraintToLP(const int ii, lp_lib::lprec* pLinProg, const time_struct &tt, int *col_ind, double *row_val) const 
{
  dv_constraint *pC=_pConstraints[ii];

  double coeff;
  int    i=0;
  int    retval;
  double RHS;
  bool   group_has_dv;
  int    DV_ind;
  
  int    nn=(int)((tt.model_time+TIME_CORRECTION)/1.0);//Options.timestep; //TMP DEBUG
 
  expressionStruct *pE= pC->pExpression;
  
  RHS=0.0;
  for (int j = 0; j < pE->nGroups; j++)
  {
    coeff=1.0;
    group_has_dv=false;
    DV_ind=DOESNT_EXIST;
    for (int k = 0; k < pE->nTermsPerGrp[j]; k++)
    {
      if (pE->pTerms[j][k]->type == TERM_DV) {
        DV_ind=pE->pTerms[j][k]->DV_ind;
        group_has_dv=true;
      }
      else if (!(pE->pTerms[j][k]->is_nested))
      {
        if (pE->pTerms[j][k]->reciprocal == true) {
          //TODO: SHOULD CHECK FOR DIVIDE BY ZERO, BUT DO WHAT?
          coeff /=  EvaluateTerm(pE->pTerms[j], k, nn)/(pE->pTerms[j][k]->mult);
        } 
        else {
          coeff *= (pE->pTerms[j][k]->mult) * EvaluateTerm(pE->pTerms[j], k, nn);
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
  int constr_type;
  int nSlack=1;
  coeff=1.0;
  if      (pE->compare==COMPARE_IS_EQUAL)   {constr_type=ROWTYPE_EQ; coeff=-1.0; nSlack=2;}
  else if (pE->compare==COMPARE_LESSTHAN)   {constr_type=ROWTYPE_LE; coeff=-1.0;}
  else if (pE->compare==COMPARE_GREATERTHAN){constr_type=ROWTYPE_GE; coeff=+1.0;}

  //Add slack variables to goals
  if (pC->is_goal) {
    col_ind[i]=pC->slack_ind1;
    row_val[i]=coeff;
    i++;
    if (nSlack == 2) {
      col_ind[i]=pC->slack_ind2;
      row_val[i]=-coeff;
      i++;
    }
  }

  retval = lp_lib::add_constraintex(pLinProg,i,row_val,col_ind,constr_type,RHS);
  ExitGracefullyIf(retval==0,"AddConstraintToLP::Error adding user-specified constraint/goal",RUNTIME_ERR);
}
#endif 
//////////////////////////////////////////////////////////////////
/// evaluates term in constraint/goal expression to numerical value (unless term is DV) 
/// \params pTerms [in] - pointer to array of terms in expression group (e.g., (A*B*C(D,E) stored as [A,B,C,D,E])
/// \param k [in] index of term to be evaluated (e.g., 2 would correspond to C)
/// \param nn [in] current timestep index
///  this is called recursively for nested terms in expression group
//
double CDemandOptimizer::EvaluateTerm(expressionTerm **pTerms,const int k, const int nn) const
{
  int nshift;
  double x,y;
  int kk;
  int p;
  const expressionTerm *pT=pTerms[k];
  
  if (pT->type == TERM_DV) 
  {
    //do nothing - can't evaluate
    return -1;
  }
  else if (pT->type == TERM_TS)
  {
    nshift=pT->timeshift;
    return pT->pTS->GetValue(nn+nshift);
  }
  else if (pT->type == TERM_LT) 
  {
    kk=pT->nested_ind1;
    x=EvaluateTerm(pTerms,kk,nn);
    return pT->pLT->GetValue(x); //decision variable can't be in lookup table 
  }
  else if (pT->type == TERM_CONST) 
  {
    return pT->value;
  }
  else if (pT->type == TERM_HISTORY) //e.g., !Q100[-3]
  { 
    nshift=pT->timeshift;
    char   tmp =pT->origexp[1];
    string tmp2=pT->origexp.substr(2);
    p=_pModel->GetSubBasinIndex(s_to_l(tmp2.c_str())); 
    if (tmp=='Q'){return _aQhist[p][nshift]; }
    if (tmp=='h'){return _ahhist[p][nshift]; }
    if (tmp=='D'){return _aDhist[p][nshift]; }
  }
  else if (pT->type == TERM_MAX) {
    x=EvaluateTerm(pTerms,pT->nested_ind1,nn);
    y=EvaluateTerm(pTerms,pT->nested_ind2,nn);
    return max(x,y);
  }
  else if (pT->type == TERM_MIN) {
    x=EvaluateTerm(pTerms,pT->nested_ind1,nn);
    y=EvaluateTerm(pTerms,pT->nested_ind2,nn);
    return min(x,y);
  } 
  else if (pT->type == TERM_CONVERT) {
    x=EvaluateTerm(pTerms,pT->nested_ind1,nn);
    return x*pT->value;
  } 
  return 0;
}