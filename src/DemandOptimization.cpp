/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2024 the Raven Development Team

  Uses lpsolve:
    https://lpsolve.sourceforge.net/5.5/
    Files from lp_solve_5.5.2.11_dev_win64.zip distribution
    under GNU Public License
  ----------------------------------------------------------------*/
#include "DemandOptimization.h"

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the Demand optimization constructor
/// \details Creates empty instance of the demand optimization class 
//
CDemandOptimizer::CDemandOptimizer(CModel *pMod)
{
  _pModel=pMod;

  _nDemands=0;
  _aDemandPenalties=NULL;
  _aDemandIDs=NULL;
  _aDelivery=NULL;
  _nReservoirs=0;
  _nEnabledSubBasins=0;

  _aResIndices=NULL;

  _nDecisionVars=0;
  _pDecisionVars=NULL;
  _nUserDecisionVars=0;

  _nSlackVars=0;
  _aSlackValues=NULL;

  _nUserConstants=0;
  _aUserConstNames=NULL;
  _aUserConstants=NULL;
  
  _nUserTimeSeries=0;
  _pUserTimeSeries=NULL;

  _nUserLookupTables=0;
  _pUserLookupTables=NULL;

  _nHistoryItems=0;
  _aQhist=NULL;
  _aDhist=NULL;
  _ahhist=NULL;

  _aCumDelDate=NULL; 
  _aCumDelivery=NULL;

  _nConstraints=0;
  _pConstraints=NULL;

}
//////////////////////////////////////////////////////////////////
/// \brief Implementation of the Demand optimization destructor
//
CDemandOptimizer::~CDemandOptimizer() 
{
  for (int i=0;i<_nDecisionVars;    i++){delete _pDecisionVars[i];    }delete [] _pDecisionVars; 
  for (int i=0;i<_nUserTimeSeries;  i++){delete _pUserTimeSeries[i];  }delete [] _pUserTimeSeries; 
  for (int i=0;i<_nUserLookupTables;i++){delete _pUserLookupTables[i];}delete [] _pUserLookupTables;
  for (int p = 0; p < _pModel->GetNumSubBasins();p++){ //TMP DEBUG - this should work with enabled subbasins only ?
    if (_aQhist!=NULL){
      delete [] _aQhist[p];
      delete [] _ahhist[p];
      delete [] _aDhist[p];
    }
  }
  delete [] _aQhist;
  delete [] _ahhist;
  delete [] _aDhist;
  delete [] _aUserConstants;
  delete [] _aUserConstNames;
  delete [] _aDelivery;
  delete [] _aDemandIDs;
  delete [] _aSlackValues;
  delete [] _aResIndices;
}
//////////////////////////////////////////////////////////////////
/// \brief gets demand index d from name
/// \params demand_tag [in] - demand ID value (e.g., '123') or Alias (e.g., 'FarmerBob' ) -WITHOUT !D or !D. start
/// \returns demand index d or DOESNT_EXIST if invalid dname 
/// must be called after _aDemandIDs is populated in CDemandOptimizer::Initialize()
//
int CDemandOptimizer::GetDemandIndexFromName(const string demand_tag) const 
{  
  for (int d = 0; d < _nDemands; d++) {
    if (s_to_i(demand_tag.c_str()) == _aDemandIDs[d]) {return d;}
    if (demand_tag == _aDemandAliases[d]) {return d;}
  } 
  return DOESNT_EXIST;
}
//////////////////////////////////////////////////////////////////
/// \brief sets history length 
/// \params n [in] - number of timesteps in history
//
void CDemandOptimizer::SetHistoryLength(const int n)
{
  _nHistoryItems=n;
}
//////////////////////////////////////////////////////////////////
/// \brief date at which cumulative delivery is rebooted to zero (usually winter date) 
/// \params julian_date [in] reboot date
/// \params demandID [in] demand identifier (e.g., 'D123')
//
void CDemandOptimizer::SetCumulativeDate(const int julian_date, const string demandID) 
{
  int d;
  d=GetDemandIndexFromName(demandID);
  if (d!=DOESNT_EXIST){
    _aCumDelDate[d]=julian_date;
  }
  else{
    WriteWarning("DemandOptimization::SetCumulativeDate: invalid demand ID",true);
  }
}
    //////////////////////////////////////////////////////////////////
/// \brief adds decision variable structure to _pDecisionVars member aarray 
/// \params pDV [in] - user-specified decision variable to be added 
//
void CDemandOptimizer::AddDecisionVar(const decision_var* pDV) 
{
  if (!DynArrayAppend((void**&)(_pDecisionVars),(void*)(pDV),_nDecisionVars)){
   ExitGracefully("CDemandOptimizer::AddDecisionVar: adding NULL DV",BAD_DATA);}
}
//////////////////////////////////////////////////////////////////
/// \brief sets bounds for user specified decision variable 
//
void CDemandOptimizer::SetDecisionVarBounds(const string name, const double& min, const double& max) 
{
  for (int i = 0; i < _nDecisionVars; i++) {
    if (_pDecisionVars[i]->name==name){
      _pDecisionVars[i]->min=min;
      _pDecisionVars[i]->max=max;
      return;
    }
  }
  string warn = "SetDecisionVarBounds: invalid decision variable name (" + name + "): must define before use";
  WriteWarning(warn.c_str(),_pModel->GetOptStruct()->noisy);
}
//////////////////////////////////////////////////////////////////
/// \brief adds user constants 
//
void CDemandOptimizer::AddUserConstant(const string name, const double& val) 
{
  string *tmp  = new string[_nUserConstants+1];
  double *tmp2 = new double[_nUserConstants+1];
  for (int i=0;i<_nUserConstants;i++){
    tmp [i] = _aUserConstNames[i];
    tmp2[i] = _aUserConstants[i];
  }
  tmp [_nUserConstants]=name;
  tmp2[_nUserConstants]=val;
  _nUserConstants++;
  delete[] _aUserConstNames;
  delete[] _aUserConstants;
  _aUserConstNames=tmp;
  _aUserConstants=tmp2;
}
//////////////////////////////////////////////////////////////////
/// \brief adds user time series 

void CDemandOptimizer::AddUserTimeSeries(const CTimeSeries *pTS)
{
  if (!DynArrayAppend((void**&)(_pUserTimeSeries),(void*)(pTS),_nUserTimeSeries)){
   ExitGracefully("CDemandOptimizer::AddDecisionVar: adding NULL time series",BAD_DATA);}
}
//////////////////////////////////////////////////////////////////
/// \brief adds user lookup table 
//
void CDemandOptimizer::AddUserLookupTable(const CLookupTable *pLUT)
{
  if (!DynArrayAppend((void**&)(_pUserLookupTables),(void*)(pLUT),_nUserLookupTables)){
   ExitGracefully("CDemandOptimizer::AddDecisionVar: adding NULL time series",BAD_DATA);}
}
  // 
//////////////////////////////////////////////////////////////////
/// \brief adds decision variable constraint OR goal to _pConstraints member array 
/// \params name [in] - constraint name
/// \params exp [in] - expression structure
/// \param soft_constraint [in] - true if this is a goal, not a constraint
/// 
dv_constraint *CDemandOptimizer::AddConstraint(string name, expressionStruct *exp, bool soft_constraint) 
{
  dv_constraint *pC=new dv_constraint();
  pC->name=name;
  pC->is_goal=soft_constraint;

  if (!DynArrayAppend((void**&)(_pConstraints),(void*)(pC),_nConstraints)){
   ExitGracefully("CDemandOptimizer::AddConstraint: adding NULL constraint",BAD_DATA);} 

  return _pConstraints[_nConstraints-1];
}
//////////////////////////////////////////////////////////////////
/// \brief Initializes Demand optimization instance
/// \notes to be called after .rvh file read and subbasin network initialization, but prior to reading .rvm file  
/// \params pModel [in] - pointer to model 
/// \params Options  [in] - model options structure 
//
void CDemandOptimizer::Initialize(CModel* pModel, const optStruct& Options) 
{
#ifndef _LPSOLVE_
  ExitGracefully("Demand optimization requires compilation with lpsolve library.",RUNTIME_ERR);
#endif
  int           p;
  string        name;
  CSubBasin    *pSB;
  decision_var *pDV;
  
  _nEnabledSubBasins=0;
  for (int pp=0;pp<pModel->GetNumSubBasins();pp++)
  {
    p=pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    if (pSB->IsEnabled()){
      if (pSB->GetReservoir()!=NULL){
        name="L"+to_string(pSB->GetID()); //special case - reach outflow to reservoir is internally used only
      }
      else{
        name="Q"+to_string(pSB->GetID());
      }
      pDV=new decision_var(name,p,DV_QOUT);

      AddDecisionVar(pDV);
    }
    _nEnabledSubBasins++;
  }

  _nReservoirs=0;
  for (int pp=0;pp<pModel->GetNumSubBasins();pp++)
  {
    p=pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    if (pSB->IsEnabled()){
      if (pSB->GetReservoir()!=NULL){
        _nReservoirs++;
        name="Q"+to_string(pSB->GetID());
        pDV=new decision_var(name,p,DV_QOUTRES);
        AddDecisionVar(pDV);
      }
    }
  }
  for (int pp=0;pp<pModel->GetNumSubBasins();pp++)
  {
    p=pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    if (pSB->IsEnabled()){
      if (pModel->GetSubBasin(p)->GetReservoir()!=NULL){
        name="h"+to_string(pModel->GetSubBasin(p)->GetID());
        pDV=new decision_var(name,p,DV_STAGE);
        AddDecisionVar(pDV);
      }
    }
  }
  int rcount=0;
  _aResIndices = new int[_nReservoirs];
  for (int pp=0;pp<pModel->GetNumSubBasins();pp++)
  {
    p=pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    _aResIndices[p]=DOESNT_EXIST;
    if (pSB->IsEnabled()){
      if (pModel->GetSubBasin(p)->GetReservoir()!=NULL){
        _aResIndices[p]=rcount;
        rcount++;
      }
    }
  }

  _nDemands=0;
  for (p=0;p<pModel->GetNumSubBasins();p++)
  {
    pSB=pModel->GetSubBasin(p);
    if (pSB->IsEnabled()){
      for (int i = 0; i < pModel->GetSubBasin(p)->GetNumWaterDemands();i++) {_nDemands++;}
    }
  }
  /* for (p = 0; p < pModel->GetNumSubBasins(); p++)
  {
    pSB=pModel->GetSubBasin(p);
    if (pSB->IsEnabled()){
      if (pModel->GetSubBasin(p)->GetReservoir()!=NULL){
        for (int i = 0; i < pModel->GetSubBasin(p)->GetReservoir()->GetNumWaterDemands();i++) {_nDemands++;}
      }
    }
  }*/

  _aDemandIDs     = new int   [_nDemands];
  _aDemandAliases = new string[_nDemands];
  int d=0;
  int ID;
  for (int pp=0;pp<pModel->GetNumSubBasins();pp++)
  {
    p=pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    if (pSB->IsEnabled()){
      for (int i = 0; i < pSB->GetNumWaterDemands();i++) 
      {
        ID = pSB->GetWaterDemandID(i);//e.g., 132
        pDV=new decision_var(name,p,DV_DELIVERY);
        pDV->dem_index=i;
        AddDecisionVar(pDV);

        _aDemandIDs[d]=ID;

        name = pSB->GetWaterDemandName(i);
        _aDemandAliases[d]=name;
        d++;
      }
    }
  }
  //TODO: Add code to get demand IDs for reservoir extractions - no slack vars for this 

  //Add slack variables for environmental flow goals 
  //----------------------------------------------------------
   int slack_ind=0;
  _nSlackVars=0;
  for (int pp=0;pp<pModel->GetNumSubBasins();pp++)
  {
    p=pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    if (pSB->IsEnabled()){
      for (int i = 0; i < pSB->GetNumWaterDemands();i++) 
      {  
        //created for all demand nodes even if environmental min flow is zero, for now
        if (pSB->GetNumWaterDemands() > 0) {
          pDV = new decision_var(name,p,DV_SLACK); 
          pDV->dem_index=slack_ind;
          AddDecisionVar(pDV);
          slack_ind++;
          _nSlackVars++;
        }
      }
    }
  }

  _aDemandPenalties = new double[_nDemands];
  _aCumDelivery     = new double[_nDemands];
  _aCumDelDate      = new int   [_nDemands];
  for (int i=0;i<_nDemands;i++){
    _aDemandPenalties[i]=0.0;
    _aCumDelivery    [i]=0.0;
    _aCumDelDate     [i]=0; //Jan 1 
  }

  //print summary
  if (true){
    cout<<"MANAGEMENT OPTIMIZATION INITIALIZE"<<endl;
    cout<<" -----------------------------------------------------------------"<<endl;
    cout<<" # Enabled subbasins: "<<_nEnabledSubBasins<<" of "<<pModel->GetNumSubBasins()<<endl;
    cout<<" # Decision Vars: "<<_nDecisionVars<<endl;
    for (int i = 0; i < _nDecisionVars; i++) {
      cout<<"    "<<i<<": "<<_pDecisionVars[i]->name<<" "<<_pDecisionVars[i]->dvar_type<<endl;
    }
    cout<<" # Demands : "<<_nDemands<<endl;
    for (int d = 0; d < _nDemands; d++) {
      cout << "    " << d << ": " << _aDemandIDs[d] << " in basin ?? "<<endl; //add alias, SBID?
    }
  }

}
void CDemandOptimizer::InitializePostRVMRead(CModel* pModel, const optStruct& Options) 
{
  //This must be done after .rvm file read \todo[move] - requires history size 
  int nSBs = _pModel->GetNumSubBasins();
  
  if (_nHistoryItems==0){_nHistoryItems=static_cast<int>(rvn_floor(1.0/Options.timestep+REAL_SMALL));} //default = 1 day

  _aQhist = new double *[nSBs];
  _ahhist = new double *[nSBs];
  _aDhist = new double *[nSBs];
  for (int p=0;p<nSBs;p++){
    _aQhist[p] = new double[_nHistoryItems];
    _ahhist[p] = new double[_nHistoryItems];
    _aDhist[p] = new double[_nHistoryItems];
    for (int i = 0; i < _nHistoryItems; i++) {
      _aQhist[p][i]=0.0;
      _ahhist[p][i]=0.0;
      _aDhist[p][i]=0.0;
    }
  }

  // create slack variables for environmental min flow constraints 
  int p;
  int slack_ind=0;
  _nSlackVars=0;
  decision_var *pDV;
  CSubBasin *pSB;
  for (int pp=0;pp<pModel->GetNumSubBasins();pp++)
  {
    p=pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    if (pSB->IsEnabled()){
      //created for all demand nodes even if environmental min flow is zero, for now
      if (pModel->GetSubBasin(p)->GetNumWaterDemands() > 0) 
      { 
        pDV =  new decision_var("E" + to_string(pModel->GetSubBasin(p)->GetID()), p, DV_SLACK); 
        pDV->dem_index=slack_ind;
        AddDecisionVar(pDV);
        slack_ind++;
        _nSlackVars++;
      }
    }
  }

  //sift through goals, update _nSlackVars, reserve memory for _aSlackValues
  // ISSUE : THIS IS TOTAL POSSIBLE NUMBER OF SLACK VARIABLES; MAY NOT BE IN PLAY - for disabled constraints, just set penalty to zero>??
  for (int j = 0; j < _nConstraints; j++) 
  {
    if (_pConstraints[j]->is_goal) {
      if (_pConstraints[j]->pExpression->compare == COMPARE_IS_EQUAL) {
        _pConstraints[j]->slack_ind1=_nSlackVars;
        _pConstraints[j]->slack_ind2=_nSlackVars+1;
        _nSlackVars+=2;
      } 
      else {
        _pConstraints[j]->slack_ind1=_nSlackVars;
        _nSlackVars++;
      }
    }
  }
  _aSlackValues = new double[_nSlackVars];

  _nUserDecisionVars=0;
  for (int ii = 0; ii < _nDecisionVars; ii++) {
    if (_pDecisionVars[ii]->dvar_type == DV_USER) {_nUserDecisionVars++;}
  }
}
//////////////////////////////////////////////////////////////////
// indexing for DV vector
// 0     ..   nSB-1 - outflow of stream reach p 
// nSB   ..   nSB+nRes-1 - reservoir outflows
// nSB+nRes.. nSB+2*Res-1 - reservoir stages 
// nSB+2*nRes .. nSB+2*nRes+nDem-1 - satisfied demand 
// nSB+2*nRes+nDem .. nSB+2*nRes+nDem+nDV - user-specified decision variables 
// remaining - slack vars - [#enviro min flow+ # user > < goals + 2* user == goals]
//  
int CDemandOptimizer::GetDVColumnInd(const dv_type typ, const int counter) const 
{
  if      (typ==DV_QOUT    ){return counter+1;} 
  else if (typ==DV_QOUTRES ){return _nEnabledSubBasins+counter+1;}
  else if (typ==DV_STAGE   ){return _nEnabledSubBasins+_nReservoirs+counter+1;}
  else if (typ==DV_DELIVERY){return _nEnabledSubBasins+2*_nReservoirs+counter+1;} 
  else if (typ==DV_USER    ){return _nEnabledSubBasins+2*_nReservoirs+_nDemands+counter+1; }
  else if (typ==DV_SLACK   ){return _nEnabledSubBasins+2*_nReservoirs+_nDemands+_nUserDecisionVars+counter+1;}
  
  return 0;
}
//////////////////////////////////////////////////////////////////
/// \brief Solves demand optimization problem 
/// \notes to be called every time step in lieu of routing mass balance 
/// \params pModel [in] - pointer to model 
/// \params Options  [in] - model options structure 
/// \params aSBrunoff [in] - array of current amount of runoff released to each subbasin (m3)  [size:nSubBasins]
/// \params t [in] - local model time 
//
void CDemandOptimizer::SolveDemandProblem(CModel *pModel, const optStruct &Options, const double *aSBrunoff,const time_struct &tt ) 
{
#ifdef _LPSOLVE_

  int    d,i,p;
  double RHS;
  int    retval;
  double demand_penalty_sum=0; //sum (P * Dstar) to be added to global obj function such that min is zero 
  double tstep=Options.timestep;

  CSubBasin *pSB;

  double t=tt.model_time;
  //int nn=tt.nn+1;//end of timestep
  int nn=static_cast<int>((t+TIME_CORRECTION)/Options.timestep)+1;//end-of timestep index
  
  int    *col_ind=new int    [_nDecisionVars]; //index of column to insert value in current row (1:nDV, not zero-indexed)
  double *row_val=new double [_nDecisionVars]; //values of row[col_ind] 

  // instantiate linear programming solver
  // ----------------------------------------------------------------
  lp_lib::lprec *pLinProg; 

  pLinProg=lp_lib::make_lp(0,_nDecisionVars);
  if (pLinProg==NULL){ExitGracefully("Error in SolveDemandProblem(): couldn construct new linear programming model",RUNTIME_ERR);}
  char tmp[200];
  
  strcpy(tmp,"RavenLPSolve");
  lp_lib::set_lp_name(pLinProg,tmp);

  for (int i=0;i<_nDecisionVars;i++)
  {
    char name[200];
    strcpy(name,_pDecisionVars[i]->name.c_str());
    lp_lib::set_col_name(pLinProg, i + 1, name);
  }
  
  //Set upper bounds of delivery (D) to demand (D*) (preferred to adding constraint)
  // ----------------------------------------------------------------
  d=0;
  for (int pp = 0; pp<pModel->GetNumSubBasins(); pp++)
  {
    p  =pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    if (pSB->IsEnabled()){
      for (int ii = 0; ii < pSB->GetNumWaterDemands();ii++)
      {
        retval=lp_lib::set_upbo(pLinProg,GetDVColumnInd(DV_DELIVERY,d), pSB->GetIrrigationDemand(ii));
        ExitGracefullyIf(retval!=0,"SolveDemandProblem::Error adding demand upper bound",RUNTIME_ERR);
        d++;
      }
    }
  }
  
  //Set bounds of user-specified decision variables 
  // ----------------------------------------------------------------
  for (int i = 0; i < _nDecisionVars; i++) {
    
    if (_pDecisionVars[i]->max < ALMOST_INF * 0.99) {
        retval=lp_lib::set_upbo(pLinProg,GetDVColumnInd(DV_USER,i), _pDecisionVars[i]->max);
        ExitGracefullyIf(retval!=0,"SolveDemandProblem::Error adding decision variable upper bound",RUNTIME_ERR);
    }
    if (_pDecisionVars[i]->min != 0.0) {
        retval=lp_lib::set_lowbo(pLinProg,GetDVColumnInd(DV_USER,i), _pDecisionVars[i]->min);
        ExitGracefullyIf(retval!=0,"SolveDemandProblem::Error adding decision variable lower bound",RUNTIME_ERR);
    }
  }


  lp_lib::set_add_rowmode(pLinProg, TRUE); //readies lp_lib to add rows and objective function
  

  // ----------------------------------------------------------------
  // create objective function : minimize(sum(penalty*undelivered demand) + sum(penalty*violations))
  // ----------------------------------------------------------------
  i=0;

  //add blank penalty to objective function for models without demand -> 0*Q_0-->min
  // ----------------------------------------------------------------
  col_ind[i]=GetDVColumnInd(DV_DELIVERY,0);
  row_val[i]=0.0;
  i++;

  // penalties for unmet demand
  // ----------------------------------------------------------------
  for (int pp = 0; pp<pModel->GetNumSubBasins(); pp++)
  {
    p  =pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    if (pSB->IsEnabled()){
      for (int ii = 0; ii < pModel->GetSubBasin(p)->GetNumWaterDemands();ii++)
      {
        demand_penalty_sum+=(_aDemandPenalties[i]*pSB->GetIrrigationDemand(ii));

        col_ind[i]=GetDVColumnInd(DV_DELIVERY,i);
        row_val[i]=_aDemandPenalties[i];
        i++;
      }
    }
  }
  // user-specified goal penalties 
  // ----------------------------------------------------------------
  for (int j = 0; j < _nConstraints; j++) 
  {
    CheckConditions(j,tt);

    if ((_pConstraints[j]->is_goal) && (_pConstraints[j]->conditions_satisfied))
    {
      if (_pConstraints[j]->pExpression->compare == COMPARE_IS_EQUAL) //two slack variables
      {
        col_ind[i]=GetDVColumnInd(DV_SLACK,i);
        row_val[i]=_pConstraints[j]->penalty_over;
        i++;
        col_ind[i]=GetDVColumnInd(DV_SLACK,i);
        row_val[i]=_pConstraints[j]->penalty_under;
        i++;  
      } 
      else // only one slack var (>= or <=) 
      {
        col_ind[i]=GetDVColumnInd(DV_SLACK,i);
        row_val[i]=_pConstraints[j]->penalty_value;
        i++;
      }
    }
  }
 
  retval = lp_lib::set_obj_fnex(pLinProg,i, row_val, col_ind);
  ExitGracefullyIf(retval==0,"SolveDemandProblem: Error specifying objective function",RUNTIME_ERR);
  
  // ----------------------------------------------------------------
  // add mass balance constraints
  // ----------------------------------------------------------------
  
  // river routing mass balance constraints 
  // ----------------------------------------------------------------
  d=0; //counter for demands
  double U_0;
  int    nInlets;
  long   SBID;
  int    p_in[10]; //assumes<10 inlets
  dv_type dvtype;
  for (int pp = 0; pp<pModel->GetNumSubBasins(); pp++)
  {
    p  =pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);

    //Q_p = (Qin +Qin )*U0 + sum(Ui*Qn-i) + Runoff - Div - Dem + Qadded + U_0*Qadded2
    //Q_p-U_0*Qin-U_0*Qin2...+Div +Dem = sum(Ui*Qn-i) + Runoff + Qadded + U_0*Qadded2
    if (pSB->IsEnabled())
    {
      U_0=pSB->GetRoutingHydrograph()[0];
      nInlets=0;
      SBID=pSB->GetID();
      for (int ppp=0;ppp<pModel->GetNumSubBasins();ppp++) //todo - pre-process this or otherwise speedup
      {
        if (pModel->GetSubBasin(ppp)->GetDownstreamID() == SBID) {
          p_in[nInlets]=pModel->GetOrderedSubBasinIndex(ppp);
          nInlets++;
        }
      }

      i=0;
      d=0;
      col_ind[i]=GetDVColumnInd(DV_QOUT,p);
      row_val[i]=1.0;
      i++;
      for (int in = 0; in < nInlets; in++) 
      {
        dvtype=DV_QOUT;
        if (pModel->GetSubBasin(p)->GetReservoir() != NULL) {dvtype=DV_QOUTRES; }

        col_ind[i]=GetDVColumnInd(dvtype,p_in[in]);
        row_val[i]=-U_0;  
        i++;
      }
      for (int ii = 0; ii < pModel->GetSubBasin(p)->GetNumWaterDemands();ii++)
      {
        col_ind[i]=GetDVColumnInd(DV_DELIVERY,d);
        row_val[i]=+1.0;
        i++; d++;
      }
      // \todo[funct]: add flow diversions

      RHS=0;
      RHS+=pSB->GetSpecifiedInflow(t)*U_0;
      RHS+=pSB->GetDownstreamInflow(t);
      for (int n = 1; n < pSB->GetInflowHistorySize(); n++) {
        RHS+=pSB->GetRoutingHydrograph()[n]*pSB->GetInflowHistory()[n-1];  // [n-1] is because this inflow history has not yet been updated
      }
      for (int n = 1; n < pSB->GetInflowHistorySize(); n++) {
        RHS += pSB->GetUnitHydrograph()[n] * pSB->GetLatHistory()[n-1]; // [n-1] is because this inflow history has not yet been updated
      }
      RHS+=aSBrunoff[p]/(tstep*SEC_PER_DAY)*pSB->GetUnitHydrograph()[0];//[m3]->[m3/s]

      retval = lp_lib::add_constraintex(pLinProg,i,row_val,col_ind,ROWTYPE_EQ,RHS);
      ExitGracefullyIf(retval==0,"SolveDemandProblem::Error adding mass balance constraint",RUNTIME_ERR);
    }
  }
  
  // Lake mass balance constraints 
  //----------------------------------------------------------------
  CReservoir *pRes;
  int    res_count=0;
  int    k;
  double Adt,ET,precip,seepage(0);
  double stage,dQdh,Qout;
  for (int pp = 0; pp<pModel->GetNumSubBasins(); pp++)
  {
    p   =pModel->GetOrderedSubBasinIndex(pp);
    pSB =pModel->GetSubBasin(p);
    pRes=pSB->GetReservoir();

    if ((pSB->IsEnabled()) && (pRes != NULL))
    {
      //dV/dh*(hnew-hold)/dt=(Qin-Qout-ET+P-resextract-seep)
      // A/dt * h -1.0 * Qin + 1.0 * Qout - 1.0*resextract = P-ET-seep
      // Qout and resextract will often be on lhs as decision vars
      // must fix area at A^n-1 and stage used for flows at h^n-1
      
      Adt=pRes->GetSurfaceArea()/Options.timestep/SEC_PER_DAY; 
      k  =pRes->GetHRUIndex();
      ET =pModel->GetHydroUnit(k)->GetForcingFunctions()->OW_PET / SEC_PER_DAY / MM_PER_METER; //average for timestep, in m/s
      ET*=pModel->GetHydroUnit(k)->GetSurfaceProps()->lake_PET_corr;

      precip=pRes->GetReservoirPrecipGains(Options.timestep)/Options.timestep/SEC_PER_DAY; //[m3]->[m3/s]
      /*if (_seepage_const>0) {seepage=_seepage_const*(_stage-_local_GW_head); //[m3/s]}*/

      i=0;
      col_ind[i]=GetDVColumnInd(DV_QOUT   ,p);
      row_val[i]=-1.0;
      i++;

      col_ind[i]=GetDVColumnInd(DV_QOUTRES,res_count);
      row_val[i]=+1.0; 
      i++;

      col_ind[i]=GetDVColumnInd(DV_STAGE  ,res_count);
      row_val[i]=Adt;
      i++;

      RHS=precip-ET-seepage; 

      retval = lp_lib::add_constraintex(pLinProg,i,row_val,col_ind,ROWTYPE_EQ,RHS);
      ExitGracefullyIf(retval==0,"SolveDemandProblem::Error adding mass balance constraint",RUNTIME_ERR);

      //Need to add second constraint such that Qout=Q(h)  must linearize storage -discharge curve
      // Qout = Qoutold + dQ/dh * (h-hold) 
      // 1.0 * Qout - dQ/dh * h = rhs = Qoutold - dQ/dh * hold  

      //TODO: SHOULD THIS BE A GOAL INSTEAD?? MAY BE OVERRIDDEN BY OTHER CONSTRAINTS 
      
      stage=pRes->GetResStage();
      dQdh =pRes->GetStageDischargeDerivative(stage,nn);
      Qout =pRes->GetOutflowRate();

      i=0;
      col_ind[i]=GetDVColumnInd(DV_QOUTRES,res_count);
      row_val[i] = 1.0;  i++;

      col_ind[i]=GetDVColumnInd(DV_STAGE  ,res_count);
      row_val[i] = -dQdh; i++;

      RHS=Qout-dQdh*stage;

      retval = lp_lib::add_constraintex(pLinProg,i,row_val,col_ind,ROWTYPE_EQ,RHS);
      ExitGracefullyIf(retval==0,"SolveDemandProblem::Error adding mass balance constraint",RUNTIME_ERR);

      res_count++;
    }
  }

  // demand goals 
  // ----------------------------------------------------------------
  d=0;
  double enviro_minQ;
  int slack_index=0;
  for (int pp = 0; pp<pModel->GetNumSubBasins(); pp++)
  {
    p  =pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    
    if (pSB->IsEnabled())
    {
      //Q>E 
      // Q + S > E 
      // important- must have higher penalty than demand undersatisfaction penalty
      
      //LATER: also add unusable flow percentage as goal  \todo[funct]
      //double ufp=pSB->GetUnusableFlowPercentage();
      // sum(D)  <= (1-ufp)*((Q+D)-E) 
      // (1-ufp) * Q - (ufp) * sum(D) + (1-ufp) * S >= (1-ufp) * E
      enviro_minQ=pSB->GetEnviroMinFlow(t);

      i=0;
      col_ind[i]=GetDVColumnInd(DV_QOUT,p);
      row_val[i]=+1.0; i++;

      RHS = enviro_minQ;

      if (enviro_minQ != 0.0) 
      {
        col_ind[i]=GetDVColumnInd(DV_SLACK,slack_index);
        row_val[i]=+1.0; i++; 

        retval = lp_lib::add_constraintex(pLinProg,i,row_val,col_ind,ROWTYPE_GE,RHS);
        ExitGracefullyIf(retval==0,"SolveDemandProblem::Error adding environmental flow goal",RUNTIME_ERR);
      } 

      slack_index++; //because slack variables are applied to all, even ones with minQ!=0
    }
  }

  // user-specified constraints 
  // ----------------------------------------------------------------
  for (int i = 0; i < _nConstraints; i++) 
  {
    if (_pConstraints[i]->conditions_satisfied)
    {
      AddConstraintToLP( i, pLinProg, tt, col_ind, row_val);
    }
  }

  lp_lib::set_add_rowmode(pLinProg, FALSE); //turn off once model constraints and obj function are added

  // ----------------------------------------------------------------
  // SOLVE OPTIMIZATION PROBLEM WITH LP_SOLVE
  // ----------------------------------------------------------------
  lp_lib::set_minim(pLinProg); 
  lp_lib::set_verbose(pLinProg,IMPORTANT);
  //lp_lib::write_LP(pLinProg);
   
  //retval=lp_lib::set_scaling(lp, CURTISREIDSCALE);  //May wish to look for scaling improvements
  //ExitGracefullyIf(retval!=0,"SolveDemandProblem: unable to set scaling algorithm",RUNTIME_ERR);

  //retval=lp_lib::set_BFP(pLinProg, "bfp_LUSOL");  //TEST THIS FOR SPEED 
  //ExitGracefullyIf(retval!=0,"SolveDemandProblem: unable to set BFP",RUNTIME_ERR);

  lp_lib::print_debugdump(pLinProg,"lp_solve_dump.txt");
  lp_lib::print_lp(pLinProg);

  retval = lp_lib::solve(pLinProg);
  if (retval!=OPTIMAL){cout<<"code: "<<retval<<endl; }
  ExitGracefullyIf(retval!=OPTIMAL,"SolveDemandProblem: non-optimal solution found. Problem is over-constrained. Remove or adjust management constraints.",RUNTIME_ERR);
  
  //This will be sum of penalties*violations
  cout<<"Objective value: "<<lp_lib::get_objective(pLinProg)+demand_penalty_sum<<endl; 

  
  //assign decision variables
  // ----------------------------------------------------------------
  double *soln   =new double [_nDecisionVars];
  lp_lib::get_variables(pLinProg,soln); 

  for (int i=0;i<_nDecisionVars;i++){
    _pDecisionVars[i]->value=soln[i];
  }
  delete [] soln;

  //set state variables corresponding to decision variables 
  // SHOULD REALLY ONLY UPDATE DELIVERY, RES_EXTRACT, AND RESERVOIR OUTFLOW - ADD SWITCH TO HAVE OPT MODEL GENERATE STATE VARS OR REGULAR RAVEN SETUP 
  // ----------------------------------------------------------------
  double value;
  dv_type typ;
  int s=0;
  for (int i=0;i<_nDecisionVars;i++)
  {
    p     = _pDecisionVars[i]->p_index;
    value = _pDecisionVars[i]->value;
    typ   = _pDecisionVars[i]->dvar_type;
    
    if      (typ == DV_QOUT)
    {
      //pModel->GetSubBasin(p)->SetQout(value);
    } 
    else if (typ == DV_STAGE)
    {
      //double h_last=pModel->GetSubBasin(p)->GetReservoir()->GetResStage();
      //pModel->GetSubBasin(p)->GetReservoir()->SetReservoirStage(value,h_last);
    } 
    else if (typ == DV_DELIVERY)
    {
      int demand_ind = _pDecisionVars[i]->dem_index; 
      _aDelivery[demand_ind]=value;
      //
      //pModel->GetSubBasin(p)->SetDelivery(demand_ind,value);  //this can be handled within CSubBasin::ApplyIrrigationDemand instead   
    } 
    else if (typ == DV_QOUTRES)
    {
      //this is byproduct of stage determination - not passed back to model 
    } 
    else if (typ == DV_SLACK) {
      //s = _pDecisionVars[i]->slack_index;
      _aSlackValues[s]=value;
      s++;
    }
    else if (typ == DV_USER)
    {
      //doesn't influence hydrologic model 
    }
  }
  
  delete [] col_ind;
  delete [] row_val;

  lp_lib::delete_lp(pLinProg);
#endif
}

//////////////////////////////////////////////////////////////////
/// \brief Write output file headers
/// \details Called prior to simulation (but after initialization) from CModel::Initialize()
/// \param &Options [in] Global model options information
//
void CDemandOptimizer::WriteOutputFileHeaders(const optStruct &Options)
{
  // DemandOptimization.csv 
  //--------------------------------------------------------
  string tmpFilename=FilenamePrepare("DemandOptimization.csv",Options);
  _DEMANDOPT.open(tmpFilename.c_str());
  if (_DEMANDOPT.fail()){
    ExitGracefully(("CModel::WriteOutputFileHeaders: unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
  }

  _DEMANDOPT<<"time [d],date,hour";
  for (int i=0;i<_nDecisionVars;i++){
    _DEMANDOPT<<","<<_pDecisionVars[i]->name;
  }
  _DEMANDOPT<<endl;

  // GoalSatisfaction.csv 
  //--------------------------------------------------------
  tmpFilename=FilenamePrepare("GoalSatisfaction.csv",Options);
  _GOALSAT.open(tmpFilename.c_str());
  if (_GOALSAT.fail()){
    ExitGracefully(("CModel::WriteOutputFileHeaders: unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
  }
  _GOALSAT<<"time [d],date,hour"; 
  for (int i = 0; i < _nDemands; i++) {
    _GOALSAT << "," << _aDemandIDs[i];
  }
  for (int i = 0; i < _nConstraints; i++) {
    if (_pConstraints[i]->is_goal) {
      _GOALSAT<<","<<_pConstraints[i]->name;
    }
  }
  _GOALSAT<<endl;
}
//////////////////////////////////////////////////////////////////
/// \brief Writes minor output to file at the end of each timestep (or multiple thereof)
/// \param &Options [in] Global model options information
/// \param &tt [in] Local (model) time *at the end of* the pertinent time step
//
void CDemandOptimizer::WriteMinorOutput(const optStruct &Options,const time_struct &tt) 
{
  string usedate=tt.date_string;                   //refers to date and time at END of time step
  string usehour=DecDaysToHours(tt.julian_day);

  // DemandOptimization.csv 
  //--------------------------------------------------------
  _DEMANDOPT<<tt.model_time <<","<<usedate<<","<<usehour; 
  for (int i=0;i<_nDecisionVars;i++){
    _DEMANDOPT<<","<<_pDecisionVars[i]->value;
  }
  _DEMANDOPT<<endl;

  // GoalSatisfaction.csv 
  //--------------------------------------------------------
  int p;
  CSubBasin *pSB;
  double pen;
  int d=0;
  _GOALSAT<<tt.model_time <<","<<usedate<<","<<usehour;

  for (int pp=0;pp<_pModel->GetNumSubBasins();pp++){
    p = _pModel->GetOrderedSubBasinIndex(pp);
    pSB=_pModel->GetSubBasin(p);
    for (int ii = 0; ii < pSB->GetNumWaterDemands(); ii++) 
    {
      //_GOALSAT << "," << _aDemandPenalties[d] * (pSB->GetIrrigationDemand(tt.model_time)-_aDelivery[d]);
      _GOALSAT << "," << (pSB->GetIrrigationDemand(tt.model_time)-_aDelivery[d]);
      d++;
    }
  }
  int k=_nDemands; //first goals are due to environmental flow constraints 
  for (int i = 0; i < _nConstraints; i++) {
    if (_pConstraints[i]->is_goal) {
      if (_pConstraints[i]->pExpression->compare != COMPARE_IS_EQUAL) {
        //pen=_pConstraints[i]->penalty_value*_aSlackValues[k];
        pen = _aSlackValues[k];
        _GOALSAT<<","<<pen;
        k++;
      } 
      else {
        //pen=_pConstraints[i]->penalty_over *_aSlackValues[k] +
        //    _pConstraints[i]->penalty_under *_aSlackValues[k+1];
        pen =_aSlackValues[k]+_aSlackValues[k+1];
        _GOALSAT<<","<<pen;
        k+=2;
      }
    }
  }
  _GOALSAT<<endl;

}
//////////////////////////////////////////////////////////////////
/// \brief Closes output file streams
/// \details after end of simulation from Main() or in ExitGracefully; All file streams are opened in WriteOutputFileHeaders() routine
//
void CDemandOptimizer::CloseOutputStreams()
{
  if (_DEMANDOPT.is_open()){_DEMANDOPT.close();}
  if (  _GOALSAT.is_open()){  _GOALSAT.close();}
}