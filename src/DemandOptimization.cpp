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
  _aDemandSBIDs=NULL;
  _aDemandAliases=NULL;
  _aDelivery=NULL;
  _nReservoirs=0;
  _nEnabledSubBasins=0;

  _aSBIndices=NULL;
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

  _aUpstreamDemands=NULL;
  _aUpCount=NULL;

  _nHistoryItems=0;
  _aQhist=NULL;
  _aDhist=NULL;
  _ahhist=NULL;

  _aCumDelDate=NULL;
  _aCumDelivery=NULL;

  _nConstraints=0;
  _pConstraints=NULL;

  _do_debug_level=0;//no debugging

}
//////////////////////////////////////////////////////////////////
/// \brief Implementation of the Demand optimization destructor
//
CDemandOptimizer::~CDemandOptimizer()
{
  for (int i=0;i<_nDecisionVars;    i++){delete _pDecisionVars[i];    }delete [] _pDecisionVars;
  for (int i=0;i<_nUserTimeSeries;  i++){delete _pUserTimeSeries[i];  }delete [] _pUserTimeSeries;
  for (int i=0;i<_nUserLookupTables;i++){delete _pUserLookupTables[i];}delete [] _pUserLookupTables;

  for (int p=0;p<_nEnabledSubBasins; p++){delete [] _aUpstreamDemands[p]; } delete [] _aUpstreamDemands;
  delete [] _aUpCount;

  for (int p = 0; p < _nEnabledSubBasins;p++){
    if (_aQhist!=NULL){
      delete [] _aQhist[p];
      delete [] _aDhist[p];
      delete [] _ahhist[p];
    }
  }
  delete [] _aQhist;
  delete [] _ahhist;
  delete [] _aDhist;
  delete [] _aUserConstants;
  delete [] _aUserConstNames;
  delete [] _aDelivery;
  delete [] _aDemandIDs;
  delete [] _aDemandSBIDs;
  delete [] _aDemandAliases;
  delete [] _aSlackValues;
  delete [] _aSBIndices;
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
/// \brief returns debug level
//
int CDemandOptimizer::GetDebugLevel() const
{
  return _do_debug_level;
}
//////////////////////////////////////////////////////////////////
// \brief returns true if time series with specified name 'TSname' exists, false otherwise
// \params TSname [in] - name to be checked
//
bool CDemandOptimizer::UserTimeSeriesExists(string TSname) const
{
  for (int i = 0; i < _nUserTimeSeries; i++) {
    if (_pUserTimeSeries[i]->GetName()==TSname){return true;}
  }
  return false;
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
/// \brief sets debug level
/// \params lev [in] -level
//
void CDemandOptimizer::SetDebugLevel(const int val)
{
  _do_debug_level=val;
}
//////////////////////////////////////////////////////////////////
/// \brief date at which cumulative delivery is rebooted to zero (usually winter date)
/// \params julian_date [in] reboot date
/// \params demandID [in] demand identifier (e.g., 'D123')
/// only called from ParseManagementFile()
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

  if (pDV->dvar_type==DV_USER){
    _nUserDecisionVars++;
  }
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
  cout <<"ADDING USER CONSTANT "<<name<<" "<<val << endl;
}
//////////////////////////////////////////////////////////////////
/// \brief adds user time series

void CDemandOptimizer::AddUserTimeSeries(const CTimeSeries *pTS)
{
  if (!DynArrayAppend((void**&)(_pUserTimeSeries),(void*)(pTS),_nUserTimeSeries)){
   ExitGracefully("CDemandOptimizer::AddUserTimeSeries: adding NULL time series",BAD_DATA);}
}
//////////////////////////////////////////////////////////////////
/// \brief adds user lookup table
//
void CDemandOptimizer::AddUserLookupTable(const CLookupTable *pLUT)
{
  if (!DynArrayAppend((void**&)(_pUserLookupTables),(void*)(pLUT),_nUserLookupTables)){
   ExitGracefully("CDemandOptimizer::AddUserLookupTable: adding NULL lookup table",BAD_DATA);}
}
//////////////////////////////////////////////////////////////////
/// \brief sets demand penalty for demand dname
//
void CDemandOptimizer::SetDemandPenalty(const string dname, const double& pen) {

  int d=GetDemandIndexFromName(dname);
  ExitGracefullyIf(pen<0,"CDemandOptimizer::SetDemandPenalty: demand penalties must be greater or equal to zero",BAD_DATA_WARN);
  if (d!=DOESNT_EXIST){
    _aDemandPenalties[d]=pen;
  }
  else {
    ExitGracefullyIf(pen<0,"CDemandOptimizer::SetDemandPenalty: invalid demand identifier in :DemandPenalty command.",BAD_DATA_WARN);
  }
}

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
  pC->pExpression=exp;

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

  //populate _aSBindices, _aResIndices
  //------------------------------------------------------------------
  _nEnabledSubBasins=0;
  _nReservoirs      =0;
  _aSBIndices  = new int[pModel->GetNumSubBasins()];
  _aResIndices = new int[pModel->GetNumSubBasins()];
  for (int pp=0;pp<pModel->GetNumSubBasins();pp++)
  {
    p=pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    _aSBIndices [p]=DOESNT_EXIST;
    _aResIndices[p]=DOESNT_EXIST;
    if (pSB->IsEnabled()){
      _aSBIndices[p]=_nEnabledSubBasins;
      _nEnabledSubBasins++;
      if (pSB->GetReservoir()!=NULL){
        _aResIndices[p]=_nReservoirs;
        _nReservoirs++;
      }
    }
  }

  //Populate ordered decision variable array _pDecisionVars[]
  // This order has to be maintained because of how the decision variables are indexed (consistent with GetDVColumnInd())
  // subbasin outflows -> reservoir outflows -> reservoir stages -> delivered demand -> user-specified DVs -> slack variables
  //------------------------------------------------------------------
  // add subbasin outflow DVs
  int SB_count=0;
  for (int pp=0;pp<pModel->GetNumSubBasins();pp++)
  {
    p=pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    if (pSB->IsEnabled()){
      if (pSB->GetReservoir()!=NULL){name="I"+to_string(pSB->GetID());} //special case - inflow to reservoir is internally used only
      else                          {name="Q"+to_string(pSB->GetID());}
      pDV=new decision_var(name,p,DV_QOUT,SB_count);
      AddDecisionVar(pDV);
      SB_count++;
    }
  }
  // add reservoir outflow DVs
  int res_count=0;
  for (int pp=0;pp<pModel->GetNumSubBasins();pp++)
  {
    p=pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    if (pSB->IsEnabled()){
      if (pSB->GetReservoir()!=NULL){
        name="Q"+to_string(pSB->GetID());
        pDV=new decision_var(name,p,DV_QOUTRES,res_count);
        AddDecisionVar(pDV);
        res_count++;
      }
    }
  }
  // add reservoir stage DVs
  res_count=0;
  for (int pp=0;pp<pModel->GetNumSubBasins();pp++)
  {
    p=pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    if (pSB->IsEnabled()){
      if (pSB->GetReservoir()!=NULL){
        name="h"+to_string(pModel->GetSubBasin(p)->GetID());
        pDV=new decision_var(name,p,DV_STAGE,res_count);
        AddDecisionVar(pDV);
        res_count++;
      }
    }
  }

  // Count demands, add demand DVs
  //------------------------------------------------------------------
  _nDemands=0;
  for (p=0;p<pModel->GetNumSubBasins();p++)
  {
    pSB=pModel->GetSubBasin(p);
    if (pSB->IsEnabled()){
      for (int i = 0; i < pModel->GetSubBasin(p)->GetNumWaterDemands();i++) {_nDemands++;}
    }
  }
  for (p = 0; p < pModel->GetNumSubBasins(); p++)
  {
    pSB=pModel->GetSubBasin(p);
    if (pSB->IsEnabled()){
      if (pModel->GetSubBasin(p)->GetReservoir()!=NULL){
        for (int i = 0; i < pModel->GetSubBasin(p)->GetReservoir()->GetNumWaterDemands();i++) {_nDemands++;}
      }
    }
  }

  _aDemandIDs     = new int   [_nDemands];
  _aDemandAliases = new string[_nDemands];
  _aDemandSBIDs   = new long  [_nDemands];
  int d=0;
  int ID;
  for (int pp=0;pp<pModel->GetNumSubBasins();pp++)
  {
    p=pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    if (pSB->IsEnabled()){
      for (int ii = 0; ii < pSB->GetNumWaterDemands();ii++)
      {
        ID   = pSB->GetWaterDemandID  (ii);//e.g., 132 (not i)
        name = pSB->GetWaterDemandName(ii);
        pDV=new decision_var(name,p,DV_DELIVERY,d);
        pDV->dem_index=ii;
        AddDecisionVar(pDV);

        _aDemandIDs    [d]=ID;
        _aDemandAliases[d]=name;
        _aDemandSBIDs  [d]=pSB->GetID();
        d++;
      }
    }
  }
  for (p = 0; p < pModel->GetNumSubBasins(); p++)
  {
    pSB=pModel->GetSubBasin(p);
    if (pSB->IsEnabled()){
      if (pModel->GetSubBasin(p)->GetReservoir()!=NULL){
        for (int ii = 0; ii < pModel->GetSubBasin(p)->GetReservoir()->GetNumWaterDemands();ii++)
        {
          ID   = pSB->GetReservoir()->GetWaterDemandID  (ii);//e.g., 132
          name = pSB->GetReservoir()->GetWaterDemandName(ii);
          pDV=new decision_var(name,p,DV_DELIVERY,d);
          pDV->dem_index=100+ii; //Valid? or should this be like 100+ii to indicate reservoir?
          AddDecisionVar(pDV);

          _aDemandIDs    [d]=ID;
          _aDemandAliases[d]=name;
          _aDemandSBIDs  [d]=pSB->GetID();
          d++;
        }
      }
    }
  }

  _aDelivery        = new double[_nDemands];
  _aDemandPenalties = new double[_nDemands];
  _aCumDelivery     = new double[_nDemands];
  _aCumDelDate      = new int   [_nDemands];
  for (int d=0;d<_nDemands;d++){
    _aDelivery       [d]=0;
    _aDemandPenalties[d]=1.0;  //Default demand penalty (make accessible parameter?)
    _aCumDelivery    [d]=0.0;
    _aCumDelDate     [d]=0; //Jan 1
  }

  IdentifyUpstreamDemands();
}
//////////////////////////////////////////////////////////////////
/// \brief adds value v to end of array a[] of original size n, expands a[] to size n+1
//
void push(int*&a, const int &v, int &n)
{
  int *tmp=new int [n+1];
  for (int i = 0; i < n; i++) { tmp[i]=a[i];}
  tmp[n]=v;
  a=tmp;
  if (n>0){delete [] a;}
  n++;
}
//////////////////////////////////////////////////////////////////
/// \brief populates _aUpstreamDemands and _aUpCount arrays
/// called by CDemandOptimizer::Initialize()
//
void  CDemandOptimizer::IdentifyUpstreamDemands()
{
  long SBID;

  _aUpCount        =new int  [_pModel->GetNumSubBasins()];
  _aUpstreamDemands=new int *[_pModel->GetNumSubBasins()];
  for (int p=0;p<_pModel->GetNumSubBasins();p++){
    _aUpstreamDemands[p]=NULL;
    _aUpCount[p]=0;
  }

  for (int d = 0; d < _nDemands; d++)
  {
    SBID=_aDemandSBIDs[d];
    while (SBID >= 0) {
      int p=_pModel->GetSubBasinByID(SBID)->GetGlobalIndex();

      if (_pModel->GetSubBasin(p)->IsEnabled()){
        push(_aUpstreamDemands[p],d,_aUpCount[p]);
        SBID=_pModel->GetSubBasin(p)->GetDownstreamID();
      }
      else {
        SBID=DOESNT_EXIST;
      }
    }
  }

  cout<<" ** ** ** ** ** ** ** ** ** ** ** "<<endl;
  cout<<"UPSTREAM DEMANDS: "<<endl;
  for (int p = 0; p < _pModel->GetNumSubBasins(); p++)
  {
    if (_pModel->GetSubBasin(p)->IsEnabled()){
      cout <<"   basin "<<p<<" SBID= "<<_pModel->GetSubBasin(p)->GetID()<<"): ";
      for (int i = 0; i < _aUpCount[p]; i++) {
        cout<<_aUpstreamDemands[p][i]<<" ";
      }
      cout<<endl;
    }
  }
  cout<<" ** ** ** ** ** ** ** ** ** ** ** "<<endl;
}
//////////////////////////////////////////////////////////////////
/// \brief Initializes Demand optimization instance
/// \notes to be called after .rvm file read
/// \params pModel [in] - pointer to model
/// \params Options  [in] - model options structure
/// note: Slack variable indices MUST be ordered in same order as used in assembly in lp_solve matrix within SolveDemandProblem()
///   1) reservoir outflow goal
///   2) user-specified goals
//
void CDemandOptimizer::InitializePostRVMRead(CModel* pModel, const optStruct& Options)
{
  if (Options.noisy){cout<<"CDemandOptimizer: Post-rvm-read initialization..."<<endl;}

  // initialize history arrays
  //------------------------------------------------------------------
  if (_nHistoryItems==0){
    _nHistoryItems=static_cast<int>(rvn_floor(1.0/Options.timestep+REAL_SMALL));//default = 1 day
  }
  _aDhist = NULL;
  _aQhist = new double *[_nEnabledSubBasins];
  _ahhist = new double *[_nEnabledSubBasins];
  _aDhist = new double *[_nEnabledSubBasins];
  ExitGracefullyIf(_aDhist==NULL,"CDemandOptimizer::InitializePostRVMRead",OUT_OF_MEMORY);
  for (int pp=0;pp<_nEnabledSubBasins;pp++){
    _aDhist[pp]=NULL;
    _aQhist[pp] = new double[_nHistoryItems];
    _ahhist[pp] = new double[_nHistoryItems];
    _aDhist[pp] = new double[_nHistoryItems];
    ExitGracefullyIf(_aDhist[pp]==NULL,"CDemandOptimizer::InitializePostRVMRead (2)",OUT_OF_MEMORY);
    for (int i = 0; i < _nHistoryItems; i++) {
      _aQhist[pp][i]=0.0;
      _ahhist[pp][i]=0.0;
      _aDhist[pp][i]=0.0;
    }
  }

  int p;
  _nSlackVars = 0;
  decision_var *pDV;
  CSubBasin *pSB;

  // create 2 slack variables for reservoir outflow by stage-discharge
  //------------------------------------------------------------------
  for (int pp=0;pp<pModel->GetNumSubBasins();pp++)
  {
    p=pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    if ((pSB->IsEnabled()) && (pSB->GetReservoir()!=NULL)){
      pDV =  new decision_var("ResSlackUp" + to_string(pModel->GetSubBasin(p)->GetID()), p, DV_SLACK,_nSlackVars);
      AddDecisionVar(pDV);
      _nSlackVars++;
      pDV =  new decision_var("ResSlackDn" + to_string(pModel->GetSubBasin(p)->GetID()), p, DV_SLACK,_nSlackVars);
      AddDecisionVar(pDV);
      _nSlackVars++;
    }
  }


  // create 1 slack variable each for environmental min flow goal and ufp
  //------------------------------------------------------------------
  for (int pp = 0; pp<pModel->GetNumSubBasins(); pp++)
  {
    p=pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    if (pSB->IsEnabled()){
      //created for all basins with environmental min flow or non-zero unusable flow pct
      if (pSB->HasEnviroMinFlow()){
        pDV =  new decision_var("E" + to_string(pModel->GetSubBasin(p)->GetID()), p, DV_SLACK,_nSlackVars);
        AddDecisionVar(pDV);
        _nSlackVars++;
      }
      /*if (pSB->GetUnusableFlowPercentage()>0) {
        pDV =  new decision_var("U" + to_string(pModel->GetSubBasin(p)->GetID()), p, DV_SLACK,_nSlackVars);
        AddDecisionVar(pDV);
        _nSlackVars++;
      }*/
    }
  }

  // Convert reservoir commands to user-specified constraints
  //------------------------------------------------------------------
  AddReservoirConstraints();

  // Sift through user-specified goals, update _nSlackVars, reserve memory for _aSlackValues
  //------------------------------------------------------------------
  // ISSUE : THIS IS TOTAL POSSIBLE NUMBER OF SLACK VARIABLES; MAY NOT BE IN PLAY DUE TO CONDITIONALS - for disabled constraints, just set penalty to zero>??
  //
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
  for (int ii = 0; ii < _nSlackVars; ii++) {
    _aSlackValues[ii]=0.0;
  }

  // Initialize user specified time series
  //------------------------------------------------------------------
  for (int i = 0; i < _nUserTimeSeries; i++)
  {
    _pUserTimeSeries[i]->Initialize(Options.julian_start_day,
                                    Options.julian_start_year,
                                    Options.duration,
                                    Options.timestep, true, Options.calendar); //is_observation=true allows for blanks in the time series
  }

  // Print summary to screen
  //------------------------------------------------------------------
  if (true)
  {
    long SBID;
    int d;

    cout<<" -----------------------------------------------------------------"<<endl;
    cout<<"  MANAGEMENT OPTIMIZATION SUMMARY                                 "<<endl;
    cout<<" -----------------------------------------------------------------"<<endl;
    cout<<" # Enabled subbasins: "<<_nEnabledSubBasins<<" of "<<pModel->GetNumSubBasins()<<endl;
    cout<<" # Enabled reservoirs/lakes: "<<_nReservoirs<<endl;
    cout<<" # Decision Vars:      "<<_nDecisionVars    <<endl;
    cout<<" # User decision vars: "<<_nUserDecisionVars<<endl;
    for (int i = 0; i < _nDecisionVars; i++)
    {
      cout<<"    DV "<<setw(4)<<i<<": "<<_pDecisionVars[i]->name<<" "<<_pDecisionVars[i]->dvar_type<<" "<<GetDVColumnInd(_pDecisionVars[i]->dvar_type,_pDecisionVars[i]->loc_index)<<" "<<_pDecisionVars[i]->loc_index<<endl;
    }
    cout<<" # Demands : "<<_nDemands<<endl;

    d=0;
    for (int pp=0;pp<pModel->GetNumSubBasins();pp++)
    {
      p   =pModel->GetOrderedSubBasinIndex(pp);
      pSB =pModel->GetSubBasin(p);
      SBID=pSB->GetID();
      if (pSB->IsEnabled()){
        for (int ii = 0; ii < pSB->GetNumWaterDemands();ii++)
        {
          cout << "    " << d << ": " << _aDemandIDs[d] << " (alias: "<<_aDemandAliases[d]<<") in basin "<< SBID <<endl;
          d++;
        }
      }
    }

    cout<<" # History Items: "     <<_nHistoryItems<<endl;
    cout<<" # Slack Vars: "        <<_nSlackVars<<endl;
    cout<<" # Constraints/Goals: " <<_nConstraints<<endl;
    for (int i = 0; i < _nConstraints; i++) {
      cout<<"    "<<i<<" : "<<_pConstraints[i]->name<<" "<<_pConstraints[i]->pExpression->origexp<<endl;
    }
    cout<<" -----------------------------------------------------------------"<<endl;
    cout<<" -----------------------------------------------------------------"<<endl;
  }
  if (Options.noisy){cout<<"   ...end Post-rvm-read initialization."<<endl;}
}


//////////////////////////////////////////////////////////////////
// converts string instring to tokenized array s (preallocated) with Len tokens
//
void TokenizeString(string instring, char **s, int &Len)
{
  // Returns first token
    static char str[256];
    strcpy(str,instring.c_str());
    char *token = strtok(str, " ");

    int i=0;
    while (token != NULL)
    {
        s[i]=token;
        printf("%s\n", token);
        token = strtok(NULL, " ");
        //cout <<"TOKENIZE: "<<s[i] << endl;
        i++;
    }
    Len=i;
}
//////////////////////////////////////////////////////////////////
// adds reservoir constraints as 'user-specified' constraints
// these are done as if the corresponding constraints were read in as expressions from the .rvm file
//
void CDemandOptimizer::AddReservoirConstraints()
{
  int p;
  CSubBasin *pSB;
  string TSname,SBIDs;
  long SBID;
  expressionStruct *exp;
  dv_constraint    *pConst=NULL;
  string expString;
  string advice;
  double stage_units_adjustment;
  int Len=0;
  char *s[MAXINPUTITEMS];

  for (int pp=0;pp<_pModel->GetNumSubBasins();pp++)
  {
    p   =_pModel->GetOrderedSubBasinIndex(pp);
    pSB =_pModel->GetSubBasin(p);
    SBID=pSB->GetID();
    SBIDs=to_string(SBID);

    if (((pSB->IsEnabled()) && (pSB->GetReservoir()!=NULL)))
    {
      int k=pSB->GetReservoir()->GetHRUIndex();
      if (k!=DOESNT_EXIST){
        stage_units_adjustment=(_pModel->GetHydroUnit(k)->GetArea()*M2_PER_KM2)/SEC_PER_DAY; //[m]/[dt (d)]->[m3/s]
      }
      else {
        advice="The reservoir in subbasin "+SBIDs+" doesnt have an :HRUID - this will negatively impact the units scaling of management optimization penalties.";
        WriteWarning(advice,false);
      }

      //Max Stage constraints
      // h <= h_max
      //------------------------------------------------------------------------
      TSname="_MaxStage_" + SBIDs;
      if (UserTimeSeriesExists(TSname))
      {
        expString=":Expression !h" + SBIDs+ " < @ts(" + TSname + ",0)";
        TokenizeString(expString, s, Len);
        exp = ParseExpression((const char**)(s), Len, 0, "internal");
        pConst=AddConstraint(TSname, exp, true);
        pConst->penalty_over=6.0*stage_units_adjustment;
        advice="A :ReservoirMaxStage time series for subbasin "+SBIDs+" has been added to the management optimization formulation";
        WriteAdvisory(advice,false);
      }

      //Min Stage constraints
      // h >= h_min
      //------------------------------------------------------------------------
      TSname="_MinStage_" + SBIDs;
      if (UserTimeSeriesExists(TSname))
      {
        expString=":Expression !h" + SBIDs+ " > @ts(" + TSname + ",0)";
        TokenizeString(expString, s, Len);
        exp = ParseExpression((const char**)(s), Len, 0, "internal");
        pConst=AddConstraint(TSname, exp, true);
        pConst->penalty_over=6.0*stage_units_adjustment;
        advice="A :ReservoirMinStage time series for subbasin "+SBIDs+" has been added to the management optimization formulation";
        WriteAdvisory(advice,false);
      }

      //Min Flow constraints
      // Q >= Q_min
      //------------------------------------------------------------------------
      TSname="_ResQmin_" + SBIDs;
      if (UserTimeSeriesExists(TSname))
      {
        expString=":Expression !Q" + SBIDs+ " > @ts(" + TSname + ",0)";
        TokenizeString(expString, s, Len);
        exp = ParseExpression((const char**)(s), Len, 0, "internal");
        pConst=AddConstraint(TSname, exp, true);
        pConst->penalty_over=5.0;
        advice="A :ReservoirMinimumFlow time series for subbasin "+SBIDs+" has been added to the management optimization formulation";
        WriteAdvisory(advice,false);
      }

      //Max Flow constraints
      // Q <= Q_max
      //------------------------------------------------------------------------
      TSname="_ResQmax_" + SBIDs;
      if (UserTimeSeriesExists(TSname))
      {
        expString = ":Expression !Q" + SBIDs + " < @ts(" + TSname + ",0)";
        TokenizeString(expString, s, Len);
        exp = ParseExpression((const char**)(s), Len, 0, "internal");
        pConst=AddConstraint(TSname, exp, true);
        pConst->penalty_over=5.0;
        advice = "A :ReservoirMaximumFlow time series for subbasin " + SBIDs + " has been added to the management optimization formulation";
        WriteAdvisory(advice,false);
      }

      // Override reservoir flow
      // Q == Q_spec
      //------------------------------------------------------------------------
      TSname="_ResFlow_" + SBIDs;
      if (UserTimeSeriesExists(TSname))
      {
        expString=":Expression !Q" + SBIDs+ " = @ts(" + TSname + ",0)";
        TokenizeString(expString, s, Len);
        exp = ParseExpression((const char**)(s), Len, 0, "internal");
        pConst=AddConstraint(TSname, exp, true);
        pConst->penalty_over=4.0;
        advice = "An :OverrideReservoirFlow time series for subbasin " + SBIDs + " has been added to the management optimization formulation";
        WriteAdvisory(advice,false);
      }

      // Maximum Q delta
      // Q^{n+1} < dQ/dt * dt +Q^{n}
      //------------------------------------------------------------------------
      TSname="_MaxQDelta_" + SBIDs;
      if (UserTimeSeriesExists(TSname))
      {
        expString = ":Expression !Q" + SBIDs + " < @ts(" + TSname + ",0) * 86400 +!Q"+SBIDs+"[-1]"; //TBD - check if this is !Q[0] or !Q[-1]
        //expString = "!q" + SBIDs + " < @ts(" + TSname + ",0)";
        TokenizeString(expString, s, Len);
        exp = ParseExpression((const char**)(s), Len, 0, "internal");
        pConst=AddConstraint(TSname, exp, true);
        pConst->penalty_over=3.0;
        advice = "An :ReservoirMaxQDelta time series for subbasin " + SBIDs + " has been added to the management optimization formulation";
        WriteAdvisory(advice,false);
      }

      // Max Q Decrease
      // Q^{n+1} > (-dQ/dt) * dt +Q^{n}
      //------------------------------------------------------------------------
      TSname="_MaxQDecrease_" + SBIDs;
      if (UserTimeSeriesExists(TSname))
      {
        expString = ":Expression !Q" + SBIDs + " > @ts(" + TSname + ",0) * -86400 +!Q"+SBIDs+"[-1]"; //TBD - check if this is !Q[0] or !Q[-1]
        //expString = "!q" + SBIDs + " < @ts(" + TSname + ",0)";
        TokenizeString(expString, s, Len);
        exp = ParseExpression((const char**)(s), Len, 0, "internal");
        pConst=AddConstraint(TSname, exp, true);
        pConst->penalty_over=3.0;
        advice = "An :ReservoirMaxQDecrease time series for subbasin " + SBIDs + " has been added to the management optimization formulation";
        WriteAdvisory(advice,false);
      }

      // Target Stage
      // h == h_target
      //------------------------------------------------------------------------
      TSname="_TargetStage_" + SBIDs;
      if (UserTimeSeriesExists(TSname))
      {
        expString=":Expression !h" + SBIDs+ " = @ts(" + TSname + ",0)";
        cout<<"TOKENIZING: "<<expString << endl;
        TokenizeString(expString, s, Len);
        cout<<"PARSING: "<<expString << endl;
        exp = ParseExpression((const char**)(s), Len, 0, "internal");
        cout<<"ADDING CONSTRAINT: "<<expString << endl;
        pConst=AddConstraint(TSname, exp, true);
        pConst->penalty_over=2.0*stage_units_adjustment;
        advice = "An :ReservoirTargetStage time series for subbasin " + SBIDs + " has been added to the management optimization formulation";
        WriteAdvisory(advice,false);
      }
    }
  }

  //TO DO: Initialize all of the time series above.
}
//////////////////////////////////////////////////////////////////
// indexing for DV vector
// 0     ..   nSB-1 - outflow of stream reach p
// nSB   ..   nSB+nRes-1 - reservoir outflows
// nSB+nRes.. nSB+2*Res-1 - reservoir stages
// nSB+2*nRes .. nSB+2*nRes+nDem-1 - satisfied demand
// nSB+2*nRes+nDem .. nSB+2*nRes+nDem+nDV - user-specified decision variables
// remaining - slack vars - [#enviro min flow+ # user <=/>= goals + 2* # user == goals + all reservoir goals ]
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
void CDemandOptimizer::UpdateHistoryArrays() {
  int pp=0;
  int p;
  for (int ppp=0;ppp<_pModel->GetNumSubBasins();ppp++)
  {
    p=_pModel->GetOrderedSubBasinIndex(ppp);
    CSubBasin *pSB=_pModel->GetSubBasin(p);
    if (pSB->IsEnabled()){
      for (int i = _nHistoryItems-1; i>0; i--) {
        _aQhist[pp][i]=_aQhist[pp][i-1];
        _aDhist[pp][i]=_aDhist[pp][i-1];
      }
      _aQhist[pp][0]=pSB->GetOutflowRate();
      if (pSB->GetReservoir()!=NULL){
        _ahhist[pp][0]=pSB->GetReservoir()->GetResStage();
      }
      else {
        _ahhist[pp][0]=0.0;
      }
      _aDhist[pp][0]=0;
      for (int ii=0;ii<pSB->GetNumWaterDemands();ii++){
        _aDhist[pp][0]+=pSB->GetDemandDelivery(ii);
      }
      pp++;
    }
  }
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

  int    d,i,p,s;
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

  double *h_iter=new double [_pModel->GetNumSubBasins()];
  int    *lprow =new int    [_pModel->GetNumSubBasins()]; //index of goal equation for non-linear reservoir stage discharge curve in subbasin p

  // update history arrays from previous timestep
  // ----------------------------------------------------------------
  UpdateHistoryArrays();

  // instantiate linear programming solver
  // ----------------------------------------------------------------
  lp_lib::lprec *pLinProg;

  pLinProg=lp_lib::make_lp(0,_nDecisionVars);
  if (pLinProg==NULL){ExitGracefully("Error in SolveDemandProblem(): couldn construct new linear programming model",RUNTIME_ERR);}

  char name[200];
  strcpy(name,"RavenLPSolve");
  lp_lib::set_lp_name(pLinProg,name);

  for (int i=0;i<_nDecisionVars;i++)
  {
    strcpy(name,_pDecisionVars[i]->name.c_str());
    lp_lib::set_col_name(pLinProg, i + 1, name);
  }

  // Set upper bounds of delivery (D) to demand (D*) (preferred to adding constraint)
  // ----------------------------------------------------------------
  d=0;
  for (int pp = 0; pp<pModel->GetNumSubBasins(); pp++)
  {
    p  =pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    if (pSB->IsEnabled()){
      for (int ii = 0; ii < pSB->GetNumWaterDemands();ii++)
      {
        retval=lp_lib::set_upbo(pLinProg,GetDVColumnInd(DV_DELIVERY,d), pSB->GetWaterDemand(ii,t));
        ExitGracefullyIf(retval!=1,"SolveDemandProblem::Error adding demand upper bound",RUNTIME_ERR);
        d++;
      }
    }
  }

  // Set lower bounds of stages to -1000m
  // ----------------------------------------------------------------
  int res_count=0;
  for (int pp = 0; pp<pModel->GetNumSubBasins(); pp++)
  {
    p  =pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    if (pSB->IsEnabled() && (pSB->GetReservoir()!=NULL)){
        retval=lp_lib::set_lowbo(pLinProg,GetDVColumnInd(DV_STAGE,res_count), -1000);
        ExitGracefullyIf(retval!=1,"SolveDemandProblem::Error adding stage lower bound",RUNTIME_ERR);
       res_count++;
    }
  }

  // Set bounds of user-specified decision variables
  // ----------------------------------------------------------------
  int userct=0;
  for (int i = 0; i < _nDecisionVars; i++) {
    if (_pDecisionVars[i]->dvar_type==DV_USER){
      if (_pDecisionVars[i]->max < ALMOST_INF * 0.99) {
        retval=lp_lib::set_upbo(pLinProg,GetDVColumnInd(DV_USER,userct), _pDecisionVars[i]->max);
        ExitGracefullyIf(retval!=1,"SolveDemandProblem::Error adding decision variable upper bound",RUNTIME_ERR);
      }
      if (_pDecisionVars[i]->min != 0.0) {
        retval=lp_lib::set_lowbo(pLinProg,GetDVColumnInd(DV_USER,userct), _pDecisionVars[i]->min);
        ExitGracefullyIf(retval!=1,"SolveDemandProblem::Error adding decision variable lower bound",RUNTIME_ERR);
      }
      userct++;
    }
  }

  lp_lib::set_add_rowmode(pLinProg, TRUE); //readies lp_lib to add rows and objective function

  // ----------------------------------------------------------------
  // create objective function : minimize(sum(penalty*undelivered demand) + sum(penalty*violations (slack vars)))
  // ----------------------------------------------------------------
  i=0;

  //add blank penalty to objective function for models without demand -> 0*Q_0-->min
  // ----------------------------------------------------------------
  col_ind[i]=GetDVColumnInd(DV_QOUT,0);
  row_val[i]=0.0;
  i++;


  d=0;
  // penalties for unmet demand
  // ----------------------------------------------------------------
  for (int pp = 0; pp<pModel->GetNumSubBasins(); pp++)
  {
    p  =pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    if (pSB->IsEnabled()){
      for (int ii = 0; ii < pSB->GetNumWaterDemands();ii++)
      {
        demand_penalty_sum+=(_aDemandPenalties[d]*pSB->GetWaterDemand(ii,t));

        col_ind[i]=GetDVColumnInd(DV_DELIVERY,d);
        row_val[i]=-_aDemandPenalties[d]; //this maximizes delivery
        i++;
        d++;
      }
    }
  }
  //ADD RESERVOIR EXTRACTION DEMAND HERE

  s=0;
  // reservoir discharge curve goal penalties
  // ----------------------------------------------------------------
  for (int pp = 0; pp<pModel->GetNumSubBasins(); pp++)
  {
    p  =pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    if ((pSB->IsEnabled()) && (pSB->GetReservoir()!=NULL))
    {
        col_ind[i]=GetDVColumnInd(DV_SLACK,s);
        row_val[i]=0.5;
        i++; s++;
        col_ind[i]=GetDVColumnInd(DV_SLACK,s);
        row_val[i]=0.5;
        i++; s++;
    }
  }

  // enviro min flow goal penalties
  // ----------------------------------------------------------------
  for (int pp = 0; pp<pModel->GetNumSubBasins(); pp++)
  {
    p  =pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    if (pSB->IsEnabled()){
      if (pSB->GetEnviroMinFlow(t)>0.0)  {
        col_ind[i]=GetDVColumnInd(DV_SLACK,s);
        row_val[i]=1.0;
        i++; s++;
      }
      /*if (pSB->GetUnusableFlowPercentage()>0) {
        col_ind[i]=GetDVColumnInd(DV_SLACK,s);
        row_val[i]=1.0;
        i++; s++;
      }*/
      /*if ((pSB->GetNumWaterDemands()>0) && ((pSB->GetUnusableFlowPercentage()>0) || (pSB->GetEnviroMinFlow(t)>0))) {
        col_ind[i]=GetDVColumnInd(DV_SLACK,s);
        row_val[i]=10.0; //Is this a fair penalty?- needs to be stronger than all upstream demand penaltys. make a parameter?
        i++; s++;
      }*/
    }
  }

  // user-specified goal penalties
  // ----------------------------------------------------------------
  for (int j = 0; j < _nConstraints; j++)
  {
    _pConstraints[j]->conditions_satisfied=CheckConditions(j,tt);

    if ((_pConstraints[j]->is_goal) && (_pConstraints[j]->conditions_satisfied))
    {
      if (_pConstraints[j]->pExpression->compare == COMPARE_IS_EQUAL) //two slack variables
      {
        col_ind[i]=GetDVColumnInd(DV_SLACK,s); //
        row_val[i]=_pConstraints[j]->penalty_over;
        i++; s++;
        col_ind[i]=GetDVColumnInd(DV_SLACK,s);
        row_val[i]=_pConstraints[j]->penalty_under;
        i++; s++;
      }
      else // only one slack var (>= or <=)
      {
        col_ind[i]=GetDVColumnInd(DV_SLACK,s);
        row_val[i]=_pConstraints[j]->penalty_over;
        i++; s++;
      }
    }
  }

  retval = lp_lib::set_obj_fnex(pLinProg,i, row_val, col_ind);
  ExitGracefullyIf(retval==0,"SolveDemandProblem: Error specifying objective function",RUNTIME_ERR);

  // ----------------------------------------------------------------
  // add mass balance constraints
  // ----------------------------------------------------------------

  // river routing mass balance constraints (WORKING)
  // ----------------------------------------------------------------
  d=0; //counter for demands
  s=0; //counter for slack vars

  double  U_0;
  int     nInlets;
  long    SBID;
  int     p_in[10]; //assumes<10 inlets
  int     p2;
  dv_type dvtype;

  for (int pp = 0; pp<pModel->GetNumSubBasins(); pp++)
  {
    p  =pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);

    //Q_p = (Qin +Qin )*U0 + sum(Ui*Qn-i) + Runoff - Div - Dem + Qadded + U_0*Qadded2
    //Q_p-U_0*Qin-U_0*Qin2...+Div +Dem = sum(Ui*Qn-i) + Runoff + Qadded + U_0*Qadded2
    if (pSB->IsEnabled())
    {
      U_0 =pSB->GetRoutingHydrograph()[0];
      SBID=pSB->GetID();

      nInlets=0;
      for (int ppp=0;ppp<pModel->GetNumSubBasins();ppp++) //todo[optimize] - pre-process this or otherwise speedup _aSBInlets[p][ppp], _aSBInletCount[p]
      {
        p2=pModel->GetOrderedSubBasinIndex(ppp);
        if (pModel->GetSubBasin(p2)->GetDownstreamID() == SBID) {
          p_in[nInlets]=p2;
          nInlets++;
        }
      }

      i=0;

      col_ind[i]=GetDVColumnInd(DV_QOUT,_aSBIndices[p]);
      row_val[i]=1.0;
      i++;

      int id;
      for (int in = 0; in < nInlets; in++)
      {
        dvtype=DV_QOUT; id=_aSBIndices[p_in[in]];
        if (pModel->GetSubBasin(p_in[in])->GetReservoir() != NULL) {
          dvtype=DV_QOUTRES; id=_aResIndices[p_in[in]];
        }
        //ExitGracefullyIf(id==DOESNT_EXIST,"SolveDemandProblem: invalid index",RUNTIME_ERR);
        col_ind[i]=GetDVColumnInd(dvtype,id);
        row_val[i]=-U_0;
        i++;
      }
      for (int dd = 0; dd<pSB->GetNumWaterDemands();dd++)
      {
        col_ind[i]=GetDVColumnInd(DV_DELIVERY,d);
        row_val[i]=+1.0;
        i++; d++;
      }
      if (pSB->GetReservoir() != NULL) {
        for (int dd = 0; dd<pSB->GetReservoir()->GetNumWaterDemands();dd++)
        {
          d++;
        }
      }
      // \todo[funct]: add flow diversions

      RHS=0;
      RHS+=pSB->GetSpecifiedInflow (t)*U_0;
      RHS+=pSB->GetDownstreamInflow(t);
      for (int n = 1; n < pSB->GetInflowHistorySize(); n++) {
        RHS+=pSB->GetRoutingHydrograph()[n]*pSB->GetInflowHistory()[n-1]; // [n-1] is because this inflow history has not yet been updated
      }
      for (int n = 1; n < pSB->GetInflowHistorySize(); n++) {
        RHS += pSB->GetUnitHydrograph()[n] * pSB->GetLatHistory()[n-1];   // [n-1] is because this inflow history has not yet been updated
      }
      RHS+=aSBrunoff[p]/(tstep*SEC_PER_DAY)*pSB->GetUnitHydrograph()[0];  // [m3]->[m3/s]

      retval = lp_lib::add_constraintex(pLinProg,i,row_val,col_ind,ROWTYPE_EQ,RHS);
      ExitGracefullyIf(retval==0,"SolveDemandProblem::Error adding mass balance constraint",RUNTIME_ERR);
    }
  }

  // Lake mass balance constraints (WORKING - EVEN NON-LINEAR!)
  //----------------------------------------------------------------
  CReservoir *pRes;
  int    k;
  double Adt,ET,precip,seepage(0);
  double h_old,dQdh,Qout_last,Qin_last;
  res_count=0;
  d=0;
  for (int pp = 0; pp<pModel->GetNumSubBasins(); pp++)
  {
    p   =pModel->GetOrderedSubBasinIndex(pp);
    pSB =pModel->GetSubBasin(p);
    pRes=pSB->GetReservoir();
    if (pSB->IsEnabled()){
      for (int dd = 0; dd<pSB->GetNumWaterDemands();dd++)
      {
        d++;
      }
    }
    if ((pSB->IsEnabled()) && (pRes != NULL))
    {
      // First constraint equation (mass balance):
      // ------------------------------------------------------------------------
      //dV/dh*(hnew-hold)/dt=(Qin-Qout-ET*A+P*A-sum demands-seep)
      // A/dt * h -1.0 * Qin + 1.0 * Qout + 1.0*{demands} = P*A-ET*A-seep +A/dt* hold
      // must fix area at A^n-1 and stage used for flows at h^n-1

      h_old    =pRes->GetResStage();
      Qout_last=pRes->GetOutflowRate();
      Qin_last =pSB->GetReservoirInflow();
      Adt      =pRes->GetSurfaceArea()/Options.timestep/SEC_PER_DAY;
      k        =pRes->GetHRUIndex();

      if (pRes->GetHRUIndex() != NULL) {
        ET   =pModel->GetHydroUnit(k)->GetForcingFunctions()->OW_PET / SEC_PER_DAY / MM_PER_METER; //average for timestep, in m/s
        ET  *=pModel->GetHydroUnit(k)->GetSurfaceProps()->lake_PET_corr;
        ET  *=pRes->GetSurfaceArea(); //m3/s
      }

      precip=pRes->GetReservoirPrecipGains(Options.timestep)/Options.timestep/SEC_PER_DAY; //[m3]->[m3/s] (confirmed: this is properly updated before MO call)
      /*if (_seepage_const>0) {seepage=_seepage_const*(_stage-_local_GW_head); //[m3/s]}*/

      i=0;
      col_ind[i]=GetDVColumnInd(DV_QOUT   ,_aSBIndices[p]);//Qin
      row_val[i]=-0.5;
      i++;

      col_ind[i]=GetDVColumnInd(DV_QOUTRES,_aResIndices[p]);//Qout
      row_val[i]=+0.5;
      i++;

      col_ind[i]=GetDVColumnInd(DV_STAGE  ,_aResIndices[p]);
      row_val[i]=Adt;
      i++;

      for (int dd = 0; dd<pSB->GetReservoir()->GetNumWaterDemands();dd++){
         col_ind[i]=GetDVColumnInd(DV_DELIVERY  ,d);
         row_val[i]=+1.0;
         i++; d++;
      }

      RHS=(precip-ET)-seepage+Adt*h_old-0.5*Qout_last+0.5*Qin_last;

      retval = lp_lib::add_constraintex(pLinProg,i,row_val,col_ind,ROWTYPE_EQ,RHS);
      ExitGracefullyIf(retval==0,"SolveDemandProblem::Error adding mass balance constraint",RUNTIME_ERR);

      // Second goal equation (relation between Qout and h):
      // may be overridden by other management expressions with larger penalties
      // ------------------------------------------------------------------------
      // linearized storage -discharge curve
      //   Qout = Qou_guess + dQ/dh * (h-h_guess)
      //   1.0 * Qout - dQ/dh * h + s1 -s2 = Q_guess - dQ/dh * h_guess

      double h_guess=h_old;
      double Q_guess;

      Q_guess=pRes->GetDischargeFromStage(h_guess,nn);
      dQdh =pRes->GetStageDischargeDerivative(h_guess,nn);

      i = 0;
      col_ind[i]=GetDVColumnInd(DV_QOUTRES,_aResIndices[p]);
      row_val[i]=1.0;
      i++;

      col_ind[i]=GetDVColumnInd(DV_STAGE  ,_aResIndices[p]);
      row_val[i]=-dQdh;
      i++;

      col_ind[i] = GetDVColumnInd(DV_SLACK, s);
      row_val[i]=-1.0;
      i++; s++;

      col_ind[i]=GetDVColumnInd(DV_SLACK  ,s);
      row_val[i]=+1.0;
      i++; s++;

      if (fabs(dQdh)>REAL_SMALL) {
        RHS=Q_guess-dQdh*h_guess;
      }
      else{ //below weir height - Q=0.0
        RHS=0.0;
      }

      retval = lp_lib::add_constraintex(pLinProg,i,row_val,col_ind,ROWTYPE_EQ,RHS);
      ExitGracefullyIf(retval==0,"SolveDemandProblem::Error adding mass balance constraint",RUNTIME_ERR);

      lprow[p]=lp_lib::get_Nrows(pLinProg);

      res_count++;
    }
  }

  // Environmental Minimum Flow Goals
  // ----------------------------------------------------------------
  d=0;
  double minQ,ufp;
  const double ENV_FLOW_PENALTY=1e6;

  for (int pp = 0; pp<pModel->GetNumSubBasins(); pp++)
  {
    p  =pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);

    if (pSB->IsEnabled())
    {
      minQ=pSB->GetEnviroMinFlow(t);
      if (minQ > REAL_SMALL)
      {
        // Q - e+ <= Qmin
        i=0;
        col_ind[i]=GetDVColumnInd(DV_QOUT,_aSBIndices[p]);
        row_val[i]=1.0;
        i++;

        col_ind[i]=GetDVColumnInd(DV_SLACK,s);
        row_val[i]=-1.0;
        i++;

        RHS=minQ;

        retval = lp_lib::add_constraintex(pLinProg,i,row_val,col_ind,ROWTYPE_LE,RHS);
        ExitGracefullyIf(retval==0,"SolveDemandProblem::Error adding environmental flow goal",RUNTIME_ERR);

        // 10^6 * e+ >= (all upstream D)

        i = 0;
        col_ind[i]=GetDVColumnInd(DV_SLACK,s);
        row_val[i]=ENV_FLOW_PENALTY;
        i++;s++;

        for (int dd = 0; dd<_aUpCount[p];dd++)
        {
          col_ind[i]=GetDVColumnInd(DV_DELIVERY,_aUpstreamDemands[p][dd]);
          row_val[i]=-1.0;
          i++;
        }

        RHS=0.0;

        retval = lp_lib::add_constraintex(pLinProg,i,row_val,col_ind,ROWTYPE_GE,RHS);
        ExitGracefullyIf(retval==0,"SolveDemandProblem::Error adding environmental flow goal (2)",RUNTIME_ERR);
      }

      //TODO - figure out how to extend this to an unusable flow percentage
      // Q + e* >= (ufp)*[sum(D)+Q] <-term in brackets is inflow upstream of demands
      ufp =pSB->GetUnusableFlowPercentage();

      // e* <= 10^6 * all upstream D

    }
  }

  // user-specified constraints / goals
  // ----------------------------------------------------------------
  for (int i = 0; i < _nConstraints; i++)
  {
    if (_pConstraints[i]->conditions_satisfied)
    {
      AddConstraintToLP( i, pLinProg, tt, col_ind, row_val);
    }
  }

  lp_lib::set_add_rowmode(pLinProg, FALSE); //turn off once model constraints and obj function are added

  for (int iter=0; iter<3; iter++)
  {
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

    if (_do_debug_level==2){
      //lp_lib::print_debugdump(pLinProg,"lp_solve_dump.txt");
      lp_lib::print_lp(pLinProg);
    }
    //ExitGracefully("Checking matrix",RUNTIME_ERR); //TMP DEBUG

    retval = lp_lib::solve(pLinProg);
    if (retval!=OPTIMAL){cout<<"code: "<<retval<<endl; }

    ExitGracefullyIf(retval!=OPTIMAL,"SolveDemandProblem: non-optimal solution found. Problem is over-constrained. Remove or adjust management constraints.",RUNTIME_ERR);

    //This will be sum of penalties*violations
    if (_do_debug_level>=1){
      cout<<"Objective value: "<<lp_lib::get_objective(pLinProg)+demand_penalty_sum<<endl;
    }

    // assign decision variables
    // ----------------------------------------------------------------
    int ctyp;
    int nrows=lp_lib::get_Nrows(pLinProg);
    double *soln   =new double [_nDecisionVars];
    double *constr =new double [nrows];
    double *resid  =new double [nrows];
    lp_lib::get_variables  (pLinProg,soln);
    lp_lib::get_constraints(pLinProg,constr); //constraint matrix * soln

    for (int i=0;i<_nDecisionVars;i++){
      _pDecisionVars[i]->value=soln[i];
    }
    for (int j=0;j<nrows;j++){
      resid[j]=constr[j]-lp_lib::get_rh(pLinProg,j+1); // LHS-RHS for each goal/constraint
      ctyp=lp_lib::get_constr_type(pLinProg,j+1);
      if (ctyp==LE){upperswap(resid[j],0.0); }
      if (ctyp==GE){lowerswap(resid[j],0.0); resid[j]*=-1; }
    }

    delete [] soln;
    delete [] constr;
    delete [] resid;

    // grab and store improved estimate of reservoir stage
    // ----------------------------------------------------------------
    double  value;
    dv_type typ;
    for (int i=0;i<_nDecisionVars;i++)
    {
      p     = _pDecisionVars[i]->p_index;
      value = _pDecisionVars[i]->value;
      typ   = _pDecisionVars[i]->dvar_type;

      if (typ == DV_STAGE) {
        h_iter[p]=value;
      }
    }

    // update lp matrix for iterative solution of stage-discharge curve
    // ----------------------------------------------------------------
    for (int pp = 0; pp<pModel->GetNumSubBasins(); pp++)
    {
      p   =pModel->GetOrderedSubBasinIndex(pp);
      pSB =pModel->GetSubBasin(p);
      pRes=pSB->GetReservoir();
      if ((pSB->IsEnabled()) && (pSB->GetReservoir()!=NULL))
      {
        double Q_guess;
        double h_guess=h_iter[p];

        Q_guess=pRes->GetDischargeFromStage(h_guess,nn);
        dQdh   =pRes->GetStageDischargeDerivative(h_guess,nn);

        int hcol=GetDVColumnInd(DV_STAGE  ,_aResIndices[p]);

        if (fabs(dQdh)>REAL_SMALL) {
          RHS=Q_guess-dQdh*h_guess;
        }
        else{ //below weir height
          RHS=0.0;
        }
        lp_lib::set_mat(pLinProg,lprow[p],hcol,-dQdh);
        lp_lib::set_rh (pLinProg,lprow[p],RHS);
      }
    }

    //cout<<"iteration: "<<iter<<endl;

  }/*end iteration loop*/

  lp_lib::delete_lp(pLinProg);
  delete [] col_ind;
  delete [] row_val;
  delete [] h_iter;
  delete [] lprow;


  //set state variables corresponding to decision variables
  // SHOULD REALLY ONLY UPDATE DELIVERY, RES_EXTRACT, AND RESERVOIR OUTFLOW - ADD SWITCH TO HAVE OPT MODEL GENERATE STATE VARS OR REGULAR RAVEN SETUP
  // ----------------------------------------------------------------
  double  value;
  dv_type typ;
  int     ii;

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
    else if (typ == DV_QOUTRES)
    {
      pModel->GetSubBasin(p)->GetReservoir()->SetOptimizedOutflow(value); //TMP DEBUG this should ONLY be in reservoirs where one or more management constraints are applied.
    }
    else if (typ == DV_DELIVERY)
    {
      d =_pDecisionVars[i]->loc_index;
      ii=_pDecisionVars[i]->dem_index;

      _aDelivery   [d] =value;
      _aCumDelivery[d]+=value*(Options.timestep*SEC_PER_DAY);

      if (tt.julian_day==_aCumDelDate[d]){_aCumDelivery[d]=0.0; }

      if (ii<100){//not a reservoir
        pModel->GetSubBasin(p)->AddToDeliveredDemand(ii,value);
      }
      else {
        pModel->GetSubBasin(p)->GetReservoir()->AddToDeliveredDemand(ii-100, value);
      }
    }
    else if (typ == DV_SLACK)
    {
      s = _pDecisionVars[i]->loc_index;
      _aSlackValues[s]=value;
    }
    else if (typ == DV_USER)
    {
      //doesn't influence hydrologic model
    }
  }


#endif
}

//////////////////////////////////////////////////////////////////
/// \brief Write output file headers
/// \details Called prior to simulation (but after initialization) from CModel::Initialize()
/// \param &Options [in] Global model options information
//
void CDemandOptimizer::WriteOutputFileHeaders(const optStruct &Options)
{
  CSubBasin *pSB;

  // DemandOptimization.csv
  //--------------------------------------------------------
  string tmpFilename=FilenamePrepare("DemandOptimization.csv",Options);
  _DEMANDOPT.open(tmpFilename.c_str());
  if (_DEMANDOPT.fail()){
    ExitGracefully(("CModel::WriteOutputFileHeaders: unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
  }

  _DEMANDOPT<<"time [d],date,hour";
  for (int i=0;i<_nDecisionVars;i++){
    if (_pDecisionVars[i]->p_index == DOESNT_EXIST) { //user specified variable
       _DEMANDOPT<<","<<_pDecisionVars[i]->name;
    }
    else{
      pSB=_pModel->GetSubBasin(_pDecisionVars[i]->p_index);
      if (pSB->IsGauged()){
        _DEMANDOPT<<","<<_pDecisionVars[i]->name;
      }
    }
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
  CSubBasin *pSB;
  string usedate=tt.date_string;                   //refers to date and time at END of time step
  string usehour=DecDaysToHours(tt.julian_day);

  if (tt.model_time==0){return;}

  // DemandOptimization.csv
  //--------------------------------------------------------
  _DEMANDOPT<<tt.model_time <<","<<usedate<<","<<usehour;
  for (int i=0;i<_nDecisionVars;i++){
    if (_pDecisionVars[i]->p_index == DOESNT_EXIST) { //user specified variable
       _DEMANDOPT<<","<<_pDecisionVars[i]->value;
    }
    else{
      pSB=_pModel->GetSubBasin(_pDecisionVars[i]->p_index);
      if (pSB->IsGauged()){ //gauged basin or not associated with basin
        _DEMANDOPT<<","<<_pDecisionVars[i]->value;
      }
    }
  }
  //_DEMANDOPT<<","<<_aDelivery[0];
  _DEMANDOPT<<endl;

  // GoalSatisfaction.csv
  //--------------------------------------------------------
  // All slack values reported should be positive
  // each entry represents the magnitude of a goal violation, expressed in m3/s
  // for environmental flow goals, for instance, this is E-Q if Q<E, 0 otherwise
  int p;
  bool include_pen=false;
  double pen;
  int d=0;
  double mult=1.0;
  _GOALSAT<<tt.model_time <<","<<usedate<<","<<usehour;

  for (int pp=0;pp<_pModel->GetNumSubBasins();pp++){
    p = _pModel->GetOrderedSubBasinIndex(pp);
    pSB=_pModel->GetSubBasin(p);
    for (int ii = 0; ii < pSB->GetNumWaterDemands(); ii++)
    {
      if (include_pen){mult=_aDemandPenalties[d];}
      _GOALSAT << "," << mult*(pSB->GetWaterDemand(ii,tt.model_time)-_aDelivery[d]);
      d++;
    }
  }
  int k=_nDemands; //first goals are due to environmental flow constraints
  for (int i = 0; i < _nConstraints; i++) {
    if (_pConstraints[i]->is_goal) {
      if (_pConstraints[i]->pExpression->compare != COMPARE_IS_EQUAL) {
        if (include_pen){mult=_pConstraints[i]->penalty_over;}
        pen = mult*_aSlackValues[k];
        _GOALSAT<<","<<pen;
        k++;
      }
      else {
        if (include_pen){mult=_pConstraints[i]->penalty_over;}
        pen =mult*_aSlackValues[k];
        if (include_pen){mult=_pConstraints[i]->penalty_under;}
        pen+=mult*_aSlackValues[k+1];
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
