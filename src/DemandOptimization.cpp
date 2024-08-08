/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2024 the Raven Development Team

  Uses lpsolve:
    https://lpsolve.sourceforge.net/5.5/
    Files from lp_solve_5.5.2.11_dev_win64.zip distribution
    under GNU Public License
  ----------------------------------------------------------------*/
#include "DemandOptimization.h"

void SummarizeExpression(const char **s, const int Len, expressionStruct* exp); //defined in DemandExpressionHandling.cpp
string DVTypeToString(dv_type t);

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the Demand optimization constructor
/// \details Creates empty instance of the demand optimization class
//
CDemandOptimizer::CDemandOptimizer(CModel *pMod)
{
  _pModel=pMod;

  _nDemands=0;
  _pDemands=NULL;

  _nDemandGroups=0;
  _pDemandGroups=NULL;

  _aDelivery=NULL;
  _aCumDelivery=NULL;
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

  _nControlVars=0;
  _pControlVars=NULL;

  _nUserTimeSeries=0;
  _pUserTimeSeries=NULL;

  _nUserLookupTables=0;
  _pUserLookupTables=NULL;

  _aUpstreamDemands=NULL;
  _aUpCount=NULL;

  _nHistoryItems=1;
  _aQhist=NULL;
  _aDhist=NULL;
  _ahhist=NULL;

  _nGoals=0;
  _pGoals=NULL;

  _demands_initialized=false;

  _do_debug_level=0;//no debugging

  _nSolverResiduals=0;
  _aSolverResiduals=NULL;
  _aSolverRowNames=NULL;

  _stage_discharge_as_goal=false; //TMP DEBUG - SHOULD BE FALSE ; should remove member variable once we know this works in general case.
}
//////////////////////////////////////////////////////////////////
/// \brief Implementation of the Demand optimization destructor
//
CDemandOptimizer::~CDemandOptimizer()
{
  if (DESTRUCTOR_DEBUG){
    cout<<"DESTROYING DEMAND OPTIMIZER"<<endl;
  }
  for (int i=0;i<_nDecisionVars;    i++){delete _pDecisionVars[i];    }delete [] _pDecisionVars;
  for (int i=0;i<_nUserTimeSeries;  i++){delete _pUserTimeSeries[i];  }delete [] _pUserTimeSeries;
  for (int i=0;i<_nUserLookupTables;i++){delete _pUserLookupTables[i];}delete [] _pUserLookupTables;
  for (int i=0;i<_nControlVars;     i++){delete _pControlVars[i];     }delete [] _pControlVars;
  for (int i=0;i<_nDemandGroups;    i++){delete _pDemandGroups[i];    }delete [] _pDemandGroups;
  for (int i=0;i<_nDemands;         i++){delete _pDemands[i];         }delete [] _pDemands;

  for (int p=0;p<_pModel->GetNumSubBasins(); p++){delete [] _aUpstreamDemands[p]; } delete [] _aUpstreamDemands;
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
  delete [] _aCumDelivery;
  delete [] _aSlackValues;
  delete [] _aSBIndices;
  delete [] _aResIndices;
  delete [] _aSolverResiduals;
  delete [] _aSolverRowNames;
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
    if      (s_to_i(demand_tag.c_str()) == _pDemands[d]->GetID()  ){return d;}
    else if (demand_tag                 == _pDemands[d]->GetName()){return d;}
  }
  return DOESNT_EXIST;
}
//////////////////////////////////////////////////////////////////
/// \brief gets demand from index 
//
CDemand* CDemandOptimizer::GetWaterDemand(const int d) {
  return _pDemands[d];
}
//////////////////////////////////////////////////////////////////
/// \brief gets number of enabled water demands
//
int CDemandOptimizer::GetNumWaterDemands() const {
  return _nDemands;
}
//////////////////////////////////////////////////////////////////
/// \brief gets demand group from index 
//
CDemandGroup* CDemandOptimizer::GetDemandGroup(const int ii) 
{
  return _pDemandGroups[ii];
}
//////////////////////////////////////////////////////////////////
/// \brief gets demand group from name, or NULL if invalid name 
//
CDemandGroup* CDemandOptimizer::GetDemandGroupFromName(const string name) 
{
  for (int ii=0;ii<_nDemandGroups;ii++){
    if (_pDemandGroups[ii]->GetName() == name) {
      return _pDemandGroups[ii];
    }
  }
  return NULL;
}
//////////////////////////////////////////////////////////////////
/// \brief gets demand group from index 
//
int CDemandOptimizer::GetNumDemandGroups() const 
{
  return _nDemandGroups;
}
//////////////////////////////////////////////////////////////////
/// \brief returns true if demands have been initialized 
//
bool CDemandOptimizer::DemandsAreInitialized() const
{
  return _demands_initialized;
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
/// \brief assigns a demand 'unrestricted' status, meaning it wont be considered when applying environmental min flow constraints
/// \params dname [in] - name of demand
//
void   CDemandOptimizer::SetDemandAsUnrestricted(const string dname) 
{
  int d = GetDemandIndexFromName(dname);
  if (d!=DOESNT_EXIST){
    _pDemands[d]->SetAsUnrestricted();
  }
  else{
    string warn = "CDemandOptimizer::SetDemandAsUnrestricted: invalid or disabled demand name " + dname +" provided in : DemandIsUnrestricted command. This will be ignored ";
    WriteWarning(warn,true); 
  }
}
//////////////////////////////////////////////////////////////////
/// \brief sets demand penalty for demand dname
//
void CDemandOptimizer::SetDemandPenalty(const string dname, const double& pen)
{
  ExitGracefullyIf(pen<0,"CDemandOptimizer::SetDemandPenalty: demand penalties must be greater or equal to zero",BAD_DATA_WARN);
  int d=GetDemandIndexFromName(dname);
  if (d!=DOESNT_EXIST){
    _pDemands[d]->SetDemandPenalty(pen);
  }
  else {
    string warn = "CDemandOptimizer::SetDemandPenalty: invalid or disabled demand name " + dname +" provided in : DemandPenalty command. This will be ignored ";
    WriteWarning(warn,true); 
  }
}
//////////////////////////////////////////////////////////////////
/// \brief date at which cumulative delivery is rebooted to zero (usually winter date)
/// \params julian_date [in] reboot date
/// \params demandID [in] demand identifier (e.g., 'D123')
/// only called from ParseManagementFile()
//
void CDemandOptimizer::SetCumulativeDate(const int julian_date, const string demandID)
{
  int d=GetDemandIndexFromName(demandID);
  if (d!=DOESNT_EXIST){
    _pDemands[d]->SetCumulDeliveryDate(julian_date);
  }
  else{
    WriteWarning("DemandOptimization::SetCumulativeDate: invalid or disabled demand ID in :DemandResetDate command",true);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief adds demand group to array of demand groups
/// \params groupname [in] - user-specified demand group name
//
void   CDemandOptimizer::AddDemandGroup(const string groupname) 
{
  CDemandGroup *pDG=new CDemandGroup(groupname,_nDemandGroups);

  if (!DynArrayAppend((void**&)(_pDemandGroups),(void*)(pDG),_nDemandGroups)){
   ExitGracefully("CDemandOptimizer::AddDemandGroup: adding NULL DG",BAD_DATA);}
}

//////////////////////////////////////////////////////////////////
/// \brief adds demand ONLY if subbasin associated with demand is enabled, otherwise tosses out.
/// \params pDem [in] - pointer to demand object
//
void  CDemandOptimizer::AddWaterDemand(CDemand* pDem) 
{
  long SBID=pDem->GetSubBasinID();
  if (_pModel->GetSubBasinByID(SBID)->IsEnabled())
  {
    if (!DynArrayAppend((void**&)(_pDemands),(void*)(pDem),_nDemands)){
     ExitGracefully("CDemandOptimizer::AddDemandGroup: adding NULL DG",BAD_DATA);}
  }
  else{
    delete pDem;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief adds decision variable structure to _pDecisionVars member aarray
/// \params pDV [in] - user-specified decision variable to be added
//
void CDemandOptimizer::AddDecisionVar(const decision_var* pDV)
{
  if (VariableNameExists(pDV->name)) {
    string warn="CDemandOptimizer::AddDecisionVar: decision variable name "+pDV->name+" is already in use.";
    ExitGracefully(warn.c_str(),BAD_DATA_WARN);
  }

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
  if (VariableNameExists(name)) {
    string warn="CDemandOptimizer::AddUserConstant: control variable name "+name+" is already in use.";
    ExitGracefully(warn.c_str(),BAD_DATA_WARN);
  }

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
/// \brief adds control variable
//
void   CDemandOptimizer::AddControlVariable(const string name, expressionStruct* pExp)
{
  if (VariableNameExists(name)) {
    string warn="CDemandOptimizer::AddControlVariable: control variable name "+name+" is already in use.";
    ExitGracefully(warn.c_str(),BAD_DATA_WARN);
  }

  control_var *pCV = new control_var();
  pCV->name=name;
  pCV->pExpression=pExp;
  pCV->current_val=RAV_BLANK_DATA;

  if (!DynArrayAppend((void**&)(_pControlVars),(void*)(pCV),_nControlVars)){
   ExitGracefully("CDemandOptimizer::AddControlVariable: adding NULL CV",BAD_DATA);}
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
/// \brief adds decision variable constraint OR goal to _pGoals member array
/// \param pGoal [in] - goal/constraint structure
///
void CDemandOptimizer::AddGoalOrConstraint(const managementGoal *pGoal)
{
  if (!DynArrayAppend((void**&)(_pGoals),(void*)(pGoal),_nGoals)){
   ExitGracefully("CDemandOptimizer::AddConstraint: adding NULL goal or constraint",BAD_DATA);}
}

//////////////////////////////////////////////////////////////////
/// \brief returns true if variable name is already in use
//
bool CDemandOptimizer::VariableNameExists(const string &name) const
{
  for (int i = 0; i < _nControlVars; i++) {
    if (_pControlVars[i]->name==name){return true;}
  }
  for (int i=0; i<_nDecisionVars; i++){
    if (_pDecisionVars[i]->name==name){return true;}
  }
  for (int i = 0; i < _nUserConstants; i++) {
    if (_aUserConstNames[i]==name){return true;}
  }
  return false;
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

  // Catalogue enabled basins/reservoirs - populate _aSBindices, _aResIndices
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

  // Populate ordered decision variable array _pDecisionVars[]
  //  This order has to be maintained because of how the decision variables are indexed (consistent with GetDVColumnInd())
  //  subbasin outflows -> reservoir outflows -> reservoir stages -> delivered demand -> user-specified DVs -> slack variables
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
  // add reservoir binary DVs
  if (!_stage_discharge_as_goal){
  res_count=0;
  for (int pp=0;pp<pModel->GetNumSubBasins();pp++)
  {
    p=pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    if (pSB->IsEnabled()){
      if (pSB->GetReservoir()!=NULL){
        name="B"+to_string(pModel->GetSubBasin(p)->GetID());
        pDV=new decision_var(name,p,DV_BINRES,res_count);
        AddDecisionVar(pDV);
        res_count++;
      }
    }
  }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief populates _aUpstreamDemands and _aUpCount arrays
/// called by CDemandOptimizer::InitializeDemands()
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
    if (!_pDemands[d]->IsUnrestricted()){
      SBID=_pDemands[d]->GetSubBasinID();
      while (SBID >= 0) {
        int p=_pModel->GetSubBasinByID(SBID)->GetGlobalIndex();

        if (_pModel->GetSubBasin(p)->IsEnabled()){
          pushIntoIntArray(_aUpstreamDemands[p],d,_aUpCount[p]);
          SBID=_pModel->GetSubBasin(p)->GetDownstreamID();
        }
        else {
          SBID=DOESNT_EXIST;
        }
      }
    }
  }

  /*cout << " ** ** ** ** ** ** ** ** ** ** ** " << endl;
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
  cout<<" ** ** ** ** ** ** ** ** ** ** ** "<<endl;*/
}
//////////////////////////////////////////////////////////////////
/// \brief Initializes Demand decision variables
/// \notes to be called during .rvm file read, PRIOR to declaring any user-specified DVs, goals, or constraints
/// \params pModel [in] - pointer to model
/// \params Options  [in] - model options structure
//
void CDemandOptimizer::InitializeDemands(CModel* pModel, const optStruct& Options)
{
  if (_demands_initialized){return;}//This routine has already been called

  if (Options.noisy){cout<<"CDemandOptimizer: Demand initialization..."<<endl;}

  // reserve memory for delivery arrays 
  //------------------------------------------------------------------
  _aDelivery        = new double[_nDemands];
  _aCumDelivery     = new double[_nDemands];
  for (int d=0;d<_nDemands;d++)
  {
    _aDelivery       [d]=0.0;
    _aCumDelivery    [d]=0.0;
  }

  // determine upstream demand connectivity
  //------------------------------------------------------------------
  IdentifyUpstreamDemands(); 

  int p;
  decision_var *pDV;

  // add delivered demand decision vars 
  //------------------------------------------------------------------
  for (int d = 0; d < _nDemands; d++) 
  {
    p=pModel->GetSubBasinByID(_pDemands[d]->GetSubBasinID())->GetGlobalIndex();
    pDV=new decision_var(_pDemands[d]->GetName(), p, DV_DELIVERY, d);
    pDV->dem_index=_pDemands[d]->GetLocalIndex();

    AddDecisionVar(pDV);
  }

  _demands_initialized=true;
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
  CSubBasin    *pSB;
  decision_var *pDV;

  // create 2 slack variables for reservoir outflow by stage-discharge
  //------------------------------------------------------------------
  for (int pp=0;pp<pModel->GetNumSubBasins();pp++)
  {
    p=pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    if ((pSB->IsEnabled()) && (pSB->GetReservoir()!=NULL)){
      pDV =  new decision_var("SR+" + to_string(pModel->GetSubBasin(p)->GetID()), p, DV_SLACK,_nSlackVars);
      AddDecisionVar(pDV);
      _nSlackVars++;
      pDV =  new decision_var("SR-" + to_string(pModel->GetSubBasin(p)->GetID()), p, DV_SLACK,_nSlackVars);
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
        pDV =  new decision_var("E+" + to_string(pModel->GetSubBasin(p)->GetID()), p, DV_SLACK,_nSlackVars);
        AddDecisionVar(pDV);
        _nSlackVars++;
      }
      /*if (pSB->GetUnusableFlowPercentage()>0) {
        pDV =  new decision_var("U+" + to_string(pModel->GetSubBasin(p)->GetID()), p, DV_SLACK,_nSlackVars);
        AddDecisionVar(pDV);
        _nSlackVars++;
      }*/
    }
  }

  // Convert reservoir commands to user-specified constraints
  //------------------------------------------------------------------
  AddReservoirConstraints();

  // Calculate penalty units correction for reservoirs 
  //------------------------------------------------------------------
  for (int j = 0; j < _nGoals; j++)
  {
    if ((_pGoals[j]->is_goal) &&  (_pGoals[j]->use_stage_units)) 
    {
      pSB =_pModel->GetSubBasin(_pGoals[j]->reservoir_index);

      if (pSB->GetReservoir()!=NULL)
      {
        int k=pSB->GetReservoir()->GetHRUIndex();
        if (k!=DOESNT_EXIST){
          _pGoals[j]->units_correction=(_pModel->GetHydroUnit(k)->GetArea()*M2_PER_KM2)/SEC_PER_DAY; //[m]/[dt (d)]->[m3/s];
        }
      }
      else {
        ExitGracefully("InitalizePostRVMRead: invalid reservoir index associated with goal",RUNTIME_ERR);
      }
    }
  }



  // set operating regime penalties to default, if not specified
  //------------------------------------------------------------------
  for (int j = 0; j < _nGoals; j++)
  {
    for (int k=0;k<_pGoals[j]->nOperRegimes;k++){
      if (_pGoals[j]->pOperRegimes[k]->penalty_over ==RAV_BLANK_DATA){_pGoals[j]->pOperRegimes[k]->penalty_over =_pGoals[j]->penalty_over; } 
      if (_pGoals[j]->pOperRegimes[k]->penalty_under==RAV_BLANK_DATA){_pGoals[j]->pOperRegimes[k]->penalty_under=_pGoals[j]->penalty_under;}
    }
  }

  // Sift through user-specified goals, update _nSlackVars, reserve memory for _aSlackValues
  //------------------------------------------------------------------
  // ISSUE : THIS IS TOTAL POSSIBLE NUMBER OF SLACK VARIABLES; MAY NOT BE IN PLAY DUE TO CONDITIONALS - for disabled constraints, just set penalty to zero>??
  //
  string nm;

  for (int j = 0; j < _nGoals; j++)
  {
    // ASSUMES ALL EXPRESSIONS IN GOAL/CONSTRAINT ARE EITHER == or >/<, NEVER MIXED.
    if (_pGoals[j]->is_goal) {
      if (_pGoals[j]->pOperRegimes[0]->pExpression->compare == COMPARE_IS_EQUAL) {
        _pGoals[j]->slack_ind1=_nSlackVars;
        _pGoals[j]->slack_ind2=_nSlackVars+1;
        pDV =  new decision_var("SL+" + to_string(j), p, DV_SLACK,_nSlackVars);
        AddDecisionVar(pDV);
        _nSlackVars++;
        pDV =  new decision_var("SL-" + to_string(j), p, DV_SLACK,_nSlackVars);
        AddDecisionVar(pDV);
        _nSlackVars++;
      }
      else {
        _pGoals[j]->slack_ind1=_nSlackVars;
        if (_pGoals[j]->pOperRegimes[0]->pExpression->compare==COMPARE_LESSTHAN){nm="-";}
        else                                                                    {nm="+";}
        pDV =  new decision_var("SL"+ nm + to_string(j), p, DV_SLACK,_nSlackVars);
        AddDecisionVar(pDV);
        _nSlackVars++;
      }
    }
    _pGoals[j]->ever_satisfied=false;
  }

  _aSlackValues = new double[_nSlackVars];
  for (int s = 0; s < _nSlackVars; s++) {
    _aSlackValues[s]=0.0;
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

  //TODO: Need to evaluate ALL conditional statements in operating regimes so that issues may be flagged PRIOR to simulation

  // Print summary to screen
  //------------------------------------------------------------------
  if (true)
  {
    cout<<" -----------------------------------------------------------------"<<endl;
    cout<<"  MANAGEMENT OPTIMIZATION SUMMARY                                 "<<endl;
    cout<<" -----------------------------------------------------------------"<<endl;
    cout<<" # Enabled subbasins: "<<_nEnabledSubBasins<<" of "<<pModel->GetNumSubBasins()<<endl;
    cout<<" # Enabled reservoirs/lakes: "<<_nReservoirs<<endl;
    cout<<" # Decision Vars:      "<<_nDecisionVars    <<endl;
    cout<<" # User decision vars: "<<_nUserDecisionVars<<endl;
    cout << "    DV " << setw(4) << "#" << ": " << setw(20) << "name" << " " << setw(12) << "type" << " " << setw(4) << "col" << " " << "loc_index"<< endl;
    for (int i = 0; i < _nDecisionVars; i++)
    {
      cout<<"    DV "<<setw(4)<<i<<": "<<setw(20)<<_pDecisionVars[i]->name<<" "<<setw(12)<<DVTypeToString(_pDecisionVars[i]->dvar_type)<<" ";
      cout<<setw(4)<<GetDVColumnInd(_pDecisionVars[i]->dvar_type,_pDecisionVars[i]->loc_index)<<" "<<_pDecisionVars[i]->loc_index<<endl;
    }

    cout<<" # Demands : "<<_nDemands<<endl;
    for (int d=0; d<_nDemands;d++){
      cout << "    " << d << ": ID="<<setw(6) << _pDemands[d]->GetID() << " (alias: "<<setw(32)<<_pDemands[d]->GetName()<<") in basin "<< _pDemands[d]->GetSubBasinID() <<endl;
    }

    string tmpstr,tmpstr2;
    cout<<" # History Items: "     <<_nHistoryItems<<endl;
    cout<<" # Slack Vars: "        <<_nSlackVars<<endl;
    cout<<" # Constraints/Goals: " <<_nGoals<<endl;
    for (int i = 0; i < _nGoals; i++) {
      if (_pGoals[i]->is_goal){tmpstr="[GOAL]      "; }
      else                          {tmpstr="[CONSTRAINT]"; }
      cout<<"    "<<i<<" "<<tmpstr<<": "<<_pGoals[i]->name<<endl;
      for (int k=0; k<_pGoals[i]->nOperRegimes; k++)
      {
        comparison ctype=_pGoals[i]->pOperRegimes[k]->pExpression->compare;
        cout<<"      +oper regime: "<<_pGoals[i]->pOperRegimes[k]->reg_name<<endl;
        cout<<"        +expression: "<<_pGoals[i]->pOperRegimes[k]->pExpression->origexp<<endl;
        cout<<"        +comp. type: "<<ctype<<endl;
        if (ctype==COMPARE_IS_EQUAL){
        cout<<"        +penalty:    ("<<_pGoals[i]->pOperRegimes[k]->penalty_under<<","<<_pGoals[i]->pOperRegimes[k]->penalty_over<<")"<<endl;
        }
        else {
        cout<<"        +penalty:    ("<<_pGoals[i]->pOperRegimes[k]->penalty_under<<")"<<endl;
        }
        for (int j = 0; j < _pGoals[i]->pOperRegimes[k]->nConditions; j++) {
          if (_pGoals[i]->pOperRegimes[k]->pConditions[j]->pExp!=NULL){tmpstr2=_pGoals[i]->pOperRegimes[k]->pConditions[j]->pExp->origexp;}
          else                                                              {tmpstr2=_pGoals[i]->pOperRegimes[k]->pConditions[j]->dv_name;}
          cout<<"        + condition: "<<tmpstr2<<endl;
        }
      }

    }
    cout<<" -----------------------------------------------------------------"<<endl;
    cout<<" -----------------------------------------------------------------"<<endl;
  }
  if (Options.noisy){cout<<"   ...end Post-rvm-read initialization."<<endl;}
}


//////////////////////////////////////////////////////////////////
// converts string instring to tokenized array s (preallocated) with Len tokens
// used within AddReservoirConstraints to emulate Parser
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
// initialization step: called from InitializePostRVMFileRead()
//
void CDemandOptimizer::AddReservoirConstraints()
{
  int               p;
  CSubBasin        *pSB;
  string            TSname,SBIDs;
  long              SBID;
  expressionStruct *exp;
  managementGoal   *pGoal=NULL;
  string            expString;
  string            advice;
  int               Len=0;
  char             *s[MAXINPUTITEMS];

  for (int pp=0;pp<_pModel->GetNumSubBasins();pp++)
  {
    p   =_pModel->GetOrderedSubBasinIndex(pp);
    pSB =_pModel->GetSubBasin(p);
    SBID=pSB->GetID();
    SBIDs=to_string(SBID);

    if (((pSB->IsEnabled()) && (pSB->GetReservoir()!=NULL)))
    {
      int k=pSB->GetReservoir()->GetHRUIndex();
      if (k==DOESNT_EXIST){
        advice="AddReservoirConstraints: The reservoir in subbasin "+SBIDs+" doesnt have an :HRUID - this will negatively impact the units scaling of management optimization penalties.";
        ExitGracefully(advice.c_str(), BAD_DATA_WARN);
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
        pGoal=new managementGoal();
        pGoal->name=TSname;
        pGoal->is_goal=true;
        pGoal->AddExpression(exp);
        pGoal->use_stage_units=true;
        pGoal->reservoir_index=p;
        pGoal->penalty_over=6.0;
        
        AddGoalOrConstraint(pGoal);
        advice="A :ReservoirMaxStage time series for subbasin "+SBIDs+" has been added to the management optimization formulation";
        WriteAdvisory(advice,false);
      }

      //Min Stage constraints
      // Q = Q_min (or zero) if h < h_min
      //------------------------------------------------------------------------
      TSname="_MinStage_" + SBIDs;
      string TSname2="_MinStageFlow_" + SBIDs;
      if (UserTimeSeriesExists(TSname))
      {
        if (UserTimeSeriesExists(TSname2)) {
          expString=":Expression !Q" + SBIDs+ " = @ts(" + TSname2 + ",0)";
        }
        else {
          expString=":Expression !Q" + SBIDs+ " = 0.0";
        }

        TokenizeString(expString, s, Len);
        exp = ParseExpression((const char**)(s), Len, 0, "internal");
        
        pGoal=new managementGoal();
        pGoal->name=TSname;
        pGoal->is_goal=true;
        pGoal->AddExpression(exp);
        pGoal->use_stage_units=true;
        pGoal->reservoir_index=p;
        pGoal->penalty_over=6.0;
        AddGoalOrConstraint(pGoal);

        exp_condition *pCond = new exp_condition();

        expString = ":Condition !h" + SBIDs + "[-1] < @ts(" + TSname + ",0)";
        TokenizeString(expString, s, Len);
        pCond->pExp=ParseExpression((const char**)(s), Len, 0, "internal");
        pGoal->AddOpCondition(pCond);

        if (GetDebugLevel()>=1){
          SummarizeExpression((const char**)(s),Len,pCond->pExp);
        }

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

        pGoal=new managementGoal();
        pGoal->name=TSname;
        pGoal->is_goal=true;
        pGoal->AddExpression(exp);
        pGoal->penalty_over=5.0;
        AddGoalOrConstraint(pGoal);

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

        pGoal=new managementGoal();
        pGoal->name=TSname;
        pGoal->is_goal=true;
        pGoal->AddExpression(exp);
        pGoal->penalty_over=5.0;
        AddGoalOrConstraint(pGoal);

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

        pGoal=new managementGoal();
        pGoal->name=TSname;
        pGoal->is_goal=true;
        pGoal->AddExpression(exp);
        pGoal->penalty_over=4.0;
        AddGoalOrConstraint(pGoal);

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

        pGoal=new managementGoal();
        pGoal->name=TSname;
        pGoal->is_goal=true;
        pGoal->AddExpression(exp);
        pGoal->penalty_over=3.0;
        AddGoalOrConstraint(pGoal);

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

        pGoal=new managementGoal();
        pGoal->name=TSname;
        pGoal->is_goal=true;
        pGoal->AddExpression(exp);
        pGoal->penalty_over=3.0;
        AddGoalOrConstraint(pGoal);

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
        TokenizeString(expString, s, Len);
        exp = ParseExpression((const char**)(s), Len, 0, "internal");

        pGoal=new managementGoal();
        pGoal->name=TSname;
        pGoal->is_goal=true;
        pGoal->AddExpression(exp);
        pGoal->use_stage_units=true;
        pGoal->reservoir_index=p;
        pGoal->penalty_over=2.0;
        AddGoalOrConstraint(pGoal);

        advice = "An :ReservoirTargetStage time series for subbasin " + SBIDs + " has been added to the management optimization formulation";
        WriteAdvisory(advice,false);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Updates control variables
/// called every time step by SolveDemandProblem() prior to solve
//
void CDemandOptimizer::UpdateControlVariables(const time_struct &tt)
{
  double t=tt.model_time;
  for (int i = 0; i < _nControlVars; i++) {
    _pControlVars[i]->current_val=EvaluateConditionExp(_pControlVars[i]->pExpression, t);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Updates history arrays
/// called every time step by SolveDemandProblem() prior to solve
//
void CDemandOptimizer::UpdateHistoryArrays()
{
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
#ifdef _LPSOLVE_
void CDemandOptimizer::IncrementAndSetRowName(lp_lib::lprec *pLinProg,int &rowcount,const string &name)
{
    rowcount++;
    string rowname=name;
    lp_lib::set_row_name(pLinProg, rowcount, &rowname[0]);
}
#endif
//////////////////////////////////////////////////////////////////
// indexing for DV vector
// 0          .. nSB-1 - outflow of stream reach p
// nSB        .. nSB+1*nRes-1 - reservoir outflows
// nSB+nRes   .. nSB+2*nRes-1 - reservoir stages
// nSB+2*nRes .. nSB+3*nRes-1 - reservoir binary variables
// nSB+2*nRes .. nSB+2*nRes+nDem-1 - satisfied demand
// nSB+2*nRes+nDem .. nSB+2*nRes+nDem+nDV - user-specified decision variables
// remaining: slack vars - [#enviro min flow + # user <=/>= goals + 2* # user == goals + all reservoir goals ]
//
int CDemandOptimizer::GetDVColumnInd(const dv_type typ, const int counter) const
{
  int N=2;
  if (!_stage_discharge_as_goal){N=3;}

  if      (typ==DV_QOUT    ){return counter+1;}
  else if (typ==DV_QOUTRES ){return _nEnabledSubBasins+counter+1;}
  else if (typ==DV_STAGE   ){return _nEnabledSubBasins+_nReservoirs+counter+1;}
  else if (typ==DV_BINRES  ){return _nEnabledSubBasins+2*_nReservoirs+counter+1;}
  else if (typ==DV_DELIVERY){return _nEnabledSubBasins+N*_nReservoirs+counter+1;}
  else if (typ==DV_USER    ){return _nEnabledSubBasins+N*_nReservoirs+_nDemands+counter+1; }
  else if (typ==DV_SLACK   ){return _nEnabledSubBasins+N*_nReservoirs+_nDemands+_nUserDecisionVars+counter+1;}

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

  int    d,i,p,s;
  double RHS;
  int    retval;
  double demand_penalty_sum=0; //sum (P * Dstar) to be added to global obj function such that min is zero
  double tstep=Options.timestep;
  string rowname;
  int    rowcount=0;

  CSubBasin *pSB;

  double t=tt.model_time;
  //int nn=tt.nn+1;//end of timestep
  int nn=static_cast<int>((t+TIME_CORRECTION)/Options.timestep);//+1;//end-of timestep index

  int    *col_ind=new int    [_nDecisionVars]; //index of column to insert value in current row (1:nDV, not zero-indexed)
  double *row_val=new double [_nDecisionVars]; //values of row[col_ind]

  double *h_iter=new double [_pModel->GetNumSubBasins()];
  double *Q_iter=new double [_pModel->GetNumSubBasins()];
  int    *lprow =new int    [_pModel->GetNumSubBasins()]; //index of goal equation for non-linear reservoir stage discharge curve in subbasin p

  // evaluates value of all control variables for this time step
  // ----------------------------------------------------------------
  UpdateControlVariables(tt);

  // update history arrays from previous timestep
  // ----------------------------------------------------------------
  UpdateHistoryArrays();

  // Update demand magnitudes
  // ----------------------------------------------------------------
  for (int d = 0; d < _nDemands; d++) {
    _pDemands[d]->UpdateDemand(Options,tt);
  }
  
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
  for (int d = 0; d < _nDemands; d++) {
    retval=lp_lib::set_upbo(pLinProg,GetDVColumnInd(DV_DELIVERY,d), _pDemands[d]->GetDemand());
    ExitGracefullyIf(retval!=1,"SolveDemandProblem::Error adding demand upper bound",RUNTIME_ERR);
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

  // Set reservoir binary variables to integer type with bounds of [0,1]
  // ----------------------------------------------------------------
  if (!_stage_discharge_as_goal) {
    int binct=0;
    for (int i = 0; i < _nDecisionVars; i++) {
      if (_pDecisionVars[i]->dvar_type==DV_BINRES){
        retval=lp_lib::set_binary(pLinProg,GetDVColumnInd(DV_BINRES,binct),1);
        ExitGracefullyIf(retval!=1,"SolveDemandProblem::Error setting decision variable as binary",RUNTIME_ERR);
        binct++;
      }
    }
  }

  // Prepare flow diversion estimates
  // ----------------------------------------------------------------
  int     pDivert;
  double  div_Q=0;
  double  sum_diverted=0;
  double *aDivert  =new double [pModel->GetNumSubBasins()]; //< diversion TO basin p (from prev timestep)
  double *aDivGuess=new double [pModel->GetNumSubBasins()]; //< diversion FROM basin p (non-linear)
  for (int p = 0; p<pModel->GetNumSubBasins(); p++)
  {
    aDivert  [p]=0;
    aDivGuess[p]=0;
  }
  for (int pp = 0; pp<pModel->GetNumSubBasins(); pp++)
  {
    p  =pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    sum_diverted=0;
    for(int i=0;i<pSB->GetNumDiversions();i++)
    {
      div_Q=pSB->GetDiversionFlow(i,pSB->GetLastOutflowRate(),Options,tt,pDivert);//[m3/s]
      if(pDivert!=DOESNT_EXIST)
      {
        aDivert[pDivert]+=div_Q;
      }
      sum_diverted+=div_Q;
    }
    aDivGuess[p]=sum_diverted; // a guess because GetDiversionFlow may be non-linear
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
  for (int d = 0; d < _nDemands; d++) 
  {
    demand_penalty_sum+=(_pDemands[d]->GetPenalty()*_pDemands[d]->GetDemand());

    col_ind[i]=GetDVColumnInd(DV_DELIVERY,d);
    row_val[i]=-_pDemands[d]->GetPenalty(); //this maximizes delivery
    i++;
  }

  s=0;
  // reservoir discharge curve goal penalties (2*nReservoirs slack vars)
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

  // enviro min flow goal penalties (_nEnviroFlowGoals slack vars)
  // ----------------------------------------------------------------
  _nEnviroFlowGoals=0;
  for (int pp = 0; pp<pModel->GetNumSubBasins(); pp++)
  {
    p  =pModel->GetOrderedSubBasinIndex(pp);
    pSB=pModel->GetSubBasin(p);
    if (pSB->IsEnabled()){
      if (pSB->GetEnviroMinFlow(t)>REAL_SMALL)  {
        col_ind[i]=GetDVColumnInd(DV_SLACK,s);
        row_val[i]=1.0;
        i++; s++;
        _nEnviroFlowGoals++;
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
  bool op_is_active;
  double penalty_over,penalty_under;

  for (int j = 0; j < _nGoals; j++)
  {
    _pGoals[j]->conditions_satisfied=false;
    _pGoals[j]->active_regime=DOESNT_EXIST;
    
    // ASSUMES ALL OP REGIME EXPRESSIONS IN GOAL/CONSTRAINT ARE EITHER == or >/<, NEVER MIXED.
    // TODO - really need to also determine wheter constraint is even valid (e.g., due to BLANK_DATA in @ts). right now this is withheld until AddConstraintToLP called
    for (int k=0;k<_pGoals[j]->nOperRegimes;k++)
    {
      op_is_active=CheckGoalConditions(j,k,tt,Options);
      //op_is_active=op_is_active && CheckGoalValidity(j,k,tt,Options);
      if (op_is_active){
        _pGoals[j]->conditions_satisfied=true;
        _pGoals[j]->active_regime=k;
        _pGoals[j]->ever_satisfied=true;

        if (_pGoals[j]->is_goal)
        {
          penalty_over =_pGoals[j]->pOperRegimes[k]->penalty_over;
          penalty_under=_pGoals[j]->pOperRegimes[k]->penalty_under;

          penalty_under*=_pGoals[j]->units_correction;
          penalty_over *=_pGoals[j]->units_correction;

          col_ind[i]=GetDVColumnInd(DV_SLACK,s);
          row_val[i]=penalty_over;
          i++;
          if (_pGoals[j]->pOperRegimes[k]->pExpression->compare == COMPARE_IS_EQUAL) //two slack variables
          {
            col_ind[i]=GetDVColumnInd(DV_SLACK,s+1);
            row_val[i]=penalty_under;
            i++;
          }
        }

        k=_pGoals[j]->nOperRegimes+1; //successfully found operating regime - stop searching through remaining operating regimes
      }
    }

    //ensure slack counters are properly incremented
    s++;
    if (_pGoals[j]->pOperRegimes[0]->pExpression->compare == COMPARE_IS_EQUAL) //two slack variables
    {
      s++;
    }
  }

  retval = lp_lib::set_obj_fnex(pLinProg,i, row_val, col_ind);
  ExitGracefullyIf(retval==0,"SolveDemandProblem: Error specifying objective function",RUNTIME_ERR);

  lp_lib::set_row_name(pLinProg, rowcount, "Objective F");

  IncrementAndSetRowName(pLinProg,rowcount,"BlankGoal");

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

    //Q_p = (Qin +Qin2 )*U0 + sum(Ui*Qn-i) + Runoff - Div_out(Q) + Div_in - Delivered + Qadded + U_0*Qadded2
    //Q_p-U_0*Qin-U_0*Qin2... +Delivered = sum(Ui*Qn-i) + Runoff - Div_out(Q) + Div_in + Qadded + U_0*Qadded2
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
        row_val[i]=1.0;
        i++; d++;
      }
      if (pSB->GetReservoir() != NULL) {
        for (int dd = 0; dd<pSB->GetReservoir()->GetNumWaterDemands();dd++)
        {
          d++;
        }
      }

      RHS=0;
      RHS+=pSB->GetSpecifiedInflow (t)*U_0;
      RHS+=pSB->GetDownstreamInflow(t);
      RHS-=aDivGuess[p]; //diverted from (guess based upon last Q)
      RHS+=aDivert[p];   //diverted to (delayed by one day, based upon yesterday's Q)
      for (int n = 1; n < pSB->GetInflowHistorySize(); n++) {
        RHS+=pSB->GetRoutingHydrograph()[n]*pSB->GetInflowHistory()[n-1]; // [n-1] is because this inflow history has not yet been updated
        RHS+=pSB->GetRoutingHydrograph()[n]*pSB->GetSpecifiedInflow(t-n*Options.timestep);
      }
      for (int n = 1; n < pSB->GetLatHistorySize(); n++) {
        RHS += pSB->GetUnitHydrograph()[n] * pSB->GetLatHistory()[n-1];   // [n-1] is because this inflow history has not yet been updated
      }
      RHS+=aSBrunoff[p]/(tstep*SEC_PER_DAY)*pSB->GetUnitHydrograph()[0];  // [m3]->[m3/s]

      retval = lp_lib::add_constraintex(pLinProg,i,row_val,col_ind,ROWTYPE_EQ,RHS);
      ExitGracefullyIf(retval==0,"SolveDemandProblem::Error adding mass balance constraint",RUNTIME_ERR);

      IncrementAndSetRowName(pLinProg,rowcount,"reach_MB_"+to_string(pSB->GetID()));
    }
  }

  //AddLakeMBEquations(pLinProg,s);
  //AddLakeOutflowEquations(pLinProg,s);
  
  // Determine which lakes/reservoirs should have their stage/discharge curve constraints disabled
  // depends upon whether overriding goals are active
  //------------------------------------------------------------------
  bool *aDisableSDCurve= new bool [pModel->GetNumSubBasins()];
  for (p = 0; p<pModel->GetNumSubBasins(); p++)
  {
    aDisableSDCurve[p]=false;
  }
  for (int j = 0; j < _nGoals; j++)
  {
    if ((_pGoals[j]->reservoir_index != DOESNT_EXIST) && 
        (_pGoals[j]->active_regime   != DOESNT_EXIST) && 
        (_pGoals[j]->overrides_SDcurve)) 
    {
      p=_pGoals[j]->reservoir_index;
      aDisableSDCurve[p]=true;
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

      if (pRes->GetHRUIndex() != DOESNT_EXIST) {
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

      IncrementAndSetRowName(pLinProg,rowcount,"reser_MB_"+to_string(pSB->GetID()));

      if (!aDisableSDCurve[p])
      {
        // Second goal equation (relation between Qout and h):
        // may be overridden by other management expressions with larger penalties
        // ------------------------------------------------------------------------
        // linearized storage -discharge curve
        //   Qout = Qout_guess + dQ/dh * (h-h_guess)
        //   1.0 * Qout - dQ/dh * h + s1 -s2 = Q_guess - dQ/dh * h_guess
        double h_guess=h_old;
        double Q_guess;
        double h_sill;

        Q_guess=pRes->GetDischargeFromStage(h_guess,nn);
        dQdh   =pRes->GetStageDischargeDerivative(h_guess,nn);
        h_sill =pRes->GetSillElevation(nn);
        if (_stage_discharge_as_goal)
        {
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

          IncrementAndSetRowName(pLinProg,rowcount,"reserv_Q_"+to_string(pSB->GetID()));

          lprow[p]=lp_lib::get_Nrows(pLinProg);
        }
        else 
        {
          // Stage discharge relation AS IF-ELSE CONSTRAINT : (Formulation from B. Tolson)
          // ------------------------------------------------------------------------
          // linearized storage -discharge curve
          //   Qout = Qout_guess + dQ/dh * (h-h_guess) if h>h_sill
          //   Qout = 0                                if h<h_sill
          // handled via 6 constraint equations 

          const double LARGE_NUMBER=1e7;

          //Eqns 1 & 2 : h^n+1 - M*(1-b) <= h_sill 
          //             h^n+1 + M*(  b) >= h_sill
          col_ind[0]=GetDVColumnInd(DV_STAGE  ,_aResIndices[p]); row_val[0]=1.0;
          col_ind[1]=GetDVColumnInd(DV_BINRES ,_aResIndices[p]); row_val[1]=+LARGE_NUMBER;
          RHS       =h_sill + LARGE_NUMBER;

          retval = lp_lib::add_constraintex(pLinProg,2,row_val,col_ind,ROWTYPE_LE,RHS);
          ExitGracefullyIf(retval==0,"SolveDemandProblem::Error adding stage discharge constraint A",RUNTIME_ERR);

          IncrementAndSetRowName(pLinProg,rowcount,"reserv_Q_A"+to_string(pSB->GetID()));

          col_ind[0]=GetDVColumnInd(DV_STAGE  ,_aResIndices[p]); row_val[0]=1.0;
          col_ind[1]=GetDVColumnInd(DV_BINRES ,_aResIndices[p]); row_val[1]=+LARGE_NUMBER;
          RHS       =h_sill;

          retval = lp_lib::add_constraintex(pLinProg,2,row_val,col_ind,ROWTYPE_GE,RHS);
          ExitGracefullyIf(retval==0,"SolveDemandProblem::Error adding stage discharge constraint B",RUNTIME_ERR);

          IncrementAndSetRowName(pLinProg,rowcount,"reserv_Q_B"+to_string(pSB->GetID()));

          //Eqns 3 & 4 : Q^n+1 >= -M(1-b)
          //             Q^n+1 <=  M(1-b)

          col_ind[0]=GetDVColumnInd(DV_QOUTRES,_aResIndices[p]); row_val[0]=1.0;
          col_ind[1]=GetDVColumnInd(DV_BINRES ,_aResIndices[p]); row_val[1]=-LARGE_NUMBER;
          RHS       =-LARGE_NUMBER;

          retval = lp_lib::add_constraintex(pLinProg,2,row_val,col_ind,ROWTYPE_GE,RHS);
          ExitGracefullyIf(retval==0,"SolveDemandProblem::Error adding stage discharge constraint C",RUNTIME_ERR);

          IncrementAndSetRowName(pLinProg,rowcount,"reserv_Q_C"+to_string(pSB->GetID()));

          col_ind[0]=GetDVColumnInd(DV_QOUTRES,_aResIndices[p]); row_val[0]=1.0;
          col_ind[1]=GetDVColumnInd(DV_BINRES ,_aResIndices[p]); row_val[1]=+LARGE_NUMBER;
          RHS       =+LARGE_NUMBER;

          retval = lp_lib::add_constraintex(pLinProg,2,row_val,col_ind,ROWTYPE_LE,RHS);
          ExitGracefullyIf(retval==0,"SolveDemandProblem::Error adding stage discharge constraint D",RUNTIME_ERR);

          IncrementAndSetRowName(pLinProg,rowcount,"reserv_Q_D"+to_string(pSB->GetID()));

          //Eqns 5 & 6 : Q^n+1 >= Qout_guess + dQ/dh * (h^n+1-h_guess) - M(b)
          //             Q^n+1 <= Qout_guess + dQ/dh * (h^n+1-h_guess) + M(b)
          //             Q^n+1 - dQdh* h^n+1 + Mb >= [Qout_guess - dQ/dh * h_guess]
          //             Q^n+1 - dQdh* h^n+1 - Mb <= [Qout_guess - dQ/dh * h_guess]

          col_ind[0]=GetDVColumnInd(DV_QOUTRES,_aResIndices[p]); row_val[0]=1.0;
          col_ind[1]=GetDVColumnInd(DV_STAGE  ,_aResIndices[p]); row_val[1]=-dQdh;
          col_ind[2]=GetDVColumnInd(DV_BINRES ,_aResIndices[p]); row_val[2]=+LARGE_NUMBER;
          RHS       =Q_guess-dQdh*h_guess;

          retval = lp_lib::add_constraintex(pLinProg,3,row_val,col_ind,ROWTYPE_GE,RHS);
          ExitGracefullyIf(retval==0,"SolveDemandProblem::Error adding stage discharge constraint E",RUNTIME_ERR);

          lprow[p]=lp_lib::get_Nrows(pLinProg);

          IncrementAndSetRowName(pLinProg,rowcount,"reserv_Q_E"+to_string(pSB->GetID()));
        
          col_ind[0]=GetDVColumnInd(DV_QOUTRES,_aResIndices[p]); row_val[0]=1.0;
          col_ind[1]=GetDVColumnInd(DV_STAGE  ,_aResIndices[p]); row_val[1]=-dQdh;
          col_ind[2]=GetDVColumnInd(DV_BINRES ,_aResIndices[p]); row_val[2]=-LARGE_NUMBER;
          RHS       =Q_guess-dQdh*h_guess;

          retval = lp_lib::add_constraintex(pLinProg,3,row_val,col_ind,ROWTYPE_LE,RHS);
          ExitGracefullyIf(retval==0,"SolveDemandProblem::Error adding stage discharge constraint F",RUNTIME_ERR);

          IncrementAndSetRowName(pLinProg,rowcount,"reserv_Q_F"+to_string(pSB->GetID()));
        }
      }
      else /* if (aDisableSDCurve[p])*/ //keep same rows - make inert equations 
      {
        if (_stage_discharge_as_goal)
        {
          col_ind[0]=GetDVColumnInd(DV_SLACK  ,s);row_val[0]=1.0; RHS=0;

          retval = lp_lib::add_constraintex(pLinProg,1,row_val,col_ind,ROWTYPE_EQ,RHS);
          IncrementAndSetRowName(pLinProg,rowcount,"reserv_Q_"+to_string(pSB->GetID()));
        }
        else {
          col_ind[0]=GetDVColumnInd(DV_BINRES ,_aResIndices[p]); row_val[0]=1.0;  RHS=0;

          retval = lp_lib::add_constraintex(pLinProg,1,row_val,col_ind,ROWTYPE_LE,RHS);
          IncrementAndSetRowName(pLinProg,rowcount,"reserv_Q_A"+to_string(pSB->GetID()));
          retval = lp_lib::add_constraintex(pLinProg,1,row_val,col_ind,ROWTYPE_LE,RHS);
          IncrementAndSetRowName(pLinProg,rowcount,"reserv_Q_B"+to_string(pSB->GetID()));
          retval = lp_lib::add_constraintex(pLinProg,1,row_val,col_ind,ROWTYPE_LE,RHS);
          IncrementAndSetRowName(pLinProg,rowcount,"reserv_Q_C"+to_string(pSB->GetID()));
          retval = lp_lib::add_constraintex(pLinProg,1,row_val,col_ind,ROWTYPE_LE,RHS);
          IncrementAndSetRowName(pLinProg,rowcount,"reserv_Q_D"+to_string(pSB->GetID()));
          retval = lp_lib::add_constraintex(pLinProg,1,row_val,col_ind,ROWTYPE_LE,RHS);
          IncrementAndSetRowName(pLinProg,rowcount,"reserv_Q_E"+to_string(pSB->GetID()));
          retval = lp_lib::add_constraintex(pLinProg,1,row_val,col_ind,ROWTYPE_LE,RHS);
          IncrementAndSetRowName(pLinProg,rowcount,"reserv_Q_F"+to_string(pSB->GetID()));
        }
        s+=2; //to ensure EnvMin counter is working
      }
      res_count++;
    }
  }

  // Environmental Minimum Flow Goals
  // ----------------------------------------------------------------
  d=0;
  double minQ,ufp;
  const double ENV_FLOW_PENALTY=1e6;
  //s=_firstEnvFlowSlackIndex;
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
        col_ind[0]=GetDVColumnInd(DV_QOUT,_aSBIndices[p]);    row_val[0]=+1.0;
        col_ind[1]=GetDVColumnInd(DV_SLACK,s);                row_val[1]=-1.0;
        RHS       =minQ;

        retval = lp_lib::add_constraintex(pLinProg,2,row_val,col_ind,ROWTYPE_LE,RHS);
        ExitGracefullyIf(retval==0,"SolveDemandProblem::Error adding environmental flow goal A",RUNTIME_ERR);

        IncrementAndSetRowName(pLinProg,rowcount,"envMin_A"+to_string(pSB->GetID()));

        // 10^6 * e+ >= (all upstream D)
        col_ind[0]=GetDVColumnInd(DV_SLACK,s);                row_val[0]=ENV_FLOW_PENALTY;
        s++;

        i=1;
        for (int dd = 0; dd<_aUpCount[p];dd++)
        {
          col_ind[i]=GetDVColumnInd(DV_DELIVERY,_aUpstreamDemands[p][dd]);
          row_val[i]=-1.0;
          i++;
        }

        RHS=0.0;

        retval = lp_lib::add_constraintex(pLinProg,i,row_val,col_ind,ROWTYPE_GE,RHS);
        ExitGracefullyIf(retval==0,"SolveDemandProblem::Error adding environmental flow goal B",RUNTIME_ERR);

        IncrementAndSetRowName(pLinProg,rowcount,"envMin_B"+to_string(pSB->GetID()));
      }

      //TODO - figure out how to extend this to an unusable flow percentage
      // Q + e* >= (ufp)*[sum(D)+Q] <-term in brackets is inflow upstream of demands
      ufp =pSB->GetUnusableFlowPercentage();
      if (ufp != 0.0) {
        // e* <= 10^6 * all upstream D
      }
    }
  }

  // user-specified constraints / goals
  // ----------------------------------------------------------------
  for (int i = 0; i < _nGoals; i++)
  {
    AddConstraintToLP( i, _pGoals[i]->active_regime, pLinProg, tt, col_ind, row_val);
    IncrementAndSetRowName(pLinProg,rowcount,_pGoals[i]->name);
  }

  // ----------------------------------------------------------------
  // ITERATIVELY SOLVE OPTIMIZATION PROBLEM WITH LP_SOLVE
  // ----------------------------------------------------------------  
  const int NUM_ITERATIONS=7;

  int     ctyp;
  int     nrows  =lp_lib::get_Nrows(pLinProg);
  double *soln   =new double [_nDecisionVars];
  double *constr =new double [nrows];
  if (_aSolverResiduals==NULL){ //TODO: problem: this needs to be max nrows 
    _aSolverResiduals=new double [nrows]; _nSolverResiduals=nrows;
    _aSolverRowNames =new string [nrows]; 
  }

  string dumpfile     =Options.output_dir+"/lp_solve_dump.txt";

  lp_lib::set_add_rowmode(pLinProg, FALSE);    //must be turned off once model goals/constraints and obj function are added
  lp_lib::set_minim      (pLinProg);           //ensures this is treated as a minimization probleme         
  lp_lib::set_verbose    (pLinProg,IMPORTANT);
  //lp_lib::set_scaling  (pLinProg, CURTISREIDSCALE);  //May wish to look for scaling improvements
  //lp_lib::set_BFP      (pLinProg, "bfp_LUSOL");      //TEST THIS FOR SPEED

  for (int iter=0; iter<NUM_ITERATIONS; iter++)
  {
    if (_do_debug_level==2)
    {
      //lp_lib::print_debugdump(pLinProg,&dumpfile[0]);
      //lp_lib::set_outputfile (pLinProg,&dumpfile[0]);
      //lp_lib::print_lp(pLinProg);
      WriteLPSubMatrix(pLinProg,"lp_matrix.csv",Options);
    }
    //

    // ----------------------------------------------------------------
    // Solve the LP problem!
    // ----------------------------------------------------------------
    retval = lp_lib::solve(pLinProg);

    // handle solve error
    // ----------------------------------------------------------------
    if (retval!=OPTIMAL)
    {
      cout<<"lp_lib::solve error code: "<<retval<<endl; 
      WriteLPSubMatrix(pLinProg,"overconstrained_lp_matrix.csv",Options);
      ExitGracefully("SolveDemandProblem: non-optimal solution found. Problem is over-constrained. Remove or adjust management constraints.",RUNTIME_ERR);
    }
    
    // assign decision variables
    // ----------------------------------------------------------------

    lp_lib::get_variables  (pLinProg,soln);
    lp_lib::get_constraints(pLinProg,constr); //constraint matrix * soln

    for (int i=0;i<_nDecisionVars;i++){
      _pDecisionVars[i]->value=soln[i];
    }

    // Extract solver residuals-  should be near zero everywhere
    // ----------------------------------------------------------------
    for (int j=0;j<nrows;j++)
    {
      if (j<_nSolverResiduals){ //tmp debug to avoid memory leak
        _aSolverRowNames [j]=to_string(lp_lib::get_row_name(pLinProg,j+1)); 
        _aSolverResiduals[j]=constr[j]-lp_lib::get_rh(pLinProg,j+1); // LHS-RHS for each goal/constraint
        ctyp=lp_lib::get_constr_type(pLinProg,j+1);
        if (ctyp==LE){upperswap(_aSolverResiduals[j],0.0); }
        if (ctyp==GE){lowerswap(_aSolverResiduals[j],0.0); _aSolverResiduals[j]*=-1; }
      }
    }

    // grab and store improved estimate of reservoir stage + flow diversion
    // ----------------------------------------------------------------
    double  value;
    dv_type typ;
    for (int i=0;i<_nDecisionVars;i++)
    {
      p     = _pDecisionVars[i]->p_index;
      value = _pDecisionVars[i]->value;
      typ   = _pDecisionVars[i]->dvar_type;

      if      (typ == DV_STAGE) {h_iter[p]=value;}
      else if (typ == DV_QOUT ) {Q_iter[p]=value;}
    }

    // update lp matrix for iterative solution of stage-discharge curve and flow diversions
    // ----------------------------------------------------------------
    int    hcol;
    double Q_guess,h_guess;
    for (int pp = 0; pp<pModel->GetNumSubBasins(); pp++)
    {
      p   =pModel->GetOrderedSubBasinIndex(pp);
      pSB =pModel->GetSubBasin(p);
      pRes=pSB->GetReservoir();
      if ((pSB->IsEnabled()) && (pRes!=NULL) && (!aDisableSDCurve[p]))
      {
        h_guess=h_iter[p];
        Q_guess=pRes->GetDischargeFromStage      (h_guess,nn);
        dQdh   =pRes->GetStageDischargeDerivative(h_guess,nn);

        hcol   =GetDVColumnInd(DV_STAGE  ,_aResIndices[p]);

        if (_stage_discharge_as_goal){
          if (fabs(dQdh)>REAL_SMALL) {RHS=Q_guess-dQdh*h_guess;}
          else                       {RHS=0.0;                 }//below weir height
          lp_lib::set_mat(pLinProg,lprow[p],hcol,-dQdh);
          lp_lib::set_rh (pLinProg,lprow[p],RHS);
        }
        else{
          //edits 2 consecutive equations 
          RHS=Q_guess-dQdh*h_guess;
          lp_lib::set_mat(pLinProg,lprow[p]  ,hcol,-dQdh);
          lp_lib::set_rh (pLinProg,lprow[p]  ,RHS);
          lp_lib::set_mat(pLinProg,lprow[p]+1,hcol,-dQdh);
          lp_lib::set_rh (pLinProg,lprow[p]+1,RHS);
        }
      }
    }

    // update lp matrix for iterative solution of flow diversions
    // ----------------------------------------------------------------
    int    junk;
    double div_Q,sum_diverted;
    for (int pp = 0; pp<pModel->GetNumSubBasins(); pp++)
    {
      p   =pModel->GetOrderedSubBasinIndex(pp);
      pSB =pModel->GetSubBasin(p);
      if (pSB->GetNumDiversions()>0)
      {
        div_Q=0;
        sum_diverted=0;
        Q_guess=Q_iter[p];
        for(int i=0;i<pSB->GetNumDiversions();i++)
        {
          div_Q=pSB->GetDiversionFlow(i,Q_guess,Options,tt,junk);
          sum_diverted+=div_Q;
        }

        RHS=lp_lib::get_rh(pLinProg,lprow[p]);//TMP DEBUG - NOT SAME LPROW AS RESERVOIRS !!!

        RHS+=(sum_diverted-aDivGuess[p]);
        aDivGuess[p]=sum_diverted;

        lp_lib::set_rh (pLinProg,lprow[p],RHS);
      }
    }
  }/*end iteration loop*/

  //Report post-iteration objective function: This will be sum of penalties*violations
  if (_do_debug_level>=2)
  {
    cout<<"Objective value: "<<lp_lib::get_objective(pLinProg)+demand_penalty_sum<<" solver rows:"<<_nSolverResiduals<<endl;
  }

  lp_lib::delete_lp(pLinProg);
  delete [] soln;
  delete [] constr;
  delete [] col_ind;
  delete [] row_val;
  delete [] h_iter;
  delete [] Q_iter;
  delete [] lprow;
  delete [] aDivert;
  delete [] aDivGuess;
  delete [] aDisableSDCurve;

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

    if      (typ == DV_QOUT) {}//do nothing; recalculated by base Raven routing
    else if (typ == DV_STAGE){}//do nothing; recalculated by base Raven reservoir routing
    else if (typ == DV_USER) {}//doesn't influence hydrologic model
    else if (typ == DV_QOUTRES)
    {
      pModel->GetSubBasin(p)->GetReservoir()->SetOptimizedOutflow(value); //TMP DEBUG this should ONLY be in reservoirs where one or more management constraints are applied?
    }
    else if (typ == DV_DELIVERY)
    {
      d =_pDecisionVars[i]->loc_index;
      ii=_pDecisionVars[i]->dem_index;

      _aDelivery   [d] =value;
      _aCumDelivery[d]+=value*(Options.timestep*SEC_PER_DAY);

      if (tt.julian_day==_pDemands[d]->GetCumulDeliveryDate()){
        _aCumDelivery[d]=0.0; 
      }

      if (!_pDemands[d]->IsReservoirDemand()){
        pModel->GetSubBasin(p)->AddToDeliveredDemand(ii,value);
      }
      else {
        pModel->GetSubBasin(p)->GetReservoir()->AddToDeliveredDemand(ii, value);
      }
    }
    else if (typ == DV_SLACK)
    { 
      s = _pDecisionVars[i]->loc_index;
      _aSlackValues[s]=value;
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

  // ManagementOptimization.csv
  //--------------------------------------------------------
  string tmpFilename=FilenamePrepare("ManagementOptimization.csv",Options);
  _MANOPT.open(tmpFilename.c_str());
  if (_MANOPT.fail()){
    ExitGracefully(("CModel::WriteOutputFileHeaders: unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
  }

  _MANOPT<<"time,date,hour";
  for (int i=0;i<_nDecisionVars;i++){
    if ((_pDecisionVars[i]->p_index == DOESNT_EXIST) &&
        ((_pDecisionVars[i]->dvar_type!=DV_SLACK) || (_do_debug_level>0))) { //user specified non-slack variable (slack included if debug>0)
       _MANOPT<<","<<_pDecisionVars[i]->name;
    }
    else{
      pSB=_pModel->GetSubBasin(_pDecisionVars[i]->p_index);
      if (pSB->IsGauged()){
        _MANOPT<<","<<_pDecisionVars[i]->name;
      }
    }
  }
  _MANOPT<<endl;

  // GoalSatisfaction.csv
  //--------------------------------------------------------
  tmpFilename=FilenamePrepare("GoalSatisfaction.csv",Options);
  _GOALSAT.open(tmpFilename.c_str());
  if (_GOALSAT.fail()){
    ExitGracefully(("CModel::WriteOutputFileHeaders: unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
  }
  _GOALSAT<<"time,date,hour";
  for (int i = 0; i < _nGoals; i++) {
    if (_pGoals[i]->is_goal) {
      _GOALSAT<<","<<_pGoals[i]->name+" penalty";
      _GOALSAT<<","<<_pGoals[i]->name+" regime";
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

  // ManagementOptimization.csv
  //--------------------------------------------------------
  _MANOPT<<tt.model_time <<","<<usedate<<","<<usehour;
  for (int i=0;i<_nDecisionVars;i++){
    if ((_pDecisionVars[i]->p_index == DOESNT_EXIST) &&
        ((_pDecisionVars[i]->dvar_type!=DV_SLACK) || (_do_debug_level>0))) { //user specified non-slack variable (slack included if debug>0)
       _MANOPT<<","<<_pDecisionVars[i]->value;
    }
    else{
      pSB=_pModel->GetSubBasin(_pDecisionVars[i]->p_index);
      if (pSB->IsGauged()){ //gauged basin or not associated with basin
        _MANOPT<<","<<_pDecisionVars[i]->value;
      }
    }
  }
  _MANOPT<<endl;

  // GoalSatisfaction.csv
  //--------------------------------------------------------
  // All slack values reported should be positive
  // each entry represents the magnitude of a goal violation, expressed in m3/s
  // for environmental flow goals, for instance, this is E-Q if Q<E, 0 otherwise
  bool   include_pen=false;
  double pen;
  double mult=1.0;
  
  _GOALSAT<<tt.model_time <<","<<usedate<<","<<usehour;

  //first slack terms are due to environmental flow constraints (_nEnviroFlowGoals) and reservoir outflow targets (2*_nReservoirs)
  //this is first index of user-specified slack variable
  
  int s=_nEnviroFlowGoals+2*_nReservoirs; 
  if (!_stage_discharge_as_goal) {
    //s=_nEnviroFlowGoals; //No slack variables for reservoirs (currently retained regardless of stage discharge status)
  }
  for (int i = 0; i < _nGoals; i++) 
  {
    mult=1.0;
    if (_pGoals[i]->is_goal) 
    {
      if (include_pen){mult=_pGoals[i]->penalty_over*_pGoals[i]->units_correction;}
      if (_pGoals[i]->pOperRegimes[0]->pExpression->compare != COMPARE_LESSTHAN){mult*=-1.0; }

      pen = mult*_aSlackValues[s];
      s++;
      if (_pGoals[i]->pOperRegimes[0]->pExpression->compare == COMPARE_IS_EQUAL) 
      {
        mult=1.0;
        if (include_pen){mult=_pGoals[i]->penalty_under*_pGoals[i]->units_correction;}
        pen+=mult*_aSlackValues[s];
        s++;
      }
      _GOALSAT<<","<<pen;

      //print active operating regime
      if (!_pGoals[i]->conditions_satisfied)
      {
        _GOALSAT<<",INACTIVE";
      }
      else {
        string active_regime;
        active_regime=_pGoals[i]->pOperRegimes[_pGoals[i]->active_regime]->reg_name;
        _GOALSAT<<","<<active_regime;
      }
    }
  }
  _GOALSAT<<endl;

  // OptimizationResidual.csv
  //--------------------------------------------------------
  if (_do_debug_level > 0) 
  {
    //Header must be written during simulation because row names don't exist until lp solver build
    if (fabs(tt.model_time - Options.timestep)<PRETTY_SMALL)
    {
      string tmpFilename=FilenamePrepare("OptimizationResidual.csv",Options);
      _RESID.open(tmpFilename.c_str());
      if (_RESID.fail()){
        ExitGracefully(("CModel::WriteOutputFileHeaders: unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
      }
      _RESID<<"time,date,hour";
      for (int i = 0; i < _nSolverResiduals; i++) {
        _RESID<<","<<_aSolverRowNames[i];
      }
      _RESID<<endl;
    }

    _RESID<<tt.model_time <<","<<usedate<<","<<usehour;
    for (int i = 0; i < _nSolverResiduals; i++) {
      _RESID<<","<<_aSolverResiduals[i];
    }
    _RESID<<endl;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Closes output file streams
/// \details after end of simulation from Main() or in ExitGracefully; All file streams are opened in WriteOutputFileHeaders() routine
//
void CDemandOptimizer::CloseOutputStreams()
{
  if (_MANOPT.is_open()){_MANOPT.close();}
  if (  _GOALSAT.is_open()){  _GOALSAT.close();}
  if (    _RESID.is_open()){    _RESID.close();}
}
//////////////////////////////////////////////////////////////////
/// \brief writes lp matrix to file in useful .csv format
/// \details used in debugging
/// file should not include full path and should end in .csv
//
void CDemandOptimizer::WriteLPSubMatrix(lp_lib::lprec* pLinProg,string file, const optStruct &Options) const
{
  ofstream LPMAT;
  string filename=Options.output_dir+file;
  LPMAT.open(filename);
  if(LPMAT.fail()) {
    ExitGracefully(("CConstituentModel::WriteOutputFileHeaders: Unable to open output file "+filename+" for writing.").c_str(),FILE_OPEN_ERR);
  }
  int nrows,ncols;
  nrows=lp_lib::get_Nrows(pLinProg);
  ncols=lp_lib::get_Ncolumns(pLinProg);
  
  LPMAT<<"rowname,";
  for (int i = 0; i < ncols; i++) {
    LPMAT<<lp_lib::get_col_name(pLinProg,i+1)<<",";
  }
  LPMAT<<"compare,RHS"<<endl;

  string comp;
  REAL *val=new REAL [ncols+1];
  for (int j=0; j<=nrows; j++)
  {
    LPMAT<<lp_lib::get_row_name(pLinProg,j)<<",";
    lp_lib::get_row(pLinProg,j,val);
    for (int i = 0; i < ncols; i++) {
      LPMAT<<val[i+1]<<",";
    }
    int ctype=lp_lib::get_constr_type(pLinProg,j);
    if      (ctype==1){comp="<="; }
    else if (ctype==2){comp=">="; }
    else if (ctype==3){comp="=="; }

    LPMAT<<comp<<",";
    LPMAT<<lp_lib::get_rh(pLinProg,j)<<endl;
  }
  LPMAT.close();
  delete [] val;
}

//////////////////////////////////////////////////////////////////
/// \brief Provides end warnings
/// \details after end of simulation from Main()
//
void CDemandOptimizer::Closure(const optStruct &Options)
{
  string warn;
  for (int j = 0; j < _nGoals; j++) {
    int nConditions=0;
    for (int k = 0; k < _pGoals[j]->nOperRegimes; k++) {
      nConditions+=_pGoals[j]->pOperRegimes[k]->nConditions;
    }

    if (!(_pGoals[j]->ever_satisfied) && (nConditions>0)) {
      warn="CDemandOptimizer::Closure: constraint/goal ["+to_string(j)+"] "+_pGoals[j]->name+" was never activated during simulation, as none of its conditions were ever met.";
      WriteWarning(warn,Options.noisy);
    }
  }

}
