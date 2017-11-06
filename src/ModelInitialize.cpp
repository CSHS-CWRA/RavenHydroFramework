/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Model.h"

/*****************************************************************
   Model Initialization Routines
------------------------------------------------------------------
   called prior to simulation
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Initializes model prior to simulation
/// \details Perform all operations required and initial check on validity of model before simulation begins;
/// -initializes mass balance arrays to zero; identifies model UTM zone;
/// -calls initialization routines of all subasins, HRUs, processes, and gauges
/// -generates gauge weights
/// -calculates stream network topology, routing orders
/// -calculates initial water/energy/mass storage
/// -calls routines to write initial conditions to minor output
///
/// \param &Options [in] Global model options information
//
void CModel::Initialize(const optStruct &Options)
{
  int g,i,j,k,p;

  // Quality control
  //--------------------------------------------------------------
  ExitGracefullyIf(_nSubBasins<1,
                   "CModel::Initialize: Must have at least one SubBasin",BAD_DATA);
  ExitGracefullyIf(_nHydroUnits<1,
                   "CModel::Initialize: Must have at least one hydrologic unit",BAD_DATA);
  ExitGracefullyIf(_nGauges<1 && _nForcingGrids<1,
                   "CModel::Initialize: Must have at least one meteorological gauge station or forcing grid",BAD_DATA);
  ExitGracefullyIf(_nProcesses==0,
                   "CModel::Initialize: must have at least one hydrological process included in model",BAD_DATA);

  //Ensure Basins & HRU IDs are unique
  for (p=0;p<_nSubBasins;p++){
    for (int pp=0;pp<p;pp++){
      if (_pSubBasins[p]->GetID()==_pSubBasins[pp]->GetID()){
        ExitGracefully("CModel::Initialize: non-unique basin identifier found",BAD_DATA);}}}

  for (k=0;k<_nSubBasins;k++){
    for (int kk=0;kk<p;kk++){
      if ((k!=kk) && (_pHydroUnits[k]->GetID()==_pHydroUnits[kk]->GetID())){
        ExitGracefully("CModel::Initialize: non-unique HRU identifier found",BAD_DATA);}}}

  if ((_nSnowLayers==0) && (StateVarExists(SNOW))){_nSnowLayers=1;}




  // initialize process algorithms, initialize water/energy balance arrays to zero
  //--------------------------------------------------------------
  _nTotalConnections=0;
  for (int j=0; j<_nProcesses;j++){
    if((_pProcesses[j]->GetProcessType()!=PRECIPITATION) && 
       (_pProcesses[j]->GetProcessType()!=LAT_ADVECTION))
       {_pProcesses[j]->Initialize();} //precip already initialized in ParseInput.cpp
    _nTotalConnections+=_pProcesses[j]->GetNumConnections();
  }
  if (_pTransModel->GetNumConstituents()>0){
    _pTransModel->CalculateLateralConnections();
    for (int j=0; j<_nProcesses;j++){
      if (_pProcesses[j]->GetProcessType()==LAT_ADVECTION){
        _pProcesses[j]->Initialize(); //must be re-initialized after all lateral processes are initialized above
        _nTotalConnections+=_pProcesses[j]->GetNumConnections();
      } 
    }
  }

  _aCumulativeBal = new double * [_nHydroUnits];
  _aFlowBal       = new double * [_nHydroUnits];
  for (k=0; k<_nHydroUnits;k++)
  {
    _aCumulativeBal[k]= NULL;
    _aCumulativeBal[k]= new double [_nTotalConnections];
    ExitGracefullyIf(_aCumulativeBal[k]==NULL,"CModel::Initialize (aCumulativeBal)",OUT_OF_MEMORY);
    for (int js=0;js<_nTotalConnections;js++){_aCumulativeBal[k][js]=0.0;}

    _aFlowBal[k]= NULL;
    _aFlowBal[k]= new double [_nTotalConnections];
    ExitGracefullyIf(_aFlowBal[k]==NULL,"CModel::Initialize (aFlowBal)",OUT_OF_MEMORY);
    for (int js=0;js<_nTotalConnections;js++){_aFlowBal[k][js]=0.0;}
  }

  _nTotalLatConnections=0;
  for (int j=0; j<_nProcesses;j++){
    _nTotalLatConnections+=_pProcesses[j]->GetNumLatConnections();
  }
  if(_nTotalLatConnections>0){
    _aFlowLatBal=NULL;
    _aCumulativeLatBal = new double[_nTotalLatConnections];
    _aFlowLatBal       = new double[_nTotalLatConnections];
    ExitGracefullyIf(_aFlowLatBal==NULL,"CModel::Initialize (_aFlowLatBal)",OUT_OF_MEMORY);
    
    for(int jss=0;jss<_nTotalLatConnections;jss++){
      _aCumulativeLatBal[jss]=0.0;
      _aFlowLatBal      [jss]=0.0;
    }
  }
  _CumulInput   =_CumulOutput  =0.0;
  _CumEnergyGain=_CumEnergyLoss=0.0;

  //Identify model UTM_zone for interpolation
  //--------------------------------------------------------------
  double cen_long(0),area_tot(0);//longiturde of area-weighted watershed centroid, total wshed area
  CHydroUnit *pHRU;
  for (k=0; k<_nHydroUnits;k++)
  {
    pHRU=_pHydroUnits[k];
    area_tot+=pHRU->GetArea();
    cen_long+=pHRU->GetCentroid().longitude/_nHydroUnits*(pHRU->GetArea());
  }
  cen_long/=area_tot;
  _UTM_zone =(int)(floor ((cen_long + 180.0) / 6) + 1);

  //Initialize HRUs, gauges and transient parameters
  //--------------------------------------------------------------
  for (k=0;k<_nHydroUnits; k++){_pHydroUnits [k]->Initialize(_UTM_zone);}
  for (g=0;g<_nGauges;     g++){_pGauges     [g]->Initialize(Options,_UTM_zone);}
  for (j=0;j<_nTransParams;j++){_pTransParams[j]->Initialize(Options);}
  // Forcing grids are not "Initialized" here because the derived data have to be populated everytime a new chunk is read

  //Generate Gauge Weights from Interpolation
  //--------------------------------------------------------------
  if (!Options.silent){cout <<"  Generating Gauge Interpolation Weights..."<<endl;}
  GenerateGaugeWeights(Options);

  //Initialize SubBasins, calculate routing orders, topology
  //--------------------------------------------------------------
  if (!Options.silent){cout<<"  Calculating basin & watershed areas..."<<endl;}
  _WatershedArea=0.0;
  for (p=0;p<_nSubBasins;p++){_WatershedArea+=_pSubBasins[p]->CalculateBasinArea();}

  if (!Options.silent){cout<<"  Calculating routing network topology..."<<endl;}
  InitializeRoutingNetwork(); //calculate proper routing orders

  if (!Options.silent){cout<<"  Initializing Basins, calculating watershed area, setting initial flow conditions..."<<endl;}
  InitializeBasinFlows(Options);

  //Calculate initial system water storage
  //--------------------------------------------------------------
  if (!Options.silent){cout<<"  Calculating initial system water storage..."<<endl;}
  _initWater=0.0;
  double S=0;
  for (i=0;i<_nStateVars;i++)
  {
    if (CStateVariable::IsWaterStorage(_aStateVarType[i])){
      S=GetAvgStateVar(i);
      _initWater+=S;
    }
  }
  _initWater+=GetTotalChannelStorage();
  _initWater+=GetTotalReservoirStorage();
  _initWater+=GetTotalRivuletStorage();
  // \todo [fix]: this fixes a mass balance bug in reservoir simulations, but there is certainly a more proper way to do it
  // I think somehow this is being double counted in the delta V calculations in the first timestep
  for(int p=0;p<_nSubBasins;p++){
    if(_pSubBasins[p]->GetReservoir()!=NULL){
      _initWater+=_pSubBasins[p]->GetIntegratedReservoirInflow(Options.timestep)/2.0/_WatershedArea*MM_PER_METER/M2_PER_KM2;
      _initWater-=_pSubBasins[p]->GetIntegratedOutflow(Options.timestep)/2.0/_WatershedArea*MM_PER_METER/M2_PER_KM2;
    }
  }

  //Initialize Transport
  //--------------------------------------------------------------
  if (_pTransModel->GetNumConstituents()>0){
    if (!Options.silent){cout<<"  Initializing Transport Model..."<<endl;}
    _pTransModel->Initialize();
  }

  // precalculate whether individual processes should apply (for speed)
  //--------------------------------------------------------------
  _aShouldApplyProcess = new bool *[_nProcesses];
  for (int j=0; j<_nProcesses;j++){
    _aShouldApplyProcess[j]=NULL;
    _aShouldApplyProcess[j] = new bool [_nHydroUnits];
    ExitGracefullyIf(_aShouldApplyProcess[j]==NULL,"CModel::Initialize (_aShouldApplyProcess)",OUT_OF_MEMORY);
    for (k=0; k<_nHydroUnits;k++){
      _aShouldApplyProcess[j][k] = _pProcesses[j]->ShouldApply(_pHydroUnits[k]);
    }
  }

  //Initialize NetCDF Output File IDs
  //--------------------------------------------------------------
  /* initialize all potential NetCDF file IDs with -9 == "not existing and hence not opened" */
  _HYDRO_ncid    = -9;   // output file ID for Hydrographs.nc         (-9 --> not opened)
  _STORAGE_ncid  = -9;   // output file ID for WatershedStorage.nc    (-9 --> not opened)
  _FORCINGS_ncid = -9;   // output file ID for ForcingFunctions.nc    (-9 --> not opened)

  //Write Output File Headers
  //--------------------------------------------------------------
  for (int c=0;c<_nCustomOutputs;c++){
    _pCustomOutputs[c]->InitializeCustomOutput(Options);
  }
  if (!Options.silent){cout<<"  Writing Output File Headers..."<<endl;}
  WriteOutputFileHeaders(Options);

  quickSort(_aOutputTimes,0,_nOutputTimes-1);

  //Prepare Output Time Series
  //--------------------------------------------------------------
  InitializeObservations(Options);

  //General QA/QC
  //--------------------------------------------------------------
  ExitGracefullyIf((GetNumGauges()<2) && (Options.orocorr_temp==OROCORR_UBCWM2),
                   "CModel::Initialize: at least 2 gauges necessary to use :OroTempCorrect method OROCORR_UBCWM2", BAD_DATA);
  for (int kk = 0; kk < _nHRUGroups; kk++){
    if (_pHRUGroups[kk]->GetNumHRUs() == 0){
      string warn = "CModel::Initialize: HRU Group " + _pHRUGroups[kk]->GetName() + " is empty.";
      WriteWarning(warn,Options.noisy);
    }
  }
  for(int i=0; i<_nObservedTS; i++){
    if(!strcmp(_pObservedTS[i]->GetName().c_str(),"RESERVOIR_STAGE"))
    {
      long SBID=s_to_l(_pObservedTS[i]->GetTag().c_str());
      if(GetSubBasinByID(SBID)->GetReservoir()==NULL){
        string warn="Observations supplied for non-existent reservoir in subbasin "+to_string(SBID);
        ExitGracefully(warn.c_str(),BAD_DATA);
      }
    }
  }
  bool wetlandsinmodel=false;
  for(k=0;k<_nHydroUnits;k++) {
    if (_pHydroUnits[k]->GetHRUType()==HRU_WETLAND){wetlandsinmodel=true;break;}
  }
  if((wetlandsinmodel) && (GetStateVarIndex(DEPRESSION)==DOESNT_EXIST)) {
    ExitGracefully("CModel::Initialize: At least one WETLAND-type soil profile is included but no DEPRESSION storage processes included in hydrologic process list.",BAD_DATA_WARN);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Initializes observation time series
/// \details Initializes _pObservedTS, _pModeledTS and _pObsWeightTS
///     Matches observation weights to the corresponding observations
///     Called from CModel::Initialize
///
/// \param &Options [in] Global model options information
//
void CModel::InitializeObservations(const optStruct &Options)
{
  int nModeledValues =(int)(ceil((Options.duration+TIME_CORRECTION)/Options.timestep)+1);
  _pModeledTS=new CTimeSeries * [_nObservedTS];
  _aObsIndex =new int           [_nObservedTS];
  CTimeSeriesABC** tmp = new CTimeSeriesABC *[_nObservedTS];
  for (int i = 0; i < _nObservedTS; i++)
  {
    _pModeledTS[i] = new CTimeSeries("MODELED" + _pObservedTS[i]->GetName(), _pObservedTS[i]->GetTag(),"",Options.julian_start_day,Options.julian_start_year,Options.timestep,nModeledValues,true);
    _pObservedTS[i]->Initialize(Options.julian_start_day, Options.julian_start_year, Options.duration, max(Options.timestep,_pObservedTS[i]->GetInterval()),true);
    _pModeledTS[i]->InitializeResample(_pObservedTS[i]->GetNumSampledValues(),_pObservedTS[i]->GetSampledInterval());
    _aObsIndex[i]=0;

    //Match weights with observations based on Name, tag and numValues
    tmp[i] = NULL;
    for (int n = 0; n < _nObsWeightTS; n++){
      if ( _pObsWeightTS[n]!=NULL
           && _pObsWeightTS[n]->GetName()     == _pObservedTS[i]->GetName()
           && _pObsWeightTS[n]->GetTag()      == _pObservedTS[i]->GetTag()
           && _pObsWeightTS[n]->GetNumValues()== _pObservedTS[n]->GetNumValues()  )
      {
        tmp[i] = _pObsWeightTS[n];
        _pObsWeightTS[n] = NULL;
        tmp[i]->Initialize(Options.julian_start_day, Options.julian_start_year, Options.duration,Options.timestep,true);
      }
    }
  }

  //clean up and warn about unmatched weights
  for (int n = 0; n < _nObsWeightTS; n++){
    if (_pObsWeightTS[n] != NULL){
      WriteWarning("Observation Weight "+_pObsWeightTS[n]->GetName()+" "+_pObsWeightTS[n]->GetTag()+" not matched to observation time series", Options.noisy);
      delete _pObsWeightTS[n]; _pObsWeightTS[n]=NULL;
    }
  }

  delete[] _pObsWeightTS;
  _pObsWeightTS = tmp;
}

//////////////////////////////////////////////////////////////////
/// \brief Initializes routing network
/// \details Calculates sub basin routing order - generates _aOrderedSBind array
/// which stores upstream -> downstream order \n
///    - outlet basins have order of zero,
///    - basins that drain to zero order basins have order of 1,
///    - basins that drain to order 1 basins have order of 2,
///    - etc.
/// \remark Called prior to simulation from CModel::Initialize
//
void CModel::InitializeRoutingNetwork()
{
  int p,pp,pTo,ord;
  const int max_iter=100;
  bool noisy=false;//useful for debugging

  _aDownstreamInds=NULL;
  _aSubBasinOrder =new int [_nSubBasins];
  _aDownstreamInds=new int [_nSubBasins];
  ExitGracefullyIf(_aDownstreamInds==NULL,"CModel::InitializeRoutingNetwork(1)",OUT_OF_MEMORY);

  //check for bad downstream IDs, populate downstream_ind array
  //----------------------------------------------------------------------
  for (p=0;p<_nSubBasins;p++)
  {
    _aSubBasinOrder[p]=0;

    pp=GetSubBasinIndex(_pSubBasins[p]->GetDownstreamID());
    ExitGracefullyIf(pp==INDEX_NOT_FOUND,
                     "CModel::InitializeRoutingNetwork: downstream basin ID not found",BAD_DATA);//Warning only?
    ExitGracefullyIf(pp==p,
                     "CModel::InitializeRoutingNetwork: subbasin empties into itself: circular reference!",BAD_DATA);
    _aDownstreamInds[p]=pp;
  }

  //iterative identification of subbasin orders
  //----------------------------------------------------------------------
  //here, order goes from 0 (trunk) to _maxSubBasinOrder (leaf)
  //This is NOT Strahler ordering!!!
  //JRC:  there might be a faster way to do this
  int last_ordersum;
  int iter(0),ordersum(0);
  do
  {
    last_ordersum=ordersum;
    for (p=0;p<_nSubBasins;p++)
    {
      pTo=_aDownstreamInds[p];
      if (pTo==DOESNT_EXIST){_aSubBasinOrder[p]=0;}//no downstream basin
      else                                                                {_aSubBasinOrder[p]=_aSubBasinOrder[pTo]+1;}
    }
    ordersum=0;
    for (p=0;p<_nSubBasins;p++)
    {
      ordersum+=_aSubBasinOrder[p];
      upperswap(_maxSubBasinOrder,_aSubBasinOrder[p]);
    }
    iter++;
  } while ((ordersum>last_ordersum) && (iter<max_iter));

  ExitGracefullyIf(iter>=max_iter,
                   "CModel::InitializeRoutingNetwork: exceeded maximum iterations. Circular reference in basin connections?",BAD_DATA);

  //this while loop should go through at most _maxSubBasinOrder times
  if (noisy){cout <<"      "<<iter<<" routing order iteration(s) completed"<<endl;}
  if (noisy){cout <<"      maximum subasin order: "<<_maxSubBasinOrder<<endl;}

  //starts at high order (leaf) basins, works its way down
  //generates _aOrderedSBind list, used in solver to order operations from
  //upstream to downstream
  //----------------------------------------------------------------------
  pp=0;
  int zerocount(0);
  _aOrderedSBind=new int [_nSubBasins];
  ExitGracefullyIf(_aOrderedSBind==NULL,"CModel::InitializeRoutingNetwork(2)",OUT_OF_MEMORY);
  for (ord=_maxSubBasinOrder;ord>=0;ord--)
  {
    if (noisy){cout<<"      order["<<ord<<"]:";}
    for (p=0;p<_nSubBasins;p++)
    {
      if (_aSubBasinOrder[p]==ord)
      {
        ExitGracefullyIf(pp>=_nSubBasins,"InitializeRoutingNetwork: fatal error",RUNTIME_ERR);
        ExitGracefullyIf(pp<0           ,"InitializeRoutingNetwork: fatal error",RUNTIME_ERR);
        _aOrderedSBind[pp]=p; pp++;
        if (noisy){cout<<" "<<p+1;}

        pTo=_aDownstreamInds[p];
        if (pTo==DOESNT_EXIST){zerocount++;}
      }
    }
    if (noisy){cout<<endl;}
  }
  if (noisy){cout <<"      number of zero-order outlets: "<<zerocount<<endl;}

  for (p = 0; p < _nSubBasins; p++)
  {
    if (_aSubBasinOrder[p] != _maxSubBasinOrder){
      _pSubBasins[p]->SetAsNonHeadwater();
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns ordered basin index, given a sub basin index
///
/// \param pp [in] Integer sub basin index
/// \return Integer index of ordered sub basin
//
int CModel::GetOrderedSubBasinIndex (const int pp) const
{
  ExitGracefullyIf((pp<0) || (pp>_nSubBasins),
                   "CModel::GetOrderedSubBasinIndex: invalid subbasin index",RUNTIME_ERR);
  return _aOrderedSBind[pp];
}

//////////////////////////////////////////////////////////////////
/// \brief Initializes basin flows
/// \details Calculates flow rates in all basins, propagates downstream;
///     Calculates all basin drainage areas and total watershed area
///     Called from CModel::Initialize AFTER initalize routing network
///     (i.e., directed stream network has already been constructed)
///
/// \param &Options [in] Global model options information
//
void CModel::InitializeBasinFlows(const optStruct &Options)
{
  int p;

  double *aSBArea=new double [_nSubBasins];//[km2] total drainage area of basin outlet
  double *aSBQin =new double [_nSubBasins];//[m3/s] avg. inflow to basin from upstream (0 for leaf)
  double *aSBQlat=NULL;
  aSBQlat=new double [_nSubBasins];//[m3/s] avg. lateral inflow within basin
  ExitGracefullyIf(aSBQlat==NULL,"CModel::InitializeBasinFlows",OUT_OF_MEMORY);

  //Estimate initial lateral runoff / flows in each basin (does not override .rvc values, if available)
  //-----------------------------------------------------------------
  double runoff_est=0.0;
  for (p=0;p<_nSubBasins;p++)
  {
    aSBArea[p]=_pSubBasins[p]->GetBasinArea();

    //runoff_est=EstimateInitialRunoff(p,Options);//[mm/d]

    if (CGlobalParams::GetParams()->avg_annual_runoff > 0){
      runoff_est= CGlobalParams::GetParams()->avg_annual_runoff/DAYS_PER_YEAR;
    }
    aSBQlat[p]=runoff_est/MM_PER_METER*(aSBArea[p]*M2_PER_KM2)/SEC_PER_DAY;//[m3/s]

    aSBQin [p]=_pSubBasins[p]->GetSpecifiedInflow(0.0);//initial conditions //[m3/s]

    //cout<<p<<": "<< aSBQin[p]<<" "<<runoff_est<<" "<<aSBQlat[p]<<" "<<aSBArea[p]<<endl;
  }

  //Propagate flows and drainage areas downstream
  //-----------------------------------------------------------------
  int pTo;
  bool warn=false;
  bool warn2=false;
  for (int pp=0;pp<_nSubBasins;pp++)
  { //calculated in order from upstream to downstream
    p=this->GetOrderedSubBasinIndex(pp);

    pTo=_aDownstreamInds[p];
    if (pTo!=DOESNT_EXIST){
      aSBQin [pTo]+=aSBQlat[p]+aSBQin[p];
      aSBArea[pTo]+=aSBArea[p];//now aSBArea==drainage area
    }
    //cout<<p<<": "<< aSBQin[p]<<" "<<aSBQlat[p]<<" "<<aSBArea[p]<<endl;
    if (_pSubBasins[p]->GetReferenceFlow() == AUTO_COMPUTE){warn =true;}
    if (_pSubBasins[p]->GetOutflowRate()   == AUTO_COMPUTE){warn2=true;}
    _pSubBasins[p]->Initialize(aSBQin[p],aSBQlat[p],aSBArea[p],Options);
  }
  if ((warn) && (_nSubBasins>1)){
    WriteWarning("CModel::InitializeBasinFlows: one or more subbasin reference discharges were autogenerated from annual average runoff", Options.noisy);
  }
  if ((warn2) && (_nSubBasins>1)){
    WriteWarning("CModel::InitializeBasinFlows: one or more subbasin initial outflows were autogenerated from annual average runoff", Options.noisy);
  }

  /*cout<<"Routing Diagnostics"<<endl;
    cout<<"ord pp pTo Qin Qlat Area"<<endl;
    for (p=0;p<_nSubBasins;p++){
    cout<< _aSubBasinOrder[p]<<" "<< _aOrderedSBind[p]<<" "<<aDownstreamInds[p]<<" ";
    cout <<aSBQin [p]<<" "<< aSBQlat[p]<<" "<<aSBArea[p]<<endl;
    }*/
  delete [] aSBQin;
  delete [] aSBArea;
  delete [] aSBQlat;
}

//////////////////////////////////////////////////////////////////
/// \brief Generates gauge weights
/// \details Populates an array _aGaugeWeights with interpolation weightings for distribution of gauge station data to HRUs
/// \remark Called after initialize routing orders
///
/// \param &Options [in] Global model options information
//
void CModel::GenerateGaugeWeights(const optStruct &Options)
{
  int k,g;
  //allocate memory
  _aGaugeWeights=new double *[_nHydroUnits];
  for (k=0;k<_nHydroUnits;k++){
    _aGaugeWeights[k]=NULL;
    _aGaugeWeights[k]=new double [_nGauges];
    ExitGracefullyIf(_aGaugeWeights[k]==NULL,"GenerateGaugeWeights",OUT_OF_MEMORY);
    for (g=0;g<_nGauges;g++){
      _aGaugeWeights[k][g]=0.0;
    }
  }
  location xyh,xyg;

  switch(Options.interpolation)
  {
  case(INTERP_NEAREST_NEIGHBOR)://---------------------------------------------
  {
    //w=1.0 for nearest gauge, 0.0 for all others
    double distmin,dist;
    int    g_min=0;
    for (k=0;k<_nHydroUnits;k++)
    {
      xyh=_pHydroUnits[k]->GetCentroid();
      g_min=0;
      distmin=ALMOST_INF;
      for (g=0;g<_nGauges;g++)
      {
        xyg=_pGauges[g]->GetLocation();

        dist=pow(xyh.UTM_x-xyg.UTM_x,2)+pow(xyh.UTM_y-xyg.UTM_y,2);
        if (dist<distmin){distmin=dist;g_min=g;}
        _aGaugeWeights[k][g]=0.0;
      }
      _aGaugeWeights[k][g_min]=1.0;
    }
    break;
  }
  case(INTERP_AVERAGE_ALL):                   //---------------------------------------------
  {
    for (k=0;k<_nHydroUnits;k++){
      for (g=0;g<_nGauges;g++){
        _aGaugeWeights[k][g]=1.0/(double)(_nGauges);
      }
    }
    break;
  }
  case(INTERP_INVERSE_DISTANCE):                      //---------------------------------------------
  {
    //wt_i = (1/r_i^2) / (sum{1/r_j^2})
    double dist;
    double denomsum;
    const double IDW_POWER=2.0;
    int atop_gauge(-1);
    for (k=0;k<_nHydroUnits;k++)
    {
      xyh=_pHydroUnits[k]->GetCentroid();
      atop_gauge=-1;
      denomsum=0;
      for (g=0;g<_nGauges;g++)
      {
        xyg=_pGauges[g]->GetLocation();
        dist=sqrt(pow(xyh.UTM_x-xyg.UTM_x,2)+pow(xyh.UTM_y-xyg.UTM_y,2));
        denomsum+=pow(dist,-IDW_POWER);
        if (dist<REAL_SMALL){atop_gauge=g;}//handles limiting case where weight= large number/large number
      }

      for (g=0;g<_nGauges;g++)
      {
        xyg=_pGauges[g]->GetLocation();
        dist=sqrt(pow(xyh.UTM_x-xyg.UTM_x,2)+pow(xyh.UTM_y-xyg.UTM_y,2));

        if (atop_gauge!=-1){_aGaugeWeights[k][g]=0.0;_aGaugeWeights[k][atop_gauge]=1.0;}
        else               {_aGaugeWeights[k][g]=pow(dist,-IDW_POWER)/denomsum;}
      }
    }

    break;
  }
  case(INTERP_INVERSE_DISTANCE_ELEVATION):                    //---------------------------------------------
  {
    //wt_i = (1/r_i^2) / (sum{1/r_j^2})
    double dist;
    double elevh,elevg;
    double denomsum;
    const double IDW_POWER=2.0;
    int atop_gauge(-1);
    for(k=0; k<_nHydroUnits; k++)
    {
      elevh=_pHydroUnits[k]->GetElevation();
      atop_gauge=-1;
      denomsum=0;
      for(g=0; g<_nGauges; g++)
      {
        elevg=_pGauges[g]->GetElevation();
        dist=abs(elevh-elevg);
        denomsum+=pow(dist,-IDW_POWER);
        if(dist<REAL_SMALL){ atop_gauge=g; }//handles limiting case where weight= large number/large number
      }

      for(g=0; g<_nGauges; g++)
      {
        elevg=_pGauges[g]->GetElevation();
        dist=abs(elevh-elevg);

        if(atop_gauge!=-1){ _aGaugeWeights[k][g]=0.0; _aGaugeWeights[k][atop_gauge]=1.0; }
        else              { _aGaugeWeights[k][g]=pow(dist,-IDW_POWER)/denomsum; }
      }
    }

    break;
  }
  case (INTERP_FROM_FILE):                    //---------------------------------------------
  {
    //format:
    //:GaugeWeightTable
    //  nGauges nHydroUnits
    //  v11 v12 v13 v14 ... v_1,nGauges
    //  ...
    //  vN1 vN2 vN3 vN4 ... v_N,nGauges
    //:EndGaugeWeightTable
    //ExitGracefullyIf no gauge file
    int   Len,line(0);
    char *s[MAXINPUTITEMS];
    ifstream INPUT;
    INPUT.open(Options.interp_file.c_str());
    if (INPUT.fail())
    {
      INPUT.close();
      string errString = "GenerateGaugeWeights:: Cannot find gauge weighting file "+Options.interp_file;
      ExitGracefully(errString.c_str(),BAD_DATA);
    }
    else
    {
      CParser *p=new CParser(INPUT,Options.interp_file,line);
      bool done(false);
      while (!done)
      {
        p->Tokenize(s,Len);
        if (IsComment(s[0],Len)){}
        else if (!strcmp(s[0],":GaugeWeightTable")){}
        else if (Len>=2){
          ExitGracefullyIf(s_to_i(s[0])!=_nGauges,
                           "GenerateGaugeWeights: the gauge weighting file has an improper number of gauges specified",BAD_DATA);
          ExitGracefullyIf(s_to_i(s[1])!=_nHydroUnits,
                           "GenerateGaugeWeights: the gauge weighting file has an improper number of HRUs specified",BAD_DATA);
          done=true;
        }
      }
      int junk;
      p->Parse2DArray_dbl(_aGaugeWeights,_nHydroUnits,_nGauges,junk);

      for (k=0;k<_nHydroUnits;k++){
        double sum=0;
        for (g=0;g<_nGauges;g++){
          sum+=_aGaugeWeights[k][g];
        }
        if(fabs(sum-1.0)>1e-4){
          ExitGracefully("GenerateGaugeWeights: INTERP_FROM_FILE: user-specified weights for gauge don't add up to 1.0",BAD_DATA);
        }
      }
      INPUT.close();
      delete p;
    }
    break;
  }
  default:
  {
    ExitGracefully("CModel::GenerateGaugeWeights: Invalid interpolation method",BAD_DATA);
  }
  }

  //check quality - weights for each HRU should add to 1
  double sum;
  for (k=0;k<_nHydroUnits;k++){
    sum=0.0;
    for (g=0;g<_nGauges;g++){
      sum+=_aGaugeWeights[k][g];
    }

    ExitGracefullyIf((fabs(sum-1.0)>REAL_SMALL) && (INTERP_FROM_FILE) && (_nGauges>1),
                     "GenerateGaugeWeights: Bad weighting scheme- weights for each HRU must sum to 1",BAD_DATA);
    ExitGracefullyIf((fabs(sum-1.0)>REAL_SMALL) && !(INTERP_FROM_FILE) && (_nGauges>1),
                     "GenerateGaugeWeights: Bad weighting scheme- weights for each HRU must sum to 1",RUNTIME_ERR);
  }
}
