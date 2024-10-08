/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2023 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Model.h"
#include "IrregularTimeSeries.h"
#include "HeatConduction.h"

/*****************************************************************
   Model Initialization Routines
------------------------------------------------------------------
   called prior to simulation
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Initializes model prior to simulation
/// \details Perform all operations required and initial check on validity of model before simulation begins;
///
///  -initializes mass balance arrays to zero; identifies model UTM zone;
///  -calls initialization routines of all subasins, HRUs, processes, and gauges
///  -generates gauge weights
///  -calculates stream network topology, routing orders
///  -calculates initial water/energy/mass storage
///  -calls routines to write initial conditions to minor output
///  -handles disabling of HRUs and subbasins in model for submodel calibration
///
/// \param &Options [in] Global model options information
//
void CModel::Initialize(const optStruct &Options)
{
  int g,i,j,k,kk,p,pp;

  // Quality control
  //--------------------------------------------------------------
  ExitGracefullyIf(_nSubBasins<1,
                    "CModel::Initialize: Must have at least one SubBasin",BAD_DATA);
  ExitGracefullyIf(_nHydroUnits<1,
                    "CModel::Initialize: Must have at least one hydrologic unit",BAD_DATA);
  ExitGracefullyIf(_nGauges<1 && _nForcingGrids<1 && (!Options.in_bmi_mode),
                    "CModel::Initialize: Must have at least one meteorological gauge station or forcing grid",BAD_DATA);
  ExitGracefullyIf(_nProcesses==0,
                    "CModel::Initialize: must have at least one hydrological process included in model",BAD_DATA);

  //Ensure Basins & HRU IDs are unique
  for (p=0;p<_nSubBasins;p++){
    for (pp=0;pp<p;pp++){
      if (_pSubBasins[p]->GetID()==_pSubBasins[pp]->GetID()){
        ExitGracefully("CModel::Initialize: non-unique (repeated) basin identifier found",BAD_DATA);}}}

  for (k=0;k<_nHydroUnits;k++){
    for (kk=0;kk<p;kk++){
      if ((k!=kk) && (_pHydroUnits[k]->GetHRUID()==_pHydroUnits[kk]->GetHRUID())){
        ExitGracefully("CModel::Initialize: non-unique (repeated) HRU identifier found",BAD_DATA);}}}

  // initialize process algorithms, initialize water/energy balance arrays to zero
  //--------------------------------------------------------------
  _nTotalConnections=0;
  for (j=0; j<_nProcesses;j++){
    if((_pProcesses[j]->GetProcessType()!=PRECIPITATION) &&
       (_pProcesses[j]->GetProcessType()!=LAT_ADVECTION))
       {_pProcesses[j]->Initialize();} //precip already initialized in ParseInput.cpp
    _nTotalConnections+=_pProcesses[j]->GetNumConnections();
  }
  if (_pTransModel->GetNumConstituents()>0){
    _pTransModel->CalculateLateralConnections();
    for (j=0; j<_nProcesses;j++){
      if (_pProcesses[j]->GetProcessType()==LAT_ADVECTION){
        _pProcesses[j]->Initialize(); //must be re-initialized after all lateral processes are initialized above
        _nTotalConnections+=_pProcesses[j]->GetNumConnections();
      }
    }
  }
  // reserve memory for mass balance arrays
  //--------------------------------------------------------------
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
  for (j=0; j<_nProcesses;j++){
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

  // Identify model UTM_zone for interpolation
  //--------------------------------------------------------------
  double cen_long(0),area_tot(0);//longiturde of area-weighted watershed centroid, total wshed area
  CHydroUnit *pHRU;
  for (k=0; k<_nHydroUnits;k++)
  {
    pHRU=_pHydroUnits[k];
    area_tot+=pHRU->GetArea();
    cen_long+=pHRU->GetCentroid().longitude*(pHRU->GetArea());
  }
  cen_long/=area_tot;
  _UTM_zone =(int)(floor ((cen_long + 180.0) / 6) + 1);

  // Initialize HRUs, gauges and transient parameters
  //--------------------------------------------------------------
  for (k=0;k<_nHydroUnits; k++ ){_pHydroUnits [k ]->Initialize(_UTM_zone);}
  for (g=0;g<_nGauges;     g++ ){_pGauges     [g ]->Initialize(Options,_UTM_zone);}
  for (j=0;j<_nTransParams;j++ ){_pTransParams[j ]->Initialize(this,Options);}
  for (kk=0;kk<_nHRUGroups;kk++){_pHRUGroups  [kk]->Initialize(); } //disables HRUs
  for (pp=0;pp<_nSBGroups; pp++){_pSBGroups   [pp]->Initialize(); } //disables SBs and HRUs

  InitializeParameterOverrides();

  // Forcing grids are not "Initialized" here because the derived data have to be populated everytime a new chunk is read

  // QA/QC Check for partial or full disabling of basin HRUs (after HRU group initialize, must be before area calculation)
  //---------------------------------------------------------------
  string disbasins="\n     ";
  bool anydisabled=false;
  int nDisBasins=0;
  for(p=0;p<_nSubBasins;p++){
    bool one_enabled =false;
    bool one_disabled=false;
    for(k=0; k<_pSubBasins[p]->GetNumHRUs();k++){
      if(_pSubBasins[p]->GetHRU(k)->IsEnabled()){one_enabled=true;}
      else                                      {one_disabled=true;}
    }
    if(one_enabled && one_disabled){
      string warn="CModel::Initialize: only some of the HRUs in subbasin "+to_string(_pSubBasins[p]->GetID())+" are disabled. It is suggested to disable all or none of the HRUs in one subbasin to avoid uninterpretable results";
      ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    }
    if(one_disabled){
      disbasins=disbasins+to_string(_pSubBasins[p]->GetID())+", "; anydisabled=true;nDisBasins++;
      if (nDisBasins % 40 == 39) {disbasins = disbasins + "\n     ";}
    }
    if((one_disabled) && (!one_enabled)){
      _pSubBasins[p]->Disable(); //all HRUs in basin disabled -disable subbasin!
    }
  }

  if(anydisabled){
    WriteWarning("CModel::Initialize: the following subbasins have one or more member HRUs disabled: "+disbasins,Options.noisy);
  }

  // Identify reservoir-linked HRUs and link them
  //--------------------------------------------------------------
  CReservoir *pRes;
  for(p=0;p<_nSubBasins;p++){
    pRes=_pSubBasins[p]->GetReservoir();
    if(pRes!=NULL){
      if(pRes->GetHRUIndex()!=DOESNT_EXIST){
        _pHydroUnits[pRes->GetHRUIndex()]->LinkToReservoir(_pSubBasins[p]->GetID());
      }
    }
  }

  // Generate Gauge Weights from Interpolation (except if in BMI mode and no RVT file is specified)
  //--------------------------------------------------------------
  if ((!Options.in_bmi_mode) || (strcmp(Options.rvt_filename.c_str(), "") != 0)) {

    if (!Options.silent){cout <<"  Generating Gauge Interpolation Weights..."<<endl;}
    forcing_type f_gauge=F_PRECIP;
    //search for the 'other' forcing
    if      (Options.air_pressure  ==AIRPRESS_DATA){f_gauge=F_AIR_PRES;}
    else if (Options.SW_radiation  ==SW_RAD_DATA  ){f_gauge=F_SW_RADIA;}
    else if (Options.evaporation   ==PET_DATA     ){f_gauge=F_PET; }
    else if (Options.rel_humidity  ==RELHUM_DATA  ){f_gauge=F_REL_HUMIDITY;}
    else if (Options.SW_radia_net  ==NETSWRAD_DATA){f_gauge=F_SW_RADIA_NET;}
    else if (Options.LW_radiation  ==LW_RAD_DATA  ){f_gauge=F_LW_RADIA_NET;}
    else if (Options.ow_evaporation==PET_DATA     ){f_gauge=F_OW_PET;}
    else if (Options.wind_velocity ==WINDVEL_DATA ){f_gauge=F_WIND_VEL;}
    if (Options.noisy){cout<<"     Gauge weights determined from "<<ForcingToString(f_gauge)<<" gauges"<<endl; }

    if(Options.write_interp_wts)
    { //creates and/or deletes existing InterpolationWeights.csv file
      ofstream WTS;
      string tmpFilename=FilenamePrepare("InterpolationWeights.csv",Options);
      WTS.open(tmpFilename.c_str());
      if(WTS.fail()) {
        ExitGracefully(("CModel::GenerateGaugeWeights: unable to open output file "+tmpFilename+" for writing.").c_str(),FILE_OPEN_ERR);
      }
      WTS.close();
    }

    GenerateGaugeWeights(_aGaugeWeights ,f_gauge   ,Options);//'other' forcings
    GenerateGaugeWeights(_aGaugeWtPrecip,F_PRECIP  ,Options);
    GenerateGaugeWeights(_aGaugeWtTemp  ,F_TEMP_AVE,Options);

  }

  //Initialize SubBasins, calculate routing orders, topology
  //--------------------------------------------------------------
  if (!Options.silent){cout<<"  Calculating basin & watershed areas..."<<endl;}
  _WatershedArea=0.0;
  for (p=0;p<_nSubBasins;p++){_WatershedArea+=_pSubBasins[p]->CalculateBasinArea();}

  if (!Options.silent){cout<<"  Calculating routing network topology..."<<endl;}
  InitializeRoutingNetwork(); //calculate proper routing orders

  if (!Options.silent){cout<<"  Initializing Basins, calculating watershed area, setting initial flow conditions..."<<endl;}
  InitializeBasins(Options,false);

  // Initialize Data Assimilation
  //--------------------------------------------------------------
  InitializeDataAssimilation(Options);

  // Initialize Groundwater
  //--------------------------------------------------------------
  if (Options.modeltype == MODELTYPE_COUPLED) {
    _pGWModel->Initialize(Options);
    // Ensure duration matches between models
    if (_pGWModel->GetTotalTSteps() != (int)(Options.duration/Options.timestep)) {
      ExitGracefully("CModel::Initialize: groundwater model total time steps not equal to Raven time steps", BAD_DATA);
    }
    // Eventually - more checks on time step compatability
  }

  // Initialize Transport
  //--------------------------------------------------------------
  if (_pTransModel->GetNumConstituents()>0){
    if (!Options.silent){cout<<"  Initializing Transport Model..."<<endl;}
    _pTransModel->Initialize(Options);
  }
  for (j = 0; j < _nProcesses; j++) {
    if (_pProcesses[j]->GetProcessType() == HEATCONDUCTION) {
      CmvHeatConduction *pHC=static_cast<CmvHeatConduction *>(_pProcesses[j]);
      pHC->StoreNumberOfHRUs(GetNumHRUs());
    }
  }

  // Precalculate whether individual processes should apply (for speed)
  //--------------------------------------------------------------
  _aShouldApplyProcess = new bool *[_nProcesses];
  for (j=0; j<_nProcesses;j++){
    _aShouldApplyProcess[j]=NULL;
    _aShouldApplyProcess[j] = new bool [_nHydroUnits];
    ExitGracefullyIf(_aShouldApplyProcess[j]==NULL,"CModel::Initialize (_aShouldApplyProcess)",OUT_OF_MEMORY);
    for (k=0; k<_nHydroUnits;k++){
      _aShouldApplyProcess[j][k] = _pProcesses[j]->ShouldApply(_pHydroUnits[k]);
      if (!_pHydroUnits[k]->IsEnabled()){ _aShouldApplyProcess[j][k] =false;}
    }
  }

  // Initialize NetCDF Output File IDs
  //--------------------------------------------------------------
  /* initialize all potential NetCDF file IDs with -9 == "not existing and hence not opened" */
  _HYDRO_ncid    = -9;   // output file ID for Hydrographs.nc         (-9 --> not opened)
  _STORAGE_ncid  = -9;   // output file ID for WatershedStorage.nc    (-9 --> not opened)
  _FORCINGS_ncid = -9;   // output file ID for ForcingFunctions.nc    (-9 --> not opened)
  _RESSTAGE_ncid = -9;   // output file ID for ReservoirStages.nc     (-9 --> not opened)
  _RESMB_ncid    = -9;

  //Initialize Custom Output
  //--------------------------------------------------------------
  for (int c=0;c<_nCustomOutputs;c++){
    _pCustomOutputs[c]->InitializeCustomOutput(Options);
  }

  quickSort(_aOutputTimes,0,_nOutputTimes-1);

  //Prepare Output Time Series
  //--------------------------------------------------------------
  InitializeObservations(Options);

  // Generate default diagnostic period - entire simulation
  //--------------------------------------------------------------
  CDiagPeriod *pDP=new CDiagPeriod("ALL","0001-01-01","9999-12-31",COMPARE_GREATERTHAN,-ALMOST_INF,Options);
  AddDiagnosticPeriod(pDP);

  //General QA/QC
  //--------------------------------------------------------------
  ExitGracefullyIf((GetNumGauges()<2) && (Options.orocorr_temp==OROCORR_UBCWM2),
                   "CModel::Initialize: at least 2 gauges necessary to use :OroTempCorrect method OROCORR_UBCWM2", BAD_DATA);

  //--Warn about interception handling
  if(!StateVarExists(CANOPY_SNOW) && (Options.interception_factor!=PRECIP_ICEPT_NONE)){
    WriteAdvisory("Because no processes with CANOPY_SNOW variable have been specified, all snow interception will be directly moved to the atmosphere as if it had sublimated.",Options.noisy);
  }
  if(!StateVarExists(CANOPY) && (Options.interception_factor!=PRECIP_ICEPT_NONE)) {
    WriteAdvisory("Because no processes with CANOPY variable have been specified, all rain interception will be directly moved to the atmosphere as if it had evaporated.",Options.noisy);
  }
  //--Check for empty HRU groups
  for (kk = 0; kk < _nHRUGroups; kk++){
    if (_pHRUGroups[kk]->GetNumHRUs() == 0){
      string warn = "CModel::Initialize: HRU Group " + _pHRUGroups[kk]->GetName() + " is empty.";
      WriteAdvisory(warn,Options.noisy);
    }
  }
  //--Check for repeated HRU groups
  for (kk = 0; kk < _nHRUGroups; kk++){
    for (int kkk = 0; kkk < kk; kkk++){
      if(StringToUppercase(_pHRUGroups[kk]->GetName()) == StringToUppercase(_pHRUGroups[kkk]->GetName())){
        string warn = "CModel::Initialize: Name of HRU Group " + _pHRUGroups[kk]->GetName() + " is repeated.";
        WriteWarning(warn,Options.noisy);
      }
    }
  }
  //--check for stage or area observations not linked to valid reservoir
  for(i=0; i<_nObservedTS; i++){
    if(!strcmp(_pObservedTS[i]->GetName().c_str(),"RESERVOIR_STAGE"))
    {
      long SBID=_pObservedTS[i]->GetLocID();
      if(GetSubBasinByID(SBID)->GetReservoir()==NULL){
        string warn="Stage observations supplied for non-existent reservoir in subbasin "+to_string(SBID);
        ExitGracefully(warn.c_str(),BAD_DATA_WARN);
      }
    }
    else if(!strcmp(_pObservedTS[i]->GetName().c_str(),"LAKE_AREA"))
    {
      long SBID=_pObservedTS[i]->GetLocID();
      if(GetSubBasinByID(SBID)->GetReservoir()==NULL){
        string warn="Lake area bservations supplied for non-existent reservoir in subbasin "+to_string(SBID);
        ExitGracefully(warn.c_str(),BAD_DATA_WARN);
      }
    }
  }
  //--check for inflow observations not linked to valid reservoir
  for (i = 0; i<_nObservedTS; i++) {
	  if (!strcmp(_pObservedTS[i]->GetName().c_str(), "RESERVOIR_INFLOW"))
	  {
		  long SBID =_pObservedTS[i]->GetLocID();
		  if (GetSubBasinByID(SBID)->GetReservoir() == NULL) {
			  string warn = "Inflow observations supplied for non-existent reservoir in subbasin " + to_string(SBID);
			  ExitGracefully(warn.c_str(), BAD_DATA_WARN);
		  }
	  }
  }
  //--check for net inflow observations not linked to valid reservoir
  for (i = 0; i<_nObservedTS; i++) {
	  if (!strcmp(_pObservedTS[i]->GetName().c_str(), "RESERVOIR_NETINFLOW"))
	  {
		  long SBID = _pObservedTS[i]->GetLocID();
		  if (GetSubBasinByID(SBID)->GetReservoir() == NULL) {
			  string warn = "Net inflow Observations supplied for non-existent reservoir in subbasin " + to_string(SBID);
			  ExitGracefully(warn.c_str(), BAD_DATA);
		  }
	  }
  }
  //--check for wetlands in model without depression storage
  bool wetlandsinmodel=false;
  for(k=0;k<_nHydroUnits;k++) {
    if (_pHydroUnits[k]->GetHRUType()==HRU_WETLAND){wetlandsinmodel=true;break;}
  }
  if((wetlandsinmodel) && (GetStateVarIndex(DEPRESSION)==DOESNT_EXIST)) {
    ExitGracefully("CModel::Initialize: At least one WETLAND-type soil profile is included but no DEPRESSION storage processes included in hydrologic process list.",BAD_DATA_WARN);
  }
  //--Check for non-standard calendar with observations
  if((Options.calendar!=CALENDAR_PROLEPTIC_GREGORIAN) && (_nObservedTS!=0)){
    WriteWarning("CModelInitialize: if a non-standard calendar is used, all observation data must be converted into the same calendar format",Options.noisy);
  }
  if(Options.calendar!=CALENDAR_PROLEPTIC_GREGORIAN){
    WriteWarning("CModelInitialize: if a non-standard calendar is used, care must be taken with forcing data. All gauge forcing data must use the same calendar convention. Mixed calendars are only supported for NetCDF forcing inputs.",Options.noisy);
  }
  //--Check for non-midnight start time
  if((Options.julian_start_day-floor(Options.julian_start_day+TIME_CORRECTION)>REAL_SMALL)){
    WriteAdvisory("CModelInitialize: if non-midnight start time is used and temperature forcings begin after midnight, average/max/min daily temperature forcings will be incorrect on the first (and often last) day of simulation.",Options.noisy);
    WriteAdvisory("                  It is suggested to either (a) use forcings which start on or prior to 00:00 of the model start date or (b) avoid algorithms which require these daily temperature inputs",Options.noisy);
  }
  //--Check for non-midnight start time
  double rem_tsteps= (1.0 - floor(Options.julian_start_day + TIME_CORRECTION))/ Options.timestep;
  if ((floor(rem_tsteps + TIME_CORRECTION)-rem_tsteps)>REAL_SMALL){
    WriteWarning("CModelInitialize: the model time step and model start time is such that midnight does not correspond to a time step ending. This will cause issues with use of daily temperature forcings (and potentially other errors) throughout the simulation.", Options.noisy);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Initializes SB demand members AFTER RVM FILE READ
/// \param &Options [in] Global model options information
//
void CModel::InitializePostRVM(const optStruct& Options)
{
  for (int p=0;p<_nSubBasins;p++){
    _pSubBasins[p]->InitializePostRVM(Options);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates initial total system water storage, updates _initWater
///
/// \param &Options [in] Global model options information
//
void CModel::CalculateInitialWaterStorage(const optStruct &Options)
{
  if (!Options.silent){cout<<"  Calculating initial system water storage..."<<endl;}
  _initWater=0.0;
  double S=0;
  for (int i=0;i<_nStateVars;i++)
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
      _initWater-=_pSubBasins[p]->GetIntegratedOutflow        (Options.timestep)/2.0/_WatershedArea*MM_PER_METER/M2_PER_KM2;
    }
  }
  if(_pTransModel!=NULL) {
    for(int c=0;c<_pTransModel->GetNumConstituents();c++) {
      _pTransModel->GetConstituentModel(c)->CalculateInitialMassStorage(Options);
    }
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
    _pModeledTS[i] = new CTimeSeries("MODELED" + _pObservedTS[i]->GetName(),
                                      _pObservedTS[i]->GetLocID(),"",
                                      Options.julian_start_day,
                                      Options.julian_start_year,
                                      Options.timestep,nModeledValues,true);

    _pObservedTS[i]->Initialize(Options.julian_start_day, Options.julian_start_year, Options.duration, Options.timestep,true,Options.calendar);
    _pModeledTS [i]->InitializeResample(nModeledValues,Options.timestep);
    _aObsIndex  [i]=0;

    //Match weights with observations based on Name, tag and numValues
    tmp[i] = NULL;
    for (int n = 0; n < _nObsWeightTS; n++){

      if (    _pObsWeightTS[n]!=NULL
           && _pObsWeightTS[n]->GetName()     == _pObservedTS[i]->GetName()
           && _pObsWeightTS[n]->GetLocID()    == _pObservedTS[i]->GetLocID()
           && _pObsWeightTS[n]->GetNumValues()== _pObservedTS[i]->GetNumValues()  )
      {
        tmp[i] = _pObsWeightTS[n];
        _pObsWeightTS[n] = NULL;
        tmp[i]->Initialize(Options.julian_start_day, Options.julian_start_year, Options.duration,Options.timestep,true,Options.calendar);
      }
    }
  }

  //clean up and warn about unmatched weights
  for (int n = 0; n < _nObsWeightTS; n++){
    if (_pObsWeightTS[n] != NULL){
      WriteWarning("Observation Weight "+_pObsWeightTS[n]->GetName()+" "+to_string(_pObsWeightTS[n]->GetLocID())+" not matched to observation time series. Please supply corresponding observation series.", Options.noisy);
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

  bool noisy=false;//useful for debugging

  _aDownstreamInds=NULL;
  _aSubBasinOrder =new int [_nSubBasins];
  _aDownstreamInds=new int [_nSubBasins];
  int *aInflowCount=new int [_nSubBasins];
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
    aInflowCount    [p]=0;
  }
  // count number of inflows to each basin
  //----------------------------------------------------------------------
  for (p=0;p<_nSubBasins;p++)
  {
    pTo=_aDownstreamInds[p];
    if (pTo!=DOESNT_EXIST){aInflowCount[pTo]++;}
  }
  // Iterative identification of subbasin orders
  //----------------------------------------------------------------------
  // here, order goes from 0 (trunk) to _maxSubBasinOrder (furthest leaf). This is NOT Strahler ordering!!!
  const int MAX_ITER=1000;
  int last_ordersum;
  int iter(0),ordersum(0);
  do
  {
    last_ordersum=ordersum;
    for (p=0;p<_nSubBasins;p++)
    {
      pTo=_aDownstreamInds[p];
      if (pTo==DOESNT_EXIST){_aSubBasinOrder[p]=0;}//no downstream basin
      else                  {_aSubBasinOrder[p]=_aSubBasinOrder[pTo]+1;}
    }
    ordersum=0;
    for (p=0;p<_nSubBasins;p++)
    {
      ordersum+=_aSubBasinOrder[p];
      upperswap(_maxSubBasinOrder,_aSubBasinOrder[p]);
    }
    iter++;
  } while ((ordersum>last_ordersum) && (iter<MAX_ITER));

  ExitGracefullyIf(iter>=MAX_ITER,
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
    // identify headwater basins - default is headwater
    if (aInflowCount[p] > 0) {
      _pSubBasins[p]->SetAsNonHeadwater();
    }

    // set subbasin pointers to their downstream subbasins
    pTo=_aDownstreamInds[p];
    if (pTo!=DOESNT_EXIST){
      _pSubBasins[p]->SetDownstreamBasin(_pSubBasins[pTo]);
    }

    for(int i=0;i<_pSubBasins[p]->GetNumDiversions();i++) {
      int pdown=_pSubBasins[p]->GetDiversionTargetIndex(i);
      if((pdown!=DOESNT_EXIST) && (pdown>0) && (pdown<_nSubBasins)) {
        _pSubBasins[pdown]->SetAsNonHeadwater();
      }
    }
  }
  delete [] aInflowCount;
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
void CModel::InitializeBasins(const optStruct &Options,const bool re_init)
{
  int p;

  double *aSBArea=new double [_nSubBasins];//[km2] total drainage area of basin outlet
  double *aSBQin =new double [_nSubBasins];//[m3/s] avg. inflow to basin from upstream (0 for leaf)
  double *aSBQlat=NULL;
  aSBQlat=new double [_nSubBasins];//[m3/s] avg. lateral inflow within basin
  ExitGracefullyIf(aSBQlat==NULL,"CModel::InitializeBasins",OUT_OF_MEMORY);

  //Estimate initial lateral runoff / flows in each basin (does not override .rvc values, if available)
  //-----------------------------------------------------------------
  double runoff_est=0.0;
  for (p=0;p<_nSubBasins;p++)
  {
    aSBArea[p]=_pSubBasins[p]->GetBasinArea();

    if (this->_pGlobalParams->GetParams()->avg_annual_runoff > 0){
      runoff_est= this->_pGlobalParams->GetParams()->avg_annual_runoff/DAYS_PER_YEAR;
    }
    aSBQlat[p]=runoff_est/MM_PER_METER*(aSBArea[p]*M2_PER_KM2)/SEC_PER_DAY;//[m3/s]
    aSBQlat[p]+=_pSubBasins[p]->GetDownstreamInflow(0.0);//[m3/s]

    aSBQin [p]=_pSubBasins[p]->GetSpecifiedInflow(0.0);//initial conditions //[m3/s]
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
    if ((pTo!=DOESNT_EXIST) && (_pSubBasins[pTo]->IsEnabled())){
      aSBQin [pTo]+=aSBQlat[p]+aSBQin[p];
      aSBArea[pTo]+=aSBArea[p];//now aSBArea==drainage area
    }

    if (_pSubBasins[p]->GetReferenceFlow() == AUTO_COMPUTE){warn =true;}
    if (_pSubBasins[p]->GetOutflowRate()   == AUTO_COMPUTE){warn2=true;}
    if(!re_init) {
      _pSubBasins[p]->Initialize          (aSBQin[p],aSBQlat[p],aSBArea[p],Options);
    }
    else {
      _pSubBasins[p]->InitializeFlowStates(aSBQin[p],aSBQlat[p],Options); //only re-calculates initial flows
    }
  }
  if ((warn) && (_nSubBasins>1)){
    WriteAdvisory("CModel::InitializeBasins: one or more subbasin reference discharges were autogenerated from annual average runoff", Options.noisy);
  }
  if ((warn2) && (_nSubBasins>1)){
    //no longer makes sense in that rvc file is called after model initialization
    //WriteAdvisory("CModel::InitializeBasins: one or more subbasin initial outflows were autogenerated from annual average runoff", Options.noisy);
  }

  /*cout<<"Routing Diagnostics"<<endl;
    cout<<"ord pp pTo Qin Qlat Area"<<endl;
    for (p=0;p<_nSubBasins;p++){
    cout<< _aSubBasinOrder[p]<<" "<< _aOrderedSBind[p]<<" "<<aDownstreamInds[p]<<" ";
    cout <<aSBQin [p]<<" "<< aSBQlat[p]<<" "<<aSBArea[p]<<endl;
    }*/ //Keep for debugging
  delete [] aSBQin;
  delete [] aSBArea;
  delete [] aSBQlat;
}

//////////////////////////////////////////////////////////////////
/// \brief clears all time series data for re-read of .rvt file
/// \remark Called only in ensemble mode
///
/// \param &Options [in] Global model options information
//
void CModel::ClearTimeSeriesData(const optStruct& Options)
{
  if(DESTRUCTOR_DEBUG) { cout<<"DELETING RVT DATA"<<endl; }
  int c,f,g,i,j,k,p;
  for (g=0;g<_nGauges;       g++){delete _pGauges       [g];} delete [] _pGauges;       _pGauges=NULL; _nGauges=0;
  for (f=0;f<_nForcingGrids; f++){delete _pForcingGrids [f];} delete [] _pForcingGrids; _pForcingGrids=NULL; _nForcingGrids=0;
  for (i=0;i<_nObservedTS;   i++){delete _pObservedTS   [i];} delete [] _pObservedTS;   _pObservedTS=NULL;
  if (_pModeledTS != NULL){
    for (i = 0; i < _nObservedTS; i++){ delete _pModeledTS[i]; } delete[] _pModeledTS;    _pModeledTS = NULL;
  }
  _nObservedTS=0;
  for (i=0;i<_nObsWeightTS;  i++){delete _pObsWeightTS  [i];} delete [] _pObsWeightTS;  _pObsWeightTS=NULL; _nObsWeightTS;

  if (_aGaugeWeights!=NULL){
    for (k=0;k<_nHydroUnits;   k++){delete [] _aGaugeWeights[k]; } delete [] _aGaugeWeights;  _aGaugeWeights=NULL;
  }
  if (_aGaugeWtPrecip!=NULL){
    for (k=0;k<_nHydroUnits;   k++){delete [] _aGaugeWtPrecip[k]; } delete [] _aGaugeWtPrecip;  _aGaugeWtPrecip=NULL;
  }
  if (_aGaugeWtTemp!=NULL){
    for (k=0;k<_nHydroUnits;   k++){delete [] _aGaugeWtTemp[k]; } delete [] _aGaugeWtTemp;  _aGaugeWtTemp=NULL;
  }
  for (j=0;j<_nTransParams;j++) {delete _pTransParams[j];} delete [] _pTransParams; _pTransParams=NULL; _nTransParams=0;
  for (j=0;j<_nClassChanges;j++){delete _pClassChanges[j];} delete [] _pClassChanges; _pClassChanges=NULL; _nClassChanges=0;

  delete [] _aObsIndex;      _aObsIndex=NULL;

  for(p=0;p<_nSubBasins;p++) {
    _pSubBasins[p]->ClearTimeSeriesData(Options);
  }
  for(c=0;c<_pTransModel->GetNumConstituents();c++) {
    _pTransModel->GetConstituentModel(c)->ClearTimeSeriesData(Options);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief initializes all paramter override structures
///
//
void CModel::InitializeParameterOverrides()
{
  for (int i=0;i<_nParamOverrides;i++)
  {
    //determines memory address location of the appropriate global variable
    // this is the only time we mess with addresses like this in the code - merely for speed
    string name=_pParamOverrides[i]->param_name;
    if      (name == "PET_BLEND_WTS") { //special case
      _pParamOverrides[i]->pxAddress=_PETBlends_wts;
      for (int j=0;j<_PETBlends_N;j++){
        _pParamOverrides[i]->aRevertValues[j] =_PETBlends_wts[j];
      }
    }
    else if (name == "POTMELT_BLEND_WTS") {
      _pParamOverrides[i]->pxAddress=_PotMeltBlends_wts;
      for (int j=0;j<_PotMeltBlends_N;j++){
        _pParamOverrides[i]->aRevertValues[j] =_PotMeltBlends_wts[j];
      }
    }
    else
    {
      double *address = this->_pGlobalParams->GetAddress(name);
      if (address==NULL){
        ExitGracefully("CModel::InitializeParameterOverrides() : Invalid or unsupported global parameter name in :LocalParameterOverride command",BAD_DATA);
      }
      _pParamOverrides[i]->pxAddress=address;

      _pParamOverrides[i]->aRevertValues[0] = this->_pGlobalParams->GetParameter(name);
    }

    _pParamOverrides[i]->aHRUIsOverridden=new bool [_nHydroUnits];
    ExitGracefullyIf(_pParamOverrides[i]->aHRUIsOverridden==NULL,"InitializeParameterOverrides()",OUT_OF_MEMORY);

    // Search for SB Group - this is done here so that group can be read in .rvp without knowledge of SB groups
    int pp=DOESNT_EXIST;
    string sbg_name=_pParamOverrides[i]->SBGroup_name;
    for (int j=0;j<_nSBGroups;j++){
      if (_pSBGroups[j]->GetName()==sbg_name){pp=j; break;}
    }
    if (pp==DOESNT_EXIST){
      ExitGracefully("CModel::InitializeParameterOverrides() : Invalid subbasin group in :LocalParameterOverride command",BAD_DATA);
    }

    // determine HRUs in which this applies
    for (int k=0;k<_nHydroUnits;k++){
      long SBID=_pSubBasins[_pHydroUnits[k]->GetSubBasinIndex()]->GetID();
       _pParamOverrides[i]->aHRUIsOverridden[k]=(_pSBGroups[pp]->IsInGroup(SBID));
    }
  }
  if (_nParamOverrides > 0) {
    string warn="One or more global parameter overrides are used in this model. These can have high computational overhead and should be used sparingly";
    WriteWarning(warn.c_str(),true);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Generates gauge weights
/// \details Populates an array aWts with interpolation weightings for distribution of gauge station data to HRUs
/// \remark Called after initialize routing orders
///
/// \param aWts [out] weights matrix [nHRUs][nGauges]
/// \param forcing [int] forcing type (F_PRECIP or F_TEMP)
/// \param &Options [in] Global model options information
//
void CModel::GenerateGaugeWeights(double **&aWts, const forcing_type forcing, const optStruct &Options)
{
  int k,g;
  bool *has_data=NULL;
  location xyh,xyg;

  //allocate memory
  aWts=NULL;
  aWts=new double *[_nHydroUnits];
  ExitGracefullyIf(aWts==NULL,"GenerateGaugeWeights",OUT_OF_MEMORY);

  for (k=0;k<_nHydroUnits;k++){
    aWts[k]=NULL;
    aWts[k]=new double [_nGauges];
    ExitGracefullyIf(aWts[k]==NULL,"GenerateGaugeWeights(2)",OUT_OF_MEMORY);
    for (g=0;g<_nGauges;g++){
      aWts[k][g]=0.0;
    }
  }

  int nGaugesWithData=0;
  has_data=new bool [_nGauges];
  ExitGracefullyIf(has_data==NULL,"GenerateGaugeWeights(3)",OUT_OF_MEMORY);
  for(g=0;g<_nGauges;g++){
    has_data[g]=_pGauges[g]->TimeSeriesExists(forcing);
    if(has_data[g]){ nGaugesWithData++; }
  }

  //handle the case that weights are allowed to sum to zero -netCDF is available
  //still need to allocate all zeros
  if (ForcingGridIsAvailable(forcing)){ return; }
  if ((forcing==F_TEMP_AVE) && (ForcingGridIsAvailable(F_TEMP_DAILY_MIN))){return;} //this is also acceptable
  if ((forcing==F_TEMP_AVE) && (ForcingGridIsAvailable(F_TEMP_DAILY_AVE))){return;} //this is also acceptable
  if ((forcing==F_PRECIP  ) && (ForcingGridIsAvailable(F_RAINFALL))){return;} //this is also acceptable

  string warn="GenerateGaugeWeights: no gauges present with the following data: "+ForcingToString(forcing);
  ExitGracefullyIf(nGaugesWithData==0,warn.c_str(),BAD_DATA_WARN);

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
        if(has_data[g]){
          xyg=_pGauges[g]->GetLocation();
          dist=pow(xyh.UTM_x-xyg.UTM_x,2)+pow(xyh.UTM_y-xyg.UTM_y,2);
          if(dist<distmin){ distmin=dist;g_min=g; }
        }
        aWts[k][g]=0.0;
      }
      aWts[k][g_min]=1.0;

    }
    break;
  }
  case(INTERP_AVERAGE_ALL):                   //---------------------------------------------
  {
    for (k=0;k<_nHydroUnits;k++){
      for (g=0;g<_nGauges;g++){
        if(has_data[g]){aWts[k][g]=1.0/(double)(nGaugesWithData);}
        else           {aWts[k][g]=0.0;}
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
    int atop_gauge(DOESNT_EXIST);
    for (k=0;k<_nHydroUnits;k++)
    {
      xyh=_pHydroUnits[k]->GetCentroid();
      atop_gauge=DOESNT_EXIST;
      denomsum=0;
      for (g=0;g<_nGauges;g++)
      {
        if(has_data[g]){
          xyg=_pGauges[g]->GetLocation();
          dist=sqrt(pow(xyh.UTM_x-xyg.UTM_x,2)+pow(xyh.UTM_y-xyg.UTM_y,2));
          denomsum+=pow(dist,-IDW_POWER);
          if(dist<REAL_SMALL){ atop_gauge=g; }//handles limiting case where weight= large number/large number
        }
      }

      for (g=0;g<_nGauges;g++)
      {
        aWts[k][g]=0.0;
        if(has_data[g]){
          xyg=_pGauges[g]->GetLocation();
          dist=sqrt(pow(xyh.UTM_x-xyg.UTM_x,2)+pow(xyh.UTM_y-xyg.UTM_y,2));

          if(atop_gauge!=DOESNT_EXIST){ aWts[k][g]=0.0;aWts[k][atop_gauge]=1.0; }
          else                        { aWts[k][g]=pow(dist,-IDW_POWER)/denomsum;         }
        }
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
    int atop_gauge(DOESNT_EXIST);
    for(k=0; k<_nHydroUnits; k++)
    {
      elevh=_pHydroUnits[k]->GetElevation();
      atop_gauge=DOESNT_EXIST;
      denomsum=0;
      for(g=0; g<_nGauges; g++)
      {
        if(has_data[g]){
          elevg=_pGauges[g]->GetElevation();
          dist=abs(elevh-elevg);
          denomsum+=pow(dist,-IDW_POWER);
          if(dist<REAL_SMALL){ atop_gauge=g; }//handles limiting case where weight= large number/large number
        }
      }

      for(g=0; g<_nGauges; g++)
      {
        aWts[k][g]=0.0;
        if(has_data[g]){
          elevg=_pGauges[g]->GetElevation();
          dist=abs(elevh-elevg);
          if(atop_gauge!=DOESNT_EXIST){ aWts[k][g]=0.0; aWts[k][atop_gauge]=1.0; }
          else                        { aWts[k][g]=pow(dist,-IDW_POWER)/denomsum; }
        }
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
      p->Parse2DArray_dbl(aWts,_nHydroUnits,_nGauges,junk);

      for (k=0;k<_nHydroUnits;k++){
        double sum=0;
        for (g=0;g<_nGauges;g++){
          sum+=aWts[k][g];
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

  //Override weights where specified
  for (k = 0; k < _nHydroUnits; k++) {
    if (_pHydroUnits[k]->GetSpecifiedGaugeIndex() != DOESNT_EXIST) {
      for (g=0;g<_nGauges;g++){
        aWts[k][g]=0.0;
      }
      g=_pHydroUnits[k]->GetSpecifiedGaugeIndex();
      aWts[k][g]=1.0;
    }
  }

  //check quality - weights for each HRU should add to 1
  double sum;
  for (k=0;k<_nHydroUnits;k++){
    sum=0.0;
    for (g=0;g<_nGauges;g++){sum+=aWts[k][g];}

    ExitGracefullyIf((fabs(sum-1.0)>REAL_SMALL) && (INTERP_FROM_FILE) && (_nGauges>1),
                     "GenerateGaugeWeights: Bad weighting scheme- weights for each HRU must sum to 1",BAD_DATA);
    ExitGracefullyIf((fabs(sum-1.0)>REAL_SMALL) && !(INTERP_FROM_FILE) && (_nGauges>1),
                     "GenerateGaugeWeights: Bad weighting scheme- weights for each HRU must sum to 1",RUNTIME_ERR);
  }

  if(Options.write_interp_wts)
  {
    ofstream WTS;
    string tmpFilename=FilenamePrepare("InterpolationWeights.csv",Options);
    WTS.open(tmpFilename.c_str(),ios::app);
    WTS<<"Weights for "<<ForcingToString(forcing)<<"--------------------------------------"<<endl;
    WTS<<"HRU index (k), HRU ID";
    for(g=0; g<_nGauges;g++) { WTS<<","<<_pGauges[g]->GetName(); }
    WTS<<endl;
    for(k=0;k<_nHydroUnits;k++) {
      WTS<<k<<","<<_pHydroUnits[k]->GetHRUID();

      for(g=0;g<_nGauges;g++) {WTS<<","<<aWts[k][g]; }
      WTS<<endl;
    }
    WTS.close();
  }

  delete[] has_data;
}
