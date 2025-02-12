/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2024 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Model.h"
#include "EnergyTransport.h"

/*****************************************************************
   Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the Model constructor
///
/// \param SM [in] Input soil model object
/// \param nsoillayers [in] Integer number of soil layers
//
CModel::CModel(const int        nsoillayers,
               const optStruct &Options)
{
  int i;
  _nSubBasins=0;      _pSubBasins=NULL;
  _nHydroUnits=0;     _pHydroUnits=NULL;
  _nHRUGroups=0;      _pHRUGroups=NULL;
  _nSBGroups=0;       _pSBGroups=NULL;
  _nGauges=0;         _pGauges=NULL;
  _nForcingGrids=0;   _pForcingGrids=NULL;
  _nProcesses=0;      _pProcesses=NULL;
  _nCustomOutputs=0;  _pCustomOutputs=NULL;
  _nTransParams=0;    _pTransParams=NULL;
  _nClassChanges=0;   _pClassChanges=NULL;
  _nParamOverrides=0; _pParamOverrides=NULL;
  _nObservedTS=0;     _pObservedTS=NULL; _pModeledTS=NULL; _aObsIndex=NULL;
  _nObsWeightTS =0;   _pObsWeightTS=NULL;
  _nDiagnostics=0;    _pDiagnostics=NULL;
  _nDiagPeriods=0;    _pDiagPeriods=NULL;
  _nAggDiagnostics=0; _pAggDiagnostics=NULL;
  _nPerturbations=0;  _pPerturbations=NULL;

  _nLandUseClasses=0;       _pLandUseClasses=NULL;
  _nAllSoilClasses = 0;     _pAllSoilClasses = NULL;
  _numVegClasses = 0;       _pAllVegClasses  = NULL;
  _nAllTerrainClasses = 0;  _pAllTerrainClasses = NULL;
  _nAllSoilProfiles = 0;    _pAllSoilProfiles = NULL;
  _nAllChannelXSects = 0;   _pAllChannelXSects = NULL;
  _nConvVariables = 0;

  _nTotalConnections=0;
  _nTotalLatConnections=0;

  _WatershedArea=0;

  _aSubBasinOrder =NULL; _maxSubBasinOrder=0;
  _aOrderedSBind  =NULL;
  _aDownstreamInds=NULL;

  _aDAscale       =NULL; //Initialized in InitializeDataAssimilation
  _aDAscale_last  =NULL;
  _aDAQadjust     =NULL;
  _aDADrainSum    =NULL;
  _aDADownSum     =NULL;
  _aDAlength      =NULL;
  _aDAtimesince   =NULL;
  _aDAoverride    =NULL;
  _aDAobsQ        =NULL;

  _pOptStruct = &Options;
  _pGlobalParams = new CGlobalParams();

  _HYDRO_ncid   =-9;
  _STORAGE_ncid =-9;
  _FORCINGS_ncid=-9;
  _RESSTAGE_ncid=-9;
  _RESMB_ncid   =-9;

  _PETBlends_N=0;
  _PETBlends_type=NULL;
  _PETBlends_wts=NULL;
  _PotMeltBlends_N=0;
  _PotMeltBlends_type=NULL;
  _PotMeltBlends_wts=NULL;

  ExitGracefullyIf(nsoillayers<1,
                   "CModel constructor::improper number of soil layers. SoilModel not specified?",BAD_DATA);

  //Initialize Lookup table for state variable indices
  for (i=0;i< MAX_STATE_VAR_TYPES;i++){
    for (int m=0;m<MAX_SV_LAYERS;m++){
      _aStateVarIndices[i][m]=DOESNT_EXIST;
    }
  }
  //determine first group of state variables based upon soil model
  //SW, atmosphere, atmos_precip always present, one for each soil layer, and 1 for GW (unless lumped)
  _nStateVars=5+nsoillayers;

  _aStateVarType =new sv_type [_nStateVars];
  _aStateVarLayer=new int     [_nStateVars];
  _nSoilVars     =nsoillayers;

  _aStateVarType[0]=SURFACE_WATER; _aStateVarLayer[0]=DOESNT_EXIST; _aStateVarIndices[(int)(SURFACE_WATER)][0]=0;
  _aStateVarType[1]=ATMOSPHERE;    _aStateVarLayer[1]=DOESNT_EXIST; _aStateVarIndices[(int)(ATMOSPHERE   )][0]=1;
  _aStateVarType[2]=ATMOS_PRECIP;  _aStateVarLayer[2]=DOESNT_EXIST; _aStateVarIndices[(int)(ATMOS_PRECIP )][0]=2;
  _aStateVarType[3]=PONDED_WATER;  _aStateVarLayer[3]=DOESNT_EXIST; _aStateVarIndices[(int)(PONDED_WATER )][0]=3;
  _aStateVarType[4]=RUNOFF;        _aStateVarLayer[4]=RUNOFF;       _aStateVarIndices[(int)(RUNOFF       )][0]=4;

  int count=0;
  for (i=5;i<5+_nSoilVars;i++)
  {
    _aStateVarType [i]=SOIL;
    _aStateVarLayer[i]=count;
    _aStateVarIndices[(int)(SOIL)][count]=i;
    count++;
  }
  _lake_sv=0; //by default, rain on lake goes direct to surface storage [0]

  _aGaugeWeights    =NULL; //Initialized in Initialize
  _aGaugeWtTemp     =NULL;
  _aGaugeWtPrecip   =NULL;
  _aCumulativeBal   =NULL;
  _aFlowBal         =NULL;
  _aCumulativeLatBal=NULL;
  _aFlowLatBal      =NULL;
  _CumulInput       =0.0;
  _CumulOutput      =0.0;
  _initWater        =0.0;

  _UTM_zone=-1;

  _nOutputTimes=0;   _aOutputTimes=NULL;
  _currOutputTimeInd=0;
  _pOutputGroup=NULL;

  _aShouldApplyProcess=NULL; //Initialized in Initialize

  _pTransModel=new CTransportModel(this);
  _pGWModel = NULL; //GW MIGRATE -should initialize with empty GW model
  _pDO =NULL;

#ifdef _MODFLOW_USG_
  _pGWModel = new CGroundwaterModel(this);
#endif

  _pEnsemble = NULL;
  _pStateVar = NULL;
}

/////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CModel::~CModel()
{
  if (DESTRUCTOR_DEBUG){cout<<"DELETING MODEL"<<endl;}
  int c,f,g,i,j,k,kk,p;

  CloseOutputStreams();

  for (p=0;p<_nSubBasins;    p++){delete _pSubBasins    [p];} delete [] _pSubBasins;    _pSubBasins=NULL;
  for (k=0;k<_nHydroUnits;   k++){delete _pHydroUnits   [k];} delete [] _pHydroUnits;   _pHydroUnits=NULL;
  for (g=0;g<_nGauges;       g++){delete _pGauges       [g];} delete [] _pGauges;       _pGauges=NULL;
  for (f=0;f<_nForcingGrids; f++){delete _pForcingGrids [f];} delete [] _pForcingGrids; _pForcingGrids=NULL;
  for (j=0;j<_nProcesses;    j++){delete _pProcesses    [j];} delete [] _pProcesses;    _pProcesses=NULL;
  for (c=0;c<_nCustomOutputs;c++){delete _pCustomOutputs[c];} delete [] _pCustomOutputs;_pCustomOutputs=NULL;
  for (i=0;i<_nObservedTS;   i++){delete _pObservedTS   [i];} delete [] _pObservedTS;   _pObservedTS=NULL;
  if (_pModeledTS != NULL){
    for (i = 0; i < _nObservedTS; i++){ delete _pModeledTS[i]; } delete[] _pModeledTS;    _pModeledTS = NULL;
  }
  for (i=0;i<_nObsWeightTS;  i++){delete _pObsWeightTS  [i];} delete [] _pObsWeightTS;  _pObsWeightTS=NULL;
  for (j=0;j<_nDiagnostics;  j++){delete _pDiagnostics  [j];} delete [] _pDiagnostics;  _pDiagnostics=NULL;
  for (j=0;j<_nDiagPeriods;  j++){delete _pDiagPeriods  [j];} delete [] _pDiagPeriods;  _pDiagPeriods=NULL;
  for (j=0;j<_nAggDiagnostics; j++){delete _pAggDiagnostics[j];} delete [] _pAggDiagnostics; _pAggDiagnostics=NULL;

  if (_aCumulativeBal!=NULL){
    for (k=0;k<_nHydroUnits;   k++){delete [] _aCumulativeBal[k];} delete [] _aCumulativeBal; _aCumulativeBal=NULL;
  }
  if (_aFlowBal!=NULL){
    for (k=0;k<_nHydroUnits;   k++){delete [] _aFlowBal[k];      } delete [] _aFlowBal;       _aFlowBal=NULL;
  }
  if (_aCumulativeLatBal!=NULL){delete [] _aCumulativeLatBal; _aCumulativeLatBal=NULL;}
  if (_aFlowLatBal      !=NULL){delete [] _aFlowLatBal;       _aFlowLatBal=NULL;}
  if (_aGaugeWeights!=NULL){
    for (k=0;k<_nHydroUnits;   k++){delete [] _aGaugeWeights[k]; } delete [] _aGaugeWeights;  _aGaugeWeights=NULL;
  }
  if (_aGaugeWtPrecip!=NULL){
    for (k=0;k<_nHydroUnits;   k++){delete [] _aGaugeWtPrecip[k]; } delete [] _aGaugeWtPrecip;  _aGaugeWtPrecip=NULL;
  }
  if (_aGaugeWtTemp!=NULL){
    for (k=0;k<_nHydroUnits;   k++){delete [] _aGaugeWtTemp[k]; } delete [] _aGaugeWtTemp;  _aGaugeWtTemp=NULL;
  }
  if (_aShouldApplyProcess!=NULL){
    for (k=0;k<_nProcesses;   k++){delete [] _aShouldApplyProcess[k]; } delete [] _aShouldApplyProcess;  _aShouldApplyProcess=NULL;
  }
  for (kk=0;kk<_nHRUGroups;kk++)  {delete _pHRUGroups[kk];    } delete [] _pHRUGroups;      _pHRUGroups  =NULL;
  for (kk=0;kk<_nSBGroups;kk++ )  {delete _pSBGroups[kk];     } delete [] _pSBGroups;       _pSBGroups  =NULL;
  for (j=0;j<_nTransParams;j++)   {delete _pTransParams[j];   } delete [] _pTransParams;    _pTransParams=NULL;
  for (j=0;j<_nClassChanges;j++)  {delete _pClassChanges[j];  } delete [] _pClassChanges;   _pClassChanges=NULL;
  for (j=0;j<_nParamOverrides;j++){delete _pParamOverrides[j];} delete [] _pParamOverrides; _pParamOverrides=NULL;

  for (i=0;i<_nPerturbations;   i++)
  {
    delete [] _pPerturbations[i]->eps;
    delete _pPerturbations   [i];
  }
  delete [] _pPerturbations;

  delete [] _aStateVarType;  _aStateVarType=NULL;
  delete [] _aStateVarLayer; _aStateVarLayer=NULL;
  delete [] _aSubBasinOrder; _aSubBasinOrder=NULL;
  delete [] _aOrderedSBind;  _aOrderedSBind=NULL;
  delete [] _aDownstreamInds;_aDownstreamInds=NULL;
  delete [] _aOutputTimes;   _aOutputTimes=NULL;
  delete [] _aObsIndex;      _aObsIndex=NULL;

  delete [] _aDAscale;       _aDAscale=NULL;
  delete [] _aDAscale_last;  _aDAscale_last=NULL;
  delete [] _aDAQadjust;     _aDAQadjust=NULL;
  delete [] _aDADrainSum;    _aDADrainSum=NULL;
  delete [] _aDADownSum;     _aDADownSum=NULL;
  delete [] _aDAlength;      _aDAlength=NULL;
  delete [] _aDAtimesince;   _aDAtimesince=NULL;
  delete [] _aDAoverride;    _aDAoverride=NULL;
  delete [] _aDAobsQ;        _aDAobsQ=NULL;

  this->DestroyAllLanduseClasses();
  this->DestroyAllSoilClasses();
  this->DestroyAllVegClasses();
  this->DestroyAllTerrainClasses();
  this->DestroyAllSoilProfiles();
  this->DestroyAllChannelXSections();

  delete _pTransModel;
  delete _pEnsemble;
  delete _pGWModel;
  delete _pStateVar;
  delete _pDO;

  delete [] _PETBlends_type;
  delete [] _PETBlends_wts;
  delete [] _PotMeltBlends_type;
  delete [] _PotMeltBlends_wts;
}

/*****************************************************************
   Accessors
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Returns pointer to global parameters object
///
/// \return Pointer to global parameters object
//
CGlobalParams* CModel::GetGlobalParams() const {
  return _pGlobalParams;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of sub basins in model
///
/// \return Integer number of sub basins
//
int CModel::GetNumSubBasins   () const{return _nSubBasins;}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of SB groups
///
/// \return Integer number of SB groups
//
int CModel::GetNumSubBasinGroups() const { return _nSBGroups; }

//////////////////////////////////////////////////////////////////
/// \brief Returns number of HRUs in model
///
/// \return Integer number of HRUs
//
int CModel::GetNumHRUs        () const{return _nHydroUnits;}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of HRU groups
///
/// \return Integer number of HRU groups
//
int CModel::GetNumHRUGroups   () const{return _nHRUGroups;}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of gauges in model
///
/// \return Integer number of gauges in model
//
int CModel::GetNumGauges      () const{return _nGauges;}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of gridded forcings in model
///
/// \return Integer number of gridded forcings in model
//
int CModel::GetNumForcingGrids () const{return _nForcingGrids;}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of state variables per HRU in model
///
/// \return Integer number of state variables per HRU in model
//
int CModel::GetNumStateVars   () const{return _nStateVars;}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of soil layers
///
/// \return Integer number of soil layers
//
int CModel::GetNumSoilLayers  () const{return _nSoilVars;}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of hydrologic processes simulated by model
///
/// \return Integer number of hydrologic processes simulated by model
//
int CModel::GetNumProcesses   () const{return _nProcesses;}

//////////////////////////////////////////////////////////////////
/// \brief Returns total modeled watershed area
///
/// \return total modeled watershed area [km2]
//
double CModel::GetWatershedArea () const{return _WatershedArea;}


//////////////////////////////////////////////////////////////////
/// \brief Returns number of observation time series
/// \return number of observation time series
//
int    CModel::GetNumObservedTS() const { return _nObservedTS; }

//////////////////////////////////////////////////////////////////
/// \brief Returns observation time series i
/// \param i [in] index of observation time series
/// \return pointer to observation time series i
//
const CTimeSeriesABC *CModel::GetObservedTS(const int i) const {
  return _pObservedTS[i];
}
//////////////////////////////////////////////////////////////////
/// \brief Returns simulated equivalent of observation time series i
/// \param i [in] index of observation time series
/// \return pointer to simulated equivalent of observation time series i
//
const CTimeSeriesABC* CModel::GetSimulatedTS(const int i) const {
  return _pModeledTS[i];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns specific hydrologic process denoted by parameter
///
/// \param j [in] Process index
/// \return pointer to hydrologic process corresponding to passed index j
//
CHydroProcessABC *CModel::GetProcess(const int j) const
{

#ifdef _STRICTCHECK_
  ExitGracefullyIf((j<0) || (j>=_nProcesses),"CModel GetProcess::improper index",BAD_DATA);
#endif
  return _pProcesses[j];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns specific gauge denoted by index
///
/// \param g [in] Gauge index
/// \return pointer to gauge corresponding to passed index g
//
CGauge *CModel::GetGauge(const int g) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((g<0) || (g>=_nGauges),"CModel GetGauge::improper index",BAD_DATA);
#endif
  return _pGauges[g];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns specific forcing grid denoted by index
///
/// \param f [in] Forcing Grid index
/// \return pointer to gauge corresponding to passed index g
//
CForcingGrid *CModel::GetForcingGrid(const forcing_type &ftyp) const
{
  int f=GetForcingGridIndexFromType(ftyp);
#ifdef _STRICTCHECK_
  if((f<0) || (f>=_nForcingGrids)) { cout<<"Invalid forcing type: "<<ForcingToString(ftyp)<<" ("<<f<<")"<<endl; }
  ExitGracefullyIf((f<0) || (f>=_nForcingGrids),"CModel GetForcingGrid::improper index",RUNTIME_ERR);
  ExitGracefullyIf(_pForcingGrids[f]==NULL,"CModel GetForcingGrid:: NULL forcing grid",RUNTIME_ERR);
#endif
  return _pForcingGrids[f];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns specific HRU denoted by index k
///
/// \param k [in] HRU index
/// \return pointer to HRU corresponding to passed index k
//
CHydroUnit *CModel::GetHydroUnit(const int k) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((k<0) || (k>=_nHydroUnits),"CModel GetHydroUnit::improper index",BAD_DATA);
#endif
  return _pHydroUnits[k];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns specific HRU with HRU identifier HRUID
///
/// \param HRUID [in] HRU identifier
/// \return pointer to HRU corresponding to passed ID HRUID, NULL if no such HRU exists
//
CHydroUnit *CModel::GetHRUByID(const long long int HRUID) const
{
  static int last_k=0;
  //smart find
  int k;
  for (int i=0;i<_nHydroUnits;i++){
    k=NearSearchIndex(i,last_k,_nHydroUnits);
    if (HRUID==_pHydroUnits[k]->GetHRUID()){ last_k=k; return _pHydroUnits[k];}
  }
  return NULL;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns specific HRU group denoted by parameter kk
///
/// \param kk [in] HRU group index
/// \return pointer to HRU group corresponding to passed index kk
//
CHRUGroup  *CModel::GetHRUGroup(const int kk) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((kk<0) || (kk>=_nHRUGroups),"CModel GetHRUGroup::improper index",BAD_DATA);
#endif
  return _pHRUGroups[kk];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns specific HRU group denoted by string parameter
///
/// \param name [in] String name of HRU group
/// \return pointer to HRU group corresponding to passed name, or NULL if this group doesn't exist
//
CHRUGroup  *CModel::GetHRUGroup(const string name) const
{
  for (int kk=0;kk<_nHRUGroups;kk++){
    if (!name.compare(_pHRUGroups[kk]->GetName())){
      return _pHRUGroups[kk];
    }
  }
  return NULL;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns true if HRU with global index k is in specified HRU Group
///
/// \param k [in] HRU global index
/// \param HRUGroupName [in] String name of HRU group
/// \return true if HRU k is in HRU Group specified by HRUGroupName
//
bool CModel::IsInHRUGroup(const int k, const string HRUGroupName) const
{
  CHRUGroup *pGrp=NULL;
  pGrp=GetHRUGroup(HRUGroupName);
  if (pGrp == NULL){ return false; }//throw warning?

  int kk = pGrp->GetGlobalIndex();
  for (int k_loc=0; k_loc<_pHRUGroups[kk]->GetNumHRUs(); k_loc++)
  {
    if (_pHRUGroups[kk]->GetHRU(k_loc)->GetGlobalIndex()==k){return true;}
  }
  return false;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns specific Sub basin denoted by index parameter
///
/// \param p [in] Sub basin index
/// \return pointer to the Sub basin object corresponding to passed index
//
CSubBasin  *CModel::GetSubBasin(const int p) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((p<0) || (p>=_nSubBasins),"CModel GetSubBasin::improper index",BAD_DATA);
#endif
  return _pSubBasins[p];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns index of subbasin downstream from subbasin referred to by index
///
/// \param p [in] List index for accessing subbasin
/// \return downstream subbasin index, if input index is valid; -1 if there is no downstream basin
//
int         CModel::GetDownstreamBasin(const int p) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((p<0) || (p>=_nSubBasins),"GetDownstreamBasin: Invalid index",BAD_DATA);
#endif
  return _aDownstreamInds[p];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns subbasin object corresponding to passed subbasin ID
///
/// \param SBID [in] long long integer sub basin ID
/// \return pointer to Sub basin object corresponding to passed ID, if ID is valid
//
CSubBasin  *CModel::GetSubBasinByID(const long long SBID) const
{
  static int last_p=0;
  int p;
  if (SBID < 0) { return NULL; }
  //smart find
  for (int i=0;i<_nSubBasins;i++){
    p = NearSearchIndex(i, last_p, _nSubBasins);
    if (_pSubBasins[p]->GetID()==SBID){last_p=p; return _pSubBasins[p];}
  }
  return NULL;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns sub basin index corresponding to passed subbasin ID
///
/// \param SBID [in] Integer subbasin ID
/// \return SubBasin index corresponding to passed ID, if ID is valid
//
int         CModel::GetSubBasinIndex(const long long SBID) const
{
  static int last_p = 0;
  int p;
  if (SBID<0){return DOESNT_EXIST;}
  //smart find
  for (int i = 0; i < _nSubBasins; i++) {
    p = NearSearchIndex(i, last_p, _nSubBasins);
    if (_pSubBasins[p]->GetID() == SBID) { last_p = p; return p; }
  }
  return INDEX_NOT_FOUND;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns array of pointers to subbasins upstream of subbasin SBID, including that subbasin
///
/// \param SBID [in] long long int subbasin ID
/// \param nUpstream [out] size of array of pointers of subbasins
/// \return array of pointers to subbasins upstream of subbasin SBID, including that subbasin
//
const CSubBasin **CModel::GetUpstreamSubbasins(const long long SBID,int &nUpstream) const
{
  static const CSubBasin **pSBs=new const CSubBasin *[_nSubBasins];

  bool *isUpstr=new bool [_nSubBasins];
  for(int p=0;p<_nSubBasins;p++) { isUpstr[p]=false; }

  int p=GetSubBasinIndex(SBID);
  if((p==DOESNT_EXIST) || (p==INDEX_NOT_FOUND)) {
    string warn="CModel::GetUpstreamSubbasins: invalid subbasin ID "+to_string(SBID)+" (:ReservoirDownstreamDemand command ? )";
    ExitGracefully(warn.c_str(),BAD_DATA);return NULL;
  }
  isUpstr[p]=true;

  const int MAX_ITER=1000;
  int numUpstr=0;
  int numUpstrOld=1;
  int iter=0;
  int down_p;
  do
  {
    numUpstrOld=numUpstr;
    for(p=0;p<_nSubBasins;p++) {
      down_p=GetSubBasinIndex(_pSubBasins[p]->GetDownstreamID());
      if(down_p!=DOESNT_EXIST) {
        if(isUpstr[down_p]==true) { isUpstr[p]=true;}
      }
    }
    numUpstr=0;
    for(p=0;p<_nSubBasins;p++) {
      if(isUpstr[p]==true) { numUpstr++; }
    }
    iter++;
  } while ((iter<MAX_ITER) && (numUpstr!=numUpstrOld));
  //cout<<"upstream basin calculations iterations = "<<iter<<" "<<numUpstr<<" basins found upstream of basin "<<SBID<<endl;
  nUpstream=numUpstr;
  int count=0;
  for(p=0;p<_nSubBasins;p++) {
    if (isUpstr[p]==true){pSBs[count]=_pSubBasins[p];count++; }
  }
  delete [] isUpstr;
  return pSBs;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns true if subbasin with ID SBID is upstream of (or is) basin with subbasin SBIDdown
/// \notes recursive call, keeps marching downstream until outlet or SBIDdown is encounterd
/// \param SBID [in] ID of subbasin being queried
/// \param SBIDdown [in] subbasin ID basis of query
/// \return true if subbasin with ID SBID is upstream of (or is) basin with subbasin SBIDdown
//
bool  CModel::IsSubBasinUpstream(const long long SBID,const long long SBIDdown) const
{
  if      (SBID==DOESNT_EXIST    ) { return false;} //end of the recursion line
  else if (SBIDdown==SBID)         { return true; } //a subbasin is upstream of itself (even handles loops on bad networks)
  else if (SBIDdown==DOESNT_EXIST) { return true; } //everything is upstream of an outlet
  else if (GetSubBasinByID(SBID)->GetDownstreamID()==SBIDdown){return true;} //directly upstream
  else {
    return IsSubBasinUpstream(GetSubBasinByID(SBID)->GetDownstreamID(),SBIDdown);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns specific subbasin group denoted by parameter pp
///
/// \param pp [in] subbasin group index
/// \return pointer to subbasin group corresponding to passed index pp
//
CSubbasinGroup  *CModel::GetSubBasinGroup(const int pp) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((pp<0) || (pp>=_nSBGroups),"CModel GetSubBasinGroup::improper index",BAD_DATA);
#endif
  return _pSBGroups[pp];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns specific subbasin group denoted by string parameter
///
/// \param name [in] String name of subbasin group
/// \return pointer to subbasin group corresponding to passed name, or NULL if this group doesn't exist
//
CSubbasinGroup  *CModel::GetSubBasinGroup(const string name) const
{
  for(int pp=0;pp<_nSBGroups;pp++) {
    if(!name.compare(_pSBGroups[pp]->GetName())) {
      return _pSBGroups[pp];
    }
  }
  return NULL;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns true if subbasin with subbasin ID SBID is in specified subbasin Group
///
/// \param SBID [in] subbasin identifier
/// \param SBGroupName [in] String name of subbasin group
/// \return true if subbasin SBID is in subbasin Group specified by SBGroupName
//
bool CModel::IsInSubBasinGroup(const long long SBID,const string SBGroupName) const
{
  CSubbasinGroup *pGrp=NULL;
  pGrp=GetSubBasinGroup(SBGroupName);
  if(pGrp == NULL) { return false; }//throw warning?

  int pp = pGrp->GetGlobalIndex();
  for(int p_loc=0; p_loc<_pSBGroups[pp]->GetNumSubbasins(); p_loc++)
  {
    if(_pSBGroups[pp]->GetSubBasin(p_loc)->GetID()==SBID) { return true; }
  }
  return false;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns hydrologic process type corresponding to passed index
///
/// \param j [in] Integer index
/// \return Process type corresponding to process with passed index
//
process_type CModel::GetProcessType(const int j) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((j<0) || (j>=_nProcesses),"CModel GetProcessType::improper index",BAD_DATA);
#endif
  return _pProcesses[j]->GetProcessType();
}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of connections of hydrological process associated with passed index
///
/// \param j [in] Integer index corresponding to a hydrological process
/// \return Number of connections associated with hydrological process symbolized by index
/// \note should be called only by solver
//
int         CModel::GetNumConnections (const int j) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((j<0) || (j>=_nProcesses),"CModel GetNumConnections::improper index",BAD_DATA);
#endif
  return _pProcesses[j]->GetNumConnections();
}
//////////////////////////////////////////////////////////////////
/// \brief Returns number of forcing perturbations in model
//
int         CModel::GetNumForcingPerturbations() const {return _nPerturbations;}

//////////////////////////////////////////////////////////////////
/// \brief Returns state variable type corresponding to passed state variable array index
///
/// \param i [in] state variable array index (>=0, <nStateVariables)
/// \return State variable type corresponding to passed state variable array index
//
sv_type     CModel::GetStateVarType(const int i) const
{
#ifdef _STRICTCHECK_
  string warn="CModel GetStateVarType::improper index ("+to_string(i)+")";
  ExitGracefullyIf((i<0) || (i>=_nStateVars),warn.c_str(),BAD_DATA);
#endif
  return _aStateVarType[i];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns index of state variable type passed
///
/// \param type [in] State variable type
/// \return Index which corresponds to the state variable type passed, if it exists; DOESNT_EXIST (-1) otherwise
/// \note  should only be used for state variable types without multiple levels; issues for soils, e.g.
//
int         CModel::GetStateVarIndex(sv_type type) const
{
  return _aStateVarIndices[(int)(type)][0];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns index of state variable type passed (for repeated state variables)
///
/// \param type [in] State variable type
/// \param layer [in] Integer identifier of the layer of interest (or, possibly DOESNT_EXIST if variable doesn't have layers)
/// \return Index which corresponds to the state variable type passed, if it exists; DOESNT_EXIST (-1) otherwise
//
int         CModel::GetStateVarIndex(sv_type type, int layer) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((layer!=DOESNT_EXIST) && ((layer<0) || (layer>=MAX_SV_LAYERS)),
                   "CModel GetStateVarIndex::improper layer",BAD_DATA);
#endif
  if (layer==DOESNT_EXIST){return _aStateVarIndices[(int)(type)][0];    }
  else                    {return _aStateVarIndices[(int)(type)][layer];}
}

//////////////////////////////////////////////////////////////////
/// \brief Uses state variable type index to access index of layer to which it corresponds
///
/// \param ii [in] Index referencing a type of state variable
/// \return Index which corresponds to soil layer, or 0 for a non-layered state variable
//
int         CModel::GetStateVarLayer(const int ii) const
{
  int count=0;
  for (int i=0;i<ii;i++){
    if (_aStateVarType[i]==_aStateVarType[ii]){count++;}
  }
  return count;
}

//////////////////////////////////////////////////////////////////
/// \brief Checks if state variable passed exists in model
///
/// \param typ [in] State variable type
/// \return Boolean indicating whether state variable exists in model
//
bool        CModel::StateVarExists(sv_type typ) const
{
  return (GetStateVarIndex(typ)!=DOESNT_EXIST);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns lake storage variable index
///
/// \return Integer index of lake storage variable
//
int         CModel::GetLakeStorageIndex() const{return _lake_sv;}

//////////////////////////////////////////////////////////////////
/// \brief Returns gauge index of gauge with specified name
///
/// \return Integer index of gauge
/// \param name [in] specified name
//
int  CModel::GetGaugeIndexFromName (const string name) const
{
  for (int g=0;g<_nGauges;g++){
    if (name==_pGauges[g]->GetName()){return g;}
  }
  return DOESNT_EXIST;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns forcing grid index of forcing grid with specified type
///
/// \return Integer index of forcing grid
/// \param name [in] specified type
//
int  CModel::GetForcingGridIndexFromType (const forcing_type &typ) const
{
  for (int f=0;f<_nForcingGrids;f++){
    if (typ==_pForcingGrids[f]->GetForcingType()){return f;}
  }
  return DOESNT_EXIST;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns current mass/energy flux (mm/d, MJ/m2/d, mg/m2/d) between two storage compartments iFrom and iTo
/// \details required for advective transport processes
///
/// \param k [in] HRU index
/// \param js [in] index of process connection (i.e., j*)
/// \param &Options [in] Global model options information
//
double CModel::GetFlux(const int k, const int js, const optStruct &Options) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((k<0) || (k>=_nHydroUnits),"CModel::GetFlux: bad HRU index",RUNTIME_ERR);
  ExitGracefullyIf((js<0) || (js>=_nTotalConnections),"CModel::GetFlux: bad connection index",RUNTIME_ERR);
#endif
  return _aFlowBal[k][js]/Options.timestep;
}
//////////////////////////////////////////////////////////////////
/// \brief returns concentration or temperature within hru k with storage index i
///
/// \param k [in] HRU index
/// \param i [in] mass/energy state variable index
//
double CModel::GetConcentration(const int k, const int i) const
{
  return _pTransModel->GetConcentration(k,i);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns current mass/energy flow (mm-m2/d, MJ/d, mg/d) between two storage compartments iFrom and iTo
/// \details required for advective transport processes
///
/// \param k [in] HRU index
/// \param qs [in] global index of process connection (i.e., q*)
/// \param &Options [in] Global model options information
//
double CModel::GetLatFlow(const int qs,const optStruct &Options) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((qs<0) || (qs>=_nTotalConnections),"CModel::GetFlux: bad connection index",RUNTIME_ERR);
#endif
  return _aFlowLatBal[qs]/Options.timestep;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns cumulative flux to/from storage unit i
///
/// \param k [in] index of HRU
/// \param i [in] index of storage compartment
/// \param to [in] true if evaluating cumulative flux to storage compartment, false for 'from'
/// \return cumulative flux to storage compartment i in hru K
//
double CModel::GetCumulativeFlux(const int k, const int i, const bool to) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((k<0) || (k>=_nHydroUnits),"CModel::GetCumulativeFlux: bad HRU index",RUNTIME_ERR);
  ExitGracefullyIf((i<0) || (i>=_nStateVars),"CModel::GetCumulativeFlux: bad state var index",RUNTIME_ERR);
#endif
  int js=0;
  double sum=0;
  double area=_pHydroUnits[k]->GetArea();
  int jss=0;
  for(int j = 0; j < _nProcesses; j++)
  {
    for(int q = 0; q < _pProcesses[j]->GetNumConnections(); q++)//each process may have multiple connections
    {
      if(( to) && (_pProcesses[j]->GetToIndices()[q]   == i)){ sum+=_aCumulativeBal[k][js]; }
      if((!to) && (_pProcesses[j]->GetFromIndices()[q] == i)){ sum+=_aCumulativeBal[k][js]; }
      js++;
    }

    if(_pProcesses[j]->GetNumLatConnections()>0)
    {
      CLateralExchangeProcessABC *pProc=(CLateralExchangeProcessABC*)_pProcesses[j];
      for(int q = 0; q < _pProcesses[j]->GetNumLatConnections(); q++)//each process may have multiple connections
      {
        if(( to) && (pProc->GetToHRUIndices()[q]  ==k) && (pProc->GetLateralToIndices()[q]  ==i)){sum+=_aCumulativeLatBal[jss]/area; }
        if((!to) && (pProc->GetFromHRUIndices()[q]==k) && (pProc->GetLateralFromIndices()[q]==i)){sum+=_aCumulativeLatBal[jss]/area; }
        jss++;
      }
    }
  }

  return sum;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns cumulative gross flux between unit iFrom and iTo in HRU k
///
/// \param k [in] index of HRU
/// \param iFrom [in] index of storage compartment
/// \param iTo [in] index of storage compartment
/// \return cumulative gross flux between unit iFrom and iTo in HRU k
//  does not address fluxes due to lateral flux
double CModel::GetCumulFluxBetween(const int k,const int iFrom,const int iTo) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((k<0) || (k>=_nHydroUnits),"CModel::GetCumulativeFlux: bad HRU index",RUNTIME_ERR);
  ExitGracefullyIf((iFrom<0) || (iTo>=_nStateVars),"CModel::GetCumulativeFlux: bad state var index",RUNTIME_ERR);
#endif
  int q,js=0;
  double sum=0;
  const int *iFromp;
  const int *iTop;
  int nConn;
  for(int j = 0; j < _nProcesses; j++)
  {
    iFromp=_pProcesses[j]->GetFromIndices();
    iTop  =_pProcesses[j]->GetToIndices();
    nConn =_pProcesses[j]->GetNumConnections();
    for (q = 0; q < nConn; q++)//each process may have multiple connections
    {
      if( (iTop  [q]== iTo) && (iFromp[q]== iFrom)){ sum+=_aCumulativeBal[k][js]; }
      if( (iFromp[q]== iTo) && (iTop  [q]== iFrom)){ sum-=_aCumulativeBal[k][js]; }
      js++;
    }
  }

  return sum;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns area-weighted average of specified cumulative flux over watershed
///
/// \param i [in] index of storage compartment
/// \param to [in] true if evaluating cumulative flux to storage compartment, false for 'from'
/// \return Area-weighted average of cumulative flux to storage compartment i
//
double CModel::GetAvgCumulFlux(const int i,const bool to) const
{
  //Area-weighted average
  double sum=0.0;
  for(int k=0;k<_nHydroUnits;k++)
  {
    if(_pHydroUnits[k]->IsEnabled())
    {
      sum +=GetCumulativeFlux(k,i,to)*_pHydroUnits[k]->GetArea();
    }
  }
  return sum/_WatershedArea;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns area-weighted average of  cumulative flux between two compartments over watershed
///
/// \param iFrom [in] index of 'from' storage compartment
/// \param iTo [in] index of 'to' storage compartment
/// \return Area-weighted average of cumulative flux between two compartments over watershed
//
double CModel::GetAvgCumulFluxBet(const int iFrom,const int iTo) const
{
  //Area-weighted average
  double sum=0.0;
  for(int k=0;k<_nHydroUnits;k++)
  {
    if(_pHydroUnits[k]->IsEnabled())
    {
      sum +=GetCumulFluxBetween(k,iFrom,iTo)*_pHydroUnits[k]->GetArea();
    }
  }
  return sum/_WatershedArea;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns options structure model
/// \return pointer to transport model
//
const optStruct  *CModel::GetOptStruct() const{ return _pOptStruct; }

//////////////////////////////////////////////////////////////////
/// \brief Returns transport model
/// \return pointer to transport model
//
CTransportModel  *CModel::GetTransportModel() const{return _pTransModel;}

//////////////////////////////////////////////////////////////////
/// \brief Returns groundwater model
/// \return pointer to groundwater model
//
CGroundwaterModel*CModel::GetGroundwaterModel() const{return _pGWModel;}

//////////////////////////////////////////////////////////////////
/// \brief Returns ensemble setup
/// \return pointer to ensemble setup
//
CEnsemble  *CModel::GetEnsemble() const { return _pEnsemble; }

//////////////////////////////////////////////////////////////////
/// \brief Returns demand optimizer
/// \return pointer to demand optimizer
//
CDemandOptimizer  *CModel::GetManagementOptimizer() const { return _pDO; }

/*****************************************************************
   Watershed Diagnostic Functions
    -aggregate data from subbasins and HRUs
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Returns area-weighted average total precipitation+irrigation rate at all HRUs [mm/d]
///
/// \return Area-weighted average of total precipitation rate [mm/d] over all HRUs
//
double CModel::GetAveragePrecip() const
{
  double sum(0);
  for (int k=0;k<_nHydroUnits;k++)
  {
    if(_pHydroUnits[k]->IsEnabled())
    {
      sum+=(_pHydroUnits[k]->GetForcingFunctions()->precip+
            _pHydroUnits[k]->GetForcingFunctions()->irrigation)*_pHydroUnits[k]->GetArea();
    }
  }
  return sum/_WatershedArea;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns current area-weighted average snowfall rate [mm/d] at all HRUs
///
/// \return Current area-weighted average snowfall rate [mm/d] at all HRUs
//
double CModel::GetAverageSnowfall() const
{
  double sum(0);
  const force_struct *f;
  for (int k=0;k<_nHydroUnits;k++)
  {
    if(_pHydroUnits[k]->IsEnabled())
    {
      f=_pHydroUnits[k]->GetForcingFunctions();
      sum+=(f->precip*f->snow_frac)*_pHydroUnits[k]->GetArea();
    }
  }
  return sum/_WatershedArea;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns current area-weighted average forcing functions at all HRUs
///
/// \return Current area-weighted average forcing functions at all HRUs
//
force_struct CModel::GetAverageForcings() const
{
  static force_struct Fave;
  const force_struct *pF_hru;
  double area_wt;
  ZeroOutForcings(Fave);

  for (int k=0;k<_nHydroUnits;k++)
  {
    if(_pHydroUnits[k]->IsEnabled())
    {
      pF_hru =_pHydroUnits[k]->GetForcingFunctions();
      area_wt=_pHydroUnits[k]->GetArea()/_WatershedArea;

      Fave.precip          +=area_wt*pF_hru->precip;
      Fave.precip_daily_ave+=area_wt*pF_hru->precip_daily_ave;
      Fave.precip_5day     +=area_wt*pF_hru->precip_5day;
      Fave.snow_frac       +=area_wt*pF_hru->snow_frac;

      Fave.temp_ave       +=area_wt*pF_hru->temp_ave;
      Fave.temp_daily_min +=area_wt*pF_hru->temp_daily_min;
      Fave.temp_daily_max +=area_wt*pF_hru->temp_daily_max;
      Fave.temp_daily_ave +=area_wt*pF_hru->temp_daily_ave;
      Fave.temp_month_max +=area_wt*pF_hru->temp_month_max;
      Fave.temp_month_min +=area_wt*pF_hru->temp_month_min;
      Fave.temp_month_ave +=area_wt*pF_hru->temp_month_ave;

      Fave.temp_ave_unc   +=area_wt*pF_hru->temp_ave_unc;
      Fave.temp_min_unc   +=area_wt*pF_hru->temp_min_unc;
      Fave.temp_max_unc   +=area_wt*pF_hru->temp_max_unc;

      Fave.air_dens       +=area_wt*pF_hru->air_dens;
      Fave.air_pres       +=area_wt*pF_hru->air_pres;
      Fave.rel_humidity   +=area_wt*pF_hru->rel_humidity;

      Fave.cloud_cover    +=area_wt*pF_hru->cloud_cover;
      Fave.ET_radia       +=area_wt*pF_hru->ET_radia;
      Fave.ET_radia_flat  +=area_wt*pF_hru->ET_radia_flat;
      Fave.SW_radia       +=area_wt*pF_hru->SW_radia;
      Fave.SW_radia_unc   +=area_wt*pF_hru->SW_radia_unc;
      Fave.SW_radia_net   +=area_wt*pF_hru->SW_radia_net;
      Fave.SW_radia_subcan+=area_wt*pF_hru->SW_radia_subcan;
      Fave.SW_subcan_net  +=area_wt*pF_hru->SW_subcan_net;
      Fave.LW_incoming    +=area_wt*pF_hru->LW_incoming;
      Fave.LW_radia_net   +=area_wt*pF_hru->LW_radia_net;
      Fave.day_length     +=area_wt*pF_hru->day_length;
      Fave.day_angle      +=area_wt*pF_hru->day_angle;   //not really necc.

      Fave.wind_vel       +=area_wt*pF_hru->wind_vel;

      Fave.PET            +=area_wt*pF_hru->PET;
      Fave.OW_PET         +=area_wt*pF_hru->OW_PET;
      Fave.PET_month_ave  +=area_wt*pF_hru->PET_month_ave;

      Fave.potential_melt +=area_wt*pF_hru->potential_melt;

      Fave.recharge       +=area_wt*pF_hru->recharge;
      Fave.precip_temp    +=area_wt*pF_hru->precip_temp;
      Fave.precip_conc    +=area_wt*pF_hru->precip_conc;

      Fave.subdaily_corr  +=area_wt*pF_hru->subdaily_corr;
    }
  }
  return Fave;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns area weighted average of state variables across all modeled HRUs
///
/// \param i [in] State variable index
/// \return Area weighted average of state variables across all modeled HRUs
//
double CModel::GetAvgStateVar(const int i) const
{
  //Area-weighted average
#ifdef _STRICTCHECK_
  ExitGracefullyIf((i<0) || (i>=_nStateVars),"CModel GetAvgStateVar::improper index",BAD_DATA);
#endif
  double sum(0.0);
  for (int k=0;k<_nHydroUnits;k++)
  {
    if(_pHydroUnits[k]->IsEnabled())
    {
      sum+=(_pHydroUnits[k]->GetStateVarValue(i)*_pHydroUnits[k]->GetArea());
    }
  }
  return sum/_WatershedArea;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns area weighted average of concentrations across all modeled HRUs
///
/// \param i [in] State variable index
/// \return Area weighted average of concentrations across all modeled HRUs
//
double CModel::GetAvgConcentration(const int i) const
{
  //Area-weighted average
#ifdef _STRICTCHECK_
  ExitGracefullyIf((i<0) || (i>=_nStateVars),"CModel GetAvgStateVar::improper index",BAD_DATA);
#endif
  double sum(0.0);
  for(int k=0;k<_nHydroUnits;k++)
  {
    if(_pHydroUnits[k]->IsEnabled())
    {
      sum+=(_pTransModel->GetConcentration(k,i)*_pHydroUnits[k]->GetArea());
    }
  }
  return sum/_WatershedArea;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns area-weighted average of specified forcing function over watershed
///
/// \param &ftype [in] enum identifier of forcing function to assess
/// \return Area-weighted average of specified forcing function
//
double CModel::GetAvgForcing (const forcing_type &ftype) const
{
  //Area-weighted average
  double sum=0.0;
  for (int k=0;k<_nHydroUnits;k++)
  {
    if (_pHydroUnits[k]->IsEnabled())
    {
      sum    +=_pHydroUnits[k]->GetForcing(ftype)*_pHydroUnits[k]->GetArea();
    }
  }
  return sum/_WatershedArea;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns total channel storage [mm]
/// \return Total channel storage in all of watershed [mm]
//
double CModel::GetTotalChannelStorage() const
{
  double sum(0);

  for (int p=0;p<_nSubBasins;p++)
  {
    sum+=_pSubBasins[p]->GetChannelStorage();            //[m3]
  }
  return sum/(_WatershedArea*M2_PER_KM2)*MM_PER_METER;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns total reservoir storage [mm]
/// \return Total reservoir storage in all of watershed [mm]
//
double CModel::GetTotalReservoirStorage() const
{
  double sum(0);

  for (int p=0;p<_nSubBasins;p++)
  {
    sum+=_pSubBasins[p]->GetReservoirStorage();          //[m3]
  }
  return sum/(_WatershedArea*M2_PER_KM2)*MM_PER_METER;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns total rivulet storage distributed over watershed [mm]
/// \return Total rivulet storage distributed over watershed [mm]
//
double CModel::GetTotalRivuletStorage() const
{
  double sum(0);

  for (int p=0;p<_nSubBasins;p++)
  {
    sum+=_pSubBasins[p]->GetRivuletStorage();//[m3]
  }

  return sum/(_WatershedArea*M2_PER_KM2)*MM_PER_METER;
}


/*****************************************************************
   Model Creation Functions
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Adds additional HRU to model
///
/// \param *pHRU [in] (valid) pointer to HRU to be added
//
void CModel::AddHRU(CHydroUnit *pHRU)
{
  if (!DynArrayAppend((void**&)(_pHydroUnits),(void*)(pHRU),_nHydroUnits)){
    ExitGracefully("CModel::AddHRU: adding NULL HRU",BAD_DATA);}
}

//////////////////////////////////////////////////////////////////
/// \brief Adds HRU group
///
/// \param *pHRUGroup [in] (valid) pointer to HRU group to be added
//
void CModel::AddHRUGroup(CHRUGroup *pHRUGroup)
{
  for(int kk=0;kk<_nHRUGroups;kk++){
    if(pHRUGroup->GetName()==_pHRUGroups[kk]->GetName()){
      WriteWarning("CModel::AddHRUGroups: cannot add two HRU groups with the same name. Group "+pHRUGroup->GetName()+ " is duplicated in input.",true);
      break;
    }
  }
  if (!DynArrayAppend((void**&)(_pHRUGroups),(void*)(pHRUGroup),_nHRUGroups)){
    ExitGracefully("CModel::AddHRUGroup: adding NULL HRU Group",BAD_DATA);}
}

//////////////////////////////////////////////////////////////////
/// \brief Adds sub basin
///
/// \param *pSB [in] (valid) pointer to Sub basin to be added
//
void CModel::AddSubBasin(CSubBasin *pSB)
{
  if (!DynArrayAppend((void**&)(_pSubBasins),(void*)(pSB),_nSubBasins)){
    ExitGracefully("CModel::AddSubBasin: adding NULL HRU",BAD_DATA);}
}

//////////////////////////////////////////////////////////////////
/// \brief Adds SubBasin group
///
/// \param *pSBGroup [in] (valid) pointer to SubBasin group to be added
//
void CModel::AddSubBasinGroup(CSubbasinGroup *pSBGroup)
{
  for(int pp=0;pp<_nSBGroups;pp++) {
    if(pSBGroup->GetName()==_pSBGroups[pp]->GetName()) {
      WriteWarning("CModel::AddSubBasinGroup: cannot add two Subbasin groups with the same name. Group "+pSBGroup->GetName()+ " is duplicated in input.",true);
      break;
    }
  }
  if(!DynArrayAppend((void**&)(_pSBGroups),(void*)(pSBGroup),_nSBGroups)) {
    ExitGracefully("CModel::AddSubBasinGroup: adding NULL SubBasin Group",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Adds gauge to model
///
/// \param *pGage [in] (valid) pointer to Gauge to be added to model
//
void CModel::AddGauge  (CGauge *pGage)
{
  if (!DynArrayAppend((void**&)(_pGauges),(void*)(pGage),_nGauges)){
    ExitGracefully("CModel::AddGauge: adding NULL Gauge",BAD_DATA);}
}

//////////////////////////////////////////////////////////////////
/// \brief Adds gridded forcing to model
///
/// \param *pGrid [in] (valid) pointer to ForcingGrid to be added to model
//
void CModel::AddForcingGrid  (CForcingGrid *pGrid, forcing_type typ)
{
  if(GetForcingGridIndexFromType(typ) == DOESNT_EXIST) {
    if(!DynArrayAppend((void**&)(_pForcingGrids),(void*)(pGrid),_nForcingGrids)){
      ExitGracefully("CModel::AddForcingGrid: adding NULL ForcingGrid",BAD_DATA);
    }
  }
  else { //overwrite grid
    int f= GetForcingGridIndexFromType(typ);
    if((_pForcingGrids[f])!=(pGrid)) {
      delete _pForcingGrids[f]; //JRC: the big offender - usually we are overwriting with the same forcing grid instance
    }
    _pForcingGrids[f]=pGrid;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Adds transient parameter to model
///
/// \param *pTP [in] (valid) pointer to transient parameter to be added to model
//
void CModel::AddTransientParameter(CTransientParam   *pTP)
{
  if (!DynArrayAppend((void**&)(_pTransParams),(void*)(pTP),_nTransParams)){
    ExitGracefully("CModel::AddTransientParameter: adding NULL transient parameter",BAD_DATA);}
}

//////////////////////////////////////////////////////////////////
/// \brief Adds parameter override to model
///
/// \param *pTP [in] (valid) pointer to transient parameter to be added to model
//
void CModel::AddParameterOverride(param_override   *pPO)
{
  if (!DynArrayAppend((void**&)(_pParamOverrides),(void*)(pPO),_nParamOverrides)){
    ExitGracefully("CModel::AddParameterOverride: adding NULL override",BAD_DATA);}
}
//////////////////////////////////////////////////////////////////
/// \brief Adds class change to model
///
/// \param *pTP [in] (valid) pointer to transient parameter to be added to model
//
void CModel::AddPropertyClassChange(const string      HRUgroup,
                                    const class_type  tclass,
                                    const string      new_class,
                                    const time_struct &tt,
                                    const optStruct   &Options)
{
  class_change *pCC=NULL;
  pCC=new class_change;
  pCC->HRU_groupID=DOESNT_EXIST;
  for (int kk = 0; kk < _nHRUGroups; kk++){
    if (!strcmp(_pHRUGroups[kk]->GetName().c_str(), HRUgroup.c_str()))
    {
      pCC->HRU_groupID=kk;
    }
  }
  if (pCC->HRU_groupID == DOESNT_EXIST){
    string warning = "CModel::AddPropertyClassChange: invalid HRU Group name: " + HRUgroup+ ". HRU group names should be defined in .rvi file using :DefineHRUGroups command. ";
    ExitGracefullyIf(pCC->HRU_groupID == DOESNT_EXIST,warning.c_str(),BAD_DATA_WARN); return;
  }
  pCC->newclass=new_class;
  if ((tclass == CLASS_LANDUSE) && (StringToLUClass(new_class) == NULL)){
    ExitGracefully("CModel::AddPropertyClassChange: invalid land use class specified",BAD_DATA_WARN);return;
  }
  if ((tclass == CLASS_VEGETATION) && (StringToVegClass(new_class) == NULL)){
    ExitGracefully("CModel::AddPropertyClassChange: invalid vegetation class specified",BAD_DATA_WARN);return;
  }
  if ((tclass == CLASS_HRUTYPE) && (StringToHRUType(new_class) == HRU_INVALID_TYPE)){
    ExitGracefully("CModel::AddPropertyClassChange: invalid HRU type specified",BAD_DATA_WARN);return;
  }

  pCC->tclass=tclass;
  if ((tclass != CLASS_VEGETATION) && (tclass != CLASS_LANDUSE) && (tclass!=CLASS_HRUTYPE)){
    ExitGracefully("CModel::AddPropertyClassChange: only vegetation, land use, and HRU type classes may be changed during the course of simulation",BAD_DATA_WARN);return;
  }

  //convert time to model time
  pCC->modeltime= TimeDifference(Options.julian_start_day,Options.julian_start_year,tt.julian_day, tt.year, Options.calendar);

  if (pCC->modeltime>Options.duration){
    string warn;
    warn="Property Class change dated "+tt.date_string+" occurs after model simulation is done; it will not effect results.";
    WriteWarning(warn,Options.noisy);
  }
  if(pCC->modeltime<0) {
    string warn;
    warn="Property Class change dated "+tt.date_string+" occurs before model simulation. All such changes will be processed before time zero, and these changes MUST be input in chronological order.";
    WriteAdvisory(warn,Options.noisy);
  }

  //cout << "PROPERTY CLASS CHANGE " << pCC->HRU_groupID << " " << pCC->tclass << " "<<pCC->modeltime<<endl;
  if (!DynArrayAppend((void**&)(_pClassChanges),(void*)(pCC),_nClassChanges)){
    ExitGracefully("CModel::AddPropertyClassChange: adding NULL property class change",BAD_DATA);}
}
//////////////////////////////////////////////////////////////////
/// \brief Adds observed time series to model
///
/// \param *pTS [in] (valid) pointer to observed time series to be added to model
//
void CModel::AddObservedTimeSeries(CTimeSeriesABC *pTS)
{
  if (!DynArrayAppend((void**&)(_pObservedTS),(void*)(pTS),_nObservedTS)){
    ExitGracefully("CModel::AddObservedTimeSeries: adding NULL observation time series",BAD_DATA);}
}
//////////////////////////////////////////////////////////////////
/// \brief Adds observation weighting time series to model
///
/// \param *pTS [in] (valid) pointer to observation weights time series to be added to model
//
void CModel::AddObservedWeightsTS(CTimeSeriesABC   *pTS)
{
  if (!DynArrayAppend((void**&)(_pObsWeightTS),(void*)(pTS),_nObsWeightTS)){
    ExitGracefully("CModel::AddObservedWeightsTS: adding NULL observed weights time series",BAD_DATA);}
}
//////////////////////////////////////////////////////////////////
/// \brief Adds diagnostic to model
///
/// \param *pDiag [in] (valid) pointer to diagnostic to be added to model
//
void CModel::AddDiagnostic(CDiagnostic       *pDiag)
{
  if (!DynArrayAppend((void**&)(_pDiagnostics),(void*)(pDiag),_nDiagnostics)){
    ExitGracefully("CModel::AddDiagnostic: adding NULL diagnostic",BAD_DATA);}
}
//////////////////////////////////////////////////////////////////
/// \brief Adds aggregate diagnostic to model
///
/// \param *pDiag [in] (valid) pointer to diagnostic to be added to model
//
void CModel::AddAggregateDiagnostic(agg_stat stat, string datatype, int group_ind)
{
  agg_diag *agg=new agg_diag();
  agg->aggtype=stat;
  agg->datatype=datatype;
  agg->kk=group_ind;

  if (!DynArrayAppend((void**&)(_pAggDiagnostics),(void*)(agg),_nAggDiagnostics)){
    ExitGracefully("CModel::AddAggregateDiagnostic: adding NULL diagnostic",BAD_DATA);}
}
//////////////////////////////////////////////////////////////////
/// \brief Adds diagnostic period to model
///
/// \param *pDiag [in] (valid) pointer to diagnostic period to be added to model
//
void CModel::AddDiagnosticPeriod(CDiagPeriod       *pDiagPer)
{
  if(!DynArrayAppend((void**&)(_pDiagPeriods),(void*)(pDiagPer),_nDiagPeriods)) {
    ExitGracefully("CModel::AddDiagnosticPeriod: adding NULL diagnostic period",BAD_DATA);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Adds model output time to model
///
/// \param timestamp [in] ISO-format timestamp string YYYY-MM-DD hh:mm:ss
/// \note must be called after model starttime & duration are set
//
void CModel::AddModelOutputTime   (const time_struct &tt_out, const optStruct &Options)
{

  //local time
  double t_loc = TimeDifference(Options.julian_start_day,Options.julian_start_year,tt_out.julian_day,tt_out.year,Options.calendar);

  ExitGracefullyIf(t_loc<0,
                   "AddModelOutputTime: Cannot have model output time prior to start of simulation",BAD_DATA_WARN);
  if (t_loc>Options.duration){
    WriteWarning("AddModelOutputTime: model output time specified after end of simulation. It will be ignored",Options.noisy);
    return;
  }

  //add to end of array
  double *tmpOut=new double [_nOutputTimes+1];
  for (int i=0;i<_nOutputTimes;i++){
    tmpOut[i]=_aOutputTimes[i];
  }
  tmpOut[_nOutputTimes]=t_loc;
  delete [] _aOutputTimes;
  _aOutputTimes=tmpOut;
  _nOutputTimes++;
}

//////////////////////////////////////////////////////////////////
/// \brief Adds a hydrological process to system
/// \details Adds a hydrological process that moves water or energy from state
/// variable (e.g., an array of storage units) i[] to state variables j[]
///
/// \param *pHydroProc [in] (valid) pointer to Hydrological process to be added
//
void CModel::AddProcess(CHydroProcessABC *pHydroProc)
{
  ExitGracefullyIf(pHydroProc==NULL,"CModel AddProcess::NULL process"  ,BAD_DATA);
  for (int q=0;q<pHydroProc->GetNumConnections();q++)
  {
    int i=pHydroProc->GetFromIndices()[q];
    int j=pHydroProc->GetToIndices  ()[q];
    ExitGracefullyIf((i<0) && (i>=_nStateVars),"CModel AddProcess::improper storage index",BAD_DATA);
    ExitGracefullyIf((j<0) && (j>=_nStateVars),"CModel AddProcess::improper storage index",BAD_DATA);
  }
  if (!DynArrayAppend((void**&)(_pProcesses),(void*)(pHydroProc),_nProcesses)){
    ExitGracefully("CModel::AddProcess: adding NULL hydrological process",BAD_DATA);}
}

//////////////////////////////////////////////////////////////////
/// \brief Adds state variables during model construction
/// \details Called during model construction (only from parse of RVI file) to dynamically
/// generate a list of state variables in the model
///
/// \note Soil structure (or any SV with multiple layers) needs to be created in Model Constructor beforehand (is this true for other multilayer SVs??)
///
/// \param *aSV [in] array of state variables types to be added
/// \param *aLev [in] array of state variable levels to be added
/// \param nSV [in] Integer number of SVs to be added (size of aSV[] and aLev[])
//
void  CModel::AddStateVariables( const sv_type *aSV,
                                 const int     *aLev,
                                 const int      nSV)
{
  int i,ii;
  bool found;
  ExitGracefullyIf(nSV > MAX_STATE_VARS,
    "CModel::AddStateVariables: exceeded maximum number of state variables manipulable by one process", RUNTIME_ERR);
  for (ii=0;ii<nSV;ii++)
  {
    found=false;
    for (i=0;i<_nStateVars;i++){
      if ((_aStateVarType [i]==aSV [ii]) &&
          (_aStateVarLayer[i]==aLev[ii])){found=true;}
    }
    //found=(_aStateVarIndices[(int)(aSV[ii])][aLev[ii]]!=DOESNT_EXIST);
    if (!found)
    {
      sv_type *tmpSV=new sv_type[_nStateVars+1];
      int     *tmpLy=NULL;
      tmpLy=new int    [_nStateVars+1];
      ExitGracefullyIf(tmpLy==NULL,"CModel::AddStateVariables",OUT_OF_MEMORY);
      for (i=0;i<_nStateVars;i++)//copy old arrays
      {
        tmpSV[i]=_aStateVarType  [i];
        tmpLy[i]=_aStateVarLayer [i];
      }
      //set values of new array items
      tmpSV[_nStateVars]=aSV[ii];
      tmpLy[_nStateVars]=aLev[ii];
      //add index to state variable lookup table
      ExitGracefullyIf((int)(aSV[ii])> MAX_STATE_VAR_TYPES,
                       "CModel::AddStateVariables: bad type specified",RUNTIME_ERR);
      ExitGracefullyIf((aLev[ii]<-1) || (aLev[ii]>=MAX_SV_LAYERS),
                       "CModel::AddStateVariables: bad layer index specified",RUNTIME_ERR);
      if (aLev[ii]==DOESNT_EXIST){_aStateVarIndices[(int)(aSV[ii])][0       ]=_nStateVars;}
      else                       {_aStateVarIndices[(int)(aSV[ii])][aLev[ii]]=_nStateVars;}
      //delete & replace old arrays
      delete [] _aStateVarType;  _aStateVarType=NULL;
      delete [] _aStateVarLayer; _aStateVarLayer=NULL;
      _aStateVarType =tmpSV;
      _aStateVarLayer=tmpLy;
      _nStateVars++;
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Adds custom output object
///
/// \param *pCO [in] (valid) pointer to Custom output object to be added
//
void CModel::AddCustomOutput(CCustomOutput *pCO)
{
  if (!DynArrayAppend((void**&)(_pCustomOutputs),(void*)(pCO),_nCustomOutputs)){
    ExitGracefully("CModel::AddCustomOutput: adding NULL custom output",BAD_DATA);}
}

//////////////////////////////////////////////////////////////////
/// \brief Adds custom table object
///
/// \param *pTab [in] (valid) pointer to Custom table object to be added
//
void CModel::AddCustomTable(CCustomTable *pTab)
{
  if (!DynArrayAppend((void**&)(_pCustomTables),(void*)(pTab),_nCustomTables)){
    ExitGracefully("CModel::AddCustomTable: adding NULL custom table",BAD_DATA);}
}
//////////////////////////////////////////////////////////////////
/// \brief adds additional forcing perturbation
/// \param type [in] forcing type to be perturbed
/// \param distrib [in] sampling distribution type
/// \param distpars [in] array of distribution parameters
/// \param group_index [in] HRU group ID (or DOESNT_EXIST if all HRUs should be perturbed)
/// \param adj [in] type of perturbation (e.g., additive or multiplicative)
//
void CModel::AddForcingPerturbation(forcing_type type, disttype distrib, double* distpars, int group_index, adjustment adj, int nStepsPerDay)
{
  force_perturb *pFP=new force_perturb();
  pFP->forcing     =type;
  pFP->distribution=distrib;
  pFP->kk          =group_index;
  pFP->adj_type    =adj;
  for (int i = 0; i < 3; i++) {
    pFP->distpar[i]=distpars[i];
  }
  pFP->eps=new double [nStepsPerDay];
  for (int n = 0; n < nStepsPerDay; n++) {
    pFP->eps[n]=0;
  }

  if (!DynArrayAppend((void**&)(_pPerturbations),(void*)(pFP),_nPerturbations)){
    ExitGracefully("CModel::AddForcingPerturbation: creating NULL perturbation",BAD_DATA_WARN);
  }

  if ((type==F_PRECIP) || (type==F_SNOWFALL) || (type==F_RAINFALL) || (type==F_TEMP_AVE)){}
  else {
    ExitGracefully("CModel::AddForcingPerturbation: only PRECIP, RAINFALL, SNOWFALL, and TEMP_AVE are supported for forcing perturbation.",BAD_DATA_WARN);
  }
}
/*****************************************************************
   Other Manipulator Functions
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Sets lake storage state variable index based on SV type and layer
///
/// \param *lak_sv [out] SV type
///     \param lev [in] Layer of interest
//
void CModel::SetLakeStorage   (const sv_type      lak_sv, const int lev)
{
  _lake_sv=GetStateVarIndex(lak_sv,lev);
  ExitGracefullyIf(_lake_sv==DOESNT_EXIST,
                   "CModel::SetLakeStorage: non-existent state variable",BAD_DATA);
}

//////////////////////////////////////////////////////////////////
/// \brief Sets HRU group for which storage output files are generated
///
/// \param *pOut [in] (assumed valid) pointer to HRU group
//
void CModel::SetOutputGroup(const CHRUGroup *pOut){
  _pOutputGroup=pOut;
}

//////////////////////////////////////////////////////////////////
/// \brief Sets Ensemble mode for model
///
/// \param *pEnsemble [in] (assumed valid) pointer to Ensemble instance
//
void CModel::SetEnsembleMode(CEnsemble *pEnsemble) {
  _pEnsemble=pEnsemble;
}
//////////////////////////////////////////////////////////////////
/// \brief Sets demand optimizer
///
/// \param *pDO [in] (assumed valid) pointer to optimizer instance
//
void CModel::AddDemandOptimization(CDemandOptimizer *pDO) {
  _pDO=pDO;
}
//////////////////////////////////////////////////////////////////
/// \brief Sets PET blend values
//
void    CModel::SetPETBlendValues(const int N, const evap_method* aEv, const double* wts) {
  _PETBlends_N=N;
  _PETBlends_type=new evap_method [N];
  _PETBlends_wts =new double      [N];
  for (int i = 0; i < N; i++) {
    _PETBlends_type[i]=aEv[i];
    _PETBlends_wts [i]=wts[i];
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Sets PotMelt blend values
//
void    CModel::SetPotMeltBlendValues(const int N, const potmelt_method* aPM, const double* wts)
{
  _PotMeltBlends_N=N;
  _PotMeltBlends_type=new potmelt_method [N];
  _PotMeltBlends_wts =new double         [N];
  for (int i = 0; i < N; i++) {
    _PotMeltBlends_type[i]=aPM[i];
    _PotMeltBlends_wts [i]=wts[i];
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Gets basin blended forcing number of weighted groups
/// \param label [in] String property identifier
/// \return number of blended weights of basin corresponding to label
//

int CModel::GetBlendedForcingsNumWeights(const string label)
{
    string label_n = StringToUppercase(label);
    if      (!label_n.compare("PET_BLEND_WTS"    )) { return _PETBlends_N; }
    else if (!label_n.compare("POTMELT_BLEND_WTS")) { return _PotMeltBlends_N; }
    else {
      ExitGracefully("CModel::GetBlendedForcingsNumWeights: unrecognized blend type",BAD_DATA_WARN);
      return 0;
    }
}
//////////////////////////////////////////////////////////////////
/// \brief Deletes all custom outputs (For FEWS override)
//
void CModel::DeleteCustomOutputs()
{
  for(int i=0;i<_nCustomOutputs;i++) {
    delete _pCustomOutputs[i];
  }
  delete [] _pCustomOutputs;
  _nCustomOutputs=0;
}

//////////////////////////////////////////////////////////////////
/// \brief sets number of layers used to simulate snowpack
/// \param nLayers # of snow layers
//
void  CModel::SetNumSnowLayers     (const int          nLayers)
{
  ExitGracefullyIf(nLayers<0,
                   "CModel::SetNumSnowLayers: cannot set negative number of snow layers",BAD_DATA);
  sv_type *aSV=new sv_type [nLayers];
  int     *aLev=new int [nLayers];
  for (int m=0;m<nLayers;m++){
    aSV[m]=SNOW;
    aLev[m]=m;
  }
  AddStateVariables(aSV,aLev,nLayers);
  delete [] aSV;
  delete [] aLev;
}

bool IsContinuousFlowObs(const CTimeSeriesABC *pObs,long long SBID);

//////////////////////////////////////////////////////////////////
/// \brief overrides streamflow with observed streamflow
/// \param SBID [in] valid subbasin identifier of basin with observations at outflow
//
void CModel::OverrideStreamflow   (const long long SBID)
{
  for (int i=0;i<_nObservedTS; i++)
  {
    if (IsContinuousFlowObs(_pObservedTS[i],SBID))
    {
      //check for blanks in observation TS
      bool bad=false;
      for (int n=0;n<_pObservedTS[i]->GetNumValues();n++){
        if (_pObservedTS[i]->GetValue(n)==RAV_BLANK_DATA){bad=true;break;}
      }
      if (bad){
        WriteWarning("CModel::OverrideStreamflow::cannot override streamflow if there are blanks in observation data",GetOptStruct()->noisy);
        return;
      }

      long long downID=GetSubBasinByID(SBID)->GetDownstreamID();
      if (downID!=DOESNT_EXIST)
      {
        //Copy time series of observed flows to new time series
        string name="Inflow_Hydrograph_"+to_string(SBID);
        CTimeSeries *pObs=dynamic_cast<CTimeSeries *>(_pObservedTS[i]);
        CTimeSeries *pTS =new CTimeSeries(name,*pObs);//copy time series
        pTS->SetLocID(SBID);

        //need to shift everything by one interval, since HYDROGRAPHS are stored as period-ending
        pTS->ShiftInTime(-(pTS->GetInterval()),*_pOptStruct);

        //add as inflow hydrograph to downstream
        GetSubBasinByID(downID)->AddInflowHydrograph(pTS);
        GetSubBasinByID(SBID)->SetDownstreamID(DOESNT_EXIST);
        GetSubBasinByID(SBID)->SetDownstreamBasin(NULL);
        return;
      }
      else{
        WriteWarning("CModel::OverrideStreamflow: overriding streamflow at an outlet subbasin has no impact on model operation",GetOptStruct()->noisy);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns the LU class corresponding to passed string
/// \details Converts string (e.g., "AGRICULTURAL" in HRU file) to LU class
///  can accept either lultclass index or lultclass tag
///  if string is invalid, returns NULL
/// \param s [in] LU class identifier (tag or index)
/// \return Pointer to LU class corresponding to identifier string s
//
CLandUseClass *CModel::StringToLUClass(const string s)
{
  string s_sup = StringToUppercase(s);
  for (int c=0;c<_nLandUseClasses;c++)
  {
    if (!s_sup.compare(StringToUppercase(_pLandUseClasses[c]->GetLanduseName()))) {
      return this->_pLandUseClasses[c];
    }
    else if (s_to_i(s.c_str())==(c+1)) {
      return this->_pLandUseClasses[c];
    }
  }

  return NULL;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns the land use  class corresponding to the passed index
///  if index is invalid, returns NULL
/// \param c [in] LandUse class index
/// \return Reference to land use class corresponding to index c
//
CLandUseClass *CModel::GetLanduseClass(int c) {
  if ((c<0) || (c >= this->_nLandUseClasses)){return NULL;}
  return this->_pLandUseClasses[c];
}

//////////////////////////////////////////////////////////////////
/// \brief Summarize LU class information to screen
//
void CModel::SummarizeLUClassesToScreen()
{
  cout<<"==================="<<endl;
  cout<<"Land Use Class Summary:"<<_nLandUseClasses<<" LU/LT classes in database"<<endl;
  for (int c=0; c<_nLandUseClasses;c++){
    cout<<"-LULT. class \""<<_pLandUseClasses[c]->GetLanduseName()<<"\" "<<endl;
    cout<<"    impermeable: "<<_pLandUseClasses[c]->GetSurfaceStruct()->impermeable_frac*100<<" %"<<endl;
    cout<<"       forested: "<<_pLandUseClasses[c]->GetSurfaceStruct()->forest_coverage*100<<" %"<<endl;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Destroy all LU classes
//
void CModel::DestroyAllLanduseClasses()
{
  if (DESTRUCTOR_DEBUG){cout <<"DESTROYING ALL LULT CLASSES"<<endl;}

  // the classes may have been already destroyed or not created
  if (_nLandUseClasses == 0) {
    if (DESTRUCTOR_DEBUG) {cout << "  No LULT classes to destroy" << endl;}
    return;
  }

  // each class must be destroyed individually, then the array
  for (int c=0; c<_nLandUseClasses;c++){
    delete _pLandUseClasses[c];
  }
  delete [] _pLandUseClasses;

  // the static variables must be reset to avoid dangling pointers and attempts to re-delete
  _pLandUseClasses = NULL;
  _nLandUseClasses = 0;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns the soil class corresponding to passed string
/// \details Converts string (e.g., "SILT" in HRU file) to soilclass
///  can accept either soilclass index or soilclass _tag
///  if string is invalid, returns NULL
/// \param s [in] Soil class identifier (_tag or index)
/// \return Reference to soil class corresponding to identifier s
//
CSoilClass *CModel::StringToSoilClass(const string s)
{
  string sup=StringToUppercase(s);
  for (int c=0;c<_nAllSoilClasses;c++)
  {
    if (!sup.compare(StringToUppercase(_pAllSoilClasses[c]->GetTag()))){return _pAllSoilClasses[c];}
    else if (s_to_i(s.c_str())==(c+1))                                 {return _pAllSoilClasses[c];}
  }
  return NULL;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns the soil class corresponding to the passed index
///  if index is invalid, returns NULL
/// \param c [in] Soil class index
/// \return Reference to soil class corresponding to index c
//
const CSoilClass *CModel::GetSoilClass(int c)
{
  if ((c < 0) || (c >= this->_nAllSoilClasses)){return NULL;}
  return this->_pAllSoilClasses[c];
}

//////////////////////////////////////////////////////////////////
/// \brief Return number of soil classes
/// \return Number of soil classes
//
int CModel::GetNumSoilClasses(){
  return this->_nAllSoilClasses;
}

//////////////////////////////////////////////////////////////////
/// \brief Add a soil class to the model
/// \param pLUClass [in] Pointer to soil class to add
//
void CModel::AddSoilClass(CSoilClass *pSoilClass) {
  if (!DynArrayAppend((void**&)(_pAllSoilClasses), (void*)pSoilClass, _nAllSoilClasses)) {
    ExitGracefully("CModel::AddSoilClass: adding NULL soil class", BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Summarize soil class information to screen
//
void CModel::SummarizeSoilClassesToScreen()
{
  cout<<"==================="<<endl;
  cout<<"Soil Class Summary:"<<this->_nAllSoilClasses<<" soils in database"<<endl;
  for (int c=0; c < this->_nAllSoilClasses;c++){
    cout<<"-Soil class \""<<_pAllSoilClasses[c]->GetTag()<<"\" "<<endl;
    cout<<"       %sand: "<<_pAllSoilClasses[c]->GetSoilStruct()->sand_con<<endl;
    cout<<"       %clay: "<<_pAllSoilClasses[c]->GetSoilStruct()->clay_con<<endl;
    cout<<"    %organic: "<<_pAllSoilClasses[c]->GetSoilStruct()->org_con<<endl;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Destroy all soil classes
//
void CModel::DestroyAllSoilClasses()
{
  if (DESTRUCTOR_DEBUG){cout <<"DESTROYING ALL SOIL CLASSES"<<endl;}

  // the classes may have been already destroyed
  if (_nAllSoilClasses == 0) {
    if (DESTRUCTOR_DEBUG){cout <<"  NO SOIL CLASSES TO DESTROY"<<endl;}
    return;
  }

  // each class must be destroyed individually, then the array
  for (int c=0; c<_nAllSoilClasses;c++){
    delete _pAllSoilClasses[c];
  }
  delete [] _pAllSoilClasses;

  // the static variables must be reset to avoid dangling pointers and attempts to re-delete
  _pAllSoilClasses = NULL;
  _nAllSoilClasses = 0;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns the vegetation class corresponding to passed string
/// \details Converts string (e.g., "BROADLEAF" in HRU file) to vegetation class
///  can accept either vegclass index or vegclass tag
///  if string is invalid, returns NULL
/// \param s [in] Vegetation class identifier (tag or index)
/// \return Reference to vegetation class corresponding to identifier s
//
CVegetationClass *CModel::StringToVegClass(const string s)
{
  string sup;
  sup = StringToUppercase(s);
  for (int c=0; c<this->_numVegClasses; c++)
  {
    if (!sup.compare(StringToUppercase(this->_pAllVegClasses[c]->GetVegetationName()))){return this->_pAllVegClasses[c];}
    else if (s_to_i(sup.c_str()) == (c+1))                                             {return this->_pAllVegClasses[c];}
  }
  return NULL;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns the vegetation class corresponding to the passed index
///  if index is invalid, returns NULL
/// \param c [in] Soil class index
/// \return Reference to vegetation class corresponding to index c
//
const CVegetationClass *CModel::GetVegClass(int c)
{
  if ((c<0) || (c >= this->_numVegClasses)) { return NULL; }
  return this->_pAllVegClasses[c];
}

//////////////////////////////////////////////////////////////////
/// \brief Return number of vegetation classes
/// \return Number of vegetation classes
//
int CModel::GetNumVegClasses(){
  return this->_numVegClasses;
}

//////////////////////////////////////////////////////////////////
/// \brief Add a vegetation to the model
/// \param pLUClass [in] Pointer to vegetation class to add
//
void CModel::AddVegClass(CVegetationClass *pVegClass) {
  if (!DynArrayAppend((void**&)(this->_pAllVegClasses),
                      (void*)pVegClass,
                      this->_numVegClasses)) {
    ExitGracefully("CModel::AddVegClass: adding NULL vegetation class", BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Summarize vegetation class information to screen
//
void CModel::SummarizeVegClassesToScreen()
{
  cout<<"==================="<<endl;
  cout<<"Vegetation Class Summary:"<<this->_numVegClasses<<" vegetation classes in database"<<endl;
  for (int c = 0; c < this->_numVegClasses; c++){
    cout<<"-Veg. class \""<<this->_pAllVegClasses[c]->GetVegetationName()<<"\" "<<endl;
    cout<<"    max. height: "<<this->_pAllVegClasses[c]->GetVegetationStruct()->max_height<<" m"<<endl;
    cout<<"       max. LAI: "<<this->_pAllVegClasses[c]->GetVegetationStruct()->max_LAI  <<endl;
    cout<<"   max. conduct: "<<this->_pAllVegClasses[c]->GetVegetationStruct()->max_leaf_cond<<" mm/s"<<endl;
    cout<<"         albedo: "<<this->_pAllVegClasses[c]->GetVegetationStruct()->albedo<<endl;
  }
}


//////////////////////////////////////////////////////////////////
/// \brief Destroy all vegetation classes
//
void CModel::DestroyAllVegClasses()
{
  if (DESTRUCTOR_DEBUG){cout <<"DESTROYING ALL VEGETATION CLASSES"<<endl;}

  // the classes may have been already destroyed or not created
  if (this->_numVegClasses == 0) {
    if (DESTRUCTOR_DEBUG){cout <<"  NO VEGETATION CLASSES TO DESTROY"<<endl;}
    return;
  }

  // each class must be destroyed individually, then the array
  for (int c=0; c<this->_numVegClasses;c++){
    delete this->_pAllVegClasses[c];
  }
  delete [] this->_pAllVegClasses;

  // the static variables must be reset to avoid dangling pointers and attempts to re-delete
  this->_pAllVegClasses = NULL;
  this->_numVegClasses = 0;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns the terrain class corresponding to passed string
/// \details Converts string (e.g., "HUMMOCKY" in HRU file) to Terrain class
///  can accept either terrainclass index or terrainclass tag
///  if string is invalid, returns NULL
/// \param s [in] terrain class identifier (tag or index)
/// \return Reference to terrain class corresponding to identifier s
//
CTerrainClass *CModel::StringToTerrainClass(const string s)
{
  string sup = StringToUppercase(s);
  for (int c=0; c < this->_nAllTerrainClasses; c++)
  {
    if (!sup.compare(StringToUppercase(this->_pAllTerrainClasses[c]->GetTag()))){return this->_pAllTerrainClasses[c];}
    else if (s_to_i(s.c_str())==(c+1))                                   {return this->_pAllTerrainClasses[c];}
  }
  return NULL;
}

//////////////////////////////////////////////////////////////////
/// \brief Return number of terrain classes
/// \return Number of terrain classes
//
int CModel::GetNumTerrainClasses(){
  return this->_nAllTerrainClasses;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns the terrain class corresponding to the passed index
///  if index is invalid, returns NULL
/// \param c [in] Soil class index
/// \return Reference to terrain class corresponding to index c
//
const CTerrainClass *CModel::GetTerrainClass(int c)
{
  if ((c<0) || (c>=this->_nAllTerrainClasses)){return NULL;}
  return this->_pAllTerrainClasses[c];
}

//////////////////////////////////////////////////////////////////
/// \brief Add a terrain class to the model
/// \param pTerrainClass [in] Pointer to terrain class to add
//
void CModel::AddTerrainClass(CTerrainClass *pTerrainClass){
  if (!DynArrayAppend((void**&)(_pAllTerrainClasses),
                      (void*)pTerrainClass,
                      _nAllTerrainClasses)) {
    ExitGracefully("CModel::AddTerrainClass: creating NULL terrain class", BAD_DATA);
  };
}

//////////////////////////////////////////////////////////////////
/// \brief Summarize terrain class information to screen
//
void CModel::SummarizeTerrainClassesToScreen()
{
  cout<<"==================="<<endl;
  cout<<"Terrain Class Summary:"<<this->_nAllTerrainClasses<<" terrain classes in database"<<endl;
  for (int c=0; c<this->_nAllTerrainClasses; c++){
    cout<<"-Terrain. class \""<<this->_pAllTerrainClasses[c]->GetTag()<<"\" "<<endl;
    cout<<"    drainage density: "<<this->_pAllTerrainClasses[c]->GetTerrainStruct()->drainage_density<<endl;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Destroy all terrain classes
//
void CModel::DestroyAllTerrainClasses()
{
  if (DESTRUCTOR_DEBUG){cout <<"DESTROYING ALL TERRAIN CLASSES"<<endl;}

  // the classes may have been already destroyed or not created
  if (this->_nAllTerrainClasses == 0) {
     if (DESTRUCTOR_DEBUG) {cout <<"  No terrain classes to destroy"<<endl;}
    return;
  }

  // each class must be destroyed individually, then the array
  for (int c=0; c<this->_nAllTerrainClasses;c++){
    delete this->_pAllTerrainClasses[c];
  }
  delete [] this->_pAllTerrainClasses;

  // the static variables must be reset to avoid dangling pointers and attempts to re-delete
  this->_pAllTerrainClasses = NULL;
  this->_nAllTerrainClasses = 0;
}

//////////////////////////////////////////////////////////////////
/// \brief Add a terrain class to the model
/// \param pTerrainClass [in] Pointer to terrain class to add
//
void CModel::AddSoilProfile(CSoilProfile *pSoilProfile){
  if (!DynArrayAppend((void**&)(this->_pAllSoilProfiles),
                      (void*)pSoilProfile,
                      this->_nAllSoilProfiles)) {
    ExitGracefully("CModel::AddSoilProfile: creating NULL terrain class", BAD_DATA);
  };
}

///////////////////////////////////////////////////////////////////
/// \brief Summarize soil profile information to screen
//
void CModel::SummarizeSoilProfilesToScreen()
{
  cout << "===================" << endl;
  cout << "Soil Profile Summary:" << this->_nAllSoilProfiles << " soils in database" << endl;
  for (int p=0; p<this->_nAllSoilProfiles; p++)
  {
    cout << "-Soil profile \"" << this->_pAllSoilProfiles[p]->GetTag() << "\" " << endl;
    cout << "    #horizons:" << this->_pAllSoilProfiles[p]->GetNumHorizons() << endl;
    for (int m=0; m < this->_pAllSoilProfiles[p]->GetNumHorizons(); m++) {
      cout << "      -layer #" << m+1 << ": " << this->_pAllSoilProfiles[p]->GetSoilTag(m) << " (thickness: " << this->_pAllSoilProfiles[p]->GetThickness(m) << " m)" << endl;
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Converts string to soil profile
/// \details Converts string (e.g., "ALL_SILT" in HRU file) to soil profile
///  can accept either soilprofile index or soilprofile tag
///  if string is invalid, returns NULL
/// \param s [in] String identifier of soil profile
/// \return Pointer to Soil profile to which passed string corresponds
//
CSoilProfile *CModel::StringToSoilProfile(const string s)
{
  string sup = StringToUppercase(s);
  for (int p=0; p<this->_nAllSoilProfiles; p++)
  {
    if (!sup.compare(StringToUppercase(this->_pAllSoilProfiles[p]->GetTag()))) {
      return this->_pAllSoilProfiles[p];
    }
    else if (s_to_i(s.c_str())==(p+1)) {
      return this->_pAllSoilProfiles[p];
    }
  }
  return NULL;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of soil profiles in model
/// \return Number of soil profiles in model
//
int CModel::GetNumSoilProfiles()
{
  return (this->_nAllSoilProfiles);
}

//////////////////////////////////////////////////////////////////
/// \brief Destroy all soil profiles in model
//
void CModel::DestroyAllSoilProfiles()
{
  if (DESTRUCTOR_DEBUG){cout <<"DESTROYING ALL SOIL PROFILES"<<endl;}

  // the classes may have been already destroyed or not created
  if (this->_nAllSoilProfiles == 0) {
    if (DESTRUCTOR_DEBUG){cout <<"  No soil profiles to destroy"<<endl; }
    return;
  }

  // each class must be destroyed individually, then the array
  for (int p=0; p<this->_nAllSoilProfiles; p++){
    delete this->_pAllSoilProfiles[p];
  }
  delete [] this->_pAllSoilProfiles;

  // reset the static variables
  this->_pAllSoilProfiles = NULL;
  this->_nAllSoilProfiles = 0;
}

//////////////////////////////////////////////////////////////////
/// \brief Converts string (e.g., "X2305" in Basin file) to channel profile
/// \param s [in] String to be converted to a channel profile
/// \return Pointer to channel profile with name equivalent to string, or NULL if string is invalid
//
CChannelXSect* CModel::StringToChannelXSect(const string s)
{
  string sup = StringToUppercase(s);
  for (int p=0; p<_nAllChannelXSects; p++)
  {
    if (!sup.compare(StringToUppercase(_pAllChannelXSects[p]->GetName()))){
      return _pAllChannelXSects[p];
    }
  }
  return NULL;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of channel cross-sections in model
/// \return Number of channel cross-sections in model
//
int CModel::GetNumChannelXSects() {
  return _nAllChannelXSects;
}

//////////////////////////////////////////////////////////////////
/// \brief Add a cross section to the model
/// \param pXSect [in] Pointer to cross section to add
//
void CModel::AddChannelXSect(CChannelXSect *pXSect)
{
  if (!DynArrayAppend((void**&)(_pAllChannelXSects),(void*)pXSect,_nAllChannelXSects)) {
    ExitGracefully("CModel::AddChannelXSect: creating NULL cross section", BAD_DATA);
  };
}

//////////////////////////////////////////////////////////////////
/// \brief Summarize profile information to screen
//
void CModel::SummarizeChannelXSectToScreen()
{
  cout << "===================" << endl;
  cout << "Channel Profile Summary:" << this->_nAllChannelXSects << " profiles in database" << endl;
  for (int p=0; p<this->_nAllChannelXSects; p++)
  {
    cout << "-Channel profile \"" << this->_pAllChannelXSects[p]->GetName() << "\" " << endl;
    cout << "           slope: " << this->_pAllChannelXSects[p]->GetBedslope()  << endl;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Deletes all channel profiles in model
//
void CModel::DestroyAllChannelXSections()
{
  if (DESTRUCTOR_DEBUG){cout <<"DESTROYING ALL CHANNEL PROFILES"<<endl;}

  // the classes may have been already destroyed or not created
  if (this->_nAllChannelXSects == 0) {
    if (DESTRUCTOR_DEBUG){cout <<"  NO CHANNEL PROFILES TO DESTROY"<<endl;}
    return;
  }

  // each class must be destroyed individually, then the array
  for (int p=0; p < this->_nAllChannelXSects; p++){
    delete this->_pAllChannelXSects[p];
  }
  delete [] this->_pAllChannelXSects;

  // reset the static variables
  this->_pAllChannelXSects = NULL;
  this->_nAllChannelXSects = 0;
}

//////////////////////////////////////////////////////////////////
/// \brief Check for duplicate channel names
//
void CModel::CheckForChannelXSectsDuplicates(const optStruct &Options)
{
  for(int p=0; p<this->_nAllChannelXSects; p++)
  {
    for(int pp=0; pp<p;pp++)
    {
      if(this->_pAllChannelXSects[p]->GetName() == this->_pAllChannelXSects[pp]->GetName()) {
        string warn = " CModel::CheckForChannelXSectsDuplicates: found duplicated channel name: " + this->_pAllChannelXSects[p]->GetName();
        WriteWarning(warn.c_str(), Options.noisy);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Write rating curves to file rating_curves.csv
//
void CModel::WriteRatingCurves(const optStruct& Options) const
{
  ofstream CURVES;
  string tmpFilename = FilenamePrepare("rating_curves.csv", Options);
  CURVES.open(tmpFilename.c_str());

  if (CURVES.fail()){
    ExitGracefully("CModel::WriteRatingCurves: Unable to open output file rating_curves.csv for writing.", FILE_OPEN_ERR);
  }
  int i;
  for (int p=0; p<this->_nAllChannelXSects; p++)
  {
    const CChannelXSect *pP = this->_pAllChannelXSects[p];
    CURVES<<pP->GetName() <<"----------------"<<endl;
    CURVES<<"Flow Rate [m3/s],";    for(i=0;i<pP->GetNPoints();i++){CURVES<<pP->GetAQAt(i)       <<",";}CURVES<<endl;
    CURVES<<"Stage Height [m],";    for(i=0;i<pP->GetNPoints();i++){CURVES<<pP->GetAStageAt(i)   <<",";}CURVES<<endl;
    CURVES<<"Top Width [m],";       for(i=0;i<pP->GetNPoints();i++){CURVES<<pP->GetATopWidthAt(i)<<",";}CURVES<<endl;
    CURVES<<"X-sect area [m2],";    for(i=0;i<pP->GetNPoints();i++){CURVES<<pP->GetAXAreaAt(i)   <<",";}CURVES<<endl;
    CURVES<<"Wetted Perimeter [m],";for(i=0;i<pP->GetNPoints();i++){CURVES<<pP->GetAPerimAt(i)   <<",";}CURVES<<endl;
  }
  CURVES.close();
}

//////////////////////////////////////////////////////////////////
/// \brief Returns the number of convolution variables
/// \return Number of convolution variables
//
int CModel::GetNumConvolutionVariables() const {
  return _nConvVariables;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns the number of convolution variables
//
void CModel::IncrementConvolutionCount(){
  _nConvVariables++;
}

//////////////////////////////////////////////////////////////////
/// \brief Gets the state variables
/// \return Pointer to state variables
//
CStateVariable* CModel::GetStateVarInfo() const {
  return _pStateVar;
}

//////////////////////////////////////////////////////////////////
/// \brief Sets the state variables
/// \param pStateVar [in] Pointer to state variables
//
void CModel::SetStateVarInfo(CStateVariable *pStateVar) {
  _pStateVar = pStateVar;
}

//////////////////////////////////////////////////////////////////
/// \brief Gets the number of lateral flow processes
/// \return Number of state variables
//
int CModel::GetNumLatFlowProcesses() {
  return this->_nLatFlowProcesses;
}

//////////////////////////////////////////////////////////////////
/// \brief Increments the number of lateral flow processes
//
void CModel::CountOneMoreLatFlowProcess() {
  this->_nLatFlowProcesses++;
}

/*****************************************************************
   Routines called repeatedly during model simulation
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Increments water/mass/energy balance
/// \details add cumulative amount of water/mass/energy for each process connection
/// \note  called from main solver
///
/// \param q_star [in] Integer index of hydrologic process connection
/// \param k [in] Integer index of HRU
/// \param moved [in] Double amount balance is to be incremented (mm or mg/m2 or MJ/m2)
//
void CModel::IncrementBalance(const int    q_star,
                              const int    k,
                              const double moved){
#ifdef _STRICTCHECK_
  ExitGracefullyIf(q_star>_nTotalConnections,"CModel::IncrementBalance: bad index",RUNTIME_ERR);
#endif
  _aCumulativeBal[k][q_star]+=moved;
  _aFlowBal      [k][q_star]= moved;
}
//////////////////////////////////////////////////////////////////
/// \brief Increments lateral flow water/energy balance
/// \details add cumulative amount of water/energy for each lateral process connection
/// \note  called from main solver
///
/// \param jss [in] Integer index of hydrologic lateral process connection
/// \param moved [in] Double amount balance is to be incremented [mm-m2] or [MJ] or [mg]
//
void CModel::IncrementLatBalance( const int jss,
                                  const double moved)//[mm-m2] or [MJ]
{
  _aCumulativeLatBal[jss]+=moved;
  _aFlowLatBal      [jss]=moved;
}
//////////////////////////////////////////////////////////////////
/// \brief Increments cumulative mass & energy added to system (precipitation/ustream basin flows, etc.)
/// \details Increment cumulative precipitation based on average preciptation and corresponding timestep [mm]
///
/// \param &Options [in] Global model options information
/// \param &tt [in] current time
//
void CModel::IncrementCumulInput(const optStruct &Options, const time_struct &tt)
{
  double area;
  area = _WatershedArea * M2_PER_KM2;

  _CumulInput+=GetAveragePrecip()*Options.timestep;

  _CumulInput+=GetAverageForcings().recharge*Options.timestep;

  for (int p=0;p<_nSubBasins;p++){
    _CumulInput+=_pSubBasins[p]->GetIntegratedSpecInflow(tt.model_time,Options.timestep)/area*MM_PER_METER;//converted to [mm] over  basin
  }

  //add from groundwater
  if(Options.modeltype!=MODELTYPE_SURFACE) {
    int iGW=GetStateVarIndex(GROUNDWATER);
    for(int k=0;k<_nHydroUnits;k++) {
      double GW=_pHydroUnits[k]->GetStateVarValue(iGW);
      area=_pHydroUnits[k]->GetArea();
      if (GW<0){_CumulInput-=GW*area/_WatershedArea;} //negative recharge
    }
    /*for(int p=0;p<_nSubBasins;p++) {
      _CumulInput+=max(pGW2River->CalcRiverFluxBySB(p),0.0)*Options.timestep*area/MM_PER_METER; ///m3/d baseflow
    }*/
  }

  _pTransModel->IncrementCumulInput(Options,tt);
}

//////////////////////////////////////////////////////////////////
/// \brief Increments cumulative outflow from system for mass balance diagnostics
/// \details Increment cumulative outflow according to timestep, flow, and area of basin
/// \details Called every timestep after MassEnergyBalance
/// \param &Options [in] Global model options information
//
void CModel::IncrementCumOutflow(const optStruct &Options, const time_struct &tt)
{
  int pDivert;
  double Qdiv;
  double area=(_WatershedArea*M2_PER_KM2);
  for(int p=0;p<_nSubBasins;p++)
  {
    bool outflowdisabled;
    CSubBasin *pSBdown=NULL;
    if(_aDownstreamInds[p]>=0) { pSBdown=_pSubBasins[_aDownstreamInds[p]]; }

    outflowdisabled=((pSBdown==NULL) || (!pSBdown->IsEnabled()));

    if(_pSubBasins[p]->IsEnabled())
    {
      if((_aSubBasinOrder[p]==0) || (outflowdisabled))//outlet does not drain into another subbasin
      {
      _CumulOutput+=_pSubBasins[p]->GetIntegratedOutflow(Options.timestep)/area*MM_PER_METER;//converted to [mm] over entire watershed
      }
      _CumulOutput+=_pSubBasins[p]->GetReservoirLosses (Options.timestep)/area*MM_PER_METER;
      _CumulOutput+=_pSubBasins[p]->GetIrrigationLosses(Options.timestep)/area*MM_PER_METER;
      _CumulOutput+=_pSubBasins[p]->GetDiversionLosses (Options.timestep)/area*MM_PER_METER; //includes all diverted water, whether or not it stays in watershed

      //Diversion losses that were actually just redirected to other basins
      for(int i=0;i<_pSubBasins[p]->GetNumDiversions();i++)
      {
        Qdiv=_pSubBasins[p]->GetDiversionFlow(i,_pSubBasins[p]->GetLastChannelOutflowRate(),Options,tt,pDivert)*Options.timestep*SEC_PER_DAY;
        if(pDivert!=DOESNT_EXIST) {
          _CumulOutput-=Qdiv/area*MM_PER_METER; //water that is not diverted out of the watershed
        }
      }
    }
  }
  //lost to groundwater
  if(Options.modeltype!=MODELTYPE_SURFACE) {
    int iGW=GetStateVarIndex(GROUNDWATER);
    for(int k=0;k<_nHydroUnits;k++) {
      double GW=_pHydroUnits[k]->GetStateVarValue(iGW);
      area=_pHydroUnits[k]->GetArea();
      if(GW>0) { _CumulOutput+=GW*area/_WatershedArea; }
    }
    /*for(int p=0;p<_nSubBasins;p++) {
      _CumulOutput-=min(pGW2River->CalcRiverFluxBySB(p),0.0)*Options.timestep*area/MM_PER_METER; ///mm baseflow
    }*/
  }
  _pTransModel->IncrementCumulOutput(Options);
}

//////////////////////////////////////////////////////////////////
/// \brief Updates values of user-specified transient parameters, updates changes to land use class
///
/// \param &Options [in] Global model options information
/// \param &tt [in] Current time structure
//
void CModel::UpdateTransientParams(const optStruct   &Options,
                                   const time_struct &tt)
{
  //--update parameters linked to time series-----------------------------------------------
  int nn=(int)((tt.model_time+REAL_SMALL)/Options.timestep);//current timestep index
  for (int j=0;j<_nTransParams;j++)
  {
    class_type ctype=_pTransParams[j]->GetParameterClassType();
    string     pname=_pTransParams[j]->GetParameterName();
    string     cname=_pTransParams[j]->GetParameterClass();
    double     value=_pTransParams[j]->GetTimeSeries()->GetSampledValue(nn);

    UpdateParameter(ctype,pname,cname,value);
  }

  //--update land use and HRU types-----------------------------------------------
  int k;
  for (int j = 0; j<_nClassChanges; j++)
  {
    if( ((_pClassChanges[j]->modeltime > tt.model_time - TIME_CORRECTION) &&
         (_pClassChanges[j]->modeltime < tt.model_time + Options.timestep)) ||
	    ((tt.model_time == 0.0) && (_pClassChanges[j]->modeltime < 0.0)) )
    {//change happens this time step

      //cout<<"updating classes on "<<tt.date_string<< endl;
      int kk   =_pClassChanges[j]->HRU_groupID;
      for(int k_loc = 0; k_loc <_pHRUGroups[kk]->GetNumHRUs();k_loc++)
      {
        k=_pHRUGroups[kk]->GetHRU(k_loc)->GetGlobalIndex();

        if      (_pClassChanges[j]->tclass == CLASS_LANDUSE)
        {
          CLandUseClass *lult_class = StringToLUClass(_pClassChanges[j]->newclass);
          _pHydroUnits[k]->ChangeLandUse(lult_class);
        }
        else if (_pClassChanges[j]->tclass == CLASS_VEGETATION)
        {
          CVegetationClass *veg_class = StringToVegClass(_pClassChanges[j]->newclass);
          _pHydroUnits[k]->ChangeVegetation(veg_class);
        }
        else if (_pClassChanges[j]->tclass == CLASS_HRUTYPE)
        {
          HRU_type typ=StringToHRUType(_pClassChanges[j]->newclass);
          _pHydroUnits[k]->ChangeHRUType(typ);
        }

        for(int jj=0; jj<_nProcesses;jj++)// kt
        {
          _aShouldApplyProcess[jj][k] = _pProcesses[jj]->ShouldApply(_pHydroUnits[k]);
        }
      }
    }
  }

}
//////////////////////////////////////////////////////////////////
/// \brief Determines parameter class (e.g., CLASS_GLOBAL or CLASS_SOIL) from parameter name and class name through slow search
///
/// \param param_str [in] string for parameter name (e.g., BASEFLOW_COEFF)
/// \param class_name [in] class name (e.g., FOREST) or SBID for subbasin parameters
/// \return parameter class enum, or CLASS_UNKNOWN if the parameter is not found
//
class_type CModel::ParamNameToParamClass(const string param_str, const string class_name) const
{
  double pval;
  class_type pclass=CLASS_UNKNOWN;
  soil_struct    S;
  veg_struct     V;
  surface_struct L;
  global_struct  G;
  pval=CSoilClass::GetSoilProperty(S,param_str,false);
  if (pval!=INDEX_NOT_FOUND){return CLASS_SOIL;}
  pval=CVegetationClass::GetVegetationProperty(V,param_str,false);
  if (pval!=INDEX_NOT_FOUND){return CLASS_VEGETATION;}
  pval = CLandUseClass::GetSurfaceProperty(L, param_str, false);
  if (pval!=INDEX_NOT_FOUND){return CLASS_LANDUSE;}
  pval = _pGlobalParams->GetGlobalProperty(G, param_str, false);
  if (pval!=INDEX_NOT_FOUND){return CLASS_GLOBAL;}
  if(_nGauges>0) {
    pval=_pGauges[0]->GetGaugeProperty(param_str);
    if(pval!=INDEX_NOT_FOUND) {return CLASS_GAUGE; }
  }

  CSubBasin *pSB=new CSubBasin(0,"",this,-1,NULL,0,AUTO_COMPUTE,false,false); //fake basin (basins not yet created)
  pval=pSB->GetBasinProperties(param_str);
  if (pval!=INDEX_NOT_FOUND){pclass=CLASS_SUBBASIN;}
  delete pSB;

  return pclass;
}
//////////////////////////////////////////////////////////////////
/// \brief adds land use class
//
void  CModel::AddLandUseClass(CLandUseClass* pLU)
{
  if (!DynArrayAppend((void**&)(_pLandUseClasses),(void*)(pLU),_nLandUseClasses)) {
    ExitGracefully("CLandUseClass::Constructor: creating NULL land use class",BAD_DATA);};
}
//////////////////////////////////////////////////////////////////
/// \brief Returns the LU class corresponding to passed string
/// \details Converts string (e.g., "AGRICULTURAL" in HRU file) to LU class
///  if string is invalid, returns NULL
/// \param s [in] LU class name
/// \return Pointer to LU class corresponding to identifier string s
//
const CLandUseClass *CModel::StringToLUClass(const string s) const
{
  int c=GetLandClassIndex(s);
  if (c==DOESNT_EXIST){return NULL;}
  else                {return _pLandUseClasses[c]; }
}
const int CModel::GetLandClassIndex(const string s) const
{
  string sup=StringToUppercase(s);
  for (int c=0;c<_nLandUseClasses;c++)
  {
    if (!sup.compare(StringToUppercase(_pLandUseClasses[c]->GetLanduseName()))){return c;}
  }
  return DOESNT_EXIST;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns the land use  class corresponding to the passed index
///  if index is invalid, returns NULL
/// \param c [in] land class index
/// \return Reference to land use class corresponding to index c
//
const CLandUseClass *CModel::GetLanduseClass(const int c) const
{
  if ((c<0) || (c>=_nLandUseClasses)){return NULL;}
  return _pLandUseClasses[c];
}
//////////////////////////////////////////////////////////////////
/// \brief Returns number of land use classes
//
int CModel::GetNumLanduseClasses() const {
  return _nLandUseClasses;
}
//////////////////////////////////////////////////////////////////
/// \brief Updates model parameter during course of simulation
///
/// \param &ctype [in] parameter class type
/// \param &pname [in] valid parameter name
/// \param &cname [in] valid parameter class name (or SBID as string for CLASS_SUBBASIN or gauge ID as string for CLASS_GAUGE)
/// \param &value [in] updated parameter value (or -1.2345, which is ignored)
//
void CModel::UpdateParameter(const class_type &ctype,const string pname,const string cname,const double &value)
{
  if (value==RAV_BLANK_DATA) { return; }

  if (ctype==CLASS_SOIL)
  {
    StringToSoilClass(cname)->SetSoilProperty(pname, value);
  }
  else if(ctype==CLASS_TRANSPORT)
  {
    ExitGracefully("UpdateParameter::CLASS_TRANSPORT", STUB);
    //CSoilClass::StringToSoilClass(cname)->SetTransportProperty(constit,pname,value);
  }
  else if(ctype==CLASS_VEGETATION)
  {
    StringToVegClass(cname)->SetVegetationProperty(pname, value);
  }
  else if(ctype==CLASS_TERRAIN)
  {
    StringToTerrainClass(cname)->SetTerrainProperty(pname, value);
  }
  else if(ctype==CLASS_LANDUSE)
  {
    StringToLUClass(cname)->SetSurfaceProperty(pname, value);
  }
  else if(ctype==CLASS_GLOBAL)
  {
    _pGlobalParams->SetGlobalProperty(pname, value);
  }
  else if(ctype==CLASS_GAUGE)
  {
    int g=GetGaugeIndexFromName(cname);
    if(g!=DOESNT_EXIST) {
      _pGauges[g]->SetGaugeProperty(pname, value);
    }
    else {
      WriteWarning("CModel::UpdateParameter: Unrecognized/invalid gauge ID ("+cname+") in input",false);
    }
  }
  else if(ctype==CLASS_SUBBASIN)
  {
    if (_nSubBasins == 0) {
      ExitGracefully("UpdateParameter: should not call before subbasins are created",RUNTIME_ERR);
    }
    long long SBID=s_to_ll(cname.c_str()); //class name should be SBID in this case
    if( (strlen(cname.c_str())>8) && //also accept SUBBASIN32 instead of 32
        (!strcmp(cname.substr(0,8).c_str(),"SUBBASIN")) ) {
      SBID=s_to_ll(cname.substr(8,strlen(cname.c_str())-8).c_str());
    }

    CSubBasin *pSB=GetSubBasinByID(SBID);
    if (pSB != NULL) {
      pSB->SetBasinProperties(pname,value);
    }
    else {
      WriteWarning("CModel::UpdateParameter: Unrecognized/invalid subbasin ID ("+to_string(SBID)+") in input",false);
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief overrides global parameters in subbasin groups
/// \notes called only from solver within HRU loop
///
/// \param k [in] global HRU index
/// \param revert [in] true if reverting to base value, false if changing to overridden value
//
void CModel::ApplyLocalParamOverrrides(const int k, const bool revert)
{
  for (int i=0;i<_nParamOverrides;i++)
  {
    if (_pParamOverrides[i]->aHRUIsOverridden[k]) {
      if (!revert){
        for (int j=0; j<_pParamOverrides[i]->nVals; j++){
          _pParamOverrides[i]->pxAddress[j] = _pParamOverrides[i]->aValues[j];

        }
      }
      else{
        for (int j=0; j<_pParamOverrides[i]->nVals; j++){
          _pParamOverrides[i]->pxAddress[j] = _pParamOverrides[i]->aRevertValues[j];
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Recalculates HRU derived parameters
/// \details Recalculate HRU derived parameters that are based upon time-of-year/day and SVs (storage, temp)
/// \remark tt is model time, not global time
///
/// \param &Options [in] Global model options information
/// \param &tt [in] Current time structure
//
void CModel::RecalculateHRUDerivedParams(const optStruct    &Options,
                                         const time_struct  &tt)
{
  for (int k=0;k<_nHydroUnits;k++)
  {
    if(_pHydroUnits[k]->IsEnabled())
    {
      _pHydroUnits[k]->RecalculateDerivedParams(Options,tt);
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief called at start of time step if needed - generates random values for perturbation of forcings
/// \param &Options [out] Global model options information
/// \params tt [in] time structure
//
void CModel::PrepareForcingPerturbation(const optStruct &Options, const time_struct &tt)
{

  for(int i=0;i<_nPerturbations;i++)
  {
    int    nStepsPerDay = (int)(rvn_round(1.0/Options.timestep));
    double partday      = Options.julian_start_day-floor(Options.julian_start_day+TIME_CORRECTION);
    int    nn           = (int)(rvn_round((tt.model_time+partday-floor(tt.model_time+partday+TIME_CORRECTION))/Options.timestep));
    bool   start_of_day = ((nn==0) || tt.day_changed); //nn==0 corresponds to midnight

    if (start_of_day)  { //get all random perturbation samples for the day
      for (int n=0;n<nStepsPerDay;n++){
        _pPerturbations[i]->eps[n]=SampleFromDistribution(_pPerturbations[i]->distribution,_pPerturbations[i]->distpar);
      }
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief called within update forcings - actually applies forcing function changes
/// \param ftype [in] forcing function type
/// \param F [out] forcing structure modified by this routine
/// \param k [in] HRU index of forcing
/// \param &Options [in] Global model options information
/// \params tt [in] time structure
//
void CModel::ApplyForcingPerturbation(const forcing_type ftype, force_struct &F, const int k, const optStruct& Options, const time_struct& tt)
{
  int kk=DOESNT_EXIST;
  for(int i=0;i<_nPerturbations;i++)
  {
    if (ftype==_pPerturbations[i]->forcing)
    {
      int    nStepsPerDay = (int)(rvn_round(1.0/Options.timestep));
      double partday      = Options.julian_start_day-floor(Options.julian_start_day+TIME_CORRECTION);
      int    nn           = (int)(rvn_round((tt.model_time+partday-floor(tt.model_time+partday+TIME_CORRECTION))/Options.timestep));
      bool   start_of_day = ((nn==0) || tt.day_changed); //nn==0 corresponds to midnight

      kk=_pPerturbations[i]->kk;

      if((kk==DOESNT_EXIST) || (GetHRUGroup(kk)->IsInGroup(k)))
      {
        _pHydroUnits[k]->AdjustHRUForcing(ftype, F,_pPerturbations[i]->eps[nn], _pPerturbations[i]->adj_type);
        if (start_of_day){
          _pHydroUnits[k]->AdjustDailyHRUForcings(ftype,F,_pPerturbations[i]->eps,_pPerturbations[i]->adj_type,nStepsPerDay);
        }
      }
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Updates values stored in modeled time series of observation data
/// modifies _pModeledTS[] time series and _aObsIndex array
/// \param &Options [in] Global model options information
/// \param &tt [in] Current time structure
//
void CModel::UpdateDiagnostics(const optStruct   &Options,
                               const time_struct &tt)
{
  //if (_nDiagnostics==0){return;}

  int n=(int)(floor((tt.model_time+TIME_CORRECTION)/Options.timestep));//current timestep index

  double     value, obsTime;
  string     datatype;
  sv_type    svtyp;
  int        layer_ind;
  CSubBasin *pBasin=NULL;

  for (int i=0;i<_nObservedTS;i++)
  {
    datatype = _pObservedTS[i]->GetName();
    svtyp    = _pStateVar->StringToSVType(datatype, layer_ind, false);

    if      (datatype=="HYDROGRAPH")//===============================================
    {
      pBasin=GetSubBasinByID (_pObservedTS[i]->GetLocID());

      if ((Options.ave_hydrograph) && (tt.model_time!=0)){
        value=pBasin->GetIntegratedOutflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);
      }
      else{
        value=pBasin->GetOutflowRate();
      }
    }
    else if (datatype == "RESERVOIR_STAGE")//=======================================
    {
      pBasin=GetSubBasinByID(_pObservedTS[i]->GetLocID());
      value = pBasin->GetReservoir()->GetResStage();
    }
    else if (datatype == "RESERVOIR_INFLOW")//======================================
    {
      pBasin =  GetSubBasinByID(_pObservedTS[i]->GetLocID());
      value = pBasin->GetIntegratedReservoirInflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);
    }
    else if (datatype == "RESERVOIR_NETINFLOW")//===================================
    {
      pBasin = GetSubBasinByID(_pObservedTS[i]->GetLocID());
      CReservoir *pRes= pBasin->GetReservoir();
      double avg_area=0.0;
      if (pRes->GetHRUIndex()!=DOESNT_EXIST){ avg_area = _pHydroUnits[pRes->GetHRUIndex()]->GetArea(); }

      double tem_precip1 = pBasin->GetAvgForcing(F_PRECIP) / (Options.timestep*SEC_PER_DAY)*avg_area*M2_PER_KM2/MM_PER_METER;
      double losses      = pRes->GetReservoirEvapLosses        (Options.timestep) / (Options.timestep*SEC_PER_DAY);
      losses            += pRes->GetReservoirGWLosses          (Options.timestep) / (Options.timestep*SEC_PER_DAY);
      value              = pBasin->GetIntegratedReservoirInflow(Options.timestep) / (Options.timestep*SEC_PER_DAY) + tem_precip1 - losses;
    }
    else if(datatype == "STREAM_CONCENTRATION")//=======================================
    {
      int c=_pObservedTS[i]->GetConstitInd();
      int p=GetSubBasinIndex(_pObservedTS[i]->GetLocID());
      if (c==DOESNT_EXIST){value=RAV_BLANK_DATA;}
      else                {value = _pTransModel->GetConstituentModel2(c)->GetOutflowConcentration(p);}
    }
    else if(datatype == "STREAM_TEMPERATURE")//=======================================
    {
      int c=_pObservedTS[i]->GetConstitInd();
      int p=GetSubBasinIndex(_pObservedTS[i]->GetLocID());
      if (c==DOESNT_EXIST){value=RAV_BLANK_DATA;}
      else                {value = _pTransModel->GetConstituentModel2(c)->GetOutflowConcentration(p);}
    }
    else if(datatype == "WATER_LEVEL")//=======================================
    {
      pBasin=GetSubBasinByID(_pObservedTS[i]->GetLocID());
      value = pBasin->GetWaterLevel();
    }
    else if (datatype == "LAKE_AREA")//========================================
    {
      pBasin=GetSubBasinByID(_pObservedTS[i]->GetLocID());
      value = pBasin->GetReservoir()->GetSurfaceArea();
    }
    else if (svtyp!=UNRECOGNIZED_SVTYPE)//==========================================
    { //State variable
      CHydroUnit *pHRU=NULL;
      pHRU=GetHRUByID((long long int)(_pObservedTS[i]->GetLocID()));
      string error="CModel::UpdateDiagnostics: Invalid HRU ID specified in observed state variable time series "+_pObservedTS[i]->GetName();
      ExitGracefullyIf(pHRU==NULL,error.c_str(),BAD_DATA);
      int sv_index=this->GetStateVarIndex(svtyp,layer_ind);
      value       =pHRU->GetStateVarValue(sv_index);
    }
    else// ERROR ===================================================================
    {
      if (tt.model_time ==0){
        string warn="CModel::UpdateDiagnostics: invalid tag ("+datatype+")used for specifying Observation type";
        WriteWarning(warn.c_str(),BAD_DATA);
      }
      value=0;
    }

    _pModeledTS[i]->SetValue(n,value);
    _pModeledTS[i]->SetSampledValue(n,value); //Handles blank value issue in final time  step

    obsTime =_pObservedTS[i]->GetSampledTime(_aObsIndex[i]); // time of the next observation
    //if(_pObservedTS[i]->GetType()==CTimeSeriesABC::TS_IRREGULAR) {obsTime+=Options.timestep;}//JRC: Not the most elegant fix, but Diagnostics are updated prior to Solver (no longer)

    //if model time is such that next unprocessed observation has occurred, update modeled
    while ((tt.model_time >= obsTime+ _pObservedTS[i]->GetSampledInterval()) &&
			     (_aObsIndex[i]<_pObservedTS[i]->GetNumSampledValues()))
		{
      value=RAV_BLANK_DATA;
      // only set values within diagnostic evaluation times. The rest stay as BLANK_DATA
      if ((obsTime >= Options.diag_start_time) && (obsTime <= Options.diag_end_time))
      {
        value= _pModeledTS[i]->GetModelledValue(obsTime,_pObservedTS[i]->GetType());
      }
      _pModeledTS[i]->SetSampledValue(_aObsIndex[i],value); //JRC : \todo[fix]: I'm not sure if this makes sense for irregular time series
      _aObsIndex[i]++;
      obsTime =_pObservedTS[i]->GetSampledTime(_aObsIndex[i]);
      //if(_pObservedTS[i]->GetType()==CTimeSeriesABC::TS_IRREGULAR) {obsTime+=Options.timestep;}
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Apply hydrological process to model
/// \details Method returns rate of mass/energy transfers rates_of_change [mm/d, mg/m2/d, or MJ/m2/d] from a set
/// of state variables iFrom[] to a set of state variables iTo[] (e.g., expected water movement [mm/d]
// from storage unit iFrom[0] to storage unit iTo[0]) in the given HRU using state_var[] as the expected
/// value of all state variables over the timestep.
///
/// \param j        [in] Integer process indentifier
/// \param *state_var [in] Array of state variables for HRU
/// \param *pHRU    [in] Pointer to HRU
/// \param &Options [in] Global model options information
/// \param &tt      [in] Time structure
/// \param *iFrom   [in] Array (size: nConnections)  of Indices of state variable losing mass or energy
/// \param *iTo     [in] Array (size: nConnections)  of indices of state variable gaining mass or energy
/// \param &nConnections [out] Number of connections between storage units/state vars
/// \param *rates_of_change [out] Double array (size: nConnections) of loss/gain rates of water [mm/d], mass [mg/m2/d], and/or energy [MJ/m2/d]
/// \return returns false if this process doesn't apply to this HRU, true otherwise
//
bool CModel::ApplyProcess ( const int          j,                    //process identifier
                            const double      *state_var,            //array of state variables for HRU
                            const CHydroUnit  *pHRU,                 //pointer to HRU
                            const optStruct   &Options,
                            const time_struct &tt,
                                  int         *iFrom,                //indices of state variable losing water or heat
                                  int         *iTo,                  //indices of state variable gaining water or heat
                                  int         &nConnections,         //number of connections between storage units/state vars
                                  double      *rates_of_change) const//loss/gain rates of water [mm/d] and energy [MJ/m2/d]
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((j<0) && (j>=_nProcesses),"CModel ApplyProcess::improper index",BAD_DATA);
#endif
  CHydroProcessABC *pProc=_pProcesses[j];

  nConnections=pProc->GetNumConnections();//total connections: nConnections
  int k=pHRU->GetGlobalIndex();
  if (!_aShouldApplyProcess[j][k]){return false;}

  for (int q=0;q<nConnections;q++)
  {
    iFrom[q]=pProc->GetFromIndices()[q];
    iTo  [q]=pProc->GetToIndices  ()[q];
    rates_of_change[q]=0.0;
  }

  pProc->GetRatesOfChange(state_var,pHRU,Options,tt,rates_of_change);

  //Special frozen flow handling - constrains flows when water is partially/wholly frozen
  //------------------------------------------------------------------------
  if (_pTransModel->GetEnthalpyModel()!=NULL)
  { //simulating enthalpy, and therefore frozen water compartments
    //if ((tt.model_time==0.0) && (j==0)){WriteWarning("JAMES: Temporarily disabled frozen ground feedback",Options.noisy); }
  /*
    double Fi,liq_stor;
    sv_type typ;
    for(int q=0;q<nConnections;q++)
    {
      typ=_aStateVarType[iFrom[q]];
      if(CStateVariable::IsWaterStorage(typ) && (typ!=ATMOS_PRECIP) && (typ!=SNOW) && (typ!=CANOPY_SNOW) && (typ!=SURFACE_WATER))
      {
        //precip is special case (frozen snow can fall)
        //snow is special case (already assume that only liquid water is lost)
        Fi=_pTransModel->GetEnthalpyModel()->GetIceContent(state_var,iFrom[q]);

        liq_stor=(1.0-Fi)*state_var[iFrom[q]]; //[mm]

        rates_of_change[q]=min(rates_of_change[q],liq_stor/Options.timestep);
      }
    }
    */
  }

  //Apply constraints
  //------------------------------------------------------------------------
  pProc->ApplyConstraints(state_var,pHRU,Options,tt,rates_of_change);

  return true;
}
//////////////////////////////////////////////////////////////////
/// \brief Apply lateral exchange hydrological process to model
/// \details Method returns rate of mass/energy transfers rates_of_change [mm/d, mg/m2/d, or MJ/m2/d] from a set
/// of state variables iFrom[] in HRUs kFrom[] to a set of state variables iTo[] in HRUs kTo[] (e.g., expected water movement [mm/d]
// from storage unit iFrom[0] in HRU kFrom[0] to storage unit iTo[0] in HRU kFrom[0]) in the given HRU using state_vars[][] as the expected
/// value of all state variables over the timestep.
///
/// \param j        [in] Integer process indentifier
/// \param **state_vars [in] Array of state variables for all HRUs (size: [nHRUs][nStateVars])
/// \param *pHRU    [in] Pointer to HRU
/// \param &Options [in] Global model options information
/// \param &tt      [in] Time structure
/// \param *kFrom   [in] Array (size: nLatConnections)  of Indices of HRU losing mass or energy
/// \param *kTo     [in] Array (size: nLatConnections)  of indices of HRU gaining mass or energy
/// \param *iFrom   [in] Array (size: nLatConnections)  of Indices of state variable losing mass or energy
/// \param *iTo     [in] Array (size: nLatConnections)  of indices of state variable gaining mass or energy
/// \param &nConnections [out] Number of connections between storage units/state vars
/// \param *exchange_rates [out] Double array (size: nConnections) of loss/gain rates of water [mm-m2/d], mass [mg/d], and/or energy [MJ/d]
/// \return returns false if this process doesn't apply to this HRU, true otherwise
//
bool CModel::ApplyLateralProcess( const int          j,
                                  const double* const* state_vars,
                                  const optStruct   &Options,
                                  const time_struct &tt,
                                        int         *kFrom,
                                        int         *kTo,
                                        int         *iFrom,
                                        int         *iTo,
                                        int         &nLatConnections,
                                        double      *exchange_rates) const
{
  CLateralExchangeProcessABC *pLatProc;

  ExitGracefullyIf((j<0) && (j>=_nProcesses),"CModel ApplyProcess::improper index",BAD_DATA);

  nLatConnections=_pProcesses[j]->GetNumLatConnections();

  if(nLatConnections==0){return false;}

  pLatProc=(CLateralExchangeProcessABC*)_pProcesses[j]; // Cast
  if (!_aShouldApplyProcess[j][pLatProc->GetFromHRUIndices()[0]]){return false;} //JRC: is the From/0 appropriate?


  for (int q=0;q<nLatConnections;q++)
  {
    iFrom[q]=pLatProc->GetLateralFromIndices()[q];
    iTo  [q]=pLatProc->GetLateralToIndices()[q];
    kFrom[q]=pLatProc->GetFromHRUIndices()[q];
    kTo  [q]=pLatProc->GetToHRUIndices()[q];
    exchange_rates[q]=0.0;
  }
  for(int q=0;q<nLatConnections;q++)
  {
    if(!_pHydroUnits[kFrom[q]]->IsEnabled()) { return false; }
    if(!_pHydroUnits[kTo  [q]]->IsEnabled()) { return false; } //ALL participating HRUs must be enabled to apply
  }
  pLatProc->GetLateralExchange(state_vars,_pHydroUnits,Options,tt,exchange_rates);

  return true;
}
