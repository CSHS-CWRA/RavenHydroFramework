/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright © 2008-2014 the Raven Development Team
----------------------------------------------------------------*/
#include "Model.h"

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
CModel::CModel(const soil_model SM, 
               const int        nsoillayers,
               const optStruct &Options)
{
  int i;
  _nSubBasins=0;    _pSubBasins=NULL;
  _nHydroUnits=0;   _pHydroUnits=NULL;
  _nHRUGroups=0;    _pHRUGroups=NULL; 
  _nGauges=0;       _pGauges=NULL;
  _nProcesses=0;    _pProcesses=NULL;
  _nCustomOutputs=0;_pCustomOutputs=NULL;
  _nTransParams=0;  _pTransParams=NULL;
  _nObservedTS=0;   _pObservedTS=NULL; _pModeledTS=NULL; _aObsIndex=NULL;
  _nObsWeightTS =0; _pObsWeightTS=NULL;
  _nDiagnostics=0;  _pDiagnostics=NULL;
  
  _nTotalConnections=0;

  _WatershedArea=0;
  
  _aSubBasinOrder =NULL; _maxSubBasinOrder=0;
  _aOrderedSBind  =NULL;
  _aDownstreamInds=NULL;

  _pOptStruct = &Options;

  ExitGracefullyIf(nsoillayers<1,
    "CModel constructor::improper number of soil layers. SoilModel not specified?",BAD_DATA);

  //Initialize Lookup table for state variable indices
  for (int i=0;i<MAX_STATE_VARS;i++){
    for (int m=0;m<MAX_SV_LAYERS;m++){
      _aStateVarIndices[i][m]=DOESNT_EXIST;
    }
  }
  //determine first group of state variables based upon soil model
  //SW, atmosphere, atmos_precip always present, one for each soil layer, and 1 for GW (unless lumped)
  _nStateVars=4+nsoillayers;

  _aStateVarType =new sv_type [_nStateVars];
  _aStateVarLayer=new int     [_nStateVars];
  _nSoilVars     =nsoillayers;
  _nAquiferLayers=0;  //initially 
  _nSnowLayers   =0;  //initially    

  _aStateVarType[0]=SURFACE_WATER; _aStateVarLayer[0]=DOESNT_EXIST; _aStateVarIndices[(int)(SURFACE_WATER)][0]=0;
  _aStateVarType[1]=ATMOSPHERE;    _aStateVarLayer[1]=DOESNT_EXIST; _aStateVarIndices[(int)(ATMOSPHERE   )][0]=1;
  _aStateVarType[2]=ATMOS_PRECIP;  _aStateVarLayer[2]=DOESNT_EXIST; _aStateVarIndices[(int)(ATMOS_PRECIP )][0]=2;
  _aStateVarType[3]=PONDED_WATER;  _aStateVarLayer[3]=DOESNT_EXIST; _aStateVarIndices[(int)(PONDED_WATER )][0]=3;

  if (SM==SOIL_LUMPED){
    _aStateVarType[4]=LUMPED_LANDFORM;  _aStateVarLayer[4]=DOESNT_EXIST; _aStateVarIndices[(int)(LUMPED_LANDFORM )][0]=4;
  }
  else
  {
    int count=0;
    for (i=4;i<4+_nSoilVars;i++)
    {
      _aStateVarType [i]=SOIL;
      _aStateVarLayer[i]=count;
      _aStateVarIndices[(int)(SOIL)][count]=i;
      count++;
    }
  }
  _lake_sv=0; //by default, rain on lake goes direct to surface storage [0]
  

  CHydroProcessABC::SetModel(this);

  _aGaugeWeights =NULL;	//Initialized in Initialize

  _aCumulativeBal=NULL;
  _aFlowBal      =NULL;
  _CumulInput    =0.0;
  _CumulOutput   =0.0;
  _CumEnergyGain =0.0;
  _CumEnergyLoss =0.0;
  _initWater     =0.0;

  _UTM_zone=-1;

  _nOutputTimes=0;   _aOutputTimes=NULL;
  _currOutputTimeInd=0;
  _pOutputGroup=NULL;

  _aShouldApplyProcess=NULL; //Initialized in Initialize

  _pTransModel=new CTransportModel(this);
}

/////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CModel::~CModel()
{
  if (DESTRUCTOR_DEBUG){cout<<"DELETING MODEL"<<endl;}
  int c,g,i,j,k,kk,p;

  CloseOutputStreams();

  for (p=0;p<_nSubBasins;    p++){delete _pSubBasins    [p];} delete [] _pSubBasins;    _pSubBasins=NULL;
  for (k=0;k<_nHydroUnits;   k++){delete _pHydroUnits   [k];} delete [] _pHydroUnits;   _pHydroUnits=NULL;
  for (g=0;g<_nGauges;       g++){delete _pGauges       [g];} delete [] _pGauges;       _pGauges=NULL;
  for (j=0;j<_nProcesses;    j++){delete _pProcesses    [j];} delete [] _pProcesses;    _pProcesses=NULL;
  for (c=0;c<_nCustomOutputs;c++){delete _pCustomOutputs[c];} delete [] _pCustomOutputs;_pCustomOutputs=NULL;
  for (i=0;i<_nObservedTS;   i++){delete _pObservedTS   [i];} delete [] _pObservedTS;   _pObservedTS=NULL;
  if (_pModeledTS != NULL){
    for (i = 0; i < _nObservedTS; i++){ delete _pModeledTS[i]; } delete[] _pModeledTS;    _pModeledTS = NULL;
  }
	for (i=0;i<_nObsWeightTS;  i++){delete _pObsWeightTS  [i];} delete [] _pObsWeightTS;  _pObsWeightTS=NULL;
  for (j=0;j<_nDiagnostics;  j++){delete _pDiagnostics  [j];} delete [] _pDiagnostics;  _pDiagnostics=NULL;

  if (_aCumulativeBal!=NULL){
    for (k=0;k<_nHydroUnits;   k++){delete [] _aCumulativeBal[k];} delete [] _aCumulativeBal; _aCumulativeBal=NULL;
  }
  if (_aFlowBal!=NULL){
    for (k=0;k<_nHydroUnits;   k++){delete [] _aFlowBal[k];      } delete [] _aFlowBal;       _aFlowBal=NULL;
  }
  if (_aGaugeWeights!=NULL){
    for (k=0;k<_nHydroUnits;   k++){delete [] _aGaugeWeights[k]; } delete [] _aGaugeWeights;  _aGaugeWeights=NULL;
  }
  if (_aShouldApplyProcess!=NULL){
    for (k=0;k<_nProcesses;   k++){delete [] _aShouldApplyProcess[k]; } delete [] _aShouldApplyProcess;  _aShouldApplyProcess=NULL;
  }
  for (kk=0;kk<_nHRUGroups;kk++){delete _pHRUGroups[kk]; } delete [] _pHRUGroups;   _pHRUGroups  =NULL;
  for (j=0;j<_nTransParams;j++) {delete _pTransParams[j];} delete [] _pTransParams; _pTransParams=NULL;

  delete [] _aStateVarType;  _aStateVarType=NULL;
  delete [] _aStateVarLayer; _aStateVarLayer=NULL;
  delete [] _aSubBasinOrder; _aSubBasinOrder=NULL;
  delete [] _aOrderedSBind;  _aOrderedSBind=NULL;
  delete [] _aDownstreamInds;_aDownstreamInds=NULL;
  delete [] _aOutputTimes;   _aOutputTimes=NULL;
  delete [] _aObsIndex;      _aObsIndex=NULL;
  
  
  CSoilClass::      DestroyAllSoilClasses();
  CVegetationClass::DestroyAllVegClasses();
  CLandUseClass::   DestroyAllLUClasses();
  CTerrainClass::   DestroyAllTerrainClasses();
  CSoilProfile::    DestroyAllSoilProfiles();
  CAquiferStack::   DestroyAllAqStacks();
  CChannelXSect::   DestroyAllChannelXSections();

  delete _pTransModel;
}
/*****************************************************************
   Accessors
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Returns number of sub basins in model
///
/// \return Integer number of sub basins
//
int CModel::GetNumSubBasins   () const{return _nSubBasins;}
  
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
/// \brief Returns number of aquifer layers
///
/// \return Integer number of aquifer layers
//
int CModel::GetNumAquiferLayers  () const{return _nAquiferLayers;}

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
/// \param g [in] Process index
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
CHydroUnit *CModel::GetHRUByID(const int HRUID) const
{
  for (int k=0;k<_nHydroUnits;k++){
    if (HRUID==_pHydroUnits[k]->GetID()){return _pHydroUnits[k];}
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
/// \param SBID [in] Integer sub basin ID
/// \return pointer to Sub basin object corresponding to passed ID, if ID is valid
//
CSubBasin  *CModel::GetSubBasinByID(const long SBID) const
{ //could be quite slow...
  if (SBID<0){return NULL;}
  for (int p=0;p<_nSubBasins;p++){
    if (_pSubBasins[p]->GetID()==SBID){return _pSubBasins[p];}
  }
  return NULL;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns sub basin index corresponding to passed subbasin ID
///
/// \param SBID [in] Integer subbasin ID
/// \return Sub basin index corresponding to passed ID, if ID is valid
//
int         CModel::GetSubBasinIndex(const long SBID) const
{  //could be quite slow...
  if (SBID<0){return DOESNT_EXIST;}
  for (int p=0;p<_nSubBasins;p++){
    if (_pSubBasins[p]->GetID()==SBID){return p;}
  }
  return INDEX_NOT_FOUND;
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
int  CModel::GetGaugeIndexFromName (const string name) const{
  for (int g=0;g<_nGauges;g++){
    if (name==_pGauges[g]->GetName()){return g;}
  }
  return DOESNT_EXIST;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns current mass/energy flux (mm/d, MJ/m2/d, mg/m2/d) between two storage compartments iFrom and iTo
/// \details required for advective transport processes
///
/// \param k [in] HRU index
/// \param iFrom_test [in] index of "from" storage state variable
/// \param iTo_test [in] index of "to" storage state variable
/// \param &Options [in] Global model options information
/// \todo [optimize]
//
double CModel::GetFlux(const int k, const int iFrom_test, const int iTo_test, const optStruct &Options) const
{
  int q,j,js(0);
  double flow=0;
  const int *iFrom,*iTo;
  int nConnections;

  ExitGracefullyIf((k<0) || (k>=_nHydroUnits),"CModel::GetFlux: bad HRU index",RUNTIME_ERR);

  //Not Optimized
  //much more effective to store aFlowBal as aFlowRate[k][iFrom][iTo]?
  if (iFrom_test==iTo_test){return 0.0;}
  for (j=0;j<_nProcesses;j++)
  {
    iTo  =_pProcesses[j]->GetToIndices();
    iFrom=_pProcesses[j]->GetFromIndices();
    nConnections=_pProcesses[j]->GetNumConnections();
    for (q=0;q<nConnections;q++)
    {
      if ((iTo_test==iTo[q]) && (iFrom_test==iFrom[q]))
      {
        flow+=_aFlowBal[k][js]/Options.timestep;
      }
      if ((iTo_test==iFrom[q]) && (iFrom_test==iTo[q]))
      {
        flow-=_aFlowBal[k][js]/Options.timestep;
      }
      js++;
    }//end for q=0 to nConnections
  }//end for j=0 to nProcesses

  return flow;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns current mass/energy flux (mm/d, MJ/m2/d, mg/m2/d) between two storage compartments iFrom and iTo
/// \details required for advective transport processes
///
/// \param k [in] HRU index
/// \param js [in] index of process connection (i.e., j*)
/// \param &Options [in] Global model options information
/// \todo [optimize]
//
double CModel::GetFlux(const int k, const int js, const optStruct &Options) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((k<0) || (k>=_nHydroUnits),"CModel::GetFlux: bad HRU index",RUNTIME_ERR);
  ExitGracefullyIf((js<0) || (js>=_nTotalConnections),"CModel::GetFlux: bad connection index",RUNTIME_ERR);
#endif
  return _aFlowBal[k][js]/Options.timestep;
}

/*****************************************************************  
   Watershed Diagnostic Functions
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Returns area-weighted average total precipitation rate at all HRUs [mm/d]
/// 
/// \return Area-weighted average of total precipitation rate [mm/d] over all HRUs
//
double CModel::GetAveragePrecip() const
{
  double sum(0);
  for (int k=0;k<_nHydroUnits;k++)
  {
    sum+=_pHydroUnits[k]->GetForcingFunctions()->precip*_pHydroUnits[k]->GetArea();
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
    f=_pHydroUnits[k]->GetForcingFunctions();
    sum+=(f->precip*f->snow_frac)*_pHydroUnits[k]->GetArea();
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
    pF_hru =_pHydroUnits[k]->GetForcingFunctions();
    area_wt=_pHydroUnits[k]->GetArea()/_WatershedArea;

    Fave.precip					+=area_wt*pF_hru->precip;
    Fave.precip_daily_ave+= area_wt*pF_hru->precip_daily_ave;
    Fave.precip_5day		+=area_wt*pF_hru->precip_5day;
    Fave.snow_frac			+=area_wt*pF_hru->snow_frac;   

    Fave.temp_ave				+=area_wt*pF_hru->temp_ave;			
    Fave.temp_daily_min	+=area_wt*pF_hru->temp_daily_min;		
    Fave.temp_daily_max	+=area_wt*pF_hru->temp_daily_max;			
    Fave.temp_daily_ave	+=area_wt*pF_hru->temp_daily_ave;
    Fave.temp_month_max	+=area_wt*pF_hru->temp_month_max;
    Fave.temp_month_min	+=area_wt*pF_hru->temp_month_min;
    Fave.temp_month_ave +=area_wt*pF_hru->temp_month_ave;

    Fave.temp_ave_unc  	+=area_wt*pF_hru->temp_ave_unc;		
    Fave.temp_min_unc	  +=area_wt*pF_hru->temp_min_unc;			
    Fave.temp_max_unc	  +=area_wt*pF_hru->temp_max_unc;

    Fave.air_dens				+=area_wt*pF_hru->air_dens;      
    Fave.air_pres				+=area_wt*pF_hru->air_pres;      
    Fave.rel_humidity	  +=area_wt*pF_hru->rel_humidity;      
    
    Fave.cloud_cover    +=area_wt*pF_hru->cloud_cover;
    Fave.ET_radia       +=area_wt*pF_hru->ET_radia;    
    Fave.SW_radia				+=area_wt*pF_hru->SW_radia;   
    Fave.SW_radia_unc   +=area_wt*pF_hru->SW_radia_unc;
    Fave.SW_radia_net   +=area_wt*pF_hru->SW_radia_net;
    Fave.LW_radia				+=area_wt*pF_hru->LW_radia;      
    Fave.day_length			+=area_wt*pF_hru->day_length;    
    Fave.day_angle      +=area_wt*pF_hru->day_angle;   //not really necc. 

    Fave.wind_vel				+=area_wt*pF_hru->wind_vel;      

    Fave.PET						+=area_wt*pF_hru->PET;           
    Fave.OW_PET					+=area_wt*pF_hru->OW_PET; 
    Fave.PET_month_ave  +=area_wt*pF_hru->PET_month_ave; 

    Fave.potential_melt +=area_wt*pF_hru->potential_melt;

    Fave.subdaily_corr  +=area_wt*pF_hru->subdaily_corr; 
  }
  return Fave;

  /*LATER?:
  double sum(0);
  int f=GetForcingIndex(force);
  for (int k=0;k<_nHydroUnits;k++)
  {
    sum+=_pHydroUnits[k]->GetForcingFunction(f)*_pHydroUnits[k]->GetArea();
  }
    sum/=watershed_area;*/
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
    sum+=(_pHydroUnits[k]->GetStateVarValue(i)*_pHydroUnits[k]->GetArea()); 
  }
  return sum/_WatershedArea;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns area-weighted average of specified forcing function over watershed
/// 
/// \param &forcing_string [in] string identifier of forcing function to assess
/// \return Area-weighted average of specified forcing function
//
double CModel::GetAvgForcing (const string &forcing_string) const
{
  //Area-weighted average
  double sum=0.0;
  for (int k=0;k<_nHydroUnits;k++)
  {
    sum    +=_pHydroUnits[k]->GetForcing(forcing_string)*_pHydroUnits[k]->GetArea();
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
    sum+=_pSubBasins[p]->GetChannelStorage();		 //[m3] 
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
/// \brief Adds observed time series to model
/// 
/// \param *pTS [in] (valid) pointer to observed time series to be added to model
//
void              CModel::AddObservedTimeSeries(CTimeSeriesABC *pTS)
{
  if (!DynArrayAppend((void**&)(_pObservedTS),(void*)(pTS),_nObservedTS)){
   ExitGracefully("CModel::AddObservedTimeSeries: adding NULL observation time series",BAD_DATA);} 
}
//////////////////////////////////////////////////////////////////
/// \brief Adds observation weighting time series to model
/// 
/// \param *pTS [in] (valid) pointer to observation weights time series to be added to model
//
void              CModel::AddObservedWeightsTS(CTimeSeriesABC   *pTS)
{
  if (!DynArrayAppend((void**&)(_pObsWeightTS),(void*)(pTS),_nObsWeightTS)){
   ExitGracefully("CModel::AddObservedWeightsTS: adding NULL observed weights time series",BAD_DATA);} 
}
//////////////////////////////////////////////////////////////////
/// \brief Adds diagnostic to model
/// 
/// \param *pDiag [in] (valid) pointer to diagnostic to be added to model
//
void              CModel::AddDiagnostic        (CDiagnostic       *pDiag)
{
  if (!DynArrayAppend((void**&)(_pDiagnostics),(void*)(pDiag),_nDiagnostics)){
   ExitGracefully("CModel::AddDiagnostic: adding NULL diagnostic",BAD_DATA);} 
}

//////////////////////////////////////////////////////////////////
/// \brief Adds model output time to model
/// 
/// \param timestamp [in] ISO-format timestamp string YYYY-MM-DD hh:mm:ss
/// \note must be called after model starttime & duration are set
//
void  CModel::AddModelOutputTime   (const time_struct &tt_out, const optStruct &Options)
{
  double t_loc;//local time

  t_loc=tt_out.julian_day-Options.julian_start_day;
  int leap;
  int yr=Options.julian_start_year;
  while (yr<tt_out.year)
  {
    leap=0;
    if (IsLeapYear(yr)){leap=1;}
    t_loc+=(365+leap);
    yr++;
  }
  ExitGracefullyIf(t_loc<0,
    "AddModelOutputTime: Cannot have model output time prior to start of simulation",BAD_DATA_WARN);
  if (t_loc>Options.duration){
    WriteWarning("AddModelOutputTime: model output time specified after end of simulation. It will be ignored",Options.noisy);
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
  ExitGracefullyIf(pHydroProc==NULL                 ,"CModel AddProcess::NULL process"  ,BAD_DATA);
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
      ExitGracefullyIf((int)(aSV[ii])>MAX_STATE_VARS,
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
      ExitGracefullyIf(_nStateVars>MAX_STATE_VARS,
        "CModel::AddStateVariables: exceeded maximum number of state variables in model",RUNTIME_ERR);
    }  
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Adds aquifer state variables during model construction
/// \details Should be called once during parse of .rvi file when :AquiferLayers command is specified
///
/// \param nLayers [in] number of aquifer layers in model
//
void CModel::AddAquiferStateVars(const int nLayers)
{
  _nAquiferLayers=nLayers;
  sv_type *aSV =new sv_type [_nAquiferLayers];
  int     *aLev=new int       [_nAquiferLayers];
  for (int i=0;i<_nAquiferLayers;i++)
  {
    aSV [i]=GROUNDWATER;
    aLev[i]=i;
  }
  AddStateVariables(aSV,aLev,_nAquiferLayers);
  delete [] aSV;
  delete [] aLev;
}
/*****************************************************************
   Other Manipulator Functions
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Sets lake storage state variable index based on SV type and layer
///
/// \param *lak_sv [out] SV type
///	\param lev [in] Layer of interest
//
void CModel::SetLakeStorage   (const sv_type      lak_sv, const int lev)
{
  _lake_sv=GetStateVarIndex(lak_sv,lev);
  ExitGracefullyIf(_lake_sv==DOESNT_EXIST,
    "CModel::SetLakeStorage: non-existent state variable",BAD_DATA);
}

//////////////////////////////////////////////////////////////////
/// \brief Sets state variable defined by SV and lev to be an aggregated variable 
///
/// \param SV [in] SV type to be aggregated
///	\param lev [in] Layer of SV to be aggregated
/// \param group_name [in] HRU group over which aggregation takes place
//
void CModel::SetAggregatedVariable(const sv_type SV, const int lev, const string group_name)
{
  int i=GetStateVarIndex(SV,lev);
  ExitGracefullyIf(i==DOESNT_EXIST,
    "CModel::SetAggregatedVariable: non-existent state variable",BAD_DATA);
  for (int kk=0;kk<_nHRUGroups;kk++){
    if (!_pHRUGroups[kk]->GetName().compare(group_name)){_pHRUGroups[kk]->SetAsAggregator(i);}
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
/// \brief Sets HRU group for which storage output files are generated
///
/// \param *pOut [in] (assumed valid) pointer to HRU group
//

void CModel::SetOutputGroup(const CHRUGroup *pOut){_pOutputGroup=pOut;}

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
  ExitGracefullyIf(_nGauges<1,
    "CModel::Initialize: Must have at least one meteorological gauge station",BAD_DATA);
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

  // initialize water/energy balance arrays to zero
  //-------------------------------------------------------------- 
  _nTotalConnections=0;
  for (int j=0; j<_nProcesses;j++){
    if (_pProcesses[j]->GetProcessType()!=PRECIPITATION){_pProcesses[j]->Initialize();} //precip already initialized in ParseInput.cpp
    _nTotalConnections+=_pProcesses[j]->GetNumConnections();
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
  _CumulInput	  =_CumulOutput  =0.0;
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
  _initWater+=GetTotalRivuletStorage();

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
	for (k=0; k<_nHydroUnits;k++)
		_aShouldApplyProcess[j][k] = _pProcesses[j]->ShouldApply(_pHydroUnits[k]);
  }

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
}

//////////////////////////////////////////////////////////////////
/// \brief Initializes observation time series
/// \details Initializes _pObservedTS, _pModeledTS and _pObsWeightTS
///	Matches observation weights to the corresponding observations
///	Called from CModel::Initialize  
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
      else								  {_aSubBasinOrder[p]=_aSubBasinOrder[pTo]+1;}
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
///	Calculates all basin drainage areas and total watershed area
///	Called from CModel::Initialize AFTER initalize routing network 
///	(i.e., directed stream network has already been constructed)
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
  
  //Estimate initial lateral runoff / flows in each basin 
  //-----------------------------------------------------------------
  double runoff_est;
  for (p=0;p<_nSubBasins;p++)
  {    
    aSBArea[p]=_pSubBasins[p]->GetBasinArea();

    runoff_est=EstimateInitialRunoff(p,Options);//[mm/d]

    aSBQlat[p]=runoff_est/MM_PER_METER*(aSBArea[p]*M2_PER_KM2)/SEC_PER_DAY;//[m3/s]
    
    aSBQin [p]=_pSubBasins[p]->GetSpecifiedInflow(0.0);//initial conditions //[m3/s]

    //cout<<p<<": "<< aSBQin[p]<<" "<<runoff_est<<" "<<aSBQlat[p]<<" "<<aSBArea[p]<<endl;
  }

  //Propagate flows and drainage areas downstream
  //-----------------------------------------------------------------
  int pTo;
  for (int pp=0;pp<_nSubBasins;pp++)
  { //calculated in order from upstream to downstream
    p=this->GetOrderedSubBasinIndex(pp);
    
    pTo=_aDownstreamInds[p]; 
    if (pTo!=DOESNT_EXIST){
      aSBQin [pTo]+=aSBQlat[p]+aSBQin[p];
      aSBArea[pTo]+=aSBArea[p];//now aSBArea==drainage area
    }
    //cout<<p<<": "<< aSBQin[p]<<" "<<aSBQlat[p]<<" "<<aSBArea[p]<<endl;
    _pSubBasins[p]->Initialize(aSBQin[p],aSBQlat[p],aSBArea[p],Options);   
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
/// \brief Estimates initial runoff [mm/d]
/// \details Returns esimate of initial runoff rate within basin. called from cModel::InitializeBasinFlows
/// \remark Called after initialize routing orders
///
/// \param p [in] Sub basin index
/// \param &Options [in] Global model options information
/// \return Double value of initial runoff [mm/d]
//
double CModel::EstimateInitialRunoff(const int p, const optStruct &Options) const
{
  //Method #1: use average annual rainfall falling on subbasin area over entire simulation
  double precip=0;
  double gauge_precip;
  for (int g=0;g<_nGauges;g++)
  {
    gauge_precip=_pGauges[g]->GetForcingValue(F_PRECIP,0,Options.duration);//[mm/d]

    for (int k=0;k<_nHydroUnits;k++)
    {
      if (_pHydroUnits[k]->GetSubBasinIndex()==p){
        precip+=(_aGaugeWeights[k][g]*gauge_precip*_pHydroUnits[k]->GetArea());
      }
    } 
  }

  precip/=_pSubBasins[p]->GetBasinArea(); //area-weighted average//[mm/d]
  
  //steady-state inflow to basin p from rainfall in that basin 
  double RUNOFF_FACT=0.4;//may wish to change runoff factor - depends upon aridity of basin, avg annual PET, etc.
  
  return RUNOFF_FACT*precip;//[mm/d]
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
    case(INTERP_AVERAGE_ALL):			//---------------------------------------------	
      {
        for (k=0;k<_nHydroUnits;k++){
          for (g=0;g<_nGauges;g++){
            _aGaugeWeights[k][g]=1.0/(double)(_nGauges);
          }
        }
        break;
      }
    case(INTERP_INVERSE_DISTANCE):			//---------------------------------------------	
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
    case(INTERP_INVERSE_DISTANCE_ELEVATION):			//---------------------------------------------	
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
    case (INTERP_FROM_FILE):			//---------------------------------------------
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
    ExitGracefullyIf((fabs(sum-1.0)>REAL_SMALL) && (INTERP_FROM_FILE),
      "GenerateGaugeWeights: Bad weighting scheme- weights for each HRU must sum to 1",BAD_DATA);
    ExitGracefullyIf((fabs(sum-1.0)>REAL_SMALL) && !(INTERP_FROM_FILE),
      "GenerateGaugeWeights: Bad weighting scheme- weights for each HRU must sum to 1",RUNTIME_ERR);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Increments water/energy balance
/// \details add cumulative amount of water/energy for each process connection
/// \note  called from main solver
///
/// \param j_star [in] Integer index of hydrologic process connection
/// \param k [in] Integer index of HRU
/// \param moved [in] Double amount balance is to be incremented
//
void CModel::IncrementBalance(const int    j_star,
                              const int    k, 
                              const double moved){
  ExitGracefullyIf(j_star>_nTotalConnections,"CModel::IncrementBalance: bad index",RUNTIME_ERR);
  _aCumulativeBal[k][j_star]+=moved;
  _aFlowBal      [k][j_star]= moved;
}

//////////////////////////////////////////////////////////////////
/// \brief Increments cumulative mass & energy added to system (precipitation/ustream basin flows, etc.)
/// \details Increment cumulative precipitation based on average preciptation and corresponding timestep [mm]
///
/// \param &Options [in] Global model options information
//
void CModel::IncrementCumulInput(const optStruct &Options, const time_struct &tt)
{
  double area;
  _CumulInput+=GetAveragePrecip()*Options.timestep;

  for (int p=0;p<_nSubBasins;p++){
    area = _WatershedArea*M2_PER_KM2;
    _CumulInput+=_pSubBasins[p]->GetIntegratedSpecInflow(tt.model_time,Options.timestep)/area*MM_PER_METER;//converted to [mm] over  basin
  }

  _pTransModel->IncrementCumulInput(Options,tt);
}

//////////////////////////////////////////////////////////////////
/// \brief Increments cumulative outflow from system for mass balance diagnostics
/// \details Increment cumulative outflow according to timestep, flow, and area of basin
///
/// \param &Options [in] Global model options information
//
void CModel::IncrementCumOutflow(const optStruct &Options)
{
  double area;
  for (int p=0;p<_nSubBasins;p++)
  {
    if (_aSubBasinOrder[p]==0)//outlet does not drain into another subbasin
    { 
      area=(_pSubBasins[p]->GetDrainageArea()*M2_PER_KM2);
      _CumulOutput+=_pSubBasins[p]->GetIntegratedOutflow(Options.timestep)/area*MM_PER_METER;//converted to [mm] over entire drainage basin
    }
  }
  _pTransModel->IncrementCumulOutput(Options);
}

//////////////////////////////////////////////////////////////////
/// \brief Updates values of user-specified transient parameters
///
/// \param &Options [in] Global model options information
/// \param &tt [in] Current time structure
//
void CModel::UpdateTransientParams(const optStruct   &Options,
                                   const time_struct &tt)
{
  int nn=(int)((tt.model_time+REAL_SMALL)/Options.timestep);//current timestep index
  for (int j=0;j<_nTransParams;j++)
  {
    class_type ctype=_pTransParams[j]->GetParameterClassType();
    string     pname=_pTransParams[j]->GetParameterName();
    string     cname=_pTransParams[j]->GetParameterClass();
    double     value=_pTransParams[j]->GetTimeSeries()->GetSampledValue(nn);
    
    if (ctype==CLASS_SOIL)
    {
      CSoilClass::StringToSoilClass(cname)->SetSoilProperty(pname,value);
    }
    else if (ctype==CLASS_VEGETATION)
    {
      CVegetationClass::StringToVegClass(cname)->SetVegetationProperty(pname,value);
    }
    else if (ctype==CLASS_TERRAIN)
    {
      CTerrainClass::StringToTerrainClass(cname)->SetTerrainProperty(pname,value);
    }
    else if (ctype==CLASS_LANDUSE)
    {
      CLandUseClass::StringToLUClass(cname)->SetSurfaceProperty(pname,value);
    }
    else if (ctype==CLASS_GLOBAL)
    {
      CGlobalParams::SetGlobalProperty(pname,value);
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
    _pHydroUnits[k]->RecalculateDerivedParams(Options,tt);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Updates values stored in modeled time series of observation data
///
/// \param &Options [in] Global model options information
/// \param &tt [in] Current time structure
//
void CModel::UpdateDiagnostics(const optStruct   &Options,
                               const time_struct &tt)
{
  if (_nDiagnostics==0){return;}


  int n=(int)(floor((tt.model_time+TIME_CORRECTION)/Options.timestep));//current timestep index
  double value, obsTime;
  int layer_ind;
  for (int i=0;i<_nObservedTS;i++)
  {
    string datatype=_pObservedTS[i]->GetName();
    sv_type svtyp=CStateVariable::StringToSVType (datatype,layer_ind,false);

    if (datatype=="HYDROGRAPH")
    {
      CSubBasin *pBasin;
      pBasin=GetSubBasinByID (s_to_l(_pObservedTS[i]->GetTag().c_str()));
      if ((Options.ave_hydrograph) && (tt.model_time!=0))
      {
         value=pBasin->GetIntegratedOutflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);
      }
      else
      {
         value=pBasin->GetOutflowRate();
      }
    }
    else if (svtyp!=DOESNT_EXIST)
    {
      CHydroUnit *pHRU;
      pHRU=GetHRUByID(s_to_i(_pObservedTS[i]->GetTag().c_str()));
      int sv_index=this->GetStateVarIndex(svtyp,layer_ind); 
      value=pHRU->GetStateVarValue(sv_index); 
    }
    else
    {
      if (tt.model_time ==0){
        WriteWarning("CModel::UpdateDiagnostics: invalid tag used for specifying Observation type",BAD_DATA);
      }
      value=0;
    }
    _pModeledTS[i]->SetValue (n,value);

		obsTime =_pObservedTS[i]->GetSampledTime(_aObsIndex[i]); // time of the next observation
    while((tt.model_time >= obsTime+_pObservedTS[i]->GetSampledInterval()) && 
          (_aObsIndex[i]<_pObservedTS[i]->GetNumSampledValues()))
		{
			value= _pModeledTS[i]->GetModelledValue(obsTime,_pObservedTS[i]->GetType());
			_pModeledTS[i]->SetSampledValue(_aObsIndex[i],value);
			_aObsIndex[i]++;
			obsTime =_pObservedTS[i]->GetSampledTime(_aObsIndex[i]);
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
/// \param j [in] Integer process indentifier
/// \param *state_var [in] Array of state variables for HRU
/// \param *pHRU [in] Pointer to HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Time structure
/// \param *iFrom [in] Array (size: nConnections)  of Indices of state variable losing mass or energy
/// \param *iTo [in] Array (size: nConnections)  of indices of state variable gaining mass or energy
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
  ExitGracefullyIf((j<0) && (j>=_nProcesses),"CModel ApplyProcess::improper index",BAD_DATA);
  
  nConnections=_pProcesses[j]->GetNumConnections();//total connections: nConnections+nCascades
	if (!_aShouldApplyProcess[j][pHRU->GetGlobalIndex()]){return false;}
  
  for (int q=0;q<nConnections;q++)
  {
    iFrom[q]=_pProcesses[j]->GetFromIndices()[q];
    iTo  [q]=_pProcesses[j]->GetToIndices  ()[q];
    rates_of_change[q]=0.0;
  }

  _pProcesses[j]->GetRatesOfChange(state_var,pHRU,Options,tt,rates_of_change);

  //special cascade handling (prior to applying constraints)
  //------------------------------------------------------------------------
  if (_pProcesses[j]->HasCascade())
  {
    int nCascades;
    static double max_state_var[MAX_STATE_VARS];
    nCascades=_pProcesses[j]->GetNumCascades();
    for (int i=0;i<_nStateVars;i++){//should only calculate for participating compartments
      max_state_var[i]=pHRU->GetStateVarMax(i,state_var,Options);
    }
    for (int q=0;q<nCascades;q++)
    {
      iFrom[nConnections-nCascades+q]=_pProcesses[j]->GetCascadeFromIndex();
      iTo  [nConnections-nCascades+q]=_pProcesses[j]->GetCascadeToIndices()[q];
    }    
    _pProcesses[j]->Cascade(rates_of_change,state_var,&max_state_var[0],Options.timestep);
  }

  _pProcesses[j]->ApplyConstraints(state_var,pHRU,Options,tt,rates_of_change);
	return true;
}
