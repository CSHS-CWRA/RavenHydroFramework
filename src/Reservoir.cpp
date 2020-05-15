/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2020 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Reservoir.h"

//////////////////////////////////////////////////////////////////
/// \brief Base Constructor for reservoir called by all other constructors
/// \param Name [in] Nickname for reservoir
/// \param SubID [in] subbasin ID
/// \param typ [in] reservoir type 
/// \details needed because versions of c++ prior to v11 don't necessaarily support delegating constructors
//
void CReservoir::BaseConstructor(const string Name,const long SubID)
{
  _name=Name;
  _SBID=SubID;

  _stage     =0.0;
  _stage_last=0.0;
  _min_stage =0.0;
  _max_stage =0.0; 
  _Qout      =0.0;
  _Qout_last =0.0;
  _MB_losses =0.0;
  _AET		   =0.0;
  _GW_seepage=0.0;

  _pHRU=NULL;
  _pExtractTS=NULL;
  _pWeirHeightTS=NULL;
  _pMaxStageTS=NULL;
  _pOverrideQ=NULL;
  _pMinStageTS=NULL;     
  _pMinStageFlowTS=NULL;  
  _pTargetStageTS=NULL;
  _pMaxQIncreaseTS=NULL;
  _pMaxQDecreaseTS=NULL;
  _pDroughtLineTS=NULL;
  _pQmaxTS=NULL;
  _pQminTS=NULL;
  _pQdownTS=NULL;
  _QdownRange=0.0;
  _pQdownSB=NULL;
  _nDemands=0;
  _aDemands=NULL;
  _demand_mult=1.0;

  _pDZTR=NULL;

  _minStageDominant=false;

  _crest_ht=0.0;
  _max_capacity=0.0;

  _aQ_back=NULL;
  _nDates=0;
  _aDates=NULL;

  _Np=0;
  _aStage  =NULL;
  _aQ      =NULL;
  _aQunder =NULL;
  _aArea   =NULL;
  _aVolume =NULL;

  _seepage_const=0;
  _local_GW_head=0.0;
}

//////////////////////////////////////////////////////////////////
/// \brief Base Constructor for reservoir 
/// \param Name [in] Nickname for reservoir
/// \param SubID [in] subbasin ID
/// \param typ [in] reservoir type 
//
CReservoir::CReservoir(const string Name,const long SubID)
{
  BaseConstructor(Name,SubID);
}

//////////////////////////////////////////////////////////////////
/// \brief Constructor for reservoir using power law rating curves
/// \param Name [in] Nickname for reservoir
/// \param SubID [in] subbasin ID
/// \param a_V [in] power law coefficient for volume rating curve
/// \param b_V [in] power law exponent for volume rating curve
//
CReservoir::CReservoir(const string Name, const long SubID, 
                       //min_stage
                       const double a_V, const double b_V,
                       const double a_Q, const double b_Q,
                       const double a_A,const double b_A)
{
   BaseConstructor(Name,SubID);

  _Np=30;
  _aStage =new double [_Np];
  _aQ     =new double [_Np];
  _aQunder=new double [_Np];
  _aArea  =new double [_Np];
  _aVolume=new double [_Np];
  ExitGracefullyIf(_aVolume==NULL,"CReservoir::constructor",OUT_OF_MEMORY);

  double ht;
  _min_stage=0.0;
  _max_stage=10; //reasonable default?
  for (int i=0;i<_Np;i++)
  {
    ht  =_min_stage+(_max_stage-_min_stage)*(double)(i)/(double)(_Np-1);
    _aStage [i]=ht;
    _aQ     [i]=a_Q*pow(ht,b_Q);
    _aArea  [i]=a_A*pow(ht,b_Q);
    _aVolume[i]=a_V*pow(ht,b_Q);
    _aQunder[i]=0.0;
  }
  _max_capacity=_aVolume[_Np-1];
}
//////////////////////////////////////////////////////////////////
/// \brief Constructor for reservoir using lookup table rating curves
/// \param Name [in] Nickname for reservoir
/// \param SubID [in] subbasin ID
/// \param a_ht[] [in] array of reservoir stages [size: nPoints]
/// \param a_Q[] [in] array of reservoir discharges [size: nPoints]
/// \param a_A[] [in] array of reservoir volumes [size: nPoints]
/// \param a_V[] [in] array of reservoir areas [size: nPoints]
//
CReservoir::CReservoir(const string Name, const long SubID, 
                       const double *a_ht,
                       const double *a_Q, const double *a_Qund,const double *a_A, const double *a_V,
                       const int     nPoints)
{
   BaseConstructor(Name,SubID);

  _min_stage =ALMOST_INF;
  _max_stage=-ALMOST_INF; 

  _Np=nPoints;
  ExitGracefullyIf(_Np<2,"CReservoir::constructor: must have more than 1 data point in stage relations",BAD_DATA_WARN);

  _aStage =new double [_Np];
  _aQ     =new double [_Np];
  _aQunder=new double [_Np];
  _aArea  =new double [_Np];
  _aVolume=new double [_Np];
  ExitGracefullyIf(_aVolume==NULL,"CReservoir::constructor (2)",OUT_OF_MEMORY);

  string warn;
  for (int i=0;i<_Np;i++)
  {
    _aStage [i]=a_ht[i];
    lowerswap(_min_stage,_aStage[i]);
    upperswap(_max_stage,_aStage[i]);
    _aQ     [i]=a_Q[i];
    _aArea  [i]=a_A[i];
    _aVolume[i]=a_V[i];
    if(a_Qund==NULL){ _aQunder[i]=0.0; }
    else            { _aQunder[i]=a_Qund[i];   }
    if ((i > 0) && ((_aStage[i]-_aStage[i-1])<0)){
      warn = "CReservoir::constructor: stage relations must be specified in order of increasing stage. [bad reservoir: " + _name + " "+to_string(SubID)+"]";
      ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    }
    if ((i > 0) && ((_aVolume[i] - _aVolume[i-1]) <= -REAL_SMALL)){
      warn = "CReservoir::constructor: volume-stage relationships must be monotonically increasing for all stages. [bad reservoir: " + _name + " "+to_string(SubID)+"]";
      ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    }
    if ((i > 0) && ((_aQ[i] - _aQ[i-1]) < -REAL_SMALL)){
      warn = "CReservoir::constructor: stage-discharge relationships must be increasing or flat for all stages. [bad reservoir: " + _name + " "+to_string(SubID)+ "]";
      ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    }
    if ((i > 0) && ((_aQunder[i] - _aQunder[i-1]) < -REAL_SMALL)){
      warn = "CReservoir::constructor: stage-discharge (underflow) relationships must be increasing or flat for all stages. [bad reservoir: " + _name + " "+to_string(SubID)+ "]";
      ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    }
  }
  _max_capacity=_aVolume[_Np-1];
}

//////////////////////////////////////////////////////////////////
/// \brief Constructor for reservoir using VARYING lookup table rating curves
/// \param Name [in] Nickname for reservoir
/// \param SubID [in] subbasin ID
/// \param nDates [in] # of flow rating curves
/// \param aDates [in] array of julian days (integer <366) at which rating curve changes
/// \param a_ht[] [in] array of reservoir stages [size: nPoints]
/// \param a_QQ[][] [in] 2-D array of reservoir discharges [size: nDates * nPoints]
/// \param a_A[] [in] array of reservoir volumes [size: nPoints]
/// \param a_V[] [in] array of reservoir areas [size: nPoints]
//
CReservoir::CReservoir(const string Name, const long SubID, 
                       const int my_nDates, const int *my_aDates, 
                       const double *a_ht,
                       double      **a_QQ, 
                       const double *a_Qund,
                       const double *a_A, 
                       const double *a_V,
                       const int     nPoints)
{
  BaseConstructor(Name,SubID);

  _nDates=my_nDates;
  _aDates = new int[_nDates];
  for (int v = 0; v<_nDates; v++){
    _aDates[v] = my_aDates[v]-1; //Julian Days in Raven from 0 to 365, not 1 to 366
  }

  _min_stage =ALMOST_INF;
  _max_stage=-ALMOST_INF; //reasonable default?

  _Np=nPoints;
  ExitGracefullyIf(_Np<2,"CReservoir::constructor: must have more than 1 data point in stage relations",BAD_DATA_WARN);

  _aStage =new double [_Np];
  _aQ     =new double [_Np];
  _aQunder=new double [_Np];
  _aArea  =new double [_Np];
  _aQ_back = new double *[_nDates];
  for (int v = 0; v<_nDates; v++){ _aQ_back[v] = new double[_Np]; }
  _aVolume=new double [_Np];
  ExitGracefullyIf(_aVolume==NULL,"CReservoir::constructor (2)",OUT_OF_MEMORY);

  string warn;
  for (int i=0;i<_Np;i++)
  {
    _aStage [i]=a_ht[i];
    lowerswap(_min_stage,_aStage[i]);
    upperswap(_max_stage,_aStage[i]);
    _aQ     [i]=a_QQ[0][i];
    _aQunder[i]=0.0;if(a_Qund!=NULL){ _aQunder[i]=a_Qund[i]; }
    _aArea  [i]=a_A[i];
    _aVolume[i]=a_V[i];
    for (int v = 0; v < _nDates; v++){
      _aQ_back[v][i] = a_QQ[v][i];
      if ((i > 0) && ((_aQ_back[v][i] - _aQ_back[v][i-1]) < -REAL_SMALL)){
        warn = "CReservoir::constructor: stage-discharge relationships must be increasing or flat for all stages. [bad varying reservoir: " + _name + " "+to_string(SubID)+ "]";
        ExitGracefully(warn.c_str(),BAD_DATA_WARN);
      }
    }
    if ((i > 0) && ((_aVolume[i] - _aVolume[i-1]) <= -REAL_SMALL)){
      warn = "CReservoir::constructor: volume-stage relationships must be monotonically increasing for all stages. [bad reservoir: " + _name + " "+to_string(SubID)+"]";
      ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    }
    if ((i > 0) && ((_aQ[i] - _aQ[i-1]) < -REAL_SMALL)){
      warn = "CReservoir::constructor: stage-discharge relationships must be increasing or flat for all stages. [bad reservoir: " + _name + " "+to_string(SubID)+ "]";
      ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    }
    if ((i > 0) && ((_aQunder[i] - _aQunder[i-1]) < -REAL_SMALL)){
      warn = "CReservoir::constructor: stage-discharge relationships (underflow) must be increasing or flat for all stages. [bad reservoir: " + _name + " "+to_string(SubID)+ "]";
      ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    }
  }
  _max_capacity=_aVolume[_Np-1];
}
//////////////////////////////////////////////////////////////////
/// \brief Constructor for Prismatic lake reservoir controlled by weir coefficient
/// \param Name [in] Nickname for reservoir
/// \param SubID [in] subbasin ID
/// \param weircoeff [in] weir coefficient, <1.0
/// \param crestw [in] width of crest, [m]
/// \param A  [in] area of reservoir [m2]
/// \param depth [in] maximum depth of prismatic reservoir
//
CReservoir::CReservoir(const string Name,
                       const long   SubID,
                       const double weircoeff, 
                       const double crestw,
                       const double crestht, 
                       const double A,
                       const double depth)
{
   BaseConstructor(Name,SubID);

  _crest_ht=crestht;
  _min_stage =-depth+_crest_ht;
  _max_stage=6.0+_crest_ht; //reasonable default?
  ExitGracefullyIf(depth<0,"CReservoir::Constructor (Lake): cannot have negative maximum lake depth",BAD_DATA_WARN);
  
  _Np=30;
  _aStage =new double [_Np];
  _aQ     =new double [_Np];
  _aQunder=new double [_Np];
  _aArea  =new double [_Np];
  _aVolume=new double [_Np];
  ExitGracefullyIf(_aVolume==NULL,"CReservoir::constructor (4)",OUT_OF_MEMORY);

  double dh;
  string warn;
  dh=(_max_stage-_min_stage)/(_Np-1);
  for (int i=0;i<_Np;i++) // - Edited by KL to reference crest height properly
  {
    _aStage [i]=_min_stage+i*dh;
    if((_min_stage+(i-1)*dh-_crest_ht)*(_min_stage+i*dh-_crest_ht)<0.0){_aStage[i]=_crest_ht;} //ensures zero point is included
    if(_aStage[i]<_crest_ht){_aQ[i]=0.0;}
    else{
      _aQ[i]=2.0/3.0*sqrt(2*GRAVITY)*crestw*pow((_aStage[i]-_crest_ht),1.5); //Overflow weir equation
    }
    _aQunder[i]=0.0;
    _aArea  [i]=A;
    _aVolume[i]=A*(_aStage[i]-_min_stage);
  }
  _max_capacity=_aVolume[_Np-1];
}

//////////////////////////////////////////////////////////////////
/// \brief Default destructor
//
CReservoir::~CReservoir()
{
  delete [] _aQ;      _aQ    =NULL;
  delete [] _aQunder; _aQunder=NULL; 
  delete [] _aArea;   _aArea =NULL;
  delete [] _aVolume; _aVolume=NULL;
  for (int v = 0; v<_nDates; v++){ delete[] _aQ_back[v]; } delete [] _aQ_back; _aQ_back=NULL;
  delete[] _aDates; _aDates=NULL;

  delete _pExtractTS;_pExtractTS=NULL;
  delete _pWeirHeightTS;_pWeirHeightTS=NULL;
  delete _pMaxStageTS;_pMaxStageTS=NULL;
  delete _pOverrideQ;_pOverrideQ=NULL;
  delete _pMinStageTS;_pMinStageTS=NULL;     
  delete _pMinStageFlowTS;_pMinStageFlowTS=NULL;  
  delete _pTargetStageTS;_pTargetStageTS=NULL;
  delete _pMaxQIncreaseTS;_pMaxQIncreaseTS=NULL; 
  delete _pMaxQDecreaseTS;_pMaxQDecreaseTS=NULL; 
  delete _pDroughtLineTS; _pDroughtLineTS=NULL;
  delete _pQminTS; _pQminTS=NULL;
  delete _pQmaxTS; _pQmaxTS=NULL;
  delete _pQdownTS; _pQdownTS=NULL;
}

//////////////////////////////////////////////////////////////////
/// \returns Subbasin ID
//
long  CReservoir::GetSubbasinID        () const{return _SBID;}

//////////////////////////////////////////////////////////////////
/// \returns reservoir storage [m3]
//
double  CReservoir::GetStorage           () const{return GetVolume(_stage);}

//////////////////////////////////////////////////////////////////
/// \returns current outflow rate [m3/s]
//
double  CReservoir::GetOutflowRate        () const{return _Qout;}

//////////////////////////////////////////////////////////////////
/// \returns current stage [m]
//
double  CReservoir::GetResStage           () const{return _stage;}

//////////////////////////////////////////////////////////////////
/// \returns previous outflow rate [m3/s]
//
double CReservoir::GetOldOutflowRate() const { return _Qout_last; }

//////////////////////////////////////////////////////////////////
/// \returns previous storage [m3]
//
double CReservoir::GetOldStorage() const{return GetVolume(_stage_last);}

//////////////////////////////////////////////////////////////////
/// \returns evaporative losses integrated over previous timestep [m3]
//
double  CReservoir::GetReservoirLosses(const double &tstep) const
{
  return _MB_losses;
}
//////////////////////////////////////////////////////////////////
/// \returns evaporative losses only integrated over previous timestep [m3]
//
double  CReservoir::GetReservoirEvapLosses(const double &tstep) const
{
	return _AET;
}
//////////////////////////////////////////////////////////////////
/// \returns seepage losses only integrated over previous timestep [m3]
//
double  CReservoir::GetReservoirGWLosses(const double &tstep) const
{
  return _GW_seepage;
}
//////////////////////////////////////////////////////////////////
/// \returns outflow integrated over timestep [m3/d]
//
double  CReservoir::GetIntegratedOutflow        (const double &tstep) const
{
  return 0.5*(_Qout+_Qout_last)*(tstep*SEC_PER_DAY); //integrated
}
//////////////////////////////////////////////////////////////////
/// \returns global HRU index (or DOESNT_EXIST) if no HRU is linked to reservoir
//
int CReservoir::GetHRUIndex() const
{
  if(_pHRU==NULL){return DOESNT_EXIST;}
  return _pHRU->GetGlobalIndex();
}
//////////////////////////////////////////////////////////////////
/// \brief initializes reservoir variables
/// \param Options [in] model options structure
//
void CReservoir::Initialize(const optStruct &Options)
{
  double model_start_day=Options.julian_start_day;
  int    model_start_yr =Options.julian_start_year;
  double model_duration =Options.duration;
  double timestep       =Options.timestep;
  if(_pExtractTS!=NULL)
  {
    _pExtractTS->Initialize(model_start_day,model_start_yr,model_duration,timestep,false,Options.calendar);
  }
  if(_pWeirHeightTS!=NULL)
  {
    _pWeirHeightTS->Initialize(model_start_day,model_start_yr,model_duration,timestep,false,Options.calendar);
  }
  if(_pMaxStageTS!=NULL)
  {
    _pMaxStageTS->Initialize(model_start_day,model_start_yr,model_duration,timestep,false,Options.calendar);
  }
  if(_pOverrideQ!=NULL)
  {
    _pOverrideQ->Initialize(model_start_day,model_start_yr,model_duration,timestep,false,Options.calendar);
  }
  if(_pMinStageTS!=NULL)
  {
    _pMinStageTS->Initialize(model_start_day,model_start_yr,model_duration,timestep,false,Options.calendar);
  }
  if(_pMinStageFlowTS!=NULL)
  {
    _pMinStageFlowTS->Initialize(model_start_day,model_start_yr,model_duration,timestep,false,Options.calendar);
  }
  if(_pTargetStageTS!=NULL)
  {
    _pTargetStageTS->Initialize(model_start_day,model_start_yr,model_duration,timestep,false,Options.calendar);
  }
  if(_pMaxQIncreaseTS!=NULL)
  {
    _pMaxQIncreaseTS->Initialize(model_start_day,model_start_yr,model_duration,timestep,false,Options.calendar);
  }
  if(_pMaxQDecreaseTS!=NULL)
  {
    _pMaxQDecreaseTS->Initialize(model_start_day,model_start_yr,model_duration,timestep,false,Options.calendar);
  }
  if(_pDroughtLineTS!=NULL)
  {
    _pDroughtLineTS->Initialize(model_start_day,model_start_yr,model_duration,timestep,false,Options.calendar);
  }
  if(_pQminTS!=NULL)
  {
    _pQminTS->Initialize(model_start_day,model_start_yr,model_duration,timestep,false,Options.calendar);
  }
  if(_pQmaxTS!=NULL)
  {
    _pQmaxTS->Initialize(model_start_day,model_start_yr,model_duration,timestep,false,Options.calendar);
  }
  if(_pQdownTS!=NULL)
  {
    _pQdownTS->Initialize(model_start_day,model_start_yr,model_duration,timestep,false,Options.calendar);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Adds extraction history
/// \param *pOutflow Outflow time series to be added
//
void    CReservoir::AddExtractionTimeSeries (CTimeSeries *pOutflow)
{
  ExitGracefullyIf(_pExtractTS!=NULL,
                   "CReservoir::AddExtractionTimeSeries: only one extraction hydrograph may be specified per reservoir",BAD_DATA_WARN);
  _pExtractTS=pOutflow;
}
//////////////////////////////////////////////////////////////////
/// \brief Adds weir height time series [m]
/// \param *pWH weir height time series to be added
//
void    CReservoir::AddWeirHeightTS (CTimeSeries *pWH)
{
  ExitGracefullyIf(_pWeirHeightTS!=NULL,
                   "CReservoir::AddWeirHeightTS: only one weir height time series may be specified per reservoir",BAD_DATA_WARN);
  _pWeirHeightTS=pWH;
}
//////////////////////////////////////////////////////////////////
/// \brief Adds maximum stage time series
/// \param *pMS maximum stage time series [m]
//
void    CReservoir::AddMaxStageTimeSeries(CTimeSeries *pMS)
{
  ExitGracefullyIf(_pMaxStageTS!=NULL,
                   "CReservoir::AddWeirHeightTS: only one weir height time series may be specified per reservoir",BAD_DATA_WARN);
  _pMaxStageTS=pMS;
}
//////////////////////////////////////////////////////////////////
/// \brief Adds override flow time series 
/// \param *pQ override flow time series [m3/s]
//
void    CReservoir::AddOverrideQTimeSeries(CTimeSeries *pQ){
  ExitGracefullyIf(_pOverrideQ!=NULL,
                   "CReservoir::AddOverrideQTimeSeries: only one overridden flow time series may be specified per reservoir",BAD_DATA_WARN);
  _pOverrideQ=pQ;
}
//////////////////////////////////////////////////////////////////
/// \brief Adds minimum stage time series
/// \param *pMS minimum stage time series [m]
//
void    CReservoir::AddMinStageTimeSeries(CTimeSeries *pMS){
  ExitGracefullyIf(_pMinStageTS!=NULL,
                   "CReservoir::AddMinStageTimeSeries: only one minimum stage time series may be specified per reservoir",BAD_DATA_WARN);
  _pMinStageTS=pMS;
}
//////////////////////////////////////////////////////////////////
/// \brief Adds minimum stage flow time series 
/// \param *pQ minimum stage flow time series [m3/s]
//
void    CReservoir::AddMinStageFlowTimeSeries(CTimeSeries *pQ){
  ExitGracefullyIf(_pMinStageFlowTS!=NULL,
                   "CReservoir::AddMinStageFlowTimeSeries: only one minimum stage flow time series may be specified per reservoir",BAD_DATA_WARN);
  _pMinStageFlowTS=pQ;
}
//////////////////////////////////////////////////////////////////
/// \brief Adds target stage time series 
/// \param *pTS target stage flow time series [m3/s]
//
void    CReservoir::AddTargetStageTimeSeries(CTimeSeries *pTS){
  ExitGracefullyIf(_pTargetStageTS!=NULL,
                   "CReservoir::AddTargetStageTimeSeries: only one target stage time series may be specified per reservoir",BAD_DATA_WARN);
  _pTargetStageTS=pTS;
}
//////////////////////////////////////////////////////////////////
/// \brief Adds drought line time series 
/// \param *pTS drought linee flow time series [m3/s]
//
void    CReservoir::AddDroughtLineTimeSeries(CTimeSeries *pTS){
  ExitGracefullyIf(_pDroughtLineTS!=NULL,
                   "CReservoir::AddDroughtLineTimeSeries: only one target stage time series may be specified per reservoir",BAD_DATA_WARN);
  _pDroughtLineTS=pTS;
}
//////////////////////////////////////////////////////////////////
/// \brief Adds maximum flow increase time series 
/// \param *pQ maximum flow increase time series [m3/s/d]
//
void    CReservoir::AddMaxQIncreaseTimeSeries(CTimeSeries *pQdelta){
  ExitGracefullyIf(_pMaxQIncreaseTS!=NULL,
                   "CReservoir::AddMaxQIncreaseTimeSeries: only one maximum flow increase time series may be specified per reservoir",BAD_DATA_WARN);
  _pMaxQIncreaseTS=pQdelta;
}
//////////////////////////////////////////////////////////////////
/// \brief Adds maximum flow decrease time series 
/// \param *pQ maximum flow decrease time series [m3/s/d]
//
void    CReservoir::AddMaxQDecreaseTimeSeries(CTimeSeries *pQdelta) {
  ExitGracefullyIf(_pMaxQDecreaseTS!=NULL,
    "CReservoir::AddMaxQDecreaseTimeSeries: only one maximum flow decrease time series may be specified per reservoir",BAD_DATA_WARN);
  _pMaxQDecreaseTS=pQdelta;
}
//////////////////////////////////////////////////////////////////
/// \brief Adds minimum flow  time series 
/// \param *pQmin minimum flow  time series [m3/s]
//
void    CReservoir::AddMinQTimeSeries(CTimeSeries *pQmin) {
  ExitGracefullyIf(_pQminTS!=NULL,
    "CReservoir::AddMinQTimeSeries: only one minimum flow time series may be specified per reservoir",BAD_DATA_WARN);
  _pQminTS=pQmin;
}
//////////////////////////////////////////////////////////////////
/// \brief Adds maximum flow  time series 
/// \param *pQmax maximum flow  time series [m3/s]
//
void    CReservoir::AddMaxQTimeSeries(CTimeSeries *pQmax) {
  ExitGracefullyIf(_pQmaxTS!=NULL,
    "CReservoir::AddMinQTimeSeries: only one minimum flow time series may be specified per reservoir",BAD_DATA_WARN);
  _pQmaxTS=pQmax;
}
//////////////////////////////////////////////////////////////////
/// \brief Adds downstream target flow  time series 
/// \param *pQ downstream target flow  time series [m3/s]
/// \param SBIDdown downstream target flow basin ID
/// \param Qrange downstream target flow range [m3/s]
//
void    CReservoir::AddDownstreamTargetQ(CTimeSeries *pQ,const CSubBasin *pSB, const double &Qrange) 
{
  ExitGracefullyIf(_pQdownTS!=NULL,
    "CReservoir::AddDownstreamTargetQ: only one downstream flow time series may be specified per reservoir",BAD_DATA_WARN);
  _pQdownTS=pQ;
  _pQdownSB=pSB;
  _QdownRange=Qrange;
}
//////////////////////////////////////////////////////////////////
/// \brief Adds reservoir downstream demand
/// \param SBID subbasin ID of demand location or AUTO_COMPUTE_LONG
/// \param pct percentage of flow demand to be satisfied by reservoir as fraction [0..1] or AUTO_COMPUTE
//
void  CReservoir::AddDownstreamDemand(const CSubBasin *pSB,const double pct, const int julian_start, const int julian_end) 
{
  down_demand *pDemand;
  pDemand=new down_demand;
  pDemand->pDownSB     =pSB;
  pDemand->percent     =pct;
  pDemand->julian_start=julian_start;
  pDemand->julian_end  =julian_end;
  if(!DynArrayAppend((void**&)(_aDemands),(void*)(pDemand),_nDemands)) {
    ExitGracefully("CReservoir::AddDownstreamDemand: adding NULL source",RUNTIME_ERR);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief links reservoir to HRU
/// \param *pHRUpointer HRU to link to
//
void  CReservoir::SetHRU(const CHydroUnit *pHRUpointer)
{
  _pHRU=pHRUpointer;
}
//////////////////////////////////////////////////////////////////
/// \brief sets max capacity
/// \param capacity, in m3
//
void  CReservoir::SetMaxCapacity(const double &max_cap)
{
  _max_capacity=max_cap;
}
//////////////////////////////////////////////////////////////////
/// \brief gets max capacity
/// \returns capacity, in m3
//
double  CReservoir::GetMaxCapacity() const
{
  return _max_capacity;
}
//////////////////////////////////////////////////////////////////
/// \brief gets demand multiplier
/// \returns reservoir demand multiplier
//
double CReservoir::GetDemandMultiplier() const
{
  return _demand_mult;
}
//////////////////////////////////////////////////////////////////
/// \brief gets current constraint name
/// \returns current constraint applied to estimate stage/flow
//
string CReservoir::GetCurrentConstraint() const {
  switch(_constraint)
  {
  case(RC_MAX_STAGE):         {return "MAX_STAGE"; break;}
  case(RC_MIN_STAGE):         {return "RC_MIN_STAGE"; break;}
  case(RC_NATURAL):           {return "RC_NATURAL"; break;}
  case(RC_TARGET):            {return "RC_TARGET"; break;}
  case(RC_DOWNSTREAM_FLOW):   {return "RC_DOWNSTREAM_FLOW"; break;}
  case(RC_MAX_FLOW_INCREASE): {return "RC_MAX_FLOW_INCREASE"; break;}
  case(RC_MIN_FLOW):          {return "RC_MIN_FLOW"; break;}
  case(RC_MAX_FLOW):          {return "RC_MAX_FLOW"; break;}
  case(RC_OVERRIDE_FLOW):     {return "RC_OVERRIDE_FLOW"; break;}
  case(RC_DRY_RESERVOIR):     {return "RC_DRY_RESERVOIR"; break;}
  default:                    {return ""; break;}
  };
}
//////////////////////////////////////////////////////////////////
/// \brief overrides volume stage curve for Lake-type reservoirs (if known)
/// \param a_ht[] - array of stage values
/// \param a_V[] array of volumes 
/// \param nPoints - number of points in array
//

void CReservoir::SetVolumeStageCurve(const double *a_ht,const double *a_V,const int nPoints) 
{
  for(int i=0;i<_Np;i++)
  {
    _aVolume[i]=Interpolate2(_aStage[i],a_ht,a_V,nPoints,false);
    if((i > 0) && ((_aVolume[i] - _aVolume[i-1]) <= -REAL_SMALL)) {
      string warn = "CReservoir::SetVolumeStageCurve: volume-stage relationships must be monotonically increasing for all stages. [bad reservoir: " + _name + " "+to_string(_SBID)+"]";
      ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief overrides area stage curve for Lake-type reservoirs (if known)
/// \param a_ht[] - array of stage values
/// \param a_A[] array of areas 
/// \param nPoints - number of points in array
//

void CReservoir::SetAreaStageCurve(const double *a_ht,const double *a_A,const int nPoints)
{
  for(int i=0;i<_Np;i++)
  {
    _aArea[i]=Interpolate2(_aStage[i],a_ht,a_A,nPoints,false);
    if((i > 0) && ((_aArea[i] - _aArea[i-1]) <= -REAL_SMALL)) {
      string warn = "CReservoir::SetAreaStageCurve: area-stage relationships must be monotonically increasing for all stages. [bad reservoir: " + _name + " "+to_string(_SBID)+"]";
      ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief sets minmum stage constraint to dominant  
//
void CReservoir::SetMinStageDominant() 
{ 
  _minStageDominant=true; 
}
//////////////////////////////////////////////////////////////////
/// \brief sets minmum stage constraint to dominant  
//
void CReservoir::SetDemandMultiplier(const double &value)
{
  _demand_mult=value;
}
//////////////////////////////////////////////////////////////////
/// \brief sets parameters for Dynamically zoned target release model  
//
void CReservoir::SetDZTRModel(const double Qmc,const double Smax,
                              const double Sci[12],const double Sni[12],const double Smi[12],
                              const double Qci[12],const double Qni[12],const double Qmi[12])
{
  _pDZTR=new DZTRmodel();
  _pDZTR->Qmc=Qmc;
  _pDZTR->Vmax=Smax;
  
  for(int i=0;i<12;i++)
  {
    _pDZTR->Vci[i]=Sci[i];     _pDZTR->Vni[i]=Sni[i];     _pDZTR->Vmi[i]=Smi[i];
    _pDZTR->Qci[i]=Qci[i];     _pDZTR->Qni[i]=Qni[i];     _pDZTR->Qmi[i]=Qmi[i];

    if((Sci[i]>Sni[i]) || (Sni[i]>Smi[i]) || (Smi[i]>Smax)) {
      WriteWarning("CReservoir::SetDZTRModel: storage ordering is off. Vci<Vni<Vmi<Vmax",false);
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief gets flows from DZTR model of Yassin et al., 2019
/// Yassin et al., Representation and improved parameterization of reservoir operation in hydrological and land-surface models
/// Hydrol. Earth Syst. Sci., 23, 3735–3764, 2019 https://doi.org/10.5194/hess-23-3735-2019
//
double CReservoir::GetDZTROutflow(const double &V, const double &Qin, const time_struct &tt, const optStruct &Options) const
{
  double Vci=InterpolateMo(_pDZTR->Vci,tt,Options);
  double Vni=InterpolateMo(_pDZTR->Vni,tt,Options);
  double Vmi=InterpolateMo(_pDZTR->Vmi,tt,Options);
  double Qci=InterpolateMo(_pDZTR->Qci,tt,Options);
  double Qni=InterpolateMo(_pDZTR->Qni,tt,Options);
  double Qmi=InterpolateMo(_pDZTR->Qmi,tt,Options);
  double Vmin=0.1*_pDZTR->Vmax;
  double Qmc=_pDZTR->Qmc;
  double tstep=Options.timestep*SEC_PER_DAY;
  
  if      (V<Vmin){return 0.0;}
  else if (V<Vci ){return min(Qci,(V-Vmin)/tstep); }
  else if (V<Vni ){return Qci+(Qni-Qci)*(V-Vci)/(Vni-Vci);}
  else if (V<Vmi ){return Qni+max((Qin-Qni),(Qmi-Qni))*(V-Vni)/(Vmi-Vni); }
  else            {return min(Qmc,max(Qmi,(V-Vmi)/tstep));}

}
//////////////////////////////////////////////////////////////////
/// \brief sets reservoir groundwater parameters
/// \param coeff - groundwater exchange coefficient [m3/s/m]
/// \param h_ref - regional groundwater head [same vertical coordinate system as stage, m] 
//
void CReservoir::SetGWParameters(const double &coeff,const double &h_ref)
{
  _local_GW_head=h_ref;
  _seepage_const=coeff;
}
//////////////////////////////////////////////////////////////////
/// \brief sets all discharges in stage-discharge curve to zero (for overriding with observations)
//
void  CReservoir::DisableOutflow()
{
  for (int i=0;i<_Np;i++){_aQ[i]=0.0;_aQunder[i]=0.0;}
}

//////////////////////////////////////////////////////////////////
/// \brief updates state variable "stage" at end of computational time step
/// \param new_stage [in] calculated stage at end of time step
//
void  CReservoir::UpdateStage(const double &new_stage,const double &res_outflow,const res_constraint &constr,const optStruct &Options,const time_struct &tt)
{
  _stage_last=_stage;
  _stage     =new_stage;
  _constraint=constr;
  _Qout_last =_Qout;
  _Qout      =res_outflow;
}
//////////////////////////////////////////////////////////////////
/// \brief updates current mass balance (called at end of time step)
/// \param tt [in] current time step information
/// \param tstep [in] time step, in days
//
void CReservoir::UpdateMassBalance(const time_struct &tt,const double &tstep)
{
  _MB_losses=0.0;
  if(_pHRU!=NULL) {
    double Evap=_pHRU->GetForcingFunctions()->OW_PET;//mm/d
    if(_pHRU->GetSurfaceProps()->lake_PET_corr>=0.0) {
      Evap*=_pHRU->GetSurfaceProps()->lake_PET_corr;
    }
    _AET      = Evap* 0.5*(GetArea(_stage)+GetArea(_stage_last)) / MM_PER_METER*tstep; //m3
    _MB_losses+=_AET;
  }
  if(_seepage_const>0) {
    _GW_seepage=_seepage_const*(0.5*(_stage+_stage_last)-_local_GW_head)*SEC_PER_DAY*tstep;
    _MB_losses+=_GW_seepage;
  }
  if(_pExtractTS!=NULL){
    int nn        =(int)((tt.model_time+TIME_CORRECTION)/tstep);//current timestep index
    _MB_losses+=_pExtractTS->GetSampledValue(nn)*SEC_PER_DAY*tstep;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief updates rating curves based upon the time
/// \notes can later support quite generic temporal changes to treatmetn of outflow-volume-area-stage relations
/// \param tt [in] current model time
/// \param Options [in] model options structure
//
void CReservoir::UpdateFlowRules(const time_struct &tt, const optStruct &Options)
{
  if (_nDates == 0){return;}
  int vv=_nDates-1;
  for (int v = 0; v < _nDates; v++){
    if (tt.julian_day >= _aDates[v]){vv=v; }
  }
  
  for (int i = 0; i < _Np; i++){
    _aQ[i] = _aQ_back[vv][i];
  }
  return;
}
//////////////////////////////////////////////////////////////////
/// \brief initialize stage,volume, area to specified initial inflow
/// \param initQ [in] initial inflow
/// \note - ignores extraction and PET ; won't work for initQ=0, which could be non-uniquely linked to stage
//
void  CReservoir::SetInitialFlow(const double &initQ,const double &initQlast,const double &t)
{
  if(initQ!=initQlast) {//reading from .rvc file
    _Qout=initQ;
    _Qout_last=initQlast;
    return;
  }
  const double RES_TOLERANCE=0.001; //[m]
  const int RES_MAXITER=20;
  double dh=0.0001;
  double h_guess=0.1;

  double weir_adj=0.0;
  if(_pWeirHeightTS!=NULL){ 
    weir_adj=_pWeirHeightTS->GetValue(t); 
  }

  for (int i=0;i<_Np;i++)
  {
    if (_aQ[i]>0.0){h_guess=_min_stage+(double)(i)/(double)(_Np)*(_max_stage-_min_stage);break;}//initialize to just topping over crest
  }
  int    iter=0;
  double change=0;
  double Q,dQdh;
  do //Newton's method with discrete approximation of dQ/dh
  {
    //if(_pDZTR==NULL) {
      Q    =GetWeirOutflow(h_guess,   weir_adj);//[m3/s]
      dQdh=(GetWeirOutflow(h_guess+dh,weir_adj)-Q)/dh;
    //}
    //else if(_pDZTR!=NULL) {
    //  Q =GetDZTROutflow(GetVolume(h_guess),initQ,tt,Options);
    //  dQdh=(GetDZTROutflow(GetVolume(h_guess+dh),initQ,tt,Options)-Q)/dh;
    //}

    change=-(Q-initQ)/dQdh;//[m]
    if (dh==0.0){change=1e-7;}
    h_guess+=change;

    //cout <<iter<<": "<<h_guess<<" "<<Q<<" "<<dQdh<<" "<<change<<endl;
    iter++;
  } while ((iter<RES_MAXITER) && (fabs(change)>RES_TOLERANCE));

  if (iter==RES_MAXITER){
    string warn="CReservoir::SetInitialFlow did not converge after "+to_string(RES_MAXITER)+"  iterations for basin "+to_string(_SBID);
    WriteWarning(warn,false);
  }

  _stage=h_guess;
  _Qout=GetWeirOutflow(_stage,weir_adj);

  //if(_pDZTR==NULL) {
  //  _Qout=GetWeirOutflow(_stage,weir_adj);
  //}
  //else if(_pDZTR!=NULL) {
  //  _Qout=GetDZTROutflow(GetVolume(_stage),initQ,tt,Options);
  //}

  _stage_last=_stage;
  _Qout_last =_Qout;
}
//////////////////////////////////////////////////////////////////
/// \brief sets initial stage
/// \param ht [in] reservoir stage elevation
//
void  CReservoir::SetInitialStage(const double &ht,const double &ht_last)
{
  _stage_last=ht_last;
  _stage     =ht;
}
//////////////////////////////////////////////////////////////////
/// \brief sets minimum stage
/// \param min_z [in] minimum elevation with respect to arbitrary datum
//
void  CReservoir::SetMinStage(const double &min_z)
{
  _min_stage=min_z;
}
//////////////////////////////////////////////////////////////////
/// \brief Routes water through reservoir
/// \param Qin_old [in] inflow at start of timestep
/// \param Qin_new [in] inflow at end of timestep
/// \param tstep [in] numerical timestep
/// \param res_ouflow [out] outflow at end of timestep
/// \returns estimate of new stage at end of timestep
//
double  CReservoir::RouteWater(const double &Qin_old, const double &Qin_new, const optStruct &Options, const time_struct &tt, double &res_outflow,res_constraint &constraint) const
{
  const double RES_TOLERANCE=0.0001; //[m]
  const int    RES_MAXITER  =100;
  
  double tstep      =Options.timestep;
  double stage_new  =0.0;

  double stage_limit=ALMOST_INF;
  double weir_adj   =0.0;
  double Qoverride  =RAV_BLANK_DATA;
  double min_stage  =-ALMOST_INF;
  double Qminstage  =0.0;
  double Qmin       =0.0;
  double Qmax       =ALMOST_INF;
  double Qtarget    =RAV_BLANK_DATA;
  double htarget    =RAV_BLANK_DATA;
  double Qdelta     =ALMOST_INF;
  double Qdelta_dec =ALMOST_INF;
  double Qshift     =RAV_BLANK_DATA;

  int nn=(int)((tt.model_time+TIME_CORRECTION)/tstep);//current timestep index

  if(_pWeirHeightTS!=NULL)  { weir_adj   =_pWeirHeightTS->  GetSampledValue(nn);}
  if(_pMaxStageTS!=NULL)    { stage_limit=_pMaxStageTS->    GetSampledValue(nn);}
  if(_pOverrideQ!=NULL)     { Qoverride  =_pOverrideQ->     GetSampledValue(nn);}
  if(_pMinStageTS!=NULL)    { min_stage  =_pMinStageTS->    GetSampledValue(nn);}    
  if(_pMinStageFlowTS!=NULL){ Qminstage  =_pMinStageFlowTS->GetSampledValue(nn);}  
  if(_pTargetStageTS!=NULL) { htarget    =_pTargetStageTS-> GetSampledValue(nn);}
  if(_pMaxQIncreaseTS!=NULL){ Qdelta     =_pMaxQIncreaseTS->GetSampledValue(nn);}
  if(_pMaxQDecreaseTS!=NULL){ Qdelta_dec =_pMaxQDecreaseTS->GetSampledValue(nn);}
  if(_pQminTS!=NULL)        { Qmin       =_pQminTS->        GetSampledValue(nn);}
  if(_pQmaxTS!=NULL)        { Qmax       =_pQmaxTS->        GetSampledValue(nn);}

  // Downstream flow targets 
  if(_pQdownTS!=NULL)       { 
    double Qdown_targ      =_pQdownTS->GetSampledValue(nn);
    double Qdown_act       =_pQdownSB->GetOutflowRate();
    //Qshift=max(min((Qdown_targ-Qdown_act)/(_QdownRange*0.5),1.0),-1.0)*(Qdown_targ-Qdown_act);

    //smoothly shift towards target - far away moves towards range edge, in zone moves towards target 
    double alpha=exp(-fabs(Qdown_targ-Qdown_act)/(_QdownRange*0.5));
    if (_QdownRange==0){alpha=0.0;}
    if(Qdown_act<Qdown_targ) {Qshift=0.8*((Qdown_targ-(1.0-alpha)*(_QdownRange*0.5))-Qdown_act);}
    else                     {Qshift=0.8*((Qdown_targ+(1.0-alpha)*(_QdownRange*0.5))-Qdown_act);}
  }

  // Downstream irrigation demand 
  for(int i=0;i<_nDemands;i++) {
    if(IsInDateRange(tt.julian_day,_aDemands[i]->julian_start,_aDemands[i]->julian_end)){
      Qmin+=(_aDemands[i]->pDownSB->GetIrrigationDemand(tt.model_time)*_aDemands[i]->percent);
      Qmin+= _aDemands[i]->pDownSB->GetEnviroMinFlow   (tt.model_time)*1.0; //assume 100% of environmental min flow must be met 
    }
  }


  // ======================================================================================
  // Mass balance on reservoir over time step:
  //
  // dV(h)/dt=(Q_in(n)+Q_in(n+1))/2-(Q_out(n)+Q_out(n+1)(h))/2-ET*(A(n)+A(n+1))/2 -> solve for h_(n+1)
  //   non-linear w.r.t., h : rewritten as f(h)-gamma=0 for Newton's method solution
  //
  // ======================================================================================
  double dh        =0.001; //[m]
  double V_old     =GetVolume(_stage);
  double A_old     =GetArea(_stage);
  double h_guess   =_stage;
  int    iter      =0;
  double change    =0;
  double lastchange=-1;
  double f,dfdh,out,out2;
  double ET      (0.0);               //[m/s]
  double ext_old (0.0),ext_new (0.0); //[m3/s]
  double seep_old(0.0),seep_new(0.0);
  double outflow_nat,stage_nat;       //[m],[m3/s]

  if(_pHRU!=NULL)
  {
    ET=_pHRU->GetForcingFunctions()->OW_PET/SEC_PER_DAY/MM_PER_METER; //average for timestep, in m/s
    if(_pHRU->GetSurfaceProps()->lake_PET_corr>=0.0) {
      ET*=_pHRU->GetSurfaceProps()->lake_PET_corr;
    }
  }
  if(_seepage_const>0) {
    seep_old=_seepage_const*(_stage-_local_GW_head); //m3/s
  }
  if(_pExtractTS!=NULL)
  {
    ext_old=_pExtractTS->GetSampledValue(nn);
    ext_new=_pExtractTS->GetSampledValue(nn); //steady rate over time step
  }

  double gamma=V_old+((Qin_old+Qin_new)-_Qout-ET*A_old-seep_old-(ext_old+ext_new))/2.0*(tstep*SEC_PER_DAY);//[m3]
  if(gamma<0)
  {//reservoir dried out; no solution available. (f is always >0, so gamma must be as well)
    string warn="CReservoir::RouteWater: basin "+to_string(_SBID)+ " dried out on " +tt.date_string;
    WriteWarning(warn,false);
    return _min_stage;
  }

  //double hg[RES_MAXITER],ff[RES_MAXITER];//retain for debugging
  double relax=1.0;
  do //Newton's method with discrete approximation of df/dh
  {
    if     (_pDZTR==NULL) {
      out =GetWeirOutflow(h_guess,   weir_adj);//[m3/s]
      out2=GetWeirOutflow(h_guess+dh,weir_adj);//[m3/s]
    }
    else if(_pDZTR!=NULL) {
      out =GetDZTROutflow(GetVolume(h_guess   ),Qin_old,tt,Options);
      out2=GetDZTROutflow(GetVolume(h_guess+dh),Qin_old,tt,Options);
    }
    out +=ET*GetArea(h_guess   )+_seepage_const*(h_guess   -_local_GW_head);//[m3/s]
    out2+=ET*GetArea(h_guess+dh)+_seepage_const*(h_guess+dh-_local_GW_head);//[m3/s]

    f   = (GetVolume(h_guess   )+out /2.0*(tstep*SEC_PER_DAY)); //[m3]
    dfdh=((GetVolume(h_guess+dh)+out2/2.0*(tstep*SEC_PER_DAY))-f)/dh; //[m3/m]

    //hg[iter]=relax*h_guess; ff[iter]=f-gamma;//retain for debugging

    change=-(f-gamma)/dfdh;//[m]
    if(dfdh==0) { change=1e-7; }

    if(iter>3) { relax *=0.98; }
    h_guess+=relax*change;
    lastchange=change;
    iter++;
  } while((iter<RES_MAXITER) && (fabs(change/relax)>RES_TOLERANCE));

  stage_new=h_guess;

  if(iter==RES_MAXITER) {
    string warn="CReservoir::RouteWater did not converge after "+to_string(RES_MAXITER)+"  iterations for basin "+to_string(_SBID)+" on "+tt.date_string;;
    WriteWarning(warn,false);
    /*for (int i = 0; i < RES_MAXITER; i++){
      string warn = to_string(hg[i]) + " " + to_string(ff[i])+ " "+to_string(gamma);WriteWarning(warn,false);
    }*/
  }

  //standard case - outflow determined through stage-discharge curve
  //---------------------------------------------------------------------------------------------
  if(_pDZTR==NULL) {
    res_outflow=GetWeirOutflow(stage_new,weir_adj);
    constraint =RC_NATURAL;
  }
  else if(_pDZTR!=NULL) {
    res_outflow=GetDZTROutflow(GetVolume(stage_new),Qin_old,tt,Options);
    constraint =RC_DZTR;
  }

  outflow_nat=res_outflow; //saved for special max stage constraint
  stage_nat  =stage_new;

  //special correction - minimum stage reached or target flow- flow overriden (but forced override takes priority)
  //---------------------------------------------------------------------------------------------
  double w=CGlobalParams::GetParams()->reservoir_relax;
  if(htarget!=RAV_BLANK_DATA) {
    double V_targ=GetVolume(htarget);
    double A_targ=GetArea(htarget);
    double seep_targ = _seepage_const*(htarget-_local_GW_head);
    Qtarget = -2 * (V_targ - V_old) / (tstep*SEC_PER_DAY) + (-_Qout + (Qin_old + Qin_new) - ET*(A_old + A_targ) - (seep_old+seep_targ) - (ext_old + ext_new));//[m3/s]
    Qtarget=max(Qtarget,Qminstage);
    if(Qtarget>_Qout) {
      Qtarget = w*Qtarget+(1-w)*_Qout; //softer move towards goal - helps with stage undershoot -you can always remove more...
    }
    constraint=RC_TARGET; //target stage
  }

  if((Qshift!=RAV_BLANK_DATA) && (Qoverride==RAV_BLANK_DATA)) {
    //or (maybe) Qtarget+=Qshift if Qtarget!=RAV_BLANK_DATA)
    Qtarget=res_outflow+Qshift;
    constraint=RC_DOWNSTREAM_FLOW; //Downstream flow correction
  }

  if(Qoverride==RAV_BLANK_DATA)
  {
    if(stage_new<min_stage)
    {
      Qoverride=Qminstage;
      constraint=RC_MIN_STAGE;
    }
    else if(Qtarget!=RAV_BLANK_DATA)
    {
      if((Qtarget-_Qout)/tstep>Qdelta) {
        Qtarget=_Qout+Qdelta*tstep; //maximum flow increase
        constraint=RC_MAX_FLOW_INCREASE; 
      }
      else if((Qtarget-_Qout)/tstep<-Qdelta_dec) {
        Qtarget=_Qout-Qdelta_dec*tstep; //maximum flow decrease
        constraint=RC_MAX_FLOW_DECREASE; //max flow decrease
      }
      Qoverride=(Qtarget+_Qout)/2.0;//converts from end of time step to average over timestep
    }
  }
        
  //special correction - flow overridden or minimum/maximum flow violated - minimum/maximum flow takes priority
  //---------------------------------------------------------------------------------------------
  if ((Qoverride!=RAV_BLANK_DATA) || (res_outflow<Qmin) || (res_outflow>Qmax))
  {         
    if(Qoverride!=RAV_BLANK_DATA) {
      if((constraint!=RC_MAX_FLOW_INCREASE) && 
         (constraint!=RC_MAX_FLOW_DECREASE) && 
         (constraint!=RC_MIN_STAGE)) { constraint=RC_OVERRIDE_FLOW; } //Specified override flow
      res_outflow=max(2*Qoverride-_Qout,Qminstage); //Qoverride is avg over dt, res_outflow is end of dt
    }

    if((constraint==RC_MIN_STAGE) && (_minStageDominant)) { 
      //In this case, min stage constraint overrides min flow constraint; nothing done 
    }
    else
    {
      if (res_outflow<Qmin){ //overwrites any other specified or target flow
        res_outflow=Qmin;
        constraint=RC_MIN_FLOW; //minimum flow
      }
      
      if(res_outflow>Qmax) { //overwrites any other specified or target flow
        res_outflow=Qmax;
        constraint=RC_MAX_FLOW; //maximum flow
      }

      double A_guess=A_old;
      double A_last,V_new;
      double seep_guess = seep_old;
      do {
        V_new= V_old+((Qin_old+Qin_new)-(_Qout+res_outflow)-ET*(A_old+A_guess)-(seep_old+seep_guess)-(ext_old+ext_new))/2.0*(tstep*SEC_PER_DAY); 
        if(V_new<0) {
          V_new=0;
          constraint=RC_DRY_RESERVOIR; //drying out reservoir
          res_outflow = -2 * (V_new - V_old) / (tstep*SEC_PER_DAY) + (-_Qout + (Qin_old + Qin_new) - ET*(A_old + 0.0) - (ext_old + ext_new));//[m3/s] //dry it out
        }
        stage_new=Interpolate2(V_new,_aVolume,_aStage,_Np,false);
        A_last     = A_guess;
        A_guess    = GetArea(stage_new);
        seep_guess = _seepage_const*(stage_new-_local_GW_head);
      } while(fabs(1.0-A_guess/A_last)>0.00001); //0.1% area error - done in one iter for constant area case
    }
  }


  //special correction : exceeded limiting stage; fix stage and re-calculate reservoir outflow
  //---------------------------------------------------------------------------------------------
  if(stage_new>stage_limit) {
    constraint=RC_MAX_STAGE; //max stage exceedance
    stage_new=stage_limit;
    double V_limit =GetVolume(stage_limit);
    double A_limit =GetArea(stage_limit);
    double seep_lim=_seepage_const*(stage_limit-_local_GW_head);
    res_outflow = -2 * (V_limit - V_old) / (tstep*SEC_PER_DAY) + (-_Qout + (Qin_old + Qin_new) - ET*(A_old + A_limit) - (seep_old+seep_lim) - (ext_old + ext_new));//[m3/s]
  }

  //other option - returns to stage-discharge curve once max stage exceeded
  if((Options.res_overflowmode==OVERFLOW_NATURAL) &&  (stage_nat>stage_limit)) {
    res_outflow=outflow_nat;
    stage_new  =stage_nat;
  }
  // ======================================================================================
  
  return stage_new;
}

//////////////////////////////////////////////////////////////////
/// \brief writes state variables to solution file
//
void CReservoir::WriteToSolutionFile (ofstream &OUT) const
{
  OUT<<"    :ResFlow, "<<_Qout<<","<<_Qout_last<<endl;
  OUT<<"    :ResStage, "<<_stage<<","<<_stage_last<<endl;
}
//////////////////////////////////////////////////////////////////
/// \brief interpolates value from rating curve
/// \param x [in] interpolation location
/// \param xmin [in] lowest point in x range
/// \param xmax [in] highest point in x range
/// \param y [in] array (size:N)  of values corresponding to N evenly spaced points in x
/// \param N size of array y
/// \returns y value corresponding to interpolation point
/// \note assumes regular spacing between min and max x value
//
double Interpolate(double x, double xmin, double xmax, double *y, int N, bool extrapbottom)
{
  //assumes N *evenly-spaced* points in x from xmin to xmax
  double dx=((xmax-xmin)/(N-1));
  if      (x<=xmin){
    if (extrapbottom){return y[0]+(y[1]-y[0])/dx*(x-xmin);}
    return y[0];
  }
  else if (x>=xmax){
    return y[N-1]+(y[N-1]-y[N-2])/dx*(x-xmax);//extrapolation-may wish to revisit
    //return y[N-1];
  }
  else
  {
    double val=(x-xmin)/dx;
    int i=(int)(floor(val));
    return y[i]+(y[i+1]-y[i])*(val-floor(val));
  }
}


//////////////////////////////////////////////////////////////////
/// \brief interpolates the volume from the volume-stage rating curve
/// \param ht [in] reservoir stage
/// \returns reservoir volume [m3] corresponding to stage ht
/// \note assumes regular spacing between min and max stage
//
double     CReservoir::GetVolume(const double &ht) const
{
  return Interpolate2(ht,_aStage,_aVolume,_Np,true);
}
//////////////////////////////////////////////////////////////////
/// \brief interpolates the surface area from the area-stage rating curve
/// \param ht [in] reservoir stage
/// \returns reservoir surface area [m2] corresponding to stage ht
/// \note assumes regular spacing between min and max stage
//
double     CReservoir::GetArea  (const double &ht) const
{
  return Interpolate2(ht,_aStage,_aArea,_Np,false);
}
//////////////////////////////////////////////////////////////////
/// \brief interpolates the discharge from the outflow-stage rating curve
/// \param ht [in] reservoir stage
/// \returns reservoir outflow [m3/s] corresponding to stage ht
/// \note assumes regular spacing between min and max stage
//
double     CReservoir::GetWeirOutflow(const double &ht, const double &adj) const
{
  double underflow=Interpolate2(ht,_aStage,_aQunder,_Np,false); //no adjustments
  return Interpolate2(ht-adj,_aStage,_aQ,_Np,false)+underflow;
}


