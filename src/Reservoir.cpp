/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2018 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Reservoir.h"

double Interpolate2(double x, double *xx, double *y, int N, bool extrapbottom);

//////////////////////////////////////////////////////////////////
/// \brief Base Constructor for reservoir called by all other constructors
/// \param Name [in] Nickname for reservoir
/// \param SubID [in] subbasin ID
/// \param typ [in] reservoir type 
/// \details needed because versions of c++ prior to v11 don't necessaarily support delegating constructors
//
void CReservoir::BaseConstructor(const string Name,const long SubID,const res_type typ)
{
  _name=Name;
  _SBID=SubID;
  _type=typ;

  _stage     =0.0;
  _stage_last=0.0;
  _min_stage =0.0;
  _max_stage =0.0; 
  _Qout      =0.0;
  _Qout_last =0.0;
  _MB_losses =0.0;
  _AET		   =0.0;

  _pHRU=NULL;
  _pExtractTS=NULL;
  _pWeirHeightTS=NULL;
  _pMaxStageTS=NULL;
  _pOverrideQ=NULL;
  _pMinStageTS=NULL;     
  _pMinStageFlowTS=NULL;  
  _pTargetStageTS=NULL;
  _pMaxQIncreaseTS=NULL; 
  _pDroughtLineTS=NULL;

  _crest_ht=0.0;

  _aQ_back=NULL;
  _nDates=0;
  _aDates=NULL;

  _Np=0;
  _aStage  =NULL;
  _aQ      =NULL;
  _aQunder =NULL;
  _aArea   =NULL;
  _aVolume =NULL;
}

//////////////////////////////////////////////////////////////////
/// \brief Base Constructor for reservoir 
/// \param Name [in] Nickname for reservoir
/// \param SubID [in] subbasin ID
/// \param typ [in] reservoir type 
//
CReservoir::CReservoir(const string Name,const long SubID,const res_type typ)
{
  BaseConstructor(Name,SubID,typ);
}

//////////////////////////////////////////////////////////////////
/// \brief Constructor for reservoir using power law rating curves
/// \param Name [in] Nickname for reservoir
/// \param SubID [in] subbasin ID
/// \param a_V [in] power law coefficient for volume rating curve
/// \param b_V [in] power law exponent for volume rating curve
//
CReservoir::CReservoir(const string Name, const long SubID, const res_type typ,
                       //min_stage
                       const double a_V, const double b_V,
                       const double a_Q, const double b_Q,
                       const double a_A,const double b_A)
{
   BaseConstructor(Name,SubID,typ);

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
CReservoir::CReservoir(const string Name, const long SubID, const res_type typ,
                       const double *a_ht,
                       const double *a_Q, const double *a_Qund,const double *a_A, const double *a_V,
                       const int     nPoints)
{
   BaseConstructor(Name,SubID,typ);

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
CReservoir::CReservoir(const string Name, const long SubID, const res_type typ,
                       const int my_nDates, const int *my_aDates, 
                       const double *a_ht,
                       double      **a_QQ, 
                       const double *a_Qund,
                       const double *a_A, 
                       const double *a_V,
                       const int     nPoints)
{
  BaseConstructor(Name,SubID,typ);

  _nDates=my_nDates;
  _aDates = new int[_nDates];
  for (int v = 0; v<_nDates; v++){
    _aDates[v] = my_aDates[v]-1; //Julian Days in Raven from 0 to 365, not 1 to 365
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
                       const res_type typ,
                       const double weircoeff, 
                       const double crestw,
                       const double crestht, 
                       const double A,
                       const double depth)
{
   BaseConstructor(Name,SubID,typ);

  _type=RESROUTE_STANDARD;

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
  delete _pDroughtLineTS; _pDroughtLineTS=NULL;
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
  if (_pExtractTS!=NULL)
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
  if(_pDroughtLineTS!=NULL)
  {
    _pDroughtLineTS->Initialize(model_start_day,model_start_yr,model_duration,timestep,false,Options.calendar);
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
/// \brief links reservoir to HRU
/// \param *pHRUpointer HRU to link to
//
void  CReservoir::SetHRU(const CHydroUnit *pHRUpointer)
{
  _pHRU=pHRUpointer;
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
void  CReservoir::UpdateStage(const double &new_stage,const double &res_outflow,const optStruct &Options,const time_struct &tt)
{
  _stage_last=_stage;
  _stage     =new_stage;

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
  if(_pHRU!=NULL){
    _AET      =_pHRU->GetSurfaceProps()->lake_PET_corr*_pHRU->GetForcingFunctions()->PET * 0.5*( GetArea(_stage)+GetArea(_stage_last)) / MM_PER_METER*tstep;
    _MB_losses+=_AET;
  }

  if(_pExtractTS!=NULL){
    int nn        =(int)((tt.model_time+REAL_SMALL)/tstep);//current timestep index
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
void  CReservoir::SetInitialFlow(const double &initQ,const double &t)
{
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
    Q   = GetOutflow(h_guess   ,weir_adj);
    dQdh=(GetOutflow(h_guess+dh,weir_adj)-Q)/dh;
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
  _Qout=GetOutflow(_stage,weir_adj);

  _stage_last=_stage;
  _Qout_last =_Qout;
}
//////////////////////////////////////////////////////////////////
/// \brief sets initial stage
/// \param ht [in] reservoir stage elevation
//
void  CReservoir::SetInitialStage(const double &ht)
{
  _stage_last=ht;
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
double  CReservoir::RouteWater(const double &Qin_old, const double &Qin_new, const double &tstep, const time_struct &tt, double &res_outflow) const
{
  double stage_new=0;
  const double RES_TOLERANCE=0.0001; //[m]
  const int RES_MAXITER=100;
  double stage_limit=ALMOST_INF;
  double weir_adj=0.0;
  double Qoverride=RAV_BLANK_DATA;
  double min_stage=-ALMOST_INF;
  double Qminstage=0.0;
  double Qtarget=RAV_BLANK_DATA;
  double htarget=RAV_BLANK_DATA;
  double Qdelta=ALMOST_INF;

  int nn=(int)((tt.model_time+REAL_SMALL)/tstep);//current timestep index

  if(_pWeirHeightTS!=NULL)  { weir_adj   =_pWeirHeightTS->  GetSampledValue(nn);}
  if(_pMaxStageTS!=NULL)    { stage_limit=_pMaxStageTS->    GetSampledValue(nn);}
  if(_pOverrideQ!=NULL)     { Qoverride  =_pOverrideQ->     GetSampledValue(nn);}
  if(_pMinStageTS!=NULL)    { min_stage  =_pMinStageTS->    GetSampledValue(nn);}    
  if(_pMinStageFlowTS!=NULL){ Qminstage  =_pMinStageFlowTS->GetSampledValue(nn);}  
  if(_pTargetStageTS!=NULL) { htarget    =_pTargetStageTS-> GetSampledValue(nn);}
  if(_pMaxQIncreaseTS!=NULL){ Qdelta     =_pMaxQIncreaseTS->GetSampledValue(nn);}

  //=============================================================================
  if (_type==RESROUTE_NONE)
  {
    //find stage such that Qout=Qin_new
    /// \todo [funct]: add RESROUTE_NONE functionality
    ExitGracefully("RESROUTE_NONE",STUB);
    stage_new=0.0;
    res_outflow=0.0;
  }
  //=============================================================================
  else if (_type==RESROUTE_STANDARD)
  {
    // Mass balance on reservoir over time step:
    //
    // dV(h)/dt=(Q_in(n)+Q_in(n+1))/2-(Q_out(n)+Q_out(n+1)(h))/2-ET*(A(n)+A(n+1))/2 -> solve for h_(n+1)
    //   non-linear w.r.t., h : rewritten as f(h)-gamma=0 for Newton's method solution
    //

    double dh=0.001; //[m]

    double V_old  =GetVolume(_stage);
    double A_old  =GetArea(_stage);
    double h_guess=_stage;
    int    iter   =0;
    double change =0;
    double lastchange=-1;
    double f,dfdh,out,out2;
    double ET(0.0);  //m/s
    double ext_old(0.0),ext_new(0.0); //m3/s
    

    if (_pHRU!=NULL)
    {
      ET=_pHRU->GetForcingFunctions()->OW_PET/SEC_PER_DAY/MM_PER_METER; //average for timestep, in m/s
      if(_pHRU->GetSurfaceProps()->lake_PET_corr>=0.0){
        ET*=_pHRU->GetSurfaceProps()->lake_PET_corr;
      }
    }
    if (_pExtractTS!=NULL)
    {
      ext_old=_pExtractTS->GetSampledValue(nn);
      ext_new=_pExtractTS->GetSampledValue(nn); //steady rate over time step
    }

    double gamma=V_old+((Qin_old+Qin_new)-_Qout-ET*A_old-(ext_old+ext_new))/2.0*(tstep*SEC_PER_DAY);//[m3]
    if (gamma<0)
    {//reservoir dried out; no solution available. (f is always >0, so gamma must be as well)
      string warn="CReservoir::RouteWater: basin "+to_string(_SBID)+ " dried out on " +tt.date_string;
      WriteWarning(warn,false);
      return _min_stage;
    }

    //double hg[RES_MAXITER],ff[RES_MAXITER];//retain for debugging
    double relax=1.0;
    do //Newton's method with discrete approximation of df/dh
    {
      out =GetOutflow(h_guess   ,weir_adj)+ET*GetArea(h_guess   );//[m3/s]
      out2=GetOutflow(h_guess+dh,weir_adj)+ET*GetArea(h_guess+dh);//[m3/s]

      f   = (GetVolume(h_guess   )+out /2.0*(tstep*SEC_PER_DAY)); //[m3]
      dfdh=((GetVolume(h_guess+dh)+out2/2.0*(tstep*SEC_PER_DAY))-f)/dh; //[m3/m]

      //hg[iter]=relax*h_guess; ff[iter]=f-gamma;//retain for debugging

      change=-(f-gamma)/dfdh;//[m]
      if (dfdh==0){change=1e-7;}

      if(iter>3){relax *=0.98;}
      h_guess+=relax*change;
      lastchange=change;
      iter++;
    } while ((iter<RES_MAXITER) && (fabs(change/relax)>RES_TOLERANCE));

    stage_new=h_guess;

    if (iter==RES_MAXITER){
      string warn="CReservoir::RouteWater did not converge after "+to_string(RES_MAXITER)+"  iterations for basin "+to_string(_SBID)+" on "+tt.date_string;;
      WriteWarning(warn,false);
      /*for (int i = 0; i < RES_MAXITER; i++){
        string warn = to_string(hg[i]) + " " + to_string(ff[i])+ " "+to_string(gamma);WriteWarning(warn,false);
        }*/
    }

    //standard case - outflow determined through stage-discharge curve
    res_outflow=GetOutflow(stage_new,weir_adj);

    //special correction - minimum stage reached or target flow- flow overriden (but forced override takes priority)
    if(htarget!=RAV_BLANK_DATA){
			double V_targ=GetVolume(htarget);
      double A_targ=GetArea(htarget);
      Qtarget = -2 * (V_targ - V_old) / (tstep*SEC_PER_DAY) + (-_Qout + (Qin_old + Qin_new) - ET*(A_old + A_targ) - (ext_old + ext_new));//[m3/s]
      Qtarget=(Qtarget+_Qout)/2.0; //should be over entire timestep
		  Qtarget=max(Qtarget,Qminstage);
    }

    if (Qoverride==RAV_BLANK_DATA)
    {
      if(stage_new<min_stage)
      {
        Qoverride=Qminstage;
      }
      else if(Qtarget!=RAV_BLANK_DATA)
      {
        if ((Qtarget-_Qout)/tstep>Qdelta){
          Qtarget=_Qout+Qdelta*tstep; //maximum flow change
        }
        Qoverride=Qtarget;
      }
    }
    
    //special correction - flow overridden
    if(Qoverride!=RAV_BLANK_DATA)
    {
      // \todo[funct] may wish to ensure that V_new is not negative
      res_outflow=2*Qoverride-_Qout;
      double A_guess=A_old;
      double A_last,V_new;
      do{
        V_new= V_old+((Qin_old+Qin_new)-(_Qout+Qoverride)-ET*(A_old+A_guess)-(ext_old+ext_new))/2.0*(tstep*SEC_PER_DAY); 
        stage_new=Interpolate2(V_new,_aVolume,_aStage,_Np,false);
        A_last=A_guess;
        A_guess=GetArea(stage_new);
      } while(fabs(1.0-A_guess/A_last)>0.00001); //0.1% area error - done in one iter for constant area case
    }
    
    //special correction : exceeded limiting stage; fix stage and re-calculate reservoir outflow
    if(stage_new>stage_limit){
      stage_new=stage_limit;
      double V_limit=GetVolume(stage_limit);
      double A_limit=GetArea  (stage_limit);
	    res_outflow = -2 * (V_limit - V_old) / (tstep*SEC_PER_DAY) + (-_Qout + (Qin_old + Qin_new) - ET*(A_old + A_limit) - (ext_old + ext_new));//[m3/s]
    }

  }
  //=============================================================================
  
  return stage_new;
}

//////////////////////////////////////////////////////////////////
/// \brief writes state variables to solution file
//
void CReservoir::WriteToSolutionFile (ofstream &OUT) const
{
  OUT<<"    :ResStage, "<<_stage<<","<<_stage_last<<endl;
}

//////////////////////////////////////////////////////////////////
/// \brief parses the :Reservoir Command
/// \note the first line (:Reservoir [name]) has already been read in
/// \param p [in] parser (likely points to .rvh file)
/// \param name [in] name of reservoir
/// \param Options [in] model options structure
//
CReservoir *CReservoir::Parse(CParser *p, string name, int &HRUID,  const optStruct &Options)
{

  /*Examples:

    :Reservoir ExampleReservoir
      :SubBasinID 23
      :HRUID 234
      :Type RESROUTE_STANDARD
      :VolumeStageRelation POWER_LAW
        0.1 2.3
      :EndVolumeStageRelation
      :OutflowStageRelation POWER_LAW
        0.1 2.3
      :EndOutflowStageRelation
      :AreaStageRelation LINEAR
        0.33
      :EndAreaStageRelation
    :EndReservoir

    :Reservoir ExampleReservoir
      :SubBasinID 23
      :HRUID 234
      :Type RESROUTE_STANDARD
      :StageRelations
        21 # number of points
        0.09 0 0 0.0  (h [m], Q [m3/s], V [m3], A [m2])
        0.1 2 43 0.2
        ...
        3.0 20000 3500 200
      :EndStageRelations
    :EndReservoir

    :Reservoir ExampleReservoir
      :SubBasinID 23
      :HRUID 234
      :Type RESROUTE_TIMEVARYING
      :VaryingStageRelations
        21 # number of points
        [jul day1] [jul day2] [jul day3] ...
        0.09 0 0 0.0 0.0 (h [m], Q [m3/s], A [m2], V [m3])
        0.1 2 43 0.2 0.3
        ...
        3.0 20000 3500 200 250
      :EndVaryingStageRelations
    :EndReservoir

    :Reservoir ExampleReservoir # 'lake type format'
      :SubBasinID 23
      :HRUID 234
      :Type RESROUTE_STANDARD
      :WeirCoefficient 0.6
      :CrestWidth 10
      :MaxDepth 5.5 # relative to minimum weir crest elevation
      :LakeArea 10.5 
      :AbsoluteCrestHeight 365 # absolute minimum weir crest height m.a.s.l. (optional) 
    :EndReservoir


  */
  char *s[MAXINPUTITEMS];
  int    Len;

  long SBID=DOESNT_EXIST;

  double a_V(1000.0),b_V(1.0);
  double a_Q(10.0),b_Q(1.0);
  double a_A(1000.0),b_A(0.0);
  double *aQ(NULL),*aQ_ht(NULL); int NQ(0);
  double *aV(NULL),*aV_ht(NULL); int NV(0);
  double *aA(NULL),*aA_ht(NULL); int NA(0);
  double *aQund(NULL);
  int nDates(0);
  double **aQQ=NULL;
  int *aDates(NULL);
  double cwidth(-1),max_depth(-1),weircoeff(-1),lakearea(-1), crestht(0.0);

  curve_function type;
  type=CURVE_POWERLAW;
  res_type restype=RESROUTE_STANDARD;

  CReservoir *pRes=NULL;
  HRUID=DOESNT_EXIST;

  while (!p->Tokenize(s, Len))
  {
    if (Options.noisy){ cout << "-->reading line " << p->GetLineNumber() << ": "; }
    if (Len == 0){}//Do nothing
    else if(IsComment(s[0],Len)){if (Options.noisy){ cout << "#" << endl; }}
    else if (!strcmp(s[0], ":SubBasinID"))
    { 
      if (Options.noisy){ cout << ":SubBasinID" << endl; }
      SBID = s_to_l(s[1]);
    }
    else if (!strcmp(s[0], ":HRUID"))
    {
      if (Options.noisy){ cout << ":HRUID" << endl; }
      HRUID = s_to_i(s[1]);
    }
    else if (!strcmp(s[0], ":Type"))
    {
      if (Options.noisy){ cout << ":Type" << endl; }
      if (!strcmp(s[1], "RESROUTE_STANDARD")){ restype = RESROUTE_STANDARD; }
      else if (!strcmp(s[1], "RESROUTE_NONE")){ restype = RESROUTE_NONE; }
    }
    else if(!strcmp(s[0],":CrestWidth"))
    {
      if (Options.noisy){ cout << ":CrestWidth" << endl; }
      cwidth=s_to_d(s[1]);
      type=CURVE_LAKE;
    }
    else if(!strcmp(s[0],":WeirCoefficient"))
    {
      if (Options.noisy){ cout << ":WeirCoefficient" << endl; }
      weircoeff=s_to_d(s[1]);
      type=CURVE_LAKE;
    }
    else if(!strcmp(s[0],":MaxDepth"))
    {
      if (Options.noisy){ cout << ":MaxDepth" << endl; }
      max_depth=s_to_d(s[1]);
      type=CURVE_LAKE;
    }
    else if(!strcmp(s[0],":LakeArea"))
    {
      if (Options.noisy){ cout << ":LakeArea" << endl; }
      lakearea=s_to_d(s[1]);
      type=CURVE_LAKE;
    }
    else if(!strcmp(s[0],":AbsoluteCrestHeight"))
    {
      if (Options.noisy){ cout << ":AbsoluteCrestHeight" << endl; }
      crestht=s_to_d(s[1]);
      type=CURVE_LAKE;
    }
    //----------------------------------------------------------------------------------------------
    else if (!strcmp(s[0], ":VolumeStageRelation"))
    {
      if (Options.noisy){ cout << ":VolumeStageRelation" << endl; }
      if (Len >= 2){
        if (!strcmp(s[1], "POWER_LAW"))
        {
          type = CURVE_POWERLAW;
          p->Tokenize(s, Len);
          if (Len >= 2){ a_V = s_to_d(s[0]); b_V = s_to_d(s[1]); }
          p->Tokenize(s, Len); //:EndVolumeStageRelation
        }
        else if (!strcmp(s[1], "LINEAR"))
        {
          p->Tokenize(s, Len);
          if (Len >= 1){ a_V = s_to_d(s[0]); b_V = 1.0; }
          p->Tokenize(s, Len); //:EndVolumeStageRelation
        }
        else if (!strcmp(s[1], "LOOKUP_TABLE"))
        {
          type = CURVE_DATA;
          p->Tokenize(s, Len);
          if (Len >= 1){ NV = s_to_i(s[0]); }
          aV = new double[NV];
          aV_ht = new double[NV];
          for (int i = 0; i < NV; i++){
            p->Tokenize(s, Len);
            if(IsComment(s[0],Len)){i--;}
            else{
              aV_ht[i] = s_to_d(s[0]);
              aV[i] = s_to_d(s[1]);
            }
          }
          p->Tokenize(s, Len); //:EndVolumeStageRelation
        }
      }
      else{
        //write warning
      }
    }
    //----------------------------------------------------------------------------------------------
    else if (!strcmp(s[0], ":AreaStageRelation"))
    {
      if (Options.noisy){ cout << ":AreaStageRelation" << endl; }
      if (Len >= 2){
        if (!strcmp(s[1], "POWER_LAW"))
        {
          type = CURVE_POWERLAW;
          p->Tokenize(s, Len);
          if (Len >= 2){ a_A = s_to_d(s[0]); b_A = s_to_d(s[1]); }
          p->Tokenize(s, Len); //:EndAreaStageRelation
        }
          else if (!strcmp(s[1], "CONSTANT"))
        {
          type = CURVE_LINEAR;
          p->Tokenize(s, Len);
          if (Len >= 1){ a_A = s_to_d(s[0]); b_A = 0.0; }
          p->Tokenize(s, Len); //:EndAreaStageRelation
        }
        else if (!strcmp(s[1], "LINEAR"))
        {
          type = CURVE_LINEAR;
          p->Tokenize(s, Len);
          if (Len >= 1){ a_A = s_to_d(s[0]); b_A = 1.0; }
          p->Tokenize(s, Len); //:EndAreaStageRelation
        }
        else if (!strcmp(s[1], "LOOKUP_TABLE"))
        {
          type = CURVE_DATA;
          p->Tokenize(s, Len);
          if (Len >= 1){ NA = s_to_i(s[0]); }
          aA = new double[NA];
          aA_ht = new double[NA];
          for (int i = 0; i < NA; i++){
            p->Tokenize(s, Len);
            if(IsComment(s[0],Len)){i--;}
            else{
              aA_ht[i] = s_to_d(s[0]);
              aA[i] = s_to_d(s[1]);
            }
          }
          p->Tokenize(s, Len); //:EndAreaStageRelation
        }
      }
      else{
        //write warning
      }
    }
    //----------------------------------------------------------------------------------------------
    else if (!strcmp(s[0], ":OutflowStageRelation"))
    {
      if (Options.noisy){ cout << ":OutflowStageRelation" << endl; }
      if (Len >= 2){
        if (!strcmp(s[1], "POWER_LAW"))
        {
          type = CURVE_POWERLAW;
          p->Tokenize(s, Len);
          if (Len >= 2){ a_Q = s_to_d(s[0]); b_Q = s_to_d(s[1]); }
          p->Tokenize(s, Len); //:EndOutflowStageRelation
        }
        else if (!strcmp(s[1], "LINEAR"))
        {
          type = CURVE_LINEAR;
          p->Tokenize(s, Len);
          if (Len >= 1){ a_Q = s_to_d(s[0]); b_Q = 1.0; }
          p->Tokenize(s, Len); //:EndOutflowStageRelation
        }
        else if (!strcmp(s[0], "LOOKUP_TABLE"))
        {
          type = CURVE_DATA;
          p->Tokenize(s, Len);
          if (Len >= 1){ NQ = s_to_i(s[0]); }
          aQ = new double[NQ];
          aQ_ht = new double[NQ];
          for (int i = 0; i < NQ; i++){
            p->Tokenize(s, Len);
            if(IsComment(s[0],Len)){i--;}
            else{
              aQ_ht[i] = s_to_d(s[0]);
              aQ[i] = s_to_d(s[1]);
            }
          }
          p->Tokenize(s, Len); //:EndOutflowStageRelation
        }
      }
      else{
        //write warning
      }
    }
    //----------------------------------------------------------------------------------------------
    else if (!strcmp(s[0], ":StageRelations"))
    {
      if (Options.noisy){ cout << ":StageRelations" << endl; }
      type = CURVE_DATA;
      p->Tokenize(s, Len);
      if (Len >= 1){ NQ = s_to_i(s[0]); }

      aQ_ht = new double[NQ];
      aQ = new double[NQ];
      aV = new double[NQ];
      aA = new double[NQ];
      aQund=new double[NQ];
      for (int i = 0; i < NQ; i++){
        p->Tokenize(s, Len);
        if(IsComment(s[0],Len)){i--;}
        else{
          if(Len>=4){
            aQ_ht[i] = s_to_d(s[0]);
            aQ[i] = s_to_d(s[1]);
            aV[i] = s_to_d(s[2]);
            aA[i] = s_to_d(s[3]);
            aQund[i]=0.0;
            if(Len>=5){
              aQund[i]=s_to_d(s[4]);
            }
          }
          else{
            WriteWarning("Incorrect line length (<4) in :Reservoir :StageRelations command",Options.noisy);
          }
        }
      }
      p->Tokenize(s, Len); //:EndStageRelations
    }
    //----------------------------------------------------------------------------------------------
    else if(!strcmp(s[0],":EndStageRelations"))
    {
    }
    //----------------------------------------------------------------------------------------------
    else if (!strcmp(s[0], ":VaryingStageRelations"))
    {
      bool hasQund=false;
      int shift=0;
      if (Options.noisy){ cout << ":VaryingStageRelations" << endl; }
      type = CURVE_VARYING;
      p->Tokenize(s, Len);
      if (Len >= 1){ NQ = s_to_i(s[0]); } //# of flows
      if(Len>=2){hasQund=true;shift=1;} //TMP DEBUG - shoudld have flag for extra column
      p->Tokenize(s, Len);
      nDates = Len;
      
      aDates = new int[nDates];
      for (int v = 0; v < nDates; v++){ aDates[v] = s_to_i(s[v]); }

      aQ_ht = new double[NQ];
      aQQ = new double *[nDates];
      for (int v = 0; v < nDates; v++){
        aQQ[v] = new double[NQ];
      }
      aV = new double[NQ];
      aA = new double[NQ];
      aQund=new double[NQ];
      for (int i = 0; i < NQ; i++){
        p->Tokenize(s, Len);
        if(IsComment(s[0],Len)){i--;}
        else{
          if(Len==1){//presumably :EndVaryingStageRelations
            ExitGracefully("CReservoir::Parse: improper number of rows in :VaryingStageRelations command", BAD_DATA);
          }
          else if (Len < 2 + nDates){
            ExitGracefully("CReservoir::Parse: improper number of columns in :VaryingStageRelations command", BAD_DATA);
          }
          aQ_ht[i] = s_to_d(s[0]); //stage
          aV[i] = s_to_d(s[1]); //volume
          aA[i] = s_to_d(s[2]); //area
          if(hasQund){ aQund[i] = s_to_d(s[3]); } //Qund
          else       { aQund[i]=0.0; }
          for (int v = 0; v < nDates; v++){
            aQQ[v][i] = s_to_d(s[3 + v +shift]); //flows for each date
          }
        }
      }
      p->Tokenize(s, Len); //:EndVaryingStageRelations
    }
    //----------------------------------------------------------------------------------------------
    else if (!strcmp(s[0], ":EndReservoir")){
      if (Options.noisy){ cout << ":EndReservoir" << endl; }
      break;
    }
    else{
      WriteWarning("Reservoir::Parse: unrecognized command (" + to_string(s[0]) + ") in :Reservoir-:EndReservoir Block",Options.noisy);
    }
  }

  if ((type==CURVE_POWERLAW) || (type==CURVE_LINEAR))
  {
    pRes=new CReservoir(name,SBID,restype,a_V,b_V,a_Q,b_Q,a_A,b_A);
  }
  else if (type==CURVE_DATA)
  {
    pRes=new CReservoir(name,SBID,restype,aQ_ht,aQ,aQund,aA,aV,NQ);//presumes aQ_ht=aV_ht=aA_ht; NA=NV=NQ
  }
  else if (type==CURVE_VARYING)
  {
    pRes=new CReservoir(name,SBID,restype,nDates,aDates,aQ_ht,aQQ,aQund,aA,aV,NQ);//presumes aQ_ht=aV_ht=aA_ht; NA=NV=NQ
  }
  else if(type==CURVE_LAKE)
  {
    ExitGracefullyIf(weircoeff==-1,"CReservoir::Parse: :WeirCoefficient must be specified for lake-type reservoirs",BAD_DATA_WARN);
    ExitGracefullyIf(cwidth==-1,"CReservoir::Parse: :CrestWidth must be specified for lake-type reservoirs",BAD_DATA_WARN);
    ExitGracefullyIf(lakearea==-1,"CReservoir::Parse: :LakeArea  must be specified for lake-type reservoirs",BAD_DATA_WARN);
    ExitGracefullyIf(max_depth==-1,"CReservoir::Parse: :LakeDepth must be specified for lake-type reservoirs",BAD_DATA_WARN);
    pRes=new CReservoir(name,SBID,restype,weircoeff,cwidth,crestht,lakearea,max_depth);
  }
  else{
    ExitGracefully("CReservoir::Parse: only currently supporting linear, powerlaw, or data reservoir rules",STUB);
  }

  for (int i = 0; i < nDates; i++){ delete[] aQQ[i]; }delete [] aQQ;
  delete [] aQ;
  delete [] aQ_ht;
  delete [] aQund;
  delete [] aV;
  delete [] aV_ht;
  delete [] aA;
  delete [] aA_ht;
  return pRes;
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
/// \brief interpolates value from rating curve
/// \param x [in] interpolation location
/// \param xx [in] array of vertices ordinates of interpolant
/// \param y [in] array (size:N)  of values corresponding to N evenly spaced points in x
/// \param N size of array y
/// \returns y value corresponding to interpolation point
/// \note does not assume regular spacing between min and max x value
//
double Interpolate2(double x, double *xx, double *y, int N, bool extrapbottom)
{
  static int ilast=0;
  if      (x<=xx[0])
  {
    if (extrapbottom){return y[0]+(y[1]-y[0])/(xx[1]-xx[0])*(x-xx[0]);}
    return y[0];
  }
  else if (x>=xx[N-1])
  {
    return y[N-1]+(y[N-1]-y[N-2])/(xx[N-1]-xx[N-2])*(x-xx[N-1]);//extrapolation-may wish to revisit
    //return y[N-1];
  }
  else
  {
    //int i=0; while ((x>xx[i+1]) && (i<(N-2))){i++;}//Dumb Search
    int i=SmartIntervalSearch(x,xx,N,ilast);
    if(i==DOESNT_EXIST){return 0.0;}
    ExitGracefullyIf(i==DOESNT_EXIST,"Interpolate2::mis-ordered list or infinite x",RUNTIME_ERR);
    ilast=i;
    return y[i]+(y[i+1]-y[i])/(xx[i+1]-xx[i])*(x-xx[i]);
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
double     CReservoir::GetOutflow(const double &ht, const double &adj) const
{
  double underflow=Interpolate2(ht,_aStage,_aQunder,_Np,false); //no adjustments
  return Interpolate2(ht-adj,_aStage,_aQ,_Np,false)+underflow;
}


