/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2024 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Reservoir.h"
#include "Model.h"     // needed to define CModel

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

  _pDownSB=NULL;

  _lakebed_thick=1.0;
  _lakebed_cond =0.0;
  _lake_convcoeff=2.0;

  _stage     =0.0;
  _stage_last=0.0;
  _min_stage =0.0;
  _max_stage =0.0;
  _Qout      =0.0;
  _Qout_last =0.0;
  _MB_losses =0.0;
  _AET		   =0.0;
  _Precip    =0.0;
  _GW_seepage=0.0;
  _aQstruct=NULL;
  _aQstruct_last=NULL;

  _pDemandTS=NULL;
  _nDemandTS=0;

  _pHRU=NULL;

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

  _Qoptimized=RAV_BLANK_DATA;
  _aQdelivered=NULL;

  _pDZTR=NULL;

  _minStageDominant=false;

  _crest_ht=0.0;
  _crest_width=DOESNT_EXIST;

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

  _nControlStructures=0;
  _pControlStructures=NULL;

  _seepage_const=0;
  _local_GW_head=0.0;

  _assimilate_stage=false;
  _assim_blank=true;
  _pObsStage=NULL;
  _DAscale=1.0;
  _DAscale_last=1.0;

  _dry_timesteps=0;
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
/// \param SBID [in] subbasin ID
/// \param a_V [in] volume of reservoir/lake when stage at absolute crest height [m3]
/// \param b_V [in] power law exponent for volume rating curve [-]
/// \param a_Q [in] power law coefficient for discharge rating curve [m3/s*m^-b_Q]
/// \param b_Q [in] power law exponent for discharge rating curve [-]
/// \param a_A[in] surface area of reservoir/lake when stage at absolute crest height [m2]
/// \param b_A [in] power law exponent for area rating curve [-] (0 for prismatic reservoir)

//
CReservoir::CReservoir(const string Name, const long SBID,
                       const double a_V,  const double b_V,
                       const double a_Q,  const double b_Q,
                       const double a_A,  const double b_A,
                       const double crestht,const double depth)
{
  BaseConstructor(Name,SBID);

  _Np=100;
  _aStage =new double[_Np];
  _aQ     =new double[_Np];
  _aQunder=new double[_Np];
  _aArea  =new double[_Np];
  _aVolume=new double[_Np];
  ExitGracefullyIf(_aVolume==NULL,"CReservoir::constructor",OUT_OF_MEMORY);

  double aA=a_A;
  double bA=b_A;
  if(aA==AUTO_COMPUTE) {aA=a_V/depth*b_V;}
  if(bA==AUTO_COMPUTE) {bA=b_V-1.0;      }

  double ht,dh;
  _crest_ht  =crestht;       // zero by default
  _min_stage =_crest_ht-depth;
  _max_stage =_crest_ht+6.0; // a postive value relative to _crest_ht
  dh=(_max_stage-_min_stage)/(double)(_Np-1);

  for(int i=0;i<_Np;i++)
  {
    ht  =_min_stage+(double)(i)*dh;
    if (((ht+dh)-_crest_ht)*(ht-_crest_ht)<0.0){ht=_crest_ht;}//ensures zero discharge point included
    double ht_above_bottom=ht-_min_stage;
    _aStage [i]=ht;
    _aQ     [i]=a_Q*pow(max(ht-_crest_ht,0.0),b_Q);
    _aArea  [i]=aA *pow(ht_above_bottom/depth,bA );
    _aVolume[i]=a_V*pow(ht_above_bottom/depth,b_V);
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

    if (_aQ[i]==0){_crest_ht=_aStage[i];}

    // QA/QC:
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
    if (_aQ[i]==0){_crest_ht=_aStage[i]; } 
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
/// \brief Constructor for prismatic lake reservoir controlled by weir coefficient
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

  _crest_width=crestw;
  _crest_ht   =crestht;
  _min_stage  =-depth+_crest_ht;
  _max_stage  =5.0+_crest_ht; //reasonable default?
  ExitGracefullyIf(depth <=0, "CReservoir::Constructor (Lake): cannot have negative maximum lake depth",BAD_DATA_WARN);
  ExitGracefullyIf(A     <=0, "CReservoir::Constructor (Lake): cannot have negative or zero lake area", BAD_DATA_WARN);
  ExitGracefullyIf(crestw<=0, "CReservoir::Constructor (Lake): cannot have negative or zero crest width", BAD_DATA_WARN);
  ExitGracefullyIf(weircoeff <= 0, "CReservoir::Constructor (Lake): cannot have negative or zero weir discharge coefficient", BAD_DATA_WARN);
  if (weircoeff>1.0){
    WriteWarning( "CReservoir::Constructor (Lake): weir discharge coefficient should be less than 1.0", true);
  }
  _Np=102;
  _aStage =new double [_Np];
  _aQ     =new double [_Np];
  _aQunder=new double [_Np];
  _aArea  =new double [_Np];
  _aVolume=new double [_Np];
  ExitGracefullyIf(_aVolume==NULL,"CReservoir::constructor (4)",OUT_OF_MEMORY);

  double dh;
  string warn;
  dh=(_max_stage-_crest_ht)/(_Np-2); //spacing = 0.05m
  _aStage [0]=_min_stage; //first point is for dry reservoir, linear interpolation of volume between this and crest height
  _aQ     [0]=0.0;
  _aQunder[0]=0.0;
  _aArea  [0]=A;
  _aVolume[0]=0.0;
  for (int i=1;i<_Np;i++) // - Edited by KL to reference crest height properly
  {
    _aStage [i]=_crest_ht+(i-1)*dh;
    _aQ     [i]=2.0/3.0*weircoeff*sqrt(2*GRAVITY)*crestw*pow((_aStage[i]-_crest_ht),1.5); //Overflow weir equation
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
  delete [] _aQ;          _aQ     =NULL;
  delete [] _aQunder;     _aQunder=NULL;
  delete [] _aArea;       _aArea  =NULL;
  delete [] _aVolume;     _aVolume=NULL;
  for (int v = 0; v<_nDates; v++){ delete[] _aQ_back[v]; } delete [] _aQ_back; _aQ_back=NULL;
  delete [] _aDates;      _aDates =NULL;

  for (int i=0; i<_nDemandTS;i++){delete _pDemandTS[i]; } delete [] _pDemandTS; _pDemandTS=NULL;

  delete _pWeirHeightTS;  _pWeirHeightTS=NULL;
  delete _pMaxStageTS;    _pMaxStageTS=NULL;
  delete _pOverrideQ;     _pOverrideQ=NULL;
  delete _pMinStageTS;    _pMinStageTS=NULL;
  delete _pMinStageFlowTS;_pMinStageFlowTS=NULL;
  delete _pTargetStageTS; _pTargetStageTS=NULL;
  delete _pMaxQIncreaseTS;_pMaxQIncreaseTS=NULL;
  delete _pMaxQDecreaseTS;_pMaxQDecreaseTS=NULL;
  delete _pDroughtLineTS; _pDroughtLineTS=NULL;
  delete _pQminTS;        _pQminTS=NULL;
  delete _pQmaxTS;        _pQmaxTS=NULL;
  delete _pQdownTS;       _pQdownTS=NULL;

  delete [] _aQstruct;      _aQstruct=NULL;
  delete [] _aQstruct_last; _aQstruct_last=NULL;
  delete [] _aQdelivered;   _aQdelivered=NULL;

  for (int i=0;i<_nControlStructures;i++){delete _pControlStructures[i];} delete [] _pControlStructures; _pControlStructures=NULL;
}

//////////////////////////////////////////////////////////////////
/// \returns Subbasin ID
//
long  CReservoir::GetSubbasinID          () const { return _SBID; }

//////////////////////////////////////////////////////////////////
/// \returns reservoir storage [m3]
//
double  CReservoir::GetStorage           () const { return GetVolume(_stage); }

//////////////////////////////////////////////////////////////////
/// \returns current stage [m]
//
double  CReservoir::GetResStage          () const { return _stage; }

//////////////////////////////////////////////////////////////////
/// \returns previous stage [m]
//
double  CReservoir::GetOldStage          () const { return _stage_last; }

//////////////////////////////////////////////////////////////////
/// \returns min stage [m]
//
double  CReservoir::GetMinStage          (const int nn) const {
  if(_pMinStageTS    !=NULL){ return _pMinStageTS->GetSampledValue(nn);}
  return _min_stage;
}
//////////////////////////////////////////////////////////////////
/// \returns max stage [m]
//
double  CReservoir::GetMaxStage          (const int nn) const {
  if(_pMaxStageTS    !=NULL){ return _pMaxStageTS->GetSampledValue(nn);}
  return ALMOST_INF;
}
//////////////////////////////////////////////////////////////////
/// \returns sill elevation [m]
/// supports time-variable discharge curves and weir height adjustments 
//
double  CReservoir::GetSillElevation(const int nn) const 
{
  double weir_adj=0.0;
  if (_pWeirHeightTS!=NULL){
    weir_adj=_pWeirHeightTS->GetSampledValue(nn);
  }
  return _crest_ht+weir_adj;
}
//////////////////////////////////////////////////////////////////
/// \returns current surface area [m2]
//
double  CReservoir::GetSurfaceArea       () const { return GetArea(_stage); }

//////////////////////////////////////////////////////////////////
/// \returns start-of-timestep surface area [m2]
//
double  CReservoir::GetOldSurfaceArea    () const { return GetArea(_stage_last); }

//////////////////////////////////////////////////////////////////
/// \returns lakebed thickness [m]
//
double  CReservoir::GetLakebedThickness() const { return _lakebed_thick; }

//////////////////////////////////////////////////////////////////
/// \returns lakebed thermal conductivity [MJ/m/K/d]
//
double  CReservoir::GetLakebedConductivity() const { return _lakebed_cond; }

//////////////////////////////////////////////////////////////////
/// \returns lake thermal convection coefficient [MJ/m2/d/K]
//
double  CReservoir::GetLakeConvectionCoeff() const { return _lake_convcoeff; }

//////////////////////////////////////////////////////////////////
/// \returns current outflow rate [m3/s]
//
double  CReservoir::GetOutflowRate       () const { return _DAscale*_Qout; }

//////////////////////////////////////////////////////////////////
/// \returns previous outflow rate [m3/s]
//
double CReservoir::GetOldOutflowRate     () const { return _DAscale_last*_Qout_last; }

//////////////////////////////////////////////////////////////////
/// \returns current outflow rate from control structure i[m3/s]
//
double  CReservoir::GetControlOutflow    (const int i) const { return _aQstruct[i]; }

//////////////////////////////////////////////////////////////////
/// \returns outflow integrated over timestep [m3/d]
//
double  CReservoir::GetIntegratedOutflow(const double& tstep) const
{
  return 0.5*(_DAscale*_Qout+_DAscale_last*_Qout_last)*(tstep*SEC_PER_DAY); //integrated
}

//////////////////////////////////////////////////////////////////
/// \returns control structure outflow integrated over timestep [m3/d]
//
double  CReservoir::GetIntegratedControlOutflow(const int i, const double& tstep) const
{
  return 0.5*(_aQstruct[i]+_aQstruct_last[i])*(tstep*SEC_PER_DAY); //integrated
}
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
/// \returns precip gains integrated over previous timestep [m3]
//
double  CReservoir::GetReservoirPrecipGains(const double& tstep) const
{
  return _Precip;
}
//////////////////////////////////////////////////////////////////
/// \returns seepage losses only integrated over previous timestep [m3]
//
double  CReservoir::GetReservoirGWLosses(const double &tstep) const
{
  return _GW_seepage;
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
/// \returns number of control structures
//
int CReservoir::GetNumControlStructures() const
{
  return _nControlStructures;
}
//////////////////////////////////////////////////////////////////
/// \returns crest width, in meters
//
double CReservoir::GetCrestWidth() const
{
  return _crest_width;
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
/// \returns multiplies discharge curve by correction factor
//
void CReservoir::MultiplyFlow(const double &mult)
{
  for(int i=0;i<_Np;i++) {
      _aQ[i]*=mult;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief returns current active regime name
/// \param i [in] control structure index (must be in range from 0...nControlStructures)
/// \param tt [in] current time structure
/// returns current active regime name (or "NONE" if none active)
//
string CReservoir::GetRegimeName(const int i,const time_struct &tt) const
{
  return _pControlStructures[i]->GetCurrentRegimeName(tt);
}

//////////////////////////////////////////////////////////////////
/// \brief returns subbasin ID of recipient of control structure i's outflow
/// \param i [in] control structure index (must be in range from 0...nControlStructures)
//
long   CReservoir::GetControlFlowTarget(const int i) const
{
  return _pControlStructures[i]->GetTargetBasinID();
}
//////////////////////////////////////////////////////////////////
/// \brief returns name of  control structure i
/// \param i [in] control structure index (must be in range from 0...nControlStructures)
//
string CReservoir::GetControlName(const int i) const
{
  return _pControlStructures[i]->GetName();
}

double CReservoir::GetStageDischargeDerivative(const double &stage, const int nn) const
{
  double dh=0.001;
  double weir_adj=0.0;
  if(_pWeirHeightTS!=NULL) {
    weir_adj=_pWeirHeightTS->GetSampledValue(nn);
  }

  return (GetWeirOutflow(stage+dh,weir_adj)-GetWeirOutflow(stage,weir_adj))/dh;
}
double CReservoir::GetDischargeFromStage(const double &stage, const int nn) const
{
  double weir_adj=0.0;
  if(_pWeirHeightTS!=NULL) {
    weir_adj=_pWeirHeightTS->GetSampledValue(nn);
  }

  return GetWeirOutflow(stage,weir_adj);
}

 int  CReservoir::GetNumDryTimesteps    () const{return _dry_timesteps;}

//////////////////////////////////////////////////////////////////
/// \brief number of water/irrigation demands
/// \return number of water/irrigation demands
//
int    CReservoir::GetNumWaterDemands() const {
  return _nDemandTS;
}
//////////////////////////////////////////////////////////////////
/// \brief returns water/irrigation demand integer ID
/// \return water/irrigation demand integer ID
//
int    CReservoir::GetWaterDemandID(const int ii) const
{
  #ifdef _STRICTCHECK_
  ExitGracefullyIf(ii < 0 || ii >= _nDemandTS, "CReservoir::GetWaterDemandID: invalid index",RUNTIME_ERR);
#endif
  return _pDemandTS[ii]->GetIDTag(); //stores demand ID
}
//////////////////////////////////////////////////////////////////
/// \brief returns water/irrigation demand name/alias
/// \return water/irrigation demand name/alias
//
string CReservoir::GetWaterDemandName      (const int ii) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf(ii < 0 || ii >= _nWaterDemands, "CSubBasin::GetWaterDemandName: invalid index",RUNTIME_ERR);
#endif
  return _pDemandTS[ii]->GetName();
}
//////////////////////////////////////////////////////////////////
/// \brief Returns specified irrigation/water use demand from reservoir at time t
/// \param &t [in] Model time at which the demand from reservoir is to be determined
/// \return specified demand from reservoir at time t
//
double CReservoir::GetWaterDemand           (const int ii,const double &t) const
{
  double Qirr;
  Qirr=_pDemandTS[ii]->GetValue(t);
  if (Qirr==RAV_BLANK_DATA){Qirr=0.0;}
  return Qirr;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns instantaneous ACTUAL irrigation use at end of current timestep
/// \return actual delivery from reservoir [m3/s]
//
double CReservoir::GetDemandDelivery(const int ii) const
{
  return _aQdelivered[ii];
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

  for (int i=0; i<_nDemandTS; i++)
  {
    _pDemandTS[i]->Initialize(model_start_day, model_start_yr, model_duration, timestep, false, Options.calendar);
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
  _aQdelivered = new double [_nDemandTS];
  for (int ii = 0; ii < _nDemandTS; ii++) {
    _aQdelivered[ii]=0.0;
  }
  _aQstruct     =new double [_nControlStructures];
  _aQstruct_last=new double [_nControlStructures];
  for (int i = 0; i < _nControlStructures; i++) {
    _aQstruct[i]=_aQstruct_last[i]=0.0;
  }
  _dry_timesteps=0;
}

//////////////////////////////////////////////////////////////////
/// \brief Adds extraction history
/// \param *pDemand demand time series to be added
//
void    CReservoir::AddDemandTimeSeries (CTimeSeries *pDemand)
{
  if (!DynArrayAppend((void**&)(_pDemandTS),(void*)(pDemand),_nDemandTS)){
    ExitGracefully("CReservoir::AddDemandTimeSeries : trying to add NULL time series",BAD_DATA);}
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
/// \brief Adds control structure outlet to Reservoir
/// \param pCS [in] pointer to valid control structure object
//
void  CReservoir::AddControlStructure(const CControlStructure* pCS)
{
  if (!DynArrayAppend((void**&)_pControlStructures, (void*)pCS, _nControlStructures)) {
    ExitGracefully("CReservoir::AddControlStructure: adding null control structure",BAD_DATA_WARN);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief links reservoir to HRU
/// \param *pHRUpointer HRU to link to
//
void  CReservoir::SetHRU(const CHydroUnit *pHRUpointer)
{
  _pHRU=pHRUpointer;
  ExitGracefullyIf(_pHRU->GetArea()<=0,"CReservoir::SetHRU: reservoir cannot be linked to a zero-area HRU.",BAD_DATA_WARN);
}
//////////////////////////////////////////////////////////////////
/// sets crest width, in meters
//
void CReservoir::SetCrestWidth(const double& width)
{
  if(_crest_width==DOESNT_EXIST) {
    WriteWarning("CReservoir::SetCrestWidth: - trying to change crest width of non-lake reservoir",true);
    return;
  }
  else {
    MultiplyFlow(width/_crest_width);
    _crest_width=width;
  }
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
/// \brief sets total incoming precip volume, in m3, for timestep
/// \param precip_m3: precip, in m3
//
void  CReservoir::SetPrecip(const double& precip_m3)
{
  _Precip=precip_m3;
}
//////////////////////////////////////////////////////////////////
/// \brief sets subbasin to which flows from main outlet go
//
void  CReservoir::SetDownstreamBasin(const CSubBasin* pDownBasin) {
  _pDownSB=pDownBasin;
}
//////////////////////////////////////////////////////////////////
/// \brief sets lakebed thickness
//
void CReservoir::SetLakebedThickness(const double& thick){
  _lakebed_thick=thick;
}
//////////////////////////////////////////////////////////////////
/// \brief sets lakebed thermal conductivity
//
void CReservoir::SetLakebedConductivity(const double& cond) {
  _lakebed_cond=cond;
}
//////////////////////////////////////////////////////////////////
/// \brief sets lake convection coeff
//
void CReservoir::SetLakeConvectionCoeff(const double& conv) {
  _lake_convcoeff=conv;
}
//////////////////////////////////////////////////////////////////
/// \brief gets current constraint name
/// \returns current constraint applied to estimate stage/flow
//
string CReservoir::GetCurrentConstraint() const {
  switch(_constraint)
  {
  case(RC_MAX_STAGE):         {return "MAX_STAGE"; }
  case(RC_MIN_STAGE):         {return "RC_MIN_STAGE"; }
  case(RC_NATURAL):           {return "RC_NATURAL"; }
  case(RC_TARGET):            {return "RC_TARGET"; }
  case(RC_DOWNSTREAM_FLOW):   {return "RC_DOWNSTREAM_FLOW"; }
  case(RC_MAX_FLOW_INCREASE): {return "RC_MAX_FLOW_INCREASE"; }
  case(RC_MAX_FLOW_DECREASE): {return "RC_MAX_FLOW_DECREASE"; }
  case(RC_MIN_FLOW):          {return "RC_MIN_FLOW"; }
  case(RC_MAX_FLOW):          {return "RC_MAX_FLOW"; }
  case(RC_OVERRIDE_FLOW):     {return "RC_OVERRIDE_FLOW"; }
  case(RC_DRY_RESERVOIR):     {return "RC_DRY_RESERVOIR"; }
  case(RC_MANAGEMENT):        {return "RC_MANAGEMENT"; }
  case(RC_DZTR):              {return "RC_DZTR";}
  default:                    {return ""; }
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
    _aVolume[i]=InterpolateCurve(_aStage[i],a_ht,a_V,nPoints,false);
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
    _aArea[i]=InterpolateCurve(_aStage[i],a_ht,a_A,nPoints,false);
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
/// \brief sets data assimilation scale factors (read from .rvc file)
//
void CReservoir::SetDataAssimFactors(const double& da_scale,const double& da_scale_last)
{
  _DAscale     =da_scale;
  _DAscale_last=da_scale_last;
}

//////////////////////////////////////////////////////////////////
/// \brief enables lake stage assimilation
/// \param pObs -time series of observed lake stage (can have NULL entries)
//
void CReservoir::TurnOnAssimilation(CTimeSeriesABC *pObs)
{
  _assimilate_stage=true;
  _pObsStage=pObs;
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
/// Hydrol. Earth Syst. Sci., 23, 3735-3764, 2019 https://doi.org/10.5194/hess-23-3735-2019
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
/// \brief scales all internal flows by scale factor (for assimilation/nudging)
/// \remark Messes with mass balance something fierce!
///
/// \return mass added to system [m3]
//
double CReservoir::ScaleFlow(const double& scale,const bool overriding, const double& tstep, const double &t)
{
  double va=0.0; //volume added
  double sf=(scale-1.0)/scale;

  _DAscale=scale;

  //Estimate volume added through scaling
  va+=0.5*(_Qout_last+_Qout)*sf*tstep*SEC_PER_DAY;

  return va;
}

//////////////////////////////////////////////////////////////////
/// \brief updates state variable "stage" and other reservoir auxiliary variables at end of computational time step
/// \param new_stage [in] calculated stage at end of time step
//
void  CReservoir::UpdateStage(const double &new_stage,const double &res_outflow,const res_constraint &constr,const double *aQstruct_new,const optStruct &Options,const time_struct &tt)
{
  _stage_last=_stage;
  _stage     =new_stage;

  _constraint=constr;

  if (_constraint==RC_DRY_RESERVOIR){_dry_timesteps++;}

  _Qout_last =_Qout;
  _Qout      =res_outflow;
  for (int i = 0; i < _nControlStructures; i++) {
    _aQstruct_last[i]=_aQstruct[i];
    _aQstruct[i]=aQstruct_new[i];
  }

  _DAscale_last=_DAscale;
  _DAscale   =1.0;
}
//////////////////////////////////////////////////////////////////
/// \brief returns AET [mm/d], meant to be called at end of time step
//
double CReservoir::GetAET() const
{
  if(_pHRU!=NULL) {
    double Evap=_pHRU->GetForcingFunctions()->OW_PET;//mm/d
    Evap*=_pHRU->GetSurfaceProps()->lake_PET_corr;
    return Evap*0.5*(GetArea(_stage)+GetArea(_stage_last))/(_pHRU->GetArea()*M2_PER_KM2); //normalized to HRU area
  }
  else {
    return 0.0;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief updates current mass balance (called at end of time step)
/// \param tt [in] current time step information
/// \param tstep [in] time step, in days
//
void CReservoir::UpdateMassBalance(const time_struct &tt,const double &tstep, const optStruct &Options)
{
  _MB_losses=0.0; //all losses outside the system

  _AET = 0.0;
  if (_pHRU!=NULL){
    _AET = GetAET()*(_pHRU->GetArea()*M2_PER_KM2)/ MM_PER_METER * tstep; //m3;
  }
  _MB_losses+=_AET;

  if(_seepage_const>0) {
    _GW_seepage=_seepage_const*(0.5*(_stage+_stage_last)-_local_GW_head)*SEC_PER_DAY*tstep;
    _MB_losses+=_GW_seepage;
  }

  for (int ii = 0; ii < _nDemandTS; ii++) {
    if (Options.management_optimization) {
      _MB_losses+=_aQdelivered[ii];
    }
    else{
      int nn        =(int)((tt.model_time+TIME_CORRECTION)/tstep);//current timestep index
      _MB_losses+=_pDemandTS[ii]->GetSampledValue(nn) * SEC_PER_DAY * tstep;
    }
  }


}
//////////////////////////////////////////////////////////////////
/// \brief adds delivered demand amount to water demand ii
/// \param ii [in] local demand index
/// \param Qdel [in] delivered flow [m3/s]
//
void CReservoir::AddToDeliveredDemand(const int ii, const double &Qdel)
{
  _aQdelivered[ii]=Qdel;
}
//////////////////////////////////////////////////////////////////
/// \brief updates rating curves based upon the time and assimilates lake stage
/// \notes can later support quite generic temporal changes to treatmetn of outflow-volume-area-stage relations
/// \param tt [in] current model time
/// \param Options [in] model options structure
//
void CReservoir::UpdateReservoir(const time_struct &tt, const optStruct &Options)
{
  // update flow rules (for time-variable stage-discharge curve)-------
  // also updates crest height
  if (_nDates != 0){
    int vv=_nDates-1;
    for (int v = 0; v < _nDates; v++){
      if (tt.julian_day >= _aDates[v]){vv=v; }
    }
    for (int i = 0; i < _Np; i++){
      _aQ[i] = _aQ_back[vv][i];
      if (_aQ[i]==0.0){_crest_ht=_aStage[i];}
    }
  }

  // Assimilate lake stage-------------------------------
  if(_assimilate_stage)
  {
    _assim_blank=true;
    if(tt.model_time>Options.assimilation_start-Options.timestep/2.0)
    {
      int nn=(int)((tt.model_time+TIME_CORRECTION)/Options.timestep)+1;//end-of timestep index

      double weir_adj=0.0;
      if(_pWeirHeightTS!=NULL) {
        weir_adj=_pWeirHeightTS->GetSampledValue(nn);
      }

      double obs_stage=_pObsStage->GetSampledValue(nn);
      if(obs_stage!=RAV_BLANK_DATA) {
        _stage=obs_stage;
        _Qout =GetWeirOutflow(_stage,weir_adj);//[m3/s]
        _assim_blank=false;
      }
    }
    //Calculate change in reservoir mass
    // \todo[funct] - correct mass balance for assimlating lake stage
  }

  // reboot delivered amount to zero each tstep ------------------
  for (int ii=0; ii<_nDemandTS; ii++){
    _aQdelivered[ii]=0.0;
  }
  return;
}
//////////////////////////////////////////////////////////////////
/// \brief initialize stage,volume, area to specified initial inflow
/// \param initQ [in] initial inflow
/// \note - ignores extraction and PET ; won't work for initQ=0, which could be non-uniquely linked to stage
//
void  CReservoir::SetInitialFlow(const double &initQ,const double &initQlast,const time_struct &tt,const optStruct& Options)
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
    weir_adj=_pWeirHeightTS->GetValue(tt.model_time);
  }

  for (int i=0;i<_Np;i++)
  {
    if (_aQ[i]>0.0){h_guess=_min_stage+(double)(i)/(double)(_Np)*(_max_stage-_min_stage);break;}//initialize to just topping over crest
  }
  int    iter=0;
  double change=0;
  double Q,dQdh,Qi;

  do //Newton's method with discrete approximation of dQ/dh
  {
    if(_pDZTR==NULL) {
      Q    =GetWeirOutflow(h_guess,   weir_adj);//[m3/s]
      dQdh=(GetWeirOutflow(h_guess+dh,weir_adj)-Q)/dh;
      for (int i = 0; i < _nControlStructures; i++) {
        Qi=_pControlStructures[i]->GetOutflow(h_guess,h_guess,_aQstruct_last[i],tt);
        Q+=Qi;
        dQdh+=(_pControlStructures[i]->GetOutflow(h_guess+dh,h_guess+dh,_aQstruct_last[i],tt)-Qi)/dh;
      }
    }
    else if(_pDZTR!=NULL) {
      Q =GetDZTROutflow(GetVolume(h_guess),initQ,tt,Options);
      dQdh=(GetDZTROutflow(GetVolume(h_guess+dh),initQ,tt,Options)-Q)/dh;
    }

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
  //if(_pDZTR==NULL) {
  _Qout=GetWeirOutflow(_stage,weir_adj);
  for (int i = 0; i < _nControlStructures; i++) {
    _aQstruct[i]=_pControlStructures[i]->GetOutflow(_stage,_stage,_aQstruct_last[i],tt);
  }
  //}
  //else if(_pDZTR!=NULL) {
  //  _Qout=GetDZTROutflow(GetVolume(_stage),initQ,tt,Options);
  //}

  _stage_last=_stage;
  _Qout_last =_Qout;
}

//////////////////////////////////////////////////////////////////
/// \brief sets flow rates of control structure i
/// \param i [in] control structure index (assumed between 0...nControlStructures)
/// \param Q [in] initial flow rate [m3/s]
/// \param Qlast [in] flow rate [m3/s] at start of time step
/// assumes we are reading from .rvc file
//
void CReservoir::SetControlFlow(const int i, const double& Q, const double& Qlast)
{
  _aQstruct[i]=Q;
  _aQstruct_last[i]=Qlast;
}

//////////////////////////////////////////////////////////////////
/// \brief sets reservoir outflow rate as generated by management optimization scheme
/// \param Qopt [in] flow rate [m3/s] at end of time step
//
void CReservoir::SetOptimizedOutflow(const double& Qopt)
{
  _Qoptimized=Qopt;
}

//////////////////////////////////////////////////////////////////
/// \brief sets initial stage
/// \param ht [in] reservoir stage elevation
//
void  CReservoir::SetReservoirStage(const double &ht,const double &ht_last)
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
double  CReservoir::RouteWater(const double &Qin_old,
                               const double &Qin_new,
                               const CModelABC* pModel,
                               const optStruct &Options,
                               const time_struct &tt,
                               double &res_outflow,
                               res_constraint &constraint,
                               double *aQstruct) const
{
  if ((_assimilate_stage) && (!_assim_blank))
  {
    res_outflow=_Qout;
    return _stage;
  }

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

  if(_pWeirHeightTS  !=NULL){ weir_adj   =_pWeirHeightTS->  GetSampledValue(nn);}
  if(_pMaxStageTS    !=NULL){ stage_limit=_pMaxStageTS->    GetSampledValue(nn);}
  if(_pOverrideQ     !=NULL){ Qoverride  =_pOverrideQ->     GetSampledValue(nn);}
  if(_pMinStageTS    !=NULL){ min_stage  =_pMinStageTS->    GetSampledValue(nn);}
  if(_pMinStageFlowTS!=NULL){ Qminstage  =_pMinStageFlowTS->GetSampledValue(nn);}
  if(_pTargetStageTS !=NULL){ htarget    =_pTargetStageTS-> GetSampledValue(nn);}
  if(_pMaxQIncreaseTS!=NULL){ Qdelta     =_pMaxQIncreaseTS->GetSampledValue(nn);}
  if(_pMaxQDecreaseTS!=NULL){ Qdelta_dec =_pMaxQDecreaseTS->GetSampledValue(nn);}
  if(_pQminTS        !=NULL){ Qmin       =_pQminTS->        GetSampledValue(nn);}
  if(_pQmaxTS        !=NULL){ Qmax       =_pQmaxTS->        GetSampledValue(nn);}

  //Management optimization - overrides flow
  if ((Options.management_optimization) && (_Qoptimized != RAV_BLANK_DATA)) {
    Qoverride=_Qoptimized;
  }

  // Downstream flow targets
  if(_pQdownTS!=NULL)
  {
    double Qdown_targ      =_pQdownTS->GetSampledValue(nn);
    double Qdown_act       =_pQdownSB->GetOutflowRate();
    //Qshift=max(min((Qdown_targ-Qdown_act)/(_QdownRange*0.5),1.0),-1.0)*(Qdown_targ-Qdown_act);

    //smoothly shift towards target - far away moves towards range edge, in zone moves towards target
    double alpha=exp(-fabs(Qdown_targ-Qdown_act)/(_QdownRange*0.5));
    if (_QdownRange==0){alpha=0.0;}
    if(Qdown_act<Qdown_targ) {Qshift=0.8*((Qdown_targ-(1.0-alpha)*(_QdownRange*0.5))-Qdown_act);}
    else                     {Qshift=0.8*((Qdown_targ+(1.0-alpha)*(_QdownRange*0.5))-Qdown_act);}
  }

  // Downstream water demand sets minimum flow
  for(int i=0;i<_nDemands;i++) {
    if(IsInDateRange(tt.julian_day,_aDemands[i]->julian_start,_aDemands[i]->julian_end)){
      Qmin+=(_aDemands[i]->pDownSB->GetTotalWaterDemand(tt.model_time)*_aDemands[i]->percent);
      Qmin+= _aDemands[i]->pDownSB->GetEnviroMinFlow   (tt.model_time)*1.0; //assume 100% of environmental min flow must be met
    }
  }

  // ======================================================================================
  // Mass balance on reservoir over time step:
  //
  // dV(h)/dt=(Q_in(n)+Q_in(n+1))/2-(Q_out(n)+Q_out(n+1)(h))/2+Precip-ET*(A(n)+A(n+1))/2 -> solve for h_(n+1)
  //   non-linear w.r.t., h : rewritten as f(h)-gamma=0 for Newton's method solution
  //
  // ======================================================================================
  double dh        =0.001; //[m]
  double V_old     =GetVolume(_stage);
  double A_old     =GetArea(_stage);
  double h_guess   =_stage;
  int    iter      =0;
  double change    =0;
  double f,dfdh,out,out2;
  double ET      (0.0);               //[m/s]
  double precip  (0.0);               //[m3/s]
  double ext_old (0.0),ext_new (0.0); //[m3/s]
  double seep_old(0.0);
  double outflow_nat,stage_nat;       //[m],[m3/s]

  //if HRU is NULL, precip shows up as runoff from HRUs, ET calculated from open water evap from HRUs (or neglected if no such process exists)
  //Otherwise, these processes are handled as seen here
  if(_pHRU!=NULL)
  {
    ET=_pHRU->GetForcingFunctions()->OW_PET/SEC_PER_DAY/MM_PER_METER; //average for timestep, in m/s
    ET*=_pHRU->GetSurfaceProps()->lake_PET_corr;
    precip=_Precip/Options.timestep/SEC_PER_DAY; //[m3]->[m3/s]
  }
  if(_seepage_const>0) {
    seep_old=_seepage_const*(_stage-_local_GW_head); //[m3/s]
  }
  ext_new=ext_old=0.0;
  for (int ii=0; ii<_nDemandTS;ii++){
    if (Options.management_optimization){
      ext_new+=_aQdelivered[ii];
      ext_old+=_aQdelivered[ii];
    }
    else{
      ext_old=_pDemandTS[ii]->GetSampledValue(nn);
      ext_new=_pDemandTS[ii]->GetSampledValue(nn); //steady rate over time step
    }
  }

  double gamma=V_old+((Qin_old+Qin_new)-_Qout+2.0*precip-ET*A_old-seep_old-(ext_old+ext_new))/2.0*(tstep*SEC_PER_DAY);//[m3]
  if(gamma<0)
  {//reservoir dried out; no solution available. (f is always >0, so gamma must be as well)
   //only remaining filling action is via seepage, which is likely not enough, and Q_out_new can't be negative
   // string warn="CReservoir::RouteWater: basin "+to_string(_SBID)+ " dried out on " +tt.date_string;
   //WriteWarning(warn,false);

    constraint=RC_DRY_RESERVOIR;
    res_outflow=0.0;
    return _min_stage;
  }

  //double hg[RES_MAXITER],ff[RES_MAXITER],fff[RES_MAXITER];//retain for debugging
  double relax=1.0;
  do //Newton's method with discrete approximation of df/dh
  {
    out=out2=0.0;
    if     (_pDZTR==NULL) {
      out =GetWeirOutflow(h_guess,   weir_adj);//[m3/s]
      out2=GetWeirOutflow(h_guess+dh,weir_adj);//[m3/s]
      for(int i=0; i<_nControlStructures; i++) {
        out +=_pControlStructures[i]->GetOutflow(h_guess   ,_stage_last, _aQstruct_last[i],tt);
        out2+=_pControlStructures[i]->GetOutflow(h_guess+dh,_stage_last, _aQstruct_last[i],tt);
      }
    }
    else if(_pDZTR!=NULL) {
      out =GetDZTROutflow(GetVolume(h_guess   ),Qin_old,tt,Options);
      out2=GetDZTROutflow(GetVolume(h_guess+dh),Qin_old,tt,Options);
    }
    out +=ET*GetArea(h_guess   )+_seepage_const*(h_guess   -_local_GW_head);//[m3/s]
    out2+=ET*GetArea(h_guess+dh)+_seepage_const*(h_guess+dh-_local_GW_head);//[m3/s]

    f   = (GetVolume(h_guess   )+out /2.0*(tstep*SEC_PER_DAY)); //[m3]
    dfdh=((GetVolume(h_guess+dh)+out2/2.0*(tstep*SEC_PER_DAY))-f)/dh; //[m3/m]

    //hg[iter]=relax*h_guess; ff[iter]=f-gamma; fff[iter]=f;//retain for debugging

    change=-(f-gamma)/dfdh;//[m]
    if(dfdh==0) { change=1e-7; }

    if(iter>3) { relax *=0.98; }
    h_guess+=relax*change;
    iter++;
  } while((iter<RES_MAXITER) && (fabs(change/relax)>RES_TOLERANCE));

  stage_new=h_guess;

  if(iter==RES_MAXITER) {
    string warn="CReservoir::RouteWater did not converge after "+to_string(RES_MAXITER)+"  iterations for basin "+to_string(_SBID)+" on "+tt.date_string;;
    WriteWarning(warn,false);
    /*for(int i = 0; i < iter; i++) { //retain for debugging
      string warn = to_string(this->GetSubbasinID())+"["+to_string(i)+"] "+to_string(hg[i]) + " " + to_string(ff[i])+ " "+ to_string(fff[i])+ " "+to_string(gamma);WriteWarning(warn,false);
    }*/
  }

  //standard case - outflow determined through stage-discharge curve or DZTR
  //---------------------------------------------------------------------------------------------
  if(_pDZTR==NULL) {
    res_outflow=GetWeirOutflow(stage_new,weir_adj);
    for(int i=0; i<_nControlStructures; i++) {
      aQstruct[i]=_pControlStructures[i]->GetOutflow(stage_new,_stage_last,_aQstruct_last[i],tt);
      res_outflow +=aQstruct[i];
    }
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
  double w = pModel->GetGlobalParams()->GetParams()->reservoir_relax;
  if(htarget!=RAV_BLANK_DATA) {
    double V_targ=GetVolume(htarget);
    double A_targ=GetArea(htarget);
    double seep_targ = _seepage_const*(htarget-_local_GW_head);
    Qtarget = -2 * (V_targ - V_old) / (tstep*SEC_PER_DAY) + (-_Qout + (Qin_old + Qin_new) +2.0*precip - ET*(A_old + A_targ) - (seep_old+seep_targ) - (ext_old + ext_new));//[m3/s]
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
      if (Options.management_optimization)
      {
        constraint=RC_MANAGEMENT;
        res_outflow=Qoverride;
      }
      else{
        res_outflow=max(2*Qoverride-_Qout,Qminstage); //Qoverride is avg over dt, res_outflow is end of dt
      }
    }

    if((constraint==RC_MIN_STAGE) && (_minStageDominant)) {
      //In this case, min stage constraint overrides min flow constraint; nothing done
    }
    else //default override of outflow - calculate resultant volume
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
        V_new= V_old+((Qin_old+Qin_new)-(_Qout+res_outflow)+2.0*precip-ET*(A_old+A_guess)-(seep_old+seep_guess)-(ext_old+ext_new))/2.0*(tstep*SEC_PER_DAY);
        if(V_new<0) {
          V_new=0;
          constraint=RC_DRY_RESERVOIR; //drying out reservoir
          res_outflow = -2.0 * (V_new - V_old) / (tstep*SEC_PER_DAY) + (-_Qout + 2.0*precip + (Qin_old + Qin_new) - ET*(A_old + 0.0) - (ext_old + ext_new));//[m3/s] //dry it out
        }
        stage_new=InterpolateCurve(V_new,_aVolume,_aStage,_Np,false);
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
    res_outflow = -2.0 * (V_limit - V_old) / (tstep*SEC_PER_DAY) + (-_Qout + (Qin_old + Qin_new) + 2.0*precip - ET*(A_old + A_limit) - (seep_old+seep_lim) - (ext_old + ext_new));//[m3/s]
  }

  //other option - returns to stage-discharge curve once max stage exceeded
  if((Options.res_overflowmode==OVERFLOW_NATURAL) &&  (stage_nat>stage_limit)) {
    res_outflow=outflow_nat;
    stage_new  =stage_nat;
  }
  // ======================================================================================

  double total_outflow=res_outflow;
  int downID=DOESNT_EXIST;
  if (_pDownSB!=NULL){downID=_pDownSB->GetID();}
  for(int i=0; i<_nControlStructures; i++)
  {
    if(outflow_nat!=0) {
      aQstruct[i]*=(total_outflow/outflow_nat); //correct for applied reservoir rules; assumes fractions unchanged -doesn't really work for stage constraints.
    }
    if (_pControlStructures[i]->GetTargetBasinID()!=downID){
      res_outflow-=aQstruct[i]; //adjust main outflow for diversions elsewhere
    }
  }

  return stage_new;
}

//////////////////////////////////////////////////////////////////
/// \brief writes state variables to solution file
//
void CReservoir::WriteToSolutionFile (ofstream &RVC) const
{
  RVC<<"    :ResFlow, "<<_Qout<<","<<_Qout_last<<endl;
  RVC<<"    :ResStage, "<<_stage<<","<<_stage_last<<endl;
  if(_DAscale_last!=1.0) {
    RVC<<"    :ResDAscale, "<<_DAscale<<","<<_DAscale_last<<endl;
  }
  for (int i = 0; i < _nControlStructures; i++) {
    RVC<<"    :ControlFlow, "<<i<<","<<_aQstruct[i]<<","<<_aQstruct_last[i]<<endl;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief interpolates the volume from the volume-stage rating curve
/// \param ht [in] reservoir stage
/// \returns reservoir volume [m3] corresponding to stage ht
//
double     CReservoir::GetVolume(const double &ht) const
{
  return InterpolateCurve(ht,_aStage,_aVolume,_Np,true);
}
//////////////////////////////////////////////////////////////////
/// \brief interpolates the surface area from the area-stage rating curve
/// \param ht [in] reservoir stage
/// \returns reservoir surface area [m2] corresponding to stage ht
//
double     CReservoir::GetArea  (const double &ht) const
{
  return InterpolateCurve(ht,_aStage,_aArea,_Np,false);
}
//////////////////////////////////////////////////////////////////
/// \brief interpolates the discharge from the outflow-stage rating curve
/// \param ht [in] reservoir stage
/// \returns reservoir outflow [m3/s] corresponding to stage ht
/// \note assumes regular spacing between min and max stage
//
double     CReservoir::GetWeirOutflow(const double &ht, const double &adj) const
{
  double underflow=InterpolateCurve(ht,_aStage,_aQunder,_Np,false); //no adjustments
  return InterpolateCurve(ht-adj,_aStage,_aQ,_Np,false)+underflow;
}
//////////////////////////////////////////////////////////////////
/// \brief clears all time series data for re-read of .rvt file
/// \remark Called only in ensemble mode
///
/// \param &Options [in] Global model options information
//
void CReservoir::ClearTimeSeriesData(const optStruct& Options)
{
  for (int i=0; i<_nDemandTS;i++){delete _pDemandTS[i]; } delete [] _pDemandTS; _pDemandTS=NULL;
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
