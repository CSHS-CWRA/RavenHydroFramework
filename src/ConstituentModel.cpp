/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2021 the Raven Development Team
----------------------------------------------------------------
Constituent Transport/Tracer Model class
coordinates information about constituent storage
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "HydroProcessABC.h"
#include "Model.h"
#include "HeatConduction.h"
#include "Transport.h"
#include "EnergyTransport.h"

string FilenamePrepare(string filebase,const optStruct &Options); //Defined in StandardOutput.cpp
bool IsContinuousConcObs(const CTimeSeriesABC *pObs,const long SBID,const int c); //Defined in StandardOutput.cpp
void WriteNetCDFGlobalAttributes(const int out_ncid,const optStruct& Options,const string descript);
int  NetCDFAddMetadata     (const int fileid,const int time_dimid,string shortname,string longname,string units);
int  NetCDFAddMetadata2D   (const int fileid,const int time_dimid,int nbasins_dimid,string shortname,string longname,string units);
void WriteNetCDFBasinList  (const int ncid,const int varid,const CModel* pModel,const optStruct& Options);
void AddSingleValueToNetCDF(const int out_ncid,const string &label,const size_t time_index,const double &value);
//////////////////////////////////////////////////////////////////
/// \brief Implentation of the Transport constructor
/// \param pModel [in] Model object
//
CConstituentModel::CConstituentModel(CModel *pMod,CTransportModel *pTMod, string name,constit_type typ,bool is_passive,const int c)
{
  _constit_index=c;

  _pModel=pMod; 
  _pTransModel=pTMod;

  _pConstitParams=NULL;

  _pSources=NULL;
  _nSources=0;

  _nSpecFlowConcs=0;
  _pSpecFlowConcs=NULL;
  _nMassLoadingTS=0;
  _pMassLoadingTS=NULL;

  _aSourceIndices=NULL;

  _aMinHist =NULL;
  _aMlatHist=NULL;
  _aMout    =NULL;

  //Initialize constituent members
  _name=name;
  _type=typ;
  _is_passive=false;
  _initial_mass=0;
  _cumul_input =0;
  _cumul_output=0;

  //Parameter initialization
  _pConstitParams= new transport_params;
  InitializeConstitParams(_pConstitParams);

  _CONC_ncid    = -9;
  _POLLUT_ncid  = -9;
}
//////////////////////////////////////////////////////////////////
/// \brief Implentation of the Transport destructor
//
CConstituentModel::~CConstituentModel()
{
  delete[] _pConstitParams; 
  delete[] _pSpecFlowConcs; _pSpecFlowConcs=NULL;
  delete[] _pMassLoadingTS; _pMassLoadingTS=NULL;
  for(int i=0;i<_nSources;i++) { delete _pSources[i]; } delete[] _pSources;
  for(int i=0;i<_nSources;i++) { delete _aSourceIndices[i]; } delete[] _aSourceIndices;
  DeleteRoutingVars();
}
//////////////////////////////////////////////////////////////////
/// \brief returns constituent type
//
constit_type CConstituentModel::GetType() const 
{
  return _type;
}
//////////////////////////////////////////////////////////////////
/// \brief returns constituent name
//
string CConstituentModel::GetName() const
{
  return _name;
}
//////////////////////////////////////////////////////////////////
/// \brief returns constituent c in model
//
const transport_params *CConstituentModel::GetConstituentParams() const
{
  return _pConstitParams;
}

//////////////////////////////////////////////////////////////////
/// \brief returns true if constituent is passive 
//
bool  CConstituentModel::IsPassive() const {
  return _is_passive;
}
//////////////////////////////////////////////////////////////////
/// \brief returns advection correction factor for constituent c
/// being transported from storage compartment iFromWater to storage compartment iToWater
/// \param pHRU [in] pointer to HRU of interest
/// \param iFromWater [in] index of "from" water storage state variable
/// \param iToWater [in] index of "to" water storage state variable
/// \param Cs [in] RAW constituent concentration (mg/L or MJ/L or L/L, not mg/L or C or o/oo)
//
double CConstituentModel::GetAdvectionCorrection(const CHydroUnit* pHRU,const int iFromWater,const int iToWater,const double& Cs) const
{
  return 1.0;
}

//////////////////////////////////////////////////////////////////
/// \brief returns decay coefficient for constituent c
/// \param c [in] constituent index
/// \param iStorWater [in] index of water storage state variable
//
double CConstituentModel::GetDecayCoefficient(const CHydroUnit *pHRU,const int iStorWater) const
{
  return _pConstitParams->decay_coeff;
}

//////////////////////////////////////////////////////////////////
/// \brief returns advection correction factor for nutrient constituent
/// being transported from storage compartment iFromWater to storage compartment iToWater
/// \param pHRU [in] pointer to HRU of interest
/// \param iFromWater [in] index of "from" water storage state variable
/// \param iToWater [in] index of "to" water storage state variable
/// \param Cs [in] RAW constituent concentration (mg/L)
//
double CNutrientModel::GetAdvectionCorrection(const CHydroUnit *pHRU,const int iFromWater,const int iToWater, const double &Cs) const
{
  sv_type fromType,toType;
  fromType=_pModel->GetStateVarType(iFromWater);
  toType  =_pModel->GetStateVarType(iToWater);

  if(IsPassive()) { return 1.0; }
  if(fromType==SOIL)
  {
    int    m=_pModel->GetStateVarLayer(iFromWater);
#ifdef _STRICTCHECK_
    ExitGracefullyIf(m==DOESNT_EXIST,"GetAdvectionCorrection:invalid iFromWater",RUNTIME_ERR);
#endif

    if(toType!=ATMOSPHERE) { return 1.0/pHRU->GetSoilProps(m)->retardation[_constit_index]; }
    else                   { return pHRU->GetVegetationProps()->uptake_moderator[_constit_index]; }
  }
  return 1.0;
}

//////////////////////////////////////////////////////////////////
/// \brief returns decay coefficient for constituent c
/// \param c [in] constituent index
/// \param iStorWater [in] index of water storage state variable
//
double CNutrientModel::GetDecayCoefficient(const CHydroUnit *pHRU,const int iStorWater) const
{
  sv_type storType;
  double decay_coeff = CConstituentModel::GetDecayCoefficient(pHRU,iStorWater);

  storType=_pModel->GetStateVarType(iStorWater);

  //add special decay_coefficients from other processes
  if(storType == SOIL)
  {
    int m = _pModel->GetStateVarLayer(iStorWater);
#ifdef _STRICTCHECK_
    ExitGracefullyIf(m==DOESNT_EXIST,"GetDecayCoefficient:invalid _iFromWater",RUNTIME_ERR);
    ExitGracefullyIf((c<0) || (c>MAX_CONSTITUENTS),"GetDecayCoefficient:invalid constit",RUNTIME_ERR);
#endif

    decay_coeff += pHRU->GetSoilProps(m)->mineraliz_rate[_constit_index];
    decay_coeff += pHRU->GetSoilProps(m)->loss_rate[_constit_index];//e.g., denitrification
  }
  return decay_coeff;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns total rivulet storage distributed over watershed [mg] or [MJ]
/// \return Total rivulet storage distributed over watershed  [mg] or [MJ]
//
double CConstituentModel::GetTotalRivuletConstituentStorage() const
{
  double sum(0);
  for(int p=0;p<_pModel->GetNumSubBasins();p++)
  {
    sum+=_rivulet_storage[p]; //[mg] or [MJ]
  }
  return sum; //[mg] or [MJ]
}
//////////////////////////////////////////////////////////////////
/// \brief Returns total channel storage distributed over watershed [mg] or [MJ]
/// \return Total channel storage distributed over watershed  [mg] or [MJ]
//
double CConstituentModel::GetTotalChannelConstituentStorage() const
{
  double sum(0);
  for(int p=0;p<_pModel->GetNumSubBasins();p++)
  {
    sum+=_channel_storage[p]; //[mg] or [MJ]
  }
  return sum; //[mg] or [MJ]
}
//////////////////////////////////////////////////////////////////
/// \brief returns net mass lost while transporting along reach p [mg]
/// \param p [in] subbasin index
/// \returns net mass lost while transporting along reach p [mg]
//
double CConstituentModel::GetNetReachLosses(const int p) const
{
  return 0.0;//default -assumes conservative transport in streams (more interesting for child classes)
}
//////////////////////////////////////////////////////////////////
/// \brief Test for whether a dirichlet condition applies to a certain compartment, place, and time
/// \note called within solver to track mass balance
/// \returns true if dirichlet source applies
/// \returns Cs, source concentration [mg/L] or [C] for enthalpy
/// \param i_stor [in] storage index of water compartment
/// \param c [in] constituent index
/// \param k [in] global HRU index
/// \param tt [in] current time structure
/// \param Cs [out] Dirichlet source concentration
//
bool  CConstituentModel::IsDirichlet(const int i_stor,const int k,const time_struct &tt,double &Cs) const
{
  Cs=0.0;
  int i_source=_aSourceIndices[i_stor][k];
  if(i_source==DOESNT_EXIST)          { return false; }
  if(!_pSources[i_source]->dirichlet) { return false; }
  Cs = _pSources[i_source]->concentration;

  if(Cs != DOESNT_EXIST) { return true; }
  else {//time series
    Cs = _pSources[i_source]->pTS->GetValue(tt.model_time);
    return true;
  } 
}
//////////////////////////////////////////////////////////////////
/// \brief returns specified mass flux for given constitutent and water storage unit at time tt
/// \returns source flux in mg/m2/d
/// \param i_stor [in] storage index of water compartment
/// \param c [in] constituent index
/// \param k [in] global HRU index
/// \param tt [in] current time structure
//
double  CConstituentModel::GetSpecifiedMassFlux(const int i_stor,const int k,const time_struct &tt) const
{
  int i_source=_aSourceIndices[i_stor][k];
  if(i_source == DOESNT_EXIST) { return 0.0; }
  if(_pSources[i_source]->dirichlet) { return 0.0; }

  double flux;
  flux=_pSources[i_source]->concentration; //'concentration' stores mass flux [mg/m2/d] if this is a flux-source
  if(flux == DOESNT_EXIST) {//get from time series
    flux=_pSources[i_source]->pTS->GetValue(tt.model_time);
  }
  return flux;
}


//////////////////////////////////////////////////////////////////
/// \brief adds dirichlet source
/// \param const_name [in] constituent name
/// \param i_stor [in] global index of water storage state variable
/// \param kk [in] HRU group index (or -1 if this applies to all HRUs)
/// \param Cs [in] Dirichlet source concentration [mg/L] or temperature [C]
//
void   CConstituentModel::AddDirichletCompartment(const int i_stor,const int kk,const double Cs)
{
  static constit_source *pLast;
  constit_source *pSource=new constit_source();

  pSource->constit_index=_constit_index; // REFACTOR - this is really no longer needed
  pSource->concentration=Cs;
  pSource->flux         =0.0;
  pSource->dirichlet    =true;
  pSource->i_stor       =i_stor;
  pSource->kk           =kk;
  pSource->pTS          =NULL;

  ExitGracefullyIf(pSource->i_stor       ==DOESNT_EXIST,"AddDirichletCompartment: invalid storage compartment index",BAD_DATA_WARN);

  if(!DynArrayAppend((void**&)(_pSources),(void*)(pSource),_nSources)) {
    ExitGracefully("CConstituentModel::AddDirichletCompartment: adding NULL source",BAD_DATA);
  }

  pLast=pSource;//so source is not deleted upon leaving this routine
}
//////////////////////////////////////////////////////////////////
/// \brief adds dirichlet source time series
/// \param const_name [in] constituent name
/// \param i_stor [in] global index of water storage state variable
/// \param kk [in] HRU group index (or -1 if this applies to all HRUs)
/// \param pTS [in] Time series of Dirichlet source concentration [mg/l] or temperature [degC]
//
void   CConstituentModel::AddDirichletTimeSeries(const int i_stor,const int kk,const CTimeSeries *pTS)
{
  static constit_source *pLast;
  constit_source *pSource=new constit_source();

  pSource->constit_index=_constit_index;
  pSource->dirichlet    =true;
  pSource->concentration=DOESNT_EXIST;
  pSource->flux         =0.0;
  pSource->i_stor       =i_stor;
  pSource->kk           =kk;
  pSource->pTS          =pTS;

  ExitGracefullyIf(pSource->i_stor       ==DOESNT_EXIST,"AddDirichletTimeSeries: invalid storage compartment index",BAD_DATA_WARN);

  if(!DynArrayAppend((void**&)(_pSources),(void*)(pSource),_nSources)) {
    ExitGracefully("CConstituentModel::AddDirichletCompartment: adding NULL source",BAD_DATA);
  }

  pLast=pSource; //so source is not deleted upon leaving this routine
}

//////////////////////////////////////////////////////////////////
/// \brief adds influx (Neumann) source
/// \param const_name [in] constituent name
/// \param i_stor [in] global index of water storage state variable
/// \param kk [in] HRU group index (or DOESNT_EXIST if this applies to all HRUs)
/// \param flux [in] specified fixed mass influx rate [mg/m2/d] or energy influx rate [MJ/m2/d]
//
void   CConstituentModel::AddInfluxSource(const int i_stor,const int kk,const double flux)
{
  static constit_source *pLast;
  constit_source *pSource=new constit_source();

  pSource->constit_index=_constit_index;
  pSource->dirichlet    =false;
  pSource->concentration=0.0;
  pSource->flux         =flux;
  pSource->i_stor       =i_stor;
  pSource->kk           =kk;
  pSource->pTS          =NULL;

  ExitGracefullyIf(pSource->i_stor       ==DOESNT_EXIST,"AddInfluxSource: invalid storage compartment index",BAD_DATA_WARN);

  if(!DynArrayAppend((void**&)(_pSources),(void*)(pSource),_nSources)) {
    ExitGracefully("CConstituentModel::AddInfluxSource: adding NULL source",BAD_DATA);
  }

  pLast=pSource;//so source is not deleted upon leaving this routine
}
//////////////////////////////////////////////////////////////////
/// \brief adds influx (Neumann) source time series
/// \param const_name [in] constituent name
/// \param i_stor [in] global index of water storage state variable
/// \param kk [in] HRU group index (or DOESNT_EXIST if this applies to all HRUs)
/// \param pTS [in] Time series of Neumann source flux rate [mg/m2/d] or [MJ/m2/d]
//
void   CConstituentModel::AddInfluxTimeSeries(const int i_stor,const int kk,const CTimeSeries *pTS)
{
  static constit_source *pLast;
  constit_source *pSource=new constit_source();

  pSource->constit_index=_constit_index;
  pSource->dirichlet    =false;
  pSource->concentration=0.0;
  pSource->flux         =DOESNT_EXIST;
  pSource->i_stor       =i_stor;
  pSource->kk           =kk;
  pSource->pTS          =pTS;

  ExitGracefullyIf(pSource->i_stor       ==DOESNT_EXIST,"AddDirichletTimeSeries: invalid storage compartment index",BAD_DATA_WARN);

  if(!DynArrayAppend((void**&)(_pSources),(void*)(pSource),_nSources)) {
    ExitGracefully("CConstituentModel::AddInfluxTimeSeries: adding NULL source",BAD_DATA);
  }

  pLast=pSource; //so source is not deleted upon leaving this routine
}
//////////////////////////////////////////////////////////////////
/// \brief adds specified in-reach flow concentration/temperature time series
/// \param pTS [in] Time series of concentrations [mg/L] or temperatures [C]
/// assumes SBID and c tags are already applied to time series
//
void   CConstituentModel::AddInflowConcTimeSeries(const CTimeSeries *pTS) {

  if(!DynArrayAppend((void**&)(_pSpecFlowConcs),(void*)(pTS),_nSpecFlowConcs)) {
    ExitGracefully("CConstituentModel::AddSpecifConcTimeSeries: adding NULL source",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief adds specified in-reach mass loading time series
/// \param pTS [in] Time series of mass loading [mg/d] or heat loading ? [MJ/d]
/// assumes SBID and c tags are already applied to time series
//
void   CConstituentModel::AddMassLoadingTimeSeries(const CTimeSeries* pTS) 
{
  if(!DynArrayAppend((void**&)(_pMassLoadingTS),(void*)(pTS),_nMassLoadingTS)) {
    ExitGracefully("CConstituentModel::AddMassLoadingTimeSeries: adding NULL source",BAD_DATA);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Set subbasin initial channel reach mass 
/// \param p [in] global subbasin index
/// \param aMout [in] total rivulet mass [mg] or [MJ]
//
void   CConstituentModel::SetChannelMass(const int p,const double mass)
{
  _channel_storage[p]=mass;
}
//////////////////////////////////////////////////////////////////
/// \brief Set subbasin initial rivulet mass 
/// \param p [in] global subbasin index
/// \param aMout [in] total rivulet mass [mg] or [MJ]
//
void   CConstituentModel::SetRivuletMass(const int p,const double mass)
{
  _rivulet_storage[p]=mass;
}
//////////////////////////////////////////////////////////////////
/// \brief Set subbasin mass outflow conditions
/// \param p [in] global subbasin index
/// \param nsegs [in] number of segments in subbasin reach
/// \param aMout [in] array of mass/energy flows [mg/d] or [MJ/d]
/// \param MoutLast [in] mass flow from end of reach from last timestep [mg/d] or [MJ/d]
//
void   CConstituentModel::SetMoutArray(const int p,const int nsegs,const double *aMout,const double MoutLast)
{
  if(nsegs!=_pModel->GetSubBasin(p)->GetNumSegments()) {
    WriteWarning("Number of reach segments in state file and input file are inconsistent. Unable to read in-reach mass flow initial conditions",false);
    return;
  }
  for(int i=0;i<nsegs;i++) { _aMout[p][i]=aMout[i]; }
  _aMout_last[p]=MoutLast;
}
//////////////////////////////////////////////////////////////////
/// \brief Set subbasin mass lateral inflow history conditions
/// \param p [in] global subbasin index
/// \param histsize [in] size of reservoir mass lateral inflow history
/// \param aMlat [in] array of mass lateral inflows [mg/d] or [MJ/d]
/// \param MlatLast [in] mass lateral inflow from last timestep [mg/d] or [MJ/d]
//
void   CConstituentModel::SetMlatHist(const int p,const int histsize,const double *aMlat,const double MlatLast)
{
  if(histsize!=_pModel->GetSubBasin(p)->GetLatHistorySize()) {
    WriteWarning("Size of lateral inflow history in state file and input file are inconsistent. Unable to read in-reach mass lateral flow initial conditions",false);
    return;
  }
  for(int i=0;i<histsize;i++) { _aMlatHist[p][i]=aMlat[i]; }
  _aMlat_last[p]=MlatLast;
}
//////////////////////////////////////////////////////////////////
/// \brief Set subbasin mass inflow history conditions
/// \param p [in] global subbasin index
/// \param histsize [in] size of reservoir mass inflow history
/// \param aMin [in] array of mass inflows [mg/d] or [MJ/d]
//
void   CConstituentModel::SetMinHist(const int p,const int histsize,const double *aMin)
{
  if(histsize!=_pModel->GetSubBasin(p)->GetLatHistorySize()) {
    WriteWarning("Size of mass inflow history in state file and input file are inconsistent. Unable to read in-reach mass inflow initial conditions",false);
    return;
  }
  for(int i=0;i<histsize;i++) { _aMinHist[p][i]=aMin[i]; }
}
//////////////////////////////////////////////////////////////////
/// \brief Set reservoir initial mass conditions
/// \param p [in] global subbasin index
/// \param Mout [in] reservoir mass  [mg]
/// \param MoutLast [in] reservoir mass from previous timestep [mg]
//
void   CConstituentModel::SetInitialReservoirMass(const int p,const double res_mass,const double res_mass_last)
{
  _aMres     [p]=res_mass;
  _aMres_last[p]=res_mass_last;
}
//////////////////////////////////////////////////////////////////
/// \brief Set reservoir rainfall input rate for this time step (called in Solver)
/// \param p [in] global subbasin index
/// \param precip_load_rate [in] reservoir mass/enthalpy input rate  [mg/d] or [MJ/d]
//
void   CConstituentModel::SetReservoirPrecipLoad(const int p,const double precip_load_rate) 
{
  _aMresRain[p]=precip_load_rate;
}
//////////////////////////////////////////////////////////////////
/// \brief Set reservoir initial mass flow conditions
/// \param p [in] global subbasin index
/// \param Mout [in] reservoir mass outflow [mg/d]
/// \param MoutLast [in] reservoir mass outflow from previous timestep [mg/d]
//
void   CConstituentModel::SetReservoirMassOutflow(const int p,const double Mout,const double MoutLast)
{
  _aMout_res     [p]=Mout;
  _aMout_res_last[p]=MoutLast;
}
//////////////////////////////////////////////////////////////////
/// \brief Set transport parameter value for specified constituent
//
/*void   CConstituentModel::SetGlobalParameter(const string const_name,const string param_name,const double &value,bool noisy)
{
  int constit_ind=GetConstituentIndex(const_name);
  if(constit_ind==DOESNT_EXIST) {
    WriteWarning("CConstituentModel::SetGlobalParameter: Unrecognized constituent name",noisy);
  }

  string ustr=StringToUppercase(param_name);
  if(ustr=="DECAY_COEFF") {
    _pConstitParams[constit_ind]->decay_coeff=value;
    ExitGracefullyIf(value<0.0," CConstituentModel::SetGlobalParameter: decay coefficient cannot be negative",BAD_DATA_WARN);
  }
  //else if (ustr=="OTHER PARAM"){
  //...
  //}
  else {
    WriteWarning("CConstituentModel::SetGlobalParameter: Unrecognized parameter name",noisy);
  }
}
*/
//////////////////////////////////////////////////////////////////
/// \brief Initialization of all transport variables
/// \note determines initial conditions for all constituents, initializes routing variables
/// \note  called after all constituents have been added by CModel::Initialize
//
void CConstituentModel::Initialize(const optStruct &Options)
{
  int i,m,i_stor;
  double watstor;
  double area=_pModel->GetWatershedArea();

  // Calculate initial mass and zero out zero volume concentrations
  //--------------------------------------------------------------------
  _cumul_input=0;
  _cumul_output=0;
  _initial_mass=0;
  
  for(int ii=0;ii<_pTransModel->GetNumWaterCompartments();ii++)
  {
    m=(_constit_index*_pTransModel->GetNumWaterCompartments())+ii;
    i=_pModel->GetStateVarIndex(CONSTITUENT,m);

    // zero out initial mass in all storage units with zero volume (to avoid absurdly high concentrations/temperatures)
    i_stor=_pTransModel->GetStorWaterIndex(ii);
    for(int k = 0; k < _pModel->GetNumHRUs(); k++)
    {
      watstor = _pModel->GetHydroUnit(k)->GetStateVarValue(i_stor);
      if(watstor<1e-9) { _pModel->GetHydroUnit(k)->SetStateVarValue(i,0.0); }
    }

    //update initial model mass (initialized from .rvc file)
    _initial_mass+=_pModel->GetAvgStateVar(i)*(area*M2_PER_KM2); //mg 
  }

  /// \todo [funct]: calculate initial mass from Dirichlet cells?

  // populate array of source indices
  //--------------------------------------------------------------------
  _aSourceIndices=NULL;
  _aSourceIndices=new int *[_pModel->GetNumStateVars()];
  ExitGracefullyIf(_aSourceIndices==NULL,"CTransport::Initialize",OUT_OF_MEMORY);
  for(i_stor=0;i_stor<_pModel->GetNumStateVars();i_stor++) {
    _aSourceIndices[i_stor]=NULL;
    _aSourceIndices[i_stor]=new int[_pModel->GetNumHRUs()];
    ExitGracefullyIf(_aSourceIndices==NULL,"CTransport::Initialize(2)",OUT_OF_MEMORY);
    for(int k=0; k<_pModel->GetNumHRUs();k++) {
      _aSourceIndices[i_stor][k]=DOESNT_EXIST; //default assumption - no source in compartment
    }
    for(i=0;i<_nSources;i++) {
      if(_pSources[i]->i_stor==i_stor)
      {
        int kk=_pSources[i]->kk;
        if(kk==DOESNT_EXIST) {
          for(int k=0; k<_pModel->GetNumHRUs();k++) {
            _aSourceIndices[i_stor][k]=i; 
          }
        }
        else {
          for(int k=0; k<_pModel->GetNumHRUs();k++) {
            if(_pModel->GetHRUGroup(kk)->IsInGroup(k)) {
              _aSourceIndices[i_stor][k]=i;
            }
          }
        }
      }
    }
  }

  // Initialize routing variables
  //--------------------------------------------------------------------
  InitializeRoutingVars();

  // Initialize time series
  //--------------------------------------------------------------------  
  for(i=0; i<_nMassLoadingTS; i++) {
    _pMassLoadingTS[i]->Initialize(Options.julian_start_day,Options.julian_start_year,Options.duration,Options.timestep,true,Options.calendar);
  }
  for(i=0; i<_nSpecFlowConcs; i++) {
    _pSpecFlowConcs[i]->Initialize(Options.julian_start_day,Options.julian_start_year,Options.duration,Options.timestep,true,Options.calendar);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Initialization of transport parameter structure
/// \param pP [in] valid pointer to transport_params structure pP
//
void CConstituentModel::InitializeConstitParams(transport_params *pP)
{
  pP->decay_coeff = 0.0;
}
//////////////////////////////////////////////////////////////////
/// \brief Preparation of all transport variables
/// \note  called after all waterbearing state variables & processes have been generated, but before constiuents created
/// \param Options [in] Global model options structure
//
void CConstituentModel::Prepare(const optStruct &Options)
{

}
//////////////////////////////////////////////////////////////////
/// \brief Calculate concentration for reporting (converts to mg/L)
/// \param M [in] mass [mg/m2] 
/// \param V [in] volume [mm]
//
double CConstituentModel::CalculateConcentration(const double &M, const double &V) const
{
  if(fabs(V)>1e-6) {
      return M/V*MM_PER_METER/LITER_PER_M3; //[mg/mm/m2]->[mg/L]
  }
  return 0.0;// JRC: should this default to zero? or NA?
}

//////////////////////////////////////////////////////////////////
/// \brief Increment cumulative input of mass to watershed.
/// \note called within solver to track mass balance
/// \param Options [in] Global model options structure
/// \param tt [in] current model time structure
//
void CConstituentModel::IncrementCumulInput(const optStruct &Options,const time_struct &tt)
{
  _cumul_input+=GetMassAddedFromInflowSources(tt.model_time,Options.timestep);

  //Most of this is handled using the CONSTITUENT_SRC storage compartment
}
//////////////////////////////////////////////////////////////////
/// \brief Increment cumulative output of mass from watershed exiting via stream channel
/// \note called within solver to track mass balance
/// \param Options [in] Global model options structure
//
void CConstituentModel::IncrementCumulOutput(const optStruct &Options)
{
  for(int p=0;p<_pModel->GetNumSubBasins();p++)
  {
    if(_pModel->GetSubBasin(p)->GetDownstreamID()==DOESNT_EXIST)//outlet does not drain into another subbasin
    {
      _cumul_output+=GetIntegratedMassOutflow(p,Options.timestep);
    }
    _cumul_output+=GetNetReachLosses(p);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief calculates mass outflow using convoluition
/// \details This parent version is simple, child versions may calculate decay, etc. in channel
/// \param aRouteHydro [in] array of routing hydrograph [size: nMinHist]
/// \param aMout_new [out] - array of mass outflows for all segments (but effectively assumes nSegments==1)
//
void  CConstituentModel::ApplyConvolutionRouting(const int p,const double *aRouteHydro,const double *aQinHist,const double *aMinHist,const int nSegments,const int nMinHist,const double &tstep,double *aMout_new) const
{
  aMout_new[nSegments-1]=0.0;
  for(int n=0;n<nMinHist;n++) {
    aMout_new[nSegments-1]+=aRouteHydro[n]*aMinHist[n];
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Write transport output file headers
/// \details Called prior to simulation (but after initialization) from CModel::Initialize()
/// \param &Options [in] Global model options information
//
void CConstituentModel::WriteOutputFileHeaders(const optStruct &Options) 
{
  string filename;

  int iCumPrecip=_pModel->GetStateVarIndex(ATMOS_PRECIP);

  string kg,mgL,kgd;

  //units names
  kg="[kg]"; kgd="[kg/d]"; mgL="[mg/l]";
  if(Options.write_constitmass) {
    kg="[mg/m2]"; kgd="[mg/m2/d]"; mgL="[mg/m2]";
  }
  if(_type==TRACER) {
    kg="[-]"; kgd="[-]"; mgL="[-]";
  }
  else if(_type==ENTHALPY) {
    kg="[MJ]"; kgd="[MJ/d]"; mgL="["+to_string(DEG_SYMBOL)+"C]";
    if(Options.write_constitmass) {
      kg="[MJ/m2]"; kgd="[MJ/m2/d]"; mgL="[MJ/m2]";
    }
  }
  else if(_type==ISOTOPE) {
    kg="[kg-L]"; kgd="[kg-L/d]"; mgL="[o/oo]";
    if(Options.write_constitmass) {
      kg="[1/m2]"; kgd="[1/m2/d]"; mgL="[1/m2]"; //very hard to interpret these units
    }
  }

  //Concentrations or Temperatures file
  //--------------------------------------------------------------------
  if(Options.write_constitmass) {
    if(_type!=ENTHALPY) { filename=_name+"Mass.csv"; }
    else                { filename="Enthalpy.csv"; }
  }
  else {
    if(_type!=ENTHALPY) { filename=_name+"Concentrations.csv";}
    else                { filename="Temperatures.csv"; }
  }
  filename=FilenamePrepare(filename,Options);

  _OUTPUT.open(filename.c_str());
  if(_OUTPUT.fail()) {
    ExitGracefully(("CConstituentModel::WriteOutputFileHeaders: Unable to open output file "+filename+" for writing.").c_str(),FILE_OPEN_ERR);
  }

  //    Header content ---------------------------
  if(_type!=ENTHALPY) {
    _OUTPUT<<"time[d],date,hour,influx"<<kgd<<",Channel Storage"<<kg<<",Rivulet Storage"<<kg;
  }
  else {
    _OUTPUT<<"time[d],date,hour,air temp."<<mgL<<",influx"<<kgd<<",Channel Storage"<<kg<<",Rivulet Storage"<<kg;
  }

  for(int ii=0;ii<_pTransModel->GetNumWaterCompartments();ii++)
  {
    if(_pTransModel->GetStorWaterIndex(ii)!=iCumPrecip)
    {
      _OUTPUT<<","<<
        CStateVariable::GetStateVarLongName(_pModel->GetStateVarType(_pTransModel->GetStorWaterIndex(ii)),_pModel->GetStateVarLayer(_pTransModel->GetStorWaterIndex(ii)))<<" "<<mgL;
    }
  }

  if(_type!=ENTHALPY) {
    _OUTPUT<<", Total Mass "<<kg<<", Cum. Loading "<<kg<<", Cum. Mass Lost "<<kg<<", MB Error "<<kg;
    _OUTPUT<<", atmos "<<kg<<", sink"<<kg<<endl;//TMP DEBUG
  }
  else {
    _OUTPUT<<", Total Energy "<<kg<<", Cum. Loading "<<kg<<", Cum. Energy Lost "<<kg<<", EB Error "<<kg;
    _OUTPUT<<", sink "<<kg<<", source "<<kg<<", atmos "<<kg<<", latent "<<kg;//TMP DEBUG
    _OUTPUT<<", sensible (r) "<<kg<<", latent (r) "<<kg<<", GW mixing (r) "<<kg<<", radiant (r) "<<kg<<", friction (r)"<<kg<< endl;//TMP DEBUG
  }

  //Pollutograph / stream temperatures file
  //--------------------------------------------------------------------
  if(_type!=ENTHALPY) {filename=_name+"Pollutographs.csv";}
  else                {filename="StreamTemperatures.csv"; }

  filename=FilenamePrepare(filename,Options);
  _POLLUT.open(filename.c_str());
  if(_POLLUT.fail()) {
    ExitGracefully(("CConstituentModel::WriteOutputFileHeaders: Unable to open output file "+filename+" for writing.").c_str(),FILE_OPEN_ERR);
  }

  _POLLUT<<"time[d],date,hour";
  if(_type==ENTHALPY) {
    _POLLUT<<",air temp.";
  }
  const CSubBasin *pBasin;
  for(int p=0;p<_pModel->GetNumSubBasins();p++) {
    pBasin=_pModel->GetSubBasin(p);
    if(pBasin->IsGauged() && (pBasin->IsEnabled())) {
      string name;
      if(pBasin->GetName()=="")   { _POLLUT<<",ID="<<pBasin->GetID()  <<" "<<mgL; }
      else                        { _POLLUT<<","   <<pBasin->GetName()<<" "<<mgL; }
      if(_type==ENTHALPY) {
        if(pBasin->GetName()=="") { _POLLUT<<",ID="<<pBasin->GetID()  <<" pct froz."; }
        else                      { _POLLUT<<","   <<pBasin->GetName()<<" pct froz."; }
      }
      for(int i = 0; i < _pModel->GetNumObservedTS(); i++) {
        if(IsContinuousConcObs(_pModel->GetObservedTS(i),pBasin->GetID(),_constit_index))
        {
          if(pBasin->GetName()=="") { _POLLUT<<",ID="<<pBasin->GetID()  <<" (observed) "<<mgL; }
          else                      { _POLLUT<<","   <<pBasin->GetName()<<" (observed) "<<mgL; }
        }
      }
    }
  }
  _POLLUT<<endl;
}

//////////////////////////////////////////////////////////////////
/// \brief Write transport output file headers in .tb0 format
/// \details Called prior to simulation (but after initialization) from CModel::Initialize()
/// \param &Options [in] Global model options information
//
void CConstituentModel::WriteEnsimOutputFileHeaders(const optStruct &Options)
{
  string filename;

  int iCumPrecip=_pModel->GetStateVarIndex(ATMOS_PRECIP);

  string kg,mgL,kgd;

  int i;
  time_struct tt,tt2;
  JulianConvert(0.0,Options.julian_start_day,Options.julian_start_year,Options.calendar,tt);
  JulianConvert(Options.timestep,Options.julian_start_day,Options.julian_start_year,Options.calendar,tt2);//end of the timestep

  //units names
  kg="kg"; kgd="kg/d"; mgL="mg/l";
  //kg="mg/m2"; kgd="mg/m2/d"; mgL="mg/m2";//TMP DEBUG OUTPUT OVERRIDE
  if(_type==TRACER) {
    kg="none"; kgd="none"; mgL="none";
  }

  //Concentrations file
  //--------------------------------------------------------------------
  filename=_name+"Concentrations.tb0";
  filename=FilenamePrepare(filename,Options);

  _OUTPUT.open(filename.c_str());
  if(_OUTPUT.fail()) {
    ExitGracefully(("CTransportModel::WriteEnsimOutputFileHeaders: Unable to open output file "+filename+" for writing.").c_str(),FILE_OPEN_ERR);
  }
  _OUTPUT<<"#########################################################################"<<endl;
  _OUTPUT<<":FileType tb0 ASCII EnSim 1.0"<<endl;
  _OUTPUT<<"#"<<endl;
  _OUTPUT<<":Application   Raven"<<endl;
  if(!Options.benchmarking) {
    _OUTPUT<<":Version       "<<Options.version<<endl;
    _OUTPUT<<":CreationDate  "<<GetCurrentMachineTime()<<endl;
  }
  _OUTPUT<<"#"<<endl;
  _OUTPUT<<"#------------------------------------------------------------------------"<<endl;
  _OUTPUT<<"#"<<endl;
  _OUTPUT<<":RunName       "<<Options.run_name<<endl;
  _OUTPUT<<":Format        Instantaneous" << endl;
  _OUTPUT<<"#"<<endl;

  if(Options.suppressICs) {
    _OUTPUT<< ":StartTime " << tt2.date_string << " " << DecDaysToHours(tt2.julian_day) << endl;
  }
  else {
    _OUTPUT<< ":StartTime " << tt.date_string << " " << DecDaysToHours(tt.julian_day) << endl;
  }

  if(Options.timestep!=1.0) { _OUTPUT<<":DeltaT " <<DecDaysToHours(Options.timestep)<<endl; }
  else { _OUTPUT<<":DeltaT 24:00:00.00"  <<endl; }
  _OUTPUT<<"#"<<endl;

  _OUTPUT<<":ColumnMetaData"<<endl;
  _OUTPUT<<"  :ColumnName influx \"Channel storage\" \"Rivulet storage\"";
  for(i=0;i<_pModel->GetNumStateVars();i++) {
    if((CStateVariable::IsWaterStorage(_pModel->GetStateVarType(i))) && (i!=iCumPrecip)) {
      _OUTPUT<<" \""<<CStateVariable::GetStateVarLongName(_pModel->GetStateVarType(i),_pModel->GetStateVarLayer(i))<<"\"";
    }
  }
  _OUTPUT<<" \"Total Mass\" \"Cum. Loading\" \"Cum. Mass lost\" \"MB error\""<<endl;

  _OUTPUT<<"  :ColumnUnits "<<kgd<<" "<<kg <<" "<< kg;
  for(i=0;i<_pModel->GetNumStateVars();i++) {
    if((CStateVariable::IsWaterStorage(_pModel->GetStateVarType(i))) && (i!=iCumPrecip)) {
      _OUTPUT<<" mgL";
    }
  }
  _OUTPUT<<" "<<kg<<" "<<kg<<" "<<kg<<" "<<kg<<endl;

  _OUTPUT<<"  :ColumnType float float float";
  for(i=0;i<_pModel->GetNumStateVars();i++) {
    if((CStateVariable::IsWaterStorage(_pModel->GetStateVarType(i))) && (i!=iCumPrecip)) {
      _OUTPUT<<" float";
    }
  }
  _OUTPUT<<" float float float float"<<endl;

  _OUTPUT<<"  :ColumnFormat -1 0 0";
  for(i=0;i<_pModel->GetNumStateVars();i++) {
    if((CStateVariable::IsWaterStorage(_pModel->GetStateVarType(i))) && (i!=iCumPrecip)) {
      _OUTPUT<<" 0";
    }
  }
  _OUTPUT<<" 0 0 0 0"<<endl;

  _OUTPUT<<":EndColumnMetaData"<<endl;
  _OUTPUT<<":EndHeader"<<endl;

  //Pollutograph file
  //--------------------------------------------------------------------
  filename=_name+"Pollutographs.tb0";
  filename=FilenamePrepare(filename,Options);

  _POLLUT.open(filename.c_str());
  if(_OUTPUT.fail()) {
    ExitGracefully(("CTransportModel::WriteEnsimOutputFileHeaders: Unable to open output file "+filename+" for writing.").c_str(),FILE_OPEN_ERR);
  }
  _POLLUT<<"#########################################################################"<<endl;
  _POLLUT<<":FileType tb0 ASCII EnSim 1.0"<<endl;
  _POLLUT<<"#"<<endl;
  _POLLUT<<":Application   Raven"<<endl;
  if(!Options.benchmarking) {
    _POLLUT<<":Version       "<<Options.version<<endl;
    _POLLUT<<":CreationDate  "<<GetCurrentMachineTime()<<endl;
  }
  _POLLUT<<"#"<<endl;
  _POLLUT<<"#------------------------------------------------------------------------"<<endl;
  _POLLUT<<"#"<<endl;
  _POLLUT<<":RunName       "<<Options.run_name<<endl;
  _POLLUT<<":Format        Instantaneous" << endl;
  _POLLUT<<"#"<<endl;

  if(Options.suppressICs) {
    _POLLUT<< ":StartTime " << tt2.date_string << " " << DecDaysToHours(tt2.julian_day) << endl;
  }
  else {
    _POLLUT<< ":StartTime " << tt.date_string << " " << DecDaysToHours(tt.julian_day) << endl;
  }

  if(Options.timestep!=1.0) { _POLLUT<<":DeltaT " <<DecDaysToHours(Options.timestep)<<endl; }
  else { _POLLUT<<":DeltaT 24:00:00.00"  <<endl; }
  _POLLUT<<"#"<<endl;

  const CSubBasin *pBasin;

  _POLLUT<<":ColumnMetaData"<<endl;
  _POLLUT<<"  :ColumnName";
  for(int p=0;p<_pModel->GetNumSubBasins();p++) {
    pBasin=_pModel->GetSubBasin(p);
    if(pBasin->IsGauged() && (pBasin->IsEnabled())) {
      if(pBasin->GetName()=="") { _POLLUT<<" ID="<<pBasin->GetID(); }
      else { _POLLUT<<" "   <<pBasin->GetName(); }
    }
  }
  _POLLUT<<endl;

  _POLLUT<<"  :ColumnUnits";
  for(int p=0;p<_pModel->GetNumSubBasins();p++) {
    pBasin=_pModel->GetSubBasin(p);
    if(pBasin->IsGauged() && (pBasin->IsEnabled())) { _POLLUT<<" "<<mgL; }
  }
  _POLLUT<<endl;

  _POLLUT<<"  :ColumnType";
  for(int p=0;p<_pModel->GetNumSubBasins();p++) {
    pBasin=_pModel->GetSubBasin(p);
    if(pBasin->IsGauged() && (pBasin->IsEnabled())) { _POLLUT<<" float"; }
  }
  _POLLUT<<endl;

  _POLLUT << "  :ColumnFormat";
  for(int p=0;p<_pModel->GetNumSubBasins();p++) {
    pBasin=_pModel->GetSubBasin(p);
    if(pBasin->IsGauged() && (pBasin->IsEnabled())) { _POLLUT<<" 0"; }
  }
  _POLLUT<< endl;

  _POLLUT<<":EndColumnMetaData"<<endl;
  _POLLUT<<":EndHeader"<<endl;
  
}

//////////////////////////////////////////////////////////////////
/// \brief Write transport output file headers in NetCDF format
/// \details Called prior to simulation (but after initialization) from CModel::Initialize()
/// \param &Options [in] Global model options information
//
void CConstituentModel::WriteNetCDFOutputFileHeaders(const optStruct& Options)
{
#ifdef _RVNETCDF_
  string kg,mgL,kgd;
  string filename;

  int iCumPrecip=_pModel->GetStateVarIndex(ATMOS_PRECIP);

  //units names
  kg="[kg]"; kgd="[kg/d]"; mgL="[mg/l]";
  if(Options.write_constitmass) {
    kg="[mg/m2]"; kgd="[mg/m2/d]"; mgL="[mg/m2]";
  }
  if(_type==TRACER) {
    kg="[-]"; kgd="[-]"; mgL="[-]";
  }
  else if(_type==ENTHALPY) {
    kg="[MJ]"; kgd="[MJ/d]"; mgL="["+to_string(DEG_SYMBOL)+"C]";
    if(Options.write_constitmass) {
      kg="[MJ/m2]"; kgd="[MJ/m2/d]"; mgL="[MJ/m2]";
    }
  }
  else if(_type==ISOTOPE) {
    kg="[kg-L]"; kgd="[kg-L/d]"; mgL="[o/oo]";
    if(Options.write_constitmass) {
      kg="[1/m2]"; kgd="[1/m2/d]"; mgL="[1/m2]"; //very hard to interpret these units
    }
  }

  time_struct tt;                       // start time structure
  int         time_dimid,varid_time;    // dimension ID (holds number of time steps) and variable ID (holds time values) for time
  const int   ndims1 = 1;
  const int   ndims2 = 2;
  int         dimids1[ndims1];          // array which will contain all dimension ids for a variable
  int         retval,ncid;

  //converts start day into "hours since YYYY-MM-DD HH:MM:SS"  (model start time)
  char  starttime[200]; // start time string in format 'hours since YYY-MM-DD HH:MM:SS'
  JulianConvert(0.0,Options.julian_start_day,Options.julian_start_year,Options.calendar,tt);
  strcpy(starttime,"hours since ");
  strcat(starttime,tt.date_string.c_str());
  strcat(starttime," ");
  strcat(starttime,DecDaysToHours(tt.julian_day,true).c_str());
  if(Options.time_zone!=0) { strcat(starttime,TimeZoneToString(Options.time_zone).c_str()); }

  // initialize all potential file IDs with -9 == "not existing and hence not opened"
  _CONC_ncid    = -9;
  _POLLUT_ncid  = -9;

  //Concentrations or Temperatures file
  //--------------------------------------------------------------------
  if(Options.write_constitmass) {
    if(_type!=ENTHALPY) { filename=_name+"Mass.nc"; }
    else                { filename="Enthalpy.nc"; }
  }
  else {
    if(_type!=ENTHALPY) { filename=_name+"Concentrations.nc"; }
    else                { filename="Temperatures.nc"; }
  }
  filename=FilenamePrepare(filename,Options);
  retval = nc_create(filename.c_str(),NC_CLOBBER|NC_NETCDF4,&ncid);  HandleNetCDFErrors(retval);
  _CONC_ncid = ncid;

  // ----------------------------------------------------------
  // global attributes
  // ----------------------------------------------------------
  WriteNetCDFGlobalAttributes(_CONC_ncid,Options,"Concentration/Temperature Output");

  // ----------------------------------------------------------
  // time vector
  // ----------------------------------------------------------
  // Define the DIMENSIONS. NetCDF will hand back an ID
  retval = nc_def_dim(_CONC_ncid,"time",NC_UNLIMITED,&time_dimid);  HandleNetCDFErrors(retval);

  /// Define the time variable.
  dimids1[0] = time_dimid;
  retval = nc_def_var(_CONC_ncid,"time",NC_DOUBLE,ndims1,dimids1,&varid_time); HandleNetCDFErrors(retval);
  retval = nc_put_att_text(_CONC_ncid,varid_time,"units",strlen(starttime),starttime);   HandleNetCDFErrors(retval);
  retval = nc_put_att_text(_CONC_ncid,varid_time,"calendar",strlen("gregorian"),"gregorian"); HandleNetCDFErrors(retval);


  int varid;
  varid= NetCDFAddMetadata(_CONC_ncid,time_dimid,"influx","influx",kgd);
  varid= NetCDFAddMetadata(_CONC_ncid,time_dimid,"channel_storage","Channel Storage",kg);
  varid= NetCDFAddMetadata(_CONC_ncid,time_dimid,"rivulet_storage","Rivulet Storage",kg);
  if(_type==ENTHALPY) {
    varid= NetCDFAddMetadata(_CONC_ncid,time_dimid,"air_temp","Air Temp",mgL);
  }
  for(int ii=0;ii<_pTransModel->GetNumWaterCompartments();ii++)
  {
    if(_pTransModel->GetStorWaterIndex(ii)!=iCumPrecip)
    {
      string sname=CStateVariable::SVTypeToString    (_pModel->GetStateVarType(_pTransModel->GetStorWaterIndex(ii)),_pModel->GetStateVarLayer(_pTransModel->GetStorWaterIndex(ii)));
      string name=CStateVariable::GetStateVarLongName(_pModel->GetStateVarType(_pTransModel->GetStorWaterIndex(ii)),_pModel->GetStateVarLayer(_pTransModel->GetStorWaterIndex(ii)));
      varid= NetCDFAddMetadata(_CONC_ncid,time_dimid,sname,name,mgL);
    }
  }
  varid= NetCDFAddMetadata(_CONC_ncid,time_dimid,"total","total",kg);
  varid= NetCDFAddMetadata(_CONC_ncid,time_dimid,"cum_loading","cumulative loading",kg);
  varid= NetCDFAddMetadata(_CONC_ncid,time_dimid,"cum_loss","cumulative loss",kg);
  varid= NetCDFAddMetadata(_CONC_ncid,time_dimid,"MB_error","mass/energy balance error",kg);

  // End define mode. This tells netCDF we are done defining metadata.
  retval = nc_enddef(_CONC_ncid);  HandleNetCDFErrors(retval);

  //Pollutograph / stream temperatures file
  //--------------------------------------------------------------------
  if(_type!=ENTHALPY) { filename=_name+"Pollutographs.nc"; }
  else                { filename="StreamTemperatures.nc"; }

  filename=FilenamePrepare(filename,Options);
  retval = nc_create(filename.c_str(),NC_CLOBBER|NC_NETCDF4,&ncid);  HandleNetCDFErrors(retval);
  _POLLUT_ncid = ncid;

  // ----------------------------------------------------------
  // global attributes
  // ----------------------------------------------------------
  WriteNetCDFGlobalAttributes(_POLLUT_ncid,Options,"Pollutograph/StreamTemperature Output");

  // ----------------------------------------------------------
  // time vector
  // ----------------------------------------------------------
  // Define the DIMENSIONS. NetCDF will hand back an ID
  retval = nc_def_dim(_POLLUT_ncid,"time",NC_UNLIMITED,&time_dimid);  HandleNetCDFErrors(retval);

  /// Define the time variable.
  dimids1[0] = time_dimid;
  retval = nc_def_var     (_POLLUT_ncid,"time",NC_DOUBLE,ndims1,dimids1,&varid_time); HandleNetCDFErrors(retval);
  retval = nc_put_att_text(_POLLUT_ncid,varid_time,"units",strlen(starttime),starttime);   HandleNetCDFErrors(retval);
  retval = nc_put_att_text(_POLLUT_ncid,varid_time,"calendar",strlen("gregorian"),"gregorian"); HandleNetCDFErrors(retval);

  // (a) count number of simulated outflows "nSim"
  string      tmp,tmp2,tmp3;
  int nbasins_dimid;
  int varid_bsim;
  int varid_Csim,varid_Cobs,varid_pctfroz;
  int nSim = 0;
  for(int p=0;p<_pModel->GetNumSubBasins();p++) {
    if(_pModel->GetSubBasin(p)->IsGauged()  && (_pModel->GetSubBasin(p)->IsEnabled())) { nSim++; }
  }

  if(nSim > 0)
  {
    // (b) create dimension "nbasins"
    retval = nc_def_dim(_POLLUT_ncid,"nbasins",nSim,&nbasins_dimid);                             HandleNetCDFErrors(retval);

    // (c) create variable  and set attributes for"basin_name"
    dimids1[0] = nbasins_dimid;
    retval = nc_def_var(_POLLUT_ncid,"basin_name",NC_STRING,ndims1,dimids1,&varid_bsim);       HandleNetCDFErrors(retval);
    tmp ="Name/ID of sub-basins with simulated concentrations";
    tmp2="timeseries_id";
    tmp3="1";
    retval = nc_put_att_text(_POLLUT_ncid,varid_bsim,"long_name",tmp.length(),tmp.c_str());    HandleNetCDFErrors(retval);
    retval = nc_put_att_text(_POLLUT_ncid,varid_bsim,"cf_role",tmp2.length(),tmp2.c_str());    HandleNetCDFErrors(retval);
    retval = nc_put_att_text(_POLLUT_ncid,varid_bsim,"units",tmp3.length(),tmp3.c_str());    HandleNetCDFErrors(retval);

    // (d) create 2D pollutograph arrays [nbasins x ntime]
    if(_type!=ENTHALPY) {
      varid_Csim    = NetCDFAddMetadata2D(_POLLUT_ncid,time_dimid,nbasins_dimid,"C_sim","Simulated concentrations",mgL);
      varid_Cobs    = NetCDFAddMetadata2D(_POLLUT_ncid,time_dimid,nbasins_dimid,"C_obs","Observed concentrations",mgL);
    }
    else {
      varid_Csim    = NetCDFAddMetadata2D(_POLLUT_ncid,time_dimid,nbasins_dimid,"T_sim","Simulated temperatures",mgL);
      varid_Cobs    = NetCDFAddMetadata2D(_POLLUT_ncid,time_dimid,nbasins_dimid,"T_obs","Observed temperatures",mgL);
      varid_pctfroz = NetCDFAddMetadata2D(_POLLUT_ncid,time_dimid,nbasins_dimid,"pct_froz","Percent frozen","0..1");
    }
  }// end if nSim>0

  if(_type==ENTHALPY) {
    varid= NetCDFAddMetadata(_POLLUT_ncid,time_dimid,"air_temp","Air Temp",mgL);
  }

  // End define mode. This tells netCDF we are done defining metadata.
  retval = nc_enddef(_POLLUT_ncid);  HandleNetCDFErrors(retval);

  // write values to NetCDF
  if(nSim > 0)
  {
    WriteNetCDFBasinList(_POLLUT_ncid,varid_bsim,_pModel,Options);
  }
#endif
}


//////////////////////////////////////////////////////////////////
/// \brief Writes minor transport output to file at the end of each timestep (or multiple thereof)
/// \note only thing this modifies should be output streams; called from CModel::WriteMinorOutput()
/// \param &Options [in] Global model options information
/// \param &tt [in] Local (model) time at the end of the pertinent time step
//
void CConstituentModel::WriteMinorOutput(const optStruct &Options,const time_struct &tt) 
{
  double currentMass,initMass; //[kg] or [MJ]
  double CumInflux,CumOutflux;
  double M; //[mg/m2] or [MJ/m2]
  double V; //[mm] 
  double concentration; //[mg/L] or [C]
  int    iCumPrecip;
  double convert;
  CEnthalpyModel *pEnthalpyModel=NULL; //Necessary evil to reduce code duplication rather than having an entirely separate WriteMinorOutput for Enthalpy

  string thisdate=tt.date_string;
  string thishour=DecDaysToHours(tt.julian_day);

  double area=_pModel->GetWatershedArea(); //[km2]

  iCumPrecip=_pModel->GetStateVarIndex(ATMOS_PRECIP);

  if((Options.suppressICs) && (tt.model_time==0.0)) { return; }

  convert=1.0/MG_PER_KG; //[mg->kg]
  if(_type==ENTHALPY)               { convert=1.0;                   } //[MJ]->[MJ]
  if(Options.write_constitmass)     { convert=1.0/(area*M2_PER_KM2); } //[mg->mg/m2] [MJ->MJ/m2]

  // Concentrations.csv or Temperatures.csv
  //----------------------------------------------------------------
  double atmos_prec  =0;//[mg] or [MJ] energy/mass advected in with snow/rain (calculated below)
  double influx      =0;//GetAverageInflux(c)*(area*M2_PER_KM2)*Options.timestep;//[mg] or [MJ] 
  // \todo [funct]: create GetAverageInflux() routine -
  //for Enthalpy, influx includes radiative heating/cooling of surface, sensible exchange, NOT latent heat flux 
  // should also include external sources from streams / diversions / etc.
  //for contaminant mass, influx includes distributed sources of contaminants
  //double net_influx =//GetAverageInflux(c)*(area*M2_PER_KM2)*Options.timestep; // [MJ] 
 
  double latent_flux=0;
  double Q_sens(0.0),Q_lat(0.0),Q_GW(0.0),Q_rad(0.0),Q_fric(0.0);
  double Qs,Ql,Qg,Qr,Qf;
  if(_type==ENTHALPY) 
  {
    pEnthalpyModel=(CEnthalpyModel*)(this);
    latent_flux =pEnthalpyModel->GetAvgLatentHeatFlux()*(area*M2_PER_KM2)*Options.timestep; // [MJ] (loss term)t)
    
    for(int p=0;p<_pModel->GetNumSubBasins();p++) {
      pEnthalpyModel->GetEnergyLossesFromReach(p,Qs,Ql,Qg,Qr,Qf);
      Q_sens+=Qs;Q_lat+=Ql;Q_GW+=Qg;Q_rad+=Qr; Q_fric+=Qf;
    }
  }
  double channel_stor=GetTotalChannelConstituentStorage();//[mg] or [MJ]
  double rivulet_stor=GetTotalRivuletConstituentStorage();//[mg] or [MJ]
  double sink        = _pModel->GetAvgStateVar(_pModel->GetStateVarIndex(CONSTITUENT_SINK,_constit_index))*(area*M2_PER_KM2);//[mg]  or [MJ] 
  double source      =-_pModel->GetAvgStateVar(_pModel->GetStateVarIndex(CONSTITUENT_SRC, _constit_index))*(area*M2_PER_KM2);//[mg]  or [MJ] 

  _OUTPUT<<tt.model_time <<","<<thisdate<<","<<thishour;
  if(_type==ENTHALPY) {
    _OUTPUT<<","<<_pModel->GetAvgForcing("TEMP_AVE");
  }

  if(tt.model_time!=0.0) { _OUTPUT<<","<<influx*convert;}
  else                   { _OUTPUT<<",---";             }
  _OUTPUT<<","<<channel_stor*convert;
  _OUTPUT<<","<<rivulet_stor*convert;

  currentMass=0.0;

  for(int ii=0;ii<_pTransModel->GetNumWaterCompartments();ii++)
  {
    //Get constituent concentration
    int m=(_constit_index*_pTransModel->GetNumWaterCompartments())+ii;

    M=_pModel->GetAvgStateVar(_pModel->GetStateVarIndex(CONSTITUENT,m)); //mass- mg/m2 or enthalpy - MJ/m2
    V=_pModel->GetAvgStateVar(_pTransModel->GetStorWaterIndex(ii)); //mm

    concentration=CalculateConcentration(M,V); // [mg/L] or [C] or [o/oo]
    if(Options.write_constitmass) { concentration=M; }//[mg/m2] or [MJ/m2] or [1/m2]

    if(_pTransModel->GetStorWaterIndex(ii)!=iCumPrecip)
    {
      _OUTPUT<<","<<concentration;     //print column entry 
      currentMass+=M*(area*M2_PER_KM2); //[mg]  or [MJ]  //increment total mass in system
    }
    else {
      atmos_prec-=M*(area*M2_PER_KM2);  //[mg]  or [MJ] //this M is always negative (loss from atmos.)
    }
  }
  currentMass+=channel_stor+rivulet_stor;  //[mg]  or [MJ] 

  initMass  =_initial_mass;

  CumInflux =source+atmos_prec+_cumul_input;        //[mg] or [MJ] 
  CumOutflux=sink  +_cumul_output;//outflow from system [mg] or [MJ] 

  //_OUTPUT<<setprecision(12);
  _OUTPUT<<","<<currentMass*convert;
  _OUTPUT<<","<<CumInflux  *convert;
  _OUTPUT<<","<<CumOutflux *convert;
  _OUTPUT<<","<<((currentMass-initMass)+(CumOutflux-CumInflux))*convert;
  if(_type==ENTHALPY)
  {
    _OUTPUT<<","<<sink       *convert;  // sink (incudes latent heat, dirichlet, latent heat)
    _OUTPUT<<","<<source     *convert;  // source (includes net surface flux,dirichlet) 
    _OUTPUT<<","<<atmos_prec *convert;  // atmospheric inputs 
    _OUTPUT<<","<<latent_flux*convert;  // latent heat
    _OUTPUT<<","<<Q_sens     *convert;  // reach sensible
    _OUTPUT<<","<<Q_lat      *convert;  // reach latent
    _OUTPUT<<","<<Q_GW       *convert;  // reach groundwater
    _OUTPUT<<","<<Q_rad      *convert;  // reach radiation
    _OUTPUT<<","<<Q_fric     *convert;  // reach friction
  }
  _OUTPUT<<endl;

  // Pollutographs.csv or StreamTemperatures.csv
  //----------------------------------------------------------------
  _POLLUT<<tt.model_time<<","<<thisdate<<","<<thishour;
  if(_type==ENTHALPY) {
    _POLLUT<<","<<_pModel->GetAvgForcing("TEMP_AVE"); //Air temperature
  }
  for(int p=0;p<_pModel->GetNumSubBasins();p++) {
    CSubBasin *pBasin=_pModel->GetSubBasin(p);
    if(pBasin->IsGauged() && (pBasin->IsEnabled()))
    {
      _POLLUT<<","<<GetOutflowConcentration(p);
      if(_type==ENTHALPY) {
        _POLLUT<<","<<pEnthalpyModel->GetOutflowIceFraction(p); 
      }
      for(int i = 0; i < _pModel->GetNumObservedTS(); i++) {
        if(IsContinuousConcObs(_pModel->GetObservedTS(i),pBasin->GetID(),_constit_index))
        {
           double val = _pModel->GetObservedTS(i)->GetAvgValue(tt.model_time,Options.timestep); 
           if((val != RAV_BLANK_DATA) && (tt.model_time>0)) { _POLLUT << "," << val; }
           else                                             { _POLLUT << ",";        }
        }
      }
    }
  }
  _POLLUT<<endl;
  
}

//////////////////////////////////////////////////////////////////
/// \brief Writes minor transport output to ensim formatted file at the end of each timestep (or multiple thereof)
/// \note only thing this modifies should be output streams; called from CModel::WriteMinorOutput()
/// \param &Options [in] Global model options information
/// \param &tt [in] Local (model) time at the end of the pertinent time step
/// \todo [reorg] merge with WriteMinorOutput - too much complex code repetition here when only difference is (1) delimeter and (2) timestep info included in the .csv file
//
void CConstituentModel::WriteEnsimMinorOutput(const optStruct &Options,const time_struct &tt) 
{
  double currentMass,CumInflux,CumOutflux,initMass; //[kg]
  double M; //[mg/m2]
  double V; //[mm]
  double concentration; //[mg/L]
  int    iCumPrecip;

  string thisdate=tt.date_string;
  string thishour=DecDaysToHours(tt.julian_day);

  double area=_pModel->GetWatershedArea(); //[km2]

  iCumPrecip=_pModel->GetStateVarIndex(ATMOS_PRECIP);

  double convert;//
  convert=1.0/MG_PER_KG; //[mg->kg]
  //convert=1.0/(area*M2_PER_KM2); //[mg->mg/m2]//TMP DEBUG OUTPUT OVERRIDE

  // Concentrations.tb0
  //----------------------------------------------------------------
  double influx      =0;//GetAverageInflux()*(area*M2_PER_KM2);//[mg/d]
  double channel_stor=GetTotalChannelConstituentStorage();//[mg]
  double rivulet_stor=GetTotalRivuletConstituentStorage();//[mg]

  if(tt.model_time!=0) { _OUTPUT<<" "<<influx*convert; }
  else { _OUTPUT<<" 0.0"; }
  _OUTPUT<<" "<<channel_stor*convert<<" "<<rivulet_stor*convert;

  currentMass=0.0;
  double atmos_prec=0;
  for(int ii=0;ii<_pTransModel->GetNumWaterCompartments();ii++)
  {
    //Get constituent concentration
    int m=ii+_constit_index*_pTransModel->GetNumWaterCompartments();
    M=_pModel->GetAvgStateVar(_pModel->GetStateVarIndex(CONSTITUENT,m)); //mg/m2
    V=_pModel->GetAvgStateVar(_pTransModel->GetStorWaterIndex(ii)); //mm
    concentration=CalculateConcentration(M,V);
    if(Options.write_constitmass) { concentration=M;}//[mg/m2]

    if(_pTransModel->GetStorWaterIndex(ii)!=iCumPrecip)
    {
      _OUTPUT<<" "<<concentration;   //print column entry
      currentMass+=M*(area*M2_PER_KM2); //mg          //increment total mass in system
    }
    else {
      atmos_prec+=M*(area*M2_PER_KM2); //mg
    }
  }
  currentMass+=channel_stor+rivulet_stor;


  CumInflux =-_pModel->GetAvgStateVar(_pModel->GetStateVarIndex(CONSTITUENT_SRC,_constit_index))*(area*M2_PER_KM2);//mg
  CumInflux +=-atmos_prec+_cumul_input;      //mg
  CumOutflux=_cumul_output;     //mg
  CumOutflux+=_pModel->GetAvgStateVar(_pModel->GetStateVarIndex(CONSTITUENT_SINK,_constit_index))*(area*M2_PER_KM2);//mg
  initMass  =_initial_mass;     //mg

  _OUTPUT<<" "<<currentMass*convert; //kg
  _OUTPUT<<" "<<CumInflux  *convert; //kg
  _OUTPUT<<" "<<CumOutflux *convert; //kg
  _OUTPUT<<" "<<((currentMass-initMass)+(CumOutflux-CumInflux))*convert; //kg
  _OUTPUT<<endl;

  // Pollutographs.tb0
  //----------------------------------------------------------------
  for(int p=0;p<_pModel->GetNumSubBasins();p++) {
    CSubBasin *pBasin=_pModel->GetSubBasin(p);
    if(pBasin->IsGauged() && pBasin->IsEnabled())
    {
      _POLLUT<<" "<<GetOutflowConcentration(p);
    }
  }
  _POLLUT<<endl;
}
//////////////////////////////////////////////////////////////////
/// \brief Writes minor transport output to NetCDF file at the end of each timestep (or multiple thereof)
/// \param &Options [in] Global model options information
/// \param &tt [in] Local (model) time at the end of the pertinent time step
//
void CConstituentModel::WriteNetCDFMinorOutput(const optStruct& Options,const time_struct& tt)
{
#ifdef _RVNETCDF_
  int    retval;                   // error value for NetCDF routines
  int    time_id;                  // variable id in NetCDF for time
  size_t time_index[1],count1[1];  // determines where and how much will be written to NetCDF; 1D variable (pre, time)
  size_t start2[2],count2[2];  // determines where and how much will be written to NetCDF; 2D variable (qsim, qobs, qin)
  double current_time[1];          // current time in hours since start time
  count1[0] = 1;                    // writes exactly one time step

  size_t time_ind2;
  current_time[0] = tt.model_time*HR_PER_DAY;
  current_time[0]=RoundToNearestMinute(current_time[0]);
  time_index[0] = int(rvn_round(tt.model_time/Options.timestep));   // element of NetCDF array that will be written
  time_ind2      =int(rvn_round(tt.model_time/Options.timestep));

  //====================================================================
  //  Concentrations.nc / Temperatures.nc
  //====================================================================
  double currentMass,initMass; //[kg] or [MJ]
  double CumInflux,CumOutflux;
  double M; //[mg/m2] or [MJ/m2]
  double V; //[mm] 
  double concentration; //[mg/L] or [C]
  double area      =_pModel->GetWatershedArea(); //[km2]
  int    iCumPrecip=_pModel->GetStateVarIndex(ATMOS_PRECIP);

  CEnthalpyModel *pEnthalpyModel=NULL; //Necessary evil to reduce code duplication rather than having an entirely separate WriteMinorOutput for Enthalpy
  if (_type==ENTHALPY){
    CEnthalpyModel *pEnthalpyModel=(CEnthalpyModel*)(this);
  }

  if((Options.suppressICs) && (tt.model_time==0.0)) { return; }
  
  double convert;
  convert=1.0/MG_PER_KG; //[mg->kg]
  if(_type==ENTHALPY)               { convert=1.0;                   } //[MJ]->[MJ]
  if(Options.write_constitmass)     { convert=1.0/(area*M2_PER_KM2); } //[mg->mg/m2] [MJ->MJ/m2]

  double atmos_prec  =0;//[mg] or [MJ] energy/mass advected in with snow/rain (calculated below)
  double influx      =0;//GetAverageInflux(c)*(area*M2_PER_KM2)*Options.timestep;//[mg] or [MJ] 
  double channel_stor=GetTotalChannelConstituentStorage();//[mg] or [MJ]
  double rivulet_stor=GetTotalRivuletConstituentStorage();//[mg] or [MJ]
  double sink        = _pModel->GetAvgStateVar(_pModel->GetStateVarIndex(CONSTITUENT_SINK,_constit_index))*(area*M2_PER_KM2);//[mg]  or [MJ] 
  double source      =-_pModel->GetAvgStateVar(_pModel->GetStateVarIndex(CONSTITUENT_SRC, _constit_index))*(area*M2_PER_KM2);//[mg]  or [MJ] 

  if(_type==ENTHALPY) {
    AddSingleValueToNetCDF(_CONC_ncid,"air_temp"  ,time_ind2,_pModel->GetAvgForcing("TEMP_AVE"));
  }

  double inf=influx*convert;
  if (tt.model_time==0.0){inf=NETCDF_BLANK_VALUE;}
  AddSingleValueToNetCDF(_CONC_ncid,"influx"           ,time_ind2,inf);
  AddSingleValueToNetCDF(_CONC_ncid,"channel_storage"  ,time_ind2,channel_stor*convert);
  AddSingleValueToNetCDF(_CONC_ncid,"rivulet_storage"  ,time_ind2,rivulet_stor*convert);

  currentMass=0.0;
  for(int ii=0;ii<_pTransModel->GetNumWaterCompartments();ii++)
  {
    //Get constituent concentration
    int m=(_constit_index*_pTransModel->GetNumWaterCompartments())+ii;
    M=_pModel->GetAvgStateVar(_pModel->GetStateVarIndex(CONSTITUENT,m)); //mass- mg/m2 or enthalpy - MJ/m2
    V=_pModel->GetAvgStateVar(_pTransModel->GetStorWaterIndex(ii)); //mm

    concentration=CalculateConcentration(M,V); // [mg/L] or [C] or [o/oo]
    if(Options.write_constitmass) { concentration=M; }//[mg/m2] or [MJ/m2] or [1/m2]

    if(_pTransModel->GetStorWaterIndex(ii)!=iCumPrecip)
    {
      string name=CStateVariable::SVTypeToString(_pModel->GetStateVarType(_pTransModel->GetStorWaterIndex(ii)),_pModel->GetStateVarLayer(_pTransModel->GetStorWaterIndex(ii)));
      AddSingleValueToNetCDF(_CONC_ncid, name.c_str(),time_ind2,concentration);
      currentMass+=M*(area*M2_PER_KM2); //[mg]  or [MJ]  //increment total mass in system
    }
    else {
      atmos_prec-=M*(area*M2_PER_KM2);  //[mg]  or [MJ] //this M is always negative (loss from atmos.)
    }
  }
  currentMass+=channel_stor+rivulet_stor;  //[mg]  or [MJ] 
  initMass  =_initial_mass;
  CumInflux =source+atmos_prec+_cumul_input;        //[mg] or [MJ] 
  CumOutflux=sink  +_cumul_output;//outflow from system [mg] or [MJ] 

  AddSingleValueToNetCDF(_CONC_ncid,"total"      ,time_ind2,currentMass*convert);
  AddSingleValueToNetCDF(_CONC_ncid,"cum_loading",time_ind2,CumInflux  *convert);
  AddSingleValueToNetCDF(_CONC_ncid,"cum_loss"   ,time_ind2,CumOutflux *convert);
  AddSingleValueToNetCDF(_CONC_ncid,"MB_error"   ,time_ind2,((currentMass-initMass)+(CumOutflux-CumInflux))*convert);

  //====================================================================
  //  Pollutographs.nc / StreamTemperatures.nc
  //====================================================================
  // (a) count how many values need to be written for q_obs, q_sim, q_in
  int nSim=0; //total # of sub-basins with simulated outflows
  for(int p=0;p<_pModel->GetNumSubBasins();p++) {
    if(_pModel->GetSubBasin(p)->IsGauged()  && (_pModel->GetSubBasin(p)->IsEnabled())) { nSim++; }
  }

  // (b) allocate memory if necessary
  double* C_obs  =NULL;   
  double* C_sim  =NULL;   
  double* pctfroz=NULL;    
  if(nSim>0) {
    C_obs  =new double[nSim];
    C_sim  =new double[nSim];
    pctfroz=new double[nSim];
  }

  // write new time step
  retval = nc_inq_varid(_POLLUT_ncid,"time",&time_id);                            HandleNetCDFErrors(retval);
  retval = nc_put_vara_double(_POLLUT_ncid,time_id,time_index,count1,&current_time[0]); HandleNetCDFErrors(retval);

  // write air temperature values
  if(_type==ENTHALPY) {
    int    temp_id;               // variable id in NetCDF for air temp
    double current_temp[1];
    current_temp[0]=_pModel->GetAvgForcing("TEMP_AVE");

    retval = nc_inq_varid(_POLLUT_ncid,"air_temp",&temp_id);                                HandleNetCDFErrors(retval);
    retval = nc_put_vara_double(_POLLUT_ncid,temp_id,time_index,count1,&current_temp[0]);   HandleNetCDFErrors(retval);
  }

  // get simulated conc/obs conc/pct froz values
  int iSim=0;
  for(int p=0;p<_pModel->GetNumSubBasins();p++) {
    CSubBasin* pBasin=_pModel->GetSubBasin(p);
    if(pBasin->IsGauged() && (pBasin->IsEnabled()))
    {
      C_sim[iSim]=GetOutflowConcentration(p);
      if(_type==ENTHALPY) {
        pctfroz[iSim]=pEnthalpyModel->GetOutflowIceFraction(p);
      }
      C_obs[iSim]=NETCDF_BLANK_VALUE;
      for(int i = 0; i < _pModel->GetNumObservedTS(); i++) {
        if(IsContinuousConcObs(_pModel->GetObservedTS(i),pBasin->GetID(),_constit_index))
        {
          double val = _pModel->GetObservedTS(i)->GetAvgValue(tt.model_time,Options.timestep);
          if((val != RAV_BLANK_DATA) && (tt.model_time>0)) { C_obs[iSim]=val; }
        }
      }
      iSim++;
    }
  }
  // write simulated conc/obs conc/pct froz values to file
  int    Csim_id;               // variable id in NetCDF for simulated concentration
  int    Cobs_id;               // variable id in NetCDF for observed concentration
  if(nSim > 0) {
    start2[0] = int(rvn_round(tt.model_time/Options.timestep)); // element of NetCDF array that will be written
    start2[1] = 0;                                              // element of NetCDF array that will be written
    count2[0] = 1;      // writes exactly one time step
    count2[1] = nSim;   // writes exactly nSim elements

    if(_type!=ENTHALPY) {
      retval = nc_inq_varid(_POLLUT_ncid,"C_sim",&Csim_id);                              HandleNetCDFErrors(retval);
      retval = nc_inq_varid(_POLLUT_ncid,"C_obs",&Cobs_id);                              HandleNetCDFErrors(retval);
      retval = nc_put_vara_double(_POLLUT_ncid,Csim_id,start2,count2,&C_sim[0]);         HandleNetCDFErrors(retval);
      retval = nc_put_vara_double(_POLLUT_ncid,Cobs_id,start2,count2,&C_obs[0]);         HandleNetCDFErrors(retval);
    }
    else if(_type==ENTHALPY) 
    {
      retval = nc_inq_varid(_POLLUT_ncid,"T_sim",&Csim_id);                              HandleNetCDFErrors(retval);
      retval = nc_inq_varid(_POLLUT_ncid,"T_obs",&Cobs_id);                              HandleNetCDFErrors(retval);
      retval = nc_put_vara_double(_POLLUT_ncid,Csim_id,start2,count2,&C_sim[0]);         HandleNetCDFErrors(retval);
      retval = nc_put_vara_double(_POLLUT_ncid,Cobs_id,start2,count2,&C_obs[0]);         HandleNetCDFErrors(retval);
      int pctfroz_id;
      retval = nc_inq_varid(_POLLUT_ncid,"pct_froz",&pctfroz_id);                        HandleNetCDFErrors(retval);
      retval = nc_put_vara_double(_POLLUT_ncid,pctfroz_id,start2,count2,&pctfroz[0]);    HandleNetCDFErrors(retval);
    }
  }

  delete[] C_obs;
  delete[] C_sim;
  delete[] pctfroz;

#endif
}
//////////////////////////////////////////////////////////////////
/// \brief Close transport output files
//
void CConstituentModel::CloseOutputFiles() 
{
  _OUTPUT.close();
  _POLLUT.close();

  #ifdef _RVNETCDF_
  int    retval;      // error value for NetCDF routines
  if (_CONC_ncid != -9)    {retval = nc_close(_CONC_ncid);    HandleNetCDFErrors(retval); }
  _CONC_ncid    = -9;
  if (_POLLUT_ncid != -9)  {retval = nc_close(_POLLUT_ncid);  HandleNetCDFErrors(retval); }
  _POLLUT_ncid  = -9;
  #endif
}

//////////////////////////////////////////////////////////////////
/// \brief Writes state variables to RVC file
//
void CConstituentModel::WriteMajorOutput(ofstream &RVC) const 
{
  RVC<<":BasinTransportVariables "<<_name<<endl;
  for(int p=0;p<_pModel->GetNumSubBasins();p++) 
  {
    RVC<<"  :BasinIndex "<<_pModel->GetSubBasin(p)->GetID()<<endl;
    int nSegs    =_pModel->GetSubBasin(p)->GetNumSegments();
    int nMlatHist=_pModel->GetSubBasin(p)->GetLatHistorySize();
    int nMinHist =_pModel->GetSubBasin(p)->GetInflowHistorySize();
    RVC<<"    :ChannelMass, "<<_channel_storage[p]<<endl;
    RVC<<"    :RivuletMass, "<<_rivulet_storage[p]<<endl;
    RVC<<"    :Mout,"<<nSegs;    for(int i=0;i<nSegs;    i++) { RVC<<","<<_aMout    [p][i]; }RVC<<","<<_aMout_last[p]<<endl;
    RVC<<"    :Mlat,"<<nMlatHist;for(int i=0;i<nMlatHist;i++) { RVC<<","<<_aMlatHist[p][i]; }RVC<<","<<_aMlat_last[p]<<endl;
    RVC<<"    :Min ,"<<nMinHist; for(int i=0;i<nMinHist; i++) { RVC<<","<<_aMinHist [p][i]; }RVC<<endl;
    if(_pModel->GetSubBasin(p)->GetReservoir()!=NULL) {
      RVC<<"    :ResMassOut, "<<_aMout_res[p]<<","<<_aMout_res_last[p]<<endl;
      RVC<<"    :ResMass, "   <<_aMres    [p]<<","<<_aMres_last    [p]<<endl;
    }
  }
  RVC<<":EndBasinTransportVariables"<<endl;
}
