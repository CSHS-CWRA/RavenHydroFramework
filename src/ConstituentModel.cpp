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

//////////////////////////////////////////////////////////////////
/// \brief Implentation of the Transport constructor
/// \param pModel [in] Model object
//
CConstituentModel::CConstituentModel(CModel *pMod,CTransportModel *pTMod, string name,constit_type typ,bool is_passive,const int c)
{
  _constit_index=c;

  _pModel=pMod; 
  _pTransModel=pTMod;

  _pConstituent=NULL;
  _pConstitParams=NULL;

  _pSources=NULL;
  _nSources=0;

  _nSpecFlowConcs=0;
  _pSpecFlowConcs=NULL;

  //_aIndexMapping=NULL; _nIndexMapping=0;
  _aSourceIndices=NULL;

  _aMinHist =NULL;
  _aMlatHist=NULL;
  _aMout    =NULL;

  _pConstituent=new constituent;

  //Initialize constituent members
  _pConstituent->name=name;
  _pConstituent->type=typ;
  _pConstituent->can_evaporate=false; //default behaviour - evaporation impossible for most contaminants
  if(_pConstituent->type==TRACER) {
    _pConstituent->can_evaporate=true;
  }
  if(_pConstituent->type==ENTHALPY) {
    _pConstituent->can_evaporate=true; //TMP DEBUG -to prevent accumulation of energy in soil/canopy stores during evaporation.
  }
  _pConstituent->is_passive=false;
  _pConstituent->initial_mass=0;
  _pConstituent->cumul_input =0;
  _pConstituent->cumul_output=0;

  //Parameter initialization
  _pConstitParams= new transport_params;
  InitializeConstitParams(_pConstitParams);
}
//////////////////////////////////////////////////////////////////
/// \brief Implentation of the Transport destructor
//
CConstituentModel::~CConstituentModel()
{
  delete[] _pConstituent;
  delete[] _pConstitParams; 
  delete[] _pSpecFlowConcs; _pSpecFlowConcs=NULL;
  for(int i=0;i<_nSources;i++) { delete _pSources[i]; } delete[] _pSources;
  delete[] _aSourceIndices;
  DeleteRoutingVars();
}

const constituent *CConstituentModel::GetConstituent() 
{
  return _pConstituent;
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
  return _pConstituent->is_passive;
}
//////////////////////////////////////////////////////////////////
/// \brief returns effective retardation factor for constituent c
/// being transported from storage compartment _iFromWater to storage compartment _iToWater
/// \param c [in] constituent index
/// \param _iFromWater [in] index of "from" water storage state variable
/// \param _iToWater [in] index of "to" water storage state variable
//
double CConstituentModel::GetRetardationFactor(const CHydroUnit *pHRU,const  int _iFromWater,const int _iToWater) const
{
  sv_type fromType,toType;
  fromType=_pModel->GetStateVarType(_iFromWater);
  toType  =_pModel->GetStateVarType(_iToWater);

  if(IsPassive()) { return ALMOST_INF; }
  if(fromType==SOIL)
  {
    int    m=_pModel->GetStateVarLayer(_iFromWater);
#ifdef _STRICTCHECK_
    ExitGracefullyIf(m==DOESNT_EXIST,"GetRetardationFactor:invalid _iFromWater",RUNTIME_ERR);
#endif

    if(toType!=ATMOSPHERE) { return pHRU->GetSoilProps(m)->retardation[_constit_index]; }
    else                   { return pHRU->GetVegetationProps()->uptake_moderator[_constit_index]; }
  }
  return 1.0;
}
//////////////////////////////////////////////////////////////////
/// \brief returns decay coefficient for constituent c
/// \param c [in] constituent index
/// \param iStorWater [in] index of water storage state variable
//
double CConstituentModel::GetDecayCoefficient(const CHydroUnit *pHRU,const int iStorWater) const
{
  sv_type storType;
  double decay_coeff = _pConstitParams->decay_coeff;

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

}//////////////////////////////////////////////////////////////////
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
void CConstituentModel::Initialize()
{
  
  //double area=pModel->GetWatershedArea();

  // Calculate initial mass and zero out zero volume concentrations
  //--------------------------------------------------------------------
  _pConstituent->cumul_input=0;
  _pConstituent->cumul_output=0;

  /// \todo [funct]: calculate initial mass from Dirichlet cells

  // populate array of source indices
  //--------------------------------------------------------------------
  // \todo [funct] will have to revise to support different sources in different HRUs (e.g., _aSourceIndices[c][i_stor][k])
  _aSourceIndices=NULL;
  _aSourceIndices=new int[_pModel->GetNumStateVars()];
  ExitGracefullyIf(_aSourceIndices==NULL,"CTransport::Initialize",OUT_OF_MEMORY);
  for(int i_stor=0;i_stor<_pModel->GetNumStateVars();i_stor++) {
    _aSourceIndices[i_stor]=DOESNT_EXIST;
    for(int i=0;i<_nSources;i++) {
      if(_pSources[i]->i_stor==i_stor)
      {
        if(_aSourceIndices[i_stor] != DOESNT_EXIST) {
          WriteWarning("CConstituentModel::Intiialize: cannot currently have more than one constitutent source per constituent/storage combination",false);
        }
        _aSourceIndices[i_stor]=i; //each
      }
    }
  }

  // Initialize routing variables
  //--------------------------------------------------------------------
  InitializeRoutingVars();

}
//////////////////////////////////////////////////////////////////
/// \brief Increment cumulative input of mass to watershed.
/// \note called within solver to track mass balance
/// \param Options [in] Global model options structure
/// \param tt [in] current model time structure
//
void CConstituentModel::IncrementCumulInput(const optStruct &Options,const time_struct &tt)
{
  _pConstituent->cumul_input+=0;// \todo [funct]: increment cumulative input [mg]
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
      _pConstituent->cumul_output+=GetIntegratedMassOutflow(p,Options.timestep);
    }
  }
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

  int i_source=_aSourceIndices[i_stor];
  if(i_source==DOESNT_EXIST) { return false; }
  if(!_pSources[i_source]->dirichlet) { return false; }
  Cs = _pSources[i_source]->concentration;

  if(_pSources[i_source]->kk==DOESNT_EXIST)
  { //Not tied to HRU Group
    if(Cs != DOESNT_EXIST) { return true; }
    else {//time series
      Cs = _pSources[i_source]->pTS->GetValue(tt.model_time);
      return true;
    }
  }
  else
  { //Tied to HRU Group - Check if we are in HRU Group
    if(_pModel->GetHRUGroup(_pSources[i_source]->kk)->IsInGroup(k))
    {
      if(Cs != DOESNT_EXIST) { return true; }
      else {//time series
        Cs = _pSources[i_source]->pTS->GetValue(tt.model_time);
        return true;
      }
    }
  }
  return false;
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
  double flux;
  bool retrieve=false;
  int i_source=_aSourceIndices[i_stor];
  if(i_source == DOESNT_EXIST) { return 0.0; }
  if(_pSources[i_source]->dirichlet) { return 0.0; }

  if(_pSources[i_source]->kk==DOESNT_EXIST)//not tied to HRU group
  {
    retrieve=true;
  }
  else { //Check if we are in HRU Group
    retrieve=_pModel->GetHRUGroup(_pSources[i_source]->kk)->IsInGroup(k);
  }
  if(retrieve)
  {
    flux=_pSources[i_source]->concentration; //'concentration' stores mass flux [mg/m2/d] if this is a flux-source
    if(flux == DOESNT_EXIST) {//get from time series
      flux=_pSources[i_source]->pTS->GetValue(tt.model_time);
    }
    return flux;
  }
  else { return 0.0; }

}
void   CConstituentModel::SetInitialMass(const double &mass) {
  _pConstituent->initial_mass=mass;
}

void  CConstituentModel::ApplyConvolutionRouting(const int p,const double *aRouteHydro,const int nSegments,const int nMinHist,const double *aMinHist,const double &tstep,double *aMout_new) const
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
  ofstream OUT;

  int iCumPrecip=_pModel->GetStateVarIndex(ATMOS_PRECIP);

  string kg,mgL,kgd;

  //units names
  kg="[kg]"; kgd="[kg/d]"; mgL="[mg/l]";
  if(Options.write_constitmass) {
    kg="[mg/m2]"; kgd="[mg/m2/d]"; mgL="[mg/m2]";
  }
  if(_pConstituent->type==TRACER) {
    kg="[-]"; kgd="[-]"; mgL="[-]";
  }
  else if(_pConstituent->type==ENTHALPY) {
    kg="[MJ]"; kgd="[MJ/d]"; mgL="[C]";
    if(Options.write_constitmass) {
      kg="[MJ/m2]"; kgd="[MJ/m2/d]"; mgL="[MJ/m2]";
    }
  }

  //Concentrations file
  //--------------------------------------------------------------------
  if(_pConstituent->type!=ENTHALPY) {
    filename=_pConstituent->name+"Concentrations.csv";
    if(Options.write_constitmass) { filename=_pConstituent->name+"Mass.csv"; }
  }
  else {
    filename="Temperatures.csv";
    if(Options.write_constitmass) { filename="Enthalpy.csv"; }
  }
  filename=FilenamePrepare(filename,Options);

  _OUTPUT.open(filename.c_str());
  if(_OUTPUT.fail()) {
    ExitGracefully(("CConstituentModel::WriteOutputFileHeaders: Unable to open output file "+filename+" for writing.").c_str(),FILE_OPEN_ERR);
  }

  //    Header content ---------------------------
  _OUTPUT<<"time[d],date,hour,influx"<<kgd<<",Channel Storage"<<kg<<",Rivulet Storage"<<kg;

  for(int ii=0;ii<_pTransModel->GetNumWaterCompartments();ii++)
  {
    if(_pTransModel->GetStorWaterIndex(ii)!=iCumPrecip)
    {
      _OUTPUT<<","<<
        CStateVariable::GetStateVarLongName(_pModel->GetStateVarType(_pTransModel->GetStorWaterIndex(ii)),_pModel->GetStateVarLayer(_pTransModel->GetStorWaterIndex(ii)))<<" "<<mgL;
    }
  }

  if(_pConstituent->type!=ENTHALPY) {
    _OUTPUT<<", Total Mass "<<kg<<", Cum. Loading "<<kg<<", Cum. Mass Lost "<<kg<<", MB Error "<<kg;
    _OUTPUT<<", atmos "<<kg<<", sink"<<kg<<endl;//TMP DEBUG
  }
  else {
    _OUTPUT<<", Total Energy "<<kg<<", Cum. Loading "<<kg<<", Cum. Energy Lost "<<kg<<", EB Error "<<kg;
    _OUTPUT<<", sink "<<kg<<", source "<<kg<<", atmos "<<", latent heat "<<endl;//TMP DEBUG
  }

  //Pollutograph / stream temperatures file
  //--------------------------------------------------------------------
  if(_pConstituent->type!=ENTHALPY) {
    filename=_pConstituent->name+"Pollutographs.csv";
  }
  else {
    filename="StreamTemperatures.csv";
  }

  filename=FilenamePrepare(filename,Options);
  _POLLUT.open(filename.c_str());
  if(_POLLUT.fail()) {
    ExitGracefully(("CConstituentModel::WriteOutputFileHeaders: Unable to open output file "+filename+" for writing.").c_str(),FILE_OPEN_ERR);
  }

  _POLLUT<<"time[d],date,hour";
  const CSubBasin *pBasin;
  for(int p=0;p<_pModel->GetNumSubBasins();p++) {
    pBasin=_pModel->GetSubBasin(p);
    if(pBasin->IsGauged() && (pBasin->IsEnabled())) {
      string name;
      if(pBasin->GetName()=="")   { _POLLUT<<",ID="<<pBasin->GetID()  <<" "<<mgL; }
      else                        { _POLLUT<<","   <<pBasin->GetName()<<" "<<mgL; }
      if(_pConstituent->type==ENTHALPY) {
        if(pBasin->GetName()=="") { _POLLUT<<",ID="<<pBasin->GetID()  <<" pct froz."; }
        else                      { _POLLUT<<","   <<pBasin->GetName()<<" pct froz."; }
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
  // \todo[funct] - support enthalpy/temperature in ensim output 
  string filename;
  ofstream OUT;

  int iCumPrecip=_pModel->GetStateVarIndex(ATMOS_PRECIP);

  string kg,mgL,kgd;

  int i;
  time_struct tt,tt2;
  JulianConvert(0.0,Options.julian_start_day,Options.julian_start_year,Options.calendar,tt);
  JulianConvert(Options.timestep,Options.julian_start_day,Options.julian_start_year,Options.calendar,tt2);//end of the timestep


  //units names
  kg="kg"; kgd="kg/d"; mgL="mg/l";
  //kg="mg/m2"; kgd="mg/m2/d"; mgL="mg/m2";//TMP DEBUG OUTPUT OVERRIDE
  if(_pConstituent->type==TRACER) {
    kg="none"; kgd="none"; mgL="none";
  }

  //Concentrations file
  //--------------------------------------------------------------------
  filename=_pConstituent->name+"Concentrations.tb0";
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
    _OUTPUT<<":CreationDate  "<<GetCurrentTime()<<endl;
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
  filename=_pConstituent->name+"Pollutographs.tb0";
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
    _POLLUT<<":CreationDate  "<<GetCurrentTime()<<endl;
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
/// \brief Writes minor transport output to file at the end of each timestep (or multiple thereof)
/// \note only thing this modifies should be output streams; called from CModel::WriteMinorOutput()
/// \param &Options [in] Global model options information
/// \param &tt [in] Local (model) time at the end of the pertinent time step
//
void CConstituentModel::WriteMinorOutput(const optStruct &Options,const time_struct &tt) 
{
  double currentMass,CumInflux,CumOutflux,initMass; //[kg] or [MJ]
  double M; //[mg/m2] or [MJ/m2]
  double V; //[mm] 
  double concentration; //[mg/L] or [C]
  int    iCumPrecip;
  double convert;
  CEnthalpyModel *pEnthalpyModel=NULL;

  string thisdate=tt.date_string;
  string thishour=DecDaysToHours(tt.julian_day);

  double area=_pModel->GetWatershedArea(); //[km2]

  iCumPrecip=_pModel->GetStateVarIndex(ATMOS_PRECIP);

  if((Options.suppressICs) && (tt.model_time==0.0)) { return; }

  convert=1.0/MG_PER_KG; //[mg->kg]
  if(_pConstituent->type==ENTHALPY) { convert=1.0; } //[MJ]->[MJ]
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
  if(_pConstituent->type==ENTHALPY) {
    pEnthalpyModel=(CEnthalpyModel*)(this);
    latent_flux =pEnthalpyModel->GetAvgLatentHeatFlux()*(area*M2_PER_KM2)*Options.timestep; // [MJ] (loss term)t)
  }
  //this->GetEnergyLossesFromReach

  double channel_stor=GetTotalChannelConstituentStorage();//[mg] or [MJ]
  double rivulet_stor=GetTotalRivuletConstituentStorage();//[mg] or [MJ]
  double sink        = _pModel->GetAvgStateVar(_pModel->GetStateVarIndex(CONSTITUENT_SINK,_constit_index))*(area*M2_PER_KM2);//[mg]  or [MJ] 
  double source      =-_pModel->GetAvgStateVar(_pModel->GetStateVarIndex(CONSTITUENT_SRC,_constit_index))*(area*M2_PER_KM2);//[mg]  or [MJ] 

  source+=GetMassAddedFromInflowSources(tt.model_time,Options.timestep);

  _OUTPUT<<tt.model_time <<","<<thisdate<<","<<thishour;

  if(tt.model_time!=0.0) {
    _OUTPUT<<","<<influx*convert;
  }
  else {
    _OUTPUT<<",---";
  }
  _OUTPUT<<","<<channel_stor*convert;
  _OUTPUT<<","<<rivulet_stor*convert;

  currentMass=0.0;

  for(int ii=0;ii<_pTransModel->GetNumWaterCompartments();ii++)
  {
    //Get constituent concentration
    int m=(_constit_index*_pTransModel->GetNumWaterCompartments())+ii;

    //concentration=_pTransModel->GetConcentration(k,_pModel->GetStateVarIndex(CONSTITUENT,m));
    //if(Options.write_constitmass) { concentration=_pModel->GetAvgStateVar(_pModel->GetStateVarIndex(CONSTITUENT,m)); }

    M=_pModel->GetAvgStateVar(_pModel->GetStateVarIndex(CONSTITUENT,m)); //mass- mg/m2 or enthalpy - MJ/m2
    V=_pModel->GetAvgStateVar(_pTransModel->GetStorWaterIndex(ii)); //mm
    
    if(_pConstituent->type!=ENTHALPY)
    {
      if(fabs(V)<=1e-6) { concentration=0.0; }
      else              { concentration=(M/V)*(MM_PER_METER/LITER_PER_M3); }//[mg/mm/m2]->[mg/L]
    }
    else {
      if(fabs(V)<=1e-6) { concentration=0.0; } // JRC: should this default to zero? or NA?
      else              { concentration=ConvertVolumetricEnthalpyToTemperature(M/V*MM_PER_METER); } //[MJ/m3]->[C]
    }

    if(Options.write_constitmass) { concentration=M; }//[mg/m2] or [MJ/m2]

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

  initMass  =_pConstituent->initial_mass;

  CumInflux =source+atmos_prec;                     //[mg] or [MJ] 
  CumOutflux=sink  +_pConstituent->cumul_output;//outflow from system [mg] or [MJ] 

  _OUTPUT<<","<<currentMass*convert;
  _OUTPUT<<","<<CumInflux  *convert;
  _OUTPUT<<","<<CumOutflux *convert;
  _OUTPUT<<","<<((currentMass-initMass)+(CumOutflux-CumInflux))*convert;
  if(_pConstituent->type==ENTHALPY)
  {
    _OUTPUT<<","<<sink       *convert;  // sink (incudes latent heat, dirichlet, latent heat)
    _OUTPUT<<","<<source     *convert;  // source (includes net surface flux,dirichlet) 
    _OUTPUT<<","<<atmos_prec *convert;  // atmospheric inputs 
    _OUTPUT<<","<<latent_flux*convert;  // latent heat
  }
  _OUTPUT<<endl;

  // Pollutographs.csv or StreamTemperatures.csv
  //----------------------------------------------------------------
  _POLLUT<<tt.model_time<<","<<thisdate<<","<<thishour;
  for(int p=0;p<_pModel->GetNumSubBasins();p++) {
    CSubBasin *pBasin=_pModel->GetSubBasin(p);
    if(pBasin->IsGauged() && (pBasin->IsEnabled()))
    {
      _POLLUT<<","<<GetOutflowConcentration(p);
      if(_pConstituent->type==ENTHALPY) {
        _POLLUT<<","<<pEnthalpyModel->GetOutflowIceFraction(p); 
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
    if(fabs(V)<=1e-6)             { concentration=0.0; }
    else                          { concentration=(M/V)*(MM_PER_METER/LITER_PER_M3); }//[mg/mm/m2]->[mg/L]
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
  CumInflux +=-atmos_prec;                    //mg
  CumOutflux=_pConstituent->cumul_output;     //mg
  CumOutflux+=_pModel->GetAvgStateVar(_pModel->GetStateVarIndex(CONSTITUENT_SINK,_constit_index))*(area*M2_PER_KM2);//mg
  initMass  =_pConstituent->initial_mass;     //mg

  _OUTPUT<<" "<<currentMass*convert; //kg
  _OUTPUT<<" "<<CumInflux*convert;   //kg
  _OUTPUT<<" "<<CumOutflux*convert;  //kg
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
/// \brief Close transport output files
//
void CConstituentModel::CloseOutputFiles() 
{
  _OUTPUT.close();
  _POLLUT.close();
}
