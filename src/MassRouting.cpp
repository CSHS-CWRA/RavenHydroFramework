/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2022 the Raven Development Team
------------------------------------------------------------------
routing of mass/energy in catchment and channel
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "SubBasin.h"
#include "Model.h"
#include "Transport.h"
//////////////////////////////////////////////////////////////////
/// \brief Initialize mass routing variables
///
void CConstituentModel::InitializeRoutingVars()
{
  int nSB=_pModel->GetNumSubBasins();
  _nMinHist       =new int     [nSB];
  _nMlatHist      =new int     [nSB];
  _aMinHist       =new double *[nSB];
  _aMlatHist      =new double *[nSB];
  _aMout          =new double *[nSB];
  _aMout_last     =new double  [nSB];
  _aMres          =new double  [nSB];
  _aMres_last     =new double  [nSB];
  _aMsed          =new double  [nSB];
  _aMsed_last     =new double  [nSB];
  _aMlat_last     =new double  [nSB];
  _aMout_res      =new double  [nSB];
  _aMout_res_last =new double  [nSB];
  _aMresRain      =new double  [nSB];
  _channel_storage=new double  [nSB];
  _rivulet_storage=new double  [nSB];
  _aMlocal        =new double  [nSB];
  _aMlocLast      =new double  [nSB];

  for(int p=0;p<nSB;p++)
  {
    int nSegments=_pModel->GetSubBasin(p)->GetNumSegments();

    _aMout[p]=NULL;
    _nMlatHist[p]=_pModel->GetSubBasin(p)->GetLatHistorySize();
    _nMinHist [p]=_pModel->GetSubBasin(p)->GetInflowHistorySize(); //TMP DEBUG - to change
    _aMinHist [p]=new double[_nMinHist[p]];
    _aMlatHist[p]=new double[_nMlatHist[p]];
    _aMout    [p]=new double[nSegments];
    ExitGracefullyIf(_aMout[p]==NULL,"CConstituentModel::InitializeRoutingVars(2)",OUT_OF_MEMORY);
    for(int i=0; i<_nMinHist[p]; i++) { _aMinHist [p][i]=0.0; }
    for(int i=0; i<_nMlatHist[p];i++) { _aMlatHist[p][i]=0.0; }
    for(int i=0; i<nSegments;    i++) { _aMout    [p][i]=0.0; }
    _aMout_last     [p]=0.0;
    _aMres          [p]=0.0;
    _aMres_last     [p]=0.0;
    _aMsed          [p]=0.0;
    _aMsed_last     [p]=0.0;
    _aMlat_last     [p]=0.0;
    _aMout_res      [p]=0.0;
    _aMout_res_last [p]=0.0;
    _aMresRain      [p]=0.0;
    _channel_storage[p]=0.0;
    _rivulet_storage[p]=0.0;
    _aMlocal        [p]=0.0;
    _aMlocLast      [p]=0.0;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Delete mass routing variables
///
void CConstituentModel::DeleteRoutingVars()
{
  int nSB=_pModel->GetNumSubBasins();
  if(_aMinHist!=NULL) {
    for(int p=0;p<nSB;p++)
    {
      delete[] _aMinHist[p];
      delete[] _aMlatHist[p];
      delete[] _aMout[p];
    }
    delete[] _aMinHist;        _aMinHist  =NULL;
    delete[] _aMlatHist;       _aMlatHist =NULL;
    delete[] _aMout;           _aMout     =NULL;
    delete[] _aMres;           _aMres     =NULL;
    delete[] _aMres_last;      _aMres_last=NULL;
    delete[] _aMsed;           _aMsed = NULL;
    delete[] _aMsed_last;      _aMsed_last = NULL;
    delete[] _aMlat_last;      _aMlat_last=NULL;
    delete[] _aMout_res;       _aMout_res =NULL;
    delete[] _aMout_res_last;  _aMout_res_last =NULL;
    delete[] _aMresRain;       _aMresRain =NULL;
    delete[] _channel_storage; _channel_storage=NULL;
    delete[] _rivulet_storage; _rivulet_storage=NULL;
    delete[] _aMout_last;      _aMout_last     =NULL;
    delete[] _nMlatHist;       _nMlatHist=NULL;
    delete[] _nMinHist;        _nMinHist=NULL;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Set upstream mass loading to primary channel and updates mass loading history
///
/// \param p    subbasin index
/// \param Minnew rates of upstream mass/energy loading for the current timestep [mg/d]/[MJ/d]
//
void   CConstituentModel::SetMassInflows(const int p,const double Minnew)
{
  for(int n=_nMinHist[p]-1;n>0;n--) {
    _aMinHist[p][n]=_aMinHist[p][n-1];
  }
  _aMinHist[p][0]=Minnew;
}
//////////////////////////////////////////////////////////////////
/// \brief Updates aMinnew, array of mass loadings, to handle fixed concentration/temperature or specified mass inflow conditions
///
/// \param p    subbasin index
/// \param t    end of timestep, in model time
/// \param Minnew rate of upstream mass/energy loading for the current timestep [mg/d]/[MJ/d]
//
void   CConstituentModel::ApplySpecifiedMassInflows(const int p,const double t,double &Minnew)
{
  double C;
  long long SBID=_pModel->GetSubBasin(p)->GetID();
  double       Q=_pModel->GetSubBasin(p)->GetOutflowRate()*SEC_PER_DAY; //[m3/d] Flow at end of time step

  //Handle additional mass/energy inflows
  for(int i=0; i<_nMassLoadingTS; i++) {
    if(_pMassLoadingTS[i]->GetLocID()==SBID) {
      Minnew+=_pMassLoadingTS[i]->GetValue(t)*MG_PER_KG; //[mg/d]
    }
  }

  //Handle specified concentrations/temperatures of streamflow
  for(int i=0; i<_nSpecFlowConcs; i++)
  {
    if(_pSpecFlowConcs[i]->GetLocID()==SBID) {
      C=_pSpecFlowConcs[i]->GetValue(t); //mg/L or C
      if(_type==ENTHALPY) { C=ConvertTemperatureToVolumetricEnthalpy(C,0.0); } //MJ/m3
      else                { C*=LITER_PER_M3; } //mg/m3
      ExitGracefullyIf(_type==ISOTOPE,"ApplySpecifiedMassInflows: cannot handle isotopes",BAD_DATA);

      Minnew=Q*C;  //[m3/d]*[mg/m3] or [m3/d]*[MJ/m3]
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief determines net mass/energy added in this timestep due to inflow source terms
///
/// \param c    constituent index
/// \param t    end of timestep, in model time
/// \param tstep [in]
/// \returns net mass/energy added in this timestep due to inflow source terms [mg] or [MJ]
//
double   CConstituentModel::GetMassAddedFromInflowSources(const double &t,const double &tstep) const
{
  double C,Q,Qold;
  long long SBID,SBID_down;
  double mass=0; //[mg] or [MJ]

  // remove mass which is overridden by specified mass inflow
  for(int p=0;p<_pModel->GetNumSubBasins();p++)
  {
    SBID_down=_pModel->GetSubBasin(p)->GetDownstreamID();
    for(int i=0; i<_nSpecFlowConcs; i++) {
      if(_pSpecFlowConcs[i]->GetLocID()==SBID_down) {
        mass-=GetIntegratedMassOutflow(p,tstep);
      }
    }
  }

  // Handle specified concentrations/temperatures of streamflow
  for(int i=0; i<_nSpecFlowConcs; i++)
  {
    SBID  =_pSpecFlowConcs[i]->GetLocID();

    Q=_pModel->GetSubBasinByID(SBID)->GetOutflowRate()*SEC_PER_DAY; //[m3/d] Flow at end of time step
    C=_pSpecFlowConcs[i]->GetValue(t+tstep);                             //mg/L or C
    if(_type==ENTHALPY) { C=ConvertTemperatureToVolumetricEnthalpy(C,0.0); } //MJ/m3
    else                { C*=LITER_PER_M3; } //mg/m3

    //TMP DEBUG - DOES NOT HANDLE ISOTOPES PROPERLY !!!
    ExitGracefullyIf(_type==ISOTOPE,"GetMassAddedFromInflowSources: cannot handle isotopes",BAD_DATA);

    mass+=0.5*Q*C*tstep;  //[mg] or [MJ]

    if(t>0) {
      Qold=_pModel->GetSubBasinByID(SBID)->GetLastOutflowRate()*SEC_PER_DAY;//[m3/d] Flow at start of time step
      C=_pSpecFlowConcs[i]->GetValue(t);                             //mg/L or C
      if(_type==ENTHALPY) { C=ConvertTemperatureToVolumetricEnthalpy(C,0.0); } //MJ/m3
      else                { C*=LITER_PER_M3; } //mg/m3
      mass+=0.5*Qold*C*tstep;
    }
  }
  // Handle external additions of mass
  for(int i=0; i<_nMassLoadingTS; i++)
  {
    mass+=0.5*_pMassLoadingTS[i]->GetValue(t+tstep)*MG_PER_KG*tstep;  //[mg]
    if(t>0) {
    mass+=0.5*_pMassLoadingTS[i]->GetValue(t      )*MG_PER_KG*tstep;  //[mg]
    }
  }

  return mass;
}
//////////////////////////////////////////////////////////////////
/// \brief Sets lateral mass/energy loading to primary channel and updates mass/energy loading history
/// \param p  [in]   subbasin index
/// \param Mlat [in] lateral mass loading rate for the current timestep [mg/d][MJ/d]
//
void   CConstituentModel::SetLateralInfluxes(const int p,const double Mlat)
{
  for(int n=_nMlatHist[p]-1;n>0;n--) {
    _aMlatHist[p][n]=_aMlatHist[p][n-1];
  }
  _aMlatHist[p][0]=Mlat;

}
//////////////////////////////////////////////////////////////////
/// \brief Gathers and prepares data prior to calling RouteMass
/// \param p [in]    subbasin index
//
void  CConstituentModel::PrepareForRouting(const int p) {
  //does nothing - interesting in child classes
}
void  CConstituentModel::PrepareForInCatchmentRouting(const int p) {
  //does nothing - interesting in child classes
}
void CConstituentModel::InCatchmentRoute(const int p, double &Mlat_new, const optStruct &Options)
{
  PrepareForInCatchmentRouting(p);

  const double * aUnitHydro =_pModel->GetSubBasin(p)->GetUnitHydrograph();
  const double * aQlatHist  =_pModel->GetSubBasin(p)->GetLatHistory();
  Mlat_new=ApplyInCatchmentRouting(p,aUnitHydro,aQlatHist,_aMlatHist[p],_nMlatHist[p], Options.timestep);

  _aMlocLast[p] = _aMlocal[p];
  _aMlocal  [p] = Mlat_new;
}
    //////////////////////////////////////////////////////////////////
/// \brief Creates aMoutnew [m^3/s], an array of point measurements for outflow at downstream end of each river segment
/// \details Array represents measurements at the end of the current timestep. if (catchment_routing==distributed),
/// lateral flow Qlat (e.g., runoff, interflow, or baseflow) is routed as part of channel routing routine otherwise,
/// lateral flow is all routed to most downstream outflow segment . Assumes uniform time step
///
/// \param p           [in]  subbasin index
/// \param aMout_new[] [out] Array of mass outflows at downstream end of each segment at end of current timestep [mg/d] [size: nsegs]
/// \param Mlat_new    [out] in-catchment routing component of aMout_new[nSegs-1] [mg/d]
/// \param Res_mass    [out] calculated reservoir mass at end of time step
/// \param ResSedMass  [out] calculated reservoir mass in sediment at end of time step
/// \param &Options    [in]  Global model options information
/// \param tt          [in]  Time structure
//
void   CConstituentModel::RouteMass(const int          p,          // SB index
                                    double            *aMout_new,  // [mg/d][size: nsegs ]
                                    double            &Mlat_new,   // [mg/d]
                                    double            &Res_mass,   // [mg]
                                    double            &ResSedMass, // [mg]
                                    const optStruct   &Options,
                                    const time_struct &tt) const
{
  //==============================================================
  // route from catchment
  //==============================================================
  //const double * aUnitHydro =_pModel->GetSubBasin(p)->GetUnitHydrograph();
  //const double * aQlatHist  =_pModel->GetSubBasin(p)->GetLatHistory();
  //Mlat_new=ApplyInCatchmentRouting(p,aUnitHydro,aQlatHist,_aMlatHist[p],_nMlatHist[p], Options.timestep);

  //==============================================================
  // route along channel
  //==============================================================
  const double * aRouteHydro=_pModel->GetSubBasin(p)->GetRoutingHydrograph();
  const double * aQinHist   =_pModel->GetSubBasin(p)->GetInflowHistory();

  int nSegments   =_pModel->GetSubBasin(p)->GetNumSegments();
  double seg_fraction=1.0/(double)(nSegments);

  if((Options.routing==ROUTE_PLUG_FLOW) || (Options.routing==ROUTE_DIFFUSIVE_WAVE))
  {
    // (sometimes fancy) convolution
    ApplyConvolutionRouting(p,aRouteHydro,aQinHist,_aMinHist[p],nSegments,_nMinHist[p],Options.timestep,aMout_new);
  }
  //==============================================================
  else if(Options.routing==ROUTE_NONE)
  {//In channel routing instantaneous
    aMout_new[nSegments-1]=_aMinHist[p][0];//spits out the mass that just entered
  }
  //==============================================================
  else {
    ExitGracefully("Unrecognized or unsupported constiuent routing method (:Routing command must be ROUTE_NONE, ROUTE_PLUG_FLOW, or ROUTE_DIFFUSIVE_WAVE to support transport)",STUB);
  }

  //all fluxes from catchment are routed directly to basin outlet (unless handled in convolution routing as source term)
  if ((!_lateral_via_convol) || (_pModel->GetSubBasin(p)->IsHeadwater())){
    aMout_new[nSegments-1]+=Mlat_new;
  }

  //Mass removal at basin outlet (Mass removal point)
  //-----------------------------------------------------------------
  //_M_loss_rate = reduction_rate * aMout_new[nSegments - 1];
  //aMout_new[nSegments - 1] -= _M_loss_rate;

  //Reservoir Routing - calculate Res_mass, Res_sed_mass
  //-----------------------------------------------------------------
  Res_mass  =0.0;
  ResSedMass=0.0;
  CReservoir *pRes=_pModel->GetSubBasin(p)->GetReservoir();
  if(pRes!=NULL)
  {
    RouteMassInReservoir(p,aMout_new,Res_mass,ResSedMass,Options,tt);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Calculates Res_mass and ResSedMass in basin with reservoir at end of time step
///
/// \param p           [in]  subbasin index
/// \param aMout_new[] [in] Array of mass outflows at downstream end of each segment at end of current timestep [mg/d] [size: nsegs]
/// \param Res_mass    [out] calculated reservoir mass at end of time step
/// \param ResSedMass [out] calculated reservoir mass in sediment at end of time step
/// \param &Options    [in]  Global model options information
/// \param tt          [in]  Time structure
//
void   CConstituentModel::RouteMassInReservoir(const int          p,          // SB index
                                               const double      *aMout_new,  // [mg/d][size: nsegs ]
                                                     double      &Res_mass,   // [mg]
                                                     double      &ResSedMass, // [mg]
                                               const optStruct   &Options,
                                               const time_struct &tt) const
{
  CReservoir* pRes = _pModel->GetSubBasin(p)->GetReservoir();
  int nSegments    = _pModel->GetSubBasin(p)->GetNumSegments();

  double V_new=pRes->GetStorage       ();
  double V_old=pRes->GetOldStorage    ();
  double Q_new=pRes->GetOutflowRate   ()*SEC_PER_DAY;
  double Q_old=pRes->GetOldOutflowRate()*SEC_PER_DAY;

  double decay_coeff=0.0;
  CHydroUnit*   pHRU=_pModel->GetHydroUnit(pRes->GetHRUIndex());
  if(pHRU!=NULL) {
    int     iSW = _pModel->GetStateVarIndex(SURFACE_WATER);
    int     ii  = _pTransModel->GetWaterStorIndexFromSVIndex(iSW);
    decay_coeff = _pTransModel->GetGeochemParam(PAR_DECAY_COEFF,_constit_index,ii,DOESNT_EXIST,pHRU);
  }

  //Explicit solution of Crank-nicolson problem
  //dM/dt=QC_in-QC_out-lambda*M
  //dM/dt=0.5(QC_in^(n+1)-QC_in^(n))-0.5(Q_outC^(n+1)-Q_outC^(n))+M_rain-0.5*lambda*(C^(n+1)+C^n)/V
  double Min= 0.5 * Options.timestep * (aMout_new[nSegments - 1] + _aMout[p][nSegments - 1]);
  double tmp;
  tmp=_aMres[p]*(1.0-0.5*Options.timestep*(Q_old/V_old+decay_coeff*Res_mass));
  tmp+=Min;
  tmp+=_aMresRain[p]*Options.timestep;
  tmp/=(1.0+0.5*Options.timestep*(Q_new/V_new+decay_coeff));
  Res_mass=tmp;

  if((V_old< REAL_SMALL) || (V_new< REAL_SMALL)) { Res_mass= ResSedMass = 0.0; } //handles dried out reservoir/lake
}
//////////////////////////////////////////////////////////////////
/// \brief Sets mass outflow from primary channel and updates flow history
/// \details Also recalculates channel storage (uglier than desired)
/// \remark Called *after* RouteMass routine is called in Solver.cpp, to reset
/// the mass outflow rates for this basin
///
/// \param **aMoutnew     [in] Array of new mass outflows [mg/d] [size:  nsegments x _nConstituents]
/// \param ResMass       [in] new reservoir mass [mg]
/// \param ResMass       [in] new reservoir sediment mass [mg]
/// \param MassOutflow   [out] new mass outflow [mg/d] from last segment or reservoir
/// \param &Options       [in] Global model options information
/// \param initialize     [in] Flag to indicate if flows are to only be initialized
//
void   CConstituentModel::UpdateMassOutflows( const int     p,
                                              const double *aMoutnew,
                                              const double &Mlat_new,
                                              const double &ResMass,
                                              const double &ResSedMass,
                                                    double &MassOutflow,
                                              const optStruct &Options,
                                              const time_struct &tt,
                                              const bool    initialize)
{
  double tstep=Options.timestep;

  CSubBasin *pBasin=_pModel->GetSubBasin(p);

  //Update mass flows
  //------------------------------------------------------
  _aMout_last[p]=_aMout[p][pBasin->GetNumSegments()-1];
  for(int seg=0;seg<pBasin->GetNumSegments();seg++) {
    _aMout[p][seg]=aMoutnew[seg];
  }
  MassOutflow=_aMout[p][pBasin->GetNumSegments()-1];// is now the new mass outflow from the channel


  //Update reservoir concentrations
  //-----------------------------------------------------
  CReservoir *pRes=_pModel->GetSubBasin(p)->GetReservoir();
  if(pRes!=NULL)
  {
    _aMres_last    [p]=_aMres[p];
    _aMres         [p]=ResMass;

    _aMout_res_last[p]=_aMout_res[p];
    _aMout_res     [p]=((pRes->GetOutflowRate()*SEC_PER_DAY)/pRes->GetStorage())*_aMres[p]; //mg/d or MJ/d

    _aMsed_last    [p]=_aMsed[p];
    _aMsed         [p]=ResSedMass;

    MassOutflow=_aMout_res[p];

    // \todo[funct] - properly handle reservoir extraction in mass balance
    // assume inflow is clean, outflow has concentration Cres
    /*_MBres_losses[p]=0.0;
    CTimeSeries *pExpRes->GetExtractionTS();
    if(pEx!=NULL){
    int nn        =(int)((tt.model_time+TIME_CORRECTION)/tstep);//current timestep index
    _MBres_losses[p]+=0.5*max(pEx->GetSampledValue(nn)*SEC_PER_DAY,0.0)*(_aMres[p]/pRes->GetStorage()+_aMres_last[p]/pRes->GetOldStorage())*tstep;
    }*/
  }

  if(initialize) { return; }//entering initial conditions

  //Update channel storage
  //------------------------------------------------------
  double dt=tstep;
  double dM=0.0;

  //mass change from linearly varying upstream inflow over time step
  dM+=0.5*(_aMinHist[p][0]+_aMinHist[p][1])*dt;

  //mass change from linearly varying downstream outflow over time step
  dM-=0.5*(_aMout[p][pBasin->GetNumSegments()-1]+_aMout_last[p])*dt;

  //energy change from loss of mass/energy in-catchment over time step
  dM-=GetCatchmentTransitLosses(p);

  //energy change from loss of mass/energy along reach over time step
  dM-=GetNetReachLosses(p);

  //mass change from lateral inflows
  if (!_lateral_via_convol){//elsewise, in get net reach losses
    dM+=0.5*(Mlat_new+_aMlat_last[p])*dt;
  }

  _channel_storage[p]+=dM;//[mg] or [MJ]

  //Update rivulet storage
  //------------------------------------------------------
  dM=0.0;
  dM+=_aMlatHist[p][0]*dt;//inflow from land surface (integrated over time step)
  dM-=0.5*(Mlat_new+_aMlat_last[p])*dt;//outflow to channel (integrated over time step)

  _rivulet_storage[p]+=dM;//[mg] or [MJ]

  _aMlat_last[p]=Mlat_new;
}
//////////////////////////////////////////////////////////////////
/// \brief returns outflow concentration of constituent c (or temperature) in subbasin p at current point in time
/// \notes only used for reporting; calculations exclusively in terms of mass/energy
///
/// \param p [in] subbasin index
//
double CConstituentModel::GetOutflowConcentration(const int p) const
{
  CReservoir *pRes=_pModel->GetSubBasin(p)->GetReservoir();

  if(pRes==NULL)
  {
    double mg_per_d=_aMout[p][_pModel->GetSubBasin(p)->GetNumSegments()-1];
    double flow    =_pModel->GetSubBasin(p)->GetOutflowRate()*SEC_PER_DAY; //[m3/d]
    if(flow<=0) { return 0.0; }
    double C=mg_per_d/flow/LITER_PER_M3; //[mg/d]/[m3/d]/[L/m3]->[mg/L]
    return C;
  }
  else {
    return _aMres[p]/pRes->GetStorage()/LITER_PER_M3; //[mg/L]
  }
}

//////////////////////////////////////////////////////////////////
/// \brief returns integrated mass outflow of constituent c from subbasin p over current timestep
///
/// \param p [in] subbasin index
/// \param c [in] constituent index
/// \returns integrated mass outflow of constituent c, in mg
//
double CConstituentModel::GetIntegratedMassOutflow(const int p,const double &tstep) const
{
  CReservoir *pRes=_pModel->GetSubBasin(p)->GetReservoir();

  if(pRes==NULL)
  {
    double mgdnew=_aMout[p][_pModel->GetSubBasin(p)->GetNumSegments()-1];
    double mgdold=_aMout_last[p];//[mg/d] or [MJ/d]

    return 0.5*(mgdnew+mgdold)*tstep; //[mg] or [MJ]
  }
  else
  {
    return 0.5*(_aMout_res[p]+_aMout_res_last[p])*tstep; //integrated - [mg]
  }
}
