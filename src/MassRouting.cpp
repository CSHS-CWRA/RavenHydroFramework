/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2020 the Raven Development Team
  ------------------------------------------------------------------
  routing of mass/energy in catchment and channel
  ----------------------------------------------------------------*/

#include "SubBasin.h"
#include "Model.h"
#include "Transport.h"
//////////////////////////////////////////////////////////////////
/// \brief Initialize mass routing variables
///
void CTransportModel::InitializeRoutingVars()
{
  int nSB=pModel->GetNumSubBasins();
  _aMinHist       =new double **[nSB];
  _aMlatHist      =new double **[nSB];
  _aMout          =new double **[nSB];
  _aMout_last     =new double  *[nSB];
  _aMres          =new double  *[nSB];
  _aMres_last     =new double  *[nSB];
  _aMlat_last     =new double  *[nSB];
  _aMout_res      =new double  *[nSB];
  _aMout_res_last =new double  *[nSB];
  _channel_storage=new double  *[nSB];
  _rivulet_storage=new double  *[nSB];
  _aEnthalpySource=new double  *[nSB];
  _aEnthalpyBeta  =new double   [nSB];

  for (int p=0;p<nSB;p++)
  {
    _rivulet_storage[p]=NULL;
    _aMinHist  [p]=new double *[_nConstituents];
    _aMlatHist [p]=new double *[_nConstituents];
    _aMout     [p]=new double *[_nConstituents];
    _aMout_last[p]=new double  [_nConstituents];
    _aMres     [p]=new double  [_nConstituents];
    _aMres_last[p]=new double  [_nConstituents];
    _aMlat_last[p]=new double  [_nConstituents];
    _aMout_res [p]=new double  [_nConstituents];
    _aMout_res_last [p]=new double [_nConstituents];
    _channel_storage[p]=new double [_nConstituents];
    _rivulet_storage[p]=new double [_nConstituents];
    
    ExitGracefullyIf(_rivulet_storage[p]==NULL,"CTransportModel::InitializeRoutingVars(1)",OUT_OF_MEMORY);
    int nMinHist =pModel->GetSubBasin(p)->GetInflowHistorySize();
    int nMlatHist=pModel->GetSubBasin(p)->GetLatHistorySize();
    int nSegments=pModel->GetSubBasin(p)->GetNumSegments();
    for (int c=0;c<_nConstituents;c++)
    {
      _aMout    [p][c]=NULL;
      _aMinHist [p][c]=new double [nMinHist];
      _aMlatHist[p][c]=new double [nMlatHist];
      _aMout    [p][c]=new double [nSegments];
      ExitGracefullyIf(_aMout[p][c]==NULL,"CTransportModel::InitializeRoutingVars(2)",OUT_OF_MEMORY);
      for (int i=0; i<nMinHist;i++ ){_aMinHist [p][c][i]=0.0;}
      for (int i=0; i<nMlatHist;i++){_aMlatHist[p][c][i]=0.0;}
      for (int i=0; i<nSegments;i++){_aMout    [p][c][i]=0.0;}
      _aMout_last[p][c]=0.0;
      _aMres          [p][c]=0.0;
      _aMres_last     [p][c]=0.0;
      _aMlat_last     [p][c]=0.0;
      _aMout_res      [p][c]=0.0;
      _aMout_res_last [p][c]=0.0;
      _channel_storage[p][c]=0.0;
      _rivulet_storage[p][c]=0.0;
    }
    _aEnthalpySource[p]=new double [nMinHist];
    for(int i=0; i<nMinHist;i++) { _aEnthalpySource[p][i]=0.0; }
    _aEnthalpyBeta  [p]=0.0;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Delete mass routing variables
///
void CTransportModel::DeleteRoutingVars()
{
  int nSB=pModel->GetNumSubBasins();
  if (_aMinHist!=NULL){
    for (int p=0;p<nSB;p++)
    {
      for (int c=0;c<_nConstituents;c++){
        delete [] _aMinHist [p][c];
        delete [] _aMlatHist[p][c];
        delete [] _aMout    [p][c];
      }
      delete [] _aMinHist  [p];
      delete [] _aMlatHist [p];
      delete [] _aMout     [p];
      delete [] _aMres     [p];
      delete [] _aMres_last[p];
      delete [] _aMlat_last[p];
      delete [] _aMout_res [p];
      delete [] _aMout_res_last [p];
      delete [] _channel_storage[p];
      delete [] _rivulet_storage[p];
      delete [] _aMout_last[p];
      delete [] _aEnthalpySource[p];
    }
    delete [] _aMinHist;        _aMinHist  =NULL;
    delete [] _aMlatHist;       _aMlatHist =NULL;
    delete [] _aMout;           _aMout     =NULL;
    delete [] _aMres;           _aMres     =NULL;
    delete [] _aMres_last;      _aMres_last=NULL;
    delete [] _aMlat_last;      _aMlat_last=NULL;
    delete [] _aMout_res;       _aMout_res =NULL;
    delete [] _aMout_res_last;  _aMout_res_last =NULL;
    delete [] _channel_storage; _channel_storage=NULL;
    delete [] _rivulet_storage; _rivulet_storage=NULL;
    delete [] _aMout_last;      _aMout_last     =NULL;
    delete [] _aEnthalpySource; _aEnthalpySource=NULL;
    delete [] _aEnthalpyBeta;   _aEnthalpyBeta  =NULL;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Set upstream mass loading to primary channel and updates mass loading history
///
/// \param p    subbasin index
/// \param aMin array of rates of upstream mass/energy loading for the current timestep [mg/d]/[MJ/d] [size: _nConstituents]
//
void   CTransportModel::SetMassInflows    (const int p, const double *aMin)
{
  int nMinHist=pModel->GetSubBasin(p)->GetInflowHistorySize();
  for (int c=0;c<_nConstituents;c++){
    for (int n=nMinHist-1;n>0;n--){
      _aMinHist[p][c][n]=_aMinHist[p][c][n-1];
    }
    _aMinHist[p][c][0]=aMin[c];
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Updates aMinnew, array of mass loadings, to handle fixed concentration/temperature or specified mass inflow conditions
///
/// \param p    subbasin index
/// \param t    end of timestep, in model time
/// \param aMinnew array of rates of upstream mass/energy loading for the current timestep [mg/d]/[MJ/d] [size: _nConstituents]
//
void   CTransportModel::ApplySpecifiedMassInflows(const int p,const double t,double *aMinnew) 
{
  double C;
  long tag1;
  int  tag2;
  long   SBID=pModel->GetSubBasin(p)->GetID();
  double    Q=pModel->GetSubBasin(p)->GetOutflowRate()*SEC_PER_DAY; //[m3/d] Flow at end of time step
  
  /*for(int c=0;c<_nConstituents;c++)
  {
    pTransportSubmodel[c]->ApplySpecifiedMassInflows(p,t,aMinnew[c]);
  }*/

  //Handle additional mass/energy inflows 
  /*for(int i=0; i<_nMassInflowTS; i++) {
    tag1=s_to_l(_pMassInflowTS[i]->GetTag().c_str());
    tag2=s_to_i(_pMassInflowTS[i]->GetTag2().c_str());
    for(int c=0;c<_nConstituents;c++) {
      if(tag1==SBID) {
        if(tag2==c) {
          aMinnew[c]+=_pMassInflowTS[i]->GetValue(t); //mg/d or MJ/d  
          // \todo[funct] - add specified mass inflow time series
        }
      }
    }
  }*/

  //Handle specified concentrations/temperatures of streamflow
  for(int i=0; i<_nSpecFlowConcs; i++) 
  {
    tag1=s_to_l(_pSpecFlowConcs[i]->GetTag().c_str());
    tag2=s_to_i(_pSpecFlowConcs[i]->GetTag2().c_str());
    for(int c=0;c<_nConstituents;c++) 
    {
      if ((tag1==SBID) && (tag2==c)) {
        C=_pSpecFlowConcs[i]->GetValue(t); //mg/L or C
        if(_pConstituents[c]->type==ENTHALPY) { C=ConvertTemperatureToVolumetricEnthalpy(C,0.0); } //MJ/m3 
        else                                  { C*=LITER_PER_M3; } //mg/m3
        aMinnew[c]=Q*C;  //[m3/d]*[mg/m3] or [m3/d]*[MJ/m3] 
      }
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
double   CTransportModel::GetMassAddedFromInflowSources(const int c, const double &t, const double &tstep) const 
{
  double C,Q,Qold;
  long SBID,SBID_down;
  long tag1; int tag2;
  double mass=0; //[mg] or [MJ]
  
  // Handle additional mass/energy inflows 
  /*for(int i=0; i<_nMassInflowTS; i++) {
    tag1=s_to_l(_pMassInflowTS[i]->GetTag ().c_str());
    tag2=s_to_i(_pMassInflowTS[i]->GetTag2().c_str());
    for(int c=0;c<_nConstituents;c++) {
      if ((tag1==SBID) && (tag2==c)) {
          mass+=0.5*(_pMassInflowTS[i]->GetValue(t)+_pMassInflowTS[i]->GetValue(t-tstep)); //mg/d or MJ/d  
        }
      }
    }
  }*/

  // remove mass which is overridden by specified mass inflow
  for(int p=0;p<pModel->GetNumSubBasins();p++) 
  {
    SBID_down=pModel->GetSubBasin(p)->GetDownstreamID();
    for(int i=0; i<_nSpecFlowConcs; i++) {
      tag1=s_to_l(_pSpecFlowConcs[i]->GetTag ().c_str());
      tag2=s_to_i(_pSpecFlowConcs[i]->GetTag2().c_str());
      if ((tag2==c) && (tag1==SBID_down)) {
        mass-=GetIntegratedMassOutflow(p,c,tstep);
      }
    }
  }
  // Handle specified concentrations/temperatures of streamflow
  int spec_c;
  for(int i=0; i<_nSpecFlowConcs; i++) 
  {
    spec_c=s_to_i(_pSpecFlowConcs[i]->GetTag2().c_str());
    SBID  =s_to_l(_pSpecFlowConcs[i]->GetTag ().c_str());
    if(spec_c==c) 
    {
      Q=pModel->GetSubBasinByID(SBID)->GetOutflowRate()*SEC_PER_DAY; //[m3/d] Flow at end of time step
      C=_pSpecFlowConcs[i]->GetValue(t);                             //mg/L or C
      if(_pConstituents[c]->type==ENTHALPY) { C=ConvertTemperatureToVolumetricEnthalpy(C,0.0); } //MJ/m3 
      else                                  { C*=LITER_PER_M3;                                 } //mg/m3
      mass+=0.5*Q*C*tstep;  //[mg] or [MJ]

      Qold=pModel->GetSubBasinByID(SBID)->GetLastOutflowRate()*SEC_PER_DAY;//[m3/d] Flow at start of time step
      C=_pSpecFlowConcs[i]->GetValue(t-tstep);                             //mg/L or C
      if(_pConstituents[c]->type==ENTHALPY) { C=ConvertTemperatureToVolumetricEnthalpy(C,0.0); } //MJ/m3 
      else                                  { C*=LITER_PER_M3;                                 } //mg/m3
      mass+=0.5*Qold*C*tstep;
    }    
  }
  return mass;
}
//////////////////////////////////////////////////////////////////
/// \brief Sets lateral mass/energy loading to primary channel and updates mass/energy loading history
///
/// \param p    subbasin index
/// \param aMin array of lateral mass loading rates for the current timestep [mg/d][MJ/d] [size: _nConstituents]
//
void   CTransportModel::SetLateralInfluxes(const int p, const double *aMlat)
{
  int nMlatHist=pModel->GetSubBasin(p)->GetLatHistorySize();
  for (int c=0;c<_nConstituents;c++){
    for (int n=nMlatHist-1;n>0;n--){
      _aMlatHist[p][c][n]=_aMlatHist[p][c][n-1];
    }
    _aMlatHist[p][c][0]=aMlat[c];
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Creates aMoutnew [m^3/s], an array of point measurements for outflow at downstream end of each river segment
/// \details Array represents measurements at the end of the current timestep. if (catchment_routing==distributed),
/// lateral flow Qlat (e.g., runoff, interflow, or baseflow) is routed as part of channel routing routine otherwise,
/// lateral flow is all routed to most downstream outflow segment . Assumes uniform time step
///
/// \param p    subbasin index
/// \param aMout_new[][] Array of mass outflows at downstream end of each segment at end of current timestep [mg/d] [size: nsegs x _nConstituents]
/// \param &Options [in] Global model options information
//
void   CTransportModel::RouteMass(const int          p,         // SB index
                                  double           **aMout_new, // [mg/d][size: nsegs x _nConstituents]
                                  double            *aRes_mass, // [mg][size: _nConstituents]
                                  const optStruct   &Options,   
                                  const time_struct &tt) const
{
  int c;//,seg;
  const double * aUnitHydro =pModel->GetSubBasin(p)->GetUnitHydrograph();
  const double * aRouteHydro=pModel->GetSubBasin(p)->GetRoutingHydrograph();

  int nMlatHist   =pModel->GetSubBasin(p)->GetLatHistorySize();
  int nMinHist    =pModel->GetSubBasin(p)->GetInflowHistorySize();
  int nSegments   =pModel->GetSubBasin(p)->GetNumSegments();
  double seg_fraction=1.0/(double)(nSegments);

  double *Mlat_new=new double [_nConstituents];

  //==============================================================
  // route from catchment
  //==============================================================
  for (int c=0;c<_nConstituents;c++)
  {
    Mlat_new[c]=0.0;
    for (int n=0;n<nMlatHist;n++){
      Mlat_new[c]+=aUnitHydro[n]*_aMlatHist[p][c][n];
    }
  }
  //==============================================================
  // route along channel
  //==============================================================
  if ((Options.routing==ROUTE_PLUG_FLOW) || (Options.routing==ROUTE_DIFFUSIVE_WAVE))
  {     //Simple convolution
    for (c=0;c<_nConstituents;c++)
    {
      if(_pConstituents[c]->type!=ENTHALPY) {
        aMout_new[nSegments-1][c]=0.0;
        for(int n=0;n<nMinHist;n++) {
          aMout_new[nSegments-1][c]+=aRouteHydro[n]*_aMinHist[p][c][n];
        }
      }
      else 
      {
        //support energy exchange along reach for enthalpy transport
        double term1, term2=0.0;
        double tstep=Options.timestep;
        double beta=_aEnthalpyBeta[p];

        //beta=0.0; //TMP DEBUG
        aMout_new[nSegments-1][c]=0.0;
        //aMinHist and source term both have same indexing; 
        for(int i=0;i<nMinHist;i++) 
        {
          term1=_aMinHist[p][c][i]*exp(-beta*i*tstep);
          term2=0.0;
          for(int j=0;j<i;j++) {
            term2+=_aEnthalpySource[p][i-j]/beta*exp(-beta*(i-j)*tstep)*(1.0-exp(-beta*tstep));
          }
          //term2=0;//TMP DEBUG
          aMout_new[nSegments-1][c]+=aRouteHydro[i]*(term1+term2);
        }
      }
    }
  }
  //==============================================================
  else if (Options.routing==ROUTE_NONE)
  {//In channel routing instantaneous
    for (c=0;c<_nConstituents;c++){
      aMout_new[nSegments-1][c]=_aMinHist[p][c][0];//spits out the mass that just entered
    }
  }
  //==============================================================
  /*else if(Options.routing==ROUTE_HYDROLOGIC){

    ExitGracefully("Transport with ROUTE_HYDROLOGIC",STUB);
  }*/
  //==============================================================
  else{
    ExitGracefully("Unrecognized or unsupported constiuent routing method (:Routing command must be ROUTE_NONE, ROUTE_PLUG_FLOW, or ROUTE_DIFFUSIVE_WAVE to support transport)",STUB);
  }

  if (Options.distrib_lat_inflow==false)
  {//all fluxes from catchment are routed directly to basin outlet
    for (c=0;c<_nConstituents;c++){
      aMout_new[nSegments-1][c]+=Mlat_new[c];
    }
  }
  else{//only last segments worth
    for (c=0;c<_nConstituents;c++){
      aMout_new[nSegments-1][c]+=Mlat_new[c]*seg_fraction;
    }
  }

  //Reservoir Routing
  //-----------------------------------------------------------------
  CReservoir *pRes=pModel->GetSubBasin(p)->GetReservoir();
  if (pRes!=NULL)
  {
    double decay_coeff(0.0);
    int iSW=pModel->GetStateVarIndex(SURFACE_WATER);
    CHydroUnit *pHRU=pModel->GetHydroUnit(pRes->GetHRUIndex());
    
    double tmp;
    double V_new=pRes->GetStorage();
    double V_old=pRes->GetOldStorage();
    double Q_new=pRes->GetOutflowRate()   *SEC_PER_DAY;
    double Q_old=pRes->GetOldOutflowRate()*SEC_PER_DAY;

    for(c=0;c<_nConstituents;c++)
    {
      decay_coeff=0.0;
      if(pHRU!=NULL){ decay_coeff =GetDecayCoefficient(c,pHRU,iSW)*aRes_mass[c]; }

      //Explicit solution of Crank-nicolson problem
      //dM/dt=QC_in-QC_out-lambda*M
      //dM/dt=0.5(QC_in^(n+1)-QC_in^(n))-0.5(Q_outC^(n+1)-Q_outC^(n))-0.5*lambda*(C^(n+1)+C^n)/V

      tmp=_aMres[p][c]*(1.0-0.5*Options.timestep*(Q_old/V_old+decay_coeff));
      tmp+=+0.5*Options.timestep*(aMout_new[nSegments-1][c]+_aMout[p][c][nSegments-1]);//inflow is outflow from channel
      tmp/=(1.0+0.5*Options.timestep*(Q_new/V_new+decay_coeff));
      aRes_mass[c]=tmp;

      if((V_old<=0.0) || (V_new<=0.0)){aRes_mass[c]=0.0;} //handles dried out reservoir/lake
    }
  }
  else
  {
    for(c=0;c<_nConstituents;c++){aRes_mass[c]=0.0;}
  }

  delete [] Mlat_new;

  return;
}
//////////////////////////////////////////////////////////////////
/// \brief Sets mass outflow from primary channel and updates flow history
/// \details Also recalculates channel storage (uglier than desired)
/// \remark Called *after* RouteMass routine is called in Solver.cpp, to reset
/// the mass outflow rates for this basin
///
/// \param **aMoutnew     [in] Array of new mass outflows [mg/d] [size:  nsegments x _nConstituents]
/// \param *aResMass      [in] Array of new reservoir masses [mg] [size:  _nConstituents]
/// \param *aMassOutflow [out] Array of new mass outflows [mg/d] from last segment or reservoir [size:  _nConstituents]
/// \param &Options       [in] Global model options information
/// \param initialize     [in] Flag to indicate if flows are to only be initialized
//
void   CTransportModel::UpdateMassOutflows(const int p,  double **aMoutnew, 
                                                          double *aResMass, 
                                                          double *aMassOutflow,
                                                 const optStruct &Options,
                                               const time_struct &tt,
                                                            bool initialize)
{
  double tstep=Options.timestep;

  CSubBasin *pBasin=pModel->GetSubBasin(p);

  for (int c=0;c<_nConstituents;c++)
  {
    //Update mass flows
    //------------------------------------------------------
    _aMout_last[p][c]=_aMout[p][c][pBasin->GetNumSegments()-1];
    for (int seg=0;seg<pBasin->GetNumSegments();seg++){
      _aMout[p][c][seg]=aMoutnew[seg][c];
    }
    aMassOutflow[c]=_aMout[p][c][pBasin->GetNumSegments()-1];// is now the new mass outflow from the channel

    //Update reservoir concentrations
    //-----------------------------------------------------
    CReservoir *pRes=pModel->GetSubBasin(p)->GetReservoir();
    if(pRes!=NULL)
    {
      _aMres_last    [p][c]=_aMres[p][c];
      _aMres         [p][c]=aResMass[c];
      _aMout_res_last[p][c]=_aMout_res[p][c];
      _aMout_res     [p][c]=(pRes->GetOutflowRate()*SEC_PER_DAY/pRes->GetStorage())*_aMres[p][c]; //mg/d

      aMassOutflow[c]=_aMout_res[p][c];

      // \todo[funct] - properly handle reservoir extraction in mass balance
      // assume inflow is clean, outflow has concentration Cres
      /*_MBres_losses[p][c]=0.0;
      CTimeSeries *pExpRes->GetExtractionTS();
      if(pEx!=NULL){
        int nn        =(int)((tt.model_time+TIME_CORRECTION)/tstep);//current timestep index
        _MBres_losses[p][c]+=0.5*max(pEx->GetSampledValue(nn)*SEC_PER_DAY,0.0)*(_aMres[p][c]/pRes->GetStorage()+_aMres_last[p][c]/pRes->GetOldStorage())*tstep;
      }*/
    }
    
    if (initialize){return;}//entering initial conditions

    //Update channel storage
    //------------------------------------------------------
    double dt=tstep;
    double dM=0.0;
    double Mlat_new(0.0);
    for (int n=0;n<pBasin->GetLatHistorySize();n++){
      Mlat_new+=pBasin->GetUnitHydrograph()[n]*_aMlatHist[p][c][n];
    }

    //mass change from linearly varying upstream inflow over time step
    dM+=0.5*(_aMinHist[p][c][0]+_aMinHist[p][c][1])*dt;

    //mass change from linearly varying downstream outflow over time step
    dM-=0.5*(_aMout[p][c][pBasin->GetNumSegments()-1]+_aMout_last[p][c])*dt;

    //energy change from loss of energy along reach over time step
    //dM-=GetEnergyLossesFromReach(); 

    //mass change from lateral inflows
    dM+=0.5*(Mlat_new+_aMlat_last[p][c])*dt;

    _channel_storage[p][c]+=dM;//[mg]

    //Update rivulet storage
    //------------------------------------------------------
    dM=0.0;
    dM+=_aMlatHist[p][c][0]*dt;//inflow from land surface (integrated over time step)
    dM-=0.5*(Mlat_new+_aMlat_last[p][c])*dt;//outflow to channel (integrated over time step)

    _rivulet_storage[p][c]+=dM;//[mg] or [MJ]

    _aMlat_last[p][c]=Mlat_new;

    //Update energy source terms for enthalpy transport
    //------------------------------------------------------
    if(_pConstituents[c]->type==ENTHALPY) {
      UpdateReachEnergySourceTerms(p);
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief returns outflow concentration of constituent c (or temperature) in subbasin p at current point in time
/// \notes only used for reporting; calculations exclusively in terms of mass/energy
/// 
/// \param p [in] subbasin index 
/// \param c [in] constituent index
//
double CTransportModel::GetOutflowConcentration(const int p, const int c) const
{
  CReservoir *pRes=pModel->GetSubBasin(p)->GetReservoir();

  //Chemical concentration [mg/L] --------------------------------------------
  if(_pConstituents[c]->type!=ENTHALPY) {
    if(pRes==NULL)
    {
      double mg_per_d=_aMout[p][c][pModel->GetSubBasin(p)->GetNumSegments()-1];
      double flow=pModel->GetSubBasin(p)->GetOutflowRate();
      if(flow<=0) { return 0.0; }
      double C=mg_per_d/flow/SEC_PER_DAY/LITER_PER_M3; //[mg/L]
      return C;
    }
    else {
      return _aMres[p][c]/pRes->GetStorage()/LITER_PER_M3; //[mg/L]
    }
  }
  //Temperatures [C] ---------------------------------------------------------
  else { 
    if(pRes==NULL)
    {
      double MJ_per_d=_aMout[p][c][pModel->GetSubBasin(p)->GetNumSegments()-1];
      double flow    =pModel->GetSubBasin(p)->GetOutflowRate();
      if(flow<=0) { return 0.0; }
      //hv=MJ/d / M3/d = [MJ/m3]
      double T=ConvertVolumetricEnthalpyToTemperature(MJ_per_d/(flow*SEC_PER_DAY)); //[C]

      return max(T,0.0); //Flowing water can't go below zero. HOWEVER, this may not conserve energy balance, as it is possible for the enthalpy to dip down below frozen conditions

      //if (T<0){cout<<"     Weird: "<<MJ_per_d<<" "<<flow<<" "<<T<<" ***************"<<endl;}
      //else {   cout<<" NOT Weird: "<<MJ_per_d<<" "<<flow<<" "<<T<<endl;}
      //return T;
    }
    else {
      return ConvertVolumetricEnthalpyToTemperature(_aMres[p][c]/pRes->GetStorage()); //[C]
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief returns integrated mass outflow of constituent c from subbasin p over current timestep
///
/// \param p [in] subbasin index 
/// \param c [in] constituent index
/// \returns integrated mass outflow of constituent c, in mg
//
double CTransportModel::GetIntegratedMassOutflow(const int p, const int c,const double &tstep) const
{
  CReservoir *pRes=pModel->GetSubBasin(p)->GetReservoir();
  
  if(pRes==NULL)
  {
    double mgdnew=_aMout     [p][c][pModel->GetSubBasin(p)->GetNumSegments()-1];
    double mgdold=_aMout_last[p][c];//[mg/d] or [MJ/d]

    return 0.5*(mgdnew+mgdold)*tstep; //[mg] or [MJ]
  }
  else
  {
    return 0.5*(_aMout_res[p][c]+_aMout_res_last[p][c])*tstep; //integrated - [mg]
  }
}
