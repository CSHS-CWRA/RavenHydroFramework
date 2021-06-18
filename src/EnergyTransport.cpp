/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2021 the Raven Development Team
----------------------------------------------------------------
CEnthalpyModel routines related to energy/enthalpy transport
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "EnergyTransport.h"

string FilenamePrepare(string filebase,const optStruct& Options); //Defined in StandardOutput.cpp

//////////////////////////////////////////////////////////////////
/// \brief enthalpy model constructor
//
CEnthalpyModel::CEnthalpyModel(CModel *pMod,CTransportModel *pTMod,string name,const int c) 
               :CConstituentModel(pMod,pTMod,name,ENTHALPY,false,c) 
{
  _aEnthalpyBeta=NULL;
  _aEnthalpySource=NULL;
}

//////////////////////////////////////////////////////////////////
/// \brief enthalpy model destructor
//
CEnthalpyModel::~CEnthalpyModel() 
{
  delete [] _aEnthalpyBeta;
  for(int p=0;p<_pModel->GetNumSubBasins();p++) { delete[] _aEnthalpySource; } delete [] _aEnthalpySource;
}

//////////////////////////////////////////////////////////////////
/// \brief Calculate temperature for reporting (converts to degrees C) (equivalent to 'concentration' for enthalpy)
/// \param M [in] energy in [MJ/m2]
/// \param V [in] volume in mm
//
double CEnthalpyModel::CalculateConcentration(const double &M,const double &V) const 
{
  if(fabs(V)>1e-6) {     
    return ConvertVolumetricEnthalpyToTemperature(M/V*MM_PER_METER);  //[MJ/m3]->[C]
  }
  return 0.0;// JRC: should this default to zero? or NA?
}

//////////////////////////////////////////////////////////////////
/// \brief returns outflow temperature from reach p (equivalent to 'concentration' for enthalpy)
/// \param p [in] subbasin index
/// \returns outflow temperature [C]
//
double CEnthalpyModel::GetOutflowConcentration(const int p) const
{
  CReservoir *pRes=_pModel->GetSubBasin(p)->GetReservoir();

  if(pRes==NULL)
  {
    double MJ_per_d=_aMout[p][_pModel->GetSubBasin(p)->GetNumSegments()-1];
    double flow    =_pModel->GetSubBasin(p)->GetOutflowRate();
    double hv      =MJ_per_d/(flow*SEC_PER_DAY); //[MJ/d] / [M3/d] = [MJ/m3]

    if(flow<=0) { return 0.0; }

    return ConvertVolumetricEnthalpyToTemperature(hv); //[C]

    //if (T<0){cout<<"     Weird: "<<MJ_per_d<<" "<<flow<<" "<<T<<" ***************"<<endl;}
    //else {   cout<<" NOT Weird: "<<MJ_per_d<<" "<<flow<<" "<<T<<endl;}
    //return T;
  }
  else {
    return ConvertVolumetricEnthalpyToTemperature(_aMres[p]/pRes->GetStorage()); 
  }
}

//////////////////////////////////////////////////////////////////
/// \brief returns ice fraction of routed water in subbasin p at current point in time
/// \notes only used for reporting; calculations exclusively in terms of mass/energy
/// \param p [in] subbasin index 
//
double CEnthalpyModel::GetOutflowIceFraction(const int p) const
{
  CReservoir* pRes=_pModel->GetSubBasin(p)->GetReservoir();

  if(pRes==NULL)
  {
    double MJ_per_d=_aMout[p][_pModel->GetSubBasin(p)->GetNumSegments()-1];
    double flow    =_pModel->GetSubBasin(p)->GetOutflowRate();
    double hv      =MJ_per_d/(flow*SEC_PER_DAY); //[MJ/d] / [M3/d] = [MJ/m3]
    if(flow<=0) { return 0.0; }

    return ConvertVolumetricEnthalpyToIceContent(hv); 
    
  }
  else {
     return ConvertVolumetricEnthalpyToIceContent(_aMres[p]/pRes->GetStorage()); //[C]
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns water temperature (in degC) of storage unit indexed by state variable index iWater
/// \param state_vars - vector of current state variables 
/// \param iWater - state variable index of water compartment of interest
//
double CEnthalpyModel::GetWaterTemperature(const double *state_vars,const int iWater) const
{
  if(iWater==DOESNT_EXIST) { return 0.0; }
 
  double Hv,stor,hv(0.0);
  int    m,iEnth;
  m    =_pTransModel->GetLayerIndex(_constit_index,iWater); //layer index of water compartment
  iEnth=_pModel->GetStateVarIndex(CONSTITUENT,m); //global index of water compartment enthalpy

  Hv   =state_vars[iEnth]; //enthalpy, [MJ/m2]
  stor =state_vars[iWater]; //water storage [mm]
  if(stor>PRETTY_SMALL) {
    hv   =Hv/(stor/MM_PER_METER); //volumetric specific enthalpy [MJ/m3]
  }
  return ConvertVolumetricEnthalpyToTemperature(hv);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns ice content [0..1] of storage unit indexed by state variable index iWater
/// \param state_vars - vector of current state variables 
/// \param iWater - state variable index of water compartment of interest
//
double CEnthalpyModel::GetIceContent(const double* state_vars,const int iWater) const
{
  if(iWater==DOESNT_EXIST) { return 0.0; }

  double Hv,stor,hv(0.0);
  int    m,iEnth;
  m    =_pTransModel->GetLayerIndex(_constit_index,iWater); //layer index of water compartment
  iEnth=_pModel->GetStateVarIndex(CONSTITUENT,m); //global index of water compartment enthalpy

  Hv   =state_vars[iEnth]; //enthalpy, [MJ/m2]
  stor =state_vars[iWater]; //water storage [mm]
  if(stor>PRETTY_SMALL) {
    hv   =Hv/(stor/MM_PER_METER); //volumetric specific enthalpy [MJ/m3]
  }
  return ConvertVolumetricEnthalpyToIceContent(hv);
}

//////////////////////////////////////////////////////////////////
/// \brief Calculate volumetric enthalpy for Dirichlet condition in storage compartment 
/// \param pHRU - pointer to the HRU of interest
/// \param Cs - dirichlet temperature (or special flag -9999)
//  \returns volumetric enthalpy, h, in MJ/mm-m2
//
double CEnthalpyModel::GetDirichletEnthalpy(const CHydroUnit* pHRU,const double& Cs) const
{
  if(Cs!=DIRICHLET_TEMP) 
  {
    double pctFroz=0.0; //assumes liquid water for flows (reasonable)
    return ConvertTemperatureToVolumetricEnthalpy(Cs,pctFroz)/MM_PER_METER;    //[C]->[MJ/m3]->[MJ/mm-m2]
  }
  else // Special temperature condition - forces temp=air temp 
  { 
    double tmp,snowfrac(0.0);
    double Tair    =pHRU->GetForcingFunctions()->temp_ave;
    if(Tair<FREEZING_TEMP) {
      snowfrac=pHRU->GetForcingFunctions()->snow_frac;
    }
    tmp=ConvertTemperatureToVolumetricEnthalpy(Tair,snowfrac)/MM_PER_METER;   //[C]->[MJ/m3]->[MJ/mm-m2]
    return max(tmp,0.0); //precip @ 0 degrees
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns watershed-wide latent heat flux determined from AET
//
double CEnthalpyModel::GetAvgLatentHeatFlux() const
{
  double area=_pModel->GetWatershedArea()*M2_PER_KM2;

  int    iAET=_pModel->GetStateVarIndex(AET);
  double  AET=_pModel->GetAvgStateVar(iAET)/MM_PER_METER;

  //int   iSubl=pModel->GetStateVarIndex(SUBLIMATED); // \todo[funct] support actual sublimation as state variable
  //double Subl=pModel->GetAvgStateVar  (iSubl)/MM_PER_METER;

  return AET*LH_VAPOR*DENSITY_WATER*area; //+Subl*LH_SUBLIM*DENSITY_WATER*area;
}

//////////////////////////////////////////////////////////////////
/// \brief returns heat per unit area generated from friction in reach [MJ/m2/d] 
/// \param Q     [in] flow rate [m3/s]
/// \param slope [in] bed slope [m/m]
/// \param perim [in] wetted perimeter [m]
//
double CEnthalpyModel::GetReachFrictionHeat(const double &Q,const double &slope,const double &perim) const
{
  //from Theurer et al 1984, US Fish and Wildlife Serviec FWS/OBS-85/15, as reported in MacDonald, Boon, and Byrne 2014
  if (perim<PRETTY_SMALL){return 0.0;}
  return 9805*Q/perim*slope*WATT_TO_MJ_PER_D;
}

//////////////////////////////////////////////////////////////////
/// \brief Updates source terms for energy balance on subbasin reaches each time step
/// \details both _aEnthalpyBeta and _aEnthalpySource (and its history) are generated here
/// \param p  subbasin index
/// \notes this is the meat behind the in-stream thermal routing algorithm
//
void   CEnthalpyModel::UpdateReachEnergySourceTerms(const int p)
{
  double tstep=_pModel->GetOptStruct()->timestep;

  const CSubBasin *pBasin=_pModel->GetSubBasin(p);

  if(pBasin->IsHeadwater()) { return; } //no reach, no need

  const CHydroUnit  *pHRU=_pModel->GetHydroUnit(pBasin->GetReachHRUIndex());
  int                iAET=_pModel->GetStateVarIndex(AET);

  int    iGW      =_pModel->GetStateVarIndex(SOIL,2);
  double temp_GW  =0.0;//pBasin->GetAvgConcentration(iGW); 
  
  double temp_air =pHRU->GetForcingFunctions()->temp_ave;           //[C]
  double SW       =pHRU->GetForcingFunctions()->SW_radia_net;       //[MJ/m2/d]
  double LW       =pHRU->GetForcingFunctions()->LW_radia_net;       //[MJ/m2/d]
  double LW_in    =pHRU->GetForcingFunctions()->LW_incoming;        //[MJ/m2/d]
  double AET      =pHRU->GetStateVarValue(iAET)/MM_PER_METER/tstep; //[m/d]

  double hstar    =pBasin->GetConvectionCoeff(); //[MJ/m2/d/K]  
  double qmix     =pBasin->GetHyporheicFlux();   //[m/d]
  double bed_ratio=pBasin->GetTopWidth()/max(pBasin->GetWettedPerimeter(),0.001);

  double dbar     =max(pBasin->GetRiverDepth(),0.001);

  double Qf       =GetReachFrictionHeat(pBasin->GetOutflowRate(),pBasin->GetBedslope(),pBasin->GetWettedPerimeter());//[MJ/m2/d]  

  double S(0.0);                            //source term [MJ/m3/d]
  S+=(SW+LW)/dbar;                          //net incoming energy term
  S-=AET*DENSITY_WATER*LH_VAPOR/dbar;       //latent heat flux term
  S+=Qf/dbar;                               //friction-generated heating term
  S+=hstar*temp_air/dbar;                   //convection with atmosphere
  S+=qmix*HCP_WATER*bed_ratio*temp_GW/dbar; //hyporheic mixing

  _aEnthalpyBeta[p]=(hstar/dbar + qmix/dbar*HCP_WATER*bed_ratio)/HCP_WATER;

  int nMinHist=_pModel->GetSubBasin(p)->GetInflowHistorySize();
  for(int n=nMinHist-1;n>0;n--) {
    _aEnthalpySource[p][n]=_aEnthalpySource[p][n-1];
  }
  _aEnthalpySource[p][0]=S;

  //cout<<" S terms: "<<(SW+LW)<<"|"<<AET*DENSITY_WATER*LH_VAPOR<<"|"<<Qf<<"|"<<hstar*temp_air<<"|"<<qmix*HCP_WATER*DENSITY_WATER*bed_ratio*temp_GW<<endl;
  //cout<<"ENTH BETA: "<<_aEnthalpyBeta[p]<<" "<<hstar/HCP_WATER<<" "<< qmix*bed_ratio<<endl;
}


//////////////////////////////////////////////////////////////////
/// \brief Calculates individual energy gain terms from reach over current time step
/// \param p    subbasin index
/// \param Q_sens [out] energy gain from convective/sensible heating [MJ/d]
/// \param Q_lat  [out] energy gain from latent heat exchange [MJ/d]
/// \param Q_GW   [out] energy gain from groundwater mixing [MJ/d]
/// \param Q_rad  [out] energy gain from radiation of surface [MJ/d]
/// \param Q_fric [out] energy gain from friction [MJ/d]
/// \returns total energy lost from reach over current time step [MJ]
//
double CEnthalpyModel::GetEnergyLossesFromReach(const int p,double &Q_sens,double &Q_lat,double &Q_GW,double &Q_rad,double &Q_fric) const
{
 
  double tstep=_pModel->GetOptStruct()->timestep;

  const CSubBasin *pBasin=_pModel->GetSubBasin(p);
  
  if(pBasin->IsHeadwater())
  {
    Q_sens=Q_lat=Q_GW=Q_rad=Q_fric=0.0;
    return 0.0;
  }

  const CHydroUnit* pHRU=_pModel->GetHydroUnit(pBasin->GetReachHRUIndex());
  int                iAET=_pModel->GetStateVarIndex(AET);

  int    iGW      =_pModel->GetStateVarIndex(SOIL,2);
  double temp_GW  =0.0;//pBasin->GetAvgConcentration(iGW); 

  double temp_air =pHRU->GetForcingFunctions()->temp_ave;           //[C]
  double SW       =pHRU->GetForcingFunctions()->SW_radia_net;       //[MJ/m2/d]
  double LW       =pHRU->GetForcingFunctions()->LW_radia_net;       //[MJ/m2/d]
  //double LW       =STEFAN_BOLTZ*EMISS*pow(T-Ta,4); T=GetOutflowConcentration(p);
  double LW_in    =pHRU->GetForcingFunctions()->LW_incoming;        //[MJ/m2/d]
  
  double AET      =pHRU->GetStateVarValue(iAET)/MM_PER_METER/tstep; //[m/d]

  double hstar    =pBasin->GetConvectionCoeff(); //[MJ/m2/d/K]  
  double qmix     =pBasin->GetHyporheicFlux();   //[m/d]
  double bed_ratio=pBasin->GetTopWidth()/max(pBasin->GetWettedPerimeter(),0.001);

  double dbar     =max(pBasin->GetRiverDepth(),0.001);

  double Qf       =GetReachFrictionHeat(pBasin->GetOutflowRate(),pBasin->GetBedslope(),pBasin->GetWettedPerimeter());//[MJ/m2/d]  

  int nMinHist              =pBasin->GetInflowHistorySize();
  const double * aRouteHydro=pBasin->GetRoutingHydrograph();
  const double * aQin       =pBasin->GetInflowHistory();

  double *zk=new double[nMinHist]; //percentage of flows in zone k
  double *Ik=new double[nMinHist];
  
  //calculate mean residence time
  double tr_mean=0.0;
  for(int i=0;i<nMinHist;i++) {
    tr_mean+=aRouteHydro[i]*i*tstep;
  }

  //calculate flow linked to each zone k, zk[k]
  double zsum=0;
  for(int k=1;k<nMinHist;k++) {
    zk[k]=0.0;
    for(int i=k;i<nMinHist;i++) {
      zk[k]+=aQin[k]*aRouteHydro[i]*(tr_mean/i/tstep)*SEC_PER_DAY; // [m3/d]
    }
    zsum+=zk[k];
  }

  //calculate integral term I_m^n [degC*d] (see paper)
  double hin;
  double beta=_aEnthalpyBeta[p];
  beta=max(beta,1e-9); //to avoid divide by zero error
  double gamma=(1.0-exp(-beta*tstep));
  Ik[0]=0.0;
  for(int m=1;m<nMinHist;m++)
  {
    hin=(_aMinHist[p][m-1]/aQin[m-1]);
    if(_aMinHist[p][m-1]<PRETTY_SMALL) { hin=0.0; }

    Ik[m]=-hin/beta*gamma*exp(-beta*(m-1)*tstep);
    for(int j=1;j<m;j++)
    {
      Ik[m]+=gamma*gamma/beta/beta*_aEnthalpySource[p][j-1]*exp(-beta*(m-j-1)*tstep);
    }
    Ik[m]+=tstep/beta     *_aEnthalpySource[p][0];
    Ik[m]+=gamma/beta/beta*_aEnthalpySource[p][0];

    Ik[m]*=1.0/HCP_WATER/DENSITY_WATER;
  }

  Q_sens = Q_lat = Q_GW = Q_rad = Q_fric =0.0; //incoming radiation [MJ/d]

  double kprime=qmix/dbar*HCP_WATER*bed_ratio;

  for(int m=1;m<nMinHist;m++)
  {
    Q_sens+=zk[m]*(hstar/dbar)*(temp_air*tstep-Ik[m]);//[m3/d]*[MJ/m2/d/K]*[1/m]*[K]*[d]=[MJ/d]

    Q_GW  +=zk[m]*kprime      *(temp_GW *tstep-Ik[m]);

    Q_rad +=zk[m]*(SW+LW)/dbar*tstep;                //[m3/d]*[MJ/m2/d]*[1/m]*[d]=[MJ/d]

    Q_lat -=zk[m]*AET*DENSITY_WATER*LH_VAPOR/dbar*tstep; //[m3/d]*[m/d]*[kg/m3]*[MJ/kg]*[1/m]*[d] = [MJ/d]

    Q_fric+=zk[m]*Qf/dbar*tstep;

   // cout<<" Losses: "<<k<<": "<<Q_sens<<" "<<Q_GW<<" "<<Q_rad<<" "<<Q_lat<<" "<<Q_fric<< " "<<zk[k]*kprime*(temp_GW*tstep-Ik[k])<<" "<<Ik[k]<<" "<<zk[k]<<endl;
   // cout<<"mean res time: "<<tr_mean<<endl;
   // cout<<" "<<hstar<<" "<< dbar<<" "<<temp_air*tstep<<" "<<zk[k]<<endl;
   //  cout<<"z["<<k<<"]: "<<zk[k]/SEC_PER_DAY<<"("<<zk[k]<<") "<<aQin[k]<<" "<<aRouteHydro[k]<<endl;
  }

  delete[] Ik;
  delete[] zk;

  return -(Q_sens+Q_GW+Q_rad+tstep*Q_lat+Q_fric)*tstep; //total energy lost [MJ] //TMP DEBUG - WHY IS THE tstep needed HERE before Q_lat!!?? Is it the water lost from the system??
}

//////////////////////////////////////////////////////////////////
/// \brief returns total energy lost from subbasin reach over current time step [MJ]
/// \param p    subbasin index
/// \returns total energy lost from subbasin reach over current time step [MJ]
//
double CEnthalpyModel::GetNetReachLosses(const int p) const 
{
  double Q_sens,Q_lat,Q_GW,Q_rad,Q_fric;
  return GetEnergyLossesFromReach(p,Q_sens,Q_lat,Q_GW,Q_rad,Q_fric);
}

//////////////////////////////////////////////////////////////////
/// \brief initializes enthalpy model
//
void CEnthalpyModel::Initialize(const optStruct& Options)
{
  // Initialize base class members
  //--------------------------------------------------------------------
  CConstituentModel::Initialize(Options);

  // Allocate memory
  //--------------------------------------------------------------------
  int nSB=_pModel->GetNumSubBasins();
  _aEnthalpySource=new double  *[nSB];
  _aEnthalpyBeta  =new double[nSB];
  for(int p=0;p<nSB;p++) {
    int nMinHist =_pModel->GetSubBasin(p)->GetInflowHistorySize();
    _aEnthalpySource[p]=new double[nMinHist];
    for(int i=0; i<nMinHist;i++) { _aEnthalpySource[p][i]=0.0; }
    _aEnthalpyBeta[p]=0.0;
  }

  // initialize stream temperatures if init_stream_temp is given
  //--------------------------------------------------------------------
  double hv;
  if(CGlobalParams::GetParams()->init_stream_temp>0.0)
  {
    const CSubBasin *pBasin;

    hv=ConvertTemperatureToVolumetricEnthalpy(CGlobalParams::GetParams()->init_stream_temp,0.0);
    for(int p=0;p<_pModel->GetNumSubBasins();p++) 
    {
      pBasin=_pModel->GetSubBasin(p);
      const double *aQin=pBasin->GetInflowHistory();
      for(int i=0; i<pBasin->GetNumSegments(); i++)
      {
        _aMout[p][i]=pBasin->GetOutflowRate()*SEC_PER_DAY*hv; //not really - requires outflow rate from all segments. Don't have access to this.
      }
      _aMout_last[p]=_aMout[p][0];
      for(int i=0; i<pBasin->GetInflowHistorySize(); i++)
      {
        _aMinHist[p][i]=aQin[i]*SEC_PER_DAY*hv;
      }

      //SOME ISSUE HERE WITH THIS INTIALIZATION - WHY?
      //NO PROBLEMS WITH ABOVE MASS FLOW RATE INITIALIZATION, ONLY EXTRA MASS BALANCE TERMS
      _channel_storage[p]=pBasin->GetChannelStorage()*hv/1.717;
      _initial_mass+=_channel_storage[p];
//      cout<<" init storage!"<<_initial_mass;
    }
  }

  // QA/QC - check for invalid reach HRU index
  //--------------------------------------------------------------------
  const CSubBasin* pBasin;
  for(int p=0;p<_pModel->GetNumSubBasins();p++)
  {
    pBasin=_pModel->GetSubBasin(p);
    if((!pBasin->IsHeadwater()) && (pBasin->GetReachHRUIndex()==DOESNT_EXIST)) {
      ExitGracefully("CEnthalpyModel::Initialize: non-headwater subbasin missing reach HRU index for temperature simulation",BAD_DATA);
    }
  }
  //  int  m=_pTransModel->GetLayerIndexFromName2("!TEMPERATURE_SOIL",2);
  //  ExitGracefullyIf(m==DOESNT_EXIST,"UpdateReachEnergySourceTerms: must have 3 layers",BAD_DATA_WARN);

  // update beta & source matrices
  //--------------------------------------------------------------------
  for(int p=0;p<nSB;p++) {
    UpdateReachEnergySourceTerms(p);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Applies special convolution with source/sink terms - analytical solution to Lagrangian heat transport problem
/// \param p      [in]  subbasin index
/// \param aRouteHydro [in] routing unit hydrograph [size: nMinHist] [-] 
/// \param aQinHist [in] reach inflow history [size: nMinHist] [m3/s]
/// \param aMinHist [in] energy inflow history [size: nMinHist] [MJ/d]
/// \param nSegments [in] number of reach segments
/// \param nMinHist [in] size of inflow history array
/// \param tstep [in] model time step [d]
/// \param aMout_new [out] array of model energy outflows [size: nSegments]
/// \brief calculates aMout_new, the energy outlow from the reach [MJ/d]
//
void    CEnthalpyModel::ApplyConvolutionRouting(const int p,const double *aRouteHydro,const double *aQinHist,const double *aMinHist,
                                                const int nSegments,const int nMinHist,const double &tstep,double *aMout_new) const
{
  double term1,term2=0.0;
  double beta=_aEnthalpyBeta[p];

  beta=max(beta,1e-9); //to avoid divide by zero error

  aMout_new[nSegments-1]=0.0;
  //aMinHist and source term both have same indexing;
  for(int i=0;i<nMinHist;i++)
  {
    term1=aMinHist[i]*exp(-beta*i*tstep);
    term2=0.0;
    for(int j=1;j<=i;j++) {
      term2+=_aEnthalpySource[p][i-j]/beta*exp(-beta*(i-j)*tstep)*(1.0-exp(-beta*tstep));
    }
    term2*=(aQinHist[i]*SEC_PER_DAY);
    aMout_new[nSegments-1]+=aRouteHydro[i]*(term1+term2);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Applies special convolution with source/sink terms - analytical solution to Lagrangian heat transport problem
/// \param p      [in]  subbasin index
/// \param aMoutnew [in] array of new mass outflows
/// \param ResMass [in] new reservoir enthalpy [MJ]
/// \param MAssOutflows [in] new mass outflows
/// \param Options [in] options structure
/// \param tt [in] time structure
/// \param initialize [in] true if initializing
//
void   CEnthalpyModel::UpdateMassOutflows(const int p,
                                          double *aMoutnew,
                                          double  &ResMass,
                                          double  &MassOutflow,
                                          const   optStruct &Options,
                                          const   time_struct &tt,
                                          bool    initialize)
{
  CConstituentModel::UpdateMassOutflows(p,aMoutnew,ResMass,MassOutflow,Options,tt,initialize);//call base class
  UpdateReachEnergySourceTerms(p);
}
//////////////////////////////////////////////////////////////////
/// \brief Write transport output file headers in .tb0 format
/// \details Called prior to simulation (but after initialization) from CModel::Initialize()
/// \param &Options [in] Global model options information
//
void CEnthalpyModel::WriteEnsimOutputFileHeaders(const optStruct &Options)
{
  this->WriteOutputFileHeaders(Options); //Ensim format not supported by enthalpy model
}

//////////////////////////////////////////////////////////////////
/// \brief Writes minor transport output to ensim formatted file at the end of each timestep (or multiple thereof)
/// \note only thing this modifies should be output streams; called from CModel::WriteMinorOutput()
/// \param &Options [in] Global model options information
/// \param &tt [in] Local (model) time at the end of the pertinent time step
//
void CEnthalpyModel::WriteEnsimMinorOutput(const optStruct &Options,const time_struct &tt)
{
  this->WriteMinorOutput(Options,tt);
}

//////////////////////////////////////////////////////////////////
/// \brief Write transport output file headers in .tb0 format
/// \details Called prior to simulation (but after initialization) from CModel::Initialize()
/// \param &Options [in] Global model options information
//
void CEnthalpyModel::WriteOutputFileHeaders(const optStruct& Options) 
{
  //StreamTemperatures.csv and Temperatures.csv (similar to concentration output files)
  CConstituentModel::WriteOutputFileHeaders(Options);

  //StreamReachEnergyBalances.csv
  //--------------------------------------------------------------------
  string filename ="StreamReachEnergyBalances.csv"; 
  filename=FilenamePrepare(filename,Options);

  _STREAMOUT.open(filename.c_str());
  if(_STREAMOUT.fail()) {
    ExitGracefully(("CEnthalpyModel::WriteOutputFileHeaders: Unable to open output file "+filename+" for writing.").c_str(),FILE_OPEN_ERR);
  }

  _STREAMOUT<<"time[d],date,hour,air temp.["+to_string(DEG_SYMBOL)+"C]"<<",";
  for(int p=0;p<_pModel->GetNumSubBasins();p++)
  {
    if(_pModel->GetSubBasin(p)->IsGauged()) {
      _STREAMOUT<<"Ein[MJ/d],Eout[MJ/d],";
      _STREAMOUT<<"Q_sens[MJ/d],Q_lat[MJ/d],Q_GW[MJ/d],Q_rad[MJ/d],Q_fric[MJ/d],channel storage[MJ],";
    }
  }
  _STREAMOUT<<endl;
}

//////////////////////////////////////////////////////////////////
/// \brief Writes minor transport output at the end of each timestep (or multiple thereof)
/// \param &Options [in] Global model options information
/// \param &tt [in] Local (model) time at the end of the pertinent time step
//
void CEnthalpyModel::WriteMinorOutput(const optStruct& Options,const time_struct& tt)
{
  //StreamTemperatures.csv and Temperatures.csv (similar to concentration output files)
  CConstituentModel::WriteMinorOutput(Options,tt);

  if(tt.model_time==0.0) { return; }

  int    p;
  double Q_sens,Q_lat,Q_GW,Q_rad,Q_fric;
  
  string thisdate=tt.date_string;
  string thishour=DecDaysToHours(tt.julian_day);

  _STREAMOUT<<tt.model_time <<","<<thisdate<<","<<thishour;
  _STREAMOUT<<","<<_pModel->GetAvgForcing("TEMP_AVE")<<",";

  for(p=0;p<_pModel->GetNumSubBasins();p++)
  {
    if(_pModel->GetSubBasin(p)->IsGauged()) {
      GetEnergyLossesFromReach(p,Q_sens,Q_lat,Q_GW,Q_rad,Q_fric);
      _STREAMOUT<<0.5*(_aMinHist[p][0]+_aMinHist[p][1])                                       <<",";
      _STREAMOUT<<0.5*(_aMout_last[p] +_aMout[p][_pModel->GetSubBasin(p)->GetNumSegments()-1])<<",";
      _STREAMOUT<<Q_sens<<","<<Q_lat<<","<<Q_GW<<","<<Q_rad<<","<<Q_fric<<",";
      _STREAMOUT<<_channel_storage[p]<<","; 
    }
  }
  _STREAMOUT<<endl;

}

//////////////////////////////////////////////////////////////////
/// \brief Writes minor transport output at the end of each timestep (or multiple thereof)
/// \param &Options [in] Global model options information
/// \param &tt [in] Local (model) time at the end of the pertinent time step
//
void CEnthalpyModel::CloseOutputFiles() 
{
  CConstituentModel::CloseOutputFiles();
  _STREAMOUT.close();
}

//////////////////////////////////////////////////////////////////
/// \brief for unit testing
//
void TestEnthalpyTempConvert()
{
  ofstream TEST;
  TEST.open("EnthalpyTest.csv");
  TEST<<"h,T,Fi,h_reverse"<<endl;
  for(double h=-500; h<0.0; h+=10) {
    double T=ConvertVolumetricEnthalpyToTemperature(h);
    double Fi=ConvertVolumetricEnthalpyToIceContent(h);
    double hr=ConvertTemperatureToVolumetricEnthalpy(T,Fi);
    double dTdh=TemperatureEnthalpyDerivative(h);
    TEST<<h<<","<<T<<","<<Fi<<","<<dTdh<<","<<hr<<endl;
  }
  for(double h=0; h<=150; h+=15) {
    double T=ConvertVolumetricEnthalpyToTemperature(h);
    double Fi=ConvertVolumetricEnthalpyToIceContent(h);
    double hr=ConvertTemperatureToVolumetricEnthalpy(T,Fi);
    double dTdh=TemperatureEnthalpyDerivative(h);
    TEST<<h<<","<<T<<","<<Fi<<","<<dTdh<<","<<hr<<endl;
  }
  TEST.close();

  ExitGracefully("TestEnthalpyTempConvert",SIMULATION_DONE);

}
