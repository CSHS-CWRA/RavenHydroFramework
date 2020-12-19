/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2021 the Raven Development Team
----------------------------------------------------------------
CEnthalpyModel routines related to energy/enthalpy transport
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "EnergyTransport.h"

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
/// \brief Returns ice content [0..1] of storage unit indexed by state variable index iWater
/// \param state_vars - vector of current state variables 
/// \param iWater - state variable index of water compartment of interest
//
double CEnthalpyModel::GetIceContent(const double *state_vars,const int iWater) const
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
/// \brief returns ice fraction of routed water in subbasin p at current point in time
/// \notes only used for reporting; calculations exclusively in terms of mass/energy
/// \param p [in] subbasin index 
//
double CEnthalpyModel::GetOutflowIceFraction(const int p) const
{
  CReservoir *pRes=_pModel->GetSubBasin(p)->GetReservoir();

  double ice_content;
  if(pRes==NULL)
  {
    double MJ_per_d=_aMout[p][_pModel->GetSubBasin(p)->GetNumSegments()-1];
    double flow=_pModel->GetSubBasin(p)->GetOutflowRate();
    //hv=MJ/d / M3/d = [MJ/m3]
    ice_content=ConvertVolumetricEnthalpyToIceContent(MJ_per_d/flow/SEC_PER_DAY); //[C]
    if(flow<=0) { ice_content= 0.0; }
  }
  else {
    ice_content=ConvertVolumetricEnthalpyToIceContent(_aMres[p]/pRes->GetStorage()); //[C]
  }
  return ice_content;
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
  return 9805*Q/perim*slope*WATT_TO_MJ_PER_D;
}

//////////////////////////////////////////////////////////////////
/// \brief Updates source terms for energy balance on subbasin reaches each time step
/// \details both _aEnthalpyBeta and _aEnthalpySource (and its history) are generated here
/// \param p    subbasin index
//
void   CEnthalpyModel::UpdateReachEnergySourceTerms(const int p)
{
  double tstep=_pModel->GetOptStruct()->timestep;

  const CSubBasin *pBasin=_pModel->GetSubBasin(p);

  if(pBasin->IsHeadwater()) { return; } //no reach, no need

  int                   k=pBasin->GetReachHRUIndex();
  ExitGracefullyIf(k==DOESNT_EXIST,"CEnthalpyModel::UpdateReachEnergySourceTerms: subbasin missing reach HRU for temperature simulation",BAD_DATA);
  const CHydroUnit  *pHRU=_pModel->GetHydroUnit(k);
  int                iAET=_pModel->GetStateVarIndex(AET);

  double SW      =pHRU->GetForcingFunctions()->SW_radia_net;       //[MJ/m2/d]
  double LW      =pHRU->GetForcingFunctions()->LW_radia_net;       //[MJ/m2/d]
  double AET     =pHRU->GetStateVarValue(iAET)/MM_PER_METER/tstep; //[m/d]
  double temp_air=pHRU->GetForcingFunctions()->temp_ave;           //[C]
  double temp_GW  =0.0; // JRC: should I get this as the average groundwater temperature in all subbasin HRUs? 

  double hstar    =pHRU->GetSurfaceProps()->convection_coeff; //[MJ/m2/d/K]  - TMP_DEBUG - THIS SHOULD BE A REACH PROPERTY, NOT SURFACE PROP.
  double qmix     =pBasin->GetHyporheicFlux(); //[m/d]
  double bed_ratio=pBasin->GetTopWidth()/pBasin->GetWettedPerimeter();
  double dbar     =pBasin->GetRiverDepth();

  double Qf       =GetReachFrictionHeat(pBasin->GetOutflowRate(),pBasin->GetBedslope(),pBasin->GetWettedPerimeter());//[MJ/m2/d]  

  double S(0.0);                         //source term [MJ/m3/d]
  S+=(SW+LW)/dbar;                       //net incoming energy term
  S-=AET*DENSITY_WATER*LH_VAPOR/dbar;    //latent heat flux term
  S+=Qf/dbar;                            //friction-generated heating term
  S+=hstar*temp_air/dbar;                //convection with atmosphere
  S+=qmix*HCP_WATER*bed_ratio*temp_GW/dbar; //hyporheic mixing

  //cout<<" S terms: "<<(SW+LW)<<"|"<<AET*DENSITY_WATER*LH_VAPOR<<"|"<<Qf<<"|"<<hstar*temp_air<<"|"<<qmix*HCP_WATER*DENSITY_WATER*bed_ratio*temp_GW<<endl;

  _aEnthalpyBeta[p]=(hstar/dbar + qmix/dbar*HCP_WATER*bed_ratio)/HCP_WATER;

  //cout<<_aEnthalpyBeta[p]<<" "<<hstar/HCP_WATER<<" "<< qmix*bed_ratio<<endl;

  int nMinHist=_pModel->GetSubBasin(p)->GetInflowHistorySize();
  for(int n=nMinHist-1;n>0;n--) {
    _aEnthalpySource[p][n]=_aEnthalpySource[p][n-1];
  }
  _aEnthalpySource[p][0]=S;

  //cout<<"EnthalpySource: p:"<<p<<" "<<S<<" "<<_aEnthalpyBeta[p]<<"/d"<<endl;
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
 
  //return 0.0;

  double tstep=_pModel->GetOptStruct()->timestep;

  const CSubBasin *pBasin=_pModel->GetSubBasin(p);
  int                   k=pBasin->GetReachHRUIndex();
  const CHydroUnit  *pHRU=_pModel->GetHydroUnit(k);
  int                iAET=_pModel->GetStateVarIndex(AET);

  double SW       =pHRU->GetForcingFunctions()->SW_radia_net;       //[MJ/m2/d]
  double LW       =pHRU->GetForcingFunctions()->LW_radia_net;       //[MJ/m2/d]
  double AET      =pHRU->GetStateVarValue(iAET)/MM_PER_METER/tstep; //[m/d]
  double temp_air =pHRU->GetForcingFunctions()->temp_ave;           //[C]
  double temp_GW  =0.0; // TMP_DEBUG JRC: should I get this as the average groundwater temperature in all subbasin HRUs? 

  double hstar    =pHRU->GetSurfaceProps()->convection_coeff; //[MJ/m2/d/K]  - TMP_DEBUG - THIS SHOULD BE A REACH PROPERTY, NOT SURFACE PROP.
  double qmix     =pBasin->GetHyporheicFlux(); //[m/d]
  double bed_ratio=pBasin->GetTopWidth()/pBasin->GetWettedPerimeter();
  double dbar     =pBasin->GetRiverDepth();

  double Qf       =GetReachFrictionHeat(pBasin->GetOutflowRate(),pBasin->GetBedslope(),pBasin->GetWettedPerimeter());//[MJ/m2/d]  

  int nMinHist              =pBasin->GetInflowHistorySize();
  const double * aRouteHydro=pBasin->GetRoutingHydrograph();
  const double * aQin       =pBasin->GetInflowHistory();

  double *zk=new double[nMinHist]; //percentage of flows in zone k
  double *Ik=new double[nMinHist];
  //double sum=0.0;

  //calculate mean residence time
  double tr_mean=0.0;
  for(int i=0;i<nMinHist;i++) {
    tr_mean+=aRouteHydro[i]*i*tstep;
  }

  //calculate flow linked to each zone k, zk[k]
  for(int k=1;k<nMinHist;k++) {
    for(int i=k;i<nMinHist;i++) {
      zk[k]=aQin[k]*aRouteHydro[i]*(tr_mean/i*tstep)*SEC_PER_DAY; // [m3/d]
    }
  }

  //calculate integral term Ik [degC*d]
  double hin;
  double beta=_aEnthalpyBeta[p];
  double betaterm=(1.0-exp(-beta*tstep));
  Ik[0]=0.0;
  for(int k=1;k<nMinHist;k++)
  {
    hin=(_aMinHist[p][k]/aQin[k]);
    Ik[k]=hin*betaterm*exp(-beta*(k-1)*tstep);
    for(int j=1;j<k;j++)
    {
      Ik[k]+=tstep/beta*_aEnthalpySource[p][j-1];
      if(j>2) {
        Ik[k]+=(betaterm*betaterm/beta/beta)*_aEnthalpySource[p][j-1]*exp(-beta*(j-2))*tstep;
      }
    }
    Ik[k]+=betaterm/beta/beta*_aEnthalpySource[p][0];
  }

  Q_sens = Q_lat = Q_GW = Q_rad = Q_fric =0.0; //incoming radiation [MJ/d]

  double kprime=qmix/dbar*HCP_WATER*DENSITY_WATER*bed_ratio;

  for(int k=1;k<nMinHist;k++)
  {
    Q_sens+=zk[k]*(hstar/dbar)*(temp_air*tstep-Ik[k]);//[m3/d]*[MJ/m2/d/K]*[1/m]*[K]*[d]=[MJ/d]

    Q_GW  +=zk[k]*kprime*(temp_GW*tstep-Ik[k]);

    Q_rad +=zk[k]*(SW+LW)/dbar*tstep;                //[m3/d]*[MJ/m2/d]*[1/m]*[d]=[MJ/d]

    Q_lat -=zk[k]*AET*DENSITY_WATER*LH_VAPOR/dbar*tstep;

    Q_fric+=zk[k]*Qf/dbar*tstep;
  }

  return -(Q_sens+Q_GW+Q_rad+Q_lat+Q_fric)*tstep; //total energy lost [MJ]

  delete[] Ik;
  delete[] zk;
}


void CEnthalpyModel::Initialize()
{

  CConstituentModel::Initialize();//call base class

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
    hv=ConvertTemperatureToVolumetricEnthalpy(CGlobalParams::GetParams()->init_stream_temp,0.0);
    for(int p=0;p<_pModel->GetNumSubBasins();p++) {
      for(int i=0; i<_pModel->GetSubBasin(p)->GetNumSegments(); i++) {
        _aMout[p][i]=_pModel->GetSubBasin(p)->GetOutflowRate()*hv; //not really - requires outflow rate from all segments. Don't have access to this.
      }
      _aMout_last[p]=_aMout[p][0];
    }
  }
}


void    CEnthalpyModel::ApplyConvolutionRouting(const int p,const double *aRouteHydro,const int nSegments,const int nMinHist,const double *aMinHist,const double &tstep,double *aMout_new) const
{
  //support energy exchange along reach for enthalpy transport
  double term1,term2=0.0;
  double beta=_aEnthalpyBeta[p];

  //beta=0.0; //TMP DEBUG
  aMout_new[nSegments-1]=0.0;
  //aMinHist and source term both have same indexing;
  for(int i=0;i<nMinHist;i++)
  {
    term1=aMinHist[i]*exp(-beta*i*tstep);
    term2=0.0;
    for(int j=0;j<i;j++) {
      term2+=_aEnthalpySource[p][i-j]/beta*exp(-beta*(i-j)*tstep)*(1.0-exp(-beta*tstep));
    }
    //term2=0;//TMP DEBUG
    aMout_new[nSegments-1]+=aRouteHydro[i]*(term1+term2);
  }
}
void   CEnthalpyModel::UpdateMassOutflows(const int p,
                                          double *aMoutnew,
                                          double  &ResMass,
                                          double  &MassOutflow,
                                          const   optStruct &Options,
                                          const   time_struct &tt,
                                          bool    initialize)
{
  CConstituentModel::UpdateMassOutflows(p,aMoutnew,ResMass,MassOutflow,Options,tt,initialize);
  UpdateReachEnergySourceTerms(p);
}
void TestEnthalpyTempConvert()
{
  ofstream OUT;
  OUT.open("EnthalpyTest.csv");
  OUT<<"h,T,Fi,h_reverse"<<endl;
  for(double h=-500; h<0.0; h+=10) {
    double T=ConvertVolumetricEnthalpyToTemperature(h);
    double Fi=ConvertVolumetricEnthalpyToIceContent(h);
    double hr=ConvertTemperatureToVolumetricEnthalpy(T,Fi);
    double dTdh=TemperatureEnthalpyDerivative(h);
    OUT<<h<<","<<T<<","<<Fi<<","<<dTdh<<","<<hr<<endl;
  }
  for(double h=0; h<=150; h+=15) {
    double T=ConvertVolumetricEnthalpyToTemperature(h);
    double Fi=ConvertVolumetricEnthalpyToIceContent(h);
    double hr=ConvertTemperatureToVolumetricEnthalpy(T,Fi);
    double dTdh=TemperatureEnthalpyDerivative(h);
    OUT<<h<<","<<T<<","<<Fi<<","<<dTdh<<","<<hr<<endl;
  }
  OUT.close();

  ExitGracefully("TestEnthalpyTempConvert",SIMULATION_DONE);

}
