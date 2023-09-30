/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2022 the Raven Development Team
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
  _aBedTemp=NULL;
  _anyGaugedLakes=false;
}

//////////////////////////////////////////////////////////////////
/// \brief enthalpy model destructor
//
CEnthalpyModel::~CEnthalpyModel()
{
  delete [] _aEnthalpyBeta;
  delete [] _aBedTemp;
  for(int p=0;p<_pModel->GetNumSubBasins();p++) { delete[] _aEnthalpySource; } delete [] _aEnthalpySource;
}

//////////////////////////////////////////////////////////////////
/// \brief Calculate temperature for reporting (converts to degrees C) (equivalent to 'concentration' for enthalpy)
/// \param M [in] energy in [MJ/m2]
/// \param V [in] volume in mm
//
double CEnthalpyModel::CalculateReportingConcentration(const double &M,const double &V) const
{
  if(fabs(V)>1e-6) {
    return ConvertVolumetricEnthalpyToTemperature(M/V*MM_PER_METER);  //[MJ/m3]->[C]
  }
  return 0.0;// JRC: should this default to zero? or NA?
}
//////////////////////////////////////////////////////////////////
/// \brief Converts basic energy units [deg C] to [mJ/mm/m2]
/// \param T [in] temperature in deg C
/// \returns concentration, in mJ/mm/m2
//
double CEnthalpyModel::ConvertConcentration(const double &T) const
{
    //specified in C, convert to MJ/m2
    double pctfroz=0.0;
    if (T<0.0){pctfroz=1.0;} //treats all 0 degree water as unfrozen
    return ConvertTemperatureToVolumetricEnthalpy(T,pctfroz)/MM_PER_METER; //[C]->[MJ/m3]*[M/mm]=[MJ/mm-m2]
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
  }
  else {
    return ConvertVolumetricEnthalpyToTemperature(_aMres[p]/pRes->GetStorage());
  }
}
//////////////////////////////////////////////////////////////////
/// \brief returns bed temperature from reach p
/// \param p [in] subbasin index
/// \returns outflow temperature [C]
//
double CEnthalpyModel::GetBedTemperature(const int p) const {
  return _aBedTemp[p];
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
/// \param T - dirichlet temperature (or special flag -9999)
//  \returns volumetric enthalpy, h, in MJ/mm-m2
//
double CEnthalpyModel::GetDirichletEnthalpy(const CHydroUnit* pHRU,const double& T) const
{
  if(T!=DIRICHLET_TEMP)
  {
    double pctFroz=0.0; //assumes liquid water for flows (reasonable)
    return ConvertTemperatureToVolumetricEnthalpy(T,pctFroz)/MM_PER_METER;    //[C]->[MJ/m3]->[MJ/mm-m2]
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
/// \param top_width [in] top_width [m]
/// from Theurer et al 1984, US Fish and Wildlife Serviec FWS/OBS-85/15, as reported in MacDonald, Boon, and Byrne 2014
//
double CEnthalpyModel::GetReachFrictionHeat(const double &Q,const double &slope,const double &top_width) const
{
  if (top_width<PRETTY_SMALL){return 0.0;}
  return GRAVITY*DENSITY_WATER*Q/top_width*slope*WATT_TO_MJ_PER_D;
}

//////////////////////////////////////////////////////////////////
/// \brief Set bed temperature
/// \param p     [in] subbasin index
/// \param val   [in] bed temperature (C)
//
void CEnthalpyModel::SetBedTemperature(const int p,const double &val)
{
  _aBedTemp[p]=val;
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
  int    mHypo           =2; //_pModel->GetHyporheicIndex();//TMP DEBUG - this might not be the compartment corresponding to "GW"
  int    iGW             =_pModel->GetStateVarIndex(SOIL,mHypo);

  double temp_air =pHRU->GetForcingFunctions()->temp_ave;           //[C]
  double SW       =pHRU->GetForcingFunctions()->SW_radia_net;       //[MJ/m2/d]
  double LW       =pHRU->GetForcingFunctions()->LW_radia_net;       //[MJ/m2/d]
  double LW_in    =pHRU->GetForcingFunctions()->LW_incoming;        //[MJ/m2/d]
  double AET      =pHRU->GetStateVarValue(iAET)/MM_PER_METER/tstep; //[m/d]

  double temp_bed =_aBedTemp[p];

  double temp_GW  =pBasin->GetAvgConcentration(iGW);
  double hstar    =pBasin->GetConvectionCoeff(); //[MJ/m2/d/K]
  double qmix     =pBasin->GetHyporheicFlux();   //[m/d]
  double bed_ratio=pBasin->GetTopWidth()/max(pBasin->GetWettedPerimeter(),0.001);
  double cel_corr =1.0;//pBasin->GetReferenceCelerity()/pBasin->GetCelerity();
  double dbar     =max(pBasin->GetRiverDepth(),0.001);
  double dbed     =pBasin->GetRiverbedThickness(); //[m]
  double kbed     =pBasin->GetRiverbedConductivity()/0.5/dbed; //[MJ/m2/d/K]

  double Qf       =GetReachFrictionHeat(pBasin->GetOutflowRate(),pBasin->GetBedslope(),pBasin->GetTopWidth());//[MJ/m2/d]

  double S(0.0);                            //source term [MJ/m3/d]
  S+=(SW+LW)/dbar;                          //net incoming energy term
  S-=AET*DENSITY_WATER*LH_VAPOR/dbar;       //latent heat flux term
  S+=Qf/dbar;                               //friction-generated heating term
  S+=hstar*temp_air/dbar;                   //convection with atmosphere
  S+=kbed* temp_bed/dbar;                   //conductive exchange with bed
  S+=qmix*HCP_WATER*bed_ratio*temp_GW/dbar; //hyporheic mixing

  _aEnthalpyBeta[p]=(hstar/dbar + qmix/dbar*HCP_WATER*bed_ratio + kbed/dbar)/HCP_WATER*cel_corr;

  int nMinHist=_pModel->GetSubBasin(p)->GetInflowHistorySize();
  for(int n=nMinHist-1;n>0;n--) {
    _aEnthalpySource[p][n]=_aEnthalpySource[p][n-1];
  }
  _aEnthalpySource[p][0]=S*cel_corr;
}


//////////////////////////////////////////////////////////////////
/// \brief Calculates individual energy gain terms from reach over current time step
/// \param p    subbasin index
/// \param Q_sens [out] energy gain from convective/sensible heating [MJ/d]
/// \param Q_sens [out] energy gain from conductive exchange with bed [MJ/d]
/// \param Q_lat  [out] energy gain from latent heat exchange [MJ/d]
/// \param Q_GW   [out] energy gain from groundwater mixing [MJ/d]
/// \param Q_rad  [out] energy gain from radiation of surface [MJ/d]
/// \param Q_fric [out] energy gain from friction [MJ/d]
/// \returns total energy lost from reach over current time step [MJ]
//
double CEnthalpyModel::GetEnergyLossesFromReach(const int p,double &Q_sens,double &Q_cond,double &Q_lat,double &Q_GW,double &Q_rad,double &Q_fric) const
{

  double tstep=_pModel->GetOptStruct()->timestep;

  const CSubBasin *pBasin=_pModel->GetSubBasin(p);

  if(pBasin->IsHeadwater())
  {
    Q_sens=Q_lat=Q_GW=Q_rad=Q_fric=0.0;
    return 0.0;
  }

  const CHydroUnit* pHRU=_pModel->GetHydroUnit(pBasin->GetReachHRUIndex());
  int               iAET=_pModel->GetStateVarIndex(AET);
  int              mHypo=2; //_pModel->GetHyporheicIndex();//TMP DEBUG - this might not be the compartment corresponding to "GW"
  int                iGW=_pModel->GetStateVarIndex(SOIL,mHypo);

  double temp_air =pHRU->GetForcingFunctions()->temp_ave;           //[C]
  double SW       =pHRU->GetForcingFunctions()->SW_radia_net;       //[MJ/m2/d]
  double LW       =pHRU->GetForcingFunctions()->LW_radia_net;       //[MJ/m2/d]
  //double LW       =STEFAN_BOLTZ*EMISS*pow(T-Ta,4); T=GetOutflowConcentration(p);
  double LW_in    =pHRU->GetForcingFunctions()->LW_incoming;        //[MJ/m2/d]
  double AET      =pHRU->GetStateVarValue(iAET)/MM_PER_METER/tstep; //[m/d]

  double temp_bed =_aBedTemp[p];

  double temp_GW  =pBasin->GetAvgConcentration(iGW);
  double hstar    =pBasin->GetConvectionCoeff(); //[MJ/m2/d/K]
  double qmix     =pBasin->GetHyporheicFlux();   //[m/d]
  double bed_ratio=pBasin->GetTopWidth()/max(pBasin->GetWettedPerimeter(),0.001);
  double cel_corr =1.0;//pBasin->GetReferenceCelerity()/pBasin->GetCelerity();
  double dbar     =pBasin->GetRiverDepth();
  double dbed     =pBasin->GetRiverbedThickness(); //[m]
  double kbed     =pBasin->GetRiverbedConductivity()/0.5/dbed; //[MJ/m2/d/K]

  double Qf       =GetReachFrictionHeat(pBasin->GetOutflowRate(),pBasin->GetBedslope(),pBasin->GetTopWidth());//[MJ/m2/d]

  dbar=max(dbar,0.001);

  int nMinHist              =pBasin->GetInflowHistorySize();
  const double * aRouteHydro=pBasin->GetRoutingHydrograph();
  const double * aQin       =pBasin->GetInflowHistory();

  double *zk=new double[nMinHist]; //percentage of flows in zone k
  double *Ik=new double[nMinHist]; //Integral term [degC-d]

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
    zsum+=zk[k]*cel_corr;
  }

  //calculate integral term I_m^n [degC*d] (see paper)
  double hin,beta,gamma;
  beta =_aEnthalpyBeta[p];//[1/d]
  beta =max(beta,1e-9); //to avoid divide by zero error
  gamma=(1.0-exp(-beta*tstep));
  Ik[0]=0.0;
  for(int m=1;m<nMinHist;m++)
  {
    hin=(_aMinHist[p][m-1]/(aQin[m-1]*SEC_PER_DAY)); //[MJ/d]/[m3/d]->[MJ/m3]
    if(_aMinHist[p][m-1]<PRETTY_SMALL) { hin=0.0; }

    Ik[m]=hin/beta*gamma*exp(-beta*(m-1)*tstep); //[MJ/m3]*[d]
    for(int j=1;j<m;j++)
    {
      Ik[m]+=gamma*gamma/beta/beta*_aEnthalpySource[p][j]*exp(-beta*(m-j-1)*tstep);
    }
    Ik[m]+=(tstep/beta-gamma/beta/beta)*_aEnthalpySource[p][0];

    Ik[m]*=1.0/HCP_WATER/tstep/cel_corr;// [MJ-d/m3] /[MJ/m3/K] /[d] = [degC];
    //double tmp=(tstep/beta-gamma/beta/beta)* _aEnthalpySource[p][0]/HCP_WATER;
    //cout<<"Tk["<<m<<"]="<<Ik[m]<<" Tin="<<hin/HCP_WATER<<" src term="<<tmp<<endl;
  }

  Q_sens = Q_cond = Q_lat = Q_GW = Q_rad = Q_fric =0.0; //incoming radiation [MJ/d]

  double kprime=qmix/dbar*HCP_WATER*bed_ratio;

  for(int m=1;m<nMinHist;m++)
  {
    Q_sens+=zk[m]*(hstar/dbar)*(temp_air-Ik[m])*tstep;   //[m3/d]*[MJ/m2/d/K]*[1/m]*[K]*[d]=[MJ/d]

    Q_cond+=zk[m]*(kbed /dbar)*(temp_bed-Ik[m])*tstep;   //[m3/d]*[MJ/m2/d/K]*[1/m]*[K]*[d]=[MJ/d]

    Q_GW  +=zk[m]*kprime      *(temp_GW-Ik[m])*tstep;

    Q_rad +=zk[m]*(SW+LW)/dbar*tstep;                    //[m3/d]*[MJ/m2/d]*[1/m]*[d]=[MJ/d]

    Q_lat -=zk[m]*AET*DENSITY_WATER*LH_VAPOR/dbar*tstep; //[m3/d]*[m/d]*[kg/m3]*[MJ/kg]*[1/m]*[d] = [MJ/d]

    Q_fric+=zk[m]*Qf/dbar*tstep;

   // cout<<"Ik term: "<<temp_air<<" "<<Ik[m]/tstep<<" beta: "<<beta<<" "<<gamma<<endl;
   // cout<<" Losses: "<<k<<": "<<Q_sens<<" "<<Q_GW<<" "<<Q_rad<<" "<<Q_lat<<" "<<Q_fric<< " "<<zk[k]*kprime*(temp_GW*tstep-Ik[k])<<" "<<Ik[k]<<" "<<zk[k]<<endl;
   // cout<<"mean res time: "<<tr_mean<<endl;
   // cout<<" "<<hstar<<" "<< dbar<<" "<<temp_air*tstep<<" "<<zk[k]<<endl;
   // cout<<"z["<<m<<"]: "<<zk[m]/SEC_PER_DAY<<"("<<zk[m]<<") "<<aQin[m]<<" "<<aRouteHydro[m]<<endl;
  }

  delete[] Ik;
  delete[] zk;

  return -(Q_sens+Q_cond+Q_GW+Q_rad+tstep*Q_lat+Q_fric)*tstep; //total energy lost [MJ] //TMP DEBUG - WHY IS THE tstep needed HERE before Q_lat!!?? Is it the water lost from the system??
}
//////////////////////////////////////////////////////////////////
/// \brief Calculates individual energy gain terms from lake/reservoir over current time step
/// \param p    subbasin index
/// \param Q_sens [out] energy gain from convective/sensible heating [MJ/d]
/// \param Q_cond [out] energy gain from conductive exchange with bed [MJ/d]
/// \param Q_lat  [out] energy gain from latent heat exchange [MJ/d]
/// \param Q_rad  [out] energy gain from radiation of surface [MJ/d]
/// \param Q_rain [out] energy gain from precip inputs [MJ/d]
/// \returns total energy lost from reach over current time step [MJ]
//
double CEnthalpyModel::GetEnergyLossesFromLake(const int p, double &Q_sens, double &Q_cond, double &Q_lat, double &Q_rad, double &Q_rain) const
{
  double tstep = _pModel->GetOptStruct()->timestep;

  CSubBasin  *pBasin = _pModel->GetSubBasin(p);
  CReservoir *pRes   = _pModel->GetSubBasin(p)->GetReservoir();
  if (pRes==NULL){return 0.0;}

  double V_new = pRes->GetStorage();
  double V_old = pRes->GetOldStorage();
  double Q_new = pRes->GetOutflowRate() * SEC_PER_DAY;
  double Q_old = pRes->GetOldOutflowRate() * SEC_PER_DAY;
  double A_new = pRes->GetSurfaceArea();
  double A_old = pRes->GetOldSurfaceArea();
  double A_avg = 0.5 * (A_new + A_old);

  CHydroUnit*   pHRU=_pModel->GetHydroUnit(pRes->GetHRUIndex());

  double Acorr=pHRU->GetArea()*M2_PER_KM2/A_avg; //handles the fact that GetAET() returns mm/d normalized by HRU area, not actual area

  double SW(0), LW(0), LW_in(0), temp_air(0), AET(0);
  double hstar(0),ksed(0),Vsed=0.001;
  if(pHRU!=NULL) { //otherwise only simulate advective mixing+ rain input
    temp_air =pHRU->GetForcingFunctions()->temp_ave;           //[C]
    SW       =pHRU->GetForcingFunctions()->SW_radia_net;       //[MJ/m2/d]
    LW       =pHRU->GetForcingFunctions()->LW_radia_net;       //[MJ/m2/d]
    LW_in    =pHRU->GetForcingFunctions()->LW_incoming;        //[MJ/m2/d]

    AET      =pRes->GetAET()*Acorr/ MM_PER_METER ;                   //[m/d] //*pHRU->GetArea()/A_avg

    hstar = pRes->GetLakeConvectionCoeff();                    //[MJ/m2/d/K]
    Vsed  = pRes->GetLakebedThickness() * pHRU->GetArea()*M2_PER_KM2;
    ksed  = pRes->GetLakebedConductivity() / 0.5 / pRes->GetLakebedThickness();//[MJ/m2/d/K]
  }
  double T_new =ConvertVolumetricEnthalpyToTemperature(_aMres[p]      / V_new);
  double T_old =ConvertVolumetricEnthalpyToTemperature(_aMres_last[p] / V_old);
  double Ts_new=ConvertVolumetricEnthalpyToTemperature(_aMsed[p]      / Vsed );
  double Ts_old=ConvertVolumetricEnthalpyToTemperature(_aMsed_last[p] / Vsed );

  Q_sens=hstar* A_avg * (temp_air -0.5*(T_new+T_old));
  Q_cond=ksed * A_avg * (0.5*(Ts_new+Ts_old)- 0.5 * (T_new + T_old));
  Q_rad =(SW+LW)*A_avg;
  Q_lat =-(AET * DENSITY_WATER * LH_VAPOR)*A_avg;
  Q_rain=_aMresRain[p];

  return -(Q_sens + Q_cond  + Q_rad +  Q_lat + Q_rain) * tstep; //[MJ]
}
//////////////////////////////////////////////////////////////////
/// \brief Calculates Res_mass (total reservoir energ) and ResSedMass (total sediment energy) in basin with reservoir at end of time step
///
/// \param p           [in]  subbasin index
/// \param aMout_new[] [in]  Array of energy outflows at downstream end of each segment at end of current timestep [MJ/d] [size: nsegs]
/// \param Res_mass    [out] calculated reservoir mass at end of time step
/// \param ResSedMass  [out] calculated reservoir mass in sediment at end of time step
/// \param &Options    [in]  Global model options information
/// \param tt          [in]  Time structure
//
void   CEnthalpyModel::RouteMassInReservoir   (const int          p,          // SB index
                                               const double      *aMout_new,  // [MJ/d][size: nsegs ]
                                                     double      &Res_mass,   // [MJ]
                                                     double      &ResSedMass, // [MJ]
                                               const optStruct   &Options,
                                               const time_struct &tt) const
{
  CSubBasin *pBasin= _pModel->GetSubBasin(p);
  CReservoir* pRes = _pModel->GetSubBasin(p)->GetReservoir();
  int nSegments    = _pModel->GetSubBasin(p)->GetNumSegments();

  double V_new=pRes->GetStorage       ();
  double V_old=pRes->GetOldStorage    ();
  double Q_new=pRes->GetOutflowRate   ()*SEC_PER_DAY;
  double Q_old=pRes->GetOldOutflowRate()*SEC_PER_DAY;
  double A_new=pRes->GetSurfaceArea();
  double A_old=pRes->GetOldSurfaceArea();
  double A_avg=0.5*(A_new+A_old);
  double tstep=Options.timestep;

  if((V_old<=0.0) || (V_new<=0.0)) { Res_mass=ResSedMass=0.0; return;} //handles dried out reservoir/lake

  CHydroUnit*   pHRU=_pModel->GetHydroUnit(pRes->GetHRUIndex());

  double Acorr=pHRU->GetArea()*M2_PER_KM2/A_avg; //handles the fact that GetAET() returns mm/d normalized by HRU area, not actual area

  double SW(0), LW(0), LW_in(0), temp_air(0), AET(0);
  double hstar(0), ksed(0), Vsed=0.001;
  if(pHRU!=NULL) { //otherwise only simulate advective mixing+ rain input
    temp_air =pHRU->GetForcingFunctions()->temp_ave;           //[C]
    SW       =pHRU->GetForcingFunctions()->SW_radia_net;       //[MJ/m2/d]
    LW       =pHRU->GetForcingFunctions()->LW_radia_net;       //[MJ/m2/d]
    LW_in    =pHRU->GetForcingFunctions()->LW_incoming;        //[MJ/m2/d]

    AET      =pRes->GetAET()*Acorr/ MM_PER_METER ;             //[m/d]

    Vsed     =pRes->GetLakebedThickness() * pHRU->GetArea()*M2_PER_KM2;
    hstar    =pRes->GetLakeConvectionCoeff(); //[MJ/m2/d/K]
    ksed     =pRes->GetLakebedConductivity()/0.5/ pRes->GetLakebedThickness();//[MJ/m2/d/K]
  }

  // N-R solution of Crank-nicolson problem as set of two non-linear algebraic equations
  // -----------------------------------------------------------------------------------------
  //dE   /dt=Qh_in-Qh_out+(Rnet*A+P*hrain*A)-ET*rho*LH*A+k*(Tair-T(E))*A+ksed*(Tsed-T(E))*A     //[MJ/d]
  //dEsed/dt=                                                           -ksed*(Tsed-T(E))*A     //[MJ/d]
  //
  //Allocate memory
  double  *B=new double  [2];
  double  *R=new double  [2];
  double **A=new double *[2];
  double **J=new double *[2];
  for (int i=0;i<2;i++){
    A[i]=new double [2];
    J[i]=new double [2];
    B[i]=R[i]=0.0;
    for (int j = 0; j < 2; j++) {A[i][j]=J[i][j]=0.0;}
  }
  double T_old    =ConvertVolumetricEnthalpyToTemperature(_aMres_last[p]/V_old);
  double T_sed_old=ConvertVolumetricEnthalpyToTemperature(_aMsed_last[p]/Vsed);

  B[0] = 0.5*(aMout_new[nSegments-1]+_aMout[p][nSegments-1])*tstep; // inflow [MJ]
  B[0]+=_aMresRain[p]*tstep;                                       // rainfall inputs [MJ]
  B[0]+= A_avg *(SW+LW)*tstep;                                     // net radiation [MJ]
  B[0]-= A_avg *(AET*DENSITY_WATER*LH_VAPOR)*tstep;                // latent heat [MJ]
  B[0]+= A_avg *(hstar*temp_air)*tstep;                            // sensible heat exhange [MJ]

  B[0]+=(1.0-0.5*tstep/V_old*Q_old)*_aMres_last[p];
  B[0]-=0.5*tstep*hstar*A_old*T_old;
  B[0]-=0.5*tstep*ksed *A_avg*T_old;
  B[0]+=0.5*tstep*ksed *A_avg*T_sed_old;

  B[1] =0.5*tstep*ksed *A_avg*T_old;
  B[1]+=_aMsed_last[p]-0.5*tstep*ksed *A_avg*T_sed_old;

  int iter=0;
  double change=ALMOST_INF;
  double tolerance=1e-2; //[deg C]
  double T_guess,Ts_guess;
  double dE,dEs;
  double E_guess  =_aMres_last[p];
  double Es_guess =_aMsed_last[p];

  while (change > tolerance)
  {
    T_guess = ConvertVolumetricEnthalpyToTemperature(E_guess/V_new);
    Ts_guess= ConvertVolumetricEnthalpyToTemperature(Es_guess/Vsed);

    A[0][0] = 1.0 + 0.5 * tstep  * Q_new/V_new;
    A[1][1] = 1.0;
    A[1][0] = 0.0;
    A[0][1] = 0.0;
    if (E_guess!=0.0){
      A[0][0]+= 0.5*tstep*hstar*A_new*T_guess/E_guess;
      A[0][0]+= 0.5*tstep*ksed *A_avg*T_guess/E_guess;  //[MJ/m2/d/K]*[m2]/[MJ/m3/K] = [m3/d]
      A[1][0]-= 0.5*tstep*ksed *A_avg*T_guess/E_guess;
    }
    if (Es_guess!=0.0){
      A[0][1]-= 0.5*tstep*ksed *A_avg*Ts_guess/Es_guess;
      A[1][1]+= 0.5*tstep*ksed *A_avg*Ts_guess/Es_guess;
    }

    //J_ij=A_ij+E*dA_ij/dE = Jacobian
    J[0][0]=A[0][0]+0.5*tstep*(hstar*A_new+ksed*A_avg) * (TemperatureEnthalpyDerivative(E_guess /V_new)/ V_new);
    J[0][1]=A[0][1]-0.5*tstep*(            ksed*A_avg) * (TemperatureEnthalpyDerivative(Es_guess/Vsed )/ Vsed );
    J[1][0]=A[1][0]-0.5*tstep*(            ksed*A_avg) * (TemperatureEnthalpyDerivative(E_guess /V_new)/ V_new);
    J[1][1]=A[1][1]+0.5*tstep*(            ksed*A_avg) * (TemperatureEnthalpyDerivative(Es_guess/Vsed )/ Vsed );

    R[0]   =-A[0][0]*E_guess-A[0][1]*Es_guess+B[0];
    R[1]   =-A[1][0]*E_guess-A[1][1]*Es_guess+B[1];

    //invert 2x2 matrix analytically
    double den= (J[0][1]*J[1][0]-J[0][0]*J[1][1]);
    dE  = (J[0][1]*R[1]-J[1][1]*R[0])/den;
    dEs = (J[1][0]*R[0]-J[0][0]*R[1])/den;

    change =sqrt(dE*dE+dEs*dEs)/ V_new / HCP_WATER; //convert to approx temp difference (for >0C water)

    E_guess  +=dE;
    Es_guess +=dEs;
    iter++;
    //if (iter>2){cout<<"iter: "<<iter<<" change: "<<change<<" tol: "<<tolerance<<endl;}
  }
  Res_mass=E_guess;
  ResSedMass=Es_guess;

  delete [] B;
  delete [] R;
  for (int i=0;i<2;i++){delete [] A[i]; delete[] J[i]; } delete[] A; delete[] J;
}
//////////////////////////////////////////////////////////////////
/// \brief returns total energy lost from subbasin reach over current time step [MJ]
/// \param p    subbasin index
/// \returns total energy lost from subbasin reach over current time step [MJ]
//
double CEnthalpyModel::GetNetReachLosses(const int p) const
{
  double Q_sens,Q_cond,Q_lat,Q_GW,Q_rad,Q_fric;
  return GetEnergyLossesFromReach(p,Q_sens,Q_cond,Q_lat,Q_GW,Q_rad,Q_fric);
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
  _aEnthalpySource=new double*[nSB];
  _aEnthalpyBeta  =new double [nSB];
  _aBedTemp       =new double [nSB];
  for(int p=0;p<nSB;p++) {
    int nMinHist =_pModel->GetSubBasin(p)->GetInflowHistorySize();
    _aEnthalpySource[p]=new double[nMinHist];
    for(int i=0; i<nMinHist;i++) { _aEnthalpySource[p][i]=0.0; }
    _aEnthalpyBeta  [p]=0.0;
    _aBedTemp       [p]=10.0;
  }

  // initialize stream temperatures if init_stream_temp is given
  //--------------------------------------------------------------------
  double hv;
  if(_pModel->GetGlobalParams()->GetParams()->init_stream_temp>0.0)
  {
    const CSubBasin *pBasin;

    hv=ConvertTemperatureToVolumetricEnthalpy(_pModel->GetGlobalParams()->GetParams()->init_stream_temp,0.0);
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

      if (pBasin->GetReservoir()!=NULL)
      {
        _aMres     [p]=hv*pBasin->GetReservoir()->GetStorage();
        _aMres_last[p]=_aMres[p];

        _aMsed     [p]=hv*_pModel->GetHydroUnit(pBasin->GetReservoir()->GetHRUIndex())->GetArea()*M2_PER_KM2*pBasin->GetReservoir()->GetLakebedThickness();
        _aMsed_last[p]=_aMsed[p];
      }

      //SOME ISSUE HERE WITH THIS INTIALIZATION - WHY?
      //NO PROBLEMS WITH ABOVE MASS FLOW RATE INITIALIZATION, ONLY EXTRA MASS BALANCE TERMS
      _channel_storage[p]=pBasin->GetChannelStorage()*hv/1.717;
      _initial_mass+=_channel_storage[p];
//      cout<<" init storage!"<<_initial_mass;

      _aBedTemp[p] = _pModel->GetGlobalParams()->GetParams()->init_stream_temp;
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
                                          double  &ResSedMass,
                                          double  &MassOutflow,
                                          const   optStruct &Options,
                                          const   time_struct &tt,
                                          bool    initialize)
{
  CConstituentModel::UpdateMassOutflows(p,aMoutnew,ResMass,ResSedMass,MassOutflow,Options,tt,initialize);//call base class
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

  string name;
  _STREAMOUT<<"time[d],date,hour,air temp.["+to_string(DEG_SYMBOL)+"C]"<<", net rad [MJ/m2/d],";
  for(int p=0;p<_pModel->GetNumSubBasins();p++)
  {
    CSubBasin *pSB=_pModel->GetSubBasin(p);
    if(pSB->GetName()=="") { name=to_string(pSB->GetID()); }
    else                   { name=pSB->GetName(); }

    if((pSB->IsGauged()) && (pSB->IsEnabled()))
    {
      _STREAMOUT<<name<<" Ein[MJ/m2/d],";
      _STREAMOUT<<name<<" Eout[MJ/m2/d],";
      _STREAMOUT<<name<<" Q_sens[MJ/m2/d],";
      _STREAMOUT<<name<<" Q_cond[MJ/m2/d],";
      _STREAMOUT<<name<<" Q_lat[MJ/m2/d],";
      _STREAMOUT<<name<<" Q_GW[MJ/m2/d],";
      _STREAMOUT<<name<<" Q_rad[MJ/m2/d],";
      _STREAMOUT<<name<<" Q_fric[MJ/m2/d],";
      _STREAMOUT<<name<<" channel storage[MJ/m2],";
    }
  }
  _STREAMOUT<<endl;

  //LakeEnergyBalances.csv
  //--------------------------------------------------------------------
  _anyGaugedLakes=false;
  for (int p = 0; p < _pModel->GetNumSubBasins(); p++)
  {
    CSubBasin* pSB = _pModel->GetSubBasin(p);
    if ((pSB->IsGauged()) && (pSB->IsEnabled()) &&
        (pSB->GetReservoir()!=NULL) && (pSB->GetReservoir()->GetHRUIndex() != DOESNT_EXIST)) { _anyGaugedLakes = true; }
  }

  if (_anyGaugedLakes)
  {
    string filename = "LakeEnergyBalances.csv";
    filename = FilenamePrepare(filename, Options);

    _LAKEOUT.open(filename.c_str());
    if (_LAKEOUT.fail()) {
      ExitGracefully(("CEnthalpyModel::WriteOutputFileHeaders: Unable to open output file " + filename + " for writing.").c_str(), FILE_OPEN_ERR);
    }

    string name;
    _LAKEOUT << "time[d],date,hour,air temp.[" + to_string(DEG_SYMBOL) + "C]" << ", net rad [MJ/m2/d],";
    for (int p = 0; p < _pModel->GetNumSubBasins(); p++)
    {
      CSubBasin *pSB = _pModel->GetSubBasin(p);
      if (pSB->GetName() == "") { name = to_string(pSB->GetID()); }
      else                      { name = pSB->GetName(); }

      if ((pSB->IsGauged()) && (pSB->IsEnabled()))
      {
        CReservoir* pRes = _pModel->GetSubBasin(p)->GetReservoir();
        if ((pRes != NULL) && (pRes->GetHRUIndex() != DOESNT_EXIST))
        {
          _LAKEOUT << name << " Ein [MJ/m2/d],";
          _LAKEOUT << name << " Eout [MJ/m2/d],";
          _LAKEOUT << name << " Q_rain [MJ/m2/d],";
          _LAKEOUT << name << " Q_sens [MJ/m2/d],";
          _LAKEOUT << name << " Q_cond [MJ/m2/d],";
          _LAKEOUT << name << " Q_lat [MJ/m2/d],";
          _LAKEOUT << name << " Q_rad [MJ/m2/d],";
          _LAKEOUT << name << " lake storage [MJ/m2],";
          _LAKEOUT << name << " lake temp [C],";
          _LAKEOUT << name << " lake sed temp [C],";
          _LAKEOUT << name << " inflow temp [C],";
          _LAKEOUT << name << " pct froz [0..1],";
        }
      }
    }
    _LAKEOUT << endl;
  }
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

  //StreamReachEnergyBalances.csv
  //--------------------------------------------------------------------
  int    p;
  double Q_sens,Q_cond,Q_lat,Q_GW,Q_rad,Q_fric,Q_rain;
  double Ein,Eout;

  string thisdate=tt.date_string;
  string thishour=DecDaysToHours(tt.julian_day);

  _STREAMOUT<<tt.model_time <<","<<thisdate<<","<<thishour<<",";
  _STREAMOUT << _pModel->GetAvgForcing(F_TEMP_AVE) << ",";
  _STREAMOUT << _pModel->GetAvgForcing(F_SW_RADIA_NET) + _pModel->GetAvgForcing(F_LW_RADIA_NET) << ",";

  double mult=1.0; //convert everything to MJ/m2/d
  for(p=0;p<_pModel->GetNumSubBasins();p++)
  {
    mult = 1.0 / _pModel->GetSubBasin(p)->GetReachLength() / _pModel->GetSubBasin(p)->GetTopWidth();
    if ((_pModel->GetSubBasin(p)->IsGauged()) && (_pModel->GetSubBasin(p)->IsEnabled()))
    {
      Ein = 0.5 * (_aMinHist[p][0] + _aMinHist[p][1]);
      Eout = 0.5 * (_aMout_last[p] + _aMout[p][_pModel->GetSubBasin(p)->GetNumSegments() - 1]);
      GetEnergyLossesFromReach(p, Q_sens, Q_cond, Q_lat, Q_GW, Q_rad, Q_fric);

      _STREAMOUT << mult * Ein << "," << mult * Eout << ",";
      _STREAMOUT << mult * Q_sens << "," << mult * Q_cond << "," << mult * Q_lat << ",";
      _STREAMOUT << mult * Q_GW << "," << mult * Q_rad << ",";
      _STREAMOUT << mult * Q_fric << "," << mult * _channel_storage[p] << ",";
    }
  }
  _STREAMOUT<<endl;

  //LakeEnergyBalances.csv
  //--------------------------------------------------------------------
  if (_anyGaugedLakes)
  {
    _LAKEOUT << tt.model_time << "," << thisdate << "," << thishour << ",";
    _LAKEOUT << _pModel->GetAvgForcing(F_TEMP_AVE) << ",";
    _LAKEOUT << _pModel->GetAvgForcing(F_SW_RADIA_NET) + _pModel->GetAvgForcing(F_LW_RADIA_NET) << ",";


    for (p = 0; p < _pModel->GetNumSubBasins(); p++)
    {
      if ((_pModel->GetSubBasin(p)->IsGauged()) && (_pModel->GetSubBasin(p)->IsEnabled()))
      {
        CReservoir* pRes = _pModel->GetSubBasin(p)->GetReservoir();

        if ((pRes != NULL) && (pRes->GetHRUIndex() != DOESNT_EXIST))
        {
          double HRUarea= _pModel->GetHydroUnit(pRes->GetHRUIndex())->GetArea()*M2_PER_KM2;
          double V_sed  = pRes->GetLakebedThickness() * HRUarea;
          double V_new  = pRes->GetStorage();

          mult = 1.0 / HRUarea;

          Ein  = 0.5 * (_aMout_last[p] + _aMout[p][_pModel->GetSubBasin(p)->GetNumSegments() - 1]);
          Eout = 0.5 * (_aMout_res [p] + _aMout_res_last[p]);

          double Qin=_pModel->GetSubBasin(p)->GetOutflowRate()*SEC_PER_DAY; //[m3/d]

          GetEnergyLossesFromLake(p,Q_sens,Q_cond,Q_lat,Q_rad,Q_rain);


          double lakeTemp  =ConvertVolumetricEnthalpyToTemperature(_aMres[p] / V_new);
          double sedTemp   =ConvertVolumetricEnthalpyToTemperature(_aMsed[p] / V_sed);
          double inflowTemp=ConvertVolumetricEnthalpyToTemperature(Ein / Qin);
          double pctFroz   =ConvertVolumetricEnthalpyToIceContent (_aMres[p] / V_new);

          if (V_new==0){lakeTemp  =-60; pctFroz=-1;}
          if (V_sed==0){sedTemp   =-60;}
          if (Qin  ==0){inflowTemp=-60;}

          _LAKEOUT << mult * Ein    << "," << mult * Eout      << ",";
          _LAKEOUT << mult * Q_rain << "," << mult * Q_sens    << ",";
          _LAKEOUT << mult * Q_cond << "," << mult * Q_lat     << ",";
          _LAKEOUT << mult * Q_rad  << "," << mult * _aMres[p] << ",";
          _LAKEOUT << lakeTemp      << "," << sedTemp          << ",";
          _LAKEOUT << inflowTemp    << "," << pctFroz          << ",";
        }
      }
    }
    _LAKEOUT << endl;
  }
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
    double T =ConvertVolumetricEnthalpyToTemperature(h);
    double Fi=ConvertVolumetricEnthalpyToIceContent(h);
    double hr=ConvertTemperatureToVolumetricEnthalpy(T,Fi);
    double dTdh=TemperatureEnthalpyDerivative(h);
    TEST<<h<<","<<T<<","<<Fi<<","<<dTdh<<","<<hr<<endl;
  }
  for(double h=0; h<=150; h+=15) {
    double T =ConvertVolumetricEnthalpyToTemperature(h);
    double Fi=ConvertVolumetricEnthalpyToIceContent(h);
    double hr=ConvertTemperatureToVolumetricEnthalpy(T,Fi);
    double dTdh=TemperatureEnthalpyDerivative(h);
    TEST<<h<<","<<T<<","<<Fi<<","<<dTdh<<","<<hr<<endl;
  }
  TEST.close();

  ExitGracefully("TestEnthalpyTempConvert",SIMULATION_DONE);

}
