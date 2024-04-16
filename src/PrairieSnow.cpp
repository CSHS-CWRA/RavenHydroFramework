/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2018 the Raven Development Team
------------------------------------------------------------------
  Prairie Blowing Snow Model
  Translated from Matt McDonald's MESH version
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "HydroUnits.h"
#include "SubBasin.h"
#include "Model.h"
#include "PrairieSnow.h"

//////////////////////////////////////////////////////////////////
/// \brief calculates total volumetric sensible heat energy of snowpack in [MJ/m2]
/// E = c_p*(rho/A)*T
/// \todo move to snow params
//
double CalculateSnowEnergyContent(const double &SWE,            //[mm]
                                  const double &snow_depth,     //[mm]
                                  const double &snow_liq,       //[mm]
                                  const double &snow_temp)      //[mm]
{
  double c_p=HCP_ICE*(SWE/snow_depth)+HCP_WATER*(snow_liq/snow_depth);  //[MJ/m3/K]

  return c_p*(snow_temp+FREEZING_TEMP)*(snow_depth/MM_PER_METER); //[MJ/m3/K]*[K]*[m]=MJ/m2
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the prairie blowing snow constructor
//
CmvPrairieBlowingSnow::CmvPrairieBlowingSnow(pbsm_type sub_type, CModel *pModel)
  :CLateralExchangeProcessABC(BLOWING_SNOW, pModel)
{
  _type=sub_type;
  _nDriftConnections=0;

  int iSnowAge      =pModel->GetStateVarIndex(SNOW_AGE);
  int iSWE          =pModel->GetStateVarIndex(SNOW);
  int iSnowLiq      =pModel->GetStateVarIndex(SNOW_LIQ);
  int iSnowDepth    =pModel->GetStateVarIndex(SNOW_DEPTH);
  int iSnowDriftTemp=pModel->GetStateVarIndex(SNODRIFT_TEMP);
  int iAtmosphere   =pModel->GetStateVarIndex(ATMOSPHERE);
  int iSnowDrift    =pModel->GetStateVarIndex(SNOW_DRIFT);

  CHydroProcessABC::DynamicSpecifyConnections(6);//nConnections=6

  iFrom[0]=iSnowAge;       iTo[0]=iSnowAge;      //rates[0]: SNOW_AGE->SNOW_AGE
  iFrom[1]=iSnowDepth;     iTo[1]=iSnowDepth;    //rates[1]: SNOW_DEPTH->SNOW_DEPTH (proxy for density)
  iFrom[2]=iSnowDriftTemp; iTo[2]=iSnowDriftTemp;//rates[2]: SNODRIFT_TEMP->SNODRIFT_TEMP
  iFrom[3]=iSnowLiq;       iTo[3]=iAtmosphere;   //rates[3]: SNOW_LIQ ->ATMOSPHERE (sublimation)
  iFrom[4]=iSWE;           iTo[4]=iAtmosphere;   //rates[4]: SNOW->ATMOSPHERE (sublimation)
  iFrom[5]=iSWE;           iTo[5]=iSnowDrift;    //rates[5]: SNOW->SNOW_DRIFT (blowing snow)
}
///////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvPrairieBlowingSnow::~CmvPrairieBlowingSnow(){}

///////////////////////////////////////////////////////////////////
/// \brief Function to initialize CmvPrairieBlowingSnow objects
//
void CmvPrairieBlowingSnow::Initialize()
{
  int nConn;
  int nHRUs_sub;

  int iDrift    =pModel->GetStateVarIndex(SNOW_DRIFT);
  int iSWE      =pModel->GetStateVarIndex(SNOW);
  int iSnowDepth=pModel->GetStateVarIndex(SNOW_DEPTH);
  int iSnowTemp =pModel->GetStateVarIndex(SNOW_TEMP);

  //sift through all HRUs, determine total number of connections
  //--------------------------------------------------
  nConn=0;
  for(int p=0;p<_pModel->GetNumSubBasins();p++)
  {
    nHRUs_sub=_pModel->GetSubBasin(p)->GetNumHRUs();
    nConn+=(nHRUs_sub*nHRUs_sub); //all snow can drift from all HRUs to all HRUs
  }
  nConn+=2*_pModel->GetNumHRUs();//handles temperature and snow depth (HRU to self)
  DynamicSpecifyLatConnections(nConn);

  //specify from/to HRU indices (k) and state variable indices (i) for all connections
  //--------------------------------------------------
  int kto,kfrom;
  int q=0; //connections counter
  for(int p=0;p<_pModel->GetNumSubBasins();p++)
  {
    for(int ks=0; ks<_pModel->GetSubBasin(p)->GetNumHRUs(); ks++)
    {
      kfrom=_pModel->GetSubBasin(p)->GetHRU(ks)->GetGlobalIndex();
      for(int kss=0; kss<_pModel->GetSubBasin(p)->GetNumHRUs(); kss++)
      {
        kto=_pModel->GetSubBasin(p)->GetHRU(kss)->GetGlobalIndex();
        _kFrom   [q]=kfrom;
        _kTo     [q]=kto;
        _iFromLat[q]=iDrift;
        _iToLat  [q]=iSWE;
        q++;
      }
    }
  }
  _nDriftConnections=q;
  for(int k=0;k<_pModel->GetNumHRUs();k++)
  {
    _kFrom   [q]=_kTo   [q]=k;
    _iFromLat[q]=_iToLat[q]=iSnowDepth;
    q++;
  }
  for(int k=0;k<_pModel->GetNumHRUs();k++)
  {
    _kFrom[q]=_kTo[q]=k;
    _iFromLat[q]=_iToLat[q]=iSnowTemp;
    q++;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for baseflow algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by baseflow algorithm (size of aP[] and aPC[])
//
void CmvPrairieBlowingSnow::GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const
{
    nP=4;
    aP[0]="VEG_DIAM";   aPC[0]=CLASS_VEGETATION;
    aP[1]="VEG_DENS";   aPC[1]=CLASS_VEGETATION;
    aP[2]="VEG_MBETA";  aPC[2]=CLASS_VEGETATION;
    aP[3]="FETCH";      aPC[3]=CLASS_LANDUSE;

    aP[nP]="MAX_HEIGHT";    aPC[nP]=CLASS_VEGETATION; nP++;
    aP[nP]="RELATIVE_HT";   aPC[nP]=CLASS_VEGETATION; nP++;
    aP[nP]="MAX_LAI";       aPC[nP]=CLASS_VEGETATION; nP++;
    aP[nP]="RELATIVE_LAI";  aPC[nP]=CLASS_VEGETATION; nP++;
    aP[nP]="MAX_LEAF_COND"; aPC[nP]=CLASS_VEGETATION; nP++;
    aP[nP]="FOREST_SPARSENESS";    aPC[nP]=CLASS_LANDUSE; nP++;
    aP[nP]="ROUGHNESS";     aPC[nP]=CLASS_LANDUSE; nP++;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns participating state variable list
///
/// \param btype [in] Baseflow algorithm type
/// \param *aSV [out] Array of state variable types needed by baseflow algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by baseflow algorithm (size of aSV[] and aLev[] arrays)
//
void CmvPrairieBlowingSnow::GetParticipatingStateVarList(pbsm_type  stype,sv_type *aSV, int *aLev, int &nSV)
{
  nSV=7;
  aSV[0]=SNOW_AGE;      aLev[0]=DOESNT_EXIST;
  aSV[1]=SNOW_DEPTH;    aLev[1]=DOESNT_EXIST;
  aSV[2]=SNODRIFT_TEMP; aLev[2]=DOESNT_EXIST;
  aSV[3]=SNOW_LIQ;      aLev[3]=DOESNT_EXIST;
  aSV[4]=SNOW;          aLev[4]=DOESNT_EXIST;
  aSV[5]=SNOW_DRIFT;    aLev[5]=DOESNT_EXIST;
  aSV[6]=ATMOSPHERE;    aLev[6]=DOESNT_EXIST;
}
//////////////////////////////////////////////////////////////////
/// \brief Calculate the probability of blowing snow occurence according to Li and Pomeroy (2000).
/// \notes Based upon MESH PBSM model coded by Matthew MacDonald
/// \param snow_depth  [in] snow depth [m]
/// \param T           [in] air temperature [C]
/// \param snowfall    [in] snowfall [mm/d]
/// \param Uten_Prob   [in] probability of blowing snow at ten meters height being less than Uten_prob
//
/// \param wind_thresh [out] threshold wind speed [m/s]
/// \param snow_age    [out] snow age [d]
/// \return probability of blowing snow occurence
//
double CmvPrairieBlowingSnow::ProbabilityThreshold( const double &snow_depth, //snow depth, [mm]
                                                    const double &T,          //air temperature [deg C]
                                                    const double &snowfall,   //[mm/d]
                                                    const double &Uten_Prob,  //[m/s]
                                                          double &wind_thresh,//Threshold wind speed [m/s]
                                                          double &snow_age,   //snow age [d]
                                                    const double &tstep) const
{
  double Mean(0.0),Variance(7.0);
  double Probability(0.0);
  bool   snow_is_dry=(snow_age<REAL_SMALL);

  wind_thresh=9.43+0.180*T+0.00330*T*T;    //[m/s] (overriden for wet snow)
  Mean       =11.0+0.365*T+0.00706*T*T+0.91*log(snow_age*HR_PER_DAY);
  Variance   =4.23+0.145*T+0.00196*T*T;

  if(snow_depth<=0.0) //no snow available
  {
    snow_age=0.0;
    return 0.0;
  }

  if(T<FREEZING_TEMP)
  {
    if     (snowfall>=0.0){snow_age= tstep;}// with concurrent snowfall: new dry snow
    else if(snow_is_dry  ){snow_age+=tstep;}// without concurrent snowfall: old dry snow
  }
  else if((T>=FREEZING_TEMP) || (!snow_is_dry)) //wet snow
  {
    snow_age=0.0;
    wind_thresh = 9.9;
    Mean=21.0;
    Variance=7.0;
  }

  if(Uten_Prob>3.0)// wind<3 m/s too weak for dry snow transport
  {
    double wind=0.0;
    double dw=0.1;
    while(wind<=Uten_Prob)
    {
      wind+=dw;
      Probability+=(1.0/(Variance*sqrt(2.0*PI)))*(exp(-0.5*pow((wind - Mean)/Variance,2.0)))*dw;//Ugly/slow way to do this - should be able to invert probability formula
    }
  }
  //cout<<"Probability "<<Probability<<" "<<Uten_Prob<<" "<<wind_thresh<<endl;
  return Probability;
}

//////////////////////////////////////////////////////////////////
/// \brief calculates drift and sublimation rates from blowing snow
// Single column calculations for blowing snow transport and sublimation.
// Ported over from FORTRAN MESH code by Matthew MacDonald (PBSMRates.F)
// Equation numbers refer to JW Pomeroy thesis (1988).
///
/// \param E_StubHt [in] [m] stubble height
/// \param Uthr [in] [m/s] threshold wind speed
/// \param T [in] air temperature
/// \param u [in] [-] wind speed
/// \param rel_hum [in] [-] relative humidity
/// \param fetch [in] [m] fetch distance
/// \param veg_dens [in] [count/m2] Vegetation density
/// \param veg_diam [in] [m] Vegetation diameter
/// \param mBeta [in] [-] unitless parameter
///
/// \param DriftH [out] [kg/m/s]
/// \param SublH [out] [kg/m^2/s]
//
// Li L, Pomeroy JW. 1997. Probability of occurrence of blowing snow. Journal of Geophysical Research 102: 21955-21964.
// Raupach MR, Gillette DA, Leys JF. 1993. The effect of roughness elements on wind erosion threshold. Journal of Geophysical Research 98: 3023-3029.
// Pomeroy JW. 1988. Wind transport of snow . Ph.D. Thesis, University of Saskatchewan.
/// JRC: Verified via comparison to CRHM code- consistent.
//
void CmvPrairieBlowingSnow::PBSMrates(const double E_StubHt, // stubble height [m]
                                      const double Uthr,     // threshold wind speed [m/s]
                                      const double T,        // air temperature [deg C]
                                      const double u,        // wind speed [m/s]
                                      const double rel_hum,  // relative humidity [0..1]
                                      const double Fetch,    // [m] (param)
                                      const double veg_dens, // [count/m2] Vegetation density
                                      const double veg_diam, // [m] Vegetation diameter
                                      const double mBeta,    // [-] unitless parameter
                                            double &DriftH,  //[kg/m/s]
                                            double &SublH) const  //[kg/m2/s] //per half hour?

{
  const double REF_FETCH=300; //XD [m]
  const double ZD=0.3; //[m]

  const double C1=2.8;
  const double C2=1.6;
  const double C3=4.2;

  //Compute stubble coefficients
  //double z_stb=0.0048*E_StubHt*100.0;  // Lettau, used for susp Z0
  double z_stb = 0.5*veg_dens*veg_diam*E_StubHt;  // [-] Essery et al (1999) from Lettau (1969)

  double SBsalt=0.0; // Sublimation in saltation layer [kg/m2/s]
  double Qsalt=0.0;  // Blowing snow flux in saltation layer  [kg/m/s]
  double SBsum=0.0;  // Sublimation above saltation layer [kg/m2/s]
  double Qsum=0.0;   // Blowing snow flux above saltation layer [kg/m/s]

  DriftH=0.0;
  SublH=0.0;
  if(u>Uthr)
  {
    double Usthr=0.03697*Uthr;           //{Eq. 6.3    }
    double Ustar=0.02264*pow(u,1.295);   //{Eq. 6.2 rev}

    //Raupach
    double RaupachTerm=1.0; //from Raupach1993
    double Sigma =(PI*veg_diam)/(4.0*E_StubHt);     // [-] Raupach Eq. 4
    double Lambda=veg_dens*veg_diam*E_StubHt;       // [-] Raupach Eq. 1 (frontal area index) (Eqn 4 MacDonald et al, 2009)
    if(E_StubHt>0.0001) {
      RaupachTerm=1.0/((1.0-Sigma*Lambda)*(1.0+mBeta*Lambda));
    }

    double Nsalt=2.0*DENSITY_AIR/(C2*C3*Ustar)*(RaupachTerm-(Usthr*Usthr)/(Ustar*Ustar));// [kg/m3] Pomeroy1988 Eq. 4.14 updated
    if(Nsalt<=0.0) {SublH=DriftH=0.0;return;}

    //-------------------------------------------------------------------------
    // calculate sublimation & drift rate in the saltation layer
    //-------------------------------------------------------------------------

    Qsalt  =C1*DENSITY_AIR*Usthr/(GRAVITY*C3*Ustar)*(Ustar*Ustar*RaupachTerm-Usthr*Usthr);// (should be [kg/m/s]; is [kg/m2]) Pomeroy1988 Eq. 4.20 (Eqn 2 MacDonaldEtAl2009)
    //UNITS DONT WORK OUT IN MESH CODE - DIVISION BY C3*Ustar not in Pomeroy1988

    double Mpr, alpha, rel_hum_z, Vsalt,Hsalt;
    Mpr=0.0001;                           // mean particle radius [m]
    alpha=5.0;                            // particle size distribution shape parameter
    Hsalt=C2/(2.0*GRAVITY)*Ustar*Ustar;   // [m] maximum saltation height {Pomeroy 1988 Eq. 4.13}
    rel_hum_z=(rel_hum-1.0)*(1.019+0.027*log(Hsalt));      // Pomeroy1988 Eq. 6.20
    upperswap(rel_hum_z,-0.01);
    Vsalt=0.6325*Ustar+2.3*Usthr;         // Pomeroy1988 Eq. 6.25

    SBsalt=SublimRateCoefficient(Mpr,alpha,Vsalt,rel_hum_z,T)*Nsalt*Hsalt;  // [kg/m2/s] Pomeroy1988 Eq. 6.11,6.13

    //-------------------------------------------------------------------------
    // calculate integrated mass blowing snow and sublimation flux in the suspended layers (height r to Bound)
    //-------------------------------------------------------------------------

    // Loop to find the first suspended drift density level, z from the reference level z_ref
    //-------------------------------------------------------------------------
    double Nz,z,z_ref;
    double dz=0.0001;
    z=z_ref=(0.05628*Ustar); //reference height [m] // Pomeroy1988 Eq. 5.27;
    while(z<=0.15)
    {
      Nz=0.8*exp(-1.55*(pow(z_ref,-0.544)-pow(z,-0.544))); // [kg/m3] Suspended level drift density (Pomeroy1988 Eq. 5.26)
      z+=dz;
      if(Nz<=Nsalt){ break; }  //drift density is less than or equal to Nsalt.
    }

    // find height of fully-developed boundary layer for turbulent diffusion
    //-------------------------------------------------------------------------
    double Bound; //[m]
    double Bd=1.0;//initial guess [m]
    double term=162.926/(Ustar*Ustar);

    Bound=ZD+VON_KARMAN*VON_KARMAN*(Fetch-REF_FETCH)*pow(log(Bd*term)*log(ZD*term),-0.5);// Pomeroy1988 Eq. 6.6
    while(fabs(Bound-Bd)>0.001)
    {
      Bd=Bound;
      Bound=ZD+VON_KARMAN*VON_KARMAN*(Fetch-REF_FETCH)*pow(log(Bd*term)*log(ZD*term),-0.5);// Pomeroy1988 Eq. 6.9
    }
    if(Fetch<REF_FETCH){ Bound=ZD; }

    // Calculate the suspended mass flux up to 5 metres
    // and the total sublimation rate to the top of the boundary layer
    // at increments of 1 mm to z=50cm & increments of 10 cm to z=Bound
    //-------------------------------------------------------------------------
    dz=0.001;
    z+=dz;
    double Uz,Vsusp;
    while(z<=Bound)
    {
      Nz    =0.8*exp(-1.55*pow(z_ref,-0.544)-pow(z,-0.544));
      Uz    =(Ustar*pow(1.2/(1.2+Nz),0.5)/VON_KARMAN)*log(z/((0.00613*(Ustar*Ustar))+z_stb));// Pomeroy1988 Eq. 4.17r,  Eq. 5.17a

      if(Uz>0.0)
      {
        Mpr= (4.6E-5)*pow(z,-0.258);                   // mean particle radius Pomeroy1988  Eq. 6.15
        upperswap(Mpr,30e-6);

        alpha=4.08+12.6*z;                             // Pomeroy1988 Eq. 6.14
        lowerswap(alpha,25);

        rel_hum_z=(rel_hum-1)*(1.019+0.027*log(z));    // Pomeroy1988 Eq. 6.20
        upperswap(rel_hum_z,-0.01);

        Vsusp=(1.1E7)*pow(Mpr,1.8)+0.0106*pow(Uz,1.36);// Pomeroy1988 Eq. 5.18

        SBsum+=SublimRateCoefficient(Mpr,alpha,Vsusp,rel_hum_z,T)*Nz*dz;  // Pomeroy1988 Eq. 6.12/6.13

        if(z<5.0) { Qsum+=(Nz*Uz)*dz; }               // Pomeroy1988 Eq. 5.4

        if(Nz<=0.00001) { //drift density small enough to finish
          SublH =-min(SBsum+SBsalt,0.0); //[kg/m^2/s]
          DriftH=(Qsalt+Qsum);           //[kg/m/s]
          return;
        }
        else{
          if(((z-dz)>=0.5) && (z<0.6)) { dz=0.1;  z=0.5; } //change step size to dz=0.1 once z>0.5
        }
      }
      else {
        SublH =0.0; //[kg/m^2/s]
        DriftH=0.0; //[kg/m/s]
        return;
      }

      z+=dz;

    }/*End while H<Hbound*/
  }/*end if u>Uthr*/

  SublH=-min(SBsum+SBsalt,0.0); // [kg/m^2/s]
  DriftH=(Qsum+Qsalt);          // [kg/m/s]
  //DriftH=0.0;//TMP DEBUG
  return;
}

//////////////////////////////////////////////////////////////////
/// \brief provides change in snow age, depth, SWE, liquid content, drift amount, and  drift temp due to sublimation/blowing snow
/// \param *storage [in] Array of state variable values for this HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of loss from baseflow [mm/d]
/// \param Options [in] model options structure
//
// Based upon MESH PBSM model coded by Matthew MacDonald (routine PBSMrun)
// Single column calculations for blowing snow transport and sublimation. Based on
// JW Pomeroy thesis(1988; UofS),Pomeroy et al. (1993; JH),and Pomeroy and Li(2000; JGR).
//
void CmvPrairieBlowingSnow::GetRatesOfChange(const double              *state_vars,
                                            const CHydroUnit  *pHRU,
                                            const optStruct   &Options,
                                            const time_struct &tt,
                                            double      *rates) const
{
  double Drift=0.0;   // [kg/m2] drift losses over time step
  double Subl=0.0;    // [kg/m2] sublimation losses over time step

  int iSWE      =pModel->GetStateVarIndex(SNOW);
  int iSnowLiq  =pModel->GetStateVarIndex(SNOW_LIQ);
  int iSnowDepth=pModel->GetStateVarIndex(SNOW_DEPTH);
  int iSnowTemp =pModel->GetStateVarIndex(SNOW_TEMP);
  int iSnowAge  =pModel->GetStateVarIndex(SNOW_AGE);
  int iSnowDriftTemp=pModel->GetStateVarIndex(SNODRIFT_TEMP);

  double SWE       =state_vars[iSWE];      // [mm]
  double snow_liq  =state_vars[iSnowLiq];  // [mm]
  double snow_depth=state_vars[iSnowDepth];// [mm] depth of snowpack
  double snow_temp =state_vars[iSnowTemp]; // [C]
  double snow_age  =state_vars[iSnowAge];  // [d]

  double snow_dens =(SWE/snow_depth)*DENSITY_ICE; //kg/m3

  //Get forcings
  const force_struct *F=pHRU->GetForcingFunctions();
  double u_meas =F->wind_vel;                              // [m/s] wind velocity @ 2m

  u_meas=10;// m/s TMP DEBUG:  override

  //Get parameters
  double fetch   =pHRU->GetSurfaceProps()->fetch;          // [m] fetch distance
  double veg_dens=pHRU->GetVegetationProps()->veg_dens;    // [count/m2] vegetation density
  double veg_diam=pHRU->GetVegetationProps()->veg_diam;    // [m] vegetation diameter (why is this m2 in MESH documentation?)
  double mBeta   =pHRU->GetVegetationProps()->veg_mBeta;   // [-] ?? ~170 or 32 if FCS>FGS
  double veg_ht  =pHRU->GetVegVarProps()->height;          // [m] vegetation height
  double meas_ht =pHRU->GetVegVarProps()->reference_height;// [m] wind speed measurement height
  double z0_mom  =pHRU->GetVegVarProps()->roughness;       // [m] momentum roughness height
  double zero_pl =pHRU->GetVegVarProps()->zero_pln_disp;   // [m] zero plane displacement

  veg_ht=1.0;//TMP DEBUG
  fetch=200;
  z0_mom=0.05;
  zero_pl=0.5;

  //cout<<"PBSM Params: "<<fetch<<" "<<veg_dens<<" "<<mBeta<<" "<<veg_ht<<" "<<ZREFM<<" "<<z0_mom<<endl;
  if(snow_depth>REAL_SMALL)
  {
    //===============================================================================
    //Set values for mB for partitioning shear stress over vegetation
    //(different for vegetation categories; see MacDonald,Pomeroy & Pietroniro(2009,Hydrol. Proc.))

    double E_StubHt;   //height of vegetation above snowpack [m]
    double z0;         //roughness length for momentum over snow/vegetation [m]
    double u10;        //vel @ 10m [m/s]
    double Ustar;      //friction velocity [m/s]
    double Uten_Prob;

    // HRU-level snow transport & sublimation calculations depths(m),SWE(mm; kg/m^2)
    E_StubHt=veg_ht-(snow_depth/MM_PER_METER);
    upperswap(E_StubHt,0.0001);

    //z0=exp(zOMLCS[i])

    z0=E_StubHt*2/3;
    if(snow_depth>0.0){z0=z0_mom;}

    z0=max(z0,meas_ht*0.5);//added by JRC

    u10=u_meas*log(10.0/z0)/log((meas_ht)/z0); //assumes z0<10, z0<ZREFM (assumes no zero plane displacement!)
    //u10=F->wind_vel*max(log((10.0-zero_pl)/z0)/log((meas_ht-zero_pl)/z0),0.0); //JRC preferred: meas_ht mus be larger than zdp+z0!

    //cout<<" u10: "<<u10<<" "<<z0<<" "<<meas_ht<<endl;
    //ExitGracefullyIf(u10<0,"Blowing Snow: u10<0: bad reference heights",RUNTIME_ERR);

    Ustar=0.02264*pow(u10,1.295); //Eq. 6.2 rev. Pomeroy1988,friction velocity over fallow

    //Calculate Uten_prob
    Uten_Prob=u10;
    if(E_StubHt>0.01)
    {
      double znod,Lambda,Ustn;
      znod=pow(Ustar,2)/163.3+0.5*E_StubHt*veg_dens*veg_diam;//> Eq. 29,Snowcover Accumulation,Relocation & Management book(1995)
      Lambda=veg_dens*veg_diam*E_StubHt;                     //> (frontal area index) Raupach Eq. 1
      Ustn=Ustar*sqrt((mBeta*Lambda)/(1.0+mBeta*Lambda));
      Uten_Prob=(log(10.0/znod))/VON_KARMAN*sqrt(Ustar-Ustn);
      //cout<<" Uten_Prob: "<<znod<<" "<<Ustar<<" "<<Ustn<<endl;
    }

    //cout<<" Uten_Prob: "<<u10<<" "<<z0<<" "<<meas_ht<<" "<<Ustar<<endl;

    // Calculate probability of blowing snow occurence (also determines Uthresh, updates snow age)
    double Prob,Uthresh(0.0);
    double DriftH(0.0),SublH(0.0);
    double snowfall=F->precip*(F->snow_frac);
    const double MIN_PROB=1e-3;

    Prob=ProbabilityThreshold(snow_depth,F->temp_ave,snowfall,Uten_Prob,Uthresh,snow_age,Options.timestep);

    //Uthresh*=0.8; //JRC: From MESH code - why?

    if(Prob>MIN_PROB)
    {
      // Single column calculations of blowing snow transport & sublimation
      PBSMrates(E_StubHt,Uthresh,F->temp_ave,u_meas,F->rel_humidity,fetch,veg_dens,veg_diam,mBeta,DriftH,SublH);// calculates DriftH, SublH

      DriftH*=Options.timestep*SEC_PER_DAY; //rates converted to incremental drift (NOT IN ORIGINAL MESH CODE//)
      SublH *=Options.timestep*SEC_PER_DAY;
      //
      Drift=DriftH*Prob/fetch; //[kg/m2]
      Subl =SublH *Prob;       //[kg/m2]
    }//end if (Prob>MIN_PROB)

  }//end  if (snow_depth>REAL_SMALL)

  //===============================================================================
  //
  //> Recalculate subarea snow properties after snow transport
  double HTCS;
  double old_energy_content(0),new_energy_content(0);
  if((Drift+Subl)>0.0)//snow mass loss is sum of transport + sublimation
  {
    if(snow_depth>0.0)
    {
      double snow_mass =(SWE/MM_PER_METER)*snow_dens;  //kg/m2
      if((Drift+Subl)>snow_mass) {// corrects for insufficient snowpack to support calculated drift/subl amounts
        Subl =snow_mass* Subl/(Subl+Drift);
        Drift=snow_mass*Drift/(Subl+Drift);
      }

      old_energy_content=CalculateSnowEnergyContent(SWE,snow_depth,snow_liq,snow_temp);

      snow_depth=max(0.0,snow_depth-((Drift+Subl)/snow_dens)*MM_PER_METER);//[mm]

      if (snow_depth==0){ snow_liq=0.0;/*Subl+=snow_liq*DENSITY_WATER/MM_PER_METER;*/} //JRC : don't like this - should be separate loss from snow_liq->atm

      SWE=snow_depth*(snow_dens/DENSITY_ICE); //no change to density

      new_energy_content=CalculateSnowEnergyContent(SWE,snow_depth,snow_liq,snow_temp);
    }
  }//end if ((Drift+Subl)>0.0)

  //assumes no change in density
  //update variables
  rates[0]=(snow_age  -state_vars[iSnowAge      ])/Options.timestep; //snow age->snow_age
  rates[1]=(snow_depth-state_vars[iSnowDepth    ])/Options.timestep; //snow_depth->snow_depth
  rates[2]=(snow_temp -state_vars[iSnowDriftTemp])/Options.timestep; //drift temp->drift temp

  rates[3]=(snow_liq*MM_PER_METER/DENSITY_WATER -state_vars[iSnowLiq ])/Options.timestep; //snow_liq->atmosphere (via sublimation)

  double denom=Subl+Drift;
  if(Subl+Drift==0){ denom=1.0; }
  rates[4]=- Subl/denom*(SWE -state_vars[iSWE])/Options.timestep; //snow->atmosphere (via sublimation)
  rates[5]=-Drift/denom*(SWE -state_vars[iSWE])/Options.timestep; //snow->drifting snow (via drift)

  rates[4]=0.0; //TMP DEBUG TO DETERMINE SUBL AMT
  rates[5]=-(SWE -state_vars[iSWE])/Options.timestep; //snow->drifting snow (via drift)

  HTCS=(new_energy_content-old_energy_content)/Options.timestep; //[MJ/m2/d] energy flux from delta snow depth

  return;
}

void CmvPrairieBlowingSnow::ApplyConstraints(const double      *state_vars,
                                             const CHydroUnit  *pHRU,
                                             const optStruct   &Options,
                                             const time_struct &tt,
                                             double      *rates) const
{
  //nothing for now - constraints handled in GetRatesOfChange()
}

//////////////////////////////////////////////////////////////////
/// \brief calculates sublimation/ saltation rate coefficient
/// \param Mpr       [in] [m] Mean snow particle radius
/// \param alpha     [in] [-] gamma shape parameter for blowing snow particle distribution
/// \param Vsalt     [in] [m/s] ventilation/saltation/sublimation velocity
/// \param rel_hum_z [in] [-] undersaturation of relative humidity at height z (<0)
/// \param T         [in] [C] air temperature
//
/// returns sublimation/saltation rate coefficient, [1/s]
/// ***JRC: VERIFIED AGAINST CRHM Classpbsm_M::Pbsm() -issue with Lamb calculation***
//
double CmvPrairieBlowingSnow::SublimRateCoefficient(const double &Mpr,
                                                   const double &alpha,
                                                   const double &Vsalt,
                                                   const double &rel_hum_z,
                                                   const double &T) const//[C]
{
  const double MMM  =18.01;
  const double RR   =8313.0;
  const double QSTAR=120.0;
  const double LATH =2.838e6;
  const double KIN_VISC=1.88E-5; // [m2/s] kinematic viscosity of atmos.

  double Es          =PA_PER_KPA*GetSaturatedVaporPressure(T); //[Pa]
  double sat_vap_dens=(Es*MMM)/(RR*(T+ZERO_CELSIUS));// [g/m3]?
  double Diff        =2.06e-5*pow((T+ZERO_CELSIUS)/ZERO_CELSIUS,1.75);// diffus. of w.vap. atmos. (m^2/s)
  double Lamb        =0.00063*(T+ZERO_CELSIUS+0.0673);                // therm. cond. of atm. (J/(msK))
  //double Lamb      =0.000076843*(T+ZERO_CELSIUS) + 0.003130762; //from CRHM

  double Htran, Reyn,Nuss, A,B,C,DmDt,Mpm;
  Htran=0.9*PI*(Mpr*Mpr)*QSTAR;
  Reyn =(2.0*Mpr*Vsalt)/KIN_VISC;    // [-] Pomeroy1988 Eq. 6.22
  Nuss =1.79+0.606*sqrt(Reyn);       // [-] Pomeroy1988Eq. 6.21
  A    =Lamb*(T+ZERO_CELSIUS)*Nuss;
  B    =LATH*MMM/(RR*(T+ZERO_CELSIUS))-1.0;
  C    =1.0/(Diff*sat_vap_dens*Nuss);
  DmDt =((2.0*PI*Mpr*rel_hum_z)-(Htran*B/A))/((LATH*B/A)+C);

  Mpm  =4.0/3.0*PI*DENSITY_ICE*(Mpr*Mpr*Mpr)*(1.0+3.0/alpha+2.0/(alpha*alpha)); // [kg] mean particle mass {Pomeroy Eq. 6.16} {Gamma Dist. Corr.}

  return DmDt/Mpm;
}

//////////////////////////////////////////////////////////////////
/// \brief weighted average of two variables v1 and v2 with weights
//
double wt_average(const double &v1,const double &v2,const double &w1,const double &w2)
{
  if((w1+w2)==0){ return 0.5*(v1+v2); }//equal zero weights
  return (w1*v1+w2*v2)/(w1+w2);
}

//////////////////////////////////////////////////////////////////
/// \brief returns lateral blowing snow exchange rates ([mm-m2/day])
/// \param **state_vars [in] 2D array of current state variables [nHRUs][nSVs]
/// \param **pHRUs [in] array of pointers to HRUs
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *exchange_rates [out] Rate of loss from "from" compartment [mm-m2/day]
//
void CmvPrairieBlowingSnow::GetLateralExchange( const double * const     *state_vars, //array of all SVs for all HRUs, [k][i]
                                                const CHydroUnit * const *pHRUs,
                                                const optStruct          &Options,
                                                const time_struct        &tt,
                                                      double             *exchange_rates) const
{
  const double DRIFT_DENS=300.;//[kg/m3]
  const double DRIFT_HCP =HCP_ICE*(DRIFT_DENS/DENSITY_ICE);//[MJ/m3/K]

  int nHRUs     =_pModel->GetNumHRUs();
  int nSubBasins=_pModel->GetNumSubBasins();

  //model inputs/outputs
  int iDrift     =pModel->GetStateVarIndex(SNOW_DRIFT);
  //int iDriftTemp =pModel->GetStateVarIndex(SNODRIFT_TEMP);
  //int iSnowTemp  =pModel->GetStateVarIndex(SNOW_TEMP);
  int iSnowDepth =pModel->GetStateVarIndex(SNOW_DEPTH);
  int iSWE       =pModel->GetStateVarIndex(SNOW);
  int iSnowLiq   =pModel->GetStateVarIndex(SNOW_LIQ);

  int k;
  const CSubBasin *pBasin;
  double RemainingDrift;  //[mm SWE]
  double TotalDrift;      //[mm SWE]
  double drift_tempSB;    //drifting snow temperature in basin
  double deltaDrift;
  double drift;
  double area_frac;       //fractional coverage of HRU in subbasin //FARE
  double SWE,snow_depth,snow_temp,snow_liq;
  int q=0; //connections counter
  for(int p=0;p<nSubBasins;p++)
  {
    pBasin=_pModel->GetSubBasin(p);

    RemainingDrift=0.0;
    TotalDrift    =0.0;
    drift_tempSB  =FREEZING_TEMP;

    //Determine total drifting snow in each subbasin (RemainingDrift) and its temperature (drift_tempSB)
    //-----------------------------------------------------------------------------------------------------
    for(int nn=0;nn<pBasin->GetNumHRUs();nn++)
    {
      k         =pBasin->GetHRU(nn)->GetGlobalIndex();
      area_frac =(pHRUs[k]->GetArea()/pBasin->GetBasinArea());//fractional coverage of HRU in subbasin

      drift       =state_vars[k][iDrift];
      deltaDrift  =drift*area_frac; //[mm SWE] average over basin
      drift_tempSB=wt_average(drift_tempSB,drift,TotalDrift,deltaDrift); //[deg C] - calculated drifting snow temp
      TotalDrift  +=deltaDrift; // total snow drift in subbasin
    }
    RemainingDrift=TotalDrift;


    double HCPS,added,dist_k,area_k;
    int nnn;
    for(int nn=0;nn<pBasin->GetNumHRUs();nn++)
    {
      //nnn=_BlowingSnowSortOrder[p][nn];
      nnn=nn;//Assumes HRUs are pre-sorted (not needed if JRC mods below are used)

      k         =pBasin->GetHRU(nnn)->GetGlobalIndex();
      area_k    =pHRUs[k]->GetArea();
      dist_k    =pHRUs[k]->GetSurfaceProps()->bsnow_distrib; //input parameter - sum of distrib must be 1 for each subbasin
      area_frac =(area_k/pBasin->GetBasinArea());
      SWE       =state_vars[k][iSWE];
      snow_depth=state_vars[k][iSnowDepth];
      snow_liq  =state_vars[k][iSnowLiq];
      snow_temp =state_vars[k][iSnowDepth];

      //if(nn==0) { added=max((state_vars[k][iDrift]   )*dist_k,0.0); } //First HRU in subbasin //added SWE - all in first basin falls on first basin?
      //else      { added=max((RemainingDrift/area_frac)*dist_k,0.0); } //Not first HRU
      added=max((TotalDrift/area_frac)*dist_k,0.0); //JRC: above should be this if all drift is to fall

      //Redistribute subbasin snow drift and calculate modified snowpack properties in HRU
      //--------------------------------------------------------------------------------------------------
      HCPS=HCP_ICE*(SWE/snow_depth)+HCP_WATER*(snow_liq/snow_depth);

      double delta_snow_temp =wt_average(snow_temp,drift_tempSB,snow_depth*HCPS,added*(DENSITY_WATER/DRIFT_DENS)*DRIFT_HCP)-snow_temp;
      double delta_snow_depth=added*(DENSITY_WATER/DRIFT_DENS);
      //double delta_SWE       =added;
      RemainingDrift-=added*area_frac;   // remove drift used from total available [mm SWE]

      //double new_energy,old_energy;
      // old_energy=CalculateSnowEnergyContent(SWE          ,snow_depth                 ,snow_liq,snow_temp);
      // new_energy=CalculateSnowEnergyContent(SWE+delta_SWE,snow_depth+delta_snow_depth,snow_liq,snow_temp);
      // HTCS[k]=(new_energy-old_energy)/Options.timestep; //Not currently used

      for(int nn=0;nn<pBasin->GetNumHRUs();nn++)
      {
        int kk=pBasin->GetHRU(nn)->GetGlobalIndex();
        exchange_rates[q]=(state_vars[kk][iDrift]*area_k/Options.timestep)*dist_k; //rate of loss from kk Drift to k SWE, mm-m2/d - should sub to added
        q++;
      }
      exchange_rates[_nDriftConnections      +k]=delta_snow_temp *area_k/Options.timestep; //self-transfer C-m2/d
      exchange_rates[_nDriftConnections+nHRUs+k]=delta_snow_depth*area_k/Options.timestep; //self-transfer
    }//end HRU loop
  }//end subbasin loop
}

//////////////////////////////////////////////////////////////////
/// \brief converts specific humidity (kg/kg) to relative humidity
// (NOT USED)
double SpecHumToRelHum(const double specific_hum,const double T,const double air_press)
{
  double a0,a1,a2,a3,a4,a5,a6;
  // convert specific (kg/kg) to relative humidity (0.xx)
  if(T>FREEZING_TEMP)  {
    // coefficients with respect to watewr
    a0=6.107799961;
    a1=4.436518521E-1;
    a2=1.428945805E-2;
    a3=2.650648471E-4;
    a4=3.031240396E-6;
    a5=2.034080948E-8;
    a6=6.136820929E-11;
  }
  else{
    //coefficients with respect to ice
    a0=6.109177956;
    a1=5.034698970E-1;
    a2=1.886013408E-2;
    a3=4.176223716E-4;
    a4=5.824720280E-6;
    a5=4.838803174E-8;
    a6=1.838826904E-10;
  }
  double eT=a0+T*(a1+T*(a2+T*(a3+T*(a4+T*(a5+T*a6)))));
  double rel_hum=28.9644/(18.01534/specific_hum+28.9644-18.01534)*air_press/100/eT;
  return rel_hum;
}
