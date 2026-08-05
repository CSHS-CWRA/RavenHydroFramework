/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2026 the Raven Development Team
  ------------------------------------------------------------------
  Frozen ground
  ----------------------------------------------------------------*/

#include "FrozenGround.h"
double CalculateThermalConductivity(const double &poro,const double&kappa_soil,const double &sat,const double &Fice);
double CalculateHeatCapacity       (const double &poro,const double&hcp_soil,  const double &sat,const double &Fice);
double CalculateDensity            (const double &poro,const double&rho_soil,  const double &sat,const double &Fice);
/*****************************************************************
   FrozenGround Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the standard constructor
/// \param groundfreeze_type [in] Model of ground freeze selected
/// \param pModel [in] model
//
CmvFrozenGround::CmvFrozenGround(groundfreeze_type type, CModelABC *pModel)
  :CHydroProcessABC(GROUND_FREEZING,pModel)
{
  _type =type;
  if (_type == FREEZE_STEFAN) {
    CHydroProcessABC::DynamicSpecifyConnections(1);
    iFrom[0]=THAW_DEPTH; iTo[0]=THAW_DEPTH;
  }
  else if (_type == FREEZE_RANKINEN) {
    CHydroProcessABC::DynamicSpecifyConnections(2);
    iFrom[0]=SOIL_TEMP;     iTo[0]=SOIL_TEMP;
    iFrom[1]=SOIL_PCT_FROZ; iTo[1]=SOIL_PCT_FROZ;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of default destructor
//
CmvFrozenGround::~CmvFrozenGround(){}

//////////////////////////////////////////////////////////////////
/// \brief does nothing
//
void CmvFrozenGround::Initialize()
{
}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for capillary rise algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by capillary rise  algorithm (size of aP[] and aPC[])
//
void CmvFrozenGround::GetParticipatingParamList(string  *aP , class_type *aPC , int &nP) const
{
  nP=0;
  if (_type==FREEZE_STEFAN)
  { 
    nP=2;
    aP [0]="THERMAL_CONDUCTIVITY";     aPC [0]=CLASS_SOIL;
    aP [1]="POROSITY";                 aPC [1]=CLASS_SOIL;
  }
  else if(_type==FREEZE_RANKINEN)
  { 
    nP=0;
    aP [nP]="THERMAL_CONDUCTIVITY";     aPC [nP]=CLASS_SOIL; nP++;
    aP [nP]="HEAT_CAPACITY";            aPC [nP]=CLASS_SOIL; nP++;
    aP [nP]="BULK_DENSITY";             aPC [nP]=CLASS_SOIL; nP++;
    aP [nP]="POROSITY";                 aPC [nP]=CLASS_SOIL; nP++;
    aP [nP]="SNOW_DAMPEN_COEFF";        aPC [nP]=CLASS_LANDUSE; nP++;
  }
  else if (_type == FREEZE_THERMAL) 
  {
    ExitGracefully("FREEZE_THERMAL",STUB);
  }
  else
  {
    ExitGracefully("CmvFrozenGround::GetParticipatingParamList: undefined frozen ground algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
/// \details User specifies from and to compartments, levels not known before construction
///
/// \param groundfreeze_type [in] algorithm used
/// \param *aSV  [out] Reference to state variable types needed by algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV  [out] Number of state variables required by algorithm (size of aSV[] and aLev[] arrays)
//
void CmvFrozenGround::GetParticipatingStateVarList(groundfreeze_type  gftype,sv_type *aSV, int *aLev, int &nSV)
{
  if (gftype==FREEZE_STEFAN)
  { 
    nSV=0;
    aSV[nSV]=THAW_DEPTH; aLev[nSV]=0; nSV++;
    aSV[nSV]=SOIL;       aLev[nSV]=0; nSV++; //requires at least top layer of soil moisture
  }
  else if(gftype==FREEZE_RANKINEN)
  { 
    nSV=0;
    aSV[nSV]=SOIL_TEMP;     aLev[nSV]=0; nSV++;
    aSV[nSV]=SOIL;          aLev[nSV]=0; nSV++; //requires at least top layer of soil moisture
    aSV[nSV]=SOIL_PCT_FROZ; aLev[nSV]=DOESNT_EXIST; nSV++;
    aSV[nSV]=SNOW;          aLev[nSV]=0; nSV++;
 }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rate of movement of frost table [mm/d]
///
/// \param *state_vars [in] array of state variables
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of loss of water from soil to another soil layer [mm/day]
//
void   CmvFrozenGround::GetRatesOfChange( const double      *state_vars,
                                          const CHydroUnit  *pHRU,
                                          const optStruct   &Options,
                                          const time_struct &tt,
                                                double      *rates) const
{
  if (pHRU->GetHRUType()!=HRU_STANDARD){return;}//Lakes & glaciers

  //--Calculate Rate of change of frost table depth------------
  //-----------------------------------------------------------------
  if (_type==FREEZE_STEFAN)
  { 
    const  soil_struct *pSoil=NULL;
    double thick, stor, max_stor, kappa_s, poro;
    double kappa_m;
    int    iSoil;

    double T_surf=pHRU->GetForcingFunctions()->temp_ave; 
    double z=state_vars[iFrom[0]];

    int m=0;
    double total_thick(0),kappaden(0),ice_sat(0);
    while (total_thick<=z)
    {
      iSoil    = pModel->GetStateVarIndex(SOIL,m);
      if (iSoil==DOESNT_EXIST){break;} //fully frozen

      pSoil    = pHRU->GetSoilProps(m);
      max_stor = pHRU->GetSoilCapacity(m);
      thick    = min(pHRU->GetSoilThickness(m),z-total_thick);
      poro     = pSoil->porosity;
      kappa_s  = pSoil->thermal_cond; //[MJ/m/d/K]

      stor = state_vars[iSoil];  
      stor=min(max(stor,0.0),max_stor); //correct for potentially invalid storage

      total_thick+=thick;
      
      kappa_m  = (1-poro)*kappa_s+(stor/max_stor)*poro*(TC_ICE)+(1-stor/max_stor)*poro*TC_AIR;

      kappaden+=thick/kappa_m; //denominator of harmonic mean

      ice_sat+=poro*(stor/max_stor)*thick;
      m++;
    }

    double kappa=total_thick/kappaden;
    double mean_ice_sat=ice_sat/total_thick; //arithmetic men

    //snow correction to surface temperature
    double Tsnow= pHRU->GetSnowTemperature();
    double Sdepth=pHRU->GetSnowDepth();
    double SWE   =state_vars[pModel->GetStateVarIndex(SNOW)];
    double rho_sno=SWE/Sdepth*DENSITY_WATER*GPCM3_PER_KGPM3; //[g/cm3]
    double kappa_snow=pow(10.0,2.65*rho_sno-1.652)/MJ_PER_D_TO_WATT;//Sturm et al. 1997, eqn 7
    if ((z>0) && (SWE>0)){
      T_surf=Tsnow*(1.0/(1.0+kappa/kappa_snow*Sdepth/2.0/z));
    }

    //dz/dt=kappa*T_surf/z/rho/f/L

    rates[0]=kappa*T_surf/z/LATENT_HEAT/DENSITY_ICE/mean_ice_sat;

    rates[0]=min(rates[0],(total_thick-z)/Options.timestep); //ensures can't freeze below max soil depth

    //must re-set annually - this is a temporary hack
    if (T_surf>5.0){
      rates[0]=-z/Options.timestep; 
    }
    ExitGracefully("FREEZE_STEFAN: STILL IN TESTING",STUB);
  }
  //-----------------------------------------------------------------
  else if(_type==FREEZE_RANKINEN) 
  {
    int iSlTemp=pModel->GetStateVarIndex(SOIL_TEMP,0);
    int iSoil  =pModel->GetStateVarIndex(SOIL,0);
    int iPctFrz=pModel->GetStateVarIndex(SOIL_PCT_FROZ);

    double Tair=pHRU->GetForcingFunctions()->temp_ave;
    
    double soil_temp=state_vars[iSlTemp];
    double pct_froz =state_vars[iPctFrz];
    double soil_wc  =state_vars[iSoil];
    double snow_depth=pHRU->GetSnowDepth()/MM_PER_METER; //[m]

    double thickness =pHRU->GetSoilThickness(0)/MM_PER_METER; //[m]
    double kappa_s   =pHRU->GetSoilProps(0)->thermal_cond;    //[MJ/m/d/K]
    double c_s       =pHRU->GetSoilProps(0)->heat_capacity;   //[MJ/m3/K]
    double poro      =pHRU->GetSoilProps(0)->porosity;      
    double soil_dens =pHRU->GetSoilProps(0)->bulk_density;    //[kg/m3]
    double damp_coeff=pHRU->GetSurfaceProps()->snow_dampen_coeff; //[m]

    double sat =soil_wc/(poro*thickness*MM_PER_METER);
    
    double heat_cap  =CalculateHeatCapacity       (poro,c_s      ,sat,pct_froz);
    double therm_cond=CalculateThermalConductivity(poro,kappa_s  ,sat,pct_froz);
    double density   =CalculateDensity            (poro,soil_dens,sat,pct_froz);

    double C_ice=0.0; 
    const double RANGE=3.0;
    if ((soil_temp>-RANGE) && (soil_temp<0.1)){C_ice=10;}//MJ/m3/K - gross approximation

    //better?
    if ((soil_temp>-RANGE) && (soil_temp<0.1)){
    //  C_ice=DENSITY_ICE*LH_FUSION*sat*poro/(RANGE+0.1); //3.1 is temp range, in  degC - allows quick freezing of dry soils
    }
    
    double snow_factor=exp(-damp_coeff*snow_depth);

    double soil_temp_new=(soil_temp+therm_cond/thickness/thickness/(soil_dens*heat_cap+C_ice)*(Tair-soil_temp)*Options.timestep)*snow_factor;

    //[MJ/m/d/K]/[m^2]/[MJ/m3/K]*[C]*[d]=[degC]

    double pct_froz_new=min(max((soil_temp_new+RANGE)/(RANGE+0.1),0.0),1.0); //linear variation from -3 to 0.1 degrees

    rates[0]=(soil_temp_new-soil_temp)/Options.timestep;   //SOIL_TEMP->SOIL_TEMP;
    rates[1]=(pct_froz_new -pct_froz )/Options.timestep;   //SOIL_PCT_FROZ->SOIL_PCT_FROZ;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
///
/// \param *state_vars [in] state variable array for this HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of loss from soil to other soil layer [mm/day]
//
void   CmvFrozenGround::ApplyConstraints( const double     *state_vars,
                                           const CHydroUnit *pHRU,
                                           const optStruct  &Options,
                                           const time_struct &tt,
                                           double     *rates) const
{
}
