/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2020 the Raven Development Team
------------------------------------------------------------------
Soil Heat Conduction
convention assumes conduction is positive downward
----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "HeatConduction.h"

/////////////////////////////////////////////////////////////////////
/// \brief calculates volumetric heat capacity[MJ/m3/K] of soil/water/ice mix
/// \details This is just a simple weighted average of soil, ice and water capacities
///
/// \param &sat_liq [in] Saturation fraction of liquid H_{2}0 [0..1]
/// \param &sat_ice [in] Saturated fraction of ice [0..1]
/// \param *pS [in] Soil properties structure
/// \return Volumetric soil heat capacity [J/m^3/K]
//
//
double CalculateHeatCapacity(const double &poro,const double&hcp_soil,const double &sat,const double &Fice) 
{
   return  (    sat)*(    poro)*(1.0-Fice)*HCP_WATER+
           (    sat)*(    poro)*(    Fice)*HCP_ICE  +
           (1.0-sat)*(    poro)           *HCP_AIR+
                     (1.0-poro)           *hcp_soil; //[MJ/m3/K]
}

//////////////////////////////////////////////////////////////////
/// \brief calculates density [kg/m3] of soil/water/ice mix
//
double CalculateDensity(const double &poro,const double&rho_soil,const double &sat,const double &Fice) 
{
   return  (    sat)*(    poro)*(1.0-Fice)*DENSITY_WATER+
           (    sat)*(    poro)*(    Fice)*DENSITY_ICE  +
           (1.0-sat)*(    poro)           *DENSITY_AIR+
                     (1.0-poro)           *rho_soil; //[MJ/m3/K]
}

//////////////////////////////////////////////////////////////////
/// \brief calculates thermal conductivity [MJ/m2/d/K] of soil/water/ice mix
//
double CalculateThermalConductivity(const double &poro,const double&kappa_soil,const double &sat,const double &Fice)
{
  double kappa;
  //arithmetic mean of thermal conductivity
  /*kappa= (    sat)*(    poro)*(1.0-Fice)*TC_WATER+
         (    sat)*(    poro)*(    Fice)*TC_ICE  +
         (1.0-sat)*(    poro)           *TC_AIR+
                   (1.0-poro)           *kappa_soil; */ //[W/m2/K]

 //geometric mean of thermal conductivity
 kappa= pow(kappa_soil,1.0-poro)*
        pow(TC_WATER  ,sat*poro*(1-Fice))*
        pow(TC_ICE    ,sat*poro*(  Fice))*
        pow(TC_AIR    ,(1.0-sat)*poro);
 
 //Kersten number
 /// \ref adapted from Community Land Model Manual 3.0, Oleson et al., 2004 \cite Oleson2012
  /*double K=sat;
  if(Fice==0.0) {
    K = max(log(sat)+1.0,0.0);
  }
  //final thermal conductivity
  //dry soil conductivity [W/m/K]
  dry_cond=(0.135*rhob+64.7)/(2700-0.947*rhob); //CLM manual eqn. 6.62
  kappa=kappa_soil;
  if(sat>REAL_SMALL) {
    kappa=(K)*kappa+(1.0-K)*dry_cond;//CLM eqn 6.58
  }
  */
  return kappa*WATT_TO_MJ_PER_D; //[W/m2/K]->[MJ/m2/d/K]
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the heat conduction constructor
///
/// \param pTransMod [in] pointer to transport model
//
CmvHeatConduction::CmvHeatConduction(const CTransportModel *pTransMod) : CHydroProcessABC(HEATCONDUCTION)
{
  _pTransModel=pTransMod;

  int nSoils=pModel->GetNumSoilLayers();
  
  CHydroProcessABC::DynamicSpecifyConnections(nSoils+1);

  int cTemp=_pTransModel->GetConstituentIndex("TEMPERATURE");
  ExitGracefullyIf(cTemp==DOESNT_EXIST,
    "CmvHeatConduction Constructor - TEMPERATURE transport must be turned on for soil heat conduction algorithm to be used",BAD_DATA_WARN);

  int layer1, layer2;

  //top connection - atmosphere to topsoil 
  //-------------------------------------------------------------
  layer2=_pTransModel->GetLayerIndex(cTemp,pModel->GetStateVarIndex(SOIL,0)); //layer index of soil[0] enthalpy 

  iFrom[0]=pModel->GetStateVarIndex(CONSTITUENT_SRC,cTemp); //global index of energy source state variable 
  iTo  [0]=pModel->GetStateVarIndex(CONSTITUENT,    layer2);//global index of enthalpy in soil[0]

  //intermediate connections - soil to soil 
  //-------------------------------------------------------------
  for(int m=0;m<nSoils-1;m++) 
  {
    layer1=_pTransModel->GetLayerIndex(cTemp,pModel->GetStateVarIndex(SOIL,m  )); //layer index of soil[m] enthalpy 
    layer2=_pTransModel->GetLayerIndex(cTemp,pModel->GetStateVarIndex(SOIL,m+1)); //layer index of soil[m+1] enthalpy 

    iFrom[m+1]=pModel->GetStateVarIndex(CONSTITUENT,layer1);//global index of enthalpy in soil[m] 
    iTo  [m+1]=pModel->GetStateVarIndex(CONSTITUENT,layer2);//global index of enthalpy in soil[m+1]
  }
  //bottom connection -geothermal 
  //-------------------------------------------------------------
  layer1=_pTransModel->GetLayerIndex(cTemp,pModel->GetStateVarIndex(SOIL,nSoils-1)); //layer index of soil[m] enthalpy 

  iFrom[nSoils]=pModel->GetStateVarIndex(CONSTITUENT,layer1);     //global index of enthalpy in soil[nSoils-1] 
  iTo  [nSoils]=pModel->GetStateVarIndex(CONSTITUENT_SRC,cTemp);  //global index of energy source state variable 

  //JRC: If we only track water enthalpy, we may need to include the sink/source term to the soil as it re-equilibrates
  /*for(int m=0;m<nSoils;m++)
  {
    layer1=_pTransModel->GetLayerIndex(cTemp,pModel->GetStateVarIndex(SOIL,m)); //layer index of soil[m] enthalpy 

    iFrom[nSoils+m+1]=pModel->GetStateVarIndex(CONSTITUENT,   layer1);//global index of enthalpy in soil[m] 
    iTo  [nSoils+m+1]=pModel->GetStateVarIndex(CONSTITUENT_SRC,cTemp);//global index of energy source state variable 
  }*/
}
//////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
//
CmvHeatConduction::~CmvHeatConduction()
{
  //no actions required
}

//inherited functions
//////////////////////////////////////////////////////////////////
/// \brief Initializes heat conduction process
///
//
void CmvHeatConduction::Initialize() 
{
  //does nothing for now
}
//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for baseflow algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by baseflow algorithm (size of aP[] and aPC[])
//
void        CmvHeatConduction::GetParticipatingParamList(string  *aP,class_type *aPC,int &nP) const
{
  nP=0;
  aP[nP]="THERMAL_COND";     aPC[nP]=CLASS_SOIL; nP++;
  aP[nP]="HEAT_CAPACITY";    aPC[nP]=CLASS_SOIL; nP++;
  aP[nP]="CONVECTION_COEFF"; aPC[nP]=CLASS_LANDUSE; nP++;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating state variable list
///
/// \param btype [in] Baseflow algorithm type
/// \param *aSV [out] Array of state variable types needed by baseflow algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by baseflow algorithm (size of aSV[] and aLev[] arrays)
//
void CmvHeatConduction::GetParticipatingStateVarList(sv_type *aSV,int *aLev,int &nSV)
{
  //NOT NEEDED. Assumes that :Transport command takes care of this
  //actual state variables just enthalpy of soil compartments
  nSV=0;
}

//////////////////////////////////////////////////////////////////
/// \brief Solves Tridiagonal matrix with diagonals e,f, and g using Thomas Algorithm. RHS=b, returns solution vector x
///
void ThomasAlgorithm(Ironclad1DArray  e,Ironclad1DArray f,
                     Ironclad1DArray  g,Ironclad1DArray b,
                     Writeable1DArray x,const int     size)
{
  int i;
  double *f1=new double[size];
  double *b1=new double[size];

  //calculate f prime vector
  f1[0]=f[0];
  for(i=1; i<size; i++) {
    f1[i]=f[i]-(e[i]*g[i-1]/f1[i-1]);
  }
  //calculate b prime vector
  b1[0]=b[0];
  for(i=1; i<size; i++) {
    b1[i]=b[i]-(e[i]*b1[i-1]/f1[i-1]);
  }
  //solve for x using back substitution
  x[size-1]=b1[size-1]/f1[size-1];
  for(i=size-2; i>=0; i--) {
    x[i]=(b1[i]-g[i]*x[i+1])/f1[i];
  }
  delete[] f1;
  delete[] b1;
}

//////////////////////////////////////////////////////////////////
/// \brief Finds thermal rate of change
/// \param *state_vars [in] Array of state variable values for this HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rates of loss from each soil unit to the next down [MJ/m2d]
///
//
void CmvHeatConduction::GetRatesOfChange(const double      *state_vars,
                                         const CHydroUnit  *pHRU,
                                         const optStruct   &Options,
                                         const time_struct &tt,
                                               double      *rates) const
{
  double tstep=Options.timestep;
  double hv; //enthalpy [MJ/m3]
  double poro, sat;
  int    iSoil,iSoilEnthalpy;
  double hcond = pHRU->GetSurfaceProps()->convection_coeff; //[MJ/m2/d/K] ;

  int    nSoils=pModel->GetNumSoilLayers();
  
  // Get Soil properties 
  double *k   =new double [nSoils]; //TODO: static storage 
  double *cp  =new double [nSoils];
  double *rho =new double [nSoils];
  double *T   =new double [nSoils];
  double *Tnew=new double [nSoils];
  double *dz = new double [nSoils];
  double *Fi = new double [nSoils];
  for(int m=0;m<nSoils;m++)
  {
    iSoil=pModel->GetStateVarIndex(SOIL,m);
    iSoilEnthalpy=iTo[m];
    poro  =pHRU->GetSoilProps(m)->porosity;
    sat   =state_vars[iSoil]/pHRU->GetSoilCapacity(m);

    hv=0.0;
    if (state_vars[iSoil]>1e-9){hv=state_vars[iSoilEnthalpy]/(state_vars[iSoil]/MM_PER_METER);}//[MJ/m3] 

    T  [m]=ConvertVolumetricEnthalpyToTemperature(hv);
    Fi [m]=ConvertVolumetricEnthalpyToIceContent (hv);

    cp [m]=CalculateHeatCapacity       (poro,pHRU->GetSoilProps(m)->heat_capacity,sat,Fi[m]);
    k  [m]=CalculateThermalConductivity(poro,pHRU->GetSoilProps(m)->thermal_cond ,sat,Fi[m]);
    rho[m]=CalculateDensity            (poro,pHRU->GetSoilProps(m)->bulk_density ,sat,Fi[m]);
    
    //cp [m]+= pHRU->GetSoilProps(m)->heat_capacity*pHRU->GetSoilProps(m)->bulk_density /sat / poro / DENSITY_WATER; //JRC: this can be unstable for small saturations. 

    dz [m]=pHRU->GetSoilThickness(m);
  }

  // Crank Nicolson approach using Thomas Algorithm
  double *AA =new double[nSoils]; //TODO: Static storage
  double *BB =new double[nSoils];
  double *CC =new double[nSoils];
  double *RHS=new double[nSoils];
  double *x  =new double[nSoils];
  double a,b,c,d;

  double Tair =pHRU->GetForcingFunctions()->temp_ave;

  //top condition - heat flux = h*(T-Tair)
  //bottom condition - flux = 0
  for(int m=0;m<nSoils;m++) 
  {
    if(m>0)       { a=0.5*(k[m]+k[m-1])/2.0/(0.5*(dz[m]+dz[m-1]))/dz[m]; } //arithmetic mean (for now)
    else          { a=hcond/2.0/dz[0];}
    if(m<nSoils-1){ c=0.5*(k[m]+k[m+1])/2.0/(0.5*(dz[m]+dz[m+1]))/dz[m]; }
    else          { c=0.0; }
    b=-(a+c);
    d=cp[m]*rho[m]/tstep;
    AA [m]=-a;
    BB [m]=(b+d);
    CC [m]=-c;
    if      (m==0       ) { RHS[m]=a*2*Tair+(d-b)*T[m]+c*T[m+1];}
    else if (m==nSoils-1) { RHS[m]=a*T[m-1]+(d-b)*T[m]         ;}
    else                  { RHS[m]=a*T[m-1]+(d-b)*T[m]+c*T[m+1];}
  }

  //solve system of equations
  ThomasAlgorithm(AA,BB,CC,RHS,Tnew,nSoils);
  
  //extract energy fluxes from solution
  //Need also to address heat exchanged with soil (sink/source term) 
  for(int m=0;m<nSoils;m++)
  {
    //-dHs/dt = - H*rho_b*c_p*dT/dt
    //rates[m+nSoils] = - dz[m]*  pHRU->GetSoilProps(m)->heat_capacity*pHRU->GetSoilProps(m)->bulk_density*(Tnew[m]-T[m])/Options.timestep; //JRC
  }
  for(int m=0;m<nSoils;m++)
  {
    T[m]=0.5*(T[m]+Tnew[m]); //T is now mean temp over time step
    if(m>0)        { rates[m]=(k[m]+k[m-1])/2.0 * 1.0/(0.5*(dz[m]+dz[m-1]))*(T[m]-T[m-1]); } // [MJ/m2/d]
    else           { rates[m]=hcond*(Tair-T[m]); } //top condition [MJ/m2/d]
  }
  rates[nSoils]=0.0; //geothermal (eventually) [MJ/m2/d]

  delete [] AA; delete [] BB; delete [] CC;  delete [] RHS; 
  delete [] k;  delete [] cp; delete [] rho; delete [] Fi;
  delete [] T;  delete [] Tnew;
}
//////////////////////////////////////////////////////////////////
/// \brief applies constraints to state variables 
/// \param *state_vars [in] Array of state variable values for this HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rates of loss from each soil unit to the next down [MJ/m2d]
///
//
void CmvHeatConduction::ApplyConstraints(const double      *state_vars,
                                         const CHydroUnit  *pHRU,
                                         const optStruct   &Options,
                                         const time_struct &tt,
                                               double      *rates) const
{
  //DOES NOTHING - constraints handled within GetRatesOfChange
}

