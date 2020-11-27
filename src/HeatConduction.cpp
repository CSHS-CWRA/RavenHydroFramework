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
  double kappa;//[MJ/m/d/K]
  //arithmetic mean of thermal conductivity
  /*kappa= (    sat)*(    poro)*(1.0-Fice)*TC_WATER+
           (    sat)*(    poro)*(    Fice)*TC_ICE  +
           (1.0-sat)*(    poro)           *TC_AIR+
                     (1.0-poro)           *kappa_soil; */ //[MJ/m/d/K]

 //geometric mean of thermal conductivity
 kappa= pow(kappa_soil,1.0-poro)*
        pow(TC_WATER  ,sat*poro*(1-Fice))*
        pow(TC_ICE    ,sat*poro*(  Fice))*
        pow(TC_AIR    ,(1.0-sat)*poro); //[MJ/m/d/K]
 
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
    kappa=(K)*kappa+(1.0-K)*dry_cond*WATT_TO_MJ_PER_D;//CLM eqn 6.58
  }
  */
  return kappa;
}
int CmvHeatConduction::_nHRUs=-1;
void CmvHeatConduction::StoreNumberOfHRUs(const int nHRUs){_nHRUs=nHRUs;}
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

  //Should at least confirm TEMPERATURE is transported??
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
/// \brief Inverts Tridiagonal matrix with diagonals cc,a, and b using algorithm from
/// assumes cc[0] and b[N] are zero, and all three indices correspond to matrix row - performs a shift to standard tridiagonal form
/// from Usmani,R. A. (1994). "Inversion of a tridiagonal jacobi matrix". Linear Algebra and its Applications. 212-213: 413–414.
/// \param cc   [in ] left diagonal  (c[i]=A[i][i-1]) [length=size]
/// \param a    [in ] diagonal       (a[i]=A[i][i]) [length=size]
/// \param b    [in ] right diagonal (b[i]=A[i][i+1] [length=size]
/// \param Ainv [out] resultant inverse matrix (assumes memory is already allocated in array of doubles)
/// \param size [in]  N, size of NxN square matrix 
///
void InvertTridiagonal(const double *cc,const double *a,const double *b,
                       double **Ainv,const int     size)
{
  int i,j,k;
  int N=size;
  double th_im1,th_jm1,ph_ip1,ph_jp1;
  double bprod,cprod;
  double *theta=new double [N];
  double *phi  =new double [N];
  double *c    =new double [N];
  for (i=0;i<N;i++){c[i]=cc[i+1]; } //requires a shift of the off-diagonal 

  theta[0]=a[0];
  theta[1]=a[1]*theta[0]-b[0]*c[0];
  for(i=2;i<size;i++) { theta[i]=a[i]*theta[i-1]-b[i-1]*c[i-1]*theta[i-2]; }
  phi[N-1]=a[N-1];
  phi[N-2]=a[N-2]*phi[N-1]-b[N-2]*c[N-2];
  for(i=N-3;i>=0;i--) { phi  [i]=a[i]*phi  [i+1]-b[i  ]*c[i  ]*phi  [i+2]; }
  for(i=0;i<N;i++) {
    if(i==0)  { th_im1=1.0;       }
    else      { th_im1=theta[i-1];}
    if(i==N-1){ ph_ip1=1.0;       }
    else      { ph_ip1=phi  [i+1];}
    for(j=0;j<N;j++) {
      if(j==0)  { th_jm1=1.0;       }
      else      { th_jm1=theta[j-1];}
      if(j==N-1){ ph_jp1=1.0;       }
      else      { ph_jp1=phi  [j+1];}
    
      if(i<j) {
        bprod=1.0; for(k=i;k<=j-1;k++) { bprod*=b[k]; }
        Ainv[i][j]=pow(-1.0,i+j+2)*bprod*th_im1*ph_jp1/theta[N-1];
      }
      else if(i==j) {
        Ainv[i][j]=th_im1*ph_jp1/theta[N-1];
      }
      else {
        cprod=1.0; for(k=j;k<=i-1;k++) { cprod*=c[k]; }
        Ainv[i][j]=pow(-1.0,i+j+2)*cprod*th_jm1*ph_ip1/theta[N-1];
      }
    }
  }
  delete [] theta;
  delete [] phi;
  delete [] c;
}
void MatVecMult(const double * const *A,const double *x,double *B,int N) {
  int i,j;
  for(i=0;i<N;i++) {
    B[i]=0.0;
    for(j=0;j<N;j++){
      B[i]+=A[j][i]*x[j];
    }
  }
}
void TestInversion() {
  int N=10;
  int i,j;
  double *a=new double [N];
  double *b=new double [N];
  double *c=new double [N];
  double **Ainv=new double *[N];
  double **A   =new double *[N];
  double *x=new double [N];
  double *xx=new double [N];
  double *B=new double [N];
  for(i=0;i<N;i++) {
    Ainv[i]=new double [N];
    A   [i]=new double [N];
    for(j=0;j<N;j++) { Ainv[i][j]=A[i][j]=0.0; }
    a[i]=2;
    b[i]=-1;
    c[i]=-1.3*x[i];
    B[i]=0;
    x[i]=(double)(i)/double(N); //arbitrary
    c[i]=-1.3*x[i];
    A[i][i]=a[i];
    if(i>0  ) { A[i][i-1]=c[i]; }
    if(i<N-1) { A[i][i+1]=b[i]; }
  }
  InvertTridiagonal(c,a,b,Ainv,N);

  MatVecMult(A,x,B,N);

  MatVecMult(Ainv,B,xx,N);

  for(i=0;i<N;i++) {
    for(j=0;j<N;j++) {
      cout<<A[i][j]<<" ";
    }
    cout<<" * "<<x[i];
    cout<<" = "<<B[i]<<endl;
  }
  cout<<"-------------------------------"<<endl;
  for(i=0;i<N;i++) {
    for(j=0;j<N;j++) {
      cout<<Ainv[i][j]<<" ";
    }
    cout<<" * "<<B[i];
    cout<<" = "<<xx[i]<<endl;
  }
  double err(0.0);
  for(i=0;i<N;i++) { err+=xx[i]-x[i]; cout<<"compare: "<<x[i]<< " "<<xx[i]<<endl;}
  cout<<"Inversion error: "<<err<<endl;

  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;
  delete [] B;
  for(i=0;i<N;i++) { delete[] Ainv[i]; } delete [] Ainv;
  for(i=0;i<N;i++) { delete[] A[i]; } delete[] A;
  ExitGracefully("Unit Testing Done",BAD_DATA);
}


//////////////////////////////////////////////////////////////////
/// \brief Generates Jacobian for Newton solution to non-linear heat conduction problem
/// \param *z [in] centroid of compartments
/// \param *poro
/// \param *H [in] vector of current energy content of each compartment
/// \param tstep [in] time step (could be local timestep)
/// \param J  [out] 3xN matrix storing diagonals of tridiagonal Jacobian [presumed memory is pre-allocated]
//
void CmvHeatConduction::GenerateJacobianMatrix(const double  *z,
                                               const double  *eta, const double * poro,const double *kappa_s,
                                               const double  *sat,  const double *satn,
                                               const double  *Vold, const double *Vnew, //[mm]
                                               const double  *hold, const double *h,
                                               const double  &tstep,
                                                     double **J, 
                                                     double  *f,
                                               const int      N) const 
{
  //    eta=thick[i]*(1-poro[i])*soil_hcp[i];
  double dTdHn;
  double kappal,kappar,kappaln,kapparn,kappaln_d,kapparn_d,kappaln_di,kapparn_di;
  double ai,bi;
  double Bn;
  double dh=0.01;
  double Fice,Ficen,Ficen_d;
  double *kap    =new double [N];
  double *kapn   =new double [N];
  double *kapn_d =new double [N];
  double *T      =new double [N];
  double *Tn     =new double [N];
  double *kapl   =new double [N];
  double *kapr   =new double [N];
  double *kapl_d =new double [N];
  double *kapr_d =new double [N]; // \todo [optimize] - static storage
  for(int i=0;i<N;i++)
  {
    Fice      =ConvertVolumetricEnthalpyToIceContent(hold[i]);
    Ficen     =ConvertVolumetricEnthalpyToIceContent(h   [i]);
    Ficen_d   =ConvertVolumetricEnthalpyToIceContent(h[i]+dh);

    kap    [i]=CalculateThermalConductivity(poro[i],kappa_s[i], sat[i],Fice);
    kapn   [i]=CalculateThermalConductivity(poro[i],kappa_s[i],satn[i],Ficen);
    kapn_d [i]=CalculateThermalConductivity(poro[i],kappa_s[i],satn[i],Ficen_d);

    T      [i]=ConvertVolumetricEnthalpyToTemperature(hold[i]);
    Tn     [i]=ConvertVolumetricEnthalpyToTemperature(   h[i]);
  }
  for(int i=0;i<N;i++) 
  {
    if(i==0) {
      ai=0.5*tstep/(z[i+1]-z[i]);
      bi=0;
    }
    else if(i==N-1) {
      ai=0;
      bi=0.5*tstep/(z[i]-z[i-1]);
    }
    else {
      ai=0.5*tstep/(z[i+1]-z[i]);
      bi=0.5*tstep/(z[i]-z[i-1]);
    }
    if(i!=0) {
      kappal    =pow(kap   [i]*kap   [i-1],0.5); // geometric mean
      kappaln   =pow(kapn  [i]*kapn  [i-1],0.5);
      kappaln_d =pow(kapn  [i]*kapn_d[i-1],0.5);
      kappaln_di=pow(kapn_d[i]*kapn  [i-1],0.5);
    }
    else {kappal=kappaln=kappaln_d=kappaln_di=0.0;}
    if(i!=N-1) {
      kappar    =pow(kap   [i]*kap   [i+1],0.5); // geometric mean
      kapparn   =pow(kapn  [i]*kapn  [i+1],0.5);
      kapparn_d =pow(kapn  [i]*kapn_d[i+1],0.5);
      kapparn_di=pow(kapn_d[i]*kapn  [i+1],0.5);
    }
    else { kappar=kapparn=kapparn_d=kapparn_di=0.0; }

    dTdHn=TemperatureEnthalpyDerivative(h[i]);

    Bn=hold[i]*Vold[i]-ai*kappar*(T[i]-T[i+1])-bi*kappal*(T[i]-T[i-1])+eta[i]*T[i];

    f[i]=h[i]*Vnew[i]+ai*kapparn*(Tn[i]-Tn[i+1])+bi*kappaln*(Tn[i]-Tn[i-1])+eta[i]*Tn[i]-Bn;

    //left diagonal
    if(i!=0) {
      J[0][i]=bi*(Tn[i]-Tn[i-1])*(kappaln_d-kappaln)/dh-bi*kappaln*TemperatureEnthalpyDerivative(h[i-1]);
    }
    else{J[0][i]=0.0; }

    //central diagonal 
    J[1][i] =Vnew[i];
    J[1][i]+=ai*(Tn[i]-Tn[i+1])*(kapparn_di-kapparn)/dh+ai*kapparn*dTdHn;
    J[1][i]+=bi*(Tn[i]-Tn[i-1])*(kappaln_di-kappaln)/dh+bi*kappaln*dTdHn;
    J[1][i]+=eta[i]*dTdHn;

    //right diagonal
    if(i!=N-1) {
      J[2][i]=ai*(Tn[i]-Tn[i+1])*(kapparn_d-kapparn)/dh-ai*kapparn*TemperatureEnthalpyDerivative(h[i+1]);
    }
    else {J[2][i]=0.0; }
  }
  delete [] kap;
  delete [] kapn;
  delete [] kapn_d;
  delete [] T;
  delete [] Tn;
  delete [] kapl;
  delete [] kapr;
  delete [] kapr_d;
  delete [] kapl_d;
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
  int    iSoil,iSoilEnthalpy;
  double hcond    = pHRU->GetSurfaceProps()->convection_coeff; //[MJ/m2/d/K] 
  double geo_grad = pHRU->GetSurfaceProps()->geothermal_grad; //[MJ/m2]
  int    nSoils=pModel->GetNumSoilLayers();
  
  int N=nSoils;//+nSnowLayers
  
  static double **v_old;
  static double *dz,*z;
  static double *poro,*sat,*satn;
  static double *eta,*Vold,*Vnew;
  static double *Tnew,*Told;
  static double *kappa_s,*kap,*kapn;//[MJ/m/d/K]
  static double *hold,*hguess,*delta_h,*f;
  static double **J,**Jinv;
  
  // Allocate Memory for static arrays
  if(tt.model_time==0.0) 
  {
    v_old=new double *[_nHRUs];
    for(int k=0;k<_nHRUs;k++) {
      v_old[k]=new double [N];
      for(int i=0; i<N;i++) {
        v_old[k][i]=0.0;
      }
    }
    dz     =new double[N];    z      =new double[N];
    sat    =new double[N];    satn   =new double[N];    poro   =new double[N];
    kappa_s=new double[N];    eta    =new double[N];
    Vold   =new double[N];    Vnew   =new double[N];
    Tnew   =new double[N];    Told   =new double[N];
    kap    =new double[N];    kapn   =new double[N];
    hguess =new double[N];    delta_h=new double[N];    hold   =new double[N];
    f      =new double[N];
    J      =new double *[3];
    for(int i=0;i<3;i++) {    J[i]=new double[N]; }
    Jinv  =new double *[N];
    for(int i=0;i<N;i++) { Jinv[i]=new double[N]; }
  }

  int k=pHRU->GetGlobalIndex();

  // Get Soil properties 
  for(int m=0;m<nSoils;m++)
  {
    iSoil=pModel->GetStateVarIndex(SOIL,m);
    iSoilEnthalpy=iTo[m];

    poro   [m]=pHRU->GetSoilProps(m)->porosity;
    kappa_s[m]=pHRU->GetSoilProps(m)->thermal_cond;
    dz     [m]=pHRU->GetSoilThickness(m)/MM_PER_METER;//dz in meters
    eta    [m]=dz[m]*(1.0-poro[m])*pHRU->GetSoilProps(m)->heat_capacity;

    if(m==0) { z[m]=0.5*dz[m]; }
    else     { z[m]=0.5*(dz[m]+dz[m-1])+z[m-1]; }

    if(tt.model_time==0.0) {
      v_old[k][m]=state_vars[iSoil]/MM_PER_METER; 
    }

    Vold[m]=v_old[k][m]; //[m]
    Vnew[m]=state_vars[iSoil]/MM_PER_METER;
    sat [m]=Vnew[m]*MM_PER_METER/pHRU->GetSoilCapacity(m);
    satn[m]=Vold[m]*MM_PER_METER/pHRU->GetSoilCapacity(m);

    //update vold for next time step
    v_old[k][m]=Vnew[m];

    hold[m]=0.0;
    if(Vold[m]>1e-6) { hold[m]=state_vars[iSoilEnthalpy]/Vold[m]; }//[MJ/m3] 
    hguess[m]=hold[m]; 
  }

  // Crank Nicolson approach using Newton Raphson
  const double TOLERANCE=1e-3; //~0.2e-3 degrees C
  const int    MAX_ITER=30;
  double       err=ALMOST_INF;
  int          i,iter;

  //double hguess_n=hguess;
  //for (n=0;n<nDivs;n++){
  iter=1;
  while((err>TOLERANCE) && (iter<MAX_ITER))
  {
    GenerateJacobianMatrix(z,eta,poro,kappa_s,sat,satn,Vold,Vnew,hold,hguess,tstep,J,f,N);

    /*double satint =sat +ddt*(n  )/dt*(satn-sat );
    double satintn=sat +ddt*(n+1)/dt*(satn-sat );
    double Vint   =Vold+ddt*(n  )/dt*(Vnew-Vold);
    double Vintn  =Vold+ddt*(n+1)/dt*(Vnew-Vold);
    //GenerateJacobianMatrix(z,eta,poro,kappa_s,satint,satintn,Vint,Vintn,hguess_l,hguess,tstep/nDivs,J,f,N);*/
    
    InvertTridiagonal(J[0],J[1],J[2],Jinv,N);

    MatVecMult(Jinv,f,delta_h,N);

    /* TESTING JACOBIAN INVERSION -UNIFORMLY WORKS  
    double *ff=new double [N];
    double ** JJ=new double *[N];
    for(int i=0;i<N;i++) { JJ[i]=new double[N];for(int j=0;j<N;j++) { JJ[i][j]=0.0; } }
    for(int i=0;i<N;i++) { JJ[i][i]=J[1][i]; if(i>0) { JJ[i][i-1]=J[0][i]; } if(i<N-1) { JJ[i][i+1]=J[2][i]; }
    }
    MatVecMult(JJ,delta_h,ff,N); 
    double err(0.0);
    for(i=0;i<N;i++) { err+=ff[i]-f[i];}
    if(fabs(err)>0.0001) {
      cout<<"Inversion error: "<<err<<endl;
    }
    delete [] ff;
    for(int i=0;i<N;i++) { delete[] JJ[i]; }delete [] JJ;*/

    err=0.0;
    for(i=0;i<N;i++) {
      hguess[i]=hguess[i]-delta_h[i];
//      if(iter>15) {cout<<std::setprecision(5)<<hguess[i]<<"|";}
      upperswap(err,fabs(delta_h[i]));
    }
//    if(iter>15) { cout<<endl; }
    //cout<<endl<<"iteration "<<iter<<" err: "<<err<<endl;
    iter++;
  }
  //if (iter==MAX_ITER){cout<<"HIT MAX ITERATIONS "<<err<<endl;}
  
  //Post-processing: extract energy fluxes from solution
  //Need also to address heat exchanged with soil (sink/source term) 
  double Ficen,Ficeo;
  for(int m=0;m<N;m++)
  {
    Tnew[m]=ConvertVolumetricEnthalpyToTemperature(hguess[m]);
    Told[m]=ConvertVolumetricEnthalpyToTemperature(hold[m]);
    Ficen  =ConvertVolumetricEnthalpyToIceContent (hguess[m]);
    Ficeo  =ConvertVolumetricEnthalpyToIceContent (hold[m]);
    kap [m]=CalculateThermalConductivity(poro[m],kappa_s[m],sat [m],Ficeo);
    kapn[m]=CalculateThermalConductivity(poro[m],kappa_s[m],satn[m],Ficen);
  }
  
  /*
  //cout<<"HEATCOND: ";
  for(int m=0;m<N;m++)
  {
    //cout <<"KAPPA: "<<kap[m]<<" "<<poro[m]<<" "<<sat[m]<<" "<<Ficeo<< " "<<kappa_s[m]<<endl;
    //cout <<Tnew[m]<<" ";
    //cout<<kap[m]<<" ";
    //cout<<ConvertVolumetricEnthalpyToIceContent(hguess[m])<<" ";
    //cout<<eta[m]<<" ";
  }
  //ExitGracefully("stub",STUB);
  //cout<<endl;
  */

  double Tair=pHRU->GetForcingFunctions()->temp_ave;
  double kappal,kappaln;
  for(int m=0;m<N;m++)
  {
    if(m!=0) {
      kappal    =pow(kap [m]*kap [m-1],0.5); // geometric mean
      kappaln   =pow(kapn[m]*kapn[m-1],0.5);
    }
    else { kappal=kappaln=0.0; }
    if(m>0)        { 
      rates[m] =-0.5*kappal  * 1.0/(0.5*(dz[m]+dz[m-1]))*(Told[m]-Told[m-1]);// [MJ/m2/d]
      rates[m]+=-0.5*kappaln * 1.0/(0.5*(dz[m]+dz[m-1]))*(Tnew[m]-Tnew[m-1]);
    } 
    else { rates[m]=hcond*(Tair-0.5*(Told[m]+Tnew[m])); } //top condition [MJ/m2/d]
  }
  for(int m=0;m<N;m++)
  {
    //-dHs/dt = - eta*dT/dt  = gain of heat by soil, lost by water if Tnew>Told
    rates[m+N] = -eta[m]*(Tnew[m]-Told[m])/Options.timestep; 
  }
  rates[N]=kap[N]*geo_grad; //geothermal [MJ/m2/d]

  // Delete static arrays
  if(tt.model_time>=Options.duration-Options.timestep-TIME_CORRECTION)
  {
    for(int k=0;k<_nHRUs;k++) {delete [] v_old[k];} delete [] v_old;
    delete[] dz;  delete[] z;    delete[] poro;
    delete[] sat; delete[] satn; delete[] kappa_s;
    delete[] eta; delete[] Vold; delete[] Vnew;
    delete[] hold;
    delete[] Tnew;   delete[] Told;    delete[] kap; delete[] kapn;
    delete[] hguess; delete[] delta_h; delete[] f;
    for(int i=0;i<3;i++) { delete[] J[i]; } delete[] J;
    for(int i=0;i<N;i++) { delete[]  Jinv[i]; } delete[] Jinv;
  }
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

