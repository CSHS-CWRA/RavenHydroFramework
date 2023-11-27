/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2023 the Raven Development Team
------------------------------------------------------------------
Soil Heat Conduction
convention assumes conduction is positive downward
----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "HeatConduction.h"
double MeanConvectiveFlux(const double &hn,const double &Vw,const double &alpha,const double &Ta,const double &dt);

/////////////////////////////////////////////////////////////////////
/// \brief calculates volumetric heat capacity[MJ/m3/K] of soil/water/ice mix
/// \details This is just a simple weighted average of soil, ice and water capacities
///
/// \param &poro [in] soil porosity [0..1]
/// \param &sat [in] Saturated fraction  [0..1]
/// \param &Fice [in] Ice fraction  [0..1]
/// \param &hcp_soil [in] heat capacity of soil [MJ/m3/K]
/// \return Volumetric soil heat capacity [MJ/m3/K]
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

 //fixed value (as if liquid water saturated)
 //kappa= pow(kappa_soil,1.0-poro)*pow(TC_WATER,poro);
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
void CmvHeatConduction::StoreNumberOfHRUs(const int nHRUs){
  _nHRUs=nHRUs;
}
//////////////////////////////////////////////////////////////////
/// \brief Implementation of the heat conduction constructor
///
/// \param pTransMod [in] pointer to transport model
//
CmvHeatConduction::CmvHeatConduction(const CTransportModel *pTransMod,
                                     CModelABC             *pModel)
  :CHydroProcessABC(HEATCONDUCTION, pModel)
{
  //g_disable_freezing=true; TMP DEBUG ONLY

  _pTransModel=pTransMod;

  _nHRUs=-1;

  int nSoils=pModel->GetNumSoilLayers();

  CHydroProcessABC::DynamicSpecifyConnections(2*nSoils+2);

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
  for(int m=0;m<nSoils;m++)
  {
    layer1=_pTransModel->GetLayerIndex(cTemp,pModel->GetStateVarIndex(SOIL,m)); //layer index of soil[m] enthalpy

    iFrom[nSoils+1+m]=pModel->GetStateVarIndex(CONSTITUENT,layer1);//global index of enthalpy in water of soil[m]
    iTo  [nSoils+1+m]=pModel->GetStateVarIndex(SOIL_TEMP,m);//actually soil enthalpy[m]
  }

  //top connection - atmosphere to topsoil
  //-------------------------------------------------------------
  iFrom[2*nSoils+1]=pModel->GetStateVarIndex(CONSTITUENT_SRC,cTemp); //global index of energy source state variable
  iTo  [2*nSoils+1]=pModel->GetStateVarIndex(SOIL_TEMP,0);


  //JRC: If we only track water enthalpy, we may need to include the sink/source term to the soil as it re-equilibrates
  for(int m=0;m<nSoils;m++)
  {
    //layer1=_pTransModel->GetLayerIndex(cTemp,pModel->GetStateVarIndex(SOIL,m)); //layer index of soil[m] enthalpy

    //iFrom[2*nSoils+m+1]=pModel->GetStateVarIndex(CONSTITUENT,   layer1);//global index of enthalpy in soil[m]
    //iTo  [2*nSoils+m+1]=pModel->GetStateVarIndex(CONSTITUENT_SRC,cTemp);//global index of energy source state variable
  }
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
  aP[nP]="GEOTHERMAL_GRAD";  aPC[nP]=CLASS_LANDUSE; nP++;
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
  //Transport command takes care of water enthalpy storage
  int nSoils=pModel->GetNumSoilLayers();
  nSV=0;
  for(int m=0;m<nSoils;m++) {
    aSV[m]=SOIL_TEMP; aLev[m]=m; nSV++;
  }
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
/// assumes b[N] is zero, and all three indices correspond to matrix row - performs a shift to standard tridiagonal form
/// from Usmani,R. A. (1994). "Inversion of a tridiagonal jacobi matrix". Linear Algebra and its Applications. 212-213: 413-414.
/// \param cc   [in ] shifted left diagonal  (cc[i+1]=A[i][i-1]) [length=size]
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
/// \param *z [in] centroid of compartments [m] [size:N]
/// \param eta [in] array of  heat capacity term in each layer (thickness*(1-poro)*soil_hcp) [MJ/m2/K][size:N]
/// \param poro [in] array of porosities in each layer [size:N]
/// \param kappa_s [in] array of soil thermal conductivities in each layer [size:N]
/// \param satn,sat [in] arrays of new and previous saturations of each compartment [size:N]
/// \param Vold,Vnew [in] [mm] arrays of new and previous water volumes of each compartment [size:N]
/// \param hold, h [in] [MJ/m3] arrays of current and previous energy content of each compartment [size:N]
/// \param tstep [in] time step (could be local timestep)
/// \param J  [out] 3xN matrix storing diagonals of tridiagonal Jacobian [presumed memory is pre-allocated]
/// \param f [out] array of functions f terms [size:N]
/// \returns boolean indicating true if matrix is non-NULL and therefore invertible
//
bool CmvHeatConduction::GenerateJacobianMatrix( const double  *z,
                                                const double  *eta, const double * poro,const double *kappa_s,
                                                const double  *sat,  const double *satn,
                                                const double  *Vold, const double *Vnew, //[mm]
                                                const double  *hold, const double *h,
                                                const double  &tstep,
                                                      double **J,
                                                      double  *f,
                                                const int      N) const
{
  double dTdHn;
  double kappal,kappar,kappaln,kapparn,kappaln_d,kapparn_d,kappaln_di,kapparn_di;
  double ai,bi,Bn;
  double Fice,Ficen,Ficen_d;
  double dh=0.01;
  double sum=ALMOST_INF;
  bool zerorow=false;

  double *kap    =new double [N];
  double *kapn   =new double [N];
  double *kapn_d =new double [N];
  double *T      =new double [N];
  double *Tn     =new double [N]; // \todo [optimize] - static storage

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
      kappal    =pow(kap   [i]*kap   [i-1],0.5); // geometric mean (should really be pow(kap[i]^(dzi/(dzi+dz(i-1))*kap[i-1]^(dz(i-1)/(dzi+dz(i-1))
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

    Bn  =hold[i]*Vold[i]-ai*kappar *(T [i]-T [i+1])-bi*kappal *(T [i]-T [i-1])+eta[i]*T [i];

    f[i]=h   [i]*Vnew[i]+ai*kapparn*(Tn[i]-Tn[i+1])+bi*kappaln*(Tn[i]-Tn[i-1])+eta[i]*Tn[i]-Bn;

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
    /*if(g_debug_vars[4]==1) {
      cout<<"J["<<i<<"]: "<<J[0][i]<<" "<<J[1][i]<<" "<<J[2][i]<<endl;
    }*/
    if ((J[0][i]+J[1][i]+J[2][i])<PRETTY_SMALL){zerorow=true;}
    sum+=J[0][i]+J[1][i]+J[2][i];
  }
  /*if(g_debug_vars[4]==1) {
    for(int i=0;i<N;i++) {
      cout<<" eta: "<<eta[i]<<" "<<kapparn<<" "<<Vnew[i]<<" "<<h[i]<<" "<<hold[i]<<dTdHn<<endl;
    }
    if(zerorow){cout<<"zero row"<<endl; }
    if (sum==0){cout<<"zero sum"<<endl; }
  }*/
  delete [] kap;
  delete [] kapn;
  delete [] kapn_d;
  delete [] T;
  delete [] Tn;
  return (sum!=0.0) && (!zerorow); //if sum==0, NULL Jacobian - no temperature gradient and no volume, can't be inverted
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

  if (pHRU->GetHRUType()!=HRU_STANDARD){return;}

  double tstep=Options.timestep;
  int    iSoil,iSoilWaterEnthalpy,iSoilEnthalpy;
  double hcond    = pHRU->GetSurfaceProps()->convection_coeff; //[MJ/m2/d/K]
  double geo_grad = pHRU->GetSurfaceProps()->geothermal_grad; //[MJ/m2]
  int    nSoils=pModel->GetNumSoilLayers();

  //g_disable_freezing=true; //TMP DEBUG - this should make the problem easier by destroying nonlinearity. It doesn't?
  g_min_storage=0.0;

  int N=nSoils;//+nSnowLayers

  static double **v_old;//[m]
  static double *dz,*z;//[m]
  static double *poro,*sat,*satn;//[-]
  static double *eta;//[MJ/m2/K]
  static double *Vold,*Vnew; //[m]
  static double *Tnew,*Told;
  static double *kappa_s,*kap,*kapn;//[MJ/m/d/K]
  static double *hold,*hguess,*delta_h,*f;
  static double **J,**Jinv;

  // Allocate Memory for static arrays
  //-----------------------------------------------------------------------
  //if (pModel->FirstHRUFirstTimestep())
  if ((tt.model_time==0.0) && (pHRU->GetGlobalIndex()==0)) //TMP DEBUG - this is not valid if HRU 0 is disabled!
  {
    Jinv=NULL;
    v_old=new double *[_nHRUs];
    for(int k=0;k<_nHRUs;k++) {
      v_old[k]=new double [N]; for(int i=0; i<N;i++) {v_old[k][i]=0.0;}
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
    if(Jinv==NULL) {ExitGracefully("CmvHeatConduction::GetRatesOfChange", OUT_OF_MEMORY);}
    for(int i=0;i<N;i++) { Jinv[i]=new double[N]; }
  }

  int k=pHRU->GetGlobalIndex();

  // Get Soil properties
  //-----------------------------------------------------------------------
  for(int m=0;m<nSoils;m++)
  {
    iSoil=pModel->GetStateVarIndex(SOIL,m);
    iSoilWaterEnthalpy=iTo[m];

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
    if(Vold[m]>1e-6) { hold[m]=state_vars[iSoilWaterEnthalpy]/Vold[m]; }//[MJ/m3]
    hguess[m]=hold[m];
  }

  // Crank Nicolson approach using Newton Raphson
  //-----------------------------------------------------------------------
  //const double TOLERANCE=1e-3; //~0.2e-3 degrees C
  const double TOLERANCE=0.5e-2;
  const int    MAX_ITER=50;
  double       err=ALMOST_INF;
  int          i,iter;
  double relax;

  //double hguess_n=hguess; //SUBTIMESTEPPING
  //for (n=0;n<nDivs;n++){
  iter=0;
  relax=1.0;
  bool zfx=false;
  while((err>TOLERANCE) && (iter<MAX_ITER))
  {
    //if(pHRU->GetGlobalIndex()==0) { g_debug_vars[4]=1; cout<<"<<!:"<<iter<<hold[0]<<" "<<hold[1]<<" "<<hold[2]<<endl;ExitGracefullyIf(isnan(hold[0]),"DONE",RUNTIME_ERR);}
    //else                          { g_debug_vars[4]=0; }

    /*double satint =sat +ddt*(n  )/dt*(satn-sat ); //SUBTIMESTEPPING
      double satintn=sat +ddt*(n+1)/dt*(satn-sat );
      double Vint   =Vold+ddt*(n  )/dt*(Vnew-Vold);
      double Vintn  =Vold+ddt*(n+1)/dt*(Vnew-Vold);
      //if (GenerateJacobianMatrix(z,eta,poro,kappa_s,satint,satintn,Vint,Vintn,hguess_l,hguess,tstep/nDivs,J,f,N))*/

    if(GenerateJacobianMatrix(z,eta,poro,kappa_s,sat,satn,Vold,Vnew,hold,hguess,tstep,J,f,N))
    {
      InvertTridiagonal(J[0],J[1],J[2],Jinv,N);

      MatVecMult(Jinv,f,delta_h,N);
    }
    else { //handles zero volume cases
      for(i=0;i<N;i++) { delta_h[i]=0.0; }
      zfx=true;
    }
    //if (iter>10){relax=0.97;}
    err=0.0;
    for(i=0;i<N;i++) {
      hguess[i]-=relax*delta_h[i];
      upperswap(err,fabs(delta_h[i]));
      if(iter>30) {
        //cout<<delta_h[i]<<" |";
//        cout<<hguess[i]<<" |";
        //cout<<ConvertVolumetricEnthalpyToTemperature(hguess[i])<<" |";
      }
    }
    if(iter>30) {
 //     cout<<"[err:"<<err<<"]"<<endl;
    }
    iter++;
  }
  if (iter==MAX_ITER){
    /*cout<<"HIT MAX ITERATIONS err: "<<err<<" V: "<<Vold[0]<<" "<<Vnew[0]<<endl;
    for(i=0;i<N;i++) {
      cout<<"  p:"<<z[i]<<" "<<eta[i]<<" "<<poro[i]<<" "<<kappa_s[i]<<" sat: "<<sat[i]<<" "<<satn[i]<<" h:"<<hold[i]<<" "<<hguess[i]<<endl;
    }
    cout<<" In HRU "<<pHRU->GetID()<<endl;
    */
    for(i=0;i<N;i++) {
      hguess[i]=hold[i];
    }
    zfx=true;//ExitGracefully("CmvHeatConduction::Hit maximum iterations. Reduce model timestep",RUNTIME_ERR);
  }

  // Post-processing: extract energy fluxes from solution
  // Need also to address heat exchanged with soil (sink/source term)
  //-----------------------------------------------------------------------
  double Ficen,Ficeo;
  for(int m=0;m<N;m++)
  {
    Told[m]=ConvertVolumetricEnthalpyToTemperature(hold[m]);
    Ficeo  =ConvertVolumetricEnthalpyToIceContent (hold[m]);
    kap[m] =CalculateThermalConductivity(poro[m],kappa_s[m],sat[m],Ficeo);

    Tnew[m]=ConvertVolumetricEnthalpyToTemperature(hguess[m]);
    if ((Vnew[m]<1e-6) && (false))
    { //special handling of tiny water volumes
      iSoilEnthalpy=pModel->GetStateVarIndex(SOIL_TEMP,m);
      eta[m]=dz[m]*(1.0-poro[m])*pHRU->GetSoilProps(m)->heat_capacity; //TMP DEBUG - only needed here if zeroed out above
      double Hsoil=state_vars[iSoilEnthalpy]+Vnew[m]*hguess[m]; //previous soil enthalpy plus this small amount of energy added to water/soil mix
      Tnew[m]=Hsoil/eta[m];
      hguess[m]=ConvertTemperatureToVolumetricEnthalpy(Tnew[m],0.0);
    }
    Ficen  =ConvertVolumetricEnthalpyToIceContent(hguess[m]);
    kapn[m]=CalculateThermalConductivity(poro[m],kappa_s[m],satn[m],Ficen);
  }
  /*if (pHRU->GetGlobalIndex()==396){
    cout<<"[k:"<<pHRU->GetGlobalIndex()<<"] ";
    cout<<"T: "; for(int m=0;m<N;m++) {cout<<Tnew[m]<<" ";}
    cout<<" ho:";for(int m=0;m<N;m++) {cout<<hold[m]<<" "; }
    cout<<"  h:";for(int m=0;m<N;m++) {cout<<hguess[m]<<" ";}
    cout<<" vn:";for(int m=0;m<N;m++) {cout<<Vnew[m]*MM_PER_METER<<" "; }
    cout<<"[iter: "<<iter<<", err:"<<err<<"]"<<endl;
  }*/

  /*
  for(int m=0;m<N;m++)
  {
    if(pHRU->GetGlobalIndex()==0) {
      //cout <<"KAPPA: "<<kap[m]<<" "<<poro[m]<<" "<<sat[m]<<" "<<Ficeo<< " "<<kappa_s[m]<<endl;
      cout <<"Temps: "<<Tnew[m]<<" "<<hguess[m]<<" "<<Told[m]<<endl;
      //cout<<ConvertVolumetricEnthalpyToIceContent(hguess[m])<<" ";
      //cout<<eta[m]<<" ";
    }
  }
  //cout<<endl;
  */

  double Tair=pHRU->GetForcingFunctions()->temp_ave;
  double kappal,kappaln;
  double Vsoil,wfrac,condflux;

  // top convection condition [MJ/m2/d]
  //-----------------------------------------------------------------------
  Vsoil   =dz[0]*(1.0-poro[0]);
  wfrac   =Vnew[0]/(Vsoil+Vnew[0]);
  //wfrac =Vnew[0]/(eta[0]*TemperatureEnthalpyDerivative(hold[0])+Vnew[0]);//pct of conductive flux going to water
  //wfrac =1.0;
  condflux=MeanConvectiveFlux(hold[0],Vnew[0],hcond,Tair,Options.timestep);

  if(pHRU->GetSnowDepth()>100) {
 //   condflux=0;//for now
  }
  rates[0    ]=(    wfrac)*condflux;
  //rates[2*N+1]=(1.0-wfrac)*condflux;

  //conductive heat fluxes
  //-----------------------------------------------------------------------
  if(!zfx) { //if the problem was properly solved,
    for(int m=1;m<N;m++)
    {
      kappal    =pow(kap [m]*kap [m-1],0.5); // geometric mean
      kappaln   =pow(kapn[m]*kapn[m-1],0.5);
      rates[m] =-0.5*kappal  * 1.0/(0.5*(dz[m]+dz[m-1]))*(Told[m]-Told[m-1]);// [MJ/m2/d]
      rates[m]+=-0.5*kappaln * 1.0/(0.5*(dz[m]+dz[m-1]))*(Tnew[m]-Tnew[m-1]);
    }
  }
  // geothermal heat flux
  //-----------------------------------------------------------------------
  wfrac   =Vnew[N-1]/(dz[N-1]*(1.0-poro[N-1])+Vnew[N-1]);
  //  wfrac   =Vnew[0]/(eta[0]*TemperatureEnthalpyDerivative(hold[0])+Vnew[0]);//pct of conductive flux going to water (can go to zero now for Vnew=deriv=0.0);
  rates[N]=wfrac*kap[N]*geo_grad; // [MJ/m2/d]

  if(!zfx) {
    for(int m=0;m<N;m++)
    {
      //-dHs/dt = - eta*dT_w/dt  = gain of heat by soil, lost by water if Tnew>Told
      rates[m+N+1] = -eta[m]*(Tnew[m]-Told[m])/Options.timestep; // energy gained by soil [MJ/m2/d]
      //cout<<" lost: "<<rates[m+N+1]<<" "<<eta[m]<<" "<<Tnew[m]<<" "<<Told[m]<<endl;
    }
  }

  // Delete static arrays
  //-----------------------------------------------------------------------
  //if (pModel->LastHRULastTimeStep())
  if((tt.model_time>=Options.duration-Options.timestep-TIME_CORRECTION) &&
    (pHRU->GetGlobalIndex()==((CModel*)(pModel))->GetNumHRUs()-1))
  {
    for(k=0;k<_nHRUs;k++) {delete [] v_old[k];} delete [] v_old;

    delete[] dz;  delete[] z;    delete[] poro;
    delete[] sat; delete[] satn; delete[] kappa_s;
    delete[] eta; delete[] Vold; delete[] Vnew;
    delete[] hold;
    delete[] Tnew;   delete[] Told;    delete[] kap; delete[] kapn;
    delete[] hguess; delete[] delta_h; delete[] f;
    for(i=0;i<3;i++) { delete[] J[i]; } delete[] J;
    for(i=0;i<N;i++) { delete[]  Jinv[i]; } delete[] Jinv;
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
//////////////////////////////////////////////////////////////////
/// \brief calculates mean convective flux over timestep when freezing may impact temperature
/// \param hn [in] volumetric enthalpy at start of timestep [MJ/m3]
/// \param Vw [in] water volume [m]
/// \param alpha [in] convection coefficient [MJ/m2/d/K]
/// \param Ta [in] air temperature for timestep [C]
/// \param dt [out] timestep [d]
/// \return mean convective flux [MJ/m2/d]
//
double MeanConvectiveFlux(const double &hn,const double &Vw,const double &alpha,const double &Ta, const double &dt)
{
  double ht,hlim;
  double hnew=hn;
  double ha;
  double tstar;
  double gamma_w=alpha/HCP_WATER/Vw;
  double gamma_i=alpha/HCP_ICE  /Vw;

  ht=-LH_FUSION*DENSITY_WATER;
  //ht=(g_freeze_temp*SPH_ICE-LH_FUSION)*DENSITY_WATER
  ha=HCP_WATER*Ta;  if (Ta<0.0){ha=HCP_ICE*Ta+ht;}

  if(Vw<1e-6) { return (ha-hn)*Vw/dt; }

  if (Ta>0.0) {
    if(hn>0.0)      { //unfrozen (case v)
      hnew=hn+(ha-hn)*(1.0-exp(-gamma_w*dt));
    }
    else if(hn<ht)  { //frozen and warming (case i) -working
      hlim=(Ta*HCP_ICE+ht);
      tstar=-gamma_i*log((hlim-ht)/(hlim-hn));
      if (dt<tstar) {hnew=hlim+(hn-hlim)*exp(-gamma_i*dt); }
      else          {hnew=ht+(dt-tstar)*gamma_i*ha;        }
    }
    else            { //partially frozen and warming (case ii) -working
      tstar=-hn/gamma_w/ha;
      if(dt<tstar)  { hnew=hn+gamma_w*ha*dt;                  }
      else          { hnew=ha*(1.0-exp(-gamma_w*(dt-tstar))); }
    }
  }
  else
  {
    if(hn<ht)       { //frozen (case v) -working
      hnew=hn+(ha-hn)*(1.0-exp(-gamma_w*dt));
    }
    else if(hn>0.0) { //unfrozen and cooling (case iii) -working
      hlim=HCP_WATER*Ta;
      tstar=-1/gamma_w*log((hlim)/(hlim-hn));
      if(dt<tstar) { hnew=hn+(hlim-hn)*(1.0-exp(-gamma_i*dt)); }
      else          {hnew=(dt-tstar)*gamma_i*(ha-ht);               }
    }
    else            { //partially frozen and cooling (case iv) -working
      tstar=(ht-hn)/gamma_i/(ha-ht);
      if(dt<tstar)  { hnew=hn+gamma_i*(ha-ht)*dt;                  }
      else          { hnew=ht+(ha-ht)*(1.0-exp(-gamma_i*(dt-tstar))); }
    }
  }

  return (hnew-hn)*Vw/dt;
}
//////////////////////////////////////////////////////////////////
/// \brief Unit test for analytical convection solution
//
void TestConvectionSolution() {
  ofstream CONV;
  double Vw=1.0;
  double alpha=3.0;
  double ht=-LH_FUSION*DENSITY_WATER;
  double h,h0;
  CONV.open("AnalyticConvTest.csv");
  CONV<<"t,case i,T,case ii,T,case iii,T,case iv,T,"<<endl;
  for(double t=0.001;t<6;t+=0.2) {
    CONV<<t<<",";
    h0=4.0;
    h=Vw*t*MeanConvectiveFlux(h0,Vw,alpha,3.0,t)+h0;
    CONV<<h<<","<<ConvertVolumetricEnthalpyToTemperature(h)<<",";
    h0=-30;
    h=Vw*t*MeanConvectiveFlux(h0,Vw,alpha,3.0,t)+h0;;
    CONV<<h<<","<<ConvertVolumetricEnthalpyToTemperature(h)<<",";
    h0=ht-8;
    h=Vw*t*MeanConvectiveFlux(h0,Vw,alpha,3.0,t)+h0;;
    CONV<<h<<","<<ConvertVolumetricEnthalpyToTemperature(h)<<",";

    h0=4;
    h=Vw*t*MeanConvectiveFlux(h0,Vw,alpha,-3.0,t)+h0;;
    CONV<<h<<","<<ConvertVolumetricEnthalpyToTemperature(h)<<",";
    h0=ht+30;
    h=Vw*t*MeanConvectiveFlux(h0,Vw,alpha,-3.0,t)+h0;;
    CONV<<h<<","<<ConvertVolumetricEnthalpyToTemperature(h)<<",";
    h0=ht-8;
    h=Vw*t*MeanConvectiveFlux(h0,Vw,alpha,-3.0,t)+h0;;
    CONV<<h<<","<<ConvertVolumetricEnthalpyToTemperature(h)<<",";
    CONV<<endl;
  }
  CONV.close();
  ExitGracefully("TestConvectionSolution",SIMULATION_DONE);
}
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
