/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2023 the Raven Development Team
  ------------------------------------------------------------------
  Convolution (routing water through UH)
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "Convolution.h"

//////////////////////////////////////////////////////////////////
/// \brief static variable initialization
//
//int  CmvConvolution::_nStores=20; //TESTING THIS - NOT FOR RELEASE
//bool CmvConvolution::_smartmode=true;

int  CmvConvolution::_nStores=MAX_CONVOL_STORES;
bool CmvConvolution::_smartmode=false;

/*****************************************************************
   Convolution Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of convolution constructor
/// \param absttype [in] Selected model of abstraction
/// \param pModel [in] pointer to model
/// \param conv_index [in] index of this convolution process
//
CmvConvolution::CmvConvolution(convolution_type type,
                               const int        to_index,
                               CModelABC       *pModel,
                               int              conv_index)
  :CHydroProcessABC(CONVOLVE, pModel)
{
  _type = type;
  _iTarget = to_index;

  if (!_smartmode){_nStores=MAX_CONVOL_STORES;}

  int N=_nStores; //shorthand

  CHydroProcessABC::DynamicSpecifyConnections(2*N);
  for (int i=0; i<N; i++)
  {
    //for applying convolution
    iFrom[i] = pModel->GetStateVarIndex(CONV_STOR, i+conv_index*N);
    iTo  [i] = _iTarget;

    //for shifting storage history
    if (i < N-1){
      iFrom[N+i] = pModel->GetStateVarIndex(CONV_STOR, i  +conv_index*N);
      iTo  [N+i] = pModel->GetStateVarIndex(CONV_STOR, i+1+conv_index*N);
    }
  }
  //updating total convolution storage
  iFrom[2*N-1] = pModel->GetStateVarIndex(CONVOLUTION, conv_index);
  iTo  [2*N-1] = pModel->GetStateVarIndex(CONVOLUTION, conv_index);
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvConvolution::~CmvConvolution(){}

//////////////////////////////////////////////////////////////////
/// \brief Initializes convolution object
//
void   CmvConvolution::Initialize(){}

//////////////////////////////////////////////////////////////////
/// \brief unit S-hydrograph (cumulative hydrograph) for all options
//
double CmvConvolution::LocalCumulDist(const double& t, const CHydroUnit* pHRU) const
{
  if (_type==CONVOL_GR4J_1)
  {
    double x4=pHRU->GetSurfaceProps()->GR4J_x4;
    return min(pow(t/x4,2.5),1.0);
  }
  else if (_type==CONVOL_GR4J_2)
  {
    double x4=pHRU->GetSurfaceProps()->GR4J_x4;
    if      (t/x4<1.0){return 0.5*pow(t/x4,2.5);}
    else if (t/x4<2.0){return 1.0-0.5*pow(2-t/x4,2.5);}
    else              {return 1.0;}
  }
  else if(_type==CONVOL_GAMMA)
  {
    double beta =pHRU->GetSurfaceProps()->gamma_scale;
    double alpha=pHRU->GetSurfaceProps()->gamma_shape;
    return GammaCumDist(t,alpha,beta);
  }
  else if(_type==CONVOL_GAMMA_2){
    double beta =pHRU->GetSurfaceProps()->gamma_scale2;
    double alpha=pHRU->GetSurfaceProps()->gamma_shape2;
    return GammaCumDist(t,alpha,beta);
  }
  return 1.0;
}

void   CmvConvolution::GenerateUnitHydrograph(const CHydroUnit *pHRU, const optStruct &Options, double *aUnitHydro, int *aInterval, int &N) const
{
  //generates unit hydrograph based upon HRU parameters
  //(called every timestep because it is potentially different in every HRU)
  double tstep=Options.timestep;
  double max_time(0);

  if (_type==CONVOL_GR4J_1)
  {
    double x4=pHRU->GetSurfaceProps()->GR4J_x4;
    max_time=x4;
  }
  else if (_type==CONVOL_GR4J_2){
    double x4=pHRU->GetSurfaceProps()->GR4J_x4;
    max_time=2*x4;
  }
  else if(_type==CONVOL_GAMMA){
    double alpha=pHRU->GetSurfaceProps()->gamma_shape;
    double beta =pHRU->GetSurfaceProps()->gamma_scale;
    max_time=4.5*pow(alpha,0.6)/beta; /*empirical estimate of tail, done in-house @UW*/
  }
  else if(_type==CONVOL_GAMMA_2){
    double alpha=pHRU->GetSurfaceProps()->gamma_shape2;
    double beta =pHRU->GetSurfaceProps()->gamma_scale2;
    max_time=4.5*pow(alpha,0.6)/beta;
  }
  ExitGracefullyIf(max_time<=0.0,"CmvConvolution::GenerateUnitHydrograph: negative or zero duration of unit hydrograph (bad GR4J_x4, or GAMMA_SHAPE/GAMMA_SCALE)",BAD_DATA);

  int NN;

  double sum=0;;
  if (!_smartmode)
  {
    max_time=min(MAX_CONVOL_STORES*tstep,max_time); //must truncate
    N =(int)(ceil(max_time/tstep));
    if (N==0){N=1;}
    NN=N;
    for (int n=0;n<N;n++)
    {
      aUnitHydro[n]=LocalCumulDist((n+1)*tstep,pHRU)-sum;
      sum+=aUnitHydro[n];

      aInterval[n]=1;
    }
  }
  else /*smart allocation - simplifies history storage*/
  {
    N=_nStores;
    NN=(int)(ceil(max_time/tstep));

    int    i=0;
    int    nlast=-1;
    double Flast=0.0;
    double Fn;

    for (int n = 0; n < NN; n++)
    {
      Fn=LocalCumulDist((n+1)*tstep,pHRU);
      if ((i == 0) && (n>0) && (Fn > 1.0 / _nStores / 4.0)) { //handles lag effect
        aInterval [i]=( n-nlast-1);
        aUnitHydro[i]=0;
        i++;
        aInterval[i]=1;
        aUnitHydro[i]=Fn-Flast;
        Flast=Fn;
        nlast=n;
        i++;
      }
      else if ((Fn >= Flast + 1.0 / _nStores) || (Fn>0.99)) { //CDF has increased enough to start new interval
        aInterval [i]=( n-nlast);
        aUnitHydro[i]=(Fn-Flast);
        sum+=aUnitHydro[i];
        Flast=Fn;
        nlast=n;
        i++;
      }
    }
    for (int ii = i; ii < N; ii++) {
      aUnitHydro[ii]=0.0;
      aInterval [ii]=1;
    }
    N=i;
  }

  ExitGracefullyIf(sum==0.0,"CmvConvolution::GenerateUnitHydrograph: bad unit hydrograph constructed",RUNTIME_ERR);

  /*static bool first=true;
  if (first) {
    ofstream CONV;
    CONV.open("convolution_test.csv");
    CONV<<"n,interval,orig,smart"<<endl;
    sum=0;
    int i=0;
    int last=0;
    for (int n=0;n<NN;n++)
    {
      double tmp=LocalCumulDist((n+1)*tstep,pHRU);
      CONV<<n<<","<<i<<","<<tmp-sum<<",";
      sum+=(tmp-sum);

      CONV<<aUnitHydro[i]/aInterval[i]<<endl;
      if (n+1>=last+aInterval[i]){
        last=n+1;
        i++;
      }
    }
    CONV.close();
    first=false;
  }*/
  //correct to ensure that sum _aUnitHydro[m]=1.0
  for (int n=0;n<N;n++){aUnitHydro[n]/=sum;}

}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for abstraction algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by abstraction algorithm (size of aP[] and aPC[])
//
void CmvConvolution::GetParticipatingParamList(string  *aP , class_type *aPC , int &nP) const
{
  if ((_type==CONVOL_GR4J_1) || (_type==CONVOL_GR4J_2))
  {
    nP=1;
    aP[0]="GR4J_X4";                    aPC[0]=CLASS_LANDUSE;
  }
  else if (_type==CONVOL_GAMMA)
  {
    nP=2;
    aP[0]="GAMMA_SHAPE";                aPC[0]=CLASS_LANDUSE;
    aP[1]="GAMMA_SCALE";                aPC[1]=CLASS_LANDUSE;
  }
  else if(_type==CONVOL_GAMMA_2)
  {
    nP=2;
    aP[0]="GAMMA_SHAPE2";               aPC[0]=CLASS_LANDUSE;
    aP[1]="GAMMA_SCALE2";               aPC[1]=CLASS_LANDUSE;
  }
  else
  {
    ExitGracefully("CmvConvolution::GetParticipatingParamList: undefined convolution algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variable list
///
/// \param absttype [in] Selected abstraction model
/// \param *aSV [out] Reference to array of state variables needed by abstraction algorithm
/// \param *aLev [out] Array of levels of multilevel state variables (or DOESNT_EXIST of single level)
/// \param &nSV [out] Number of participating state variables (length of aSV and aLev arrays)
/// \param conv_index [in] Index of this convolution process in model (from CModel::_nConvVariables)
//
void CmvConvolution::GetParticipatingStateVarList(convolution_type absttype, sv_type *aSV, int *aLev, int &nSV, int conv_index)
{
  nSV = _nStores+1;

  for (int i=0;i<_nStores;i++)
  {
    aSV [i]=CONV_STOR;
    aLev[i]=(conv_index*_nStores)+i;
  }
  aSV [_nStores]=CONVOLUTION;
  aLev[_nStores]=conv_index;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rates of water loss to abstraction
//
/// \param *state_vars [in] Array of state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time at which abstraction is to be calculated
/// \param *rates [out] rates[0]= rate of abstraction [mm/d]
//
void   CmvConvolution::GetRatesOfChange( const double      *state_vars,
                                         const CHydroUnit  *pHRU,
                                         const optStruct   &Options,
                                         const time_struct &tt,
                                               double      *rates) const
{
  int i;
  double TS_old;
  double tstep=Options.timestep;
  static double S         [MAX_CONVOL_STORES];
  static double aUnitHydro[MAX_CONVOL_STORES];
  static int    aInterval [MAX_CONVOL_STORES];

  int N =0;
  GenerateUnitHydrograph(pHRU,Options,&aUnitHydro[0],&aInterval[0],N); //THIS IS SLOW!! - create aUnitHydro as process array

  //Calculate S[0] as change in convolution total storage
  TS_old=state_vars[iFrom[2*_nStores-1]]; //total storage after water added to convol stores earlier in process list
  double sum(0.0);
  int NN=0;
  for (i=0;i<N;i++){
    S[i]=state_vars[iFrom[i]];
    sum+=S[i];
    NN+=aInterval[i];
  }
  S[0]+=(TS_old-sum); //amount of water added this time step to convol stores

  //outflow from convolution storage
  double orig_storage;
  double sumrem;
  if (!_smartmode){
    //usually, we just would do Q=sum UH(n)*R(t-n), but we store depleted R rather than original inputs
    sumrem=1.0;
    for (i=0; i<N; i++)
    {
      if (sumrem<REAL_SMALL){orig_storage=0.0;}
      else                  {orig_storage=S[i]/sumrem;}
      rates[i]=aUnitHydro[i]*orig_storage/tstep; //water released to target store from internal store
      S[i]-=rates[i]*tstep;
      sumrem-=aUnitHydro[i];
    }//sumrem should ==0 at end, so should S[N-1]
    //challengee here - original storage is not necessarily S/sum_0^i-1 UH because of partial shifts of mass
  }
  if (_smartmode){
    int sumn=0;
    sumrem=1.0;
    i=0;
    double rel;
    for (int n = 0; n < NN; n++) {
      if (i<_nStores){
        if (sumrem<REAL_SMALL){orig_storage=0.0;}
        else                  {orig_storage=S[i]/aInterval[i];}
        rel=aUnitHydro[i]*orig_storage/tstep;
        rates[i]+=rel;//water released to target store from internal store
        S[i]-=rel*tstep;

        sumrem-=aUnitHydro[i]/aInterval[i];
        if ((n+1 - sumn)>=aInterval[i]) {
          sumn=n+1;
          i++;
        }
      }
      else {
        rates[i]=S[i]/tstep;
        S[i]=0;
      }
    }
  }
  //cout<<"convend: "<<sumrem<<" "<<S[N-1]<<endl;

  //time shift convolution history
  double move;
  for (i=N-2; i>=0; i--)
  {
    move=S[i]/aInterval[i];
    rates[_nStores+i]=move/tstep;
    S[i  ]-=move;
    S[i+1]+=move;
  }

  //update total convolution storage:
  double TS_new=0;
  for (i=1; i<N; i++){TS_new+=S[i];}
  rates[2*_nStores-1]=(TS_new-TS_old)/tstep;
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow cannot drain "from" compartment over timestep.
///
/// \param *state_vars [in] Array of state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time at which abstraction is to be calculated
/// \param *rates [out] rates[0]= rate of abstraction [mm/d]
//
void   CmvConvolution::ApplyConstraints(const double             *state_vars,
                                        const CHydroUnit *pHRU,
                                        const optStruct      &Options,
                                        const time_struct &tt,
                                        double     *rates) const
{
  //no constraints.
}
