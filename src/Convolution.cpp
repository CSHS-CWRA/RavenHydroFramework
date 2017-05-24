/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ------------------------------------------------------------------
  Convolution (routing water through UH)
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "Convolution.h"

/*****************************************************************
   Convolution Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of convolution constructor
/// \param absttype [in] Selected model of abstraction
//
CmvConvolution::CmvConvolution(convolution_type type,
                               const int        to_index)
  :CHydroProcessABC(CONVOLVE)
{
  _type=type;
  _nConv++; //starts at -1
  _iTarget=to_index;

  int N=MAX_CONVOL_STORES;

  CHydroProcessABC::DynamicSpecifyConnections(2*MAX_CONVOL_STORES);
  for (int i=0;i<N;i++)
  {
    //for applying convolution
    iFrom[i]=pModel->GetStateVarIndex  (CONV_STOR,i+_nConv*MAX_CONVOL_STORES);
    iTo  [i]=_iTarget;

    //for shifting storage history
    if (i<N-1){
      iFrom[N+i]=pModel->GetStateVarIndex(CONV_STOR,i  +_nConv*MAX_CONVOL_STORES);
      iTo  [N+i]=pModel->GetStateVarIndex(CONV_STOR,i+1+_nConv*MAX_CONVOL_STORES);
    }
  }
  //updating total convolution storage
  iFrom[2*N-1]=pModel->GetStateVarIndex(CONVOLUTION,_nConv);
  iTo  [2*N-1]=pModel->GetStateVarIndex(CONVOLUTION,_nConv);

}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvConvolution::~CmvConvolution(){}
//////////////////////////////////////////////////////////////////
/// \brief static variable initialization
//
int CmvConvolution::_nConv=-1;
//////////////////////////////////////////////////////////////////
/// \brief Initializes convolution object
//
void   CmvConvolution::Initialize()
{


}
//////////////////////////////////////////////////////////////////
/// \brief unit S-hydrograph (cumulative hydrograph) for GR4J #2
//
double GR4J_SH2(const double &t, const double &x4)
{
  if (t/x4<1.0)     {return 0.5*pow(t/x4,2.5);}
  else if (t/x4<2.0){return 1.0-0.5*pow(2-t/x4,2.5);}
  else              {return 1.0;}
}
void   CmvConvolution::GenerateUnitHydrograph(const CHydroUnit *pHRU, const optStruct &Options, double *aUnitHydro, int &N) const
{
  //generates unit hydrograph based upon HRU parameters (called every timestep
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

  N =(int)(ceil(max_time/tstep));
  if (N>MAX_CONVOL_STORES) { printf("N = %i    MAX_CONVOL_STORES = %i ",N,MAX_CONVOL_STORES ); }
  ExitGracefullyIf(N>MAX_CONVOL_STORES,"CmvConvolution::GenerateUnitHydrograph: unit hydrograph duration for convolution too long",BAD_DATA);

  if (N==0){N=1;}
  ExitGracefullyIf(max_time<=0.0,"CmvConvolution::GenerateUnitHydrograph: negative or zero duration of unit hydrograph (bad GR4J_x4)",BAD_DATA);

  if (_type==CONVOL_TRIANGLE){
    for (int n=0; n<N; n++)
    {
      // aUnitHydro[i]=
    }
  }
  else if (_type==CONVOL_GR4J_1)
  {
    double x4=pHRU->GetSurfaceProps()->GR4J_x4;
    for (int n=0; n<N; n++)
    {
      aUnitHydro[n]=min(pow((n+1)*tstep/x4,2.5),1.0)-min(pow(n*tstep/x4,2.5),1.0);
    }
  }
  else if (_type==CONVOL_GR4J_2){
    double x4=pHRU->GetSurfaceProps()->GR4J_x4;
    for (int n=0; n<N; n++)
    {
      aUnitHydro[n]=GR4J_SH2((n+1)*tstep,x4)-GR4J_SH2(n*tstep,x4);
    }
  }
  //cout<<"** ";for (int i=0;i<N;i++){cout<<aUnitHydro[i]<<" ";} cout<<endl;
  //correct to ensure that sum _aUnitHydro[m]=1.0
  double sum=0.0;
  for (int n=0;n<N;n++){sum+=aUnitHydro[n];}
  ExitGracefullyIf(sum==0.0,"CmvConvolution::GenerateUnitHydrograph: bad unit hydrograph constructed",RUNTIME_ERR);
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
  else
  {
    ExitGracefully("CmvConvolution::GetParticipatingParamList: undefined abstraction algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variable list
///
/// \param absttype [in] Selected abstraction model
/// \param *aSV [out] Reference to array of state variables needed by abstraction algorithm
/// \param *aLev [out] Array of levels of multilevel state variables (or DOESNT_EXIST of single level)
/// \param &nSV [out] Number of participating state variables (length of aSV and aLev arrays)
//
void CmvConvolution::GetParticipatingStateVarList(convolution_type absttype, sv_type *aSV, int *aLev, int &nSV)
{
  nSV=MAX_CONVOL_STORES+1;
  for (int i=0;i<MAX_CONVOL_STORES;i++)
  {
    aSV [i]=CONV_STOR;
    aLev[i]=i+(_nConv+1)*MAX_CONVOL_STORES;
  }
  aSV [MAX_CONVOL_STORES]=CONVOLUTION;
  aLev[MAX_CONVOL_STORES]=(_nConv+1);
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
void   CmvConvolution::GetRatesOfChange( const double                   *state_vars,
                                         const CHydroUnit       *pHRU,
                                         const optStruct        &Options,
                                         const time_struct &tt,
                                         double     *rates) const
{
  int i;
  double TS_old;
  double tstep=Options.timestep;
  static double S[MAX_CONVOL_STORES];
  static double aUnitHydro[MAX_CONVOL_STORES];

  int N =0;
  GenerateUnitHydrograph(pHRU,Options,&aUnitHydro[0],N);

  //Calculate S[0] from Convolution total storage
  TS_old=state_vars[iFrom[2*MAX_CONVOL_STORES-1]]; //total storage after water added to convol stores earlier in process list
  double sum(0.0);
  for (i=1;i<N;i++){
    S[i]=state_vars[iFrom[i]];
    sum+=S[i];
  }
  S[0]=TS_old-sum; //amount of water added this time step to convol stores
  //cout<<TS_old<<" "<<sum<<endl;
  //cout<<"** ";for (int i=0;i<N;i++){cout<<S[i]<<" ";} cout<<endl;

  //outflow from convolution storage
  double orig_storage;
  double sumrem=1.0;
  for (int i=0; i<N; i++)
  {
    if (sumrem<REAL_SMALL){orig_storage=0.0;}
    else                  {orig_storage=S[i]/sumrem;}
    rates[i]=aUnitHydro[i]*orig_storage/tstep;
    S[i]-=rates[i]*tstep;
    sumrem-=aUnitHydro[i];
  }//sumrem should ==0 at end, so should S[N-1]
  //cout<<sumrem<<endl;
  //cout<<"** ";for (int i=0;i<N;i++){cout<<S[i]<<" ";} cout<<endl;

  //time shift convolution history
  double move;
  for (int i=N-2; i>=0; i--)
  {
    move=S[i];
    rates[MAX_CONVOL_STORES+i]=move/tstep;
    S[i  ]-=move;
    S[i+1]+=move;
  }
  //cout<<"*> ";for (int i=0;i<N;i++){cout<<S[i]<<" ";} cout<<endl;

  //update total convolution storage:
  double TS_new=0;
  for (int i=1; i<N; i++){TS_new+=S[i];}
  rates[2*MAX_CONVOL_STORES-1]=(TS_new-TS_old)/tstep;
  //cout<<TS_old<<" "<<TS_new<<endl;
  //cout<<"* ----------------------"<<endl;

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
//no constraints?
}
