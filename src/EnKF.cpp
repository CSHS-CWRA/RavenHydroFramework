/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2022 the Raven Development Team
----------------------------------------------------------------*/

#include "EnKF.h"
#include "Matrix.h"

bool IsContinuousFlowObs(const CTimeSeriesABC* pObs,long SBID);
bool ParseInitialConditionsFile(CModel*& pModel,const optStruct& Options);
bool ParseTimeSeriesFile(CModel*& pModel,const optStruct& Options);

string FilenamePrepare(string filebase,const optStruct& Options); //Defined in StandardOutput.cpp

//step 1) get perturbations working - run in ensemble mode [DONE]
//step 2) Get order working - must write soln at t0 and reboot e>N ensembles at t0 with this soln (to avoid duplicating perturbations) [DONE]
//step 3) get assimilation matrix calcs running [DONE]
//step 4) build .R visualization post processor [DONE]
//step 5) Implement observation perturbations (eQ treatment) [DONE]
//step 6) Test 
//step 7) Documentation 

// in .rvi file:
//  :EnsembleMode ENSEMBLE_ENKF [double #members]
//  :AssimilationStartTime 2022-10-02 00:13:00
//  :WarmEnsemble fc2 # used if ICs should be read from ensemble solution.rvc file
//  :SilentMode # useful
// 
// in .rve file:
//  :DataHorizon 1 # no. of timesteps (1 for standard EnKF or huge if all data since sim start is used; 2+ for variational approach)
//
//  :ForecastRVTFilename ./meteo/model_forecast.rvt 
// 
//  :ForcingPerturbation RAINFALL GAMMA_DIST [dist param1] [dist param2] {HRU_Group}
//  :ForcingPerturbation TEMP_AVE NORMAL_DIST [mean=0] [std dev] {HRU_Group}
//  :ForcingPerturbation RECHARGE GAMMA_DIST [dist param1] [dist param2] {HRU_Group}
// 
//  :AssimilatedState STREAMFLOW {SubBasinGroup} # only STREAMFLOW to be supported initially
//  :AssimilatedState SNOW       {HRUGroup}
//     
//  //plus, uses streamflow locations tagged for assimilation 
//  :AssimilateStreamflow [SBID1] # (optionally in .rvt file)
//  :AssimilateStreamflow [SBID2]
//   later:
//  :AssimilateHRUState SNOW [HRU1] 
// 
//////////////////////////////////////////////////////////////////
/// \brief EnKF Ensemble Construcutor 
/// \param num_members [in] EnKF ensemble size 
/// \param &Options [out] Global model options information
//
CEnKFEnsemble::CEnKFEnsemble(const int num_members,const optStruct &Options)
  :CEnsemble(num_members*2 ,Options) //*2 for pre and post-assimilation runs 
{
  _type=ENSEMBLE_ENKF;
  _nEnKFMembers=num_members;

  _pPerturbations=NULL; //Specified during .rve file parse
  _nPerturbations=0;
  _pObsPerturbations=NULL;
  _nObsPerturbations=0;
  _aAssimStates=NULL;
  _aAssimLayers=NULL;
  _aAssimGroupID=NULL;
  _nAssimStates=0;
  _t_assim_start=-1.0;

  _state_matrix=NULL; //Built in ::Initialize
  _nStateVars=0;

  _warm_ensemble=false;
  _warm_runname="";
  _forecast_rvt="";

  _obs_matrix=NULL;
  _output_matrix=NULL;
  _noise_matrix=NULL;
  _nObsDatapoints=0;

  _window_size=1;
  _nTimeSteps=0;

  _truncate_hind=true; //default behaviour - hindcasts stop at t0
  _duration_stored=Options.duration;
}
//////////////////////////////////////////////////////////////////
/// \brief EnKF Ensemble Destrucutor 
//
CEnKFEnsemble::~CEnKFEnsemble()
{
  for(int e=0;e<_nEnKFMembers;e++) {
    delete [] _state_matrix [e]; _state_matrix [e]=NULL;
    delete [] _obs_matrix   [e]; _obs_matrix   [e]=NULL;
    delete [] _output_matrix[e]; _output_matrix[e]=NULL;
    delete [] _noise_matrix [e]; _noise_matrix [e]=NULL;
  }
  delete [] _state_matrix;
  delete [] _obs_matrix;
  delete [] _output_matrix;
  delete [] _noise_matrix;

  for (int i=0;i<_nPerturbations;   i++){delete _pPerturbations[i];   } delete [] _pPerturbations;
  for (int i=0;i<_nObsPerturbations;i++){delete _pObsPerturbations[i];} delete [] _pObsPerturbations;
  
  delete [] _aAssimStates;
  delete [] _aAssimLayers;
  delete [] _aAssimGroupID;
  delete [] _aObsIndices;
}
//////////////////////////////////////////////////////////////////
/// \brief adds additional forcing perturbation
/// \param type [in] forcing type to be perturbed
/// \param distrib [in] sampling distribution type 
/// \param distpars [in] array of distribution parameters
/// \param group_index [in] HRU group ID (or DOESNT_EXIST if all HRUs should be perturbed)
/// \param adj [in] type of perturbation (e.g., additive or multiplicative)
//
void CEnKFEnsemble::AddPerturbation(forcing_type type, disttype distrib, double* distpars, int group_index, adjustment adj) 
{
  force_perturb *pFP=new force_perturb();
  pFP->forcing=type;
  pFP->distribution=distrib;
  for (int i = 0; i < 3; i++) {
    pFP->distpar[i]=distpars[i];
  }
  pFP->kk=group_index;
  pFP->adj_type=adj;

  if (!DynArrayAppend((void**&)(_pPerturbations),(void*)(pFP),_nPerturbations)){
    ExitGracefully("CEnKFEnsemble::AddPerturbation: creating NULL perturbation",BAD_DATA_WARN);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief adds additional state observation perturbation - applied to ALL observations of this type
/// \param type [in] state variable type to be perturbed
/// \param distrib [in] sampling distribution type 
/// \param distpars [in] array of distribution parameters
/// \param adj [in] type of perturbation (e.g., additive or multiplicative)
//
void CEnKFEnsemble::AddObsPerturbation(sv_type type, disttype distrib, double* distpars, adjustment adj) 
{
  obs_perturb *pOP=new obs_perturb();
  pOP->state=type;
  pOP->distribution=distrib;
  for (int i = 0; i < 3; i++) {
    pOP->distpar[i]=distpars[i];
  }
  pOP->kk=DOESNT_EXIST;
  //pOP->kk=group_index;
  pOP->adj_type=adj;

  if (!DynArrayAppend((void**&)(_pObsPerturbations),(void*)(pOP),_nObsPerturbations)){
    ExitGracefully("CEnKFEnsemble::AddObsPerturbation: creating NULL observation perturbation",BAD_DATA_WARN);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief adds additional state variable to 
/// \param sv [in] state variable type updated in assimilation
/// \param assim_groupID [in] HRU or subbasin group ID (or DOESNT_EXIST if all should be used)
//
void CEnKFEnsemble::AddAssimilationState(sv_type sv, int lay, int assim_groupID) 
{  
  _nAssimStates++;
  sv_type *tmp =new sv_type[_nAssimStates];
  int     *tmp2=new     int[_nAssimStates];
  int     *tmp3=new     int[_nAssimStates];
  for (int i = 0; i < _nAssimStates - 1; i++) {
    tmp [i]=_aAssimStates[i];
    tmp2[i]=_aAssimGroupID[i];
    tmp3[i]=_aAssimLayers[i];
  }
  delete [] _aAssimStates;
  delete [] _aAssimGroupID;
  delete [] _aAssimLayers;
  _aAssimStates=tmp;
  _aAssimGroupID=tmp2;
  _aAssimLayers=tmp3;
  _aAssimStates [_nAssimStates-1]=sv;
  _aAssimGroupID[_nAssimStates-1]=assim_groupID;
  _aAssimLayers [_nAssimStates-1]=lay;
}
//////////////////////////////////////////////////////////////////
/// \brief turns on use of warm ensemble
//
void CEnKFEnsemble::SetToWarmEnsemble(string runname){_warm_ensemble=true;_warm_runname=runname;}

//////////////////////////////////////////////////////////////////
/// \brief truncates hindcasts 
//
void CEnKFEnsemble::DontTruncateHindcasts() { _truncate_hind=false; }

//////////////////////////////////////////////////////////////////
/// \brief set value of data horizon
//
void CEnKFEnsemble::SetWindowSize(const int nTimesteps){ _window_size=nTimesteps;}

//////////////////////////////////////////////////////////////////
/// \brief set forecast RVT filename
//
void CEnKFEnsemble::SetForecastRVTFile(string filename){_forecast_rvt=filename;}

//////////////////////////////////////////////////////////////////
/// \brief Accessor - gets ensemble start time (in local model time)
/// \param e [in] ensemble index
/// \return ensemble start time (in local model time)
//
double CEnKFEnsemble::GetStartTime(const int e) const 
{
  if(e<_nEnKFMembers) { return 0.0; } //uncorrected members
  return  _t_assim_start; //members post-state update - only forecasting
}

//////////////////////////////////////////////////////////////////
/// \brief initializes EnKF assimilation run
/// \param &Options [out] Global model options information
//
void CEnKFEnsemble::Initialize(const CModel* pModel,const optStruct &Options)
{
  CEnsemble::Initialize(pModel,Options);

  // QA/QC
  //-----------------------------------------------
  if(_nMembers==0) {
    ExitGracefully("CEnKFEnsemble::Initialize: number of EnKF ensemble members must be >0",BAD_DATA);
  }
  if(_nPerturbations==0) {
    ExitGracefully("CEnKFEnsemble::Initialize: at least one forcing perturbation must be set using :ForcingPerturbation command",BAD_DATA);
  }
  if(_nAssimStates==0) {
    ExitGracefully("CEnKFEnsemble::Initialize: at least one assimilated state must be set in :ForcingPerturbation command",BAD_DATA);
  }
  if(Options.assimilation_start<0) {
    ExitGracefully("CEnKFEnsemble::Initialize: must set :AssimilationStartTime in .rvi file",BAD_DATA_WARN);
  }
  if(Options.assimilation_start>Options.duration) {
    ExitGracefully("CEnKFEnsemble::Initialize: assimilation start time after simulation has ended.",BAD_DATA_WARN);
  }
  if(Options.assimilation_start<0.0) {
    ExitGracefully("CEnKFEnsemble::Initialize: assimilation start time before start time of simulation.",BAD_DATA_WARN);
  }
  _t_assim_start=Options.assimilation_start;

  // determine total number of state variables
  //-----------------------------------------------
  int kk;
  _nStateVars=0;
  for(int i=0;i<_nAssimStates;i++) {
    kk=_aAssimGroupID[i];
    if(_aAssimStates[i]==STREAMFLOW) {
      for(int pp=0;pp<pModel->GetSubBasinGroup(kk)->GetNumSubbasins();pp++)
      {
        CSubBasin* pBasin=pModel->GetSubBasinGroup(kk)->GetSubBasin(pp);
        _nStateVars+=pBasin->GetInflowHistorySize();
        _nStateVars+=pBasin->GetNumSegments();
        _nStateVars++; //Qlast
      }
    }
    else {
      _nStateVars+=pModel->GetHRUGroup(kk)->GetNumHRUs();
    }
  }

  //create empty _state_matrix
  //-----------------------------------------------
  _state_matrix=new double* [_nEnKFMembers];
  ExitGracefullyIf(_state_matrix==NULL,"CEnKFEnsemble::Initialize",OUT_OF_MEMORY);
  for(int e = 0; e < _nEnKFMembers; e++)
  {
    _state_matrix[e]=NULL;
    _state_matrix[e]=new double [_nStateVars]; // gets populated later upon selecting states
    ExitGracefullyIf(_state_matrix[e]==NULL,"CEnKFEnsemble::Initialize(2)",OUT_OF_MEMORY);
    for(int i=0;i<_nStateVars;i++) {
      _state_matrix[e][i]=0.0;
    }
  }
  cout<<"ENKF found "<<_nStateVars<<" state variables for assimilation. State matrix built."<<endl;

  //determine total number of observations 
  //-----------------------------------------------
  _nTimeSteps=(int)(rvn_floor((_t_assim_start+TIME_CORRECTION)/Options.timestep));
  _window_size=min(_window_size,_nTimeSteps);

  bool good;
  double obsval=0;
  _nObsDatapoints=0;
  _nObs=0;
  _aObsIndices=new int [pModel->GetNumObservedTS()];
  for(int i=0; i<pModel->GetNumObservedTS();i++) 
  {
    const CTimeSeriesABC *pTSObs=pModel->GetObservedTS(i);
    long SBID=pTSObs->GetLocID();
    CSubBasin *pSB=pModel->GetSubBasinByID(SBID);
    good= ((pSB!=NULL) && (pSB->UseInFlowAssimilation()));
    
    if ((good) && (IsContinuousFlowObs(pTSObs,SBID)))
    {
      _aObsIndices[_nObs]=i;
      _nObs++;
      for(int nn=_nTimeSteps-_window_size;nn<_nTimeSteps;nn++) {
        obsval=pTSObs->GetSampledValue(nn);
        if(obsval!=RAV_BLANK_DATA) { _nObsDatapoints++; }
      }
    }
  }
  ExitGracefullyIf(_nObs==0          ,"EnKF::Initialize : no observations available for assimilation. Use :AssimilateStreamflow command to indicate enabled observation time series",BAD_DATA);
  ExitGracefullyIf(_nObsDatapoints==0,"EnKF::Initialize : no observations available for assimilation. Blank data during simulation period?",BAD_DATA);
  cout<<"ENKF found "<<_nObsDatapoints<<" valid observation data points for assimilation. Observation matrices built."<<endl;

  //allocate _output_matrix and _noise_matrix
  //-----------------------------------------------
  _obs_matrix   =new double *[_nEnKFMembers];
  _output_matrix=new double *[_nEnKFMembers];
  _noise_matrix =new double *[_nEnKFMembers];
  for(int e=0;e<_nEnKFMembers;e++) {
    _noise_matrix[e]=NULL;
    _obs_matrix   [e]=new double [_nObsDatapoints];
    _output_matrix[e]=new double [_nObsDatapoints];
    _noise_matrix [e]=new double [_nObsDatapoints];
    ExitGracefullyIf(_output_matrix[e]==NULL,"EnKF Initialization: output matrix",OUT_OF_MEMORY);
  }
  
  //populate _output_matrix, _obs_matrix
  //-----------------------------------------------
  int j;
  double eps;
  for(int e=0;e<_nEnKFMembers;e++) 
  {
    j=0;
    for(int ii=0;ii<_nObs;ii++) 
    {
      const CTimeSeriesABC* pTSObs=pModel->GetObservedTS(_aObsIndices[ii]);

      //get perturbation info from observation time series name. if not found, obs data is not perturbed
      // \todo[clean] - this is pretty ugly 
      // \todo[funct] - support other observational types, e.g., lake level
      sv_type sv; int lay;
      obs_perturb *pPerturb=NULL;
      if (!strcmp(pTSObs->GetName().c_str(),"HYDROGRAPH")){sv=STREAMFLOW;}
      else { sv=CStateVariable::StringToSVType(pTSObs->GetName(),lay,true);}
      for (int p = 0; p < _nObsPerturbations; p++) {
        if (_pObsPerturbations[p]->state==sv){pPerturb=_pObsPerturbations[p];}
      }

      for(int nn=_nTimeSteps-_window_size;nn<_nTimeSteps;nn++)
      {
        obsval=pTSObs->GetSampledValue(nn);
        if(obsval!=RAV_BLANK_DATA) {
          _noise_matrix[e][j]=0;
          if (pPerturb!=NULL){
            eps=SampleFromDistribution(pPerturb->distribution,pPerturb->distpar);
            if (pPerturb->adj_type==ADJ_ADDITIVE){
              _noise_matrix[e][j]=eps;
            }
            else if (pPerturb->adj_type == ADJ_MULTIPLICATIVE) {
              _noise_matrix[e][j]=(eps*obsval)-obsval;
            }
          }
          _obs_matrix[e][j]=obsval+_noise_matrix[e][j];
          j++;
        }
      }
    }
  
    for(j=0;j<_nObsDatapoints;j++) {
       _output_matrix[e][j]=0; //filled during simulation
    }
  }

  // Create and open EnKFOutput file
  //-----------------------------------------------
  string filename= FilenamePrepare("EnKFOutput.csv",Options);
  _ENKFOUT.open(filename.c_str());
}
//////////////////////////////////////////////////////////////////
/// \brief called at start of each time step
/// \param pModel [out] pointer to global model instance
/// \param &Options [out] Global model options information
/// \params tt [in] time structure
/// \param e [out] ensembe member index
//
void CEnKFEnsemble::StartTimeStepOps(CModel* pModel,optStruct& Options,const time_struct &tt,const int e)
{
  //Write EnKF-adjusted and unadjusted solutions immediately after state update
  //Update model state matrix
  //WANTED TO MOVE TO CLOSETIMESTEP OPS SO WRITTEN EVEN IF end==t_assim_start 
  //YOU CANT DO THIS BECAUSE ASSIMILATION HANDLED AFTER SIMULATION IS PAST T_0
  if(fabs(tt.model_time-_t_assim_start)<0.1*Options.timestep) 
  {
    pModel->WriteMajorOutput(Options,tt,"solution_t0",false);
  
    if(e<_nEnKFMembers) {
      AddToStateMatrix(pModel,Options,e);
    }

    if (_truncate_hind){
      _disable_output=true; //so output is not written for single time step after assimilation
    }
  }

  if(tt.model_time+TIME_CORRECTION<_t_assim_start)  //only perturb and build output matrix prior to t0
  {
    //Apply forcing perturbations for the current time step
    //----------------------------------------------------------
    int kk=DOESNT_EXIST;
    double eps=0;
    for(int i=0;i<_nPerturbations;i++) {
      kk=_pPerturbations[i]->kk;
      eps=SampleFromDistribution(_pPerturbations[i]->distribution,_pPerturbations[i]->distpar);
      for(int k=0;k<pModel->GetNumHRUs();k++) {
        if((kk==DOESNT_EXIST) || (pModel->GetHRUGroup(kk)->IsInGroup(k))) {
          pModel->GetHydroUnit(k)->AdjustHRUForcing(_pPerturbations[i]->forcing,eps,_pPerturbations[i]->adj_type);
        }
      }
    }

    //Build output Matrix
    //----------------------------------------------------------
    int curr_nn=(int)((tt.model_time+TIME_CORRECTION)/Options.timestep);//current timestep index

    if (curr_nn<_nTimeSteps-_window_size){return;} // haven't reached the data availability period yet

    int j=0;
    double obsval;
    for(int ii=0;ii<_nObs;ii++) { //messy that we have to go through whole loop like this
      const CTimeSeriesABC* pTSObs=pModel->GetObservedTS (_aObsIndices[ii]);
      const CTimeSeriesABC* pTSMod=pModel->GetSimulatedTS(_aObsIndices[ii]);
      for(int nn=_nTimeSteps-_window_size;nn<_nTimeSteps;nn++) {
        obsval=pTSObs->GetSampledValue(nn);
        if(obsval!=RAV_BLANK_DATA) {
          if(nn==curr_nn) {
            _output_matrix[e][j]=pTSMod->GetSampledValue(nn); break;
          }
          j++;
        }
      }
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief called at end of each time step
/// \param pModel [out] pointer to global model instance
/// \param &Options [out] Global model options information
/// \params tt [in] time structure
/// \param e [out] ensembe member index
//
void CEnKFEnsemble::CloseTimeStepOps(CModel* pModel,optStruct& Options,const time_struct& tt,const int e)
{
  
}
//////////////////////////////////////////////////////////////////
/// \brief the EnKF assimilation matrix calculations for updating states
/// Based upon Mandel, J., Efficient Implementation of the Ensemble Kalman Filter, Report, Univ of Colorado, 2006.
//
void CEnKFEnsemble::AssimilationCalcs()
{
  //return;
  int i,j,k;

  //some shorthand to clean up local variable names
  int N       =_nEnKFMembers;   //# of enkf members
  int M       =_nStateVars;     //# of assimilated state vars
  int Nobs    =_nObsDatapoints; //# of observed state variables 

  double** X  =_state_matrix;   //state matrix is [NxM]
  double** eQT=_noise_matrix; //outerr matrix is [NxNobs]

  double** HA;   //output matrix difference from ensemble mean [NobsxN] 
  double** HAT;  //HA transpose [NxNobs]
  double** A;    //prediction ensemble variation matrix [MxN] 
  double** P;    //inverted matrix term 
  double** eQ;   
  double** tmp;
  double** MM;
  double** Z;
  double** X_deltaT;
  double** X_delta;
  double* ans =new double[Nobs];
  double* diff=new double[Nobs];

  AllocateMatrix(M,N,A);     
  AllocateMatrix(Nobs,N,HA); 
  AllocateMatrix(N,Nobs,HAT);
  AllocateMatrix(Nobs,Nobs,P);
  AllocateMatrix(Nobs,N,eQ); 
  AllocateMatrix(Nobs,Nobs,tmp);
  AllocateMatrix(Nobs,N,MM); 
  AllocateMatrix(N,N,Z); 
  AllocateMatrix(M,N,X_deltaT);
  AllocateMatrix(N,M,X_delta);

  //HA=O_sim-1/N*(O_sim*e_N1)*e_1N
  double outMean=0;
  for(j = 0; j < Nobs; j++) {
    outMean=0;
    for(i = 0; i < N; i++) {outMean+=_output_matrix[i][j]/N;}
    for(i = 0; i < N; i++) {
      HA[j][i]=_output_matrix[i][j]-outMean; // Nobs x N
    }
  }
  
  //A =X -1/N*(X *e_N1)*e_1N 
  double Xmean=0;
  for(k = 0; k < M; k++) {
    Xmean=0;
    for(i = 0; i < N; i++) {Xmean+=X[i][k]/N;}
    for(i = 0; i < N; i++) {
      A[k][i]=X[i][k]-Xmean; // M x N
    }
  }
  //cout<<"  -calculated HA, A matrices"<<endl;

  //P=1/(N-1)*HA*(HA)'+1/(N-1)*(eQ*eQ');
  // 1st term is covariance matrix of simulated output states
  // 2nd term is the covariance matrix (R) of ths measurement error
  TransposeMat(HA,Nobs,N,HAT);
  MatMult(HA,HAT,Nobs,N,Nobs,P);        //P is Nobs x Nobs
  TransposeMat(eQT,N,Nobs,eQ);          //eQ is Nobs x N
  MatMult(eQ,eQT,Nobs,N,Nobs,tmp);      //tmp is Nobs x Nobs
  MatAdd(P,tmp,Nobs,Nobs,P);            //P is Nobs x Nobs
  ScalarMatMult(P,1.0/(N-1),Nobs,Nobs,P);

  //Z=HA'*(inv(P)*(O_obs-O_sim));
  double svd_tol=1e-8;
  for(i=0; i<N; i++) 
  {
    for(j=0;j<Nobs;j++) { diff[j]=_obs_matrix[i][j]-_output_matrix[i][j]; }
    SVD(P,diff,ans,Nobs,svd_tol);
    for(j=0;j<Nobs;j++) { MM[j][i]=ans[j]; } //MM is Nobs x N
  }
  MatMult(HAT,MM,N,Nobs,N,Z);           //Z is N x N

  //X_delta = 1/(N-1)*A*Z;
  MatMult(A,Z,M,N,N,X_deltaT);          //X_deltaT is M x N
  ScalarMatMult(X_deltaT,1.0/(N-1),M,N,X_deltaT);
  TransposeMat (X_deltaT,M,N,X_delta);  //X_delta is N x M

  //update state matrix Xa=X+X_delta
  MatAdd(_state_matrix,X_delta,N,M,_state_matrix);   

  /*cout<<"X_delta[][]"<<endl;
  for(int i=0;i<N;i++) {
    for(int j=0;j<M;j++) {
      cout<<X_delta[i][j]<<"| ";
    }
    cout<<endl;
  }*/
  //clean up
  delete[] ans; 
  delete[] diff;
  DeleteMatrix(M,N,A);
  DeleteMatrix(Nobs,N,HA);
  DeleteMatrix(N,Nobs,HAT);
  DeleteMatrix(Nobs,Nobs,P);
  DeleteMatrix(Nobs,N,eQ);
  DeleteMatrix(Nobs,Nobs,tmp);
  DeleteMatrix(Nobs,N,MM);
  DeleteMatrix(N,N,Z);
  DeleteMatrix(M,N,X_deltaT);
  DeleteMatrix(N,M,X_delta);

}
//////////////////////////////////////////////////////////////////
/// \brief updates model - called PRIOR to each model ensemble run 
/// \param pModel [out] pointer to global model instance
/// \param &Options [out] Global model options information
//
void CEnKFEnsemble::UpdateStates(CModel* pModel,optStruct& Options,const int e)
{
  int kk,N,ii;
  ii=0;
  for(int i=0;i<_nAssimStates;i++) 
  {
    kk=_aAssimGroupID[i];
    if(_aAssimStates[i]==STREAMFLOW) 
    {
       //cout<<"UPDATE STATES : # subbasins: "<<pModel->GetSubBasinGroup(kk)->GetNumSubbasins()<<endl;
       for(int pp=0;pp<pModel->GetSubBasinGroup(kk)->GetNumSubbasins();pp++) 
       {

         CSubBasin *pBasin=pModel->GetSubBasinGroup(kk)->GetSubBasin(pp);
         //cout<<"   updating basin "<<pBasin->GetID()<<endl;
         N=pBasin->GetInflowHistorySize();
         const double *U=pBasin->GetRoutingHydrograph();
         const double *Qi=pBasin->GetInflowHistory();
         double *aQi=new double [N];
         double sum=0;
         for(int j=0;j<N;j++) {
           sum+=U[j];
           //aQi[j]=sum*Qi[j]+_state_matrix[e-_nEnKFMembers][ii+j]; //only adjusts portion of flows not already lost from reach (sum*Qi)
           aQi[j]=max(_state_matrix[e-_nEnKFMembers][ii+j],0.0); ii++;
         }
         pBasin->SetQinHist(N,aQi);
         delete [] aQi;

         N=pBasin->GetNumSegments();
         double* aQo=new double[N];
         for(int j=0;j<N;j++) {
           aQo[j]=max(_state_matrix[e-_nEnKFMembers][ii+j],0.0); ii++;
         }
         double QoLast=max(_state_matrix[e-_nEnKFMembers][ii],0.0); ii++;
         pBasin->SetQoutArray(N,aQo,QoLast);
         delete [] aQo;

         ExitGracefullyIf(ii>_nStateVars,"UpdateStates: bad calc",RUNTIME_ERR);
       }
    }
    else {
       for (int k=0;k<pModel->GetHRUGroup(kk)->GetNumHRUs();k++)
       {
          CHydroUnit *pHRU=pModel->GetHRUGroup(kk)->GetHRU(k);
          int iii=pModel->GetStateVarIndex(_aAssimStates[i],_aAssimLayers[0]); 
          pHRU->SetStateVarValue(iii,_state_matrix[e-_nEnKFMembers][ii]); ii++;
       }
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief updates model - called PRIOR to each model ensemble run 
/// \param pModel [out] pointer to global model instance
/// \param &Options [out] Global model options information
//
void CEnKFEnsemble::AddToStateMatrix(CModel* pModel,optStruct& Options,const int e)
{
  //only retrieve states from open loop, not EnKF forecasts
  //This routine should closely echo UpdateStates, except retrieving, not asssigning
  int kk;
  int ii=0;
  int N;
  for(int i=0;i<_nAssimStates;i++)
  {
    kk=_aAssimGroupID[i];
    if(_aAssimStates[i]==STREAMFLOW)
    {
      //cout<<"RETRIEVE STATES : # subbasins: "<<pModel->GetSubBasinGroup(kk)->GetNumSubbasins()<<endl;
      for(int pp=0;pp<pModel->GetSubBasinGroup(kk)->GetNumSubbasins();pp++)
      {
        CSubBasin* pBasin=pModel->GetSubBasinGroup(kk)->GetSubBasin(pp);

        const double* aQi=pBasin->GetInflowHistory();
        const double*   U=pBasin->GetRoutingHydrograph();
        const double*  Qi=pBasin->GetInflowHistory();
        double sum=0;
        for(int j=0;j<pBasin->GetInflowHistorySize();j++) {
          sum+=U[j];
          //_state_matrix[e][ii+j]=(1-sum)*aQi[j]; ii++;
          _state_matrix[e][ii+j]=aQi[j]; ii++;
        }

        N=pBasin->GetNumSegments();
        ExitGracefullyIf(N>1,"EnKF does not currently support nSegments>1 in subbasin reach",BAD_DATA);
        
        _state_matrix[e][ii]=pBasin->GetOutflowRate();      ii++;
        _state_matrix[e][ii]=pBasin->GetLastOutflowRate();  ii++;

        ExitGracefullyIf(ii>_nStateVars,"AddToStateMatrix: bad calc",RUNTIME_ERR);
      }
    }
    else {
      for (int k=0;k<pModel->GetHRUGroup(kk)->GetNumHRUs();k++)
      {
        CHydroUnit *pHRU=pModel->GetHRUGroup(kk)->GetHRU(k);
        int iii=pModel->GetStateVarIndex(_aAssimStates[i],_aAssimLayers[i]); 
        _state_matrix[e][ii]=pHRU->GetStateVarValue(iii); ii++;
      }
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief updates model - called PRIOR to each model ensemble run 
/// \param pModel [out] pointer to global model instance
/// \param &Options [out] Global model options information
//
void CEnKFEnsemble::UpdateModel(CModel *pModel,optStruct &Options,const int e)
{
  string solfile;


  ExitGracefullyIf(e>=_nMembers,"CEnKFEnsemble::UpdateMode: invalid ensemble member index",RUNTIME_ERR);

  //- update output file/ run names ----------------------------
  Options.output_dir=_aOutputDirs[e];
  Options.run_name  =_aRunNames  [e];

  //if :WarmEnsemble, re-read initial conditions from state-adjusted solution_t0.rvc
  // otherwise uses the base model .rvc file for all members
  if(e<_nEnKFMembers) {
    if(_warm_ensemble) {
      if(_warm_runname=="") { solfile="solution_t0.rvc"; }
      else                  { solfile=_warm_runname+"_"+"solution_t0.rvc"; }
      Options.rvc_filename=_aOutputDirs[e+_nEnKFMembers]+solfile;
    }
    ParseInitialConditionsFile(pModel,Options);
  }

  //truncates all hindcasts so they don't run all the way to forecast horizon
  // they must run at least one time step past t0
  if(_truncate_hind) {
    if(e<_nEnKFMembers) { Options.duration=_t_assim_start+ Options.timestep; }
    else                { Options.duration=max(_duration_stored,_t_assim_start + Options.timestep); }
  }
  
  //- read forecast .rvt file if required
  if((e==_nEnKFMembers) && (_forecast_rvt!="")) {
    pModel->ClearTimeSeriesData(Options); //wipes out everything read pre-forecast
    Options.rvt_filename=_forecast_rvt;
    ParseTimeSeriesFile(pModel,Options);
  }

  if(e>=_nEnKFMembers) {
    //For EnKF forecasts, read solution from hindcast ensembles 
    if(Options.run_name=="") { solfile="solution_t0.rvc"; }
    else { solfile=Options.run_name+"_"+"solution_t0.rvc"; }
    Options.rvc_filename=_aOutputDirs[e-_nEnKFMembers]+solfile;
    ParseInitialConditionsFile(pModel,Options);
  }

  if(e>=_nEnKFMembers) {
    //Update state vector in model 
    UpdateStates(pModel,Options,e);
  }
  //The model is run following this routine call...
}
//////////////////////////////////////////////////////////////////
/// \brief called AFTER each model ensemble run 
/// \param pModel [out] pointer to global model instance
/// \param &Options [out] Global model options information
/// \param e [out] ensembe member index
//
void CEnKFEnsemble::FinishEnsembleRun(CModel *pModel,optStruct &Options,const int e)
{
  if(e==_nEnKFMembers-1) //Just prior to first forecast run - generates updated state vector 
  {
    _ENKFOUT<<"PRE-ASSIMILATION STATE MATRIX:"<<endl;
    for(int e=0;e<_nEnKFMembers;e++) {
      for(int i=0;i<_nStateVars;i++) {
        _ENKFOUT<<_state_matrix[e][i]<<",";
      }
      _ENKFOUT<<endl;
    }
    cout<<"PERFORMING ASSIMILATION CALCULATIONS..."<<endl;
    AssimilationCalcs();

    _ENKFOUT<<"POST-ASSIMILATION STATE MATRIX:"<<endl;
    for(int e=0;e<_nEnKFMembers;e++) {
      for(int i=0;i<_nStateVars;i++) {
        _ENKFOUT<<_state_matrix[e][i]<<",";
      }
      _ENKFOUT<<endl;
    }
  }
  if(e==_nMembers-1)
  {
    _ENKFOUT.close();
  }
}
