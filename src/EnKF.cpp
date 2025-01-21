/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2024 the Raven Development Team
----------------------------------------------------------------*/

#include "EnKF.h"
#include "Matrix.h"

bool IsContinuousFlowObs(const CTimeSeriesABC* pObs,long SBID);
bool ParseInitialConditions(CModel*& pModel,const optStruct& Options);
bool ParseTimeSeriesFile(CModel*& pModel,const optStruct& Options);

// in .rvi file:
//  :EnsembleMode ENSEMBLE_ENKF [double #members]
//  :WarmEnsemble fc2 # used if ICs should be read from ensemble fc2_solution_enKF.rvc file
//  :SilentMode # useful
//
// in .rve file:
//  :EnKFMode ENKF_CLOSED_LOOP or ENKF_SPINUP or ENKF_CLOSED_FORECAST...
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
  :CEnsemble(num_members,Options)
{
  _type        =ENSEMBLE_ENKF;
  _EnKF_mode   =ENKF_UNSPECIFIED;
  _nEnKFMembers=num_members;

  _pObsPerturbations=NULL;
  _nObsPerturbations=0;
  _aAssimStates =NULL;
  _aAssimLayers =NULL;
  _aAssimGroupID=NULL;
  _nAssimStates=0;

  _state_matrix =NULL; //Built in ::Initialize
  _nStateVars=0;
  _state_names=NULL;

  _warm_runname="";
  _extra_rvt="";
  _orig_rvc_file=Options.rvc_filename;

  _obs_matrix   =NULL;
  _output_matrix=NULL;
  _noise_matrix =NULL;
  _nObsDatapoints=0;

  _window_size=1;
  _nTimeSteps =0;
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
  delete [] _state_names;

  for (int i=0;i<_nObsPerturbations;i++){delete _pObsPerturbations[i];} delete [] _pObsPerturbations;

  delete [] _aAssimStates;
  delete [] _aAssimLayers;
  delete [] _aAssimGroupID;
  delete [] _aObsIndices;
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
  pOP->state       =type;
  pOP->distribution=distrib;
  pOP->kk          =DOESNT_EXIST;
  //pOP->kk        =group_index;
  pOP->adj_type    =adj;
  for (int i = 0; i < 3; i++) {
    pOP->distpar[i]=distpars[i];
  }
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
    tmp [i]=_aAssimStates [i];
    tmp2[i]=_aAssimGroupID[i];
    tmp3[i]=_aAssimLayers [i];
  }
  delete [] _aAssimStates;
  delete [] _aAssimGroupID;
  delete [] _aAssimLayers;
  _aAssimStates =tmp;
  _aAssimGroupID=tmp2;
  _aAssimLayers =tmp3;
  _aAssimStates [_nAssimStates-1]=sv;
  _aAssimGroupID[_nAssimStates-1]=assim_groupID;
  _aAssimLayers [_nAssimStates-1]=lay;
}
//////////////////////////////////////////////////////////////////
/// \brief sets runnname of solution file used in open- or closed-loop simulations
//
void CEnKFEnsemble::SetWarmRunname(string runname)
{
  _warm_runname=runname;
}
//////////////////////////////////////////////////////////////////
/// \brief sets EnKF mode - spinup/open loop/closed loop
//
void CEnKFEnsemble::SetEnKFMode(EnKF_mode mode)
{
  _EnKF_mode=mode;
}

//////////////////////////////////////////////////////////////////
/// \brief set value of data horizon
//
void CEnKFEnsemble::SetWindowSize(const int nTimesteps){ _window_size=nTimesteps;}

//////////////////////////////////////////////////////////////////
/// \brief set extra RVT filename
//
void CEnKFEnsemble::SetExtraRVTFile(string filename){_extra_rvt=filename;}

//////////////////////////////////////////////////////////////////
/// \brief Accessor - gets ensemble start time (in local model time)
/// \param e [in] ensemble index
/// \return ensemble start time (in local model time)
//
double CEnKFEnsemble::GetStartTime(const int e) const
{
  return 0.0;
}

//////////////////////////////////////////////////////////////////
/// \brief Accessor - gets ensemble mode
/// \return ensemble mode
//
EnKF_mode CEnKFEnsemble::GetEnKFMode() const
{
  return _EnKF_mode;
}

//////////////////////////////////////////////////////////////////
/// \brief initializes EnKF assimilation run
/// \param &Options [out] Global model options information
//
void CEnKFEnsemble::Initialize(const CModel* pModel,const optStruct &Options)
{
  CEnsemble::Initialize(pModel,Options);
  CStateVariable *pStateVar = pModel->GetStateVarInfo();

  // QA/QC
  //-----------------------------------------------
  if(_nMembers==0) {
    ExitGracefully("CEnKFEnsemble::Initialize: number of EnKF ensemble members must be >0",BAD_DATA);
  }
  if ((pModel->GetNumForcingPerturbations() == 0) && (_EnKF_mode != ENKF_FORECAST) && (_EnKF_mode != ENKF_OPEN_FORECAST)) {
    ExitGracefully("CEnKFEnsemble::Initialize: at least one forcing perturbation must be set using :ForcingPerturbation command",BAD_DATA);
  }
  if ((_nAssimStates==0) && (_EnKF_mode!=ENKF_FORECAST) && (_EnKF_mode!=ENKF_OPEN_FORECAST)) {
    ExitGracefully("CEnKFEnsemble::Initialize: at least one assimilated state must be set in :ForcingPerturbation command",BAD_DATA);
  }
  if (Options.assimilate_flow == true) {
    WriteWarning("CEnKFEnsemble::Initialize: It is not recommended to use EnKF data assimilation simultaneously with direct insertion streamflow assimilation (enabled with the :AssimilateStreamflow command in the .rvi file). The two assimilation strategies may interfere with each other, particularly when the assimilated observations are shared by both.",Options.noisy);
  }
  if (Options.assimilate_stage == true) {
    WriteWarning("CEnKFEnsemble::Initialize: It is not recommended to use EnKF data assimilation simultaneously with direct insertion reservoir stage assimilation (enabled with the :AssimilateReservoirStage command in the .rvi file). The two assimilation strategies may interfere with each other, particularly when the assimilated observations are shared by both.",Options.noisy);
  }

  //TO ADD - complain if .rvc file is anything other than runname_solution.rvc

  // determine total number of state variables
  // should closely echo contents of AddToStateMatrix/UpdateFromStateMatrix
  //-----------------------------------------------
  int kk;

  _nStateVars=0;
  int ii=0;
  for(int i=0;i<_nAssimStates;i++)
  {
    kk=_aAssimGroupID[i];
    if(_aAssimStates[i]==STREAMFLOW) {
      for(int pp=0;pp<pModel->GetSubBasinGroup(kk)->GetNumSubbasins();pp++)
      {
        CSubBasin* pBasin=pModel->GetSubBasinGroup(kk)->GetSubBasin(pp);
        _nStateVars+=pBasin->GetInflowHistorySize();
        _nStateVars+=pBasin->GetNumSegments();
        _nStateVars++; //Qlast

        if (pBasin->GetReservoir() != NULL){_nStateVars+=2; }//resflow + resflow_last
      }
    }
    else if(_aAssimStates[i]==RESERVOIR_STAGE) {
      for(int pp=0;pp<pModel->GetSubBasinGroup(kk)->GetNumSubbasins();pp++)
      {
        CSubBasin* pBasin=pModel->GetSubBasinGroup(kk)->GetSubBasin(pp);
        if (pBasin->GetReservoir() != NULL){_nStateVars+=2; }
      }
    }
    else { //HRU state variable
      _nStateVars+=pModel->GetHRUGroup(kk)->GetNumHRUs();
    }
  }

  // get names of state variables
  //-----------------------------------------------
  _state_names=new string [ _nStateVars];
  for(int i=0;i<_nAssimStates;i++)
  {
    kk=_aAssimGroupID[i];
    if(_aAssimStates[i]==STREAMFLOW) {
      for(int pp=0;pp<pModel->GetSubBasinGroup(kk)->GetNumSubbasins();pp++)
      {
        CSubBasin* pBasin=pModel->GetSubBasinGroup(kk)->GetSubBasin(pp);

        for (int n=0;n<pBasin->GetInflowHistorySize();n++){_state_names[ii]="inflow_" +to_string(pp)+"_"+to_string(n); ii++; }
        for (int n=0;n<pBasin->GetNumSegments();      n++){_state_names[ii]="outflow_"+to_string(pp)+"_"+to_string(n); ii++; }
        _state_names[ii]="outflowlast_"+to_string(pp); ii++;

        if (pBasin->GetReservoir() != NULL)
        {
          _state_names[ii]="resflow_"    +to_string(pp); ii++;
          _state_names[ii]="resflowlast_"+to_string(pp); ii++;
        }
      }
    }
    else if(_aAssimStates[i]==RESERVOIR_STAGE) {
      for(int pp=0;pp<pModel->GetSubBasinGroup(kk)->GetNumSubbasins();pp++)
      {
        CSubBasin* pBasin=pModel->GetSubBasinGroup(kk)->GetSubBasin(pp);
        if (pBasin->GetReservoir() != NULL)
        {
          _state_names[ii]="resstage_"    +to_string(pp); ii++;
          _state_names[ii]="resstagelast_"+to_string(pp); ii++;
        }
      }
    }
    else { //HRU state variable
      for (int n=0;n<pModel->GetHRUGroup(kk)->GetNumHRUs();n++)
      {
        long long int k=pModel->GetHRUGroup(kk)->GetHRU(n)->GetHRUID();
        string svname = pModel->GetStateVarInfo()->SVTypeToString(_aAssimStates[i], _aAssimLayers[i]);
        _state_names[ii]=svname+"_" +to_string(k); ii++;
      }
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
  cout<<"ENKF: Found "<<_nStateVars<<" state variables for assimilation. State matrix built."<<endl<<endl;

  //determine total number of observations
  //-----------------------------------------------
  _nTimeSteps=(int)(rvn_floor((Options.duration+TIME_CORRECTION)/Options.timestep));
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
  if ((_nObs==0) && (_EnKF_mode!=ENKF_FORECAST) && (_EnKF_mode!=ENKF_OPEN_FORECAST)) {
    ExitGracefully("EnKF::Initialize : no observations available for assimilation. Use :AssimilateStreamflow command to indicate enabled observation time series",BAD_DATA);
  }

  if ((_nObsDatapoints==0) && (_EnKF_mode!=ENKF_FORECAST) && (_EnKF_mode!=ENKF_OPEN_FORECAST)) {//DELTARES_REV - if no obs found, just skip assimilation
    WriteWarning("EnKF::Initialize : no observations were available for assimilation. EnKF assimilation updates will not be performed.",Options.noisy);
  }
  else if (_nObsDatapoints>0) {
    cout<<"ENKF: Found "<<_nObsDatapoints<<" valid observation data points for assimilation. Observation matrices built."<<endl;
  }

  //allocate _output_matrix and _noise_matrix
  //-----------------------------------------------
  _obs_matrix   =new double *[_nEnKFMembers];
  _output_matrix=new double *[_nEnKFMembers];
  _noise_matrix =new double *[_nEnKFMembers];
  for(int e=0;e<_nEnKFMembers;e++) {
    _noise_matrix [e]=NULL;
    _obs_matrix   [e]=new double [_nObsDatapoints];
    _output_matrix[e]=new double [_nObsDatapoints];
    _noise_matrix [e]=new double [_nObsDatapoints];
    ExitGracefullyIf(_noise_matrix[e]==NULL,"EnKF Initialization: noise matrix",OUT_OF_MEMORY);
  }

  //populate _output_matrix, _obs_matrix
  //-----------------------------------------------
  int j;
  double eps;
  for(int e=0;e<_nEnKFMembers;e++)
  {
    j=0;
    for(ii=0;ii<_nObs;ii++)
    {
      const CTimeSeriesABC* pTSObs=pModel->GetObservedTS(_aObsIndices[ii]);

      //get perturbation info from observation time series name. if not found, obs data is not perturbed
      sv_type sv; int lay;
      obs_perturb *pPerturb=NULL;
      if      (!strcmp(pTSObs->GetName().c_str(),"HYDROGRAPH")     ){ sv = STREAMFLOW;}
      else if (!strcmp(pTSObs->GetName().c_str(),"RESERVOIR_STAGE")){ sv = RESERVOIR_STAGE;}
      else                                                          { sv = pStateVar->StringToSVType(pTSObs->GetName(),lay,true);}

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
            if      (pPerturb->adj_type == ADJ_ADDITIVE      ){ _noise_matrix[e][j]=eps;}
            else if (pPerturb->adj_type == ADJ_MULTIPLICATIVE){ _noise_matrix[e][j]=(eps*obsval)-obsval; }
          }
          _obs_matrix[e][j]=obsval+_noise_matrix[e][j];
          j++;
        }
      }
    }

    for(j=0;j<_nObsDatapoints;j++) {
       _output_matrix[e][j]=0; //filled during simulation in CloseTimeStepOps
    }
  }

  // Create and open EnKFOutput file
  //-----------------------------------------------
  string filename= FilenamePrepare("EnKFOutput.csv",Options);
  _ENKFOUT.open(filename.c_str());
  if(_ENKFOUT.fail()) {
    ExitGracefully(("CEnKFEnsemble::Initialize: Unable to open output file "+filename+" for writing.").c_str(),FILE_OPEN_ERR);
  }

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
  //Apply forcing perturbations for the current time step
  //----------------------------------------------------------
  if ((_EnKF_mode!=ENKF_FORECAST) && (_EnKF_mode!=ENKF_OPEN_FORECAST))
  {
    pModel->PrepareForcingPerturbation(Options, tt);
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
  if (fabs(tt.model_time-Options.duration)<0.25*Options.timestep) { //t=t_end
    AddToStateMatrix(pModel,Options,e);
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
          _output_matrix[e][j]=pTSMod->GetSampledValue(nn);
          //cout<<"GRABBING OUTPUT ens: "<<e<<" nn="<<curr_nn<<" obs: "<<_output_matrix[e][j]<<" "<<obsval<<" duration "<<Options.duration<<endl;
          break;
        }
        j++;
      }
    }
  }
  if (j!=_nObsDatapoints){
    ExitGracefully("Indexing issue in EnKF::CloseTimeStepOps: This shouldnt happen.",RUNTIME_ERR);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief the EnKF assimilation matrix calculations for updating states
/// Based upon Mandel, J., Efficient Implementation of the Ensemble Kalman Filter, Report, Univ of Colorado, 2006.
//
void CEnKFEnsemble::AssimilationCalcs()
{
  if (_nObsDatapoints==0){return; } //skips assimilation if no observations available

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

  /*cout << "_output_matrix[][]:" << endl;
  for(int i=0;i<N;i++) {
    for(int j=0;j<Nobs;j++) {cout<<_output_matrix[i][j]<<"| ";}cout<<endl;
  }
  cout<<"_obs_matrix[][]:"<<endl;
  for(int i=0;i<N;i++) {
    for(int j=0;j<Nobs;j++) {cout<<_obs_matrix[i][j]<<"| ";}cout<<endl;
  } */


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

  /*cout << "X_delta[][]" << endl;
  for(int i=0;i<N;i++) {
    for(int j=0;j<M;j++) {cout<<X_delta[i][j]<<"| ";}cout<<endl;
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
/// \brief updates model states - called in FinishEnsembleRun right after assimilation
/// \param pModel [out] pointer to global model instance
/// \param &Options [out] Global model options information
//
void CEnKFEnsemble::UpdateFromStateMatrix(CModel* pModel,optStruct& Options,const int e)
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

         N=pBasin->GetInflowHistorySize();
         double *aQi=new double [N];
         for(int j=0;j<N;j++) {
           aQi[j]=max(_state_matrix[e][ii+j],0.0);
         }
         ii+=N;
         pBasin->SetQinHist(N,aQi);

         N=pBasin->GetNumSegments();
         double* aQo=new double[N];
         for(int j=0;j<N;j++) {
           aQo[j]=max(_state_matrix[e][ii+j],0.0);
         }
         ii+=N;

         double QoLast=max(_state_matrix[e][ii],0.0);
         ii++;

         pBasin->SetQoutArray(N,aQo,QoLast);

         if (pBasin->GetReservoir() != NULL)
         {
           time_struct tt_dummy;
           pBasin->GetReservoir()->SetInitialFlow(_state_matrix[e][ii],_state_matrix[e][ii+1],tt_dummy,Options);
           ii+=2;
         }

         delete [] aQi;
         delete [] aQo;

         ExitGracefullyIf(ii>_nStateVars,"UpdateFromStateMatrix: bad calc",RUNTIME_ERR);
       }
    }
    else if(_aAssimStates[i]==RESERVOIR_STAGE)
    {
      for(int pp=0;pp<pModel->GetSubBasinGroup(kk)->GetNumSubbasins();pp++)
      {
         CSubBasin *pBasin=pModel->GetSubBasinGroup(kk)->GetSubBasin(pp);
         if (pBasin->GetReservoir() != NULL)
         {
           pBasin->GetReservoir()->SetReservoirStage(_state_matrix[e][ii],_state_matrix[e][ii+1]);
           ii+=2;
         }
       }
    }
    else //Assimilated states are HRU state variables
    {
       for (int k=0;k<pModel->GetHRUGroup(kk)->GetNumHRUs();k++)
       {
          CHydroUnit *pHRU=pModel->GetHRUGroup(kk)->GetHRU(k);
          int iii=pModel->GetStateVarIndex(_aAssimStates[i],_aAssimLayers[i]);
          pHRU->SetStateVarValue(iii,_state_matrix[e][ii]); ii++;
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
  //This routine should closely echo UpdateFromStateMatrix, except retrieving, not asssigning
  int kk,N;
  int ii=0;
  for(int i=0;i<_nAssimStates;i++)
  {
    kk=_aAssimGroupID[i];
    if(_aAssimStates[i]==STREAMFLOW)
    {
      //cout<<"RETRIEVE STATES : # subbasins: "<<pModel->GetSubBasinGroup(kk)->GetNumSubbasins()<<endl;
      for(int pp=0;pp<pModel->GetSubBasinGroup(kk)->GetNumSubbasins();pp++)
      {
        CSubBasin* pBasin=pModel->GetSubBasinGroup(kk)->GetSubBasin(pp);

        N=pBasin->GetInflowHistorySize();
        const double* aQi=pBasin->GetInflowHistory();
        for(int j=0;j<N;j++) {
          _state_matrix[e][ii+j]=aQi[j];
        }
        ii+=N;

        N=pBasin->GetNumSegments();
        const double* aQo=pBasin->GetOutflowArray();
        for(int j=0;j<N;j++) {
          _state_matrix[e][ii+j]=aQo[j];
        }
        ii+=N;

        _state_matrix[e][ii]=pBasin->GetLastOutflowRate();
        ii++;

        if (pBasin->GetReservoir() != NULL)
        {
          _state_matrix[e][ii]=pBasin->GetReservoir()->GetOutflowRate();
          ii++;
          _state_matrix[e][ii]=pBasin->GetReservoir()->GetOldOutflowRate();
          ii++;
        }
        ExitGracefullyIf(ii>_nStateVars,"AddToStateMatrix: bad calc",RUNTIME_ERR);
      }
    }
    else if(_aAssimStates[i]==RESERVOIR_STAGE)
    {
      for(int pp=0;pp<pModel->GetSubBasinGroup(kk)->GetNumSubbasins();pp++)
      {
        CSubBasin* pBasin=pModel->GetSubBasinGroup(kk)->GetSubBasin(pp);
        if (pBasin->GetReservoir() != NULL)
        {
          _state_matrix[e][ii]=pBasin->GetReservoir()->GetResStage();
          ii++;
          _state_matrix[e][ii]=pBasin->GetReservoir()->GetOldStage();
          ii++;
        }
      }
    }
    else //Assimilated states are HRU state variables
    {
      for (int k=0;k<pModel->GetHRUGroup(kk)->GetNumHRUs();k++)
      {
        CHydroUnit *pHRU=pModel->GetHRUGroup(kk)->GetHRU(k);
        int iii=pModel->GetStateVarIndex(_aAssimStates[i],_aAssimLayers[i]);
        _state_matrix[e][ii]=pHRU->GetStateVarValue(iii);
        ii++;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief updates model - called PRIOR to t=0 for each model ensemble member
/// \param pModel [out] pointer to global model instance
/// \param &Options [out] Global model options information
//
void CEnKFEnsemble::UpdateModel(CModel *pModel,optStruct &Options,const int e)
{
  CEnsemble::UpdateModel(pModel,Options,e);
  string solfile;

  ExitGracefullyIf(e>=_nMembers,"CEnKFEnsemble::UpdateMode: invalid ensemble member index",RUNTIME_ERR);

  //- update output file/ run names ----------------------------
  Options.output_dir=_aOutputDirs[e];
  Options.run_name  =_aRunNames  [e]; //typically fixed for EnKF

  //if closed-loop, re-read initial conditions from state-adjusted solution_EnKF.rvc
  //if open-loop, read initial conditions from unadjusted solution.rvc files
  // otherwise uses the single base model .rvc file for all members

  if((_EnKF_mode == ENKF_CLOSED_LOOP) || (_EnKF_mode == ENKF_FORECAST)) {
    if(_warm_runname=="") { solfile="solution_EnKF.rvc"; }
    else                  { solfile=_warm_runname+"_"+"solution_EnKF.rvc"; }
    Options.rvc_filename=_aOutputDirs[e]+solfile;
    if (_aSolutionFiles[e] != ""){Options.rvc_filename=_aSolutionFiles[e]; } //Allows storage elsewhere -requires user to know whether these are EnKF or not
  }
  else if ((_EnKF_mode==ENKF_OPEN_LOOP) || (_EnKF_mode == ENKF_OPEN_FORECAST)) {
    if(_warm_runname=="") { solfile="solution.rvc"; }
    else                  { solfile=_warm_runname+"_"+"solution.rvc"; }
    Options.rvc_filename=_aOutputDirs[e]+solfile;
    if (_aSolutionFiles[e] != ""){Options.rvc_filename=_aSolutionFiles[e]; } //Allows storage elsewhere -requires user to know whether these are EnKF or not
  }
  else {//if(_EnKF_mode==ENKF_SPINUP)
    Options.rvc_filename=_orig_rvc_file;
  }

  ParseInitialConditions(pModel,Options);
  pModel->CalculateInitialWaterStorage(Options);

  // read ensemble-member specific time series (e.g., upstream flows in model cascade), if present
  if (_extra_rvt != "") {
    string rvt_name=_extra_rvt;
    SubstringReplace(rvt_name,"*",to_string(e+1));//replaces wildcard
    Options.rvt_filename=rvt_name;
    ParseTimeSeriesFile(pModel,Options);
  }

  //The model is run following this routine call...
}

//////////////////////////////////////////////////////////////////
/// \brief called AFTER each model ensemble run
/// \param pModel [out] pointer to global model instance
/// \param &Options [out] Global model options information
/// \param e [out] ensembe member index
//
void CEnKFEnsemble::FinishEnsembleRun(CModel *pModel,optStruct &Options,const time_struct &tt,const int e)
{
  if ((_EnKF_mode==ENKF_OPEN_LOOP) || (_EnKF_mode==ENKF_FORECAST) || (_EnKF_mode==ENKF_OPEN_FORECAST)){
    return; //EnKF output file not populated, solution_EnKF.rvc file not generated
  }

  //grabs states and stores them in state matrix
  AddToStateMatrix(pModel,Options,e);

  if(e==_nEnKFMembers-1) //After all ensemble members have run
  {
    _ENKFOUT<<"PRE-ASSIMILATION STATE MATRIX:"<<endl;
    _ENKFOUT<<"member"<<",";
    for (int i = 0; i<_nStateVars; i++) {
      _ENKFOUT<<_state_names[i]<<",";
    }
    _ENKFOUT<<endl;
    for(int ee=0;ee<_nEnKFMembers;ee++) {
      _ENKFOUT<<ee+1<<",";
      for(int i=0;i<_nStateVars;i++) {
        _ENKFOUT<<_state_matrix[ee][i]<<",";
      }
      _ENKFOUT<<endl;
    }

    cout<<endl<<"ENKF: Performing Assimilation Calculations with "<<_nObsDatapoints<<" observation datapoints..."<<endl;
    AssimilationCalcs();

    _ENKFOUT<<"POST-ASSIMILATION STATE MATRIX:"<<endl;
    _ENKFOUT<<"member"<<",";
    for (int i = 0; i<_nStateVars; i++) {
      _ENKFOUT<<_state_names[i]<<",";
    }
    _ENKFOUT<<endl;
    for(int ee=0;ee<_nEnKFMembers;ee++) {
      _ENKFOUT<<ee+1<<",";
      for(int i=0;i<_nStateVars;i++) {
        _ENKFOUT<<_state_matrix[ee][i]<<",";
      }
      _ENKFOUT<<endl;
    }

    _ENKFOUT<<"NOISE MATRIX:"<<endl;
    for(int ee=0;ee<_nEnKFMembers;ee++) {
      _ENKFOUT<<ee+1<<",";
      for(int i=0;i<_nObsDatapoints;i++) {
        _ENKFOUT<<_noise_matrix[ee][i]<<",";
      }
      _ENKFOUT<<endl;
    }
    _ENKFOUT<<"OBSERVATION MATRIX:"<<endl;
    for(int ee=0;ee<_nEnKFMembers;ee++) {
      _ENKFOUT<<ee+1<<",";
      for(int i=0;i<_nObsDatapoints;i++) {
        _ENKFOUT<<_obs_matrix[ee][i]<<",";
      }
      _ENKFOUT<<endl;
    }
    _ENKFOUT<<"SIMULATED OUTPUT MATRIX:"<<endl;
    for(int ee=0;ee<_nEnKFMembers;ee++) {
      _ENKFOUT<<ee+1<<",";
      for(int i=0;i<_nObsDatapoints;i++) {
        _ENKFOUT<<_output_matrix[ee][i]<<",";
      }
      _ENKFOUT<<endl;
    }
    _ENKFOUT.close();

    //Write EnKF-updated solution files
    // requires full re-reading of ensemble member solution files
    string solfile;
    cout<<"ENKF: Writing  Solution Files..."<<endl;
    for(int ee=0;ee<_nEnKFMembers;ee++)
    {
      Options.output_dir=_aOutputDirs[ee];
      Options.run_name  =_aRunNames  [ee];

      if(Options.run_name=="") { solfile="solution.rvc"; }
      else                     { solfile=Options.run_name+"_"+"solution.rvc"; }
      Options.rvc_filename=_aOutputDirs[ee]+solfile;

      ParseInitialConditions(pModel,Options);

      //Update state vector in model
      UpdateFromStateMatrix(pModel,Options,ee);

      pModel->WriteMajorOutput(Options,tt,"solution_EnKF",false);
    }
  }
}
