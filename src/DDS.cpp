/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2023 the Raven Development Team
----------------------------------------------------------------*/
#include "ModelEnsemble.h"

//external function declarations
double UniformRandom();
double GaussRandom();
bool   ParseInitialConditions(CModel *&pModel,const optStruct &Options);

//////////////////////////////////////////////////////////////////
/// \brief DDS Ensemble Construcutor
/// \param num_members [in] number of DDS iterations
/// \param &Options [out] Global model options information
//
CDDSEnsemble::CDDSEnsemble(const int num_members,const optStruct &Options)
  :CEnsemble(num_members,Options)
{
  _type=ENSEMBLE_DDS;
  _nParamDists=0;
  _pParamDists=NULL;
  _BestParams=NULL;
  _TestParams=NULL;
  _Fbest=ALMOST_INF;
  _r_val=0.2;

  _calib_SBID=DOESNT_EXIST;
  _calib_Obj=DIAG_NASH_SUTCLIFFE;
  _calib_Period="ALL";
}

//////////////////////////////////////////////////////////////////
/// \brief DDS Ensemble Destrucutor
//
CDDSEnsemble::~CDDSEnsemble()
{
  for(int i=0;i<_nParamDists;i++) {
    delete _pParamDists[i];
  }
  delete[] _pParamDists; _nParamDists=0;
  delete[] _BestParams;
  delete[] _TestParams;
}

//////////////////////////////////////////////////////////////////
/// \brief Adds parameter distribution to DDS setup
/// \param dist [in] pointer to parameter distribution structure for specified parameter
//
void CDDSEnsemble::AddParamDist(const param_dist *dist)
{
  if(!DynArrayAppend((void**&)(_pParamDists),(void*)(dist),_nParamDists)) {
    ExitGracefully("CDDSEnsemble::AddParamDist: adding NULL distribution",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief set DDS perturbation value
//
void CDDSEnsemble::SetPerturbationValue(const double &perturb) {
  _r_val=perturb;
  if((_r_val<0) || (_r_val>1.0)) {
    ExitGracefully("CDDSEnsemble::SetPerturbationValue",BAD_DATA_WARN);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief set DDS calibration target
//
void CDDSEnsemble::SetCalibrationTarget(const long long SBID,const diag_type object_diag,const string period)
{
  _calib_SBID  =SBID;
  _calib_Obj   =object_diag;
  _calib_Period=period;
}

//////////////////////////////////////////////////////////////////
/// \brief initializes DDS caliobration run
/// \param &Options [out] Global model options information
//
void CDDSEnsemble::Initialize(const CModel* pModel,const optStruct &Options)
{
  // QA/QC
  //-----------------------------------------------
  if(_calib_SBID==DOESNT_EXIST) {
    ExitGracefully("CDDSEnsemble::Initialize: DDS calibration target/objective function must be set",BAD_DATA);
  }
  if(_nMembers==0) {
    ExitGracefully("CDDSEnsemble::Initialize: number of DDS iterations/ensemble members must be >0",BAD_DATA);
  }
  if(_nParamDists==0) {
    ExitGracefully("CDDSEnsemble::Initialize: at least one parameter distribution must be set using :ParameterDistributions command",BAD_DATA);
  }

  // Initialize best and test parames to default
  //-----------------------------------------------
  _BestParams=new double[_nParamDists];
  _TestParams=new double[_nParamDists];
  for(int i=0; i<_nParamDists;i++)
  {
    _BestParams[i]=_TestParams[i]=_pParamDists[i]->default_val;
  }

  // Create and open DDSOutput file
  //-----------------------------------------------
  string filename=Options.main_output_dir+"DDSOutput.csv";
  _DDSOUT.open(filename.c_str());
}

/**********************************************************************
 PerturbParam
-----------------------------------------------------------------------
Generates a neighboring decision variable value for a single
decision variable value being perturbed by the DDS optimization algorithm.
New DV value respects the upper and lower DV bounds.

Translated by James Craig July 2006 from Bryan Tolson's Fortran DDS code

returns  new decision variable value (within specified min and max)
**********************************************************************/
double CDDSEnsemble::PerturbParam(const double &x_best, //current best decision variable (DV) value
                                  const double &lowerbound,
                                  const double &upperbound)
{
  double x_new;
  double x_min = lowerbound;
  double x_max = upperbound;

  x_new=x_best+GaussRandom()*_r_val*(x_max-x_min);

  // need if statements to check within DV bounds.  If not, bounds are reflecting.
  if(x_new<x_min)
  {
    x_new=x_min+(x_min-x_new); //reflect
    if(x_new>x_max) { x_new=x_min; }
  }
  else if(x_new>x_max)
  {
    x_new=x_max-(x_new-x_max);//reflect
    if(x_new<x_min) { x_new=x_max; }
  }
  return x_new;
}

//////////////////////////////////////////////////////////////////
/// \brief updates model - called PRIOR to each model ensemble run
/// \param pModel [out] pointer to global model instance
/// \param &Options [out] Global model options information
//
void CDDSEnsemble::UpdateModel(CModel *pModel,optStruct &Options,const int e)
{
  CEnsemble::UpdateModel(pModel,Options,e);
  ExitGracefullyIf(e>=_nMembers,"CDDSEnsemble::UpdateMode: invalid ensemble member index",RUNTIME_ERR);

  //int iters_remaining=_nMembers-e;
  double u;

  //- update output file/ run names ----------------------------
  Options.output_dir=_aOutputDirs[e];
  Options.run_name  =_aRunNames[e];

  //- Update parameter values ----------------------------------
  // Determine variable selected as neighbour
  double Pn=1.0-log(double(e))/log(double(_nMembers));
  int dvn_count=0;

  //- define TestParams initially as best current solution------
  for(int k=0;k<_nParamDists;k++)
  {
    _TestParams[k]=_BestParams[k];
  }

  //- perturb TestParams --------------------------------------
  for(int k=0;k<_nParamDists;k++)
  {
    u=UniformRandom();
    if(u<Pn) {
      dvn_count++;
      _TestParams[k]=PerturbParam(_BestParams[k],_pParamDists[k]->distpar[0],_pParamDists[k]->distpar[1]);
    }
  }
  if(dvn_count==0) {
    u=UniformRandom();
    int dv=(int)(ceil((double)(_nParamDists)*u))-1; // index for one DV
    _TestParams[dv]=PerturbParam(_BestParams[dv],_pParamDists[dv]->distpar[0],_pParamDists[dv]->distpar[1]);
  }

  //- update parameters in model -----------------------------
  for(int k=0;k<_nParamDists;k++)
  {
    pModel->UpdateParameter(_pParamDists[k]->param_class,
                            _pParamDists[k]->param_name,
                            _pParamDists[k]->class_group,
                            _TestParams[k]);
  }
  //- Re-read initial conditions to update state variables----
  if(!ParseInitialConditions(pModel,Options)) {
    ExitGracefully("Cannot find or read .rvc file",BAD_DATA);
  }
  pModel->CalculateInitialWaterStorage(Options);

  //The model is run following this routine call...
}
//////////////////////////////////////////////////////////////////
/// \brief called AFTER each model ensemble run
/// \param pModel [out] pointer to global model instance
/// \param &Options [out] Global model options information
/// \param e [out] ensembe member index
//
void CDDSEnsemble::FinishEnsembleRun(CModel *pModel,optStruct &Options,const time_struct &tt,const int e)
{
  double Ftest=pModel->GetObjFuncVal(_calib_SBID,_calib_Obj,_calib_Period);

  // update current (best) solution - optimization is minimization
  //----------------------------------------------
  if(Ftest<=_Fbest)
  {
    _Fbest = Ftest;
    for(int k=0;k<_nParamDists;k++) { _BestParams[k]=_TestParams[k]; }

    //write results
    _DDSOUT<<e+1<<", "<<_Fbest<<","<<endl;
  }

  //Write objective function to screen
  //----------------------------------------------
  cout<<"DDS Obj. Function: "<<Ftest <<" [best: "<<_Fbest<<"]"<<endl;
  //for(int i=0;i<_nParamDists;i++) {cout<<"P["<<i<<"]: "<<_TestParams[i]<<", ";}cout<<endl;

  // write best parameter vector
  //----------------------------------------------
  if(e==_nMembers-1)
  {
    _DDSOUT<<"Best parameter vector:"<<endl;
    for(int k=0;k<_nParamDists;k++) {_DDSOUT<<_pParamDists[k]->param_name<<"("<<_pParamDists[k]->class_group <<"), "<<_BestParams[k]<<endl; }
    _DDSOUT.close();
  }
}
