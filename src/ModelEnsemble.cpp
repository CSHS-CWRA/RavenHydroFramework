/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2019 the Raven Development Team
----------------------------------------------------------------*/
#include "ModelEnsemble.h"
//#include <random>
//see http://anadoxin.org/blog/c-shooting-yourself-in-the-foot-4.html

bool ParseInitialConditionsFile(CModel *&pModel,const optStruct &Options);

double UniformRandom(){return rand()/RAND_MAX; }
double GaussRandom(){return 1.0;}//TMP DEBUG

//////////////////////////////////////////////////////////////////
/// \brief Ensemble Default Constructor
///
/// \param num_members [in] number of ensemble members/realizations
/// \param &Options [out] Global model options information
//
CEnsemble::CEnsemble(const int num_members,const optStruct &Options)
{
  _type=ENSEMBLE_NONE;
  _nMembers=num_members;

  _aOutputDirs=NULL;
  _aOutputDirs=new string [_nMembers];
  ExitGracefullyIf(_aOutputDirs==NULL,"CEnsemble constructor",OUT_OF_MEMORY);

  _aRunNames  =NULL;
  _aRunNames=new string[_nMembers];
  ExitGracefullyIf(_aRunNames==NULL,"CEnsemble constructor (2)",OUT_OF_MEMORY);

  for(int e=0;e<_nMembers;e++) 
  {
    _aOutputDirs[e]=Options.main_output_dir;
    _aRunNames  [e]=Options.run_name;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Ensemble Default Destructor
//
CEnsemble::~CEnsemble()
{
  delete [] _aOutputDirs;
  delete [] _aRunNames;
}
//////////////////////////////////////////////////////////////////
/// \brief Accessor - gets number of ensemble members
/// \return number of ensemble members
//
int CEnsemble::GetNumMembers() {
  return _nMembers;
}
//////////////////////////////////////////////////////////////////
/// \brief Accessor - gets ensemble type
/// \return ensemble type
//
ensemble_type CEnsemble::GetType() {
  return _type;
}

//Manipulator Functions
//////////////////////////////////////////////////////////////////
/// \brief sets random seed
/// \param [in] random seed
//
void CEnsemble::SetRandomSeed(const int seed) {

}
//////////////////////////////////////////////////////////////////
/// \brief sets output directory for ensemble member output
/// \param OutDirString [in] string with or without '*' random card. If * is present, will be replaced with ensemble ID
//
void CEnsemble::SetOutputDirectory(const string OutDirString) 
{
  size_t pos=OutDirString.find("*");

  for(int e=0;e<_nMembers;e++)
  {
    _aOutputDirs[e]=OutDirString;
    if(pos!=string::npos){ _aOutputDirs[e].replace(pos,1,to_string(e+1)); }
    _aOutputDirs[e]=_aOutputDirs[e]+"/"; //add trailing backslash
  }
}
//////////////////////////////////////////////////////////////////
/// \brief sets run name for ensemble members
/// \param RunNameString [in] string with or without '*' random card. If * is present, will be replaced with ensemble ID
//
void CEnsemble::SetRunNames(const string RunNameString) 
{
  size_t pos=RunNameString.find("*");

  for(int e=0;e<_nMembers;e++)
  {
    _aRunNames[e]=RunNameString;
    if(pos!=string::npos) { _aRunNames[e].replace(pos,1,to_string(e+1)); }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief initializes ensemble
/// \param &Options [out] Global model options information
//
void CEnsemble::Initialize(const optStruct &Options) {
  //default does nothing - abstract base class
}

//////////////////////////////////////////////////////////////////
/// \brief updates model - called prior to each model ensemble run 
/// \param pModel [out] pointer to global model instance
/// \param &Options [out] Global model options information
//
void CEnsemble::UpdateModel(CModel *pModel,optStruct &Options, int e) {
  //default does nothing - abstract base class
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// MONTE CARLO CLASS

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////
/// \brief Monte Carlo Ensemble Construcutor 
/// \param num_members [in] number of MC ensemble realizations
/// \param &Options [out] Global model options information
//
CMonteCarloEnsemble::CMonteCarloEnsemble(const int num_members,const optStruct &Options) 
                    :CEnsemble(num_members,Options)
{
  _type=ENSEMBLE_MONTECARLO;
  _nParamDists=0;
  _pParamDists=NULL;
}
//////////////////////////////////////////////////////////////////
/// \brief Monte Carlo Ensemble Destrucutor 
//
CMonteCarloEnsemble::~CMonteCarloEnsemble() 
{
  for(int i=0;i<_nParamDists;i++) {
    delete _pParamDists[i];
  }
  delete [] _pParamDists; _nParamDists=0;
}
//////////////////////////////////////////////////////////////////
/// \brief Adds parameter distribution to MC setup
/// \param dist [in] pointer to parameter distribution structure for specified parameter
//
void CMonteCarloEnsemble::AddParamDist(const param_dist *dist) 
{
  if(!DynArrayAppend((void**&)(_pParamDists),(void*)(dist),_nParamDists)) {
    ExitGracefully("CMonteCarloEnsemble::AddParamDist: adding NULL distribution",BAD_DATA);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief initializes ensemble
/// \param &Options [out] Global model options information
//
void CMonteCarloEnsemble::Initialize(const optStruct &Options)
{
  CEnsemble::Initialize(Options);
  ofstream MCOUT;
  MCOUT.open(Options.main_output_dir+"MonteCarloOutput.csv");
  if(MCOUT.fail()) {
    string warn="CMonteCarloEnsemble::Initialize: unable to open file "+Options.main_output_dir+"MonteCarloOutput.csv for writing";
    ExitGracefully(warn.c_str(),BAD_DATA);
  }
  MCOUT<<"eID,";
  for (int i=0;i<_nParamDists;i++){MCOUT<<_pParamDists[i]->param_name+" ("+_pParamDists[i]->class_group+"),";}
  MCOUT<<endl;
  MCOUT.close();
}
//////////////////////////////////////////////////////////////////
/// \brief updates model - called prior to each model ensemble run 
/// \param pModel [out] pointer to global model instance
/// \param &Options [out] Global model options information
//
void CMonteCarloEnsemble::UpdateModel(CModel *pModel,optStruct &Options, const int e) 
{
  ExitGracefullyIf(e>=_nMembers,"CMonteCarloEnsemble::UpdateMode: invalid ensemble member index",RUNTIME_ERR);

  //- update output file/ run names ----------------------------
  Options.output_dir=_aOutputDirs[e];
  Options.run_name  =_aRunNames[e];

  //- Update parameter values ----------------------------------
  ofstream MCOUT;
  MCOUT.open(Options.main_output_dir+"MonteCarloOutput.csv",ios::app);
  MCOUT<<e+1<<", ";

  double val;
  for(int i=0;i<_nParamDists;i++) 
  {
    val=SampleFromDistribution(_pParamDists[i]);
    pModel->UpdateParameter(_pParamDists[i]->param_class,
                            _pParamDists[i]->param_name,
                            _pParamDists[i]->class_group,
                            val);
    MCOUT<<to_string(val)<<", ";
  //  cout<<"RAND PARAM: "<<val<<" between "<<_pParamDists[i]->distpar[0]<<" and "<< _pParamDists[i]->distpar[1]<<endl;
  }
  MCOUT<<endl;

  //- Re-read initial conditions to update state variables----
  if(!ParseInitialConditionsFile(pModel,Options)) {
    ExitGracefully("Cannot find or read .rvc file",BAD_DATA);}

  MCOUT.close();
}
//////////////////////////////////////////////////////////////////
/// \brief updates model - called prior to each model ensemble run 
/// \param pModel [out] pointer to global model instance
/// \param &Options [out] Global model options information
//
double SampleFromDistribution(param_dist *dist) 
{
  double value;

  //std::mt19937 generator;

  if(dist->distribution==DIST_UNIFORM) 
  {
    //std::uniform_real_distribution<double> distribution(dist->distpar[0],dist->distpar[1]);
    //value=distribution(generator);
    value=dist->distpar[0]+(dist->distpar[1]-dist->distpar[0])*UniformRandom();
  }
  else if(dist->distribution==DIST_NORMAL) 
  {  
    //std::normal_distribution<double> distribution(dist->distpar[0],dist->distpar[1]);
    value=0;
    //value = distribution(generator);
  }
  return value;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// DDS ENSEMBLE CLASS

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////
/// \brief Monte Carlo Ensemble Construcutor 
/// \param num_members [in] number of MC ensemble realizations
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
}
//////////////////////////////////////////////////////////////////
/// \brief Monte Carlo Ensemble Destrucutor 
//
CDDSEnsemble::~CDDSEnsemble()
{
  for(int i=0;i<_nParamDists;i++) {
    delete _pParamDists[i];
  }
  delete [] _pParamDists; _nParamDists=0;
  delete [] _BestParams;
  delete [] _TestParams;
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

void CDDSEnsemble::Initialize(const optStruct &Options) {

  _BestParams=new double [_nParamDists];
  _TestParams=new double [_nParamDists];
  for(int i=0; i<_nParamDists;i++)
  {
    _BestParams[i]=_TestParams[i]=_pParamDists[i]->default_val;
  }
}

void CDDSEnsemble::SetPerturbationValue(const double &perturb) {
  _r_val=perturb;
  if((_r_val<0) || (_r_val>1.0)) {
    ExitGracefully("CDDSEnsemble::SetPerturbationValue",BAD_DATA_WARN);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief updates model - called prior to each model ensemble run 
/// \param pModel [out] pointer to global model instance
/// \param &Options [out] Global model options information
//
void CDDSEnsemble::UpdateModel(CModel *pModel,optStruct &Options,const int e)
{
  ExitGracefullyIf(e>=_nMembers,"CDDSEnsemble::UpdateMode: invalid ensemble member index",RUNTIME_ERR);

  int iters_remaining=_nMembers-e;

  //- update output file/ run names ----------------------------
  Options.output_dir=_aOutputDirs[e];
  Options.run_name  =_aRunNames[e];

  //- Update parameter values ----------------------------------
  // Determine variable selected as neighbour 
  double Pn=1.0-log(double(e))/log(double(_nMembers));
  int dvn_count=0;
  
  // define TestParams initially as best current solution
  for(int k=0;k<_nParamDists;k++)
  {
    _TestParams[k]=_BestParams[k];
  }

  for(int k=0;k<_nParamDists;k++)
  {
    if(UniformRandom()<Pn) {
      dvn_count++;
      _TestParams[k]=PerturbParam(_BestParams[k],_pParamDists[k]->distpar[1],_pParamDists[k]->distpar[1]);
    }
  }
  if(dvn_count==0) {
    int dv=(int)(ceil((double)(_nParamDists)*UniformRandom()))-1; // index for one DV 
    _TestParams[dv]=PerturbParam(_BestParams[dv],_pParamDists[dv]->distpar[1],_pParamDists[dv]->distpar[1]);
  }

  //- Re-read initial conditions to update state variables----
  if(!ParseInitialConditionsFile(pModel,Options)) {
    ExitGracefully("Cannot find or read .rvc file",BAD_DATA);
  }


}
void CDDSEnsemble::FinishEnsembleRun(CModel *pModel,optStruct &Options,const int e)
{
  double Ftest=-1;//=pModel->GetObjFuncVal(); //TMP DEBUG - must extract objective function!
  
  if(Ftest<=_Fbest) // update current (best) solution
  {
    _Fbest = Ftest;
    for(int k=0;k<_nParamDists;k++) { _BestParams[k]=_TestParams[k]; }
    
    //write results
    ofstream DDSOUT;
    DDSOUT.open(Options.main_output_dir+"DDSOutput.csv",ios::app);
    DDSOUT<<e+1<<", "<<_Fbest<<",";
    DDSOUT<<endl;
    DDSOUT.close();
  }
  if(e==_nMembers-1) {
    //write best parameter vector?

  }

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
                                  const double &upperbound, 
                                  const double &lowerbound)
{

  double x_new;

  double x_max = upperbound;
  double x_min = lowerbound;

  x_new=x_best+GaussRandom()*_r_val*(x_max-x_min);

  // need if statements to check within DV bounds.  If not, bounds are reflecting.
  if(x_new<x_min)
  {
    x_new=x_min+(x_min-x_new); //reflect
                               // if reflection goes past x_max then value should be x_min since 
                               // without reflection the approach goes way past lower bound.  
                               // This keeps x close to lower bound when x_best is close to lower bound
    if(x_new>x_max) { x_new=x_min; }
  }
  else if(x_new>x_max)
  {
    x_new=x_max-(x_new-x_max);//reflect

    if(x_new<x_min) { x_new=x_max; }
  }

  return x_new;
}