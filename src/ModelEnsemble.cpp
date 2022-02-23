/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2021 the Raven Development Team
----------------------------------------------------------------*/
#include "ModelEnsemble.h"
//#include <random>
//see http://anadoxin.org/blog/c-shooting-yourself-in-the-foot-4.html

bool ParseInitialConditionsFile(CModel *&pModel,const optStruct &Options);

//////////////////////////////////////////////////////////////////
/// \brief returns uniformly distributed random variable between 0 and 1
/// \return uniformly distributed random variable between 0 and 1
//
double UniformRandom()
{
  return (double)(rand())/RAND_MAX; 
}
//////////////////////////////////////////////////////////////////
/// \brief returns normally distributed random variable with mean of 0, variance=1
/// \notes uses Box-Muller transform
/// \return normally distributed random variable with mean of 0, variance=1
//
double GaussRandom()
{
  double u1=(double)(rand())/RAND_MAX;
  double u2=(double)(rand())/RAND_MAX;
  return sqrt(-2.0*log(u1))*cos(2.0*PI*u2);
}

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

  for(int e=0;e<_nMembers;e++) //default - everything runs in same directory
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
int CEnsemble::GetNumMembers() const {
  return _nMembers;
}
//////////////////////////////////////////////////////////////////
/// \brief Accessor - gets ensemble type
/// \return ensemble type
//
ensemble_type CEnsemble::GetType() const {
  return _type;
}
//////////////////////////////////////////////////////////////////
/// \brief Accessor - gets ensemble type
/// \return ensemble type
//
double CEnsemble::GetStartTime(const int e) const {
  return 0.0;
}

//Manipulator Functions
//////////////////////////////////////////////////////////////////
/// \brief sets random seed
/// \param [in] random seed
//
void CEnsemble::SetRandomSeed(const unsigned int seed) 
{
  srand(seed);
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
void CEnsemble::Initialize(const CModel* pModel,const optStruct &Options) {
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
void CMonteCarloEnsemble::Initialize(const CModel* pModel,const optStruct &Options)
{
  CEnsemble::Initialize(pModel,Options);
  ofstream MCOUT;
  string filename=Options.main_output_dir+"MonteCarloOutput.csv";
  MCOUT.open(filename.c_str());
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
  string filename=Options.main_output_dir+"MonteCarloOutput.csv";
  MCOUT.open(filename.c_str(),ios::app);
  MCOUT<<e+1<<", ";

  double val;
  for(int i=0;i<_nParamDists;i++) 
  {
    val=SampleFromDistribution(_pParamDists[i]->distribution,_pParamDists[i]->distpar);
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
double SampleFromDistribution(disttype distribution, double distpar[3]) 
{
  double value=0;

  //std::mt19937 generator;

  if(distribution==DIST_UNIFORM) 
  {
    //std::uniform_real_distribution<double> distribution(dist->distpar[0],dist->distpar[1]);
    //value=distribution(generator);
    value=distpar[0]+(distpar[1]-distpar[0])*UniformRandom();
  }
  else if(distribution==DIST_NORMAL) 
  {  
    //std::normal_distribution<double> distribution(dist->distpar[0],dist->distpar[1]);
    //value = distribution(generator);
    value=distpar[0]+(distpar[1]*GaussRandom());
  }
  return value;
}

