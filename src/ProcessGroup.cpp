/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2023 the Raven Development Team
  ----------------------------------------------------------------*/
#include "ProcessGroup.h"
#include "Model.h"

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the constructor for an abstract group of hydrological processes
//
CProcessGroup::CProcessGroup(string name, CModelABC *pModel)
  :CHydroProcessABC(PROCESS_GROUP, pModel)
{
  _pSubProcesses=NULL;
  _nSubProcesses=0;
  _aWeights=NULL;
  _name=name;
}
//////////////////////////////////////////////////////////////////
/// \brief Implementation of the CProcessGroup destructor
//
CProcessGroup::~CProcessGroup()
{
  delete[] _pSubProcesses; _pSubProcesses=NULL;
  delete[]_aWeights; _aWeights=NULL;
  _nSubProcesses=0;
}
//////////////////////////////////////////////////////////////////
/// \brief initialization: determines connection characteristics from subprocesses
/// \remark Called before solution
//
void CProcessGroup::Initialize()
{
  //populates iFrom, iTo, once group is populated
  int nConn;
  int N=0;
  for(int j=0;j<_nSubProcesses;j++)
  {
    nConn=_pSubProcesses[j]->GetNumConnections();
    N+=nConn;
  }
  DynamicSpecifyConnections(N);

  int q1=0;
  for(int j=0;j<_nSubProcesses;j++)
  {
    nConn=_pSubProcesses[j]->GetNumConnections();
    for(int qq=0;qq<nConn;qq++){
      iFrom[q1+qq]=_pSubProcesses[j]->GetFromIndices()[qq];
      iTo  [q1+qq]=_pSubProcesses[j]->GetToIndices()[qq];
    }
    q1+=nConn;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Returns rates of change in all state variables modeled over time step
/// \param *state_var [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Rate of loss of water from snow melt/refreeze combined [mm/d]
//
void CProcessGroup::GetRatesOfChange(const double      *state_vars,
                                      const CHydroUnit  *pHRU,
                                      const optStruct   &Options,
                                      const time_struct &tt,
                                      double      *rates) const
{
  double *loc_rates;
  loc_rates=new double[_nConnections];//has to be smaller than this
  for(int q=0;q<_nConnections;q++){ loc_rates[q]=0.0; }

  int q1=0;
  int nConn;
  for(int j=0;j<_nSubProcesses;j++)
  {
    nConn=_pSubProcesses[j]->GetNumConnections();
    _pSubProcesses[j]->GetRatesOfChange(state_vars,pHRU,Options,tt,loc_rates);
    _pSubProcesses[j]->ApplyConstraints(state_vars,pHRU,Options,tt,loc_rates);
    for(int qq=0;qq<nConn;qq++){
      rates[q1+qq]=_aWeights[j]*loc_rates[qq];
    }

    q1+=nConn;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Rate of change in state variables due to snow balance calculations
//
void CProcessGroup::ApplyConstraints( const double      *state_vars,
                                      const CHydroUnit  *pHRU,
                                      const optStruct   &Options,
                                      const time_struct &tt,
                                      double      *rates) const
{
  return;
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
///
/// \param bal_type [in] Snow energy/mass balance model selected
/// \param *aSV [out] Array of state variable types needed by energy balance algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by energy balance algorithm (size of aSV[] and aLev[] arrays)
//
void CProcessGroup::GetParticipatingStateVarList(sv_type *aSV,int *aLev,int &nSV)
{
  //this is not necessary; groups don't add new state vars - this has to be called for each process anyhow
}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for snow balance algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by snow balance algorithm (size of aP[] and aPC[])
//
void CProcessGroup::GetParticipatingParamList(string *aP,class_type *aPC,int &nP) const
{
  //generated by collating all participating parameters of subprocesses
  int nPP;
  const int MAX_PARAMS=100;
  string     *aPP =new string [MAX_PARAMS];
  class_type *aPCP=new class_type [MAX_PARAMS];
  nP=0;
  int pp=0;
  for(int i=0;i<_nSubProcesses;i++)
  {
    _pSubProcesses[i]->GetParticipatingParamList(aPP,aPCP,nPP);
    nP+=nPP;
    for(int p=0;p<nPP;p++) {
      aP [pp]=aPP[p];
      aPC[pp]=aPCP[p];
      pp++;
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief returns number of subprocesses in group
//
int CProcessGroup::GetGroupSize() const
{
  return _nSubProcesses;
}
//////////////////////////////////////////////////////////////////
/// \brief adds hydrologic process to list of subprocesses
///
/// \param *pProc [in] pointer to hydrologic process to be added
//
void CProcessGroup::AddProcess(CHydroProcessABC *pProc)
{
  if (!DynArrayAppend((void**&)(_pSubProcesses),(void*)(pProc),_nSubProcesses)){
    ExitGracefully("CModel::AddHRU: adding NULL HRU",BAD_DATA);}
  delete[] _aWeights;
  _aWeights=new double[_nSubProcesses];
  for(int i=0;i<_nSubProcesses;i++){ _aWeights[i]=1.0; } //defaults to 1.0
}
//////////////////////////////////////////////////////////////////
/// \brief sets process weights for averaging of processes
///
/// \param *aWts [in] array of process group weights
/// \param nVal [in] number of weights
//
void CProcessGroup::SetWeights(const double *aWts, const int nVal)
{
  if(nVal!=_nSubProcesses){
    WriteWarning("CProcessGroup::SetWeights: Incorrect number of weights for process group. Weights will be ignored",false);
    return;
  }
  double sum=0.0;
  for(int q=0; q<_nSubProcesses;q++){
    _aWeights[q]=aWts[q];
    sum+=_aWeights[q];
  }
  if(fabs(sum-1.0)>0.01) {
    for(int q=0; q<_nSubProcesses;q++) {_aWeights[q]/=sum; } //clean up weights that almost add to 1 (roundoff error)
  }
  if(fabs(sum-1.0)>REAL_SMALL){ WriteWarning("CProcessGroup::SetWeights: Process group weights do not sum to 1.0",false); }
}
//////////////////////////////////////////////////////////////////
/// \brief calculates N process weights from N-1 numbers ranging from 0 to 1
///
/// \param *aVals [in] array of weight seeds between 0 and 1
/// \param nVal [in] number of weight seeds
//
void CProcessGroup::CalcWeights(const double *aVals,const int nVal)
{
  if(nVal!=_nSubProcesses-1) {
    WriteWarning("CProcessGroup::SetWeights: Incorrect number of numbers for process group. Weights will be ignored",false);
    return;
  }
  CalcWeightsFromUniformNums(aVals, _aWeights, _nSubProcesses);
}
