/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ----------------------------------------------------------------*/
#include "HydroUnits.h"

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the HRU group constructor
/// \details Creates empty HRU group. Sets HRU group name, initializes HRU array to NULL and nHRUs to 0
/// \param tag [in] Name of HRU group
//
CHRUGroup::CHRUGroup(string tag, int global_ind)
{
  _name=tag;
  _nHRUs=0;  _pHRUs=NULL;
  _aAggregateSV=NULL;
  _aAggregateSV=new bool [MAX_STATE_VARS];
  for (int i=0;i<MAX_STATE_VARS;i++){
    _aAggregateSV[i]=false;
  }
  _global_kk=global_ind;
  _disabled=false;
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the HRU Group destructor
/// \details Deletes pointers to HRUs only
//
CHRUGroup::~CHRUGroup()
{
  delete [] _pHRUs; _pHRUs=NULL; //deletes pointers only
  delete [] _aAggregateSV;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of HRUs in instantiated HRU group
/// \return Number of HRUs in group
//
int CHRUGroup::GetNumHRUs        () const{return _nHRUs;}

//////////////////////////////////////////////////////////////////
/// \brief Returns name of HRU group
/// \return Name of HRU group
//
string CHRUGroup::GetName () const {return _name;}

//////////////////////////////////////////////////////////////////
/// \brief Returns unique HRU group global model index (kk)
///
/// \return Integer index of HRU in global model array
//
int    CHRUGroup::GetGlobalIndex       () const {return _global_kk;}

//////////////////////////////////////////////////////////////////
/// \brief Returns true if HRU with global index k_global is in group
/// \param k_global [in] Global model index of an HRU
///
/// \return true if HRU with global index k_global is in group
//
bool  CHRUGroup::IsInGroup          (const int k_global) const
{
  for (int k=0;k<_nHRUs; k++){
    if (_pHRUs[k]->GetGlobalIndex()==k_global){return true;}
  }
  return false;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns HRU corresponding to index k in group
/// \param k [in] Index referring to kth element of the HRU Group
/// \return HRU corresponding to index k
//
CHydroUnit *CHRUGroup::GetHRU(const int k) const
{
  ExitGracefullyIf((k<0) || (k>=_nHRUs),"CHRUGroup GetHydroUnit::improper index",BAD_DATA);
  return _pHRUs[k];
}

//////////////////////////////////////////////////////////////////
/// \brief Add an HRU to HRU group by dynamically appending to array pHRUs
//
void CHRUGroup::AddHRU(CHydroUnit *pHRU)
{
  if (!DynArrayAppend((void**&)(_pHRUs),(void*)(pHRU),_nHRUs)){
   ExitGracefully("CHRUGroup::AddHRU: adding NULL HRU",BAD_DATA);} 
} 
//////////////////////////////////////////////////////////////////
/// \brief initializes HRU Groups
//
void CHRUGroup::Initialize()
{
  if(_disabled)
  {
    for(int k=0;k<_nHRUs;k++){
      _pHRUs[k]->Disable();
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief disables HRU Group
//
void CHRUGroup::DisableGroup()
{
  _disabled=true;
}
//////////////////////////////////////////////////////////////////
/// \brief returns true if state variable i is aggregated across this HRU group
/// \returns true if state variable i is aggregated across this HRU group
//
bool  CHRUGroup::IsAggregatorGroup   (const int i) const
{
  return _aAggregateSV[i];
}
//////////////////////////////////////////////////////////////////
/// \brief ensures that state variable i is aggregated across this HRU group
//
void  CHRUGroup::SetAsAggregator    (const int i){
  _aAggregateSV[i]=true;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns average value of a state variable specified by index i over the total area covered by the HRU group
/// \param i [in] Index corresponding to the state variable whose average will be calculated
/// \return Average of state variable with index i across all HRUs in group, per unit area coverage of group
//
double CHRUGroup::GetAvgStateVar (const int i) const
{
  double sum=0.0;
  double areasum=0.0;
  double area;
  for (int k=0;k<_nHRUs;k++)
  {
    area    =_pHRUs[k]->GetArea();
    sum    +=_pHRUs[k]->GetStateVarValue(i)*area;
    areasum+=area;
  }
  return sum/areasum;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns average value of a state variable specified by index i over the total area covered by the HRU group
/// \param i [in] Index corresponding to the state variable whose average will be calculated
/// \return Average of state variable with index i across all HRUs in group, per unit area coverage of group
//
double CHRUGroup::GetAvgConcentration(const int i) const
{
  double sum=0.0;
  double areasum=0.0;
  double area;
  for(int k=0;k<_nHRUs;k++)
  {
    area    =_pHRUs[k]->GetArea();
    sum    +=_pHRUs[k]->GetConcentration(i)*area;
    areasum+=area;
  }
  return sum/areasum;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns average value of a forcing function specified by forcing_string over the total area covered by the HRU group
/// \param &forcing_string [in] Index corresponding to the state variable whose average will be calculated
/// \return Average of forcing function with identifier forcing_string, across all HRUs in group, per unit area coverage of group
//
double CHRUGroup::GetAvgForcing (const string &forcing_string) const
{
  double sum=0.0;
  double areasum=0.0;
  double area;
  for (int k=0;k<_nHRUs;k++)
  {
    if(_pHRUs[k]->IsEnabled()){
      area    =_pHRUs[k]->GetArea();
      sum    +=_pHRUs[k]->GetForcing(forcing_string)*area;
      areasum+=area;
    }
  }
  return sum/areasum;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns area-weighted average of specified cumulative flux over HRU group
///
/// \param i [in] index of storage compartment
/// \param to [in] true if evaluating cumulative flux to storage compartment, false for 'from'
/// \return Area-weighted average of cumulative flux to storage compartment i
//
double CHRUGroup::GetAvgCumulFlux (const int i, const bool to) const
{
  //Area-weighted average
  double sum=0.0;
  double areasum=0.0;
  double area;
  for (int k=0;k<_nHRUs;k++)
  {
    if(_pHRUs[k]->IsEnabled()){
      area    =_pHRUs[k]->GetArea();
      sum    +=_pHRUs[k]->GetCumulFlux(i,to)*area;
      areasum+=area;
    }
  }
  return sum/areasum;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns area-weighted average of  cumulative flux between two compartments over HRU Group
///
/// \param iFrom [in] index of 'from' storage compartment
/// \param iTo [in] index of 'to' storage compartment
/// \return Area-weighted average of cumulative flux between two compartments over HRU Group
//
double CHRUGroup::GetAvgCumulFluxBet (const int iFrom, const int iTo) const
{
  //Area-weighted average
  double sum=0.0;
  double areasum=0.0;
  double area;
  for (int k=0;k<_nHRUs;k++)
  {
    if(_pHRUs[k]->IsEnabled()){
      area    =_pHRUs[k]->GetArea();
      sum    +=_pHRUs[k]->GetCumulFluxBet(iFrom,iTo)*area;
      areasum+=area;
    }
  }
  return sum/areasum;
}
