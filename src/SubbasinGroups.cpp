/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2020 the Raven Development Team
  ----------------------------------------------------------------*/
#include "SubBasin.h"

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the Subbasin group constructor
/// \details Creates empty Subbasin group. Sets subbasin group name, initializes subbasin array to NULL and _nSubBasins to 0
/// \param tag [in] Name of HRU group
//
CSubbasinGroup::CSubbasinGroup(string tag, int global_ind)
{
  _name=tag;
  _nSubbasins=0;  _pSubbasins=NULL;
  _global_pp=global_ind;
  _disabled=false;
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the Group destructor
/// \details Deletes pointers to basins only
//
CSubbasinGroup::~CSubbasinGroup()
{
  delete [] _pSubbasins; _pSubbasins=NULL; //deletes pointers only
}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of Subbasins in instantiated subbasin group
/// \return Number of Subbasins in group
//
int CSubbasinGroup::GetNumSubbasins        () const{return _nSubbasins;}

//////////////////////////////////////////////////////////////////
/// \brief Returns name of Subbasin group
/// \return Name of group
//
string CSubbasinGroup::GetName () const {return _name;}

//////////////////////////////////////////////////////////////////
/// \brief Returns unique subbasin group global model index (pp)
///
/// \return Integer index of Subbasin in global model array
//
int    CSubbasinGroup::GetGlobalIndex       () const {return _global_pp;}

//////////////////////////////////////////////////////////////////
/// \brief Returns true if subbasins with ID SBID is in group
/// \param p_global [in] Global model index of a Subbasins
///
/// \return true if Subbasin with subbasin ID SBID is in group
//
bool  CSubbasinGroup::IsInGroup          (const long SBID) const
{
  for (int p=0;p<_nSubbasins; p++){
    if (_pSubbasins[p]->GetID()==SBID){return true;}
  }
  return false;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns Subbasin corresponding to index p in group
/// \param p [in] Index referring to pth element of the Subbasin Group
/// \return Subbasin corresponding to index p
//
CSubBasin *CSubbasinGroup::GetSubBasin(const int p) const
{
  ExitGracefullyIf((p<0) || (p>=_nSubbasins),"CSubbasinGroup GetSubbasin::improper index",BAD_DATA);
  return _pSubbasins[p];
}
//////////////////////////////////////////////////////////////////
/// \return true if subbasin group is disabled
//
bool CSubbasinGroup::IsDisabled() const
{
  return _disabled;
}
//////////////////////////////////////////////////////////////////
/// \brief Add an Subbasin to group by dynamically appending to array
//
void CSubbasinGroup::AddSubbasin(CSubBasin *pSB)
{
  if (!DynArrayAppend((void**&)(_pSubbasins),(void*)(pSB),_nSubbasins)){
   ExitGracefully("CSubbasinGroup::AddSubbasin: adding NULL subbasin",BAD_DATA);} 
} 
//////////////////////////////////////////////////////////////////
/// \brief initializes Subbasin Groups
//
void CSubbasinGroup::Initialize()
{
  if(_disabled)
  {
    for(int p=0;p<_nSubbasins;p++){
      _pSubbasins[p]->Disable(); //this disables constituent HRUs
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief disables SubBasin Group
//
void CSubbasinGroup::DisableGroup()
{
  _disabled=true; //propagates once initialized
  for(int p=0;p<_nSubbasins;p++) {
    _pSubbasins[p]->Disable(); //this disables constituent HRUs
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Returns average value of a state variable specified by index i over the total area covered by the subbasin group
/// \param i [in] Index corresponding to the state variable whose average will be calculated
/// \return Average of state variable with index i across all subbasins in group, per unit area coverage of group
//
double CSubbasinGroup::GetAvgStateVar (const int i) const
{
  double sum=0.0;
  double areasum=0.0;
  double area;
  for (int p=0;p<_nSubbasins;p++)
  {
    area    =_pSubbasins[p]->GetBasinArea();
    sum    +=_pSubbasins[p]->GetAvgStateVar(i)*area;
    areasum+=area;
  }
  return sum/areasum;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns average value of a concentration/temperature specified by index i over the total area covered by the subbasin group
/// \param i [in] Index corresponding to the constituent state variable whose average will be calculated
/// \return Average of concentration/temperature with index i across all subbasins in group, per unit area coverage of group
//
double CSubbasinGroup::GetAvgConcentration(const int i) const
{
  double sum=0.0;
  double areasum=0.0;
  double area;
  for(int p=0;p<_nSubbasins;p++)
  {
    area    =_pSubbasins[p]->GetBasinArea();
    sum    +=_pSubbasins[p]->GetAvgConcentration(i)*area;
    areasum+=area;
  }
  return sum/areasum;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns average value of a forcing function specified by forcing_string over the total area covered by the HRU group
/// \param &forcing_string [in] Index corresponding to the state variable whose average will be calculated
/// \return Average of forcing function with identifier forcing_string, across all HRUs in group, per unit area coverage of group
//
double CSubbasinGroup::GetAvgForcing (const forcing_type &ftype) const
{
  double sum=0.0;
  double areasum=0.0;
  double area;
  for(int p=0;p<_nSubbasins;p++)
  {
    area    =_pSubbasins[p]->GetBasinArea();
    sum    +=_pSubbasins[p]->GetAvgForcing(ftype)*area;
    areasum+=area;
  }
  return sum/areasum;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns area-weighted average of specified cumulative flux over subbasin group
///
/// \param i [in] index of storage compartment
/// \param to [in] true if evaluating cumulative flux to storage compartment, false for 'from'
/// \return Area-weighted average of cumulative flux to storage compartment i
//
double CSubbasinGroup::GetAvgCumulFlux (const int i, const bool to) const
{
  //Area-weighted average
  double sum=0.0;
  double areasum=0.0;
  double area;
  for(int p=0;p<_nSubbasins;p++)
  {
    area    =_pSubbasins[p]->GetBasinArea();
    sum    +=_pSubbasins[p]->GetAvgCumulFlux(i,to)*area;
    areasum+=area;
  }
  return sum/areasum;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns area-weighted average of  cumulative flux between two compartments over subbasin Group
///
/// \param iFrom [in] index of 'from' storage compartment
/// \param iTo [in] index of 'to' storage compartment
/// \return Area-weighted average of cumulative flux between two compartments over subbasin Group
//
double CSubbasinGroup::GetAvgCumulFluxBet (const int iFrom, const int iTo) const
{
  //Area-weighted average
  double sum=0.0;
  double areasum=0.0;
  double area;
  for(int p=0;p<_nSubbasins;p++)
  {
    area    =_pSubbasins[p]->GetBasinArea();
    sum    +=_pSubbasins[p]->GetAvgCumulFluxBet(iFrom,iTo)*area;
    areasum+=area;
  }
  return sum/areasum;
}
