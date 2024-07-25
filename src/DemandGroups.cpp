/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2024 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Demands.h"

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the demand group constructor
/// \details Creates empty demand group. Sets demand group name, initializes demandID array to NULL and _nDemands to 0
/// \param tag [in] Name of demand group
//
CDemandGroup::CDemandGroup(string tag, int global_ind)
{
  _name=tag;
  _nDemands=0;
  _pDemands=NULL;
  _global_ii=global_ind;
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the Group destructor
//
CDemandGroup::~CDemandGroup()
{
  delete [] _pDemands; _pDemands=NULL; //deletes pointers only
}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of demands in instantiated demand group
/// \return Number of demands in group
//
int CDemandGroup::GetNumDemands        () const{return _nDemands;}

//////////////////////////////////////////////////////////////////
/// \brief Returns name of demand group
/// \return Name of group
//
string CDemandGroup::GetName () const {return _name;}

//////////////////////////////////////////////////////////////////
/// \brief Returns unique demand group global model index (pp)
///
/// \return Integer index of demand in global model array
//
int    CDemandGroup::GetGlobalIndex       () const {return _global_ii;}

//////////////////////////////////////////////////////////////////
/// \brief Returns true if demands with ID demandID is in group
/// \param demandID [in] demand id
///
/// \return true if demand with  ID demandID is in group
//
bool  CDemandGroup::IsInGroup          (const int demandID) const
{
  for (int p=0;p<_nDemands; p++){
    if (_pDemands[p]->GetID() == demandID) { return true; }
  }
  return false;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns demand corresponding to index ii in group
/// \param ii [in] Index referring to pth element of the demand Group
/// \return demand corresponding to index ii
//
CDemand* CDemandGroup::GetDemand(const int ii) const
{
  return _pDemands[ii];
}
//////////////////////////////////////////////////////////////////
/// \brief Returns demand name corresponding to index ii in group
/// \param ii [in] Index referring to pth element of the demand Group
/// \return demand name corresponding to index ii
//
/*string CDemandGroup::GetDemandName(const int ii) const
{
  ExitGracefullyIf((ii<0) || (ii>=_nDemands),"CDemandGroup GetDemandName::improper index",BAD_DATA);
  return _aDemandNames[ii];
}*/

//////////////////////////////////////////////////////////////////
/// \brief Add an demand to group by dynamically appending to array
//
void CDemandGroup::AddDemand(CDemand *pDem)
{
  if (!DynArrayAppend((void**&)(_pDemands),(void*)(pDem),_nDemands)){
     ExitGracefully("CDemandGroup::AddDemand: adding NULL demand",BAD_DATA);}
}
//////////////////////////////////////////////////////////////////
/// \brief initializes demand Groups
//
void CDemandGroup::Initialize()
{
}
