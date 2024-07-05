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
  _nDemands=0;  _aDemandIDs=NULL;
  _global_ii=global_ind;
  _disabled=false;
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the Group destructor
//
CDemandGroup::~CDemandGroup()
{
  delete [] _aDemandIDs; _aDemandIDs=NULL; //deletes pointers only
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
    if (_aDemandIDs[p]==demandID){return true;}
  }
  return false;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns demand corresponding to index ii in group
/// \param ii [in] Index referring to pth element of the demand Group
/// \return demand corresponding to index ii
//
int CDemandGroup::GetDemandID(const int ii) const
{
  ExitGracefullyIf((ii<0) || (ii>=_nDemands),"CDemandGroup GetDemandID::improper index",BAD_DATA);
  return _aDemandIDs[ii];
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
/// \return true if subbasin group is disabled
//
bool CDemandGroup::IsDisabled() const
{
  return _disabled;
}
//////////////////////////////////////////////////////////////////
/// \brief Add an Subbasin to group by dynamically appending to array
//
void CDemandGroup::AddDemand(int demandID)
{
  ExitGracefully("STUB: AddDemand",STUB);
}
//////////////////////////////////////////////////////////////////
/// \brief initializes demand Groups
//
void CDemandGroup::Initialize()
{
  if(_disabled)
  {
    for(int p=0;p<_nDemands;p++){
     // _aDemands[p]->Disable(); //this disables constituent demands
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief disables demand Group
//
void CDemandGroup::DisableGroup()
{
  _disabled=true; //propagates once initialized
  for(int p=0;p<_nDemands;p++) {
    // _aDemands[p]->Disable(); //this disables constituent demands
  }
}
