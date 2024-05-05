/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2024 the Raven Development Team
  ----------------------------------------------------------------*/
#ifndef DEMANDS_H
#define DEMANDS_H

#include "RavenInclude.h"

class CDemandGroup
{
private:/*------------------------------------------------------*/

  string          _name;        ///< Demand group name
  int             _nDemands ;   ///< Number of Demands present in group
  int             *_aDemandIDs; ///< Array of demand indices [size:_nDemands]
  //CDemand        **_pDemands;   ///< Array of pointers to demand objects [size: _nDemands]
  int             _global_ii;   ///< index of group in master Demand Group array (in CModel)
  bool            _disabled;    ///< true if all Subbasins in group are disabled

public:/*-------------------------------------------------------*/
  //Constructors:
  CDemandGroup(string tag,int global_ind);
  ~CDemandGroup();

  void AddDemand(int demandID);
  void DisableGroup();
  void Initialize();

  string            GetName            () const;
  int               GetNumDemands      () const;
  int               GetGlobalIndex     () const;
  int               GetDemandID        (const int ii_local) const;
  //CDemand        *GetDemand          (const int ii_local) cosnst;
  bool              IsInGroup          (const int demandID) const;
  bool              IsDisabled         () const;
};

#endif