/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2024 the Raven Development Team
  ----------------------------------------------------------------*/
#ifndef DEMANDS_H
#define DEMANDS_H

#include "RavenInclude.h"
#include "TimeSeries.h"

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for irrigation or other water demand
//
class CDemand
{
private:/*------------------------------------------------------*/
  int       _ID;           ///< unique integer identifier
  string    _name;         ///< unique name/alias identifier 

  long      _SBID;         ///< subbasin ID 
  int       _loc_index;    ///< local demand index (counter in each subbasin) 
  bool      _is_reservoir; ///< true if withdrawal is from reservoir 

  double    _penalty;      ///< penalty for not satisfying demand [s/m3]

  bool      _unrestricted; ///< true if demand is unrestricted by downstream flow regulations [default: false]

  int       _cumDelivDate; ///< julian date to calculate cumulative deliveries from {default: Jan 1}
  //double _annualLicense; ///< annual maximum allocation

  //int    _irrigHRUGroup; ///< index kk of HRU group on which withdrawn water is applied
  //double _returnPct;     ///< percentage of flow which returns to stream

  CTimeSeries *_pDemandTS; ///< pointer to time series of demands (or NULL, if calculated elsewise)

  double   _currentDemand; ///< (time-variable) demand in current time step [m3/s]

public:/*-------------------------------------------------------*/
  CDemand(int ID, string name, long SBID, bool is_res);
  CDemand(int ID, string name, long SBID, bool is_res, CTimeSeries *pTS);
  ~CDemand();

  int     GetID() const;
  string  GetName() const;
  long    GetSubBasinID() const;
  int     GetLocalIndex() const;
  double  GetPenalty() const;
  bool    IsUnrestricted() const;
  int     GetCumulDeliveryDate() const;
  bool    IsReservoirDemand() const;

  double  GetDemand() const;

  //Manipulators
  void    SetLocalIndex(const int ii);
  void    SetDemandPenalty(const double &P);
  void    SetCumulDeliveryDate(const int date);
  void    SetAsUnrestricted();

  //Called during simulation
  void    UpdateDemand(const optStruct &Options,const time_struct &tt);
};
///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for irrigation or other water demand
//
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
  //string          GetDemandName      (const int ii_local) const;
  //CDemand        *GetDemand          (const int ii_local) const;
  bool              IsInGroup          (const int demandID) const;
  bool              IsDisabled         () const;
};

#endif