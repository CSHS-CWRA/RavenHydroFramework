/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2024 the Raven Development Team
  ----------------------------------------------------------------*/
#ifndef DEMANDS_H
#define DEMANDS_H

#include "RavenInclude.h"
#include "TimeSeries.h"
#include "Expression.h"
#include "Model.h"

enum demand_type{
  DEMAND_TIME_SERIES,
  DEMAND_PCT,
  DEMAND_LOOKUP,
  DEMAND_EXPRESSION //could even link to soil moisture or AET
};
enum return_type{
  RETURN_TIME_SERIES,
  RETURN_MAX,
  RETURN_EXPRESSION
};

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for irrigation or other water demand
//
class CDemand
{
private:/*------------------------------------------------------*/
  long long            _ID;           ///< unique long long integer identifier
  string               _name;         ///< unique name/alias identifier

  const CModel        *_pModel;       ///< pointer to model object

  long long            _SBID;         ///< subbasin ID
  int                  _loc_index;    ///< local demand index ii (counter in each subbasin or reservoir)
  int                  _global_index; ///< global demand index d (from list of all demands in water management model)
  bool                 _is_reservoir; ///< true if withdrawal is from reservoir

  double               _penalty;      ///< penalty for not satisfying demand [s/m3]

  bool                _unrestricted; ///< true if demand is unrestricted by downstream flow regulations [default: false]

  int                 _cumDelivDate;  ///< julian date to calculate cumulative deliveries from {default: Jan 1}

  //Demand characterization
  demand_type         _demType;       ///< demand type
  double              _multiplier;    ///< multiplies time series or any other means of calculating demand
  CTimeSeries        *_pDemandTS;     ///< pointer to time series of demands (or NULL, if calculated elsewise)
  //CLookupTable     *_pDemandLT;     ///< pointer to demand lookup table D(Qin+Runoff) (or NULL)
  expressionStruct   *_pDemandExp;    ///< return expression (or NULL, if calculated elsewise)
  double              _demandFract;   ///< [0..1] percentage of flow demanded
  //double            _annualLicense; ///< annual maximum allocation
  //
  //Return flow variables
  return_type        _retType;       ///< return type
  long long          _targetSBID;    ///< subbasin id of return destination (defaults to _SBID, is -1 for irrigation to HRU group)
  int                _irrigHRUGroup; ///< index kk of HRU group on which withdrawn water is applied
  double             _returnPct;     ///< percentage of delivered demand which returns to stream or land
  CTimeSeries       *_pReturnTS;     ///< pointer to time series of target return flows (or NULL, if calculated elsewise)
  expressionStruct  *_pReturnExp;    ///< return expression (or NULL, if calculated elsewise)

  double             _currentDemand;   ///< (time-variable) demand in current time step [m3/s]
  double             _currentRetTarget;///< (time-variable) return flow target in current time step [m3/s]

public:/*-------------------------------------------------------*/
  CDemand(long long ID, string name, long long SBID, bool is_res, const CModel *pMod);
  CDemand(long long ID, string name, long long SBID, bool is_res, CTimeSeries *pTS, const CModel *pMod);
  ~CDemand();

  long long GetDemandID() const;
  string    GetName() const;
  long long GetSubBasinID() const;
  int       GetGlobalIndex() const;
  int       GetLocalIndex() const;
  double    GetPenalty() const;
  bool      IsUnrestricted() const;
  int       GetCumulDeliveryDate() const;
  bool      IsReservoirDemand() const;
  bool      HasReturnFlow() const;
  long long GetTargetSBID() const;
  double    GetReturnFlowFraction() const;

  double    GetDemand() const;
  double    GetReturnFlowTarget() const;

  //Manipulators
  void      SetLocalIndex(const int ii);
  void      SetGlobalIndex(const int d);
  void      SetDemandPenalty(const double &P);
  void      SetCumulDeliveryDate(const int date);
  void      SetAsUnrestricted();
  void      SetMultiplier(const double &M);
  void      SetTargetSBID(const long long ID);
  void      SetDemandFraction(const double &val);
  void      SetReturnFraction(const double &val);
  void      SetDemandTimeSeries(CTimeSeries *pTS);
  void      SetReturnTimeSeries(CTimeSeries *pTS);
  void      SetDemandExpression(expressionStruct *pExp);
  void      SetReturnExpression(expressionStruct *pExp);
  void      SetIrrigationGroup(const int kk);

  void      Initialize(const optStruct &Options);

  //Called during simulation
  void      UpdateDemand(const optStruct &Options,const time_struct &tt);
};
///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for irrigation or other water demand
//
class CDemandGroup
{
private:/*------------------------------------------------------*/

  string          _name;        ///< Demand group name
  int             _nDemands ;   ///< Number of Demands present in group
  CDemand       **_pDemands;    ///< Array of pointers to demand objects [size: _nDemands]
  int             _global_ii;   ///< index of group in master Demand Group array (in CModel)

public:/*-------------------------------------------------------*/
  //Constructors:
  CDemandGroup(string tag,int global_ind);
  ~CDemandGroup();

  void AddDemand(CDemand *pDem);
  void Initialize();

  string            GetName            () const;
  int               GetNumDemands      () const;
  int               GetGlobalIndex     () const;
  CDemand          *GetDemand          (const int ii_local) const;
  bool              IsInGroup          (const long long demandID) const;
};

#endif
