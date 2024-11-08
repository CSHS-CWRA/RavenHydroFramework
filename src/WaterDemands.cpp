/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2024 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Demands.h"

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the water demand constructor
/// \details Creates empty instance of the water demand class
//
CDemand::CDemand(int ID, string name, long SBID, bool is_res)
{
  _ID=ID;
  _name=name;
  _SBID=SBID;
  _is_reservoir=is_res;

  _demType=DEMAND_TIME_SERIES;
  _retType=RETURN_MAX;

  _loc_index=DOESNT_EXIST;
  _global_index=DOESNT_EXIST;
  _unrestricted=0;
  _cumDelivDate=0; //Jan-1

  _penalty=1.0;
  _multiplier=1.0;
  _demandFract=0.0;//default

  _pDemandTS=NULL;
  _pReturnTS=NULL;
  _pDemandExp=NULL;
  _pReturnExp=NULL;
  _targetSBID=_SBID;           //default return is to same location
  _irrigHRUGroup=DOESNT_EXIST;
  _returnPct=0.0;              //default is no return

  _currentDemand=0;
  _currentRetTarget=0;
}
CDemand::CDemand(int ID, string name, long SBID, bool is_res, CTimeSeries* pTS):CDemand(ID,name,SBID,is_res)
{
  _pDemandTS=pTS;
  _demType=DEMAND_TIME_SERIES;
}
CDemand::~CDemand()
{
  delete _pDemandTS;
  delete _pReturnTS;
  delete _pDemandExp;
  delete _pReturnExp;
}
//////////////////////////////////////////////////////////////////
/// \brief Water demand accessors
//
int     CDemand::GetID() const
{
  return _ID;
}
string  CDemand::GetName() const
{
  return _name;
}
int     CDemand::GetGlobalIndex() const
{
  ExitGracefullyIf(_global_index==DOESNT_EXIST,"CDemand::GetLocalIndex(): didnt set local index",RUNTIME_ERR);
  return _global_index; //d
}
int     CDemand::GetLocalIndex() const
{
  ExitGracefullyIf(_loc_index==DOESNT_EXIST,"CDemand::GetLocalIndex(): didnt set local index",RUNTIME_ERR);
  return _loc_index;
}
long    CDemand::GetSubBasinID() const
{
  return _SBID;
}
double  CDemand::GetPenalty() const
{
  return _penalty;
}
bool    CDemand::IsUnrestricted() const
{
  return _unrestricted;
}
int     CDemand::GetCumulDeliveryDate() const
{
  return _cumDelivDate;
}
bool    CDemand::IsReservoirDemand() const
{
  return _is_reservoir;
}
double  CDemand::GetDemand() const
{
  return _currentDemand;
}
double  CDemand::GetReturnFlowTarget() const
{
  return _currentRetTarget;
}
double  CDemand::GetReturnFlowFraction() const
{
  return _returnPct;
}
// returns true if this demand has a return flow
bool    CDemand::HasReturnFlow() const
{
  return (_returnPct*_multiplier>0.0);
}
long    CDemand::GetTargetSBID() const {
  return _targetSBID;
}
//////////////////////////////////////////////////////////////////
/// \brief Water demand manipulators
//
void    CDemand::SetLocalIndex(const int ii)
{
  _loc_index=ii;
}
void    CDemand::SetGlobalIndex(const int d)
{
  _global_index=d;
}
void    CDemand::SetDemandPenalty(const double& P)
{
  _penalty=P;
}
void    CDemand::SetCumulDeliveryDate(const int date)
{
  _cumDelivDate=date;
}
void    CDemand::SetAsUnrestricted()
{
  _unrestricted=true;
}
void    CDemand::SetMultiplier(const double& M)
{
  _multiplier*=M; //allows for multiple calls to set multiplier
}
void    CDemand::SetDemandTimeSeries(CTimeSeries* pTS) {
  //todo: throw warning if _pDemandTS exists (AND DELETE PREVIOUS TS - OR COULD CAUSE PROBLEMS WITH ENSEMBLES)
  _demType=DEMAND_TIME_SERIES;
  _pDemandTS=pTS;
}
void    CDemand::SetReturnTimeSeries(CTimeSeries* pTS) {
  //todo: throw warning if _pReturnTS exists
  if (_returnPct==0.0){_returnPct=1.0;} // problem: if multiplier intentionally set to zero, this overrides it
  _retType=RETURN_TIME_SERIES;
  _pReturnTS=pTS;
}
void    CDemand::SetReturnFraction(const double &val)
{
  ExitGracefullyIf((val>1.0) || (val<0.0),"CDemand::SetReturnFraction: demand must be set between 0 and 1",BAD_DATA_WARN);
  _returnPct=val;
}
void    CDemand::SetDemandFraction(const double &val)
{
  ExitGracefullyIf((val>1.0) || (val<0.0),"CDemand::SetDemandFraction: demand must be set between 0 and 1",BAD_DATA_WARN);
  _demType=DEMAND_PCT;
  _demandFract=val;
}
void    CDemand::SetDemandExpression(expressionStruct* pExp)
{
  //TODO: check for valid LHS ?
  _demType=DEMAND_EXPRESSION;
  _pDemandExp=pExp;
}
void    CDemand::SetReturnExpression(expressionStruct* pExp)
{
  //TODO: check for valid LHS ?
  _retType=RETURN_EXPRESSION;
  _pReturnExp=pExp;
}
void    CDemand::SetTargetSBID(const long ID)
{
  _targetSBID=ID;
}
void    CDemand::SetIrrigationGroup(const int kk)
{
  _targetSBID=DOESNT_EXIST;
  _irrigHRUGroup=kk;
}
//////////////////////////////////////////////////////////////////
/// \brief initialized demand and return time series prior to simulation
/// \param &Options [in] Global model options information
//
void    CDemand::Initialize(const optStruct &Options)
{
  if (_pDemandTS != NULL) {
    _pDemandTS->Initialize(Options.julian_start_day,Options.julian_start_year,Options.duration,Options.timestep,false,Options.calendar);
  }
  if (_pReturnTS != NULL) {
    _pReturnTS->Initialize(Options.julian_start_day,Options.julian_start_year,Options.duration,Options.timestep,false,Options.calendar);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief re-calculates current demand magnitude (_currentDemand) - called at start of time step
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time structure
//
void    CDemand::UpdateDemand(const optStruct &Options,const time_struct& tt)
{
  // Calculate Delivery Targets
  //-------------------------------------------------------------------
  if      (_demType==DEMAND_TIME_SERIES)
  {
    //int nn=tt.nn;
    int nn=(int)((tt.model_time+TIME_CORRECTION)/Options.timestep);//current timestep index

    double Qirr=_pDemandTS->GetSampledValue(nn);
    if (Qirr==RAV_BLANK_DATA){Qirr=0.0;}

    _currentDemand=_multiplier*Qirr;
    //cout<<"UPDATE DEMAND : "<<_currentDemand<<" "<<_multiplier<<" "<<Qirr<< " " <<nn<<endl;
  }
  else if (_demType == DEMAND_PCT)
  {
    //pct of inflow from end of previous time step
    //double Q=pModel->GetSubBasin(p)->GetInflow(); // outflow from last time step
    //_currentDemand=demPct * Q;
    ExitGracefully("UpdateDEmand::DEMAND_PCT",STUB);
  }
  else if (_demType == DEMAND_EXPRESSION)
  {
    //double val=_pModel->GetDO->EvaluateExpression(pExp,tt.model_time);
    //if (fabs(val-RAV_BLANK_DATA)<REAL_SMALL){_currentDemand=0.0;}
    //else                                    {_currentDemand=val;}
    ExitGracefully("UpdateDEmand::DEMAND_EXPRESSION",STUB);
  }
  /*
  else if (_demType == DEM_SOIL_MOISTURE) {

  }*/

  // Calculate Return Flow Targets
  //-------------------------------------------------------------------
  _currentRetTarget=0;
  if (HasReturnFlow())
  {
    if      (_retType==RETURN_TIME_SERIES)
    {
      //int nn=tt.nn;
      int nn=(int)((tt.model_time+TIME_CORRECTION)/Options.timestep);//current timestep index

      double Qret=_pReturnTS->GetSampledValue(nn);
      if (Qret==RAV_BLANK_DATA){Qret=0.0;}

      _currentRetTarget=_multiplier*Qret;
    }
    else if (_retType == RETURN_MAX) {
      _currentRetTarget = 1e8;
    }
    else if (_retType == RETURN_EXPRESSION)
    {
      ExitGracefully("UpdateDEmand::RETURN_EXPRESSION",STUB);
    }
    _currentRetTarget=min(_currentRetTarget,_currentDemand); //this is likely unnecessary
  }
}
