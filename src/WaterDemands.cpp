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

  _loc_index=DOESNT_EXIST;
  _unrestricted=0;
  _cumDelivDate=0; //Jan-1
  
  _penalty=1.0;
  _multiplier=1.0;

  _pDemandTS=NULL;

  _currentDemand=0;
}
CDemand::CDemand(int ID, string name, long SBID, bool is_res, CTimeSeries* pTS):CDemand(ID,name,SBID,is_res) 
{
  _pDemandTS=pTS;
}
CDemand::~CDemand() 
{
  //delete _pDemandTS; //only copy of time series here - actual is tied to subbasin 
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
//////////////////////////////////////////////////////////////////
/// \brief Water demand manipulators
//
void    CDemand::SetLocalIndex(const int ii) 
{
  _loc_index=ii;
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
//////////////////////////////////////////////////////////////////
/// \brief re-calculates current demand magnitude (_currentDemand) - called at start of time step
//
void    CDemand::UpdateDemand(const optStruct &Options,const time_struct& tt) 
{
  //if (_demand_type==DEM_TIME_SERIES)
  if (_pDemandTS!=NULL)
  {
    //int nn=tt.nn;
    int nn=(int)((tt.model_time+TIME_CORRECTION)/Options.timestep);//current timestep index
    
    double Qirr=_pDemandTS->GetSampledValue(nn);
    if (Qirr==RAV_BLANK_DATA){Qirr=0.0;}

    _currentDemand=_multiplier*Qirr;
  }
  /*else if (_demand_type == DEM_AET_PERCENTAGE) {

  }
  else if (_demand_type == DEM_SOIL_MOISTURE) {

  }*/
}