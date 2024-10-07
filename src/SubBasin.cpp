/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2024 the Raven Development Team
  ----------------------------------------------------------------*/
#include "SubBasin.h"

/*****************************************************************
   Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the subbasin constructor
///
/// \param Identifier [in] Subbasin ID
/// \param Name       [in] Subbasin name
/// \param *pMod      [in] Pointer to surface water model
/// \param down_ID    [in] Index of downstream SB, if < 0, downstream outlet
/// \param *pChan     [in] Channel
/// \param reach_len  [in] Reach length [m]
/// \param Qreference [in] Reference flow [m^3/s] [or AUTO_COMPUTE]
/// \param gaged      [in] If true, hydrographs are generated
/// \param is_conduit [in] if true, this is a conduit
//
CSubBasin::CSubBasin( const long           Identifier,
                      const string         Name,
                      const CModelABC     *pMod,
                      const long           down_ID,     //index of downstream SB, if <0, downstream outlet
                      const CChannelXSect *pChan,       //Channel
                      const double         reach_len,   //reach length [m]
                      const double         Qreference,  //reference flow, m3/s [or AUTO_COMPUTE]
                      const bool           gaged,
                      const bool           is_conduit)
{
  ExitGracefullyIf(pMod==NULL,
    "CSubBasin:Constructor: NULL model",BAD_DATA);
  ExitGracefullyIf(((Qreference<=0.0) && (Qreference!=AUTO_COMPUTE)),
    "CSubBasin::Constructor: Reference flow must be non-zero and positive (or _AUTO)",BAD_DATA);

  _pModel            =pMod;

  _ID                =Identifier;
  _name              =Name;
  _is_conduit        =is_conduit;
  _global_p          =DOESNT_EXIST;

  _basin_area        =0.0;
  _drainage_area     =0.0;
  _reach_length      =reach_len;
  _reach_length2     =reach_len;
  _is_headwater      =true;
  _reach_HRUindex    =DOESNT_EXIST; //default
  _hyporheic_flux    =0.0; //default
  _convect_coeff     =2.0; //default
  _sens_exch_coeff   =0.0;
  _GW_exch_coeff     =0.0;
  _bed_conductivity  =0.0;
  _bed_thickness     =0.5; //m

  _t_conc            =AUTO_COMPUTE;
  _t_peak            =AUTO_COMPUTE;
  _t_lag             =AUTO_COMPUTE;
  _gamma_shape       =3.0;
  _gamma_scale       =1.0;
  _reservoir_constant=AUTO_COMPUTE;
  _num_reservoirs    =1;

  _rain_corr         =1.0;
  _snow_corr         =1.0;
  _temperature_corr = 0.0;
  _unusable_flow_pct =0.0;

  _res_disabled      =false;
  _assimilate        =false;

  // estimate reach length and _nSegments if needed
  //-----------------------------------------------------------------------
  double max_len = _pModel->GetGlobalParams()->GetParams()->max_reach_seglength*M_PER_KM;

  if((_reach_length==AUTO_COMPUTE) && (max_len/M_PER_KM<0.99*DEFAULT_MAX_REACHLENGTH))
  {
    string warn="Basin "+to_string(_ID)+" with _AUTO reach length will be assigned a single reach segment; to override, specify a specific reach length";
    WriteWarning(warn.c_str(),BAD_DATA_WARN);
    _nSegments = 1;
  }
  else{
    _nSegments = max((int)(ceil(_reach_length/max_len)),1);
  }

  string warn="CSubBasin:Constructor: exceeded maximum river segments in basin " +to_string(_ID)+": MAX_REACH_SEGLENGTH may be too small";
  if(_nSegments>MAX_RIVER_SEGS){ WriteWarning(warn,true); _nSegments=MAX_RIVER_SEGS-1; }

  _downstream_ID     =down_ID;
  if (_downstream_ID<0){_downstream_ID=DOESNT_EXIST;}//outflow basin
  _gauged            =gaged;
  _disabled          =false;

  _pDownSB           =NULL;
  _pChannel          =pChan; //Can be NULL
  _pReservoir        =NULL;
  _nHydroUnits       =0;
  _pHydroUnits       =NULL;

  //Initialized in Initialize
  _aQout=NULL;
  _aQout=new double [_nSegments];
  ExitGracefullyIf(_aQout==NULL,"CSubBasin Constructor",OUT_OF_MEMORY);
  for (int seg=0;seg<_nSegments;seg++){
    _aQout[seg]=AUTO_COMPUTE; //can be overridden by initial conditions
  }
  _QoutLast=AUTO_COMPUTE; //can be overridden by initial conditions
  _QlatLast=AUTO_COMPUTE; //can be overridden by initial conditions
  _Qlocal=0;
  _QlocLast=0;
  _channel_storage=0.0;   //calculated from initial flows
  _rivulet_storage=0.0;   //calculated from initial flows

  _Qirr=0.0;
  _QirrLast=0.0;
  _Qdelivered=0.0;
  _Qreturn=0.0;
  _QdivLast=0.0;
  _Qdiverted=0;
  _aQdelivered=NULL;
  _aQreturned=NULL;

  //Below are initialized in GenerateCatchmentHydrograph, GenerateRoutingHydrograph
  _aQlatHist     =NULL;  _nQlatHist     =0;
  _aQinHist      =NULL;  _nQinHist      =0;
  _aUnitHydro    =NULL;
  _aRouteHydro   =NULL;
  _c_hist        =NULL;

  //Below are modified using AddDiversion()
  _nDiversions   =0;
  _pDiversions   =NULL;

  //initialized in AddInflowHydrograph, similar
  _pInflowHydro  =NULL;
  _pInflowHydro2 =NULL;
  _pEnviroMinFlow=NULL;

  _pWaterDemands =NULL;
  _nWaterDemands =0;

  _Q_ref      =Qreference;
  _c_ref      =AUTO_COMPUTE;
  _w_ref      =AUTO_COMPUTE;
  _mannings_n =AUTO_COMPUTE;
  _slope      =AUTO_COMPUTE;
  _diffusivity=AUTO_COMPUTE;
  if(pChan!=NULL) {
    _mannings_n=_pChannel->GetMinMannings();
    _slope     =_pChannel->GetBedslope();
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
/// \details Point array references to null
//
CSubBasin::~CSubBasin()
{
  if (DESTRUCTOR_DEBUG){cout<<"  DELETING SUBBASIN"<<endl;}
  delete [] _pHydroUnits;_pHydroUnits=NULL; //just deletes pointer array, not hydrounits
  delete [] _aQout;      _aQout      =NULL;
  delete [] _aQlatHist;  _aQlatHist  =NULL;
  delete [] _aQinHist;   _aQinHist   =NULL;
  delete [] _aUnitHydro; _aUnitHydro =NULL;
  delete [] _aRouteHydro;_aRouteHydro=NULL;
  delete [] _c_hist;     _c_hist     =NULL;
  delete [] _pDiversions;_pDiversions=NULL; _nDiversions=0;
  delete _pInflowHydro;  _pInflowHydro=NULL;
  delete _pInflowHydro2; _pInflowHydro2=NULL;
  if (_pWaterDemands!=NULL){
    for (int ii = 0; ii < _nWaterDemands; ii++) {delete _pWaterDemands[ii];} delete [] _pWaterDemands;  _pWaterDemands=NULL;
  }
  delete [] _aQdelivered;
  delete [] _aQreturned;
  delete _pEnviroMinFlow;_pEnviroMinFlow=NULL;
  delete _pReservoir;
}
/*****************************************************************
   Accessors
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Returns subbasin ID
/// \return subbasin ID
//
long                    CSubBasin::GetID               () const {return _ID;            }

//////////////////////////////////////////////////////////////////
/// \brief Returns global index p
/// \return global index p
//
int                     CSubBasin::GetGlobalIndex      () const {return _global_p;      }

//////////////////////////////////////////////////////////////////
/// \brief Returns subbasin name
/// \return subbasin name
//
string                  CSubBasin::GetName             () const {return _name;          }

//////////////////////////////////////////////////////////////////
/// \brief Returns subbasin physical area [km^2]
/// \return subbasin area [km^2]
//
double                  CSubBasin::GetBasinArea        () const {return _basin_area;    }//[km2]

//////////////////////////////////////////////////////////////////
/// \brief Returns SB drainage area (includes basin and upstream contributing areas) [km^2]
/// \return SB drainage area [km^2]
//
double                  CSubBasin::GetDrainageArea     () const {return _drainage_area; }//[km2]

//////////////////////////////////////////////////////////////////
/// \brief Returns ID of downstream subbasin
/// \return ID of downstream subbasin (or DOESNT_EXIST if there is none)
//
long                    CSubBasin::GetDownstreamID     () const {return _downstream_ID; }

//////////////////////////////////////////////////////////////////
/// \brief Returns reach length [m]
/// \return Reach length [m]
//
double                   CSubBasin::GetReachLength      () const {return _reach_length;  }//[m]

//////////////////////////////////////////////////////////////////
/// \brief Returns True if subbasin is gauged, false otherwise
/// \return Boolean indicating if subbasin is gauged
//
bool                     CSubBasin::IsGauged            () const {return _gauged;        }

//////////////////////////////////////////////////////////////////
/// \brief Returns True if subbasin is headwater, false otherwise
/// \return Boolean indicating if subbasin is headwater basin
//
bool                     CSubBasin::IsHeadwater         () const {return _is_headwater;  }

//////////////////////////////////////////////////////////////////
/// \brief Returns number of river segments used in routing
/// \return Number of river segments used in routing
//
int                 CSubBasin::GetNumSegments          () const {return _nSegments;     }

//////////////////////////////////////////////////////////////////
/// \brief Returns true if basin is enabled
/// \return true if basin is enabled
//
bool                 CSubBasin::IsEnabled              () const {return !_disabled;     }

//////////////////////////////////////////////////////////////////
/// \brief Returns rainfall correction factor for basin
/// \return rainfall correction factor for basin
//
double            CSubBasin::GetRainCorrection         () const {return _rain_corr;     }

//////////////////////////////////////////////////////////////////
/// \brief Returns snowfall correction factor for basin
/// \return snowfall correction factor for basin
//
double            CSubBasin::GetSnowCorrection         () const{ return _snow_corr;     }

//////////////////////////////////////////////////////////////////
/// \brief Returns temperature correction factor for basin
/// \return temperature correction factor for basin
//
double            CSubBasin::GetTemperatureCorrection  () const { return _temperature_corr; }

//////////////////////////////////////////////////////////////////
/// \brief returns Unit Hydrograph as array pointer
/// \return Unit Hydrograph as array pointer
//
const double        *CSubBasin::GetUnitHydrograph    () const{return _aUnitHydro;}

//////////////////////////////////////////////////////////////////
/// \brief returns routing Hydrograph as array pointer
/// \return routing Hydrograph as array pointer
//
const double        *CSubBasin::GetRoutingHydrograph () const{return _aRouteHydro;}

//////////////////////////////////////////////////////////////////
/// \brief returns historical inflow hydrograph as array pointer
/// \return historical inflow hydrograph as array pointer
//
const double        *CSubBasin::GetInflowHistory     () const{return _aQinHist;}

//////////////////////////////////////////////////////////////////
/// \brief returns outflow array
/// \return outflow array as array pointer
//
const double        *CSubBasin::GetOutflowArray     () const {return _aQout;}

//////////////////////////////////////////////////////////////////
/// \brief returns historical inflow hydrograph as array pointer
/// \return historical inflow hydrograph as array pointer
//
const double        *CSubBasin::GetLatHistory       () const{return _aQlatHist;}

//////////////////////////////////////////////////////////////////
/// \brief returns number of timesteps stored in unit hydrograph history
/// \return number of timesteps stored in unit hydrograph history
//
int                  CSubBasin::GetLatHistorySize    () const{return _nQlatHist;}

//////////////////////////////////////////////////////////////////
/// \brief returns number of timesteps stored in routing hydrograph history
/// \return number of timesteps stored in routing hydrograph history
//
int                  CSubBasin::GetInflowHistorySize () const{return _nQinHist;}

//////////////////////////////////////////////////////////////////
/// \brief returns number of timesteps stored in routing hydrograph history
/// \return number of timesteps stored in routing hydrograph history
//
int                  CSubBasin::GetOutflowArraySize () const{return _nSegments;}

//////////////////////////////////////////////////////////////////
/// \brief returns true if subbasin outflow data should be used in assimilation
/// \return true if subbasin outflow data should be used in assimilation
//
bool                 CSubBasin::UseInFlowAssimilation() const { return _assimilate; }

//////////////////////////////////////////////////////////////////
/// \brief Returns Number of HRUs in SB
/// \return Number of HRUs in SB
//
int                CSubBasin::GetNumHRUs             () const{return _nHydroUnits;   }

//////////////////////////////////////////////////////////////////
/// \brief Returns HRU corresponding to index ks
/// \param ks [in] Index correspondng to the HRU of interest (local index, not global)
/// \return HRU corresponding to index ks
//
const CHydroUnit*CSubBasin::GetHRU(const int ks) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((ks<0) && (ks>=_nHydroUnits),"CSubBasin:GetHRU::improper index",BAD_DATA);
#endif
  return _pHydroUnits[ks];
}
//////////////////////////////////////////////////////////////////
/// \brief Returns pointer to reservoir associated with Subbasin (or NULL)
/// \return  pointer to reservoir associated with Subbasin (or NULL if no reservoir)
//
CReservoir    *CSubBasin::GetReservoir () const
{
  return _pReservoir;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns area-weighted average value of state variable with index i over all HRUs
/// \param i [in] Index corresponding to a state variable
/// \return area-weighted average value of state variable with index i over all HRUs
//
double CSubBasin::GetAvgStateVar (const int i) const
{
#ifdef  _STRICTCHECK_
  ExitGracefullyIf((i<0) && (i>=_pModel->GetNumStateVars()),"CSubBasin:GetAverageStateVar::improper index",BAD_DATA);
#endif
  double sum=0.0;
  for (int k=0;k<_nHydroUnits;k++)
  {
    if(_pHydroUnits[k]->IsEnabled())
    {
      sum+=_pHydroUnits[k]->GetStateVarValue(i)*_pHydroUnits[k]->GetArea();
    }
  }
  if (_basin_area==0.0){return 0.0;}
  return sum/_basin_area;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns area-weighted average value of concentration/temperature with index i over all HRUs
/// \param i [in] Index corresponding to a state variable
/// \return area-weighted average value of concentration/temperature with index i over all HRUs
//
double CSubBasin::GetAvgConcentration(const int i) const
{
#ifdef  _STRICTCHECK_
  ExitGracefullyIf((i<0) && (i>=_pModel->GetNumStateVars()),"CSubBasin:GetAvgConcentration::improper index",BAD_DATA);
#endif
  double sum=0.0;
  for(int k=0;k<_nHydroUnits;k++)
  {
    if(_pHydroUnits[k]->IsEnabled())
    {
      sum+=_pHydroUnits[k]->GetConcentration(i)*_pHydroUnits[k]->GetArea();
    }
  }
  if (_basin_area==0.0){return 0.0;}
  return sum/_basin_area;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns area-weighted average value of forcing function over entire subbasin
/// \param &ftype [in] forcing type identifier corresponding to a forcing function
/// \return area-weighted average value of forcing function over entire subbasin
//
double CSubBasin::GetAvgForcing (const forcing_type &ftype) const
{
  double sum=0.0;
  for (int k=0;k<_nHydroUnits;k++)
  {
    if(_pHydroUnits[k]->IsEnabled())
    {
      sum    +=_pHydroUnits[k]->GetForcing(ftype)*_pHydroUnits[k]->GetArea();
    }
  }
  if (_basin_area==0.0){return 0.0;}
  return sum/_basin_area;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns area-weighted average of specified cumulative flux over HRU group
///
/// \param i [in] index of storage compartment
/// \param to [in] true if evaluating cumulative flux to storage compartment, false for 'from'
/// \return Area-weighted average of cumulative flux to storage compartment i
//
double CSubBasin::GetAvgCumulFlux(const int i,const bool to) const
{
  //Area-weighted average
  double sum=0.0;
  for(int k=0;k<_nHydroUnits;k++)
  {
    if(_pHydroUnits[k]->IsEnabled())
    {
      sum    +=_pHydroUnits[k]->GetCumulFlux(i,to)*_pHydroUnits[k]->GetArea();
    }
  }
  if (_basin_area==0.0){return 0.0;}
  return sum/_basin_area;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns area-weighted average of  cumulative flux between two compartments over subbasin
///
/// \param iFrom [in] index of 'from' storage compartment
/// \param iTo [in] index of 'to' storage compartment
/// \return Area-weighted average of cumulative flux between two compartments over subbasin
//
double CSubBasin::GetAvgCumulFluxBet(const int iFrom,const int iTo) const
{
  //Area-weighted average
  double sum=0.0;
  for(int k=0;k<_nHydroUnits;k++)
  {
    if(_pHydroUnits[k]->IsEnabled())
    {
      sum    +=_pHydroUnits[k]->GetCumulFluxBet(iFrom,iTo)*_pHydroUnits[k]->GetArea();
    }
  }
  if (_basin_area==0.0){return 0.0;}
  return sum/_basin_area;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns specified inflow to subbasin at time t
/// \param &t [in] Model time at which the inflow to SB is to be determined
/// \return specified inflow to subbasin at time t
//
double CSubBasin::GetSpecifiedInflow(const double &t) const
{
  if (_pInflowHydro==NULL){return 0.0;}
  return _pInflowHydro->GetValue(t);
}
//////////////////////////////////////////////////////////////////
/// \brief Returns specified downstream inflow to subbasin at time t
/// \param &t [in] Model time at which the inflow to SB is to be determined
/// \return specified inflow to subbasin at time t
//
double CSubBasin::GetDownstreamInflow(const double &t) const
{
  if (_pInflowHydro2==NULL){return 0.0;}
  return _pInflowHydro2->GetValue(t);
}
//////////////////////////////////////////////////////////////////
/// \brief Returns specified irrigation/water use demand from subbasin at time t
/// \param &t [in] Model time at which the demand from SB is to be determined
/// \return specified demand from subbasin at time t
//
double CSubBasin::GetTotalWaterDemand() const
{
  double sum=0.0;
  for (int ii=0;ii<_nWaterDemands;ii++)
  {
    sum+=_pWaterDemands[ii]->GetDemand();
  }
  return sum;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns specified irrigation/water use demand from subbasin
/// \param &t [in] Model time at which the demand from SB is to be determined
/// \return specified demand from subbasin at current time
//
double CSubBasin::GetWaterDemand(const int ii) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf(ii < 0 || ii >= _nWaterDemands, "CSubBasin::GetWaterDemand: invalid index",RUNTIME_ERR);
#endif

  return _pWaterDemands[ii]->GetDemand();
}
//////////////////////////////////////////////////////////////////
/// \brief Returns instantaneous ACTUAL water use at end of current timestep
/// \return actual delivery from subbasin [m3/s]
//
double CSubBasin::GetDemandDelivery() const
{
  double sum=0;
  for (int ii=0;ii<_nWaterDemands;ii++){
    sum+=_aQdelivered[ii];
  }
  return sum;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns calculated irrigation/water use return flow to subbasin at time t
/// \param &t [in] Model time at which the return flow from SB is to be determined
/// \return specified return to subbasin at time t [m3/s]
//
double CSubBasin::GetTotalReturnFlow() const
{
  return _Qreturn;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns instantaneous ACTUAL water delivery at end of current timestep
/// \return actual delivery from subbasin [m3/s]
//
double CSubBasin::GetDemandDelivery(const int ii) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((ii<0) || (ii>=_nWaterDemands),"CSubBasin::GetDemandDelivery: invalid demand index",RUNTIME_ERR);
#endif
  return _aQdelivered[ii];
}
//////////////////////////////////////////////////////////////////
/// \brief Returns instantaneous ACTUAL return flow from water use ii at end of current timestep (may be going to another basin)
/// \return actual return from subbasin demand ii [m3/s]
//
double CSubBasin::GetReturnFlow(const int ii) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((ii<0) || (ii>=_nWaterDemands),"CSubBasin::GetDemandDelivery: invalid demand index",RUNTIME_ERR);
#endif
  return _aQreturned[ii];
}
//////////////////////////////////////////////////////////////////
/// \brief Returns cumulative downstream specified water demand, including from this subbasin
/// \param &t [in] Model time at which the demand from SB is to be determined
/// \return specified cumulative downstream demand (in [m3/s]) from subbasin at time t
//
double CSubBasin::GetDownstreamIrrDemand(const double &t) const
{
  double Qirr;
  Qirr=GetTotalWaterDemand();

  if (_downstream_ID==DOESNT_EXIST){ return Qirr; }
  else                             { return _pDownSB->GetDownstreamIrrDemand(t)+Qirr; }
}
//////////////////////////////////////////////////////////////////
/// \brief Returns specified environmental minimum flow at subbasin outlet at time t
/// \param &t [in] Model time at which the environmental minimum flow is to be determined
/// \return specified environmental minimum flow in subbasin at time t
//
double CSubBasin::GetEnviroMinFlow(const double &t) const
{
  if(_pEnviroMinFlow==NULL) { return 0.0; }
  double Qmin=_pEnviroMinFlow->GetValue(t); if(Qmin==RAV_BLANK_DATA) { Qmin=0; }
  return Qmin;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns subbasin index p of diversion target basin
/// \return subbasin index p of diversion target basin
//
int CSubBasin::GetDiversionTargetIndex(const int i) const
{
#ifdef _STRICTCHECK_
  if ((i<0) || (i>=_nDiversions)){ExitGracefully("GetDiversionTargetIndex: bad index",RUNTIME_ERR);}
#endif
  return _pDiversions[i]->target_p;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns true if subbasin has irrigation demand
/// \return true if subbasin has irrigation demand
//
bool CSubBasin::HasIrrigationDemand() const
{
  int nResDemands=0;
  if (_pReservoir!=NULL){nResDemands=_pReservoir->GetNumWaterDemands();}
  return ((_nWaterDemands+nResDemands)>0);
}
//////////////////////////////////////////////////////////////////
/// \brief Returns true if subbasin has irrigation demand
/// \return true if subbasin has irrigation demand
//
bool CSubBasin::HasEnviroMinFlow() const
{
  return (_pEnviroMinFlow!=NULL);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns downstream outflow due to irrigation demand after flow constraints applied
/// \param &t [in] Model time at which the outflow from SB is to be determined
/// \param &Q [in] Estimate of subbasin outflow prior to applying demand [m3/s]
/// \return irrigation demand outflow from subbasin at time t [m3/s]
//
double CSubBasin::ApplyIrrigationDemand(const double &t,const double &Q,const bool optimized) const
{
  if (_nWaterDemands==0){return 0.0;}

  //if using Management optimization, delivery determined by optimizer
  if (optimized) {
    //cout<<setprecision(5)<< _Qdelivered<<" of "<<GetTotalWaterDemand()<<" applied ("<<setprecision(3)<<_Qdelivered/GetTotalWaterDemand(t)*100<<"%) ("<<_Qdelivered/Q*100<<"% of flow)" << endl;
    return min(Q,_Qdelivered); //\todo[fix] - this min shouldnt be required.
    return _Qdelivered; //total delivered as calculated from management optimization
  }

  double Qdelivered;
  double Qdemand=GetTotalWaterDemand();
  double Qmin   =GetEnviroMinFlow(t);

  Qdelivered=min(max((1.0-_unusable_flow_pct)*(Q-Qmin),0.0),Qdemand);

  return Qdelivered;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns number of flow diversions
/// \return Returns number of flow diversions
//
int CSubBasin::GetNumDiversions() const {
  return _nDiversions;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns diverted flow rate to target subbasin with index pDivert
/// \param &i [in] diversion index
/// \param &Q [in] Estimate of subbasin outflow prior to applying diversion [m3/s]
/// \param &Options [in] model options structure
/// \param &tt [in] current model time structure
/// \param &pDivert [out] index of subbasin to divert to
/// \return irrigation demand outflow from subbasin at time t [m3/s]
//
double CSubBasin::GetDiversionFlow(const int i, const double &Q, const optStruct &Options, const time_struct &tt, int &p_Divert) const
{
  if ((i<0) || (i>=_nDiversions)){ExitGracefully("CSubBasin::GetDiversionFlow",RUNTIME_ERR); }
  if (_pDiversions[i]==NULL) { p_Divert=DOESNT_EXIST; return 0.0; } //no diversion

  p_Divert=_pDiversions[i]->target_p;

  if (IsInDateRange(tt.julian_day,_pDiversions[i]->julian_start,_pDiversions[i]->julian_end))
  {
    if(_pDiversions[i]->aQdivert==NULL) {
      if(Q>_pDiversions[i]->min_flow)
      {
        return (Q-_pDiversions[i]->min_flow)*_pDiversions[i]->percentage;
      }
    }
    else { //uses lookup table
      return InterpolateCurve(Q,_pDiversions[i]->aQsource,_pDiversions[i]->aQdivert,_pDiversions[i]->nPoints,false);
    }
  }

  return 0.0;
}
//////////////////////////////////////////////////////////////////
/// \brief add diversion
/// \param jul_start [in] julian start date of diversion (must be 0-365)
/// \param jul_end [in] julian end date of diversion (must be 0-365, >start date)
/// \param target_p [in] index (NOT SBID) of target subbasin
/// \param min_flow [in] minimum flow for demand to be applied
/// \param pct [in] percent of flow diverted to other subbasin (0-1.0)
/// \return irrigation demand outflow from subbasin at time t [m3/s]
//
void  CSubBasin::AddFlowDiversion(const int jul_start,const int jul_end,const int target_p,const double min_flow,const double pct) {

  if((jul_start<0) || (jul_start>366)) {
    ExitGracefully("CSubBasin::AddFlowDiversion: invalid julian start date. must be between 0 and 366",BAD_DATA_WARN);
  }
  if((jul_end<0) || (jul_end>366)) {
    ExitGracefully("CSubBasin::AddFlowDiversion: invalid julian start date. must be between 0 and 366",BAD_DATA_WARN);
  }
  if(min_flow<0) {
    WriteWarning("CSubBasin::AddFlowDiversion: negative minimum flow. will be corrected to 0.0",true);
  }
  if((pct<0) || (pct>100)){
    WriteWarning("CSubBasin::AddFlowDiversion: invalid diversion fraction. Must be between 0.0 (0%) and 1.0 (100%)",true);
  }
  if(jul_start>jul_end) {
    WriteWarning("CSubBasin::AddFlowDiversion: julian start date after julian end date. No diversion will occur",true);
  }
  if(target_p==DOESNT_EXIST) {
    WriteWarning("CSubBasin::AddFlowDiversion: target diversion subbasin is -1. Watershed mass balance reporting does not yet incorporate this loss, and will incorrectly indicate a mass balance error.",true);
  }
  diversion *div=new diversion();
  div->julian_end=jul_end;
  div->julian_start=jul_start;
  div->min_flow=max(min_flow,0.0);
  div->percentage=min(max(pct,0.0),1.0);
  div->target_p=target_p;

  if(!DynArrayAppend((void**&)(_pDiversions),(void*)(div),_nDiversions)) {
    ExitGracefully("CSubBasin::AddFlowDiversion: adding NULL diversion",BAD_DATA);
  }
}
void  CSubBasin::AddFlowDiversion(const int jul_start,const int jul_end,const int target_p,const double *aQ1,const double *aQ2,const int NQ) {
  if((jul_start<0) || (jul_start>366)) {
    ExitGracefully("CSubBasin::AddFlowDiversion: invalid julian start date. must be between 0 and 366",BAD_DATA_WARN);
  }
  if((jul_end<0) || (jul_end>366)) {
    ExitGracefully("CSubBasin::AddFlowDiversion: invalid julian start date. must be between 0 and 366",BAD_DATA_WARN);
  }
  if(jul_start>jul_end) {
    WriteWarning("CSubBasin::AddFlowDiversion: julian start date after julian end date. No diversion will occur",true);
  }
  if(target_p==DOESNT_EXIST) {
    WriteWarning("CSubBasin::AddFlowDiversion: target diversion subbasin is -1. Watershed mass balance reporting does not yet incorporate this loss, and will incorrectly indicate a mass balance error.",true);
  }
  diversion *div=new diversion();
  div->julian_end=jul_end;
  div->julian_start=jul_start;
  div->target_p=target_p;

  div->min_flow=RAV_BLANK_DATA;
  div->percentage=RAV_BLANK_DATA;
  div->nPoints=NQ;
  div->aQsource=new double [NQ];
  div->aQdivert=new double [NQ];
  ExitGracefullyIf(div->aQdivert==NULL,"CSubBasin::AddFlowDiversion",OUT_OF_MEMORY);
  for(int i=0;i<NQ;i++) {
    div->aQsource[i]=aQ1[i];
    div->aQdivert[i]=aQ2[i];
  }

  if(!DynArrayAppend((void**&)(_pDiversions),(void*)(div),_nDiversions)) {
    ExitGracefully("CSubBasin::AddFlowDiversion: adding NULL diversion",BAD_DATA);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Returns reach HRU index
/// \return reach HRU index
//
int   CSubBasin::GetReachHRUIndex() const {
  return _reach_HRUindex;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns reach convection coefficient
/// \return reach convection coefficient
//
double CSubBasin::GetConvectionCoeff() const {
  return _convect_coeff;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns reach hyporheix gross flux
/// \return reach hyporheic flux
//
double CSubBasin::GetHyporheicFlux() const {
  return _hyporheic_flux;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns reach bed conductivity
/// \return reach hyporheic bed conductivity
//
double CSubBasin::GetRiverbedConductivity() const {
  return _bed_conductivity;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns reach bed thickness
/// \return reach bed thickness
//
double CSubBasin::GetRiverbedThickness() const {
  return _bed_thickness;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns channel storage [m^3]
/// \note Should only be called after _aQinHist has been updated by calling UpdateInflow
/// \return Channel storage [m^3]
//
double CSubBasin::GetChannelStorage () const
{
  return _channel_storage;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns channel storage [m^3]
/// \note Should only be called after _aQinHist has been updated by calling UpdateInflow
/// \return Channel storage [m^3]
//
double CSubBasin::GetReservoirStorage () const
{
  if(_pReservoir==NULL){return 0.0;}
  return _pReservoir->GetStorage();
}
//////////////////////////////////////////////////////////////////
/// \brief Returns rivulet storage [m^3]
/// \details Returns uphill tributary storage [m^3] obtained from convoltuion
/// of inflow history with unit hydrograph
/// \note Should only be called after _aQlatHist has been updated by calling UpdateLateralInflow
///
/// \return Current Rivulet (non-main channel) surface water storage
//
double CSubBasin::GetRivuletStorage () const
{
  return _rivulet_storage;
}
//////////////////////////////////////////////////////////////////
/// \brief number of water/irrigation demands
/// \return number of water/irrigation demands
//
int CSubBasin::GetNumWaterDemands() const
{
  return _nWaterDemands;
}
//////////////////////////////////////////////////////////////////
/// \brief returns water/irrigation demand integer ID
/// \return water/irrigation demand integer ID
//
CDemand *CSubBasin::GetWaterDemandObj (const int ii) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf(ii < 0 || ii >= _nWaterDemands, "CSubBasin::GetWaterDemandObj: invalid index",RUNTIME_ERR);
#endif
  return _pWaterDemands[ii];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns Outflow from subbasin at start of current timestep (during solution) or end of completed timestep [m^3/s]
/// \return Outflow at start of current timestep (during solution) or end of completed timestep [m^3/s]
//
double    CSubBasin::GetOutflowRate   () const
{
  if (_pReservoir!=NULL){return _pReservoir->GetOutflowRate();}
  return _aQout[_nSegments-1]; //[m3/s](from start of time step until after solver is called)
}

//////////////////////////////////////////////////////////////////
/// \brief Returns Outflow from channel reach (even if reservoir present) before diversions at start of current timestep (during solution) or end of completed timestep [m^3/s]
/// \note used for calculating flow diversions
/// \return Channel outflow at start of current timestep (during solution) or end of completed timestep [m^3/s]
//
double    CSubBasin::GetChannelOutflowRate   () const
{
  return _aQout[_nSegments-1]+_Qdiverted; //[m3/s](from start of time step until after solver is called)
}
//////////////////////////////////////////////////////////////////
/// \brief Returns Outflow from channel reach (even if reservoir present) before diversions at start of current timestep (during solution) or end of completed timestep [m^3/s]
/// \note used for calculating flow diversions
/// \return Channel outflow at start of current timestep (during solution) or end of completed timestep [m^3/s]
//
double    CSubBasin::GetLastChannelOutflowRate   () const
{
  return _QoutLast+_QdivLast; //[m3/s](from start of time step until after solver is called)
}



//////////////////////////////////////////////////////////////////
/// \brief Returns Outflow *from channel* at start of previous timestep (during solution) or start of completed timestep [m^3/s]
/// \return Outflow at start of previous timestep (during solution) or start of completed timestep [m^3/s]
//
double    CSubBasin::GetLastOutflowRate() const
{
  return _QoutLast;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns local (in-catchment) contribution to outflow at start of current timestep (during solution) or end of completed timestep [m^3/s]
/// \return local outflow at start of current timestep (during solution) or end of completed timestep [m^3/s]
//
double CSubBasin::GetLocalOutflowRate() const
{
  return _Qlocal;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns Outflow at start of current timestep (during solution) or end of completed timestep [m^3/s]
/// \return Outflow at start of current timestep (during solution) or end of completed timestep [m^3/s]
//
double  CSubBasin::GetReservoirInflow() const
{
  if (_pReservoir==NULL){return 0.0;}
  return _aQout[_nSegments-1]; //[m3/s](from start of time step until after solver is called)
}
//////////////////////////////////////////////////////////////////
/// \brief Returns Reservoir losses integrated over specified timestep [m^3]
/// \return Reservoir losses over specified timestep  [m^3]
//
double CSubBasin::GetReservoirLosses(const double &tstep) const
{
  if(_pReservoir==NULL){ return 0.0; }
  return _pReservoir->GetReservoirLosses(tstep);
}
//////////////////////////////////////////////////////////////////
/// \brief Returns Irrigation losses integrated over specified timestep [m^3]
/// \return Irrigation losses over specified timestep  [m^3]
//
double CSubBasin::GetIrrigationLosses(const double &tstep) const
{
  return 0.5*(_QirrLast+_Qirr)*tstep*SEC_PER_DAY;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns Irrigation losses integrated over specified timestep [m^3]
/// \return Irrigation losses over specified timestep  [m^3]
//
double CSubBasin::GetDiversionLosses(const double &tstep) const
{
  return 0.5*(_Qdiverted+_QdivLast)*tstep*SEC_PER_DAY;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns total volume lost from main reach over timestep [m^3]
/// \note Should be called only at end of completed tstep
/// used in mass balance to estimate water loss through streamflow
/// \return Total volume lost from main reach over timestep [m^3]
//
double CSubBasin::GetIntegratedOutflow(const double &tstep) const//[m3]
{
  if (_pReservoir!=NULL){return _pReservoir->GetIntegratedOutflow(tstep);}

  return 0.5*(_aQout[_nSegments-1]+_QoutLast)*(tstep*SEC_PER_DAY); //integrated
}
//////////////////////////////////////////////////////////////////
/// \brief Returns total volume lost from main reach to reservoir over timestep [m^3]
/// \note Should be called only at end of completed tstep
/// \return Total volume lost from main reach over timestep [m^3]
//
double CSubBasin::GetIntegratedReservoirInflow(const double &tstep) const
{
  if (_pReservoir==NULL){return 0.0;}
  return 0.5*(_aQout[_nSegments-1]+_QoutLast)*(tstep*SEC_PER_DAY); //integrated
}
//////////////////////////////////////////////////////////////////
/// \brief Returns total volume gained by main reach from unmodeled sources over timestep [m^3]
/// \note Should be called only at end of completed tstep
/// \return Total volume gained by main reach from unmodeled sources over timestep [m^3]
//
double CSubBasin::GetIntegratedSpecInflow(const double &t, const double &tstep) const//[m3]
{
  //used in mass balance to estimate water gain from unmodeled upstream sources
  double sum=0.0;
  sum+=0.5*(GetSpecifiedInflow(t) +GetSpecifiedInflow (t+tstep))*(tstep*SEC_PER_DAY); //integrated
  sum+=GetDownstreamInflow(t)*(tstep*SEC_PER_DAY);                   //integrated -period starting
  return sum;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns total volume lost by main reach sourced from in-catcment routing over timestep [m^3]
/// \note Should be called only at end of completed tstep
/// \return Total volume lost by main reach sourced from in-catcment routing over timestep [m^3]
//
double CSubBasin::GetIntegratedLocalOutflow(const double &tstep) const//[m3]
{
  return 0.5*(_Qlocal+_QlocLast)*(tstep*SEC_PER_DAY); //integrated
}
//////////////////////////////////////////////////////////////////
/// \brief Returns reference discharge [m^3/s]
/// \return  reference discharge [m^3/s] or AUTO_COMPUTE if not yet calculated
//
double CSubBasin::GetReferenceFlow() const
{
  return _Q_ref;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns x-sect area at reference discharge [m^2]
/// \return  x-sect area at reference discharge [m^2] or AUTO_COMPUTE if not yet calculated
//
double CSubBasin::GetReferenceXSectArea() const
{
  if (_Q_ref==AUTO_COMPUTE){return AUTO_COMPUTE;}
  if(_pChannel==NULL)      {return ALMOST_INF;  }
  return _pChannel->GetArea(_Q_ref,_slope,_mannings_n);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns reference celerity [m/s]
/// \return  reference celerity [m/s] or AUTO_COMPUTE if not yet calculated
//
double CSubBasin::GetReferenceCelerity() const {
  return _c_ref;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns current celerity [m/s]
/// \return  reference celerity [m/s]
///
double CSubBasin::GetCelerity() const
{
  if(_pChannel!=NULL) {
    return _pChannel->GetCelerity(_aQout[_nSegments-1],_slope,_mannings_n);
  }
  else {
    return 0.0;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Returns reference diffusivity [m2/s]
/// \return  reference diffusivity [m2/s] or AUTO_COMPUTE if not yet calculated
//
double CSubBasin::GetDiffusivity() const
{
  if((_pChannel!=NULL) && (_diffusivity==AUTO_COMPUTE)){
    return _pChannel->GetDiffusivity(_Q_ref,_slope,_mannings_n);
  }
  else {
    return _diffusivity;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Returns channel depth, in meters
/// \return  channel depth [m], with minimum depth of 1cm
//
double CSubBasin::GetRiverDepth() const
{
  if (_pChannel==NULL){return ALMOST_INF;}
  const double MIN_CHANNEL_DEPTH=0.01;
  return max(_pChannel->GetDepth(_aQout[_nSegments-1],_slope,_mannings_n),MIN_CHANNEL_DEPTH);
}
//////////////////////////////////////////////////////////////////
/// \brief Returns channel area, in m2
/// \return  channel area [m], with minimum area of 0.01 m2
//
double CSubBasin::GetXSectArea() const
{
  if (_pChannel==NULL){return ALMOST_INF;}
  const double MIN_CHANNEL_AREA=0.01;
  return max(_pChannel->GetArea(_aQout[_nSegments-1],_slope,_mannings_n),MIN_CHANNEL_AREA);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns water surface elevation, in meters
/// \return  reference water surface elevation [m]
//
double CSubBasin::GetWaterLevel() const
{
  if (_pChannel==NULL){return ALMOST_INF;}
  const double MIN_CHANNEL_DEPTH=0.01;
  return _pChannel->GetStageElev(_aQout[_nSegments-1],_slope,_mannings_n);
}
//////////////////////////////////////////////////////////////////
/// \brief Returns reach bedslope
/// \return  reference reach bedslope [m/m]
//
double CSubBasin::GetBedslope() const{
  if(_pChannel==NULL) { return 0; }
  return _pChannel->GetBedslope();
}

//////////////////////////////////////////////////////////////////
/// \brief Returns current wetted perimeter of stream
/// \return  wetted perimeter [m]
//
double CSubBasin::GetWettedPerimeter() const
{
  if(_pChannel==NULL) { return ALMOST_INF; }
  return _pChannel->GetWettedPerim(_aQout[_nSegments-1],_slope,_mannings_n);
}
//////////////////////////////////////////////////////////////////
/// \brief Returns top surface width of stream
/// \return current top width [m]
//
double CSubBasin::GetTopWidth() const {
  if(_pChannel==NULL) { return ALMOST_INF; }
  return _pChannel->GetTopWidth(_aQout[_nSegments-1],_slope,_mannings_n);
}
//////////////////////////////////////////////////////////////////
/// \brief gets unusable flow percentatge
//
double CSubBasin::GetUnusableFlowPercentage() const {
  return _unusable_flow_pct;
}

/*****************************************************************
   Manipulators
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Adds an HRU to HRU coverage of subbasin
/// \param *pHRU [in] HRU to be added
//
void CSubBasin::AddHRU(CHydroUnit *pHRU)
{
  if (!DynArrayAppend((void**&)(_pHydroUnits),(void*)(pHRU),_nHydroUnits)){
    ExitGracefully("CSubBasin::AddHRU: adding NULL HRU",BAD_DATA);}
}

//////////////////////////////////////////////////////////////////
/// \brief Adds a reservoir to the outlet of the basin
/// \param *pRes [in] pointer to reservoir to be added
//
void CSubBasin::AddReservoir(CReservoir *pRes)
{
  ExitGracefullyIf(_pReservoir!=NULL,
                   "CSubBasin::AddReservoir: only one reservoir may be specified per basin",BAD_DATA);
  if(_res_disabled) {
    delete pRes;
    _reach_length=_reach_length2;
    return;
  }
  _pReservoir=pRes;
}

void CSubBasin::SetGlobalIndex(const int p)
{
  _global_p=p;
}
//////////////////////////////////////////////////////////////////
/// \brief Sets basin properties
/// \param label [in] String property identifier
/// \param &value [in] Property set value
/// \return True if procedure was successful
//
bool CSubBasin::SetBasinProperties(const string label,
                                   const double &value)
{
  string label_n = StringToUppercase(label);
  if      (!label_n.compare("TIME_CONC"     ))  {_t_conc=value;}
  else if (!label_n.compare("TIME_TO_PEAK"  ))  {_t_peak=value;}
  else if (!label_n.compare("TIME_LAG"      ))  {_t_lag=value;}
  else if (!label_n.compare("RES_CONSTANT"  ))  {_reservoir_constant=value;}
  else if (!label_n.compare("NUM_RESERVOIRS"))  {_num_reservoirs=(int)(value);}
  else if (!label_n.compare("GAMMA_SHAPE"   ))  {_gamma_shape=value;}
  else if (!label_n.compare("GAMMA_SCALE"   ))  {_gamma_scale=value;}

  else if (!label_n.compare("Q_REFERENCE"   ))  {_Q_ref=value;}
  else if (!label_n.compare("MANNINGS_N"    ))  {_mannings_n=value;}
  else if (!label_n.compare("SLOPE"         ))  {_slope=value;}
  else if (!label_n.compare("DIFFUSIVITY"   ))  {_diffusivity=value; }
  else if (!label_n.compare("CELERITY"      ))  {_c_ref=value; }

  else if (!label_n.compare("RAIN_CORR"     ))  {_rain_corr=value;}
  else if (!label_n.compare("SNOW_CORR"     ))  {_snow_corr=value;}
  else if (!label_n.compare("TEMP_CORR"     ))  {_temperature_corr=value; }

  else if (!label_n.compare("REACH_HRU_ID"  ))  { _reach_HRUindex=(int)(value); }
  else if (!label_n.compare("HYPORHEIC_FLUX"))  { _hyporheic_flux=value; }
  else if (!label_n.compare("CONVECT_COEFF" ))  { _convect_coeff=value; }
  else if (!label_n.compare("SENS_EXCH_COEFF")) { _sens_exch_coeff=value; }
  else if (!label_n.compare("GW_EXCH_COEFF" ))  { _GW_exch_coeff=value; }

  else if (!label_n.compare("RIVERBED_CONDUCTIVITY")){ _bed_conductivity = value; }
  else if (!label_n.compare("RIVERBED_THICKNESS"   )){ _bed_thickness = value; }
  else if (!label_n.compare("RESERVOIR_DISABLED"  )) { _res_disabled=(bool)(value); }
  else if (!label_n.compare("CORR_REACH_LENGTH"   )) { _reach_length2=value; }
  else if (!label_n.compare("LAKEBED_CONDUCTIVITY")) {
    if (_pReservoir != NULL) {_pReservoir->SetLakebedConductivity(value); }
  }
  else if (!label_n.compare("LAKEBED_THICKNESS")) {
    if (_pReservoir != NULL) {_pReservoir->SetLakebedThickness(value);}
  }
  else if (!label_n.compare("LAKE_CONVECT_COEFF")) {
    if (_pReservoir != NULL) {_pReservoir->SetLakeConvectionCoeff(value);}
  }
  else if (!label_n.compare("RESERVOIR_CREST_WIDTH")) {
    if (_pReservoir != NULL) { _pReservoir->SetCrestWidth(value);}
  }
  else{
    return false;//bad string
  }
  return true;
}
//////////////////////////////////////////////////////////////////
/// \brief Gets basin properties
/// \param label [in] String property identifier
/// \return value of basin property corresponding to label
//
double CSubBasin::GetBasinProperties(const string label) const
{
  string label_n = StringToUppercase(label);
  if      (!label_n.compare("TIME_CONC"     ))  {return _t_conc;}
  else if (!label_n.compare("TIME_TO_PEAK"  ))  {return _t_peak;}
  else if (!label_n.compare("TIME_LAG"      ))  {return _t_lag;}
  else if (!label_n.compare("RES_CONSTANT"  ))  {return _reservoir_constant;}
  else if (!label_n.compare("NUM_RESERVOIRS"))  {return (double)(_num_reservoirs);}
  else if (!label_n.compare("GAMMA_SHAPE"   ))  {return _gamma_shape;}
  else if (!label_n.compare("GAMMA_SCALE"   ))  {return _gamma_scale;}

  else if (!label_n.compare("Q_REFERENCE"   ))  {return _Q_ref;}
  else if (!label_n.compare("MANNINGS_N"    ))  {return _mannings_n;}
  else if (!label_n.compare("SLOPE"         ))  {return _slope;}
  else if (!label_n.compare("DIFFUSIVITY"   ))  {return _diffusivity; }
  else if (!label_n.compare("CELERITY"      ))  { return _c_ref; }

  else if (!label_n.compare("RAIN_CORR"     ))  { return _rain_corr;}
  else if (!label_n.compare("SNOW_CORR"     ))  { return _snow_corr;}
  else if (!label_n.compare("TEMP_CORR"     ))  { return _temperature_corr; }

  else if (!label_n.compare("REACH_HRU_ID"  ))  { return (double)(_reach_HRUindex); }
  else if (!label_n.compare("HYPORHEIC_FLUX"))  { return _hyporheic_flux; }
  else if (!label_n.compare("CONVECT_COEFF"))   { return _convect_coeff; }
  else if (!label_n.compare("SENS_EXCH_COEFF")) { return _sens_exch_coeff; }
  else if (!label_n.compare("GW_EXCH_COEFF"))   { return _GW_exch_coeff; }

  else if (!label_n.compare("RESERVOIR_DISABLED")) { return (double)(_res_disabled); }
  else if (!label_n.compare("CORR_REACH_LENGTH"))  { return _reach_length2; }

  else if (!label_n.compare("RESERVOIR_CREST_WIDTH")) {
    if(_pReservoir!=NULL) {
      return _pReservoir->GetCrestWidth();
    }
    else{return 0.0;}
  }
  else{
    return INDEX_NOT_FOUND;//bad string
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Sets basin headwater status (called by model during subbasin initialization)
//
void CSubBasin::SetAsNonHeadwater()
{
  _is_headwater=false;
}
//////////////////////////////////////////////////////////////////
/// \brief Adds inflow hydrograph
/// \param *pInflow pointer to Inflow time series to be added
//
void    CSubBasin::AddInflowHydrograph (CTimeSeries *pInflow)
{
  ExitGracefullyIf((_pInflowHydro!=NULL) && (g_current_e==DOESNT_EXIST),
                   "CSubBasin::AddInflowHydrograph: only one inflow hydrograph may be specified per basin",BAD_DATA);
  delete _pInflowHydro; // must delete previous if this is for an ensemble run
  _pInflowHydro=pInflow;
}
//////////////////////////////////////////////////////////////////
/// \brief Adds inflow hydrograph
/// \param *pInflow pointer to Inflow time series to be added
//
void    CSubBasin::AddDownstreamInflow (CTimeSeries *pInflow)
{
  ExitGracefullyIf((_pInflowHydro2!=NULL)  && (g_current_e==DOESNT_EXIST),
                   "CSubBasin::AddDownstreamInflow: only one inflow hydrograph may be specified per basin",BAD_DATA);
  delete _pInflowHydro; // must delete previous if this is for an ensemble run
  _pInflowHydro2=pInflow;
}
//////////////////////////////////////////////////////////////////
/// \brief Adds water demand
/// \param *pDemand pointer to demand object to be added
//
void    CSubBasin::AddWaterDemand(CDemand *pDemand)
{
  if (!DynArrayAppend((void**&)(_pWaterDemands),(void*)(pDemand),_nWaterDemands)){
    ExitGracefully("CSubBasin::AddIrrigationDemand: trying to add NULL regime",BAD_DATA);}
}
//////////////////////////////////////////////////////////////////
/// \brief Adds minimum enviro flow time series
/// \param *pInflow pointer to minimum enviro flow time series to be added
//
void    CSubBasin::AddEnviroMinFlow(CTimeSeries *pMinFlow)
{
  ExitGracefullyIf(_pEnviroMinFlow!=NULL,
    "CSubBasin::AddEnviroMinFlow: only one environmental minimum flow time series may be specified per basin",BAD_DATA);
  _pEnviroMinFlow=pMinFlow;
}

/////////////////////////////////////////////////////////////////
/// \brief Sets channel storage, usually upon read of state file
/// \param &V [in] channel storage [m3]
//
void CSubBasin::SetChannelStorage   (const double &V){
  _channel_storage=V;
}

/////////////////////////////////////////////////////////////////
/// \brief Sets rivulet storage, usually upon read of state file
/// \param &V [in] rivulet storage [m3]
//
void CSubBasin::SetRivuletStorage   (const double &V){
  _rivulet_storage=V;
}

/////////////////////////////////////////////////////////////////
/// \brief Sets qout storage array, usually upon read of state file
/// \param &N [in] size of aQo array (=_nSegments) or DOESNT_EXIST if _nSegments is unknown
/// \param *aQo array (size N) outflow from each reach segment [m3/s]
/// \param QoLast: most recent outflow from reach [m3/s]
//
void CSubBasin::SetQoutArray(const int N, const double *aQo, const double QoLast)
{
  if (N!=_nSegments){
    string warn="Number of reach segments in state file and input file are inconsistent in basin "+to_string(_ID)+". Unable to read in-reach flow initial conditions";
    WriteWarning(warn.c_str(),false);
  }
  else if (N==_nSegments){
    for (int i=0;i<_nSegments;i++){_aQout[i]=aQo[i];}
    _QoutLast=QoLast;
  }
}
/////////////////////////////////////////////////////////////////
/// \brief Sets qout storage array, usually upon read of state file
/// \param Q flow [m3/s]
//
void CSubBasin::SetQout(const double &Q)
{
  for (int i=0;i<_nSegments;i++){_aQout[i]=Q;}
  _QoutLast=Q;
}
/////////////////////////////////////////////////////////////////
/// \brief Sets Qlat storage array, usually upon read of state file
/// \param &N [in] size of aQl array (=_nQlatHist)
/// \param *aQl array (size N) history of lateral loss rates from HRUs, water still in subbasin [m3/s]
/// \param QoLast: most recent lateral inflow to reach [m3/s]
//
void CSubBasin::SetQlatHist(const int N, const double *aQl, const double QlLast)
{
  if (N != _nQlatHist) {
    WriteWarning("CSubBasin::SetQlatHist: size of lateral flow history differs between current model and initial conditions file. Array will be truncated",false);
  }
  if(N==0) { return; }
  for (int i=0;i<min(_nQlatHist,N);i++){_aQlatHist[i]=aQl[i];}
  _QlatLast=QlLast;
}

/////////////////////////////////////////////////////////////////
/// \brief Sets Qin storage array, usually upon read of state file
/// \param &N [in] size of aQi array (=_nQinHist)
/// \param *aQi array (size N) recent history of inflow to reach [m3/s]
/// \param QoLast: most recent lateral inflow to reach [m3/s]
//
void CSubBasin::SetQinHist          (const int N, const double *aQi)
{
  if (N != _nQinHist) {
    WriteWarning("CSubBasin::SetQinHist: size of inflow history differs between current model and initial conditions file. Array will be truncated",false);
  }
  if(N==0) { return; }
  for (int i=0;i<min(_nQinHist,N);i++){_aQinHist[i]=aQi[i];}
}

//////////////////////////////////////////////////////////////////
/// \brief Sets inflow into primary channel and updates flow history
/// \remark Called *before* RouteWater routine in Solver.cpp, updating the
/// input flow rates to the valid values so that _aQinHist[0]
/// the *current* inflow rate needed for RouteWater and the history is properly stored
/// \note This water will eventually reach the reach outflow, but possibly
/// after a lag or convolution (due to channel/rivulet storage within
/// the basin)
///
/// \param &Qin [in] New flow [m^3/s]
//
void CSubBasin::UpdateInflow    (const double &Qin)//[m3/s]
{
  for (int n=_nQinHist-1;n>0;n--){
    _aQinHist[n]=_aQinHist[n-1];
  }
  _aQinHist[0]=Qin;
}

//////////////////////////////////////////////////////////////////
/// \brief Sets lateral inflow into primary channel and updates flow history
/// \remark Called *before* RouteWater routine in Solver.cpp, updating the
/// input flow rate to the valid values so that _aQlatHist[0] are the
/// *current* inflow rates needed for RouteWater and the history is properly stored
///
/// \param &Qlat [in] New lateral inflow [m^3/s]
//
void CSubBasin::UpdateLateralInflow    (const double &Qlat)//[m3/s]
{
  for (int n=_nQlatHist-1;n>0;n--){
    _aQlatHist[n]=_aQlatHist[n-1];
  }
  _aQlatHist[0]=Qlat;
}

//////////////////////////////////////////////////////////////////
/// \brief sets unusable flow percentatge
//
void CSubBasin::SetUnusableFlowPercentage(const double &val)
{
  ExitGracefullyIf((val<0) || (val>1.0),
    "CSubBasin:SetUnusableFlowPercentage: invalid value for :UnusableFlowPercentage (must be between zero and one).",BAD_DATA_WARN);
  _unusable_flow_pct=val;
}

//////////////////////////////////////////////////////////////////
/// \brief scales all internal flows by scale factor (for assimilation/nudging)
/// \remark Messes with mass balance something fierce!
/// \param &scale [in] Flow scaling factor
/// \param &scale_last [in] True if Qlast should be scaled (overriding); false for no-data scaling
/// \param &tstep [in] time step [d]
/// \return volume added to system [m3]
//
double CSubBasin::ScaleAllFlows(const double &scale, const bool overriding, const double &tstep, const double &t)
{
  double va=0.0; //volume added [m3]
  double sf=(scale-1.0)/scale;

  if(!overriding)
  {
    for(int n=0;n<_nQlatHist;n++) {
      _aQlatHist[n]*=scale;  upperswap(_aQlatHist[n],0.0);
      va+=_aQlatHist[n]*sf*tstep*SEC_PER_DAY;
    }
    for(int n=0;n<_nQinHist; n++) {
      _aQinHist[n]*=scale; upperswap(_aQinHist[n],0.0);
      va+=_aQinHist[n]*sf*tstep*SEC_PER_DAY;
    }
    for(int i=0;i<_nSegments;i++) {
      _aQout[i]*=scale; upperswap(_aQout[i],0.0);
      va+=0.5*(_aQout[i]+_aQout[i+1])*sf*tstep*SEC_PER_DAY;
    }
  }

  if((overriding) && (_pReservoir==NULL)) {
    for (int i=0;i<_nSegments;i++)
    {
      _aQout[i]*=scale;  upperswap(_aQout[i],0.0);
      va+=0.5*(_aQout[i]+_aQout[i+1])*sf*tstep*SEC_PER_DAY;
    }
  }

  if(_pReservoir!=NULL) {
    va+=_pReservoir->ScaleFlow(scale,overriding,tstep,t);
  }

  //Estivate volume added through scaling
  _channel_storage*=scale; va+=_channel_storage*sf;
  _rivulet_storage*=scale; va+=_rivulet_storage*sf;
  return va;
}
/////////////////////////////////////////////////////////////////
/// \brief Sets (i.e., overrides initial) Downstream ID (use sparingly!)
/// \param down_SBID [in] ID of downstream subbasin
//
void CSubBasin::SetDownstreamID(const long down_SBID)
{
  _downstream_ID=down_SBID;
}
/////////////////////////////////////////////////////////////////
/// \brief Sets Downstream subbasin (performed by model during routing initialization)
/// \param pSB [in] pointer to downstream subbasin
//
void CSubBasin::SetDownstreamBasin(const CSubBasin* pSB)
{
  _pDownSB=pSB;
  if (_pReservoir!=NULL){_pReservoir->SetDownstreamBasin(pSB); }
}
/////////////////////////////////////////////////////////////////
/// \brief increment quantity of total delivered demand from management
/// *called by DemandOptimization routines only*
/// this is zeroed out in UpdateSubBasin() every time step prior to management optimization call
/// \param ii [in] demand index
/// \param Qdel [in] demand delivery [m3/s]
//
void CSubBasin::AddToDeliveredDemand(const int ii, const double &Qdel)
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((ii<0) || (ii>=_nWaterDemands),"CSubBasin::AddToDeliveredDemand: invalid demand index",RUNTIME_ERR);
#endif
  _aQdelivered[ii]=Qdel;
  _Qdelivered+=Qdel;
}
/////////////////////////////////////////////////////////////////
/// \brief update demand magnitudes, called in solver at start of every time step
/// \param &Options [in] Global model options information
/// \param &tt [in] time structure at start of current time step
//
void CSubBasin::UpdateDemands(const optStruct& Options, const time_struct& tt)
{
  for (int ii = 0; ii < _nWaterDemands; ii++) {
    _pWaterDemands[ii]->UpdateDemand(Options,tt);
  }
  if (_pReservoir != NULL) {
    _pReservoir->UpdateDemands(Options,tt);
  }
}

/////////////////////////////////////////////////////////////////
/// \brief increment quantity of total delivered demand from management
/// *called by DemandOptimization routines only*
/// this is zeroed out in UpdateSubBasin() every time step prior to management optimization call
/// \param Qret [in] return flow [m3/s]
//
void CSubBasin::AddToReturnFlow(const double &Qret)
{
  _Qreturn+=Qret;
}
void CSubBasin::RecordReturnFlow(const int ii, const double& Qret)
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((ii<0) || (ii>=_nWaterDemands),"CSubBasin::RecordReturnFlow: invalid demand index",RUNTIME_ERR);
#endif
  _aQreturned[ii]=Qret;
}
/////////////////////////////////////////////////////////////////
/// \brief Sets whether basin is gauged
/// \param isgauged [in] boolean whether subbasin is gauged
//
void CSubBasin::SetGauged(const bool isgauged) {
    _gauged = isgauged;
}
/////////////////////////////////////////////////////////////////
/// \brief Sets basin as disabled
//
void CSubBasin::Disable(){
  _disabled=true;
  for(int k=0; k<_nHydroUnits; k++){ _pHydroUnits[k]->Disable(); }
}
/////////////////////////////////////////////////////////////////
/// \brief Sets basin as enabled
//
void CSubBasin::Enable(){
  _disabled=false;
  for(int k=0; k<_nHydroUnits; k++){ _pHydroUnits[k]->Enable(); }
}
/////////////////////////////////////////////////////////////////
/// \brief Sets assimilation to true
//
void CSubBasin::IncludeInAssimilation() {
  _assimilate=true;
}
//////////////////////////////////////////////////////////////////
/// \brief Calculates subbasin area as a sum of HRU areas
/// \remark Called in CModel::Initialize
/// \note After CModel::Initialize, GetBasinArea should be used
/// \return subbasin Area [km^2]
//
double    CSubBasin::CalculateBasinArea()
{
  _basin_area=0.0;
  for (int k=0;k<_nHydroUnits;k++)
  {
    ExitGracefullyIf(_pHydroUnits[k]->GetArea()<=0.0,
      "CSubBasin::CalculateBasinArea: one or more HRUs has a negative or zero area",BAD_DATA);
    if(_pHydroUnits[k]->IsEnabled())
    {
      _basin_area+=_pHydroUnits[k]->GetArea();
    }
  }
  ExitGracefullyIf(_basin_area<0.0,
                   "CSubBasin::CalculateBasinArea: negative subbasin area!", BAD_DATA);

  return _basin_area;
}

//////////////////////////////////////////////////////////////////
/// \brief Initializes SB attributes
/// \details Calculates total subbasin area; initializes _basin_area, _drainage_area,
/// sets initial flow conditions based on annual avg. rainfall and upstream flows,
/// reserves memory for Qin, Qlat history storage, generates Unit hydrographs
//
/// \note Should be called only after basins are ordered from upstream to
/// downstream, in order from upstream to down
///
/// \param &Qin_avg [in] Average inflow from upstream [m^3/s]
/// \param &Qlat_avg [in] Average lateral inflow [m^3/s]
/// \param &total_drain_area [in] Total drainage area [km^2]
/// \param &Options [in] Global model options information
//
void CSubBasin::Initialize(const double    &Qin_avg,          //[m3/s] from upstream
                           const double    &Qlat_avg,         //[m3/s] from lateral inflow
                           const double    &total_drain_area, //[km2]
                           const optStruct &Options)
{
  string warn;
  if ((_nHydroUnits==0) && (!_is_conduit)){ //allowed if conduit
    warn="CSubBasin::Initialize: subbasin "+to_string(_ID)+" has no constituent HRUs and therefore zero area";
    ExitGracefully(warn.c_str(),BAD_DATA_WARN);
  }
  if ((_nHydroUnits>0) && (_is_conduit)){
    warn="CSubBasin::Initialize: conduit "+to_string(_ID)+" has constituent HRUs. HRUs should not be linked to conduits";
    ExitGracefully(warn.c_str(),BAD_DATA_WARN);
  }
  if((_pChannel==NULL) && (Options.routing!=ROUTE_NONE) && (Options.routing!=ROUTE_EXTERNAL)){
    warn="CSubBasin::Initialize: channel profile for basin "+to_string(_ID)+" may only be 'NONE' if Routing=ROUTE_NONE or ROUTE_EXTERNAL";
    ExitGracefully(warn.c_str(),BAD_DATA);
  }
  if (_pInflowHydro != NULL){_is_headwater=false;}

  _drainage_area=total_drain_area;

  if(!_disabled){
    //Estimate reference flow in non-headwater basins from annual average runoff
    //------------------------------------------------------------------------
    if(_Q_ref==AUTO_COMPUTE)
    {
      if(((Qin_avg + Qlat_avg) <= 0) && (!_is_headwater) && (!_disabled)){//reference flow only matters for non-headwater basins
        string warn="CSubBasin::Initialize: negative or zero average flow specified in initialization (basin "+to_string(_ID)+")";
        ExitGracefully(warn.c_str(),BAD_DATA);
      }
      double mult=_pModel->GetGlobalParams()->GetParams()->reference_flow_mult;
      ResetReferenceFlow(mult*(Qin_avg+Qlat_avg)); //VERY APPROXIMATE - much better to specify!

      //string advice="Reference flow in basin " +to_string(_ID)+" was estimated from :AvgAnnualRunoff to be "+to_string(_Q_ref) +" m3/s. (this will not be used in headwater basins)";
      //WriteAdvisory(advice,false);
    }
    else{
      ResetReferenceFlow(_Q_ref);
    }

    if(_pChannel!=NULL) {
      _pChannel->CheckReferenceFlow(_Q_ref,_slope,_mannings_n,_ID);
    }

    //Estimate reach length from area if not provided
    //------------------------------------------------------------------------
    if (_reach_length==AUTO_COMPUTE)
    {
      _reach_length =pow(_basin_area,0.67)*M_PER_KM;//[m] // \ref from Grand river data, JRC 2010
      string advice="Reach length in basin " +to_string(_ID)+" was estimated from basin area to be "+to_string(_reach_length/M_PER_KM) +" km. (this will not be used in headwater basins)" ;
      WriteAdvisory(advice,false);
    }

    // Estimate reach routing params: t_peak, t_conc, gamma shape/scale
    //------------------------------------------------------------------------
    ///< \ref from Williams (1922), as cited in Handbook of Hydrology, eqn. 9.4.3 \cite williams1922
    if (_t_conc==AUTO_COMPUTE)
    {
      //_t_conc=14.6*_reach_length/M_PER_KM*pow(_basin_area,-0.1)*pow( [[AVERAGE VALLEY SLOPE???]],-0.2)/MIN_PER_DAY;
      _t_conc=0.76*pow(max(_basin_area,0.001),0.38)/HR_PER_DAY;// [d] \ref Austrailian Rainfall and runoff guidelines [McDermott and Pilgrim (1982)]
      _t_conc *= _pModel->GetGlobalParams()->GetParams()->TOC_multiplier;
      //if (Options.catchment_routing!=ROUTE_NONE){
      //  WriteAdvisory("Time of concentration has been estimated as "+to_string(_t_conc)+" days for basin "+to_string(_ID),false);
      //}
    }
    if(Options.catchment_routing==ROUTE_GAMMA_CONVOLUTION ){_t_peak=AUTO_COMPUTE;}
    if(Options.catchment_routing==ROUTE_DUMP              ){_t_peak=AUTO_COMPUTE;}

    if (_t_peak     ==AUTO_COMPUTE){
      _t_peak=0.3*_t_conc;
      _t_peak *= _pModel->GetGlobalParams()->GetParams()->TIME_TO_PEAK_multiplier;
    }
    if (_t_lag      ==AUTO_COMPUTE){_t_lag =0.0;}
    if (_gamma_shape==AUTO_COMPUTE){
      _gamma_shape=3.0;
      _gamma_shape *= _pModel->GetGlobalParams()->GetParams()->GAMMA_SHAPE_multiplier;
    }
    if (_gamma_scale==AUTO_COMPUTE){
      _gamma_scale=max((_gamma_shape-1.0)/_t_peak,0.01); //only really should be used if _gamma_shape>1.0
      _gamma_scale *= _pModel->GetGlobalParams()->GetParams()->GAMMA_SCALE_multiplier;
    }

    if (_reservoir_constant==AUTO_COMPUTE){
      _reservoir_constant=-log(_t_conc/(1+_t_conc));
    }
    ExitGracefullyIf(_t_conc<_t_peak,"CSubBasin::Initialize: time of concentration must be greater than time to peak",BAD_DATA);
    ExitGracefullyIf(_t_peak<=0,     "CSubBasin::Initialize: time to peak must be greater than zero",BAD_DATA);
    ExitGracefullyIf(_t_conc<=0,     "CSubBasin::Initialize: time of concentration must be greater than zero",BAD_DATA);
    ExitGracefullyIf(_gamma_shape<=0,"CSubBasin::Initialize: gamma shape parameter must be greater than zero",BAD_DATA);
    ExitGracefullyIf(_gamma_scale<=0,"CSubBasin::Initialize: gamma scale parameter must be greater than zero",BAD_DATA);

    _diffusivity=GetDiffusivity();

    //generate catchment & routing hydrograph weights
    //reserves memory and populates _aQinHist,_aQlatHist,_aRouteHydro
    //------------------------------------------------------------------------
    GenerateCatchmentHydrograph(Qlat_avg,Options);

    GenerateRoutingHydrograph  (Qin_avg, Options);

    //Initialize flow states aQo, channel_storage, rivulet storage
    //------------------------------------------------------------------------
    InitializeFlowStates(Qin_avg,Qlat_avg,Options);

  } //end if disabled

  //initialize reservoir & time series members
  //------------------------------------------------------------------------
  if (_pReservoir!=NULL){
    _pReservoir->Initialize(Options);
  }
  if (_pInflowHydro != NULL){
    _pInflowHydro->Initialize(Options.julian_start_day,Options.julian_start_year,Options.duration,Options.timestep,false,Options.calendar);
  }
  if (_pInflowHydro2 != NULL){
    _pInflowHydro2->Initialize(Options.julian_start_day,Options.julian_start_year,Options.duration,Options.timestep,false,Options.calendar);
  }
  if(_pEnviroMinFlow != NULL) {
    _pEnviroMinFlow->Initialize(Options.julian_start_day,Options.julian_start_year,Options.duration,Options.timestep,false,Options.calendar);
  }
  //must be initialized AFTER RVM FILE READ



  //QA/QC check of Muskingum parameters, if necessary
  //------------------------------------------------------------------------
  if ((Options.routing==ROUTE_MUSKINGUM) || (Options.routing==ROUTE_MUSKINGUM_CUNGE))
  {
    double dx=_reach_length/(double)(_nSegments);
    double K=GetMuskingumK(dx);
    double X=GetMuskingumX(dx);
    if ((Options.timestep<2*K*X) || (Options.timestep>2*K*(1-X)))
    {
      string warn;
      if (Options.timestep<2*K*X){
        warn ="CSubBasin::Initialize: inappropriate global time step for Muskingum routing in subbasin "+_name+" timestep too small, must decrease size of reach segments";
      }
      else{
        warn="CSubBasin::Initialize: inappropriate global time step for Muskingum routing in subbasin "+_name+": local timestepping will be used";
      }
      WriteWarning(warn,Options.noisy);
      if (Options.noisy){
        cout<<warn<<endl;
        cout<<"   For Muskingum stability, the following condition should hold:"<<endl;
        cout<<"   2KX < dt < 2K(1-X): "<< 2*K*X<<" < "<<Options.timestep<< " < " << 2*K*(1-X)<<endl;
        cout<<"   where K="<<K<<" X="<<X<<" dt="<<Options.timestep<<endl<<endl;
      }
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Initializes SB demand members AFTER RVM FILE READ
/// \param &Options [in] Global model options information
//
void CSubBasin::InitializePostRVM(const optStruct &Options)
{
  if (_nWaterDemands>0){
    _aQdelivered=new double [_nWaterDemands];
    _aQreturned =new double [_nWaterDemands];
    for (int ii=0;ii<_nWaterDemands;ii++){
      _pWaterDemands[ii]->Initialize(Options);
      _aQdelivered  [ii]=0.0;
      _aQreturned   [ii]=0.0;
    }
  }
  if (_pReservoir!=NULL){_pReservoir->InitializePostRVM(Options); }
}
//////////////////////////////////////////////////////////////////
/// \brief Initializes channel flows, channel storage, rivulet storage from estimated mean inflows
/// \details called from base initialization, but also when re-booting .rvc file
///
/// \param &Qin_avg [in] Average inflow from upstream [m^3/s]
/// \param &Qlat_avg [in] Average lateral inflow [m^3/s]
/// \param &Options [in] Global model options information
//
void CSubBasin::InitializeFlowStates(const double& Qin_avg,const double& Qlat_avg,const optStruct &Options)
{
  if(!_disabled) {
    int seg;
    //Set initial conditions for flow history variables (may later be overwritten by .rvc file)
    //------------------------------------------------------------------------
    for(seg=0;seg<_nSegments;seg++) {
      _aQout[seg]=Qin_avg+Qlat_avg*(double)(seg+1)/(double)(_nSegments);
    }
    _QoutLast    =_aQout[_nSegments-1];
    _QlatLast    =Qlat_avg;

    //Calculate Initial Channel Storage from flowrate
    //------------------------------------------------------------------------
    _channel_storage=0.0;
    if((Options.routing!=ROUTE_NONE) && (Options.routing!=ROUTE_EXTERNAL))
    {
      for(seg=0;seg<_nSegments;seg++)
      {
        _channel_storage+=_pChannel->GetArea(_aQout[seg],_slope,_mannings_n)*(_reach_length/_nSegments); //[m3]
      }
    }

    //initialize rivulet storage
    //------------------------------------------------------------------------
    double sum(0.0);
    for(int n=0;n<_nQlatHist;n++)
    {
      sum+=(n+1)*_aUnitHydro[n];
    }
    _rivulet_storage=sum*Qlat_avg*(Options.timestep*SEC_PER_DAY);//[m3];
  } //end if disabled
}
/////////////////////////////////////////////////////////////////
/// \brief Resets reference flow
/// \details Resets celerity and channel top width at reference flow based on reference flow
/// \param &Qreference [in] Reference flow
//
void CSubBasin::ResetReferenceFlow(const double &Qreference)
{
  _Q_ref=Qreference;
  if ((_Q_ref!=AUTO_COMPUTE) && (_pChannel!=NULL))
  {
    if ((_Q_ref <= 0.0) && (!_is_headwater)){
      cout << "_Q_ref=" << _Q_ref << endl;
      ExitGracefully("CSubBasin::ResetReferenceFlow: invalid (negative or zero) reference flow rate in non-headwater basin", BAD_DATA);
    }
    if(_c_ref==AUTO_COMPUTE) {
      _c_ref=_pChannel->GetCelerity(_Q_ref,_slope,_mannings_n);
      ExitGracefullyIf(_c_ref<=0.0,"CSubBasin::ResetReferenceFlow: negative or zero celerity",RUNTIME_ERR);
    }
    _w_ref=_pChannel->GetTopWidth(_Q_ref,_slope,_mannings_n);
  }
  else
  {
    _c_ref=AUTO_COMPUTE;
    _w_ref=AUTO_COMPUTE;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Generates routing (channel) unit hydrograph for subbasin
/// \details generates channel transfer function, reserves memory for and
/// populates _aQinHist,_aRouteHydro
/// \remark Called once for each basin from CSubBasin::initialize
///
/// \param &Qin_avg [in] estimated Average inflow [m^3/s]
/// \param &Options [in & out] Global model options information
//
void CSubBasin::GenerateRoutingHydrograph(const double &Qin_avg,
                                          const optStruct &Options)
{
  int n;
  double tstep=Options.timestep;
  double travel_time; //duration for wave to traverse reach
  int OldnQinHist=_nQinHist;

  travel_time=_reach_length/_c_ref/SEC_PER_DAY; //[day]

  if      (Options.routing==ROUTE_PLUG_FLOW)
  {
    _nQinHist=(int)(ceil(travel_time/tstep))+2;
  }
  else if (Options.routing==ROUTE_DIFFUSIVE_WAVE)
  {
    _nQinHist=(int)(ceil(2*travel_time/tstep))+2;//really a function of reach length and diffusivity
  }
  else if(Options.routing==ROUTE_DIFFUSIVE_VARY)
  {
    _nQinHist=(int)(ceil(2.5*travel_time/tstep))+2;//really a function of reach length and diffusivity
  }
  else if ((Options.routing==ROUTE_HYDROLOGIC) || (Options.routing==ROUTE_TVD))
  {
    _nQinHist=2;
  }
  else
  {
    //_nQinHist=(int)(ceil(2*travel_time/tstep))+2;
    _nQinHist=20;
  }

  //reserve memory, initialize
  _aQinHist   =new double [_nQinHist+1 ];
  for (n=0;n<_nQinHist;n++){_aQinHist[n]=Qin_avg; }

  _aRouteHydro=new double [_nQinHist+1 ];
  for (n=0;n<_nQinHist;n++){_aRouteHydro[n]=0.0;}

  double sum;
  //---------------------------------------------------------------
  if (Options.routing==ROUTE_PLUG_FLOW)
  {
    double ts;
    for (n=0;n<_nQinHist-1;n++)
    {
      ts=(travel_time-n*tstep)/tstep;
      if ((ts>=0.0) && (ts<1.0)){
        _aRouteHydro[n  ]=1-ts;
        _aRouteHydro[n+1]=ts;
      }
    }
  }
  //---------------------------------------------------------------
  else if (Options.routing==ROUTE_DIFFUSIVE_WAVE)
  {
    ExitGracefullyIf(_nSegments>1,
                     "ROUTE_DIFFUSIVE_WAVE only valid for single-segment rivers",BAD_DATA);
    sum=0.0;
    double cc         =_c_ref*SEC_PER_DAY; //[m/day]
	double diffusivity=_diffusivity*SEC_PER_DAY; //[m2/d]

    for (n=0;n<_nQinHist;n++)
    {
      _aRouteHydro[n  ]=max(ADRCumDist(n*tstep,_reach_length,cc,diffusivity)-sum,0.0);
      sum+=_aRouteHydro[n];
    }
    _aRouteHydro[_nQinHist-1]=0.0;//must truncate infinite distrib.

    if (travel_time<tstep){ //very sharp ADR CDF or reach length==0 -override
      _aRouteHydro[0]=1.0-travel_time/tstep;
      _aRouteHydro[1]=travel_time/tstep;
      for (n=2;n<_nQinHist;n++){_aRouteHydro[n]=0.0;}
    }
  }
  //---------------------------------------------------------------
  else if(Options.routing==ROUTE_DIFFUSIVE_VARY)
  {
    //intialized here to dump hydrograph. Updated each timestep within the UpdateSubBasin() routine.
    for(n=0;n<_nQinHist;n++) { _aRouteHydro[n]=0.0; } //should not be used
    _aRouteHydro[0]=1.0;

    double cc=_c_ref*SEC_PER_DAY; //[m/day]
    //double diffusivity=_pChannel->GetDiffusivity(_Q_ref,_slope,_mannings_n)*SEC_PER_DAY;// m2/d
    _c_hist=new double [_nQinHist];
    for(n=0;n<_nQinHist;n++) {_c_hist[n]=cc; }
  }
  //---------------------------------------------------------------
  else
  {
    for (n=0;n<_nQinHist;n++){_aRouteHydro[n  ]=0.0;} //should not be used
    _aRouteHydro[0]=1.0;
  }

  //correct to ensure that sum _aRouteHydro[m]=1.0
  sum=0.0;
  for (n=0;n<_nQinHist;n++){sum+=_aRouteHydro[n];}

  if(sum<0.2) {
    cout<<"diff: "<<_pChannel->GetDiffusivity(_Q_ref,_slope,_mannings_n)<<" cel"<<_c_ref<<endl;
    cout<<"---";for(n=0;n<_nQinHist;n++) { cout<< _aRouteHydro[n]<<" "; }cout<<" SUM: "<<sum<<endl;
    string warning="CSubBasin::GenerateRoutingHydrograph: bad routing hydrograph constructed in subbasin "+to_string(_ID);
    ExitGracefully(warning.c_str(),RUNTIME_ERR); //for very diffusive channels - reach length too long
  }
  for (n=0;n<_nQinHist;n++){_aRouteHydro[n]/=sum;}
}

//////////////////////////////////////////////////////////////////
/// \brief Generates catchment hydrograph for subbasin
/// \details generates catchment hydrographs, reserves memory for and
/// populates aQiHist,_aUnitHydro
/// \remark Called once for each basin from CSubBasin::initialize
///
/// \param &Qlat_avg [in] Estimate of (steady state) inflow [m^3/s]
/// \param &Options [in & out] Global model options information
//
void CSubBasin::GenerateCatchmentHydrograph(const double    &Qlat_avg,
                                            const optStruct &Options)
{
  int n;
  double tstep=Options.timestep;
  double t;

  int OldnQlatHist=_nQlatHist;

  if (Options.catchment_routing==ROUTE_TRI_CONVOLUTION)
  {
    _nQlatHist=(int)(ceil((_t_conc)/tstep))+3;
  }
  else if (Options.catchment_routing==ROUTE_GAMMA_CONVOLUTION)
  {
    _nQlatHist=(int)(ceil(4.5*pow(_gamma_shape,0.6)/_gamma_scale/tstep))+1; //empirical estimate of tail of gamma distribution
  }
  else if (Options.catchment_routing==ROUTE_DELAYED_FIRST_ORDER)
  {
    _nQlatHist=2;
  }
  else if (Options.catchment_routing==ROUTE_DUMP)
  {
    _nQlatHist=3;
  }
  else if (Options.catchment_routing==ROUTE_RESERVOIR_SERIES)
  {
    _nQlatHist=(int)(ceil((4.0/_reservoir_constant)/tstep))*_num_reservoirs+2; //approximate
  }

  //add additional history required for lag
  _nQlatHist+=(int)(ceil(_t_lag/tstep));

  if(_nQlatHist<=0) {
    string warn="GenerateCatchmentHydrograph: negative or zero inflow history. Something is wrong with catchment routing parameters in basin "+to_string(_ID);
    ExitGracefully(warn.c_str(),RUNTIME_ERR);
  }

  //reserve memory, initialize
  _aQlatHist  =new double [_nQlatHist];
  for (n=0;n<_nQlatHist;n++){_aQlatHist [n]=Qlat_avg;}//set to initial (steady-state) conditions

  _aUnitHydro =new double [_nQlatHist];
  for (n=0;n<_nQlatHist;n++){_aUnitHydro[n]=0.0;}

  //generate unit hydrograph
  double sum;
  //---------------------------------------------------------------
  if (Options.catchment_routing==ROUTE_GAMMA_CONVOLUTION)
  {
    sum=0;
    double alpha= _gamma_shape;
    double beta = _gamma_scale;
    for (n=0;n<_nQlatHist;n++)
    {
      t=n*tstep-_t_lag;
      _aUnitHydro[n]=GammaCumDist(t+tstep,alpha, beta)-sum;

      //ExitGracefullyIf(!_finite(_aUnitHydro[n]),
      ExitGracefullyIf(_aUnitHydro[n]>ALMOST_INF,
        "GenerateCatchmentHydrograph: issues with gamma distribution. Time to peak may be too small relative to timestep",RUNTIME_ERR);

      sum+=_aUnitHydro[n];
    }
    _aUnitHydro[_nQlatHist-1]=0.0;//must truncate infinite distribution
  }
  //---------------------------------------------------------------
  else if (Options.catchment_routing==ROUTE_TRI_CONVOLUTION)
  {
    sum=0.0;
    for (n=0;n<_nQlatHist;n++)
    {
      t=n*tstep-_t_lag;
      _aUnitHydro[n]=TriCumDist(t+tstep,_t_conc,_t_peak)-sum;
      sum+=_aUnitHydro[n];
    }
  }
  //---------------------------------------------------------------
  else if (Options.catchment_routing==ROUTE_RESERVOIR_SERIES)
  {
    sum=0.0;
    for (n=0;n<_nQlatHist;n++)
    {
      t=n*tstep-_t_lag;
      _aUnitHydro[n]=NashCumDist(t+tstep,_reservoir_constant,_num_reservoirs)-sum;
      sum+=_aUnitHydro[n];
    }
  }
  //---------------------------------------------------------------
  else if (Options.catchment_routing==ROUTE_DUMP)
  {
    for (n=0;n<_nQlatHist;n++){_aUnitHydro[n]=0.0;}
    _aUnitHydro[0]=1.0;
  }
  //---------------------------------------------------------------

  //correct to ensure that area under _aUnitHydro[n]=1.0
  //---------------------------------------------------------------
  sum=0.0;
  for (n=0;n<_nQlatHist;n++){sum+=_aUnitHydro[n];}
  ExitGracefullyIf(sum==0.0,"CSubBasin::GenerateCatchmentHydrograph: bad unit hydrograph constructed",RUNTIME_ERR);
  if(fabs(sum-1.0)>0.05){ WriteWarning("CSubBasin::GenerateCatchmentHydrograph: unit hydrograph truncated",Options.noisy); }
  for (n=0;n<_nQlatHist;n++){_aUnitHydro[n]/=sum;}
}

//////////////////////////////////////////////////////////////////
/// \brief updates flow algorithms based upon the time
/// \notes can later support quite generic temporal changes to treatment of flow phenom (e.g., changes to unit hydrograph)
/// \param tt [in] current model time
/// \param Options [in] model options structure
//
void  CSubBasin::UpdateSubBasin(const time_struct &tt, const optStruct &Options)
{
  _Qdelivered=0; //reset each time step before demand optmization applied
  for (int ii = 0; ii < _nWaterDemands; ii++) {
    _aQdelivered[ii]=0.0;
    _aQreturned [ii]=0.0;
  }
  _Qreturn=0;
  if (Options.routing==ROUTE_DIFFUSIVE_VARY){UpdateRoutingHydro(Options.timestep); }

  if (_pReservoir != NULL){ _pReservoir->UpdateReservoir(tt,Options); }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets outflow from primary channel and updates flow history
/// \details Also recalculates channel storage (uglier than desired), resets _QlatLast and _QoutLast
/// \remark Called *after* RouteWater routine is called in Solver.cpp, to reset
/// the outflow rates for this basin
///
/// \param *aQo [in] Array of new outflows [m3/s]
/// \param irr_Q [in] actual delivered water demand [m3/s]
/// \param Qdiv [in] diverted flow [m3/s]
/// \param &res_ht [in] new reservoir stage [m]
/// \param res_outflow [in] reservoir outflow [m3/s]
/// \param constraint [in] current constraint applied to reservoir flow [-]
/// \param res_Qstruct [in] pointer to structure for multiple-control reservoir
/// \param &Options [in] Global model options information
/// \param &tt [in] time structure at start of current time step
/// \param initialize [in] Flag to indicate if flows are to only be initialized
//
void CSubBasin::UpdateOutflows   (const double         *aQo,    //[m3/s]
                                  const double         &irr_Q,  //[m3/s]
                                  const double         &Qdiv,   //[m3/s]
                                  const double         &res_ht, //[m]
                                  const double         &res_outflow, //[m3/s]
                                  const res_constraint &constraint,
                                  const double         *res_Qstruct,
                                  const optStruct      &Options,
                                  const time_struct    &tt,
                                  const bool            initialize)
{
  double tstep=Options.timestep;

  //Update flows
  //------------------------------------------------------
  _QoutLast=_aQout[_nSegments-1];
  for (int seg=0;seg<_nSegments;seg++){
    _aQout[seg]=aQo[seg];
  }
  double irrQ_act = min(irr_Q, _aQout[_nSegments-1]);
  _aQout[_nSegments-1]-=irrQ_act;

  double divQ_act = min(Qdiv, _aQout[_nSegments-1]);
  _aQout[_nSegments-1]-=divQ_act;

   //_aQout[num_segments-1] is now the new outflow from the channel

  _QdivLast=_Qdiverted;
  _Qdiverted=divQ_act;
  _QirrLast =_Qirr;
  _Qirr     =irrQ_act;

  if (!Options.management_optimization){
    double total=GetTotalWaterDemand();
    for (int ii = 0; ii < _nWaterDemands; ii++) {
      double demand=_pWaterDemands[ii]->GetDemand();
      _aQdelivered[ii]=_Qirr*(demand/total); //_aQirr may be less than total
      _aQreturned [ii]=0.0; //only valid in man opt mode
    }
  }

  if (_pReservoir!=NULL){
    _pReservoir->UpdateStage(res_ht,res_outflow,constraint,res_Qstruct,Options,tt);
    _pReservoir->UpdateMassBalance(tt,Options.timestep,Options);
  }

  if (initialize){return;}//entering initial conditions

  //Update channel storage
  //------------------------------------------------------
  double dt=tstep*SEC_PER_DAY;
  double dV=0.0;
  double Qlat_new(0.0);
  for (int n=0;n<_nQlatHist;n++){
    Qlat_new+=_aUnitHydro[n]*_aQlatHist[n];
  }
  _QlocLast=_Qlocal;
  _Qlocal=Qlat_new;

  //volume change from linearly varying upstream inflow over this time step
  dV+=0.5*(_aQinHist[0]+_aQinHist[1])*dt;

  //volume change from linearly varying downstream outflow over this time step
  dV-=0.5*(_aQout[_nSegments-1]+_QoutLast)*dt;

  //volume change from lateral inflows over previous time step (corrects for the fact that _aQout is channel flow plus that added from lateral inflow)
  dV+=0.5*(Qlat_new+_QlatLast)*dt;

  //volume change from irrigation
  dV-=0.5*(_Qirr+_QirrLast)*dt;

  //volume change from diversion
  dV-=0.5*(_Qdiverted+_QdivLast)*dt;

  _channel_storage+=dV;//[m3]

  //Update rivulet storage
  //------------------------------------------------------
  dV=0.0;

  //inflow from ground surface (integrated over this time step)
  dV+=_aQlatHist[0]*dt;

  //outflow to channel (integrated over this time step)
  dV-=0.5*(Qlat_new+_QlatLast)*dt;

  _rivulet_storage+=dV;//[m3]

  /* KEEP HERE FOR TESTING
     cout.precision(4);
     cout<<"   UH: ";
     for (int i=0;i<_nQlatHist;i++){cout<<_aUnitHydro[i]<<"  ";}cout<<endl;
     cout<<"   Ql: ";
     for (int i=0;i<_nQlatHist;i++){cout<<_aQlatHist[i]<<"  ";}cout<<endl;
     cout<<"   dV: "<<dV/dt<<" ("<<_aQlatHist[0]<<" - "<< 0.5*(Qlat_new+_QlatLast)<<") [Qlatlast="<<QlatLast<<"]"<<endl;*/

  _QlatLast=Qlat_new;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns Muskingum parameter, K
/// \details Calculates K based upon routing step size, dx, and the channel rating curve properties
/// \param [in] &dx Routing step size [m]
/// \return Muskingum parameter K [d]
//
double  CSubBasin::GetMuskingumK(const double &dx) const
{
  return dx/(_c_ref*SEC_PER_DAY);//[m]/([m/s]*[s/d])
}
//////////////////////////////////////////////////////////////////
/// \brief Returns Muskingum parameter, X
/// \details Calculates X based upon routing step size, dx, and the channel rating curve properties
/// \param [in] &dx Routing step size [m]
/// \return Muskingum parameter X [m]
/// \notes calculated based upon analogy to diffusion equation
//
double  CSubBasin::GetMuskingumX(const double &dx) const
{
  ///< X ranges from 0 to .5 \cite Maidment1993
  ExitGracefullyIf(_pChannel==NULL,"CSubBasin::GetMuskingumX",BAD_DATA);
  double bedslope=_slope;
  if(_slope==AUTO_COMPUTE){bedslope=_pChannel->GetBedslope(); }//overridden by channel

  return max(0.0,0.5*(1.0-_Q_ref/bedslope/_w_ref/_c_ref/dx));//[m3/s]/([m]*[m/s]*[m])
}

//This routine only used in TVD scheme
double CSubBasin::GetReachSegVolume(const double &Qin,  // [m3/s]
                                    const double &Qout, //[m3/s]
                                    const double &dx) const  //[m]
{
  //muskingum-type storage
  double bedslope=_slope;
  if(_slope==AUTO_COMPUTE){bedslope=_pChannel->GetBedslope(); }//overridden by channel

  /*
	// Old version
  double c_avg=_pChannel->GetCelerity(0.5*(Qin+Qout),_slope,_mannings_n);//[m/s]
  double w_avg=_pChannel->GetTopWidth(0.5*(Qin+Qout),_slope,_mannings_n);//[m]

  double K=dx/(c_avg*SEC_PER_DAY); //[d]
  double X =max(0.0,0.5*(1.0-Qin /bedslope/w_avg /c_avg /dx));//[-]
  return (K*X*Qin+K*(1-X)*Qout)*SEC_PER_DAY; //[m3]
  */

	double c_in =_pChannel->GetCelerity(Qin, _slope,_mannings_n);//[m/s]
  double w_in =_pChannel->GetTopWidth(Qin, _slope,_mannings_n);//[m]
  double c_out=_pChannel->GetCelerity(Qout,_slope,_mannings_n);//[m/s]
  double w_out=_pChannel->GetTopWidth(Qout,_slope,_mannings_n);//[m]

  double Kin =dx/(c_in *SEC_PER_DAY); //[d]
  double Kout=dx/(c_out*SEC_PER_DAY); //[d]
  double Xin =max(0.0,0.5*(1.0-Qin /bedslope/w_in /c_in /dx));//[-]
  double Xout=max(0.0,0.5*(1.0-Qout/bedslope/w_out/c_out/dx));//[-]

  //cout<<"K,X: "<<Kin<<" "<<Kout<<" "<<Xin<<" "<<Xout<<" "<<Qin<<" "<<Qout<<" "<<dx<<endl;
  return (Kin*Xin*Qin+Kout*(1-Xout)*Qout)*SEC_PER_DAY; //[m3]

  //level pool-type storage - works great
  //return dx*_pChannel->GetArea(Qout);
  //return dx*0.5*(_pChannel->GetArea(Qin)+_pChannel->GetArea(Qout));  //oscillatory
}


double CSubBasin::thQ(double In_old,double In_new,double Out_old,double Out_new,double th_in,double dx,double tstep) const
{
  //Schwanenberge and Montero, eqn 21
  if(fabs(Out_old-Out_new)<REAL_SMALL){return 0.5;}
  double thQ;
  thQ=(GetReachSegVolume(In_new,Out_new,dx)-GetReachSegVolume(In_old,Out_old,dx))/tstep/SEC_PER_DAY/(Out_old-Out_new);
  thQ+=(-(1.0-th_in)*In_old-(th_in)*In_new+Out_old)/(Out_old-Out_new);
  return thQ;
}

double CSubBasin::TVDTheta(double In_old,double In_new,double Out_old,double Out_new,double th_in, double dx, double tstep) const
{//Schwanenberge and Montero, eqn 22
  double th;
  if     (In_new> Out_old){
    th=thQ(In_old,In_new,Out_old,min(In_new,In_old),th_in,dx,tstep);
    //cout<<"TH!  "<<th<<" "<<min(1.0,max(0.5,th))<<endl;
    //return 1.0;
    return min(1.0,max(0.5,th));
  }
  else if(In_new< Out_old){
    th=thQ(In_old,In_new,Out_old,max(In_new,In_old),th_in,dx,tstep);
    return min(1.0,max(0.5,th));
  }
  return 0.5;
}

//////////////////////////////////////////////////////////////////
/// \brief updates c history, _aRouteHydro for ROUTE_DIFFUSIVE_VARY method
//
void CSubBasin::UpdateRoutingHydro(const double &tstep)
{
  int n;
  double Q=_aQinHist[0]; //Initial guess -we may have to revise this
  //could also be flow-weighted avg flow
  //amount outflowed = aQin[0]*_aRouteHydro[0]+a_Qin[1]*(_aRouteHydro[0]+_aRouteHydro[1])
  //amount still in channel= aQin[0]*(1-A[0])+aQin[1]*(1-A[0]-A[1])...
  double sum=0,sumwt=0,sum2=0;
  for(n=0; n<_nQinHist;n++) {
    sumwt+=_aRouteHydro[n];
    sum+=_aQinHist[n]*(1.0-sumwt);
    sum2+=(1.0-sumwt);
  }
  //Q=sum/sum2; //alternate weighted channel flowrate (ideally would be using current _aRouteHydro, not possible)

  double cc   =_pChannel->GetCelerity   (     Q,_slope,_mannings_n)*SEC_PER_DAY;// [m/day]
  double D_ref=_pChannel->GetDiffusivity(_Q_ref,_slope,_mannings_n)*SEC_PER_DAY;// [m2/d]

  if(_ID==1) {
    cout<<setprecision(3);
    cout<<"Q:  "<<Q<<endl;
    cout<<"c hist: [m/s]";
    for(n=0;n<_nQinHist;n++) {
      cout<<_c_hist[n]/SEC_PER_DAY<<" ";
    }
    cout<<" c_ref: "<<_c_ref<<" "<<_Q_ref<<endl;
  }
  for(n=_nQinHist; n>0; n--) {
    _c_hist[n]=_c_hist[n-1];
  }
  _c_hist[0]=cc;

  sum=0;
  double alpha=D_ref/2;
  for(n=0;n<_nQinHist;n++) {
    //may have to shift _c_hist
    _aRouteHydro[n]=max(TimeVaryingADRCumDist((n)*tstep,_reach_length,_c_hist,_nQinHist,alpha,tstep)-sum,0.0);
    sum+=_aRouteHydro[n];
  }
  for(n=0;n<_nQinHist;n++) {
    _aRouteHydro[n]/=sum;
  }

  //cout<<" SUM: "<<sum<<endl;
  //if reference flow is too fast, fancy stuff goes out the window?
  double travel_time=_reach_length/_c_ref;
  if(travel_time<tstep) { //very sharp ADR CDF or reach length==0 -override
    _aRouteHydro[0]=1.0-travel_time/tstep;
    _aRouteHydro[1]=travel_time/tstep;
    for(int n=2;n<_nQinHist;n++) { _aRouteHydro[n]=0.0; }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Creates aQout_new [m^3/s], an array of point measurements for outflow at downstream end of each river segment
/// \details Array represents measurements at the end of the current timestep. if (catchment_routing==distributed),
/// lateral flow Qlat (e.g., runoff, interflow, or baseflow) is routed as part of channel routing routine otherwise,
/// lateral flow is all routed to most downstream outflow segment . Assumes uniform time step
///
/// \param [out] *aQout_new Array of outflows at downstream end of each segment at end of current timestep [m^3/s]
/// \param &Options [in] Global model options information
/// \param tt [in] - time structure
//
void CSubBasin::RouteWater(double *aQout_new,//[m3/s][size:_nSegments]
                           const optStruct &Options,
                           const time_struct &tt) const
{
  int    seg,n;
  double tstep;       //[d] time step
  double dx;          //[m]
  double Qlat_new;    //[m3/s] flow at end of timestep to reach from catchment storage
  double seg_fraction;//[0..1] % of reach length assoc with segment

  tstep       =Options.timestep;
  dx          =_reach_length/(double)(_nSegments);
  seg_fraction=1.0/(double)(_nSegments);

  //==============================================================
  // route from catchment
  //==============================================================
  Qlat_new=0.0;
  for (n=0;n<_nQlatHist;n++){
    Qlat_new+=_aUnitHydro[n]*_aQlatHist[n];
  }

  //==============================================================
  // route in channel
  //==============================================================
  routing_method route_method;
  route_method=Options.routing;
  if ((_is_headwater) && (route_method!=ROUTE_EXTERNAL)){route_method=ROUTE_NONE;}

  if ((route_method==ROUTE_MUSKINGUM) ||
      (route_method==ROUTE_MUSKINGUM_CUNGE))
  {
    double K,X,c1,c2,c3,c4,denom,cunge;
    K=GetMuskingumK(dx);
    X=GetMuskingumX(dx);

    //may need local time-stepping (2*KX<dt<2*K(1-X))
    //Local time stepping works great if dt is too large, but no fix (other than adding more segments) if dt is too small!
    double dt,Qin,Qin_new;
    dt=min(K,tstep);
    //dt=tstep;

    static double aQoutStored[MAX_RIVER_SEGS];
    for (seg=0;seg<_nSegments;seg++){aQoutStored[seg]=_aQout[seg];}
    //cout<<"check: "<< 2*K*X<<" < "<<dt<< " < " << 2*K*(1-X)<<" K="<<K<<" X="<<X<<" dt="<<dt<<endl;
    for (double t=0;t<tstep;t+=dt)//Local time-stepping
    {
      if (dt>(tstep-t)){dt=tstep-t;}
      cunge=0;
      if (route_method==ROUTE_MUSKINGUM_CUNGE){cunge=1;}

      //Standard Muskingum/Muskingum Cunge
      denom=(2*K*(1.0-X)+dt);
      c1 = ( dt-2*K*(  X)) / denom;
      c2 = ( dt+2*K*(  X)) / denom;
      c3 = (-dt+2*K*(1-X)) / denom;
      c4 = ( dt          ) / denom;

      Qin    =_aQinHist[1]+((t   )/tstep)*(_aQinHist[0]-_aQinHist[1]);
      Qin_new=_aQinHist[1]+((t+dt)/tstep)*(_aQinHist[0]-_aQinHist[1]);

      for (seg=0;seg<_nSegments;seg++)//move downstream
      {
        aQout_new[seg] = c1*Qin_new + c2*Qin + c3*_aQout[seg] + cunge*c4*(Qlat_new*seg_fraction);
        Qin    =_aQout    [seg];
        Qin_new=aQout_new[seg];

        _aQout[seg]=aQout_new[seg];//only matters for dt<tstep
      }

    } //Local time-stepping
    for (seg=0;seg<_nSegments;seg++){_aQout[seg]=aQoutStored[seg];}
  }
  //==============================================================
  else if (route_method==ROUTE_STORAGECOEFF)
  {
    ///< Variable storage routing method \ref (Williams, 1969/ SWAT) \cite williams1969flood
    ///< As interpreted from Williams, 1969/Kim and Lee, HP 2010 \cite williams1969flood
    double ttime;                                       //[d] reach segment travel time
    double storage_coeff;
    double c1,c2,c3;

    double Qin_new=_aQinHist[0];
    double Qin    =_aQinHist[1];
    for (seg=0;seg<_nSegments;seg++)
    {
      ttime=GetMuskingumK(dx)*seg_fraction;// K=representative travel time for the reach segment (else storagecoeff changes and mass flows not preserved)
      storage_coeff = min(1.0 / (ttime/tstep + 0.5),1.0);

      c1=storage_coeff/2;
      c2=storage_coeff/2;
      c3=(1.0-storage_coeff);

      //new outflow proportional to old outflow without correction
      double corr=1.0;
      aQout_new[seg] = c1*Qin + c2*Qin_new + c3*(_aQout[seg]-corr*_Qlocal);
      Qin    =_aQout   [seg];
      Qin_new=aQout_new[seg];
    }
  }
  //==============================================================
  else if (route_method==ROUTE_HYDROLOGIC)
  { ///< basic hydrologic routing using Newton's algorithm
    ///<  ONE SEGMENT ONLY FOR NOW

    ///dV(Q)/dt=(Q_in(n)+Q_in(n+1))/2-(Q_out(n)+Q_out(n+1)(h))/2 -> solve for >Qout_new
    /// rewritten as f(Q)-gamma=0 for Newton's method solution

    const double ROUTE_MAXITER=20;
    const double ROUTE_TOLERANCE=0.0001;//[m3/s]

    double Qout_old =_aQout   [_nSegments-1]-_Qlocal;
    double Qin_new  =_aQinHist[0];
    double Qin_old  =_aQinHist[1];
    double V_old=_pChannel->GetArea(Qout_old,_slope,_mannings_n)*_reach_length;

    int    iter=0;
    double change=0;
    double gamma=V_old+(Qin_old+Qin_new-Qout_old)/2.0*(tstep*SEC_PER_DAY);//[m3]

    double dQ=0.1; //[m3/s]
    double f,dfdQ;
    double Q_guess=Qout_old;
    double relax=1.0;//0.99;

    //double Qg[ROUTE_MAXITER],ff[ROUTE_MAXITER]; //for debugging; retain
    do //Newton's method with discrete approximation of df/dQ
    {
      f   =(_pChannel->GetArea(Q_guess   ,_slope,_mannings_n)*_reach_length+(Q_guess   )/2.0*(tstep*SEC_PER_DAY));
      dfdQ=(_pChannel->GetArea(Q_guess+dQ,_slope,_mannings_n)*_reach_length+(Q_guess+dQ)/2.0*(tstep*SEC_PER_DAY)-f)/dQ;
      change=-(f-gamma)/dfdQ;//[m3/s]
      if (dfdQ==0.0){change=1e-7;}

      //Qg[iter]=Q_guess; ff[iter]=f-gamma;//for debugging; retain

      Q_guess+=relax*change; //relaxation required for some tough cases
      if (Q_guess<0){Q_guess=0; change=0.0;}

      //if (Q_guess < 0){ cout << iter << ":" << Q_guess << " " << f << " " << dfdQ << " " << change << endl; }
      iter++;
      if (iter > 3){ relax = 0.9; }
      if (iter > 10){ relax = 0.7; }
    } while ((iter<ROUTE_MAXITER) && (fabs(change)>ROUTE_TOLERANCE));
    //the above is really fast for most trivial changes in flow, since f is locally linear, 1 to 3 iterations sufficient

    aQout_new[_nSegments-1]=Q_guess;
    if (Q_guess<0){ cout << "Negative flow in basin "<<to_string(this->_ID)<<" Qoutold: "<<Qout_old<<" Qin: "<<Qin_new<<" "<<Qin_old<<endl; }

    if (iter==ROUTE_MAXITER){
      string warn="CSubBasin::RouteWater did not converge after "+to_string(ROUTE_MAXITER)+"  iterations for basin "+to_string(_ID)+
                  " flow: " +to_string(Q_guess)+" stage: "+to_string(_pChannel->GetStageElev(Q_guess,_slope,_mannings_n));
      WriteWarning(warn,false);
      /*for (int i = 0; i < 20; i++){
        string warn = to_string(Qg[i]) + " " + to_string(ff[i]);
        WriteWarning(warn,false);
        }*/
    }
  }
  //==============================================================
  else if(route_method==ROUTE_TVD)
  { //Total variation diminishing method of Schwanenberg and Montero, 2016 -Journal of Hydrology 539 p188-195
    ///< basic hydrologic routing using Newton's algorithm
    ///<  ONE SEGMENT ONLY FOR NOW

    ///dV(Q)/dt=(1-th_I)Q_in(n)+(th_I)Q_in(n+1)-(1-th_Q)Q_out(n)+(th_Q)Q_out(n+1)(h) -> solve for >Qout_new
    /// rewritten as f(Q)-gamma=0 for Newton's method solution

    const double ROUTE_MAXITER=20;
    const double ROUTE_TOLERANCE=0.0001;//[m3/s]

    double Qout_old =_aQout   [_nSegments-1]-_Qlocal;
    double Qin_new  =_aQinHist[0];
		double Qin_old  =_aQinHist[1];

    int    iter=0;
    double change=0;

    double f,df;
    double Q_guess=Qin_new;
    double relax=1.0;//0.99;
    double th_in(0.5); //default for single reach
    double th_out;
    double dx=_reach_length;
    double dQ=0.0001;
    do //Newton's method
    {
      th_out=TVDTheta(Qin_old,Qin_new,Qout_old,Q_guess   ,th_in,dx,tstep);
      cout<<"th_out: "<<th_out<<" "<<Qin_old<<" "<<Qin_new<<" "<<Qout_old<<" "<<Q_guess<<" "<<endl;
      f   =(GetReachSegVolume(Qin_new,Q_guess,dx)-GetReachSegVolume(Qin_old,Qout_old,dx))/(tstep*SEC_PER_DAY);
      f-=(1.0-th_in )*Qin_old +(th_in )*Qin_new;
      f+=(1.0-th_out)*Qout_old+(th_out)*Q_guess;

      th_out=TVDTheta(Qin_old,Qin_new,Qout_old,Q_guess+dQ,th_in,dx,tstep);
      df   =(GetReachSegVolume(Qin_new,Q_guess+dQ,dx)-GetReachSegVolume(Qin_old,Qout_old,dx))/(tstep*SEC_PER_DAY);
      df-=(1.0-th_in )*Qin_old +(th_in )*Qin_new;
      df+=(1.0-th_out)*Qout_old+(th_out)*(Q_guess+dQ);
      df-=f;

      change=-f/(df/dQ);//[m3/s]
      if (df==0.0){change=1e-7;}

      Q_guess+=relax*change; //relaxation required for some tough cases

      /*double minQ=min(min(Qin_old,Qin_new),Qout_old);
      if(Q_guess<minQ){ Q_guess+=0.98*(minQ-Q_guess); }*///cheat -introduces MB errors

      cout<<"Q_guess ["<<iter<<"]: "<< Q_guess<< " "<<f<<endl;
      if (Q_guess<0){Q_guess=REAL_SMALL; change=0.0;}

      iter++;
      if (iter > 3){ relax = 0.9; }
      if (iter > 10){ relax = 0.7; }
    } while ((iter<ROUTE_MAXITER) && (fabs(change)>ROUTE_TOLERANCE));
    //the above is really fast for most trivial changes in flow, since f is locally linear, 1 to 3 iterations sufficient

    aQout_new[_nSegments-1]=Q_guess;
    if (Q_guess<0){ cout << "Negative flow in basin "<<to_string(this->_ID)<<" Qoutold: "<<Qout_old<<" Qin: "<<Qin_new<<" "<<Qin_old<<endl; }

    if (iter==ROUTE_MAXITER){
      string warn="CSubBasin::RouteWater:TVD did not converge after "+to_string(ROUTE_MAXITER)+"  iterations for basin "
                   +to_string(_ID)+ " flow: " +to_string(Q_guess)+" stage: "+to_string(_pChannel->GetStageElev(Q_guess,_slope,_mannings_n));
      WriteWarning(warn,false);
      /*for (int i = 0; i < 20; i++){
        string warn = to_string(Qg[i]) + " " + to_string(ff[i]);
        WriteWarning(warn,false);
      }*/
    }

  }
  //==============================================================
  else if ((route_method==ROUTE_PLUG_FLOW)
           || (route_method==ROUTE_DIFFUSIVE_WAVE))
  {
    //Simple convolution - segmentation unused
    aQout_new[_nSegments-1]=0.0;
    for (n=0;n<_nQinHist;n++){
      aQout_new[_nSegments-1]+=_aRouteHydro[n]*_aQinHist[n];
    }
  }
  //==============================================================
  else if (route_method==ROUTE_DIFFUSIVE_VARY)
  {
    aQout_new[_nSegments-1]=0.0;
    for(n=0;n<_nQinHist;n++) {
      aQout_new[_nSegments-1]+=_aRouteHydro[n]*_aQinHist[n];
    }
  }
  //==============================================================
  else if (route_method==ROUTE_EXTERNAL)
  {//In channel routing skipped, overridden flow from .rvl file used
    for (seg=0;seg<_nSegments;seg++){
      aQout_new[seg]=0.0;
    }
    aQout_new[_nSegments-1]=_aQout[_nSegments-1]; //flow is same as overriden flow from start of timestep
  }
  //==============================================================
  else if (route_method==ROUTE_NONE)
  {//In channel routing skipped, inflow directly routed out - segmentation unused
    for (seg=0;seg<_nSegments;seg++){
      aQout_new[seg]=0.0;
    }
    aQout_new[_nSegments-1]=_aQinHist[0];
  }
  //==============================================================
  else{
    ExitGracefully("Unrecognized routing method",STUB);
  }

  //all fluxes from catchment are routed directly to basin outlet
  aQout_new[_nSegments-1]+=Qlat_new;


  /*if (aQout_new[_nSegments-1]<-REAL_SMALL){
    cout<<"Lateral inflow:" <<Qlat_new<<endl;
    cout<<"Outflow: "<<aQout_new[_nSegments-1]<<endl;
    ExitGracefully("RouteWater: negative flow!",RUNTIME_ERR);
    }*/

  return;
}

/////////////////////////////////////////////////////////////////
/// \brief Write subbasin state variable data to solution file
/// \param &OUT [out] Output file stream to which data will be written
//
void CSubBasin::WriteToSolutionFile (ofstream &RVC) const
{
  RVC<<_name<<endl;
  RVC<<"    :ChannelStorage, "<<_channel_storage<<endl;
  RVC<<"    :RivuletStorage, "<<_rivulet_storage<<endl;
  RVC<<"    :Qout,"<<_nSegments  ;for (int i=0;i<_nSegments;i++){RVC<<","<<_aQout    [i]; if (i%200==199){RVC<<endl<<"    ";}}RVC<<","<< _QoutLast << endl;
  RVC<<"    :Qlat,"<<_nQlatHist  ;for (int i=0;i<_nQlatHist;i++){RVC<<","<<_aQlatHist[i]; if (i%200==199){RVC<<endl<<"    ";}}RVC<<","<< _QlatLast << endl;
  RVC<<"    :Qin ,"<<_nQinHist   ;for (int i=0;i<_nQinHist; i++){RVC<<","<<_aQinHist [i]; if (i%200==199){RVC<<endl<<"    ";}}RVC<<endl;
  if (_pReservoir!=NULL){
    _pReservoir->WriteToSolutionFile(RVC);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief clears all time series data for re-read of .rvt file
/// \remark Called only in ensemble mode
///
/// \param &Options [in] Global model options information
//
void CSubBasin::ClearTimeSeriesData(const optStruct& Options)
{
  delete _pInflowHydro;  _pInflowHydro=NULL;
  delete _pInflowHydro2; _pInflowHydro2=NULL;
  //delete [] _aQdelivered;
  //\todo[fix]: THERE WILL BE PROBLEMS IF DEMAND TIME SERIES IS IN ENSEMBLE SIMULATION
  delete _pEnviroMinFlow;_pEnviroMinFlow=NULL;
  _pReservoir->ClearTimeSeriesData(Options);
}
