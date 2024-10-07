/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2024 the Raven Development Team
  ----------------------------------------------------------------*/
#ifndef SUBBASIN_H
#define SUBBASIN_H

#include "RavenInclude.h"
#include "HydroUnits.h"
#include "SoilAndLandClasses.h"
#include "ChannelXSect.h"
#include "TimeSeries.h"
#include "Reservoir.h"
class CReservoir;
class CDemand;
class CChannelXSect;  // defined in ChannelXSect.h
enum res_constraint;

///////////////////////////////////////////////////////////////////
/// \brief flow diversion data strucure
/// \details relates diversion quantity to flow rate in subbasin
struct diversion
{
  int     julian_start;   //Julian start date of flow diversion (e.g., 90 for Apr 1)
  int     julian_end;     //Julian end date of flow diversion (e.g., 242 for Aug 31)
  int     target_p;       //INDEX (not subbasin ID) of target basin (or -1 if diverted from system)

  double  min_flow;       //minimum flow above which flow is diverted [m3/s] (Q>0)
  double  percentage;     //percentage of flow above minimum which is diverted

  int     nPoints;        //number of points in flow-diversion lookup table
  double *aQsource;       //array of discharges [m3/s] in flow diversion lookup table
  double *aQdivert;       //array of diversion flow rates [m3/s] correspionding to discharges in flow diversion lookup table

  diversion()  {aQsource=NULL; aQdivert=NULL;julian_start=julian_end=0; target_p=DOESNT_EXIST;min_flow=0;percentage=1.0;nPoints=0;}
  ~diversion() {delete [] aQsource; delete [] aQdivert;}
};

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction class for contiguous watershed section with a primary channel, contains a collection of HRUs
/// \details Used primarily to route water
class CSubBasin
{
private:/*------------------------------------------------------*/

  long                     _ID;   ///< unique ID of subbasin (must be positive)
  string                 _name;   ///< name
  bool               _disabled;   ///< true if disabled
  int                _global_p;   ///< p index in _pModel subbasin array

  const CModelABC     *_pModel;   ///< Pointer to model

  //basin properties
  double           _basin_area;   ///< contributing surface area for subbasin [km2]
  double        _drainage_area;   ///< total upstream drainage area [km2] (includes subbasin area)
  long          _downstream_ID;   ///< ID of downstream subbasin; if <0, then this outflows outside the model domain
  bool                 _gauged;   ///< if true, hydrographs are generated for downstream flows
  bool           _is_headwater;   ///< true if no subbasins drain into this one and _pInflowHydro==NULL
  bool             _is_conduit;   ///< true if basin is 'conduit' basin, and shouldn't have HRUs
  bool             _assimilate;   ///< true if observed data in this basin should be assimilated.
  const CSubBasin    *_pDownSB;   ///< pointer to downstream subbasin instance

  //catchment routing properties
  double               _t_conc;   ///< basin time of concentration [d]
  double               _t_peak;   ///< basin time to peak [d] (<=_t_conc)
  double                _t_lag;   ///< basin time lag [d]
  double   _reservoir_constant;   ///< linear basin/catchment routing constant [1/d]
  int          _num_reservoirs;   ///< number of linear reservoirs used for in-catchment routing
  double          _gamma_shape;   ///< shape parameter of gamma unit hydrograph
  double          _gamma_scale;   ///< scale parameter of gamma unit hydrograph
  int          _reach_HRUindex;   ///< HRU *index* k (not ID) associated with reach. Used for reach-specific forcings.
  double      _sens_exch_coeff;   ///< in-catchment sensible exchange coefficient [1/d]
  double        _GW_exch_coeff;   ///< in-catchment groundwater temperature exchange coefficient [1/d]

  double       _hyporheic_flux;   ///< gross exchange flux with groundwater [m/d]
  double        _convect_coeff;   ///< convection coefficient [MJ/m2/d/K]
  double     _bed_conductivity;   ///< bed thermal conductivity [MJ/m/d/K]
  double        _bed_thickness;   ///< mean thickness of stream bed [m]

  //River/stream  channel data:
  const CChannelXSect*_pChannel;  ///< Main channel

  double         _reach_length;   ///< length of subbasin reach [m]
  double                _Q_ref;   ///< reference flow rate [m3/s]
  double                _c_ref;   ///< celerity at reference flow rate (if ==AUTO_COMPUTE, or user-specified otherwiae)[m/s]
  double                _w_ref;   ///< channel top width at reference flow rate [m]
  double           _mannings_n;   ///< manning's n for channel (or uses channel one, if =AUTO_COMPUTE)
  double                _slope;   ///< channel slope (rise over run) (or uses channel one, if =AUTO_COMPUTE)
  double          _diffusivity;   ///< channel diffusivity (specified or calculated, if =AUTO_COMPUTE) [m2/s]

  //Other params
  double            _rain_corr;   ///< correction factor for rainfall [-]
  double            _snow_corr;   ///< correction factor for snowfall [-]
  double            _temperature_corr;   ///< correction factor (additive) for temperature [-]

  int               _nSegments;   ///< Number of river segments used in routing(>=1)

  //Reservoir
  CReservoir      *_pReservoir;   ///< Reservoir object (or NULL, if no reservoir)
  bool           _res_disabled;   ///< true if lake or reservoir should be disabled/never created
  double        _reach_length2;   ///< corrected reach length [m] if reservoir/lake is disabled

  //state variables:
  double               *_aQout;   ///< downstream river (out)flow [m3/s] at start of time step at end of each channel segment [size=_nSegments]
  double           *_aQlatHist;   ///< history of lateral runoff into surface water [m3/s][size:_nQlatHist] - uniform (time-averaged) over timesteps
  //                              ///  if Ql=Ql(t), aQlatHist[0]=Qlat(t to t+dt), aQlatHist[1]=Qlat(t-dt to t)...
  int               _nQlatHist;   ///< size of _aQlatHist array
  double      _channel_storage;   ///< water storage in channel [m3]
  double      _rivulet_storage;   ///< water storage in rivulets [m3]
  double             _QoutLast;   ///< Qout from downstream channel segment [m3/s] at start of previous timestep- needed for reporting integrated outflow
  double             _QlatLast;   ///< Qlat (after convolution) at start of previous timestep [m3/s]
  //double           _aQchanLast;   ///< channel outflow
  double               _Qlocal;   ///< local contribution to current subbasin outflow [m3/s] (just from in-catchment routing)
  double             _QlocLast;   ///< last local contribution [m3/s]

  double                 _Qirr;   ///< Qirr (delivered water demand) at end of timestep [m3/s] (for MB accounting)
  double             _QirrLast;   ///< Qirr (delivered water demand) at start of timestep [m3/s] (for MB accounting)
  double             _QdivLast;   ///< diverted flow at start of timestep [m3/s] (for MB accounting)
  double            _Qdiverted;   ///< diverted flow at end of timestep [m3/s] (for MB accounting)
  double           _Qdelivered;   ///< delivered flow at end of time step if management optimization is used [m3/s]
  double         *_aQdelivered;   ///< delivered flow at end of time step for each water demand [m3/s] [size: _nIrrigDemands]
  double          *_aQreturned;   ///< return flow associated with deliveries (may be redirected to other basin) [m3/s]
  double              _Qreturn;   ///< return flow at end of time step if management optimization is used [m3/s] (may be sourced from other basin)

  //Hydrograph Memory
  double            *_aQinHist;   ///< history of inflow from upstream into primary channel [m3/s][size:nQinHist] (aQinHist[n] = Qin(t-ndt))
  //                              ///  _aQinHist[0]=Qin(t), _aQinHist[1]=Qin(t-dt), _aQinHist[2]=Qin(t-2dt)...
  int                _nQinHist;   ///< size of _aQinHist array
  double              *_c_hist;   ///< reach celerity history [size: _nQinHist] (used for ROUTE_DIFFUSIVE_VARY only)

  //characteristic weighted hydrographs
  double          *_aUnitHydro;   ///< [size:_nQlatHist] catchment unit hydrograph (time step-dependent). area under = 1.0.
  double         *_aRouteHydro;   ///< [size:_nQinHist ] routing unit hydrograph. area under = 1.0.

  //HRUs
  int             _nHydroUnits;   ///< constituent HRUs with different hydrological characteristics
  CHydroUnit    **_pHydroUnits;   ///< [size:nHydroUnits] Array of pointers to constituent HRUs

  //Treatment Plant/Irrigation/Other incoming hydrograph
  CTimeSeries   *_pInflowHydro;   ///< pointer to time series of inflows; NULL if no specified input - Inflow at upstream entrance of basin
  CTimeSeries  *_pInflowHydro2;   ///< pointer to time series of inflows/extractions ; at downstream end of basin reach
  CTimeSeries *_pEnviroMinFlow;   ///< pointer to time series of environmental minimum flow targets that force reduced irrigation demand (=0 by default)

  int           _nWaterDemands;   ///< number of water demand objects
  CDemand     **_pWaterDemands;   ///< pointer to array of demand objects; demand applied at downstream end of basin reach [size: _nIrrigDemand]

  int             _nDiversions;   ///< number of flow diversions from basin
  diversion     **_pDiversions;   ///< array of pointers to flow diversion structures

  double    _unusable_flow_pct;   ///< for irrigation, fraction of flow above enviro min. flow that is not allowed to be used [0..1] (0 default)

  //Methods implemented in SubBasin.cpp
  double               GetMuskingumK(const double &dx) const;
  double               GetMuskingumX(const double &dx) const;
  void     GenerateRoutingHydrograph(const double &Qin_avg, const optStruct &Options);
  void   GenerateCatchmentHydrograph(const double &Qlat_avg,const optStruct &Options);
  double           GetReachSegVolume(const double &Qin, const double &Qout, const double &dx) const;

  double                         thQ(double In_old,double In_new,double Out_old,double Out_new,double th_in,double dx,double tstep) const;
  double                    TVDTheta(double In_old,double In_new,double Out_old,double Out_new,double th_in,double dx,double tstep) const;

  void            UpdateRoutingHydro(const double &tstep);
public:/*-------------------------------------------------------*/
  //Constructors:
  CSubBasin(const long           ID,
            const string         Name,
            const CModelABC     *pMod,
            const long           downstream_ID, //index of downstream SB, if <0, downstream outlet
            const CChannelXSect *pChan,         //Channel
            const double         reach_len,     //reach length [m]
            const double         Qreference,    //reference flow [m3/s]
            const bool           gauged,        //if true, hydrographs are generated
            const bool           is_conduit);   //to tag subbasin as conduit
  ~CSubBasin();

  //Accessor functions
  long                 GetID                () const;
  int                  GetGlobalIndex       () const; //[p]
  string               GetName              () const;
  double               GetBasinArea         () const;
  double               GetDrainageArea      () const;
  double               GetAvgStateVar       (const int i) const;
  double               GetAvgConcentration  (const int i) const;
  double               GetAvgForcing        (const forcing_type &ftype) const;
  double               GetAvgCumulFlux      (const int i, const bool to) const;
  double               GetAvgCumulFluxBet   (const int iFrom, const int iTo) const;
  double               GetReferenceFlow     () const;
  double               GetReferenceCelerity () const;
  double               GetReferenceXSectArea() const;
  double               GetCelerity          () const;
  double               GetDiffusivity       () const;
  long                 GetDownstreamID      () const;
  int                  GetNumHRUs           () const;
  const CHydroUnit    *GetHRU               (const int k) const;
  bool                 IsGauged             () const;
  double               GetReachLength       () const;
  int                  GetNumSegments       () const;
  bool                 IsEnabled            () const;
  bool                 IsHeadwater          () const;
  double               GetRainCorrection    () const;
  double               GetSnowCorrection    () const;
  double               GetTemperatureCorrection () const;
  int                  GetReachHRUIndex     () const;
  double               GetRiverDepth        () const;
  double               GetXSectArea         () const;
  double               GetWaterLevel        () const;
  double               GetHyporheicFlux     () const;
  double               GetConvectionCoeff   () const;
  double               GetRiverbedConductivity() const;
  double               GetRiverbedThickness () const;
  double               GetBedslope          () const;
  double               GetWettedPerimeter   () const;
  double               GetTopWidth          () const;
  bool                 UseInFlowAssimilation() const;
  int                  GetNumWaterDemands   () const;
  double               GetUnusableFlowPercentage() const;

  const double   *GetUnitHydrograph        () const;
  const double   *GetRoutingHydrograph     () const;
  const double   *GetInflowHistory         () const;
  const double   *GetLatHistory            () const;
  const double   *GetOutflowArray          () const;
  int             GetLatHistorySize        () const;
  int             GetInflowHistorySize     () const;
  int             GetOutflowArraySize      () const;
  int             GetNumDiversions         () const;

  double          GetOutflowRate           () const;                   //[m3/s] from final reach segment OR reservoir, point in time
  double          GetChannelOutflowRate    () const;                   //[m3/s] from final reach segment (NOT reservoir), point in time, BEFORE diversions included
  double          GetLastChannelOutflowRate() const;                   //[m3/s] from final reach segment (NOT reservoir), point in time, BEFORE diversions included
  double          GetLastOutflowRate       () const;                   //[m3/s] from final segment, previous timestep
  double          GetLocalOutflowRate      () const;                   //[m3/s] local contribution to outflow
  double          GetIntegratedOutflow     (const double &tstep) const;//[m3] from final segment integrated over timestep
  double          GetIntegratedSpecInflow  (const double &t,
                                            const double &tstep) const;//[m3] from specified inflows integrated over timestep
  double          GetIntegratedLocalOutflow(const double &tstep) const;//[m3] from in-catchment local contributions, integrated over timestep
  double          GetReservoirInflow       () const;                   //[m3/s] from final segment upstream of reservoir, point in time
  double          GetReservoirLosses       (const double &tstep) const;//[m3] from reservoir integrated over timestep
  double      GetIntegratedReservoirInflow (const double &tstep) const;//[m3] from final segment upstream of reservoir integrated over timestep
  double          GetIrrigationLosses      (const double &tstep) const;//[m3] from actual irrigation (not just demand)
  double          GetDiversionLosses       (const double &tstep) const;//[m3] from flow diversions

  double          GetRivuletStorage        () const;                   //[m3] volume en route to outflow
  double          GetChannelStorage        () const;                   //[m3] volume in channel
  double          GetReservoirStorage      () const;                   //[m3] volume in reservoir

  double          GetSpecifiedInflow       (const double &t) const;    //[m3/s] to upstream end of channel at point in time
  double          GetDownstreamInflow      (const double &t) const;    //[m3/s] to downstream end of channel at point in time
  double          GetTotalWaterDemand      () const;                   //[m3/s] total from downstream end of channel at point in time
  double          GetWaterDemand           (const int ii) const;       //[m3/s] iith demand from downstream end of channel at point in time
  double          GetTotalReturnFlow       () const;                   //[m3/s] total instantaneous return flow TO this basin
  double          GetReturnFlow            (const int ii) const;       //[m3/s] return flow associated with iith demand (may be redirected OUT of basin)
  double          GetDownstreamIrrDemand   (const double &t) const;    //[m3/s] cumulative downstream irrigation demand, including from this subbasin
  double          GetDemandDelivery        () const;                   //[m3/s] instantaneous total delivery rate Qirr
  double          GetDemandDelivery        (const int ii) const;       //[m3/s] instantaneous delivery rate to demand ii
  double          GetEnviroMinFlow         (const double &t) const;    //[m3/s] environmental minimum flow target from downstream outlet
  bool            HasEnviroMinFlow         () const;                   // true if basin has enviro min flow time series
  bool            HasIrrigationDemand      () const;                   // true if basin has specified irrigation demand
  int             GetDiversionTargetIndex  (const int i) const;        // returns subbasin index p of diversion i

  CReservoir     *GetReservoir             () const;
  CDemand        *GetWaterDemandObj        (const int ii) const;       //returns pointer to ith water demand

  double          GetBasinProperties       (const string label) const;

  //Manipulator functions
  //called during model construction/assembly:
  void            AddHRU                   (CHydroUnit *pHRU);
  void            AddReservoir             (CReservoir *pReservoir);
  void            SetGlobalIndex           (const int p);
  bool            SetBasinProperties       (const string label,const double &value);
  void            SetAsNonHeadwater        ();
  double          CalculateBasinArea       ();
  void            Initialize               (const double    &Qin_avg,          //[m3/s]
                                            const double    &Qlat_avg,         //[m3/s]
                                            const double    &total_drain_area, //[km2]
                                            const optStruct &Options);
  void            InitializePostRVM        (const optStruct &Options);
  void            InitializeFlowStates     (const double& Qin_avg,const double& Qlat_avg,const optStruct &Options);
  void            AddInflowHydrograph      (CTimeSeries *pInflow);
  void            AddDownstreamInflow      (CTimeSeries *pInflow);
  void            AddEnviroMinFlow         (CTimeSeries *pMinFlow);
  void            AddWaterDemand           (CDemand     *pDemand);
  void            AddFlowDiversion         (const int jul_start, const int jul_end, const int target_p, const double min_flow, const double pct);
  void            AddFlowDiversion         (const int jul_start, const int jul_end, const int target_p, const double *aQ1, const double *aQ2, const int NQ);
  void            ClearTimeSeriesData      (const optStruct& Options);

  // reservoir manipulators
  void            ResetReferenceFlow       (const double &Qreference);

  void            SetChannelStorage        (const double &V);
  void            SetRivuletStorage        (const double &V);
  void            SetQout                  (const double &Q);
  void            SetQoutArray             (const int N, const double *aQo, const double QoLast);
  void            SetQlatHist              (const int N, const double *aQl, const double QlLast);
  void            SetQinHist               (const int N, const double *aQi);
  void            SetDownstreamID          (const long down_SBID);
  void            SetDownstreamBasin       (const CSubBasin *pSB);
  void            SetGauged                (const bool isgauged);
  void            Disable                  ();
  void            Enable                   ();
  double          ScaleAllFlows            (const double &scale_factor, const bool scale_last, const double &tstep, const double &t);
  void            SetUnusableFlowPercentage(const double &val);
  void            IncludeInAssimilation    ();

  //called during model operation:
  void            UpdateDemands            (const optStruct &Options,const time_struct &tt);
  void            AddToDeliveredDemand     (const int ii, const double &Qdel);
  void            RecordReturnFlow         (const int ii, const double& Qret);
  void            AddToReturnFlow          (const double &Qret);
  void            UpdateInflow             (const double &Qin );//[m3/s]
  void            UpdateLateralInflow      (const double &Qlat);//[m3/s]
  void            UpdateSubBasin           (const time_struct &tt, const optStruct &Options);
  void            UpdateOutflows           (const double *Qout_new,
                                            const double &Qirr,
                                            const double &Qdiv,
                                            const double &res_ht,
                                            const double &res_outflow,
                                            const res_constraint &constraint,
                                            const double    *res_qstruct,
                                            const optStruct &Options,
                                            const time_struct &tt,
                                            const bool    initialize);//[m3/s]

  double          ApplyIrrigationDemand    (const double &t,const double &Q,const bool optimized) const; //[m3/s]
  double          GetDiversionFlow         (const int i, const double &Q, const optStruct &Options, const time_struct &tt, int &pDivert) const;

  void            RouteWater               (      double      *Qout_new,
                                            const optStruct   &Options,
                                            const time_struct &tt) const;

  void            WriteToSolutionFile      (ofstream &OUT) const;
};

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction class for collection of Subbasins
//
class CSubbasinGroup
{
private:/*------------------------------------------------------*/

  string          _name;       ///< Subbasin group name
  int             _nSubbasins; ///< Number of HRUs present in group
  CSubBasin     **_pSubbasins; ///< Array of pointers to member HRUs (size = nHRUs)
  int             _global_pp;  ///< index of group in master Subbasin Group array (in CModel)
  bool            _disabled;   ///< true if all Subbasins in group are disabled

public:/*-------------------------------------------------------*/
  //Constructors:
  CSubbasinGroup(string tag,int global_ind);
  ~CSubbasinGroup();

  void AddSubbasin(CSubBasin *pSubBasin);
  void DisableGroup();
  void Initialize();

  string            GetName            () const;
  int               GetNumSubbasins    () const;
  int               GetGlobalIndex     () const;
  CSubBasin        *GetSubBasin        (const int p_local) const;
  double            GetAvgStateVar     (const int i) const;
  double            GetAvgConcentration(const int i) const;
  double            GetAvgForcing      (const forcing_type &ftype) const;
  double            GetAvgCumulFlux    (const int i,const bool to) const;
  double            GetAvgCumulFluxBet (const int iFrom,const int iTo) const;
  bool              IsInGroup          (const long SBID) const;
  bool              IsDisabled         () const;
};

#endif
