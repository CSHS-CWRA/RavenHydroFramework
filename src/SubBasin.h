/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2018 the Raven Development Team
  ----------------------------------------------------------------*/
#ifndef SUBBASIN_H
#define SUBBASIN_H

#include "RavenInclude.h"
#include "HydroUnits.h"
#include "SoilAndLandClasses.h"
#include "ChannelXSect.h"
#include "TimeSeries.h"
#include "Reservoir.h"

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction class for contiguous watershed section with a primary channel, contains a collection of HRUs
/// \details Used primarily to route water
class CSubBasin
{
private:/*------------------------------------------------------*/

  long                     _ID;   ///< unique ID of subbasin (must be positive)
  string                 _name;   ///< name
  bool               _disabled;   ///< true if disabled

  const CModelABC     *_pModel;   ///< Pointer to model

  //basin properties
  double           _basin_area;   ///< contributing surface area for subbasin [km2]
  double        _drainage_area;   ///< total upstream drainage area [km2] (includes subbasin area)
  double         _avg_ann_flow;   ///< average annual flowrate [m3/s] (averaged along reach)
  long          _downstream_ID;   ///< ID of downstream subbasin; if <0, then this outflows outside the model domain
  bool                 _gauged;   ///< if true, hydrographs are generated for downstream flows
  bool           _is_headwater;   ///< true if no subbasins drain into this one and _pInflowHydro==NULL

  //catchment routing properties
  double               _t_conc;   ///< basin time of concentration [d]
  double               _t_peak;   ///< basin time to peak [d] (<=_t_conc)
  double                _t_lag;   ///< basin time lag [d]
  double   _reservoir_constant;   ///< linear basin/catchment routing constant [1/d]
  int          _num_reservoirs;   ///< number of linear reservoirs used for in-catchment routing
  double          _gamma_shape;   ///< shape parameter of gamma unit hydrograph 
  double          _gamma_scale;   ///< scale parameter of gamma unit hydrograph

  //River/stream  channel data:
  const CChannelXSect*_pChannel;  ///< Main channel

  double         _reach_length;   ///< length of subbasin reach [m]
  double                _Q_ref;   ///< reference flow rate [m3/s]
  double                _c_ref;   ///< celerity at reference flow rate [m/s]
  double                _w_ref;   ///< channel top width at reference flow rate [m]
  double           _mannings_n;   ///< manning's n for channel (or uses channel one, if =AUTO_COMPUTE)
  double                _slope;   ///< channel slope (rise over run) (or uses channel one, if =AUTO_COMPUTE)

  //Other params
  double            _rain_corr;   ///< correction factor for rainfall [-]
  double            _snow_corr;   ///< correction factor for snowfall [-]

  int               _nSegments;   ///< Number of river segments used in routing(>=1)

  //Reservoir
  CReservoir      *_pReservoir;   ///< Reservoir object (or NULL, if no reservoir)

  //state variables:
  double               *_aQout;   ///< downstream river (out)flow [m3/s] at start of time step at end of each channel segment [size=_nSegments]
  double           *_aQlatHist;   ///< history of lateral runoff into surface water [m3/s][size:_nQlatHist] - uniform (time-averaged) over timesteps
  //                              ///  if Ql=Ql(t), aQlatHist[0]=Qlat(t to t+dt), aQlatHist[1]=Qlat(t-dt to t)...
  int               _nQlatHist;   ///< size of _aQlatHist array
  double      _channel_storage;   ///< water storage in channel [m3]
  double      _rivulet_storage;   ///< water storage in rivulets [m3]
  double             _QoutLast;   ///< Qout from downstream channel segment [m3/s] at start of previous timestep- needed for reporting integrated outflow
  double             _QlatLast;   ///< Qlat (after convolution) at start of previous timestep [m3/s]

  //Hydrograph Memory
  double            *_aQinHist;   ///< history of inflow from upstream into primary channel [m3/s][size:nQinHist] (aQinHist[n] = Qin(t-ndt))
  //                              ///  _aQinHist[0]=Qin(t), _aQinHist[1]=Qin(t-dt), _aQinHist[2]=Qin(t-2dt)...
  int                _nQinHist;   ///< size of _aQinHist array

  //characteristic weighted hydrographs
  double          *_aUnitHydro;   ///< [size:_nQlatHist] catchment unit hydrograph (time step-dependent). area under = 1.0.
  double         *_aRouteHydro;   ///< [size:_nQinHist ] routing unit hydrograph. area under = 1.0.

  //HRUs
  int             _nHydroUnits;   ///< constituent HRUs with different hydrological characteristics
  CHydroUnit    **_pHydroUnits;   ///< [size:nHydroUnits] Array of pointers to constituent HRUs

  //Treatment Plant/Other incoming hydrograph
  CTimeSeries   *_pInflowHydro;   ///< pointer to time series of inflows; NULL if no specified input - Inflow at upstream entrance of basin
  CTimeSeries   *_pInflowHydro2;  ///< pointer to time series of inflows/extractions ; at downstream end of basin reach

  //Methods implemented in SubBasin.cpp
  double               GetMuskingumK(const double &dx) const;
  double               GetMuskingumX(const double &dx) const;
  void     GenerateRoutingHydrograph(const double &Qin_avg, const optStruct &Options);
  void   GenerateCatchmentHydrograph(const double &Qlat_avg,const optStruct &Options);
  double           GetReachSegVolume(const double &Qin, const double &Qout, const double &dx) const;

  double                         thQ(double In_old,double In_new,double Out_old,double Out_new,double th_in,double dx,double tstep) const;
  double                    TVDTheta(double In_old,double In_new,double Out_old,double Out_new,double th_in,double dx,double tstep) const;

public:/*-------------------------------------------------------*/
  //Constructors:
  CSubBasin(const long           ID,
            const string         Name,
            const CModelABC     *pMod,
            const long           downstream_ID, //index of downstream SB, if <0, downstream outlet
            const CChannelXSect *pChan,         //Channel
            const double         reach_len,     //reach length [m]
            const double         Qreference,    //reference flow [m3/s]
            const bool           gauged);       //if true, hydrographs are generated
  ~CSubBasin();

  //Accessor functions
  long                 GetID                () const;
  string               GetName              () const;
  double               GetBasinArea         () const;
  double               GetDrainageArea      () const;
  double               GetAvgAnnualFlow     () const;
  double               GetAvgStateVar       (const int i) const;
  double               GetAvgForcing        (const string &forcing_string) const;
  double               GetAvgCumulFlux      (const int i, const bool to) const;
  double               GetAvgCumulFluxBet   (const int iFrom, const int iTo) const;
  double               GetReferenceFlow     () const;
  long                 GetDownstreamID      () const;
  int                  GetNumHRUs           () const;
  const CHydroUnit    *GetHRU               (const int k) const;
  bool                 IsGauged             () const;
  double               GetReachLength       () const;
  int                  GetNumSegments       () const;
  bool                 IsEnabled            () const;
  double               GetRainCorrection    () const;
  double               GetSnowCorrection    () const;

  const double        *GetUnitHydrograph    () const;
  const double        *GetRoutingHydrograph () const;
  int                  GetLatHistorySize    () const;
  int                  GetInflowHistorySize () const;

  double               GetOutflowRate       () const;                   //[m3/s] from final segment, point in time
  double               GetIntegratedOutflow (const double &tstep) const;//[m3] from final segment integrated over timestep
  double             GetIntegratedSpecInflow(const double &t,
                                             const double &tstep) const;//[m3] from specified inflows integrated over timestep
  double               GetReservoirInflow   () const;                   //[m3/s] from final segment upstream of reservoir, point in time
  double               GetReservoirLosses   (const double &tstep) const;//[m3] from reservoir integrated over timestep
  double             GetReservoirEvapLosses (const double &tstep) const;//[m3] from reservoir integrated over timestep AET only
  double        GetIntegratedReservoirInflow(const double &tstep) const;//[m3] from final segment upstream of reservoir integrated over timestep
  double               GetRivuletStorage    () const;                   //[m3] volume en route to outflow
  double               GetChannelStorage    () const;                   //[m3] volume in channel
  double               GetReservoirStorage  () const;                   //[m3] volume in reservoir
  double               GetSpecifiedInflow   (const double &t) const;    //[m3/s] to upstream end of channel at point in time
  double               GetDownstreamInflow  (const double &t) const;    //[m3/s] to downstream end of channel at point in time

  CReservoir          *GetReservoir         () const;

  //Manipulator functions
  //called during model construction/assembly:
  void            AddHRU              (CHydroUnit *pHRU);
  void            AddReservoir        (CReservoir *pReservoir);
  bool            SetBasinProperties  (const string label,
                                       const double &value);
  void            SetAsNonHeadwater   ();
  double          CalculateBasinArea  ();
  void            Initialize          (const double    &Qin_avg,          //[m3/s]
                                       const double    &Qlat_avg,         //[m3/s]
                                       const double    &total_drain_area, //[km2]
                                       const optStruct &Options);
  void            AddInflowHydrograph (CTimeSeries *pInflow);
  void            AddDownstreamInflow (CTimeSeries *pInflow);
  void            AddReservoirExtract (CTimeSeries *pOutflow);
  void            AddWeirHeightTS     (CTimeSeries *pWeirHt);
  void            AddMaxStageTS       (CTimeSeries *pMaxStage);
  void            AddOverrideFlowTS   (CTimeSeries *pQoverride);
  void           AddMinStageTimeSeries(CTimeSeries *pMS);
  void       AddMinStageFlowTimeSeries(CTimeSeries *pQ);
  void        AddTargetStageTimeSeries(CTimeSeries *pTS);
  void       AddMaxQIncreaseTimeSeries(CTimeSeries *pQdelta);
  void            ResetReferenceFlow  (const double    &Qreference);
  void            SetReservoirFlow    (const double &Q,const double &t);
  void            SetInitialReservoirStage   (const double &h);
  void            SetChannelStorage   (const double &V);
  void            SetRivuletStorage   (const double &V);
  void            SetQoutArray        (const int N, const double *aQo, const double QoLast);
  void            SetQlatHist         (const int N, const double *aQl, const double QlLast);
  void            SetQinHist          (const int N, const double *aQi);
  void            SetDownstreamID     (const long down_SBID);
  void            Disable             ();
  void            Enable              ();
  double          ScaleAllFlows       (const double &scale_factor, const double &tstep);

  //called during model operation:
  void            SetInflow           (const double &Qin );//[m3/s]
  void            UpdateFlowRules     (const time_struct &tt, const optStruct &Options);
  void            UpdateOutflows      (const double *Qout_new,
                                       const double &res_ht,
                                       const double &res_outflow,
                                       const optStruct &Options,
                                       const time_struct &tt,
                                       bool initialize);//[m3/s]
  void            SetLateralInflow    (const double &Qlat);//[m3/s]

  void            RouteWater          (      double      *Qout_new,
                                             double      &res_ht,
                                             double      &res_outflow,
                                             const optStruct   &Options,
                                             const time_struct &tt) const;
  double          ChannelLosses       (const double      &reach_volume,
                                       const double      &PET,
                                       const optStruct   &Options) const;

  void            WriteMinorOutput    (const time_struct &tt) const;
  void            WriteToSolutionFile (ofstream &OUT) const;
};

#endif
