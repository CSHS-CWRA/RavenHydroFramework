/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2024 the Raven Development Team
  ----------------------------------------------------------------
  Reservoir.h
  ------------------------------------------------------------------
  defines abstract reservoir
  ----------------------------------------------------------------*/
#ifndef RESERVOIR_H
#define RESERVOIR_H

#include "RavenInclude.h"
#include "ParseLib.h"
#include "HydroUnits.h"
#include "TimeSeries.h"
#include "SubBasin.h"
#include "ControlStructures.h"

class CSubBasin;
enum curve_function{
  CURVE_LINEAR,      ///< y =a*x
  CURVE_POWERLAW,    ///< y=a*x^b
  CURVE_DATA,        ///< y=interp(xi,yi)
  CURVE_VARYING,     ///< y=interp(xi,yi,t)
  CURVE_LAKE         ///< y=interp(xi,yi) (curve from weir)
};
struct DZTRmodel //Dynamically zoned target release model structure
{
  double Qmc;        ///< maximum channel capacity [m3/s]
  double Vmax;       ///< maximum reservoir storage [m3]
  double Qci[12];    ///< critical release target [m3/s]
  double Qni[12];    ///< normal release target [m3/s]
  double Qmi[12];    ///< maximum release target [m3/s]
  double Vci[12];    ///< critical release storage [m3]
  double Vni[12];    ///< normal release storage [m3]
  double Vmi[12];    ///< critical release storage [m3]
};
struct down_demand {
  const CSubBasin *pDownSB;   ///< pointer to subbasin of downstream demand location
  double percent;             ///< percentage of demand met by reservoir
  int    julian_start;        ///< julian start day for commencement of demand (beginning @ 0)
  int    julian_end;          ///< julian end day of demand (wraps, such that if julian_end < julian_start, demand in winter)
};
class CSubBasin;
class CControlStructure;
/*****************************************************************
   Class CReservoir
------------------------------------------------------------------
   Data Abstraction for reservoir
******************************************************************/
class CReservoir
{
private:/*-------------------------------------------------------*/
  string       _name;                ///< reservoir name
  long         _SBID;                ///< subbasin ID
  double       _max_capacity;        ///< maximum reservoir capacity [m3]

  double       _lakebed_thick;       ///< lakebed thickness [m]
  double       _lakebed_cond;        ///< lakebed thermal conductivity [MJ/m/K/d]
  double       _lake_convcoeff;      ///< lake surface thermal convection coefficient [MJ/m2/d/K]

  const CHydroUnit  *_pHRU;          ///< (potentially zero-area) HRU used for Precip/ET calculation (or NULL for no ET)

  CTimeSeries **_pDemandTS;          ///< array of Time Series of water demand [m3/s] (or NULL for zero extraction) [size: _nDemandTS]
  int           _nDemandTS;          ///< number of reservoir demand time series

  CTimeSeries *_pWeirHeightTS;       ///< Time series of weir heights [m] (or NULL for fixed weir height)
  CTimeSeries *_pMaxStageTS;         ///< Time series of rule curve upper stage constraint [m] (or NULL for no maximum stage)
  CTimeSeries *_pOverrideQ;          ///< Time series of overrride flows [m3/s] (or NULL for no override)
  CTimeSeries *_pMinStageTS;         ///< Time series of minimum stage [m] (or NULL for no minimum)
  CTimeSeries *_pMinStageFlowTS;     ///< Time series of flows when below minimum stage [m3/s] (or NULL for no minimum)
  CTimeSeries *_pTargetStageTS;      ///< Time series of target stage [m] (or NULL for no target)
  CTimeSeries *_pMaxQIncreaseTS;     ///< Time series maximum rate of flow increase [m3/s/d] (or NULL for no maximum)
  CTimeSeries *_pMaxQDecreaseTS;     ///< Time series maximum rate of flow decrease [m3/s/d] (or NULL for no maximum)
  CTimeSeries *_pDroughtLineTS;      ///< Time series of drought line stage [m] (or NULL if none exists)
  CTimeSeries *_pQmaxTS;             ///< Time series of maximum flow constraint [m3/s] (or NULL for no maximum)
  CTimeSeries *_pQminTS;             ///< Time series of minimum flow constraint [m3/s] (or NULL for no minimum)
  CTimeSeries *_pQdownTS;            ///< Time series of downstream flow soft target [m3/s] (or NULL for none)
  double       _QdownRange;          ///< range of acceptable target flows from Qdown [m3/s] (or zero for hard target)

  double       _Qoptimized;          ///< proscribed outflow provided by management optimization routine [m3/s] (or RAV_BLANK_DATA for none)

  const CSubBasin *_pQdownSB;        ///< pointer to downstream SubBasin for diversions (or NULL for none)
  const CSubBasin *_pDownSB;         ///< pointer to downstream subbasin

  down_demand**_aDemands;            ///< array of pointers to downstream demand information used to determine Qmin [size:_nDemands]
  int          _nDemands;            ///< size of downstream demand location array

  DZTRmodel    *_pDZTR;              ///< pointer to DZTR model, if used (default=NULL)

  bool          _minStageDominant;   ///< true if minimum stage dominates minflow/overrideflow constraints (false by default)
  double        _demand_mult;        ///< reservoir demand multiplier that indicates percentage of requested downstream irrigation demand
                                     ///< satisfied from this reservoir.

  bool          _assimilate_stage;   ///< true if assimilating lake stage for this reservoir
  const CTimeSeriesABC *_pObsStage;  ///< observed lake stage
  bool          _assim_blank;        ///< true if observed stage is blank this time step

  double        _DAscale;            //< outflow scale factor - used for reporting overriden flows
  double        _DAscale_last;       //< outflow scale factor for previous time step

  int           _dry_timesteps;      //< number of time steps this reservoir dried out  during simulation

  //state variables :
  double       _stage;               ///< current stage [m] (actual state variable)
  double       _stage_last;          ///< stage at beginning of current time step [m]
  double       _Qout;                ///< outflow corresponding to current stage [m3/s]
  double       _Qout_last;           ///< outflow at beginning of current time step [m3/s]
  double       _MB_losses;           ///< losses over current time step [m3]
  double       _AET;                 ///< losses through AET only [m3]
  double       _Precip;              ///< gains through precipitation [m3]
  double       _GW_seepage;          ///< losses to GW only [m3] (negative for GW gains)
  double      *_aQstruct;            ///< array of flows from control structures at end of time step [m3/s]
  double      *_aQstruct_last;       ///< array of flows from control structures at start of time step [m3/s]
  double      *_aQdelivered;         ///< amount of water demand delivered for each demand [m3/s] (for management optimization)

  res_constraint _constraint;        ///< current constraint type

  //rating curve characteristics :
  double       _min_stage;           ///< reference elevation [m] (below which, no volume; flow can be zero)
  double       _max_stage;           ///< maximum reference elevation [m]

  double       _crest_ht;            ///< absolute crest height, m.a.s.l (0 by default)
  double       _crest_width;         ///< crest width [m]

  int          _Np;                   ///< number of points on rating curves
  double      *_aStage;              ///< Rating curve for stage elevation [m]
  double      *_aQ;                  ///< Rating curve for overflow (e.g., weir) flow rates [m3/s]
  double      *_aQunder;             ///< Rating curve for underflow/orifice flow (if specified; 0 by default) [m3/s]
  double      *_aArea;               ///< Rating curve for surface area [m2]
  double      *_aVolume;             ///< Rating curve for storage volume [m3]

  double     **_aQ_back;             ///< Rating curve for flow rates for different times [m3/s]
  int         *_aDates;              ///< Array of Julian days at which aQ changes [days after Jan 1]
  int          _nDates;              ///< size of aDates

  int                 _nControlStructures; ///< number of outflow control structures
  CControlStructure **_pControlStructures; ///< Array of pointer to outflow control structures

  //GW information :
  double       _seepage_const;       ///< seepage constant [m3/s/m] for groundwater losses Q=k*(h-h_loc)
  double       _local_GW_head;       ///< local head [m] (same  for groundwater losses Q=k*(h-h_loc)

  void       BaseConstructor(const string Name,const long SubID); //because some versions of c++ don't like delegating constructors

  double     GetVolume (const double &ht) const;
  double     GetArea   (const double &ht) const;
  double     GetWeirOutflow(const double &ht, const double &adj) const;

  double     GetDZTROutflow(const double &V,const double &Qin,const time_struct &tt,const optStruct &Options) const;

  void       MultiplyFlow(const double &mult);

public:/*-------------------------------------------------------*/
  //Constructors:
  CReservoir(const string Name, const long SubID);
  CReservoir(const string name, const long SubID, //Power law constructor
             const double a_v, const double b_v,
             const double a_Q, const double b_Q,
             const double a_A, const double b_A,
             const double crestht, const double depth);
  CReservoir(const string name, const long SubID,
             const double *a_ht,
             const double *a_Q, const double *aQ_und,const double *a_A, const double *a_V,
             const int     nPoints);
  CReservoir(const string Name, const long SubID,
             const int nDates, const int *aDates, const double *a_ht,
             double **a_QQ, const double *aQ_und, const double *a_A, const double *a_V,
             const int     nPoints);
  CReservoir(const string Name, const long SubID, const double weircoeff, //Lake constructor
             const double crestw, const double crestht, const double A, const double depth);
  ~CReservoir();

  //Accessors
  long              GetSubbasinID            () const;

  double            GetStorage               () const; //[m3]
  double            GetOldStorage            () const; //[m3]
  double            GetOutflowRate           () const; //[m3/s]
  double            GetOldOutflowRate        () const; //[m3/s]
  double            GetIntegratedOutflow     (const double &tstep) const; //[m3]
  double            GetControlOutflow        (const int i) const; //[m3]
  double            GetIntegratedControlOutflow(const int i, const double& tstep) const; //[m3]
  double            GetReservoirLosses       (const double &tstep) const; //[m3]
  double            GetReservoirEvapLosses   (const double &tstep) const; //[m3]
  double            GetReservoirGWLosses     (const double &tstep) const; //[m3]
  double            GetReservoirPrecipGains  (const double &tstep) const; //[m3]
  double            GetResStage              () const; //[m]
  double            GetOldStage              () const; //[m]
  double            GetSurfaceArea           () const; //[m2]
  double            GetOldSurfaceArea        () const; //[m2]
  double            GetLakebedThickness      () const; //[m]
  double            GetLakebedConductivity   () const; //[MJ/m/K/d]
  double            GetLakeConvectionCoeff   () const; //[MJ/m2/d/K]
  double            GetDischargeFromStage    (const double &stage, const int nn) const; //[m3/s]
  double            GetStageDischargeDerivative(const double &stage,const int nn) const; //[m3/s/d]
  double            GetMinStage              (const int nn) const;//[m]
  double            GetMaxStage              (const int nn) const;//[m]
  double            GetSillElevation         (const int nn) const;//[m]

  int               GetNumWaterDemands       () const;
  int               GetWaterDemandID         (const int ii) const;
  string            GetWaterDemandName       (const int ii) const;
  double            GetWaterDemand           (const int ii,const double &t) const;  //[m3/s] iith demand from reservoir at point in time
  double            GetDemandDelivery        (const int ii) const;       //[m3/s] instantaneous delivery rate to demand ii

  int               GetHRUIndex              () const;
  double            GetMaxCapacity           () const; //[m3]
  string            GetCurrentConstraint     () const;
  double            GetDemandMultiplier      () const;
  double            GetCrestWidth            () const;
  int               GetNumControlStructures  () const;
  double            GetAET                   () const; //[mm/d]
  string            GetRegimeName            (const int i, const time_struct &tt) const;
  long              GetControlFlowTarget     (const int i) const;
  string            GetControlName           (const int i) const;

  int               GetNumDryTimesteps       () const;

  //Manipulators
  void              SetMinStage              (const double &min_z);
  void              SetMaxCapacity           (const double &max_cap);
  void              Initialize               (const optStruct &Options);
  void              SetInitialFlow           (const double &Q,const double &Qlast,const time_struct &tt, const optStruct &Options);
  void              SetReservoirStage        (const double &ht, const double &ht_last);
  void              SetControlFlow           (const int i, const double &Q, const double &Qlast);
  void              SetOptimizedOutflow      (const double &Qout);
  void              SetVolumeStageCurve      (const double *a_ht,const double *a_V,const int nPoints);
  void              SetAreaStageCurve        (const double *a_ht,const double *a_A,const int nPoints);
  void              SetGWParameters          (const double &coeff, const double &h_ref);
  void              SetCrestWidth            (const double &width);
  void              SetDataAssimFactors      (const double &da_scale, const double &da_scale_last);
  void              TurnOnAssimilation       (CTimeSeriesABC *pObs);
  void              SetPrecip                (const double &precip_m3);
  void              SetDownstreamBasin       (const CSubBasin *pDownBasin);
  void              SetLakebedThickness      (const double &thick);
  void              SetLakebedConductivity   (const double &cond);
  void              SetLakeConvectionCoeff   (const double &conv);

  void              AddDemandTimeSeries      (CTimeSeries *pOutflow);
  void              AddWeirHeightTS          (CTimeSeries *pWeirHt);
  void              AddMaxStageTimeSeries    (CTimeSeries *pMS);
  void              AddOverrideQTimeSeries   (CTimeSeries *pQ);
  void              AddMinStageTimeSeries    (CTimeSeries *pMS);
  void              AddMinStageFlowTimeSeries(CTimeSeries *pQ);
  void              AddTargetStageTimeSeries (CTimeSeries *pTS);
  void              AddMaxQIncreaseTimeSeries(CTimeSeries *pQdelta);
  void              AddMaxQDecreaseTimeSeries(CTimeSeries *pQdelta);
  void              AddDroughtLineTimeSeries (CTimeSeries *pTS);
  void              AddMinQTimeSeries        (CTimeSeries *pQmin);
  void              AddMaxQTimeSeries        (CTimeSeries *pQmax);
  void              AddDownstreamTargetQ     (CTimeSeries *pQ, const CSubBasin *pSB, const double &range);

  void              AddDownstreamDemand      (const CSubBasin *pSB,const double pct, const int julian_start, const int julian_end);
  void              AddControlStructure      (const CControlStructure *pCS);

  void              SetDZTRModel             (const double Qmc, const double Smax,
                                              const double Sci[12], const double Sni[12],const double Smi[12],
                                              const double Qci[12], const double Qni[12],const double Qmi[12]);
  void              SetMinStageDominant      ();
  void              SetDemandMultiplier      (const double &value);

  void              SetHRU                   (const CHydroUnit *pHRU);
  void              DisableOutflow           ();
  void              ClearTimeSeriesData      (const optStruct& Options);

  //Called during simulation:
  double            RouteWater               (const double      &Qin_old,
                                              const double      &Qin_new,
                                              const CModelABC*  pModel,
                                              const optStruct   &Options,
                                              const time_struct &tt,
                                                    double      &res_outflow,
                                                 res_constraint &constraint,
                                                    double      *aQstruct) const;
  void              UpdateStage              (const double      &new_stage,
                                              const double      &new_ouflow,
                                              const res_constraint &constraint,
                                              const double      *new_struct_flows,
                                              const optStruct   &Options,
                                              const time_struct &tt);
  void              WriteToSolutionFile      (ofstream &OUT) const;
  void              UpdateReservoir          (const time_struct &tt, const optStruct &Options);
  void              UpdateMassBalance        (const time_struct &tt, const double &tstep, const optStruct &Options);
  double            ScaleFlow                (const double &scale, const bool overriding,const double &tstep,const double &t);
  void              AddToDeliveredDemand     (const int ii, const double &Q);
};
#endif
