/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2019 the Raven Development Team
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
class CSubBasin;
enum curve_function{
  CURVE_LINEAR,      ///< y =a*x
  CURVE_POWERLAW,    ///< y=a*x^b
  CURVE_DATA,        ///< y=interp(xi,yi)
  CURVE_VARYING,     ///< y=interp(xi,yi,t)
  CURVE_LAKE         ///< y=interp(xi,yi) (curve from weir)
};
enum res_type{
  RESROUTE_STANDARD,
  RESROUTE_NONE
};
enum res_constraint {
  RC_MAX_STAGE,
  RC_MIN_STAGE,
  RC_NATURAL, 
  RC_TARGET, 
  RC_DOWNSTREAM_FLOW,
  RC_MAX_FLOW_INCREASE,
  RC_MAX_FLOW_DECREASE,
  RC_MIN_FLOW,
  RC_MAX_FLOW,
  RC_OVERRIDE_FLOW,
  RC_DRY_RESERVOIR
};
struct down_demand {
  const CSubBasin *pDownSB;   ///< pointer to subbasin of downstream demand location
  double percent;             ///< percentage of demand met by reservoir
};
class CSubBasin;
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
  res_type     _type;                ///< reservoir type (default: RESROUTE_STANDARD)
  double       _max_capacity;        ///< maximum reservoir capacity [m3]

  const CHydroUnit  *_pHRU;          ///< (potentially zero-area) HRU used for Precip/ET calculation (or NULL for no ET)

  CTimeSeries *_pExtractTS;          ///< Time Series of extraction [m3/s] (or NULL for zero extraction)

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
  
  const CSubBasin *_pQdownSB;        ///< pointer to downstream SubBasin (or NULL for none)

  down_demand**_aDemands;            ///< array of pointers to downstream demand information used to determine Qmin [size:_nDemands]
  int          _nDemands;            ///< size of downstream demand location array  

  //state variables :
  double       _stage;               ///< current stage [m] (actual state variable)
  double       _stage_last;          ///< stage at beginning of current time step [m]
  double       _Qout;                ///< outflow corresponding to current stage [m3/s]
  double       _Qout_last;           ///< outflow at beginning of current time step [m3/s]
  double       _MB_losses;           ///< losses over current time step [m3]
  double       _AET;                 ///< losses through AET only [m3]
  double       _GW_seepage;          ///< losses to GW only [m3] (negative for GW gains)

  res_constraint _constraint;          ///< current constraint type


  //rating curve characteristics :
  double       _min_stage;           ///< reference elevation [m] (below which, no volume; flow can be zero)
  double       _max_stage;           ///< maximum reference elevation [m]

  double       _crest_ht;            ///< absolute crest height, m.a.s.l (0 by default)

  int          _Np;                   ///< number of points on rating curves
  double      *_aStage;              ///< Rating curve for stage elevation [m]
  double      *_aQ;                  ///< Rating curve for overflow (e.g., weir) flow rates [m3/s]
  double      *_aQunder;             ///< Rating curve for underflow/orifice flow (if specified; 0 by default) [m3/s]
  double      *_aArea;               ///< Rating curve for surface area [m2]
  double      *_aVolume;             ///< Rating curve for storage volume [m3]

  double     **_aQ_back;             ///< Rating curve for flow rates for different times [m3/s]
  int         *_aDates;              ///< Array of Julian days at which aQ changes [days after Jan 1]
  int          _nDates;              ///< size of aDates

  //GW information :
  double       _seepage_const;       ///< seepage constant [m3/s/m] for groundwater losses Q=k*(h-h_loc)  
  double       _local_GW_head;       ///< local head [m] (same  for groundwater losses Q=k*(h-h_loc)

  void       BaseConstructor(const string Name,const long SubID,const res_type typ); //because some versions of c++ don't like delegating constructors

  double     GetVolume (const double &ht) const;
  double     GetArea   (const double &ht) const;
  double     GetOutflow(const double &ht, const double &adj) const;

public:/*-------------------------------------------------------*/
  //Constructors:
  CReservoir(const string Name,const long SubID,const res_type typ);
  CReservoir(const string name, const long SubID, const res_type typ,
             const double a_v, const double b_v,
             const double a_Q, const double b_Q,
             const double a_A, const double b_A);
  CReservoir(const string name, const long SubID, const res_type typ,
             const double *a_ht,
             const double *a_Q, const double *aQ_und,const double *a_A, const double *a_V,
             const int     nPoints);
  CReservoir(const string Name, const long SubID, const res_type typ,
             const int nDates, const int *aDates, const double *a_ht,
             double **a_QQ, const double *aQ_und, const double *a_A, const double *a_V,
             const int     nPoints);
  CReservoir(const string Name, const long SubID, const res_type typ, const double weircoeff, //Lake constructor
             const double crestw, const double crestht, const double A, const double depth);
  ~CReservoir();

  //Accessors
  long              GetSubbasinID            () const;
  double            GetStorage               () const; //[m3]
  double            GetOutflowRate           () const; //[m3/d]
  double            GetReservoirLosses       (const double &tstep) const; //[m3]
  double            GetReservoirEvapLosses   (const double &tstep) const; //[m3]
  double            GetReservoirGWLosses     (const double &tstep) const; //[m3]
  double            GetIntegratedOutflow     (const double &tstep) const; //[m3]
  double            GetResStage              () const; //[m]
  double            GetOldOutflowRate        () const; //[m3/d]
  double            GetOldStorage            () const; //[m3]
  int               GetHRUIndex              () const;
  double            GetMaxCapacity           () const; //[m3]
  string            GetCurrentConstraint     () const;

  //Manipulators
  void              SetMinStage              (const double &min_z);
  void              SetMaxCapacity           (const double &max_cap);
  void              Initialize               (const optStruct &Options);
  void              SetInitialFlow           (const double &Q,const double &Qlast,const double &t);
  void              SetInitialStage          (const double &ht, const double &ht_last);
  void              SetVolumeStageCurve      (const double *a_ht,const double *a_V,const int nPoints);
  void              SetGWParameters          (const double &coeff, const double &h_ref);

  void              AddExtractionTimeSeries  (CTimeSeries *pOutflow);
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

  void              AddDownstreamDemand      (const CSubBasin *pSB,const double pct);

  void              SetHRU                   (const CHydroUnit *pHRU);
  void              DisableOutflow           ();

  //Called during simulation:
  double            RouteWater               (const double      &Qin_old,
                                              const double      &Qin_new,
                                              const double      &tstep,
                                              const time_struct &tt,
                                                    double      &res_outflow,
                                                 res_constraint &constraint) const;
  void              UpdateStage              (const double      &new_stage,
                                              const double      &new_ouflow, 
                                              const res_constraint &constraint,
                                              const optStruct   &Options,
                                              const time_struct &tt);
  void              WriteToSolutionFile      (ofstream &OUT) const;
  void              UpdateFlowRules          (const time_struct &tt, const optStruct &Options);
  void              UpdateMassBalance        (const time_struct &tt, const double &tstep);

  static CReservoir  *Parse(CParser *p, string name, int &HRUID, const optStruct &Options);
};
#endif
