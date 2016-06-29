/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright © 2008-2014 the Raven Development Team
------------------------------------------------------------------
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

enum curve_function{
  CURVE_LINEAR,      ///< y =a*x
  CURVE_POWERLAW,    ///< y=a*x^b
  CURVE_DATA,        ///< y=interp(xi,yi)
  CURVE_VARYING      ///< y=interp(xi,yi,t)
};
enum res_type{
  RESROUTE_STANDARD,
  RESROUTE_NONE
};
/*****************************************************************
   Class CReservoir
------------------------------------------------------------------
   Data Abstraction for reservoir 
******************************************************************/
class CReservoir
{
  private:/*-------------------------------------------------------*/
    string       name;                ///< reservoir name
    long         SBID;                ///< subbasin ID
    res_type     type;                ///< reservoir type (default: RESROUTE_STANDARD)

    const CHydroUnit  *pHRU;          ///< (potentially zero-area) HRU used for Precip/ET calculation (or NULL for no ET)
          CTimeSeries *pExtractTS;    ///< Time Series of extraction [m3/s] (or NULL for zero extraction)
    
    double       stage;               ///< current stage [m] (actual state variable)
    double       stage_last;          ///< stage at beginning of current time step [m]
    double       Qout;                ///< outflow corresponding to current stage [m3/s] 
    double       Qout_last;           ///< outflow at beginning of current time step [m3/s]
     
    double       min_stage;           ///< reference elevation [m] (below which, no volume; flow can be zero)
    double       max_stage;           ///< maximum reference elevation [m] 

    int          N;                   ///< number of points on rating curves
    double      *aStage;              ///< Rating curve for stage elevation [m]
    double      *aQ;                  ///< Rating curve for flow rates [m3/s]
    double      *aArea;               ///< Rating curve for surface area [m2]
    double      *aVolume;             ///< Rating curve for storage volume [m3]

    double     **aQ_back;             ///< Rating curve for flow rates for different times [m3/s]
    int         *aDates;              ///< Array of Julian days at which aQ changes [days after Jan 1]
    int          nDates;              ///< size of aDates

    double     GetVolume(const double &ht) const;
    double     GetArea  (const double &ht) const;
    double     GetOutflow(const double &ht) const; 

  public:/*-------------------------------------------------------*/
    //Constructors:
    CReservoir(const string name, const long SubID, const res_type typ,
               const double a_v, const double b_v, 
               const double a_Q, const double b_Q, 
               const double a_A, const double b_A);
    CReservoir(const string name, const long SubID, const res_type typ,
               const double *a_ht, 
               const double *a_Q, const double *a_A, const double *a_V, 
               const int     nPoints);
    CReservoir(const string Name, const long SubID, const res_type typ,
               const int nDates, const int *aDates, const double *a_ht, 
                     double **a_QQ, const double *a_A, const double *a_V, 
               const int     nPoints);
    ~CReservoir();

    //Accessors
    long              GetSubbasinID        () const;
    double            GetStorage           () const;//[m3]
    double            GetOutflowRate       () const;//[m3/d] 
    double            GetEvapLosses        () const;//[m3/d]
    double            GetIntegratedOutflow (const double &tstep) const;//[m3]
    double            GetStage             () const;//[m]

    //Manipulators
    void              SetMinStage            (const double &min_z);
    void              Initialize             (const optStruct &Options);
    void              UpdateStage            (const double &new_stage);
    void              UpdateFlowRules        (const time_struct &tt, const optStruct &Options);
    void              SetInitialFlow         (const double &Q);
    void	            AddExtractionTimeSeries(CTimeSeries *pOutflow);
    void              SetHRU                 (const CHydroUnit *pHRU);

    //
    double            RouteWater(const double      &Qin_old, 
                                 const double      &Qin_new, 
                                 const double      &tstep,
                                 const time_struct &tt) const;
    void              WriteToSolutionFile (ofstream &OUT) const;

    static CReservoir  *Parse(CParser *p, string name, int &HRUID, const optStruct &Options);
};
#endif