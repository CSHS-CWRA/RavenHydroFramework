/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ----------------------------------------------------------------*/
#ifndef _DIAGNOSTICS_H
#define _DIAGNOSTICS_H

#include "RavenInclude.h"
#include "TimeSeries.h"

enum diag_type {
  DIAG_NASH_SUTCLIFFE,
  DIAG_DAILY_NSE,
  DIAG_RMSE,
  DIAG_PCT_BIAS,
  DIAG_ABSERR,
  DIAG_ABSMAX,
  DIAG_PDIFF,
  DIAG_TMVOL,
  DIAG_RCOEF,
  DIAG_NSC,
  DIAG_RSR,
  DIAG_R2,
  DIAG_CUMUL_FLOW,
  DIAG_LOG_NASH,
  DIAG_KLING_GUPTA,
  DIAG_NASH_SUTCLIFFE_DER,
  DIAG_RMSE_DER,
  DIAG_KLING_GUPTA_DER,
  DIAG_NASH_SUTCLIFFE_RUN,
  DIAG_MBF
};


///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for time series comparison diagnostics
class CDiagnostic
{
private:/*------------------------------------------------------*/

  diag_type   _type;    ///< diagnostic type
  int         _width;   ///< moving window width (in timesteps)

public:/*------------------------------------------------------*/

	CDiagnostic(const diag_type  typ);
  CDiagnostic(const diag_type  typ, const int wid);
  ~CDiagnostic();

  string    GetName() const;
  diag_type GetType() const;

  double CalculateDiagnostic(CTimeSeriesABC  *pTSmod, 
                             CTimeSeriesABC  *pTSObs, 
                             CTimeSeriesABC  *pTSWeights, 
                             const double    &starttime,
                             const double    &endtime,
                             const optStruct &Options) const;
};

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for diagnostic period
class CDiagPeriod
{
private:/*------------------------------------------------------*/

  string    _name;    ///< period name (e.g., "CALIBRATION")
  double    _t_start; ///< starttime (in local model time)
  double    _t_end;   ///< endtime (in local model time)

public:/*------------------------------------------------------*/

  CDiagPeriod(string name, string startdate, string enddate, const optStruct &Options);
  ~CDiagPeriod();

  string   GetName()      const;
  double   GetStartTime() const;
  double   GetEndTime()   const;
};
#endif
