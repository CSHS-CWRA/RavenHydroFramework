/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------*/
#ifndef TIMESERIES_H
#define TIMESERIES_H

#include "TimeSeriesABC.h"
#include "RavenInclude.h"
#include "ParseLib.h"
#include "Forcings.h"

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for a continuous time series recorded by a gauge
/// \details Data Abstraction for time series conceptualized as a set of step
///   functions starting at time[i] ending at time[i+1] with value val[i].
///   time series will return val[N-1] for all times past t[N-1]
/// \note The series should be ordered, i.e., t[i+1]>t[i] always
///  Time units should be in days
/// \remark Used for forcing functions: precip, temp, pressure, etc.
/// \remark Storage, even with resampling to model time step is not bad.
/// \remark For hourly time step data & model, 1 gauge = ~160KB per year of simulation+~140KB per year of data
//
class CTimeSeries: public CTimeSeriesABC
{
private:/*------------------------------------------------------*/

  double _start_day; ///< Day corresponding to local TS time 0.0 (beginning of time series)
  int   _start_year; ///< Year corresponding to local TS time 0.0 (beginning of time series)

  double  _interval; ///< uniform interval between data points (in days)
  double     *_aVal; ///< Array of magnitude of pulse (variable units)
  int      _nPulses; ///< number of pulses (total duration=(nPulses-1)*_interval)

  double *_aSampVal; ///< Array of resampled time series values every timestep for model duration
  int     _nSampVal; ///< size of aSampVal (~model_duration/timestep)
  double  _sampInterval; ///< timestep of resampled timeseries

  bool   _sub_daily; ///< true if smallest time interval is sub-daily

  double    _t_corr; ///< number of days between model start date and gauge start date (positive if data exists before model start date)
  ///< \brief correction from model time (t) to time series/local time

  bool       _pulse; ///< flag determining whether this is a pulse-based or piecewise-linear time series
  ///< \remark forcing functions are all pulse-based

  int     GetTimeIndex(const double &t_loc) const;

  void        Resample(const double &tstep,          //days
                       const double &model_duration);//days

  CTimeSeries(const CTimeSeries &t); //suppresses default copy constructor

public:/*-------------------------------------------------------*/
  //Constructors:
  CTimeSeries(string name,
              string tag,
              double one_value);
  CTimeSeries(string name,
              string tag,
              string filename,
              double start_day,
              int start_yr,
              double interval, //in days
              double *aValues,
              const int NumValues,
              const bool is_pulse_type);
  CTimeSeries(string name,
              string tag,
              string filename,
              double start_day,
              int start_yr,
              double interval, //in days
              const int NumValues,
              const bool is_pulse_type);
  CTimeSeries(string name,
              const CTimeSeries &t);
  ~CTimeSeries();

  void Initialize(const double model_start_day, //jul day
                  const    int model_start_year,//year
                  const double model_duration,  //days
                  const double timestep,        //days
                  const   bool is_observation);

  void        InitializeResample(const int nSampVal, const double sampInterval);

  bool   IsDaily      () const;
  int    GetStartYear () const;
  double GetStartDay  () const;
  double GetInterval  () const;

  double GetValue     (const double &t) const;
  double GetAvgValue  (const double &t, const double &tstep) const;
  double GetMinValue  (const double &t, const double &tstep) const;
  double GetMaxValue  (const double &t, const double &tstep) const;
  int    GetNumValues () const;
  double GetTime      (const int n) const;
  double GetValue     (const int n) const;

  double GetSampledValue(const int nn) const; //nn is timestep number
  double GetSampledTime(const int nn) const; //nn is timestep number
  double GetSampledInterval() const;
  int    GetNumSampledValues() const;

  double GetDailyAvg    (const int model_day) const;
  double GetDailyMin    (const int model_day) const;
  double GetDailyMax    (const int model_day) const;

  bool   IsPulseType()  const;

  static CTimeSeries  *Sum          (CTimeSeries *pTS1, CTimeSeries *pTS2, string name);
  static CTimeSeries  *Parse        (CParser *p, bool is_pulse, string name, string tag, const optStruct &Options, bool shift_to_per_ending=false);
  static CTimeSeries **ParseMultiple(CParser *p, int &nTS, forcing_type *aType, bool is_pulse);
  static CTimeSeries **ParseEnsimTb0(string filename, int &nTS, forcing_type *aType);

  void   Multiply       (const double &factor);

  double GetModelledValue(const double &t,const ts_type type) const;
  void   SetValue(const int n, const double &val);
  void   SetSampledValue(const int nn, const double &val);
};

#endif
