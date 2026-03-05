/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2026 the Raven Development Team
  ----------------------------------------------------------------*/
#ifndef IRREGULARTIMESERIES_H
#define IRREGULARTIMESERIES_H

#include "TimeSeriesABC.h"
#include "RavenInclude.h"
#include "ParseLib.h"

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for discontinuous irregularly spaced time series
/// \details Data Abstraction for time series conceptualized as instantaneous values
/// \note The series should be ordered, i.e., t[i+1]>t[i] always
///  Time units should be in days
//
class CIrregularTimeSeries: public CTimeSeriesABC
{
private:/*------------------------------------------------------*/
  int        *_aYears;    ///< Array of years
  double     *_aDays;     ///< Array of julian days
  double     *_aTimes ;   ///< Array of model times
  double     *_aVal;      ///< Array of values
  int         _nVals;     ///< number of values
  int         _firstIndex;///< index of first value after model start time (or DOESNT_EXIST if all observations prior to model start)

  void Resample(const double &tstep, const double &model_duration);

public:/*-------------------------------------------------------*/
  //Constructors:
  CIrregularTimeSeries(string name,
                       long long loc_ID,
                       string filename,
                       double *aValues,
                       double *aDays,
                       int    *aYears,
                       const int NumValues);
  CIrregularTimeSeries(string name,
                       const CIrregularTimeSeries &t);
  ~CIrregularTimeSeries();

  void Initialize(const double model_start_day, //jul day
                  const    int model_start_year,//year
                  const double model_duration,   //days
                  const double timestep,         //days
                  const   bool is_observation,
                  const    int calendar);

  int    GetNumValues() const;
  double GetTime      (const int n) const;
  double GetValue     (const int n) const;
  double GetDay       (const int n) const;
  int    GetYear      (const int n) const;
  double GetInterval() const;

  double GetAvgValue  (const double &t, const double &tstep) const;
  double GetMinValue  (const double &t, const double &tstep) const;
  double GetMaxValue  (const double &t, const double &tstep) const;

  int    GetTimeIndexFromModelTime(const double &t_mod) const;

  static CIrregularTimeSeries  *Parse (CParser *p, const string name, const long long loc_ID, const int nMeasurements);
};

#endif
