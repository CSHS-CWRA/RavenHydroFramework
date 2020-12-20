/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2020 the Raven Development Team
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
  int         _indexCorr; ///< index of first value within model time
  int         _nSampVal;  ///< number of values within model time

public:/*-------------------------------------------------------*/
  //Constructors:
  CIrregularTimeSeries(string name,
                       long   loc_ID,
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

  double GetSampledValue(const int nn) const; //nn is timestep number
  double GetSampledTime(const int nn) const; //nn is timestep number
  double GetSampledInterval() const;
  int    GetNumSampledValues() const;

  int    GetTimeIndexFromModelTime(const double &t_mod) const;

  static CIrregularTimeSeries  *Parse (CParser *p, const string name, const long loc_ID, const int nMeasurements);
};

#endif
