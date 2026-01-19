/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2026 the Raven Development Team
  ----------------------------------------------------------------*/
#ifndef TIMESERIESABC_H
#define TIMESERIESABC_H

#include "RavenInclude.h"
#include "ParseLib.h"

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for time series (abstract base class)

class CTimeSeriesABC
{
public:
  enum ts_type {TS_REGULAR, TS_IRREGULAR};

  double   *_aSampVal;     ///< Array of resampled time series values every timestep for model duration
  int       _nSampVal;     ///< size of aSampVal (~model_duration/timestep)
  double    _sampInterval; ///< timestep of resampled timeseries

private:/*------------------------------------------------------*/
  ts_type   _type;         ///< type - regular or irregular
  string    _name;         ///< name of time series (used only for error messages)
  long long _loc_ID;       ///< location ID (stores additional info, like HRU or SBID for observation data)
  int       _constit_ind;  ///< constituent index, if a concentration/temperature observation
  string    _srcfile;      ///< original source file
  long long _demand_ID;    ///< integer ID tag for (e.g.,) demand ID

  CTimeSeriesABC(const CTimeSeriesABC &t); //suppresses default copy constructor

public:/*-------------------------------------------------------*/
  //Constructors:
  CTimeSeriesABC(ts_type type,
                 string  name,
                 long long loc_ID,
                 string filename);
  CTimeSeriesABC(string name,
                 const CTimeSeriesABC &t);
  virtual ~CTimeSeriesABC();

  virtual void Initialize(const double model_start_day, //jul day
                          const    int model_start_year,//year
                          const double model_duration,  //days
                          const double timestep,        //days
                          const   bool is_observation,
                          const    int calendar) =0;

  virtual void InitializeResample(const int nSampVal, const double &sampInterval); //virtual for constant time series case

  ts_type   GetType      () const;
  string    GetName      () const;
  long long GetLocID     () const;
  int       GetConstitInd() const;
  string    GetSourceFile() const;
  long long GetDemandID  () const;

  void      SetLocID     (long long ID);
  void      SetConstitInd(const int c);
  void      SetDemandID  (long long demandID);

  virtual double GetInterval() const=0;
  virtual double GetTime      (const int n) const=0;
  virtual double GetValue     (const int n) const=0;
  virtual int    GetNumValues ()            const=0;
  virtual double GetAvgValue  (const double &t, const double &tstep) const=0;
  virtual double GetMinValue  (const double &t, const double &tstep) const=0;
  virtual double GetMaxValue  (const double &t, const double &tstep) const=0;

  virtual double GetSampledValue(const int nn) const;
  virtual double GetSampledTime(const int nn) const;
  virtual double GetSampledInterval() const;
  virtual int    GetNumSampledValues() const;

  virtual int    GetTimeIndexFromModelTime(const double &t_mod) const=0;

  void   SetSampledValue (const int nn,const double &val);
};

#endif
