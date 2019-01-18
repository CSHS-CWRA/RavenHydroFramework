/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
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
  enum ts_type {ts_regular, ts_irregular};

private:/*------------------------------------------------------*/
  ts_type     _type;
  string      _name; ///< name of time series (used only for error messages)
  string       _tag; ///< data tag (stores additional info, like HRU or SB ID for observation data)
  string   _srcfile; ///< original source file

  CTimeSeriesABC(const CTimeSeriesABC &t); //suppresses default copy constructor

public:/*-------------------------------------------------------*/
  //Constructors:
  CTimeSeriesABC(ts_type type,
                 string  name,
                 string  tag,
                 string filename="");
  CTimeSeriesABC(string name,
                 const CTimeSeriesABC &t);
  virtual ~CTimeSeriesABC();

  virtual void Initialize(const double model_start_day, //jul day
                          const    int model_start_year,//year
                          const double model_duration,  //days
                          const double timestep,        //days
                          const   bool is_observation,
                          const    int calendar) =0;

  ts_type GetType      () const;
  string  GetName      () const;
  string  GetTag       () const;
  string  GetSourceFile() const;

  void    SetTag       (string tag){_tag=tag;}

  virtual double GetInterval() const=0;
  virtual double GetTime      (const int n) const=0;
  virtual double GetValue     (const int n) const=0;
  virtual int    GetNumValues ()            const=0;
  virtual double GetAvgValue  (const double &t, const double &tstep) const=0;
  virtual double GetMinValue  (const double &t, const double &tstep) const=0;
  virtual double GetMaxValue  (const double &t, const double &tstep) const=0;

  virtual double GetSampledValue(const int nn) const=0;
  virtual double GetSampledTime(const int nn) const=0;
  virtual double GetSampledInterval() const=0;
  virtual int    GetNumSampledValues() const=0;

};

#endif
