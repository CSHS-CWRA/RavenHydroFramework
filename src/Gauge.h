/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------*/

#ifndef GAUGE_H
#define GAUGE_H

#include "RavenInclude.h"
#include "TimeSeries.h"
#include "Forcings.h"

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction class for atmospheric measurement locations
//
class CGauge
{
private:/*------------------------------------------------------*/

  string        _name;           ///< Gauge name
  location      _Loc;            ///< Gauge Location
  double        _elevation;      ///< Gauge Elevation [masl]
  double        _meas_ht;        ///< Gauge Measurement height above land surface [m]

  CTimeSeries **_pTimeSeries;    ///< Array of pointers to time series
  int           _nTimeSeries;    ///< number of time series linked to gauge
  int           _aTSindex[MAX_FORCING_TYPES]; ///< lookup array for time series by forcing type

  double        _aAveTemp[12];   ///< representative average monthly temperatures [C]
  double        _aMinTemp[12];   ///< representative minimum monthly temperatures [C]
  double        _aMaxTemp[12];   ///< representative maximum monthly temperatures [C]
  double        _aAvePET [12];   ///< representative average monthly PET [mm/d] (or monthly PET factor [mm/d/K], if MONTHLY_FACTOR is used)

  double        _rainfall_corr;  ///< correction factor for rainfall (stored with gauge, used elsewhere)
  double        _snowfall_corr;  ///< correction factor for snowfall (stored with gauge, used elsewhere)
  double        _cloud_min_temp; ///< minimum temperature threshold used to determine cloud_cover factor
  double        _cloud_max_temp; ///< maximum temparature threshold used to determine cloud_cover factor

  double       DailyTempCorrection              (const double t) const; ///< Daily temperature correction [C]
  void         GenerateMinMaxAveTempFromSubdaily(const optStruct &Options);
  void         GenerateAveSubdailyTempFromMinMax(const optStruct &Options);
  void         GenerateMinMaxSubdailyTempFromAve(const optStruct &Options);

  CTimeSeries *GetTimeSeries(const forcing_type ftype) const;

  void         WarnAboutForcing(bool is_needed, const forcing_type ftype) const;

public:/*-------------------------------------------------------*/
  //Constructors:
  CGauge(      string gauge_name,
               const double latit,
               const double longit,
               const double elev);
  ~CGauge();

  void   Initialize(const optStruct  &Options,
                    const int         UTM_zone);

  //Accessor functions
  string   GetName            () const;
  location GetLocation        () const;
  double   GetElevation       () const;
  double   GetMeasurementHt   () const;
  double   GetRainfallCorr    () const;
  double   GetSnowfallCorr    () const;
  double   GetCloudMinRange   () const;
  double   GetCloudMaxRange   () const;
  double   GetGaugeProperty   (const string &pname) const;
  bool     TimeSeriesExists   (const forcing_type ftype) const;

  //special accessors (built from multiple time series):
  double   GetAverageSnowFrac (const double &t, const double &tstep) const; //snow & rain
  double   GetAverageSnowFrac (const int nn) const;

  double   GetForcingValue    (const forcing_type ftype, const double &t, const double &tstep) const;
  double   GetForcingValue    (const forcing_type ftype, const int nn) const;

  double   GetMonthlyAveTemp  (const int month) const;
  double   GetMonthlyMinTemp  (const int month) const;
  double   GetMonthlyMaxTemp  (const int month) const;
  double   GetMonthlyAvePET   (const int month) const;
  //...

  //Manipulator functions (used in model creation/parsing)
  void     SetLatitude        (const double &l);
  void     SetLongitude       (const double &l);
  void     SetElevation       (const double &e);
  void     SetMeasurementHt   (const double &ht);
  bool     SetGaugeProperty   (const string prop_tag, const double &value);

  void     AddTimeSeries        (CTimeSeries *pTS, const forcing_type ftype);

  void     SetMonthlyAveTemps   (const double temps[12]);
  void     SetMonthlyMinTemps   (const double temps[12]);
  void     SetMonthlyMaxTemps   (const double temps[12]);
  void     SetMonthlyPET        (const double PET[12]);
};
#endif
