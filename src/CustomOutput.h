/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2023 the Raven Development Team
  ----------------------------------------------------------------
  Custom output generation
  ----------------------------------------------------------------*/
#ifndef _CUSTOM_OUTPUT_H
#define _CUSTOM_OUTPUT_H

#include "RavenInclude.h"
#include "ModelABC.h"
#include "Model.h"
#include "Forcings.h"

const int MAX_HISTOGRAM_BINS=40; ///< Maximum allowable number of histogram bins

class CModel; //required for compilation

//////////////////////////////////////////////////////////////////
/// \brief Methods of aggregating geographical coverage
//
enum spatial_agg
{
  BY_HRU,        ///< Aggregate by HRU
  BY_BASIN,      ///< Aggregate by drainage basin
  BY_WSHED,      ///< Aggregate by watershed
  BY_HRU_GROUP,  ///< Aggregate by HRU group
  BY_SB_GROUP,   ///< Aggregate by SubBasin group
  BY_SELECT_HRUS ///< Aggregate by HRU, but only include one HRU Group
};

///////////////////////////////////////////////////////////////////
/// \brief Methods of aggregating time
//
enum time_agg
{
  DAILY,       ///< Aggregate by day
  MONTHLY,     ///< Aggregate by month
  YEARLY,      ///< Aggregate by year
  WATER_YEARLY,///< Aggregate by water year
  EVERY_NDAYS, ///< Aggregate by N-day intervals
  EVERY_TSTEP  ///< Aggregate by time-step
  //ENTIRE_SIM   ///< Aggregate over entire simulation (not currently implemented)
};

///////////////////////////////////////////////////////////////////////
/// \brief
/// \docminor This enumerated list and its values need to be described
enum diagnostic
{
  VAR_STATE_VAR,        ///< track state variable
  VAR_FORCING_FUNCTION, ///< track forcing function
  VAR_HYD_COND,         ///< track hydraulic conductivity
  VAR_FROM_FLUX,        ///< track gross flux from specific state variable
  VAR_TO_FLUX,          ///< track gross flux to specific state variable
  VAR_BETWEEN_FLUX      ///< track net flux between specific state variables
};

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for custom model output generator
class CCustomOutput
{
private:/*------------------------------------------------------*/

  ofstream     _CUSTOM;     ///< output file stream

  int          _netcdf_ID;  ///< netCDF file identifier

  diagnostic   _var;        ///< output variable identifier
  sv_type      _svtype;     ///< state variable output type (if output var is a SV)
  int          _svind;      ///< state variable index (if output var is a SV or flux)
  int          _svind2;     ///< target state variable index (if output var is a flux between two compartments)
  string       _force_str;  ///< forcing function name (if output var is a forcing function)
  forcing_type _ftype;      ///< forcing function type (if output var is a forcing function)

  agg_stat     _aggstat;    ///< time aggregation statistic(average, max, min, etc.) (spatial average is always used)
  time_agg     _timeAgg;    ///< how aggregated (monthly, daily, hourly, etc.)
  spatial_agg  _spaceAgg;   ///< how aggregated (by HRU, by Basin, etc.)

  double       _hist_min;   ///< histogram min
  double       _hist_max;   ///< Histogram max
  int          _nBins;      ///< Histogram # of bins

  string       _filename;   ///< custom output filename (relative path, with extension)
  string       _filename_user; ///< custom output filename provided by user

  string       _varName;    ///< forcing variable or state variable name
  string       _varUnits;   ///< forcing variable or state variable units
  string       _timeAggStr; ///< temporal aggregation type string
  string       _statStr;    ///< statistic type string
  string       _spaceAggStr;///< spatial aggregation type string

  double     **data;        ///< stores accumulated data for each HRU,Basin, or WShed (size:[num_store][num_data])
  int          num_data;    ///< number of data points
  int          num_store;   ///< number of data items needed for each HRU, Basin or WShed
                            //(e.g., =2 if max and min are both tracked)
  int         _time_index;  ///< index tracking current output line (e.g., 3=3 years/months/days passed, dependent upon _timeAgg

  int          count;       ///< counts accumulated data (# of timesteps since last output dump)

  int          kk_only;     ///< index of HRUGroup for which output is generated when spaceAgg==BY_SELECT_HRUS

  const CModel *pModel;     ///< Reference to model

  void DetermineCustomFilename(const optStruct& Options); 

public:/*------------------------------------------------------*/

  CCustomOutput(const diagnostic    variable,
                const sv_type       sv,
                const int           sv_index,
                const int           sv_index2,
                const string        force_string,
                const agg_stat      stat,
                const time_agg      time_aggregation,
                const spatial_agg   space_aggregation,
                const string        filename_spec,
                const int           kk,
                const CModel       *pMod,
                const optStruct                 &Options);
  ~CCustomOutput();

  void             CloseFiles(const optStruct &Options);

  void     SetHistogramParams(const double min,const double max, const int numBins);

  void InitializeCustomOutput(const optStruct &Options);

  void      WriteFileHeader  (const optStruct &Options);
  void      WriteCustomOutput(const time_struct &tt, const optStruct &Options);

  static CCustomOutput *ParseCustomOutputCommand(char *s[MAXINPUTITEMS],const int Len,CModel *&pModel,const optStruct &Options);
  
private:
  void     WriteCSVFileHeader();
  void   WriteEnSimFileHeader(const optStruct &Options);
  void  WriteNetCDFFileHeader(const optStruct &Options);
};
#endif
