/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2023 the Raven Development Team
  ----------------------------------------------------------------*/
#ifndef FORCINGGRID_H
#define FORCINGGRID_H

#include "RavenInclude.h"
#include "ParseLib.h"
#include "Forcings.h"
#include "Model.h"

#ifdef _RVNETCDF_
#include <netcdf.h>
#endif

///////////////////////////////////////////////////////////////////
/// \brief   Data abstraction for gridded, 3D forcings
/// \details Data Abstraction for gridded, 3D forcing data.
///          Dimensions are (x,y,t). Only grid cells specified in
///          GridWeights will be used. Grid cells are consecutively
///          numbered line by line, i.e.\n
///          [ [    1,   2,   3, ...,   NC, ], \n
///            [ NC+1,NC+2,NC+3, ..., 2 NC, ], \n
///            [                 ...      ] ]
/// \note    The series axis should be ordered, i.e., t[i+1]>t[i] always
///          Time units should be in days
/// \remark Used for forcing functions: precip, temp, pressure, etc.
//
class CForcingGrid //: public CForcingGridABC
{

private:/*------------------------------------------------------*/
  forcing_type _ForcingType;                 ///< Forcing type, e.g. PRECIP or TEMP
  string       _filename;                    ///< Name of NetCDF file
  string       _varname;                     ///< Name of forcing variable in NetCDF file

  string       _DimNames[3];                 ///< Names of all three dimensions as in NetCDF file
  ///                                        ///< [ dim_cols, dim_row, dim_time]
  int          _GridDims[3];                 ///< Length of dimensions of grid [ size = (x,y,t) = (NC, NR, _ChunkSize) ]

  int          _nHydroUnits;                 ///< number of HRUs (important for weights)

  int          _dim_order;                   ///< code (1-6) for different dimension orders  
  //                                         ///< (x,y,t) = 1, (y,x,t) = 2, (x,t,y) = 3,
  //                                         ///< (t,x,y) = 4, (y,t,x) = 5, (t,y,x) = 6

  int          _ChunkSize;                   ///< number of time points read before upper limit of
  ///                                        ///< allowed storage is reached
  int          _nChunk;                      ///< number of chunks (blocks) which can be read
  int          _iChunk;                      ///< current chunk read and stored in _aVal

  double       _start_day;                   ///< Day corresponding to local TS time 0.0 (beginning of time series)
  int          _start_year;                  ///< Year corresponding to local TS time 0.0 (beginning of time series)
  string       _tag;                         ///< data tag (additional information for data)
  double       _interval;                    ///< uniform interval between data points (in days); delta t
  int          _steps_per_day;               ///< number of data intervals per day (pre-calculated for speed) =1.0/_interval
  bool         _is_derived;                  ///< true if forcing grid is derived from input forcings (e.g. t_ave from t_min and t_max)
  ///                                        ///< false if forcing grid is truely read from NetCDF file (e.g. t_min or t_max)
  
  double     **_aVal;                        ///< Array of magnitudes of pulses (variable units)
  ///                                        ///< [size _ChunkSize, _nNonZeroWeightedGridCells]
  ///                                        ///< time steps are in model resolution (means original input data are
  ///                                        ///< already aggregated to match model resolution)

  double     **_GridWeight;                  ///< Sparse array of weights for each HRU for a list of cells    
  //                                         ///< Dimensions : [_nHydroUnits][_nWeights[k]] (variable)
  //                                         ///< _GridWeight[k][i] is fraction of forcing for HRU k is from grid cell _GridWtCellIDs[k][i]
  //                                         ///< where grid cell index l=_GridWtCellIDs[k][i] is derived by l = j * dim_cols + i
  ///                                        ///< where i and j are the zero-indexed column and row of cell l respectively and
  ///                                        ///< dim_cols is the total number of columns.
  ///                                        ///< Following contraint must be satisfied:
  ///                                        ///<      sum(_GridWeight[k][i], {i=0,_nWeights[k]-1}) = 1.0 for all HRUs k=0,...,_nHydroUnits-1
  int        **_GridWtCellIDs;               ///< cell IDs for all non-zero grid weights for HRU k size=[_nHydroUnits][_nWeights[k]] (variable)
  int         *_CellIDToIdx;                 ///< local cell index ic (ranging from 0 to _nNonZeroWeightedGridCells-1)) corresponding to cell ID [size: ncells]
  int         *_nWeights;                    ///< number of weights for each HRU k (size=_nHydroUnits) (each entry greater or equal to 1)
  int          _nNonZeroWeightedGridCells;   ///< Number of non-zero weighted grid cells:
  //                                         ///< This is effectively the number of data points which is stored from the original data. 
  int         *_IdxNonZeroGridCells;         ///< indexes of non-zero weighted grid cells [size = _nNonZeroWeightedGridCells]
  int          _WinLength[3];                ///< length of data grid window in each dimension (x,y,t - defaults to _GridDims)
  int          _WinStart [3];                ///< data grid window starting point (x,y,t - defaults to 0, 0, chunksize)
  
  int          _nPulses;                     ///< number of pulses (total duration=(nPulses-1)*_interval)
  bool         _pulse;                       ///< flag determining whether this is a pulse-based or
  ///                                        ///< piecewise-linear time series
  double       _t_corr;                      ///< correction time _t_corr, i.e. distance between
  ///                                        ///< current chunk start day and model start day (in days)
  bool         _deaccumulate;                ///< true if input precipitation needs to be deaccumulated from cum. mm to mm/d
  double       _TimeShift;                   ///< time shift of data (fractional day by which read data should be shifted)
  double       _LinTrans_a;                  ///< linear transformation of read data: new = a*data + b
  double       _LinTrans_b;                  ///< linear transformation of read data: new = a*data + b
  bool         _period_ending;               ///< true if data is period ending - subtracts additional interval *on top of _TimeShift*
  bool         _is_3D;                       ///< true if forcings are 3D (lat, lon, time); false if 2D (stations, time)
  double       _rainfall_corr;               ///< correction factor for rainfall (stored with gauge, used elsewhere)
  double       _snowfall_corr;               ///< correction factor for snowfall (stored with gauge, used elsewhere)
  double       _cloud_min_temp;              ///< minimum temperature threshold used to determine cloud_cover factor
  double       _cloud_max_temp;              ///< maximum temparature threshold used to determine cloud_cover factor
  double       _aAveTemp[12];                ///< representative average monthly temperatures [C]
  double       _aMinTemp[12];                ///< representative minimum monthly temperatures [C]
  double       _aMaxTemp[12];                ///< representative maximum monthly temperatures [C]
  double       _aAvePET [12];                ///< representative average monthly PET [mm/d] (or monthly PET factor [mm/d/K], if MONTHLY_FACTOR is used)

  string       _AttVarNames[4];              ///< variable attribute names for [0]: lat [1]: long [2]:elev [3]: station identifier
  double*      _aLatitude;                   ///< fixed array of cell centroid latitudes, if provided (size: _IdxNonZeroGridCells)
  double*      _aLongitude;                  ///< fixed array of cell centroid longitudes, if provided (size: _IdxNonZeroGridCells)
  double*      _aElevation;                  ///< fixed array of cell representative elevations, if provided (size: _IdxNonZeroGridCells)
  string*      _aStationIDs;                 ///< fixed array of cell station/cell IDS (size:_IdxNonZeroGridCells)

  void   CellIdxToRowCol(const int        cellid,
                         int              &row,
                         int              &column) const;             ///< returns row and column index of cell ID

  void   ReadAttGridFromNetCDF (const int ncid,const string varname,const int ncols,const int nrows,double *&values);
 // void   ReadAttGridFromNetCDF2(const int ncid,const string varname,const int nrows,const int ncols,string *values);

public:/*------------------------------------------------------*/
  //Constructors:

  // simple constructor
  CForcingGrid(
    string       ForcingType,
    string       filename,
    string       varname,
    string       DimNames[3],
    bool         is_3D
    );

  // copy constructor
  CForcingGrid( const CForcingGrid &grid );

  ~CForcingGrid();

  // Parses all information from NetCDF file and sets variables like grid dimensions and buffer size
  void ForcingGridInit( const optStruct   &Options );

  // Initialize sets the correction time _t_corr
  // (= distance between time series start day and model start day) and
  // QA/QC to check for example that modeling period is completely covered by forcing data
  void Initialize(const optStruct &Options);

  // Reallocate all arrays in class to (potentially updated) grid dimensions
  // mainly used when sub-daily grids have to be added to model
  void ReallocateArraysInForcingGrid( );

  // ReadData populates _aVal or does nothing if no new chunk need to be read (= current modeling time step is within current chunk)
  bool   ReadData(const optStruct   &Options,
                  const double global_model_time);

  // accessors
  double GetValue                   (const int ic, const int it) const;
  double GetValue_avg               (const int ic, const double &t, const int n) const;
  double GetValue_min               (const int ic, const double &t, const int n) const;
  double GetValue_max               (const int ic, const double &t, const int n) const;

  // Weighting matrix associated routines
  void   AllocateWeightArray(              const int        nHydroUnits,
                                           const int        nGridCells);                ///< allocates _GridWeight [nHydroUnits, nGridCells]
  void   SetWeightVal(                     const int        HRUID,
                                           const int        CellID,
                                           const double     weight);                    ///< sets one entry of _GridWeight[HRUID, CellID] = weight
  bool   CheckWeightArray(                 const int        nHydroUnits,
                                           const int        nGridCells,
                                           const CModel    *pModel);                    ///< checks if sum(_GridWeight[HRUID, :]) = 1.0 for all HRUIDs
  double GetGridWeight(                    const int        k,
                                           const int        CellID) const;              ///< returns weighting of HRU and CellID pair \todo[clean]: function not used
  double GetChunkIndexFromModelTimeStep(   const optStruct &Options,
                                           const double     global_model_time)  const;  ///< returns index in current chunk corresponding to model time step \todo[clean]: function not used
  double GetChunkIndexFromModelTimeStepDay(const optStruct &Options,
                                           const double     global_model_time)  const;  ///< returns index in current chunk corresponding to beginning of day of currentmodel time step \todo[clean]: function not used

  // Routines for checking content
  void   CheckValue3D(                     const double value,
                                           const double checkval,
                                           const int dim1_idx,
                                           const int dim2_idx,
                                           const int dim3_idx);                         ///< checks for missing and fill values in 3D array and exists Raven if found
  void   CheckValue2D(                     const double value,
                                           const double checkval,
                                           const int dim1_idx,
                                           const int dim2_idx);                         ///< checks for missing and fill values in 2D array and exists Raven if found

  // set class variables
  void         SetForcingType(             const forcing_type &ForcingType);               ///< set _ForcingType               of class
  void         SetFilename(                const string filename);                  ///< set _filename                  of class
  void         SetVarname(                 const string varname);                   ///< set _varname                   of class
  void         SetDimNames(                const string DimNames[3]);               ///< set _DimNames                  of class
  void         SetGridDims(                const int    GridDims[3]);               ///< set _GridDims                  of class
  void         SetToDeaccumulate();                                                 ///< set _deaccumulate              of class
  void         SetLinearTransform(         const double LinTrans_a,                 ///< set _LinTrans_a and _b         of class
                                           const double LinTrans_b);
  void         SetTimeShift(               const double TimeShift);                 ///< set _TimeShift                 of class
  void         SetAsPeriodEnding           ();                                      ///< set _period_ending
  void         SetIs3D(                    const bool   is3D);                      ///< set _is3D                      of class
  void         SetIdxNonZeroGridCells(     const int    nHydroUnits,
                                           const int    nGridCells, const optStruct &Options); 
                                                                                    ///< set _IdxNonZeroGridCells       of class
  void         CalculateChunkSize         (const optStruct &Options);
                                                                                    
  void         SetnHydroUnits(             const int    nHydroUnits);               ///< set _nHydroUnits               of class
  void         SetChunkSize(               const int    ChunkSize);                 ///< set _ChunkSize                 of class
  void         SetInterval(                const double interval);                  ///< set _interval                  of class
  void         SetSnowfallCorr(            const double snowfall_corr);             ///< set _snowfall_corr             of class
  void         SetRainfallCorr(            const double rainfall_corr);             ///< set _rainfall_corr             of class
  void         Setcloud_min_temp(          const double cloud_min_temp);            ///< set _cloud_min_temp            of class
  void         Setcloud_max_temp(          const double cloud_max_temp);            ///< set _cloud_max_temp            of class
  void         SetaAveTemp(                const double aAveTemp[12]);              ///< set _aAveTemp[12]              of class
  void         SetaMinTemp(                const double aMinTemp[12]);              ///< set _aMinTemp[12]              of class
  void         SetaMaxTemp(                const double aMaxTemp[12]);              ///< set _aMaxTemp[12]              of class
  void         SetaAvePET(                 const double aAvePET [12]);              ///< set _aAvePET [12]              of class
  void         SetValue(                   const int idx,
                                           const int t,
                                           const double aVal);                      ///< set _aVal              of class
  void         SetAttributeVarName(        const string var, const string varname); ///< set elevation, lat, or long var name
  void         SetStationElevation(        const int idx, const double &elev);      ///< set elevation of station idx
  void         SetIsDerived               (const bool is_derived);

  // get class variables
  double       GetInterval()                                      const; ///< data interval (in days)
  bool         GetIsDerived()                                     const; ///< if data are read from NetCDF (false) or derived from these data (true)
  bool         GetIs3D()                                          const; ///< true if NetCDF data are (lat,lon,time), false if data are (nstations,time)
  int          GetStartYear()                                     const; ///< start year of gridded time series data
  double       GetStartDay()                                      const; ///< start day of time series data
  int          GetCols()                                          const; ///< Number of columns (= 1st dimension of gridded data)
  int          GetRows()                                          const; ///< Number of rows    (= 2nd dimension of gridded data)
  int          GetNumValues()                                     const; ///< Number of pulses  (= 3rd dimension of gridded data)
  int          GetNumberNonZeroGridCells()                        const; ///< Number of non-zero weighted grid cells 
  int          GetChunkSize()                                     const; ///< Current chunk size
  int          GetnHydroUnits()                                   const; ///< get number of HRUs _nHydroUnits
  forcing_type GetForcingType()                                   const; ///< Type of forcing data, e.g. PRECIP, TEMP
  bool         ShouldDeaccumulate()                               const; ///< true if data must be deaccumulated
  double       GetSnowfallCorr()                                  const; ///< snowfall correction factor
  double       GetRainfallCorr()                                  const; ///< rainfall correction factor
  double       GetCloudMinRange()                                 const; ///< Minimum temperature threshold used to determine cloud_cover factor
  double       GetCloudMaxRange()                                 const; ///< Maximum temperature threshold used to determine cloud_cover factor
  double       GetMonthlyAveTemp(const int month)                 const; ///< Representative Average temperature for month
  double       GetMonthlyMinTemp(const int month)                 const; ///< Representative minimum temperature for month
  double       GetMonthlyMaxTemp(const int month)                 const; ///< Representative maximum temperature for month
  double       GetMonthlyAvePET (const int month)                 const; ///< Average PET over month
  double       DailyTempCorrection(const double t)                const; ///< Daily temperature correction [C]
  int          GetTimeIndex(const double &t, const double &tstep) const; ///< get time index corresponding to t+tstep/2
  string       GetFilename()                                      const; ///< return forcing filename

  double       GetCellLatitude       (const int l) const;        ///< returns Latitude of cell l (or 0, if not available)
  double       GetCellLongitude      (const int l) const;        ///< returns Longitude of cell l (or 0, if not available)
  double       GetRefElevation       (const int k) const;        ///< returns representative elevation of forcing data in HRU k (or 0, if not available)

  double       GetWeightedValue          (const int k, const double &t,const double &tstep) const; ///<returns weighted value in HRU k
  double       GetDailyWeightedValue     (const int k, const double &t,const double &tstep, const optStruct &Options) const; ///<returns daily weighted value in HRU k
  double       GetWeightedAverageSnowFrac(const int k, const double &t,const double &tstep,const CForcingGrid *pRain) const; ///<returns daily weighted value in HRU k
};

#endif


