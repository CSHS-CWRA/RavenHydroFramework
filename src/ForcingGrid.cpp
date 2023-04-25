/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2023 the Raven Development Team
  ----------------------------------------------------------------*/

#include "ForcingGrid.h"
#include "ParseLib.h"  // for GetFilename()
#include "Forcings.h"
#include <string.h>

/*****************************************************************
   Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

///////////////////////////////////////////////////////////////////
/// \brief Implementation of forcing grid constructor \n
///        Sets only the main four arguments and calls then the initialization
///        routine determining grid dimensions, buffersize etc.
///
/// \param ForcingType   [in] Forcing type, e.g. PRECIP or TEMP
/// \param filename      [in] Name of NetCDF file
/// \param varname       [in] Name of variable in NetCDF file
/// \param DimNames[3]   [in] Names of all three dimensions as in NetCDF file [ dim_cols, dim_row, dim_time]
/// \param is_3D         [in] true if NetCDF grid is 3D (gridded) data rather than 2D (station) data
//
CForcingGrid::CForcingGrid(string       ForcingType,
                           string       filename,
                           string       varname,
                           string       DimNames[3],
                           bool         is_3D
                          )
{
  _ForcingType = GetForcingTypeFromString(ForcingType);
  _filename		 = filename;
  _varname		 = varname;
  _DimNames[0] = DimNames[0]; _DimNames[1]  = DimNames[1]; _DimNames[2]  = DimNames[2];
  _is_3D		   = is_3D;
  _GridDims [0]=0; _GridDims [1]=0; _GridDims [2]=0;
  _WinLength[0]=0; _WinLength[1]=0; _WinLength[2]=0;
  _WinStart [0]=0; _WinStart [1]=0; _WinStart [2]=0;

  _interval		    = 1.0;
  _steps_per_day	= 1;
  
  _TimeShift		  = 0.0;
  _LinTrans_a		  = 1.0;
  _LinTrans_b		  = 0.0;
  _deaccumulate		= false;
  _period_ending	= false;

  // -------------------------------
  // Additional variables initialized and eventually overwritten by ParseTimeSeries
  // -------------------------------
  _rainfall_corr  = 1.0;
  _snowfall_corr  = 1.0;

  _cloud_min_temp =-20.0; //ensures cloud-free status always unless overriden
  _cloud_max_temp = 40.0;

  for (int i=0;i<12;i++)
  {
    _aAveTemp[i]=NOT_SPECIFIED;
    _aMinTemp[i]=NOT_SPECIFIED;
    _aMaxTemp[i]=NOT_SPECIFIED;
    _aAvePET [i]=NOT_SPECIFIED;
  }

  //initialized in ReallocateArraysInForcingGrid
  _aVal                = NULL;

  // initialized in AllocateWeightArray,SetIdxNonZeroGridCells()
  _GridWeight          = NULL;
  _GridWtCellIDs       = NULL;
  _CellIDToIdx         = NULL;
  _nWeights            = NULL;
  _IdxNonZeroGridCells = NULL;
  _nNonZeroWeightedGridCells=0;

  //initialized in CalculateChunkSize()
  _ChunkSize           =0;
  _nChunk              =1;
  _iChunk              = -1; // current chunk read (-1 = no chunk read, 0 = first chunk...)

  //initialized in SetAttributeVarName
  _aLatitude           = NULL;
  _aLongitude          = NULL;
  _aElevation          = NULL;
  _aStationIDs         = NULL;
  _AttVarNames[0]      ="NONE";
  _AttVarNames[1]      ="NONE";
  _AttVarNames[2]      ="NONE";
  _AttVarNames[3]      ="NONE";
  
}
///////////////////////////////////////////////////////////////////
/// \brief Copy constructor.
///
/// \param ForcingType   [in] an existing grid
CForcingGrid::CForcingGrid( const CForcingGrid &grid )
{
  _ForcingType                 = grid._ForcingType                     ;
  _filename                    = grid._filename                        ;
  _varname                     = grid._varname                         ;
  _deaccumulate                = grid._deaccumulate                    ;
  _is_3D                       = grid._is_3D                           ;
  _TimeShift                   = grid._TimeShift                       ;
  _LinTrans_a                  = grid._LinTrans_a                      ;
  _LinTrans_b                  = grid._LinTrans_b                      ;
  _period_ending               = grid._period_ending                   ;
  for (int ii=0; ii<3;  ii++) {_DimNames [ii]= grid._DimNames [ii]; }
  for (int ii=0; ii<3;  ii++) {_GridDims [ii]= grid._GridDims [ii]; }
  for (int ii=0; ii<3;  ii++) {_WinLength[ii]= grid._WinLength[ii]; }
  for (int ii=0; ii<3;  ii++) {_WinStart [ii]= grid._WinStart [ii]; }
  _nNonZeroWeightedGridCells   = grid._nNonZeroWeightedGridCells;
  _nHydroUnits                 = grid._nHydroUnits                     ;
  _ChunkSize                   = grid._ChunkSize                       ;
  _nChunk                      = grid._nChunk                          ;
  _iChunk                      = grid._iChunk                          ;
  _start_day                   = grid._start_day                       ;
  _start_year                  = grid._start_year                      ;
  _tag                         = grid._tag                             ;
  _interval                    = grid._interval                        ;
  _is_derived                  = true                                  ;
  _nPulses                     = grid._nPulses                         ;
  _pulse                       = grid._pulse                           ;
  _t_corr                      = grid._t_corr                          ;
  _rainfall_corr               = grid._rainfall_corr                   ;
  _snowfall_corr               = grid._snowfall_corr                   ;
  _cloud_min_temp              = grid._cloud_min_temp                  ;
  _cloud_max_temp              = grid._cloud_max_temp                  ;
  for (int ii=0; ii<12; ii++) {_aAveTemp[ii] = grid._aAveTemp[ii];}
  for (int ii=0; ii<12; ii++) {_aMinTemp[ii] = grid._aMinTemp[ii];}
  for (int ii=0; ii<12; ii++) {_aMaxTemp[ii] = grid._aMaxTemp[ii];}
  for (int ii=0; ii<12; ii++) {_aAvePET [ii] = grid._aAvePET [ii];}
  
  _aVal=NULL;
  _aVal = new double *[_ChunkSize];
  ExitGracefullyIf(_aVal==NULL,"CForcingGrid::Copy Constructor(1)",OUT_OF_MEMORY);
  for (int it=0; it<_ChunkSize; it++) {                       // loop over time points in buffer
    _aVal[it]=NULL;
    _aVal[it] = new double [_nNonZeroWeightedGridCells];
    ExitGracefullyIf(_aVal[it]==NULL,"CForcingGrid::Constructor",OUT_OF_MEMORY);
    for (int ic=0; ic<_nNonZeroWeightedGridCells;ic++){       // loop over non-zero weighted cells
      _aVal[it][ic]=grid._aVal[it][ic];                            // copy the value
    }
  }

  int ncells;
  if (_is_3D) { ncells = _GridDims[0] * _GridDims[1];}
  else        { ncells = _GridDims[0];}

 
  //cout<<"Creating new GridWeights array (Copy Constructor): "<<ForcingToString(_ForcingType)<<endl;

  AllocateWeightArray(_nHydroUnits,ncells);
  for (int k=0; k<_nHydroUnits; k++) {                       
    for(int i=0;i<grid._nWeights[k];i++) {
      SetWeightVal(k,grid._GridWtCellIDs[k][i],grid._GridWeight[k][i]); //dynamically grows sparse matrix
    }
  }
  
  _CellIDToIdx =NULL;
  _CellIDToIdx = new int [ncells];
  ExitGracefullyIf(_CellIDToIdx==NULL,"CForcingGrid::Copy Constructor(8)",OUT_OF_MEMORY);
  for(int c=0; c<ncells; c++) {             // loop over non-zero weighted cells
    _CellIDToIdx[c]=grid._CellIDToIdx[c];              // copy the value
  }

  _IdxNonZeroGridCells = NULL;
  _IdxNonZeroGridCells = new int [_nNonZeroWeightedGridCells];
  ExitGracefullyIf(_IdxNonZeroGridCells==NULL,"CForcingGrid::Copy Constructor(3)",OUT_OF_MEMORY);
  for (int ic=0; ic<_nNonZeroWeightedGridCells; ic++) {             // loop over non-zero weighted cells
    _IdxNonZeroGridCells[ic]=grid._IdxNonZeroGridCells[ic];              // copy the value
  }

  _aLatitude=NULL;_aLongitude=NULL;_aElevation=NULL;_aStationIDs=NULL;
  if(grid._aLatitude!=NULL) {
    _aLatitude=new double [_nNonZeroWeightedGridCells];
    ExitGracefullyIf(_aLatitude==NULL,"CForcingGrid::Copy Constructor(4)",OUT_OF_MEMORY);
    for(int c=0; c<_nNonZeroWeightedGridCells; c++) {_aLatitude[c]=grid._aLatitude[c];}
  }
  if(grid._aLongitude!=NULL) {
    _aLongitude=new double[_nNonZeroWeightedGridCells];
    ExitGracefullyIf(_aLongitude==NULL,"CForcingGrid::Copy Constructor(5)",OUT_OF_MEMORY);
    for(int c=0; c<_nNonZeroWeightedGridCells; c++) { _aLongitude[c]=grid._aLongitude[c]; }
  }
  if(grid._aElevation!=NULL) {
    _aElevation=new double[_nNonZeroWeightedGridCells];
    ExitGracefullyIf(_aElevation==NULL,"CForcingGrid::Copy Constructor(6)",OUT_OF_MEMORY);
    for(int c=0; c<_nNonZeroWeightedGridCells; c++) { _aElevation[c]=grid._aElevation[c]; }
  }
  if(grid._aStationIDs!=NULL) {
    _aStationIDs=new string[_nNonZeroWeightedGridCells];
    ExitGracefullyIf(_aStationIDs==NULL,"CForcingGrid::Copy Constructor(7)",OUT_OF_MEMORY);
    for(int c=0; c<_nNonZeroWeightedGridCells; c++) { _aStationIDs[c]=grid._aStationIDs[c]; }
  }
  
}

///////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
//
CForcingGrid::~CForcingGrid()
{
  if (DESTRUCTOR_DEBUG){cout<<"    DELETING GRIDDED DATA"<<endl;}
  if(_aVal!=NULL) {
    for(int it=0; it<_ChunkSize; it++) { delete[] _aVal[it];      _aVal[it]=NULL; }      delete[] _aVal;_aVal= NULL;
  }
  
  for(int k=0; k<_nHydroUnits; k++) {
    delete[] _GridWeight[k];    _GridWeight   [k]=NULL;
    delete[] _GridWtCellIDs[k]; _GridWtCellIDs[k]=NULL;
  }
  delete [] _GridWeight;            _GridWeight          = NULL;
  delete [] _GridWtCellIDs;         _GridWtCellIDs       = NULL;
  delete [] _CellIDToIdx;           _CellIDToIdx         = NULL;
  delete [] _nWeights;              _nWeights            = NULL;
  delete [] _IdxNonZeroGridCells;   _IdxNonZeroGridCells = NULL;
  delete [] _aLatitude;             _aLatitude           = NULL;
  delete [] _aLongitude;            _aLongitude          = NULL;
  delete [] _aElevation;            _aElevation          = NULL;
  delete [] _aStationIDs;           _aStationIDs         = NULL;
}


///////////////////////////////////////////////////////////////////
/// \brief Implementation of forcing grid constructor only determining grid
///        dimensions and buffersize (data are only initalized but not read from file)
///
/// \note  Needs _ForcingType, _filename, _varname, _DimNames to be set already.
///        Use for that either "SetForcingType", "SetFilename", "SetVarname", and "SetDimNames" \n
///        or "CForcingGrid()".
void CForcingGrid::ForcingGridInit(const optStruct   &Options)
{
#ifdef _RVNETCDF_
  int    ncid;                  // file unit
  int    dimid_x;               // id of x dimension (columns)
  int    dimid_y(0);            // id of y dimension (rows)
  int    dimid_t;               // id of t dimension (time)
  int    dimids_var[3];         // ids of dimensions of a NetCDF variable
  int    varid_t;               // id of time variable
  int    varid_f;               // id of forcing variable read
  int    retval;                // error value for NetCDF routines
  size_t GridDim_t;             // special type for GridDims required by nc routine
  char * unit_t;                // special type for string of variable's unit     required by nc routine
  char * unit_f;                // special type for string of variable's unit     required by nc routine
  int    calendar;              // enum int of calendar used
  char * long_name_f;           // special type for string of variable's long name required by nc routine

  int    iatt;
  char   attrib_name[500];

  int    retval1;
  size_t att_len;        // length of the attribute's text

  string dash;           // to check format of time unit string
  string colon;          // to check format of time unit string
  string unit_t_str;     // to check format of time unit string
  int    ntime;          // number of time steps
  
  if(_ForcingType==F_UNRECOGNIZED) {
    ExitGracefully("ParseTimeSeriesFile: :GriddedForcing and :StationForcing blocks requires valid :ForcingType command",BAD_DATA);
  }
  if(_filename=="NONE") {
    ExitGracefully("ParseTimeSeriesFile: :GriddedForcing and :StationForcing blocks requires valid :FileNameNC command",BAD_DATA);
  }
  if(_varname=="NONE") {
    ExitGracefully("ParseTimeSeriesFile: :GriddedForcing and :StationForcing blocks requires valid :VarNameNC command",BAD_DATA);
  }
  if(_DimNames[0]=="NONE") {
    ExitGracefully("ParseTimeSeriesFile: :GriddedForcing and :StationForcing blocks requires valid :DimNamesNC command",BAD_DATA);
  }
  if(_AttVarNames[2]=="NONE") {
    if (((_ForcingType==F_TEMP_DAILY_MIN) && (Options.orocorr_temp  ==OROCORR_NONE)) ||
        ((_ForcingType==F_TEMP_DAILY_MAX) && (Options.orocorr_temp  ==OROCORR_NONE)) ||
        ((_ForcingType==F_TEMP_AVE      ) && (Options.orocorr_temp  ==OROCORR_NONE)) ||
        ((_ForcingType==F_PRECIP        ) && (Options.orocorr_precip==OROCORR_NONE)) ||
        ((_ForcingType==F_RAINFALL      ) && (Options.orocorr_precip==OROCORR_NONE)) ||
        ((_ForcingType==F_PET           ) && (Options.orocorr_PET   ==OROCORR_NONE)))
    {
      //warning not needed
    }
    else {
      string warning="Since no elevation data are in NetCDF file "+_filename+", all orographic corrections for "+ForcingToString(_ForcingType)+" will be disabled unless :StationElevations or :StationElevationsByIdx command is present. (or :ElevationVarNameNC is after grid weights command in .rvt file)";
      WriteAdvisory(warning,Options.noisy);
    }
  }

  // open NetCDF read-only; ncid will be set
  if(Options.noisy) { cout<<"Initializing grid file "<<_filename<<endl; }

  string filename_e=_filename;
  SubstringReplace(filename_e,"*",to_string(g_current_e+1)); //replaces wildcard for ensemble runs

  retval = nc_open(filename_e.c_str(),NC_NOWRITE,&ncid);          HandleNetCDFErrors(retval);


  // Get the id of dimensions based on its name; dimid will be set
  // Find length of dimension and store it in GridDim
  if(_is_3D) {
    // dimension x = number of columns of the grid
    retval = nc_inq_dimid(ncid,_DimNames[0].c_str(),&dimid_x);   HandleNetCDFErrors(retval);
    retval = nc_inq_dimlen(ncid,dimid_x,&GridDim_t);             HandleNetCDFErrors(retval);
    _GridDims[0] = static_cast<int>(GridDim_t);  // convert returned 'size_t' to 'int'

    // dimension y = number of rows of the grid
    retval = nc_inq_dimid(ncid,_DimNames[1].c_str(),&dimid_y);   HandleNetCDFErrors(retval);
    retval = nc_inq_dimlen(ncid,dimid_y,&GridDim_t);             HandleNetCDFErrors(retval);
    _GridDims[1] = static_cast<int>(GridDim_t);  // convert returned 'size_t' to 'int'

    // dimension t = number of time steps in NetCDF file
    retval = nc_inq_dimid(ncid,_DimNames[2].c_str(),&dimid_t);   HandleNetCDFErrors(retval);
    retval = nc_inq_dimlen(ncid,dimid_t,&GridDim_t);             HandleNetCDFErrors(retval); //this should also be able to handle UNLIMITED time dimension
    _GridDims[2] = static_cast<int>(GridDim_t);  // convert returned 'size_t' to 'int'
  }
  else {
    // dimension x = number of stations in NetCDF file
    retval = nc_inq_dimid(ncid,_DimNames[0].c_str(),&dimid_x);   HandleNetCDFErrors(retval);
    retval = nc_inq_dimlen(ncid,dimid_x,&GridDim_t);             HandleNetCDFErrors(retval);
    _GridDims[0] = static_cast<int>(GridDim_t);  // convert returned 'size_t' to 'int'

    // dimension t = number of time steps in NetCDF file
    retval = nc_inq_dimid(ncid,_DimNames[1].c_str(),&dimid_t);   HandleNetCDFErrors(retval);
    retval = nc_inq_dimlen(ncid,dimid_t,&GridDim_t);             HandleNetCDFErrors(retval);
    _GridDims[1] = static_cast<int>(GridDim_t);  // convert returned 'size_t' to 'int'
  }

  // Get the id of the data variable based on its name; varid will be set
  string varname_e=_varname;
  SubstringReplace(varname_e,"*",to_string(g_current_e+1));//replaces wildcard for ensemble runs

  retval = nc_inq_varid(ncid,varname_e.c_str(),&varid_f);         HandleNetCDFErrors(retval);

  // determine in which order the dimensions are in variable
  retval = nc_inq_vardimid(ncid,varid_f,dimids_var);             HandleNetCDFErrors(retval);

  if(_is_3D) {
    //   if  (t,         y=lat=row, x=lon=col)   --> (2,1,0) --> (dimid_t, dimid_y, dimid_x)
    //   if  (t,         x=lon=col, y=lat=row)   --> (2,0,1) --> (dimid_t, dimid_x, dimid_y)
    //   if  (y=lat=row, t,         x=lon=col)   --> (1,2,0) --> (dimid_y, dimid_t, dimid_x)
    //   if  (y=lat=row, x=lon=col, t        )   --> (1,0,2) --> (dimid_y, dimid_x, dimid_t)
    //   if  (x=lon=col, t,         y=lat=row)   --> (0,2,1) --> (dimid_x, dimid_t, dimid_y)
    //   if  (x=lon=col, y=lat=row, t        )   --> (0,1,2) --> (dimid_x, dimid_y, dimid_t)
    if     ((dimids_var[0] == dimid_x) && (dimids_var[1] == dimid_y) && (dimids_var[2] == dimid_t)) {  // dimensions are (x,y,t)
      _dim_order = 1;
    }
    else if((dimids_var[0] == dimid_y) && (dimids_var[1] == dimid_x) && (dimids_var[2] == dimid_t)) {  // dimensions are (y,x,t)
      _dim_order = 2;
    }
    else if((dimids_var[0] == dimid_x) && (dimids_var[1] == dimid_t) && (dimids_var[2] == dimid_y)) {  // dimensions are (x,t,y)
      _dim_order = 3;
    }
    else if((dimids_var[0] == dimid_t) && (dimids_var[1] == dimid_x) && (dimids_var[2] == dimid_y)) {  // dimensions are (t,x,y)
      _dim_order = 4;
    }
    else if((dimids_var[0] == dimid_y) && (dimids_var[1] == dimid_t) && (dimids_var[2] == dimid_x)) {  // dimensions are (y,t,x)
      _dim_order = 5;
    }
    else if((dimids_var[0] == dimid_t) && (dimids_var[1] == dimid_y) && (dimids_var[2] == dimid_x)) {  // dimensions are (t,y,x)
      _dim_order = 6;
    }
  }
  else {
    //   if  (station,t)   --> (1,0) --> (dimid_x, dimid_t)
    //   if  (t,station)   --> (0,1) --> (dimid_t, dimid_x)
    if((dimids_var[0] == dimid_x) && (dimids_var[1] == dimid_t)) {  // dimensions are (station,t)
      _dim_order = 1;
    }
    else if((dimids_var[0] == dimid_t) && (dimids_var[1] == dimid_x)) {  // dimensions are (t,station)
      _dim_order = 2;
    }
  }
  if(Options.noisy) { printf("  Order of dimensions in NetCDF is Case %i...\n",_dim_order); }

  // start year, date
  // ----------------------------------------------------------------------------------------------
  // Start day and year are extracted from the UNITS attribute of time axis
  // can be, e.g., "days    since 1989-01-01 00:00:00" or
  //               "hours   since 1989-01-01 00:00:00" or
  //               "minutes since 1989-01-01 00:00:00" or
  //               "seconds since 1989-01-01 00:00:00"
  // then whole time variable has to be read and first time step needs to be added to get _start_year and _start_day
  // difference between t[0] and t[1] is then used to get _interval
  // -------------------------------
  if(_is_3D) {retval = nc_inq_varid(ncid,_DimNames[2].c_str(),&varid_t);   HandleNetCDFErrors(retval);}
  else       {retval = nc_inq_varid(ncid,_DimNames[1].c_str(),&varid_t);   HandleNetCDFErrors(retval);}

  // unit of time
  //------------------------------------------------------------------------------------------------
  retval = nc_inq_attlen(ncid,varid_t,"units",&att_len);                   HandleNetCDFErrors(retval);
  unit_t=new char[att_len+1];
  retval = nc_get_att_text(ncid,varid_t,"units",unit_t);                   HandleNetCDFErrors(retval);// read attribute text
  unit_t[att_len] = '\0';// add string determining character
  unit_t_str=to_string(unit_t);
  if(!IsValidNetCDFTimeString(unit_t_str))   // check that unit of time is in format "[days/minutes/...] since YYYY-MM-DD HH:MM:SS{+0000}"
  {
    cout<<"time unit string: "<<unit_t_str<<endl;
    ExitGracefully("CForcingGrid::ForcingGridInit: time unit string is not in the format '[days/hours/...] since YYYY-MM-DD HH:MM:SS +0000' !",BAD_DATA);
  }

  // calendar attribute
  // ----------------------------------------------------------------------------------------------
  calendar=GetCalendarFromNetCDF(ncid,varid_t,_filename,Options);
			 
  // time attribute: set my_time[]
  // ----------------------------------------------------------------------------------------------
  if (_is_3D) {ntime = _GridDims[2];}
  else        {ntime = _GridDims[1];}

  double *my_time=new double[ntime];
  ExitGracefullyIf(my_time==NULL,"CForcingGrid::ForcingGridInit",OUT_OF_MEMORY);
  GetTimeVectorFromNetCDF(ncid,varid_t,ntime,my_time);

  double time_zone=0;
  GetTimeInfoFromNetCDF(unit_t,calendar,my_time,ntime,_filename,_interval,_start_day,_start_year,time_zone);
  _steps_per_day=(int)(rvn_round(1.0/_interval)); //pre-calculate for speed.
  delete[] unit_t;

  /*
  printf("ForcingGrid: unit_t:          %s\n",unit_t_str.c_str());
  printf("ForcingGrid: tt.julian_day:   %f\n",tt.julian_day);
  printf("ForcingGrid: tt.day_of_month: %i\n",tt.day_of_month);
  printf("ForcingGrid: tt.month:        %i\n",tt.month);
  printf("ForcingGrid: tt.year:         %i\n",tt.year);
  printf("ForcingGrid: my_time[0]:      %f\n",my_time[0]);
  printf("ForcingGrid: _interval:       %f\n",_interval);
  printf("ForcingGrid: _start_day:      %f\n",_start_day);
  printf("ForcingGrid: _start_year:     %d\n",_start_year); 
  printf("ForcingGrid: time_zone:       %f\n",time_zone);
  printf("ForcingGrid: time shift:       %f\n",_TimeShift);
  printf("ForcingGrid: # times:       %i\n",ntime);
  time_struct tp;
  JulianConvert(0.0,_start_day-_interval,_start_year,calendar,tp);
  printf("ForcingGrid: start string:    %s %s\n",tp.date_string.c_str(),DecDaysToHours(tp.julian_day).c_str());
  */

  delete[] my_time;

  //QA/QC:
  //--------------------------------
  if(_interval<=0) {
    ExitGracefully("CForcingGrid: ForcingGridInit: negative time interval is not allowed",BAD_DATA);
  }
  if((ForcingToString(_ForcingType) == "TEMP_DAILY_AVE") && (_interval != 1.0)) {
    ExitGracefully("CForcingGrid: ForcingGridInit: Gridded forcing 'TEMP_DAILY_AVE' must have daily timestep. Please use 'TEMP_AVE' instead (see *.rvt file).",BAD_DATA);
  }

  if(_period_ending) { _TimeShift-=_interval; }

  // -------------------------------
  // add time shift to data
  //      --> only applied when _interval < 1.0 (daily)
  //      --> otherwise ignored and warning written to RavenErrors.txt
  // -------------------------------             
  if (_interval >= 1.0) {   // data are not sub-daily
    if ( ceil(_TimeShift) == _TimeShift) {  // time shift of whole days requested
      AddTime(_start_day,_start_year,_TimeShift,calendar,_start_day,_start_year) ;
    }
    else {  // sub-daily shifts (e.g. 1.25) of daily data requested
      WriteAdvisory("CForcingGrid: ForcingGridInit: time shift specified for NetCDF time series will be ignored", Options.noisy);
      WriteAdvisory("                               because inputs are daily and time shift is sub-daily", Options.noisy);
    }
  }
  else {  // data are sub-daily
    AddTime(_start_day,_start_year,_TimeShift,calendar,_start_day,_start_year) ;
  }
   
  // -------------------------------
  // set _tag
  //     ---> looks for attributes "long_name" and "units" of forcing variable
  //     --> if either attribute not available it will be set to "?"
  // -------------------------------
  retval = 0;
  iatt   = 0;
  long_name_f = (char *) malloc(3);
  strcpy(long_name_f, "?\0");
  unit_f = (char *) malloc(3);
  strcpy(unit_f, "?\0");

  while (retval==0)
  {
    // inquire attributes name
    retval = nc_inq_attname(ncid, varid_f, iatt, attrib_name);
    if (strcmp(attrib_name,"long_name") == 0)// long_name of forcing
    {
      retval1 = nc_inq_attlen (ncid, varid_f, attrib_name, &att_len);     HandleNetCDFErrors(retval1);
      long_name_f = (char *) malloc(att_len + 1);// allocate memory of char * to hold attribute's text
      retval1 = nc_get_att_text(ncid, varid_f, attrib_name, long_name_f); HandleNetCDFErrors(retval1);
      long_name_f[att_len] = '\0';// add string determining character
    }
    else if (strcmp(attrib_name,"units") == 0)// unit of forcing
    {
      retval1 = nc_inq_attlen (ncid, varid_f, attrib_name, &att_len);     HandleNetCDFErrors(retval1);
      unit_f = (char *) malloc(att_len + 1);// allocate memory of char * to hold attribute's text
      retval1 = nc_get_att_text(ncid, varid_f, attrib_name, unit_f);      HandleNetCDFErrors(retval1);
      unit_f[att_len] = '\0';// add string determining character
    }
    iatt++;
  }
  _tag = to_string(long_name_f)+" in ["+to_string(unit_f)+"]";
  free( unit_f); free(long_name_f);
  if (Options.noisy){ printf("Forcing found in NetCDF file: %s \n",_tag.c_str()); }

  // -------------------------------
  // Set number of pulses and pulse type to be consistent with class CTimeSeries
  // -------------------------------
  _nPulses      = ntime;
  _pulse        = true;
  ExitGracefullyIf(_nPulses<=0,
                   "CForcingGrid: ForcingGridInit: no time point entries in forcing grid",BAD_DATA);

  // -------------------------------
  // Close the file. This frees up any internal NetCDF resources
  // associated with the file, and flushes any buffers.
  // -------------------------------
  retval = nc_close(ncid);         HandleNetCDFErrors(retval);

  _is_derived = false;

  // -------------------------------
  // Correction time initialized with zero but will be set correctly in Initialize()
  // -------------------------------
  _t_corr=0.0;

  // -------------------------------
  // Weighting array will be allocated and initialized with AllocateWeightArray()
  // -------------------------------
  if (Options.noisy){
    printf("ForcingGrid: Interval:        %f\n",_interval);
    printf("ForcingGrid: _start_day:      %f\n",_start_day);
    printf("ForcingGrid: _start_year:     %i\n",_start_year);
    printf("ForcingGrid: duration:        %f\n",(double)(_nPulses)*_interval);
    printf("\n");
  }
#endif   // ends #ifdef _RVNETCDF_
}

///////////////////////////////////////////////////////////////////
/// \brief  Reallocate all arrays in class to (potentially updated) time grid dimensions
//          mainly used when sub-daily grids have to be added to model using ForcingCopyCreate
//
void CForcingGrid::ReallocateArraysInForcingGrid( )
{
  int ntime ;  // number of time steps
  
  ntime  = _GridDims[2]; //assumes _Is_3D
  //if (_is_3D) {ntime = _GridDims[2];}
  //else        {ntime = _GridDims[1];}
  // 
  // -------------------------------
  // Initialize data array and set all entries to NODATA value
  // -------------------------------
  _aVal = NULL;
  _aVal =  new double *[ntime];
  for (int it=0; it<ntime; it++) {                       // loop over time points in buffer
    _aVal[it]=NULL;
    _aVal[it] = new double [_nNonZeroWeightedGridCells];
    ExitGracefullyIf(_aVal[it]==NULL,"CForcingGrid::ReallocateArraysInForcingGrid",OUT_OF_MEMORY);
    for (int ic=0; ic<_nNonZeroWeightedGridCells;ic++){       // loop over non-zero weighted cells
      _aVal[it][ic]=NETCDF_BLANK_VALUE;                       // initialize
    }
  }
}

///////////////////////////////////////////////////////////////////
/// \brief  Updates class variable _aVal containing current chunk of data \\
///         chunk = 0         --> data[          0*_chunksize : 1*_chunksize-1][:][:] \\
///         chunk = 1         --> data[          1*_chunksize : 2*_chunksize-1][:][:] \\
///         chunk = 2         --> data[          2*_chunksize : 3*_chunksize-1][:][:] \\
///         ...
///         chunk = nchunks-1 --> data[(nchunks-1)*buffersize : _nPulses]      [:][:] \\
/// \return Returns 'true' if new chunk was read, otherwise 'false'.
///         Updates class variable _aVal containing current chunk of data
///
/// \param &Options          [in] Global model options information
/// \param global_model_time [in] current simulation time (in days)
//
bool CForcingGrid::ReadData(const optStruct   &Options,
                            const double       global_model_time)
{
  bool new_chunk_read=false;  // true if new chunk was read, otherwise false

#ifdef _RVNETCDF_

  // local variables
  int     ir,ic,it;
  int     ncid;          // file unit
  int     dim1;          // length of 1st dimension in NetCDF data
  int     dim2;          // length of 2nd dimension in NetCDF data
  int     dim3;          // length of 3rd dimension in NetCDF data

  int     varid_f;       // id of forcing variable read
  double  missval;       // value of "missing_value" attribute of forcing variable 
  double  fillval;       // value of "_FillValue"    attribute of forcing variable
  double  add_offset;    // value of "add_offset"    attribute of forcing variable
  double  scale_factor;  // value of "scale_factor"  attribute of forcing variable
  size_t  att_len;       // length of the attribute's text
  nc_type att_type;      // type of attribute
  int     retval;        // error value for NetCDF routines
  int     iChunkSize;    // size of current chunk; always equal _ChunkSize except for last chunk in file (might be shorter)
  int     iChunk_new;    // chunk in which current model time step falls

  // check if chunk id is valid
  // -------------------------------
  ExitGracefullyIf((_iChunk>=_nChunk) || (_iChunk <-1),"CForcingGrid: ReadData: this is not a valid chunk",BAD_DATA);

  // -------------------------------
  // only initialization required
  // -------------------------------
  if (_iChunk == -1)
  {
    Initialize(Options);

    // allocate _aVal matrix using maximum chunk size
    // -------------------------------
    _aVal = NULL;
    _aVal = new double *[_ChunkSize];
    for (it=0; it<_ChunkSize; it++) {                    // loop over time points in buffer
      _aVal[it]=NULL;
      _aVal[it] = new double [_nNonZeroWeightedGridCells];
      ExitGracefullyIf(_aVal[it]==NULL,"CForcingGrid::ReadData",OUT_OF_MEMORY);
      for (ic=0; ic<_nNonZeroWeightedGridCells;ic++){    // loop over all non-zero weighted grid cells
        _aVal[it][ic]=NETCDF_BLANK_VALUE;
      }
    }

    // set _is_derived_data to False because data are truely read from a file
    // -------------------------------
    _is_derived = false;

  } // That's all to do if chunk == -1

  // --------------------------------------------------------------------------------------------
  // Now read first proper chunk (if possible)
  // --------------------------------------------------------------------------------------------

  iChunk_new = max(int(floor((global_model_time+TIME_CORRECTION) / (_interval * _ChunkSize))),0);  

  // check if given model time step is covered by current chunk; if yes, do nothing; if no,  read next chunk
  if(_iChunk != iChunk_new)
  {  
    if(Options.noisy){ 
      cout<<endl<<" Start reading new chunk... iChunk = "<<iChunk_new<<" (var = "<<_varname.c_str()<<")"<<endl;
      time_struct tt_tmp;
      JulianConvert(global_model_time,Options.julian_start_day,Options.julian_start_year,Options.calendar,tt_tmp);
      cout<<tt_tmp.date_string<<endl;
    }
    
    _iChunk = iChunk_new;

    // check if chunk id is valid
    // -------------------------------
    ExitGracefullyIf((_iChunk>=_nChunk) || (_iChunk<-1),"CForcingGrid: ReadData: this is not a valid chunk",BAD_DATA);

    // determine chunk size 
    // -------------------------------
    iChunkSize = min(_ChunkSize,int((Options.duration - global_model_time) / _interval));
    
    // Open NetCDF file, Get the id of the forcing data, varid_f
    // -------------------------------
    string filename_e=_filename;
    SubstringReplace(filename_e,"*",to_string(g_current_e+1)); //replaces wildcard for ensemble runs

    retval = nc_open(filename_e.c_str(),NC_NOWRITE,&ncid);      HandleNetCDFErrors(retval);

    string varname_e=_varname;
    SubstringReplace(varname_e,"*",to_string(g_current_e+1)); //replaces wildcard for ensemble runs

    retval = nc_inq_varid(ncid,varname_e.c_str(),&varid_f);     HandleNetCDFErrors(retval);

    // find "_FillValue" of forcing data
    // -------------------------------
    fillval = NETCDF_BLANK_VALUE; //Default
    retval = nc_inq_att(ncid, varid_f, "_FillValue", &att_type, &att_len);      
    if (retval != NC_ENOTATT) {
      HandleNetCDFErrors(retval);
      retval = nc_get_att_double(ncid, varid_f, "_FillValue", &fillval);       HandleNetCDFErrors(retval);// read attribute value
    }

    // find "missing_value" of forcing data
    // -------------------------------
    missval = NETCDF_BLANK_VALUE; //Default
    retval = nc_inq_att(ncid, varid_f, "missing_value", &att_type, &att_len);      
    if (retval != NC_ENOTATT) {
      HandleNetCDFErrors(retval);
      retval = nc_get_att_double(ncid, varid_f, "missing_value", &missval);     HandleNetCDFErrors(retval);// read attribute value
    }

    // check for attributes "add_offset" of forcing data
    // -------------------------------
    add_offset = 0.0;
    retval = nc_inq_att(ncid, varid_f, "add_offset", &att_type, &att_len);
    if (retval != NC_ENOTATT) {
      HandleNetCDFErrors(retval);
      retval = nc_get_att_double(ncid, varid_f, "add_offset", &add_offset);       HandleNetCDFErrors(retval);// read attribute value
    }
    
    // check for attributes "scale_factor" of forcing data
    // -------------------------------
    scale_factor = 1.0;
    retval = nc_inq_att(ncid, varid_f, "scale_factor", &att_type, &att_len);
    if (retval != NC_ENOTATT) {
      HandleNetCDFErrors(retval);
      retval = nc_get_att_double(ncid, varid_f, "scale_factor", &scale_factor);       HandleNetCDFErrors(retval);// read attribute value
    }
    if (Options.noisy){ 
      cout << "iChunksize:  = " << iChunkSize   << endl;
      cout << "add_offset   = " << add_offset   << endl;
      cout << "scale_factor = " << scale_factor << endl;
    }

    // allocate aTmp matrix
    // -------------------------------
    dim1 = 1; dim2 = 1; dim3 = 1;

    if ( _is_3D ) {
      switch(_dim_order)
      {
      case(1):
        dim1 = _WinLength[0]; dim2 = _WinLength[1]; dim3 = iChunkSize;    break; // dimensions are (x,y,t)
      case(2):
        dim1 = _WinLength[1]; dim2 = _WinLength[0]; dim3 = iChunkSize;    break; // dimensions are (y,x,t)
      case(3):
        dim1 = _WinLength[0]; dim2 = iChunkSize;    dim3 = _WinLength[1]; break; // dimensions are (x,t,y)
      case(4):
        dim1 = iChunkSize;    dim2 = _WinLength[0]; dim3 = _WinLength[1]; break; // dimensions are (t,x,y)
      case(5):
        dim1 = _WinLength[1]; dim2 = iChunkSize;    dim3 = _WinLength[0]; break; // dimensions are (y,t,x)
      case(6):
        dim1 = iChunkSize;    dim2 = _WinLength[1]; dim3 = _WinLength[0]; break; // dimensions are (t,y,x)
      }
    }
    else {
      switch(_dim_order)
      {
      case(1):
        dim1 = _GridDims[0]; dim2 = iChunkSize;   dim3 = 1; break; // dimensions are (station,t)
      case(2):
        dim1 = iChunkSize;   dim2 = _GridDims[0]; dim3 = 1; break; // dimensions are (t, station)
      }
    }

    // -------------------------------
    // emulate VLA 3D array storage - store 3D array as vector using Row Major Order
    // -------------------------------
    double *aVec=NULL;
    aVec=new double[dim1*dim2*dim3];//stores actual data
    ExitGracefullyIf(aVec==NULL,"CForcingGrid::ReadData : aVec",OUT_OF_MEMORY);
    for(int i=0; i<dim1*dim2*dim3; i++) {
      aVec[i]=NETCDF_BLANK_VALUE;
    }

    double ***aTmp3D=NULL; //stores pointers to rows/columns of 3D data
    double  **aTmp2D=NULL; //stores pointers to rows/columns of 2D data
    if ( _is_3D ) {
      aTmp3D=new double **[dim1];
      ExitGracefullyIf(aTmp3D==NULL,"CForcingGrid::ReadData : aTmp3D(0)",OUT_OF_MEMORY);
      for(it=0;it<dim1;it++){
        aTmp3D[it]=NULL;
        aTmp3D[it]=new double *[dim2];
        ExitGracefullyIf(aTmp3D[it]==NULL,"CForcingGrid::ReadData : aTmp3D(1)",OUT_OF_MEMORY);
        for(ir=0;ir<dim2;ir++){
          aTmp3D[it][ir]=&aVec[it*dim2*dim3+ir*dim3]; //points to correct location in aVec data storage
        }
      }
    }
    else {
      aTmp2D=new double *[dim1];
      ExitGracefullyIf(aTmp2D==NULL,"CForcingGrid::ReadData : aTmp2D(0)",OUT_OF_MEMORY);
      for(it=0;it<dim1;it++){
        aTmp2D[it]=&aVec[it*dim2]; //points to correct location in aVec data storage
      }
    }

    // Read chunk of data.
    // -------------------------------
    if ( _is_3D ) 
    {
      int       start_point = _ChunkSize * _iChunk+(int)(_t_corr/_interval);;//JRC_TIME_FIX:
      size_t    nc_start [3];
      size_t    nc_length[3];
      ptrdiff_t nc_stride[3];

      nc_length[0] = (size_t)(dim1); nc_stride[0] = 1; 
      nc_length[1] = (size_t)(dim2); nc_stride[1] = 1;
      nc_length[2] = (size_t)(dim3); nc_stride[2] = 1;
      
      switch(_dim_order) {
      case(1): // dimensions are (x,y,t)
        nc_start[0]  = (size_t)(_WinStart[0]);  nc_start[1]  = (size_t)(_WinStart[1]);  nc_start[2]  = (size_t)(start_point);
        break;
      case(2): // dimensions are (y,x,t)
        nc_start[0]  = (size_t)(_WinStart[1]);  nc_start[1]  = (size_t)(_WinStart[0]);  nc_start[2]  = (size_t)(start_point); 
        break;
      case(3): // dimensions are (x,t,y)
        nc_start[0]  = (size_t)(_WinStart[0]);  nc_start[1]  = (size_t)(start_point);   nc_start[2]  = (size_t)(_WinStart[1]);
        break;
      case(4): // dimensions are (t,x,y)
        nc_start[0]  = (size_t)(start_point);   nc_start[1]  = (size_t)(_WinStart[0]);  nc_start[2]  = (size_t)(_WinStart[1]);
        break;
      case(5): // dimensions are (y,t,x)
        nc_start[0]  = (size_t)(_WinStart[1]);  nc_start[1]  = (size_t)(start_point);   nc_start[2]  = (size_t)(_WinStart[0]);
        break;
      case(6): // dimensions are (t,y,x)
        nc_start[0]  = (size_t)(start_point);   nc_start[1]  = (size_t)(_WinStart[1]);  nc_start[2]  = (size_t)(_WinStart[0]);
        break;
      }

      //Read giant chunk of data from NetCDF (this is the bottleneck of this code)
      retval=nc_get_vars_double(ncid,varid_f,nc_start,nc_length,nc_stride,&aTmp3D[0][0][0]);   HandleNetCDFErrors(retval);
      new_chunk_read = true;

      if (Options.noisy) {
        cout<<" CForcingGrid::ReadData - is3D"<<endl;
        printf("  Dim of chunk read: dim3 = %i   dim2 = %i   dim1 = %i\n",dim3,dim2,dim1);
        printf("  start  chunk: (%zu, %zu, %zu)\n", nc_start[0], nc_start[1], nc_start[2]);
        printf("  length chunk: (%zu, %zu, %zu)\n",nc_length[0],nc_length[1],nc_length[2]);
        printf("  stride chunk: (%zu, %zu, %zu)\n",nc_stride[0],nc_stride[1],nc_stride[2]);
      }
    }
    else //2D
    {
      int       start_point = _ChunkSize * _iChunk+(int)(_t_corr/_interval);//JRC_TIME_FIX:
      size_t    nc_start[2];
      size_t    nc_length[2];
      ptrdiff_t nc_stride[2];

      nc_length[0] = (size_t)(dim1); nc_stride[0] = 1;
      nc_length[1] = (size_t)(dim2); nc_stride[1] = 1;

      switch(_dim_order) {
        case(1): // dimensions are (station,t)
          nc_start[0]  = 0;
          nc_start[1]  = (size_t)(start_point);
          break;
        case(2): // dimensions are (t,station)
          nc_start[0]  = (size_t)(start_point);
          nc_start[1]  = 0;
          break;
      }
    
      //Read from NetCDF (this is the bottleneck of this code)
      retval=nc_get_vars_double(ncid,varid_f,nc_start,nc_length,nc_stride,&aTmp2D[0][0]);
      HandleNetCDFErrors(retval);
      new_chunk_read = true;

      if (Options.noisy) {
        cout<<" CForcingGrid::ReadData - !is3D"<<endl;
        printf("  Dim of chunk read: dim2 = %i   dim1 = %i\n",dim2,dim1);
        printf("  start  chunk: (%zu, %zu)\n", nc_start[0], nc_start[1]);
        printf("  length chunk: (%zu, %zu)\n",nc_length[0],nc_length[1]);
        printf("  stride chunk: (%zu, %zu)\n",nc_stride[0],nc_stride[1]);
      }
    }

    // Re-scale NetCDF variables based on their internal add-offset and scale_factor
    // MANDATORY to do before any value of these data are used
    // -------------------------------
    if ( _is_3D ) {
      for (it=0;it<dim1;it++){
        for (ir=0;ir<dim2;ir++){
	        for (ic=0;ic<dim3;ic++){
	          aTmp3D[it][ir][ic] = aTmp3D[it][ir][ic] * scale_factor + add_offset;
            //if  ((it==0) || (it==dim1-1)) {cout<<setprecision(4)<<setw(9)<<aTmp3D[it][ir][ic]<<" ";}
	        }
        //if ((it==0) || (it==dim1-1)) { cout << endl; }
        }
       //if  ((it==0) || (it==dim1-1)) {cout<<endl<<endl;}
      }
    }
    else {
      for (it=0;it<dim1;it++){
	      for (ir=0;ir<dim2;ir++){
	        aTmp2D[it][ir] = aTmp2D[it][ir] * scale_factor + add_offset;
	      }
      }
    }
    

    // Copy all data from aTmp array to member array _aVal.
    // -------------------------------
    double val;
    if ( _is_3D ) 
    {
      int irow,icol;
      if (_dim_order == 1) {
        for (it=0; it<iChunkSize; it++){                     // loop over time points in buffer
          for (ic=0; ic<_nNonZeroWeightedGridCells; ic++){   // loop over non-zero weighted grid cells
            CellIdxToRowCol(_IdxNonZeroGridCells[ic],irow,icol); 
            val=aTmp3D[icol-_WinStart[0]][irow-_WinStart[1]][it];
            if(val==missval) { CheckValue3D(val,missval,it,irow,icol); }
            if(val==fillval) { CheckValue3D(val,fillval,it,irow,icol); }
            _aVal[it][ic]=_LinTrans_a*val+_LinTrans_b;

          }
        }
      }
      else if (_dim_order == 2) {
        for (it=0; it<iChunkSize; it++){                     // loop over time points in buffer
          for (ic=0; ic<_nNonZeroWeightedGridCells; ic++){   // loop over non-zero weighted grid cells
		        CellIdxToRowCol(_IdxNonZeroGridCells[ic],irow,icol);
		        val=aTmp3D[irow-_WinStart[1]][icol-_WinStart[0]][it];
            if(val==missval) { CheckValue3D(val,missval,it,irow,icol); }
            if(val==fillval) { CheckValue3D(val,fillval,it,irow,icol); }
            _aVal[it][ic]=_LinTrans_a*val+_LinTrans_b;
		      }
		    }
      }
      else if (_dim_order == 3) {
        for (it=0; it<iChunkSize; it++){                     // loop over time points in buffer
          for (ic=0; ic<_nNonZeroWeightedGridCells; ic++){   // loop over non-zero weighted grid cells
            CellIdxToRowCol(_IdxNonZeroGridCells[ic],irow,icol);
            val=aTmp3D[icol-_WinStart[0]][it][irow-_WinStart[1]];
            if(val==missval) { CheckValue3D(val,missval,it,irow,icol); }
            if(val==fillval) { CheckValue3D(val,fillval,it,irow,icol); }
            _aVal[it][ic]=_LinTrans_a*val+_LinTrans_b;
          }
        }
      }
      else if (_dim_order == 4) {
        for (it=0; it<iChunkSize; it++){                     // loop over time points in buffer
          for (ic=0; ic<_nNonZeroWeightedGridCells; ic++){   // loop over non-zero weighted grid cells
            CellIdxToRowCol(_IdxNonZeroGridCells[ic],irow,icol);
            val=aTmp3D[it][icol-_WinStart[0]][irow-_WinStart[1]];
            if(!((Options.deltaresFEWS) && (it==0))) {
              if(val==missval) { CheckValue3D(val,missval,it,irow,icol); }
              if(val==fillval) { CheckValue3D(val,fillval,it,irow,icol); }
            }
            _aVal[it][ic]=_LinTrans_a*val+_LinTrans_b;
          }
        }
      }
      else if (_dim_order == 5) {
        for (it=0; it<iChunkSize; it++){                     // loop over time points in buffer
          for (ic=0; ic<_nNonZeroWeightedGridCells; ic++){   // loop over non-zero weighted grid cells
            CellIdxToRowCol(_IdxNonZeroGridCells[ic],irow,icol);
            val=aTmp3D[irow-_WinStart[1]][it][icol-_WinStart[0]];
            if(val==missval) { CheckValue3D(val,missval,it,irow,icol); }
            if(val==fillval) { CheckValue3D(val,fillval,it,irow,icol); }
            _aVal[it][ic]=_LinTrans_a*val+_LinTrans_b;
          }
        }
      }
      else if (_dim_order == 6) {
        for (it=0; it<iChunkSize; it++){                      // loop over time points in buffer
          for (ic=0; ic<_nNonZeroWeightedGridCells; ic++){    // loop over non-zero weighted grid cells
            CellIdxToRowCol(_IdxNonZeroGridCells[ic],irow,icol);
            val=aTmp3D[it][irow-_WinStart[1]][icol-_WinStart[0]];
            if(val==missval) { CheckValue3D(val,missval,it,irow,icol);}   
            if(val==fillval) { CheckValue3D(val,fillval,it,irow,icol); } 
            _aVal[it][ic]=_LinTrans_a*val+_LinTrans_b;
          }
        }
      }
    }
    else // 2D
    {
      if (_dim_order == 1) {
        for (it=0; it<iChunkSize; it++){                     // loop over time points in buffer
          for (ic=0; ic<_nNonZeroWeightedGridCells; ic++){   // loop over non-zero weighted grid cells
            val=aTmp2D[_IdxNonZeroGridCells[ic]][it];
            if(val==missval) { CheckValue2D(val,missval,_IdxNonZeroGridCells[ic],it); }   // throw error  if value to read in equals "missing_value"
            if(val==fillval) { CheckValue2D(val,fillval,_IdxNonZeroGridCells[ic],it); }   // throw error  if value to read in equals "_FillValue"
            _aVal[it][ic]=_LinTrans_a*val+_LinTrans_b;
          }
        }
      }
      else if (_dim_order == 2) {
        for (it=0; it<iChunkSize; it++){                     // loop over time points in buffer
          for (ic=0; ic<_nNonZeroWeightedGridCells; ic++){   // loop over non-zero weighted grid cells
            val=aTmp2D[it][_IdxNonZeroGridCells[ic]];
            if(val==missval)  { CheckValue2D(val,missval,it,_IdxNonZeroGridCells[ic]); }  // throw error if value to read in equals "missing_value"
            if(val==fillval)  { CheckValue2D(val,fillval,it,_IdxNonZeroGridCells[ic]); }  // throw error if value to read in equals "_FillValue"
            if(rvn_isnan(val)){ CheckValue2D(val,NAN,    it,_IdxNonZeroGridCells[ic]); }
            _aVal[it][ic]=_LinTrans_a*val+_LinTrans_b;
          }
        }
      }
    }

    //delete dynamic arrays
    // -------------------------------
    if ( _is_3D ) {for (it=0;it<dim1;it++){delete [] aTmp3D[it];} delete [] aTmp3D;}
    else          {delete [] aTmp2D;}
    delete [] aVec;

    // read attribute grids - lat, long, elevation of grid cells
    // -------------------------------
    if(_is_3D) {
      switch(_dim_order)
      {
        case(1): dim1 = _GridDims[0]; dim2 = _GridDims[1]; break; // dimensions are (x,y,t)->(x,y)
        case(2): dim1 = _GridDims[1]; dim2 = _GridDims[0]; break; // dimensions are (y,x,t)->(y,x)*
        case(3): dim1 = _GridDims[0]; dim2 = _GridDims[1]; break; // dimensions are (x,t,y)->(x,y)
        case(4): dim1 = _GridDims[0]; dim2 = _GridDims[1]; break; // dimensions are (t,x,y)->(x,y)
        case(5): dim1 = _GridDims[1]; dim2 = _GridDims[0]; break; // dimensions are (y,t,x)->(y,x)*
        case(6): dim1 = _GridDims[1]; dim2 = _GridDims[0]; break; // dimensions are (t,y,x)->(y,x)*
      }
    }
    else {
      dim1 = _GridDims[0]; dim2 = 1;
    }
    ReadAttGridFromNetCDF(ncid,_AttVarNames[0],dim1,dim2,_aLatitude);
    ReadAttGridFromNetCDF(ncid,_AttVarNames[1],dim1,dim2,_aLongitude);
    ReadAttGridFromNetCDF(ncid,_AttVarNames[2],dim1,dim2,_aElevation);
    //ReadAttGridFromNetCDF2(ncid,_AttVarNames[3],dim1,dim2,_aStationIDs);

    if (_aElevation!=NULL){
      for(int ic=0; ic<_nNonZeroWeightedGridCells; ic++) {
        ExitGracefullyIf(rvn_isnan(_aElevation[ic]),"CForcingGrid::ReadData - NaN elevation found in NetCDF elevation grid with non-zero HRU weight",BAD_DATA);
      }
    }

    // -------------------------------
    // Close NetCDF file
    // -------------------------------
    retval = nc_close(ncid);       HandleNetCDFErrors(retval);
    
  }// end if(_iChunk != iChunk_new)

#endif   // end #ifdef _RVNETCDF_

  return new_chunk_read;

}

///////////////////////////////////////////////////////////////////
/// \brief   Enables queries of time series values using model time
/// \details Calculates _t_corr, correction to global model time, checks for overlap
//           with model duration, resamples to model time step/day
/// \remark  t=0 corresponds to first day with recorded values at that gauge
/// checks if forcing data cover whole simulation period
/// \param   Options          [in] options structure
///
void CForcingGrid::Initialize( const optStruct &Options )    
{
  double model_start_day =Options.julian_start_day;    // fractional day of the year (here called Julian day) [days]
  int    model_start_year=Options.julian_start_year;   //         [year]
  double model_duration  =Options.duration;            //         [days]
  double model_timestep  =Options.timestep;            // delta t [days]
                                                    
  //_t_corr is number of days between model start date and forcing 
  // start date (positive if data exists before model start date)
  //------------------------------------------------------------------------------
  _t_corr = -TimeDifference(model_start_day,model_start_year,_start_day,_start_year,Options.calendar);
 
  //QA/QC: Check for overlap issues between time series duration and model duration
  //------------------------------------------------------------------------------
  double duration               = (double)(_nPulses)*_interval;
  double local_simulation_start = _t_corr;
  double local_simulation_end   = _t_corr + model_duration;

  if(Options.deltaresFEWS) { local_simulation_start-=_interval; }//first timestep in FEWS netCDF output is noise
  if ((local_simulation_start<0) && (local_simulation_start>-TIME_CORRECTION)){ local_simulation_start=0.0; }//correct for floating point errors

  if (Options.noisy){ cout << endl; }
  if (duration < local_simulation_start)  //out of data before simulation starts!
  {
    cout << "Initialize forcing grid  '" << _varname.c_str() << "'" << endl;
    cout << "  time series start day, year, duration :" << _start_day      << "," << _start_year      << " " << duration       << endl;
    cout << "  model start day, year, duration :"       << model_start_day << "," << model_start_year << " " << model_duration << endl;
   
    ExitGracefully("CForcingGrid::Initialize: gridded forcing data not available for entire model simulation duration", BAD_DATA);
  }
  if (duration + model_timestep < local_simulation_end)    //run out of data before simulation finishes
  {                                                        // +model_timestep is for coincdent duration & data
    cout << "Initialize forcing grid  '" << _varname.c_str() << "'" << endl;
    cout << "  time series start day, year, duration :" << _start_day      << "," << _start_year      << " " << duration       << endl;
    cout << "  model start day, year, duration :"       << model_start_day << "," << model_start_year << " " << model_duration << endl;
   
    ExitGracefully("CForcingGrid::Initialize: gridded forcing data not available at end of model simulation", BAD_DATA);
  }
  if ((local_simulation_start<0) || (_start_year>model_start_year))     //data does not begin until after simulation
  {
    cout << "Initialize forcing grid  '" << _varname.c_str() << "'" << endl;
    cout << "  time series start day, year, duration :" << _start_day      << "," << _start_year      << " " << duration       << endl;
    cout << "  model start day, year, duration :"       << model_start_day << "," << model_start_year << " " << model_duration << endl;
    
    ExitGracefully("CForcingGrid::Initialize: gridded forcing data not available at beginning of model simulation", BAD_DATA);
  }
}

///////////////////////////////////////////////////////////////////
/// \brief returns row and column index of cell ID
///
/// \param cellid  [in] int of cell ID
//
void CForcingGrid::CellIdxToRowCol(const int cellid, int &row, int &column) const
{
  int ncols = GetCols();
  row    = int(cellid / ncols);
  column = cellid % ncols;
  return;
}

///////////////////////////////////////////////////////////////////
/// \brief nothing returned only exits Raven when value of 3D array equals checkval and dumps an error message
///
/// \param     value    [in] double value to be checked
/// \param     checkval [in] value to be checked against
//
void CForcingGrid::CheckValue3D(const double value, const double checkval, const int dim1_idx, const int dim2_idx, const int dim3_idx)
{
  if (value == checkval) {
    cout << "Forcing grid  '" << _varname.c_str() << "'" << endl;
    cout << "   dimension 1:  idx = " << dim1_idx << endl;
    cout << "   dimension 2:  idx = " << dim2_idx << endl;
    cout << "   dimension 3:  idx = " << dim3_idx << endl;
    cout << "   check value:        " << checkval << endl;
    cout << "   data  value:        " << value    << endl;
    ExitGracefully("CForcingGrid::ReadData: 3D forcing data contain missing or fill values", BAD_DATA);
  }
  return;
}

///////////////////////////////////////////////////////////////////
/// \brief nothing returned only exits Raven when value of 2D array equals checkval and dumps an error message
///
/// \param     value    [in] double value to be checked
/// \param     checkval [in] value to be checked against
//
void CForcingGrid::CheckValue2D(const double value, const double checkval, const int dim1_idx, const int dim2_idx)
{
  if (value == checkval) {
    cout << "Forcing grid  '" << _varname.c_str() << "'" << endl;
    cout << "   dimension 1:  idx = " << dim1_idx << endl;
    cout << "   dimension 2:  idx = " << dim2_idx << endl;
    cout << "   check value:        " << checkval << endl;
    cout << "   data  value:        " << value    << endl;
    ExitGracefully("CForcingGrid::ReadData: 2D forcing data contain missing or fill values", BAD_DATA);
  }
  return;
}

///////////////////////////////////////////////////////////////////
/// \brief sets the _ForcingType in class CForcingGrid
///
/// \param forcingtype  [in] string of forcing type, e.g. PRECIP
//
void CForcingGrid::SetForcingType(const forcing_type &ForcingType)
{
  _ForcingType=ForcingType;
}

///////////////////////////////////////////////////////////////////
/// \brief sets the _Filename in class CForcingGrid
///
/// \param filename  [in] filename of NetCDF file containing gridded forcing data
//
void CForcingGrid::SetFilename(const string filename)
{
  _filename=filename;
}

///////////////////////////////////////////////////////////////////
/// \brief sets the _varname in class CForcingGrid
///
/// \param varname  [in] string of variable name in NetCDF file
//
void CForcingGrid::SetVarname(const string varname)
{
  _varname=varname;
}

///////////////////////////////////////////////////////////////////
/// \brief sets the _DimNames in class CForcingGrid
///
/// \param DimNames  [in] array of three strings containing the names of the dimensions NetCDF file
///                       order must be: x-dimension, y-dimension, time-dimension
//
void CForcingGrid::SetDimNames(const string DimNames[3])
{
  _DimNames[0]=DimNames[0];
  _DimNames[1]=DimNames[1];
  _DimNames[2]=DimNames[2];
}

///////////////////////////////////////////////////////////////////
/// \brief sets the _SetGridDims in class CForcingGrid
///
/// \param GridDims  [in] array of three intgers determining length of each dimension (x=#cols, y=#rows, t=#timepoints)
///                          order must be: x-dimension, y-dimension, time-dimension
//
void CForcingGrid::SetGridDims(const int GridDims[3])
{
  _GridDims[0]=GridDims[0];
  _GridDims[1]=GridDims[1];
  _GridDims[2]=GridDims[2];
}
///////////////////////////////////////////////////////////////////
/// \brief allocated memory for and populates the _aLastNonZeroWt and _IdxNonZeroGridCells arrays, Calculates _nNonZeroWeightedGridCells
/// called once after :GridWeights read in. 
///
/// \param nHydroUnits number of HRUs
/// \param nGridCells number of grid cells
//
void CForcingGrid::SetIdxNonZeroGridCells(const int nHydroUnits, const int nGridCells, const optStruct &Options)
{
  int row, col;
  int minrow=_GridDims[1];
  int maxrow=0;
  int mincol=_GridDims[0];
  int maxcol=0;

  bool *nonzero= NULL;
  nonzero =  new bool [nGridCells]; 
  for (int il=0; il<nGridCells; il++) { // loop over all cells of NetCDF
    nonzero[il] = false;
  }

  if (_GridWeight != NULL){
    for(int k=0; k<nHydroUnits; k++) {  // loop over HRUs
      for(int i=0; i<_nWeights[k]; i++) { // loop over all cells of NetCDF
        if(_GridWeight[k][i] > 0.00001) {
          nonzero[_GridWtCellIDs[k][i]] = true;
          CellIdxToRowCol(_GridWtCellIDs[k][i],row,col);
          if(row>maxrow) { maxrow=row; }
          if(col>maxcol) { maxcol=col; }
          if(row<minrow) { minrow=row; }
          if(col<mincol) { mincol=col; }

        }
      }
    }
  }
  else{
    ExitGracefully(
      "CForcingGrid: SetIdxNonZeroGridCells: _GridWeight is not allocated yet. Call AllocateWeightArray(nHRUs) first.", RUNTIME_ERR);
  }

  _WinLength[0]=maxcol-mincol+1;
  _WinLength[1]=maxrow-minrow+1;
  _WinLength[2]=_GridDims[2];
  _WinStart [0]=mincol;
  _WinStart [1]=minrow;
  _WinStart [2]=0; //temporary - this shifts over course of simualtion
  
  //To remove support for local window:
  //_WinLength[0]=_GridDims[0];_WinStart[0]=0;
  //_WinLength[1]=_GridDims[1];_WinStart[1]=0;
  
  /*cout<<"INITIALIZING GRID WINDOW"<<endl;
  cout<<" start   (col, row): ("<<_WinStart[0]              <<", "<<_WinStart[1]              <<")"<<endl;
  cout<<" end     (col, row): ("<<maxcol                    <<", "<<maxrow<<")"<<endl;
  cout<<" lengths (col, row): ("<<_WinLength[0]             <<", "<<_WinLength[1]             <<")"<<endl;
  cout<<"-----------------------------"<<endl;*/

  // count number of non-zero weighted grid cells 
  _nNonZeroWeightedGridCells=0;
  for (int il=0; il<nGridCells; il++) { // loop over all cells of NetCDF
    if ( nonzero[il] ) { _nNonZeroWeightedGridCells++; }
  }

  _IdxNonZeroGridCells = NULL;
  _IdxNonZeroGridCells = new int [_nNonZeroWeightedGridCells];
  for (int il=0; il<_nNonZeroWeightedGridCells; il++) { // loop over all cells non-zero weighted grid cells
    _IdxNonZeroGridCells[il] = -1;
  }
  _CellIDToIdx = NULL;
  _CellIDToIdx = new int[nGridCells];
  for(int c=0; c<nGridCells; c++) { // loop over all cells non-zero weighted grid cells
    _CellIDToIdx[c] = DOESNT_EXIST;
  }

  int ic = 0;
  for (int il=0; il<nGridCells; il++) { // loop over all cells of NetCDF
    if ( nonzero[il] ) {
      _IdxNonZeroGridCells[ic] = il;
      _CellIDToIdx[il]=ic;
      ic++;
    }
  }
  delete[] nonzero;

  if (Options.noisy){
    cout<<"Finished SetIdxNonZeroGridCells routine, # of non-zero weighted cells: "<<_nNonZeroWeightedGridCells<<endl;
  }
}

///////////////////////////////////////////////////////////////////
/// \brief calculates _ChunkSize and total number of chunks to read _nChunks, sets the id of the current chunk 
/// depending upon size of grid (_nNonZeroWeightedGridCells*buffersize*8byte <=  10 MB=10*1024*1024 byte)
/// needs to be called after SetIdxNonZeroGridCells()
///
/// \param nHydroUnits number of HRUs
/// \param nGridCells number of grid cells
//
void CForcingGrid::CalculateChunkSize(const optStruct& Options) 
{  
  if(_nNonZeroWeightedGridCells==0) { //Never called SetIdxNonZeroGridCells()
    ExitGracefully("CalculateChunkSize: called at wrong location",RUNTIME_ERR);
  }

  int ntime;
  if(_is_3D) { ntime = _GridDims[2]; }
  else       { ntime = _GridDims[1]; }

  int    BytesPerTimestep;      // Memory requirement for one timestep of gridded forcing file [Bytes]
  if(_is_3D) { BytesPerTimestep = 8 * _WinLength[0] * _WinLength[1]; }
  else       { BytesPerTimestep = 8 * _WinLength[0]; }

  //BytesPerTimestep = 8 * _nNonZeroWeightedGridCells; //?? amount actually stored in memory?

  int CHUNK_MEMORY=Options.NetCDF_chunk_mem*1024 * 1024;

  int tmpChunkSize;

  tmpChunkSize = (int)(max(min(CHUNK_MEMORY  / BytesPerTimestep,ntime),1));  // number of timesteps per chunk - 10MB chunks

  if(!Options.deltaresFEWS) {
    tmpChunkSize = (int)((int)(tmpChunkSize*_interval)/_interval);               // make sure chunks are complete days (have to relax for FEWS)
  }  
  
  tmpChunkSize = max(int(rvn_round(1.0/_interval)),tmpChunkSize);            // make sure  at least one day is read
                                                                                 // support larger chunk if model duration is small
  double partday=Options.julian_start_day-floor(Options.julian_start_day);
  tmpChunkSize = min(tmpChunkSize,(int) ceil(ceil(Options.duration+partday)/_interval));//ensures goes to midnight of last day 

  SetChunkSize(tmpChunkSize);

  _nChunk    = int(ceil((Options.duration/_interval)/_ChunkSize));                      // total number of chunks
  
  if (Options.noisy){
    cout<<"Finished CalculateChunkSize routine,     # of time steps per chunk:    "<<_ChunkSize<<endl;
    cout<<"                                         # of time chunks:             "<<_nChunk   <<endl;
  }
  /*cout<<"ntime: "<<ntime<<endl;*/
}
///////////////////////////////////////////////////////////////////
/// \brief sets the _nHydroUnits in class CForcingGrid
///
/// \param nHydroUnits  [in] number of HRUs
//
void CForcingGrid::SetnHydroUnits(const int nHydroUnits)
{
  _nHydroUnits=nHydroUnits;
}

///////////////////////////////////////////////////////////////////
/// \brief sets the _ChunkSize in class CForcingGrid
///
/// \param ChunkSize  [in] size of currently read chunk
//
void CForcingGrid::SetChunkSize(const int    ChunkSize)
{
  _ChunkSize=ChunkSize;
}

///////////////////////////////////////////////////////////////////
/// \brief sets the _interval in class CForcingGrid
///
/// \param interval  [in] delta t of this specific gridded forcing data
//
void CForcingGrid::SetInterval(const double interval)
{
  _interval=interval;
  _steps_per_day=(int)(rvn_round(1.0/_interval));
}

///////////////////////////////////////////////////////////////////
/// \brief sets the _snowfall_corr in class CForcingGrid
///
/// \param snowfall_corr  [in] value giving the snowfall correction (default 1.0)
//
void CForcingGrid::SetSnowfallCorr(const double snowfall_corr)
{
  _snowfall_corr = snowfall_corr;
}

///////////////////////////////////////////////////////////////////
/// \brief sets the _rainfall_corr in class CForcingGrid
///
/// \param rainfall_corr  [in] value giving the rainfall correction (default 1.0)
//
void CForcingGrid::SetRainfallCorr(const double rainfall_corr)
{
  _rainfall_corr = rainfall_corr;
}

///////////////////////////////////////////////////////////////////
/// \brief sets the _cloud_min_temp in class CForcingGrid
///
/// \param cloud_min_temp  [in] value giving the minimum temperature threshold used to determine
///                             cloud_cover factor
///                             (default -20.0, ensures cloud-free status always unless overriden)
//
void CForcingGrid::Setcloud_min_temp(const double cloud_min_temp){
  _cloud_min_temp=cloud_min_temp;
}

///////////////////////////////////////////////////////////////////
/// \brief sets the _cloud_max_temp in class CForcingGrid
///
/// \param cloud_max_temp  [in] value giving the maximum temperature threshold used to determine
///                             cloud_cover factor
///                             (default -20.0, ensures cloud-free status always unless overriden)
//
void CForcingGrid::Setcloud_max_temp(const double cloud_max_temp){
  _cloud_max_temp=cloud_max_temp;
}

///////////////////////////////////////////////////////////////////
/// \brief sets the _aAveTemp[12] in class CForcingGrid
///
/// \param aAveTemp[12]  [in] value giving representative average monthly temperatures [C]
//
void CForcingGrid::SetaAveTemp(const double aAveTemp[12]){
  for (int ii=0; ii<12; ii++) {
    _aAveTemp[ii] = aAveTemp[ii];
  }
}

///////////////////////////////////////////////////////////////////
/// \brief sets the _aMinTemp[12] in class CForcingGrid
///
/// \param aMinTemp[12]  [in] value giving representative minimum monthly temperatures [C]
//
void CForcingGrid::SetaMinTemp(const double aMinTemp[12]){
  for (int ii=0; ii<12; ii++) {
    _aMinTemp[ii] = aMinTemp[ii];
  }
}

///////////////////////////////////////////////////////////////////
/// \brief sets the _aMaxTemp[12] in class CForcingGrid
///
/// \param aMaxTemp[12]  [in] value giving representative maximum monthly temperatures [C]
//
void CForcingGrid::SetaMaxTemp(const double aMaxTemp[12]){
  for (int ii=0; ii<12; ii++) {
    _aMaxTemp[ii] = aMaxTemp[ii];
  }
}

///////////////////////////////////////////////////////////////////
/// \brief sets the _aAvePET [12] in class CForcingGrid
///
/// \param aAvePET [12]  [in] value giving representative average monthly PET [mm/d]
///                            (or monthly PET factor [mm/d/K], if MONTHLY_FACTOR is used)
//
void CForcingGrid::SetaAvePET(const double aAvePET [12]){
  for (int ii=0; ii<12; ii++) {
    _aAvePET[ii] = aAvePET[ii];
  }
}

///////////////////////////////////////////////////////////////////
/// \brief sets the _aVal in class CForcingGrid
///
/// \param ic    [in] Index of grid cell with non-zero weighting (value between 0 and _nNonZeroWeightedGridCells)
/// \param it    [in] time   index of _aVal[t][idx]
/// \param aVal  [in] value to be set
//
void CForcingGrid::SetValue( const int ic, const int it, const double aVal) {
#ifdef _STRICTCHECK_
  if(it>=_ChunkSize) {
    ExitGracefully("CForcingGrid::SetValue:invalid index",RUNTIME_ERR);}
  if(ic>=_nNonZeroWeightedGridCells) {
    ExitGracefully("CForcingGrid::SetValue:invalid index",RUNTIME_ERR);}
#endif
  _aVal[it][ic] = aVal;
}

///////////////////////////////////////////////////////////////////
/// \brief indicates if NetCDF data are 3D (lat,lon,time)   --> True or
///                                 are 2D (nstations,time) --> False
//
void CForcingGrid::SetIs3D( bool is_3D ){
  _is_3D=is_3D;
}

///////////////////////////////////////////////////////////////////
/// \brief indicates that the forcing grid stores an accumulated variable
//
void CForcingGrid::SetToDeaccumulate(){
  _deaccumulate=true;
}

///////////////////////////////////////////////////////////////////
/// \brief sets the time shift (TimeShift) in class CForcingGrid
///
/// \param TimeShift  [in] float in fractional days showing how read in data
///                        timestamps need to be shifted
//
void CForcingGrid::SetTimeShift(const double TimeShift)
{
  _TimeShift = TimeShift;
}

///////////////////////////////////////////////////////////////////
/// \brief sets as period ending
///
//
void CForcingGrid::SetAsPeriodEnding()
{
  _period_ending=true;
}

///////////////////////////////////////////////////////////////////
/// \brief sets NetCDF variable names for lat, long, or elevation
/// \param var [in] one of "Latitude","Longitude", or "Elevation"
/// \param varname [in] variable name
///  can be done anytime prior to calling ReadData()
//
void  CForcingGrid::SetAttributeVarName(const string var,const string varname) 
{
  if      (var=="Latitude" ) { _AttVarNames[0]=varname;}
  else if (var=="Longitude") { _AttVarNames[1]=varname; }
  else if (var=="Elevation") { _AttVarNames[2]=varname; }
  else if (var=="StationIDs"){ _AttVarNames[3]=varname; }
}

///////////////////////////////////////////////////////////////////
/// \brief sets linear transform parameters a and b in class CForcingGrid
///        to allow for linear transformation new = a*data + b of read-in data
///
/// \param LinTrans_a  [in] linear transformation of read data: new = a*data + b
/// \param LinTrans_b  [in] linear transformation of read data: new = a*data + b
//
void CForcingGrid::SetLinearTransform(const double LinTrans_a, const double LinTrans_b)
{
  _LinTrans_a = LinTrans_a;
  _LinTrans_b = LinTrans_b;
}

///////////////////////////////////////////////////////////////////
/// \brief allocates _GridWeight [nHydroUnits, nGridCells] in class CForcingGrid
///
/// \param nHydroUnits [in] Number of HRUs must be given
///                         (is read from :GridWeights block in *.rvi file)
/// \param nGridCells  [in] Number of Grid Cells must be given
///                         (is read from :GridWeights block in *.rvi file)
//
void   CForcingGrid::AllocateWeightArray(const int nHydroUnits, const int nGridCells)
{
  // -------------------------------
  // Initialize weighting array and set all entries to zero (default value)
  // since only non-zero entries need to be specified in :GridWeight block in *.rvi file
  // -------------------------------
  //cout<<"Creating new GridWeights array (Base Constructor): "<<ForcingToString(_ForcingType)<<endl;
  
  _GridWeight = NULL;
  _GridWeight = new double *[nHydroUnits];
  ExitGracefullyIf(_GridWeight==NULL   ,"AllocateWeightArray(1)",OUT_OF_MEMORY);
  _GridWtCellIDs=NULL;
  _GridWtCellIDs = new int *[nHydroUnits];
  ExitGracefullyIf(_GridWtCellIDs==NULL,"AllocateWeightArray(2)",OUT_OF_MEMORY);
  _nWeights=NULL;
  _nWeights = new int [nHydroUnits];
  ExitGracefullyIf(_nWeights==NULL     ,"AllocateWeightArray(3)",OUT_OF_MEMORY);
  for(int k=0; k<nHydroUnits; k++) {  
    _GridWeight   [k] =NULL;
    _GridWtCellIDs[k] =NULL;
    _nWeights     [k] =0;
  }


}

///////////////////////////////////////////////////////////////////
/// \brief sets one entry of _GridWeight[HRUID][i] = weight in class CForcingGrid; dynamically grows gridweight
///
/// \param HRUID  [in] HRU id (numbered: 0,1,2,...,nHydroUnits-1)
/// \param CellID [in] cell ID in NetCDF; cells are numbered linewise from left to right starting with 0:
///                         0 1 2 3 4 \n
///                         5 6 7 8 9 \n
///                         ...
/// \param weight [in] weight[k][l] is fraction of forcing for HRUID k is from grid cell l=(i,j)
///                    and zero-based grid cell index l is derived by l = j * NC + i
///                    where i and j are the column and row of cell l respectively and
///                    NC is the total number of columns.
///                    Following contraint must be satisfied:
///                        sum(w_kl, {l=1,NC*NR}) = 1.0 for all HRUs k
//                     This will be checked with routine CheckWeightArray().
//
void   CForcingGrid::SetWeightVal(const int HRUID,
                                  const int CellID,
                                  const double weight)
{
#ifdef _STRICTCHECK_
  if (_GridWeight == NULL){
    ExitGracefully(
      "CForcingGrid: SetWeightVal: _GridWeight is not allocated yet. Call AllocateWeightArray(nHRUs) first.",RUNTIME_ERR);
  }
  int ncells;
  if(_is_3D) { ncells = _GridDims[0] * _GridDims[1]; }
  else       { ncells = _GridDims[0]; }

  if((HRUID<0) || (HRUID>=_nHydroUnits)) {
    ExitGracefully("CForcingGrid: SetWeightVal: invalid HRU ID",BAD_DATA);}
  if((CellID<0) || (CellID>=ncells)) {
    ExitGracefully("CForcingGrid: SetWeightVal: invalid cell ID",BAD_DATA);}
#endif
  int k=HRUID;

  //entry already exists in sparse weights matrix, replace
  for(int i=0;i<_nWeights[k];i++) {
    if (_GridWtCellIDs[k][i]==CellID){
      _GridWeight[k][i]=weight; return;
    }
  }

  //dynamically grow weights matrix
  _nWeights[k]++;
  int    *tmpid=new int   [_nWeights[k]];
  double *tmpwt=new double[_nWeights[k]];
  for(int i=0;i<_nWeights[k]-1;i++) { 
    tmpwt[i]=_GridWeight   [k][i]; 
    tmpid[i]=_GridWtCellIDs[k][i];
  } 
  tmpwt[_nWeights[k]-1]=weight; //append to end
  tmpid[_nWeights[k]-1]=CellID;

  delete [] _GridWeight   [k]; _GridWeight   [k]=tmpwt;
  delete [] _GridWtCellIDs[k]; _GridWtCellIDs[k]=tmpid;
}
///////////////////////////////////////////////////////////////////
/// \brief sets one entry of _aElevation[CellID] 
//
/// \param cellID [in] cell ID/stationID in NetCDF (from 0 to ncells-1)
/// \param elev [in] elevation of cell in NetCDF

void   CForcingGrid::SetStationElevation(const int CellID,const double &elev) 
{
  int ncells;
  if(_is_3D) { ncells = _GridDims[0] * _GridDims[1]; }
  else       { ncells = _GridDims[0]; }

  if((CellID<0) || (CellID>=ncells)) {
    ExitGracefully("CForcingGrid: SetStationElevation: invalid cell/station identifier (likely in :StationElevationsByAttribute or ByIdx command)",BAD_DATA);
  }
  if(_aElevation==NULL) { 
    _aElevation=new double[_nNonZeroWeightedGridCells];
    for(int i=0;i<_nNonZeroWeightedGridCells;i++) { _aElevation[i]=0.0; }
  }
  if (_CellIDToIdx[CellID]!=DOESNT_EXIST){ // only assign elevations to cells used by model (i.e., ones with non-zero weights)
    _aElevation[_CellIDToIdx[CellID]]=elev;
  }
}

///////////////////////////////////////////////////////////////////
/// \brief checks if sum(_GridWeight[HRUID, :]) = 1.0 for all enabled HRUIDs
///
/// \param nHydroUnits [in] Number of HRUs must be given
///                         (is read from :GridWeights block in *.rvi file)
/// \param nGridCells  [in] Number of Grid Cells must be given
///                         (is read from :GridWeights block in *.rvi file)
/// \param pModel     [in] pointer to model
//
/// return true if sum for each enabled HRUID is one; otherwise false
//
bool   CForcingGrid::CheckWeightArray(const int nHydroUnits, const int nGridCells, const CModel *pModel)
{
  double sum_HRU;
  bool   check = true;

  if (_GridWeight != NULL){
     for(int k=0; k<_nHydroUnits; k++) {  // loop over HRUs
      sum_HRU = 0.0;
      for(int i=0; i<_nWeights[k]; i++) { // loop over all cells
        sum_HRU += _GridWeight[k][i];
      }
      if(fabs(sum_HRU - 1.0) < 0.05) {//repair if less than 5%
        for(int i=0; i<_nWeights[k];i++) {
          _GridWeight[k][i]/=sum_HRU;
        }
        sum_HRU = 1.0;
      }
      bool enabled=pModel->GetHydroUnit(k)->IsEnabled();
      if ((fabs(sum_HRU - 1.0) > 0.0001) && (enabled)) {
        cout<<"HRU ID = "<<pModel->GetHydroUnit(k)->GetID()<<" Sum Forcing Weights = "<<sum_HRU<<endl;
        for(int i=0; i<_nWeights[k]; i++) { cout<< _GridWeight[k][i]<<" ";} cout<<endl;
        check = false;
      }
    }
  }
  else{
    ExitGracefully(
      "CForcingGrid: CheckWeightArray: _GridWeight is not allocated yet. Call AllocateWeightArray(nHRUs) first.", BAD_DATA);
  }

  return check;
}

///////////////////////////////////////////////////////////////////
/// \brief weighting of HRU and CellID pair
/// return weighting of HRU and CellID pair (FUNCTION UNUSED)
//
double CForcingGrid::GetGridWeight(const int k,
                                   const int CellID) const
{
#ifdef _STRICTCHECK_
  int nCells   =GetRows()*GetCols(); //handles 3D or 2D
  if ((k <0)     || (k >=_nHydroUnits)){ExitGracefully("CForcingGrid::GetGridWeight: invalid k index",RUNTIME_ERR);}
  if ((CellID<0) || (CellID>=nCells  )){ExitGracefully("CForcingGrid::GetGridWeight: invalid CellID index",RUNTIME_ERR); }
  if (_GridWeight==NULL){ ExitGracefully("CForcingGrid::GetGridWeight: NULL Grid weight matrix",RUNTIME_ERR); }
#endif
  for(int i=0; i<_nWeights[k];i++) 
  {
    if(_GridWtCellIDs[k][i]==CellID) { return _GridWeight[k][i]; }
  } 
  return 0.0;
}

///////////////////////////////////////////////////////////////////
/// \brief Returns data interval
/// \return data interval (in days)
//
double  CForcingGrid::GetInterval() const{return _interval;}

///////////////////////////////////////////////////////////////////
/// \brief Returns _is_derived class variable.\n
///        If data are read from NetCDF (false) or derived from these data (true)
/// \return _is_derived class variable
//
bool   CForcingGrid::GetIsDerived() const{return _is_derived;}

///////////////////////////////////////////////////////////////////
/// \brief changes _is_derived member of ForcingGrid
//
void   CForcingGrid::SetIsDerived(const bool is_derived){_is_derived=is_derived;}

///////////////////////////////////////////////////////////////////
/// \brief Returns _is_3D class variable.\n
///        True if NetCDF data are (lat,lon,time), false if data are (nstations,time).
/// \return _is_derived class variable
//
bool   CForcingGrid::GetIs3D() const{return _is_3D;}

///////////////////////////////////////////////////////////////////
/// \brief Returns start year of time series data
/// \return Start year of time series data
//
int    CForcingGrid::GetStartYear() const{return _start_year;}

///////////////////////////////////////////////////////////////////
/// \brief Returns start day of time series data
/// \return Start day of time series data (Julian decimal date)
//
double CForcingGrid::GetStartDay() const{return _start_day;}

///////////////////////////////////////////////////////////////////
/// \brief Returns number of columns (= 1st dimension of gridded data)
/// \return Number of columns (= 1st dimension of gridded data)
//
int CForcingGrid::GetCols() const{return _GridDims[0];}

///////////////////////////////////////////////////////////////////
/// \brief Returns number of rows (= 2nd dimension of gridded data)
/// \return Number of rows (= 2nd dimension of gridded data)
//
int CForcingGrid::GetRows() const
{
  if (_is_3D) {return _GridDims[1];}
  else        {return 1;}
}

///////////////////////////////////////////////////////////////////
/// \brief Returns number of pulses (= 3rd dimension of gridded data)
/// \return Number of pulses (= 3rd dimension of gridded data)
//
int CForcingGrid::GetNumValues() const{return _nPulses;}

///////////////////////////////////////////////////////////////////
/// \brief Returns number of non-zero weighted grid cells
/// \return Number of non-zero weighted grid cells
//
int CForcingGrid::GetNumberNonZeroGridCells() const{return _nNonZeroWeightedGridCells;}

///////////////////////////////////////////////////////////////////
/// \brief Returns current chunk size
/// \return Current chunk size
//
int  CForcingGrid::GetChunkSize() const{return _ChunkSize;}

///////////////////////////////////////////////////////////////////
/// \brief Returns number of HRUs of class (_nHydroUnits)
/// \return number of HRUs of class
//
int CForcingGrid::GetnHydroUnits() const{return _nHydroUnits;}

///////////////////////////////////////////////////////////////////
/// \brief returns time index idx corresponding to t+tstep/2
/// \return time index idx corresponding to t+tstep/2
//
int CForcingGrid::GetTimeIndex(const double &t, const double &tstep) const
{
  return int((t+0.5*tstep) *rvn_round(1.0/_interval))  % _ChunkSize;
}

///////////////////////////////////////////////////////////////////
/// \brief returns weighted value of gridded forcing in HRU k over timestep starting at time t
/// \param k    [in] HRU index
/// \param t      [in] model time [days]
/// \return tstep [in] model time step [days]
/// \returns timestep weighted value of gridded forcing in HRU k
//
double CForcingGrid::GetWeightedValue(const int k,const double &t,const double &tstep) const
{
  int idx_new = GetTimeIndex(t,tstep);
  int nSteps = max(1,(int)(rvn_round(tstep/_interval)));//# of intervals in time step
  double wt,sum=0.0;
  for(int i = 0;i <_nWeights[k]; i++)
  {
    wt   = _GridWeight[k][i];
    sum += wt * GetValue_avg(_CellIDToIdx[_GridWtCellIDs[k][i]],idx_new,nSteps);
  }
  return sum;
}
///////////////////////////////////////////////////////////////////
/// \brief returns daily weighted value of gridded forcing in HRU k
/// \param k     [in] HRU index
/// \param t     [in] model time [days]
/// \param tstep [in] model time step [days]
/// \returns daily weighted value of gridded forcing in HRU k
//
double CForcingGrid::GetDailyWeightedValue(const int k,const double &t,const double &tstep, const optStruct &Options) const
{
  double time_shift=Options.julian_start_day-floor(Options.julian_start_day+TIME_CORRECTION);
  int it_new_day = GetTimeIndex(t-time_shift,tstep);//index corresponding to start of day
  double wt,sum=0;
  for(int i = 0;i <_nWeights[k]; i++)
  {
    wt   = _GridWeight[k][i];
    sum += wt * GetValue_avg(_CellIDToIdx[_GridWtCellIDs[k][i]],it_new_day,_steps_per_day);
  }
  return sum;
}
///////////////////////////////////////////////////////////////////
/// \brief returns daily weighted snowfrac value if this is a gridded snow dataset and pRain is provided
/// \param k     [in] HRU index
/// \param t     [in] model time [days]
/// \param pRain [in] pointer to gridded rain dataset
/// \param tstep [in] model time step [days]
//
double CForcingGrid::GetWeightedAverageSnowFrac(const int k,const double &t,const double &tstep,const CForcingGrid *pRain) const
{
#ifdef _STRICTCHECK_
  if ((k<0) || (k>_nHydroUnits)){ExitGracefully("CForcingGrid::GetWeightedAverageSnowFrac: invalid HRU index",RUNTIME_ERR); }
#endif 

  int nSteps = max(1,(int)(rvn_round(tstep/_interval)));//# of intervals in time step
  double wt,sum=0.0;
  double snow; double rain;
  for (int i=0;i<_nWeights[k];i++)
  {
    int ic=_CellIDToIdx[_GridWtCellIDs[k][i]];
    wt   = _GridWeight[k][i];
    snow = GetValue_avg(ic, t, nSteps);
    if(snow>0.0){
      rain=pRain->GetValue_avg(ic, t, nSteps);
      sum+= wt * snow/(snow+rain);
    }
  }
  return sum;
}

///////////////////////////////////////////////////////////////////
/// \brief returns weighted value of gridded reference elevation in HRU k
/// \param k    [in] HRU index
/// \returns reference elevation of forcing in HRU k
//
double CForcingGrid::GetRefElevation(const int k) const
{
  if(_aElevation==NULL) { return RAV_BLANK_DATA; }
  double wt,sum=0.0;
  for(int i = 0;i <_nWeights[k]; i++)
  {
    wt   = _GridWeight[k][i];
    sum += wt * _aElevation[_CellIDToIdx[_GridWtCellIDs[k][i]]];
  }
  return sum;
}
///////////////////////////////////////////////////////////////////
/// \brief returns latitude of cell l
/// \param l    [in] cell index
/// \returns latitude of cell l
//
double CForcingGrid::GetCellLatitude(const int l) const
{
  if(_aLatitude==NULL) { return 0.0; }
  if(l<0)              { return 0.0; }
  return _aLatitude[l];
}
///////////////////////////////////////////////////////////////////
/// \brief returns longitude of cell l
/// \param l    [in] cell index
/// \returns longitude of cell l
//
double CForcingGrid::GetCellLongitude(const int l) const
{
  if(_aLongitude==NULL) { return 0.0; }
  if(l<0) { return 0.0; }
  return _aLongitude[l];
}
///////////////////////////////////////////////////////////////////
/// \brief Returns magnitude of time series data point for which t is a float index
/// \param ic    [in] Index of grid cell with non-zero weighting (value between 0 and _nNonZeroWeightedGridCells)
/// \param t      [in] Time index
/// \return Magnitude of time series data point for which t is an index
//
double CForcingGrid::GetValue(const int ic, const int it) const
{
  return _aVal[it][ic];
}

///////////////////////////////////////////////////////////////////
/// \brief Returns average over n timesteps of time series data point for which t is an index
/// \param ic     [in] Index of grid cell with non-zero weighting (value between 0 and _nNonZeroWeightedGridCells)
/// \param t      [in] local time (with respect to chunk start, in days)
/// \param n      [in] Number of time steps
/// \return Magnitude of time series data point for which t is an index
//
double CForcingGrid::GetValue_avg(const int ic, const double &t, const int nsteps) const
{
  int it_start=max((int)(t),0);
  int lim=min(nsteps,_ChunkSize-it_start);
  double sum = 0.0;
  for (int it=it_start; it<it_start+lim;it++){
    sum += _aVal[it][ic];
  }
  sum /= (double)(lim);

  return sum;
}

///////////////////////////////////////////////////////////////////
/// \brief Returns minimum of n timesteps of time series data point for which t is an index
/// \param ic     [in] Index of grid cell with non-zero weighting (value between 0 and _nNonZeroWeightedGridCells)
/// \param t      [in] Time index
/// \param n      [in] Number of time steps
/// \return Magnitude of time series data point for which t is an index
//
double CForcingGrid::GetValue_min(const int ic, const double &t, const int nsteps) const
{
  double min_val = ALMOST_INF ;
  int it_start=max((int)(t),0);
  int lim=min(nsteps,_ChunkSize-it_start);
  for (int it=it_start; it<it_start+lim;it++){
    if(_aVal[it][ic] < min_val){min_val=_aVal[it][ic];}
  }
  return min_val;
}

///////////////////////////////////////////////////////////////////
/// \brief Returns maximum of n timesteps of time series data point for which t is an index
/// \param ic     [in] Index of grid cell with non-zero weighting (value between 0 and _nNonZeroWeightedGridCells)
/// \param t      [in] Time index
/// \param n      [in] Number of time steps
/// \return Magnitude of time series data point for which t is an index
//
double CForcingGrid::GetValue_max(const int ic, const double &t, const int nsteps) const
{
  double max_val = -ALMOST_INF ;
  int it_start=max((int)(t),0);
  int lim=min(nsteps,_ChunkSize-it_start);
  for (int it=it_start; it<it_start+lim;it++){
    if(_aVal[it][ic] > max_val){max_val=_aVal[it][ic];}
  }
  return max_val;
}

///////////////////////////////////////////////////////////////////
/// \brief Returns grid filename
//
string CForcingGrid::GetFilename() const {return _filename; }

///////////////////////////////////////////////////////////////////
/// \brief Returns snowfall correction factor
/// \param None
/// \return Snowfall correction factor
//
double CForcingGrid::GetSnowfallCorr() const {return _snowfall_corr; }

///////////////////////////////////////////////////////////////////
/// \brief Returns rainfall correction factor
/// \param None
/// \return Rainfall correction factor
//
double CForcingGrid::GetRainfallCorr() const {return _rainfall_corr; }

///////////////////////////////////////////////////////////////////
/// \brief Returns minimum temperature threshold used to determine cloud_cover factor
/// \param None
/// \return Minimum temperature threshold used to determine cloud_cover factor
//
double CForcingGrid::GetCloudMinRange() const {return _cloud_min_temp;}

///////////////////////////////////////////////////////////////////
/// \brief Returns maximum temperature threshold used to determine cloud_cover factor
/// \param None
/// \return Maximum temperature threshold used to determine cloud_cover factor
//
double CForcingGrid::GetCloudMaxRange() const {return _cloud_max_temp;}

//////////////////////////////////////////////////////////////////
/// \brief Returns representative average temperature over month
/// \param month [in] Month (from 1-12)
/// \return Representative Average temperature for month
//
double CForcingGrid::GetMonthlyAveTemp  (const int month) const
{
  return _aAveTemp[month-1];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns representative minimum temperature over month
/// \param month [in] Month (from 1-12)
/// \return Representative minimum temperature for month
//
double CForcingGrid::GetMonthlyMinTemp  (const int month) const
{
  return _aMinTemp[month-1];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns representative maximum temperature over month
/// \param month [in] Month (from 1-12)
/// \return Representative maximum temperature for month
//
double CForcingGrid::GetMonthlyMaxTemp  (const int month) const
{
  return _aMaxTemp[month-1];
}
//////////////////////////////////////////////////////////////////
/// \brief Returns average PET over month
/// \param month [in] Month for which the average PET is to be determined
/// \return Average PET over month
//
double CForcingGrid::GetMonthlyAvePET   (const int month) const
{
  return _aAvePET[month-1];
}

//////////////////////////////////////////////////////////////////
/// \brief Return hourly temperature correction at time t
/// \param t [in] Time at which daily temperature correction is to be determined
/// \return Daily temperature correction at time t
//
double CForcingGrid::DailyTempCorrection(const double t) const
{
  return -cos(2.0*PI*(t-PEAK_TEMP_HR/HR_PER_DAY));
}

///////////////////////////////////////////////////////////////////
/// \brief  Returns time index of currently read chunk corresponding to current model time step.
///         Returned value is double because time series could start at fractional day,
///         i.e. not at 00:00:00 hours.
///
/// \param  &Options          [in] Global model options information
/// \param  global_model_time [in] current simulation time stime
///
/// \return Index of current read chunk
///
double CForcingGrid::GetChunkIndexFromModelTimeStep(const optStruct &Options,
                                                    const double    global_model_time  // current model time step in [days]
                                                    ) const
{
  // overall un-chunked index in NetCDF file starting from very first entry in NetCDF file
  double idx_NC = (_t_corr + global_model_time) / _interval ;
  double delta  = (_t_corr + global_model_time) / _interval - floor(idx_NC); // difference from integer index

  // for example overall index = 17 and _ChunkSize = 5 [time points]  --> index in current chunk = 2
  int idx_chunk = (int)idx_NC % _ChunkSize;

  return (double)idx_chunk + delta;
}

///////////////////////////////////////////////////////////////////
/// \brief  Returns time index of currently read chunk corresponding to beginning of currently modelled day.
///         Returned value is double because time series could start at fractional day,
///         i.e. not at 00:00:00 hours.
///
/// \param  &Options          [in] Global model options information
/// \param  global_model_time [in] current simulation time stime
///
/// \return Index of current read chunk
///
double CForcingGrid::GetChunkIndexFromModelTimeStepDay(const optStruct &Options,
                                                       const double    global_model_time  // current model time step in [days]
                                                       ) const
{
  // overall un-chunked index in NetCDF file starting from very first entry in NetCDF file
  double idx_NC = (_t_corr + global_model_time) / _interval ;

  // for example overall index = 17 and _ChunkSize = 5 [time points]  --> index in current chunk = 2
  int idx_chunk = (int)idx_NC % _GridDims[2];

  // number of values per day
  int nn = (int)(1.0/_interval);

  // returned index is beginning of day not actual time step
  idx_chunk = idx_chunk - (idx_chunk % nn);

  return (double)idx_chunk; // + delta;
}

///////////////////////////////////////////////////////////////////
/// \brief  true if the variable should be deaccumulated
/// \return true if the variable should be deaccumulated
///
bool CForcingGrid::ShouldDeaccumulate()      const{return _deaccumulate;}

//////////////////////////////////////////////////////////////////
/// \brief Returns forcing type
/// \return forcing type
//
forcing_type CForcingGrid::GetForcingType  () const{return _ForcingType;}

void CForcingGrid::ReadAttGridFromNetCDF(const int ncid, const string varname, const int ncols,const int nrows,double *&values)
{
  // -------------------------------
  // Open NetCDF file, Get the lat long elev information
  // -------------------------------
#ifdef _RVNETCDF_

  if(varname!="NONE")
  {
    int       retval,varid;
    int       nCells,irow,icol;
    size_t    nc_start[2];
    size_t    nc_length[2];
    ptrdiff_t nc_stride[2];
    nc_start[0]  = 0;  nc_stride[0]=1;  nc_length[0]=ncols;
    nc_start[1]  = 0;  nc_stride[1]=1;  nc_length[1]=nrows;

    if(_is_3D) { nCells = ncols*nrows; }
    else       { nCells = ncols; }
    
    delete [] values; //re-read if buffered
    values=new double [_nNonZeroWeightedGridCells]; //allocate memory

    double *aVec=NULL;
    aVec=new double[nCells];//stores actual data
    for(int i=0; i<nCells; i++) { aVec[i]=NETCDF_BLANK_VALUE; }
    
    //get the varid for this attribute
    retval = nc_inq_varid(ncid,varname.c_str(),&varid); 
    if(retval==NC_ENOTVAR) {
      string warning="Variable \""+varname+"\" not found in NetCDF file "+_filename;
      ExitGracefully(warning.c_str(),BAD_DATA); retval=0;
    }
    HandleNetCDFErrors(retval);
    
    //get the data
    if(_is_3D) 
    {
      double  **aTmp2D=NULL; //stores pointers to rows/columns of 2D data
      aTmp2D = new double* [nrows];
      ExitGracefullyIf(aTmp2D==NULL,"CForcingGrid::ReadAttGridFromNetCDF",OUT_OF_MEMORY);
      for(int irow=0;irow<nrows;irow++) {
        aTmp2D[irow]=&aVec[irow*ncols]; //points to correct location in aVec data storage
      }

      retval=nc_get_vars_double(ncid,varid,nc_start,nc_length,nc_stride,&aTmp2D[0][0]);    HandleNetCDFErrors(retval);
      //copy matrix
      for(int ic=0; ic<_nNonZeroWeightedGridCells; ic++) {   // loop over non-zero weighted grid cells
        CellIdxToRowCol(_IdxNonZeroGridCells[ic],irow,icol);
        values[ic]=aTmp2D[irow][icol];
      }

      delete[] aTmp2D;
    }
    else 
    {
      retval=nc_get_vars_double(ncid,varid,nc_start,nc_length,nc_stride,&aVec[0]);   HandleNetCDFErrors(retval);
      //copy matrix
      for(int ic=0; ic<_nNonZeroWeightedGridCells; ic++) {   // loop over non-zero weighted grid cells
        CellIdxToRowCol(_IdxNonZeroGridCells[ic],irow,icol);
        values[ic]=aVec[icol];
      }
    }
    delete[] aVec;
    
    //TMP DEBUG - plot gridded elevations/lat/long - retain this code 
    /*bool found;
    double val;
    cout<<" #non-zero: "<<_nNonZeroWeightedGridCells<<endl;
    for (int j = 0; j < nrows; j++) {
      for (int i = 0; i < ncols; i++) {
        int l=j*ncols+i;
        found=false; val=-1;
        for(int ic=0; ic<_nNonZeroWeightedGridCells; ic++) { 
          if (_IdxNonZeroGridCells[ic]==l){ 
            found=true;
            val=values[ic];
          }
        }
        if (!found){l=-1;}
        //cout<<std::setw(3)<<l<<" "; //report cell ID
        cout<<std::setw(8)<<val<<" "; //report elevation/lat/long
      }
      cout<<endl;
    }*/
  }
#endif
}
// same as above but for string information
// Note this needs to be fixed by switching nrows/ncols as in above routine
/*void CForcingGrid::ReadAttGridFromNetCDF2(const int ncid,const string varname,const int nrows,const int ncols,string *values)
{
  // -------------------------------
  // Open NetCDF file, Get the lat long elev information
  // -------------------------------
#ifdef _RVNETCDF_

  if(varname!="NONE")
  {
    int       retval,varid;
    int       nCells,irow,icol;
    size_t    nc_start[2];
    size_t    nc_length[2];
    ptrdiff_t nc_stride[2];
    nc_start[0]  = 0; nc_stride[0]=1;
    nc_start[1]  = 0; nc_stride[1]=1;

    if(_is_3D) { nCells = nrows*ncols; }
    else { nCells = nrows; }

    delete[] values; //re-read if buffered
    values=new string[_nNonZeroWeightedGridCells]; //allocate memory

    double *aVec=NULL;
    aVec=new double[nCells];//stores actual data
    for(int i=0; i<nCells; i++) { aVec[i]=NETCDF_BLANK_VALUE; }

    //get the varid for this attribute
    retval = nc_inq_varid(ncid,varname.c_str(),&varid);
    if(retval==NC_ENOTVAR) {
      string warning="Variable "+varname+" not found in NetCDF file "+_filename;
      ExitGracefully(warning.c_str(),BAD_DATA); retval=0;
    }
    HandleNetCDFErrors(retval);

    //get the data
    if(_is_3D)
    {
      char  **aTmp2D=NULL; //stores pointers to rows/columns of 2D data
      aTmp2D=new char *[nrows];
      for(int row=0;row<nrows;row++) {
        aTmp2D[row]=&aVec[row*ncols]; //points to correct location in aVec data storage
      }

      retval=nc_get_vars_double(ncid,varid,nc_start,nc_length,nc_stride,&aTmp2D[0][0]);    HandleNetCDFErrors(retval);
      //copy matrix
      for(int ic=0; ic<_nNonZeroWeightedGridCells; ic++) {   // loop over non-zero weighted grid cells
        CellIdxToRowCol(_IdxNonZeroGridCells[ic],irow,icol);
        values[ic]=aTmp2D[irow][icol];
      }

      for(int i=0;i<nrows;i++) { delete[] aTmp2D[i]; } delete[] aTmp2D;
    }
    else
    {
      retval=nc_get_vars_text(ncid,varid,nc_start,nc_length,nc_stride,&aVec[0]);   HandleNetCDFErrors(retval);
      //copy matrix
      for(int ic=0; ic<_nNonZeroWeightedGridCells; ic++) {   // loop over non-zero weighted grid cells
        CellIdxToRowCol(_IdxNonZeroGridCells[ic],irow,icol);
        values[ic]=aVec[irow];
      }
    }
    delete[] aVec;
  }
#endif

}*/