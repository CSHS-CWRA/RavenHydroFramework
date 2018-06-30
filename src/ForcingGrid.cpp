/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2018 the Raven Development Team
  ----------------------------------------------------------------*/

// start of compilation if NetCDF library is available

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
CForcingGrid::CForcingGrid(string       ForcingType,
                           string       filename,
                           string       varname,
                           string       DimNames[3]
  )
{
  _ForcingType  = GetForcingTypeFromString(ForcingType);
  _filename     = filename;
  _varname      = varname;
  _DimNames[0]  = DimNames[0]; _DimNames[1]  = DimNames[1]; _DimNames[2]  = DimNames[2];
  _deaccumulate = false;

  // -------------------------------
  // Additional variables initialized and eventually overwritten by ParseTimeSeries
  // -------------------------------
  _rainfall_corr = 1.0;
  _snowfall_corr = 1.0;

  _cloud_min_temp=-20.0; //ensures cloud-free status always unless overriden
  _cloud_max_temp= 40.0;

  for (int i=0;i<12;i++)
  {
    _aAveTemp[i]=NOT_SPECIFIED;
    _aMinTemp[i]=NOT_SPECIFIED;
    _aMaxTemp[i]=NOT_SPECIFIED;
    _aAvePET [i]=NOT_SPECIFIED;
  }

  //initialized in ReallocateArraysInForcingGrid
  _aVal                = NULL;
  _GridWeight          = NULL;

  //initialized in SetIdxNonZeroGridCells()
  _IdxNonZeroGridCells = NULL; 
  _aFirstNonZeroWt     = NULL;
  _aLastNonZeroWt      = NULL;

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
  for (int ii=0; ii<3;  ii++) {_DimNames[ii] = grid._DimNames[ii]; }
  for (int ii=0; ii<3;  ii++) {_GridDims[ii] = grid._GridDims[ii]; }
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

  _aVal =  new double *[grid._ChunkSize];
  for (int it=0; it<grid._ChunkSize; it++) {                       // loop over time points in buffer
    _aVal[it]=NULL;
    _aVal[it] = new double [grid._nNonZeroWeightedGridCells];
    ExitGracefullyIf(_aVal[it]==NULL,"CForcingGrid::Constructor",OUT_OF_MEMORY);
    for (int ic=0; ic<grid._nNonZeroWeightedGridCells;ic++){       // loop over non-zero weighted cells
      _aVal[it][ic]=grid._aVal[it][ic];                            // copy the value
    }
  }

  _GridWeight =  new double *[grid._nHydroUnits];
  for (int ik=0; ik<grid._nHydroUnits; ik++) {                           // loop over HRUs
    _GridWeight[ik] =NULL;
    _GridWeight[ik] = new double [grid._GridDims[0]*grid._GridDims[1]];
    ExitGracefullyIf(_GridWeight[ik]==NULL,"CForcingGrid::Constructor(2)",OUT_OF_MEMORY);
    for (int ic=0; ic<grid._GridDims[0]*grid._GridDims[1]; ic++) {       // loop over cells = rows*cols
      _GridWeight[ik][ic]=grid._GridWeight[ik][ic];                      // copy the value
    }
  }

  _IdxNonZeroGridCells = NULL; 
  _IdxNonZeroGridCells = new int [grid._nNonZeroWeightedGridCells];
  for (int ic=0; ic<grid._nNonZeroWeightedGridCells; ic++) {             // loop over non-zero weighted cells
    _IdxNonZeroGridCells[ic]=grid._IdxNonZeroGridCells[ic];              // copy the value
  }

  _aFirstNonZeroWt=new int [grid._nHydroUnits];
  _aLastNonZeroWt =new int [grid._nHydroUnits];
  for(int k=0; k<_nHydroUnits;k++){
    _aFirstNonZeroWt[k]=grid._aFirstNonZeroWt[k];
    _aLastNonZeroWt [k]=grid._aLastNonZeroWt[k];
  }
}

///////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
//
CForcingGrid::~CForcingGrid()
{
  if (DESTRUCTOR_DEBUG){cout<<"    DELETING GRIDDED DATA"<<endl;}
  for(int it=0; it<_ChunkSize; it++) {delete[] _aVal[it]; _aVal[it]=NULL;} delete[] _aVal;_aVal= NULL;
  for(int k=0; k<_nHydroUnits; k++) { delete[] _GridWeight[k];  _GridWeight[k]=NULL;}    delete[] _GridWeight;_GridWeight= NULL;   
  delete [] _GridWeight;            _GridWeight          = NULL;
  delete [] _IdxNonZeroGridCells;   _IdxNonZeroGridCells = NULL;
  delete [] _aFirstNonZeroWt;       _aFirstNonZeroWt     = NULL;
  delete [] _aLastNonZeroWt;        _aLastNonZeroWt     = NULL;
}


///////////////////////////////////////////////////////////////////
/// \brief Implementation of forcing grid constructor only determining grid
///        dimensions and buffersize (data are only initalized but not read from file)
///
/// \note  Needs _ForcingType, _filename, _varname, _DimNames to be set already.
///        Use for that either "SetForcingType", "SetFilename", "SetVarname", and "SetDimNames" \n
///        or "CForcingGrid()".
void CForcingGrid::ForcingGridInit( const optStruct   &Options )
{
#ifdef _RVNETCDF_
  int    ncid;                  // file unit
  int    dimid_x;               // id of x dimension (columns)
  int    dimid_y;               // id of y dimension (rows)
  int    dimid_t;               // id of t dimension (time)
  int    dimids_var[3];         // ids of dimensions of a NetCDF variable
  int    varid_t;               // id of time variable
  int    varid_f;               // id of forcing variable read
  int    retval;                // error value for NetCDF routines
  size_t GridDim_t;             // special type for GridDims required by nc routine
  char   unit_t[200];           // special type for string of variable's unit required by nc routine
  char * unit_f;                // special type for string of variable's unit required by nc routine
  char * long_name_f;           // special type for string of variable's long name required by nc routine
  int    BytesPerTimestep;      // Memory requirement for one timestep of gridded forcing file [Bytes]

  int iatt;
  char attrib_name[500];

  int retval1;
  size_t att_len;     // length of the attribute's text

  string dash;           // to check format of time unit string
  string colon;          // to check format of time unit string
  string unit_t_str;     // to check format of time unit string

  // open NetCDF read-only; ncid will be set
  // _filename.c_str() converts filename from 'string' to 'const char *'
  if(Options.noisy){ cout<<"Initializing grid file "<<_filename<<endl; }
  retval = nc_open(_filename.c_str(), NC_NOWRITE, &ncid);
  HandleNetCDFErrors(retval);


  // Get the id of dimensions based on its name; dimid will be set
  // Find length of dimension and store it in GridDim
  // dimension x = number of columns of the grid
  retval = nc_inq_dimid (ncid, _DimNames[0].c_str(), &dimid_x);  HandleNetCDFErrors(retval);
  retval = nc_inq_dimlen(ncid, dimid_x, &GridDim_t);             HandleNetCDFErrors(retval);
  _GridDims[0] = static_cast<int>(GridDim_t);  // convert returned 'size_t' to 'int'
  
  // dimension y = number of rows of the grid
  retval = nc_inq_dimid (ncid, _DimNames[1].c_str(), &dimid_y);  HandleNetCDFErrors(retval);
  retval = nc_inq_dimlen(ncid, dimid_y, &GridDim_t);             HandleNetCDFErrors(retval);
  _GridDims[1] = static_cast<int>(GridDim_t);  // convert returned 'size_t' to 'int'

  // dimension t = number of time steps in NetCDF file
  retval = nc_inq_dimid(ncid, _DimNames[2].c_str(), &dimid_t);   HandleNetCDFErrors(retval);
  retval = nc_inq_dimlen(ncid, dimid_t, &GridDim_t);             HandleNetCDFErrors(retval);
  _GridDims[2] = static_cast<int>(GridDim_t);  // convert returned 'size_t' to 'int'

  // Get the id of the data variable based on its name; varid will be set
  retval = nc_inq_varid(ncid, _varname.c_str(), &varid_f);       HandleNetCDFErrors(retval);

  // determine in which order the dimensions are in variable
  //   if  (t,         y=lat=row, x=lon=col)   --> (2,1,0) --> (dimid_t, dimid_y, dimid_x)
  //   if  (t,         x=lon=col, y=lat=row)   --> (2,0,1) --> (dimid_t, dimid_x, dimid_y)
  //   if  (y=lat=row, t,         x=lon=col)   --> (1,2,0) --> (dimid_y, dimid_t, dimid_x)
  //   if  (y=lat=row, x=lon=col, t        )   --> (1,0,2) --> (dimid_y, dimid_x, dimid_t)
  //   if  (x=lon=col, t,         y=lat=row)   --> (0,2,1) --> (dimid_x, dimid_t, dimid_y)
  //   if  (x=lon=col, y=lat=row, t        )   --> (0,1,2) --> (dimid_x, dimid_y, dimid_t)
  retval = nc_inq_vardimid(ncid, varid_f, dimids_var);          HandleNetCDFErrors(retval);

  if    ((dimids_var[0]==dimid_x) && (dimids_var[1]==dimid_y) && (dimids_var[2]==dimid_t))  {  // dimensions are (x,y,t)
    _dim_order = 1;
  }
  else if((dimids_var[0]==dimid_y) && (dimids_var[1]==dimid_x) && (dimids_var[2]==dimid_t))  {  // dimensions are (y,x,t)
    _dim_order = 2;
  }
  else if((dimids_var[0] == dimid_x) && (dimids_var[1] == dimid_t) && (dimids_var[2] == dimid_y))  {  // dimensions are (x,t,y)
    _dim_order = 3;
  }
  else if((dimids_var[0] == dimid_t) && (dimids_var[1] == dimid_x) && (dimids_var[2] == dimid_y))  {  // dimensions are (t,x,y)
    _dim_order = 4;
  }
  else if((dimids_var[0] == dimid_y) && (dimids_var[1] == dimid_t) && (dimids_var[2] == dimid_x))  {  // dimensions are (y,t,x)
    _dim_order = 5;
  }
  else if((dimids_var[0] == dimid_t) && (dimids_var[1] == dimid_y) && (dimids_var[2] == dimid_x))  {  // dimensions are (t,y,x)
    _dim_order = 6;
  }
  if(Options.noisy){ printf("  Order of dimensions in NetCDF is Case %i...\n",_dim_order); }

  // -------------------------------
  // Start day and year are extracted from the UNITS attribute of time axis
  // can be, e.g., "days    since 1989-01-01 00:00:00" or
  //               "hours   since 1989-01-01 00:00:00" or
  //               "minutes since 1989-01-01 00:00:00" or
  //               "seconds since 1989-01-01 00:00:00" 
  // then whole time variable has to be read and first time step needs to be added to get _start_year and _start_day
  // difference between t[0] and t[1] is then used to get _interval
  // -------------------------------
  retval = nc_inq_varid(ncid,_DimNames[2].c_str(),&varid_t);// varid of time
  HandleNetCDFErrors(retval);

  retval = nc_get_att_text(ncid,varid_t,"units",unit_t);    // unit of time
  HandleNetCDFErrors(retval);
  

  int *time=NULL;
  time=new int [_GridDims[2]];
  ExitGracefullyIf(time==NULL,"CForcingGrid::ForcingGridInit",OUT_OF_MEMORY);
  
  retval = nc_get_var_int(ncid,varid_t,&time[0]);        
  HandleNetCDFErrors(retval);

  // -------------------------------
  // check if we have equal time steps
  // -------------------------------
  double delta_t = time[1] - time[0];
  for (int ii=1; ii<_GridDims[2]-1 ; ii++)
  {
    double delta_t2 = time[ii+1] - time[ii];
    if ( abs(delta_t2 - delta_t) > 0.00001 )
    {
      printf("\n\n\nCForcingGrid: ForcingGridInit: variable name: %s\n",_varname.c_str());
      ExitGracefully("CForcingGrid: ForcingGridInit: time steps are not equal in gridded input",BAD_DATA);
    }
  }

  // -------------------------------
  // delta t (in days)
  // -------------------------------
  if (strstr(unit_t, "hours")) {  // if unit of time in hours: hours since 1989-01-01 00:00:00

    // ---------------------------
    // check if format is YYYY-MM-DD HH:MM:SS
    // fill with leading zeros if necessary
    // ---------------------------
    unit_t_str = to_string(unit_t);
    dash = unit_t_str.substr(16, 1);  // first dash in date
    if ( !strstr(dash.c_str(), "-") ){
      printf("time unit string: %s\n",unit_t_str.c_str());
      ExitGracefully("CForcingGrid: ForcingGridInit: time unit string has weird format!",BAD_DATA);
    }
    dash  = unit_t_str.substr(19, 1);  // second dash in date
    if ( !strstr(dash.c_str(), "-") ){
      unit_t_str.insert(17,"0");
      //if (Options.noisy){ printf("corrected time unit string: %s\n",unit_t_str.c_str()); }
    }
    colon = unit_t_str.substr(25, 1);  // first colon in time
    if ( !strstr(colon.c_str(), ":") ){
      unit_t_str.insert(23,"0");
      //if (Options.noisy){ printf("corrected time unit string: %s\n",unit_t_str.c_str()); }
    }
    colon = unit_t_str.substr(28, 1);  // second colon in time
    if ( !strstr(colon.c_str(), ":") ){
      unit_t_str.insert(26,"0");
      //if (Options.noisy){ printf("corrected time unit string: %s\n",unit_t_str.c_str()); }
    }

    // ---------------------------
    // set class variables
    // ---------------------------
    string sDate = unit_t_str.substr(12,10);
    string sTime = unit_t_str.substr(23,8);
    time_struct tt = DateStringToTimeStruct(sDate, sTime);
    if (time[0]*_interval > 365.) {
      ExitGracefullyIf(_interval<=0,
                       "CForcingGrid: ForcingGridInit: first time step in NetCDF is more than one year apart from date written in time unit!",BAD_DATA);
    }
    _interval   = (time[1] - time[0])/24.;
    _start_day  = tt.julian_day + time[0]*(1./24.);  // I hope this is not more than one year
    _start_year = tt.year;

  /* printf("ForcingGrid: unit_t:          %s\n",unit_t_str.c_str());
     printf("ForcingGrid: sDate:           %s\n",sDate.c_str());
     printf("ForcingGrid: sTime:           %s\n",sTime.c_str());
     printf("ForcingGrid: tt.julian_day:   %f\n",tt.julian_day);
     printf("ForcingGrid: tt.day_of_month: %i\n",tt.day_of_month);
     printf("ForcingGrid: tt.month:        %i\n",tt.month);
     printf("ForcingGrid: tt.year:         %i\n",tt.year);
     printf("ForcingGrid: time[0]:         %i\n",time[0]);
     printf("ForcingGrid: _interval:       %f\n",_interval);*/

  }
  else 
  {
    if (strstr(unit_t, "days")) {  // if unit of time in days: days since 1989-01-01 00:00:00
      
      // ---------------------------
      // check if format is YYYY-MM-DD HH:MM:SS
      // fill with leading zeros if necessary
      // ---------------------------
      unit_t_str = to_string(unit_t);
      dash = unit_t_str.substr(15, 1);  // first dash in date
      if ( !strstr(dash.c_str(), "-") ){
        printf("time unit string: %s\n",unit_t_str.c_str());
        ExitGracefully("CForcingGrid: ForcingGridInit: time unit string has weird format!",BAD_DATA);
      }
      dash  = unit_t_str.substr(18, 1);  // second dash in date
      if ( !strstr(dash.c_str(), "-") ){
        unit_t_str.insert(16,"0");
        //if (Options.noisy){ printf("corrected time unit string: %s\n",unit_t_str.c_str()); }
      }
      colon = unit_t_str.substr(24, 1);  // first colon in time
      if ( !strstr(colon.c_str(), ":") ){
        unit_t_str.insert(22,"0");
        //if (Options.noisy){ printf("corrected time unit string: %s\n",unit_t_str.c_str()); }
      }
      colon = unit_t_str.substr(27, 1);  // second colon in time
      if ( !strstr(colon.c_str(), ":") ){
        unit_t_str.insert(25,"0");
        //if (Options.noisy){ printf("corrected time unit string: %s\n",unit_t_str.c_str()); }
      }

      // ---------------------------
      // set class variables
      // ---------------------------
      string sDate = unit_t_str.substr(11,10);
      string sTime = unit_t_str.substr(22,8);
      time_struct tt = DateStringToTimeStruct(sDate, sTime);
      if (time[0]*_interval > 365.) {
        ExitGracefullyIf(_interval<=0,
                         "CForcingGrid: ForcingGridInit: first time step in NetCDF is more than one year apart from date written in time unit!",BAD_DATA);
      }
      _interval   = (time[1] - time[0])/1.;
      _start_day  = tt.julian_day + time[0]*(1.);  // I hope this is not more than one year
      _start_year = tt.year; 
    }
    else {
      if (strstr(unit_t, "minutes")) {  // if unit of time in minutes: minutes since 1989-01-01 00:00:00
	  
	  // ---------------------------
	  // check if format is YYYY-MM-DD HH:MM:SS
	  // fill with leading zeros if necessary
	  // ---------------------------
	  unit_t_str = to_string(unit_t);
	  dash = unit_t_str.substr(18, 1);  // first dash in date
	  if ( !strstr(dash.c_str(), "-") ){
	    printf("time unit string: %s\n",unit_t_str.c_str());
	    ExitGracefully("CForcingGrid: ForcingGridInit: time unit string has weird format!",BAD_DATA);
	  }
	  dash  = unit_t_str.substr(21, 1);  // second dash in date
	  if ( !strstr(dash.c_str(), "-") ){
	    unit_t_str.insert(19,"0");
	    //if (Options.noisy){ printf("corrected time unit string: %s\n",unit_t_str.c_str()); }
	  }
	  colon = unit_t_str.substr(27, 1);  // first colon in time
	  if ( !strstr(colon.c_str(), ":") ){
	    unit_t_str.insert(25,"0");
	    //if (Options.noisy){ printf("corrected time unit string: %s\n",unit_t_str.c_str()); }
	  }
	  colon = unit_t_str.substr(30, 1);  // second colon in time
	  if ( !strstr(colon.c_str(), ":") ){
	    unit_t_str.insert(28,"0");
	    //if (Options.noisy){ printf("corrected time unit string: %s\n",unit_t_str.c_str()); }
	  }
	  
	  // ---------------------------
	  // set class variables
	  // ---------------------------
	  string sDate = unit_t_str.substr(14,10);
	  string sTime = unit_t_str.substr(25,8);
	  time_struct tt = DateStringToTimeStruct(sDate, sTime);
	  if (time[0]*_interval > 365.) {
	    ExitGracefullyIf(_interval<=0,
			     "CForcingGrid: ForcingGridInit: first time step in NetCDF is more than one year apart from date written in time unit!",BAD_DATA);
	  }
	  _interval   = (time[1] - time[0])/24./60.;
	  _start_day  = tt.julian_day + time[0]*(1./24./60.);  // I hope this is not more than one year
	  _start_year = tt.year;
	  
	  /* printf("ForcingGrid: unit_t:          %s\n",unit_t_str.c_str());
	     printf("ForcingGrid: sDate:           %s\n",sDate.c_str());
	     printf("ForcingGrid: sTime:           %s\n",sTime.c_str());
	     printf("ForcingGrid: tt.julian_day:   %f\n",tt.julian_day);
	     printf("ForcingGrid: tt.day_of_month: %i\n",tt.day_of_month);
	     printf("ForcingGrid: tt.month:        %i\n",tt.month);
	     printf("ForcingGrid: tt.year:         %i\n",tt.year);
	     printf("ForcingGrid: time[0]:         %i\n",time[0]);
	     printf("ForcingGrid: _interval:       %f\n",_interval);*/
	  
      }
      else {
	if (strstr(unit_t, "seconds")) {  // if unit of time in seconds: seconds since 1989-01-01 00:00:00
	  
	  // ---------------------------
	  // check if format is YYYY-MM-DD HH:MM:SS
	  // fill with leading zeros if necessary
	  // ---------------------------
	  unit_t_str = to_string(unit_t);
	  dash = unit_t_str.substr(18, 1);  // first dash in date
	  if ( !strstr(dash.c_str(), "-") ){
	    printf("time unit string: %s\n",unit_t_str.c_str());
	    ExitGracefully("CForcingGrid: ForcingGridInit: time unit string has weird format!",BAD_DATA);
	  }
	  dash  = unit_t_str.substr(21, 1);  // second dash in date
	  if ( !strstr(dash.c_str(), "-") ){
	    unit_t_str.insert(19,"0");
	    //if (Options.noisy){ printf("corrected time unit string: %s\n",unit_t_str.c_str()); }
	  }
	  colon = unit_t_str.substr(27, 1);  // first colon in time
	  if ( !strstr(colon.c_str(), ":") ){
	    unit_t_str.insert(25,"0");
	    //if (Options.noisy){ printf("corrected time unit string: %s\n",unit_t_str.c_str()); }
	  }
	  colon = unit_t_str.substr(30, 1);  // second colon in time
	  if ( !strstr(colon.c_str(), ":") ){
	    unit_t_str.insert(28,"0");
	    //if (Options.noisy){ printf("corrected time unit string: %s\n",unit_t_str.c_str()); }
	  }
	  
	  // ---------------------------
	  // set class variables
	  // ---------------------------
	  string sDate = unit_t_str.substr(14,10);
	  string sTime = unit_t_str.substr(25,8);
	  time_struct tt = DateStringToTimeStruct(sDate, sTime);
	  if (time[0]*_interval > 365.) {
	    ExitGracefullyIf(_interval<=0,
			     "CForcingGrid: ForcingGridInit: first time step in NetCDF is more than one year apart from date written in time unit!",BAD_DATA);
	  }
	  _interval   = (time[1] - time[0])/24./60./60.;
	  _start_day  = tt.julian_day + time[0]*(1./24./60./60.);  // I hope this is not more than one year
	  _start_year = tt.year;
	  
	  /* printf("ForcingGrid: unit_t:          %s\n",unit_t_str.c_str());
	     printf("ForcingGrid: sDate:           %s\n",sDate.c_str());
	     printf("ForcingGrid: sTime:           %s\n",sTime.c_str());
	     printf("ForcingGrid: tt.julian_day:   %f\n",tt.julian_day);
	     printf("ForcingGrid: tt.day_of_month: %i\n",tt.day_of_month);
	     printf("ForcingGrid: tt.month:        %i\n",tt.month);
	     printf("ForcingGrid: tt.year:         %i\n",tt.year);
	     printf("ForcingGrid: time[0]:         %i\n",time[0]);
	     printf("ForcingGrid: _interval:       %f\n",_interval);*/
	  
	}
	ExitGracefullyIf(_interval<=0,
			 "CForcingGrid: ForcingGridInit: this unit in time is not implemented yet (only days, hours, minutes, seconds)",BAD_DATA);
      }
    }
  }
  ExitGracefullyIf(_interval<=0,
                   "CForcingGrid: ForcingGridInit: negative time interval is not allowed",BAD_DATA);

  if ( ForcingToString(_ForcingType) == "TEMP_DAILY_AVE" && _interval != 1.0 )
    ExitGracefully( "Gridded forcing 'TEMP_DAILY_AVE' must have daily timestep. Please use 'TEMP_AVE' instead (see *.rvt file).",BAD_DATA);

  delete[] time;

  // -------------------------------
  // Determine chunk size (_ChunkSize),
  //           number of chunks read in total (_nChunk) and
  // Set       id of current chunk (_iChunk)
  // depending on size of grid (#cols x #rows x buffersize x 8byte <=  10 MB=10*1024*1024 byte)
  // -------------------------------
  BytesPerTimestep = 8 * _GridDims[0] * _GridDims[1];
  _ChunkSize = int(max(min( (10*1024 * 1024) / BytesPerTimestep, _GridDims[2]),1));     // number of timesteps per chunk
  _ChunkSize = max(int(round(1./_interval)),int(int(_ChunkSize*_interval)/_interval));  // make sure complete days and at least one day is read
  _nChunk    = int(ceil(1.*_GridDims[2] / _ChunkSize));                                 // total number of chunks
  _iChunk    = -1;                                                                      // current chunk read (-1 = no chunk read, 0 = first chunk...)

  // -------------------------------
  // set _tag
  //     ---> looks for attributes "long_name" and "units" of forcing variable
  //     --> if either attribute not available it will be set to "?"
  // -------------------------------
  retval = 0;
  iatt   = 0;
  long_name_f = (char *) malloc(2);
  strcpy(long_name_f, "?\0");
  unit_f = (char *) malloc(2);
  strcpy(unit_f, "?\0");

  while (retval==0) 
  {
    // inquire attributes name
    retval = nc_inq_attname(ncid, varid_f, iatt, attrib_name);

    // long_name of forcing
    if (strcmp(attrib_name,"long_name") == 0) 
    {
      retval1 = nc_inq_attlen (ncid, varid_f, attrib_name, &att_len);// inquire length of attribute's text
      HandleNetCDFErrors(retval1);
      long_name_f = (char *) malloc(att_len + 1);// allocate memory of char * to hold attribute's text
      retval1 = nc_get_att_text(ncid, varid_f, attrib_name, long_name_f);// read attribute text
      HandleNetCDFErrors(retval1);
      long_name_f[att_len] = '\0';// add string determining character
    }

    // unit of forcing
    if (strcmp(attrib_name,"units") == 0) 
    {
      retval1 = nc_inq_attlen (ncid, varid_f, attrib_name, &att_len);// inquire length of attribute's text
      HandleNetCDFErrors(retval1);
      unit_f = (char *) malloc(att_len + 1);// allocate memory of char * to hold attribute's text
      retval1 = nc_get_att_text(ncid, varid_f, attrib_name, unit_f);// read attribute text
      HandleNetCDFErrors(retval1);
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
  _nPulses      = _GridDims[2];
  _pulse        = true;
  ExitGracefullyIf(_nPulses<=0,
                   "CForcingGrid: ForcingGridInit: no time point entries in forcing grid",BAD_DATA);

  // -------------------------------
  // Close the file. This frees up any internal NetCDF resources
  // associated with the file, and flushes any buffers.
  // -------------------------------
  retval = nc_close(ncid);
  HandleNetCDFErrors(retval);

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
/// \brief  Reallocate all arrays in class to (potentially updated) grid dimensions
//          mainly used when sub-daily grids have to be added to model
//
void CForcingGrid::ReallocateArraysInForcingGrid( )
{

  // -------------------------------
  // Initialize data array and set all entries to NODATA value
  // -------------------------------
  _aVal = NULL;
  _aVal =  new double *[_GridDims[2]];
  for (int it=0; it<_GridDims[2]; it++) {                       // loop over time points in buffer
    _aVal[it]=NULL;
    _aVal[it] = new double [_nNonZeroWeightedGridCells];
    ExitGracefullyIf(_aVal[it]==NULL,"CForcingGrid::ReallocateArraysInForcingGrid",OUT_OF_MEMORY);
    for (int ic=0; ic<_nNonZeroWeightedGridCells;ic++){       // loop over non-zero weighted cells
      _aVal[it][ic]=NETCDF_BLANK_VALUE;                       // initialize
    }
  }

  // -------------------------------
  // Initialize weight array and set all entries to NODATA value
  // -------------------------------
  _GridWeight =  NULL;
  _GridWeight =  new double *[_nHydroUnits];
  for (int ik=0; ik<_nHydroUnits; ik++) {                      // loop over HRUs
    _GridWeight[ik] =NULL;
    _GridWeight[ik] = new double [_GridDims[0]*_GridDims[1]];
    ExitGracefullyIf(_GridWeight[ik]==NULL,"CForcingGrid::ReallocateArraysInForcingGrid(2)",OUT_OF_MEMORY);
    for (int ic=0; ic<_GridDims[0]*_GridDims[1]; ic++) {       // loop over cells = rows*cols
      _GridWeight[ik][ic]=NETCDF_BLANK_VALUE;
    }
  }

  // -------------------------------
  // Initialize indexes of non-zero weighted cell ids and set all entries to NODATA value
  // -------------------------------
  _IdxNonZeroGridCells = NULL;
  _IdxNonZeroGridCells =  new int [_nNonZeroWeightedGridCells];
  for (int ic=0; ic<_nNonZeroWeightedGridCells; ic++) { // loop over all cells non-zero weighted grid cells
    _IdxNonZeroGridCells[ic] = -1;
  }

  //initialize _aFirstNonZeroWt and _aLastNonZeroWt
  _aFirstNonZeroWt=new int [_nHydroUnits];
  _aLastNonZeroWt =new int [_nHydroUnits];
  for(int k=0; k<_nHydroUnits;k++){
    _aFirstNonZeroWt[k]=0;
    _aLastNonZeroWt [k]=_nNonZeroWeightedGridCells-1;
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
  int    ncid;          // file unit
  int    dim1;          // length of 1st dimension in NetCDF data
  int    dim2;          // length of 2nd dimension in NetCDF data
  int    dim3;          // length of 3rd dimension in NetCDF data

  int    varid_f;       // id of forcing variable read
  int    retval;        // error value for NetCDF routines
  int    iChunkSize;    // size of current chunk; always equal _ChunkSize except for last chunk in file (might be shorter)
  int    iChunk_new;    // chunk in which current model time step falls

  double model_start_day;
  int    model_start_yr;
  double model_duration;
  double model_timestep;

  new_chunk_read = false;

  // -------------------------------
  // check if chunk id is valid
  // -------------------------------
  ExitGracefullyIf((_iChunk>=_nChunk) || (_iChunk <-1),
                   "CForcingGrid: ReadData: this is not a valid chunk",BAD_DATA);

  // -------------------------------
  // only initialization required
  // -------------------------------
  if (_iChunk == -1) 
  {
    // -------------------------------
    // allocate array: usually array is of size [_ChunkSize, nRows, nCols]
    //                 except for last chunk it is of size [_GridDims[2]-(nchunks-1)*_ChunkSize+1, nRows, nCols]
    // -------------------------------
    if (_iChunk < _nChunk-1) {
      iChunkSize = _ChunkSize;
    }
    else {
      //
      // t1 = Options.duration - global_model_time       --> days still to simulate
      // t2 = t1 / _interval                             --> number of timesteps need to be read from forcing
      iChunkSize = int((Options.duration - global_model_time) / _interval);
    }

    if (Options.noisy){ printf("Initialize before reading first chunk  (var = %s)\n",_varname.c_str()); }
    if (Options.noisy){ printf("  Dim of chunk read: nCol = %i   nRow = %i   nTime = %i\n",_GridDims[0],_GridDims[1],iChunkSize); }

    // (1) set the correction time _t_corr, i.e. days between forcing data start and model simulation start
    // (2) check if forcing data cover whole simulation period
    model_start_day = Options.julian_start_day;
    model_start_yr  = Options.julian_start_year;
    model_duration  = Options.duration;  // model duration in days --> number of timesteps must be duration/timestep
    model_timestep  = Options.timestep;

    Initialize(model_start_day,model_start_yr,model_duration,model_timestep,Options);

    // -------------------------------
    // allocate _aVal matrix
    // -------------------------------
    _aVal = NULL;
    _aVal = new double *[iChunkSize];
    for (int it=0; it<iChunkSize; it++) {                    // loop over time points in buffer
      _aVal[it]=NULL;
      _aVal[it] = new double [_nNonZeroWeightedGridCells];
      ExitGracefullyIf(_aVal[it]==NULL,"CForcingGrid::ReadData",OUT_OF_MEMORY);
      for (int ic=0; ic<_nNonZeroWeightedGridCells;ic++){    // loop over all non-zero weighted grid cells
        _aVal[it][ic]=NETCDF_BLANK_VALUE;
      }
    }

    // -------------------------------
    // set _is_derived_data to False because data are truely read from a file
    // -------------------------------
    _is_derived = false;



    // -------------------------------
    // Print data on screen
    // -------------------------------
    // for (int it=0; it<4; it++)  // York.nc : loop over time points in buffer: 4 -->
    //   {
    //     printf("t = %i (Init)\n",it);
    //     printf("   %f %f %f %f \n",_aVal[it][0][0],_aVal[it][0][1],_aVal[it][0][2],_aVal[it][0][3]);
    //     printf("   %f %f %f %f \n",_aVal[it][1][0],_aVal[it][1][1],_aVal[it][1][2],_aVal[it][1][3]);
    //     printf("   %f %f %f %f \n",_aVal[it][2][0],_aVal[it][2][1],_aVal[it][2][2],_aVal[it][2][3]);
    //     printf("   %f %f %f %f \n",_aVal[it][3][0],_aVal[it][3][1],_aVal[it][3][2],_aVal[it][3][3]);
    //     printf("   %f %f %f %f \n",_aVal[it][4][0],_aVal[it][4][1],_aVal[it][4][2],_aVal[it][4][3]);
    //     printf("   %f %f %f %f \n",_aVal[it][5][0],_aVal[it][5][1],_aVal[it][5][2],_aVal[it][5][3]);
    //     printf("\n");
    //   }

    if (Options.noisy){ printf("  Finished initialization before reading\n"); }
    if (Options.noisy){ printf("  \n"); }

  } // That's all to do if chunk == -1

  // ------------------------------------------------------------------
  // Now read first proper chunk (if possible)
  // ------------------------------------------------------------------

  // -------------------------------
  // check if given model time step is covered by current chunk
  // if yes, do nothing
  // if no,  read next chunk
  // -------------------------------
  model_start_day = Options.julian_start_day;
  model_start_yr  = Options.julian_start_year;
  model_duration  = Options.duration;
  model_timestep  = Options.timestep;

  double length_chunk = _interval * _ChunkSize; // [days]
  iChunk_new = int(floor((_t_corr + global_model_time) / length_chunk ));

  if(_iChunk != iChunk_new) 
  {  // current model time step is not covered by current chunk --> read new chunk

    if(Options.noisy){ printf("\n"); }
    if(Options.noisy){ printf("Start reading new chunk... iChunk = %i (var = %s)\n",iChunk_new,_varname.c_str()); }
    _iChunk = iChunk_new;

    // -------------------------------
    // check if chunk id is valid
    // -------------------------------
    ExitGracefullyIf((_iChunk>=_nChunk) || (_iChunk<-1),
      "CForcingGrid: ReadData: this is not a valid chunk",BAD_DATA);

    // -------------------------------
    // allocate array: usually array is of size [_ChunkSize, nRows, nCols]
    //                 except for last chunk it is of size [_GridDims[2]-(nchunks-1)*_ChunkSize+1, nRows, nCols]
    // -------------------------------
    if(_iChunk < _nChunk-1) {
      iChunkSize = _ChunkSize;
    }
    else 
    {
      // t1 = _t_corr + Options.duration - global_model_time       --> days still to simulate (plus days at the beginning)
      // t2 = t1 / _interval                                       --> number of timesteps need to be read from forcing
      if(Options.noisy){ printf("   >>> t1 = %f\n",Options.duration - global_model_time); }
      if(Options.noisy){ printf("   >>> t2 = %i\n",int((Options.duration - global_model_time) / _interval)); }
      iChunkSize = int((_t_corr + Options.duration - global_model_time) / _interval);
    }
    if(Options.noisy){ printf("  iChunksize: %i\n",iChunkSize); }
    
    // -------------------------------
    // Open NetCDF file, Get the id of the forcing data, varid_f
    // -------------------------------
    retval = nc_open(_filename.c_str(),NC_NOWRITE,&ncid);
    HandleNetCDFErrors(retval);

    retval = nc_inq_varid(ncid,_varname.c_str(),&varid_f);
    HandleNetCDFErrors(retval);
      
    // -------------------------------
    // allocate aTmp matrix
    // -------------------------------
    dim1 = 1; dim2 = 1; dim3 = 1;

    switch(_dim_order)
    {
      case(1): 
        dim1 = _GridDims[0]; dim2 = _GridDims[1]; dim3 = iChunkSize;   break;// dimensions are (x,y,t)
      case(2):
        dim1 = _GridDims[1]; dim2 = _GridDims[0]; dim3 = iChunkSize;   break;// dimensions are (y,x,t)
      case(3):
        dim1 = _GridDims[0]; dim2 = iChunkSize;   dim3 = _GridDims[1]; break;// dimensions are (x,t,y)
      case(4):
        dim1 = iChunkSize;   dim2 = _GridDims[0]; dim3 = _GridDims[1]; break; // dimensions are (t,x,y)
      case(5):
        dim1 = _GridDims[1]; dim2 = iChunkSize;   dim3 = _GridDims[0]; break; // dimensions are (y,t,x)
      case(6):
        dim1 = iChunkSize;   dim2 = _GridDims[1]; dim3 = _GridDims[0]; break;// dimensions are (t,y,x)    
    }
 
    
    // -------------------------------
    //emulate VLA 3D array storage - store 3D array as vector using Row Major Order
    // -------------------------------
    double *aVec=new double[dim1*dim2*dim3];//stores actual data
    for(int i=0; i<dim1*dim2*dim3; i++) {
      aVec[i]=NETCDF_BLANK_VALUE;
    }
    ExitGracefullyIf(aVec==NULL,"CForcingGrid::ReadData : aVec",OUT_OF_MEMORY);

    double ***aTmp=NULL; //stores pointers to rows/columns of 3D data 
    aTmp=new double **[dim1];
    ExitGracefullyIf(aTmp==NULL,"CForcingGrid::ReadData : aTmp(0)",OUT_OF_MEMORY);
    for(int it=0;it<dim1;it++){
      aTmp[it]=NULL;
      aTmp[it]=new double *[dim2];
      ExitGracefullyIf(aTmp[it]==NULL,"CForcingGrid::ReadData : aTmp(1)",OUT_OF_MEMORY);
      for(int ir=0;ir<dim2;ir++){
        aTmp[it][ir]=&aVec[it*dim2*dim3+ir*dim3]; //points to correct location in aVec data storage
      }
    }
        
    // -------------------------------
    // Read chunk of data.
    // -------------------------------
    int       start_point = _ChunkSize * _iChunk;
    size_t    nc_start [] ={(size_t)(start_point),0,0}; 
    size_t    nc_length[] ={(size_t)(dim1),(size_t)(dim2),(size_t)(dim3)};
    ptrdiff_t nc_stride[] ={1,1,1};

    if(Options.noisy){
      printf("  Dim of chunk read: nCol = %i   nRow = %i   nTime = %i\n",dim3,dim2,dim1);
      printf("  start  chunk: (%lu, %lu, %lu)\n", nc_start[0], nc_start[1], nc_start[2]);
      printf("  length chunk: (%lu, %lu, %lu)\n",nc_length[0],nc_length[1],nc_length[2]);
      printf("  stride chunk: (%lu, %lu, %lu)\n",nc_stride[0],nc_stride[1],nc_stride[2]);
    }

    retval=nc_get_vars_double(ncid,varid_f,nc_start,nc_length,nc_stride,&aTmp[0][0][0]);
    HandleNetCDFErrors(retval);
    new_chunk_read = true;
    
    /*for(int it=0;it<dim1;it++){
      for(int ir=0;ir<dim2;ir++){
       for (int it=0;it<dim3;it++){
        cout<<aTmp[it][ir][it]<<" ";
        }
      }
      cout<<endl;
    }
     cout<<"NCGETVARS"<<endl;*/
    // -------------------------------
    // Copy all data from aTmp array to member array _aVal.
    // -------------------------------
    int irow,icol;

    if(_dim_order == 1) {
      for(int it=0; it<iChunkSize; it++){                     // loop over time points in buffer
        for(int ic=0; ic<_nNonZeroWeightedGridCells; ic++){   // loop over non-zero weighted grid cells
          CellIdxToRowCol(_IdxNonZeroGridCells[ic],irow,icol);
          _aVal[it][ic]=aTmp[icol][irow][it];
        }
      }
    }
    else if(_dim_order == 2) {
      for(int it=0; it<iChunkSize; it++){                     // loop over time points in buffer
        for(int ic=0; ic<_nNonZeroWeightedGridCells; ic++)    // loop over non-zero weighted grid cells
        {
          CellIdxToRowCol(_IdxNonZeroGridCells[ic],irow,icol);
          _aVal[it][ic]=aTmp[irow][icol][it];
        }
      }
    }
    else if(_dim_order == 3) {
      for(int it=0; it<iChunkSize; it++){                     // loop over time points in buffer
        for(int ic=0; ic<_nNonZeroWeightedGridCells; ic++){   // loop over non-zero weighted grid cells
          CellIdxToRowCol(_IdxNonZeroGridCells[ic],irow,icol);
          _aVal[it][ic]=aTmp[icol][it][irow];
        }
      }
    }
    else if(_dim_order == 4) {
      for(int it=0; it<iChunkSize; it++){                     // loop over time points in buffer
        for(int ic=0; ic<_nNonZeroWeightedGridCells; ic++){   // loop over non-zero weighted grid cells
          CellIdxToRowCol(_IdxNonZeroGridCells[ic],irow,icol);
          _aVal[it][ic]=aTmp[it][icol][irow];
        }
      }
    }
    else if(_dim_order == 5) {
      for(int it=0; it<iChunkSize; it++){                     // loop over time points in buffer
        for(int ic=0; ic<_nNonZeroWeightedGridCells; ic++){   // loop over non-zero weighted grid cells
          CellIdxToRowCol(_IdxNonZeroGridCells[ic],irow,icol);
          _aVal[it][ic]=aTmp[irow][it][icol];
        }
      }
    }
    else if(_dim_order == 6) {
      for(int it=0; it<iChunkSize; it++){                      // loop over time points in buffer
        for(int ic=0; ic<_nNonZeroWeightedGridCells; ic++){    // loop over non-zero weighted grid cells
          CellIdxToRowCol(_IdxNonZeroGridCells[ic],irow,icol);
          _aVal[it][ic]=aTmp[it][irow][icol];
        }
      }
    }

    // -------------------------------
    // Print data on screen
    // -------------------------------
    // for (int it=0; it<min(_ChunkSize,10); it++)  // York.nc : loop over time points in buffer: 4 -->
    //   {
    //     printf("t = %i (Read)\n",it);
    //     printf("   %f %f %f %f \n",_aVal[it][ 0],_aVal[it][ 1],_aVal[it][ 2],_aVal[it][ 3]);
    //     printf("   %f %f %f %f \n",_aVal[it][ 4],_aVal[it][ 5],_aVal[it][ 6],_aVal[it][ 7]);
    //     printf("   %f %f %f %f \n",_aVal[it][ 8],_aVal[it][ 9],_aVal[it][10],_aVal[it][11]);
    //     printf("   %f %f %f %f \n",_aVal[it][12],_aVal[it][13],_aVal[it][14],_aVal[it][15]);
    //     printf("\n");
    //   }

    // -------------------------------
    //delete aTmp, aVec array
    // -------------------------------
    for(int it=0;it<dim1;it++){delete [] aTmp[it];}delete [] aTmp;   
    delete [] aVec;   
    

    // -------------------------------
    // Close NetCDF file
    // -------------------------------
    retval = nc_close(ncid);
    HandleNetCDFErrors(retval);

  }// end if(_iChunk != iChunk_new) 

#endif   // end #ifdef _RVNETCDF_

  return new_chunk_read;

}

///////////////////////////////////////////////////////////////////
/// \brief   Enables queries of time series values using model time
/// \details Calculates _t_corr, correction to global model time, checks for overlap
//           with model duration, resamples to model time step/day
/// \remark  t=0 corresponds to first day with recorded values at that gauge
///
/// \param   model_start_day  [in] Julian start day of model
/// \param   model_start_year [in] start year of model
/// \param   model_duration   [in] Duration of model, in days
/// \param   timestep         [in] model timestep, in days
///
void CForcingGrid::Initialize( const double model_start_day,   // fractional day of the year (here called Julian day) [days]
                               const    int model_start_year,  //         [year]
                               const double model_duration,    //         [days]
                               const double model_timestep,    // delta t [days]
                               const optStruct &Options        // Options
  )
{
  //_t_corr is number of days between model start date and forcing chunk
  // start date (positive if data exists before model start date)

  _t_corr = -TimeDifference(model_start_day,model_start_year,_start_day,_start_year);

  //QA/QC: Check for overlap issues between time series duration and model duration
  //------------------------------------------------------------------------------
  double duration = (double)(_nPulses)*_interval;
  double local_simulation_start = (_t_corr);
  double local_simulation_end = (model_duration + _t_corr);

  if (Options.noisy){ cout << endl; }
  if (duration < local_simulation_start)  //out of data before simulation starts!
  {
    cout << "Initialize forcing grid  '" << _varname.c_str() << "'" << endl;
    cout << "  time series start day, year, duration :" << _start_day << "," << _start_year << " " << duration << endl;
    cout << "  model start day, year, duration :" << model_start_day << "," << model_start_year << " " << model_duration << endl;
    ExitGracefully(
      "CForcingGrid::Initialize: gridded forcing data not available for entire model simulation duration", BAD_DATA);
  }
  if (duration + model_timestep < local_simulation_end)    //run out of data before simulation finishes
  {                                                      //+model_timesteps is for coincdent duration & data
    cout << "Initialize forcing grid  '" << _varname.c_str() << "'" << endl;
    cout << "  time series start day, year, duration :" << _start_day << "," << _start_year << " " << duration << endl;
    cout << "  model start day, year, duration :" << model_start_day << "," << model_start_year << " " << model_duration << endl;
    ExitGracefully(
      "CForcingGrid::Initialize: gridded forcing data not available at end of model simulation", BAD_DATA);
  }
  if ((local_simulation_start<0) || (_start_year>model_start_year))     //data does not begin until after simulation
  {
    cout << "Initialize forcing grid  '" << _varname.c_str() << "'" << endl;
    cout << "  time series start day, year, duration :" << _start_day << "," << _start_year << " " << duration << endl;
    cout << "  model start day, year, duration :" << model_start_day << "," << model_start_year << " " << model_duration << endl;
    ExitGracefully(
      "CForcingGrid::Initialize: gridded forcing data not available at beginning of model simulation", BAD_DATA);
  }

  /*if (Options.noisy){ 
    cout << "Initialize forcing grid '" << _varname.c_str() << "'" << endl; 
    cout << "  time series start day, year, duration : " << _start_day << "," << _start_year << " " << duration << endl; 
    cout << "  model start day, year, duration       : " << model_start_day << "," << model_start_year << " " << model_duration << endl; 
    cout << "Finished forcing grid'" << _varname.c_str() << "'" << endl<< endl; 
  }*/

}

///////////////////////////////////////////////////////////////////
/// \brief returns row and column index of cell ID
///
/// \param cellid  [in] int of cell ID
//
void CForcingGrid::CellIdxToRowCol(const int cellid, int &row, int &column)
{
  int ncols = GetCols();
  row    = int(cellid / ncols);
  column = cellid % ncols;
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
/// \param DimNames  [in] array of three strings containing the names of the dimensions NetCDF file\n
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
/// \param GridDims  [in] array of three intgers determining length of each dimension (x=#cols, y=#rows, t=#timepoints)\n
///                          order must be: x-dimension, y-dimension, time-dimension
//
void CForcingGrid::SetGridDims(const int GridDims[3])
{
  _GridDims[0]=GridDims[0];
  _GridDims[1]=GridDims[1];
  _GridDims[2]=GridDims[2];
}

///////////////////////////////////////////////////////////////////
/// \brief sets the _nNonZeroWeightedGridCells in class CForcingGrid
///
/// \param nNonZeroWeightedGridCells  [in] number of grid cells with non-zero weighting across all HRUs
//
void CForcingGrid::SetNumberNonZeroGridCells(const int nNonZeroWeightedGridCells)
{
  _nNonZeroWeightedGridCells = nNonZeroWeightedGridCells;
}

void CForcingGrid::SetIdxNonZeroGridCells(const int nHydroUnits, const int nGridCells)
{
  bool *nonzero;

  nonzero = NULL;
  nonzero =  new bool [nGridCells]; // \todo [optimization]: may wish to declare as static
  for (int il=0; il<nGridCells; il++) { // loop over all cells of NetCDF
    nonzero[il] = false;
  }

  //printf("ForcingGrid: nNonZeroGridCells = %i\n",nNonZeroGridCells);

  _IdxNonZeroGridCells = NULL;
  _IdxNonZeroGridCells = new int [_nNonZeroWeightedGridCells];
  for (int il=0; il<_nNonZeroWeightedGridCells; il++) { // loop over all cells non-zero weighted grid cells
    _IdxNonZeroGridCells[il] = -1;
  }

  if (_GridWeight != NULL){
    for (int ik=0; ik<nHydroUnits; ik++) {  // loop over HRUs
      for (int il=0; il<nGridCells; il++) { // loop over all cells of NetCDF
        if ( _GridWeight[ik][il] > 0.00001 ) {
          nonzero[il] = true;
        }
      }
    }
  }
  else{
    ExitGracefully(
      "CForcingGrid: SetIdxNonZeroGridCells: _GridWeight is not allocated yet. Call AllocateWeightArray(nHRUs) first.", BAD_DATA);
  }

  int ic = 0;
  for (int il=0; il<nGridCells; il++) { // loop over all cells of NetCDF
    if ( nonzero[il] ) {
      _IdxNonZeroGridCells[ic] = il;
      //printf("ForcingGrid: _IdxNonZeroGridCells[%i] = %i\n",iCell,il);
      ic++;
    }
  }
  //determine _aFirstNonZeroWt and _aLastNonZeroWt
  _aFirstNonZeroWt=new int [_nHydroUnits];
  _aLastNonZeroWt =new int [_nHydroUnits];
  bool firstfound=false;
  int cell_idx;
  for(int k=0; k<_nHydroUnits;k++){
    firstfound=false;
    for(int ic = 0; ic < _nNonZeroWeightedGridCells; ic++)
    {
      cell_idx = _IdxNonZeroGridCells[ic];
      if (_GridWeight[k][cell_idx]>REAL_SMALL) 
      {
        if(!firstfound){ _aFirstNonZeroWt[k]=ic; firstfound=true; }
        _aLastNonZeroWt[k]=ic;
      }
    }
  }

  delete[] nonzero;
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

// ///////////////////////////////////////////////////////////////////
// /// \brief sets the _aVal in class CForcingGrid
// ///
// /// \param x_col  [in] column index of _aVal[t][row][col]
// /// \param y_row  [in] row    index of _aVal[t][row][col]
// /// \param t      [in] time   index of _aVal[t][row][col]
// /// \param aVal   [in] value to be set
// //
// void CForcingGrid::SetValue( const int x_col, const int y_row, const int t, const double aVal) {
//   _aVal[t][y_row][x_col] = aVal;
// }

///////////////////////////////////////////////////////////////////
/// \brief sets the _aVal in class CForcingGrid
///
/// \param ic    [in] Index of grid cell with non-zero weighting (value between 0 and _nNonZeroWeightedGridCells)
/// \param t      [in] time   index of _aVal[t][idx]
/// \param aVal   [in] value to be set
//
void CForcingGrid::SetValue( const int ic, const int t, const double aVal) {
  _aVal[t][ic] = aVal;
}
///////////////////////////////////////////////////////////////////
/// \brief indicates that the forcing grid stores an accumulated variable
//
void CForcingGrid::SetToDeaccumulate(){
  _deaccumulate=true;
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
  _GridWeight = NULL;
  _GridWeight = new double *[nHydroUnits];
  for (int ik=0; ik<nHydroUnits; ik++) {  // loop over HRUs
    _GridWeight[ik] = new double [nGridCells];
    for (int il=0; il<nGridCells; il++) { // loop over all cells of NetCDF
      _GridWeight[ik][il] = 0.0;
    }
  }
}

///////////////////////////////////////////////////////////////////
/// \brief sets one entry of _GridWeight[HRUID, CellID] = weight in class CForcingGrid
///
/// \param HRUID  [in] HRU id (numbered: 0,1,2,...,nHydroUnits-1)
/// \param CellID [in] cell ID in NetCDF; cells are numbered linewise from left to right starting with 0:
///                         0 1 2 3 4 \n
///                         5 6 7 8 9 \n
///                         ...
/// \param weight [in] weight[k][l] is fraction of forcing for HRUID k is from grid cell l=(i,j)
///                    and grid cell index l is derived by l = (i-1) * NC + j
///                    where i and j are the row and column of cell l respectively and
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
      "CForcingGrid: SetWeightVal: _GridWeight is not allocated yet. Call AllocateWeightArray(nHRUs) first.",BAD_DATA);
  }
#endif
  _GridWeight[HRUID][CellID] = weight;

}

///////////////////////////////////////////////////////////////////
/// \brief checks if sum(_GridWeight[HRUID, :]) = 1.0 for all HRUIDs
///
/// \param nHydroUnits [in] Number of HRUs must be given
///                         (is read from :GridWeights block in *.rvi file)
/// \param nGridCells  [in] Number of Grid Cells must be given
///                         (is read from :GridWeights block in *.rvi file)
//
/// return true if sum for each HRUID is one; otherwise false
//
bool   CForcingGrid::CheckWeightArray(const int nHydroUnits, const int nGridCells)
{
  double sum_HRU;
  bool   check = true;

  if (_GridWeight != NULL){
    for (int ik=0; ik<nHydroUnits; ik++) {  // loop over HRUs
      sum_HRU = 0.0;
      for (int il=0; il<nGridCells; il++) { // loop over all cells of NetCDF
        sum_HRU = sum_HRU + _GridWeight[ik][il];
      }
      if(abs(sum_HRU - 1.0) < 0.05) {//repair if less than 5%
        for(int il=0; il<nGridCells; il++) {
          _GridWeight[ik][il]/=sum_HRU;
        }    
        sum_HRU = 1.0;
      }
      if ( abs(sum_HRU - 1.0) > 0.0001 ) {
        printf("HRU ID = %i    Sum Forcing Weights = %f\n", ik,sum_HRU);
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
/// \brief determines number of non-zero weighted grid cells
///
/// \param nHydroUnits [in] Number of HRUs must be given
///                         (is read from :GridWeights block in *.rvi file)
/// \param nGridCells  [in] Number of Grid Cells must be given
///                         (is read from :GridWeights block in *.rvi file)
//
/// return number of non-zero weighted grid cells
//
int CForcingGrid::NumberNonZeroWeightedGridCells(const int nHydroUnits, const int nGridCells)
{
  bool *nonzero;
  int nNonZeroGridCells;

  nNonZeroGridCells = 0;
  nonzero =  new bool [nGridCells];  // false if grid cell has only zero weight;
  //                                  // true  if grid cell has at least one non-zero weight;

  for (int il=0; il<nGridCells; il++) { // loop over all cells of NetCDF
    nonzero[il] = false;
  }

  if (_GridWeight != NULL){
    for (int ik=0; ik<nHydroUnits; ik++) {  // loop over HRUs
      for (int il=0; il<nGridCells; il++) { // loop over all cells of NetCDF
        if ( _GridWeight[ik][il] > 0.00001 ) {
          nonzero[il] = true;
        }
      }
    }
  }
  else{
    ExitGracefully(
      "CForcingGrid: CheckWeightArray: _GridWeight is not allocated yet. Call AllocateWeightArray(nHRUs) first.", BAD_DATA);
  }

  for (int il=0; il<nGridCells; il++) { // loop over all cells of NetCDF
    if ( nonzero[il] ) {
      nNonZeroGridCells += 1;
    }
  }

  return nNonZeroGridCells;
}

///////////////////////////////////////////////////////////////////
/// \brief weighting of HRU and CellID pair
///
/// return weighting of HRU and CellID pair
//
double CForcingGrid::GetGridWeight(const int HRUID,
                                   const int CellID) const
{
  return _GridWeight[HRUID][CellID];
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
int CForcingGrid::GetRows() const{return _GridDims[1];}

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
/// \brief Returns ID of i-th grid cell with non-zero weighting
/// \return ID of i-th grid cell with non-zero weighting
//
int CForcingGrid::GetIdxNonZeroGridCell(int ic) const{return _IdxNonZeroGridCells[ic];}

///////////////////////////////////////////////////////////////////
/// \brief Returns current chunk size
/// \return Current chunk size
//
int  CForcingGrid::GetChunkSize() const{return _ChunkSize;}

///////////////////////////////////////////////////////////////////
/// \brief Returns name of forcing grid data
/// \return name of forcing data
//
forcing_type CForcingGrid::GetName()  const{return _ForcingType;}

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
  return int((_t_corr + t) * round(1.0/_interval)+0.5*tstep)  % _ChunkSize;
}

///////////////////////////////////////////////////////////////////
/// \brief returns weighted value of gridded forcing in HRU k over timestep
/// \param k    [in] HRU index
/// \param t      [in] model time [days]
/// \return tstep [in] model time step [days]
/// \returns timestep weighted value of gridded forcing in HRU k
//
double CForcingGrid::GetWeightedValue(const int k,const double &t,const double &tstep) const
{

  int idx_new = GetTimeIndex(t,tstep);
  int nSteps = (int)(max(1.0,round(tstep/_interval)));//# of intervals in time step
  double wt;
  double sum=0;
  for(int ic = _aFirstNonZeroWt[k]; ic <= _aLastNonZeroWt[k]; ic++)
  {
    wt       = _GridWeight[k][_IdxNonZeroGridCells[ic]];
    sum += wt * GetValue(ic,idx_new,nSteps);
  }
  return sum;
}
///////////////////////////////////////////////////////////////////
/// \brief returns daily weighted value of gridded forcing in HRU k
/// \param k    [in] HRU index
/// \param t      [in] model time [days]
/// \param tstep [in] model time step [days]
/// \returns daily weighted value of gridded forcing in HRU k
//
double CForcingGrid::GetDailyWeightedValue(const int k,const double &t,const double &tstep) const
{
  int idx_new_day = GetTimeIndex((double)(int(t)),tstep);//start of day
  int nStepsDaily  = (int)(round(1.0/_interval));//# of intervals in day
  double wt;
  double sum=0;
  for(int ic = _aFirstNonZeroWt[k]; ic <= _aLastNonZeroWt[k]; ic++)
  {
    wt       = _GridWeight[k][_IdxNonZeroGridCells[ic]];
    sum += wt * GetValue(ic,idx_new_day,nStepsDaily);
  }
  return sum;
}
///////////////////////////////////////////////////////////////////
/// \brief returns daily weighted snowfrac value if this is a gridded snow dataset and pRain is provided 
/// \param k    [in] HRU index
/// \param t      [in] model time [days]
/// \param pRain [in] pointer to gridded rain dataset
/// \param tstep [in] model time step [days]
//
double CForcingGrid::GetWeightedAverageSnowFrac(const int k,const double &t,const double &tstep,const CForcingGrid *pRain) const
{

  int nSteps = (int)(max(1.0,round(tstep/_interval)));//# of intervals in time step
  double wt;
  double sum=0.0;
  double snow; double rain;
  for(int ic = _aFirstNonZeroWt[k]; ic <= _aLastNonZeroWt[k]; ic++)
  {
    wt       = _GridWeight[k][_IdxNonZeroGridCells[ic]];
    snow = GetValue(ic, t, nSteps);
    if(snow>0.0){
      rain=pRain->GetValue(ic, t, nSteps);
      sum+= wt * snow/(snow+rain);
    }
  }
  return sum;
}

///////////////////////////////////////////////////////////////////
/// \brief Returns magnitude of time series data point for which t is a float index
/// \param ic    [in] Index of grid cell with non-zero weighting (value between 0 and _nNonZeroWeightedGridCells)
/// \param t      [in] Time index
/// \return Magnitude of time series data point for which t is an index
//
double CForcingGrid::GetValue(const int ic, const double &t) const 
{
  return _aVal[(int)t][ic];
}

///////////////////////////////////////////////////////////////////
/// \brief Returns average over n timesteps of time series data point for which t is an index
/// \param ic     [in] Index of grid cell with non-zero weighting (value between 0 and _nNonZeroWeightedGridCells)
/// \param t      [in] local time (with respect to chunk start, in days)
/// \param n      [in] Number of time steps
/// \return Magnitude of time series data point for which t is an index
//
double CForcingGrid::GetValue(const int ic, const double &t, const int nsteps) const 
{


  int ii_start=(int)(t);
  int lim=min(nsteps,_ChunkSize-ii_start); 
  double sum = 0.0;
  for (int ii=ii_start; ii<ii_start+lim;ii++){
    sum += _aVal[ii][ic];
  }
  sum /= (double)(lim);

  return sum;
}

///////////////////////////////////////////////////////////////////
/// \brief Returns average over n timesteps of time series data point for which t is an index
/// \param idx    [in] Index of grid cell with non-zero weighting (value between 0 and _nNonZeroWeightedGridCells)
/// \param t      [in] Time index
/// \param n      [in] Number of time steps
/// \return Magnitude of time series data point for which t is an index
//
/*double CForcingGrid::GetValue_ave(const int idx, const double &t, const int n) const 
{
  double sum = 0.0;
  for (int ii=0; ii<n; ii++) {
    sum += _aVal[min(_ChunkSize-1,(int)t+ii)][idx];
  };
  sum /= double(n);

  return sum;
}*/ //identical to GetValue()

///////////////////////////////////////////////////////////////////
/// \brief Returns minimum of n timesteps of time series data point for which t is an index
/// \param ic     [in] Index of grid cell with non-zero weighting (value between 0 and _nNonZeroWeightedGridCells)
/// \param t      [in] Time index
/// \param n      [in] Number of time steps
/// \return Magnitude of time series data point for which t is an index
//
double CForcingGrid::GetValue_min(const int ic, const double &t, const int nsteps) const 
{
  /*double mini = ALMOST_INF ;
  for (int ii=0; ii<nsteps; ii++) {
    if (_aVal[min(_ChunkSize-1,(int)t+ii)][idx] < mini) {mini = _aVal[min(_ChunkSize-1,(int)t+ii)][idx]; };
  };
  return mini;*/

  double mini = ALMOST_INF ;
  int ii_start=(int)(t);
  int lim=min(nsteps,_ChunkSize-ii_start); 
  for (int ii=ii_start; ii<ii_start+lim;ii++){
    if(_aVal[ii][ic] < mini){mini=_aVal[ii][ic];}
  }
  return mini;

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
  /*double maxi = -ALMOST_INF ;
  for (int ii=0; ii<nsteps; ii++) {
    if (_aVal[min(_ChunkSize-1,(int)t+ii)][idx] > maxi) {maxi = _aVal[min(_ChunkSize-1,(int)t+ii)][idx]; };
  };

  return maxi;*/

  double maxi = -ALMOST_INF ;
  int ii_start=(int)(t);
  int lim=min(nsteps,_ChunkSize-ii_start); 
  for (int ii=ii_start; ii<ii_start+lim;ii++){
    if(_aVal[ii][ic] > maxi){maxi=_aVal[ii][ic];}
  }
  return maxi;
}

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
  double model_timestep;

  // delta t of model
  model_timestep  = Options.timestep;

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
  double model_timestep;

  // delta t of model
  model_timestep  = Options.timestep;

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
bool CForcingGrid::ShouldDeaccumulate()                const{return _deaccumulate;}
