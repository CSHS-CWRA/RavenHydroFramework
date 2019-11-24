/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2019 the Raven Development Team
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
/// \param is_3D         [in] true if NetCDF grid is 3D (gridded) data rather than 2D (station) data
//
CForcingGrid::CForcingGrid(string       ForcingType,
                           string       filename,
                           string       varname,
                           string       DimNames[3],
                           bool         is_3D
                          )
{
  _ForcingType  = GetForcingTypeFromString(ForcingType);
  _filename     = filename;
  _varname      = varname;
  _DimNames[0]  = DimNames[0]; _DimNames[1]  = DimNames[1]; _DimNames[2]  = DimNames[2];
  _is_3D        = is_3D;

  _interval     = 1.0;
  _steps_per_day= 1;
  
  _TimeShift    = 0.0;
  _LinTrans_a   = 1.0;
  _LinTrans_b   = 0.0;
  _deaccumulate = false;
  _period_ending= false;

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

  // initialized in AllocateWeightArray,SetIdxNonZeroGridCells()
  _GridWeight          = NULL;
  _GridWtCellIDs       = NULL;
  _CellIDToIdx         = NULL;
  _nWeights            = NULL;
  _IdxNonZeroGridCells = NULL;

  //initialized in SetAttributeVarName
  _aLatitude           = NULL;
  _aLongitude          = NULL;
  _aElevation          = NULL;
  _aStationIDs         = NULL;
  _AttVarNames[0]="NONE";
  _AttVarNames[1]="NONE";
  _AttVarNames[2]="NONE";
  _AttVarNames[3]="NONE";
}

void CForcingGrid::ReadAttGridFromNetCDF(const int ncid, const string varname,const int nrows, const int ncols,double *values) 
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
    else       { nCells = nrows; }
    
    delete [] values; //re-read if buffered
    values=new double [_nNonZeroWeightedGridCells]; //allocate memory

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
      double  **aTmp2D=NULL; //stores pointers to rows/columns of 2D data
      aTmp2D=new double *[nrows];
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
      retval=nc_get_vars_double(ncid,varid,nc_start,nc_length,nc_stride,&aVec[0]);   HandleNetCDFErrors(retval);
      //copy matrix
      for(int ic=0; ic<_nNonZeroWeightedGridCells; ic++) {   // loop over non-zero weighted grid cells
        CellIdxToRowCol(_IdxNonZeroGridCells[ic],irow,icol);
        values[ic]=aVec[irow];
      }
    }
    delete[] aVec;
  }
#endif

}
// same as above but for string information
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
  
  _aVal=NULL;
  _aVal =  new double *[_ChunkSize];
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
  int    BytesPerTimestep;      // Memory requirement for one timestep of gridded forcing file [Bytes]

  int    iatt;
  char   attrib_name[500];

  int    retval1;
  size_t att_len;        // length of the attribute's text
  nc_type att_type;      // type of attribute

  string dash;           // to check format of time unit string
  string colon;          // to check format of time unit string
  string unit_t_str;     // to check format of time unit string
  int    ntime;          // number of time steps
  nc_type type;          // type of a NetCDF variable, e.g., nc_INT, nc_FLOAT, nc_DOUBLE

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
  if(_aElevation==NULL) {
    string warning="Since no elevation data are in NetCDF file "+_filename+", all orographic corrections for "+ForcingToString(_ForcingType)+" will be disabled unless :StationElevations command is present.";
    WriteAdvisory(warning,Options.noisy);
  }

  // open NetCDF read-only; ncid will be set
  if(Options.noisy) { cout<<"Initializing grid file "<<_filename<<endl; }
  retval = nc_open(_filename.c_str(),NC_NOWRITE,&ncid);          HandleNetCDFErrors(retval);


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
    retval = nc_inq_dimlen(ncid,dimid_t,&GridDim_t);             HandleNetCDFErrors(retval);
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
  retval = nc_inq_varid(ncid,_varname.c_str(),&varid_f);         HandleNetCDFErrors(retval);

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
  retval = nc_inq_attlen(ncid,varid_t,"units",&att_len);                   HandleNetCDFErrors(retval);

  unit_t=new char[att_len+1];
  retval = nc_get_att_text(ncid,varid_t,"units",unit_t);                   HandleNetCDFErrors(retval);// read attribute text
  unit_t[att_len] = '\0';// add string determining character
  
  // check that unit of time is in format "[days/minutes/...] since YYYY-MM-DD HH:MM:SS{+0000}"
  unit_t_str=to_string(unit_t);
  if(!IsValidNetCDFTimeString(unit_t_str)) 
  {
    cout<<"time unit string: "<<unit_t_str<<endl;
    ExitGracefully("CForcingGrid::ForcingGridInit: time unit string is not in the format '[days/hours/...] since YYYY-MM-DD HH:MM:SS +0000' !",BAD_DATA);
  }
  
  // calendar attribute
  // ----------------------------------------------------------------------------------------------
  retval = nc_inq_att(ncid, varid_t, "calendar", &att_type, &att_len);      
  if (retval == NC_ENOTATT) {
    // (a) if not found, set to proleptic_gregorian
    calendar = StringToCalendar("PROLEPTIC_GREGORIAN");
  }
  else 
  {  
    HandleNetCDFErrors(retval);
    // (b) if found, read and make sure '\0' is terminating character
    retval = nc_inq_attlen(ncid, varid_t, "calendar", &att_len);            HandleNetCDFErrors(retval);// inquire length of attribute's text
    
    char *calendar_t = new char [att_len + 1]; // allocate memory of char * to hold attribute's text
    retval = nc_get_att_text(ncid, varid_t, "calendar", calendar_t);        HandleNetCDFErrors(retval);// read attribute text
    calendar_t[att_len] = '\0';                // add string determining character

    calendar = StringToCalendar(calendar_t);
    delete [] calendar_t;
  }
			 
  // time attribute: set my_time[]
  // ----------------------------------------------------------------------------------------------
  if (_is_3D) {ntime = _GridDims[2];}
  else        {ntime = _GridDims[1];}

  // the final time variable will be in double precision no matter what the incoming time variable type was
  double *my_time=NULL;
  my_time=new double [ntime];
  ExitGracefullyIf(my_time==NULL,"CForcingGrid::ForcingGridInit",OUT_OF_MEMORY);

  // find out if time is an integer or double variable
  retval = nc_inq_var( ncid, varid_t,
                       NULL,     // out: name
                       &type,    // out: type: nc_INT, nc_FLOAT, nc_DOUBLE, etc.
                       NULL,     // out: dims
                       NULL,     // out: dimids
                       NULL );   HandleNetCDFErrors(retval);// out: natts

  if (type == NC_DOUBLE)    // time is given as double --> use as is
  {
    retval = nc_get_var_double(ncid,varid_t,&my_time[0]);  HandleNetCDFErrors(retval);
  }
  else if (type == NC_FLOAT)    // time is given as float (single precision) --> convert to double
  {
    float *ttime=NULL;
    ttime=new float [ntime];
    ExitGracefullyIf(ttime==NULL,"CForcingGrid::ForcingGridInit",OUT_OF_MEMORY);

    retval = nc_get_var_float(ncid,varid_t,&ttime[0]);     HandleNetCDFErrors(retval);

    // convert to double to proceed
    for (int itime=0; itime<ntime;itime++){       // loop over all time steps
      my_time[itime]=double(ttime[itime]);        // convert to double
    }
    delete[] ttime;
  }
  else if ((type == NC_INT) || (type == NC_INT64))   // time is given as integer --> convert to double
  {
    int *ttime=NULL;
    ttime=new int [ntime];
    ExitGracefullyIf(ttime==NULL,"CForcingGrid::ForcingGridInit",OUT_OF_MEMORY);

    retval = nc_get_var_int(ncid,varid_t,&ttime[0]);       HandleNetCDFErrors(retval);

    // convert to double to proceed
    for (int itime=0; itime<ntime;itime++){       // loop over all time steps
      my_time[itime]=(double) ttime[itime];        // convert to double
    }
    delete [] ttime;
  }
  else
  {
    ExitGracefully("CForcingGrid: ForcingGridInit: time variable is not of type DOUBLE, FLOAT, or INT. Contact the developers to implement the type you need.",BAD_DATA);
  }

  // -------------------------------
  // check if we have equal time steps
  // -------------------------------  
  double delta_t = my_time[1] - my_time[0];
  for (int ii=1; ii<ntime-1 ; ii++)
  {
    double delta_t2 = my_time[ii+1] - my_time[ii];
    if ( abs(delta_t2 - delta_t) > 0.00001 )
    {
      printf("\n\n\nCForcingGrid: ForcingGridInit: variable name: %s\n",_varname.c_str());
      printf("\n\n\n              expected delta_t = %15.8lf\n",delta_t);
      printf("\n\n\n              derived  delta_t = %15.8lf\n",delta_t2);
      printf("\n\n\n                       i  =%i  t[i]  =%15.8lf\n",ii,  my_time[ii]);
      printf("\n\n\n                       i+1=%i  t[i+1]=%15.8lf\n",ii+1,my_time[ii+1]);
      ExitGracefully("CForcingGrid: ForcingGridInit: time steps are not equal in gridded input",BAD_DATA);
    }
  }
  
  // -------------------------------
  // delta t and interval (in days) / steps per day
  // -------------------------------
  double time_zone=0;
  time_struct tt;
  if (strstr(unit_t, "hours")) 
  {  
    tt =TimeStructFromNetCDFString(unit_t_str,"hours",calendar,time_zone);
    _interval   = (my_time[1] - my_time[0])/HR_PER_DAY;
    ExitGracefullyIf(_interval<=0,"CForcingGrid: ForcingGridInit: Interval is negative!",BAD_DATA);
    AddTime(tt.julian_day,tt.year,my_time[0]/HR_PER_DAY,calendar,_start_day,_start_year);
  }
  else if (strstr(unit_t, "days")) 
  {  
    tt =TimeStructFromNetCDFString(unit_t_str,"days",calendar,time_zone);
    _interval   = (my_time[1] - my_time[0]);
    ExitGracefullyIf(_interval<=0, "CForcingGrid: ForcingGridInit: Interval is negative!",BAD_DATA);
    AddTime(tt.julian_day,tt.year,my_time[0],calendar,_start_day,_start_year) ;
  }
  else if (strstr(unit_t, "minutes")) 
  {  
    tt =TimeStructFromNetCDFString(unit_t_str,"minutes",calendar,time_zone);
    _interval   = (my_time[1] - my_time[0])/MIN_PER_DAY;
    ExitGracefullyIf(_interval<=0, "CForcingGrid: ForcingGridInit: Interval is negative!",BAD_DATA);
    AddTime(tt.julian_day,tt.year,my_time[0]/MIN_PER_DAY,calendar,_start_day,_start_year) ;
  }
  else if (strstr(unit_t, "seconds"))
  {  
    tt =TimeStructFromNetCDFString(unit_t_str,"seconds",calendar,time_zone);
    _interval   = (my_time[1] - my_time[0])/SEC_PER_DAY;
    ExitGracefullyIf(_interval<=0,"CForcingGrid: ForcingGridInit: Interval is negative!",BAD_DATA);
    AddTime(tt.julian_day,tt.year,my_time[0]/SEC_PER_DAY,calendar,_start_day,_start_year) ;
  }
  else{
    ExitGracefully("CForcingGrid: ForcingGridInit: this unit in time is not implemented yet (only days, hours, minutes, seconds)",BAD_DATA);
  }
  _steps_per_day=(int)(round(1.0/_interval)); //pre-calculate for speed.
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

  ExitGracefullyIf(_interval<=0,
                   "CForcingGrid: ForcingGridInit: negative time interval is not allowed",BAD_DATA);

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
  // Determine chunk size (_ChunkSize),
  //           number of chunks read in total (_nChunk) and
  // Set       id of current chunk (_iChunk)
  // depending on size of grid (#cols x #rows x buffersize x 8byte <=  10 MB=10*1024*1024 byte)
  // -------------------------------
  if ( _is_3D ) {BytesPerTimestep = 8 * _GridDims[0] * _GridDims[1];}
  else          {BytesPerTimestep = 8 * _GridDims[0];               }

  const int CHUNK_MEMORY=10; //MB

  _ChunkSize = int(max(min( (CHUNK_MEMORY *1024 * 1024) / BytesPerTimestep, ntime),1));  // number of timesteps per chunk - 10MB chunks
// _ChunkSize = max(int(round(1./_interval)),int(int(_ChunkSize*_interval)/_interval));  // make sure complete days and at least one day is read
  if(!Options.deltaresFEWS) { 
    _ChunkSize = int(int(_ChunkSize*_interval)/_interval); 
  }  // make sure chunks are complete days (have to relax for FEWS)

  _ChunkSize = max(int(round(1./_interval)),_ChunkSize);                                // make sure  at least one day is read
  //support larger chunk if model duration is small

  //_nChunk    = int(ceil(1.0*ntime / _ChunkSize));                                        // total number of chunks
  _nChunk    = int(ceil((Options.duration/_interval)/_ChunkSize)); //JRC_TIME_FIX
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

    if (strcmp(attrib_name,"long_name") == 0)// long_name of forcing
    {
      retval1 = nc_inq_attlen (ncid, varid_f, attrib_name, &att_len);// inquire length of attribute's text
      HandleNetCDFErrors(retval1);
      long_name_f = (char *) malloc(att_len + 1);// allocate memory of char * to hold attribute's text
      retval1 = nc_get_att_text(ncid, varid_f, attrib_name, long_name_f);// read attribute text
      HandleNetCDFErrors(retval1);
      long_name_f[att_len] = '\0';// add string determining character
    }
    else if (strcmp(attrib_name,"units") == 0)// unit of forcing
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
//          mainly used when sub-daily grids have to be added to model
//
void CForcingGrid::ReallocateArraysInForcingGrid( )
{
  int ntime ;  // number of time steps
  int ncells ; // number of cells (either rows*cols or nstations)

  ntime  = _GridDims[2];

  if (_is_3D) {ncells = _GridDims[0] * _GridDims[1];}
  else        {ncells = _GridDims[0]; }

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
  int     ncid;          // file unit
  int     dim1;          // length of 1st dimension in NetCDF data
  int     dim2;          // length of 2nd dimension in NetCDF data
  int     dim3;          // length of 3rd dimension in NetCDF data

  int     varid_f;       // id of forcing variable read
  double  missval;       // value of "missing_value" attribute of forcing variable 
  double  fillval;       // value of "_FillValue"    attribute of forcing variable
  size_t  att_len;       // length of the attribute's text
  nc_type att_type;      // type of attribute
  int     retval;        // error value for NetCDF routines
  int     iChunkSize;    // size of current chunk; always equal _ChunkSize except for last chunk in file (might be shorter)
  int     iChunk_new;    // chunk in which current model time step falls

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
      iChunkSize = int((Options.duration - global_model_time) / _interval);
    }

    if (Options.noisy){ 
      printf("Initialize before reading first chunk  (var = %s)\n",_varname.c_str()); 
      if (_is_3D) {printf("  Dim of chunk read: dim1 = %i   dim2 = %i   ChunkSize = %i\n",_GridDims[0],_GridDims[1],iChunkSize);}
      else        {printf("  Dim of chunk read: nStations = %i   nTime = %i\n",_GridDims[0],iChunkSize);}
    }

    // (1) set the correction time _t_corr, i.e. days between forcing data start and model simulation start
    // (2) check if forcing data cover whole simulation period
    model_start_day = Options.julian_start_day;
    model_start_yr  = Options.julian_start_year;
    model_duration  = Options.duration;  // model duration in days --> number of timesteps must be duration/timestep
    model_timestep  = Options.timestep;

    Initialize(Options);

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
    // for (int it=0; it<4; it++)  // .nc : loop over time points in buffer: 4 -->
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

  //iChunk_new = max(int(floor((_t_corr + global_model_time+TIME_CORRECTION) / (_interval * _ChunkSize))),0); //chunk index correspond to current model time
  //JRC_TIME_FIX: 
  iChunk_new = max(int(floor((global_model_time+TIME_CORRECTION) / (_interval * _ChunkSize))),0);

  if(_iChunk != iChunk_new)
  {  // current model time step is not covered by current chunk --> read new chunk
    if(Options.noisy){ printf("\n Start reading new chunk... iChunk = %i (var = %s)\n",iChunk_new,_varname.c_str()); }
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
      //iChunkSize = int((_t_corr + Options.duration - global_model_time) / _interval);
      //JRC_TIME_FIX: 
      iChunkSize = int((Options.duration - global_model_time) / _interval);
    }
    if(Options.noisy){ printf("  iChunksize: %i\n",iChunkSize); }

    // -------------------------------
    // Open NetCDF file, Get the id of the forcing data, varid_f
    // -------------------------------
    retval = nc_open(_filename.c_str(),NC_NOWRITE,&ncid);      HandleNetCDFErrors(retval);
    retval = nc_inq_varid(ncid,_varname.c_str(),&varid_f);     HandleNetCDFErrors(retval);

    // -------------------------------
    // find "_FillValue" of forcing data
    // -------------------------------
    retval = nc_inq_att(ncid, varid_f, "_FillValue", &att_type, &att_len);      
    if (retval == NC_ENOTATT) {
      // if not found, set to NETCDF_BLANK_VALUE
      fillval = NETCDF_BLANK_VALUE;
    }
    else 
    {  
      HandleNetCDFErrors(retval);
      retval = nc_get_att_double(ncid, varid_f, "_FillValue", &fillval);       HandleNetCDFErrors(retval);// read attribute value
    }

    // -------------------------------
    // find "missing_value" of forcing data
    // -------------------------------
    retval = nc_inq_att(ncid, varid_f, "missing_value", &att_type, &att_len);      
    if (retval == NC_ENOTATT) {
      // if not found, set to NETCDF_BLANK_VALUE
      missval = NETCDF_BLANK_VALUE;
    }
    else 
    {  
      HandleNetCDFErrors(retval);
      retval = nc_get_att_double(ncid, varid_f, "missing_value", &missval);     HandleNetCDFErrors(retval);// read attribute value
    }

    // -------------------------------
    // allocate aTmp matrix
    // -------------------------------
    dim1 = 1; dim2 = 1; dim3 = 1;

    if ( _is_3D ) {
      switch(_dim_order)
        {
        case(1):
          dim1 = _GridDims[0]; dim2 = _GridDims[1]; dim3 = iChunkSize;   break; // dimensions are (x,y,t)
        case(2):
          dim1 = _GridDims[1]; dim2 = _GridDims[0]; dim3 = iChunkSize;   break; // dimensions are (y,x,t)
        case(3):
          dim1 = _GridDims[0]; dim2 = iChunkSize;   dim3 = _GridDims[1]; break; // dimensions are (x,t,y)
        case(4):
          dim1 = iChunkSize;   dim2 = _GridDims[0]; dim3 = _GridDims[1]; break; // dimensions are (t,x,y)
        case(5):
          dim1 = _GridDims[1]; dim2 = iChunkSize;   dim3 = _GridDims[0]; break; // dimensions are (y,t,x)
        case(6):
          dim1 = iChunkSize;   dim2 = _GridDims[1]; dim3 = _GridDims[0]; break; // dimensions are (t,y,x)
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
    //emulate VLA 3D array storage - store 3D array as vector using Row Major Order
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
      for(int it=0;it<dim1;it++){
        aTmp3D[it]=NULL;
        aTmp3D[it]=new double *[dim2];
        ExitGracefullyIf(aTmp3D[it]==NULL,"CForcingGrid::ReadData : aTmp3D(1)",OUT_OF_MEMORY);
        for(int ir=0;ir<dim2;ir++){
          aTmp3D[it][ir]=&aVec[it*dim2*dim3+ir*dim3]; //points to correct location in aVec data storage
        }
      }
    }
    else {
      aTmp2D=new double *[dim1];
      ExitGracefullyIf(aTmp2D==NULL,"CForcingGrid::ReadData : aTmp2D(0)",OUT_OF_MEMORY);
      for(int it=0;it<dim1;it++){
        aTmp2D[it]=&aVec[it*dim2]; //points to correct location in aVec data storage
      }
    }

    // -------------------------------
    // Read chunk of data.
    // -------------------------------
    if ( _is_3D ) {

      //int       start_point = _ChunkSize * _iChunk;
      int       start_point = _ChunkSize * _iChunk+(int)(_t_corr/_interval);;//JRC_TIME_FIX:
      size_t    nc_start [3];
      size_t    nc_length[3];
      ptrdiff_t nc_stride[3];

      nc_length[0] = (size_t)(dim1); nc_stride[0] = 1;
      nc_length[1] = (size_t)(dim2); nc_stride[1] = 1;
      nc_length[2] = (size_t)(dim3); nc_stride[2] = 1;
      
      switch(_dim_order) {
      case(1): // dimensions are (x,y,t)
        nc_start[0]  = 0;
        nc_start[1]  = 0;
        nc_start[2]  = (size_t)(start_point); 
        retval=nc_get_vars_double(ncid,varid_f,nc_start,nc_length,nc_stride,&aTmp3D[0][0][0]);
        HandleNetCDFErrors(retval);
        new_chunk_read = true;
        break;
      case(2): // dimensions are (y,x,t)
        nc_start[0]  = 0;                     
        nc_start[1]  = 0;                     
        nc_start[2]  = (size_t)(start_point); 
        retval=nc_get_vars_double(ncid,varid_f,nc_start,nc_length,nc_stride,&aTmp3D[0][0][0]);
        HandleNetCDFErrors(retval);
        new_chunk_read = true;
        break;
      case(3): // dimensions are (x,t,y)
        nc_start[0]  = 0;                     
        nc_start[1]  = (size_t)(start_point);                     
        nc_start[2]  = 0; 
        retval=nc_get_vars_double(ncid,varid_f,nc_start,nc_length,nc_stride,&aTmp3D[0][0][0]);
        HandleNetCDFErrors(retval);
        new_chunk_read = true;
        break;
      case(4): // dimensions are (t,x,y)
        nc_start[0]  = (size_t)(start_point);                     
        nc_start[1]  = 0;                     
        nc_start[2]  = 0; 
        retval=nc_get_vars_double(ncid,varid_f,nc_start,nc_length,nc_stride,&aTmp3D[0][0][0]);
        HandleNetCDFErrors(retval);
        new_chunk_read = true;
        break;
      case(5): // dimensions are (y,t,x)
        nc_start[0]  = 0;                     
        nc_start[1]  = (size_t)(start_point);                     
        nc_start[2]  = 0; 
        retval=nc_get_vars_double(ncid,varid_f,nc_start,nc_length,nc_stride,&aTmp3D[0][0][0]);
        HandleNetCDFErrors(retval);
        new_chunk_read = true;
        break;
      case(6): // dimensions are (t,y,x)
        nc_start[0]  = (size_t)(start_point); 
        nc_start[1]  = 0;                     
        nc_start[2]  = 0;                     
        retval=nc_get_vars_double(ncid,varid_f,nc_start,nc_length,nc_stride,&aTmp3D[0][0][0]);
        HandleNetCDFErrors(retval);
        new_chunk_read = true;
        break;
      }

      if (Options.noisy) {
        cout<<" CForcingGrid::ReadData - is3D"<<endl;
        printf("  Dim of chunk read: dim3 = %i   dim2 = %i   dim1 = %i\n",dim3,dim2,dim1);
        printf("  start  chunk: (%zu, %zu, %zu)\n", nc_start[0], nc_start[1], nc_start[2]);
        printf("  length chunk: (%zu, %zu, %zu)\n",nc_length[0],nc_length[1],nc_length[2]);
        printf("  stride chunk: (%zu, %zu, %zu)\n",nc_stride[0],nc_stride[1],nc_stride[2]);
      }

      /*for (int it=0;it<dim1;it++){
        for (int ir=0;ir<dim2;ir++){
        for (int it=0;it<dim3;it++){
        cout<<aTmp3D[it][ir][it]<<" ";
        }
        }
        cout<<endl;
        }
        cout<<"NCGETVARS"<<endl;*/
    }
    else {

      //int       start_point = _ChunkSize * _iChunk;
      int       start_point = _ChunkSize * _iChunk+(int)(_t_corr/_interval);;//JRC_TIME_FIX:
      size_t    nc_start [2];
      size_t    nc_length[2];
      ptrdiff_t nc_stride[2];

      nc_length[0] = (size_t)(dim1); nc_stride[0] = 1;
      nc_length[1] = (size_t)(dim2); nc_stride[1] = 1;
      
      switch(_dim_order) {
      case(1): // dimensions are (station,t)
        nc_start[0]  = 0;
        nc_start[1]  = (size_t)(start_point);
        retval=nc_get_vars_double(ncid,varid_f,nc_start,nc_length,nc_stride,&aTmp2D[0][0]);
        HandleNetCDFErrors(retval);
        new_chunk_read = true;
        break;
      case(2): // dimensions are (t,station)
        nc_start[0]  = (size_t)(start_point);                     
        nc_start[1]  = 0;                     
        retval=nc_get_vars_double(ncid,varid_f,nc_start,nc_length,nc_stride,&aTmp2D[0][0]);
        HandleNetCDFErrors(retval);
        new_chunk_read = true;
        break;
      }

      if (Options.noisy) {
        cout<<" CForcingGrid::ReadData - !is3D"<<endl;
        printf("  Dim of chunk read: dim2 = %i   dim1 = %i\n",dim2,dim1);
        printf("  start  chunk: (%zu, %zu)\n", nc_start[0], nc_start[1]);
        printf("  length chunk: (%zu, %zu)\n",nc_length[0],nc_length[1]);
        printf("  stride chunk: (%zu, %zu)\n",nc_stride[0],nc_stride[1]);
      }

      // for (int it=0;it<dim1;it++){
      //        cout<<"it="<<it<<" :: ";
      //        for (int ir=0;ir<dim2;ir++){
      //          cout<<aTmp2D[it][ir]<<" ";
      //   }
      //        cout<<endl;
      // }
      // cout<<"NCGETVARS"<<endl;
    }

    double val;
    // -------------------------------
    // Copy all data from aTmp array to member array _aVal.
    // -------------------------------
    // \todo [optimize] - check values in 3D as done below for 2D (avoiding function call overhead)
    if ( _is_3D ) {
      int irow,icol;
      if (_dim_order == 1) {
        for (int it=0; it<iChunkSize; it++){                     // loop over time points in buffer
          for (int ic=0; ic<_nNonZeroWeightedGridCells; ic++){   // loop over non-zero weighted grid cells
            CellIdxToRowCol(_IdxNonZeroGridCells[ic],irow,icol);
            CheckValue3D(aTmp3D[icol][irow][it], missval, icol, irow, it);   // check if value to read in equals "missing_value"
            CheckValue3D(aTmp3D[icol][irow][it], fillval, icol, irow, it);   // check if value to read in equals "_FillValue"
            _aVal[it][ic]=_LinTrans_a*aTmp3D[icol][irow][it]+_LinTrans_b;
          }
        }
      }
      else if (_dim_order == 2) {
        for (int it=0; it<iChunkSize; it++){                     // loop over time points in buffer
          for (int ic=0; ic<_nNonZeroWeightedGridCells; ic++){   // loop over non-zero weighted grid cells
		        CellIdxToRowCol(_IdxNonZeroGridCells[ic],irow,icol);
		        CheckValue3D(aTmp3D[irow][icol][it], missval, irow, icol, it);   // check if value to read in equals "missing_value"
		        CheckValue3D(aTmp3D[irow][icol][it], fillval, irow, icol, it);   // check if value to read in equals "_FillValue"
		        _aVal[it][ic]=_LinTrans_a*aTmp3D[irow][icol][it]+_LinTrans_b;
		      }
		    }
      }
      else if (_dim_order == 3) {
        for (int it=0; it<iChunkSize; it++){                     // loop over time points in buffer
          for (int ic=0; ic<_nNonZeroWeightedGridCells; ic++){   // loop over non-zero weighted grid cells
            CellIdxToRowCol(_IdxNonZeroGridCells[ic],irow,icol);
            CheckValue3D(aTmp3D[icol][it][irow], missval, icol, it, irow);   // check if value to read in equals "missing_value"
	          CheckValue3D(aTmp3D[icol][it][irow], fillval, icol, it, irow);   // check if value to read in equals "_FillValue"
            _aVal[it][ic]=_LinTrans_a*aTmp3D[icol][it][irow]+_LinTrans_b;
          }
        }
      }
      else if (_dim_order == 4) {
        for (int it=0; it<iChunkSize; it++){                     // loop over time points in buffer
          for (int ic=0; ic<_nNonZeroWeightedGridCells; ic++){   // loop over non-zero weighted grid cells
            CellIdxToRowCol(_IdxNonZeroGridCells[ic],irow,icol);
            if (!((Options.deltaresFEWS) && (it==0))) {
              CheckValue3D(aTmp3D[it][icol][irow], missval, it, icol, irow);   // check if value to read in equals "missing_value"
              CheckValue3D(aTmp3D[it][icol][irow], fillval, it, icol, irow);   // check if value to read in equals "_FillValue"
              //cout<<aTmp3D[it][icol][irow]<<endl;
            }
            _aVal[it][ic]=_LinTrans_a*aTmp3D[it][icol][irow]+_LinTrans_b;
          }
        }
      }
      else if (_dim_order == 5) {
        for (int it=0; it<iChunkSize; it++){                     // loop over time points in buffer
          for (int ic=0; ic<_nNonZeroWeightedGridCells; ic++){   // loop over non-zero weighted grid cells
            CellIdxToRowCol(_IdxNonZeroGridCells[ic],irow,icol);
		        CheckValue3D(aTmp3D[irow][it][icol], missval, irow, it, icol);   // check if value to read in equals "missing_value"
		        CheckValue3D(aTmp3D[irow][it][icol], fillval, irow, it, icol);   // check if value to read in equals "_FillValue"
            _aVal[it][ic]=_LinTrans_a*aTmp3D[irow][it][icol]+_LinTrans_b;
          }
        }
      }
      else if (_dim_order == 6) {
        for (int it=0; it<iChunkSize; it++){                      // loop over time points in buffer
          for (int ic=0; ic<_nNonZeroWeightedGridCells; ic++){    // loop over non-zero weighted grid cells
            CellIdxToRowCol(_IdxNonZeroGridCells[ic],irow,icol);
		        CheckValue3D(aTmp3D[it][irow][icol], missval, it, irow, icol);   // check if value to read in equals "missing_value"
		        CheckValue3D(aTmp3D[it][irow][icol], fillval, it, irow, icol);   // check if value to read in equals "_FillValue"
            _aVal[it][ic]=_LinTrans_a*aTmp3D[it][irow][icol]+_LinTrans_b;
          }
        }
      }
    }
    else {

      if (_dim_order == 1) {
        for (int it=0; it<iChunkSize; it++){                     // loop over time points in buffer
          for (int ic=0; ic<_nNonZeroWeightedGridCells; ic++){   // loop over non-zero weighted grid cells
            val=aTmp2D[_IdxNonZeroGridCells[ic]][it];
            if(val==missval) { CheckValue2D(val,missval,_IdxNonZeroGridCells[ic],it); }   // throw error  if value to read in equals "missing_value"
            if(val==fillval) { CheckValue2D(val,fillval,_IdxNonZeroGridCells[ic],it); }   // throw error  if value to read in equals "_FillValue"
            _aVal[it][ic]=_LinTrans_a*val+_LinTrans_b;
          }
        }
      }
      else if (_dim_order == 2) {
        for (int it=0; it<iChunkSize; it++){                     // loop over time points in buffer
          for (int ic=0; ic<_nNonZeroWeightedGridCells; ic++){   // loop over non-zero weighted grid cells
            val=aTmp2D[it][_IdxNonZeroGridCells[ic]];

            if(val==missval) { CheckValue2D(val,missval,it,_IdxNonZeroGridCells[ic]); }  // throw error if value to read in equals "missing_value"
            if(val==fillval) { CheckValue2D(val,fillval,it,_IdxNonZeroGridCells[ic]); }  // throw error if value to read in equals "_FillValue"
            _aVal[it][ic]=_LinTrans_a*val+_LinTrans_b;
          }
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
    if ( _is_3D ) {
      for (int it=0;it<dim1;it++){delete [] aTmp3D[it];}
      delete [] aTmp3D;
    }
    else {
      delete [] aTmp2D;
    }
    delete [] aVec;

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
///
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

  _t_corr = -TimeDifference(model_start_day,model_start_year,_start_day,_start_year,Options.calendar);
 
  //QA/QC: Check for overlap issues between time series duration and model duration
  //------------------------------------------------------------------------------
  double duration = (double)(_nPulses)*_interval;
  double local_simulation_start = (_t_corr);
  double local_simulation_end = (model_duration + _t_corr);
  if((local_simulation_start<0) && (local_simulation_start>-TIME_CORRECTION)){ local_simulation_start=0.0; }

  if(Options.deltaresFEWS) { local_simulation_start-=_interval; }//first timestep in FEWS netCDF output is noise

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
void CForcingGrid::SetIdxNonZeroGridCells(const int nHydroUnits, const int nGridCells)
{
  bool *nonzero;

  nonzero = NULL;
  nonzero =  new bool [nGridCells]; // \todo [optimization]: may wish to declare as static
  for (int il=0; il<nGridCells; il++) { // loop over all cells of NetCDF
    nonzero[il] = false;
  }

  if (_GridWeight != NULL){
    for(int k=0; k<nHydroUnits; k++) {  // loop over HRUs
      for(int i=0; i<_nWeights[k]; i++) { // loop over all cells of NetCDF
        if(_GridWeight[k][i] > 0.00001) {
          nonzero[_GridWtCellIDs[k][i]] = true;
        }
      }
    }
  }
  else{
    ExitGracefully(
      "CForcingGrid: SetIdxNonZeroGridCells: _GridWeight is not allocated yet. Call AllocateWeightArray(nHRUs) first.", BAD_DATA);
  }

  // count number of non-zero weighted grid cells
  _nNonZeroWeightedGridCells=0;
  for (int il=0; il<nGridCells; il++) { // loop over all cells of NetCDF
    if ( nonzero[il] ) { _nNonZeroWeightedGridCells++; }
  }
  //printf("ForcingGrid: nNonZeroGridCells = %i\n",nNonZeroGridCells);

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
  _steps_per_day=(int)(round(1.0/_interval));
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
  for(int i=0;i<_nWeights[k]-1;i++) {
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
///
/// \param stat_ID [in] cell ID/stationID in NetCDF (from 0 to _nNonZeroWeightedGridCells-1)
/// \param weight [in] elevation of cell in NetCDF

void   CForcingGrid::SetStationElevation(const int CellID,const double &elev) 
{
  int ncells;
  if(_is_3D) { ncells = _GridDims[0] * _GridDims[1]; }
  else       { ncells = _GridDims[0]; }

  if((CellID<0) || (CellID>=ncells)) {
    ExitGracefully("CForcingGrid: SetStationElevation: invalid cell/station identifier (likely in :StationElevationsByAttribute command)",BAD_DATA);
  }
  if(_aElevation==NULL) {
    _aElevation=new double[_nNonZeroWeightedGridCells];
    for(int i=0;i<_nNonZeroWeightedGridCells;i++) { _aElevation[i]=0.0; }
  }
  
  _aElevation[_CellIDToIdx[CellID]]=elev;
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
      if(fabs(sum_HRU - 1.0) > 0.0001) {
        cout<<"HRU ID = "<<k<<" Sum Forcing Weights = "<<sum_HRU<<endl;
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
  //  \todo[optimize] - could we send in i instead?
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
  //JRC_TIME_FIX
  return int(( t) * round(1.0/_interval)+0.5*tstep)  % _ChunkSize;
  //return int((t+0.25*tstep)*_steps_per_day)  % _ChunkSize; //JRC \todo[fix] - I believe this is the correct approach

  //return int((_t_corr + t) * round(1.0/_interval)+0.5*tstep)  % _ChunkSize;
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
  int nSteps = max(1,(int)(round(tstep/_interval)));//# of intervals in time step
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

  int nSteps = (int)(max(1.0,round(tstep/_interval)));//# of intervals in time step
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
bool CForcingGrid::ShouldDeaccumulate()      const{return _deaccumulate;}

//////////////////////////////////////////////////////////////////
/// \brief Returns forcing type
/// \return forcing type
//
forcing_type CForcingGrid::GetForcingType  () const{return _ForcingType;}
