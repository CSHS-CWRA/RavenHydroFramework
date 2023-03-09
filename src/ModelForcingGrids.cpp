/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2019 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Model.h"

/*****************************************************************
   Routines for generating forcing grids from other forcing grids
   All member functions of CModel
   -GenerateAveSubdailyTempFromMinMax
   -GenerateMinMaxAveTempFromSubdaily
   -GenerateMinMaxSubdailyTempFromAve
   -GeneratePrecipFromSnowRain
   -GetAverageSnowFrac
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief if forcing grid of type typ doesnt exist, creates full copy of pGrid 
///    but with potentially new time interval specified
///    otherwise just returns existing grid of type typ
///    assume
//
CForcingGrid *CModel::ForcingCopyCreate(const CForcingGrid *pGrid,
                                        const forcing_type typ,
                                        const double &interval, 
                                        const int nVals, const optStruct &Options)
{
  static CForcingGrid *pTout;
  if (GetForcingGridIndexFromType(typ) == DOESNT_EXIST )
  { // for the first chunk, the derived grid does not exist and has to be added to the model
    // all weights, etc., are copied from the base grid 

    pTout = new CForcingGrid(*pGrid);  // copy everything from pGrid; matrices are deep copies

    //following values are overwritten:
    int    GridDims[3];
    GridDims[0] = pGrid->GetCols();
    GridDims[1] = pGrid->GetRows();
    GridDims[2] = nVals;
    pTout->SetForcingType(typ);
    pTout->SetInterval(interval);            
    pTout->SetGridDims(GridDims);

    pTout->SetChunkSize(nVals); 
    //pTout->CalculateChunkSize(Options);
    pTout->ReallocateArraysInForcingGrid();
  }
  else
  {
    // for all latter chunks the the grid already exists and values will be just overwritten
    pTout=GetForcingGrid(typ);
  }

  //int nGridHRUs=pGrid->GetnHydroUnits();
  //int nCells   =pGrid->GetRows()*pGrid->GetCols(); //handles 3D or 2D
  //cout<<"**ForcingGrid Copy Create : "<<ForcingToString(pGrid->GetForcingType())<<"-->"<<ForcingToString(pTout->GetForcingType())<< " is new?:"<<is_new<<"**"<<endl;
  
  return pTout;
}

//////////////////////////////////////////////////////////////////
/// \brief Creates all missing gridded snow/rain/precip data based on gridded information available,
///        precip data are assumed to have the same resolution and hence can be initialized together.
///
/// \param Options [in]  major options of the model
//
void CModel::GenerateGriddedPrecipVars(const optStruct &Options)
{
  CForcingGrid * pGrid_pre  = NULL;

  // see if gridded forcing is read from a NetCDF
  bool pre_gridded   = ForcingGridIsInput(F_PRECIP);
  bool rain_gridded  = ForcingGridIsInput(F_RAINFALL);
  bool snow_gridded  = ForcingGridIsInput(F_SNOWFALL);

  // find the correct grid
  if(pre_gridded)  { pGrid_pre   = GetForcingGrid((F_PRECIP)); }

  // Minimum requirements of forcing grids: must have precip or rain
  ExitGracefullyIf(!pre_gridded && !rain_gridded,"CModel::InitializeForcingGrids: No precipitation forcing found",BAD_DATA);

  if(Options.noisy) { cout<<"Generating gridded precipitation variables..."<<endl; }
  //--------------------------------------------------------------------------
  // Handle snowfall availability
  //--------------------------------------------------------------------------
  if(snow_gridded && rain_gridded && (Options.rainsnow!=RAINSNOW_DATA))
  {
    WriteWarning("CModel::GenerateGriddedPrecipVars: both snowfall and rainfall data are provided at a gauge, but :RainSnowFraction method is something other than RAINSNOW_DATA. Snow fraction will be recalculated.",Options.noisy);
  }
  //deaccumulate if necessary
  double rainfall_rate;
  if((pre_gridded) && (pGrid_pre->ShouldDeaccumulate()))
  {
    for(int it=0; it<pGrid_pre->GetChunkSize()-1; it++) {                   // loop over time points in buffer
      for(int ic=0; ic<pGrid_pre->GetNumberNonZeroGridCells(); ic++) {       // loop over non-zero grid cell indexes
        rainfall_rate=(pGrid_pre->GetValue(ic,it+1)-pGrid_pre->GetValue(ic,it))/pGrid_pre->GetInterval();
        pGrid_pre->SetValue(ic,it,rainfall_rate);   // copies precipitation values
      }
    }
    for(int ic=0; ic<pGrid_pre->GetNumberNonZeroGridCells(); ic++) {       // loop over non-zero grid cell indexes
      rainfall_rate=0;
      pGrid_pre->SetValue(ic,pGrid_pre->GetChunkSize()-1,rainfall_rate);
    }
    //pGrid_pre->SetChunkSize(pGrid_pre->GetChunkSize()-1);
  }
  if(Options.noisy) { cout<<"SNOW="<<snow_gridded<<" RAIN="<<rain_gridded<<" PRECIP="<<pre_gridded<<endl; }
  if(snow_gridded && rain_gridded && !pre_gridded) {
    GeneratePrecipFromSnowRain(Options);
  }
  if(!rain_gridded && !snow_gridded && pre_gridded) {
    if(Options.noisy) { cout<<"Generating rain grid"<<endl; }
    GenerateRainFromPrecip(Options);
  }
  if(!snow_gridded && (pre_gridded || rain_gridded)) {
    if(Options.noisy) { cout<<"Generating empty snow grid"<<endl; }
    GenerateZeroSnow(Options);
    if(!pre_gridded) { 
      if(Options.noisy) { cout<<"Generating precip grid"<<endl; }
      GeneratePrecipFromSnowRain(Options); 
    }
  }
  if(Options.noisy) { cout<<"...done generating gridded precipitation variables."<<endl; }
}

//////////////////////////////////////////////////////////////////
/// \brief Creates all missing gridded temperature data based on gridded information available,
///        e.g when only sub-daily temperature is provided estimate daily average, minimum and maximum temperature.
///        data are assumed to have the same resolution and hence can be initialized together.
///
/// \param Options [in]  major options of the model
//
void CModel::GenerateGriddedTempVars(const optStruct &Options)
{
  CForcingGrid * pGrid_tave       = NULL;
  CForcingGrid * pGrid_daily_tave = NULL;

  // see if gridded forcing is read from a NetCDF
  bool temp_ave_gridded       = ForcingGridIsInput(F_TEMP_AVE);
  bool temp_daily_min_gridded = ForcingGridIsInput(F_TEMP_DAILY_MIN);
  bool temp_daily_max_gridded = ForcingGridIsInput(F_TEMP_DAILY_MAX);
  bool temp_daily_ave_gridded = ForcingGridIsInput(F_TEMP_DAILY_AVE);

  // find the correct grid
  if(temp_ave_gridded      ) { pGrid_tave       = GetForcingGrid((F_TEMP_AVE)); }
  if(temp_daily_ave_gridded) { pGrid_daily_tave = GetForcingGrid((F_TEMP_DAILY_AVE)); }

  // Temperature min/max grids must be both available
  ExitGracefullyIf((temp_daily_min_gridded && !temp_daily_max_gridded) || (!temp_daily_min_gridded && temp_daily_max_gridded),
    "CModel::UpdateHRUForcingFunctions: Input minimum and maximum temperature have to be given either both as time series or both as gridded input.",BAD_DATA);

  if(Options.noisy) { cout<<"Generating gridded temperature variables..."<<endl; }

  //--------------------------------------------------------------------------
  // Populate Temperature forcing grid: by timestep, daily min/max/average
  //--------------------------------------------------------------------------
  string tmp[3];
  tmp[0] = "NONE";  tmp[1] = "NONE";  tmp[2] = "NONE";

  if ((temp_ave_gridded) && (fabs(pGrid_tave->GetInterval()-1.0)<REAL_SMALL)) // (A) daily temperature data provided, labeled as TEMP_AVE - >convert to temp_daily_ave_gridded
  {

    pGrid_daily_tave = ForcingCopyCreate(pGrid_tave,F_TEMP_DAILY_AVE,pGrid_tave->GetInterval(),pGrid_tave->GetChunkSize(),Options);
    int chunk_size =pGrid_tave->GetChunkSize();
    int nNonZero   =pGrid_tave->GetNumberNonZeroGridCells();
    for(int it=0; it<chunk_size; it++) {           // loop over time points in buffer
      for(int ic=0; ic<nNonZero;ic++) {             // loop over non-zero grid cell indexes
        pGrid_daily_tave->SetValue(ic,it,pGrid_tave->GetValue(ic,it)); //copy values //\todo[optimize] - is this necessary - should copy in ForcingCopyCreate
      }
    }
    AddForcingGrid(pGrid_daily_tave,F_TEMP_DAILY_AVE);
    temp_ave_gridded=false;
    temp_daily_ave_gridded=true;
    pGrid_daily_tave->SetIsDerived(false); //treat as if it were input directly
  }
  if(temp_ave_gridded)                                                // (B) Sub-daily temperature data provided
  {
    if(Options.noisy) { cout<<"Generating min/max/ave temp"<<endl; }
    GenerateMinMaxAveTempFromSubdaily(Options);                        // ---> Generate daily min, max, & avg
  }
  else if(temp_daily_min_gridded && temp_daily_max_gridded)            // (C) Daily min/max temperature data provided
  {
    if(Options.noisy) { cout<<"Generating ave/hourly temp from min/max"<<endl; }
    GenerateAveSubdailyTempFromMinMax(Options);                        // --> Generate T_daily_ave, T_ave (downscaled)
  }
  else if(temp_daily_ave_gridded)                                      // (D) only daily average data provided
  {
    if(Options.noisy) { cout<<"Generating min/max temp from ave"<<endl; }
    GenerateMinMaxSubdailyTempFromAve(Options);                        // --> Generate daily T_min, T_max, and subdaily ave (downscaled)
  }
  else
  {
    ExitGracefully("CModel::GenerateGriddedTempVars: Insufficient data to generate temperature forcing grids",BAD_DATA);
  }

  ExitGracefullyIf(GetForcingGrid(F_TEMP_DAILY_MAX)->GetInterval() != GetForcingGrid(F_TEMP_DAILY_MIN)->GetInterval(),
    "CModel::InitializeForcingGrids: Input minimum and maximum temperature must have same time resolution.",BAD_DATA);

  if(Options.noisy) { cout<<"...done generating gridded temperature variables."<<endl; }
}
//////////////////////////////////////////////////////////////////
/// \brief Generates Tave and subhourly time series from daily Tmin & Tmax time series
/// \note presumes existence of valid F_TEMP_DAILY_MIN and F_TEMP_DAILY_MAX time series
//
void CModel::GenerateAveSubdailyTempFromMinMax(const optStruct &Options)
{
  CForcingGrid *pTmin,*pTmax,*pTave,*pTave_daily;
  pTmin=GetForcingGrid(F_TEMP_DAILY_MIN);  
  pTmax=GetForcingGrid(F_TEMP_DAILY_MAX); 

  pTmin->Initialize(Options);// needed for correct mapping from time series to model time
  pTmax->Initialize(Options);

  int nValsPerDay= (int)(1.0 / pTmin->GetInterval()); //typically should be 1.0

  // ----------------------------------------------------
  // Generate daily average grid, if it doesnt exist
  // ----------------------------------------------------
  if(!ForcingGridIsInput(F_TEMP_DAILY_AVE)) 
  {
    int    nVals     = (int)ceil(pTmin->GetChunkSize() * pTmin->GetInterval());
    double Tave;
    pTave_daily = ForcingCopyCreate(pTmin,F_TEMP_DAILY_AVE,1.0,nVals,Options);

    double t=0.0;
    int chunk_size =pTave_daily->GetChunkSize();
    int nNonZero   =pTave_daily->GetNumberNonZeroGridCells();
    for(int it=0; it<chunk_size; it++) {            // loop over time points in buffer
      for(int ic=0; ic<nNonZero;ic++) {             // loop over non-zero grid cell indexes
        Tave=0.5*(pTmin->GetValue_avg(ic,t * nValsPerDay,nValsPerDay) +
                  pTmax->GetValue_avg(ic,t * nValsPerDay,nValsPerDay));
        pTave_daily->SetValue(ic,it,Tave);
      }
      t+=1.0;
    }
    AddForcingGrid(pTave_daily,F_TEMP_DAILY_AVE);
  }
  else {
    pTave_daily=GetForcingGrid(F_TEMP_DAILY_AVE);
  }

  // ----------------------------------------------------
  // Generate subdaily temperature values
  // ----------------------------------------------------
  if (Options.timestep<(1.0-TIME_CORRECTION)) //tstep < 1 day
  {
    int    nVals     = (int)ceil(pTave_daily->GetChunkSize()/Options.timestep);

    pTave = ForcingCopyCreate(pTmin,F_TEMP_AVE,Options.timestep,nVals,Options);

    // set forcing values
    // Tmax       is with input time resolution
    // Tmin       is with input time resolution
    // Tave       is with model time resolution
    // Tave_daily is with daily resolution
    double t=0.0; // model time
    double Tmax,Tmin,T1corr,T2corr,val;
    int time_idx_chunk;
    int nNonZero  =pTave->GetNumberNonZeroGridCells();
    double time_shift=Options.julian_start_day-floor(Options.julian_start_day+TIME_CORRECTION);
    for (int it=0; it<nVals; it++) {                   // loop over all time points (nVals)
      for (int ic=0; ic<nNonZero; ic++){               // loop over non-zero grid cell indexes
        time_idx_chunk = int(floor(t+TIME_CORRECTION));
        Tmin   = pTmin->GetValue_avg(ic, floor(t +time_shift+TIME_CORRECTION)*nValsPerDay, nValsPerDay);
        Tmax   = pTmax->GetValue_avg(ic, floor(t +time_shift+TIME_CORRECTION)*nValsPerDay, nValsPerDay);
        T1corr = pTave->DailyTempCorrection(t+time_shift);
        T2corr = pTave->DailyTempCorrection(t+time_shift+Options.timestep); 
        val=pTave_daily->GetValue(ic, time_idx_chunk)+0.25*(Tmax-Tmin)*(T1corr+T2corr);
        pTave->SetValue( ic, it, val);
      }
      t += Options.timestep;
    }
    AddForcingGrid(pTave,F_TEMP_AVE);
  }
  else //tstep ==  1 day
  {
    if(!ForcingGridIsInput(F_TEMP_AVE))
    {
      int    nVals     = pTave_daily->GetChunkSize();
      pTave = ForcingCopyCreate(pTave_daily,F_TEMP_AVE,1.0,nVals,Options);

      int nNonZero  =pTave->GetNumberNonZeroGridCells();
      for(int it=0; it<nVals; it++) {                                 // loop over time points in buffer
        for(int ic=0; ic<nNonZero; ic++) {                                  // loop over non-zero grid cell indexes
          pTave->SetValue(ic,it,pTave_daily->GetValue(ic,it));  // --> just copy daily average values
        }
      }
      AddForcingGrid(pTave,F_TEMP_AVE);
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Generates daily Tmin,Tmax,Tave time series grids from T (subdaily) time series grid
/// \note presumes existence of valid F_TEMP_AVE time series with subdaily timestep
//
void CModel::GenerateMinMaxAveTempFromSubdaily(const optStruct &Options)
{
  CForcingGrid *pTave,*pTmin_daily,*pTmax_daily,*pTave_daily;

  pTave=GetForcingGrid(F_TEMP_AVE);
  pTave->Initialize(Options);  // needed for correct mapping from time series to model time
 
  double interval = pTave->GetInterval();
  int    nVals    = (int)ceil(pTave->GetChunkSize()*interval); // number of daily values

  pTmin_daily = ForcingCopyCreate(pTave,F_TEMP_DAILY_MIN,1.0,nVals,Options);
  pTmax_daily = ForcingCopyCreate(pTave,F_TEMP_DAILY_MAX,1.0,nVals,Options);
  pTave_daily = ForcingCopyCreate(pTave,F_TEMP_DAILY_AVE,1.0,nVals,Options);

  double time;
  double time_shift=Options.julian_start_day-floor(Options.julian_start_day+TIME_CORRECTION);
  int    nsteps_in_day=int(1.0/interval);

  for (int it=0; it<nVals; it++) {                    // loop over time points in buffer
    time=(double)it*1.0/interval;
    time=floor(time-time_shift+TIME_CORRECTION); //model time corresponding to 00:00 on day of it

    for (int ic=0; ic<pTave->GetNumberNonZeroGridCells(); ic++){         // loop over non-zero grid cell indexes
      pTmin_daily->SetValue(ic, it, pTave->GetValue_min(ic, time, nsteps_in_day));
      pTmax_daily->SetValue(ic, it, pTave->GetValue_max(ic, time, nsteps_in_day));
      pTave_daily->SetValue(ic, it, pTave->GetValue_avg(ic, time, nsteps_in_day));
    }
  }

  AddForcingGrid(pTmin_daily,F_TEMP_DAILY_MIN);
  AddForcingGrid(pTmax_daily,F_TEMP_DAILY_MAX);
  AddForcingGrid(pTave_daily,F_TEMP_DAILY_AVE);
}

//////////////////////////////////////////////////////////////////
/// \brief Generates Tmin, Tmax and subhourly time series grids from daily average temperature time series grid
/// \note presumes existence of valid F_TEMP_DAILY_AVE
/// \note necessarily naive - it is hard to downscale with little temp data
//
void CModel::GenerateMinMaxSubdailyTempFromAve(const optStruct &Options)
{
  CForcingGrid *pTmin_daily,*pTmax_daily,*pTave_daily;

  pTave_daily=GetForcingGrid(F_TEMP_DAILY_AVE);

  pTave_daily->Initialize(Options);//  needed for correct mapping from time series to model time

  double interval = pTave_daily->GetInterval();
  int       nVals = pTave_daily->GetChunkSize(); // number of subdaily values (input resolution) - should be 1
  
  pTmin_daily = ForcingCopyCreate(pTave_daily,F_TEMP_DAILY_MIN,interval,nVals,Options);
  pTmax_daily = ForcingCopyCreate(pTave_daily,F_TEMP_DAILY_MAX,interval,nVals,Options);
  int chunksize=pTave_daily->GetChunkSize();
  int nNonZero =pTave_daily->GetNumberNonZeroGridCells();
  double daily_temp;
  for (int it=0; it<nVals; it++) {          // loop over time points in buffer
    for (int ic=0; ic<nNonZero; ic++){      // loop over non-zero grid cell indexes
      daily_temp=pTave_daily->GetValue(ic, min(chunksize,it));// should be it+0.5
      pTmin_daily->SetValue(ic, it, daily_temp-4.0); 
      pTmax_daily->SetValue(ic, it, daily_temp+4.0); 
    }
  }
  AddForcingGrid(pTmin_daily,F_TEMP_DAILY_MIN);
  AddForcingGrid(pTmax_daily,F_TEMP_DAILY_MAX);

  // ----------------------------------------------------
  // Generate subdaily averages from daily values (min, max)
  // ----------------------------------------------------
  GenerateAveSubdailyTempFromMinMax(Options);
}

////////////////////////////////////////////////////////////////// 
/// \brief Generates precipitation grid as sum of snowfall and rainfall grids
/// \note  presumes existence of valid F_SNOWFALL and F_RAINFALL time series
//
void CModel::GeneratePrecipFromSnowRain(const optStruct &Options)
{

  CForcingGrid *pPre,*pSnow,*pRain;
  pSnow=GetForcingGrid((F_SNOWFALL));
  pRain=GetForcingGrid((F_RAINFALL));

  //below needed for correct mapping from time series to model time
  pSnow->Initialize(Options);
  pRain->Initialize(Options);

  double interval_snow = pSnow->GetInterval();
  double interval_rain = pRain->GetInterval();

  ExitGracefullyIf(interval_snow != interval_rain,
                   "CModel::GeneratePrecipFromSnowRain: rainfall and snowfall must have the same time resolution!",BAD_DATA);

  int nVals = pSnow->GetChunkSize();

  pPre = ForcingCopyCreate(pSnow,F_PRECIP,interval_snow,nVals,Options);

  int nNonZero  =pPre->GetNumberNonZeroGridCells();
  for (int it=0; it<nVals; it++) {     // loop over time points in buffer
    for (int ic=0; ic<nNonZero; ic++){      // loop over non-zero grid cell indexes
      pPre->SetValue(ic, it, pSnow->GetValue(ic, it) + pRain->GetValue(ic, it));  // precipitation = sum of snowfall and rainfall
    }
  }

  AddForcingGrid(pPre,F_PRECIP);
}

//////////////////////////////////////////////////////////////////
/// \brief Generates rainfall grid as copy of precipitation grid 
/// \note  presumes existence of valid F_PRECIP time series
//
void CModel::GenerateRainFromPrecip(const optStruct &Options)
{
  CForcingGrid *pPre,*pRain;
  pPre=GetForcingGrid((F_PRECIP));

  pPre->Initialize(Options);  //needed for correct mapping from time series to model time

  int nVals = pPre->GetChunkSize();

  pRain = ForcingCopyCreate(pPre,F_RAINFALL,pPre->GetInterval(),nVals,Options);

  // set forcing values
  int nNonZero  =pRain->GetNumberNonZeroGridCells();
  for(int it=0; it<nVals; it++) {                   // loop over time points in buffer
    for(int ic=0; ic<nNonZero; ic++){                    // loop over non-zero grid cell indexes
      pRain->SetValue(ic,it,pPre->GetValue(ic,it));      // copies precipitation values
    }
  }

  AddForcingGrid(pRain,F_RAINFALL);
}

//////////////////////////////////////////////////////////////////
/// \brief Generates snowfall time series grid being constantly zero
/// \note  presumes existence of either rainfall or snowfall
//
void CModel::GenerateZeroSnow(const optStruct &Options)
{
  CForcingGrid *pPre(NULL),*pSnow;
  if      (ForcingGridIsAvailable(F_PRECIP))   { pPre=GetForcingGrid((F_PRECIP)); }
  else if (ForcingGridIsAvailable(F_RAINFALL)) { pPre=GetForcingGrid((F_RAINFALL)); }
  else {ExitGracefully("CModel:GenerateZeroSnow: missing precip or rainfall grid",RUNTIME_ERR); return;}

  pPre->Initialize(Options);//needed for correct mapping from time series to model time

  int    nVals     = pPre->GetChunkSize();

  pSnow = ForcingCopyCreate(pPre,F_SNOWFALL,pPre->GetInterval(),nVals,Options);

  int nNonZero  =pSnow->GetNumberNonZeroGridCells();
  for (int it=0; it<nVals; it++) {                   // loop over time points in buffer
    for (int ic=0; ic<nNonZero; ic++){                    // loop over non-zero grid cell indexes
      pSnow->SetValue(ic, it , 0.0);                      // fills everything with 0.0
    }
  }

  AddForcingGrid(pSnow,F_SNOWFALL);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns average fraction of snow in precipitation between time t and following n timesteps
/// \param x_col  [in] Column index
/// \param y_row  [in] Row index
/// \param t      [in] Time index
/// \param n      [in] Number of time steps
/// \return average fraction of snow in precipitation between time t and following n timesteps
//
double CModel::GetAverageSnowFrac(const int ic, const double t, const int n) const
{

  CForcingGrid *pSnow,*pRain;
  pSnow=GetForcingGrid((F_SNOWFALL));
  pRain=GetForcingGrid((F_RAINFALL));

  double snow = pSnow->GetValue_avg(ic, t, n);
  double rain = pRain->GetValue_avg(ic, t, n);

  if ((snow+rain)==0.0){return 0.0;}
  return snow/(snow+rain);

}






