/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
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
/// \brief Generates Tave and subhourly time series from daily Tmin & Tmax time series
/// \note presumes existence of valid F_TEMP_DAILY_MIN and F_TEMP_DAILY_MAX time series
//
void CModel::GenerateAveSubdailyTempFromMinMax(const optStruct &Options)
{
  CForcingGrid *pTmin,*pTmax,*pTave,*pTave_daily;
  pTmin=GetForcingGrid(F_TEMP_DAILY_MIN);  
  pTmax=GetForcingGrid(F_TEMP_DAILY_MAX); 

  double start_day=Options.julian_start_day;
  int    start_yr =Options.julian_start_year;
  double duration =Options.duration;
  double timestep =Options.timestep;

  //below needed for correct mapping from time series to model time
  pTmin->Initialize(start_day,start_yr,duration,timestep,Options);
  pTmax->Initialize(start_day,start_yr,duration,timestep,Options);

  int    nVals     = (int)ceil(pTmin->GetChunkSize() * pTmin->GetInterval());
  int    GridDims[3];
  GridDims[0] = pTmin->GetCols(); GridDims[1] = pTmin->GetRows(); GridDims[2] = nVals;

  // ----------------------------------------------------
  // Generate daily average values Tave=(Tmin+Tmax)/2
  // --> This is always a daily time series (also if TEMP_DAILY_MIN and TEMP_DAILY_MAX are subdaily)
  // ----------------------------------------------------
  if ( GetForcingGridIndexFromType(F_TEMP_DAILY_AVE) == DOESNT_EXIST )
  {
    // for the first chunk the derived grid does not exist and has to be added to the model
    pTave_daily = new CForcingGrid(* pTmin);  // copy everything from tmin; matrixes are deep copies
    pTave_daily->SetForcingType(F_TEMP_DAILY_AVE);
    pTave_daily->SetInterval(1.0);        // always daily
    pTave_daily->SetGridDims(GridDims);
    pTave_daily->SetChunkSize(nVals);     // if Tmin/Tmax are subdaily, several timesteps might be merged to one
    pTave_daily->ReallocateArraysInForcingGrid();
  }
  else
  {
    // for all latter chunks the the grid already exists and values will be just overwritten
    pTave_daily=GetForcingGrid(F_TEMP_DAILY_AVE);
  }

  // (1) set weighting
  int nGridHRUs=pTave_daily->GetnHydroUnits();
  int nCells=pTave_daily->GetRows()*pTave_daily->GetCols();
  for (int ik=0; ik<nGridHRUs; ik++) {                // loop over HRUs
    for (int ic=0; ic<nCells; ic++) {                 // loop over cells = rows*cols
      pTave_daily->SetWeightVal(ik, ic, pTmin->GetGridWeight(ik, ic));
    }
  }

  // (2) set indexes of on-zero weighted grid cells
  pTave_daily->SetIdxNonZeroGridCells(pTave_daily->GetnHydroUnits(),pTave_daily->GetRows()*pTave_daily->GetCols());

  // (3) set forcing values
  double t=0.0;

  int chunk_size=pTave_daily->GetChunkSize();
  int nNonZero  =pTave_daily->GetNumberNonZeroGridCells();
  for (int it=0; it<chunk_size; it++) {           // loop over time points in buffer
    for (int ic=0; ic<nNonZero;ic++){             // loop over non-zero grid cell indexes
      // found in Gauge.cpp: CGauge::GenerateAveSubdailyTempFromMinMax(const optStruct &Options)
      //    double t=0.0;//model time
      //    for (int n=0;n<nVals;n++)
      //      {
      //        aAvg[n]=0.5*(pTmin->GetValue(t+0.5)+pTmax->GetValue(t+0.5));
      //        t+=1.0;
      //      }
      // TODO: it --> should be it+0.5

      //double time_idx_chunk = pTmax->GetChunkIndexFromModelTimeStep(Options,t+0.5);
      int    nValsPerDay    = (int)(1.0 / pTmin->GetInterval());
      double time_idx_chunk = t * nValsPerDay;
      pTave_daily->SetValue(ic, it , 0.5*(pTmin->GetValue(ic, time_idx_chunk, nValsPerDay) + pTmax->GetValue(ic, time_idx_chunk, nValsPerDay)));
    }
    t+=1.0;
  }

  if ( GetForcingGridIndexFromType(F_TEMP_DAILY_AVE) == DOESNT_EXIST ) {
    this->AddForcingGrid(pTave_daily);
    if (Options.noisy){ printf("\n------------------------> TEMP_DAILY_AVE Added \n"); }
  }
  else {
    if (Options.noisy){ printf("\n------------------------> TEMP_DAILY_AVE Replace \n"); }
  }

  // ----------------------------------------------------
  // Generate subdaily temperature values
  // ----------------------------------------------------
  if (Options.timestep<(1.0-TIME_CORRECTION))
  {
    int    nVals     = (int)ceil(pTave_daily->GetChunkSize()/Options.timestep);
    int    GridDims[3];
    GridDims[0] = pTmin->GetCols(); GridDims[1] = pTmin->GetRows(); GridDims[2] = nVals;

    // Generate subdaily average values
    if ( GetForcingGridIndexFromType(F_TEMP_AVE) == DOESNT_EXIST )
    {
      // for the first chunk the derived grid does not exist and has to be added to the model
      pTave = new CForcingGrid(* pTave_daily);  // copy everything from tmin; matrixes are deep copies
      pTave->SetForcingType(F_TEMP_AVE);
      pTave->SetInterval(Options.timestep);  // is always model time step; no matter which _interval Tmin/Tmax had
      pTave->SetGridDims(GridDims);
      pTave->SetChunkSize(nVals);
      pTave->ReallocateArraysInForcingGrid();
    }
    else
    {
      // for all latter chunks the the grid already exists and values will be just overwritten
      pTave=GetForcingGrid((F_TEMP_AVE));
    }

    // (1) set weighting
    int nGridHRUs=pTmin->GetnHydroUnits();
    int nCells=pTmin->GetRows()*pTmin->GetCols();
    for (int ik=0; ik<nGridHRUs; ik++) {              // loop over HRUs
      for (int ic=0; ic<nCells; ic++) {               // loop over cells = rows*cols
        pTave->SetWeightVal(ik, ic, pTmin->GetGridWeight(ik, ic));
      }
    }

    // (2) set indexes of on-zero weighted grid cells
    pTave->SetIdxNonZeroGridCells(pTave->GetnHydroUnits(),pTave->GetRows()*pTave->GetCols());

    // (3) set forcing values
    // Tmax       is with input time resolution
    // Tmin       is with input time resolution
    // Tave       is with model time resolution
    // Tave_daily is with daily resolution
    double t=0.0; // model time
    double time_idx_chunk,Tmax,Tmin,T1corr,T2corr,val;
    int nNonZeroWeightedGridCells = pTave->GetNumberNonZeroGridCells();
    for (int it=0; it<GridDims[2]; it++) {                   // loop over all time points (nVals)
      for (int ic=0; ic<nNonZeroWeightedGridCells; ic++){    // loop over non-zero grid cell indexes
        time_idx_chunk = double(int((t+Options.timestep/2.0)/pTmin->GetInterval()));
        Tmax   = pTmax->GetValue(ic, time_idx_chunk);
        Tmin   = pTmin->GetValue(ic, time_idx_chunk);
        T1corr = pTave->DailyTempCorrection(t);
        T2corr = pTave->DailyTempCorrection(t+Options.timestep);
        val    = 0.5*(Tmax+Tmin)+0.5*(Tmax-Tmin)*0.5*(T1corr+T2corr);
        pTave->SetValue( ic, it, val);
      }
      t += Options.timestep;
    }

    if ( GetForcingGridIndexFromType(F_TEMP_AVE) == DOESNT_EXIST ) {
      this->AddForcingGrid(pTave);
      if (Options.noisy){ printf("\n------------------------> TEMP_AVE case 1 Added \n"); }
    }
    else {
      if (Options.noisy){ printf("\n------------------------> TEMP_AVE case 1 Replace \n"); }
    }
  }
  else
  {
    // Tmax       is with input time resolution
    // Tmin       is with input time resolution
    // Tave       is with model time resolution
    // Tave_daily is with daily resolution

    // model does not run with subdaily time step
    // --> just copy daily average values
    if ( GetForcingGridIndexFromType(F_TEMP_AVE) == DOESNT_EXIST )
    {
      // for the first chunk the derived grid does not exist and has to be added to the model
      pTave = new CForcingGrid(* pTave_daily);  // copy everything from tave; matrixes are deep copies
      pTave->SetForcingType(F_TEMP_AVE);
    }
    else
    {
      // for all latter chunks the the grid already exists and values will be just overwritten
      pTave = GetForcingGrid((F_TEMP_AVE));
    }

    // (1) set weighting
    int nGridHRUs=pTave->GetnHydroUnits();
    int nCells=pTave->GetRows()*pTave->GetCols();
    for (int ik=0; ik<nGridHRUs; ik++) {              
      for (int ic=0; ic<nCells; ic++) {               
        pTave->SetWeightVal(ik, ic, pTave_daily->GetGridWeight(ik, ic)); // --> just copy daily average values
      }
    }

    // (2) set indexes of on-zero weighted grid cells
    pTave->SetIdxNonZeroGridCells(pTave->GetnHydroUnits(),pTave->GetRows()*pTave->GetCols());

    // (3) set forcing values
    int chunk_size=pTave->GetChunkSize();
    int nNonZero  =pTave->GetNumberNonZeroGridCells();
    for (int it=0; it<chunk_size; it++) {                                 // loop over time points in buffer
      for (int ic=0; ic<nNonZero; ic++){                                  // loop over non-zero grid cell indexes
        pTave->SetValue(ic, it , pTave_daily->GetValue(ic, (double)it));  // --> just copy daily average values
      }
    }

    if ( GetForcingGridIndexFromType(F_TEMP_AVE) == DOESNT_EXIST ) {
      this->AddForcingGrid(pTave);
      if (Options.noisy){ printf("\n------------------------> TEMP_AVE case 2 Added \n"); }
    }
    else {
      if (Options.noisy){ printf("\n------------------------> TEMP_AVE case 2 Replace \n"); }
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Generates daily Tmin,Tmax,Tave time series from T (subdaily) time series
/// \note presumes existence of valid F_TEMP_AVE time series with subdaily timestep
//
void CModel::GenerateMinMaxAveTempFromSubdaily(const optStruct &Options)
{

  CForcingGrid *pTave,*pTmin_daily,*pTmax_daily,*pTave_daily;
  double interval;

  pTave=GetForcingGrid((F_TEMP_AVE));
  interval = pTave->GetInterval();

  double start_day = Options.julian_start_day; //floor(pT->GetStartDay());
  int    start_yr  = Options.julian_start_year;//pT->GetStartYear();
  double duration  = Options.duration;         //(interval*pTave->GetNumValues());
  double timestep  = Options.timestep;

  // below needed for correct mapping from time series to model time
  pTave->Initialize(start_day,start_yr,duration,timestep,Options);

  int nVals=(int)ceil(pTave->GetChunkSize()*interval); //Options.timestep);  // number of daily values
  int GridDims[3];
  GridDims[0] = pTave->GetCols(); GridDims[1] = pTave->GetRows(); GridDims[2] = nVals;

  // ----------------------------------------------------
  // Generate daily values (min, max, ave) from subdaily
  // ----------------------------------------------------
  if ( GetForcingGridIndexFromType(F_TEMP_DAILY_MIN) == DOESNT_EXIST ) {
    // for the first chunk the derived grid does not exist and has to be added to the model
    pTmin_daily = new CForcingGrid(* pTave);  // copy everything from tave; matrixes are deep copies
    pTmin_daily->SetForcingType(F_TEMP_DAILY_MIN);
    pTmin_daily->SetInterval(1.0);
    pTmin_daily->SetGridDims(GridDims);
    pTmin_daily->SetChunkSize(nVals);
    pTmin_daily->ReallocateArraysInForcingGrid();
  }
  else {
    // for all latter chunks the the grid already exists and values will be just overwritten
    pTmin_daily=GetForcingGrid(F_TEMP_DAILY_MIN);
  }

  if ( GetForcingGridIndexFromType(F_TEMP_DAILY_MAX) == DOESNT_EXIST ) {
    // for the first chunk the derived grid does not exist and has to be added to the model
    pTmax_daily = new CForcingGrid(* pTave);  // copy everything from tave; matrixes are deep copies
    pTmax_daily->SetForcingType(F_TEMP_DAILY_MAX);
    pTmax_daily->SetInterval(1.0);
    pTmax_daily->SetGridDims(GridDims);
    pTmax_daily->SetChunkSize(nVals);
    pTmax_daily->ReallocateArraysInForcingGrid();
  }
  else {
    // for all latter chunks the the grid already exists and values will be just overwritten
    pTmax_daily=GetForcingGrid(F_TEMP_DAILY_MAX);
  }

  if ( GetForcingGridIndexFromType(F_TEMP_DAILY_AVE) == DOESNT_EXIST ) {
    // for the first chunk the derived grid does not exist and has to be added to the model
    pTave_daily = new CForcingGrid(* pTave);  // copy everything from tave; matrixes are deep copies
    pTave_daily->SetForcingType(F_TEMP_DAILY_AVE);
    pTave_daily->SetInterval(1.0);
    pTave_daily->SetGridDims(GridDims);
    pTave_daily->SetChunkSize(nVals);
    pTave_daily->ReallocateArraysInForcingGrid();
  }
  else {
    // for all latter chunks the the grid already exists and values will be just overwritten
    pTave_daily=GetForcingGrid((F_TEMP_DAILY_AVE));
  }

  // (1) set weighting
  int nGridHRUs=pTave->GetnHydroUnits();
  int nCells=pTave->GetRows()*pTave->GetCols();
  double wt;
  for (int ik=0; ik<nGridHRUs; ik++) {                
    for (int ic=0; ic<nCells; ic++) {    
      wt = pTave->GetGridWeight(ik, ic);
      pTmin_daily->SetWeightVal(ik, ic, wt);
      pTmax_daily->SetWeightVal(ik, ic, wt);
      pTave_daily->SetWeightVal(ik, ic, wt);
    }
  }

  // (2) set indexes of on-zero weighted grid cells
  pTmin_daily->SetIdxNonZeroGridCells(nGridHRUs,nCells);
  pTmax_daily->SetIdxNonZeroGridCells(nGridHRUs,nCells);
  pTave_daily->SetIdxNonZeroGridCells(nGridHRUs,nCells);

  // (3) set forcing values
  for (int it=0; it<nVals; it++) {                    // loop over time points in buffer
    for (int ic=0; ic<pTave->GetNumberNonZeroGridCells(); ic++){         // loop over non-zero grid cell indexes
      pTmin_daily->SetValue(ic, it, pTave->GetValue_min(ic, (double)it*1.0/interval, int(1.0/interval)));
      pTmax_daily->SetValue(ic, it, pTave->GetValue_max(ic, (double)it*1.0/interval, int(1.0/interval)));
      pTave_daily->SetValue(ic, it, pTave->GetValue    (ic, (double)it*1.0/interval, int(1.0/interval)));
    }
  }

  if ( GetForcingGridIndexFromType(F_TEMP_DAILY_MIN) == DOESNT_EXIST ) {
    this->AddForcingGrid(pTmin_daily);
    if (Options.noisy){ printf("\n------------------------> TEMP_DAILY_MIN Added \n"); }
  }
  else {
    if (Options.noisy){ printf("\n------------------------> TEMP_DAILY_MIN Replace \n"); }
  }

  if ( GetForcingGridIndexFromType(F_TEMP_DAILY_MAX) == DOESNT_EXIST ) {
    this->AddForcingGrid(pTmax_daily);
    if (Options.noisy){ printf("\n------------------------> TEMP_DAILY_MAX Added \n"); }
  }
  else {
    if (Options.noisy){ printf("\n------------------------> TEMP_DAILY_MAX Replace \n"); }
  }

  if ( GetForcingGridIndexFromType(F_TEMP_DAILY_AVE) == DOESNT_EXIST ) {
    this->AddForcingGrid(pTave_daily);
    if (Options.noisy){ printf("\n------------------------> TEMP_DAILY_AVE Added \n"); }
  }
  else {
    if (Options.noisy){ printf("\n------------------------> TEMP_DAILY_AVE Replace \n"); }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Generates Tmin, Tmax and subhourly time series from daily average temperature time series
/// \note presumes existence of valid F_TEMP_DAILY_AVE
/// \note necessarily naive - it is hard to downscale with little temp data
//
void CModel::GenerateMinMaxSubdailyTempFromAve(const optStruct &Options)
{

  CForcingGrid *pTmin_daily,*pTmax_daily,*pTave_daily;
  double interval,wt;

  pTave_daily=GetForcingGrid((F_TEMP_DAILY_AVE));
  interval = pTave_daily->GetInterval();

  double start_day = Options.julian_start_day; //floor(pT->GetStartDay());
  int    start_yr  = Options.julian_start_year;//pT->GetStartYear();
  double duration  = Options.duration;         //(interval*pTave->GetNumValues());
  double timestep  = Options.timestep;

  // below needed for correct mapping from time series to model time
  pTave_daily->Initialize(start_day,start_yr,duration,timestep,Options);

  int nVals=(int)ceil(pTave_daily->GetChunkSize());                            // number of subdaily values (input resolution)
  int GridDims[3];
  GridDims[0] = pTave_daily->GetCols(); GridDims[1] = pTave_daily->GetRows(); GridDims[2] = nVals;

  // ----------------------------------------------------
  // Generate daily values (min, max) from daily average
  // ----------------------------------------------------
  if ( GetForcingGridIndexFromType(F_TEMP_DAILY_MIN) == DOESNT_EXIST ) {
    // for the first chunk the derived grid does not exist and has to be added to the model
    pTmin_daily = new CForcingGrid(* pTave_daily);  // copy everything from tave_daily; matrixes are deep copies
    pTmin_daily->SetForcingType(F_TEMP_DAILY_MIN);
    pTmin_daily->SetInterval(interval);  // input tmp_ave resolution //Options.timestep);
    pTmin_daily->SetGridDims(GridDims);
    pTmin_daily->SetChunkSize(nVals);
    pTmin_daily->ReallocateArraysInForcingGrid();
  }
  else {
    // for all latter chunks the the grid already exists and values will be just overwritten
    pTmin_daily=GetForcingGrid((F_TEMP_DAILY_MIN));
  }

  if ( GetForcingGridIndexFromType(F_TEMP_DAILY_MAX) == DOESNT_EXIST ) {
    // for the first chunk the derived grid does not exist and has to be added to the model
    pTmax_daily = new CForcingGrid(* pTave_daily);  // copy everything from tave_daily; matrixes are deep copies
    pTmax_daily->SetForcingType(F_TEMP_DAILY_MAX);
    pTmax_daily->SetInterval(interval); // input tmp_ave resolution //Options.timestep);
    pTmax_daily->SetGridDims(GridDims);
    pTmax_daily->SetChunkSize(nVals);
    pTmax_daily->ReallocateArraysInForcingGrid();
  }
  else {
    // for all latter chunks the the grid already exists and values will be just overwritten
    pTmax_daily=GetForcingGrid((F_TEMP_DAILY_MAX));
  }

  // (1) set weighting
  int nGridHRUs=pTave_daily->GetnHydroUnits();
  int nCells=pTave_daily->GetRows()*pTave_daily->GetCols();
  for (int ik=0; ik<nGridHRUs; ik++) {                     
    for (int ic=0; ic<nCells; ic++) {   
      wt = pTave_daily->GetGridWeight(ik, ic);
      pTmin_daily->SetWeightVal(ik, ic, wt);
      pTmax_daily->SetWeightVal(ik, ic, wt);
    }
  }

  // (2) set indexes of on-zero weighted grid cells
  pTmin_daily->SetIdxNonZeroGridCells(pTmin_daily->GetnHydroUnits(),pTmin_daily->GetRows()*pTmin_daily->GetCols());
  pTmax_daily->SetIdxNonZeroGridCells(pTmax_daily->GetnHydroUnits(),pTmax_daily->GetRows()*pTmax_daily->GetCols());

  // (3) set forcing values
  // Tmax       is with input time resolution
  // Tmin       is with input time resolution
  // Tave       is with model time resolution
  // Tave_daily is with daily resolution
  int chunksize=pTave_daily->GetChunkSize();
  int nNonZero=pTave_daily->GetNumberNonZeroGridCells();
  for (int it=0; it<nVals; it++) {          // loop over time points in buffer
    for (int ic=0; ic<nNonZero; ic++){      // loop over non-zero grid cell indexes
      pTmin_daily->SetValue(ic, it, pTave_daily->GetValue(ic, (double)(min(chunksize,it)))-4.0); // should be it+0.5
      pTmax_daily->SetValue(ic, it, pTave_daily->GetValue(ic, (double)(min(chunksize,it)))+4.0); // should be it+0.5
    }
  }

  if ( GetForcingGridIndexFromType(F_TEMP_DAILY_MIN) == DOESNT_EXIST ) {
    this->AddForcingGrid(pTmin_daily);
    if (Options.noisy){ printf("\n------------------------> TEMP_DAILY_MIN Added \n"); }
  }
  else {
    if (Options.noisy){ printf("\n------------------------> TEMP_DAILY_MIN Replace \n"); }
  }

  if ( GetForcingGridIndexFromType(F_TEMP_DAILY_MAX) == DOESNT_EXIST ) {
    this->AddForcingGrid(pTmax_daily);
    if (Options.noisy){ printf("\n------------------------> TEMP_DAILY_MAX Added \n"); }
  }
  else {
    if (Options.noisy){ printf("\n------------------------> TEMP_DAILY_MAX Replace \n"); }
  }

  // ----------------------------------------------------
  // Generate subdaily averages from daily values (min, max)
  // ----------------------------------------------------
  GenerateAveSubdailyTempFromMinMax(Options);
}

////////////////////////////////////////////////////////////////// 
/// \brief Generates precipitation as sum of snowfall and rainfall
/// \note  presumes existence of valid F_SNOWFALL and F_RAINFALL time series
//
void CModel::GeneratePrecipFromSnowRain(const optStruct &Options)
{

  CForcingGrid *pPre,*pSnow,*pRain;
  pSnow=GetForcingGrid((F_SNOWFALL));
  pRain=GetForcingGrid((F_RAINFALL));

  double start_day=Options.julian_start_day;
  int    start_yr =Options.julian_start_year;
  double duration =Options.duration;
  double timestep =Options.timestep;

  //below needed for correct mapping from time series to model time
  pSnow->Initialize(start_day,start_yr,duration,timestep,Options);
  pRain->Initialize(start_day,start_yr,duration,timestep,Options);

  double interval_snow = pSnow->GetInterval();
  double interval_rain = pRain->GetInterval();

  ExitGracefullyIf(interval_snow != interval_rain,
                   "CModel::GeneratePrecipFromSnowRain: rainfall and snowfall must have the same time resolution!",BAD_DATA);

  int    nVals     = pSnow->GetChunkSize();
  int    GridDims[3];
  GridDims[0] = pSnow->GetCols(); GridDims[1] = pSnow->GetRows(); GridDims[2] = nVals;

  // ----------------------------------------------------
  // Generate precipitation
  // ----------------------------------------------------
  if ( GetForcingGridIndexFromType(F_PRECIP) == DOESNT_EXIST ) {

    // for the first chunk the derived grid does not exist and has to be added to the model
    pPre = new CForcingGrid(* pSnow);  // copy everything from snowfall; matrixes are deep copies
    pPre->SetForcingType(F_PRECIP);
    pPre->SetInterval(pSnow->GetInterval());        // will be at same time resolution as precipitation
    pPre->SetGridDims(GridDims);
    pPre->SetChunkSize(nVals);                     // has same number of timepoints as precipitation
    pPre->ReallocateArraysInForcingGrid();
  }
  else {

    // for all latter chunks the the grid already exists and values will be just overwritten
    pPre=GetForcingGrid((F_PRECIP));
  }

  // (1) set weighting
  int nGridHRUs=pPre->GetnHydroUnits();
  int nCells=pPre->GetRows()*pPre->GetCols();
  for (int ik=0; ik<nGridHRUs; ik++) {                          // loop over HRUs
    for (int ic=0; ic<nCells; ic++) {               // loop over cells = rows*cols
      pPre->SetWeightVal(ik, ic, pSnow->GetGridWeight(ik, ic));
    }
  }

  // (2) set indexes of on-zero weighted grid cells
  pPre->SetIdxNonZeroGridCells(pPre->GetnHydroUnits(),pPre->GetRows()*pPre->GetCols());

  // (3) set forcing values
  int chunk_size=pPre->GetChunkSize();
  int nNonZero  =pPre->GetNumberNonZeroGridCells();
  for (int it=0; it<chunk_size; it++) {     // loop over time points in buffer
    for (int ic=0; ic<nNonZero; ic++){      // loop over non-zero grid cell indexes
      pPre->SetValue(ic, it, pSnow->GetValue(ic, it) + pRain->GetValue(ic, it));  // precipitation = sum of snowfall and rainfall
    }
  }

  if ( GetForcingGridIndexFromType(F_PRECIP) == DOESNT_EXIST ) {
    this->AddForcingGrid(pPre);
    if (Options.noisy){ printf("\n------------------------> PRECIP Added \n"); }
  }
  else {
    if (Options.noisy){ printf("\n------------------------> PRECIP Replace \n"); }
  }

}

//////////////////////////////////////////////////////////////////
/// \brief Generates rainfall as copy of precipitation
/// \note  presumes existence of valid F_PRECIP time series
//
void CModel::GenerateRainFromPrecip(const optStruct &Options)
{
  CForcingGrid *pPre,*pRain;
  pPre=GetForcingGrid((F_PRECIP));

  double start_day=Options.julian_start_day;
  int    start_yr =Options.julian_start_year;
  double duration =Options.duration;
  double timestep =Options.timestep;

  //below needed for correct mapping from time series to model time
  pPre->Initialize(start_day,start_yr,duration,timestep,Options);

  int    nVals     = pPre->GetChunkSize();
  int    GridDims[3];
  GridDims[0] = pPre->GetCols(); GridDims[1] = pPre->GetRows(); GridDims[2] = nVals;

  // ----------------------------------------------------
  // Generate rainfall
  // ----------------------------------------------------
  if ( GetForcingGridIndexFromType(F_RAINFALL) == DOESNT_EXIST ) {

    // for the first chunk the derived grid does not exist and has to be added to the model
    pRain = new CForcingGrid(* pPre);  // copy everything from precip; matrixes are deep copies
    pRain->SetForcingType(F_RAINFALL);
    pRain->SetInterval(pPre->GetInterval());        // will be at same time resolution as precipitation
    pRain->SetGridDims(GridDims);
    pRain->SetChunkSize(nVals);                     // has same number of timepoints as precipitation
    pRain->ReallocateArraysInForcingGrid(); //This should NOT be done (JRC)

    // copy weighting (is this even necessary? should already be there from Copy constructor \todo [optimize] JRC
    int nGridHRUs=pRain->GetnHydroUnits();
    int nCells=pRain->GetRows()*pRain->GetCols();
    for(int ik=0; ik<nGridHRUs; ik++) {              // loop over HRUs
      for(int ic=0; ic<nCells; ic++) {               // loop over cells = rows*cols
        pRain->SetWeightVal(ik,ic,pPre->GetGridWeight(ik,ic));
      }
    }

    // determine indexes of on-zero weighted grid cells //should already be there from Copy constructor \todo [optimize] JRC
    pRain->SetIdxNonZeroGridCells(pRain->GetnHydroUnits(),pRain->GetRows()*pRain->GetCols());

  }
  else {

    // for all latter chunks the the grid already exists and values will be just overwritten
    pRain=GetForcingGrid((F_RAINFALL));
  }


  // set forcing values
  int chunk_size=pRain->GetChunkSize();
  int nNonZero  =pRain->GetNumberNonZeroGridCells();
  for (int it=0; it<chunk_size; it++) {                   // loop over time points in buffer
    for (int ic=0; ic<nNonZero; ic++){                    // loop over non-zero grid cell indexes
      pRain->SetValue(ic, it , pPre->GetValue(ic, it));   // copies precipitation values
    }
  }

  if ( GetForcingGridIndexFromType(F_RAINFALL) == DOESNT_EXIST ) {
    this->AddForcingGrid(pRain);
    if (Options.noisy){ printf("\n------------------------> RAINFALL Added \n"); }
  }
  else {
    if (Options.noisy){ printf("\n------------------------> RAINFALL Replace \n"); }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Generates snowfall time series being constantly zero
/// \note  presumes existence of either rainfall or snowfall
//
void CModel::GenerateZeroSnow(const optStruct &Options)
{

  CForcingGrid *pPre(NULL),*pSnow;
  if (ForcingGridIsAvailable(F_PRECIP))   { pPre=GetForcingGrid((F_PRECIP)); }
  if (ForcingGridIsAvailable(F_RAINFALL)) { pPre=GetForcingGrid((F_RAINFALL)); }

  double start_day=Options.julian_start_day;
  int    start_yr =Options.julian_start_year;
  double duration =Options.duration;
  double timestep =Options.timestep;

  //below needed for correct mapping from time series to model time
  pPre->Initialize(start_day,start_yr,duration,timestep,Options);

  int    nVals     = pPre->GetChunkSize();
  int    GridDims[3];
  GridDims[0] = pPre->GetCols(); GridDims[1] = pPre->GetRows(); GridDims[2] = nVals;

  // ----------------------------------------------------
  // Generate snowfall
  // ----------------------------------------------------
  if ( GetForcingGridIndexFromType(F_SNOWFALL) == DOESNT_EXIST ) {

    // for the first chunk the derived grid does not exist and has to be added to the model
    pSnow = new CForcingGrid(* pPre);  // copy everything from precip; matrixes are deep copies
    pSnow->SetForcingType(F_SNOWFALL);
    pSnow->SetInterval(pPre->GetInterval());        // will be at same time resolution as precipitation
    pSnow->SetGridDims(GridDims);
    pSnow->SetChunkSize(nVals);                     // has same number of timepoints as precipitation
    pSnow->ReallocateArraysInForcingGrid();
  }
  else {

    // for all latter chunks the the grid already exists and values will be just overwritten
    pSnow=GetForcingGrid((F_SNOWFALL));
  }

  // (1) set weighting
  int nGridHRUs=pSnow->GetnHydroUnits();
  int nCells=pSnow->GetRows()*pSnow->GetCols();
  for (int ik=0; ik<nGridHRUs; ik++) {                           
    for (int ic=0; ic<nCells; ic++) {           
      pSnow->SetWeightVal(ik, ic, pPre->GetGridWeight(ik, ic));
    }
  }

  // (2) set indexes of on-zero weighted grid cells
  pSnow->SetIdxNonZeroGridCells(pSnow->GetnHydroUnits(),pSnow->GetRows()*pSnow->GetCols());

  // (3) set forcing values
  int chunk_size=pSnow->GetChunkSize();
  int nNonZero  =pSnow->GetNumberNonZeroGridCells();
  for (int it=0; it<chunk_size; it++) {                   // loop over time points in buffer
    for (int ic=0; ic<nNonZero; ic++){                    // loop over non-zero grid cell indexes
      pSnow->SetValue(ic, it , 0.0);                      // fills everything with 0.0
    }
  }

  if ( GetForcingGridIndexFromType(F_SNOWFALL) == DOESNT_EXIST ) {
    this->AddForcingGrid(pSnow);
    if (Options.noisy){ printf("\n------------------------> SNOWFALL Added \n"); }
  }
  else {
    if (Options.noisy){ printf("\n------------------------> SNOWFALL Replace \n"); }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns average fraction of snow in precipitation between time t and following n timesteps
/// \param x_col  [in] Column index
/// \param y_row  [in] Row index
/// \param t      [in] Time index
/// \param n      [in] Number of time steps
/// \return average fraction of snow in precipitation between time t and following n timesteps
//
double CModel::GetAverageSnowFrac(const int idx, const double t, const int n) const
{

  CForcingGrid *pSnow,*pRain;
  pSnow=GetForcingGrid((F_SNOWFALL));
  pRain=GetForcingGrid((F_RAINFALL));

  double snow = pSnow->GetValue(idx, t, n);
  double rain = pRain->GetValue(idx, t, n);

  if ((snow+rain)==0.0){return 0.0;}
  return snow/(snow+rain);

}






