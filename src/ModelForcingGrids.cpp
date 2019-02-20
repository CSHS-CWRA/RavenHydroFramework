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
/// \brief if forcing grid of type typ doesnt exist, creates full copy of pGrid 
///    but with potentially new time interval specified
///    otherwise just returns existing grid of type typ
//
CForcingGrid *CModel::ForcingCopyCreate(const CForcingGrid *pGrid,
                                        const forcing_type typ,
                                        const double &interval, 
                                        const int nVals)
{
  static CForcingGrid *pTout;
  if ( GetForcingGridIndexFromType(typ) == DOESNT_EXIST )
  { // for the first chunk the derived grid does not exist and has to be added to the model
    int    GridDims[3];
    GridDims[0] = pGrid->GetCols(); 
    GridDims[1] = pGrid->GetRows(); 
    GridDims[2] = nVals;

    pTout = new CForcingGrid(* pGrid);  // copy everything from pGrid; matrixes are deep copies

    pTout->SetForcingType(typ);
    pTout->SetInterval(interval);            
    pTout->SetGridDims(GridDims);
    pTout->SetChunkSize(nVals);        
    pTout->ReallocateArraysInForcingGrid();
  }
  else
  {
    // for all latter chunks the the grid already exists and values will be just overwritten
    pTout=GetForcingGrid(typ);
  }

  //JRC - is below  necessary? - shouldn't this be done in deep copy? [it seems to be necessary]
  int nGridHRUs=pGrid->GetnHydroUnits();
  int nCells=pGrid->GetRows()*pGrid->GetCols(); //handles 3D or 2D
  double wt;
  for (int k=0; k<nGridHRUs; k++) {                
    for (int ic=0; ic<nCells; ic++) {    
      wt = pGrid->GetGridWeight(k, ic);
      pTout->SetWeightVal(k, ic, wt); 
    }
  }

  pTout->SetIdxNonZeroGridCells(nGridHRUs,nCells);

  return pTout;
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

  int    nVals     = (int)ceil(pTmin->GetChunkSize() * pTmin->GetInterval());

  pTave_daily = ForcingCopyCreate(pTmin,F_TEMP_DAILY_AVE,1.0,nVals);

  double t=0.0;
  int chunk_size =pTave_daily->GetChunkSize();
  int nNonZero   =pTave_daily->GetNumberNonZeroGridCells();
  int nValsPerDay= (int)(1.0 / pTmin->GetInterval());
  for (int it=0; it<chunk_size; it++) {           // loop over time points in buffer
    for (int ic=0; ic<nNonZero;ic++){             // loop over non-zero grid cell indexes
      pTave_daily->SetValue(ic, it , 0.5*(pTmin->GetValue_avg(ic, t * nValsPerDay, nValsPerDay) + 
                                          pTmax->GetValue_avg(ic, t * nValsPerDay, nValsPerDay)));
    }
    t+=1.0;
  }

  AddForcingGrid(pTave_daily,F_TEMP_DAILY_AVE);

  // ----------------------------------------------------
  // Generate subdaily temperature values
  // ----------------------------------------------------
  if (Options.timestep<(1.0-TIME_CORRECTION)) //tstep < 1 day
  {
    int    nVals     = (int)ceil(pTave_daily->GetChunkSize()/Options.timestep);

    pTave = ForcingCopyCreate(pTmin,F_TEMP_AVE,Options.timestep,nVals);

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
    int    nVals     = pTave_daily->GetChunkSize();

    pTave = ForcingCopyCreate(pTave_daily,F_TEMP_AVE,1.0,nVals);

    int nNonZero  =pTave->GetNumberNonZeroGridCells();
    for (int it=0; it<nVals; it++) {                                 // loop over time points in buffer
      for (int ic=0; ic<nNonZero; ic++){                                  // loop over non-zero grid cell indexes
        pTave->SetValue(ic, it , pTave_daily->GetValue(ic, it));  // --> just copy daily average values
      }
    }

    AddForcingGrid(pTave,F_TEMP_AVE);
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

  pTmin_daily = ForcingCopyCreate(pTave,F_TEMP_DAILY_MIN,1.0,nVals);
  pTmax_daily = ForcingCopyCreate(pTave,F_TEMP_DAILY_MAX,1.0,nVals);
  pTave_daily = ForcingCopyCreate(pTave,F_TEMP_DAILY_AVE,1.0,nVals);

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
  int       nVals = (int)ceil(pTave_daily->GetChunkSize()); // number of subdaily values (input resolution)
  
  pTmin_daily = ForcingCopyCreate(pTave_daily,F_TEMP_DAILY_MIN,interval,nVals);
  pTmax_daily = ForcingCopyCreate(pTave_daily,F_TEMP_DAILY_MAX,interval,nVals);

  int chunksize=pTave_daily->GetChunkSize();
  int nNonZero=pTave_daily->GetNumberNonZeroGridCells();
  for (int it=0; it<nVals; it++) {          // loop over time points in buffer
    for (int ic=0; ic<nNonZero; ic++){      // loop over non-zero grid cell indexes
      pTmin_daily->SetValue(ic, it, pTave_daily->GetValue(ic, min(chunksize,it))-4.0); // should be it+0.5
      pTmax_daily->SetValue(ic, it, pTave_daily->GetValue(ic, min(chunksize,it))+4.0); // should be it+0.5
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

  pPre = ForcingCopyCreate(pSnow,F_PRECIP,interval_snow,nVals);

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

  pRain = ForcingCopyCreate(pPre,F_RAINFALL,pPre->GetInterval(),nVals);

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
  if (ForcingGridIsAvailable(F_PRECIP))   { pPre=GetForcingGrid((F_PRECIP)); }
  if (ForcingGridIsAvailable(F_RAINFALL)) { pPre=GetForcingGrid((F_RAINFALL)); }

  pPre->Initialize(Options);//needed for correct mapping from time series to model time

  int    nVals     = pPre->GetChunkSize();

  pSnow = ForcingCopyCreate(pPre,F_SNOWFALL,pPre->GetInterval(),nVals);

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






