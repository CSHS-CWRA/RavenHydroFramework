/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2026 the Raven Development Team
  ----------------------------------------------------------------*/
#include "TimeSeriesABC.h"

/*****************************************************************
   Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/
///////////////////////////////////////////////////////////////////
/// \brief Implementation of the constructor of an abstract time serie
///
/// \param type      [in] type of time series being constructed (subclass)
/// \param Name      [in] name of time series
/// \param tag       [in] data tag
/// \param filename  [in] original source file
//
CTimeSeriesABC::CTimeSeriesABC(ts_type type,
                               string  Name,
                               long long loc_ID,
                               string  filename="")
{
  _type       =type;
  _name       =Name;
  _loc_ID     =loc_ID;
  _srcfile    =filename;
  _constit_ind=DOESNT_EXIST;
  _demand_ID  =DOESNT_EXIST;

  _aSampVal =NULL; //generated in Resample() routine
  _nSampVal =0;    //generated in Resample() routine
  _sampInterval=1.0;
}

///////////////////////////////////////////////////////////////////
/// \brief Implementation of copy constructor (for which the address of a time series is passed)
/// \param &t [in] Address of a time series of which a "copy" is made
//
CTimeSeriesABC::CTimeSeriesABC(string Name,
                               const CTimeSeriesABC &t)
{
  _type=t.GetType();
  _name=Name;
  _loc_ID=t.GetLocID();
  _srcfile = "";
  _constit_ind=t.GetConstitInd();
  _demand_ID=t.GetDemandID();
  _nSampVal=t.GetNumSampledValues();
  _sampInterval=t.GetSampledInterval();
  _aSampVal=NULL;
  _aSampVal=new double [_nSampVal];
  ExitGracefullyIf(_aSampVal==NULL,"CTimeSeriesABC copy constructor",OUT_OF_MEMORY);
  for(int nn=0;nn<_nSampVal;nn++){
    _aSampVal[nn]=t.GetSampledValue(nn);
  }
}
///////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
//
CTimeSeriesABC::~CTimeSeriesABC()
{
  delete [] _aSampVal; _aSampVal=NULL;
}

/*****************************************************************
   Accessors
------------------------------------------------------------------
*****************************************************************/

///////////////////////////////////////////////////////////////////
/// \brief Returns type (subclass) of time series
/// \return type of time series
//
CTimeSeriesABC::ts_type CTimeSeriesABC::GetType() const{return _type;}

///////////////////////////////////////////////////////////////////
/// \brief Returns name of time series
/// \return name of time series
//
string CTimeSeriesABC::GetName()  const{return _name;}

//////////////////////////////////////////////////////////////////
/// \brief Returns location identifier
/// \return location identifier (e.g., HRU ID or Basin ID of observation data)
//
long long CTimeSeriesABC::GetLocID()      const{return _loc_ID;}

//////////////////////////////////////////////////////////////////
/// \brief returns demand ID
/// \return demand ID
//
long long   CTimeSeriesABC::GetDemandID() const {return _demand_ID; }

//////////////////////////////////////////////////////////////////
/// \brief Returns location identifier
/// \return constituent index (or DOESNT_EXIST, otherwise)
//
int    CTimeSeriesABC::GetConstitInd() const { return _constit_ind; }

//////////////////////////////////////////////////////////////////
/// \brief Sets location identifier
/// \param ID - HRU ID or SubBasin ID linked to observations
//
void   CTimeSeriesABC::SetLocID(long long ID) { _loc_ID=ID; }

//////////////////////////////////////////////////////////////////
/// \brief Sets demand ID
/// \param ID - long long integer tag (e.g., demand ID)
//
void   CTimeSeriesABC::SetDemandID(long long ID) { _demand_ID=ID; }

//////////////////////////////////////////////////////////////////
/// \brief Sets constituent index
/// \param c - constituent index
//
void   CTimeSeriesABC::SetConstitInd(const int c) {_constit_ind=c; }

//////////////////////////////////////////////////////////////////
/// \brief Returns source input file
/// \return source file of time series (if one exists)
//
string   CTimeSeriesABC::GetSourceFile()      const{return _srcfile;}

///////////////////////////////////////////////////////////////////
/// \brief Returns average value of time series during timestep nn of model simulation (nn=0..nSteps-1)
/// \notes must be called after resampling
///
/// \param nn [in] time step number (measured from simulation start)
/// \return Average value of time series data over time step nn
//
double CTimeSeriesABC::GetSampledValue(const int nn) const
{
  if (nn>_nSampVal-1){
    return RAV_BLANK_DATA;
  }
  return _aSampVal[nn];
}
///////////////////////////////////////////////////////////////////
/// \brief Returns the model time of the resampled time series at index nn
/// \param nn [in] Index
/// \return time of time series data point for which nn is an index
//
double CTimeSeriesABC::GetSampledTime(const int nn) const
{
  return _sampInterval*nn;
}
///////////////////////////////////////////////////////////////////
/// \brief Returns data interval of resampled timeseries
/// \notes this is usually model timestep
/// \return resampled data interval (in days)
//
double CTimeSeriesABC::GetSampledInterval() const
{
  return _sampInterval;
}
///////////////////////////////////////////////////////////////////
/// \brief Returns the number resampled values
/// \return Number of resampled values
//
int CTimeSeriesABC::GetNumSampledValues() const
{
  return _nSampVal;
}
///////////////////////////////////////////////////////////////////
/// \brief enables sampled values to be overwritten to store model-generated information
/// \notes must use with caution!! Poor encapsulation of data
///
/// \param nn [in] index of sampled value
/// \param val [in] value to overwrite
//
void CTimeSeriesABC::SetSampledValue(const int nn, const double &val)
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf(nn>=_nSampVal, "CTimeSeries::SetSampledValue: Overwriting array allocation",RUNTIME_ERR);
#endif
  _aSampVal[nn]=val;
}
//////////////////////////////////////////////////////////////////
/// \brief Initializes arrays for the resampled time series
///
/// \param nSampVal [in] Number of values in the resampled time series
/// \param sampInterval [in] Time interval of the resampled time series (in days)
//
void CTimeSeriesABC::InitializeResample(const int nSampVal, const double &sampInterval)
{
  _nSampVal=nSampVal;
  _sampInterval=sampInterval;
  ExitGracefullyIf(_nSampVal<=0,"CTimeSeriesABC::InitializeResample: bad # of samples",RUNTIME_ERR);

  if(_aSampVal!=NULL){delete[] _aSampVal;}

  _aSampVal=NULL;
  _aSampVal=new double [_nSampVal];
  ExitGracefullyIf(_aSampVal==NULL,"CTimeSeriesABC::InitializeResample",OUT_OF_MEMORY);

  //Initialize with blanks
  for (int nn=0;nn<_nSampVal;nn++){
    _aSampVal[nn] = RAV_BLANK_DATA;
  }
}
