/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------*/
#include "TimeSeries.h"
#include "ParseLib.h"
#include "Forcings.h"
/*****************************************************************
   Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the time series constructor when time series is single constant value
/// \param one_value [in] Constant value of time series
//
CTimeSeries::CTimeSeries(string Name, string tag, double one_value)
  :CTimeSeriesABC(ts_regular,Name,tag)
{
  _start_day=0.0;
  _start_year=1900;
  _nPulses  =2;
  _interval =ALMOST_INF;
  _aVal     =new double [_nPulses];
  _aVal [0] =one_value;
  _aVal [1] =one_value;
  _sub_daily=false;
  _t_corr   =0.0;
  _pulse    =true;
  _aSampVal =NULL; //generated in Resample() routine
  _nSampVal =0;    //generated in Resample() routine
}

///////////////////////////////////////////////////////////////////
/// \brief Implementation of the time series constructor for uniformly spaced data points
///
/// \param strt_day      [in] Start day of time series (Julian decimal date)
/// \param start_yr      [in] Start year of time series
/// \param *aT_start     [in] Array of start times of pulses (time since start day) [size NumPulses]
/// \param *aValues      [in] Array of magnitudes of pulses (variable units) [size NumPulses]
/// \param NumPulses     [in] Number of pulses which occur
/// \param is_pulse_type [in] Flag determining time-series type (pulse-based vs. piecewise-linear)
//
CTimeSeries::CTimeSeries(string     Name,
                         string     tag,
                         string     filename,
                         double     strt_day,
                         int        start_yr,
                         double     data_interval,
                         double    *aValues,
                         const int  NumPulses,
                         const bool is_pulse_type)
  :CTimeSeriesABC(ts_regular,Name,tag,filename)
{
  _start_day =strt_day;
  _start_year=start_yr;
  _pulse     =is_pulse_type;
  _interval  =data_interval;
  _nPulses   =NumPulses;

  ExitGracefullyIf(NumPulses<=0,
                   "CTimeSeries: Constructor: no entries in time series",BAD_DATA);
  ExitGracefullyIf(_interval<=0,
                   "CTimeSeries: Constructor: negative time interval is not allowed",BAD_DATA);

  _aVal=NULL;
  _aVal=new double [_nPulses];
  ExitGracefullyIf(_aVal==NULL,"CTimeSeries: Constructor",OUT_OF_MEMORY);

  for (int n=0; n<_nPulses;n++)
  {
    _aVal [n]=aValues [n];
  }
  _sub_daily=(_interval<(1.0-TIME_CORRECTION));//to account for potential roundoff error
  _t_corr=0.0;

  _aSampVal =NULL; //generated in Resample() routine
  _nSampVal =0;
}

///////////////////////////////////////////////////////////////////
/// \brief Implementation of the time series constructor for empty time series used to store model generated data
///
CTimeSeries::CTimeSeries(string     Name,
                         string     tag,
                         string     filename,
                         double     strt_day,
                         int        start_yr,
                         double     data_interval,
                         const int  NumPulses,
                         const bool is_pulse_type)
  :CTimeSeriesABC(ts_regular,Name,tag,filename)
{
  _start_day =strt_day;
  _start_year=start_yr;
  _pulse     =is_pulse_type;
  _interval  =data_interval;
  _nPulses   =NumPulses;

  ExitGracefullyIf(NumPulses<=0,
                   "CTimeSeries: Constructor: no entries in time series",BAD_DATA);
  ExitGracefullyIf(_interval<=0,
                   "CTimeSeries: Constructor: negative time interval is not allowed",BAD_DATA);

  _aVal=NULL;
  _aVal=new double [_nPulses];
  ExitGracefullyIf(_aVal==NULL,"CTimeSeries: Constructor",OUT_OF_MEMORY);

  for (int n=0; n<_nPulses;n++)
  {
    _aVal [n]=RAV_BLANK_DATA;
  }

  _sub_daily=(_interval<(1.0-TIME_CORRECTION));//to account for potential roundoff error
  _t_corr=0.0;

  _aSampVal =NULL; //generated in Resample() routine
  _nSampVal =0;
}

///////////////////////////////////////////////////////////////////
/// \brief Implementation of copy constructor (for which the address of a time series is passed)
/// \param &t [in] Address of a time series of which a "copy" is made
//
CTimeSeries::CTimeSeries(string Name,
                         const CTimeSeries &t)
  :CTimeSeriesABC(Name,t)
{
  _start_day =t.GetStartDay();
  _start_year=t.GetStartYear();
  _pulse     =t.IsPulseType();
  _interval  =t.GetInterval();
  _nPulses   =t.GetNumValues();
  _aVal      =NULL;
  _aVal      =new double [_nPulses];
  ExitGracefullyIf(_aVal==NULL,"CTimeSeries copy constructor",OUT_OF_MEMORY);
  for (int n=0; n<_nPulses;n++)
  {
    _aVal[n]=t.GetValue(n);
  }
  _sub_daily=t._sub_daily;
  _t_corr   =0.0;

  _aSampVal =NULL; //generated in Resample() routine
  _nSampVal =0;
}
///////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
//
CTimeSeries::~CTimeSeries()
{
  if (DESTRUCTOR_DEBUG){cout<<"    DELETING TIME SERIES"<<endl;}
  delete [] _aVal;     _aVal =NULL;
  delete [] _aSampVal; _aSampVal=NULL;
}

/*****************************************************************
   Accessors
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Returns true if smallest time interval is daily
/// \return Boolean indicating if time interval is daily
//
bool   CTimeSeries::IsDaily()      const{return !_sub_daily;}

///////////////////////////////////////////////////////////////////
/// \brief Returns number of pulses
/// \return Number of pulses
//
int    CTimeSeries::GetNumValues() const{return _nPulses;}

///////////////////////////////////////////////////////////////////
/// \brief Returns data interval
/// \return data interval (in days)
//
double  CTimeSeries::GetInterval() const{return _interval;}

///////////////////////////////////////////////////////////////////
/// \brief Returns start year of time series data
/// \return Start year of time series data
//
int    CTimeSeries::GetStartYear() const{return _start_year;}

///////////////////////////////////////////////////////////////////
/// \brief Returns start day of time series data
/// \return Start day of time series data (Julian decimal date)
//
double CTimeSeries::GetStartDay () const{return _start_day;}

///////////////////////////////////////////////////////////////////
/// \brief Returns the time for index n in local time.
/// \remark For pulse-based time series, time is time at the start of the pulse.
/// \param n [in] Index
/// \return Magnitude of time series data point for which n is an index
//
double CTimeSeries::GetTime(const int n) const{return _interval*n;}

///////////////////////////////////////////////////////////////////
/// \brief Returns magnitude of time series data point for which n is an index
/// \param n [in] Pulse index
/// \return Magnitude of time series data point for which n is an index
//
double CTimeSeries::GetValue(const int n)const{return _aVal[n];}

///////////////////////////////////////////////////////////////////
/// \brief Returns true if time series is pulse-based
/// \return True if time series is pulse-based (i.e., data points represent a constant value over the time interval)
//
bool CTimeSeries::IsPulseType()  const{return _pulse;}

///////////////////////////////////////////////////////////////////
/// \brief Enables queries of time series values using model time
/// \details Calculates _t_corr, correction to global model time, checks for overlap with model duration, resamples to model time step/day
/// \remark t=0 corresponds to first day with recorded values at that gauge
///
/// \param model_start_day [in] Julian start day of model
/// \param model_start_year [in] start year of model
/// \param model_duration [in] Duration of model, in days
/// \param timestep [in] mdoel timestep, in days
/// \param is_observation [in] - true if this is an observation time series, rather than model input
void CTimeSeries::Initialize( const double model_start_day,   //julian day
                              const    int model_start_year,  //year
                              const double model_duration,     //days
                              const double timestep,           //days
                              const bool   is_observation)
{
  //_t_corr is number of days between model start date and gauge
  //start date (positive if data exists before model start date)

  _t_corr = -TimeDifference(model_start_day,model_start_year,_start_day,_start_year);

  //QA/QC: Check for overlap issues between time series duration and model duration
  //------------------------------------------------------------------------------
  double duration = (double)(_nPulses)*_interval;
  double local_simulation_start = (_t_corr);
  double local_simulation_end = (model_duration + _t_corr);

  if (!is_observation) //time series overlap is only necessary for forcing data
  {
    if (duration < local_simulation_start)  //out of data before simulation starts!
    {
      cout << " Time series " << GetName() << endl;
      cout << "  time series start day, year, duration :" << _start_day << "," << _start_year << " " << duration << endl;
      cout << "  model start day, year, duration :" << model_start_day << "," << model_start_year << " " << model_duration << endl;
      ExitGracefully(
        "CTimeSeries::Initialize: time series forcing data not available for entire model simulation duration", BAD_DATA);
    }
    if (duration + timestep < local_simulation_end)    //run out of data before simulation finishes
    {                                                  //+timesteps is for coincdent duration & data
      cout << " Time series " << GetName() << endl;
      cout << "  time series start day, year, duration :" << _start_day << "," << _start_year << " " << duration << endl;
      cout << "  model start day, year, duration :" << model_start_day << "," << model_start_year << " " << model_duration << endl;
      ExitGracefully(
        "CTimeSeries::Initialize: time series forcing data not available at end of model simulation", BAD_DATA);
    }
    if ((local_simulation_start<0) || (_start_year>model_start_year))     //data does not begin until after simulation
    {
      cout << " Time series " << GetName() << endl;
      cout << "  time series start day, year, duration :" << _start_day << "," << _start_year << " " << duration << endl;
      cout << "  model start day, year, duration :" << model_start_day << "," << model_start_year << " " << model_duration << endl;
      ExitGracefully(
        "CTimeSeries::Initialize: time series forcing data not available at beginning of model simulation", BAD_DATA);
    }
  }

  // Resample time series
  //------------------------------------------------------------------------------
  if (is_observation){Resample(timestep, model_duration+timestep);} //extra timestep needed for last observation of continuous hydrograph
  else               {Resample(timestep, model_duration);}
}

//////////////////////////////////////////////////////////////////
/// \brief Resamples time series to regular intervals of 1 timestep (tstep) over model duration
///
/// \param &tstep [in] Model timestep (in days)
/// \param &model_duration [in] Model simulation duration (in days)
//
void CTimeSeries::Resample(const double &tstep,          //days
                           const double &model_duration)//days
{
  int nSampVal=(int)(ceil(model_duration/tstep-TIME_CORRECTION));

  if (!_pulse){nSampVal++;}
  InitializeResample(nSampVal,tstep);

  double t=0;
  for (int nn=0;nn<_nSampVal;nn++){
    if (_pulse){
      _aSampVal[nn] = GetAvgValue(t, tstep);
    }
    else{
      _aSampVal[nn] = GetValue(t);
    }
    t+=tstep;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Initializes arrays for the resampled time series
///
/// \param nSampVal [in] Number of values in the resampled time series
/// \param sampInterval [in] Time interval of the resampled time series (in days)
//
void CTimeSeries::InitializeResample(const int nSampVal, const double sampInterval)
{
  _nSampVal=nSampVal;
  _sampInterval=sampInterval;
  ExitGracefullyIf(_nSampVal<=0,"CTimeSeries::Resample: bad # of samples",RUNTIME_ERR);

  if(_aSampVal!=NULL){delete[] _aSampVal;}

  _aSampVal=NULL;
  _aSampVal=new double [_nSampVal];
  ExitGracefullyIf(_aSampVal==NULL,"CTimeSeries::Resample",OUT_OF_MEMORY);

  //Initialize with blanks
  for (int nn=0;nn<_nSampVal;nn++){
    _aSampVal[nn] = RAV_BLANK_DATA;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns index of time period for time t in terms of local time
///
/// \param &t_loc [in] Local time
/// \return Index of time period for time t; if t_loc<0, returns 0, if t_loc>_nPulses-1, returns nPulses-1
//
int CTimeSeries::GetTimeIndex(const double &t_loc) const
{
  return min((int)(max(floor(t_loc/_interval),0.0)),_nPulses-1);
}
///////////////////////////////////////////////////////////////////
/// \brief Returns point value of time series at time t
/// \param &t [in] Global time at which point value of time series is to be determined
/// \return Point value of time series data at time t
//
double CTimeSeries::GetValue(const double &t) const
{
  //takes in GLOBAL time
  static int n;
  n = 0;
  double t_loc = t + _t_corr;
  n = GetTimeIndex(t_loc);

  //if ((n < 0) ||(n>_nPulses - 1)){ return RAV_BLANK_DATA;}

  if (_pulse){return _aVal[n];}
  else      {
    if (n==_nPulses-1){return _aVal[n];}
    return _aVal[n]+(t_loc-(double)(n)*_interval)/(_interval)*(_aVal[n+1]-_aVal[n]);
  }
}

///////////////////////////////////////////////////////////////////
/// \brief Returns average value of time series over interval t to t+tstep
/// \param &t [in] Left bound of model time interval over which average value of time series is to be determined
/// \param &tstep [in] Duration of interval over which average value of time series is to be determined
/// \return Average time series value across interval t to t+tstep
/// if interval is (mostly) blank data, returns BLANK_DATA, correction for blank data slices
//
double CTimeSeries::GetAvgValue(const double &t, const double &tstep) const
{
  int n1(0), n2(0);
  double sum = 0;
  double t_loc = t + _t_corr; //time re-referenced to start of time series
  n1 = GetTimeIndex(t_loc);
  n2 = GetTimeIndex(t_loc + tstep);
  //t_loc       is now between n1*_interval and (n1+1)*_interval
  //t_loc+tstep is now between n2*_interval and (n2+1)*_interval
  double inc;
  double blank = 0;
  if (t_loc < -TIME_CORRECTION) {return RAV_BLANK_DATA; }
  if (t_loc >= _nPulses*_interval) {return RAV_BLANK_DATA; }
  if (_pulse){
    if (n1 == n2){ return _aVal[n1]; }
    else{
      sum = 0;
      inc = ((double)(n1 + 1)*_interval - t_loc);
      if (_aVal[n1] == RAV_BLANK_DATA)  { blank += inc; }
      else                                          { sum += _aVal[n1] * inc; }

      for (int n = n1 + 1; n < n2; n++){
        if (_aVal[n] == RAV_BLANK_DATA) { blank += _interval; }
        else                                        { sum += _aVal[n] * _interval; }
      }
      inc = ((t_loc + tstep) - (double)(n2)*_interval);
      if (_aVal[n2] == RAV_BLANK_DATA)  { blank += inc; }
      else                                          { sum += _aVal[n2] * inc; }
    }
  }
  else{
    ExitGracefully("CTimeSeries::GetAvgValue (non-pulse)", STUB);
  }
  if (blank / tstep > 0.001){return RAV_BLANK_DATA;}

  return sum/tstep;
}

///////////////////////////////////////////////////////////////////
/// \brief Returns minimum value of time series between time t to t+tstep
/// \param &t [in] Left endpoint of time interval over which minimum value of time series is to be determined
/// \param &tstep [in] Duration of time interval over which minimum time series value is to be determined
/// \return Minimum value of time series over interval t to t+tstep
//
double CTimeSeries::GetMinValue(const double &t, const double &tstep) const
{
  int n1(0),n2(0);
  double vmin(ALMOST_INF);
  double t_loc=t+_t_corr;
  n1=GetTimeIndex(t_loc      +TIME_CORRECTION);//to account for potential roundoff error
  n2=GetTimeIndex(t_loc+tstep-TIME_CORRECTION);

  ExitGracefullyIf(!_pulse,"CTimeSeries::GetMinValue (non-pulse)",STUB);

  if (n1==n2){return _aVal[n1];}
  for (int n=n1;n<=n2;n++){lowerswap(vmin,_aVal[n]);}
  return vmin;
}

///////////////////////////////////////////////////////////////////
/// \brief Returns maximum value of time series between time t to t+tstep
/// \param &t [in] Left endpoint of time interval over which maximum value of time series is to be determined
/// \param &tstep [in] Duration of time interval over which maximum time series value is to be determined
/// \return Maximum value of time series over interval t to t+tstep
//
double CTimeSeries::GetMaxValue(const double &t, const double &tstep) const
{
  int n1(0),n2(0);
  double vmax(-ALMOST_INF);
  double t_loc=t+_t_corr;
  n1=GetTimeIndex(t_loc      +TIME_CORRECTION);//to account for potential roundoff error
  n2=GetTimeIndex(t_loc+tstep-TIME_CORRECTION);

  ExitGracefullyIf(!_pulse,"CTimeSeries::GetMaxValue (non-pulse)",STUB);

  if (n1==n2){return _aVal[n1];}
  for (int n=n1;n<=n2;n++){upperswap(vmax,_aVal[n]);}
  return vmax;
}

///////////////////////////////////////////////////////////////////
/// \brief Multiply all time series values by factor
/// \param &factor [in] Multiplier for time series data
//
void   CTimeSeries::Multiply     (const double &factor)
{
  for (int n=0; n<_nPulses;n++)
  {
    _aVal [n]*=factor;
  }
}
///////////////////////////////////////////////////////////////////
/// \brief Returns average value of time series during timestep nn of model simulation
/// \notes must be called after resampling
///
/// \param nn [in] time step number (measured from simulation start)
/// \return Average value of time series data over time step nn
//
double CTimeSeries::GetSampledValue(const int nn) const
{
  if (nn>_nSampVal-1){
    //cout <<"CTimeSeries::GetSampledValue "<< nn << " "<<_nSampVal<<endl;
    return RAV_BLANK_DATA;
  }
  return _aSampVal[nn];
}
///////////////////////////////////////////////////////////////////
/// \brief Returns the model time of the resampled time series at index nn
/// \param nn [in] Index
/// \return time of time series data point for which nn is an index
//
double CTimeSeries::GetSampledTime(const int nn) const
{
  return _sampInterval*nn;
}
///////////////////////////////////////////////////////////////////
/// \brief Returns data interval of resampled timeseries
/// \notes this is usually model timestep
/// \return resampled data interval (in days)
//
double CTimeSeries::GetSampledInterval() const
{
  return _sampInterval;
}
///////////////////////////////////////////////////////////////////
/// \brief Returns the number resampled values
/// \return Number of resampled values
//
int CTimeSeries::GetNumSampledValues() const
{
  return _nSampVal;
}
///////////////////////////////////////////////////////////////////
/// \brief Returns average value of time series during day model_day of model simulation
/// \notes uses precalculated value if available
///
/// \param model_day [in] (measured from simulation start)
/// \return Average value of time series data over specified model day
//
double CTimeSeries::GetDailyAvg(const int model_day) const
{
  return GetAvgValue((double)(model_day),1.0);
}
///////////////////////////////////////////////////////////////////
/// \brief Returns minimum value of time series during day model_day of model simulation
/// \notes uses precalculated value if available
///
/// \param model_day [in] (measured from simulation start)
/// \return Minimum value of time series data over specified model day
//
double CTimeSeries::GetDailyMin(const int model_day) const
{
  return GetMinValue((double)(model_day),1.0);
}
///////////////////////////////////////////////////////////////////
/// \brief Returns maximum value of time series during day model_day of model simulation
/// \notes uses precalculated value if available
///
/// \param model_day [in] (measured from simulation start)
/// \return Maximum value of time series data over specified model day
//
double CTimeSeries::GetDailyMax(const int model_day) const
{
  return GetMaxValue((double)(model_day),1.0);
}

///////////////////////////////////////////////////////////////////
/// \brief Returns the value of the modelled time series at time t
/// \notes only really valid for special ts storing model results for diagnostics.
/// \param &t [in] Global time at which point value of time series is to be determined
/// \param type [in] type of observed time series
/// \return value of modelled time series
//
double CTimeSeries::GetModelledValue(const double &t,const ts_type type) const
{
  static int n;
  n=0;
  double t_loc=t+_t_corr;
  n=GetTimeIndex(t_loc);

  if (type == ts_regular){return GetAvgValue(t,_sampInterval);}
  else      {
    if (n==_nPulses-1){return _aVal[n];}
    return _aVal[n]+(t_loc-(double)(n)*_interval)/(_interval)*(_aVal[n+1]-_aVal[n]);
  }
}

///////////////////////////////////////////////////////////////////
/// \brief enables time series values to be overwritten to store model-generated information
/// \notes must use with caution!! Poor encapsulation of data
///
/// \param n [in] index of  value
/// \param val [in] value to overwrite
//
void CTimeSeries::SetValue(const int n, const double &val)
{
  ExitGracefullyIf(n>=_nPulses, "CTimeSeries::SetValue: Overwriting array allocation",BAD_DATA);
  _aVal[n]=val;
}

///////////////////////////////////////////////////////////////////
/// \brief enables sampled values to be overwritten to store model-generated information
/// \notes must use with caution!! Poor encapsulation of data
///
/// \param nn [in] index of sampled value
/// \param val [in] value to overwrite
//
void CTimeSeries::SetSampledValue(const int nn, const double &val)
{
  ExitGracefullyIf(nn>=_nSampVal, "CTimeSeries::SetSampledValue: Overwriting array allocation",BAD_DATA);
  _aSampVal[nn]=val;
}

///////////////////////////////////////////////////////////////////
/// \brief sums together two time series
/// \notes time series must have same start time, duration, and sample interval
///
///
/// \param pTS1 [in] pointer to first time series
/// \return pTS2 [in] pointer to second time series
//
CTimeSeries  *CTimeSeries::Sum(CTimeSeries *pTS1, CTimeSeries *pTS2, string name) //STATIC
{

  double start_day=pTS1->GetStartDay();
  int    start_yr =pTS1->GetStartYear();
  double interval =pTS1->GetInterval();
  int    nPulses  =pTS1->GetNumValues();
  bool   is_pulse =pTS1->IsPulseType();

  ExitGracefullyIf((start_day!=pTS2->GetStartDay()) || (start_yr!=pTS1->GetStartYear()),
                   "CTimeSeries::Sum: time series must have same start date",BAD_DATA);
  ExitGracefullyIf(interval!=pTS2->GetInterval(),
                   "CTimeSeries::Sum: time series must have same interval",BAD_DATA);
  ExitGracefullyIf(nPulses!=pTS2->GetNumValues(),
                   "CTimeSeries::Sum: time series must have same number of data points",BAD_DATA);

  double *_aVal=new double [nPulses];
  for (int n=0;n<nPulses;n++){
    _aVal[n]=pTS1->GetValue(n)+pTS2->GetValue(n);
  }
  CTimeSeries *pOut=new CTimeSeries(name,"","",start_day,start_yr,interval,_aVal,nPulses,is_pulse);
  delete [] _aVal;
  return pOut;
}

///////////////////////////////////////////////////////////////////
/// \brief Parses standard single time series format from file and creates time series object
/// \param *p [in] CParser object pointing to input file
/// \param is_pulse [in] Flag determining time-series type (piecewise-uniform [pulsed] vs. piecewise-linear)
/// \return Pointer to created time series
// Assumed format: 
//  [yyyy-mm-dd] [hh:mm:ss.0] [timestep] [nMeasurements]
//  data1 data2
//  data3 data4 data5
//  ...
//  dataN-1 dataN
// where N=nMeasurements
// OR 
// :AnnualCycle J F M A M J J A S O N D 
//
CTimeSeries *CTimeSeries::Parse(CParser *p, bool is_pulse, string name, string tag, const optStruct &Options, bool shift_to_per_ending)
{

  char   *s[MAXINPUTITEMS];
  int    Len;
  double start_day;
  int    start_yr;
  double tstep;
  int    nMeasurements;

  // for NetCDF handling
  string FileNameNC;            // filename of NetCDF
  string VarNameNC;             // name of variable in NetCDF
  string DimNamesNC_stations;   // name of station dimension (optional)
  string DimNamesNC_time;       // name of time dimension
  int    StationIdx;            // station to read (numbering starts with 1;
  //                            // mandatory if DimNamesNC_stations given)
  double TimeShift;             // time shift of data (fractional day by which
  //                            // read data should be shifted)
  double LinTrans_a;            // Linear transformation: new = a*data + b
  double LinTrans_b;            // Linear transformation: new = a*data + b

  CTimeSeries *pTimeSeries=NULL;

  p->Tokenize(s,Len);
  if (IsComment(s[0],Len)){p->Tokenize(s,Len);} //try again

  // READFROMNETCDF FORMAT =======================================================
  if(!strcmp(s[0],":ReadFromNetCDF"))
  {

    // initialize
    FileNameNC = "None";
    VarNameNC  = "None";
    StationIdx = 0;   // numbering starts with 1
    TimeShift  = 0.0; // no time shift applied
    LinTrans_a = 1.0; // no linear transformation: 1.0*x+0.0 --> a = 1.0
    LinTrans_b = 0.0; // no linear transformation: 1.0*x+0.0 --> b = 0.0 
    
    while (strcmp(s[0],":EndReadFromNetCDF")) {
      p->Tokenize(s,Len);
      if(!strcmp(s[0],":FileNameNC")) { FileNameNC = s[1]; FileNameNC = CorrectForRelativePath(FileNameNC ,Options.rvt_filename); }
      if(!strcmp(s[0],":VarNameNC"))  { VarNameNC  = s[1]; }
      if(!strcmp(s[0],":DimNamesNC")) {
        if (Len == 2) {
          DimNamesNC_stations = "None";
          DimNamesNC_time     = s[1];
        }
        else {
          if (Len == 3) {
            DimNamesNC_stations = s[1];
            DimNamesNC_time     = s[2];
          }
          else {
            ExitGracefully("CTimeSeries::Parse: ReadFromNetCDF requires 1 or 2 arguments for :DimNamesNC",BAD_DATA); return NULL;
          }
        }
      }
      if(!strcmp(s[0],":StationIdx"))      { StationIdx = s_to_i(s[1]); }
      if(!strcmp(s[0],":TimeShift"))       { TimeShift  = s_to_d(s[1]); }
      if(!strcmp(s[0],":LinearTransform")) {
        if (Len != 3) {
          ExitGracefully("CTimeSeries::Parse: ReadFromNetCDF 2 arguments for :LinearTransform (i.e. a,b in a*x+b)",BAD_DATA); return NULL;
        }
        else {
          LinTrans_a = s_to_d(s[1]);
          LinTrans_b = s_to_d(s[2]);
        }
      }
    }

    // Check that all information given and consistent
    if ( !strcmp(FileNameNC.c_str(),"None") ) {
      ExitGracefully("CTimeSeries::Parse: ReadFromNetCDF: FileNameNC needs to be given",BAD_DATA); return NULL;
    }
    if ( !strcmp(VarNameNC.c_str(),"None") ) {
      ExitGracefully("CTimeSeries::Parse: ReadFromNetCDF: VarNameNC needs to be given",BAD_DATA); return NULL;
    }
    if ( !strcmp(DimNamesNC_time.c_str(),"None") ) {
      ExitGracefully("CTimeSeries::Parse: ReadFromNetCDF: DimNamesNC needs to be given",BAD_DATA); return NULL;
    }
    if ( strcmp(DimNamesNC_stations.c_str(),"None") && (StationIdx == 0) ) {
      ExitGracefully("CTimeSeries::Parse: ReadFromNetCDF: If there are multiple stations in NetCDF file, the StationIdx argument must be given (numbering starts with 1) to identify station that should be read in",BAD_DATA); return NULL;
    }

     /* cout <<"Filename   = "<<FileNameNC<<endl;
     cout <<"VarNameNC  = "<<VarNameNC<<endl;
     cout <<"DimNamesNC = "<<DimNamesNC_stations<<endl;
     cout <<"DimNamesNC = "<<DimNamesNC_time<<endl;
     cout <<"StationIdx = "<<StationIdx<<endl;
     cout <<"TimeShift  = "<<TimeShift<<endl;  */

    // if DimNamesNC_stations == "None" then NetCDF variable is assumed to be 1D (only time)
    // if DimNamesNC_stations != "None" then StationIdx must be >= 1
    pTimeSeries = ReadTimeSeriesFromNetCDF(
                                           Options,              // model options (such as simulation period)
                                           name,                 // ForcingType
					   tag,                  // critical information about timeseries, e.g. subbasin ID or HRU ID
                                           FileNameNC,           // file name of NetCDF
                                           VarNameNC,            // name of variable in NetCDF
                                           DimNamesNC_stations,  // name of station dimension (optional; default=None)
                                           DimNamesNC_time,      // name of time dimension (mandatory)
                                           StationIdx,           // idx of station to be read
                                           //                    // (only used if DimNamesNC_stations not None)
                                           TimeShift,            // time shift of data (fractional day by which
                                           //                    // read data should be shifted)
                                           LinTrans_a,           // linear transformation: a in new = a*data+b
                                           LinTrans_b            // linear transformation: b in new = a*data+b
                                           );

    return pTimeSeries;
  }
  
  if (Len<4){p->ImproperFormat(s); cout <<"Length:" <<Len<<endl;}
  
  // ANNUALCYCLE FORMAT ==========================================================
  if(!strcmp(s[0],":AnnualCycle"))
  {
    if(Len<13){ ExitGracefully("CTimeSeries::Parse: incorrect format of AnnualCycle command",BAD_DATA); return NULL; }
    double *aVal;
    int nVals=int(Options.duration)+1;
    aVal =new double [nVals];
    time_struct tt;   
    double aMonVal[12];
    for(int i=0;i<12;i++){ aMonVal[i]=s_to_d(s[i+1]); }
    for(int n=0;n<nVals;n++){
      JulianConvert(double(n),Options.julian_start_day,Options.julian_start_year,tt);
      aVal[n]=InterpolateMo(aMonVal,tt,Options);
    }
    pTimeSeries=new CTimeSeries(name,tag,p->GetFilename(),Options.julian_start_day,Options.julian_start_year,1.0,aVal,nVals,is_pulse); 
    delete[] aVal; aVal =NULL;
    p->Tokenize(s,Len);//read closing term (e.g., ":EndData")
    if(string(s[0]).substr(0,4)!=":End"){
      ExitGracefully("CTimeSeries: Parse: no :EndData command used with :AnnualCycle command ",BAD_DATA);
    }
    return pTimeSeries;
  }

  // REGULAR FORMAT ==============================================================
  if ((string(s[0]).length()==10) &&
      ((string(s[0]).substr(4,1)=="/") ||
       (string(s[0]).substr(4,1)=="-")))
  { //in timestamp format [yyyy-mm-dd] [hh:mm:ss.0] [timestep] [nMeasurements]
    time_struct tt;

    tt=DateStringToTimeStruct(string(s[0]),string(s[1]));
    start_day=tt.julian_day;
    start_yr =tt.year;

    string tString=s[2];
    if ((tString.length()>=2) && ((tString.substr(2,1)==":") || (tString.substr(1,1)==":"))){//support for hh:mm:ss.00 format in timestep
      time_struct tt;
      tt=DateStringToTimeStruct("0000-01-01",tString);
      tstep=FixTimestep(tt.julian_day);
    }
    else{
      tstep =FixTimestep(s_to_d(s[2]));
    }

    nMeasurements=s_to_i(s[3]);
    //units =s[4]
  }
  else
  { //julian date format [nMeasurements] [start_day] [start_yr] [timestep]
    nMeasurements=s_to_i(s[0]);
    start_day    =s_to_d(s[1]);
    start_yr     =s_to_i(s[2]);
    tstep        =FixTimestep(s_to_d(s[3]));
    //units =s[4]
  }
  if (shift_to_per_ending)
  {
    start_day+=Options.timestep;
    int leap=0;
    if (IsLeapYear(start_yr)){ leap = 1; }
    if (start_day>=365+leap){start_day-=365+leap; start_yr++;}
  }

  double *aVal;
  aVal =new double [nMeasurements];
  if (aVal == NULL){
    ExitGracefully("CTimeSeries::Parse",OUT_OF_MEMORY);
  }

  int n=0;
  //cout << n << " "<<nMeasurements << " " << s[0] << " "<<Len<<" "<<strcmp(s[0],"&")<<" "<<p->Tokenize(s,Len)<<endl;
  while ((n<nMeasurements) && (!p->Tokenize(s,Len)))
  {
    if (IsComment(s[0],Len)){p->Tokenize(s,Len);}//try again
    for(int i=0;i<Len;i++){
      if (n>=nMeasurements)
      {
        cout << " n | nMeasurements " << n << " "<<nMeasurements<<endl;
        ExitGracefully("CTimeSeries::Parse: Bad number of time series points",BAD_DATA);
      }
      if (!is_numeric(s[i]))
      {
        ExitGracefully( ("Non-numeric value found in time series (line " +to_string(p->GetLineNumber())+" of file "+p->GetFilename()+")").c_str(),BAD_DATA_WARN);
      }
      aVal [n]=fast_s_to_d(s[i]);
      n++;
    }
  }

  if (n!=nMeasurements){
    cout<<p->GetFilename()<<endl;
    cout << " n | nMeasurements " << n << " "<<nMeasurements<<endl;
    ExitGracefully("CTimeSeries::Parse: Bad number of time series points",BAD_DATA);}

  p->Tokenize(s,Len);//read closing term (e.g., ":EndData")
  if(string(s[0]).substr(0,4)!=":End"){
    ExitGracefully("CTimeSeries: Parse: exceeded specified number of time series points in sequence or no :EndData command used. ",BAD_DATA);
  }

  pTimeSeries=new CTimeSeries(name,tag,p->GetFilename(),start_day,start_yr,tstep,aVal,n,is_pulse);
  delete [] aVal;  aVal =NULL;
  return pTimeSeries;
}

///////////////////////////////////////////////////////////////////
/// \brief Parses multicolumn single time series format and create array of time series objects
/// \note Creates output array and aType array; these must be
///  deleted outside of this routine
/// \param *p [in] CParser object pointing to input file
/// \param &nTS [out] number of time series parsed
/// \param *aType [out] array (size nTS)  of forcing types, one per time series parsed
/// \param is_pulse [in] Flag determining time-series type (pulse-based vs. piecewise-linear)
/// \return array (size nTS) of pointers to time series
//
CTimeSeries **CTimeSeries::ParseMultiple(CParser *p, int &nTS, forcing_type *aType, bool is_pulse)
{
  char *s[MAXINPUTITEMS];
  int Len;
  int i;
  double start_day;
  int    start_yr;
  double tstep;
  int    nMeasurements;
  bool   noisy=false;
  CTimeSeries **pTimeSeries=NULL;

  //timestamp & numdata info ----------------------------------------------
  p->Tokenize(s,Len);
  if (IsComment(s[0],Len)){p->Tokenize(s,Len);}//try again
  if (Len<4){p->ImproperFormat(s);}

  if ((string(s[0]).length()==10) &&
      ((string(s[0]).substr(4,1)=="/") ||
       (string(s[0]).substr(4,1)=="-")))
  {//in timestamp format  [yyyy-mm-dd] [hh:mm:ss.0] [timestep] [nMeasurements]
    time_struct tt;
    tt=DateStringToTimeStruct(string(s[0]),string(s[1]));
    start_day=tt.julian_day;
    start_yr =tt.year;

    string tString=s[2];
    if ((tString.length()>=2) && ((tString.substr(2,1)==":") || (tString.substr(1,1)==":"))){//support for hh:mm:ss.00 format in timestep
      time_struct tt;
      tt=DateStringToTimeStruct("0000-01-01",tString);
      tstep=FixTimestep(tt.julian_day);
    }
    else{
      tstep =FixTimestep(s_to_d(s[2]));
    }

    nMeasurements=s_to_i(s[3]);
  }
  else{//julian date format [nMeasurements] [start_day] [start_yr] [timestep]
    nMeasurements=s_to_i(s[0]);
    start_day    =s_to_d(s[1]);
    start_yr     =s_to_i(s[2]);
    tstep        =FixTimestep(s_to_d(s[3]));
  }

  //:Parameters line ----------------------------------------------
  p->Tokenize(s,Len);
  if (IsComment(s[0],Len)){p->Tokenize(s,Len);}//try again
  if (strcmp(s[0],":Parameters")){
    ExitGracefully("CTimeSeries::ParseMultiple : MultiData command improperly formatted",BAD_DATA);
  }
  nTS=Len-1;

  ExitGracefullyIf(nTS>MAX_MULTIDATA,
                   "CTimeSeries::ParseMultiple: exceeded max entries in :Multidata time series input",BAD_DATA);
  //for (i=0; i<Len;i++){cout<<s[i]<<" ";}cout<<endl;

  for (i=1;i<Len;i++){
    aType[i-1]=GetForcingTypeFromString(s[i]);
    if (noisy){cout<<"  "<<s[i]<<endl;}
    ExitGracefullyIf(i-1>=MAX_MULTIDATA,
                     "CTimeSeries::ParseMultiple : exceeded max number of time series in MultiData command",BAD_DATA);

    string warn="CTimeSeries::ParseMultiple : unrecognized time series type "+to_string(s[i])+ " in MultiData command";
    ExitGracefullyIf(aType[i-1]==F_UNRECOGNIZED,  warn.c_str(),BAD_DATA);
  }

  //:Units line ----------------------------------------------
  p->Tokenize(s,Len);
  if (IsComment(s[0],Len)){p->Tokenize(s,Len);}//try again
  if (strcmp(s[0],":Units")){
    ExitGracefully("CTimeSeries::ParseMultiple : MultiData command improperly formatted",BAD_DATA);
  }

  //Tabular Data ----------------------------------------------
  double **aVal;
  aVal=new double *[nTS];
  for (i=0;i<nTS;i++){
    aVal[i] =new double [nMeasurements];
  }
  int n=0;
  while (!p->Tokenize(s,Len))
  {
    if (!IsComment(s[0],Len))
    {
      if (!strcmp(s[0],":EndMultiData")){break;}
      else{
        if (Len!=nTS)
        {
          cout<<"line number: "<<p->GetLineNumber()<<endl;
          cout<<"measurement " <<n+1 << " of "<<nMeasurements<<endl;
          string error="CTimeSeries:ParseMultiple: wrong number of columns in :MultiData data. File: "+p->GetFilename();
          ExitGracefully(error.c_str(),BAD_DATA);
        }
        if (n>=nMeasurements)
        {
          cout<<" first val | nMeaurements: "<<s[0]<<" | "<<n<<" nMeasurements"<<endl;
          string error="CTimeSeries::ParseMultiple: Bad number of time series points. File: "+p->GetFilename();
          ExitGracefully(error.c_str(),BAD_DATA);
        }

        for(i=0;i<Len;i++){
          if (!is_numeric(s[i]))
          {
            ExitGracefully( ("Non-numeric value found in time series (line " +to_string(p->GetLineNumber())+" of file "+p->GetFilename()+")").c_str(),BAD_DATA_WARN);
          }
          aVal [i][n]=fast_s_to_d(s[i]);
          //aVal [i][n]=s_to_d(s[i]);
        }
        n++;
      }
    }
  }

  if (n!=nMeasurements){
    string error="CTimeSeries::ParseMultiple: Insufficient number of time series points. File: "+p->GetFilename();
    ExitGracefully(error.c_str(),BAD_DATA);
  }

  // finished. Now process data --------------------------------
  pTimeSeries=new CTimeSeries *[nTS];
  for (i=0;i<nTS;i++){
    pTimeSeries[i]=NULL;
    pTimeSeries[i]=new CTimeSeries(ForcingToString(aType[i]),"",p->GetFilename(),start_day,start_yr,tstep,aVal[i],nMeasurements,is_pulse);
  }

  //delete dynamic memory---------------------------------------
  for (i=0;i<nTS;i++){
    delete [] aVal[i];   aVal[i] =NULL;
  }
  delete [] aVal; aVal=NULL;
  return pTimeSeries;

}
//////////////////////////////////////////////////////////////////
/// \brief Parses set of time series from Ensim .tb0 file format
/// \note Column names must correspond to Raven forcing tags
///
/// \param filename [in] Ensim file name, with .tb0 extension
/// \param &nTS [out] number of time series parsed
/// \param *aType [out] array (size nTS)  of forcing types, one per time series parsed
/// \return array (size nTS) of pointers to time series
//
CTimeSeries **CTimeSeries::ParseEnsimTb0(string filename, int &nTS, forcing_type *aType)
{
  char *s[MAXINPUTITEMS];
  int    Len;
  int    i,line=0;
  double start_day=0;
  int    start_yr=1900;
  double tstep=1.0;
  int    nMeasurements;
  CTimeSeries **pTimeSeries=NULL;
  bool noisy=false;
  bool period_ending=false;
  time_struct tt;

  ifstream INPUT;
  INPUT.open(filename.c_str());
  if (INPUT.fail())
  {
    string errString="ERROR opening file: "+filename;
    ExitGracefully(errString.c_str(),BAD_DATA);
    return NULL;
  }

  ifstream inFile(filename.c_str());
  int linecount=(int)std::count(istreambuf_iterator<char>(inFile),
                                istreambuf_iterator<char>(), '\n');//count # of lines in file

  CParser *p=new CParser(INPUT,filename,line);

  while (!p->Tokenize(s,Len))
  {
    if (noisy){ cout << "-->reading line " << p->GetLineNumber() << ": ";}

    if      (Len==0){}//Do nothing
    else if (s[0][0]=='#')                      {if (noisy){cout<<"#"<<endl;}}//Comment, do nothing
    else if (!strcmp(s[0],":FileType")         ){if (noisy){cout<<"Filetype"<<endl;}}//do nothing
    else if (!strcmp(s[0],":Application")      ){if (noisy){cout<<"Application"<<endl;}}//do nothing
    else if (!strcmp(s[0],":Version")          ){if (noisy){cout<<"Version"<<endl;}}// do nothing
    else if (!strcmp(s[0],":WrittenBy")        ){if (noisy){cout<<"WrittenBy"<<endl;}}// do nothing
    else if (!strcmp(s[0],":CreationDate")     ){if (noisy){cout<<"CreationDate"<<endl;}}// do nothing
    else if (!strcmp(s[0],":ColumnMetaData")   ){if (noisy){cout<<"ColumnMetaData"<<endl;}}// do nothing
    else if (!strcmp(s[0],":EndColumnMetaData")){if (noisy){cout<<"EndColumnMetaData"<<endl;}}// do nothing
    else if (!strcmp(s[0],":ColumnUnits")      ){if (noisy){cout<<"ColumnUnits"<<endl;}}// do nothing
    else if (!strcmp(s[0],":ColumnType")       ){if (noisy){cout<<"ColumnType"<<endl;}}// do nothing
    else if (!strcmp(s[0],":ColumnFormat")     ){if (noisy){cout<<"ColumnFormat"<<endl;}}// do nothing
    else if (!strcmp(s[0],":Format"))
    {
      if (string(s[1]) == "PeriodEnding"){period_ending=true;}
    }
    else if (!strcmp(s[0],":StartTime")      )
    {
      if (noisy){cout<<"StartTime"<<endl;}
      if ((string(s[1]).length()==10) &&
          ((string(s[1]).substr(4,1)=="/") || (string(s[1]).substr(4,1)=="-")))
        //if (IsValidDateString(s[1]))
      {//in timestamp format
        tt=DateStringToTimeStruct(string(s[1]),string(s[2]));
        start_day=tt.julian_day;
        start_yr =tt.year;
      }
      else
      {
        string errString = "ParseEnsimTb0: Bad date format: " + string(s[1]);
        ExitGracefully(errString.c_str(),BAD_DATA);
      }
    }
    else if (!strcmp(s[0],":DeltaT"))
    {
      if (noisy){cout<<"DeltaT"<<endl;}
      string tString=s[1];
      if ((tString.length()>=2) && ((tString.substr(2,1)==":") || (tString.substr(1,1)==":"))){//support for hh:mm:ss.00 format
        tstep=DateStringToTimeStruct("0000-01-01",tString).julian_day;
      }
      else{
        string errString = "ParseEnsimTb0: Bad DeltaT format: " + string(s[1]);
        ExitGracefully(errString.c_str(),BAD_DATA);
      }
    }
    else if (!strcmp(s[0],":ColumnName"))
    {
      if (noisy){cout<<"ColumnName"<<endl;}
      nTS=Len-1;
      for (i=1;i<Len;i++){
        aType[i-1]=GetForcingTypeFromString(s[i]);
        if (noisy){cout<<"  "<<s[i]<<endl;}
        ExitGracefullyIf(i-1>=MAX_MULTIDATA,
                         "CTimeSeries::ParseMultiple : exceeded max number of time series in MultiData command",BAD_DATA);
        ExitGracefullyIf(aType[i-1]==F_UNRECOGNIZED,
                         "CTimeSeries::ParseEnsimtb0 : unrecognized time series type in MultiData command",BAD_DATA);
      }
    }
    else if (!strcmp(s[0],":EndHeader")         )
    { //Start reading data
      //need to peek to end of file to see how many items remain!
      if (noisy){cout<<"EndHeader. Reading data..."<<endl;}
      nMeasurements=linecount-p->GetLineNumber()+1;
      if (noisy){cout <<"  Number of records in file: "<<nMeasurements<<endl;}
      int n=0;
      double **aVal;
      aVal=new double *[nTS];
      for (i=0;i<nTS;i++){
        aVal[i] =new double [nMeasurements];
      }
      while ((n<nMeasurements) && (!p->Tokenize(s,Len)))
      {
        if (Len!=nTS){
          cout<<"measurement " <<n << " of "<<nMeasurements<<endl;
          ExitGracefully("CTimeSeries:ParseMultiple: wrong number of columns in Ensim tb0 data",BAD_DATA);
        }
        for(i=0;i<Len;i++){
          aVal [i][n]=s_to_d(s[i]);
          //cout<<  aVal [i][n]<<" | ";
        }
        //cout<<endl;
        n++;
      }

      // Make sure that the nMeasurements is actually equal to the number of valid records found (for n<nMeasurements)
      nMeasurements = n;
      if (noisy){cout <<"  Number of valid measurements found: "<<nMeasurements<<endl;}

      // finished. Now process data --------------------------------
      if (period_ending){
        start_day=start_day-tstep;
        if ((start_day<0) && IsLeapYear(start_yr-1)){start_day+=366;start_yr-=1;}
        else if (start_day<0)                       {start_day+=365;start_yr-=1;}
      }
      pTimeSeries=new CTimeSeries *[nTS];
      for (i=0;i<nTS;i++){
        pTimeSeries[i]=NULL;
        pTimeSeries[i]=new CTimeSeries(ForcingToString(aType[i]),"",filename,start_day,start_yr,tstep,aVal[i],nMeasurements,true);
      }

      //delete dynamic memory---------------------------------------
      for (i=0;i<nTS;i++){
        delete [] aVal[i];   aVal[i] =NULL;
      }
      delete [] aVal; aVal=NULL;
      INPUT.close();

      return pTimeSeries;
    }
  }
  INPUT.close();
  return NULL;
}

//////////////////////////////////////////////////////////////////
/// \brief Reads a time series from a NetCDF file
/// \note  Will return just a vector that needs to be converted into TimeSeries object afterwards
///
/// \param  Options             [in] global model otions such as simulation period
/// \param  name                [in] forcing type
/// \param  tag                 [in] critical information about timeseries, e.g. subbasin ID or HRU ID
/// \param  FileNameNC          [in] file name of NetCDF
/// \param  VarNameNC           [in] name of variable in NetCDF
/// \param  DimNamesNC_stations [in] name of station dimension (optional; default=None)
/// \param  DimNamesNC_time     [in] name of time dimension (mandatory)
/// \param  StationIdx          [in] idx of station to be read (only used if DimNamesNC_stations not None)
/// \param  TimeShift           [in] time shift of data (fractional day by which read data should be shifted)
/// \param  LinTrans_a,         [in] linear transformation: a in new = a*data+b
/// \param  LinTrans_b          [in] linear transformation: b in new = a*data+b
/// \return array (size nTS) of pointers to time series
//
CTimeSeries *CTimeSeries::ReadTimeSeriesFromNetCDF(const optStruct &Options, string name,
						   string tag, string FileNameNC, string VarNameNC,
                                                   string DimNamesNC_stations, string DimNamesNC_time,
                                                   int StationIdx, double TimeShift, double LinTrans_a, double LinTrans_b)
{
  CTimeSeries *pTimeSeries=NULL; // time series of data
  
#ifdef _RVNETCDF_
  int    ncid;                  // file unit
  int    retval;                // error value for NetCDF routines
  int    dimid_x;               // id of x dimension
  int    dimid_t;               // id of time dimension 
  size_t dummy;                 // special type for GridDims required by nc routine
  int    nstations;             // length of x dimension (columns)
  int    ntime;                 // length of t dimension (time)
  int    varid_t;               // id of time variable
  int    varid_f;               // id of VarNameNC variable
  char   unit_t[200];           // special type for string of variable's unit required by nc routine
  string dash;                  // to check format of time unit string
  string colon;                 // to check format of time unit string
  string unit_t_str;            // to check format of time unit string
  
  // handling variable dimension orders
  int    dimids_var[2];         // ids of dimensions of a NetCDF variable
  int    dim1;                  // length of 1st dimension in NetCDF data
  int    dim2;                  // length of 2nd dimension in NetCDF data
  int    dim_order;             // order of dimensions of variable to read

  // final variables to create time series object
  int    nMeasurements;          // number of data points read
  double start_day=0.0;          // start day of time series
  int    start_yr=1900;          // start year of time series
  double tstep=1.0;              // time difference between two data points
  bool   is_pulse;               // if data are pulses

  // -------------------------------
  // (1) open NetCDF read-only (get ncid)
  // -------------------------------
  if (Options.noisy){ cout<<"Start reading time series for "<< VarNameNC << " from NetCDF file "<< FileNameNC << endl; }
  retval = nc_open(FileNameNC.c_str(), NC_NOWRITE, &ncid);
  HandleNetCDFErrors(retval);

  // -------------------------------
  // (2) get dimension lengths
  // -------------------------------
  //     (a) time dimension 
  retval = nc_inq_dimid (ncid, DimNamesNC_time.c_str(), &dimid_t);  HandleNetCDFErrors(retval);
  retval = nc_inq_dimlen(ncid, dimid_t, &dummy);                    HandleNetCDFErrors(retval);
  ntime  = static_cast<int>(dummy);  // convert returned 'size_t' to 'int'
  //     (b) other dimension (optional)
  if ( strcmp(DimNamesNC_stations.c_str(),"None") ) {
    retval = nc_inq_dimid (ncid, DimNamesNC_stations.c_str(), &dimid_x);  HandleNetCDFErrors(retval);
    retval = nc_inq_dimlen(ncid, dimid_x, &dummy);                        HandleNetCDFErrors(retval);
    nstations  = static_cast<int>(dummy);  // convert returned 'size_t' to 'int'
  }
  else {
    nstations = 0;
  }

  // -------------------------------
  // (3) time values and unit
  // -------------------------------
  //     (a) ID of time variable
  retval = nc_inq_varid(ncid,DimNamesNC_time.c_str(),&varid_t);
  HandleNetCDFErrors(retval);
  //     (b) time unit
  retval = nc_get_att_text(ncid,varid_t,"units",unit_t);   
  HandleNetCDFErrors(retval);
  //     (c) allocate array for time values
  int *time=NULL;
  time=new int [ntime];
  ExitGracefullyIf(time==NULL,"CTimeSeries::ReadTimeSeriesFromNetCDF",OUT_OF_MEMORY);
  //     (d) time values
  retval = nc_get_var_int(ncid,varid_t,&time[0]);
  HandleNetCDFErrors(retval);

  // -------------------------------
  // (4) check if time intervals are equal and set tstep
  // -------------------------------
  double delta_t = time[1] - time[0];
  for (int ii=1; ii<ntime-1 ; ii++)
  {
    double delta_t2 = time[ii+1] - time[ii];
    if ( abs(delta_t2 - delta_t) > 0.00001 )
    {
      printf("\n\n\nCTimeSeries::ReadTimeSeriesFromNetCDF: variable name: %s\n",VarNameNC.c_str());
      ExitGracefully("CTimeSeries::ReadTimeSeriesFromNetCDFt: time steps are not equal in NetCDF time variable",BAD_DATA);
    }
  }

  // -------------------------------
  // (5) determine tstep, start_day and start_yr depending on time unit
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
      ExitGracefully("CTimeSeries::ReadTimeSeriesFromNetCDF: time unit string has weird format!",BAD_DATA);
    }
    dash  = unit_t_str.substr(19, 1);  // second dash in date
    if ( !strstr(dash.c_str(), "-") ){
      unit_t_str.insert(17,"0");
    }
    colon = unit_t_str.substr(25, 1);  // first colon in time
    if ( !strstr(colon.c_str(), ":") ){
      unit_t_str.insert(23,"0");
    }
    colon = unit_t_str.substr(28, 1);  // second colon in time
    if ( !strstr(colon.c_str(), ":") ){
      unit_t_str.insert(26,"0");
    }
    string sDate = unit_t_str.substr(12,10);
    string sTime = unit_t_str.substr(23,8);
    time_struct tt = DateStringToTimeStruct(sDate, sTime);
    tstep   = (time[1] - time[0])/24.;
    ExitGracefullyIf(tstep<=0,
                     "CTimeSeries::ReadTimeSeriesFromNetCDF: Interval is negative!",BAD_DATA);
    AddTime(tt.julian_day,tt.year,time[0]*(1./24.),start_day,start_yr) ;

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
        ExitGracefully("CTimeSeries::ReadTimeSeriesFromNetCDF: time unit string has weird format!",BAD_DATA);
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
      string sDate = unit_t_str.substr(11,10);
      string sTime = unit_t_str.substr(22,8);
      time_struct tt = DateStringToTimeStruct(sDate, sTime);
      tstep   = (time[1] - time[0])/1.;
      ExitGracefullyIf(tstep<=0, "CTimeSeries::ReadTimeSeriesFromNetCDF: Interval is negative!",BAD_DATA);
      AddTime(tt.julian_day,tt.year,time[0]*(1.),start_day,start_yr) ;
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
          ExitGracefully("CTimeSeries::ReadTimeSeriesFromNetCDF: time unit string has weird format!",BAD_DATA);
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
        string sDate = unit_t_str.substr(14,10);
        string sTime = unit_t_str.substr(25,8);
        time_struct tt = DateStringToTimeStruct(sDate, sTime);
        tstep   = (time[1] - time[0])/24./60.;
        ExitGracefullyIf(tstep<=0,
                         "CTimeSeries::ReadTimeSeriesFromNetCDF: Interval is negative!",BAD_DATA);
        AddTime(tt.julian_day,tt.year,time[0]*(1./24./60.),start_day,start_yr) ;

      }
      else
        {
          if (strstr(unit_t, "seconds")) {  // if unit of time in seconds: seconds since 1989-01-01 00:00:00

            // ---------------------------
            // check if format is YYYY-MM-DD HH:MM:SS
            // fill with leading zeros if necessary
            // ---------------------------
            unit_t_str = to_string(unit_t);
            dash = unit_t_str.substr(18, 1);  // first dash in date
            if ( !strstr(dash.c_str(), "-") ){
              printf("time unit string: %s\n",unit_t_str.c_str());
              ExitGracefully("CTimeSeries::ReadTimeSeriesFromNetCDF: time unit string has weird format!",BAD_DATA);
            }
            dash  = unit_t_str.substr(21, 1);  // second dash in date
            if ( !strstr(dash.c_str(), "-") ){
              unit_t_str.insert(19,"0");
            }
            colon = unit_t_str.substr(27, 1);  // first colon in time
            if ( !strstr(colon.c_str(), ":") ){
              unit_t_str.insert(25,"0");
            }
            colon = unit_t_str.substr(30, 1);  // second colon in time
            if ( !strstr(colon.c_str(), ":") ){
              unit_t_str.insert(28,"0");
            }
            string sDate = unit_t_str.substr(14,10);
            string sTime = unit_t_str.substr(25,8);
            time_struct tt = DateStringToTimeStruct(sDate, sTime);
            tstep   = (time[1] - time[0])/24./60./60.;
            ExitGracefullyIf(tstep<=0,"CTimeSeries::ReadTimeSeriesFromNetCDF: Interval is negative!",BAD_DATA);
            AddTime(tt.julian_day,tt.year,time[0]*(1./24./60./60.),start_day,start_yr) ;

          }
          ExitGracefullyIf(tstep<=0,"CTimeSeries::ReadTimeSeriesFromNetCDF: this unit in time is not implemented yet (only days, hours, minutes, seconds)",BAD_DATA);
        }
    }
  }

  // -------------------------------
  // (6) set n (= number of datapoints)
  // -------------------------------
  nMeasurements = ntime;

  // -------------------------------
  // (7) determine dimension order in variable
  // -------------------------------
  //     (a) Get the id of the data variable based on its name; varid will be set
  retval = nc_inq_varid(ncid, VarNameNC.c_str(), &varid_f);       HandleNetCDFErrors(retval);
  //     (b) determine in which order the dimensions are in variable
  retval = nc_inq_vardimid(ncid, varid_f, dimids_var);          HandleNetCDFErrors(retval);
  //     (c) determine order
  if ( strcmp(DimNamesNC_stations.c_str(),"None") ) {
    if      ((dimids_var[0] == dimid_x) && (dimids_var[1] == dimid_t))  {  // dimensions are (x,t)
      dim_order = 1;
      dim1 = nstations;
      dim2 = ntime;
    }
    else {                                                                 // dimensions are (t,x)
      dim_order = 2;
      dim1 = ntime;
      dim2 = nstations;
    }
  }
  else {                                                                   // dimensions are (t)
    dim_order = 1;
    dim1 = ntime;
    dim2 = 1;
  }

  // -------------------------------
  // (8) Initialize data array (data type that is needed by NetCDF library)
  // -------------------------------
  double *aVec=new double[dim1*dim2];//stores actual data
  for(int i=0; i<dim1*dim2; i++) {
    aVec[i]=NETCDF_BLANK_VALUE;
  }
  ExitGracefullyIf(aVec==NULL,"CForcingGrid::ReadData : aVec",OUT_OF_MEMORY);
    
  double **aTmp2D=NULL; //stores pointers to rows/columns of 2D data (ntime x nstations)
  double  *aTmp1D=NULL; //stores pointers to rows/columns of 1D data (ntime)

  if ( strcmp(DimNamesNC_stations.c_str(),"None") ) {
    aTmp2D=new double *[dim1];
    ExitGracefullyIf(aTmp2D==NULL,"CForcingGrid::ReadData : aTmp2D(0)",OUT_OF_MEMORY);
    for(int it=0;it<dim1;it++){
      aTmp2D[it]=&aVec[it*dim2]; //points to correct location in aVec data storage
    }
  }
  else {
    aTmp1D=new double [dim1];
    ExitGracefullyIf(aTmp1D==NULL,"CForcingGrid::ReadData : aTmp1D(0)",OUT_OF_MEMORY);
  }

  // -------------------------------
  // (9) Read data
  // -------------------------------
  if ( strcmp(DimNamesNC_stations.c_str(),"None") ) {

    size_t    nc_start [2];
    size_t    nc_length[2];
    ptrdiff_t nc_stride[2];
    
    switch(dim_order) {
    case(1): // dimensions are (x,t)
      nc_start[0]  = StationIdx-1; nc_length[0] = (size_t)(1);      nc_stride[0] = 1;
      nc_start[1]  = 0;            nc_length[1] = (size_t)(ntime);  nc_stride[1] = 1;
      retval=nc_get_vars_double(ncid,varid_f,nc_start,nc_length,nc_stride,&aTmp2D[0][0]);
      HandleNetCDFErrors(retval);
      break;
    case(2): // dimensions are (t,x)
      nc_start[0]  = 0;            nc_length[0] = (size_t)(ntime);  nc_stride[0] = 1;
      nc_start[1]  = StationIdx-1; nc_length[1] = (size_t)(1);      nc_stride[1] = 1;
      retval=nc_get_vars_double(ncid,varid_f,nc_start,nc_length,nc_stride,&aTmp2D[0][0]);
      HandleNetCDFErrors(retval);
      break;
    }
    if (Options.noisy) {
      printf("  Dim of chunk read: dim1 = %i   dim2 = %i\n",dim1,dim2);
      printf("  start  chunk: (%lu, %lu)\n", nc_start[0], nc_start[1]);
      printf("  length chunk: (%lu, %lu)\n", nc_length[0],nc_length[1]);
      printf("  stride chunk: (%lu, %lu)\n", nc_stride[0],nc_stride[1]);
      printf("  Data read: [%.4f, %.4f, %.4f ... %.4f]\n",aTmp2D[0][0], aTmp2D[0][1], aTmp2D[0][2], aTmp2D[0][ntime-1]);
    }
  }
  else {

    size_t    nc_start [1];
    size_t    nc_length[1];
    ptrdiff_t nc_stride[1];
    
    switch(dim_order) {
      case(1):
        nc_start [0] = 0;  nc_length[0] = (size_t)(ntime);  nc_stride[0] = 1;
        retval=nc_get_vars_double(ncid,varid_f,nc_start,nc_length,nc_stride,&aTmp1D[0]);
        HandleNetCDFErrors(retval);
        break;
      }
    if (Options.noisy) {
      printf("  Dim of chunk read: dim1 = %i \n",dim1);
      printf("  start  chunk: (%lu)\n", nc_start[0]);
      printf("  length chunk: (%lu)\n", nc_length[0]);
      printf("  stride chunk: (%lu)\n", nc_stride[0]);
      printf("  Data read: [%.4f, %.4f, %.4f ... %.4f]\n",aTmp1D[0], aTmp1D[1], aTmp1D[2], aTmp1D[ntime-1]);
    }
  }

  // -------------------------------
  // (10) close NetCDF
  // -------------------------------
  retval = nc_close(ncid);
  HandleNetCDFErrors(retval);

  // -------------------------------
  // (11) Initialize data array (data type that is needed by RAVEN)
  // -------------------------------
  double *aVal;
  aVal = NULL;
  aVal =  new double [ntime];
  for (int it=0; it<ntime; it++) {                   // loop over time points in buffer
      aVal[it]=NETCDF_BLANK_VALUE;                   // initialize
  }

  // -------------------------------
  // (12) Convert into RAVEN data array and apply linear transformation: new = a*data+b
  // -------------------------------
  if ( strcmp(DimNamesNC_stations.c_str(),"None") ) {

    // no switch of dim_order required because aTmp2D already points to
    // correct location in aVec data storage (see step 9)

    for (int it=0; it<ntime; it++){                     // loop over time points read
      for (int ic=0; ic<1; ic++){                       // loop over stations
        aVal[it] = LinTrans_a * aTmp2D[ic][it] + LinTrans_b;
      }
    };

  }
  else {

    // no switch of dim_order required because aTmp2D already points to
    // correct location in aVec data storage (see step 9)

    for (int it=0; it<ntime; it++){                     // loop over time points read
      aVal[it]=aTmp1D[it];
    };
  }

  if (Options.noisy) {
    printf("  aVal: [%.4f, %.4f, %.4f ... %.4f]\n",aVal[0], aVal[1], aVal[2], aVal[ntime-1]);
  }

  // -------------------------------
  // (13) add time shift to data
  //      --> only applied when tstep < 1.0 (daily)
  //      --> otherwise ignored and warning written to RavenErrors.txt
  // -------------------------------

  if (tstep >= 1.0) {   // data are not sub-daily
    if ( ceil(TimeShift) == TimeShift) {  // time shift of whole days requested
      // if (TimeShift != 0.0) { cout<<"before: start_day "<<start_day<<"  start_yr "<<start_yr<<endl; }
      AddTime(start_day,start_yr,TimeShift,start_day,start_yr) ;
      // if (TimeShift != 0.0) { cout<<"after:  start_day "<<start_day<<"  start_yr "<<start_yr<<endl; }
    }
    else {  // sub-daily shifts (e.g. 1.25) of daily data requested
      WriteAdvisory("CTimeSeries::ReadTimeSeriesFromNetCDF: time shift specified for NetCDF time series will be ignored", Options.noisy);
      WriteAdvisory("                                       because inputs are daily and time shift is sub-daily", Options.noisy);
    }
  }
  else {  // data are sub-daily
    // if (TimeShift != 0.0) { cout<<"before: start_day "<<start_day<<"  start_yr "<<start_yr<<endl; }
    AddTime(start_day,start_yr,TimeShift,start_day,start_yr) ;
    // if (TimeShift != 0.0) { cout<<"after:  start_day "<<start_day<<"  start_yr "<<start_yr<<endl; }
  }

  // -------------------------------
  // (14) convert to time series object
  // -------------------------------
  is_pulse = true;
  
  pTimeSeries=new CTimeSeries(name,tag,FileNameNC.c_str(),start_day,start_yr,tstep,aVal,nMeasurements,is_pulse);

  // -------------------------------
  // (15) delete dynamic memory
  // -------------------------------
  delete [] aVal;  aVal =NULL;  
#endif   // ends #ifdef _RVNETCDF_
  
  return pTimeSeries;
    
}
