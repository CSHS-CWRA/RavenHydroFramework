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
  _aSampVal =NULL;//generated in Resample() routine
  _nSampVal =0;//generated in Resample() routine
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
CTimeSeries::CTimeSeries(string    Name,
                         string    tag,
                         string    filename,
                         double    strt_day,
                         int   start_yr,
                         double    data_interval,
                         double   *aValues,
                         const int NumPulses,
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
CTimeSeries::CTimeSeries(string    Name,
                         string    tag,
                         string    filename,
                         double    strt_day,
                         int   start_yr,
                         double    data_interval,
                         const int NumPulses,
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
    _aVal [n]=BLANK_DATA;
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
    _aSampVal[nn] = BLANK_DATA;
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

  //if ((n < 0) ||(n>_nPulses - 1)){ return CTimeSeries::BLANK_DATA;}

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
  if (t_loc < -TIME_CORRECTION) {return CTimeSeriesABC::BLANK_DATA; }
  if (t_loc >= _nPulses*_interval) {return CTimeSeriesABC::BLANK_DATA; }
  if (_pulse){
    if (n1 == n2){ return _aVal[n1]; }
    else{
      sum = 0;
      inc = ((double)(n1 + 1)*_interval - t_loc);
      if (_aVal[n1] == CTimeSeriesABC::BLANK_DATA)  { blank += inc; }
      else                                          { sum += _aVal[n1] * inc; }

      for (int n = n1 + 1; n < n2; n++){
        if (_aVal[n] == CTimeSeriesABC::BLANK_DATA) { blank += _interval; }
        else                                        { sum += _aVal[n] * _interval; }
      }
      inc = ((t_loc + tstep) - (double)(n2)*_interval);
      if (_aVal[n2] == CTimeSeriesABC::BLANK_DATA)  { blank += inc; }
      else                                          { sum += _aVal[n2] * inc; }
    }
  }
  else{
    ExitGracefully("CTimeSeries::GetAvgValue (non-pulse)", STUB);
  }
  if (blank / tstep > 0.001){return BLANK_DATA;}

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
    return CTimeSeriesABC::BLANK_DATA;
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
//
CTimeSeries *CTimeSeries::Parse(CParser *p, bool is_pulse, string name, string tag, const optStruct &Options, bool shift_to_per_ending)
{

  char *s[MAXINPUTITEMS];
  int Len;
  double start_day;
  int    start_yr;
  double tstep;
  int    nMeasurements;

  CTimeSeries *pTimeSeries=NULL;

  p->Tokenize(s,Len);
  if (IsComment(s[0],Len)){p->Tokenize(s,Len);}//try again
  if (Len<4){p->ImproperFormat(s); cout <<"Length:" <<Len<<endl;}
  
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
      //aVal [n]=s_to_d(s[i]);
			n++;
    }
  }

  if (n!=nMeasurements){
    cout<<p->GetFilename()<<endl;
    cout << " n | nMeasurements " << n << " "<<nMeasurements<<endl;
    ExitGracefully("CTimeSeries::Parse: Bad number of time series points",BAD_DATA);}

  p->Tokenize(s,Len);//read closing term (e.g., ":EndRain")
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
