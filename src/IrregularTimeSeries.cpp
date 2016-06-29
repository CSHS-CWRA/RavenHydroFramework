/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright © 2008-2014 the Raven Development Team
----------------------------------------------------------------*/
#include "IrregularTimeSeries.h"
#include "ParseLib.h"
/*****************************************************************
   Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/
///////////////////////////////////////////////////////////////////
/// \brief Implementation of the time series constructor for non-uniformly spaced data points
///
/// \param Name          [in] name of time series
/// \param tag           [in] data tag
/// \param filename      [in] original source file  
/// \param *aValues      [in] Array of time series values [size NumValues]
/// \param *aDays        [in] Array of julian days  [size NumValues]
/// \param *aYears       [in] Array of years [size NumValues]
/// \param NumValues     [in] Number of entries in the time series
//
CIrregularTimeSeries::CIrregularTimeSeries(	string    Name,
											string    tag,
											string    filename,
											double   *aValues,
											double   *aDays,
					 						int      *aYears,
											const int NumValues)
	:CTimeSeriesABC(ts_irregular,Name,tag,filename)
{
	_nVals= NumValues;

	ExitGracefullyIf(NumValues <= 0,
		"CIrregularTimeSeries: Constructor: no entries in time series", BAD_DATA);

	_aVal = NULL;
	_aVal = new double[_nVals];

	_aDays = NULL;
	_aDays = new double[_nVals];

	_aYears = NULL;
	_aYears = new int[_nVals];
	
	ExitGracefullyIf(_aVal==NULL || _aDays==NULL || _aYears==NULL, "CIrregularTimeSeries: Constructor", OUT_OF_MEMORY);

	for (int n = 0; n<_nVals; n++)
	{
		_aVal[n]  = aValues[n];
		_aDays[n] = aDays[n];
		_aYears[n]= aYears[n];
	}

	_aTimes=NULL; //generated in initialize
	_indexCorr=-1;
	_nSampVal=0;
}

///////////////////////////////////////////////////////////////////
/// \brief Implementation of copy constructor (for which the address of a time series is passed)
/// \param &t [in] Address of a time series of which a "copy" is made
//
 CIrregularTimeSeries::CIrregularTimeSeries(string Name,
                          const CIrregularTimeSeries &t)
	  :CTimeSeriesABC(Name,t)
 {
	 _nVals = t.GetNumValues();

	 _aVal = NULL;
	 _aVal = new double[_nVals];

	_aDays = NULL;
	_aDays = new double[_nVals];

	_aYears = NULL;
	_aYears = new int[_nVals];

	ExitGracefullyIf(_aVal==NULL || _aDays==NULL || _aYears==NULL, "CIrregularTimeSeries: Constructor", OUT_OF_MEMORY);

	 for (int n = 0; n<_nVals; n++)
	 {
	   _aVal[n]   = t.GetValue(n);
	   _aDays[n]  = t.GetDay(n);
	   _aYears[n] = t.GetYear(n);
	 }

	_aTimes=NULL; //generated in initialize
	_indexCorr=-1;
	_nSampVal=-1;

 }
///////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
//
CIrregularTimeSeries::~CIrregularTimeSeries()
{
  if (DESTRUCTOR_DEBUG){cout<<"    DELETING IRREGULAR TIME SERIES"<<endl;}
  delete[] _aVal;   _aVal   = NULL;
  delete[] _aDays;  _aDays  = NULL;
  delete[] _aYears; _aYears = NULL;
  delete[] _aTimes; _aTimes = NULL;
}

///////////////////////////////////////////////////////////////////
/// \brief Enables queries of time series values using global time
/// \details Calculates _aTimes, global time for each time series value
/// \details Calculates _indexCorr, correction to global index to get only values within model time.
/// \remark t=0 corresponds to first day with recorded values at that gauge
///
/// \param global_start_day [in] Julian global start day
/// \param global_start_year [in] Global start year
/// \param model_duration [in] Duration of model, in days
/// \param timestep [in] mdoel timestep, in days
/// \param is_observation [in] - true if this is an observation time series, rather than model input
void CIrregularTimeSeries::Initialize(const double global_start_day,
                                      const    int global_start_year,
                                      const double model_duration,
                                      const double timestep,
                                      const   bool is_observation)
{
  double time;
  int yr, leap;

  _aTimes = new double[_nVals];
  for (int n = 0; n < _nVals; n++)
  {
    time = _aDays[n] - global_start_day;
    yr = _aYears[n];

    if (IsLeapYear(yr)){time-=1;} //leap day already added to julian day for first year

	  while (yr!=global_start_year)
    {
      leap = 0;
      if (IsLeapYear(yr)){ leap = 1; }
      if (yr > global_start_year)
      {
        time += (365 + leap);
        yr--;
      }
      else
      {
        time -= (365 + leap);
        yr++;
      }
    }
    _aTimes[n] = time;
    if (_indexCorr == -1 && time >= 0){ _indexCorr = n; }
    if (time >= 0 && time < model_duration){ _nSampVal++; }
  }
  if (!is_observation)
  {
    ExitGracefullyIf(_nSampVal <= 0,
      "CIrregularTimeSeries::Initialize: time series data not available during model simulation", BAD_DATA);
  }
}

///////////////////////////////////////////////////////////////////
/// \brief Returns number of values
/// \return Number of values
//
int    CIrregularTimeSeries::GetNumValues() const{ return _nVals; }


///////////////////////////////////////////////////////////////////
/// \brief Returns the time of the time series data point for which n is an index
/// \param n [in] Index
/// \return time of time series data point for which n is an index
//
double CIrregularTimeSeries::GetTime(const int n) const{return _aTimes[n];}

///////////////////////////////////////////////////////////////////
/// \brief Returns magnitude of time series data point for which n is an index
/// \param n [in] Index
/// \return Magnitude of time series data point for which n is an index
//
double CIrregularTimeSeries::GetValue(const int n)const{return _aVal[n];}

///////////////////////////////////////////////////////////////////
/// \brief Returns the julian day of index n
/// \param n [in] Index
/// \return julian day of the timeseries at index n
//
double CIrregularTimeSeries::GetDay(const int n) const{return _aDays[n];}
///////////////////////////////////////////////////////////////////
/// \brief Returns the year of timeseries ar index n
/// \param n [in] Index
/// \return year of the timeseries at index n
//
int CIrregularTimeSeries::GetYear(const int n) const{return _aYears[n];}
///////////////////////////////////////////////////////////////////
/// \brief Returns average of the values in interval t to t+tstep
/// \remark returns BLANK_DATA if no values are found in the interval
/// \param &t [in] Left bound of (global) time interval over which average value of time series is to be determined
/// \param &tstep [in] Duration of interval over which average value of time series is to be determined
/// \return Average of time series values across interval t to t+tstep
//
double CIrregularTimeSeries::GetAvgValue(const double &t, const double &tstep) const
{
	double sum = 0;
	int count = 0;
	for (int n = 0; n < _nVals; n++)
	{
	  if ( _aTimes[n]>=t && _aTimes[n]<=(t+tstep) )
	  {
		sum += _aVal[n];
		count++;
	  }
	}
	if (count > 0){ return sum / count; }
	else { return BLANK_DATA; }
}

///////////////////////////////////////////////////////////////////
/// \brief Returns minimum value of time series between time t to t+tstep
/// \remark returns BLANK_DATA if no values are found in the interval
/// \param &t [in] Left endpoint of time interval over which minimum value of time series is to be determined
/// \param &tstep [in] Duration of time interval over which minimum time series value is to be determined
/// \return Minimum value of time series over interval t to t+tstep
//
double CIrregularTimeSeries::GetMinValue(const double &t, const double &tstep) const
{
	double vmin(ALMOST_INF);
	bool blank = true;
	for (int n = 0; n < _nVals; n++)
	{
		if (_aTimes[n] >= t && _aTimes[n] <= (t + tstep))
		{
			lowerswap(vmin, _aVal[n]);
			blank = false;
		}
	}
	if (blank){ return BLANK_DATA; }

	return vmin;
}

///////////////////////////////////////////////////////////////////
/// \brief Returns maximum value of time series between time t to t+tstep
/// \remark returns BLANK_DATA if no values are found in the interval
/// \param &t [in] Left endpoint of time interval over which maximum value of time series is to be determined
/// \param &tstep [in] Duration of time interval over which maximum time series value is to be determined
/// \return Maximum value of time series over interval t to t+tstep
//
double CIrregularTimeSeries::GetMaxValue(const double &t, const double &tstep) const
{
	double vmax(-ALMOST_INF);
	bool blank = true;
	for (int n = 0; n < _nVals; n++)
	{
		if (_aTimes[n] >= t && _aTimes[n] <= (t + tstep))
		{
			upperswap(vmax, _aVal[n]);
			blank = false;
		}
	}
	if (blank){ return BLANK_DATA; }

	return vmax;
}

///////////////////////////////////////////////////////////////////
/// \brief Returns 0 since time series is assumed to be instantaneous
/// \return 0.0
//
double CIrregularTimeSeries::GetSampledInterval() const{return 0.0;}
double CIrregularTimeSeries::GetInterval() const{return 0.0;}

///////////////////////////////////////////////////////////////////
/// \brief Returns nnth value within the model time
/// \notes must be called after initializing
///  
/// \param nn [in] index
/// \return time series value
//
double CIrregularTimeSeries::GetSampledValue(const int nn) const{return _aVal[nn+_indexCorr];}

///////////////////////////////////////////////////////////////////
/// \brief Returns time of the nnth value within the model time
/// \notes must be called after initializing
///  
/// \param nn [in] index
/// \return time
//
double CIrregularTimeSeries::GetSampledTime(const int nn) const{return _aTimes[nn+_indexCorr];}

///////////////////////////////////////////////////////////////////
/// \brief Returns number of values within model time
/// \notes must be called after initializing
///  
/// \return number of values
//
int    CIrregularTimeSeries::GetNumSampledValues() const{return _nSampVal;}

///////////////////////////////////////////////////////////////////
/// \brief Parses standard single time series format from file and creates time series object
/// \param *p [in] CParser object pointing to input file 
/// \return Pointer to created time series
//
CIrregularTimeSeries  *CIrregularTimeSeries::Parse (CParser *p, const string name, const string tag, const int nMeasurements)
{
  char *s[MAXINPUTITEMS];
  int Len;
  double *aDays;
  int    *aYears;
  double *aVal;

  aVal = new double [nMeasurements];
  aDays= new double [nMeasurements];
  aYears=new int    [nMeasurements];
  ExitGracefullyIf(aVal==NULL || aDays==NULL || aYears==NULL, "CIrregularTimeSeries: Parse", OUT_OF_MEMORY);

  int n=0;
  s[0]="";
  while ((n<nMeasurements) && (strcmp(s[0],"&")) && (!p->Tokenize(s,Len)))
  {
		if (Len<3){
		  p->ImproperFormat(s); cout <<"Length:" <<Len<<endl;
		  ExitGracefully("CIrregularTimeSeries::Parse: Bad number of time series points",BAD_DATA);
		}
		ExitGracefullyIf(n>=nMeasurements,"CIrregularTimeSeries::Parse: Bad number of time series points",BAD_DATA);

		if ((string(s[0]).length()==10) &&
		   ((string(s[0]).substr(4,1)=="/") || 
			(string(s[0]).substr(4,1)=="-")))
		{ //in timestamp format [yyyy-mm-dd] [hh:mm:ss.0]
			time_struct tt;
			tt=DateStringToTimeStruct(string(s[0]),string(s[1]));
			aDays[n]=tt.julian_day;
			aYears[n] =tt.year;
		} 
		else
		{ //julian date format [nMeasurements] [start_day]
			aDays [n]=s_to_d(s[0]);
			aYears[n]=s_to_i(s[1]);
		}

		aVal [n]=s_to_d(s[2]);
		n++;
  }
  if (n!=nMeasurements){
    ExitGracefully("CIrregularTimeSeries::Parse: Bad number of time series points",BAD_DATA);}
  p->Tokenize(s,Len);//read closing term (e.g., ":EndRain")

  CIrregularTimeSeries *pTimeSeries=NULL;
  pTimeSeries=new CIrregularTimeSeries(name,tag,p->GetFilename(),aVal,aDays,aYears,nMeasurements);
  delete [] aVal;  aVal =NULL;
  delete [] aDays; aDays =NULL;
  delete [] aYears; aYears =NULL;
  return pTimeSeries;
}