/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2019 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Gauge.h"

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the Gauge constructor
/// \details Initialize all time series to NULL and average Temp and PET to NOT_SPECIFIED
///
/// \param gauge_name [in] Gauge name
/// \param latit [in] Latitude of gauge location [degrees]
/// \param longit [in] Longitude of gauge location [degrees]
/// \param elev [in] Elevation of gauge [masl]
//
CGauge::CGauge(string gauge_name,
               const double latit, //[deg]
               const double longit,//[deg]
               const double elev)  //[masl]
{
  _Loc.latitude  =latit;
  _Loc.longitude =longit;
  _Loc.UTM_x                    =0.0;
  _Loc.UTM_y                    =0.0;//set in initialize

  _elevation     =elev;
  _meas_ht       =2.0; // [m], default
  _rainfall_corr =1.0;
  _snowfall_corr =1.0;
  _cloud_min_temp=-20.0;//ensures cloud-free status always unless overriden
  _cloud_max_temp= 40.0;


  _name=gauge_name;

  _nTimeSeries=0;
  _pTimeSeries  =NULL;
  for (int i=0;i<MAX_FORCING_TYPES;i++){
    _aTSindex[i]=DOESNT_EXIST;
  }

  for (int i=0;i<12;i++)
  {
    _aAveTemp[i]=NOT_SPECIFIED;
    _aMinTemp[i]=NOT_SPECIFIED;
    _aMaxTemp[i]=NOT_SPECIFIED;
    _aAvePET [i]=NOT_SPECIFIED;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
/// \details deletes all associated time series objects
//
CGauge::~CGauge()
{
  if (DESTRUCTOR_DEBUG){cout<<"  DELETING GAUGE"<<endl;}

  for (int i=0;i<_nTimeSeries;i++){delete _pTimeSeries[i];_pTimeSeries[i]=NULL;}
  delete [] _pTimeSeries; _pTimeSeries=NULL;
}

//////////////////////////////////////////////////////////////////
/// \brief Performs all operations required before the simulation begins
/// \details Calculates UTM location of gauge, initializes all time series
/// \param &Options [in] Global model options information
/// \param UTM_zone [in] UTM zone in which the gauge lies
//
void CGauge::Initialize(const optStruct   &Options,
                        const int          UTM_zone)
{
  LatLonToUTMXY(_Loc.latitude,_Loc.longitude,
                UTM_zone,
                _Loc.UTM_x,   _Loc.UTM_y);

  //Populate Temperature time series: by timestep, daily min/max/average
  //--------------------------------------------------------------------------
  if((GetTimeSeries(F_TEMP_AVE)!=NULL) || (GetTimeSeries(F_TEMP_DAILY_AVE)!=NULL) ||
    ((GetTimeSeries(F_TEMP_DAILY_MAX)!=NULL) && (GetTimeSeries(F_TEMP_DAILY_MIN)!=NULL)))
  {
    if((GetTimeSeries(F_TEMP_AVE)!=NULL) && (GetTimeSeries(F_TEMP_AVE)->IsDaily()))
    {
      // if daily temp is specified, copy to temp_daily_ave time series
      CTimeSeries *temp_daily_ave=new CTimeSeries("TEMP_DAILY_AVE",*GetTimeSeries(F_TEMP_AVE));
      AddTimeSeries(temp_daily_ave,F_TEMP_DAILY_AVE);
    }

    if((GetTimeSeries(F_TEMP_AVE)!=NULL) && (!GetTimeSeries(F_TEMP_AVE)->IsDaily()))//Sub-daily temperature data provided
    {
      GenerateMinMaxAveTempFromSubdaily(Options);//Generate daily min,max, & avg
    }
    else if((GetTimeSeries(F_TEMP_DAILY_MAX)!=NULL) && //Min/max temperature data provided
      (GetTimeSeries(F_TEMP_DAILY_MIN)!=NULL))
    {
      GenerateAveSubdailyTempFromMinMax(Options);//Generate T_daily_ave, T_ave (downscaled)
    }
    else if(GetTimeSeries(F_TEMP_DAILY_AVE)!=NULL) //only daily average data provided
    {
      GenerateMinMaxSubdailyTempFromAve(Options);
    }
  }

  bool hasPrecip  =(GetTimeSeries(F_PRECIP)!=NULL);
  bool hasSnowfall=(GetTimeSeries(F_SNOWFALL)!=NULL);
  bool hasRainfall=(GetTimeSeries(F_RAINFALL)!=NULL);

  // Handle snowfall availability
  //--------------------------------------------------------------------------
  ExitGracefullyIf((hasPrecip || hasRainfall) && (!hasSnowfall) && (Options.rainsnow==RAINSNOW_DATA),
                   "CGauge::Initialize: snow autogeneration is off at gauge with precip data, but no snow data has been supplied.",BAD_DATA);

  if(hasSnowfall && hasRainfall && (Options.rainsnow!=RAINSNOW_DATA))
  {
    WriteWarning("Gauge:Initialize: both snowfall and rainfall data are provided at a gauge, but :RainSnowFraction method is something other than RAINSNOW_DATA. Snow fraction will be recalculated.",Options.noisy);
  }

  if (hasSnowfall && hasRainfall){
    CTimeSeries *sum;
    sum=CTimeSeries::Sum(GetTimeSeries(F_RAINFALL),GetTimeSeries(F_SNOWFALL),"PRECIP");//sum together rain & snow to get precip
    AddTimeSeries(sum,F_PRECIP);
  }
  if (!hasRainfall && !hasSnowfall && hasPrecip){
    AddTimeSeries(new CTimeSeries("RAINFALL",*GetTimeSeries(F_PRECIP)),F_RAINFALL); //if no snow or rain, copy precip to rainfall
  }
  if (!hasSnowfall && (hasRainfall || hasPrecip))
  {
    AddTimeSeries(new CTimeSeries("SNOWFALL",DOESNT_EXIST,0.0),F_SNOWFALL);//blank series, all 0.0s
  }


  double model_start_day=Options.julian_start_day;
  int    model_start_yr =Options.julian_start_year;
  double model_duration =Options.duration;
  double timestep       =Options.timestep;

  /*for (int i=0; i<MAX_FORCING_TYPES;i++){
    if (aTSindex[i]!=DOESNT_EXIST){
    cout<<" Forcing: "<<ForcingToString((forcing_type)(i))<<endl;
    CTimeSeries *pTS=pTimeSeries[aTSindex[i]];
    cout<<pTS->GetNumValues()<<" "<<pTS->GetInterval()<<" "<<pTS->GetStartDay()<<endl;
    }
    }*/
  for (int i=0; i<_nTimeSeries;i++){
    _pTimeSeries[i]->Initialize(model_start_day,model_start_yr,model_duration,timestep,false,Options.calendar);
  }

  //QA/QC: check time series for valid values
  //--------------------------------------------------------------------------
  double val;
  int index;
  int nSamples=(int)(ceil(model_duration/Options.timestep-TIME_CORRECTION));
  hasPrecip  =(GetTimeSeries(F_PRECIP)!=NULL);
  if(hasPrecip){
    //precip greater than zero, daily less than 2000 mm/d (world record = 1825 mm/d)
    index=_aTSindex[(int)(F_PRECIP)];
    
    if(index!=DOESNT_EXIST){
      for(int nn=0;nn<nSamples; nn++)
      {
        val=_pTimeSeries[index]->GetSampledValue(nn);
        if(val==RAV_BLANK_DATA){
          ExitGracefully("CGauge::Initialize: Raven cannot have blank data in precipitation time series",BAD_DATA);
        }
        if((val<-1e-6) || (val>10000)){
          cout<<GetName()<<" "<<nn<<" "<<val<<endl;
          ExitGracefully("CGauge::Initialize: negative or excessively large (>10000mm/d) precipitation intensity reported at gauge",BAD_DATA);
        }
      }
    }
  }
  //temp greater than -60C less than 60C
  bool hasAveTemp      =(GetTimeSeries(F_TEMP_DAILY_AVE)!=NULL);
  bool derivedAveTemp  =((GetTimeSeries(F_TEMP_DAILY_MIN)!=NULL) && (GetTimeSeries(F_TEMP_DAILY_MAX)!=NULL));
  if(hasAveTemp){
    index=_aTSindex[(int)(F_TEMP_DAILY_AVE)];
    if(index!=DOESNT_EXIST){
      for(int nn=0;nn<nSamples; nn++)
      {
        val=_pTimeSeries[index]->GetSampledValue(nn);
	      if (val==RAV_BLANK_DATA && !(derivedAveTemp)){
          string warning;
          warning ="CGauge::Initialize: Raven cannot have blank data in daily temperature time series (Gauge: "+_name+", n="+to_string(nn)+")";
          ExitGracefully(warning.c_str(),BAD_DATA);
        }
        else if((val<-60) || (val>60)){
          ExitGracefully("CGauge::Initialize: excessively small or large average temperature (<-60C or >60C) reported at gauge",BAD_DATA);
        }
      }
    }
  }
  //PET unreasonable
  bool blank_in_PET=false;
  index=_aTSindex[(int)(F_PET)];
  if (index!=DOESNT_EXIST){
    for (int nn=0;nn<nSamples; nn++)
    {
      val=_pTimeSeries[index]->GetSampledValue (nn);
      if(val==RAV_BLANK_DATA) {
        blank_in_PET=true;
      }
      else if (val<-REAL_SMALL){
        ExitGracefully("CGauge::Initialize: negative PET reported at gauge",BAD_DATA);
      }
    }
    if(blank_in_PET) {
      WriteWarning("CGauge::Initialize: blank PET data in PET time series; will infill via infill ET estimation method",Options.noisy);
    }
  }
  //Positive non-blank wind velocity
  index=_aTSindex[(int)(F_WIND_VEL)];
  if (index!=DOESNT_EXIST){
    for (int nn=0;nn<nSamples; nn++)
    {
      val=_pTimeSeries[index]->GetSampledValue (nn);
      if(val==RAV_BLANK_DATA) {
        ExitGracefully("CGauge::Initialize: blank wind speed reported at gauge",BAD_DATA);
      }
      else if (val<-REAL_SMALL){
        ExitGracefully("CGauge::Initialize: negative wind speed reported at gauge",BAD_DATA);
      }
    }
  }

  WarnAboutForcing(Options.SW_radiation  ==SW_RAD_DATA,   F_SW_RADIA);
  WarnAboutForcing(Options.SW_radia_net  ==NETSWRAD_DATA, F_SW_RADIA_NET);
  WarnAboutForcing(Options.LW_radiation  ==LW_RAD_DATA,   F_LW_RADIA_NET);
  WarnAboutForcing(Options.cloud_cover   ==CLOUDCOV_DATA, F_CLOUD_COVER);
  WarnAboutForcing(Options.ow_evaporation==PET_DATA,      F_OW_PET);
  WarnAboutForcing(Options.evaporation   ==PET_DATA,      F_PET);
  WarnAboutForcing(Options.pot_melt      ==POTMELT_DATA,  F_POTENTIAL_MELT);
  WarnAboutForcing(Options.rel_humidity  ==RELHUM_DATA,   F_REL_HUMIDITY);
  WarnAboutForcing(Options.air_pressure  ==AIRPRESS_DATA, F_AIR_PRES);
  WarnAboutForcing(Options.wind_velocity ==WINDVEL_DATA,  F_WIND_VEL);

  // Check for monthly values, when needed
  //--------------------------------------------------------------------------
  ExitGracefullyIf((_aAveTemp[0]==NOT_SPECIFIED) && ((Options.evaporation==PET_FROMMONTHLY)),
                   "CGauge::Initialize: monthly temps for gauge not specified, but are needed",BAD_DATA);
  ExitGracefullyIf((_aAvePET[0]==NOT_SPECIFIED) && ((Options.evaporation==PET_FROMMONTHLY) || (Options.evaporation==PET_MONTHLY_FACTOR)),
                   "CGauge::Initialize: monthly PET values for gauge not specified, but are needed",BAD_DATA);

}
//////////////////////////////////////////////////////////////////
/// \brief Gets time series representing forcing data type ftype
/// \param ftype [in] forcing data type
/// \returns time series representing forcing data type ftype, if time series exists, NULL otherwise
//
CTimeSeries *CGauge::GetTimeSeries(const forcing_type ftype) const
{
  int index=_aTSindex[(int)(ftype)];
  if (index!=DOESNT_EXIST){return _pTimeSeries[index];}
  return NULL;
}
//////////////////////////////////////////////////////////////////
/// \brief returns true if forcing exists at gauge
/// \param ftype [in] forcing data type
/// \returns true if forcing exists at gauge
//
bool     CGauge::TimeSeriesExists(const forcing_type ftype) const
{
  return (_aTSindex[(int)(ftype)]!=DOESNT_EXIST);
}
//////////////////////////////////////////////////////////////////
/// \brief Checks if warnings about forcings are needed at this gauge, warns if required
/// \param is_needed [in] true if specific forcing data is needed
/// \param ftype [in] forcing data type
//
void CGauge::WarnAboutForcing(bool is_needed, const forcing_type ftype) const
{
  string errmsg;
  string warning;

  //check for FROM_DATA without data or vice versa
  if ((is_needed) && (GetTimeSeries(ftype) == NULL)){
    errmsg = "USE_DATA is specified as the means of generating "+ForcingToString(ftype) +", but no data is provided at gauge " + this->_name;
    ExitGracefully(errmsg.c_str(),BAD_DATA_WARN);
  }
  else if ((GetTimeSeries(ftype) != NULL) && (!is_needed)){
    warning = ForcingToString(ftype) +" data supplied at gauge " + this->_name + " but will not be used due to choice of forcing generation algorithm";
    WriteWarning(warning,false);
  }
}
/*****************************************************************
   Manipulator functions
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Sets gauge latitude
/// \param &lat [in] Latitude of gauge [degrees]
//
void     CGauge::SetLatitude          (const double &lat){_Loc.latitude =lat;}

//////////////////////////////////////////////////////////////////
/// \brief Sets gauge longitude
/// \param &lon [in] Longitude of gauge [degrees]
//
void     CGauge::SetLongitude         (const double &lon){_Loc.longitude=lon;}

//////////////////////////////////////////////////////////////////
/// \brief Sets gauge elevation
/// \param &e [in] Elevation of gauge [masl]
//
void     CGauge::SetElevation         (const double &e){_elevation=e;}

//////////////////////////////////////////////////////////////////
/// \brief Sets gauge height w.r.t. land surface (default=2.0m)
/// \param &ht [in] Height of gauge [m]
//
void     CGauge::SetMeasurementHt     (const double &ht){_meas_ht=ht;}

//////////////////////////////////////////////////////////////////
/// \brief Sets gauge property
/// \param prop_tag [in] Property Identifier (string)
/// \param &value [in] Set value
/// \return True if operation was successful, false otherwise (unrecognized string)
//
bool     CGauge::SetGaugeProperty          (const string prop_tag, const double &value)
{
  string label_n = StringToUppercase(prop_tag);
  if      (!label_n.compare("RAINFALL_CORR"   ))  {_rainfall_corr=value;}
  else if (!label_n.compare("SNOWFALL_CORR"   ))  {_snowfall_corr=value;}
  else if (!label_n.compare("ELEVATION"       ))  {_elevation=value;}
  else if (!label_n.compare("CLOUD_MIN_RANGE" ))  {_cloud_min_temp=value;}
  else if (!label_n.compare("CLOUD_MAX_RANGE" ))  {_cloud_max_temp=value;}
  else if (!label_n.compare("LATITUDE"        ))  {_Loc.latitude=value;}
  else if (!label_n.compare("LONGITUDE"       ))  {_Loc.longitude=value;}
  else{
    WriteWarning("CGauge::SetGaugeProperty: unrecognized gauge property "+prop_tag,false);
    return false;//bad string
  }
  return true;
}

/*****************************************************************
   Add Precip, Temperature Time Series
*****************************************************************/
//////////////////////////////////////////////////////////////////
/// \brief Loads minimum temperature time series into gauge
/// \param *pTS [in]time series to be loaded
/// \param ftype [in] data type represented by time series
//
void CGauge::AddTimeSeries      (CTimeSeries *pTS, forcing_type ftype)
{
  ExitGracefullyIf(pTS==NULL,"AddTempTimeSeries::NULL time series added",BAD_DATA);
  int index=_aTSindex[(int)(ftype)];
  if (index!=DOESNT_EXIST)//overwriting existing time series
  {
    cout <<"ftype : "<<ftype<<" Forcing: "<<ForcingToString(ftype)<<endl;
    string warn="CGauge::AddTimeSeries: a time series of data has been overwritten at gauge "+_name;
    WriteWarning(warn,true);
    delete _pTimeSeries[index]; _pTimeSeries[index]=pTS;
  }
  else{
    if (!DynArrayAppend((void**&)(_pTimeSeries),(void*)(pTS),_nTimeSeries)){
      ExitGracefully("CGauge::AddTimeSeries: adding NULL time series",BAD_DATA);}
    _aTSindex[(int)(ftype)]=_nTimeSeries-1;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets monthly temperature array to temps[]
/// \param temps[] [in] New monthly temperature array
//
void CGauge::SetMonthlyAveTemps(const double temps[12])
{
  for (int i=0;i<12;i++){_aAveTemp[i]=temps[i];}
}

//////////////////////////////////////////////////////////////////
/// \brief Sets monthly maximum temperature array to temps[]
/// \param temps[] [in] New monthly temperature array
//
void CGauge::SetMonthlyMaxTemps(const double temps[12])
{
  for (int i=0;i<12;i++){_aMaxTemp[i]=temps[i];}
}

//////////////////////////////////////////////////////////////////
/// \brief Sets monthly minimum temperature array to temps[]
/// \param temps[] [in] New monthly temperature array
//
void CGauge::SetMonthlyMinTemps(const double temps[12])
{
  for (int i=0;i<12;i++){_aMinTemp[i]=temps[i];}
}
//////////////////////////////////////////////////////////////////
/// \brief Sets monthly PET array to PET[]
/// \param PET[] [in] New monthly PET array
//
void CGauge::SetMonthlyPET (const double PET[12])
{
  for (int i=0;i<12;i++){_aAvePET[i]=PET[i];}
}
/*****************************************************************
   Accessor Functions
   the time t is always model time (i.e., 0 at start of model)
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Returns name of gauge
/// \return name of gauge
//
string   CGauge::GetName                        () const {return _name;}

//////////////////////////////////////////////////////////////////
/// \brief Returns location of gauge
/// \return Location of gauge
//
location CGauge::GetLocation                    () const {return _Loc;}

//////////////////////////////////////////////////////////////////
/// \brief Returns elevation of gauge
/// \return Elevation of gauge [masl]
//
double   CGauge::GetElevation     () const {return _elevation;}
//////////////////////////////////////////////////////////////////
/// \brief Returns measurement height of gauge above land surface
/// \return Height of gauge [m]
//
double   CGauge::GetMeasurementHt () const {return _meas_ht;}
//----------------------------------------------------------------
double   CGauge::GetRainfallCorr  () const {return _rainfall_corr;}
//----------------------------------------------------------------
double   CGauge::GetSnowfallCorr  () const {return _snowfall_corr;}
//----------------------------------------------------------------
double   CGauge::GetCloudMinRange () const {return _cloud_min_temp;}
//----------------------------------------------------------------
double   CGauge::GetCloudMaxRange () const {return _cloud_max_temp;}

//////////////////////////////////////////////////////////////////
/// \brief Return hourly temperature correction at time t
/// \param t [in] Time at which daily temperature correction is to be determined
/// \return Daily temperature correction at time t
//
double CGauge::DailyTempCorrection(const double t) const
{
  return -cos(2.0*PI*(t-PEAK_TEMP_HR/HR_PER_DAY));
}

//////////////////////////////////////////////////////////////////
/// \brief Returns average fraction of snow in precipitation between time t and t+tstep
/// \param &t [in] Left bound of interval over which average fraction of snow in precipitation will be calculated
/// \param &tstep [in] Duration of interval over which average fraction of snow in precipitation will be calculated
/// \return average fraction of snow in precipitation betweeb time t and t+tstep
//
double CGauge::GetAverageSnowFrac       (const double &t, const double &tstep) const
{
  double rain=GetTimeSeries(F_RAINFALL)->GetAvgValue(t,tstep);
  double snow=GetTimeSeries(F_SNOWFALL)->GetAvgValue(t,tstep);
  if ((snow+rain)==0){return 0.0;}
  return snow/(snow+rain);
}
double CGauge::GetAverageSnowFrac       (const int nn) const
{
  double rain=GetTimeSeries(F_RAINFALL)->GetSampledValue(nn);
  double snow=GetTimeSeries(F_SNOWFALL)->GetSampledValue(nn);
  if ((snow+rain)==0){return 0.0;}
  return snow/(snow+rain);
}
//////////////////////////////////////////////////////////////////
/// \brief Gets gauge property
/// \param pname [in] parameter name (string)
/// \return parameter value if pname is valid, DOESNT_EXIST otherwise (unrecognized string)
//
double     CGauge::GetGaugeProperty          (const string &pname) const
{
  string label_n = StringToUppercase(pname);
  if      (!label_n.compare("RAINFALL_CORR"   ))  { return _rainfall_corr;}
  else if (!label_n.compare("SNOWFALL_CORR"   ))  { return _snowfall_corr;}
  else if (!label_n.compare("ELEVATION"       ))  { return _elevation;}
  else if (!label_n.compare("CLOUD_MIN_RANGE" ))  { return _cloud_min_temp;}
  else if (!label_n.compare("CLOUD_MAX_RANGE" ))  { return _cloud_max_temp;}
  else if (!label_n.compare("LATITUDE"        ))  { return _Loc.latitude;}
  else if (!label_n.compare("LONGITUDE"       ))  { return _Loc.longitude;}
  else{
    WriteWarning("CGauge::GetGaugeProperty: unrecognized gauge property "+pname,false);
    return DOESNT_EXIST;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Returns average forcing function value between time t and t+tstep
/// \param f [in] desired forcing function
/// \param &t [in] Left bound of interval over which average will be calculated
/// \param &tstep [in] Duration of interval over which average will be calculated
/// \return average forcing function value between time t and t+tstep
//
double CGauge::GetForcingValue         (const forcing_type ftype, const double &t, const double &tstep) const
{
  CTimeSeries *pTS;
  pTS=GetTimeSeries(ftype);
  if (pTS==NULL){return 0.0;}
  return pTS->GetAvgValue(t,tstep);
}
//////////////////////////////////////////////////////////////////
/// \brief Returns average forcing function value for time step nn
/// \param f [in] desired forcing function
/// \param nn [in] time step index
/// \return average forcing function value over time step nn
//
double CGauge::GetForcingValue         (const forcing_type ftype, const int nn) const
{
  CTimeSeries *pTS;
  pTS=GetTimeSeries(ftype);
  if (pTS==NULL){return 0.0;}

  return pTS->GetSampledValue(nn);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns representative average temperature over month
/// \param month [in] Month (from 1-12)
/// \return Representative Average temperature for month
//
double   CGauge::GetMonthlyAveTemp  (const int month) const
{
  return _aAveTemp[month-1];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns representative minimum temperature over month
/// \param month [in] Month (from 1-12)
/// \return Representative minimum temperature for month
//
double   CGauge::GetMonthlyMinTemp  (const int month) const
{
  return _aMinTemp[month-1];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns representative maximum temperature over month
/// \param month [in] Month (from 1-12)
/// \return Representative maximum temperature for month
//
double   CGauge::GetMonthlyMaxTemp  (const int month) const
{
  return _aMaxTemp[month-1];
}
//////////////////////////////////////////////////////////////////
/// \brief Returns average PET over month
/// \param month [in] Month for which the average PET is to be determined
/// \return Average PET over month
//
double   CGauge::GetMonthlyAvePET   (const int month) const
{
  return _aAvePET[month-1];
}

//////////////////////////////////////////////////////////////////
/// \brief Converts string (e.g., "X2305" in Gauge file) to gauge
/// \remarks if string is invalid, returns NULL
/// \return gauge corresponding to input string
/*const CGauge*CGauge::StringToGauge(const string s)
  {
  string sup=StringToUppercase(s);
  for (int p=0;p<NumGauges;p++)
  {
  if (!sup.compare(StringToUppercase(pAllNumGauges[p]->GetTag()))){return pAllNumGauges[p];}
  //else if (s_to_i(s.c_str())==(p+1))          {return pAllNumGauges[p];}
  }
  return NULL;
  }*/

//////////////////////////////////////////////////////////////////
/// \brief Generates daily Tmin,Tmax,Tave time series from T (subdaily) time series
/// \note presumes existence of valid F_TEMP_AVE time series with subdaily timestep
///    for incomplete days of TEMP_AVE data, min/max/average is estimated from available data
//
void CGauge::GenerateMinMaxAveTempFromSubdaily(const optStruct &Options)
{
  CTimeSeries *pT;
  pT=GetTimeSeries(F_TEMP_AVE);

  bool incomp=false;
  double start_day=Options.julian_start_day; 
  int    start_yr =Options.julian_start_year;
  double duration =Options.duration;
  double timestep =Options.timestep;

  //below needed for correct mapping from time series to model time
  pT->Initialize(start_day,start_yr,duration,timestep,false,Options.calendar);

  double time_shift=Options.julian_start_day-floor(Options.julian_start_day+TIME_CORRECTION);
  int nVals=(int)ceil(duration+time_shift);//must provide full days overlapping model duration
  double *aMin=new double [nVals];
  double *aMax=new double [nVals];
  double *aAvg=new double [nVals];
  
  double t=0.0;//model time
  for (int n=0;n<nVals;n++){
    aMin[n]=pT->GetMinValue(t-time_shift,1.0); //t-time_shift corresponds to 00:00 on day of model time t
    aMax[n]=pT->GetMaxValue(t-time_shift,1.0);
    aAvg[n]=pT->GetAvgValue(t-time_shift,1.0);
    if     ((aAvg[n]==RAV_BLANK_DATA) && ((t-time_shift    )<0       )) { aAvg[n]= pT->GetAvgValue(t,1.0-time_shift); incomp=true;}//incomplete start day
    else if((aAvg[n]==RAV_BLANK_DATA) && ((t-time_shift+1.0)>duration)) { aAvg[n]= pT->GetAvgValue(t-time_shift,duration); incomp=true;}//incomplete end day
    //
    t+=1.0;
  }
  this->AddTimeSeries(new CTimeSeries("TEMP_DAILY_MIN",DOESNT_EXIST,"",start_day-time_shift,start_yr,1.0,aMin,nVals,true),F_TEMP_DAILY_MIN);
  this->AddTimeSeries(new CTimeSeries("TEMP_DAILY_MAX",DOESNT_EXIST,"",start_day-time_shift,start_yr,1.0,aMax,nVals,true),F_TEMP_DAILY_MAX);
  this->AddTimeSeries(new CTimeSeries("TEMP_DAILY_AVE",DOESNT_EXIST,"",start_day-time_shift,start_yr,1.0,aAvg,nVals,true),F_TEMP_DAILY_AVE);
  delete [] aMin;
  delete [] aMax;
  delete [] aAvg;
  if (incomp){
    string warning;
	  warning="CGauge:GenerateMinMaxAveTempFromSubdaily: incomplete subdaily temperature data was used to estimate daily average temperature on first or last day of simulation at gauge "+_name;
    WriteWarning(warning,Options.noisy);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Generates Tave and subhourly time series from daily Tmin & Tmax time series
/// \note presumes existence of valid F_TEMP_DAILY_MIN and F_TEMP_DAILY_MAX time series
//
void CGauge::GenerateAveSubdailyTempFromMinMax(const optStruct &Options)
{
  CTimeSeries *pTmin,*pTmax,*pTdaily_ave;
  int          nVals;

  pTmin      =GetTimeSeries(F_TEMP_DAILY_MIN);
  pTmax      =GetTimeSeries(F_TEMP_DAILY_MAX);
  pTdaily_ave=GetTimeSeries(F_TEMP_DAILY_AVE);

  double start_day=Options.julian_start_day;
  int    start_yr =Options.julian_start_year;
  double duration =Options.duration;
  double timestep =Options.timestep;

  //below needed for correct mapping from time series to model time
  pTmin->Initialize(start_day,start_yr,duration,timestep,false,Options.calendar);
  pTmax->Initialize(start_day,start_yr,duration,timestep,false,Options.calendar);

  //Generate daily average values Tave=(Tmin+Tmax)/2 unless provided
  if(pTdaily_ave==NULL)
  {
    nVals=(int)ceil(duration);
    double *aAvg=new double[nVals];
    double t=0.0;//model time
    for(int n=0;n<nVals;n++)
    {
      aAvg[n]=0.5*(pTmin->GetValue(t+0.5)+pTmax->GetValue(t+0.5));
      t+=1.0;
    }
    pTdaily_ave=new CTimeSeries("TEMP_DAILY_AVE",DOESNT_EXIST,"",start_day,start_yr,1.0,aAvg,nVals,true);
    ExitGracefullyIf(pTdaily_ave==NULL,"GenerateAveSubdailyTempFromMinMax",OUT_OF_MEMORY);
    this->AddTimeSeries(pTdaily_ave,F_TEMP_DAILY_AVE);
    delete[] aAvg;
  }
  //Generate subdaily temperature values Tave
  if (Options.timestep<(1.0-TIME_CORRECTION))
  {
    nVals=(int)ceil(duration/Options.timestep);
    double *aT=NULL;
    aT=new double [nVals];
    ExitGracefullyIf(aT==NULL,"GenerateAveSubdailyTempFromMinMax",OUT_OF_MEMORY);
    double t=0.0;//model time
    for (int n=0;n<nVals;n++)
    {
      double Tmax=pTmax->GetValue(t+Options.timestep/2.0);
      double Tmin=pTmin->GetValue(t+Options.timestep/2.0);
      aT[n]=0.5*(Tmax+Tmin)+0.5*(Tmax-Tmin)*0.5*(DailyTempCorrection(t)+DailyTempCorrection(t+Options.timestep));
      t+=Options.timestep;
    }

    CTimeSeries *pNewTS=new CTimeSeries("TEMP_AVE",DOESNT_EXIST,"",start_day,start_yr,Options.timestep,aT,nVals,true);
    this->AddTimeSeries(pNewTS,F_TEMP_AVE);
    delete [] aT;
  }
  else{
    CTimeSeries *T=GetTimeSeries(F_TEMP_DAILY_AVE);
    this->AddTimeSeries(new CTimeSeries("TEMP_AVE",*T),F_TEMP_AVE); //just copy daily average values
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Generates Tmin, Tmax and subhourly time series from daily average temperature time series
/// \note presumes existence of valid F_TEMP_DAILY_AVE
/// \note necessarily naive - it is hard to downscale with little temp data
//
void CGauge::GenerateMinMaxSubdailyTempFromAve(const optStruct &Options)
{
  CTimeSeries *pT;
  pT=GetTimeSeries(F_TEMP_DAILY_AVE);

  double start_day=Options.julian_start_day; //floor(pT->GetStartDay());
  int    start_yr =Options.julian_start_year;//pT->GetStartYear();
  double duration =Options.duration;         //(interval*pTave->GetNumValues());
  double timestep =Options.timestep;

  //below needed for correct mapping from time series to model time
  pT->Initialize(start_day,start_yr,duration,timestep,false,Options.calendar);

  int nVals=(int)ceil(duration);
  double *aMin=new double [nVals];
  double *aMax=new double [nVals];
  double t=0.0;//model time
  for (int n=0;n<nVals;n++){
    aMin[n]=pT->GetValue(t+0.5)-4.0;//Options.temp_swing*0.5;
    aMax[n]=pT->GetValue(t+0.5)+4.0;//Options.temp_swing*0.5;
    t+=1.0;
  }
  this->AddTimeSeries(new CTimeSeries("TEMP_DAILY_MIN",DOESNT_EXIST,"",start_day,start_yr,1.0,aMin,nVals,true),F_TEMP_DAILY_MIN);
  this->AddTimeSeries(new CTimeSeries("TEMP_DAILY_MAX",DOESNT_EXIST,"",start_day,start_yr,1.0,aMax,nVals,true),F_TEMP_DAILY_MAX);
  delete [] aMin;
  delete [] aMax;

  GenerateAveSubdailyTempFromMinMax(Options);
}
