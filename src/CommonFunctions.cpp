//////////////////////////////////////////////////////////////////
///  Raven Library Source Code
///  Copyright (c) 2008-2023 the Raven Development Team
//////////////////////////////////////////////////////////////////

#include <time.h>
#include "RavenInclude.h"

//////////////////////////////////////////////////////////////////
/// \brief Returns a string describing the process corresponding to the enumerated process type passed
///
/// \param p [in] Type of proceses
/// \return String equivalent of the process_type identifier
//
string GetProcessName(process_type p)
{
  string name;
  switch(p)
  {
  case(NULL_PROCESS_TYPE):  {name="NULL";                     break;}
  case(FLUSH):              {name="Flush";                    break;}
  case(SPLIT):              {name="Split";                    break;}
  case(OVERFLOW_PROC):      {name="Overflow";                 break;}
  case(EXCHANGE_FLOW):      {name="Exchange Flow";            break;}

  case(PRECIPITATION):      {name="Precipitation";            break;}

  case(BASEFLOW):           {name="Baseflow";                 break;}
  case(INFILTRATION):       {name="Infiltration";             break;}
  case(INTERFLOW):          {name="Interflow";                break;}
    //case(REDISTRIBUTION):     {name="Soil Moisture Redistribution"; break;}
  case(PERCOLATION):        {name="Percolation";              break;}
  case(SOIL_EVAPORATION):   {name="Soil Evaporation";         break;}
  case(CAPILLARY_RISE):     {name="Capillary Rise";           break;}
  case(SOIL_BALANCE):       {name="Soil Balance";             break;}

  case(CANOPY_EVAPORATION): {name="Canopy Evaporation";       break;}
  case(CANOPY_SNOW_EVAPORATION):{name="Canopy Snow Sublimation";break;}
  case(CANOPY_DRIP):        {name="Canopy Drip";              break;}
  case(OPEN_WATER_EVAPORATION):{name="Open Water Evaporation";break;}
  case(LAKE_EVAPORATION):   {name="Lake Evaporation";         break;}
  case(LAKE_RELEASE):       {name="Lake Release";             break;}
  case(DEPRESSION_OVERFLOW):{name="Depression Overflow";      break;}
  case(SEEPAGE):            {name="Seepage from Depression";  break;}
  case(RECHARGE):           {name="Recharge";                 break;}
  case(DRAIN):              {name="Drain";                    break;}
  case(GWRECHARGE):         {name="Groundwater Recharge";     break;}

  case(SNOWMELT):           {name="Snow Melt";                break;}
  case(SNOWSQUEEZE):        {name="Liquid snow release";      break;}
  case(REFREEZE):           {name="Snow Refreeze";            break;}
  case(SUBLIMATION):        {name="Sublimation";              break;}
  case(SNOW_BALANCE):       {name="Snow Balance";             break;}
  case(GLACIER_MELT):       {name="Glacier Melt";             break;}
  case(GLACIER_RELEASE):    {name="Glacier Release";          break;}
  case(GLACIER_INFIL):      {name="Glacier Infiltration";     break;}
  case(SNOW_ALBEDO_EVOLVE): {name="Snow Albedo Evolution";    break;}
  case(BLOWING_SNOW):       {name="Blowing Snow";             break;}
  case(LAKE_FREEZING):      {name="Lake Freezing";            break;}
  case(GROUND_FREEZING):    {name="Ground Freezing";          break;}

  case(CROP_HEAT_UNIT_EVOLVE):{name="Crop Heat Unit Evolution";break;}
  case(ABSTRACTION):        {name="Abstraction";              break;}
  case(SNOWTEMP_EVOLVE):    {name="Snow Temp. Evolution";     break;}
  case(CONVOLVE):           {name="Convolution";              break;}

  case(PROCESS_GROUP):      {name="Process Group";            break;}

  case(ADVECTION):          {name="Advection";                break;}
  case(LAT_ADVECTION):      {name="Lateral Advection";        break;}
  case(LAT_FLUSH):          {name="Lateral Flush";            break;}
  case(LAT_EQUIL):          {name="Lateral Equilibrate";      break;}
  case(MASS_LOADING):       {name="Mass Loading";             break;}
  case(DECAY):              {name="Decay";                    break;}
  case(TRANSFORMATION):     {name="Transformation";           break;}

  case(HEATCONDUCTION):     {name="Heat Conduction";          break;}
    //..
  default:              {
    name="Unknown Hydrological Process";
    ExitGracefully("GetProcessName: Unknown Hydrological Process",RUNTIME_ERR);//STRICT
    break;
  }
  }
  return name;
}

///////////////////////////////////////////////////////////////////////////
/// \brief Determines if parameter value is to be autocalculated
/// \remark Used in override of parameter/property values with defaults (in AutoCalculate routines)
///
/// \param &val [out] Parameter value to be determined
/// \param set_val [in] Parameter value specified, or USE_TEMPLATE_VALUE if the template is to be applied
/// \param template_val [in] Parameter value if the default template is to be applied
/// \return Boolean indicating if the parameter value is to be autocalculated
//
bool SetCalculableValue(double &val, double set_val, double template_val)
{
  //  preferred approach if master parameter list is complete
 if ((template_val==NOT_NEEDED) || (template_val==NOT_NEEDED_AUTO)){
    val=template_val; //even if value specified, overriden, because it is not needed
    return false;     //autocalculation not needed
  }
  //approach needed to account for the fact that master parameter list may be incomplete
  if (template_val==NOT_NEEDED_AUTO){
    template_val=AUTO_COMPUTE;
  }

  val =set_val;
  if (val==USE_TEMPLATE_VALUE){val=template_val;} //override with default
  return (val==AUTO_COMPUTE); //returns true if the parameter value is to be autocalculated
}

///////////////////////////////////////////////////////////////////////////
/// \brief Sets the value of a non-autocalculable model parameter
/// \remark Used in override of parameter/property values with defaults (in AutoCalculate routines)
///
/// \param &val [out] Parameter value to be specified
/// \param set_val [in] Parameter value specified, or USE_TEMPLATE_VALUE if the template is to be applied
/// \param template_val [in] Parameter value if the default template is to be applied
/// \param needed [in] Indicates if parameter value is needed
/// \return Boolean indicating if the parameter value is not specified but required by model
/// \note this is only used for non-autogeneratable parameters
/// \todo [QA/QC] this routine should take in parameter name and class for better warnings to be produced
//
bool SetSpecifiedValue(double &val, double set_val, double template_val, bool needed, string pname)
{
  val =set_val;
  if (val==USE_TEMPLATE_VALUE){val=template_val;} //override with template
  if ((val==USE_TEMPLATE_VALUE) && (template_val==AUTO_COMPUTE))
  {
    string errmsg;
    errmsg = "SetSpecifiedValue::Cannot specify _AUTO as the template value for any parameter that cannot be autocomputed. A required parameter "+pname+"is missing from the .rvp file.";
    ExitGracefully(errmsg.c_str(),BAD_DATA_WARN);
  }
  if (val==AUTO_COMPUTE)
  {
    string errmsg;
    errmsg = "SetSpecifiedValue::Cannot specify _AUTO as the default value for any parameter that cannot be autocomputed. A required parameter " + pname + "is missing from the .rvp file";
    ExitGracefully(errmsg.c_str(),BAD_DATA_WARN);
  }
  return ((val==NOT_SPECIFIED) && (needed)); //returns true (bad) if the parameter value is not specified, but needed
}

///////////////////////////////////////////////////////////////////////////
/// \brief Determines default parameter value
/// \remark Default template values depend upon whether the parameter is autocalculable
///
/// \param is_template [in] true if the default value being set is for the template class
/// \param is_computable [in] true if the value is autocalculable
/// \return default parameter value for the parameter
//
double DefaultParameterValue(bool is_template, bool is_computable)
{
  if (!is_template){return USE_TEMPLATE_VALUE;} //this is the default value for all parameters
  else //default template values depend upon whether the parameter is autocalculable
  {
    //if (is_computable){return AUTO_COMPUTE;}
    //else              {return NOT_SPECIFIED;}
    if (is_computable){return NOT_NEEDED_AUTO;}
    else              {return NOT_NEEDED;}
  }
}

///////////////////////////////////////////////////////////////////////////
/// \brief Dynamically append pointer to array
/// \details Dynamically adds additional pointer (*xptr) onto array of pointers (**pArr) of
/// initial size indicated by parameter. Increments size by 1
///
/// \param **&pArr [out] Array of pointers to which *xptr will be added
/// \param *xptr [in] Pointer to be added to array
/// \param &size [in & out] Integer size of array
/// \return Boolean indicating success of method
//
bool DynArrayAppend(void **& pArr, void *xptr,int &size)
{
  void **tmp=NULL;
  if (xptr==NULL){return false;}
  if ((pArr==NULL) && (size>0)) {return false;}
  size=size+1;                                                //increment size
  tmp=new void *[size+1];                                     //allocate memory
  if (tmp==NULL){ExitGracefully("DynArrayAppend::Out of memory",OUT_OF_MEMORY);}
  for (int i=0; i<(size-1); i++){                             //copy array
#ifdef _STRICTCHECK_
    if (pArr[i]==NULL){ExitGracefully("DynArrayAppend::Bad existing array",RUNTIME_ERR);}
#endif
    tmp[i]=pArr[i];
  }
  tmp[size-1]=xptr;                                           //add new pointer
  if (size>1){delete [] pArr; pArr=NULL;}                     //delete old array of pointers
  pArr=tmp;                                                   //redirect pointer
  return true;
}

/**************************************************************************
      Threshold Smoothing functions
---------------------------------------------------------------------------
From Kavetski & Kuczera, 2007
**************************************************************************/

///////////////////////////////////////////////////////////////////////////
/// \brief Enforces positivity of input value \cite Kavetski2007WRR
/// \todo [funct] Work needed here if smoothing is to be used - should probably DELETE ENTIRELY
///
/// \param &val [in] input value on which positivity is enforced
/// \return Returns the value itself if it is greater than 0.0, otherwise returns 0.0.
//
double threshPositive(const double &val)
{
  return max(val,0.0);
  //return threshMax(val,0.0,0.01);
}

////////////////////////////////////////////////////////////////////////////
/// \brief Implementation of Chen-Harker-Kanzow-Smale max function \cite Chen2000JRSOJ
/// \remark Result is normalized: deviation from actual max does not exceed smooth_coeff (@ v1-v2=0)
/// \param &v1 [in] First input value to function
/// \param &v2 [in] Second input value to function
/// \param &smooth_coeff [in] Smoothing coefficient parameter
/// \return Double output of function
double threshMax(const double &v1, const double &v2, const double &smooth_coeff)
{
  double x=(v2-v1);
  if      (smooth_coeff==0.0)         {return max(v1,v2);}
  else if (fabs(x)/smooth_coeff>100.0){return max(v1,v2 );}//avoids numerical errors
  else{
    return v1 + 0.5*(x + sqrt(x*x + 4.0*smooth_coeff*smooth_coeff));
  }
}

////////////////////////////////////////////////////////////////////////////
/// \brief Implementation of Chen-Harker-Kanzow-Smale min function \cite Chen2000JRSOJ
/// \remark Result is normalized: deviation from actual min does not exceed smooth_coeff (@ v1-v2=0)
///
/// \param &v1 [in] First input value to function
/// \param &v2 [in] Second input value to function
/// \param &smooth_coeff [in] Smoothing coefficient parameter
/// \return Double output of function
//
double threshMin(const double &v1, const double &v2, const double &smooth_coeff)
{
  double x=(v2-v1);
  if      (smooth_coeff==0.0)         {return min(v1,v2);}
  else if (fabs(x)/smooth_coeff>100.0){return min(v1,v2);}//avoids numerical errors
  else{
    return v2 - 0.5 * (x + sqrt(x*x + 4.0*smooth_coeff*smooth_coeff));
  }
}

////////////////////////////////////////////////////////////////////////////
/// \brief Rounds time to nearest minute to prevent roundoff error in NetCDF reporting of time
/// \param &t [in] model time, in hours
/// \return time, in hours, rounded to nearest minute
//
double RoundToNearestMinute(const double &t)
{
  const double MIN_PER_HOUR=60;
  return floor(t+TIME_CORRECTION)+(rvn_round((t-floor(t+TIME_CORRECTION))*MIN_PER_HOUR))/MIN_PER_HOUR;
}

////////////////////////////////////////////////////////////////////////////
/// \brief Returns boolean value indicating if year passed is a leap year
/// \note Only valid until ~4000 AD
///
/// \param year     [in] Integer year
/// \param calendar [in] enum int of calendar used
/// \return Boolean value indicating if year passed is a leap year

bool IsLeapYear(const int year, const int calendar)
{
    bool leap = false;

    // handle default CALENDAR_PROLEPTIC_GREGORIAN first for code efficiency
    if (calendar == CALENDAR_PROLEPTIC_GREGORIAN){
      leap = (((year % 4 == 0) && (year % 100 != 0)) || (year % 400 == 0)); //valid until ~4000 AD:)
      return leap;
    }

    if(calendar==CALENDAR_365_DAY){return false;}
    if(calendar==CALENDAR_366_DAY){return true; }

    // other calendars
    if ((calendar == CALENDAR_JULIAN ||
         calendar == CALENDAR_GREGORIAN ) &&
        (year % 4 == 0)) {
      leap = true;
      if ((calendar == CALENDAR_GREGORIAN) &&
	         (year % 100 == 0) && (year % 400 != 0) && (year > 1583)) {
	      leap = false;
      }
    }
    return leap;
}

///////////////////////////////////////////////////////////////////////////
/// \brief Fills time structure tt
/// \details Converts Julian decimal date to string and returns day of month, month,
/// and year in time_struct. If dec_date >365/366, then year is incremented. Accounts for leap years
///
/// \param &model_time [in]  Time elapsed since start of simulation
/// \param start_date  [in]  double simulation start date (Julian date)
/// \param start_year  [in]  Integer simulation start year
/// \param calendar    [in]  enum int of calendar used
/// \param &tt         [out] Time structure to house date information
//
void JulianConvert(double model_time, const double start_date, const int start_year, const int calendar, time_struct &tt)
{
  int leap(0);
  string mon;
  double sum,days,ddate;
  double dday;
  int    dmonth,dyear;

  //handles daily roundoff error, (e.g., t=4.999873->t=5.0)
  if( (model_time-floor(model_time)) > (1-TIME_CORRECTION))
  {
    model_time = floor(model_time+TIME_CORRECTION);
  }
  //handles hourly roundoff error (e.g., t=5.24999873->t=5.25)
  if(model_time*HR_PER_DAY-floor(model_time*HR_PER_DAY) > (1.0-TIME_CORRECTION)*HR_PER_DAY) {
    model_time = floor(HR_PER_DAY*(model_time+TIME_CORRECTION))/HR_PER_DAY;
  }

  double dec_date=start_date+model_time; //decimal date calculated from start_date,start year

  dyear=start_year;
  ddate=dec_date;

  if (IsLeapYear(dyear,calendar)){leap=1;}
  while (ddate>=(365+leap)) //correct for years
  {
    ddate-=(365.0+leap);
    dyear++;
    leap=0;if (IsLeapYear(dyear,calendar)){leap=1;}
  }
  //ddate is now decimal julian date from Jan 1 0:00:00 of current dyear

  dmonth=1; days=31;sum=31-TIME_CORRECTION; mon="Jan";
  if      (ddate>=sum){dmonth+=1;days=28+leap;sum+=days;mon="Feb";}//Feb
  if      (ddate>=sum){dmonth+=1;days=31;     sum+=days;mon="Mar";}//Mar
  if      (ddate>=sum){dmonth+=1;days=30;     sum+=days;mon="Apr";}//Apr
  if      (ddate>=sum){dmonth+=1;days=31;     sum+=days;mon="May";}//May
  if      (ddate>=sum){dmonth+=1;days=30;     sum+=days;mon="Jun";}//Jun
  if      (ddate>=sum){dmonth+=1;days=31;     sum+=days;mon="Jul";}//Jul
  if      (ddate>=sum){dmonth+=1;days=31;     sum+=days;mon="Aug";}//Aug
  if      (ddate>=sum){dmonth+=1;days=30;     sum+=days;mon="Sep";}//Sep
  if      (ddate>=sum){dmonth+=1;days=31;     sum+=days;mon="Oct";}//Oct
  if      (ddate>=sum){dmonth+=1;days=30;     sum+=days;mon="Nov";}//Nov
  if      (ddate>=sum){dmonth+=1;days=31;     sum+=days;mon="Dec";}//Dec

  dday=ddate-sum+days; //decimal days since 0:00 on first of month


  tt.model_time=model_time;
  tt.julian_day=ddate;
  tt.day_of_month=(int)(ceil(dday+REAL_SMALL)); //real_small to handle dday=1.0
  if (tt.day_of_month==0){tt.day_of_month=1;}
  tt.month =dmonth;
  tt.year=dyear;
  //tt.nn = static_cast<int>(rvn_floor((tt.model_time+TIME_CORRECTION)/tstep)));

  tt.day_changed = false;
  if((model_time <= PRETTY_SMALL) || (tt.julian_day-floor(tt.julian_day+TIME_CORRECTION)<0.001)) { tt.day_changed = true; }

  static char out[50];
  sprintf(out,"%4.4d-%2.2i-%2.2d",dyear,tt.month,tt.day_of_month); //2006-02-28 (ISO Standard)

  tt.date_string=string(out);
  tt.leap_yr=IsLeapYear(tt.year,calendar);

  //tt.nn=(int)((tt.model_time+TIME_CORRECTION)/Options.timestep);//current timestep index - likely needs Options.timestep to be in global member

}

////////////////////////////////////////////////////////////////////////////
/// \brief Returns time-of-day string for decimal date (e.g., day=124.3-->"07:12",day=1.5  -->"12:00")
///
/// \param dec_date [in] Decimal date
/// \return String hours of day in 00:00:00.00 format
//
string DecDaysToHours(const double dec_date, const bool truncate)
{
  double hours=(dec_date-floor(dec_date))*24.0;
  double mind =(hours   -floor(hours   ))*60.0;
  double sec  =(mind    -floor(mind    ))*60.0;

  int hr=(int)(floor(hours));
  int min=(int)(floor(mind));

  if (sec>=59.995){min++;sec=0.0;} //to account for case where time is rounded up to 60 sec
  if (min==60)    {hr++; min=0;}
  if (hr==24)     {hr=0;}

  static char out[12];
  if (truncate){sprintf(out,"%2.2d:%2.2d:%2.2d",hr,min,(int)(sec));}
  else         {sprintf(out,"%2.2d:%2.2d:%05.2f",hr,min,sec);}
  return string(out);
}

////////////////////////////////////////////////////////////////////////////
/// \brief returns time struct corresponding to string in the following format
/// \param sDate    [in] date string in ISO standard format yyyy-mm-dd or yyyy/mm/dd
/// \param calendar [in] Enum int of calendar used
/// \param sTime    [in] time string in ISO standard format hh:mm:ss.00
/// \return Time structure equivalent of passed date and time
//
time_struct DateStringToTimeStruct(const string sDate, string sTime, const int calendar)
{

  static time_struct tt;
  if (sDate.length()!=(size_t)(10)){
    string errString = "DateStringToTimeStruct: Invalid date format used: "+sDate;
    ExitGracefully(errString.c_str(),BAD_DATA);
  }
  if (sTime.length()<(size_t)(7))  {
    string errString = "DateStringToTimeStruct: Invalid time format used (hourstamp): "+sTime;
    ExitGracefully(errString.c_str(),BAD_DATA);
  }

  tt.date_string=sDate;
  tt.year        =s_to_i(sDate.substr(0,4).c_str());
  tt.month       =s_to_i(sDate.substr(5,2).c_str());
  if (tt.month>12)  {
    string errString = "DateStringToTimeStruct: Invalid time format used (month>12): "+sDate;
    ExitGracefully(errString.c_str(),BAD_DATA);
  }
  tt.day_of_month=s_to_i(sDate.substr(8,2).c_str());
  tt.model_time  =0.0;//unspecified
  tt.leap_yr     =IsLeapYear(tt.year,calendar);
  tt.julian_day  =tt.day_of_month-1;

  if (tt.month>= 2){tt.julian_day+=31;}
  if (tt.month>= 3){tt.julian_day+=28;}
  if (tt.month>= 4){tt.julian_day+=31;}
  if (tt.month>= 5){tt.julian_day+=30;}
  if (tt.month>= 6){tt.julian_day+=31;}
  if (tt.month>= 7){tt.julian_day+=30;}
  if (tt.month>= 8){tt.julian_day+=31;}
  if (tt.month>= 9){tt.julian_day+=31;}
  if (tt.month>=10){tt.julian_day+=30;}
  if (tt.month>=11){tt.julian_day+=31;}
  if (tt.month==12){tt.julian_day+=30;}
  if ((tt.leap_yr  ) && (tt.month> 2)){tt.julian_day+= 1;}

  int leap=0;
  if (tt.leap_yr){leap=1;}

  if (tt.day_of_month > (DAYS_PER_MONTH[tt.month - 1]+leap)) {
    ExitGracefully("DateStringToTimeStruct: Invalid time format used - exceeded max day of month",BAD_DATA);
  }
  if (tt.day_of_month <=0 ) {
    ExitGracefully("DateStringToTimeStruct: Invalid time format used - negative or zero day of month",BAD_DATA);
  }

  int hr, min;
  double sec;

  if (sTime.substr(1,1)==":"){sTime="0"+sTime;} //for h:mm:ss.00 format to hh:mm:ss.00

  ExitGracefullyIf((sTime.substr(2,1)!=":"),"DateStringToTimeStruct: Invalid time format used",BAD_DATA);
  ExitGracefullyIf((sTime.substr(5,1)!=":"),"DateStringToTimeStruct: Invalid time format used",BAD_DATA);

  hr        =s_to_i(sTime.substr(0,2).c_str());
  min       =s_to_i(sTime.substr(3,2).c_str());
  sec       =s_to_d(sTime.substr(6,6).c_str());
  tt.julian_day+=(double)(hr )/HR_PER_DAY;
  tt.julian_day+=(double)(min)/MIN_PER_DAY;
  tt.julian_day+=(double)(sec)/SEC_PER_DAY;

  //Below reprocesses date string (optional)
  JulianConvert(0.0, tt.julian_day, tt.year, calendar, tt);

  return tt;
}
////////////////////////////////////////////////////////////////////////////
/// \brief returns true if julian date is between two julian days (days inclusive)
/// \param julian_day   [in] julian date from 0.0 to 365.0
/// \param julian_start [in] integer start day of date range (0=Jan 1, 364=Dec 31 in non-leap)
/// \param julian_end [in] integer end day of date range (0=Jan 1, 364=Dec 31 in non-leap)
//
bool        IsInDateRange(const double &julian_day,const int &julian_start,const int &julian_end)
{
  if(julian_start<julian_end) {
    return ((julian_day>=julian_start) && (julian_day<=julian_end));
  }
  else {
    return ((julian_day>=julian_start) || (julian_day<=julian_end)); //wraps around Dec 31-Jan 1
  }
}
////////////////////////////////////////////////////////////////////////////
/// \brief returns zero-indexed julian date (0..364(5)) given string in format Mmm-dd, mm-dd, or Julian date
/// \param date_str [in] date string ('Apr-23', '04-23', or '113')
/// \param calendar   [in] enumerated calendar type
/// \return zero indexed julian date
//
int GetJulianDayFromMonthYear   (const string &date_str, const int calendar)
{
  if((date_str.length()==3) && (date_str.substr(1,1)=="-"))//support for m-d format (4-2)
  {
    string fulldate_str="1999-0"+date_str.substr(0,1)+"-0"+date_str.substr(2,1); //no leap year
    time_struct tt = DateStringToTimeStruct(fulldate_str, "00:00:00", calendar);
    return (int)(tt.julian_day);
  }
  else if((date_str.length()>=2) && (date_str.substr(1,1)=="-"))//support for m-dd format (4-23)
  {
    string fulldate_str="1999-0"+date_str; //no leap year
    time_struct tt = DateStringToTimeStruct(fulldate_str, "00:00:00", calendar);
    return (int)(tt.julian_day);
  }
  else if((date_str.length()==4) && (date_str.substr(2,1)=="-"))//support for mm-d format (10-2)
  {
    string fulldate_str="1999-"+date_str.substr(0, 3) + "0" + date_str.substr(3, 1); //no leap year
    time_struct tt = DateStringToTimeStruct(fulldate_str, "00:00:00", calendar);
    return (int)(tt.julian_day);
  }
  else if((date_str.length()==5) && (date_str.substr(2,1)=="-"))//support for mm-dd format (04-23)
  {
    string fulldate_str="1999-"+date_str; //no leap year
    time_struct tt = DateStringToTimeStruct(fulldate_str, "00:00:00", calendar);
    return (int)(tt.julian_day);
  }
  else if ((date_str.length()>=3) && (date_str.substr(3,1)=="-"))//support for Mmm-dd & Mmm-d format (Apr-23, Apr-02, or Apr-2)
  {
    string months[12] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
    string fulldate_str="";
    for (int mon=0;mon<12;mon++){
      if (date_str.substr(0,3)==months[mon]){
        if (mon<9){fulldate_str="0"; }
        if (date_str.length()<=5){ //Apr-1
          fulldate_str=fulldate_str+to_string(mon+1)+"-0"+date_str.substr(4, 1);
        }
        else{ //Apr-01 or Apr-23
          fulldate_str=fulldate_str+to_string(mon+1)+"-"+date_str.substr(4, 2);
        }
      }
    }
    if (fulldate_str == "") {
      ExitGracefully("GetJulianDayFromMonthYear: invalild month indicated in Mmm-dd date format",BAD_DATA);
    }
    fulldate_str="1999-"+fulldate_str;
    time_struct tt = DateStringToTimeStruct(fulldate_str, "00:00:00", calendar);
    return (int)(tt.julian_day);
  }
  else{ //julian date format
    int val=s_to_i(date_str.c_str());
    if ((val==0) || (val>365)){
      ExitGracefully("GetJulianDayFromMonthYear: invalild format of julian date. Should be in 'Mmm-dd', 'mm-dd', or supplied as julian date 1..365",BAD_DATA);
    }
    return s_to_i(date_str.c_str())-1; //convert from integer julian days to Raven 0-indexed julian days
  }
}
////////////////////////////////////////////////////////////////////////////
/// \brief returns time struct corresponding to string in the following format
/// \param unit_t_str [in] full time string from NetCDF file (e.g., 'days since YYYY-MM-dd 00:00:00+0000')
/// \param timestr    [in] first word of string (e.g., 'days')
/// \param calendar   [in] enumerated calendar type
/// \param timezone   [out] time shift from GMT, in days
/// \return Raven Time structure equivalent of passed date and time, time shift if applicable
//
time_struct TimeStructFromNetCDFString(const string unit_t_str,const string timestr,const int calendar,double &timezone)
{
  string dash,colon,tmp;
  tmp=unit_t_str;
  timezone=0.0;
  int start=(int)strlen(timestr.c_str());
  start+=7; //first char of year YYYY (7=length(' since '))
  // ---------------------------
  // check if format is hours since YYYY-MM-DD HH:MM:SS, fill with leading zeros if necessary
  // blank between date and time (position 10) can be either ' ' or 'T'
  // Y  Y  Y  Y  -  M  M  -  d  d  _  0  0  :  0  0  :  0  0  .  0     +  0  0  0  0
  // 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
  // ---------------------------
  dash = tmp.substr(start+4,1);  // first dash in date

  if(!strstr(dash.c_str(),"-"))
  {
    printf("time unit string: %s\n",tmp.c_str());
    ExitGracefully("CommonFunctions:TimeStructFromNetCDFString: time unit string has weird format!",BAD_DATA);
  }
  if(!strstr(tmp.substr(start+7 ,1).c_str(),"-")) { tmp.insert(start+5,"0"); } // second dash in date - fixes YYYY-M-dd

  bool date_only=false;
  if ((int)(strlen(tmp.c_str()))<=(start+10+2)){date_only=true;} //00:00:00 +0000 not included in string

  if(!date_only) {
    if(!strstr(tmp.substr(start+10,1).c_str()," ") && !strstr(tmp.substr(start+10,1).c_str(),"T")) { tmp.insert(start+8,"0"); } // second dash in date - fixes YYYY-MM-d
    if( strstr(tmp.substr(start+10,1).c_str(),"T")) {tmp.replace(start+10,1," "); } // replacing 'T' with ' ' in case date looks like YYYY-MM-DDTHH:MM:SS
    if(!strstr(tmp.substr(start+13,1).c_str(),":")) { tmp.insert(start+11,"0"); } // first colon in time  - fixes 1:00:00
    if(!strstr(tmp.substr(start+16,1).c_str(),":")) { tmp.insert(start+14,"0"); } // second colon in time - fixes 11:0:00 (?)
  }
  else {
    if(strlen(tmp.c_str())==(size_t)(start+9)) { tmp.insert(start+8,"0"); } // second dash in date - fixes YYYY-MM-d
  }
  string sTime,sDate;
  sDate = tmp.substr(start,10); //YYYY-MM-DD
  if(!date_only) { sTime = tmp.substr(start+11,8); }  //HH:MM:SS
  else           { sTime = "00:00:00";             }  // assumes start of day if no timestamp given
  //cout<<"sTime"<<sTime<<" sDate"<< sDate<<" "<<unit_t_str<<" "<<date_only<<" "<<tmp<<endl;

  timezone=0;
  if((strlen(tmp.c_str())-start)==(size_t)(26)) {
    if(!strcmp(tmp.substr(start+22,1).c_str(),"+")) {
      timezone=(double)(s_to_i(tmp.substr(start+23,4).c_str()))/HR_PER_DAY/100; //time zone, in days
    }
    else if(!strcmp(tmp.substr(start+22,1).c_str(),"-")) {
      timezone=-(double)(s_to_i(tmp.substr(start+23,4).c_str()))/HR_PER_DAY/100; //time zone, in days
    }
  }
  return DateStringToTimeStruct(sDate,sTime,calendar);
}
////////////////////////////////////////////////////////////////////////////
/// \brief returns time zone string in format " +0700" or " -1200"
/// \param tz [in] time zone as integer - hours from GMT
/// \return time zone as string in format " +0600" or "" if tz is zero
//
string TimeZoneToString(const int tz) {
  string ends="";
  if(tz<0) {
    if(tz<-9) { ends=" -" +to_string(-tz)+"00"; }
    else      { ends=" -0"+to_string(-tz)+"00"; }
  }
  else if(tz>0) {
    if(tz>9) { ends=" +" +to_string(tz)+"00"; }
    else     { ends=" +0"+to_string(tz)+"00"; }
  }
  return ends;
}

////////////////////////////////////////////////////////////////////////////
/// \brief returns time struct corresponding to string in the following format
/// \param unit_t_str [in] full time string from NetCDF file (e.g., '[days/minutes/hours] since YYYY-MM-dd 00:00:00{+0000}' or '[days/minutes/hours] since YYYY-MM-dd' )
/// \return true if string is valid
//
bool IsValidNetCDFTimeString(const string time_string)
{
  int att_len=(int)strlen(time_string.c_str());
  bool isvalid = true;
  if (att_len<15) {return false;}

  size_t pos=time_string.find("since",0);
  if(pos==string::npos) { return false; } //no "since" in string

  pos+=6;
  string date_string=time_string.substr(pos,10);

  if (!strstr(date_string.substr(4,1).c_str(),"-")){isvalid=false;}//properly located dashes in date string
  if (!strstr(date_string.substr(7,1).c_str(),"-")){isvalid=false;}

  if(time_string.length()<(pos+19)) { return isvalid; } //no time stamp

  string hr_string  =time_string.substr(pos+11,8);
  //cout<<"TIME STRING: "<<time_string<<" "<<pos<<" "<<date_string<<" "<<hr_string<<endl;

  if(!strstr(hr_string.substr(2,1).c_str(),":")) { isvalid=false; }//properly located dashes in date string
  if(!strstr(hr_string.substr(5,1).c_str(),":")) { isvalid=false; }

  return isvalid;
}
///////////////////////////////////////////////////////////////////
/// \brief calculates time difference, in days, between two specified dates
/// \details positive if day 2 is after day 1
///
/// \param jul_day1 [in] Julian day of date 1 (measured from Jan 1 of year @ 00:00:00)
/// \param year1    [in] year of date 1
/// \param jul_day2 [in] Julian day of date 2 (measured from Jan 1 of year @ 00:00:00)
/// \param year1    [in] year of date 2
/// \param calendar [in] enum int of calendar used

double TimeDifference(const double jul_day1,const int year1,const double jul_day2,const int year2, const int calendar)
{
  int leap,yr;
  double diff= jul_day2 - jul_day1;
  yr=year2-1;
  while (yr >= year1)
  {
    leap=0; if (IsLeapYear(yr,calendar)){ leap = 1; }
    diff += (365+leap);
    yr--;
  }
  yr=year2;
  while (yr<year1)
  {
    leap=0; if (IsLeapYear(yr,calendar)){ leap = 1; }
    diff -= (365+leap);
    yr++;
  }
  return diff;
}

///////////////////////////////////////////////////////////////////
/// \brief adds specified number of days to julian date and returns resultant julian date
///
/// \param jul_day1    [in]  Julian day of date 1 (measured from Jan 1 of year @ 00:00:00)
/// \param year1       [in]  year of date 1
/// \param daysadded   [in]  positive or negative number of days (can be fractional days) added to date 1
/// \param &Options    [in] Global model options information
/// \param jul_day_out [out] Julian day of output date (measured from Jan 1 of year @ 00:00:00)
/// \param year_out    [out] year of output date
//
void AddTime(const double jul_day1,const int year1,const double &daysadded,const int calendar, double &jul_day_out,int &year_out)
{
  int    yr;
  double leap;
  double daysleft;

  yr=year1;
  jul_day_out=jul_day1;

  if(daysadded>=0)
  {
    daysleft=daysadded;
    do {
      leap=0; if(IsLeapYear(yr,calendar)) { leap=1; }
      if((jul_day_out+daysleft)<(365.0+leap)) {
        jul_day_out+=daysleft;
        year_out=yr;
        break;
      }
      else {
        yr++;
        daysleft-=(365.0+leap-jul_day_out);
        jul_day_out=0.0;
      }
      ExitGracefullyIf(daysleft<0.0,"Invalid input to AddTime routine (negative julian day?)",RUNTIME_ERR);
    } while(true);
  }
  else
  { //if daysadded<0
    daysleft=-daysadded;
    do {
      if((jul_day_out-daysleft)>=0.0) { //99% of cases
        jul_day_out-=daysleft;
        year_out=yr;
        return;
      }
      else {
        yr--;
        leap=0; if(IsLeapYear(yr,calendar)) { leap=1; }
        daysleft-=jul_day_out;
        if(daysleft<(365+leap)){ jul_day_out=(365+leap)-daysleft;year_out=yr;break; }
        else                   { jul_day_out=0.0;daysleft-=(365+leap); }//skip whole year
      }
      ExitGracefullyIf(daysleft<0.0,"Invalid input to AddTime routine (negative julian day?)",RUNTIME_ERR);
    } while(true);
  }

  // if calendar is STANDARD or GREGORIAN and the original time is before 4 October 1582
  // while the final time is after, one has to add additional 10 days
  // because in this calendar the day following 4 October 1582 is 15 October 1582 (there are 10 days missing)
  // --> THIS is why people usually use the Proleptic Gregorian calendar :)
  if ((calendar == CALENDAR_GREGORIAN) &&
      ((year1 == 1582 && jul_day1 <= 277) || (year1 < 1582)) &&
      ((year_out > 1582) || ((year_out == 1582) && (jul_day_out >= 278)))) {

    double tmp_day;
    int    tmp_yr;

    tmp_yr  = year_out;
    tmp_day = jul_day_out;

    AddTime(tmp_day,tmp_yr,10.0,calendar,jul_day_out,year_out);
    return;

  }
  return;

}

///////////////////////////////////////////////////////////////////
/// \brief Parse chars of calendar and return calendar integer
///
/// \param cal_chars        [in]  String conatining calendar name, e.g., "PROLEPTIC_GREGORIAN"
/// \param StringToCalendar [out] enum integer representing calendar
//
int StringToCalendar(string cal_chars)
{
  string str=StringToUppercase(cal_chars);
  if (strcmp("STANDARD", str.c_str()) == 0) {
    return CALENDAR_GREGORIAN;
  }
  else if (strcmp("GREGORIAN", str.c_str()) == 0) {
    return CALENDAR_GREGORIAN;
  }
  else if (strcmp("PROLEPTIC_GREGORIAN", str.c_str()) == 0) {
    return CALENDAR_PROLEPTIC_GREGORIAN;
  }
  else if ((strcmp("NOLEAP", str.c_str()) == 0) || (strcmp("NO_LEAP", str.c_str()) == 0)) {
    return CALENDAR_365_DAY;
  }
  else if (strcmp("365_DAY", str.c_str()) == 0) {
    return CALENDAR_365_DAY;
  }
  else if (strcmp("360_DAY", str.c_str()) == 0) {
    ExitGracefully("CommonFunctions: StringToCalendar: Raven does not support 360_DAY calendars!", BAD_DATA);
    return CALENDAR_360_DAY;
  }
  else if (strcmp("JULIAN", str.c_str()) == 0) {
    return CALENDAR_JULIAN;
  }
  else if (strcmp("ALL_LEAP", str.c_str()) == 0) {
    return CALENDAR_366_DAY;
  }
  else if (strcmp("366_DAY", str.c_str()) == 0) {
    return CALENDAR_366_DAY;
  }
  else {
    printf("Calendar used: %s", str.c_str());
    ExitGracefully("CommonFunctions: StringToCalendar: Unknown calendar specified!", BAD_DATA);
  }
  return -1;  // just to avoid compiler warning of void function return
}


////////////////////////////////////////////////////// /////////////////////
/// \brief Round the timestep to the nearest fractional day
/// \return improved timestep
double    FixTimestep(double tstep)
{
  double tmp = rvn_round(1.0/tstep);
  ExitGracefullyIf(fabs(tstep*tmp-1.0)>0.1,
                   "CommonFunctions::FixTimestep: timesteps and time intervals must evenly divide into one day",BAD_DATA);
  return 1.0/tmp;
}
////////////////////////////////////////////////////// /////////////////////
/// \brief True if string is proper iso date (e.g., yyyy-mm-dd or yyyy/mm/dd)
/// \return true if valid date string
//
bool IsValidDateString(const string sDate)
{
  return ((sDate.length()==10) &&
          ((sDate.substr(4,1)=="/") || (sDate.substr(4,1)=="-")) &&
          ((sDate.substr(7,1)=="/") || (sDate.substr(7,1)=="-")));
}
////////////////////////////////////////////////////// /////////////////////
/// \brief Get the current system date/time
/// \return "now" as an ISO formatted string
string GetCurrentMachineTime(void)
{
  // Get the current wall clock time
  time_t now;
  time(&now);
  struct tm *curTime = localtime(&now);

  // generate the ISO string
  char s[20];
  sprintf(s,"%4i-%02i-%02i %02i:%02i:%02i",
          curTime->tm_year+1900, curTime->tm_mon+1, curTime->tm_mday,
          curTime->tm_hour, curTime->tm_min, curTime->tm_sec);

  return string(s);
}

///////////////////////////////////////////////////////////////////////////
/// \brief Interpolates monthly data stored in array aVal during year based on specified time
/// \remark Model time, t[days] is specified
///
/// \param aVal     [in] Array of doubles representing monthly data
/// \param &tt      [in] Time structure which specifies interpolation
/// \param &method [in] Interpolation method
/// \return Interpolated value at time denoted by &tt
double InterpolateMo(const double       aVal[12],
                     const time_struct &tt,
                     const monthly_interp &method,
                     const optStruct &Options)
{
  double wt;
  int leap(0),mo,nextmo;
  int day, month, year;

  day    =tt.day_of_month;
  month  =tt.month;
  year   =tt.year;

  if      (method==MONTHINT_UNIFORM)//uniform over month
  {
    return aVal[month-1];
  }
  else if (method==MONTHINT_LINEAR_FOM)//linear from first of month
  {
    mo=month-1;
    nextmo=mo+1;
    if (nextmo==12){nextmo=0;}
    leap=0;if ((IsLeapYear(year,Options.calendar)) && (mo==1)){leap=1;}
    wt=1.0-(double)(day)/(DAYS_PER_MONTH[mo]+leap);
    return wt*aVal[mo]+(1-wt)*aVal[nextmo];
  }
  else if ((method==MONTHINT_LINEAR_21) ||
           (method==MONTHINT_LINEAR_MID))
    //linear from 21st of month to 21st of next month (e.g., UBC_WM) or other day
  {
    double pivot=0.0;
    if      (method==MONTHINT_LINEAR_21){pivot=21;}
    else if (method==MONTHINT_LINEAR_MID){
      pivot=0.5*DAYS_PER_MONTH[month-1];
      if ((IsLeapYear(year,Options.calendar)) && (month==2)){pivot+=0.5;}
    }

    if (day <= pivot)
    {
      mo=month-2;
      nextmo=mo+1;
      if (mo==-1){mo=11;nextmo=0;}
      leap=0;if ((IsLeapYear(year,Options.calendar)) && (mo==1)){leap=1;}
      wt=1.0-(double)(day+(DAYS_PER_MONTH[mo]+leap)-pivot)/(double)(DAYS_PER_MONTH[mo]+leap);
    }
    else{
      mo=month-1;
      nextmo=mo+1;
      int lastmo=mo-1;
      if (lastmo==-1){lastmo=11;}
      if (nextmo==12){nextmo=0;}
      leap=0;if ((IsLeapYear(year,Options.calendar)) && (lastmo==1)){leap=1;}
      wt=1.0-(double)(day-pivot)/(double)(DAYS_PER_MONTH[lastmo]+leap);
    }

    //\math \f$ wt=0.5-0.5*cos(wt*PI) \f$ ; //Useful smoothing function

    return wt*aVal[mo]+(1-wt)*aVal[nextmo];
  }
  return 0.0;
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates saturation vapor pressure [KPa] \cite Murray1966JAM
/// \remark Uses Dingman equation 7.4 \cite Dingman1994
// (Murray, Applied Meteorol 6:203, 1967)
///
/// \param &T [in] Temperature in Celsius
/// \return Saturated vapor pressure [kPa] corresponding to passed temperature
//
double GetSaturatedVaporPressure(const double &T)//[C]
{
  const double A1=0.61078;
  const double A2=17.26939;
  const double A3=237.3;
  const double A4=21.87456;
  const double A5=265.5;

  //0.61115*exp(22.452*T/(T+ZERO_CELSIUS)); //MESH

  //0.611*exp(17.27*Ta/(Ta+237.3)); //Common simplification of Dingman for T>0

  if (T>=0){return A1*exp(A2*T/(T+A3));}  // Dingman/Brook90 version (Tetens equation, 1930)
  else     {return A1*exp(A4*T/(T+A5));}
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates saturation vapor pressure slope [de/dT]
/// \remark Uses Dingman equation 7.4 \cite Dingman1994 (Murray, Applied Meteorol 6:203, 1967) \cite Murray1966JAM
///
/// \param &T [in] Temperature in Celsius
/// \param &e_sat [in] Saturated vapour pressure [kPa]
/// \return Saturated vapor pressure corresponding to passed temperature
//
double GetSatVapSlope(const double &T, const double &e_sat)
{
  const double A2=17.26939;
  const double A3=237.3;
  const double A4=21.87456;
  const double A5=265.5;
  //calculate d(sat_vap)/dT - Dingman 7.6 , SWAT 1:2.3.4
  if (T>0){return A2*A3/pow(T+A3,2)*e_sat;}
  else    {return A4*A5/pow(T+A5,2)*e_sat;}

  //from CRHM routine ClassCRHMCanopyClearingGap::delta
  //if(T>0) { return(2504.0*exp(17.27 * T/(T+237.3)) / sqrt(T+237.3)); }
  //else    { return(3549.0*exp(21.88 * T/(T+265.5)) / sqrt(T+265.5)); }

}

//////////////////////////////////////////////////////////////////
/// \brief Calculates latent heat of vaporization [MJ/kg]
/// \remark Uses Dingman equation 7-8 (Harrison, 1963) \cite Dingman1994  \cite Harrison1963HaM
///
/// \param &T [in] Temperature in Celsius
/// \return Latent heat of vaporization [MJ/kg]
//
double GetLatentHeatVaporization(const double &T)
{
  return 2.495-0.002361*T;//[MJ/kg]
}

//////////////////////////////////////////////////////////////////
/// \brief Returns latent psychrometric constant  [kPa/K] \cite Brunt1952
/// \remark Uses Dingman eqn. 7-13 \cite Dingman1994, SWAT 1:2.3.7\cite Neitsch2005 (Brunt, 1952)
///
/// \param &P [in] Atmospheric pressure [kPa]
/// \param &LH_vapor [in] Latent heat of vaporization [MJ/kg]
/// \return Psychrometric 'constant'  [kPa/K] (~0.07)
//
double GetPsychrometricConstant  (const double &P,const double &LH_vapor)
{
  return SPH_AIR/AIR_H20_MW_RAT*P/LH_vapor;//[kPa/K];
}

//////////////////////////////////////////////////////////////////
/// \brief returns wet bulb temperature, calculated using newton Raphson method
//
/// \param &P [in] Atmospheric pressure [kPa]
/// \param &Ta [in] air temperature
/// \return Psychrometric 'constant'  [kPa/K] (~0.07)
//
double GetWetBulbTemperature  (const double &P,const double &Ta,const double &rel_hum)
{
  double LH_vapor=GetLatentHeatVaporization(Ta);
  double gamma   =GetPsychrometricConstant(P,LH_vapor);
  double e_sat   =GetSaturatedVaporPressure(Ta);
  double e_a     =e_sat*rel_hum;

  double change=10;
  double tolerance=0.025;
  double Twet=Ta; //initial guess
  double F,dF,e_wet;

  while (fabs(change) > tolerance) //solve using Newton's method (usually ~2-3 iters)
  {
    e_wet=GetSaturatedVaporPressure(Twet);
    F=Ta-Twet-1/gamma*(e_wet-e_a);
    dF=-1.0-1/gamma*(GetSatVapSlope(Twet,e_wet));
    change=-F/dF;
    Twet+=change;
  }

  return Twet; //[C]
}

//////////////////////////////////////////////////////////////////
/// \brief Returns air density [kg/m3]
/// \remark From CLM Manual, pg. 12 \cite Oleson2012
///
/// \param &T [in] Temperature in Celsius
/// \param &P [in] Atmospheric pressure [kPa]
/// \return Air density [kg/m3]
//
double GetAirDensity(const double &T, const double &P)
{
  double e=GetSaturatedVaporPressure(T);
  return (P-0.378*e)/(DRY_GAS_CONST*(T+ZERO_CELSIUS))*GRAMS_PER_KG;
}

//////////////////////////////////////////////////////////////////
/// \brief Converts relative humidity to specific humidity
///
/// \param rel_hum [in] relative humidity [0-1]
/// \param air_press [in] Air pressure [kPa]
/// \param T [in] Air temperature [Celsius]
/// \return specific humidity, [kg/kg]
//
double GetSpecificHumidity(const double &rel_hum,const double &air_press,const double &T)
{
  double e_a=rel_hum*GetSaturatedVaporPressure(T);
  //simplified:
  //return AIR_H20_MW_RAT*e_a/air_press;

  return AIR_H20_MW_RAT*e_a/(air_press- e_a*(1-AIR_H20_MW_RAT));
}

//////////////////////////////////////////////////////////////////
/// \brief Converts passed temperature from Celsius to Farenheit
///
/// \param &T [in] Temperature in Celsius
/// \return Double temperature in Fahrenheit
//
double CelsiusToFarenheit(const double &T)
{
  return 1.8*T+32.0;
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates vertical transport efficiency of water vapor by turbulent eddies
/// \remark From Dingman pg. 273 \cite Dingman1994
///
/// \param &P [in] air pressure [kPa]
/// \param &ref_ht [in] Measurement height [m]
/// \param &zero_pl [in] Zero plane displacement [m]
/// \param &z0 [in] Coefficient of roughness
/// \return Vertical transport efficiency [m s^2/kg]
//
double GetVerticalTransportEfficiency(const double &P,
                                      const double &ref_ht,
                                      const double &zero_pl,
                                      const double &z0)
{
  double numer,denom;

  numer = AIR_H20_MW_RAT*DENSITY_AIR;
  denom = P*DENSITY_WATER*(1.0/pow(VON_KARMAN,2)*(pow((log((ref_ht - zero_pl)/z0)),2)));

  return numer/denom; //[m s^2 Kg^-1]
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates atmospheric conductivity [mm/s] for ET calculations
/// \ref From Dingman eqn. 7-49 \cite Dingman1994, Howell, T.A and Evett, S.R., USDA-ARS \cite Howell2004
///
/// \param &wind_vel [in] Wind velocity [m/d]
/// \param &meas_ht [in] Measurement height of wind vel [m] - must be greater than zero plane displacement
/// \param &zero_pl [in] Zero plane displacement [m]
/// \param &rough_ht [in] Roughness height [m]
/// \param &vap_rough_ht [in] Vapour roughness height [m]
/// \return Atmospheric conductivity [mm/s]
//
double CalcAtmosphericConductance(const double &wind_vel,     //[m/d]
                                  const double &meas_ht,      //[m]
                                  const double &zero_pl,      //[m]
                                  const double &rough_ht,     //[m]
                                  const double &vap_rough_ht) //[m]
{
  double atmos_cond;
  if (zero_pl==0.0){return 0.0;}

  //'6.25' from Dingman equation 7-49 is roughly 1/VK^2 (~6)
  atmos_cond=(wind_vel*MM_PER_METER*pow(VON_KARMAN,2));
  atmos_cond/=(log((meas_ht-zero_pl)/rough_ht)*log((meas_ht-zero_pl)/vap_rough_ht));

  return atmos_cond;//[mm/s]
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates dew point temperature, in celsius
///
/// \param &Ta [in] Air temperature in Celsius
/// \param &rel_hum [in] Relative humidity [0..1]
/// \ref Magnus-Tetens approximation, Murray, F. W., On the computation of saturation vapor pressure, J. Appl. Meteorol., 6, 203-204, 1967. \cite Murray1966JAM
/// \return Dew point temperature [C]
//
double GetDewPointTemp(const double &Ta,      //air temp, [C]
                       const double &rel_hum)//relative humidity [0..1]
{
  const double a=17.27;
  const double b=237.7;//[C]

  double tmp=(a*Ta/(b+Ta))+log(rel_hum);
  return b*tmp/(a-tmp);
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates dew point temperature
/// \ref from Dingman eqn D-11 \cite Dingman1994; can also be used to estimate rain temperature
///
/// \param &e [in] vapour pressure [kPa]
/// \return Dew point temperature [C]
//
double GetDewPointTemp(const double &e)
{
  double numer,denom;

  numer=  log(e)+0.4926;
  denom=  0.0708 - 0.00421*log(e);

  return numer/denom; //[C]
}
//////////////////////////////////////////////////////////////////
/// \brief converts volumetric enthalpy of water/ice only [MJ/m3 water] to temperature
///
/// \param hv [in] volumetric enthapy [MJ/m3 water]
/// \return water temperature [C]
//
double ConvertVolumetricEnthalpyToTemperature(const double &hv)
{

  if      (g_disable_freezing)         {
    if(hv/SPH_WATER/DENSITY_WATER  > 40) { return  40.0; } //upper limit - due to small volume, small energy roundoff error
    return hv/SPH_WATER/DENSITY_WATER;
  }

  double g_freeze_temp=-0.0;
  double hvt=(g_freeze_temp*SPH_ICE-LH_FUSION)*DENSITY_WATER;//transition enthalpy [MJ/m3 water]

  if(hv/SPH_WATER/DENSITY_WATER  > 40) { return  40.0; } //upper limit - due to small volume, small energy roundoff error
  if((hv-hvt)/SPH_ICE/DENSITY_ICE<-40) { return -40.0; } //lower limit

  if      (hv>0  ){ return hv/SPH_WATER/DENSITY_WATER;}
  else if (hv>hvt){ return g_freeze_temp*hv/hvt;}
  else            { return (hv-hvt)/SPH_ICE/DENSITY_ICE; }
}
//////////////////////////////////////////////////////////////////
/// \brief converts volumetric enthalpy of water/ice only [MJ/m3 water] to ice fraction
///
/// \param hv [in] volumetric enthapy [MJ/m3 water]
/// \return ice fraction [0..1]
//
double ConvertVolumetricEnthalpyToIceContent(const double &hv)
{
  if      (g_disable_freezing)            { return 0; }

  double g_freeze_temp=-0.0;
  double hvt=(g_freeze_temp*SPH_ICE-LH_FUSION)*DENSITY_WATER;//transition enthalpy [MJ/m3 water]

  if      (hv>0  ){ return 0;}
  else if (hv>hvt){ return hv/hvt;}
  else            { return 1.0; }
}

//////////////////////////////////////////////////////////////////
/// \brief calculates dT/dh corresponding to volumetric enthalpy hv
///
/// \param hv [in] volumetric enthapy [MJ/m3 water]
/// \return dT/dh [C-m3/MJ]
//
double TemperatureEnthalpyDerivative(const double &hv)
{
  if(g_disable_freezing)                  { return 1.0/SPH_WATER/DENSITY_WATER; }

  double g_freeze_temp=-0.0;
  double hvt=(g_freeze_temp*SPH_ICE-LH_FUSION)*DENSITY_WATER; //transition enthalpy [MJ/m3 water]

  //double alpha=0.05; //smoothing factor
  //return 0.5*rvn_erfc(alpha*(hv-hvt))*1.0/SPH_ICE/DENSITY_ICE+0.5*erfc(alpha*(-hv))*1.0/SPH_WATER/DENSITY_WATER;

  //derivative of temperature with respect to volumetric enthalpy
  if      (hv>0  ){ return 1.0/SPH_WATER/DENSITY_WATER;}
  else if (hv>hvt){ return g_freeze_temp/hvt;}
  else            { return 1.0/SPH_ICE/DENSITY_ICE; }   //JRC: NOTE - since all water is in SWE, not ice volume [mm], shouldnt this be DENSITY_WATER?? I think so.
}

//////////////////////////////////////////////////////////////////
/// \brief converts temperature to volumetric enthalpy of water/ice only [MJ/m3 water]
///
/// \param T [in] water temperature [C]
/// \param pctfroz [in] percent frozen [0..1]
/// \return volumetric enthapy [MJ/m3 water]
//
double ConvertTemperatureToVolumetricEnthalpy(const double &T,const double &pctfroz)
{
  if(g_disable_freezing) { return T*SPH_WATER*DENSITY_WATER; }

  double g_freeze_temp=-0.0;
  double hvt=(g_freeze_temp*SPH_ICE-LH_FUSION)*DENSITY_WATER;//transition enthalpy [MJ/m3 water]

  if     ((g_freeze_temp==0.0) && (fabs(T)<REAL_SMALL)){ return pctfroz*hvt; } //along zero line
  else if(T>0)                                         { return T*SPH_WATER*DENSITY_WATER; }
  else if((hvt<0.0) && (T>g_freeze_temp))              { return pctfroz*hvt; }
  else                                                 { return T*SPH_ICE  *DENSITY_ICE+hvt; } //JRC: NOTE - since all water is in SWE, not ice volume [mm], shouldnt this be DENSITY_WATER?? I think so.
}

//////////////////////////////////////////////////////////////////
/// \brief Converts any lowercase characters in a string to uppercase, returning the converted string
/// \param &s [in] String to be converted to uppercase
/// \return &s converted to uppercase
//
string StringToUppercase(const string &s)
{
  string ret(s.size(), char());
  for(int i = 0; i < (int)(s.size()); ++i)
  {
    if ((s[i] <= 'z' && s[i] >= 'a')){ret[i] =  s[i]-('a'-'A');}
    else                             {ret[i]  = s[i];}
  }
  return ret;
}

//////////////////////////////////////////////////////////////////
/// \brief Simple and fast atof (ascii to float) function.
/// \notes Executes about 5x faster than standard MSCRT library atof().
///
/// \notes ported 09-May-2009 from Tom Van Baak (tvb) www.LeapSecond.com
//

#define white_space(c) ((c) == ' ' || (c) == '\t')
#define valid_digit(c) ((c) >= '0' && (c) <= '9')

double fast_s_to_d (const char *p)
{
  int frac;
  double sign, value, scale;

  // Skip leading white space, if any.

  while (white_space(*p) ) {
    p += 1;
  }

  // Get sign, if any.

  sign = 1.0;
  if (*p == '-') {
    sign = -1.0;
    p += 1;

  } else if (*p == '+') {
    p += 1;
  }

  // Get digits before decimal point or exponent, if any.

  for (value = 0.0; valid_digit(*p); p += 1) {
    value = value * 10.0 + (*p - '0');
  }

  // Get digits after decimal point, if any.

  if (*p == '.') {
    double pow10 = 10.0;
    p += 1;
    while (valid_digit(*p)) {
      value += (*p - '0') / pow10;
      pow10 *= 10.0;
      p += 1;
    }
  }

  // Handle exponent, if any.

  frac = 0;
  scale = 1.0;
  if ((*p == 'e') || (*p == 'E')) {
    unsigned int expon;

    // Get sign of exponent, if any.

    p += 1;
    if (*p == '-') {
      frac = 1;
      p += 1;

    } else if (*p == '+') {
      p += 1;
    }

    // Get digits of exponent, if any.

    for (expon = 0; valid_digit(*p); p += 1) {
      expon = expon * 10 + (*p - '0');
    }
    if (expon > 308) expon = 308;

    // Calculate scaling factor.

    while (expon >= 50) { scale *= 1E50; expon -= 50; }
    while (expon >=  8) { scale *= 1E8;  expon -=  8; }
    while (expon >   0) { scale *= 10.0; expon -=  1; }
  }

  // Return signed and scaled floating point result.

  return sign * (frac ? (value / scale) : (value * scale));
}


//////////////////////////////////////////////////////////////////
/// \brief Converts any string to corresponding HRU Type
/// \param &s [in] String to be converted to uppercase
/// \return type, defaults to standard (doesn't complain)
//
HRU_type StringToHRUType(const string s)
{
  string sup;
  sup=StringToUppercase(s);

  if      (!s.compare("GLACIER" )){return HRU_GLACIER;}
  else if (!s.compare("LAKE"    )){return HRU_LAKE;}
  else if (!s.compare("WATER"   )){return HRU_WATER;}
  else if (!s.compare("ROCK"    )){return HRU_ROCK;}
  else if (!s.compare("PAVEMENT")){return HRU_ROCK; }
  else if (!s.compare("WETLAND" )){return HRU_WETLAND;}
  else if (!s.compare("STANDARD")){return HRU_STANDARD;}

#ifdef _STRICTCHECK_
  ExitGracefully("StringToHRUType: unrecognized hru type code",BAD_DATA);
#endif

  return HRU_INVALID_TYPE;
}

//////////////////////////////////////////////////////////////////
/// \brief returns true if line is empty, begins with '#' or '*'
/// \param &s [in] first string token in file line
/// \param Len length of line
/// \return true if line is empty or a comment
//
bool IsComment(const char *s, const int Len)
{
  if ((Len==0) || (s[0]=='#') || (s[0]=='*')){return true;}
  return false;
}
//////////////////////////////////////////////////////////////////
/// \brief replaces all instances of substring 'from' with string 'to' in string str
/// \param &str [in/out] string subjected to modification
/// \param &from [in] substring to be replaced
/// \param &to [in]  substring to replace it with
/// from solution by Michael Mrozek in https://stackoverflow.com/questions/3418231/replace-part-of-a-string-with-another-string
//
void SubstringReplace(string &str,const string &from,const string &to)
{
  if(from.empty()) { return; }
  size_t start_pos = 0;
  while((start_pos = str.find(from,start_pos)) != std::string::npos) {
    str.replace(start_pos,from.length(),to);
    start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
  }
}

/////////////////////////////////////////////////////////////////
/// \brief writes warning to screen and to Raven_errors.txt file
/// \param warn [in] warning message printed
//
void WriteWarning(const string warn, bool noisy)
{
  if (!g_suppress_warnings){
    ofstream WARNINGS;
    WARNINGS.open((g_output_directory+"Raven_errors.txt").c_str(),ios::app);
    if (noisy){cout<<"WARNING!: "<<warn<<endl;}
    WARNINGS<<"WARNING  : "<<warn<<endl;
    WARNINGS.close();
  }
}
/////////////////////////////////////////////////////////////////
/// \brief writes advisory to screen and to Raven_errors.txt file
/// \param warn [in] warning message printed
//
void WriteAdvisory(const string warn, bool noisy)
{
  if (!g_suppress_warnings){
    ofstream WARNINGS;
    WARNINGS.open((g_output_directory+"Raven_errors.txt").c_str(),ios::app);
    if (noisy){cout<<"ADVISORY: "<<warn<<endl;}
    WARNINGS<<"ADVISORY : "<<warn<<endl;
    WARNINGS.close();
  }
}
///////////////////////////////////////////////////////////////////
/// \brief NetCDF error handling
/// \return Error string and NetCDF exit code
//
void HandleNetCDFErrors(int error_code){

#ifdef _RVNETCDF_
  if(error_code==0){ return; }
  else{
    string warn;
    warn="NetCDF error ["+ to_string(nc_strerror(error_code))+"] occured.";
    ExitGracefully(warn.c_str(),BAD_DATA);
  }
#endif
}
///////////////////////////////////////////////////////////////////
/// \brief Return AUTO_COMPUTE tag if passed string is tagged, otherwise convert to double
/// \param s [in] Input string
/// \return AUTO_COMPUTE or USE_TEMPLATE_VALUE tag if string is tagged, otherwise double conversion of string
//
double AutoOrDouble(const string s)
{
  if (!s.compare("_AUTO"   )){return AUTO_COMPUTE;}
  if (!s.compare("AUTO"    )){return AUTO_COMPUTE;}
  if (!s.compare("_DEFAULT")){return USE_TEMPLATE_VALUE;}
  if (!s.compare("_DEF"    )){return USE_TEMPLATE_VALUE;}
  return s_to_d(s.c_str());
}
///////////////////////////////////////////////////////////////////
/// \brief Returns same number unless really small and g_suppress_zeros is true, then zeros out
/// \param d [in] Input double
/// \return d or zero, if d is near zero
//
double FormatDouble(const double &d)
{
  if((g_suppress_zeros) && (fabs(d)<REAL_SMALL)){return 0.0;}
  return d;
}

//////////////////////////////////////////////////////////////////
/// \brief Determines the complementary error of passed double &x (i.e., erfc(x))
/// \remark From Charbeneau, Groundwater Hydraulics and Pollutant Transport, 2000\cite Charbeneau2002AMR
/// \param &x [in] Double input into the function
/// \return Complementary error of &x, erfc(x)
//
double rvn_erfc(const double &x)
{
  double tmp(fabs(x)); //take abs so that we are always in positive quadrant.
  double fun;
  double f1;
  double tmp2;
  double tmp3;

  if(tmp > 3.0){
    f1  = (1.0 - 1.0/(2.0 * tmp * tmp)
           + 3.0/(4.0 * pow(tmp,4))
           - 5.0/(6.0 * pow(tmp,6)));
    fun = f1 * exp(-tmp * tmp) / (tmp * sqrt(PI));
  }
  else{
    tmp2 = 1.0 / (1.0 + (0.3275911 * tmp));
    tmp3 =   0.254829592  * tmp2       //5th order polynomial interpolation
      - (0.284496736 * tmp2 * tmp2)
      + (1.421413741 * pow(tmp2,3))
      - (1.453152027 * pow(tmp2,4))
      + (1.061405429 * pow(tmp2,5));
    fun = tmp3 * exp(-tmp * tmp);
  }
  if (tmp == x) {return fun;}
  else{return (2-fun);}
}
//////////////////////////////////////////////////////////////////
/// \brief Calculates the error function of passed value x
/// \details Uses pre-defined complementary error function to define error function
///
/// \param &x [in] Double input into the funciton
/// \return Error of &x, erf(x)
//
double rvn_erf(const double &x)
{
  return 1-rvn_erfc(x);
}
//////////////////////////////////////////////////////////////////
/// local versions of math functiosn trunc, floor, round, and isnan
/// - not available for standard BORLAND C compiler
//
double rvn_trunc(const double& x) {
#ifdef __BORLANDC__
  return pdMath_trunc(x);
#else
  return trunc(x);
#endif
}
double rvn_round(const double& x) {
#ifdef __BORLANDC__
  return pdMath_round(x);
#else
  return round(x);
#endif
}
double rvn_floor(const double& x) {
#ifdef __BORLANDC__
  return pdMath_floor(x);
#else
  return floor(x);
#endif
}
bool rvn_isnan(const double& x) {
#ifdef __BORLANDC__
  return (bool)(pdMath_isnan(x));
#else
  return std::isnan(x);
#endif
}
////////////////////////////////////////////////////////////////////
/// \brief Returns the value of the lognormal distribution function f_x(x) for the specified value x
/// \param &x [in] Double whose lognormal probability distribution value is to be returned
/// \param &mu [in] Mean of transformed function Y = ln(x)
/// \param &sig [in] Standard deviation of transformed distribution Y = ln(x)
/// \return Lognormal probability distribution value at x
//
double log_pdf(const double &x, const double &mu, const double &sig)
{
  if (x<=0){return 0.0;}
  return 1.0/x/sig/sqrt(2.0*PI)*exp(-0.5*(pow((log(x)-mu)/sig,2)));
}
/////////////////////////////////////////////////////////////////
/// \brief returns Nth recursive approximation of Lambert W_{-1} function
/// \remark From D.A. Barry et al. Used in some Green Ampt infiltration schemes \cite Barry2005AiWR.
//
double LambertN(const double &x, const int N)
{
  double sigma,tmp;
  sigma=pow(-2.0-2.0*log(-x),0.5);
  if (N<=2){tmp=-1.0-0.5*sigma*sigma-sigma/(1.0+sigma/6.3);}
  else     {tmp=LambertN(x,N-1);}
  return (1+log(-x)-log(-tmp))*tmp/(1+tmp);
}

/////////////////////////////////////////////////////////////////
/// \brief Returns gamma function of argument x
/// \note Returns 1e308 if argument is a negative integer or 0, or if argument exceeds 171.
///
/// \param x [in] Argument whose gamma function will be determined
/// \return Gamma function of argument x
//
double gamma2(double x)
{
  int i,k,m;
  double ga=0,gr,r=0,z;

  static double g[] = {
    1.0,                  0.5772156649015329, -0.6558780715202538,
    -0.420026350340952e-1, 0.1665386113822915, -0.421977345555443e-1,
    -0.9621971527877e-2,   0.7218943246663e-2, -0.11651675918591e-2,
    -0.2152416741149e-3,   0.1280502823882e-3, -0.201348547807e-4,
    -0.12504934821e-5,     0.1133027232e-5,    -0.2056338417e-6,
    0.6116095e-8,         0.50020075e-8,      -0.11812746e-8,
    0.1043427e-9,         0.77823e-11,        -0.36968e-11,
    0.51e-12,            -0.206e-13,          -0.54e-14,
    0.14e-14};

  if (x > 171.0){return 0.0;}    // This value is an overflow flag.
  if (x == (int)x)
  {
    if (x > 0.0) {
      ga = 1.0;               // use factorial
      for (i=2;i<x;i++) {ga *= i;}
    }
    else{
      ga=0.0;
      ExitGracefully("Gamma:negative integer values not allowed",RUNTIME_ERR);
    }
  }
  else
  {
    z=x;
    r=1.0;
    if (fabs(x) > 1.0) {
      z = fabs(x);
      m = (int)(z);
      r = 1.0;
      for (k=1;k<=m;k++) {r *= (z-k);}
      z -= m;
    }
    gr = g[24];
    for (k=23;k>=0;k--) {gr = gr*z+g[k];}
    ga = 1.0/(gr*z);
    if (fabs(x) > 1.0) {
      ga *= r;
      if (x < 0.0) {ga = -PI/(x*ga*sin(PI*x));}
    }
  }
  return ga;
}
/////////////////////////////////////////////////////////////////
/// \brief returns value of incomplete gamma function  with parameter a for input x
/// incomplete Gamma is g(x,a)=int_0^x of t^a-1 exp(-t) dt
/// \param &x [in] upper limit of integral
/// \param &a [in] shape parameter
/// \return Incomplete gamma function g(x,a)
//
double IncompleteGamma(const double &x, const double &a)
{
  //cumulative distribution
  /// \ref from http://algolist.manual.ru/maths/count_fast/gamma_function.php
  const int N=100;
  if (x<=0){return 0.0;}
  double num=1;
  double sum=0.0;
  double prod=1.0;
  for (int n=0;n<N;n++){
    if (n>0){num*=x;}
    prod*=(a+n);
    sum+=num/prod;
  }
  return sum*pow(x,a)*exp(-x);
}

/////////////////////////////////////////////////////////////////
/// \brief returns value of Gamma distribution for argument x, with parameters alpha and beta
/// gamma(x,a,b)=b^a/Gamma(a)*x^(a-1)*exp(-b*x)
///
/// \param &x [in] Argument x whose Gamma distribution value will be determined
/// \param &alpha [in] shape parameter
/// \param &beta [in] scale parameter
/// \return Gamma distribution value
//
double GammaDist(const double &x, const double &alpha, const double &beta)
{
  //mean=alpha/beta
  double bx=beta*x;
  return pow(bx,alpha)/x/gamma2(x)*exp(-bx);
}
//////////////////////////////////////////////////////////////////
/// \brief Calculates cumulative two parameter cumulative gamma distribution \cite Clark2008WRR
///
/// \param &t [in] time
/// \param &alpha [in] shape parameter
/// \param &beta [in] scaling parameter
/// \return integrated gamma distribution from 0..t
//
double GammaCumDist(const double &t,  const double &alpha,const double &beta)
{
  return IncompleteGamma(beta*t,alpha)/gamma2(alpha);
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates cumulative triangular distribution
/// \remark Area under=1.0-->peak=2/tp
///
/// \param &t [in] argument of triangular distribution
/// \param &tc [in] End of triangle
/// \param &tp [in] Peak of triangle
/// \return Cumulative triangular distribution value for input t
//
double TriCumDist(const double &t, const double &tc, const double &tp)
{
  double b=2.0/tc;
  double m;
  if (t<0.0){return 0.0;}
  if (t<=tp){
    m=(b/tp);
    return 0.5*m*t*t;
  }
  else if (t<=tc){
    m=-b/(tc-tp);
    return tp/tc+b*(t-tp)+0.5*m*(t-tp)*(t-tp);
  }
  else{
    return 1.0;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns the cumulative distribution function at input t for a sequence of linear reservoirs
/// \remark Basically the Gamma distribution for integer shape parameters. A common unit hydrograph format \n
///  - \math \f$ PDF(t)/UH(t)=t^{N-1}k^{N}e^{-kt} \f$
///  - \math \f$ CDF(t)/cum UH(t)=1-e^{-kt}\sum_{n=0}^{N-1}t^n/n! \f$
///
/// \param &t [in] The input value whose CDF is to be determined
/// \param &k [in] CDF Parameter (linear storage coeff)
/// \param &NR [in] Integer number of reservoirs
/// \return CDF at point t
//
double NashCumDist(const double &t, const double &k, const int &NR)
{
  if (t<0.0){return 0.0;}
  double fact=1.0;
  double prod=1.0;
  double sum =1.0;
  for (int n=1; n<NR; n++){
    fact=fact*n;
    prod=prod*(k*t);
    sum+=prod/fact;
  }
  return 1.0-exp(-k*t)*sum;
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates cumulative kinematic wave solution distribution
///
/// \param &t time [d]
/// \param &L reach length [m]
/// \param &v celerity [m/d]
/// \param &D diffusivity [m2/d]
/// \return Returns cumulative kinematic wave solution distribution
//
//int_0^time L/2/t^(3/2)/sqrt(pi*D)*exp(-(v*t-L)^2/(4*D*t)) dt
// extreme case (D->0): =1 for v*t<L, 0 otherwise
double ADRCumDist(const double &t, const double &L, const double &v, const double &D)
{
  ExitGracefullyIf(D<=0,"ADRCumDist: Invalid diffusivity",RUNTIME_ERR);

 /* double term =L/2.0*pow(PI*D,-0.5);
  double beta =L/sqrt(D);
  double alpha=v/sqrt(D);
  double dt   =t/5000.0; //t in [day] (5000.0 is # of integral divisions)
  double integ=0.0;
  for (double tt=0.5*dt;tt<t;tt+=dt)
  {
    // D: unit = m2 / day
    // [m]/[day]^1.5/[m]/[day]^0.5 * exp([m/day]*[day]-[m])^2/([m]^2)/[day])
    integ+=pow(tt,-1.5)*exp(-((alpha*tt-beta)*(alpha*tt-beta))/4.0/tt);
    // equivalent to (old) version, but more stable since alpha, beta ~1-10 whereas both v^2 and D can be very large
    //integ+=pow(tt,-1.5)*exp(-((v*tt-L)*(v*tt-L))/4.0/D/tt); //old version
  }
  return integ*term*dt;
  */
  //Analytical- Ogata Banks, 1969
  double F=0;
  if (t<=0){return 0.0;}
  if(L < 500*(D/v)) {
    F=0.5*(exp((v*L)/D)*rvn_erfc((L+v*t)/sqrt(4.0*D*t)));
  }
  F+=0.5*(rvn_erfc((L-v*t)/sqrt(4.0*D*t)));
  return F;
}
//////////////////////////////////////////////////////////////////
/// \brief Calculates time-varying cumulative kinematic wave solution distribution
///
/// \param &t time [d]
/// \param &L reach length [m]
/// \param *v array of celerities [m/d] for nv timesteps of length dt
/// \param &nv number of celerity values in array
/// \param &D diffusivity [m2/d]
/// \param &dt timestep [d]
/// \return Returns cumulative kinematic wave solution distribution for temporally variable velocity
//
double TimeVaryingADRCumDist(const double &t,const double &L,const double *v, int nv,const double &D,const double &dt)
{
  ExitGracefullyIf(D<=0   ,"TimeVaryingADRCumDist: Invalid diffusivity",RUNTIME_ERR);
  ExitGracefullyIf(t>nv*dt,"TimeVaryingADRCumDist: not enough celerity terms in array",RUNTIME_ERR);
  int i;
  i=(int)(floor((t+TIME_CORRECTION)/dt));
  double tstar=v[i]/v[0]*(t-(double)(i)*dt);
  for(int j=0;j<i;j++) {
    tstar+=dt*v[j]/v[0];
  }
  double F=0;
  if (t<=0){return 0.0;}
  if(L < 500*(D/v[0])) {
    F=0.5*(exp((v[0]*L)/D)*rvn_erfc((L+v[0]*tstar)/sqrt(4.0*D*tstar)));
  }
  F+=0.5*(rvn_erfc((L-v[0]*tstar)/sqrt(4.0*D*tstar)));
  return F;
}
//////////////////////////////////////////////////////////////////
/// \brief calculates N process weights from N-1 numbers ranging from 0 to 1
///
/// \param *aVals [in/out] array of weight seeds (uniform numbers) between 0 and 1 (size:N-1)
/// \param *aWeights [in/out] array of weights between 0 and 1 (size:N)
/// \param N [in] size of aWeights to use
/// From Mai et al., The pie sharing problem: Unbiased sampling of N+1 summative weights, Environmental Modelling and Software, 148, 105282, doi:10.1016/j.envsoft.2021.105282, 2022
//
void CalcWeightsFromUniformNums(const double* aVals, double* aWeights, const int N)
{
    double sum = 0.0;
    for (int q = 0; q < (N - 1); q++) {
        aWeights[q] = (1.0 - sum) * (1.0 - pow(1.0 - aVals[q], 1.0 / (N - q - 1)));
        sum += aWeights[q];
    }
    aWeights[N - 1] = 1.0 - sum;
}

//////////////////////////////////////////////////////////////////
/// \brief Quicksort algorithm
/// \author coded by Ayman Khedr, 3A Environmental University of Waterloo
///
/// \param arr[] [in & out] Unordered array of doubles to be sorted
/// \param left [in] Left bound of sort
/// \param right [in] Right bound of sort (should be ==[size of array-1] outside of recursion call)
//
void quickSort(double arr[], int left, int right)
{
  if (right<=left){return;}//e.g., if array size==0 or 1
  int i = left, j = right;
  double tmp;
  double pivot = arr[(left + right) / 2];

  // partition
  while (i <= j) {
    while (arr[i] < pivot){i++;}
    while (arr[j] > pivot){j--;}
    if (i <= j) {
      tmp = arr[i];
      arr[i] = arr[j];
      arr[j] = tmp;
      i++;
      j--;
    }
  };

  // recursion
  if (left <  j){quickSort(arr, left,  j);}
  if (i < right){quickSort(arr, i, right);}
}
//////////////////////////////////////////////////////////////////
/// \brief adds value v to end of integer array a[] of original size n, expands a[] to size n+1
//
void pushIntoIntArray(int*&a, const int &v, int &n)
{
  int *tmp=new int [n+1];
  for (int i = 0; i < n; i++) { tmp[i]=a[i];}
  tmp[n]=v;
  if (n>0){delete [] a;}
  a=tmp;
  n++;
}
//////////////////////////////////////////////////////////////////
// given unsorted or sorted array arr[] returns rank of each term, where 0 indicates the largest value (none smaller) and N-1 the smallest value (all are greater)
// if values are not unique, gaps will form between ranks
// NOT OPTIMIZED
void getRanks(const double* arr, const int N, int* ranks)
{
  double soft=0;
	for(int i=0;i<N;i++){
    ranks[i]=0;
		for(int j=0;j<N;j++){
      if (arr[j] >= arr[i]+soft) {ranks[i]++;}
		}
	}
}
//////////////////////////////////////////////////////////////////
/// \brief Returns index of array for near searching (ordered search around a guessed array index)
///
/// \param i [in] near search index, ranges from 0 to size-1
/// \param guess_p [in] best guess of appropriate array index p in search
/// \param size [in] size of array being searched
/// \return Returns index of array for near searching (ordered search around a guessed array index)
/// \details: NearSearch should proceed as
/// a[size];     //some kind of array of items where we are looking for the 'correct' item
/// int guess=4; //a best guess of the index of a[] which we believe is close to the 'correct' one
/// for (int i=0; i<size;i++){
///   p=NearSearchIndex(i,guess_p,size);
///   if (comparison_operator using a[p] that needs to go through minimal a[p] entries){break;}
/// }
//
int NearSearchIndex(const int i, int guess_p, const int size)
{
  int p;
  if ((guess_p >= size) || (guess_p < 0)) { guess_p = 0; }//fix bad guess
  if ((i < 0) || (i >= size)) {
    ExitGracefully("NearSearchIndex: bad index", RUNTIME_ERR);
  }
  if (i % 2 == 0) { p = guess_p - ((i + 0) / 2); }
  else { p = guess_p + ((i + 1) / 2); }
  if (p < 0) { p += size; } //valid wraparound
  if (p >= size) { p -= size; }
  return p;
}
///////////////////////////////////////////////////////////////////////////
/// \brief identifies index location of value in uneven continuous list of sorted value ranges
///
/// \param &x [in] value for which the interval index is to be found
/// \param *ax [in] array of consecutive values from ax[0] to ax[N-1] indicating interval boundaries
/// \param N [in] size of array ax
/// \param iguess [in] best guess as to which interval x is in
/// \return interval index value (index refers to lower bound of interval, i.e., i indicates x is between ax[i] and ax[i+1]
/// \note returns -1 if outside of range
//
int SmartIntervalSearch(const double &x,const double *ax,const int N,const int iguess)
{
  int i=iguess;
  if((iguess>N-1) || (iguess<0)) { i=0; }
  if((x>=ax[i]) && (x<ax[i+1])) { return i; }

  int plus,plus2;
  for(int d=1;d<(int)(rvn_trunc(N/2)+1);d++)
  {
    plus =i+d;    if(plus >N-1) { plus -=N; } //wrap
    plus2=i+d+1;  if(plus2>N-1) { plus2-=N; } //wrap
    if((x>=ax[plus]) && (x<ax[plus2])) { return plus; }
    plus =i-d;    if(plus <0) { plus +=N; } //wrap
    plus2=i-d+1;  if(plus2<0) { plus2+=N; }   //wrap
    if((x>=ax[plus]) && (x<ax[plus2])) { return plus; }
  }
  return DOESNT_EXIST;
}


//////////////////////////////////////////////////////////////////
/// \brief interpolates value from rating curve
/// \param x [in] interpolation location
/// \param xx [in] array (size:N) of vertices ordinates of interpolant
/// \param y [in] array (size:N)  of values corresponding to array points xx
/// \param N size of arrays x and y
/// \returns y value corresponding to interpolation point
/// \note does not assume regular spacing between min and max x value
/// \note if below minimum xx, either extrapolates (if extrapbottom=true), or uses minimum value
/// \note if above maximum xx, always extrapolates
//
double InterpolateCurve(const double x,const double *xx,const double *y,int N,bool extrapbottom)
{
  static int ilast=0;
  if(x<=xx[0])
  {
    if(extrapbottom) { return y[0]+(y[1]-y[0])/(xx[1]-xx[0])*(x-xx[0]); }
    return y[0];
  }
  else if(x>=xx[N-1])
  {
    return y[N-1]+(y[N-1]-y[N-2])/(xx[N-1]-xx[N-2])*(x-xx[N-1]);//extrapolation-may wish to revisit
                                                                //return y[N-1];
  }
  else
  {
    //int i=0; while ((x>xx[i+1]) && (i<(N-2))){i++;}//Dumb Search
    int i=SmartIntervalSearch(x,xx,N,ilast);
    if(i==DOESNT_EXIST) { return 0.0; }
    ExitGracefullyIf(i==DOESNT_EXIST,"InterpolateCurve::mis-ordered list or infinite x",RUNTIME_ERR);
    ilast=i;
    if (fabs(xx[i+1]-xx[i]) < REAL_SMALL) { return (y[i]+y[i+1])/2; }  // x locations too close to each other
    return y[i]+(y[i+1]-y[i])/(xx[i+1]-xx[i])*(x-xx[i]);
  }
}
