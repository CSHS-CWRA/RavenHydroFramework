/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2021 the Raven Development Team
----------------------------------------------------------------*/
#include "RavenInclude.h"

//////////////////////////////////////////////////////////////////
/// reads calendar from "calendar" attribute of NetCDF, converts to Raven enumerated type
//
/// \param ncid [in] netCDF file ID
/// \param varid_t [in] netCDF time variable ID
/// \param filename [in] filename of netCDF file
/// \param &Options [in] Global model options information
/// \return calendar (PROLEPTIC_GREGORIAN if calendar not in file)
//
 int GetCalendarFromNetCDF(const int ncid,int varid_t,const string filename,const optStruct &Options) 
{
   int calendar=CALENDAR_PROLEPTIC_GREGORIAN;

#ifdef _RVNETCDF_
  int     retval;
  size_t  att_len;        // length of the attribute's text
  nc_type att_type;      // type of attribute
  

  retval = nc_inq_att(ncid,varid_t,"calendar",&att_type,&att_len);
  if(retval != NC_ENOTATT) 
  {
    HandleNetCDFErrors(retval);
    // if found, read and make sure '\0' is terminating character
    retval = nc_inq_attlen(ncid,varid_t,"calendar",&att_len);            HandleNetCDFErrors(retval);// inquire length of attribute's text

    char *calendar_t = new char[att_len + 1]; // allocate memory of char * to hold attribute's text
    retval = nc_get_att_text(ncid,varid_t,"calendar",calendar_t);        HandleNetCDFErrors(retval);// read attribute text
    calendar_t[att_len] = '\0';                // add string determining character

    calendar = StringToCalendar(calendar_t);

    if(calendar!=Options.calendar) {
      if((calendar==CALENDAR_GREGORIAN) && (Options.calendar==CALENDAR_PROLEPTIC_GREGORIAN)) {} //basically the same, not worth warning
      else {
        string warn;
        warn="[Critical]: NetCDF file "+filename+" does not use same calendar as the simulation. Use the :Calendar command to ensure consistency";
        WriteWarning(warn,Options.noisy);
      }
    }

    delete[] calendar_t;
  }
#endif
  return calendar;
}

 //////////////////////////////////////////////////////////////////
 /// reads time informaiton from NetCDF vector of times, calculates start_time, interval, and time zone of data
 //
 /// \param unit_t [in] netCDF time unit string
 /// \param calendar [in] netCDF calendar
 /// \param time [in] (preferably uniform) array of times, in seconds, hrs, minutes, or days
 /// \param filename [in] NetCDF filename
 /// \param tstep [out] uniform timestep of data
 /// \param start_day [out] julian start day of data
 /// \param start_yr [out] julian start year of data
 /// \param time_zone [out] time zone of data
 //
 void GetTimeInfoFromNetCDF(const char *unit_t,int calendar,const double *time,const int ntime,const string filename,double &tstep,double &start_day,int &start_yr,double &time_zone) 
 {
#ifdef _RVNETCDF_
   // -------------------------------
   // check if we have equal time steps in time data vector
   // -------------------------------  
   double delta_t = time[1] - time[0];
   for(int ii=1; ii<ntime-1; ii++)
   {
     double delta_t2 = time[ii+1] - time[ii];
     if(fabs(delta_t2 - delta_t) > 0.00001)
     {
       printf("\n\n\nGetTimeInfoFromNetCDF: file name: %s\n",filename.c_str());
       printf("\n\n\n              expected delta_t = %15.8lf\n",delta_t);
       printf("\n\n\n              derived  delta_t = %15.8lf\n",delta_t2);
       printf("\n\n\n                       i  =%i  t[i]  =%15.8lf\n",ii,time[ii]);
       printf("\n\n\n                       i+1=%i  t[i+1]=%15.8lf\n",ii+1,time[ii+1]);
       ExitGracefully("GetTimeInfoFromNetCDF: time steps are not equal in NetCDF input",BAD_DATA);
     }
   }

   //-----------------------------
   // extract julian start date and year, calculate time interval in days
   //----------------------------
   string unit_t_str = to_string(unit_t);
   if     (strstr(unit_t,"hours"  )){tstep   = (time[1] - time[0])/HR_PER_DAY;}
   else if(strstr(unit_t,"days"   )){tstep   = (time[1] - time[0]);}
   else if(strstr(unit_t,"minutes")){tstep   = (time[1] - time[0])/MIN_PER_DAY;}
   else if(strstr(unit_t,"seconds")){tstep   = (time[1] - time[0])/SEC_PER_DAY;}
   else {
     ExitGracefully("GetTimeInfoFromNetCDF: this unit in time is not implemented yet (only days, hours, minutes, seconds)",BAD_DATA);
   }
   ExitGracefullyIf(tstep<=0,"GetTimeInfoFromNetCDF: Interval is negative!",BAD_DATA);

   //Get julian date/year of time[0]
   GetJulianDateFromNetCDFTime(unit_t,calendar,time[0],start_day,start_yr);
#endif

}
 //////////////////////////////////////////////////////////////////
 /// reads converts Netcdf time since a reference date to a julian date 
 //
 /// \param unit_t_str [in] netCDF time unit string - includes information about reference date
 /// \param calendar [in] netCDF calendar
 /// \param time [in] time in seconds, hrs, minutes, or days since reference date
 /// \param julian_day [out] julian start day of data
 /// \param year [out] julian start year of data
 //
 void GetJulianDateFromNetCDFTime(const string unit_t_str,const int calendar,const double &time,double &julian_day,int &year)
 {
#ifdef _RVNETCDF_
   time_struct tt_ref;   //reference time in NetCDF time string (e.g., time in days since 2001-10-01)
   double delta_t;       //timeshift in days
   double time_zone;
   string unit_string;
   if     (strstr(unit_t_str.c_str(),"hours"  )){    unit_string="hours";   delta_t=time/HR_PER_DAY; }
   else if(strstr(unit_t_str.c_str(),"days"   )){    unit_string="days";    delta_t=time;            }
   else if(strstr(unit_t_str.c_str(),"minutes")){    unit_string="minutes"; delta_t=time/MIN_PER_DAY;}
   else if(strstr(unit_t_str.c_str(),"seconds")){    unit_string="seconds"; delta_t=time/SEC_PER_DAY;}
   else {
     ExitGracefully("GetJulianDateFromNetCDFTime: this unit in time is not implemented yet (only days, hours, minutes, seconds)",BAD_DATA);
   }
   tt_ref =TimeStructFromNetCDFString(unit_t_str,unit_string,calendar,time_zone);
   AddTime(tt_ref.julian_day,tt_ref.year,delta_t,calendar,julian_day,year);
#endif
}
 //////////////////////////////////////////////////////////////////
 /// reads reads in a vector of NetCDF times in integer/float/double format
 //
 /// \param ncid [in] netCDF file ID
 /// \param varid_t [in] netCDF time index
 /// \param ntime [in] number of time entries in vector
 /// \param my_time [out] vector of times ( in hours/minutes/days/seconds since reference date), pre-allocated to size ntime
 //
void GetTimeVectorFromNetCDF(const int ncid,const int varid_t,const int ntime,double *my_time) 
{
  
#ifdef _RVNETCDF_
  int     retval;
  nc_type type;

  // find out if time is an integer or double variable
  retval = nc_inq_var(ncid,varid_t,
                        NULL,     // out: name
                        &type,    // out: type: nc_INT, nc_FLOAT, nc_DOUBLE, etc.
                        NULL,     // out: dims
                        NULL,     // out: dimids
                        NULL);   HandleNetCDFErrors(retval);// out: natts

  if(type == NC_DOUBLE)    // time is given as double --> use as is
  {
    retval = nc_get_var_double(ncid,varid_t,&my_time[0]);  HandleNetCDFErrors(retval);
  }
  else if(type == NC_FLOAT)    // time is given as float (single precision) --> convert to double
  {
    float *ttime=NULL;
    ttime=new float[ntime];
    ExitGracefullyIf(ttime==NULL,"CForcingGrid::ForcingGridInit",OUT_OF_MEMORY);
    retval = nc_get_var_float(ncid,varid_t,&ttime[0]);     HandleNetCDFErrors(retval);
    for(int itime=0; itime<ntime;itime++) {       // loop over all time steps
      my_time[itime]=double(ttime[itime]);        // convert to double
    }
    delete[] ttime;
  }
  else if((type == NC_INT) || (type == NC_INT64))   // time is given as integer --> convert to double
  {
    int *ttime=NULL;
    ttime=new int[ntime];
    ExitGracefullyIf(ttime==NULL,"CForcingGrid::ForcingGridInit",OUT_OF_MEMORY);
    retval = nc_get_var_int(ncid,varid_t,&ttime[0]);       HandleNetCDFErrors(retval);
    for(int itime=0; itime<ntime;itime++) {       // loop over all time steps
      my_time[itime]=(double)ttime[itime];        // convert to double
    }
    delete[] ttime;
  }
  else
  {
    ExitGracefully("GetTimeVectorFromNetCDF: time variable in NetCDF is not of type DOUBLE, FLOAT, or INT.",BAD_DATA);
  }
#endif
}
