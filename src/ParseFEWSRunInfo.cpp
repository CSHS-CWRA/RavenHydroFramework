/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2020 the Raven Development Team
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Model.h"


//////////////////////////////////////////////////////////////////
/// \brief Parses Deltares FEWS RunInfo file
///
/// \param *&pModel [out] Reference to model object
/// \param &Options [out] Global model options information
/// \return True if operation is successful
//
bool ParseNetCDFRunInfoFile(CModel *&pModel, optStruct &Options)
{
  if (Options.runinfo_filename==""){return true;}

#ifdef _RVNETCDF_
  int ncid;          //NetCDF fileid
  int retval;        //return value of NetCDF routines
  size_t att_len;    //character string length
  double time_zone;

  if(Options.noisy) { cout<<"Opening runinfo file "<<Options.runinfo_filename<<endl; }
  retval = nc_open(Options.runinfo_filename.c_str(),NC_NOWRITE,&ncid);          HandleNetCDFErrors(retval);

  // Ingest start time=============================================================================
  int varid_startt;
  double start_time;
  retval = nc_inq_varid   (ncid,"start_time",&varid_startt);
  if(retval != NC_ENOTATT) 
  {
    HandleNetCDFErrors(retval);
    retval = nc_inq_attlen  (ncid,varid_startt,"units",&att_len);                 HandleNetCDFErrors(retval);
    char *unit_t=new char[att_len+1];
    retval = nc_get_att_text(ncid,varid_startt,"units",unit_t);                   HandleNetCDFErrors(retval);// read attribute text
    unit_t[att_len] = '\0';// add string determining character

    if(!IsValidNetCDFTimeString(unit_t))
    {
      cout<<"time unit string: "<<unit_t<<endl;
      ExitGracefully("ParseRunInfoFile: start_time time unit string is not in the format '[days/hours/...] since YYYY-MM-DD HH:MM:SS +0000' !",BAD_DATA);
    }

    retval=nc_get_var_double(ncid,varid_startt,&start_time);
    GetJulianDateFromNetCDFTime(unit_t,Options.calendar,start_time,Options.julian_start_day,Options.julian_start_year);
    if (time_zone!=0.0){ExitGracefully("ParseRunInfoFile: does not yet support time zone shifts",BAD_DATA_WARN); }

    if(Options.noisy) { cout<<"ParseRunInfoFile: read variable start_time from NetCDF: start day: "<<Options.julian_start_day<<" yr: "<< Options.julian_start_year<<endl; }

    delete [] unit_t;
  }

  // Ingest end time=============================================================================
  int varid_endt;
  double end_time;
  double end_day;
  int end_year;
  retval = nc_inq_varid(ncid,"end_time",&varid_endt);
  if(retval == NC_ENOTATT) {
    //end_time not found. Do nothing
  }
  else {
    HandleNetCDFErrors(retval);
    retval = nc_inq_attlen(ncid,varid_endt,"units",&att_len);                 HandleNetCDFErrors(retval);
    char *unit_t=new char[att_len+1];
    retval = nc_get_att_text(ncid,varid_endt,"units",unit_t);                 HandleNetCDFErrors(retval);// read attribute text
    unit_t[att_len] = '\0';// add string determining character

    if(!IsValidNetCDFTimeString(unit_t))
    {
      cout<<"time unit string: "<<unit_t<<endl;
      ExitGracefully("ParseRunInfoFile: end_time time unit string is not in the format '[days/hours/...] since YYYY-MM-DD HH:MM:SS +0000' !",BAD_DATA);
    }
    
    retval=nc_get_var_double(ncid,varid_endt,&end_time);
    GetJulianDateFromNetCDFTime(unit_t,Options.calendar,end_time,end_day,end_year);
    if(time_zone!=0.0) { ExitGracefully("ParseRunInfoFile: does not yet support time zone shifts",BAD_DATA_WARN); }

    Options.duration=TimeDifference(Options.julian_start_day,Options.julian_start_year,end_day,end_year,Options.calendar);

    if(Options.noisy) { cout<<"ParseRunInfoFile: read variable end_time from NetCDF: duration: "<<Options.duration<< " days.  end day: "<<end_day<<" yr: "<< end_year<<endl; }

    delete[] unit_t;
  }
  // Ingest working directory=====================================================================
  int varid_workdir;
  retval = nc_inq_varid(ncid,"work_dir",&varid_workdir);
  if(retval!= NC_ENOTATT) 
  {
    HandleNetCDFErrors(retval);

    char *workdir=new char[255];
    retval=nc_get_var(ncid,varid_workdir,workdir);

    Options.output_dir=to_string(workdir);
    if(Options.noisy) { cout<<"ParseRunInfoFile: read variable work_dir from NetCDF: "<<Options.output_dir<<endl; }
    delete [] workdir;
  }
  // Ingest properties      =====================================================================
  int varid_props;
  retval = nc_inq_varid(ncid,"properties",&varid_props);
  if(retval != NC_ENOTATT) 
  {
    HandleNetCDFErrors(retval);
    
    //  RunName
    retval = nc_inq_attlen(ncid,varid_props,"RunName",&att_len);                    
    if(retval !=NC_ENOTATT)
    {
      HandleNetCDFErrors(retval);
      char *runname=new char[att_len+1];
      retval = nc_get_att_text(ncid,varid_props,"RunName",runname);                   HandleNetCDFErrors(retval);// read attribute text
      runname[att_len] = '\0';// add string determining character

      Options.run_name=runname;
      if(Options.noisy) { cout<<"ParseRunInfoFile: read properties:RunName from NetCDF: "<<Options.run_name<<endl; }
      delete[] runname;
    }
  }

  retval = nc_close(ncid);         HandleNetCDFErrors(retval);

  ExitGracefully("TESTING FEWS RUNINFO",STUB);

  return true;
#else
  ExitGracefully("ParseRunInfoFile: Deltares FEWS Runinfo file can only be read using NetCDF version of Raven",BAD_DATA_WARN);
#endif
}