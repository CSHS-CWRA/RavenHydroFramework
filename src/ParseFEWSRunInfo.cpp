/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2021 the Raven Development Team
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
  //double time_zone;

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
    //if (time_zone!=0.0){ExitGracefully("ParseRunInfoFile: does not yet support time zone shifts",BAD_DATA_WARN); }

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
    //if(time_zone!=0.0) { ExitGracefully("ParseRunInfoFile: does not yet support time zone shifts",BAD_DATA_WARN); }

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
  return false;
  ExitGracefully("ParseRunInfoFile: Deltares FEWS Runinfo file can only be read using NetCDF version of Raven",BAD_DATA_WARN);
#endif
}
//////////////////////////////////////////////////////////////////
/// \brief Parses Deltares FEWS State info file
/// assumes state info file has 
///
/// \param *&pModel [out] Reference to model object
/// \param &Options [out] Global model options information
/// \return True if operation is successful
//
bool ParseNetCDFStateFile(CModel *&pModel,optStruct &Options,const double &t)
{
  if(Options.runinfo_filename=="") { return true; }

#ifdef _RVNETCDF_
  int ncid;          //NetCDF fileid
  int retval;        //return value of NetCDF routines
  //size_t att_len;    //character string length
                     //double time_zone;
  int HRUvecID;      //attribute ID of HRU vector
  size_t nHRUsLocal;    //size of HRU Vector

  //The state info file includes:
    //1)	A vector of size[1:NumHRUs] which includes HRU IDs which reference the HRUID used in the Raven .rvh file. The index i of each entry will be the same as that used by the matrix in part 2. Suggested variable name: HRUIDs.
    //2)	One or more matrices of size[1:NumTimes,1:NumHRUs] which includes the state variable being updated,using the same indexing[i=1:NumHRUs] as the vector of HRUIDs. Blank values will not be updated. Raven model HRUs with IDs not in this list will not be updated. HRUIDs in this list but not in Raven model will throw a warning to RavenErrors.txt.
    //Ideally,the naming convention of this would be that either(1) the variable name of this 2D matrix or (2) the long_name attribute of this matrix would be equivalent to the raven state variable tag,e.g.,SNOW,PONDED_WATER,SOIL[0]. Another option would be to add a raven_name attribute.

  string statefile="state_mods.nc"; //=Options.stateinfo_filename;

  if(Options.noisy) { cout<<"Opening state info file "<<statefile<<endl; }
  retval = nc_open(statefile.c_str(),NC_NOWRITE,&ncid);          HandleNetCDFErrors(retval);

  retval = nc_inq_varid(ncid,"HRUIDs",&HRUvecID);                HandleNetCDFErrors(retval);
  if(retval==NC_ENOTVAR) {
    string warning="Variable HRUIDs not found in NetCDF state file "+statefile;
    ExitGracefully(warning.c_str(),BAD_DATA); retval=0;
  }
  HandleNetCDFErrors(retval);

  //Get size of HRUID vector
  nHRUsLocal= 0;
  int dimid=0; //how to get this?? "stations" is ID?
  size_t lenp[]={0};
  retval= nc_inq_dimlen(ncid,dimid,lenp);
  if(retval!=NC_NOERR) {
    string warning="Failure determining size of HRU ID vector in "+statefile;
    ExitGracefully(warning.c_str(),BAD_DATA); retval=0;
  }
  nHRUsLocal=lenp[0];

  //parse vector of HRU IDs
  size_t start[] ={0};
  size_t count[] ={nHRUsLocal};
  int *HRUIDs=new int [nHRUsLocal];
  retval=nc_get_vara_int(ncid,HRUvecID,start,count,HRUIDs);

  retval = nc_close(ncid);         HandleNetCDFErrors(retval);

  //Parse array of states 

  ExitGracefully("TESTING FEWS STATE FILE",STUB);

  return true;
#else
  return false;
  ExitGracefully("ParseRunInfoFile: Deltares FEWS Runinfo file can only be read using NetCDF version of Raven",BAD_DATA_WARN);
#endif
}