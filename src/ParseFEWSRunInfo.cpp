/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2025 the Raven Development Team
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Model.h"
#include "SoilAndLandClasses.h"
#include "EnKF.h"

void GetNetCDFTimeInfo(const int ncid, int &time_dimid,int &ntime, int &start_time_index, const optStruct &Options);
void GetNetCDFStationArray(const int ncid, const string filename,int &stat_dimid,int &stat_varid, long *&aStations,  string *&aStat_string,int &nStations);

//////////////////////////////////////////////////////////////////
/// \brief Parses Deltares FEWS RunInfo file
///
/// \param *&pModel [out] Reference to model object
/// \param &Options [out] Global model options information
/// \return True if operation is successful
//
bool ParseNetCDFRunInfoFile(CModel *&pModel, optStruct &Options, bool runname_overridden, bool mode_overridden)
{
  if (Options.runinfo_filename==""){return true;}

#ifdef _RVNETCDF_
  int ncid;          //NetCDF file id
  int retval;        //return value of NetCDF routines
  size_t att_len;    //character string length

  if(Options.noisy) { cout<<"Opening FEWS runinfo file "<<Options.runinfo_filename<<endl; }
  retval = nc_open(Options.runinfo_filename.c_str(),NC_NOWRITE,&ncid);
  if(retval != 0) {
    string warn = "ParseNetCDFRunInfoFile: :FEWSRunInfoFile command: Cannot find file "+ Options.runinfo_filename;
    ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    return false;
  }
  HandleNetCDFErrors(retval);

  // Ingest start time=============================================================================
  int varid_startt;
  double start_time;
  retval = nc_inq_varid   (ncid,"start_time",&varid_startt);
  if(retval != NC_ENOTVAR)
  {
    HandleNetCDFErrors(retval);
    retval = nc_inq_attlen (ncid,varid_startt,"units",&att_len);                  HandleNetCDFErrors(retval);
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
  if(retval != NC_ENOTVAR)
  {
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

  // Ingest assimilation start time=============================================================================
  int varid_asst;
  double assim_time;
  double assim_day;
  int    assim_year;
  retval = nc_inq_varid(ncid,"assimilation_time",&varid_asst);
  if(retval != NC_ENOTVAR)
  {
    HandleNetCDFErrors(retval);
    retval = nc_inq_attlen(ncid,varid_asst,"units",&att_len);                 HandleNetCDFErrors(retval);
    char *unit_t=new char[att_len+1];
    retval = nc_get_att_text(ncid,varid_asst,"units",unit_t);                 HandleNetCDFErrors(retval);// read attribute text
    unit_t[att_len] = '\0';// add string determining character

    if(!IsValidNetCDFTimeString(unit_t))
    {
      cout<<"time unit string: "<<unit_t<<endl;
      ExitGracefully("ParseRunInfoFile: end_time time unit string is not in the format '[days/hours/...] since YYYY-MM-DD HH:MM:SS +0000' !",BAD_DATA);
    }

    retval=nc_get_var_double(ncid,varid_asst,&assim_time);
    GetJulianDateFromNetCDFTime(unit_t,Options.calendar,assim_time,assim_day,assim_year);

    Options.assimilation_start=TimeDifference(Options.julian_start_day,Options.julian_start_year,assim_day,assim_year,Options.calendar); //model time

    if(Options.noisy) { cout<<"ParseRunInfoFile: read variable assimilation_time from NetCDF: day: "<<assim_day<<" yr: "<< assim_year<<endl; }

    delete[] unit_t;
  }

  // Ingest working directory=====================================================================
  /*int varid_workdir;
  retval = nc_inq_varid(ncid,"work_dir",&varid_workdir);
  if(retval!= NC_ENOTVAR)
  {
    HandleNetCDFErrors(retval);

    char *workdir=new char[255];
    retval=nc_get_var(ncid,varid_workdir,workdir);

    Options.output_dir=to_string(workdir)+"/";
    if(Options.noisy) { cout<<"ParseRunInfoFile: read variable work_dir from NetCDF: "<<Options.output_dir<<endl; }
    delete [] workdir;
  }*/

  // Ingest properties      =====================================================================
  int varid_props;
  retval = nc_inq_varid(ncid,"properties",&varid_props);
  if(retval != NC_ENOTVAR)
  {
    HandleNetCDFErrors(retval);

    // RunName
    retval = nc_inq_attlen(ncid,varid_props,"RunName",&att_len);
    if(retval !=NC_ENOTATT)
    {
      HandleNetCDFErrors(retval);
      char *runname=new char[att_len+1];
      retval = nc_get_att_text(ncid,varid_props,"RunName",runname);                   HandleNetCDFErrors(retval);// read attribute text
      runname[att_len] = '\0';// add string determining character
      if(!runname_overridden) {
        Options.run_name=runname;
      }
      else {
        WriteWarning("ParseRunInfoFile: when run_name is specified from command line, it cannot be overridden in the runinfo file. RunName attribute ignored.",Options.noisy);
      }
      if(Options.noisy) { cout<<"ParseRunInfoFile: read properties:RunName from NetCDF: "<<Options.run_name<<endl; }
      delete[] runname;
    }

    // SuppressWarnings
    retval = nc_inq_attlen(ncid,varid_props,"BlockRavenWarnings",&att_len);
    if(retval !=NC_ENOTATT)
    {
      HandleNetCDFErrors(retval);
      char* boolean=new char[att_len+1];
      retval = nc_get_att_text(ncid,varid_props,"BlockRavenWarnings",boolean);       HandleNetCDFErrors(retval);// read attribute text
      boolean[att_len] = '\0';// add string determining character

      g_suppress_warnings=(!strcmp(boolean,"true"));
      if(Options.noisy) { cout<<"ParseRunInfoFile: read properties:BlockRavenWarnings from NetCDF: "<<g_suppress_warnings<<endl; }
      delete[] boolean;
    }

    //Block custom output
    retval = nc_inq_attlen(ncid,varid_props,"BlockRavenCustomOutput",&att_len);
    if(retval !=NC_ENOTATT)
    {
      HandleNetCDFErrors(retval);
      char* boolean=new char[att_len+1];
      retval = nc_get_att_text(ncid,varid_props,"BlockRavenCustomOutput",boolean);       HandleNetCDFErrors(retval);// read attribute text
      boolean[att_len] = '\0';// add string determining character

      if (!strcmp(boolean,"true")) {
        pModel->DeleteCustomOutputs();
      }
      if(Options.noisy) { cout<<"ParseRunInfoFile: read properties:BlockRavenCustomOutput from NetCDF: "<<(!strcmp(boolean,"true")) <<endl; }
      delete[] boolean;
    }

    //NoisyMode
    retval = nc_inq_attlen(ncid, varid_props, "NoisyMode", &att_len);
    if (retval != NC_ENOTATT)
    {
      HandleNetCDFErrors(retval);
      char* boolean = new char[att_len + 1];
      retval = nc_get_att_text(ncid, varid_props, "NoisyMode", boolean);       HandleNetCDFErrors(retval);// read attribute text
      boolean[att_len] = '\0';// add string determining character

      Options.noisy = (!strcmp(boolean,"true"));

      if (Options.noisy) { cout << "ParseRunInfoFile: read properties:NoisyMode from NetCDF: " << (!strcmp(boolean,"true"))  << endl; }
      delete[] boolean;
    }

    //SilentMode
    retval = nc_inq_attlen(ncid, varid_props, "SilentMode", &att_len);
    if (retval != NC_ENOTATT)
    {
      HandleNetCDFErrors(retval);
      char* boolean = new char[att_len + 1];
      retval = nc_get_att_text(ncid, varid_props, "SilentMode", boolean);       HandleNetCDFErrors(retval);// read attribute text
      boolean[att_len] = '\0';// add string determining character

      Options.silent = (!strcmp(boolean,"true"));
      if (Options.silent) { Options.noisy = false; }

      if (Options.noisy) { cout << "ParseRunInfoFile: read properties:SilentMode from NetCDF: " << (!strcmp(boolean,"true"))  << endl; }
      delete[] boolean;
    }

    // Mode
    retval = nc_inq_attlen(ncid, varid_props, "Mode", &att_len);
    if(retval != NC_ENOTATT)
    {
      HandleNetCDFErrors(retval);
      char* mode = new char[att_len];
      retval = nc_get_att_text(ncid,varid_props,"Mode",mode);       HandleNetCDFErrors(retval);// read attribute text

      Options.run_mode=mode[0];//only uses first character

      if(Options.noisy) { cout << "ParseRunInfoFile: read properties:Mode from NetCDF: " << Options.run_mode<< endl; }
      delete[] mode;
    }

    // EnKFMode
    retval = nc_inq_attlen(ncid, varid_props, "EnKFMode", &att_len);
    if(retval != NC_ENOTATT)
    {
      HandleNetCDFErrors(retval);
      char* modestr = new char[att_len];
      retval = nc_get_att_text(ncid,varid_props,"EnKFMode",modestr);       HandleNetCDFErrors(retval);// read attribute text
      CEnsemble *pEnsemble=pModel->GetEnsemble();
      EnKF_mode mode=ENKF_SPINUP;
      if(pEnsemble->GetType()==ENSEMBLE_ENKF) {
        CEnKFEnsemble* pEnKF=((CEnKFEnsemble*)(pEnsemble));

        if      (!strcmp(modestr,"ENKF_SPINUP"       )){mode=ENKF_SPINUP;}
        else if (!strcmp(modestr,"ENKF_CLOSED_LOOP"  )){mode=ENKF_CLOSED_LOOP;}
        else if (!strcmp(modestr,"ENKF_OPEN_LOOP"    )){mode=ENKF_OPEN_LOOP;}
        else if (!strcmp(modestr,"ENKF_FORECAST"     )){mode=ENKF_FORECAST;}
        else if (!strcmp(modestr,"ENKF_OPEN_FORECAST")){mode=ENKF_OPEN_FORECAST;}
        else {
          ExitGracefully("ParseEnsembleFile: EnKFMode-  invalid mode specified",BAD_DATA_WARN);
        }
        pEnKF->SetEnKFMode(mode);
      }

      if(Options.noisy) { cout << "ParseRunInfoFile: read properties:EnKFMode from NetCDF: " << mode<< endl; }
      delete[] modestr;
    }

    // AssimilateStreamflow
    retval = nc_inq_attlen(ncid,varid_props,"AssimilateStreamflow",&att_len);
    if(retval != NC_ENOTATT)
    {
      HandleNetCDFErrors(retval);
      char* boolean = new char[att_len + 1];
      retval = nc_get_att_text(ncid,varid_props,"AssimilateStreamflow",boolean);       HandleNetCDFErrors(retval);// read attribute text
      boolean[att_len] = '\0';// add string determining character

      Options.assimilate_flow = (!strcmp(boolean,"true"));

      if(Options.noisy) { cout << "ParseRunInfoFile: read properties:AssimilateStreamflow from NetCDF: " << (!strcmp(boolean,"true"))  << endl; }
      delete[] boolean;
    }

    // AssimilateReservoirStage
    retval = nc_inq_attlen(ncid,varid_props,"AssimilateReservoirStage",&att_len);
    if(retval != NC_ENOTATT)
    {
      HandleNetCDFErrors(retval);
      char* boolean = new char[att_len + 1];
      retval = nc_get_att_text(ncid,varid_props,"AssimilateReservoirStage",boolean);       HandleNetCDFErrors(retval);// read attribute text
      boolean[att_len] = '\0';// add string determining character

      Options.assimilate_stage = (!strcmp(boolean,"true"));

      if(Options.noisy) { cout << "ParseRunInfoFile: read properties:AssimilateReservoirStage from NetCDF: " << (!strcmp(boolean,"true"))  << endl; }
      delete[] boolean;
    }


    // Other attributes with info embedded in attribute name
    int nAttributes;
    char att_name[NC_MAX_NAME];
    double att_value;
    retval = nc_inq_varnatts(ncid,varid_props,&nAttributes);
    if (retval==NC_NOERR){
      for (int i = 0; i < nAttributes; i++) {
        retval=nc_inq_attname(ncid,varid_props,i,att_name);
        if (retval==NC_NOERR){
          string att_name_s=to_string(att_name);

          // NamedConstant(s) ----------------------------------------------------
          if (att_name_s.substr(0, 14) == "NamedConstant_")
          {
            string name=att_name_s.substr(14,att_name_s.length());

            retval = nc_inq_attlen(ncid,varid_props,att_name,&att_len);         HandleNetCDFErrors(retval);
            char* my_value = new char[att_len + 1];
            retval = nc_get_att_text(ncid,varid_props,att_name,my_value);       HandleNetCDFErrors(retval);// read attribute text
            my_value[att_len] = '\0';// add string determining character
            att_value=s_to_d(my_value);
            delete [] my_value;

            //cout<<"RUN INFO NAMED CONSTANT = "<<name<<" value = "<<att_value <<endl;
            if (Options.management_optimization){
              pModel->GetManagementOptimizer()->AddUserConstant(name,att_value);
            }
          }
        }
      }
    }
  }
  else {
    WriteWarning("ParseNetCDFRunInfoFile: Properties variable not found",Options.noisy);
  }

  retval = nc_close(ncid);         HandleNetCDFErrors(retval);

  return true;
#else
  ExitGracefully("ParseRunInfoFile: Deltares FEWS Runinfo file can only be read using NetCDF version of Raven",BAD_DATA_WARN);
  return false;
#endif
}

//////////////////////////////////////////////////////////////////
/// \brief Parses Deltares FEWS State info file
///
/// \param *&pModel [out] Reference to model object
/// \param &Options [out] Global model options information
/// \return True if operation is successful
//
// The state info file includes:
// 1) 'time' and 'stations' dimensions
// 2)	A vector named "station_id" of size[1:stations] which includes HRU IDs which reference the HRUID used in the Raven .rvh file.
//  The index i of each entry will be the same as that used by the matrix in part 3.
// 3)	One or more matrices of size[1:time,1:stations] which includes the state variable being updated,using the
//  same indexing[i=1:stations] as the vector of HRUIDs. Blank values will not be updated. Raven model HRUs with IDs
//  not in this list will not be updated. HRUIDs in this list but not in Raven model will throw a warning to Raven_errors.txt.
//  The naming convention of this attribute is equivalent to the raven state variable tag,e.g.,SNOW,PONDED_WATER,SOIL[0].
//  The state variable corresponding to the model start time will be used for initialization; all other values in the time vector are ignored
//  Note : this ONLY looks for state variables that are in model, all other state variable arrays will be ignored
//
bool ParseNetCDFStateFile(CModel *&pModel,const optStruct &Options)
{
  if(Options.stateinfo_filename=="") { return true; }

#ifdef _RVNETCDF_
  int     ncid;                // NetCDF file id
  int     retval;              // return value of NetCDF routines
  int     ntime    = 0;        // size of time vector
  size_t  att_len;             // length of the attribute's text
  nc_type att_type;            // type of attribute
  int     time_dimid;          // dimension id of time variable
  int     start_time_index=0;  // index of NetCDF time array which corresponds to Raven start time

  // Open file
  //====================================================================
  string statefile=Options.stateinfo_filename;

  if(Options.noisy) { cout<<"Opening state info file "<<statefile<<endl; }
  retval = nc_open(statefile.c_str(),NC_NOWRITE,&ncid);
  if(retval != 0) {
    string warn = "ParseNetCDFStateFile: :FEWSStateInfoFile command: Cannot find file "+ statefile;
    ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    return false;
  }
  HandleNetCDFErrors(retval);

  // Get information from time vector
  //====================================================================
  GetNetCDFTimeInfo(ncid,time_dimid,ntime,start_time_index,Options);

  if ((start_time_index < 0) || (start_time_index >= ntime))
  {
    time_struct tt;
    JulianConvert(0, Options.julian_start_day, Options.julian_start_year, Options.calendar, tt);
    string warn = "ParseNetCDFStateFile: Model start time ("+tt.date_string+" "+ DecDaysToHours(Options.julian_start_day)+") does not occur during time window provided in NetCDF state update file. State update file will be ignored.";
    WriteWarning(warn,Options.noisy);
    return true;
    //ExitGracefully(warn.c_str(), BAD_DATA);
    //return false;
  }

  // Get vector of HRUIDs
  //====================================================================
  int   HRUvecID;         //attribute ID of HRU vector
  int   HRUdimID;         // nstations dimension ID
  long * HRUIDs    =NULL;  // vector of HRUIDs //\todo[funct]: support long long int
  string *junk=NULL;
  int   nHRUsLocal=0;     // size of HRU Vector

  GetNetCDFStationArray(ncid, statefile,HRUdimID,HRUvecID, HRUIDs, junk, nHRUsLocal);

  if(Options.noisy) {
    cout << "[STATEFILE]: " << " station id found:" << HRUvecID << endl;
    cout << "[STATEFILE]: " << " station dimid found:" << HRUdimID << endl;
    cout << "[STATEFILE]: " << " time vector size:" << ntime << endl;
    cout << "[STATEFILE]: " << " nHRUs=" << nHRUsLocal << endl;
  }

  //check that valid (integer) HRUIDs exist in Raven Model
  // ===============================================================================
  for (int k = 0; k < nHRUsLocal; k++) {
    if ((HRUIDs[k]!=DOESNT_EXIST) && (pModel->GetHRUByID((long long int)(HRUIDs[k])) == NULL)) {
      ExitGracefully("ParseNetCDFStateFile: invalid HRU ID found in FEWS state update file. ", BAD_DATA);
    }
  }

  //Parse arrays of states   (look for all state variables in model)
  // ===============================================================================
  string* aSV_names = new string[pModel->GetNumStateVars()];
  for (int i = 0; i < pModel->GetNumStateVars(); i++) {
    aSV_names[i] = pModel->GetStateVarInfo()->SVTypeToString(pModel->GetStateVarType(i), pModel->GetStateVarLayer(i));
  }

  int       SVID;
  double*   state_vec;

  size_t    nc_start [] = { 0, 0 };
  size_t    nc_length[] = { (size_t)(ntime), (size_t)(nHRUsLocal) };
  ptrdiff_t nc_stride[] = { 1, 1 };

  state_vec = new double[nHRUsLocal*ntime];

  for (int i = 0; i < pModel->GetNumStateVars(); i++)
  {
    //look to see if state variable name is in file
    retval = nc_inq_varid(ncid, aSV_names[i].c_str(), &SVID);
    if (retval != NC_ENOTVAR)
    {
      HandleNetCDFErrors(retval);

      // verify dimensions are ntimes*nHRUsLocal, else throw error
      // -------------------------------
      int      sv_dimid[] = { 0 , 0 };
      retval = nc_inq_var(ncid, SVID, NULL, NULL, NULL, sv_dimid, NULL);                    HandleNetCDFErrors(retval);
      if ((sv_dimid[0] != time_dimid) || (sv_dimid[1] != HRUdimID)) {
        ExitGracefully("ParseNetCDFStateFile: dimensions of state variable array are not (time,nstations)", BAD_DATA);
      }

      // find "_FillValue" of data, if present
      // -------------------------------
      double fillval = NETCDF_BLANK_VALUE; //Default
      retval = nc_inq_att(ncid, SVID, "_FillValue", &att_type, &att_len);
      if (retval != NC_ENOTATT) {
        HandleNetCDFErrors(retval);
        retval = nc_get_att_double(ncid, SVID, "_FillValue", &fillval);                     HandleNetCDFErrors(retval);
      }

      // read array of doubles, update model states
      // -------------------------------
      retval = nc_get_vars_double(ncid, SVID, nc_start, nc_length, nc_stride, state_vec);   HandleNetCDFErrors(retval);

      for (int k = 0; k < nHRUsLocal; k++)
      {
        int index = start_time_index * nHRUsLocal + k;
        if ((state_vec[index] != NETCDF_BLANK_VALUE) && (state_vec[index] != fillval) && (state_vec[index]!=RAV_BLANK_DATA))
        {
          //cout << "[STATEFILE]: OVERRIDING:" <<aSV_names[i]<<" ["<< k<<"]: "<< state_vec[index] << endl;
          if (HRUIDs[k] != DOESNT_EXIST) {
            pModel->GetHRUByID(HRUIDs[k])->SetStateVarValue(i, state_vec[index]);
          }
        }
        else {
          //cout << "[STATEFILE]: " << aSV_names[i] << " [" << k << "]: NULL"  << endl;
        }
      }

    }
    else {
      //state variable not found in NetCDF file; do nothing
    }
  }
  delete[] state_vec;
  delete[] aSV_names;
  delete[] HRUIDs;
  delete[] junk;

  //close file
  //====================================================================
  retval = nc_close(ncid);         HandleNetCDFErrors(retval);

  // check quality of initial state variables
  //======================================================================
  CHydroUnit* pHRU;
  double *v = new double[pModel->GetNumStateVars()];
  for (int k = 0; k < pModel->GetNumHRUs(); k++)
  {
    pHRU = pModel->GetHydroUnit(k);
    for (int i = 0; i < pModel->GetNumStateVars(); i++)
    {
      v[i] = pHRU->GetStateVarValue(i);
    }
    for (int i = 0; i < pModel->GetNumStateVars(); i++)
    {
      double maxv = max(pHRU->GetStateVarMax(i, v, Options, true), 0.0); //ignores all variable maximum thresholds that are dependent upon model state
      if ((v[i] - maxv) > PRETTY_SMALL)// check for capacity
      {
        string name = pModel->GetStateVarInfo()->GetStateVarLongName(pModel->GetStateVarType(i), pModel->GetStateVarLayer(i),pModel->GetTransportModel());
        string warn = "maximum state variable limit exceeded in initial conditions for " + name + " (in HRU " + to_string(pHRU->GetHRUID()) + ") in FEWS state update file";
        WriteWarning(warn, Options.noisy);

        if (!Options.keepUBCWMbugs) {
          pHRU->SetStateVarValue(i, maxv);
        }
      }
    }
  }
  delete[] v;
  return true;
#else
  ExitGracefully("ParseNetCDFStateFile: Deltares FEWS State Update file can only be read using NetCDF version of Raven",BAD_DATA_WARN);
  return false;
#endif
}

//////////////////////////////////////////////////////////////////
/// \brief Parses Deltares FEWS Flow State update file
///
/// \param *&pModel [out] Reference to model object
/// \param &Options [out] Global model options information
/// \return True if operation is successful
//
// The flow state update file includes:
// 1) 'time' and 'stations' dimensions
// 2)	A vector named "station_id" of size[1:stations] which includes Subbasin IDs which reference the Subbasin used in the Raven .rvh file.
//  The index i of each entry will be the same as that used by the matrix in part 3.
// 3)	One or more matrices of size[1:time,1:stations] which includes the flow/stage state variable being initialized,using the
//  same indexing[i=1:stations] as the vector of SBIDs. Blank values will not be updated. Raven model subbasins with IDs
//  not in this list will not be updated. SBIDs in this list but not in Raven model will throw a warning to Raven_errors.txt.
// .The state variable corresponding to the model start time will be used for initialization of flow/stage; all other values in the time vector are ignored
//  The naming convention of this attribute is one of: {QOUT, STAGE,  }
//
bool ParseNetCDFFlowStateFile(CModel*& pModel,const optStruct& Options)
{
  if(Options.flowinfo_filename=="") { return true; }

#ifdef _RVNETCDF_
  int     ncid;                // NetCDF file id
  int     retval;              // return value of NetCDF routines
  int     ntime    = 0;        // size of time vector
  size_t  att_len;             // length of the attribute's text
  nc_type att_type;            // type of attribute
  int     time_dimid;          // dimension id of time variable
  int     start_time_index=0;  // index of NetCDF time array which corresponds to Raven start time

  // Open file
  //====================================================================
  string statefile=Options.flowinfo_filename;

  if(Options.noisy) { cout<<"Opening flow state info file "<<statefile<<endl; }
  retval = nc_open(statefile.c_str(),NC_NOWRITE,&ncid);
  if(retval != 0) {
    string warn = "ParseNetCDFStateFile: :FEWSBasinStateInfoFile command: Cannot find file "+ statefile;
    ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    return false;
  }
  HandleNetCDFErrors(retval);

  // Get information from time vector
  //====================================================================
  GetNetCDFTimeInfo(ncid,time_dimid,ntime,start_time_index,Options);

  if ((start_time_index < 0) || (start_time_index >= ntime))
  {
    time_struct tt;
    JulianConvert(0, Options.julian_start_day, Options.julian_start_year, Options.calendar, tt);
    string warn = "ParseNetCDFFlowStateFile: Model start time ("+tt.date_string+" "+ DecDaysToHours(Options.julian_start_day)+") does not occur during time window provided in NetCDF state update file";
    ExitGracefully(warn.c_str(), BAD_DATA);
    return false;
  }

  // Get vector of Subbasin IDs
  //====================================================================
  int   SBvecID;         //attribute ID of SB vector
  int   SBdimID;         // nstations dimension ID
  long* SBIDs    =NULL;  // vector of Subbasin IDs (allocated in GetNetCDFStationArray)
  string *junk  =NULL;
  int   nSBsLocal=0;     // size of subbasin Vector

  GetNetCDFStationArray(ncid, statefile,SBdimID,SBvecID, SBIDs,junk, nSBsLocal);

  if(Options.noisy) {
    cout << "[FLOWFILE]: " << " station id found:" << SBvecID << endl;
    cout << "[FLOWFILE]: " << " station dimid found:" << SBdimID << endl;
    cout << "[FLOWFILE]: " << " time vector size:" << ntime << endl;
    cout << "[FLOWFILE]: " << " nSBs=" << nSBsLocal << endl;
    cout << "[FLOWFILE]: Check"<<endl;
  }

  //check that valid (integer) SBIDs exist in Raven Model
  // ===============================================================================
  for (int p = 0; p < nSBsLocal; p++) {
    if ((SBIDs[p]!=DOESNT_EXIST) && (pModel->GetSubBasinByID((int)(SBIDs[p])) == NULL)) {
      ExitGracefully("ParseNetCDFStateFile: invalid SB ID found in FEWS state update file. ", BAD_DATA);
    }
  }

  //Parse arrays of states   (look for all state variables in model)
  // ===============================================================================
  string* aSV_names = new string[3];
  aSV_names[0]="BASIN_OUTFLOW";
  aSV_names[1]="RESERVOIR_STAGE";
  aSV_names[2]="BASIN_INFLOW";

  int       SVID;
  double*   state_vec;

  size_t    nc_start [] = { 0, 0 };
  size_t    nc_length[] = { (size_t)(ntime), (size_t)(nSBsLocal) };
  ptrdiff_t nc_stride[] = { 1, 1 };

  state_vec = new double[nSBsLocal*ntime];

  for (int i = 0; i < 2; i++)
  {
    //look to see if state variable name is in file
    retval = nc_inq_varid(ncid, aSV_names[i].c_str(), &SVID);
    if (retval != NC_ENOTVAR)
    {
      HandleNetCDFErrors(retval);

      // verify dimensions are ntimes*nSBsLocal, else throw error
      // -------------------------------
      int      sv_dimid[] = { 0 , 0 };
      retval = nc_inq_var(ncid, SVID, NULL, NULL, NULL, sv_dimid, NULL);                    HandleNetCDFErrors(retval);
      if ((sv_dimid[0] != time_dimid) || (sv_dimid[1] != SBdimID)) {
        ExitGracefully("ParseNetCDFFlowStateFile: dimensions of state variable array are not (time,nstations)", BAD_DATA);
      }

      // find "_FillValue" of data, if present
      // -------------------------------
      double fillval = NETCDF_BLANK_VALUE; //Default
      retval = nc_inq_att(ncid, SVID, "_FillValue", &att_type, &att_len);
      if (retval != NC_ENOTATT) {
        HandleNetCDFErrors(retval);
        retval = nc_get_att_double(ncid, SVID, "_FillValue", &fillval);                     HandleNetCDFErrors(retval);
      }

      // read array of doubles, update model states
      // -------------------------------
      retval = nc_get_vars_double(ncid, SVID, nc_start, nc_length, nc_stride, state_vec);   HandleNetCDFErrors(retval);

      for (int p = 0; p < nSBsLocal; p++)
      {
        int index = start_time_index * nSBsLocal + p;
        if ((state_vec[index] != NETCDF_BLANK_VALUE) && (state_vec[index] != fillval))
        {
          //cout << "[FLOWFILE]: OVERRIDING" <<aSV_names[i]<<" ["<< p<<"]: "<< state_vec[index] << endl;
          if (SBIDs[p] != DOESNT_EXIST)
          {
            CSubBasin *pSB=pModel->GetSubBasinByID(SBIDs[p]);
            if      (i==0){//"BASIN_OUTFLOW"
              pSB->SetQout(state_vec[index]);
            }
            else if (i==1){//"RESERVOIR_STAGE"
              pSB->GetReservoir()->SetReservoirStage(state_vec[index],state_vec[index]);
            }
            else if (i==2){//"BASIN_INFLOW"
              int N=pSB->GetInflowHistorySize();
              double *Qin=new double [N];
              for(int i=0;i<N;i++) { Qin[i]=pSB->GetInflowHistory()[i]; }
              Qin[0]=state_vec[index]; //overwrite Qin for this time step
              pSB->SetQinHist(N,Qin);
              delete [] Qin;
            }
          }
        }
        else {
          //cout << "[FLOWFILE]: " << aSV_names[i] << " [" << p << "]: NULL"  << endl;
        }
      }

    }
    else {
      //state variable not found in NetCDF file; do nothing
    }
  }
  delete[] state_vec;
  delete[] aSV_names;
  delete[]SBIDs;
  delete[]junk;

  //close file
  //====================================================================
  retval = nc_close(ncid);         HandleNetCDFErrors(retval);

  return true;
#else
  ExitGracefully("ParseNetCDFFlowStateFile: Deltares FEWS State Update file can only be read using NetCDF version of Raven",BAD_DATA_WARN);
  return false;
#endif
}


//////////////////////////////////////////////////////////////////
/// \brief Parses Deltares FEWS Parameter update file
///
/// \param *&pModel [out] Reference to model object
/// \param &Options [out] Global model options information
/// \return True if operation is successful
//
// The parameter update file includes:
// 1)	a time vector of size ntime
// 2) Vectors of size (ntime) with parameter time series; these are named PARAM_in_CLASS where
//    PARAM is a valid parameter name and CLASS is a valid soil class, veg class, land use class, or GLOBALS
//
bool ParseNetCDFParamFile(CModel*&pModel, const optStruct &Options)
{
  if (Options.paraminfo_filename == "") { return true; }

#ifdef _RVNETCDF_
  int     ncid;          // NetCDF file id
  int     retval;        // return value of NetCDF routines

  string paramfile=Options.paraminfo_filename;

  if(Options.noisy) { cout<<"Opening FEWS parameter info file "<<paramfile<<endl; }
  retval = nc_open(paramfile.c_str(),NC_NOWRITE,&ncid);
  if(retval != 0) {
    string warn = "ParseNetCDFParamFile: :FEWSParamInfoFile command: Cannot find file "+ paramfile;
    ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    return false;
  }
  HandleNetCDFErrors(retval);

  //Get information from time vector
  //====================================================================
  int ntime;
  int time_dimid;
  int start_time_index=0;
  GetNetCDFTimeInfo(ncid,time_dimid,ntime,start_time_index,Options);

  //cout<<"[PARAMFILE]: "<<" start time index="<<start_time_index<<endl;
  //start_time_index = 1;  //TMP DEBUG!!

  if ((start_time_index < 0) || (start_time_index >= ntime))
  {
    time_struct tt;
    JulianConvert(0, Options.julian_start_day, Options.julian_start_year, Options.calendar, tt);
    string warn = "ParseNetCDFStateFile: Model start time ("+tt.date_string+" "+ DecDaysToHours(Options.julian_start_day)+") does not occur during time window provided in NetCDF state update file";
    ExitGracefully(warn.c_str(), BAD_DATA_WARN);
    return false;
  }

  // Sift through all attributes that include "_in_"
  //====================================================================
  int nvars[]={0};
  int nVars=0;
  retval=nc_inq_nvars(ncid,nvars);                    HandleNetCDFErrors(retval);
  nVars=nvars[0];

  int *varids=new int[nVars];
  retval=nc_inq_varids(ncid,nvars,varids);            HandleNetCDFErrors(retval);

  string namestr;
  string param_str;
  string class_str;
  char *varname=new char [256];
  for (int j = 0; j < nVars; j++)
  {
    retval=nc_inq_varname(ncid,varids[j],varname);    HandleNetCDFErrors(retval);
    //cout<<"[PARAMFILE]: "<<" var="<<varname<<" ";
    namestr=varname;
    if (namestr.find("_in_")!=string::npos)
    {
      int i=(int)namestr.find("_in_");
      param_str=namestr.substr(0,i);
      class_str=namestr.substr(i+4,1000);

      //figure out if soil, veg, lult, subbasin, or global parameter
      class_type pclass=CLASS_UNKNOWN;
      pclass=pModel->ParamNameToParamClass(param_str,class_str);

      if (pclass!=CLASS_UNKNOWN)
      {
         // find "_FillValue" of attribute, if present
        // -------------------------------
        double fillval;
        size_t  att_len;       // length of the attribute's text
        nc_type att_type;      // type of attribute
        fillval = NETCDF_BLANK_VALUE; //Default
        retval = nc_inq_att(ncid, varids[j], "_FillValue", &att_type, &att_len);
        if (retval != NC_ENOTATT) {
          HandleNetCDFErrors(retval);
          retval = nc_get_att_double(ncid, varids[j], "_FillValue", &fillval);       HandleNetCDFErrors(retval);// read attribute value
        }

        // check for valid class name
        bool bad=false;

        if      (pclass==CLASS_SOIL      ){if (pModel->StringToSoilClass(class_str)==NULL){bad=true;}}
        else if (pclass==CLASS_VEGETATION){if (pModel->StringToVegClass (class_str)==NULL){bad=true;}}
        else if (pclass==CLASS_LANDUSE   ){if (pModel->StringToLUClass  (class_str)==NULL){bad=true;}}

        if (bad){
          string warn="ParseNetCDFParamFile:: unrecognized soil/veg/lult class found ("+class_str+") in FEWS parameter update file";
          WriteWarning(warn.c_str(),Options.noisy);
        }
        else{
          // get parameter vector from NetCDF
          // ===============================================================================
          size_t    start[]   = {0,0};
          size_t    count[]   = {(size_t)(ntime),1};
          ptrdiff_t stridep[] = {1,1};
          float     *ip=new float [ntime];

          nc_get_vars_float(ncid,varids[j],start,count,stridep,ip);    HandleNetCDFErrors(retval);

          double pval=(double)ip[start_time_index];
          /*cout<<" * ";
          cout<<" |"<<param_str<<"| ";
          cout<<" |"<<class_str<<"| ";
          cout<<" CLASS: "<<pclass<<" ";
          cout<<" val: "<<pval<<endl;*/
          delete [] ip;

          if ((pval!=NETCDF_BLANK_VALUE) && (pval!=fillval))
          {

            if (Options.noisy) {cout<<"  ParseNetCDFParamFile: Updating "<<param_str<<" in "<<class_str<<" value: "<<pval<<endl;}
            pModel->UpdateParameter(pclass,param_str,class_str,pval);// update parameter value
          }
        }
      }
      else {
        cout<<" unrecognized: "<<param_str<<" ";
        string warn="ParseNetCDFParamFile:: unrecognized parameter found ("+param_str+") in FEWS parameter update file";
        WriteWarning(warn.c_str(),Options.noisy);
      }
    }
    //cout<<endl;
  }

  retval = nc_close(ncid);   HandleNetCDFErrors(retval);

#endif
  return true;
}



//////////////////////////////////////////////////////////////////
/// \brief Extracts time array info and timestep index within NetCDF time() array corresponding to Raven model start time
///
/// \param ncid [in] NetCDF file id
/// \param time_dimid[out] NetCDF dimension idea for time dimension
/// \param ntime[out] size of NetCDF time array
/// \param start_time_index[out] index of NetCDF time array corresponding to Raven model start time
/// \param &Options [in] Global model options information
//
void GetNetCDFTimeInfo(const int ncid, int &time_dimid, int &ntime, int &start_time_index, const optStruct& Options)
{
#ifdef _RVNETCDF_
  size_t att_len;
  int    retval;
  string unit_t;
  double ref_day, delta_t;
  int    ref_year;
  char*  unitstr;

  //get number of timesteps, time dimension
  int      time_dim[]   = { 0 };
  size_t   timedimlen[] = { 0 };
  retval = nc_inq_dimid (ncid, "time", time_dim);           HandleNetCDFErrors(retval);
  retval = nc_inq_dimlen(ncid, time_dim[0], timedimlen);    HandleNetCDFErrors(retval);
  ntime = (int)(timedimlen[0]);
  time_dimid=time_dim[0];

  retval  = nc_inq_attlen(ncid, time_dimid, "units", &att_len);       HandleNetCDFErrors(retval);
  unitstr = new char[att_len + 1];
  retval  = nc_get_att_text(ncid, time_dimid, "units", unitstr);      HandleNetCDFErrors(retval);
  unitstr[att_len] = '\0';
  unit_t  = unitstr;

  float *timearr = new float[ntime];
  retval = nc_get_var_float(ncid, time_dimid, timearr); HandleNetCDFErrors(retval);

  GetJulianDateFromNetCDFTime(unit_t, Options.calendar, timearr[0], ref_day, ref_year);

  delta_t=TimeDifference(ref_day, ref_year, Options.julian_start_day, Options.julian_start_year, Options.calendar); //gets start day with reference to reference day (start-ref) - should be postive
  start_time_index = (int)(rvn_round(delta_t / Options.timestep));

  //TMP DEBUG - only for debugging below
  /*time_struct tt;
  JulianConvert(0, Options.julian_start_day, Options.julian_start_year, Options.calendar, tt);
  cout << "[FEWS FILE]: " << " TIME ARRAY :" << endl;
  for (int i = 0; i < ntime; i++) {
  cout << "[FEWS FILE]: " << " time[" << i << "]: " << timearr[i] << endl;
  }
  cout << "[FEWS FILE]: " << " ref day, year: " << ref_day << ", " << ref_year << endl;
  cout << "[FEWS FILE]: " << " time starts "<< start_time_index <<" timesteps after model start ("<<tt.date_string<<")"<< endl;
  */ //END TMP DEBUG

  delete[] timearr;
  delete[] unitstr;
#endif
}

//////////////////////////////////////////////////////////////////
/// \brief Extracts station array info from Deltares FEWS-generated NetCDF file
///
/// \param ncid [in] NetCDF file id
/// \param filename [in] filename of NetCDF file
/// \param stat_dimid [out] dimension ID of 'stations' dimension
/// \param stat_varid [out] NetCDF variable id for variable 'station_id'
/// \param aStations [out] array of station IDs, as long (this vector is created here, and aStations should be deleted outside this routine)
/// \param nStations [out] size of aStations array
//
void GetNetCDFStationArray(const int ncid, const string filename,int &stat_dimid,int &stat_varid, long *&aStations, string *&aStat_string, int &nStations)
{
#ifdef _RVNETCDF_
  int retval;
  int charsize;

  aStations=NULL;

  // get varid of station_id
  //---------------------------------------------------
  retval = nc_inq_varid(ncid,"station_id",&stat_varid);
  if(retval==NC_ENOTVAR) {
    string warning="GetNetCDFStationArray: Variable 'station_id' not found in NetCDF forcing file "+filename;
    ExitGracefully(warning.c_str(),BAD_DATA); retval=0;
  }
  HandleNetCDFErrors(retval);

  // get character length dimension
  //---------------------------------------------------
  int      char_dimid[] = { 0 };
  size_t   chdimlens [] = { 0 };
  retval = nc_inq_dimid (ncid, "char_leng_id", char_dimid);   HandleNetCDFErrors(retval);
  retval = nc_inq_dimlen(ncid, char_dimid[0], chdimlens);     HandleNetCDFErrors(retval);
  charsize = (int)(chdimlens[0]); //usually 64

  //Get size of station vector, allocate memory
  //---------------------------------------------------
  int    station_dimid[] = { 0 };
  size_t       dimlens[] = { 0 };
  retval = nc_inq_dimid (ncid, "stations", station_dimid);    HandleNetCDFErrors(retval);
  retval = nc_inq_dimlen(ncid, station_dimid[0],dimlens);
  if(retval!=NC_NOERR) {
    string warning="GetNetCDFStationArray: Failure determining size of stations dimension in "+filename;
    ExitGracefully(warning.c_str(),BAD_DATA); retval=0;
  }
  HandleNetCDFErrors(retval);
  stat_dimid=station_dimid[0];
  nStations=(int)(dimlens[0]);
  aStations= new long [nStations];
  aStat_string = new string [nStations];

  // parse and store vector of station IDs
  // ===============================================================================
  size_t    start[]   = {0,0};
  size_t    count[]   = {(size_t)(nStations),(size_t)(charsize)};
  ptrdiff_t stridep[] = {1,1};
  char     *ID_chars;
  void     *ip;
  string    str;
  char      tmp;

  ip     = new char [nStations * charsize];
  retval = nc_get_vars(ncid, stat_varid, start, count, stridep, ip);   HandleNetCDFErrors(retval); //read string array
  ID_chars = (char*)(ip);

  //convert character array to long IDs, store in aStations[]
  for (int k = 0; k < nStations; k++)
  {
    str="";
	  aStat_string[k]="";
    for (int j = 0; j < charsize; j++) {
      tmp = ID_chars[k * charsize + j];
      if (isdigit(tmp)) {str += tmp;}
      aStat_string[k]+=tmp;
      //if (isdigit(tmp)) { cout << tmp << " |"; }
    }
    //cout << endl;
    aStations[k] = s_to_l(str.c_str());
    if (aStations[k] == 0) { aStations[k] = DOESNT_EXIST;}
    //cout << "[STATION_ID]: " << " ID["<<k<<"]="<< aStations[k] <<endl;
  }
  delete[] ip;
#endif
}
