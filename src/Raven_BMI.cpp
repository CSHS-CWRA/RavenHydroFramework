/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2023 the Raven Development Team
  ----------------------------------------------------------------*/
#include <time.h>
#include "RavenInclude.h"
#include "RavenMain.h"
#include "Raven_BMI.h"

const int GRID_SUBBASIN=1;
const int GRID_HRU     =0;

//////////////////////////////////////////////////////////////////
/// \brief RavenBMI class constructor and destructor
//
CRavenBMI::CRavenBMI()
{
  //not sure how this will work with global variable defined in RavenMain.
  pModel=NULL;
}

CRavenBMI::~CRavenBMI() {}


//////////////////////////////////////////////////////////////////
/// \brief Splits a string by whitespace into a vector of char*
///
/// \param line [in] string to be split
/// \return vector of char* containing the split string
std::vector<char *> CRavenBMI::SplitLine(std::string line)
{
  std::vector<char *> args;
  std::istringstream iss(line);
  std::string token;

  while(iss >> token)
  {
    char *arg = new char[token.size() + 1];
    copy(token.begin(), token.end(), arg);
    arg[token.size()] = '\0';
    args.push_back(arg);
  }
  args.push_back(0);
  return args;
}


//////////////////////////////////////////////////////////////////
/// @brief Sets values to Option struct based on command line arguments
///
/// @param config_key Argument key of the yaml config file
/// @param config_value Value of the argument
/// @return void - sets values in Options struct
void CRavenBMI::ProcessConfigFileArgument(std::string config_key, std::string config_value)
{
  if (config_key == "rvi_file") {
    Options.rvi_filename = config_value;
  } else if (config_key == "rvp_file") {
    Options.rvp_filename = config_value;
  } else if (config_key == "rvh_file") {
    Options.rvh_filename = config_value;
  } else if (config_key == "rvt_file") {
    Options.rvt_filename = config_value;
  } else if (config_key == "rvc_file") {
    Options.rvc_filename = config_value;
  } else if (config_key == "rvg_file") {
    Options.rvg_filename = config_value;
  } else if (config_key == "rve_file") {
    Options.rve_filename = config_value;
  } else if (config_key == "rvl_file") {
    Options.rvl_filename = config_value;
  } else if (config_key == "output_directory") {
    Options.output_dir = config_value;
  } else if (config_key == "main_output_directory") {
    Options.main_output_dir = config_value;
  } else {
    // TODO: include all remaining options
    cout << "WARNING: Unknown key '" << config_key << "' in Raven config file. Ignoring it." << endl;
  }
  return;
}


//////////////////////////////////////////////////////////////////
/// \brief Reads the BMI configuration file and checks if model options are valid
///
/// \param config_file [in] name of configuration file
/// \return void - sets values in Options struct
void CRavenBMI::ReadConfigFile(std::string config_file)
{
  ifstream CONFIG;
  int config_str_ini, config_str_end;    // used to parse the lines of the config file
  std::string config_key, config_value;  // used to parse the lines of the config file
  bool cli_args = false;                 // used to check if command line arguments were used
  std::vector<char *> args;              // used to parse the command line arguments
  char** argv;

  CONFIG.open(config_file);
  if (CONFIG.fail()) {
    cout << "Cannot find configuration file " << config_file << endl;
    return;
  }

  // read line by line and parse
  for( std::string line; getline( CONFIG, line ); )
  {
    // skip blank lines and comments
    if (line.empty() || line[0] == '#') {
      continue;
    }

    // parse line and update options struct accordingly
    config_str_end = -1;
    config_key = "";
    config_value = "";
    do {
        // find ini and end of substring
        config_str_ini = config_str_end + 1;
        config_str_end = (int)line.find(":", config_str_ini);

        // extract substring and trim leading and trailing whitespaces
        config_value = line.substr(config_str_ini, config_str_end - config_str_ini);
        config_value.erase(0, config_value.find_first_not_of(" \t"));
        config_value.erase(config_value.find_last_not_of(" \t") + 1);

        // set config_key if it is the first substring
        if (config_key == "") {
          config_key = config_value;
          config_value = "";
          continue;
        }

        // set options struct accordingly
        if (config_key == "cli_args") {
          cli_args = true;
          // "Raven.exe" is a dummy argument to mimic the command line arguments
          args = SplitLine("Raven.exe " + config_value);
          argv = args.data();
          ProcessExecutableArguments((int)(args.size())-1, argv, Options);
          break;
        } else {
          ProcessConfigFileArgument(config_key, config_value);
        }

    } while (config_str_end != -1);

    // if a "cli_args:" line was found then we can stop parsing
    if (cli_args) { break; }

  }
  CONFIG.close();

  return;
}


//////////////////////////////////////////////////////////////////
/// \brief Initialization called prior to model simulation
///
/// \param config_file [in] name of configuration file (unused)
//
void CRavenBMI::Initialize(std::string config_file)
{
  //NOTE: ENSEMBLE MODE NOT SUPPORTED WITH BMI

  ReadConfigFile(config_file);

  PrepareOutputdirectory(Options);

  for (int i=0;i<10;i++){g_debug_vars[i]=0;}

  ofstream WARNINGS;
  WARNINGS.open((Options.main_output_dir+"Raven_errors.txt").c_str());
  if (WARNINGS.fail()){
    ExitGracefully("Main::Unable to open Raven_errors.txt. Bad output directory specified?",RAVEN_OPEN_ERR);
  }
  WARNINGS.close();

  Options.in_bmi_mode  = true;  // flag needed to ignore some arguments of the .rvi file
  Options.rvt_filename = "";   // just a dummy filename to avoid errors
  Options.duration     = ALMOST_INF;  // "infinity": will run as long as "Update()" is called

  //Read input files, create model, set model options
  if (!ParseInputFiles(pModel, Options)){
    ExitGracefully("Main::Unable to read input file(s)",BAD_DATA);}

  CheckForErrorWarnings(true, pModel);

  pModel->Initialize                  (Options);
  ParseInitialConditions              (pModel, Options);
  pModel->CalculateInitialWaterStorage(Options);
  pModel->SummarizeToScreen           (Options);
  pModel->GetEnsemble()->Initialize   (pModel,Options);

  CheckForErrorWarnings(false, pModel);

  PrepareOutputdirectory(Options); //adds new output folders, if needed
  pModel->WriteOutputFileHeaders(Options);

  double t_start=0.0;
  t_start=pModel->GetEnsemble()->GetStartTime(0);

  //Write initial conditions-------------------------------------
  JulianConvert(t_start,Options.julian_start_day,Options.julian_start_year,Options.calendar,tt);
  pModel->RecalculateHRUDerivedParams(Options,tt);
  pModel->UpdateHRUForcingFunctions  (Options,tt);
  pModel->UpdateDiagnostics          (Options,tt);
  pModel->WriteMinorOutput           (Options,tt);

}
//////////////////////////////////////////////////////////////////
/// \brief run simulation for a single time step
///
//
void CRavenBMI::Update()
{
  pModel->UpdateTransientParams      (Options,tt);
  pModel->RecalculateHRUDerivedParams(Options,tt);
  pModel->UpdateHRUForcingFunctions  (Options,tt);
  pModel->UpdateDiagnostics          (Options,tt);
  pModel->PrepareAssimilation        (Options,tt);
  pModel->WriteSimpleOutput          (Options,tt);
 // pModel->GetEnsemble()->StartTimeStepOps(pModel,Options,tt,e);
  CallExternalScript                 (Options,tt);
  ParseLiveFile                      (pModel,Options,tt);

  MassEnergyBalance(pModel,Options,tt); //where the magic happens!

  pModel->IncrementCumulInput        (Options,tt);
  pModel->IncrementCumOutflow        (Options,tt);

  JulianConvert(tt.model_time+Options.timestep,Options.julian_start_day,Options.julian_start_year,Options.calendar,tt);//increments time structure

  pModel->WriteMinorOutput           (Options,tt);
  //pModel->WriteProgressOutput        (Options,clock()-t1,step,(int)ceil(Options.duration/Options.timestep));
  //pModel->GetEnsemble()->CloseTimeStepOps(pModel,Options,tt,e);

  //if ((Options.use_stopfile) && (CheckForStopfile(step,tt))) { break; }
  //step++;
}
//////////////////////////////////////////////////////////////////
/// \brief run simulation for multiple time steps, ending at time 'time'
///
/// \param time [in] model time at which simulation should stop (in days)
//
void CRavenBMI::UpdateUntil(double time)
{
  for(double t=tt.model_time; t<(time-TIME_CORRECTION); t+=Options.timestep)  // in [d]
  {
    Update();
  }
}
//////////////////////////////////////////////////////////////////
/// \brief called after simulation is done
//
void CRavenBMI::Finalize()
{
  pModel->UpdateDiagnostics (Options,tt);
  pModel->RunDiagnostics    (Options);
  pModel->WriteMajorOutput  (Options,tt,"solution",true);
  pModel->CloseOutputStreams();

  FinalizeGracefully("Successful Simulation", SIMULATION_DONE);
}

//------------------------------------------------------------------
// Model information functions.
//------------------------------------------------------------------

//////////////////////////////////////////////////////////////////
/// \brief returns name of BMI component
/// \return name of BMI component
//
std::string CRavenBMI::GetComponentName()
{
  return "Raven Hydrological Modelling Framework "+Options.version;
}
//////////////////////////////////////////////////////////////////
/// \brief returns number of accessible input datasets
/// \return number of accessible input datasets
//
int CRavenBMI::GetInputItemCount()
{
  return 2;
}
//////////////////////////////////////////////////////////////////
/// \brief returns array of accessible input dataset names
/// \return array of accessible input dataset names
//
std::vector<std::string> CRavenBMI::GetInputVarNames()
{
  vector<string> names;
  names.push_back("precipitation");
  names.push_back("temp_ave");
  return names;
}
//////////////////////////////////////////////////////////////////
/// \brief returns number of accessible output datasets
/// \return number of accessible output datasets
//
int CRavenBMI::GetOutputItemCount()
{
  return 3;
}
//////////////////////////////////////////////////////////////////
/// \brief returns array of accessible output dataset names
/// \return array of accessible output dataset names
//
std::vector<std::string> CRavenBMI::GetOutputVarNames()
{
  vector<string> names;
  names.push_back("streamflow");
  names.push_back("soil[0]");
  names.push_back("snow");
  return names;
}

//------------------------------------------------------------------
// Variable information functions
//------------------------------------------------------------------
//////////////////////////////////////////////////////////////////
/// \brief returns grid id of input or output variable
/// \param name [in] - name of  variable
/// \return grid id of input or output variable
//
int CRavenBMI::GetVarGrid(std::string name)
{
  //input variables
  if      (name=="precipitation"){return GRID_HRU;}
  else if (name=="temp_ave"     ){return GRID_HRU;}

  //output variables
  if      (name=="streamflow"   ){return GRID_SUBBASIN;}
  else if (name=="soil[0]"      ){return GRID_HRU;}
  else if (name=="snow"         ){return GRID_HRU;}

  return 0;
}
//////////////////////////////////////////////////////////////////
/// \brief returns units of input or output variable as string
/// \param name [in] - name of  variable
/// \return units of input or output variable as string
//
std::string CRavenBMI::GetVarUnits(std::string name)
{
  //input variables
  if      (name=="precipitation"){return "mm/d";}
  else if (name=="temp_ave"     ){return "C";}

  //output variables
  if      (name=="streamflow"   ){return "m3/s";}
  else if (name=="soil[0]"      ){return "mm";  }
  else if (name=="snow"         ){return "mm";  }

  return "";
}
//////////////////////////////////////////////////////////////////
/// \brief returns type of input or output variable as string
/// \param name [in] - name of  variable
/// \return type of input (now, assumes all outputs are double precision numbers)
//
std::string CRavenBMI::GetVarType(std::string name)
{
  return "double";
}
//////////////////////////////////////////////////////////////////
/// \brief returns size of input or output variable, in bytes
/// \param name [in] - name of  variable
/// \return size of individual input or output variable item, in bytes
//
int CRavenBMI::GetVarItemsize(std::string name)
{
  return sizeof(double);
}
//////////////////////////////////////////////////////////////////
/// \brief returns total size of input or output variable array, in bytes
/// \param name [in] - name of  variable
/// \return total size of input or output variable array, in bytes
//
int CRavenBMI::GetVarNbytes(std::string name)
{
  return GetVarItemsize(name)*GetGridSize(GetVarGrid(name));
}
//////////////////////////////////////////////////////////////////
/// \brief returns variable location on unstructured grid
/// \param name [in] - name of  variable
/// \return "node", always
//
std::string CRavenBMI::GetVarLocation(std::string name)
{
  return "node";
}
//////////////////////////////////////////////////////////////////
/// \brief accessors for model basic elements
//
double CRavenBMI::GetCurrentTime()   {  return tt.model_time;    }
double CRavenBMI::GetStartTime()     {  return 0.0;              }
double CRavenBMI::GetEndTime()       {  return Options.duration; }
double CRavenBMI::GetTimeStep()      {  return Options.timestep; }
std::string CRavenBMI::GetTimeUnits(){  return "d";              }

//------------------------------------------------------------------
// Variable getters
//------------------------------------------------------------------

//////////////////////////////////////////////////////////////////
/// \brief returns array of variable values for variable with supplied name
/// \param name [in] - name of  variable
/// \param dest [out] - pointer to array of variable values
//
void CRavenBMI::GetValue(std::string name, void* dest)
{
  double *out=NULL;
  int k,p;
  bool is_HRU_SV=false;
  int iSV=DOESNT_EXIST;

  //output variables
  if      (name=="streamflow")
  {
    out=new double [pModel->GetNumSubBasins()];
    for (p = 0; p < pModel->GetNumSubBasins(); p++)
    {
      if (Options.ave_hydrograph){
        out[p]=pModel->GetSubBasin(p)->GetIntegratedOutflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);
      }
      else{
        out[p]=pModel->GetSubBasin(p)->GetOutflowRate();
      }
    }
    memcpy(dest,out,pModel->GetNumSubBasins()*sizeof(double));
    delete [] out;
  }
  else if (name == "soil[0]") {is_HRU_SV=true; iSV=pModel->GetStateVarIndex(SOIL,0);}
  else if (name == "snow")    {is_HRU_SV=true; iSV=pModel->GetStateVarIndex(SNOW);  }
  else {
    throw std::logic_error("RavenBMI.GetValue: variable '" + name + "' not covered by this function.");
    return;
  }

  if ((is_HRU_SV) && (iSV!=DOESNT_EXIST))
  {
    out=new double [pModel->GetNumHRUs()];
    for (k = 0; k < pModel->GetNumHRUs(); k++) {
      out[k]=pModel->GetHydroUnit(k)->GetStateVarArray()[iSV];
    }
    memcpy(dest,out,pModel->GetNumSubBasins()*sizeof(double));
    delete [] out;
  }

}
//////////////////////////////////////////////////////////////////
/// \brief returns array of variable values for variable with supplied name at subset o flocations
/// \param name [in] - name of  variable
/// \param dest [out] - pointer to array of variable values
/// \param inds [in] - array of array indices/grid nodes from which to grab variables
/// \param count[in] - size of inds and dest arrays
//
void CRavenBMI::GetValueAtIndices(std::string name, void* dest, int* inds, int count)
{
  //output variables
  double *out=new double [count];
  bool is_HRU_SV=false;
  int iSV=DOESNT_EXIST;
  int i,k,p;

  if      (name=="streamflow")
  {
    for (i = 0; i <count; i++) {
      p=inds[i];
      out[p]=pModel->GetSubBasin(p)->GetIntegratedOutflow(Options.timestep);
    }
  }
  else if (name == "soil[0]") {is_HRU_SV=true; iSV=pModel->GetStateVarIndex(SOIL,0);}
  else if (name == "snow")    {is_HRU_SV=true; iSV=pModel->GetStateVarIndex(SNOW);  }

  if ((is_HRU_SV) && (iSV!=DOESNT_EXIST))
  {
    for (i = 0; i <count; i++) {
      k=inds[i];
      out[k]=pModel->GetHydroUnit(k)->GetStateVarArray()[iSV];
    }
  }
  memcpy(dest,out,count*sizeof(double));
  delete [] out;
}
//////////////////////////////////////////////////////////////////
/// \brief unavailable BMI functionality
/// \param name [in] - name of  variable
//
void *CRavenBMI::GetValuePtr(std::string name)
{
  double *out = new double[pModel->GetNumSubBasins()];  // allocate memory for output
  try {
    this->GetValue(name, out);                            // get the value
  } catch (std::logic_error &e) {
    delete[] out;
    throw std::logic_error(std::string("CRavenBMI.GetValuePtr: ") + e.what());
    return NULL;
  }
  return out;
}

//------------------------------------------------------------------
// Variable setters
//------------------------------------------------------------------

//////////////////////////////////////////////////////////////////
/// \brief sets values for full variable array with supplied name
/// \param name [in] - name of  variable
/// \param src  [in] - pointer to array of variable values
//
void CRavenBMI::SetValue(std::string name, void* src)
{
  //input variables
  double *input=(double*)src;
  forcing_type Ftype;
  bool is_forcing=false;

  if      (name=="precipitation"){is_forcing=true; Ftype=F_PRECIP;  }
  else if (name=="temp_ave")     {is_forcing=true; Ftype=F_TEMP_AVE;}

  if (is_forcing){
    for (int k=0;k<pModel->GetNumHRUs();k++){
      pModel->GetHydroUnit(k)->SetHRUForcing(Ftype, input[k]);
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief sets values for subset of variable array with supplied name
/// \param name [in] - name of  variable
/// \param src  [in] - pointer to array of variable values
/// \param inds [in] - array of array indices/grid nodes from which to grab variables
/// \param count[in] - size of inds and src arrays
//
void CRavenBMI::SetValueAtIndices(std::string name, int* inds, int count, void* src)
{
  //input variables
  double *input=(double*)src;
  forcing_type Ftype;
  bool is_forcing=false;

  if      (name=="precipitation"){is_forcing=true; Ftype=F_PRECIP;}
  else if (name=="temp_ave")     {is_forcing=true; Ftype=F_TEMP_AVE;}

  if (is_forcing){
    for (int i=0;i<count;i++){
      int k=inds[i];
      pModel->GetHydroUnit(k)->SetHRUForcing(Ftype, input[k]);
    }
  }
}

//------------------------------------------------------------------
// Grid information functions
//------------------------------------------------------------------

//////////////////////////////////////////////////////////////////
/// \brief returns rank of grid
/// \param grid [in]
/// \return rank of grid (always 1 for Raven unstructured 'grid' - HRUs or Basins)
//
int CRavenBMI::GetGridRank(const int grid)
{
  return 1; //vector rank (0=scalar, 1=1-D vector, 2=grid)
}
//////////////////////////////////////////////////////////////////
/// \brief returns size of unstructured grid
/// \param grid [in]
/// \return size of grid (# HRUs or #SubBasins)
//
int CRavenBMI::GetGridSize(const int grid)
{
  if      (grid==GRID_HRU     ){return pModel->GetNumHRUs();      }
  else if (grid==GRID_SUBBASIN){return pModel->GetNumSubBasins(); }

  return 0;
}
//////////////////////////////////////////////////////////////////
/// \brief returns type of grid
/// \param grid [in]
/// \returns 'unstructured', always
//
std::string CRavenBMI::GetGridType(const int grid)
{
  return "unstructured";
}
//////////////////////////////////////////////////////////////////
/// \brief returns x-coordinates of centroids of subbasins or HRUs
/// \param grid [in]
/// \param x [out] - x locations of centroids of subbasins or HRUs
//
void CRavenBMI::GetGridX(const int grid, double *x)
{
  if      (grid == GRID_HRU) {
    for (int k = 0; k < pModel->GetNumHRUs(); k++) {
      x[k]=pModel->GetHydroUnit(k)->GetCentroid().longitude;
    }
  }
  else if (grid == GRID_SUBBASIN)
  {
    throw std::logic_error("Not Implemented");
    //x[p]=pModel->GetSubBasin(p)->GetCentroid().longitude;
  }
}
void CRavenBMI::GetGridY(const int grid, double *y)
{
  if (grid == GRID_HRU) {
    for (int k = 0; k < pModel->GetNumHRUs(); k++) {
      y[k]=pModel->GetHydroUnit(k)->GetCentroid().latitude;
    }
  }
  else if (grid == GRID_SUBBASIN)
  {
    throw std::logic_error("Not Implemented");
    //y[p]=pModel->GetSubBasin(p)->GetCentroid().latitude;
  }
}
void CRavenBMI::GetGridZ(const int grid, double *z){
  if (grid == GRID_HRU) {
    for (int k = 0; k < pModel->GetNumHRUs(); k++) {
      z[k]=pModel->GetHydroUnit(k)->GetElevation();
    }
  }
  else if (grid == GRID_SUBBASIN)
  {
    throw std::logic_error("Not Implemented");
    //for (int p=0; p<pModel->GetNumSubBasins(); p++){
    //z[p]=pModel->GetSubBasin(p)->GetAvgElevation();
    //}
  }
}

int  CRavenBMI::GetGridNodeCount(const int grid)
{
  return GetGridSize(grid);
}
int  CRavenBMI::GetGridEdgeCount   (const int grid){return 0;}
int  CRavenBMI::GetGridFaceCount   (const int grid){return 0;}

//the following are not valid for unstructured grids (e.g., HRUs and basins)
void CRavenBMI::GetGridShape       (const int grid, int *shape)         {throw std::logic_error("Not Implemented");}
void CRavenBMI::GetGridSpacing     (const int grid, double *spacing)    {throw std::logic_error("Not Implemented");}
void CRavenBMI::GetGridOrigin      (const int grid, double *origin)     {throw std::logic_error("Not Implemented");}
void CRavenBMI::GetGridEdgeNodes   (const int grid, int *edge_nodes)    {throw std::logic_error("Not Implemented");}
void CRavenBMI::GetGridFaceEdges   (const int grid, int *face_edges)    {throw std::logic_error("Not Implemented");}
void CRavenBMI::GetGridFaceNodes   (const int grid, int *face_nodes)    {throw std::logic_error("Not Implemented");}
void CRavenBMI::GetGridNodesPerFace(const int grid, int *nodes_per_face){throw std::logic_error("Not Implemented");}
