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

// constants used in the config file (cfg_file) parsing
const std::string CFG_FILE_CLIARGS_TAG = "cli_args";
const std::string CFG_FILE_INPVARS_TAG = "input_vars";
const std::string CFG_FILE_OUTVARS_TAG = "output_vars";
const std::string CFG_FILE_IOVAR_FORCE_TAG = "forcing";
const std::string CFG_FILE_IOVAR_STATE_TAG = "state";
const std::string CFG_FILE_IOVAR_DERIV_TAG = "derived";

//////////////////////////////////////////////////////////////////
/// \brief RavenBMI class constructor and destructor
//
CRavenBMI::CRavenBMI()
{
  pModel=NULL;

  // vectors storing input/output variable names and ids are filled in the Initialize function
  input_var_names = std::vector<std::string>();
  input_var_ids   = std::vector<int>();
  output_var_names= std::vector<std::string>();
  output_var_ids  = std::vector<int>();
  output_var_layer_index = std::vector<int>();
  output_var_type = std::vector<std::string>();
}

CRavenBMI::~CRavenBMI() {}


//////////////////////////////////////////////////////////////////
/// \brief Splits a string by whitespace into a vector of char*
///
/// \param line [in] string to be split
/// \return vector of char* containing the split string
std::vector<char *> CRavenBMI::_SplitLineByWhitespace(std::string line)
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
/// \brief Splits a string by colon into a vector of char*
/// 
/// \param line [in] string to be split
/// \return vector of char* containing the split string
std::vector<std::string> CRavenBMI::_SplitLineByColon(std::string line)
{
  int config_str_ini, config_str_mid;
  std::vector<std::string> args;

  config_str_ini = 0;
  config_str_mid = (int)line.find(":", config_str_ini);

  // simple case: no colon found - return a list of one element (the whole string)
  if (config_str_mid == -1) {
    args.push_back(line);
    return args;
  }

  // first string is the key, add it to the list
  args.push_back(line.substr(config_str_ini, config_str_mid - config_str_ini));

  // even if there is nothing after the colon, add an empty string to the list to communicate that
  //  the line HAD a final colon
  args.push_back(line.substr(config_str_mid + 1));

  return args;
}


//////////////////////////////////////////////////////////////////
/// \brief Checks if a derived variable is valid
///
/// \param var_name [in] name of the variable
/// \return true if the variable is valid, false otherwise
bool CRavenBMI::_IsValidDerivedVariable(std::string var_name)
{
  if (var_name == "OUTFLOW") {
    return true;
  }
  return false;
}


//////////////////////////////////////////////////////////////////
/// \brief Reads the BMI configuration file and checks if model options are valid
///
/// \param config_file [in] name of configuration file
/// \return void - sets values in Options struct
void CRavenBMI::_ReadConfigFile(std::string config_file)
{
  ifstream CONFIG;
  int config_str_ini, config_str_end;    // used to parse the lines of the config file
  std::string config_key, config_value;  // used to parse the lines of the config file
  bool cli_args = false;                 // used to check if command line arguments were used
  std::vector<char *> args;              // used to parse the command line arguments
  char** argv;
  int config_value_id;
  std::vector<std::string> line_split_by_colon;
  bool listing_inp_vars = false;         // flag to check if the config file is listing input variables
  bool listing_out_vars = false;         // flag to check if the config file is listing output variables

  CONFIG.open(config_file);
  if (CONFIG.fail()) {
    throw std::logic_error("Cannot find configuration file " + config_file);
    return;
  }

  // read and parse line by line
  for( std::string line; getline( CONFIG, line ); )
  {
    // skip blank lines and comments
    if (line.empty() || line[0] == '#') {
      continue;
    }

    // all lines of the config file must have one colon
    line_split_by_colon = _SplitLineByColon(line);
    if (line_split_by_colon.size() != 2) {
      throw std::logic_error("WARNING: Invalid line in Raven config file: '" + line + "'");
    }

    // just to make the code more readable
    config_key = line_split_by_colon[0];
    config_value = line_split_by_colon[1];

    // remove leading and trailing whitespaces
    config_key.erase(0, config_key.find_first_not_of(" \t"));
    config_key.erase(config_key.find_last_not_of(" \t") + 1);
    config_value.erase(0, config_value.find_first_not_of(" \t"));
    config_value.erase(config_value.find_last_not_of(" \t") + 1);

    if (config_key == CFG_FILE_CLIARGS_TAG) {

      cli_args = true;
      // "Raven.exe" is a dummy argument to mimic the command line arguments
      args = _SplitLineByWhitespace("Raven.exe " + config_value);
      argv = args.data();
      ProcessExecutableArguments((int)(args.size())-1, argv, Options);
      listing_inp_vars = false;
      listing_out_vars = false;
      continue;
      
    } else if (config_key == CFG_FILE_INPVARS_TAG) {

      listing_inp_vars = true;
      listing_out_vars = false;
      continue;

    } else if (config_key == CFG_FILE_OUTVARS_TAG) {

      listing_out_vars = true;
      listing_inp_vars = false;
      continue;

    } else if (listing_inp_vars) {

      // all lines after the "input_vars:" line must be in the form "- [VAR_NAME]: [var_type]"
      if (config_key.substr(0, 2) != "- ") {
        throw std::logic_error("WARNING: Invalid input var line in Raven config file: " + line);
      }

      // remove the '- ' from the beginning of the line
      config_key = config_key.substr(2);

      if (config_value == CFG_FILE_IOVAR_FORCE_TAG) {
        // case: regular - forcing as input
        config_value_id = GetForcingTypeFromString(config_key);
        if (config_value_id == F_UNRECOGNIZED) {
          throw std::logic_error("WARNING: Invalid forcing type in Raven config file: '" + line + "'");
          continue;
        }
        input_var_names.push_back(config_key);
        input_var_ids.push_back(config_value_id);
        continue;
      } else if (config_value == CFG_FILE_IOVAR_STATE_TAG) {
        // case: data assimilation - state as input
        // TODO: implement
        throw std::logic_error("WARNING: State variables as input are not supported in the BMI interface yet.");
        continue;
      } else if (config_value == CFG_FILE_IOVAR_DERIV_TAG) {
        // case: data assimilation - derivative as input
        // TODO: implement
        throw std::logic_error("WARNING: Derivative variables as input are not supported in the BMI interface yet.");
        continue;
      }
      throw std::logic_error("WARNING: Invalid input var line in Raven config file: " + line);
      continue;

    } else if (listing_out_vars) {

      // all lines after the "input_vars:" line must be in the form "- [VAR_NAME]: [var_type]"
      if (config_key.substr(0, 2) != "- ") {
        throw std::logic_error("WARNING: Invalid input var line in Raven config file: " + line );
        continue;
      }

      // remove the '- ' from the beginning of the line
      config_key = config_key.substr(2);

      if (config_value == CFG_FILE_IOVAR_STATE_TAG) {
        // case: regular - state as output (evaluated after the model initialized)
        output_var_names.push_back(config_key);
        output_var_type.push_back(CFG_FILE_IOVAR_STATE_TAG);
        continue;

      } else if (config_value == CFG_FILE_IOVAR_DERIV_TAG) {
        // case: regular - special cases as output
        if (_IsValidDerivedVariable(config_key)) {
          output_var_names.push_back(config_key);
          output_var_type.push_back(CFG_FILE_IOVAR_DERIV_TAG);
        } else {
          throw std::logic_error("WARNING: Invalid derived variable in Raven config file: '" + line + "'");
        }
        continue;

      } else if (config_value == CFG_FILE_IOVAR_FORCE_TAG) {
        // case: forcing as output? does it make sense?
        throw std::logic_error("WARNING: Forcing variables as output are not supported in Raven's BMI interface.");
        continue;

      }

      // if we got here, the line is invalid
      throw std::logic_error("WARNING: Invalid output var line in Raven config file: " + line);
      continue;
    }

    // if we got here, the line is invalid
    throw std::logic_error("WARNING: Invalid line in Raven config file: " + line);

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

  _ReadConfigFile(config_file);

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

  // all the output variables must be checked
  int config_out_var_layer_index;
  int config_value_id;
  for (int i = 0; i < output_var_names.size(); i++) {
    if (output_var_type[i] == CFG_FILE_IOVAR_DERIV_TAG) {
      output_var_ids.push_back(-1);
      output_var_layer_index.push_back(-1);
      continue;
    }
    if (output_var_type[i] == CFG_FILE_IOVAR_STATE_TAG) {
      config_value_id = pModel->GetStateVarInfo()->StringToSVType(output_var_names[i],
                                                                  config_out_var_layer_index,
                                                                  true);
      output_var_ids.push_back(config_value_id);
      output_var_layer_index.push_back(config_out_var_layer_index);
      continue;
    }
    throw std::logic_error("WARNING: Output variable '" + output_var_names[i] + "' has an invalid type '" + output_var_type[i] + "'.");
  }

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
  return (input_var_names.size());
}

//////////////////////////////////////////////////////////////////
/// \brief returns array of accessible input dataset names
/// \return array of accessible input dataset names
//
std::vector<std::string> CRavenBMI::GetInputVarNames()
{
  return(input_var_names);
}

//////////////////////////////////////////////////////////////////
/// \brief returns number of accessible output datasets
/// \return number of accessible output datasets
//
int CRavenBMI::GetOutputItemCount()
{
  return (output_var_names.size());
}

//////////////////////////////////////////////////////////////////
/// \brief returns array of accessible output dataset names
/// \return array of accessible output dataset names
//
std::vector<std::string> CRavenBMI::GetOutputVarNames()
{
  return(output_var_names);
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
  // TODO - adjust it

  // special cases are handled first to make them feel... special
  if (name == "OUTFLOW") {
    return(GRID_SUBBASIN);
  }

  // the variable can be a forcing, so we need to search for it
  forcing_type Ftype = GetForcingTypeFromString(name);
  if (Ftype != F_UNRECOGNIZED) {
    return(GRID_HRU);  // is this a good assumption?
  }

  // if not a special case nor a forcing, it MUST be a state var (or the code will just break)
  int state_var_layer_idx;
  sv_type SType = pModel->GetStateVarInfo()->StringToSVType(name, state_var_layer_idx, true);
  return(GRID_HRU);    // is this a good assumption?
}

//////////////////////////////////////////////////////////////////
/// \brief returns units of input or output variable as string
/// \param name [in] - name of  variable
/// \return units of input or output variable as string
//
std::string CRavenBMI::GetVarUnits(std::string name)
{
  // special cases are handled first to make them feel... special
  if (name == "OUTFLOW") {
    return("m3/s");
  }

  // the variable can be a forcing, so we need to search for it
  forcing_type Ftype = GetForcingTypeFromString(name);
  if (Ftype != F_UNRECOGNIZED) {
    return(GetForcingTypeUnits(Ftype));
  }

  // if not a special case nor a forcing, it MUST be a state var (or the code will just break)
  int state_var_layer_idx;
  sv_type SType = pModel->GetStateVarInfo()->StringToSVType(name, state_var_layer_idx, true);
  return(CStateVariable::GetStateVarUnits(SType));
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
  sv_type Stype= (sv_type)DOESNT_EXIST;
  int Slayer=DOESNT_EXIST;

  // special case: output variables
  if (name == "OUTFLOW" )
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
    return;
  }
  else {
    // search in the output variable names
    for (int i = 0; i < output_var_names.size(); i++) {
      if (output_var_names[i] != name) {
        continue;
      }
      if (output_var_type[i] != CFG_FILE_IOVAR_STATE_TAG) {
        throw std::logic_error("CRavenBMI.GetValue: variable '" + name + "' is not a state variable.");
        return;
      }
      
      Stype = (sv_type)output_var_ids[i];
      Slayer = output_var_layer_index[i];
      is_HRU_SV=true;
      if (Slayer > -1) {
        iSV = pModel->GetStateVarIndex(Stype, Slayer);
      } else {
        iSV = pModel->GetStateVarIndex(Stype);
      }

      break;
    }
  }

  if ((is_HRU_SV) && (iSV!=DOESNT_EXIST))
  {
    out=new double [pModel->GetNumHRUs()];
    for (k = 0; k < pModel->GetNumHRUs(); k++) {
      out[k]=pModel->GetHydroUnit(k)->GetStateVarArray()[iSV];
    }
    memcpy(dest,out,pModel->GetNumSubBasins()*sizeof(double));
    delete [] out;
    return;
  }
  
  // if it got here... dude, I have a bad news for you...
  throw std::logic_error("RavenBMI.GetValue: variable '" + name + "' not covered by this function.");
  return;
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
  sv_type Stype=(sv_type)DOESNT_EXIST;
  int Slayer=DOESNT_EXIST;

  // special case: output variables
  if (name == "OUTFLOW" ) {
    for (i = 0; i <count; i++) {
      p=inds[i];
      out[p]=pModel->GetSubBasin(p)->GetIntegratedOutflow(Options.timestep);
    }
  }
  else {
    // search in the output variable names
    for (int i = 0; i < output_var_names.size(); i++) {
      if (output_var_names[i] != name) {
        continue;
      }
      if (output_var_type[i] != CFG_FILE_IOVAR_STATE_TAG) {
        throw std::logic_error("CRavenBMI.GetValue: variable '" + name + "' is not a state variable.");
        return;
      }
      
      Stype = (sv_type)output_var_ids[i];
      Slayer = output_var_layer_index[i];
      is_HRU_SV=true;
      if (Slayer > -1) {
        iSV = pModel->GetStateVarIndex(Stype, Slayer);
      } else {
        iSV = pModel->GetStateVarIndex(Stype);
      }

      break;
    }
  }

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

  Ftype = GetForcingTypeFromString(name);
  if (Ftype != F_UNRECOGNIZED) {
    is_forcing = true;
  } else {
    throw std::logic_error("CRavenBMI.SetValue: only accepts setting forcing values (for now).");
  }

  // all HRUs receive the same forcing
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

  Ftype = GetForcingTypeFromString(name);
  if (Ftype != F_UNRECOGNIZED) {
    is_forcing = true;
  } else {
    throw std::logic_error("WARNING: CRavenBMI.SetValueAtIndices: only accepts forcings now. Tried to set value for variable '"+name+"'");
  }

  // all HRUs receive the same forcing
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
