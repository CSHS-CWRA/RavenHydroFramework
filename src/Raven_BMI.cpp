/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2026 the Raven Development Team
  ----------------------------------------------------------------*/
#include <time.h>
#include "RavenInclude.h"
#include "RavenMain.h"
#include "Raven_BMI.h"

void ParseManagementFile       (CModel *&pModel, const optStruct &Options);

//////////////////////////////////////////////////////////////////
/// \brief RavenBMI class constructor and destructor
//
CRavenBMI::CRavenBMI()
{
  pModel=NULL;
  _duration=0;

  // vectors storing input/output variable names and ids are filled in the Initialize function
  _input_vars  = std::vector<rvn_var_data>();
  _output_vars = std::vector<rvn_var_data>();
}

CRavenBMI::~CRavenBMI(){}


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
/// \brief Checks if a subbasin variable is valid
///
/// \param var_name [in] name of the variable
/// \return true if the variable is valid, false otherwise
bool CRavenBMI::_IsValidSubBasinStateVariable(std::string var_name)
{
  if ((var_name == "OUTFLOW") || (var_name == "STREAMFLOW") ||
      (var_name == "RESERVOIR_STAGE")) {
    return true;
  }
  return false;
}


//////////////////////////////////////////////////////////////////
/// \brief Reads the BMI configuration file and checks if model options are valid
///
/// \param config_file [in] name of configuration file
/// \return void - sets values in Options struct
///
/// Example Config file:
/// cli_args: -t input.rvt -o ./output/ #equivalent to executable args
/// input_vars:
/// - SOIL[0]: hru_state
/// - SNOW: hru_state
/// - RAINFALL: forcing
/// output_vars:
///
void CRavenBMI::_ReadConfigFile(std::string config_file)
{
  string config_key, config_value;  // used to parse the lines of the config file
  bool cli_args = false;                 // used to check if command line arguments were used
  std::vector<char *> args;              // used to parse the command line arguments
  char** argv;

  std::vector<string> line_split_by_colon;
  bool listing_inp_vars = false;         // flag to check if the config file is listing input variables
  bool listing_out_vars = false;         // flag to check if the config file is listing output variables

  ifstream CONFIG;
  CONFIG.open(config_file);
  if (CONFIG.fail()) {
    cout<<"ERROR: Cannot find configuration file " + config_file<<endl;
    throw std::logic_error("Cannot find configuration file " + config_file);
    return;
  }

  // read and parse line by line
  for( std::string line; getline( CONFIG, line ); )
  {
    std::cout<<" BMIconfig line :" <<line<<std::endl;
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
    config_key   = line_split_by_colon[0];
    config_value = line_split_by_colon[1];

    // remove leading and trailing whitespaces
    config_key.erase(0, config_key.find_first_not_of(" \t"));
    config_key.erase(config_key.find_last_not_of(" \t") + 1);
    config_value.erase(0, config_value.find_first_not_of(" \t"));
    config_value.erase(config_value.find_last_not_of(" \t") + 1);

    if (config_key == "cli_args")
    {
      cli_args = true;
      args = _SplitLineByWhitespace("Raven.exe " + config_value);
      argv = args.data();
      ProcessExecutableArguments((int)(args.size())-1, argv, Options);

      listing_inp_vars = false;
      listing_out_vars = false;
      continue;
    }
    else if(config_key == "duration")
    {
      _duration=s_to_d(config_value.c_str());

      listing_inp_vars = false;
      listing_out_vars = false;
      continue;
    }
    else if (config_key == "input_vars")
    {
      listing_inp_vars = true;
      listing_out_vars = false;
      continue;
    }
    else if (config_key == "output_vars")
    {
      listing_out_vars = true;
      listing_inp_vars = false;
      continue;
    }
    else if ((listing_inp_vars) || (listing_out_vars))
    {
      // all lines after the "input_vars:" line must be in the form "- [VAR_NAME]: [var_type]"
      if (config_key.substr(0, 2) != "- ") {
        throw std::logic_error("WARNING: Invalid input var line in Raven config file: " + line);
      }

      // remove the '- ' from the beginning of the line
      config_key = config_key.substr(2);

      //std::cout<<" BMIconfig line :" <<line<<std::endl;

      diagnostic var_type;
      if      (config_value=="forcing"  ){var_type=VAR_FORCING_FUNCTION;}
      else if (config_value=="hru_state"){var_type=VAR_STATE_VAR;}
      else if (config_value=="subbasin_state")
      {
        if      (config_key=="STREAMFLOW"     ){var_type=VAR_STREAMFLOW; }
        else if (config_key=="RESERVOIR_STAGE"){var_type=VAR_RESERVOIR_STAGE;}
      //else if (config_key=="RUNOFF"         ){var_type=VAR_RUNOFF;} //to return Qlat
        else {
          throw std::logic_error("WARNING: Invalid subbasin state variable type in Raven config file input_vars: '" + line + "'");
          continue;
        }
      }
      else if(config_value=="flux")
      {
        // case: model interfacing
        // TODO: implement
        throw std::logic_error("WARNING: flux variables as input are not supported in the BMI interface yet.");
        continue;
      }
      else{
        throw std::logic_error("WARNING: Invalid input var line in Raven config file: " + line);
        continue;
      }

      rvn_var_data *tmp=new rvn_var_data(var_type,config_key);

      if      (listing_inp_vars){ _input_vars.push_back(*tmp);}
      else if (listing_out_vars){_output_vars.push_back(*tmp);}
    }
  }
  CONFIG.close();

  if (!cli_args){
    throw std::logic_error("RavenBMI:ReadConfigFile: missing cli_args command in configuration file");
  }
  if(_duration==0.0) {
    throw std::logic_error("RavenBMI:ReadConfigFile: missing duration command in configuration file");
  }
  return;
}
//////////////////////////////////////////////////////////////////
/// \brief populates extra data in output and input vars
/// verifies whether variables are valid raven variables
///
//
void CRavenBMI::_CheckConfigVars(std::vector<rvn_var_data>  &vars)
{
  // all the input & output variables must be checked
  for (int i = 0; i < vars.size(); i++)
  {
    if      (vars[i].type == VAR_STREAMFLOW){}
    else if (vars[i].type == VAR_RESERVOIR_STAGE){}
    else if (vars[i].type == VAR_FORCING_FUNCTION)
    {
      forcing_type ftype=GetForcingTypeFromString(vars[i].name);
      vars[i].f_type=ftype;

      if(ftype==UNRECOGNIZED_SVTYPE) {
        throw std::logic_error("WARNING: config variable '" + vars[i].name + "' has an invalid state variable type.");
        return;
      }
    }
    else if (vars[i].type == VAR_STATE_VAR)
    {
      int layer_ind;
      sv_type typ = pModel->GetStateVarInfo()->StringToSVType(vars[i].name,layer_ind,true);
      vars[i].sv_layer_ind=layer_ind;
      vars[i].sv_type     =typ;

      if (typ==UNRECOGNIZED_SVTYPE){
        throw std::logic_error("WARNING: config variable '" + vars[i].name + "' has an invalid state variable type.");
        return;
      }
    }
    else{
      throw std::logic_error("WARNING: config variable '" + vars[i].name + "' has an invalid or unsupported type '" + to_string(vars[i].type) + "'.");
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Initialization called prior to model simulation
/// this routine echoes the contents of RavenMain.cpp and must be updated accordingly when the main routine is revised
///
/// \param config_file [in] name of configuration file (unused)
/// NOTE:   ENSEMBLE MODE NOT SUPPORTED WITH BMI
//
void CRavenBMI::Initialize(std::string config_file)
{
  //cout<<"BMI - READING CONFIG FILE"<<endl;
  _ReadConfigFile(config_file);

  //cout<<"BMI - PREPARING OUTPUT DIR"<<endl;
  PrepareOutputdirectory(Options);

  for (int i=0;i<10;i++){g_debug_vars[i]=0;}

  ofstream WARNINGS;
  WARNINGS.open((Options.main_output_dir+"Raven_errors.txt").c_str());
  if (WARNINGS.fail()){
    ExitGracefully("Main::Unable to open Raven_errors.txt. Bad output directory specified?",RAVEN_OPEN_ERR);
  }
  WARNINGS.close();

  Options.in_bmi_mode  = true;  // flag needed to ignore some arguments of the .rvi file

  //cout<<"BMI - PARSING RAVEN INPUT FILES"<<endl;

  //Read input files, create model, set model options
  if (!ParseInputFiles(pModel, Options)){
    ExitGracefully("Main::Unable to read input file(s)",BAD_DATA);}

  CheckForErrorWarnings(true, pModel);

  //model duration overridden by config file
  //This must be done prior to initialize so time series sampling is done correctly
  if (_duration>Options.duration){
    throw std::logic_error("WARNING: Duration indicated in .rvi file must be greater than config file indicates");
  }
  Options.duration = _duration;

  pModel->Initialize                  (Options);

  ParseManagementFile                 (pModel,Options);
  pModel->InitializePostRVM           (Options);

  ParseInitialConditions              (pModel, Options);
  pModel->CalculateInitialWaterStorage(Options);
  pModel->SummarizeToScreen           (Options);

  CheckForErrorWarnings(false, pModel);

  //Checks whether config variables are valid in this Raven model
  _CheckConfigVars(_input_vars);
  _CheckConfigVars(_output_vars);

  PrepareOutputdirectory(Options); //adds new output folders, if needed
  pModel->WriteOutputFileHeaders(Options);

  double t_start=0.0;
  //t_start=pModel->GetEnsemble()->GetStartTime(0);

  //Write initial conditions-------------------------------------
  JulianConvert(t_start,Options.julian_start_day,Options.julian_start_year,Options.calendar,tt);
  pModel->RecalculateHRUDerivedParams(Options,tt);
  pModel->UpdateHRUForcingFunctions  (Options,tt);
  pModel->UpdateDiagnostics          (Options,tt);
  pModel->InitializePostRVC          (Options   );
  pModel->WriteMinorOutput           (Options,tt);
}

//////////////////////////////////////////////////////////////////
/// \brief run simulation for a single time step
///
//
void CRavenBMI::Update()
{
  pModel->UpdateTransientParams      (Options,tt);
  pModel->ApplyStateOverrrides       (Options,tt);
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
  pModel->UpdateDiagnostics          (Options,tt); //required to read stuff!!

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
  for(double t=tt.model_time; t<(time+0.1*Options.timestep); t+=Options.timestep)  // in [d]
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
  return (int)(_input_vars.size());
}

//////////////////////////////////////////////////////////////////
/// \brief returns array of accessible input dataset names
/// \return array of accessible input dataset names
//
std::vector<std::string> CRavenBMI::GetInputVarNames()
{
  std::vector<string> names= std::vector<string>();

  for (int i=0; i<_input_vars.size(); i++){
    names.push_back(_input_vars[i].name);
  }
  return names;
}

//////////////////////////////////////////////////////////////////
/// \brief returns number of accessible output datasets
/// \return number of accessible output datasets
//
int CRavenBMI::GetOutputItemCount()
{
  return (int)(_output_vars.size());
}

//////////////////////////////////////////////////////////////////
/// \brief returns array of accessible output dataset names
/// \return array of accessible output dataset names
//
std::vector<std::string> CRavenBMI::GetOutputVarNames()
{
  std::vector<string> names= std::vector<string>();

  for (int i=0; i<_output_vars.size(); i++){
    names.push_back(_output_vars[i].name);
  }
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

  for (int i = 0; i < _output_vars.size(); i++) {
    if (_output_vars[i].name == name)
    {
      if      (_output_vars[i].type==VAR_STREAMFLOW)      {return GRID_SUBBASIN;}
      else if (_output_vars[i].type==VAR_RESERVOIR_STAGE) {return GRID_SUBBASIN;}
      else if (_output_vars[i].type==VAR_FORCING_FUNCTION){return GRID_HRU;}
      else if (_output_vars[i].type==VAR_STATE_VAR)       {return GRID_HRU; }
    }
  }
  for (int i = 0; i < _input_vars.size(); i++) {
    if (_input_vars[i].name == name)
    {
      if      (_input_vars[i].type==VAR_STREAMFLOW)      {return GRID_SUBBASIN;}
      else if (_input_vars[i].type==VAR_RESERVOIR_STAGE) {return GRID_SUBBASIN;}
      else if (_input_vars[i].type==VAR_FORCING_FUNCTION){return GRID_HRU;}
      else if (_input_vars[i].type==VAR_STATE_VAR)       {return GRID_HRU; }
    }
  }

  throw std::logic_error("RavenBMI.GetVarUnits: variable '" + name + "' is invalid, not in config.txt, or not yet supported.");
  return 0;
}

//////////////////////////////////////////////////////////////////
/// \brief returns units of input or output variable as string
/// \param name [in] - name of  variable
/// \return units of input or output variable as string
//
std::string CRavenBMI::GetVarUnits(std::string name)
{
  //NEEDS TO BE BMI FORMAT UNITS - will work for mm but not for (e.g., MJ/m2/d which should be MJ m-2 d-1
  // todo: add bmi argument to GetForcingTypeUnits and GetStateVarUnits
  for (int i = 0; i < _output_vars.size(); i++) {
    if (_output_vars[i].name == name)
    {
      if      (_output_vars[i].type==VAR_STREAMFLOW)      {return "m3 s-1";}
      else if (_output_vars[i].type==VAR_RESERVOIR_STAGE) {return "m";}
      else if (_output_vars[i].type==VAR_FORCING_FUNCTION){return GetForcingTypeUnits(_output_vars[i].f_type);}
      else if (_output_vars[i].type==VAR_STATE_VAR)       {return CStateVariable::GetStateVarUnits(_output_vars[i].sv_type); }
    }
  }
  for (int i = 0; i < _input_vars.size(); i++) {
    if (_input_vars[i].name == name)
    {
      if      (_input_vars[i].type==VAR_STREAMFLOW)      {return "m3 s-1";}
      else if (_input_vars[i].type==VAR_RESERVOIR_STAGE) {return "m";}
      else if (_input_vars[i].type==VAR_FORCING_FUNCTION){return GetForcingTypeUnits(_input_vars[i].f_type);}
      else if (_input_vars[i].type==VAR_STATE_VAR)       {return CStateVariable::GetStateVarUnits(_input_vars[i].sv_type); }
    }
  }

  throw std::logic_error("RavenBMI.GetVarUnits: variable '" + name + "' is invalid, not in config.txt, or not yet supported.");
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
  int k,p,iSV;

  for (int i = 0; i < _output_vars.size(); i++) {
    if (_output_vars[i].name == name)
    {
      if(_output_vars[i].type== VAR_STREAMFLOW) {
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
      }
      else if(_output_vars[i].type== VAR_RESERVOIR_STAGE)
      {
        out=new double [pModel->GetNumSubBasins()];
        for (p = 0; p < pModel->GetNumSubBasins(); p++) {
          CSubBasin* pBasin = pModel->GetSubBasin(p);
          out[p] = 0.0;
          if (pBasin->GetReservoir() != NULL) {
            out[p] = pBasin->GetReservoir()->GetResStage();
          }
        }
      }
      else if (_output_vars[i].type==VAR_STATE_VAR){
        iSV = pModel->GetStateVarIndex(_output_vars[i].sv_type, _output_vars[i].sv_layer_ind);
        out=new double[pModel->GetNumHRUs()];
        for (k = 0; k < pModel->GetNumHRUs(); k++) {
          out[k]=pModel->GetHydroUnit(k)->GetStateVarArray()[iSV];
        }
      }
      else if (_output_vars[i].type==VAR_FORCING_FUNCTION)
      {
        out=new double[pModel->GetNumHRUs()];
        for (k = 0; k < pModel->GetNumHRUs(); k++) {
          out[k]=pModel->GetHydroUnit(k)->GetForcing(_output_vars[i].f_type);
        }
      }
    }
  }
  if (out!=NULL){
    memcpy(dest,out,pModel->GetNumHRUs()*sizeof(double));
    delete [] out;
  }
  else{
    throw std::logic_error("CRavenBMI.GetValue: unsupported variable name or not included in config file.");
  }

  return;
}

//////////////////////////////////////////////////////////////////
/// \brief returns array of variable values for variable with supplied name at subset of locations
/// \param name [in] - name of  variable
/// \param dest [out] - pointer to array of variable values
/// \param inds [in] - array of array indices/grid nodes from which to grab variables
/// \param count[in] - size of inds and dest arrays
//
void CRavenBMI::GetValueAtIndices(std::string name, void* dest, int* inds, int count)
{
  int iSV,k,p;

  // search in the output variable names
  for (int i = 0; i < _output_vars.size(); i++) {
    if (_output_vars[i].name == name) {
      double *out=new double[count];
      if      (_output_vars[i].type== VAR_STREAMFLOW){
        for (i = 0; i <count; i++) {
          p=inds[i];
          if (Options.ave_hydrograph){
            out[p]=pModel->GetSubBasin(p)->GetIntegratedOutflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);
          }
          else{
            out[p]=pModel->GetSubBasin(p)->GetOutflowRate();
          }
        }
      }
      else if (_output_vars[i].type== VAR_RESERVOIR_STAGE)
      {
        for (i = 0; i <count; i++) {
          p=inds[i];
          CSubBasin* pBasin = pModel->GetSubBasin(p);
          out[p] = 0.0;
          if (pBasin->GetReservoir() != NULL) {
            out[p] = pBasin->GetReservoir()->GetResStage();
          }
        }
      }
      else if (_output_vars[i].type== VAR_FORCING_FUNCTION)
      {
        for(i = 0; i <count; i++) {
          k=inds[i];
          out[k]=pModel->GetHydroUnit(k)->GetForcing(_output_vars[i].f_type);
        }
      }
      else if (_output_vars[i].type== VAR_STATE_VAR)
      {
        iSV = pModel->GetStateVarIndex(_output_vars[i].sv_type, _output_vars[i].sv_layer_ind);
        for(i = 0; i <count; i++) {
          k=inds[i];
          out[k]=pModel->GetHydroUnit(k)->GetStateVarArray()[iSV];
        }
      }
      memcpy(dest,out,count*sizeof(double));
      delete [] out;
      return;
    }
  }
  throw std::logic_error("CRavenBMI.GetValueAtIndices: unsupported variable name or not included in config file.");
  return;

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
  }
  catch(std::logic_error &e) {
    delete[] out;
    throw std::logic_error(std::string("CRavenBMI.GetValuePtr: ") + e.what());
    return NULL;
  }
  return out;
}

//------------------------------------------------------------------
// Manipulators
//------------------------------------------------------------------

//////////////////////////////////////////////////////////////////
/// \brief sets values for full variable array with supplied name
/// \param name [in] - name of  variable
/// \param src  [in] - pointer to array of variable values
//
void CRavenBMI::SetValue(std::string name, void* src)
{
  double *input=(double*)src;
  for (int i=0; i<_input_vars.size(); i++){
    if (_input_vars[i].name==name){
      if      (_input_vars[i].type==VAR_STREAMFLOW)
      {
        for(int p=0;p<pModel->GetNumSubBasins();p++) {
          pModel->GetSubBasin(p)->SetQout(input[p]);;
        }
      }
      else if (_input_vars[i].type==VAR_RESERVOIR_STAGE)
      {
        for (int p=0;p<pModel->GetNumSubBasins();p++){
          CReservoir *pRes=pModel->GetSubBasin(p)->GetReservoir();
          if (pRes!=NULL){
            pRes->SetReservoirStage(input[p],RAV_BLANK_DATA);
          }
        }
      }
      else if (_input_vars[i].type==VAR_FORCING_FUNCTION)
      {
        for (int k=0;k<pModel->GetNumHRUs();k++){
          pModel->GetHydroUnit(k)->SetHRUForcing(_input_vars[i].f_type, input[k]);
        }
      }
      else if (_input_vars[i].type==VAR_STATE_VAR)
      {
        int iSV=pModel->GetStateVarIndex(_input_vars[i].sv_type,_input_vars[i].sv_layer_ind);

        for (int k=0;k<pModel->GetNumHRUs();k++){
          pModel->GetHydroUnit(k)->SetStateVarValue(iSV, input[k]);
        }
      }
      return;
    }
  }
  throw std::logic_error("CRavenBMI.SetValue: this input variable is not accessible or was not included as such in the config.txt file.");
  src=NULL;
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
  int p,k;
  double *input=(double*)src;

  for (int i=0; i<_input_vars.size(); i++){
    if (_input_vars[i].name==name){
      if      (_input_vars[i].type==VAR_STREAMFLOW)
      {
        for(int j=0;j<count;j++) {
          p=inds[j];
          pModel->GetSubBasin(p)->SetQout(input[j]);;
        }
      }
      else if (_input_vars[i].type==VAR_RESERVOIR_STAGE)
      {
        for(int j=0;j<count;j++) {
          p=inds[j];
          CReservoir *pRes=pModel->GetSubBasin(p)->GetReservoir();
          if (pRes!=NULL){
            pRes->SetReservoirStage(input[j],RAV_BLANK_DATA);
          }
        }
      }
      else if (_input_vars[i].type==VAR_FORCING_FUNCTION)
      {
        for(int j=0;j<count;j++) {
          k=inds[j];
          pModel->GetHydroUnit(k)->SetHRUForcing(_input_vars[i].f_type, input[j]);
        }
      }
      else if (_input_vars[i].type==VAR_STATE_VAR)
      {
        int iSV=pModel->GetStateVarIndex(_input_vars[i].sv_type,_input_vars[i].sv_layer_ind);
        for(int j=0;j<count;j++) {
          k=inds[j];
          pModel->GetHydroUnit(k)->SetStateVarValue(iSV, input[j]);
        }
      }
      return;
    }
  }
  throw std::logic_error("CRavenBMI.SetValue: this input variable is not accessible or was not included as such in the config.txt file.");
  src=NULL;
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
