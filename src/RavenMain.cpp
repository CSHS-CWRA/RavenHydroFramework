/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2023 the Raven Development Team
  ----------------------------------------------------------------*/
#include <time.h>
#include "RavenInclude.h"
#include "Model.h"
#include "UnitTesting.h"

// Function Declarations------------------------------------------
//Defined in ParseInput.cpp
bool ParseInputFiles       (CModel*&pModel,        optStruct &Options);
bool ParseInitialConditions(CModel*& pModel, const optStruct &Options);

//Defined in Solvers.cpp
void MassEnergyBalance     (CModel *pModel,const optStruct   &Options,const time_struct &tt);        
void ParseLiveFile         (CModel*&pModel,const optStruct   &Options, const time_struct &tt);

//Local functions defined below main() in RavenMain.cpp
void ProcessExecutableArguments(int argc, char* argv[], optStruct   &Options);
void CheckForErrorWarnings     (bool quiet);
bool CheckForStopfile          (const int step, const time_struct &tt);
void CallExternalScript        (const optStruct &Options, const time_struct &tt);

// Main Driver Variables------------------------------------------
static optStruct   Options;
static CModel      *pModel;

// Global variables - declared as extern in RavenInclude.h--------
string g_output_directory ="";
bool   g_suppress_warnings=false;
bool   g_suppress_zeros   =false;
double g_debug_vars[10];
bool   g_disable_freezing =false;
double g_min_storage      =0.0;
int    g_current_e        =0;
static string RavenBuildDate(__DATE__);

//////////////////////////////////////////////////////////////////
//
/// \brief Primary Raven driver routine
//
/// \param argc [in] number of arguments to executable
/// \param argv[] [in] executable arguments; Raven.exe [base_filename] [-p rvp_file] [...] [-t rvt_file] [-o output_dir]
/// for using WD\output subdirectory, can use "-o .\output\"
/// \return Success of main method
//
int main(int argc, char* argv[])
{
  double      t;
  clock_t     t0, t1;          //computational time markers
  time_struct tt;
  int         nEnsembleMembers;
  
  Options.version="3.6";
#ifdef _NETCDF_ 
  Options.version+=" w/ netCDF";
#endif
  
  ProcessExecutableArguments(argc,argv,Options);
  PrepareOutputdirectory(Options);

  for (int i=0;i<10;i++){g_debug_vars[i]=0;}

  RavenUnitTesting(Options);

  if (!Options.silent){
    int year = s_to_i(RavenBuildDate.substr(RavenBuildDate.length()-4,4).c_str());
    cout <<"============================================================"<<endl;
    cout <<"                        RAVEN                               "<<endl;
    cout <<" a robust semi-distributed hydrological modelling framework "<<endl;
    cout <<"    Copyright 2008-"<<year<<", the Raven Development Team "  <<endl;
    cout <<"                    Version "<<Options.version               <<endl;
    cout <<"                BuildDate "<<RavenBuildDate                  <<endl;
    cout <<"============================================================"<<endl;
  }

  ofstream WARNINGS;
  WARNINGS.open((Options.main_output_dir+"Raven_errors.txt").c_str());
  if (WARNINGS.fail()){
    ExitGracefully("Main::Unable to open Raven_errors.txt. Bad output directory specified?",RAVEN_OPEN_ERR);
  }
  WARNINGS.close();
  
  t0=clock();

  CStateVariable::Initialize();

  //Read input files, create model, set model options
  if (!ParseInputFiles(pModel, Options)){
    ExitGracefully("Main::Unable to read input file(s)",BAD_DATA);}
   
  CheckForErrorWarnings(true);

  if (!Options.silent){
    cout <<"======================================================"<<endl;
    cout <<"Initializing Model..."<<endl;
  }
  pModel->Initialize                  (Options);
  ParseInitialConditions              (pModel, Options);
  pModel->CalculateInitialWaterStorage(Options);
  pModel->SummarizeToScreen           (Options);
  pModel->GetEnsemble()->Initialize   (pModel,Options);

  CheckForErrorWarnings(false);

  nEnsembleMembers=pModel->GetEnsemble()->GetNumMembers();
    
  for(int e=0;e<nEnsembleMembers; e++) //only run once in standard mode
  {
    
    pModel->GetEnsemble()->UpdateModel(pModel,Options,e);
    PrepareOutputdirectory(Options); //adds new output folders, if needed
    pModel->WriteOutputFileHeaders(Options);
      
    if(!Options.silent) {
      cout <<endl<<"======================================================"<<endl;
      if(nEnsembleMembers>1) { cout<<"Ensemble Member "<<e+1<<" "; g_suppress_warnings=true;}
      cout <<"Simulation Start..."<<endl;
    }
    
    double t_start=0.0;
    t_start=pModel->GetEnsemble()->GetStartTime(e);
      
    //Write initial conditions-------------------------------------
    JulianConvert(t_start,Options.julian_start_day,Options.julian_start_year,Options.calendar,tt);
    pModel->RecalculateHRUDerivedParams(Options,tt);
    pModel->UpdateHRUForcingFunctions  (Options,tt);
    pModel->UpdateDiagnostics          (Options,tt);
    pModel->WriteMinorOutput           (Options,tt);

    //Solve water/energy balance over time--------------------------------
    t1=clock();
    int step=0;
      
    for(t=t_start; t<Options.duration-TIME_CORRECTION; t+=Options.timestep)  // in [d]
    {
      pModel->UpdateTransientParams      (Options,tt);
      pModel->RecalculateHRUDerivedParams(Options,tt);
      pModel->UpdateHRUForcingFunctions  (Options,tt);
      pModel->PrepareAssimilation        (Options,tt);
      pModel->WriteSimpleOutput          (Options,tt);
      pModel->GetEnsemble()->StartTimeStepOps(pModel,Options,tt,e);
      CallExternalScript                 (Options,tt);
      ParseLiveFile                      (pModel,Options,tt);

      MassEnergyBalance(pModel,Options,tt); //where the magic happens!

      pModel->IncrementCumulInput        (Options,tt);
      pModel->IncrementCumOutflow        (Options,tt);

      JulianConvert(t+Options.timestep,Options.julian_start_day,Options.julian_start_year,Options.calendar,tt);//increments time structure

      pModel->WriteMinorOutput           (Options,tt);
      pModel->WriteProgressOutput        (Options,clock()-t1,step,(int)ceil(Options.duration/Options.timestep));
      pModel->UpdateDiagnostics          (Options,tt); //required to read stuff!! 
      pModel->GetEnsemble()->CloseTimeStepOps(pModel,Options,tt,e);

      if ((Options.use_stopfile) && (CheckForStopfile(step,tt))) { break; }
      step++;
    }

    //Finished Solving----------------------------------------------------
    pModel->UpdateDiagnostics (Options,tt);
    pModel->RunDiagnostics    (Options);
    pModel->WriteMajorOutput  (Options,tt,"solution",true);
    pModel->CloseOutputStreams();

    if(!Options.silent) 
    {
      cout <<"======================================================"<<endl;
      cout <<"...Raven Simulation Complete: "<<Options.run_name<<endl;
      cout <<"    Parsing & initialization: "<< float(t1     -t0)/CLOCKS_PER_SEC << " seconds elapsed . "<<endl;
      cout <<"                  Simulation: "<< float(clock()-t1)/CLOCKS_PER_SEC << " seconds elapsed . "<<endl;
      if(Options.output_dir!="") {
        cout <<"  Output written to "        << Options.output_dir                                       <<endl;
      }
      cout <<"======================================================"<<endl;
    }
    if (Options.benchmarking) {
      cout <<"                              "<< pModel->GetNumHRUs()*(Options.duration/Options.timestep)/(float(clock()-t1)/CLOCKS_PER_SEC)<<" HRU-time steps/second"<<endl;
    }

    pModel->GetEnsemble()->FinishEnsembleRun(pModel,Options,tt,e);
  }/* end ensemble loop*/


  ExitGracefully("Successful Simulation",SIMULATION_DONE);
  return 0;
}

//////////////////////////////////////////////////////////////////
/// \param argc [in] number of arguments to executable
/// \param argv[] [in] executable arguments; Raven.exe [filebase] [-p rvp_file] [-h hru_file] [-t rvt_file] [-c rvc_file] [-o output_dir]
/// \details initializes input files and output directory
/// \details filebase has no extension, all others require .rv* extension
/// \param Options [in] Global model options
//
void ProcessExecutableArguments(int argc, char* argv[], optStruct   &Options)
{
  int i=1;
  string word,argument;
  bool version_announce=false;
  int mode=0;
  argument="";
  //initialization:
  Options.run_name    ="";
  Options.run_mode    =' ';
  Options.rvi_filename="";
  Options.rvh_filename="";
  Options.rvp_filename="";
  Options.rvt_filename="";
  Options.rvc_filename="";
  Options.rvg_filename="";
  Options.rve_filename="";
  Options.rvl_filename="";
  Options.output_dir  ="";
  Options.main_output_dir="";
  Options.silent=false;
  Options.noisy =false;
  Options.pause =true;
  Options.forecast_shift=0.0;
  Options.warm_ensemble_run="";

  //Parse argument list
  while (i<=argc)
  {
    if (i!=argc){
      word=to_string(argv[i]);
    }
    if ((word=="-p") || (word=="-h") || (word=="-t") || (word=="-e") || (word=="-c") || (word=="-o") || 
        (word=="-s") || (word=="-r") || (word=="-n") || (word=="-l") || (word=="-m") || (word=="-v") || 
        (word=="-we")|| (word=="-tt")|| (i==argc))
    {
      if      (mode==0){
        Options.rvi_filename=argument+".rvi";
        Options.rvp_filename=argument+".rvp";
        Options.rvh_filename=argument+".rvh";
        Options.rvt_filename=argument+".rvt";
        Options.rvc_filename=argument+".rvc";
        Options.rvg_filename=argument+".rvg";
        Options.rve_filename=argument+".rve";
        Options.rvl_filename=argument+".rvl";
        argument="";
        mode=10;
      }
      else if (mode==1 ){Options.rvp_filename=argument; argument="";}
      else if (mode==2 ){Options.rvh_filename=argument; argument="";}
      else if (mode==3 ){Options.rvt_filename=argument; argument="";}
      else if (mode==4 ){Options.rvc_filename=argument; argument="";}
      else if (mode==5 ){Options.output_dir  =argument; argument="";}
      else if (mode==6 ){Options.run_name    =argument; argument="";}
      else if (mode==7 ){Options.rve_filename=argument; argument="";}
      else if (mode==8 ){Options.rvg_filename=argument; argument="";}
      else if (mode==9 ){Options.rvl_filename=argument; argument="";}
      else if (mode==11){Options.run_mode =argument[0]; argument="";}
      else if (mode==12){Options.forecast_shift=s_to_d(argument.c_str()); argument=""; }
      else if (mode==13){Options.warm_ensemble_run=argument; argument=""; }

      if      (word=="-p"){mode=1; }
      else if (word=="-h"){mode=2; }
      else if (word=="-t"){mode=3; }
      else if (word=="-c"){mode=4; }
      else if (word=="-o"){mode=5; }
      else if (word=="-s"){Options.silent=true; mode=10;}
      else if (word=="-n"){Options.noisy=true;  mode=10;}
      else if (word=="-r"){mode=6; }
      else if (word=="-e"){mode=7; }
      else if (word=="-g"){mode=8; }	  
      else if (word=="-l"){mode=9; }
      else if (word=="-m"){mode=11;}
      else if (word=="-tt"){mode=12; }
      else if (word=="-we"){mode=13; }
      else if (word=="-v"){Options.pause=false; version_announce=true; mode=10;} //For PAVICS
    }
    else{
      if (argument==""){argument+=word;}
      else             {argument+=" "+word;}
    }
    i++;
  }
  if (argc==1){//no arguments
    Options.rvi_filename="nomodel.rvi";
    Options.rvp_filename="nomodel.rvp";
    Options.rvh_filename="nomodel.rvh";
    Options.rvt_filename="nomodel.rvt";
    Options.rvc_filename="nomodel.rvc";
    Options.rvg_filename="nomodel.rvg";
    Options.rve_filename="nomodel.rve";
    Options.rvl_filename="nomodel.rvl";
  }

  // make sure that output dir has trailing '/' if not empty
  if ((Options.output_dir.compare("") != 0) && (Options.output_dir.back()!='/')){ Options.output_dir=Options.output_dir+"/"; }

  //convert to days 
  Options.forecast_shift/=HR_PER_DAY;

  char cCurrentPath[FILENAME_MAX];
  if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath))){
    ExitGracefully("RavenMain: unable to retrieve current directory.", RUNTIME_ERR);
  }
  Options.working_dir = to_string(cCurrentPath);
  Options.main_output_dir=Options.output_dir;

  if(version_announce) {
    cout<<Options.version<<endl;
    ExitGracefully("Version check",SIMULATION_DONE);
  }
}
/////////////////////////////////////////////////////////////////
/// \brief Exits gracefully from program, explaining reason for exit and destructing simulation all pertinent parameters
/// \remark Called from within code
///
/// \param statement [in] String to print to user upon exit
/// \param code [in] Code to determine why the system is exiting
//
void ExitGracefully(const char *statement,exitcode code)
{

  string typeline;
  switch (code){
  case(SIMULATION_DONE):  {typeline="============================================================";break;}
  case(RUNTIME_ERR):      {typeline="Error Type: Runtime Error";       break;}
  case(BAD_DATA):         {typeline="Error Type: Bad input data";      break;}
  case(BAD_DATA_WARN):    {typeline="Error Type: Bad input data";      break;}
  case(OUT_OF_MEMORY):    {typeline="Error Type: Out of memory";       break;}
  case(FILE_OPEN_ERR):    {typeline="Error Type: File opening error";  break;}
  case(STUB):             {typeline="Error Type: Stub function called";break;}
  default:                {typeline="Error Type: Unknown";             break;}
  }

  if (code != RAVEN_OPEN_ERR){//avoids recursion problems
    ofstream WARNINGS;
    WARNINGS.open((Options.main_output_dir+"Raven_errors.txt").c_str(),ios::app);
    if (WARNINGS.fail()){
      WARNINGS.close();
      string message="Unable to open errors file ("+Options.main_output_dir+"Raven_errors.txt)";
      ExitGracefully(message.c_str(),RAVEN_OPEN_ERR);
    }
    if (code!=SIMULATION_DONE){WARNINGS<<"ERROR : "<<statement<<endl;}
    else                      {WARNINGS<<"SIMULATION COMPLETE :)"<<endl;}
    WARNINGS.close();
  }
  if (code==BAD_DATA_WARN){return;}//just write these errors to a file if not in strict mode

  cout <<endl<<endl;
  cout <<"============== Exiting Gracefully =========================="<<endl;
  cout <<"Exiting Gracefully: "<<statement                             <<endl;
  cout << typeline                                                     <<endl;
  cout <<"============================================================"<<endl;

  delete pModel; pModel=NULL;//deletes EVERYTHING!
  CStateVariable::Destroy();

  if(Options.pause)
  {
    cout << "Press the ENTER key to continue"<<endl;
    cin.get();
  }
  exit(0);
}
/////////////////////////////////////////////////////////////////
/// \brief Checks if errors have been written to Raven_errors.txt, if so, exits gracefully
/// \note called prior to simulation initialization, after parsing everything
///
//
void CheckForErrorWarnings(bool quiet)
{
  int      Len;
  char    *s[MAXINPUTITEMS];
  bool     errors_found(false);
  bool     warnings_found(false);

  ifstream WARNINGS;
  WARNINGS.open((Options.main_output_dir+"Raven_errors.txt").c_str());
  if (WARNINGS.fail()){WARNINGS.close();return;}

  CParser *p=new CParser(WARNINGS,Options.main_output_dir+"Raven_errors.txt",0);

  while (!(p->Tokenize(s,Len)))
  {
    if(Len>0){
      if(!strcmp(s[0],"ERROR"  )){ errors_found  =true; }
      if(!strcmp(s[0],"WARNING")){ warnings_found=true; }
    }
  }
  WARNINGS.close();
  if ((warnings_found) && (!quiet)){
    cout<<"*******************************************************"<<endl<<endl;
    cout<<"WARNING: Warnings have been issued while parsing data. "<<endl;
    cout<<"         See Raven_errors.txt for details              "<<endl<<endl;
    cout<<"*******************************************************"<<endl<<endl;
  }

  if (errors_found){
    ExitGracefully("Errors found in input data. See Raven_errors.txt for details",BAD_DATA);
  }
}
/////////////////////////////////////////////////////////////////
/// \brief Checks if stopfile exists in current working directory
/// \note called during simulation to determine whether progress should be stopped
///
//
bool CheckForStopfile(const int step, const time_struct &tt)
{
  if(step%100!=0){ return false; } //only check every 100th timestep 
  ifstream STOP;
  STOP.open("stop");
  if (STOP.fail()){STOP.close(); return false;}
  else //Stopfile found
  {
    STOP.close();
    pModel->WriteMajorOutput  (Options,tt,"solution",true);
    pModel->CloseOutputStreams();
    ExitGracefully("CheckForStopfile: simulation interrupted by user using stopfile",SIMULATION_DONE);
    return true;
  }
}
/////////////////////////////////////////////////////////////////
/// \brief Calls external script to be run 
/// idea/code from Kai Tsuruta, PCIC
//
void CallExternalScript(const optStruct &Options,const time_struct &tt) 
{
  if(Options.external_script!="") {
    string script=Options.external_script;
    SubstringReplace(script,"<model_time>",to_string(tt.model_time));
    SubstringReplace(script,"<date>"      ,tt.date_string);
    SubstringReplace(script,"<version>"   ,Options.version);
    SubstringReplace(script,"<output_dir>",Options.output_dir);
    system(script.c_str()); //Calls script
  }
}
