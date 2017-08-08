/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------*/
#include <time.h>
#include "RavenInclude.h"
#include "Model.h"
#include "UnitTesting.h"

//Defined in ParseInput.cpp
bool ParseInputFiles (CModel      *&pModel,
                      optStruct    &Options);
//Defined in Solvers.cpp
void MassEnergyBalance(      CModel      *pModel,
                             const optStruct   &Options,
                             const time_struct &tt);        //time
//Defined below
void ProcessExecutableArguments(int argc, char* argv[], optStruct   &Options);
void CheckForErrorWarnings();
bool CheckForStopfile(const time_struct &tt);

// Main Driver Variables------------------------------------------
static optStruct   Options;
static CModel      *pModel;

// Global variables - declared as extern in RavenInclude.h--------
string g_output_directory="";
bool   g_suppress_warnings=false;
double g_debug_vars[5];

static string RavenBuildDate(__DATE__);

//////////////////////////////////////////////////////////////////
//
/// \brief Primary Raven driver routine
//
/// \param argc [in] number of arguments to executable
/// \param argv[] [in] executable arguments; Raven.exe [filebase] [-p rvp_file] [-h hru_file] [-t rvt_file] [-o output_dir]
/// for using WD\output subdirectory, can use "-o .\output\"
/// \return Success of main method
//
int main(int argc, char* argv[])
{
  double      t;
  string      filebase;
  clock_t     t0, t1;          //computational time markers
  time_struct tt;

  ProcessExecutableArguments(argc, argv, Options);
  PrepareOutputdirectory(Options);

  Options.pause=true;
  Options.version="2.7.1";

  for (int i=0;i<5;i++){g_debug_vars[i]=0;}

  RavenUnitTesting(Options);

  if (!Options.silent){
    int year = s_to_i(RavenBuildDate.substr(RavenBuildDate.length()-4,4).c_str());
    cout <<"=========================================================="<<endl;
    cout <<"                        RAVEN                             "<<endl;
    cout <<" a numerically robust semi-distributed hydrological model "<<endl;
    cout <<"    Copyright 2008-"<<year<<", the Raven Development Team "<<endl;
    cout <<"                    Version "<<Options.version             <<endl;
    cout <<"                BuildDate "<<RavenBuildDate                <<endl;
    cout <<"=========================================================="<<endl;
  }

  ofstream WARNINGS((Options.output_dir+"Raven_errors.txt").c_str());
  if (WARNINGS.fail()){
    ExitGracefully("Main::Unable to open Raven_errors.txt. Bad output directory specified?",RUNTIME_ERR);
  }
  WARNINGS.close();
  t0=clock();

  CStateVariable::Initialize();

  //Read input files, create model, set model options
  if (ParseInputFiles(pModel, Options))
  {
    CheckForErrorWarnings();

    if (!Options.silent){
      cout <<"======================================================"<<endl;
      cout <<"Initializing Model..."<<endl;}
    pModel->Initialize       (Options);
    pModel->SummarizeToScreen(Options);

    if (!Options.silent){
      cout <<"======================================================"<<endl;
      cout <<"Simulation Start..."<<endl;}

    //Write initial conditions-------------------------------------
    JulianConvert(0.0,Options.julian_start_day,Options.julian_start_year,tt);
    pModel->RecalculateHRUDerivedParams(Options,tt);
    pModel->UpdateHRUForcingFunctions  (Options,tt);
    pModel->UpdateDiagnostics          (Options,tt);
    pModel->WriteMinorOutput           (Options,tt);

    //Solve water/energy balance over time--------------------------------
    t1=clock();
    int step=0;

    for (t=0; t<Options.duration-TIME_CORRECTION; t+=Options.timestep)
    {
      pModel->UpdateTransientParams      (Options,tt);
      pModel->RecalculateHRUDerivedParams(Options,tt);
      pModel->UpdateHRUForcingFunctions  (Options,tt);
      pModel->UpdateDiagnostics          (Options,tt);

      MassEnergyBalance(pModel,Options,tt); //where the magic happens!

      pModel->IncrementCumulInput       (Options,tt);
      pModel->IncrementCumOutflow       (Options);

      JulianConvert(t+Options.timestep,Options.julian_start_day,Options.julian_start_year,tt);//increments time structure
      pModel->WriteMinorOutput          (Options,tt);

      if ((step%100==0) && (CheckForStopfile(tt))){break;}
      step++;
    }

    //Finished Solving----------------------------------------------------
    pModel->UpdateDiagnostics (Options,tt);
    pModel->RunDiagnostics    (Options);
    pModel->WriteMajorOutput  ("solution",Options,tt,true);
    pModel->CloseOutputStreams();

    if (!Options.silent){
      cout <<"======================================================"<<endl;
      cout <<"...Simulation Complete: "   <<Options.run_name<<endl;
      cout <<"  Parsing & initialization: "<< float(t1     -t0)/CLOCKS_PER_SEC << " seconds elapsed . "<<endl;
      cout <<"                Simulation: "<< float(clock()-t1)/CLOCKS_PER_SEC << " seconds elapsed . "<<endl;
      if (Options.output_dir!=""){
        cout <<"  Output written to "        <<Options.output_dir                                        <<endl;
      }
      cout <<"======================================================"<<endl;
    }
  }
  else
  {
    ExitGracefully("Main::Unable to read input file(s)",BAD_DATA);
  }

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
  int mode=0;
  argument="";
  //initialization:
  Options.run_name    ="";
  Options.rvi_filename="";
  Options.rvh_filename="";
  Options.rvp_filename="";
  Options.rvt_filename="";
  Options.rvc_filename="";
  Options.output_dir="";
  Options.silent=false;
  Options.noisy =false;

  //Parse argument list
  while (i<=argc)
  {
    if (i!=argc){
      word=to_string(argv[i]);
    }
    if ((word=="-p") || (word=="-h") || (word=="-t") || (word=="-c") || (word=="-o") || (word=="-s") || (word=="-r") || (i==argc))
    {
      if      (mode==0){
        Options.rvi_filename=argument+".rvi";
        Options.rvp_filename=argument+".rvp";
        Options.rvh_filename=argument+".rvh";
        Options.rvt_filename=argument+".rvt";
        Options.rvc_filename=argument+".rvc";
        argument="";
        mode=10;
      }
      else if (mode==1){Options.rvp_filename=argument; argument="";}
      else if (mode==2){Options.rvh_filename=argument; argument="";}
      else if (mode==3){Options.rvt_filename=argument; argument="";}
      else if (mode==4){Options.rvc_filename=argument; argument="";}
      else if (mode==5){Options.output_dir  =argument; argument="";}
      else if (mode==6){Options.run_name    =argument; argument="";}
      if      (word=="-p"){mode=1;}
      else if (word=="-h"){mode=2;}
      else if (word=="-t"){mode=3;}
      else if (word=="-c"){mode=4;}
      else if (word=="-o"){mode=5;}
      else if (word=="-s"){Options.silent=true;}//should be specified prior to other flags
      else if (word=="-r"){mode=6;}
    }
    else{
      if (argument==""){argument+=word;}
      else             {argument+=" "+word;}
    }
    i++;
  }
  if (argc==1){//no arguments
    Options.rvi_filename="model.rvi";
    Options.rvp_filename="model.rvp";
    Options.rvh_filename="model.rvh";
    Options.rvt_filename="model.rvt";
    Options.rvc_filename="model.rvc";
  }

  char cCurrentPath[FILENAME_MAX];
  if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath))){
    ExitGracefully("RavenMain: unable to retrieve current directory.", RUNTIME_ERR);
  }
  Options.working_dir = to_string(cCurrentPath);

  // identify executable directory
  //char basePath[255] = "";
  //_fullpath(basePath, argv[0], sizeof(basePath)); //_realpath in linux
  //cout << to_string(basePath) << endl;

  /*cout<<"PROCESSED:"<<endl;
    cout<<Options.rvi_filename<<endl;
    cout<<Options.rvp_filename<<endl;
    cout<<Options.rvh_filename<<endl;
    cout<<Options.rvt_filename<<endl;
    cout<<Options.rvc_filename<<endl;
    cout<<Options.output_dir<<endl;*/
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
  case(SIMULATION_DONE):  {typeline="===============================================";break;}
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
    WARNINGS.open((Options.output_dir+"Raven_errors.txt").c_str(),ios::app);
    if (WARNINGS.fail()){
      string message="Unable to open errors file ("+Options.output_dir+"Raven_errors.txt)";
      ExitGracefully(message.c_str(),RAVEN_OPEN_ERR);
    }
    if (code!=SIMULATION_DONE){WARNINGS<<"ERROR : "<<statement<<endl;}
    else                      {WARNINGS<<"SIMULATION COMPLETE :)"<<endl;}
    WARNINGS.close();
  }
  if (code==BAD_DATA_WARN){return;}//just write these errors to a file if not in strict mode

  cout <<endl<<endl;
  cout <<"===============Exiting Gracefully=============="<<endl;
  cout <<"Exiting Gracefully: "<<statement                <<endl;
  cout << typeline                                        <<endl;
  cout <<"==============================================="<<endl;

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
void CheckForErrorWarnings()
{

  ifstream WARNINGS;
  WARNINGS.open((Options.output_dir+"Raven_errors.txt").c_str());
  if (WARNINGS.fail()){WARNINGS.close();return;}

  CParser *p=new CParser(WARNINGS,Options.output_dir+"Raven_errors.txt",0);
  int      Len;
  char    *s[MAXINPUTITEMS];
  bool     errors_found(false);
  bool     warnings_found(false);

  while (!(p->Tokenize(s,Len)))
  {
    if  (!strcmp(s[0],"ERROR"  )){errors_found  =true;}
    if  (!strcmp(s[0],"WARNING")){warnings_found=true;}
  }
  WARNINGS.close();
  if (warnings_found){// && (!Options.silent){
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
bool CheckForStopfile(const time_struct &tt)
{
  ifstream STOP;
  STOP.open("stop");
  if (!STOP.is_open()){return false;}
  else //Stopfile found
  {
    STOP.close();
    pModel->WriteMajorOutput  ("solution",Options,tt,true);
    pModel->CloseOutputStreams();
    ExitGracefully("CheckForStopfile: simulation interrupted by user using stopfile",SIMULATION_DONE);
    return true;
  }
}

