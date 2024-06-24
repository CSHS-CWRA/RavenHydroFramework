/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2023 the Raven Development Team
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Model.h"

static CModel *pModel = NULL;  // in the standalone version, the pModel is a global variable

/////////////////////////////////////////////////////////////////
/// \brief Finalizes program gracefully, explaining reason for finalizing, and destructing simulation and all pertinent parameters
/// \remark Called from within ExitGracefully() or from the BMI Finalize() function
///
/// \param statement [in] String to print to user upon exit
/// \param code [in] Code to determine why the system is exiting
//
inline void FinalizeGracefully(const char *statement, exitcode code)
{
  // errors can occur earlier than the initialization of global pModel variable
  if (pModel == NULL) {
    cout << "============== Early Gracefully Exit ==============" << endl;
    cout << "Error Statement: " << statement << endl;
    if (code!=SIMULATION_DONE) {
    cout << "Error Code: " << code << endl;
    }
    cout << "---------------------------------------------------" << endl;
    return;
  }

  const optStruct* Options = pModel->GetOptStruct();

  string typeline;
  switch (code){
    case(SIMULATION_DONE): {typeline="============================================================";break;}
    case(RUNTIME_ERR):     {typeline="Error Type: Runtime Error";       break;}
    case(BAD_DATA):        {typeline="Error Type: Bad input data";      break;}
    case(BAD_DATA_WARN):   {typeline="Error Type: Bad input data";      break;}
    case(OUT_OF_MEMORY):   {typeline="Error Type: Out of memory";       break;}
    case(FILE_OPEN_ERR):   {typeline="Error Type: File opening error";  break;}
    case(STUB):            {typeline="Error Type: Stub function called";break;}
    default:               {typeline="Error Type: Unknown";             break;}
  }

  if (code != RAVEN_OPEN_ERR) { //avoids recursion problems
    ofstream WARNINGS;
    WARNINGS.open((Options->main_output_dir+"Raven_errors.txt").c_str(),ios::app);
    if (WARNINGS.fail()) {
      WARNINGS.close();
      string message="Unable to open errors file ("+Options->main_output_dir+"Raven_errors.txt)";
      ExitGracefully(message.c_str(),RAVEN_OPEN_ERR);
    }
    if (code!=SIMULATION_DONE) {WARNINGS<<"ERROR    : "<< statement << endl;
                                cerr    <<"ERROR    : "<< statement << endl;}
    else                       {WARNINGS<<"SIMULATION COMPLETE :)"<<endl;}

    WARNINGS.close();
  }
  if (code==BAD_DATA_WARN){return;}//just write these errors to a file if not in strict mode

  cout <<endl<<endl;
  cout <<"============== Exiting Gracefully =========================="<<endl;
  cout <<"Exiting Gracefully: "<<statement                             <<endl;
  cout << typeline                                                     <<endl;
  cout <<"============================================================"<<endl;

  delete pModel; pModel=NULL; //deletes EVERYTHING!

  if(Options->pause) {
    cout << "Press the ENTER key to continue"<<endl;
    cin.get();
  }
}


/////////////////////////////////////////////////////////////////
/// \brief Exits gracefully from program, explaining reason for exit and destructing simulation all pertinent parameters
/// \remark Called from within code
///
/// \param statement [in] String to print to user upon exit
/// \param code [in] Code to determine why the system is exiting
//
 void ExitGracefully(const char *statement, exitcode code)
{
  FinalizeGracefully(statement, code);
  if (code!=BAD_DATA_WARN){exit(0);}
}
