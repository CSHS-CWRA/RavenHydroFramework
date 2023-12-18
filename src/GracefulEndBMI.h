/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2023 the Raven Development Team
----------------------------------------------------------------*/
#include "RavenInclude.h"

/////////////////////////////////////////////////////////////////
/// \brief Finalizes program gracefully, explaining reason for finalizing, and destructing simulation and all pertinent parameters
/// \remark Called from within ExitGracefully() or from the BMI Finalize() function
///
/// \param statement [in] String to print to user upon exit
/// \param code [in] Code to determine why the system is exiting
//
inline void FinalizeGracefully(const char *statement, exitcode code)
{
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

  // error message is printed to the standard error stream
  if (code != RAVEN_OPEN_ERR) { //avoids recursion problems
    if (code!=SIMULATION_DONE) { cerr << "ERROR : " << statement << endl;}
    else                       { cout << "SIMULATION COMPLETE :)" << endl;}
  }
  if (code==BAD_DATA_WARN){
    return;
  }

  // success messages are printed to the standard output stream
  cout << endl << endl;
  cout << "============== Exiting Gracefully ==========================" << endl;
  cout << "Exiting Gracefully: " << statement                            << endl;
  cout << typeline                                                       << endl;
  cout << "============================================================" << endl;

}


/////////////////////////////////////////////////////////////////
/// \brief Exits gracefully from program, explaining reason for exit and destructing simulation all pertinent parameters
/// \remark Called from within code
///
/// \param statement [in] String to print to user upon exit
/// \param code [in] Code to determine why the system is exiting
//
inline void ExitGracefully(const char *statement, exitcode code)
{
  FinalizeGracefully(statement, code);
}
