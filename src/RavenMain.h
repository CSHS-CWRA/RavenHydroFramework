/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2023 the Raven Development Team
  ----------------------------------------------------------------*/
#ifndef RAVEN_MAIN
#define RAVEN_MAIN

#include "RavenInclude.h"
#include "Model.h"

//Defined in ParseInput.cpp
bool ParseInputFiles       (CModel*&pModel,        optStruct &Options);
bool ParseInitialConditions(CModel*& pModel, const optStruct &Options);

//Defined in Solvers.cpp
void MassEnergyBalance     (CModel *pModel,const optStruct   &Options, const time_struct &tt);
void ParseLiveFile         (CModel*&pModel,const optStruct   &Options, const time_struct &tt);

//Local functions defined below main() in RavenMain.cpp
void ProcessExecutableArguments(int argc, char* argv[], optStruct   &Options);
void CheckForErrorWarnings     (bool quiet, CModel *pModel);
bool CheckForStopfile          (const int step, const time_struct &tt, CModel *pModel);
void CallExternalScript        (const optStruct &Options, const time_struct &tt);

#endif
