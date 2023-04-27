/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2023 the Raven Development Team
  ----------------------------------------------------------------
  Class CGlobalParams
  ----------------------------------------------------------------*/
#ifndef GLOBAL_PARAMS_H
#define GLOBAL_PARAMS_H

#include "RavenInclude.h"
#include "Properties.h"

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for global model parameters
//
class CGlobalParams
{
protected:/*----------------------------------------------------*/

  static global_struct             G;   ///< global parameters

public:/*-------------------------------------------------------*/

  CGlobalParams();
  ~CGlobalParams();

  //Accessors
  static const global_struct    *GetParams();
  static double GetParameter(const string param_name);
  static void SetGlobalProperty          (const string  &param_name, const double &value);

  //routines
  static void AutoCalculateGlobalParams(const global_struct &Gtmp, const global_struct &Gdefault);


  static void InitializeGlobalParameters (global_struct &G, bool is_template);
  static void SetGlobalProperty          (global_struct &G, const string  param_name, const double value);
  static double GetGlobalProperty        (const global_struct &G, string  param_name, const bool strict=true);

  static double *GetAddress(const string param_name);

  static void SummarizeToScreen();
};

#endif
