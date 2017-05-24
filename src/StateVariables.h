/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------*/
#ifndef STATEVARIABLE_H
#define STATEVARIABLE_H

#include "RavenInclude.h"

///////////////////////////////////////////////////////////////////
/// \brief Methods class for related state variable parsing and querying routines
/// \details Implemented in StateVariables.cpp
//
class CStateVariable
{
private:/*------------------------------------------------------*/

  static int            _nAliases;         ///< total number of aliases
  static string        *_aAliases;         ///< Array of alias values [size: nAliases]
  static string        *_aAliasReferences; ///< Array of strings referenced by aliases [size: nAliases]

  static string        CheckAliasList      (const string s);

public:/*-------------------------------------------------------*/

  static string        SVStringBreak       (const string s, int &num);

  static void          Initialize          ();
  static void          Destroy             ();

  static void          AddAlias            (const string s1, const string s2);

  //static functions
  static string        GetStateVarLongName (sv_type      typ, const int layer_index);
  static string        GetStateVarUnits    (sv_type      typ);
  static sv_type       StringToSVType      (const string s, int &layer_index, bool strict);
  static string        SVTypeToString      (const sv_type typ, const int layerindex);

  static bool          IsWaterStorage      (sv_type      typ);
  static bool          IsEnergyStorage     (sv_type      typ);
};
#endif
