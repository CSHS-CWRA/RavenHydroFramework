/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2024 the Raven Development Team
  ----------------------------------------------------------------*/
#ifndef STATEVARIABLE_H
#define STATEVARIABLE_H

#include "RavenInclude.h"

class CTransportModel;

///////////////////////////////////////////////////////////////////
/// \brief Methods class for related state variable parsing and querying routines
/// \details Implemented in StateVariables.cpp
//
class CStateVariable
{
private:/*------------------------------------------------------*/

  int              _nAliases;          ///< total number of aliases
  string          *_aAliases;          ///< Array of alias values [size: nAliases]
  string          *_aAliasReferences;  ///< Array of strings referenced by aliases [size: nAliases]

  string           CheckAliasList      (const string s);
  CTransportModel *_pTransportModel;   ///< Pointer to transport model

public:/*-------------------------------------------------------*/

  CStateVariable();  ///< Constructor
  ~CStateVariable(); ///< Destructor

  static string        SVStringBreak       (const string s, int &num);

  void                 Initialize          ();
  void                 Destroy             ();

  void                 AddAlias            (const string s1, const string s2);

  void                 SetTransportModel   (CTransportModel *pTransportModel);

  //static functions
  static string        GetStateVarLongName (sv_type               typ,
                                            const int             layer_index,
                                            const CTransportModel *pTransportModel);
  static string        GetStateVarUnits    (sv_type      typ);                         // can be kept static
  sv_type              StringToSVType      (const string  s, int &layer_index, bool strict);
  string               SVTypeToString      (const sv_type typ, const int layerindex);

  static bool          IsWaterStorage      (sv_type      typ,bool conv_coverup=true);
  static bool          IsEnergyStorage     (sv_type      typ);
};
#endif
