/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ----------------------------------------------------------------*/
#ifndef MODELABC_H
#define MODELABC_H

#include "RavenInclude.h"


///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for surface water model with limited access to variables
/// \details Allows lesser classes (HydroProcesses, SubBasins, and Hydrounits)
///   to obtain information from model without changing the model
///   PURELY VIRTUAL CLASS
//
class CModelABC //version of model visible to lesser units
{
protected:/*------------------------------------------------------*/

public:/*-------------------------------------------------------*/

  //Accessor functions
  virtual int         GetNumStateVars    () const=0;
  virtual sv_type     GetStateVarType    (const int i) const=0;
  virtual int         GetStateVarIndex   (sv_type type) const=0;
  virtual int         GetStateVarIndex   (sv_type type, int layer) const=0;
  virtual int         GetStateVarLayer   (const int i) const=0;

  virtual double      GetConcentration   (const int k, const int i) const=0;
  virtual double      GetFlux            (const int k, const int js, const optStruct &Options) const=0;
  virtual double      GetLatFlow         (const int js, const optStruct &Options) const=0;
  virtual double      GetCumulativeFlux  (const int k, const int i, const bool to) const=0;
  virtual double      GetCumulFluxBetween(const int k,const int iFrom,const int iTo) const=0;

  virtual int         GetNumSoilLayers   () const=0;

  virtual double      GetAvgStateVar     (const int i) const=0;
  virtual double      GetAvgConcentration(const int i) const=0;
  virtual double      GetAvgCumulFlux    (const int i, const bool to) const=0;
  virtual double      GetAvgCumulFluxBet (const int iFrom, const int iTo) const=0;

  virtual bool        StateVarExists     (sv_type type) const=0;

  virtual bool        IsInHRUGroup       (const int k, const string HRUGroupName) const=0;

  virtual int         GetLakeStorageIndex() const=0;

  virtual const optStruct   *GetOptStruct() const = 0;
};

#endif
