/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ----------------------------------------------------------------
  class definitions:
  CmvOWEvaporation
  CmvLakeEvaporation
  ----------------------------------------------------------------*/
#ifndef OPENWATER_EVAP_H
#define OPENWATER_EVAP_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"

///////////////////////////////////////////////////////////////////
/// \brief Types of open water evaporation models
//
enum owevap_type
{
  OPEN_WATER_EVAP,        ///< evaporates at OW_PET rate
  OPEN_WATER_RIPARIAN,    ///< evaporates at OW_PET rate, but only for fraction of landscape covered in streams
  OPEN_WATER_UWFS         ///< Upscaled wetland fill and spill approach - tracks depression storage and deficit distribution
};

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for open water evaporation
/// \details Calculates loss of water from open water to atmosphere
//
class CmvOWEvaporation: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  owevap_type   _type; ///< Model of open water evaporation used

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvOWEvaporation(owevap_type ow_type, const int i_from);      //general constructor
  ~CmvOWEvaporation();

  //inherited functions
  void Initialize();
  void GetRatesOfChange(const double              *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double      *rates) const;
  void ApplyConstraints(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double      *rates) const;

  static void GetParticipatingStateVarList(owevap_type ow_type,
                                           sv_type *aSV, int *aLev, int &nSV);
  void GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const;
};

///////////////////////////////////////////////////////////////////
/// \brief Types of lake evaporation
//
enum lakeevap_type
{
  LAKE_EVAP_BASIC,                              ///< constant evaporation rate - fraction of PET
};

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for open water (lake) evaporation
/// \details Calculates loss of water from open water (lake) to atmosphere
//
class CmvLakeEvaporation: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  lakeevap_type                  type; ///< Model of lake evaporation used

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvLakeEvaporation(lakeevap_type lk_type, const int fromIndex);       //general constructor
  ~CmvLakeEvaporation();

  //inherited functions
  void Initialize();
  void GetRatesOfChange(const double              *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double      *rates) const;
  void ApplyConstraints(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double      *rates) const;

  static void GetParticipatingStateVarList(lakeevap_type lk_type,
                                           sv_type *aSV, int *aLev, int &nSV);
  void GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const;
};

#endif
