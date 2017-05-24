/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------
  class definitions:
  CmvGlacierMelt
  CmvGlacierRunoff
  ----------------------------------------------------------------*/

#ifndef GLACIER_PROCESSES_H
#define GLACIER_PROCESSES_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"

///////////////////////////////////////////////////////////////////
/// \brief Methods of modeling glacier melt
//
enum glacial_melt_type{
  GMELT_SIMPLE_MELT, ///< simple potential melt approach
  GMELT_HBV,         ///< HBV-EC Degree-day method for modellign glacier melt
  GMELT_UBC          ///< UBC method for modellign glacier melt
};

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for melt from glacier
//
class CmvGlacierMelt: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  glacial_melt_type type; ///< Specified algorithm for glacier melt simulation

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvGlacierMelt(const glacial_melt_type mtype);
  ~CmvGlacierMelt();

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
                        const time_struct &t,
                        double      *rates) const;

  void        GetParticipatingParamList   (string  *aP, class_type *aPC, int &nP) const;
  static void GetParticipatingStateVarList(glacial_melt_type  mtype,
                                           sv_type *aSV,
                                           int     *aLev,
                                           int     &nSV);
};

////////////////////////////////////////////////////////////////////
/// \brief Methods of modeling glacier water release
//
enum glacial_release_type{
  GRELEASE_LINEAR,         ///< Linear storage method for simulating glacial release
  GRELEASE_LINEAR_ANALYTIC,///< Linear storage method for simulating glacial release, analytical sol'n over timestep
  GRELEASE_HBV_EC          ///< HBV-EC method for simulating glacial release
};
////////////////////////////////////////////////////////////////////
/// \brief Data abstraction for release of glacier meltwater storage
//
class CmvGlacierRelease: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  glacial_release_type type; ///< Specified algorithm for glacial meltwater release simulation

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvGlacierRelease(const glacial_release_type r_type);
  ~CmvGlacierRelease();

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

  void        GetParticipatingParamList   (string  *aP, class_type *aPC, int &nP) const;
  static void GetParticipatingStateVarList(glacial_release_type  mtype,
                                           sv_type *aSV,
                                           int     *aLev,
                                           int     &nSV);
};

////////////////////////////////////////////////////////////////////
/// \brief Methods of modeling glacier water infiltration to deep groundwater
//
enum glacial_infil_type{
  GINFIL_UBCWM     ///< UBCWM (BC Hydro version) method for simulating infiltration of rain & meltwater
};
////////////////////////////////////////////////////////////////////
/// \brief Data abstraction for glacier water infiltration to deep groundwater
//
class CmvGlacierInfil: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  glacial_infil_type type; ///< Specified algorithm for glacial infiltration simulation

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvGlacierInfil(const glacial_infil_type i_type);
  ~CmvGlacierInfil();

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

  void        GetParticipatingParamList   (string  *aP, class_type *aPC, int &nP) const;
  static void GetParticipatingStateVarList(glacial_infil_type   mtype,
                                           sv_type *aSV,
                                           int     *aLev,
                                           int     &nSV);
};
#endif
