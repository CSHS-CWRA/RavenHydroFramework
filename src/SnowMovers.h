/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ----------------------------------------------------------------
  class definitions:
  CmvSublimation
  CmvSnowMelt
  CmvSnowRefreeze
  CmvSnowBalance
  ----------------------------------------------------------------*/

#ifndef SNOWMOVERS_H
#define SNOWMOVERS_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"

///////////////////////////////////////////////////////////////////
/// \brief Techniques for modelling sublimation
//
enum sublimation_type
{
  SUBLIM_ALL,             ///< HBV model
  SUBLIM_MAXIMUM,         ///< sublimates at PET rate (not very realistic)
  SUBLIM_SVERDRUP,        ///< Sverdrup 1946 : adapted from Gray 1973 \cite Sverdrup1946JoAS \cite gray1974
  SUBLIM_KUZMIN,          ///< Kuzmin 1972 : adapted from Gray 1973
  SUBLIM_CENTRAL_SIERRA,  ///< US Army Corps of Engineers : adapted from Gray 1973 \cite Engineers1956NPDPO
  SUBLIM_PBSM,            ///< Pomeroy Prairie Blowing Snow Model Sublimation : adapted from Pomeroy 1993 \cite Pomeroy1993JoH
  SUBLIM_WILLIAMS,        ///< Canadian National Research Council [1959] - developed from Sverdrup
  SUBLIM_KUCHMENT_GELFAN, ///< Kuchment Gelfan (1996) - essentially equivalent to Kuzmin 
  SUBLIM_BULK_AERO        ///< Bulk aerodynamic approach of Hock (2005)
};

///////////////////////////////////////////////////////////////////
/// \brief Calculates loss of water from snow to atmosphere
//
class CmvSublimation: public CHydroProcessABC
{
private:/*------------------------------------------------------*/

  sublimation_type      type; ///< sublimation algorithm type
  bool              model_depth_change; ///< True if depth change is modeled

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvSublimation(sublimation_type sub_type);
  ~CmvSublimation();

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
  static void GetParticipatingStateVarList(sublimation_type     stype,
                                           sv_type *aSV, int *aLev, int &nSV);
};

///////////////////////////////////////////////////////////////////
/// \brief Techniques for modelling water loss from snow to soil/runoff
//
enum snowmelt_type
{
  MELT_POTMELT               //< melt at potential melt rate
};

////////////////////////////////////////////////////////////////////
/// \brief Calculates rate of snowmelt
//
class CmvSnowMelt: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  snowmelt_type type; ///< Model of snow melt selected

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvSnowMelt(snowmelt_type melt_type,
              int Out_index);
  ~CmvSnowMelt();

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
  static void GetParticipatingStateVarList(snowmelt_type stype,
                                           sv_type *aSV, int *aLev, int &nSV);
};

///////////////////////////////////////////////////////////////////
/// \brief Calculates loss of water from liquid snow to ponded water/surface water as maximum holding capacity of snow changes
//
class CmvSnowSqueeze: public CHydroProcessABC
{
private:/*------------------------------------------------------*/

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvSnowSqueeze(int out_index);
  ~CmvSnowSqueeze();

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
  static void GetParticipatingStateVarList(sv_type *aSV, int *aLev, int &nSV);
};

///////////////////////////////////////////////////////////////////
/// \brief Techniques for modelling snow refreeze
//
enum refreeze_type
{
  FREEZE_DEGREE_DAY///< Degree day method
};

///////////////////////////////////////////////////////////////////
/// \brief Calculates loss of water from snow to soil/runoff
//
class CmvSnowRefreeze: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  refreeze_type type; ///< Model of snow refreeze selected

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvSnowRefreeze(refreeze_type freeze_type);
  ~CmvSnowRefreeze();

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
  static void GetParticipatingStateVarList(refreeze_type type,
                                           sv_type *aSV, int *aLev, int &nSV);
};

///////////////////////////////////////////////////////////////////
/// \brief Methods of balancing energy based on melt/refreeze
//
enum snowbal_type
{
  SNOBAL_SIMPLE_MELT,  ///< melt rate = potential melt, only SWE modified
  SNOBAL_COLD_CONTENT, ///< Cold content-based
  SNOBAL_HBV,          ///< HBV method
  SNOBAL_UBCWM,        ///< UBCWM method
  SNOBAL_CEMA_NEIGE,   ///< Cema Neige method
  SNOBAL_TWO_LAYER,    ///< Two layer cold content (converted from GJ C# code)
  SNOBAL_GAWSER,       ///< GAWSER snow melt model (modified from Object GAWSER to replicate behavior)
  SNOBAL_CRHM_EBSM,    ///< CRHM's energy balance snow model (from Marks, 1997)
  SNOBAL_HMETS         ///< HMETS single layer snow balance routine (Martel et al., 2017)
};

///////////////////////////////////////////////////////////////////
/// \brief Calculates melt,refreeze, density change, etc. based upon mass & energy balance on snowpack
//
class CmvSnowBalance: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  snowbal_type type; ///< Method of energy-balance selected

  void ColdContentBalance(const double           *state_vars,
                            const CHydroUnit       *pHRU,
                            const optStruct        &Options,
                            const time_struct      &t,
                            double           *rates) const;

  void TwoLayerBalance(const double           *state_vars,
                            const CHydroUnit       *pHRU,
                            const optStruct        &Options,
                            const time_struct      &t,
                            double           *rates) const;

  void GawserBalance    (const double           *state_vars,
                            const CHydroUnit       *pHRU,
                            const optStruct        &Options,
                            const time_struct      &t,
                            double           *rates) const;
  void CRHMSnowBalance  (const double           *state_vars,
                            const CHydroUnit       *pHRU,
                            const optStruct        &Options,
                            const time_struct      &t,
                                  double           *rates) const;

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvSnowBalance(snowbal_type  bal_type);
  CmvSnowBalance(snowbal_type bal_type, int iSnowTo);
  ~CmvSnowBalance();

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
  static void GetParticipatingStateVarList(snowbal_type stype,
                                           sv_type *aSV, int *aLev, int &nSV);
};
///////////////////////////////////////////////////////////////////
/// \brief Methods of updating snow temperature
//
enum snowtemp_evolve_type
{
  SNOTEMP_NEWTONS  ///< Cema Neige approach - linear exchange with atmosphere
};

///////////////////////////////////////////////////////////////////
/// \brief Calculates evolution of snow temperature based upon air temperature
//
class CmvSnowTempEvolve: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  snowtemp_evolve_type _type; ///< Method of snow temp evolution

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvSnowTempEvolve(snowtemp_evolve_type  ste_type);
  ~CmvSnowTempEvolve();

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
  static void GetParticipatingStateVarList(snowtemp_evolve_type ste_type,
                                           sv_type *aSV, int *aLev, int &nSV);
};
#endif
