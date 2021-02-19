/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ----------------------------------------------------------------
  class definitions:
  CmvDepressionOverflow
  CmvAbstraction
  ----------------------------------------------------------------*/

#ifndef DEPRESSION_PROCESSES_H
#define DEPRESSION_PROCESSES_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"

///////////////////////////////////////////////////////////////////
/// \brief Method of calculating overflow of depression storage into surface water
//
enum depflow_type{
  DFLOW_THRESHPOW, ///< threshold-based power law overflow
  DFLOW_LINEAR,    ///< threshold-based linear overflow
  DFLOW_WEIR       ///< threshold-based weir eqn overflow
};

///////////////////////////////////////////////////////////////////
/// \brief Calculates the loss depression storage into surface water
//
class CmvDepressionOverflow: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  depflow_type _type; ///< Model of abstaction selected

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvDepressionOverflow(const depflow_type dtype);
  ~CmvDepressionOverflow();

  //inherited functions
  void Initialize();
  void GetRatesOfChange(const double      *state_vars,
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
  static void GetParticipatingStateVarList(depflow_type  dtype,
                                           sv_type *aSV,
                                           int     *aLev,
                                           int     &nSV);
};

///////////////////////////////////////////////////////////////////
/// \brief Method of calculating seepage of depression storage into soil storage
//
enum seepage_type{
  SEEP_LINEAR    ///< linear seepage
};

///////////////////////////////////////////////////////////////////
/// \brief Calculates the seepage of depression storage into lower soils
//
class CmvSeepage: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  seepage_type _type; ///< Model of abstaction selected

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvSeepage(const seepage_type seeptype, int iToSoil);
  ~CmvSeepage();

  //inherited functions
  void Initialize();
  void GetRatesOfChange(const double      *state_vars,
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
  static void GetParticipatingStateVarList(seepage_type  dtype,
                                           sv_type *aSV,
                                           int     *aLev,
                                           int     &nSV);
};


///////////////////////////////////////////////////////////////////
/// \brief Method of calculating abstraction of rainfall/snowmelt
//
enum abstraction_type{
  ABST_PERCENTAGE, ///< Abstraction by percentage
  ABST_FILL,       ///< Fill abstraction
  ABST_SCS,        ///< Abstraction using SCS method
  ABST_PDMROF,     ///< Abstraction using PDMROF Method (Mekkonen et al., 2014)
  ABST_UWFS        ///< Abstraction using upscaled wetland fill and spill model (Taheri et al.)
};

///////////////////////////////////////////////////////////////////
/// \brief Calculates the abstraction of rainfall/snowmelt (ponded water) to depression storage
//
class CmvAbstraction: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  abstraction_type _type; ///< Model of abstaction selected

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvAbstraction(const abstraction_type absttype);
  ~CmvAbstraction();

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
  static void GetParticipatingStateVarList(abstraction_type  absttype,
                                           sv_type *aSV,
                                           int     *aLev,
                                           int     &nSV);
};
///////////////////////////////////////////////////////////////////
/// \brief Method of calculating release of lake storage into surface water
//
enum lakerel_type{
  LAKEREL_LINEAR   ///< linear loss/gain rate
};

///////////////////////////////////////////////////////////////////
/// \brief calculates release of lake storage into surface water
//
class CmvLakeRelease: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  lakerel_type _type; ///< Model of abstaction selected

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvLakeRelease(const lakerel_type absttype);
  ~CmvLakeRelease();

  //inherited functions
  void Initialize();
  void GetRatesOfChange(const double      *state_vars,
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
  static void GetParticipatingStateVarList(lakerel_type  lr_type,
                                           sv_type *aSV,
                                           int     *aLev,
                                           int     &nSV);
};
#endif
