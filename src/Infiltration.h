/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ----------------------------------------------------------------
  class definitions:
  CmvInfiltration
  ----------------------------------------------------------------*/

#ifndef INFILTRATION_H
#define INFILTRATION_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"


///////////////////////////////////////////////////////////////////
/// \brief Algorithms for modeling infiltration

//
enum infil_type{
  INF_SCS,                      ///< SCS curve number method for infiltration
  INF_SCS_NOABSTRACTION,        ///< SCS curve number method for infiltration without abstraction representation
  INF_RATIONAL,                 ///< Rational method for infiltration
  INF_ALL_INFILTRATES,          ///< Simple infiltration method where all rainfall infiltrates into soil
  INF_GREEN_AMPT,               ///< Standard Green-Ampt infiltration model
  INF_GA_SIMPLE,                ///< Simple Green-Ampt infiltration model
  INF_UPSCALED_GREEN_AMPT,	///< Upscaled Green-Ampt infiltration model from Craig (2010)
  INF_PHILIP_1957,              ///< Philip 1957 Infiltration model
  INF_VIC,                      ///< VIC (Variable Infiltration Capacity) Infiltration model
  INF_VIC_ARNO,                 ///< VIC-ARNO Infiltration model
  INF_TOPMODEL,                 ///< TOPMODEL Infiltration model
  INF_PRMS,                     ///< PRMS Infiltration model
  INF_HBV,                      ///< HBV Infiltration model
  INF_UBC,                      ///< UBC Watershed model infiltration algorithm
  INF_GR4J,                     ///< GR4J Watershed model infiltration algorithm
  INF_HMETS,                    ///< HMETS model infiltration algorithm
  INF_XINANXIANG,               ///< Xinanxiang model infiltration algorithm
  INF_PDM,                      ///< Probability Distributed Model (PDM) of Moore (1985) [HYMOD]
  INF_AWBM                      ///< Austrailian Water Balance Model of Boughton (1993)
};

////////////////////////////////////////////////////////////////////
/// \Brief Data abstraction for infiltration process
/// \details partitions ponded water (rainfall/snowmelt) to either soil storage or runoff
//
class CmvInfiltration: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  infil_type type; ///< Infiltration algorithm

  static double GreenAmptCumInf   (const double &t,         //[d], time from start of rainfall (or time step)
                                   const double &alpha,     //[mm]
                                   const double &Ks,        //Ksat [mm/d]
                                   const double &w);  //rainfall rate, [mm/d]

  double GetSCSRunoff      (const CHydroUnit  *pHRU,
                            const optStruct       &Options,
                            const time_struct &tt,
                            const double      &rainthru) const;

  void   GetUBCWMRunoff    (const double                  *state_vars,
                            const CHydroUnit  *pHRU,
                            const optStruct       &Options,
                            const time_struct &tt,
                            double      *rates,
                            const double      &rainthru) const;

  void GetGreenAmptRunoff(const double      *state_vars,
                          const CHydroUnit  *pHRU,
                          const optStruct   &Options,
                          const time_struct &tt,
                          double            *rates,
                          const double      &rainthru) const;

  double GetHeterogeneousGreenAmptRunoff(
                          const double      &rainthru,//[mm/d]
                          const soil_struct *S,
                          const double       soil_thickness,//[m]
                          const double      &soil_water_content,//[mm]
                          const time_struct &tt,//[d]
                          const optStruct   &Options) const;

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvInfiltration(const infil_type itype);
  ~CmvInfiltration();

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
  void        GetParticipatingParamList   (string  *aP , class_type *aPC , int &nP) const;
  static void GetParticipatingStateVarList(infil_type btype,sv_type *aSV, int *aLev, int &nSV);
};

#endif
