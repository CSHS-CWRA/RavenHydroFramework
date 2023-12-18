/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2018 the Raven Development Team
  ----------------------------------------------------------------
  class definitions:
  CmvPrairieBlowingSnow
  ----------------------------------------------------------------*/

#ifndef PBSM_H
#define PBSM_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"
#include "Model.h"

///////////////////////////////////////////////////////////////////
/// \brief Techniques for modelling blowing snow
//
enum pbsm_type
{
  PBSM_SUBLONLY,      ///< Pomeroy blowing snow sublimation (no lateral transfer)
  PBSM_FULL          ///< Pomeroy Prairie Blowing Snow Model with lateral transfer
};

///////////////////////////////////////////////////////////////////
/// \brief Calculates loss of water from snow to atmosphere
//
class CmvPrairieBlowingSnow: public CLateralExchangeProcessABC
{
private:/*------------------------------------------------------*/

  pbsm_type        _type;              ///< blowing snow algorithm type
  int              _nDriftConnections; ///< number of drifting snow corrections (nHRU_p*nHRU_p)*nSubBasins

  double SublimRateCoefficient(const double &Mpr,
                               const double &alpha,
                               const double &Vsalt,
                               const double &rel_hum_z,
                               const double &T) const;

  double ProbabilityThreshold(const double &snow_depth, //snow depth, [m]
                              const double &T,          //air temperature [deg C]
                              const double &snowfall,   //[mm/d]
                              const double &Uten_Prob,  //[m/s]
                                    double &wind_thresh,//Threshold wind speed [m/s]
                                    double &snow_age,   //snow age [d]
                              const double &tstep) const;

  void  PBSMrates           ( const double E_StubHt, // [m] stubble height
                              const double Uthr,     // [m/s] threshold wind speed
                              const double T,        // air temperature [deg C]
                              const double u,        // wind speed [m/s]
                              const double rel_hum,  // relative humidity [0..1]
                              const double Fetch,    // [m] (param)
                              const double veg_dens, // [count/m2] Vegetation density
                              const double veg_diam, // [m] Vegetation diameter
                              const double mBeta,    // [-] unitless parameter
                                    double &DriftH,  //[kg/m/s]
                                    double &SublH) const;  //[kg/m2/s] //per half hour?

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvPrairieBlowingSnow(pbsm_type sub_type, CModel *pModel);
  ~CmvPrairieBlowingSnow();

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
  static void GetParticipatingStateVarList(pbsm_type     stype,
                                           sv_type *aSV, int *aLev, int &nSV);

  void GetLateralExchange(const double * const *state_vars, //array of all SVs for all HRUs, [k][i]
                          const CHydroUnit * const *pHRUs,
                          const optStruct   &Options,
                          const time_struct &tt,
                                double      *exchange_rates) const;//purely virtual
};

#endif
