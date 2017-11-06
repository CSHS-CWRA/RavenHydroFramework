/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2017 the Raven Development Team
----------------------------------------------------------------
class definitions:
CmvLatAdvection
----------------------------------------------------------------*/

#ifndef LATADVECTION_H
#define LATADVECTION_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"
#include "LateralExchangeABC.h"
#include "Model.h"

////////////////////////////////////////////////////////////////////
/// \brief Calculates the advection of a substance
//
class CmvLatAdvection : public CLateralExchangeProcessABC
{
private:/*------------------------------------------------------*/
  const CTransportModel* pTransModel;

  int      _constit_ind;         ///< index of constituent being advected

  int         _nSources;         ///< number of Dirichlet source compartments
  int         *_iSource;         ///< state variable indices of Dirichlet source compartments [size: nSources]
  double *_aSourceConcs;         ///< array of source concentrations for Dirichlet sources [size: nSources]

public:/*-------------------------------------------------------*/
       //Constructors/destructors:
  CmvLatAdvection(string constit_name,CTransportModel *pTransportModel);
  ~CmvLatAdvection();

  //inherited functions
  void Initialize();
  void GetRatesOfChange(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                              double      *rates) const{}//does nothing
  void ApplyConstraints(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &t,
                              double      *rates) const{}//does nothing

  void        GetParticipatingParamList(string  *aP,class_type *aPC,int &nP) const;

  void GetLateralExchange(const double * const *state_vars, //array of all SVs for all HRUs, [k][i]
                          const CHydroUnit * const *pHRUs,
                          const optStruct   &Options,
                          const time_struct &tt,
                                double      *exchange_rates) const;

};
#endif
