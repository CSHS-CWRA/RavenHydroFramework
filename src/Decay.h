/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2018 the Raven Development Team
  ----------------------------------------------------------------
  class definitions:
  CmvDecay
  CmvTransformation
  ----------------------------------------------------------------*/

#ifndef DECAY_H
#define DECAY_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"
#include "Model.h"
enum decay_type
{
  DECAY_BASIC,    /// < basic decay, calculated as -k*C
  DECAY_ANALYTIC  /// < analytic treatment of decay over finite time step
};

////////////////////////////////////////////////////////////////////
/// \brief Calculates the decay of a constituent
//
class CmvDecay: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  const CTransportModel* _pTransModel;

  decay_type _dtype;              ///< decay algorithm type
  int _constit_ind;              ///< index of constituent which is decaying

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvDecay(string constit_name, decay_type dtyp, CTransportModel *pTransportModel);
  ~CmvDecay();

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

};

enum transformation_type
{
  TRANSFORM_LINEAR,            /// < basic linear transformation, calculated as dC/dt=-k*C; dA/dt=+k*C
  TRANSFORM_LINEAR_ANALYTIC,   /// < analytic treatment of decay over finite time step
  TRANSFORM_NONLINEAR,         /// < basic power law transformation, calculated as dC/dt=-k*C^n; dA/dt=+k*C^n
  TRANSFORM_NONLIN_ANALYTIC,   /// < analytic treatment of nonlinear transformation over finite time step
  TRANSFORM_MINERALIZATION     /// < Nitrate mineralization 
};
////////////////////////////////////////////////////////////////////
/// \brief Calculates the decay of a substance
//
class CmvTransformation: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  const CTransportModel* _pTransModel;

  transformation_type _ttype;  ///< transformation algorithm type
  int _constit_ind1;           ///< index of reactant constituent 
  int _constit_ind2;           ///< index of product constituent 

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvTransformation(string reactant_name, string product_name, transformation_type ttyp, CTransportModel *pTransportModel);
  ~CmvTransformation();

  //inherited functions
  void Initialize();
  void GetRatesOfChange(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double            *rates) const;
  void ApplyConstraints(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &t,
                        double            *rates) const;

  void        GetParticipatingParamList   (string  *aP, class_type *aPC, int &nP) const;

};
#endif
