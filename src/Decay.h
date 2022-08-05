/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2022 the Raven Development Team
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
  DECAY_LINEAR,   /// < analytic treatment of first-order decay/loss over finite time step
  DECAY_ZEROORDER,/// < analytic treatment of zeroth-order decay/loss over finite time step
  DECAY_DENITRIF, /// < adds temperature/saturation correction factor to decay rate
  DECAY_HOOKESLAW /// < exchange with atmosphere 
};

////////////////////////////////////////////////////////////////////
/// \brief Calculates the decay of a constituent
//
class CmvDecay: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  const CTransportModel* _pTransModel;

  decay_type _dtype;              ///< decay algorithm type
  int        _constit_ind;        ///< index, c, of constituent which is decaying

  int        _process_ind;        ///< index of transport process name
  int        _iWaterStore;        ///< index ii of water store (or DOESNT_EXIST if this should occur in all water compartments)

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvDecay(string constit_name, decay_type dtyp,int proc_ind,int iWatStor,CTransportModel *pTransportModel);
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
  TRANS_LINEAR,            /// < analytic treatment of decay over finite time step, dC/dt=-k*C dA/dt=+s*k*C
  TRANS_NONLINEAR         /// < basic power law transformation, analytic treatment  calculated as dC/dt=-k*C^n; dA/dt=+s*k*C^n
};
////////////////////////////////////////////////////////////////////
/// \brief Calculates the decay of a substance in a particular water store
//
class CmvTransformation: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  const CTransportModel* _pTransModel;

  transformation_type _ttype;  ///< transformation algorithm type
  int _constit_ind1;           ///< index of reactant constituent 
  int _constit_ind2;           ///< index of product constituent 
  int _process_ind;            ///< index of transport process name
  int _iWaterStore;            ///< index of water store (or DOESNT_EXIST if this should occur in all water compartments)

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvTransformation(string reactant_name, string product_name, transformation_type ttyp,int proc_ind,int iWatStor, CTransportModel *pTransportModel);
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
