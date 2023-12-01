/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2022 the Raven Development Team
  ----------------------------------------------------------------
  class definitions:
  CmvChemEquil
  ----------------------------------------------------------------*/

#ifndef CHEM_EQUIL_H
#define CHEM_EQUIL_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"
#include "Model.h"
enum chem_equil_type
{
  EQUIL_FIXED_RATIO,     /// < constitutent 2 mass is a linear function of constituent 1 mass (B=cc*A)
  EQUIL_LINEAR,          /// < rate-limited linear equilibrium dB/dt=-k*(B-cc*A)
  EQUIL_LINEAR_SORPTION  /// < linear equilibrium sorption Cs=Kd*Cw (constit1=aqueous, constit2=sorbed)
};

////////////////////////////////////////////////////////////////////
/// \brief Calculates the decay of a constituent
//
class CmvChemEquil : public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  const CTransportModel* _pTransModel;

  chem_equil_type _eq_type;            ///< decay algorithm type
  int             _constit_ind1;       ///< index, c, of references constituent
  int             _constit_ind2;       ///< index, c, of constituent at equilibrium with c1

  int             _process_ind;        ///< index of equilibrium process name
  int             _iWaterStore;        ///< index ii of water store (or DOESNT_EXIST if this should occur in all water compartments)

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvChemEquil(string constit_name1,
               string constit_name2,
               chem_equil_type dtyp,
               int proc_ind,
               int iWatStor,
               CTransportModel *pTransportModel,
               CModelABC *pModel);
  ~CmvChemEquil();

  //inherited functions
  void Initialize();
  void GetRatesOfChange(const double     * state_vars,
                        const CHydroUnit * pHRU,
                        const optStruct  & Options,
                        const time_struct& tt,
                              double     * rates) const;
  void ApplyConstraints(const double     * state_vars,
                        const CHydroUnit * pHRU,
                        const optStruct  & Options,
                        const time_struct& t,
                              double     * rates) const;

  void        GetParticipatingParamList(string* aP,class_type* aPC,int& nP) const;

};
#endif
