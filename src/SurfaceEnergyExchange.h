/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2021 the Raven Development Team
------------------------------------------------------------------
  definitions of surface energy transfer classes
----------------------------------------------------------------*/
#ifndef SURFACEEXCHANGE_H
#define SURFACEEXCHANGE_H

#include "RavenInclude.h"
#include "Model.h"
#include "Transport.h"
#include "EnergyTransport.h"

////////////////////////////////////////////////////////////////////
/// \brief Data abstraction for transfer of energy from atmosphere to land surface (canopy/soil/snow/depression)
//
class CmvPartitionEnergy : public CHydroProcessABC
{
private:/*------------------------------------------------------*/

  const CTransportModel *_pTransModel;

  double GetStorTemp(int iEnth,int iWatStor,const double &Tair,const double *state_vars) const;

public:/*-------------------------------------------------------*/
       //Constructors/destructors:
  CmvPartitionEnergy(const CTransportModel *pTransMod, CModelABC *pModel);
  ~CmvPartitionEnergy();

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
                        const time_struct &tt,
                              double      *rates) const;

  void        GetParticipatingParamList(string  *aP,class_type *aPC,int &nP) const;
  static void GetParticipatingStateVarList(sv_type *aSV,int *aLev,int &nSV);

};

#endif
