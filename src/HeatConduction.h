/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2023 the Raven Development Team
----------------------------------------------------------------
class definitions: CmvHeatConduction
----------------------------------------------------------------*/

#ifndef HEATCONDUCTION_H
#define HEATCONDUCTION_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"
#include "Transport.h"

////////////////////////////////////////////////////////////////////
/// \brief Data abstraction for loss of water from soil/groundwater to surface water
//
class CmvHeatConduction : public CHydroProcessABC
{
private:/*------------------------------------------------------*/

  const CTransportModel *_pTransModel;
  int                    _nHRUs;   //must store locally to retain start-of-timestep water storage

  bool                   _initialized;

  bool GenerateJacobianMatrix(const double  *z,
                              const double  *eta,const double * poro,const double *kappa_s,
                              const double  *sat,const double *satn,
                              const double  *Vold,const double *Vnew, //[m]
                              const double  *hold,const double *h,
                              const double  &tstep,
                                    double **J,
                                    double  *f,
                              const int      N) const;

public:/*-------------------------------------------------------*/
       //Constructors/destructors:
  CmvHeatConduction(const CTransportModel *pTransMod,
                    CModelABC             *pModel);
  ~CmvHeatConduction();

  void StoreNumberOfHRUs(const int nHRUs);

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

  void GetParticipatingParamList(string  *aP,class_type *aPC,int &nP) const;
  void GetParticipatingStateVarList(sv_type *aSV,int *aLev,int &nSV);

};

#endif
