/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2020 the Raven Development Team
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

public:/*-------------------------------------------------------*/
       //Constructors/destructors:
  CmvHeatConduction(const CTransportModel *pTransMod);
  ~CmvHeatConduction();

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

  void        GetParticipatingParamList(string  *aP,class_type *aPC,int &nP) const;
  static void GetParticipatingStateVarList(sv_type *aSV,int *aLev,int &nSV);

};

#endif