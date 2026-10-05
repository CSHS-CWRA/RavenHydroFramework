/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2026 the Raven Development Team
----------------------------------------------------------------
Irrigation input process
----------------------------------------------------------------*/
#ifndef IRRIGATION_H
#define IRRIGATION_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"

///////////////////////////////////////////////////////////////////
enum irrigation_type
{
  IRRIG_FROMFILE   ///< irrigation supplied from IRRIGATION forcing time series/grid
};

////////////////////////////////////////////////////////////////////
/// \brief Moves irrigation water (IRRIGATION forcing) from ATMOS_PRECIP through the IRRIGATION_SRC storage
//
class CmvIrrigation: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  irrigation_type _type;   ///< irrigation method

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvIrrigation(irrigation_type itype,
                int             to_index,
                CModelABC      *pModel);
  ~CmvIrrigation();

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
                        const time_struct &tt,
                        double            *rates) const;

  void        GetParticipatingParamList   (string *aP, class_type *aPC, int &nP) const{nP=0;}
  static void GetParticipatingStateVarList(irrigation_type itype, sv_type *aSV, int *aLev, int &nSV);
};
#endif
