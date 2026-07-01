/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2026 the Raven Development Team
----------------------------------------------------------------*/

#ifndef AGE_TRACERS_H
#define AGE_TRACERS_H

#include "RavenInclude.h"
#include "Transport.h"

///////////////////////////////////////////////////////////////////
/// \brief Class for coordinating age tracer simulation (child of ConstituentModel)
/// \details Implemented in AgeTracers.cpp
//
class CAgeTracer :public CConstituentModel
{
private:
   double _model_time; //local storage of tt.model_time
   double _t_diff;

public:/*-------------------------------------------------------*/
  CAgeTracer(CModel *pMod,CTransportModel *pTMod,string name,const int c);
  ~CAgeTracer();

  double CalculateReportingConcentration(const double &mass, const double &vol) const; //returns age [d]
  double ConvertConcentration(const double &delta) const;
  double GetOutflowConcentration(const int p) const;

  double GetAdvectionCorrection(const CHydroUnit* pHRU,const int iFromWater,const int iToWater,const double& mass, const double &vol, const double &Q) const;

  void   SetModelTime(const double &t);
  void   SetInitialTimeOffset(const double &tdiff);

  //Manipulators (inherited from CConstitModel)
  void   Initialize                 (const optStruct &Options);
  void   UpdateInitialConditions    (const optStruct &Options);
  void   WriteEnsimOutputFileHeaders(const optStruct &Options);
  void   WriteEnsimMinorOutput      (const optStruct &Options,const time_struct &tt);
};
#endif
