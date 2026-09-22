/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2026 the Raven Development Team
  ----------------------------------------------------------------
  struct class_change
  class  CTransientParam
  ----------------------------------------------------------------*/
#ifndef TRANSIENT_PARAM_H
#define TRANSIENT_PARAM_H

#include "RavenInclude.h"
#include "Properties.h"
#include "TimeSeries.h"
#include "SoilAndLandClasses.h"
#include "Model.h"

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for abrupt class changes
//
struct class_change
{
  int        HRU_groupID; // HRU group id (kk)
  class_type tclass;      // type of class (e.g., CLASS_LANDUSE)
  string     newclass;    // new class tag
  double     modeltime;   // modeltime of shift
};
///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for continuous class transitions
//
struct class_transition
{
  int        HRU_groupID;    // HRU group id (kk)
  class_type tclass;         // type of class (e.g., CLASS_LANDUSE)
  string     newclass;       // new class tag
  double     starttime;      // modeltime of transition start
  double     endtime;        // modeltime of transititon end
  transition_function funct; // transition function
  double     params[3];      // parameters of transition function

  surface_struct surf_params; //HRUs point to these, not class structures, during transition
  veg_struct     veg_params;
};
///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for global model parameters
//
class CTransientParam
{
protected:/*----------------------------------------------------*/

  string      param_name;   ///< Name of parameter (e.g., "HYDRAULIC_COND")
  class_type  param_type;   ///< parameter class type (e.g., SOIL, TERRAIN,LANDUSE,GLOBAL)
  string      class_name;   ///< class name (e.g., "GuelphLoam")(ignored for global params)

  CTimeSeries *pTimeSeries; ///< Time series of parameter value
  CModel      *_pModel;

public:/*-------------------------------------------------------*/

  CTransientParam(      CTimeSeries      *pTS,
                        const string      pname,
                        const class_type  ptype,
                        const string      classname,
                        CModel           *pModel);
  ~CTransientParam();

  //Accessors
  const CTimeSeries *GetTimeSeries();
  string             GetParameterName();
  string             GetParameterClass();
  class_type         GetParameterClassType();

  //routines
  void Initialize(const CModel *pModel, const optStruct &Options);

};

#endif
