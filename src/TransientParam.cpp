/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2023 the Raven Development Team
  ----------------------------------------------------------------*/
#include "TransientParam.h"
#include "SoilAndLandClasses.h"
#include "Model.h"

/*****************************************************************
   Constructor / Destructor
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the transient parameter constructor
/// \param pTS [in] time series of parameter history
/// \param pname [in] name of parameter (e.g., "HYDRAULIC_COND")
/// \param ptype [in] class of parameter(e.g., CLASS_SOIL or CLASS_GLOBAL)
/// \param classname [in] class of parameter time series (e.g., "GuelphLoam") (not used/ignored for CLASS_GLOBAL parameters)
//
CTransientParam::CTransientParam(CTimeSeries      *pTS,
                                 const string      pname,
                                 const class_type  ptype,
                                 const string      classname,
                                 CModel *pModel)
{
  pTimeSeries=pTS;
  param_name=pname;
  param_type=ptype;
  class_name="";
  if (param_type!=CLASS_GLOBAL){
    class_name=classname;
  }
  _pModel = pModel;
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
//
CTransientParam::~CTransientParam()
{
  if (DESTRUCTOR_DEBUG){cout<<"  DELETING TRANSIENT PARAMETER"<<endl;}
  delete pTimeSeries; pTimeSeries=NULL;
}

//////////////////////////////////////////////////////////////////
/// \returns associated time series
//
const CTimeSeries *CTransientParam::GetTimeSeries(){return pTimeSeries;}

//////////////////////////////////////////////////////////////////
/// \returns associated parameter name
//
string CTransientParam::GetParameterName(){return param_name;}

//////////////////////////////////////////////////////////////////
/// \returns associated parameter class
//
string CTransientParam::GetParameterClass(){return class_name;}

//////////////////////////////////////////////////////////////////
/// \returns associated parameter class type (e.g., CLASS_SOIL)
//
class_type CTransientParam::GetParameterClassType(){return param_type;}

//////////////////////////////////////////////////////////////////
/// \brief Initialization routine, called prior to simulation
//
void CTransientParam::Initialize(const CModel *pModel, const optStruct &Options)
{
  pTimeSeries->Initialize(Options.julian_start_day,
                          Options.julian_start_year,
                          Options.duration,
                          Options.timestep,
                          false,
                          Options.calendar);

  if (param_type==CLASS_SOIL)
  {
    if (_pModel->StringToSoilClass(class_name) == NULL){
      string msg="CTransientParam::Initialize: invalid soil class name: "+class_name;
      ExitGracefully(msg.c_str(),BAD_DATA);
    }
  }
  else if (param_type==CLASS_VEGETATION)
  {
    if (_pModel->StringToVegClass(class_name) == NULL){
      string msg="CTransientParam::Initialize: invalid vegetation class name: "+class_name;
      ExitGracefully(msg.c_str(),BAD_DATA);
    }
  }
  else if (param_type==CLASS_TERRAIN)
  {
    if (_pModel->StringToTerrainClass(class_name)==NULL){
      string msg="CTransientParam::Initialize: invalid terrain class name: "+class_name;
      ExitGracefully(msg.c_str(),BAD_DATA);
    }
  }
  else if (param_type == CLASS_LANDUSE)
  {
    if (_pModel->StringToLUClass(class_name) == NULL) {
      string msg = "CTransientParam::Initialize: invalid land use/land type class name: " + class_name;
      ExitGracefully(msg.c_str(), BAD_DATA);
    }
  }
  else if (param_type==CLASS_GLOBAL)
  {
    //
  }

}
