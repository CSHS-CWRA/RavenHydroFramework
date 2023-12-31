/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2022 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Properties.h"
#include "SoilAndLandClasses.h"
#include "Model.h"

/*****************************************************************
   Constructor / Destructor
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the terrain class constructor
/// \param name [in] String nickname for terrain class
//
CTerrainClass::CTerrainClass(const string name, CModel* pModel)
{
  tag = name;
  pModel->AddTerrainClass(this);
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
//
CTerrainClass::~CTerrainClass()
{
  if (DESTRUCTOR_DEBUG){cout<<"  DELETING TERRAIN CLASS "<<endl;}
}

///////////////////////////////////////////////////////////////////
/// \brief Returns pointer to terrain properties
/// \return pointer to Terrain properties associated with terrain class
//
const terrain_struct *CTerrainClass::GetTerrainStruct() const{return &T;}

///////////////////////////////////////////////////////////////////
/// \brief Return nick name identifier of terrain class
/// \return nick name identifier of terrain class
//
string                CTerrainClass::GetTag                                      () const{return tag;}
/*****************************************************************
   Static Initialization, Accessors
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Automatically calculates terrain propeties
/// \details Sets terrain properties based upon simple terrain parameters
///     Input [Ttmp] has been read from .rvp file - if parameter ==
///     AUTO_COMPLETE, then empirical relationships are used to estimate
///     parameters
///
/// \param &Ttmp [in] Input terrain class parameters (read from .rvp file)
/// \param &Tdefault [in] Default terrain class parameters
//
void CTerrainClass::AutoCalculateTerrainProps(const terrain_struct &Ttmp,
                                              const terrain_struct &Tdefault)
{
  //bool autocalc=false;

  //these parameters are required
  T.drainage_density =Ttmp.drainage_density;
  T.hillslope_length =Ttmp.hillslope_length;

  //Standard terrain properties
  //----------------------------------------------------------------------------
  /*autocalc=SetCalculableValue(T.lambda,Ttmp.lambda,Tdefault.lambda);
    if (autocalc)
    {
    T.topmodel_lambda=0.0;
    }*/

  //Model-specific terrain properties - cannot be autocomputed, must be specified by user
  //----------------------------------------------------------------------------
  SetSpecifiedValue(T.topmodel_lambda,Ttmp.topmodel_lambda,Tdefault.topmodel_lambda,false,"TOPMODEL_LAMBDA");
}

//////////////////////////////////////////////////////////////////
/// \brief Sets default terrain properties
/// \param &T [out] Terrain properties class
/// \param is_template [in] True if the default value being set is for the template class
//
void CTerrainClass::InitializeTerrainProperties(terrain_struct &T, bool is_template)//static
{
  //required parameters
  T.hillslope_length=100;//[m]
  T.drainage_density=1.0;//[?]
  T.topmodel_lambda =7.5;//[m]needs reasonable estimate

  //T.topmodel_lambda =DefaultParameterValue(is_template,true);
}
//////////////////////////////////////////////////////////////////
/// \brief Sets the value of the terrain property corresponding to param_name
/// \param param_name [in] Parameter identifier
/// \param value [in] Value of parameter to be set
//
void  CTerrainClass::SetTerrainProperty(const string &param_name,
                                        const double &value)
{
  SetTerrainProperty(T,param_name,value);
}
//////////////////////////////////////////////////////////////////
/// \brief Sets the value of the terrain property corresponding to param_name
/// \param &T [out] Terrain properties class
/// \param param_name [in] Parameter identifier
/// \param value [in] Value of parameter to be set
//
void  CTerrainClass::SetTerrainProperty(terrain_struct &T,
                                        const string param_name,
                                        const double value)
{
  string name;
  name = StringToUppercase(param_name);

  if      (!name.compare("HILLSLOPE_LENGTH"  )){T.hillslope_length=value;}
  else if (!name.compare("DRAINAGE_DENSITY"  )){T.drainage_density=value;}
  else if (!name.compare("TOPMODEL_LAMBDA"   )){T.topmodel_lambda=value; }
  else{
    WriteWarning("CTerrainClass::SetTerrainProperty: Unrecognized/invalid terrain parameter name ("+name+") in .rvp file",false);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief gets terrain property corresponding to param_name
/// \param param_name [in] Parameter identifier
/// \returns value of parameter
//
double CTerrainClass::GetTerrainProperty(string param_name) const
{
  return GetTerrainProperty(T,param_name);
}

///////////////////////////////////////////////////////////////////////////
/// \brief Returns terrain property value corresponding to param_name from structure provided
/// \param &T [in] terrain structure
/// \param param_name [in] Parameter name
/// \return terrain property corresponding to parameter name
//
double CTerrainClass::GetTerrainProperty(const terrain_struct &T, string param_name)
{
  string name;
  name = StringToUppercase(param_name);

  if      (!name.compare("HILLSLOPE_LENGTH"       )){return T.hillslope_length;}
  else if (!name.compare("DRAINAGE_DENSITY"       )){return T.drainage_density;}
  else if (!name.compare("TOPMODEL_LAMBDA"        )){return T.topmodel_lambda; }
  else{
    string msg="CTerrainClass::GetTerrainProperty: Unrecognized/invalid terrain parameter name in .rvp file: "+name;
    ExitGracefully(msg.c_str(),BAD_DATA_WARN);
    return 0.0;
  }

}
