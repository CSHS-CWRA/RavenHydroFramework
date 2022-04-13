/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2022 the Raven Development Team
  ----------------------------------------------------------------
  Class CSoilClass
  Class CVegetationClass
  Class CLandUseClass
  Class CTerrainClass
  ----------------------------------------------------------------*/
#ifndef SOIL_LAND_CLASSES_H
#define SOIL_LAND_CLASSES_H

#include "RavenInclude.h"
#include "Properties.h"

class CHydroUnit;
class CModelABC;

////////////////////////////////////////////////////////////////////
/// \brief Types of class
//
enum class_type{
  CLASS_SOIL,       ///< Soil class
  CLASS_VEGETATION, ///< Vegetation class
  CLASS_LANDUSE,    ///< Landuse class
  CLASS_TERRAIN,    ///< Terrain class
  CLASS_TRANSPORT,  ///< Transport class
  CLASS_GLOBAL,     ///< Global parameter class
  CLASS_SUBBASIN,   ///< Subbasin parameter class
  CLASS_GAUGE,      ///< Gauge parameter class
  CLASS_HRUTYPE,    ///< special flag to handle HRU type changes; not assoc. with parameters
  CLASS_UNKNOWN
};

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for soil classification
//
class CSoilClass
{
protected:/*----------------------------------------------------*/

  string              _tag;                 ///< nickname for soil, e.g., "SILTY_SAND"
  soil_struct         _Soil;                ///< corresponding properties for soil

  static CSoilClass **_pAllSoilClasses;     ///< array of pointers to all soil classes that have been created
  static int          _nAllSoilClasses;     ///< number of soil classes

public:/*-------------------------------------------------------*/
  //Constructors:
  CSoilClass(const string name, const int nConstits);
  ~CSoilClass();

  //Accessors
  string                   GetTag() const;
  const soil_struct       *GetSoilStruct() const;
  double                   GetSoilProperty(string &param_name) const;
  void                     SetSoilProperty(string param_name, const double &value);

  //routines
  void AutoCalculateSoilProps(const soil_struct &Stmp,const soil_struct &Sdefault,const int nConstit);

  static int               GetNumClasses();
  static const CSoilClass *GetSoilClass(int c);
  static       CSoilClass *StringToSoilClass(const string s);
  static void              DestroyAllSoilClasses();
  static void              SetSoilProperty         (soil_struct &S, string param_name, const double value);
  static double            GetSoilProperty         (const soil_struct &S, string param_name, const bool strict=true);
  static void              InitializeSoilProperties(soil_struct &S, bool is_template,int nConstits);

  static double            CalcSoilResistance(const double &psi,const soil_struct &S);
  static double            CalcSoilHydCond   (const double      &sat_liq,
                                              const double      &sat_ice,
                                              const soil_struct &S);
  static double            CalcSoilSaturation(const double      &psi,
                                              const soil_struct &S);
  static double            CalcSoilPotential (const double      &sat_liq,
                                              const double      &sat_ice,
                                              const soil_struct &S);

  static void              SummarizeToScreen();
};

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for vegetation classification
//
class CVegetationClass
{
protected:/*----------------------------------------------------*/

  veg_struct                V;                    ///< corresponding canopy/root properties

  static CVegetationClass **pAllVegClasses;       ///< array of pointers to all vegetation classes that have been created
  static int                NumVegClasses;        ///< Number of vegetation classes that have been created (length of pAllVegClasses array)

public:/*-------------------------------------------------------*/
  //Constructors:
  CVegetationClass(const string name);
  ~CVegetationClass();

  //Accessors
  string                   GetVegetationName() const;
  const veg_struct        *GetVegetationStruct() const;
  double                   GetParameter(const string param_name) const;//not currently used
  double                   GetVegetationProperty(string param_name) const;
  void                     SetVegetationProperty(const string &param_name, const double &value);

  //routines
  void AutoCalculateVegetationProps(const veg_struct    &Vtmp,
                                    const veg_struct    &Vdefault);

  static void                    PreInitialize();
  static int                     GetNumClasses();
  static const CVegetationClass *GetVegClass(int c);
  static       CVegetationClass *StringToVegClass(const string s);
  static void                    DestroyAllVegClasses();

  static void                    SetVegetationProperty(veg_struct &V, const string param_name, const double &value);
  static void                    SetVegTransportProperty( int          constit_ind,int          constit_ind2,
                                                          veg_struct  &V,string param_name, const double value);
  static double                  GetVegetationProperty(const veg_struct &V, string param_name, const bool strict=true);
  static double                  GetVegTransportProperty(int constit_ind, const veg_struct &V,string &param_name);

  static void                    InitializeVegetationProps(string name, veg_struct &V, bool is_template);
  //static void                    InitializeVegetationProps(double *params, bool is_template);

  static void                    SummarizeToScreen();

  //----------------------------------------------------------------
  //below routines defined in VegetationParams.cpp:
  static void   RecalculateCanopyParams(      veg_var_struct    &VV,         //transient canopy params
                                              const CHydroUnit       *pHRU,
                                              const CModelABC        *pModel,
                                              //const CVegetationClass  *pVeg,       //fixed vegetation params
                                              //const surface_struct    &G,          //ground surface
                                              //const force_struct      *F,          //forcing functions over tstep
                                              //const double             soil_def,   //soil moisture deficit [mm]
                                              //const double             snowdepth,  //snowdepth [mm]
                                              const time_struct       &tt,
                                              const optStruct         &Options);

  static void     RecalculateRootParams(        veg_var_struct  &VV,
                                                const CHydroUnit       *pHRU,
                                                const CModelABC        *pModel,
                                                //const CVegetationClass  *pVeg,          //fixed vegetation params
                                                //const surface_struct    &G,          //ground surface
                                                //const soil_struct        *S   [MAX_SOILLAYERS],
                                                //const double             thick[MAX_SOILLAYERS],//[m]
                                                //const int                nSoilLayers,
                                                const time_struct       &tt,
                                                const optStruct         &Options);
  static void      WriteVegetationProps(const veg_struct        &V,
                                        const veg_var_struct    &VV);
  static double     CalcClosedRoughness(const double            &height);
  static double CalcClosedZeroPlaneDisp(const double            &height,
                                        const double            &closed_rough);

};

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for land use classification
/// \details land use class determines ground surface properties
/// (e.g., impermeable fraction, fraction wetland,
/// depression storage, etc.)
//
class CLandUseClass
{
protected:/*----------------------------------------------------*/

  surface_struct            S;                    ///< corresponding surface properties

  static CLandUseClass    **pAllLUClasses;        ///< array of pointers to all LU classes that have been created
  static int                NumLUClasses;         ///< Number of land use classes (length of pAllLUClasses array)

public:/*-------------------------------------------------------*/
  //Constructors:
  CLandUseClass(const string name);
  ~CLandUseClass();

  //Accessors
  string                   GetLanduseName() const;
  const surface_struct    *GetSurfaceStruct() const;
  double                   GetSurfaceProperty(string param_name) const;
  void                     SetSurfaceProperty(const string &param_name, const double &value);

  //routines
  void AutoCalculateLandUseProps(const surface_struct &Stmp,
                                 const surface_struct &Sdefault);

  static int                     GetNumClasses();
  static const CLandUseClass    *GetLUClass(int c);
  static       CLandUseClass    *StringToLUClass(const string s);
  static void                    DestroyAllLUClasses();

  static void                    InitializeSurfaceProperties(string name, surface_struct &S, bool is_template);
  static void                    SetSurfaceProperty         (surface_struct &S, const string param_name, const double value);
  static double                  GetSurfaceProperty         (const surface_struct &S, string param_name, const bool strict=true);

  static void                    SummarizeToScreen();
};

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for terrain classification
/// \details terrain class determines topographic properties / relief
//
class CTerrainClass
{
protected:/*----------------------------------------------------*/

  string                    tag;                  ///< nickname for terr. class, e.g., "HILLY"
  terrain_struct            T;                    ///< corresponding terrain properties


  static CTerrainClass    **pAllTerrainClasses;   ///< array of pointers to all terrain classes that have been created
  static int                NumTerrainClasses;  ///< Number of terrain classes that have been created length of pAllTerrainClasses


public:/*-------------------------------------------------------*/
  //Constructors:
  CTerrainClass(const string name);
  ~CTerrainClass();

  //Accessors
  string                   GetTag() const;
  const terrain_struct    *GetTerrainStruct() const;
  double                   GetTerrainProperty(string param_name) const;
  void                     SetTerrainProperty(const string &param_name, const double &value);

  //routines
  void AutoCalculateTerrainProps(const terrain_struct &Ttmp, const terrain_struct &Tdefault);

  static int                     GetNumClasses();
  static const CTerrainClass    *GetTerrainClass(int c);
  static       CTerrainClass    *StringToTerrainClass(const string s);
  static void                    DestroyAllTerrainClasses();

  static void                    InitializeTerrainProperties(terrain_struct &T, bool is_template);
  static void                    SetTerrainProperty(terrain_struct &T, const string  param_name, const double value);
  static double                  GetTerrainProperty(const terrain_struct &T, string param_name);

  static void                    SummarizeToScreen();
};

#endif
