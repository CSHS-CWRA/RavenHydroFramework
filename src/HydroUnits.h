/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2023 the Raven Development Team
  ----------------------------------------------------------------*/
#ifndef HYDROUNITS_H
#define HYDROUNITS_H

#include "RavenInclude.h"
#include "ModelABC.h"
#include "Properties.h"
#include "SoilAndLandClasses.h"
#include "GlobalParams.h"
#include "SoilProfile.h"

///////////////////////////////////////////////////////////////////
/// \brief Abstract class representing single portion of watershed
/// \details Data abstraction for single contiguous portion of watershed
/// with homogeneous land use/ land type and homogeneous (or at least
/// statistically homogenous) properties
//
class CHydroUnit
{
private:/*------------------------------------------------------*/

  const CModelABC            *_pModel;  ///< Pointer to model

  int                             _ID;  ///< Unique HRU identifier
  int                       _global_k;  ///< Global model index as stored in CModel array (can/will be different from ID)
  double                        _Area;  ///< contributing drainage area for HydroUnit [km^2]
  int                    _SubbasinInd;  ///< global index p (not ID!) of subbasin this HRU is in
  location                  _Centroid;  ///< centroid of HRU
  HRU_type                   _HRUType;  ///< Standard, Lake, Rock, Glacier, etc...
  bool                      _Disabled;  ///< true if processes are not simulated for this HRU
  bool                    _res_linked;  ///> true if HRU is linked to Reservoir

  //Model State variables:
  double                  *_aStateVar;  ///< Array of *current value* of state variable i with size CModel::nStateVars [mm] for water storage, permafrost depth, snow depth, [MJ/m^2] for energy storage

  //Model Forcing functions:
  force_struct              _Forcings;  ///< *current values* of forcing functions for time step (precip, temp, etc.)

  //property structures (from soil, veg, aq, surface classes, class-specific)
  double                _AvgElevation;  ///< average elevation of HRU [masl]
  double                   _AvgAspect;  ///< average terrain aspect [rad counterclockwise from north]
  double                    _AvgSlope;  ///< average terrain slope [rad]

  double                      _LatRad;  ///< latitude of centroid [rad]
  double                       _LatEq;  ///< equivalent latitude for slope/aspect  [rad]
  double                   _SolarNoon;  ///< effective solar noon correction for slope [days]

  double                  _PrecipMult;  ///< HRU-specific precipitation correction factor
  int              _SpecifiedGaugeIdx;  ///< index of met gauge to be used regardless of interpolation scheme (or DOESNT_EXIST, by default)

  const soil_struct           *_pSoil[MAX_SOILLAYERS];     ///< array of pointers to structures with profile soil properties

  const veg_struct             *_pVeg;  ///< pointer to structure with vegetation properties
  const surface_struct     *_pSurface;  ///< pointer to structure with land use/land type properties
  const terrain_struct     *_pTerrain;  ///< pointer to structure with terrain properties

  //variable property structures (locally stored, HRU-specific)
  double                   aThickness[MAX_SOILLAYERS];     ///< soil layer thicknesses [m]
  //double                  aRootFrac[MAX_SOILLAYERS];     ///< root distribution (sum over all layers == 1)

  /// /todo allow for variable aquifer properties in each HRU (see above)
  veg_var_struct              _VegVar;  ///< Points to derived vegetation properties

public:/*-------------------------------------------------------*/
  //Constructors:

  CHydroUnit(const CModelABC        *pMod,
             const int               ID,
             const int               global_ind,
             const double            drainage_area,    //land surface area, km^2
             const int               basin_index,      //index of parent watershed
             const double            elevation,        //[ma.s.l.]
             const double            Latitude,         //decimal lat,long
             const double            Longitude,
             const double            slope,            //[rad]
             const double            aspect,           //[rad counterclockwise from north]
             const HRU_type          type,
             const CSoilProfile     *soil_profile,     //soil profile
             const CVegetationClass *veg_class,        //vegetation class
             const void             *PLACEHOLDER,      //old aquifer class
             const CTerrainClass    *terrain_class,    //terrain class
             const CLandUseClass    *lult_class);      //land use

  ~CHydroUnit();

  //Accessor functions (some inlined for speed)
  inline int             GetID           () const { return _ID;          }
  inline int             GetGlobalIndex  () const { return _global_k;    }
  inline bool            IsEnabled       () const { return !_Disabled;   }
  inline location        GetCentroid     () const { return _Centroid;    }
  inline double          GetArea         () const { return _Area;        }
  inline int             GetSubBasinIndex() const { return _SubbasinInd; }
  inline HRU_type        GetHRUType      () const { return _HRUType;     } //(standard, lake, rock, or glacier)
  inline bool            IsLake          () const { return (_HRUType==HRU_LAKE);}
  bool                   IsLinkedToReservoir() const;

  inline double          GetStateVarValue(const int i) const { return _aStateVar[i]; }
  inline double          GetConcentration(const int i) const { return _pModel->GetConcentration(_global_k,i); }
  double*                GetStateVarArray() const;

  double                 GetElevation    () const;//[masl]
  double                 GetSlope        () const;//[rad]
  double                 GetAspect       () const;//[rad]
  double                 GetLatRad       () const;//[rad]
  double                 GetLatEq        () const;//[rad]
  double                 GetSolarNoon    () const;//[days]

  double                 GetPrecipMultiplier() const;
  int                    GetSpecifiedGaugeIndex() const;

  soil_struct     const *GetSoilProps       (const int m) const;
  double                 GetSoilThickness   (const int m) const;//[mm]
  double                 GetSoilCapacity    (const int m) const;//[mm]
  double                 GetSoilTensionStorageCapacity(const int m) const;//[mm]



  veg_struct      const *GetVegetationProps () const;
  veg_var_struct  const *GetVegVarProps     () const;
  surface_struct  const *GetSurfaceProps    () const;
  terrain_struct  const *GetTerrainProps    () const;

  force_struct    const *GetForcingFunctions() const;
  double                 GetForcing         (const forcing_type &ftype) const;

  double                 GetCumulFlux       (const int i, const bool to) const;
  double                 GetCumulFluxBet    (const int iFrom, const int iTo) const;

  double                 GetStateVarMax     (const int     i,
                                             const double *curr_state_var,
                                             const optStruct &Options,
                                             const bool ignorevar=false) const;

  //specialized state variable accessors 
  double                 GetSnowAlbedo        () const;
  double                 GetSnowTemperature   () const;
  double                 GetSnowSWE           () const;
  double                 GetSnowCover         () const;
  double                 GetSurfaceTemperature() const;
  double                 GetTotalAlbedo       () const;
  double                 GetSnowDepth         () const;

  //Manipulator functions (used in parser)
  void                   Disable              ();
  void                   Enable               ();
  void                   LinkToReservoir      (const long SBID);

  //Manipulator functions (used in initialization)
  void          Initialize              (const int UTM_zone);

  //Manipulator functions (used in solution method)
  void          SetStateVarValue        (const int           i,
                                         const double       &new_value);
  void          UpdateForcingFunctions  (const force_struct &Fnew);
  void          CopyDailyForcings       (force_struct &F);
  void          SetPrecipMultiplier     (const double factor);
  void          SetSpecifiedGaugeIndex  (const int g);
  void          ChangeLandUse           (const CLandUseClass    *lult_class);
  void          ChangeVegetation        (const CVegetationClass *veg_class);
  void          ChangeHRUType           (const HRU_type typ);
  void          AdjustHRUForcing        (const forcing_type Ftyp,const double& epsilon, const adjustment adj); 

  //will be removed with landscape elements:
  void          RecalculateDerivedParams(const optStruct    &Options,
                                         const time_struct  &tt);
};

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction class for collection of HRUs
//
class CHRUGroup
{
private:/*------------------------------------------------------*/

    string            _name; ///< HRU group name
    int              _nHRUs; ///< Number of HRUs present in group
    CHydroUnit     **_pHRUs; ///< Array of pointers to member HRUs (size = nHRUs)
    int          _global_kk; ///< index of group in master HRU Group array (in CModel)
    bool          _disabled; ///< true if all HRUs in group are disabled

public:/*-------------------------------------------------------*/
  //Constructors:
  CHRUGroup(string tag, int global_ind);
  ~CHRUGroup();

  void AddHRU(CHydroUnit *pHRU);
  void DisableGroup();
  void Initialize();
  void EmptyGroup();

  string            GetName            () const;
  int               GetNumHRUs         () const;
  int               GetGlobalIndex     () const;
  CHydroUnit       *GetHRU             (const int k_local) const;
  double            GetAvgStateVar     (const int i) const;
  double            GetAvgConcentration(const int i) const;
  double            GetAvgForcing      (const forcing_type &ftype) const;
  double            GetAvgCumulFlux    (const int i, const bool to) const;
  double            GetAvgCumulFluxBet (const int iFrom, const int iTo) const;
  bool              IsInGroup          (const int k) const;
};
#endif
