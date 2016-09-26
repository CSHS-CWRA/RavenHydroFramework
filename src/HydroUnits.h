/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright © 2008-2014 the Raven Development Team
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
    int                    _SubbasinInd;  ///< global index (not ID!) of subbasin this HRU is in
    location                  _Centroid;  ///< centroid of HRU
    HRU_type                   _HRUType;  ///< Standard, Lake, Glacier, etc...
    bool                      _Disabled;  ///< true if processes are not simulated for this HRU

    //Model State variables:
    double                  *_aStateVar;  ///< Array of *current value* of state variable i with size CModel::nStateVars [mm] for water storage, permafrost depth, snow depth, [MJ/m^2] for energy storage

    //Model Forcing functions:
    force_struct              _Forcings;  ///< *current values* of forcing functions for time step (precip, temp, etc.)   
    
    //property structures (from soil, veg, aq, surface classes, class-specific)
    double                _AvgElevation;  ///< average elevation of HRU [masl]
    double                   _AvgAspect;  ///< average terrain aspect [rad]
    double                    _AvgSlope;  ///< average terrain slope [rad]

	  double					            _LatRad;  ///< latitude of centroid [rad]
    double						           _LatEq;  ///< equivalent latitude for slope/aspect  [rad]
    double					       	 _SolarNoon;  ///< effective solar noon correction for slope [days]

    double					       	_PrecipMult; ///< HRU-specific precipitation correction factor

    const soil_struct           *_pSoil[MAX_SOILLAYERS];     //< array of pointers to structures with profile soil properties
    const soil_struct        *pAquifers[MAX_AQUIFER_LAYERS]; //< array of pointers to structures with aquifer layer properties
    //const soil_struct      *pAquitard[MAX_AQUIFER_LAYERS]; //< /todo Add aquitard information
    const veg_struct             *_pVeg;  //< pointer to structure with vegetation properties
    const CVegetationClass * pVegetation;  //< pointer to corresponding vegetation class (not currently used)
    const surface_struct     *_pSurface;  //< pointer to structure with land use/land type properties
    const terrain_struct     *_pTerrain;  //< pointer to structure with terrain properties

    //variable property structures (locally stored, HRU-specific)
    double                   aThickness[MAX_SOILLAYERS];     ///< soil layer thicknesses [m]
    double                 aAqThickness[MAX_AQUIFER_LAYERS]; ///< aquifer layer thicknesses [m]
    //double             aAquitardThick[MAX_AQUIFER_LAYERS]; ///< aquitard layer thicknesses [m]
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
               const double            aspect,           //[rad from north]
               const HRU_type          type,              
               const CSoilProfile     *soil_profile,     //soil profile
               const CVegetationClass *veg_class,        //vegetation class 
               const void             *PLACEHOLDER,      //old aquifer class
               //const CAquiferStack    *aquifer_system,   //aquifer system profile
               const CTerrainClass    *terrain_class,    //terrain class
               const CLandUseClass    *lult_class);      //land use 

    ~CHydroUnit();

    //Accessor functions
    int                    GetID           () const;
    int                    GetGlobalIndex  () const;
    bool                   IsEnabled       () const;
    location               GetCentroid     () const;
    double                 GetArea         () const;
    int                    GetSubBasinIndex() const;
    HRU_type               GetType         () const;
    bool                   IsLake          () const;//TMP?

    double                 GetStateVarValue(const int i) const;
	  double*                GetStateVarArray() const;

    double                 GetElevation    () const;//[masl]
    double                 GetSlope        () const;//[rad]
    double                 GetAspect       () const;//[rad]
    double				         GetLatRad       () const;//[rad]
    double                 GetLatEq        () const;//[rad]
    double                 GetSolarNoon    () const;//[days]

	  double				         GetPrecipMultiplier () const;

    soil_struct     const *GetSoilProps    (const int m) const;
    double                 GetSoilThickness(const int m) const;//[mm]
    double                 GetSoilCapacity (const int m) const;//[mm]
    double                 GetSoilTensionStorageCapacity(const int m) const;//[mm]

    soil_struct     const *GetAquiferProps    (const int m) const;
    double                 GetAquiferThickness(const int m) const;//[mm]
    double                 GetAquiferCapacity (const int m) const;//[mm]

    veg_struct      const *GetVegetationProps() const;
    veg_var_struct  const *GetVegVarProps    () const;
    surface_struct  const *GetSurfaceProps   () const;
    terrain_struct  const *GetTerrainProps   () const;
    
    force_struct    const *GetForcingFunctions() const;
    double                 GetForcing         (const string &forcing_string) const;

    double                 GetCumulFlux       (const int i, const bool to) const;

    double                 GetStateVarMax     (const int     i,
                                               const double *curr_state_var,
                                               const optStruct &Options) const;

    //specialized state variable accessors // \todo [funct] - should likely hand in representative state_var array
    double                 GetSnowAlbedo        () const;
    double                 GetSnowTemperature   () const;
    double                 GetSnowSWE           () const;
    double                 GetSnowCover         () const;
    double                 GetSurfaceTemperature() const;
    double                 GetTotalAlbedo       () const;

    //Manipulator functions (used in parser)
    void          Disable                 ();

    //Manipulator functions (used in initialization)
    void          Initialize              (const int UTM_zone);

    //Manipulator functions (used in solution method)
    void          SetStateVarValue        (const int           i, 
                                           const double       &new_value);
    void          UpdateForcingFunctions  (const force_struct &Fnew);
	  void          CopyDailyForcings       (force_struct &F);
	  void				  SetPrecipMultiplier     (const double factor);
    void          ChangeLandUse           (const CLandUseClass    *lult_class);
	  void          ChangeVegetation        (const CVegetationClass *veg_class);

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

    string name;        ///< HRU group name
    int nHRUs;          ///< Number of HRUs present in group
    CHydroUnit **pHRUs; ///< Array of pointers to member HRUs (size = nHRUs)
    int global_kk;      ///< index of group in master HRU Group array (in CModel)

    bool *aAggregateSV; ///< array (size MAX_STATE_VARS) - true if state variable is averaged over HRU group each timestep

  public:/*-------------------------------------------------------*/
    //Constructors:
    CHRUGroup(string tag, int global_ind);
    ~CHRUGroup();

    void AddHRU(CHydroUnit *pHRU);
    //void DisableAllInGroup();

    string            GetName            () const;
    int               GetNumHRUs         () const;
    int               GetGlobalIndex     () const;
    CHydroUnit       *GetHRU             (const int k_local) const;
    double            GetAvgStateVar     (const int i) const;
    double            GetAvgForcing      (const string &forcing_string) const;
    double            GetAvgCumulFlux    (const int i, const bool to) const;
    bool              IsInGroup          (const int k) const;

    bool              IsAggregatorGroup  (const int i) const; 

    void              SetAsAggregator    (const int i);
};
#endif