/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2024 the Raven Development Team
  ----------------------------------------------------------------*/
#include "HydroUnits.h"
#include "Forcings.h"
#include "Radiation.h"

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the HydroUnit (HRU) constructor
///
/// \param *pMod [in] pointer to surface water model
/// \param identifier [in] HRU identifier/ID
/// \param drainage_area [in] Drainage area [km^2]
/// \param basin_ind [in] Index of Subbasin to which HRU belongs
/// \param elevation [in] HRU elevation [masl]
/// \param latit [in] HRU centroid latitutde [degrees]
/// \param longit [in] HRU centroid longitude [degrees]
/// \param slope [in] HRU average/representative slope [rad]
/// \param aspect [in] HRU average/representative aspect [rad]
/// \param typ [in] HRU type (e.g., standard/glacier/lake...)
/// \param  *soil_profile [in] Soil profile type
/// \param *aquifer_system [in] aquifer profile
/// \param *veg_class [in] Dominant vegetation class
/// \param *terrain_class [in] Dominant terrain class
/// \param *lult_class [in] Dominant land use class
//
CHydroUnit::CHydroUnit(const CModelABC        *pMod,
                       const long long int     identifier,
                       const int               global_ind,
                       const double            drainage_area,
                       const int               basin_ind,
                       const double            elevation,
                       const double            latit,
                       const double            longit,
                       const double            slope,//[rad]
                       const double            aspect,//[rad counterclock from N]
                       const HRU_type          typ,
                       const CSoilProfile     *soil_profile,
                       const CVegetationClass *veg_class,
                       const void             *PLACEHOLDER,
                       const CTerrainClass    *terrain_class,
                       const CLandUseClass    *lult_class)
{
  string error;
  if (drainage_area<=0.0){
    error="CHydroUnit constructor:: HRU "+to_string(identifier)+" has a negative or zero area";
    ExitGracefully(error.c_str(),BAD_DATA_WARN);
  }
  if (basin_ind<0){
    error="CHydroUnit constructor:: HRU "+to_string(identifier)+" has a bad subbasin index specified";
    ExitGracefully(error.c_str(),BAD_DATA_WARN);
  }
  if (soil_profile==NULL){
    error="CHydroUnit constructor:: HRU "+to_string(identifier)+": invalid soil profile specified";
    ExitGracefully(error.c_str(),BAD_DATA);
  }
  if (lult_class==NULL){
    error="CHydroUnit constructor:: HRU "+to_string(identifier)+": invalid land use class specified";
    ExitGracefully(error.c_str(),BAD_DATA);
  }
  if (veg_class==NULL){
    error="CHydroUnit constructor:: HRU "+to_string(identifier)+": invalid vegetation class specified";
    ExitGracefully(error.c_str(),BAD_DATA);
  }
  if (pMod==NULL){
    error="CHydroUnit constructor:: HRU "+to_string(identifier)+": invalid model pointer specified";
    ExitGracefully(error.c_str(),RUNTIME_ERR);
  }

  int i;
  _ID                   =identifier;
  _global_k             =global_ind;
  _pModel               =pMod;
  _Area                 =drainage_area;
  _SubbasinInd          =basin_ind;
  _Centroid.latitude    =latit;
  _Centroid.longitude   =longit;
  _Centroid.UTM_x       =0.0;//corrected in Initialize Routine
  _Centroid.UTM_y       =0.0;
  _Disabled             =false;
  _res_linked           =false;
  _HRUType              =typ;

  _aStateVar=new double [_pModel->GetNumStateVars()];
  for (i=0;i<_pModel->GetNumStateVars();i++){
    _aStateVar[i]=0.0;
  }

  ZeroOutForcings(_Forcings);

  _AvgElevation =elevation;
  _AvgAspect    =aspect; //counterclockwise from north
  _AvgSlope     =slope;

  _LatRad       = latit*DEGREES_TO_RADIANS;
  _LatEq        = CRadiation::CalculateEquivLatitude(_LatRad,slope,aspect);

  //calculate corrected solar noon for given slope/aspect
  _SolarNoon =CRadiation::CalculateSolarNoon(_LatRad,slope,aspect);

  _pVeg      =veg_class    ->GetVegetationStruct();
  _pSurface  =lult_class   ->GetSurfaceStruct();
  if (terrain_class!=NULL){ _pTerrain  =terrain_class->GetTerrainStruct();}
  else                    { _pTerrain=NULL;}
  _PrecipMult = 1.0;
  _SpecifiedGaugeIdx=DOESNT_EXIST;

  soil_profile->AllocateSoilLayers(_pModel->GetNumSoilLayers(),_pSoil,aThickness);
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the HRU destructor
///
CHydroUnit::~CHydroUnit()
{
  if (DESTRUCTOR_DEBUG){cout<<"    DELETING HYDROUNIT"<<endl;}
  delete [] _aStateVar; _aStateVar=NULL;
}
/*****************************************************************
   Accessors
*****************************************************************/
//////////////////////////////////////////////////////////////////
/// \brief Returns a pointer to the state variable array
///
/// \return Double*pointer to the state variable array
//
double*    CHydroUnit::GetStateVarArray                 () const
{
  return _aStateVar;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns average elevation of HRU above sea level [m]
///
/// \return Double value of average HRU elevation [masl]
//
double    CHydroUnit::GetElevation                      () const {return _AvgElevation;}

//////////////////////////////////////////////////////////////////
/// \brief Returns average terrain slope [rad]
///
/// \return Double average terrain slope [rad]
//
double    CHydroUnit::GetSlope                          () const {return _AvgSlope;}//[rad]

//////////////////////////////////////////////////////////////////
/// \brief Returns average terrain aspect [rad]
///
/// \return Double average terrain aspect [rad]
//
double    CHydroUnit::GetAspect                         () const {return _AvgAspect;}//[rad]

//////////////////////////////////////////////////////////////////
/// \brief Returns centroid latitude [rad]
///
/// \return Double centroid latitude                    (in radians)
//
double    CHydroUnit::GetLatRad                         () const {return _LatRad;}//[rad]

//////////////////////////////////////////////////////////////////
/// \brief Returns the equivalent latitude for slope [rad]
///
/// \return Double equivalent latitude for slope [rad]
//
double    CHydroUnit::GetLatEq                          () const {return _LatEq;}//[rad]

//////////////////////////////////////////////////////////////////
/// \brief Returns the effective solar noon correction for slope [days]
///
/// \return Double solar noon correction
//
double    CHydroUnit::GetSolarNoon                      () const {return _SolarNoon;}//[days]

//////////////////////////////////////////////////////////////////
/// \brief Returns the HRU-specific precipitation correction factor
///
/// \return double precipitation correction factor
//
double    CHydroUnit::GetPrecipMultiplier               () const {return _PrecipMult;}

//////////////////////////////////////////////////////////////////
/// \brief Returns the HRU-specific gauge index if interpolation is to be overriden (or DOESNT_EXIST, if base interpolation scheme is to be used)
///
/// \return int gauge ID
//
int       CHydroUnit::GetSpecifiedGaugeIndex            () const {return _SpecifiedGaugeIdx;}

//////////////////////////////////////////////////////////////////
/// \brief Checks if HRU is linked to reservoir
///
/// \return Boolean showing if HRU is linked to reservoir or not
//
bool    CHydroUnit::IsLinkedToReservoir                 () const {return _res_linked;}

//////////////////////////////////////////////////////////////////
/// \brief Returns canopy properties of vegetation in HRU
///
/// \return Canopy structure containing canopy properties of HRU
//
veg_struct   const *CHydroUnit::GetVegetationProps      () const {return _pVeg;}

//////////////////////////////////////////////////////////////////
/// \brief Returns canopy variable properties
///
/// \return Pointer to structure containing canopy variable properties
//
veg_var_struct  const *CHydroUnit::GetVegVarProps       () const {return &_VegVar;}

//////////////////////////////////////////////////////////////////
/// \brief Returns surface properties of HRU
///
/// \return Pointer to structure containing surface properties of HRU
//
surface_struct  const *CHydroUnit::GetSurfaceProps      () const {return _pSurface;}

//////////////////////////////////////////////////////////////////
/// \brief Returns terrain properties of HRU
///
/// \return Pointer to structure containing terrain properties of HRU
//
terrain_struct  const *CHydroUnit::GetTerrainProps      () const {return _pTerrain;}

//////////////////////////////////////////////////////////////////
/// \brief Returns soil properties of a layer of soil in HRU
///
/// \param m [in] Integer index specifying soil layer
/// \return Pointer to structure containing soil properties of HRU in layer m
//
soil_struct     const *CHydroUnit::GetSoilProps       (const int m) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((m<0) || (m>=_pModel->GetNumSoilLayers()),
                   "CHydroUnit GetSoilProps::improper index",BAD_DATA);
#endif
  return _pSoil[m];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns soil thickness of a layer of soil in HRU [mm]
///
/// \param m [in] Integer index specifying soil layer
/// \return soil thickness of a layer of soil in HRU [mm]
//
double  CHydroUnit::GetSoilThickness   (const int m) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((m<0) || (m>=_pModel->GetNumSoilLayers()),
                   "CHydroUnit GetSoilThickness::improper index",BAD_DATA);
#endif
  return aThickness[m]*MM_PER_METER;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns soil capacity of a layer of soil in HRU
///
/// \param m [in] Integer index specifying soil layer
/// \return maximum water storage capacity of soil layer in [mm]
//
double CHydroUnit::GetSoilCapacity   (const int m) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((m<0) || (m>=_pModel->GetNumSoilLayers()),
                   "CHydroUnit GetSoilCapacity::improper index",BAD_DATA);
#endif
  return aThickness[m]*MM_PER_METER*_pSoil[m]->cap_ratio;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns soil tension storage capacity of a layer of soil in HRU  (MC at field cap-MC at wilting point ) [mm]
///
/// \param m [in] Integer index specifying soil layer
/// \return Double indicating soil tension storage capacity of layer m [mm]
//
double CHydroUnit::GetSoilTensionStorageCapacity(const int m) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((m<0) || (m>=_pModel->GetNumSoilLayers()),
                   "CHydroUnit GetSoilTensionStorageCapacity::improper index",BAD_DATA);
#endif
  double storcap=aThickness[m]*MM_PER_METER*_pSoil[m]->cap_ratio*(_pSoil[m]->field_capacity-_pSoil[m]->sat_wilt);
  if (storcap<-REAL_SMALL){
    ExitGracefully("Wilting point saturation (SAT_WILT) is greater than field capacity (FIELD_CAPACITY) for one or more soil classes, leading to negative tension storage",BAD_DATA);
  }
  return storcap;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns pointer to HRU forcing functions
///
/// \return pointer to current forcing functions
//
force_struct     const *CHydroUnit::GetForcingFunctions () const
{
  return &_Forcings;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns forcing function value for this HRU according to string input
///
/// \param &ftype [in] enum of forcing function (e.g. F_PRECIP)
/// \return Double specified forcing function value
//
double CHydroUnit::GetForcing(const forcing_type &ftype) const
{
  return GetForcingFromType(ftype,_Forcings);
}
//////////////////////////////////////////////////////////////////
/// \brief Returns area-weighted average of specified cumulative flux over HRU group
///
/// \param i [in] index of storage compartment
/// \param to [in] true if evaluating cumulative flux to storage compartment, false for 'from'
/// \return cumulative flux to date to storage compartment i
//
double CHydroUnit::GetCumulFlux(const int i, const bool to) const
{
  return _pModel->GetCumulativeFlux(_global_k, i,to);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns cumulative flux between two compartments in HRU
///
/// \param iFrom [in] index of 'from' storage compartment
/// \param iTo [in] index of 'to' storage compartment
/// \return cumulative flux to date between two compartments within HRU
//
double CHydroUnit::GetCumulFluxBet(const int iFrom, const int iTo) const
{
  return _pModel->GetCumulFluxBetween(_global_k, iFrom,iTo);
}


/*****************************************************************
   Manipulators
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Sets state variable value
/// \remarks should only be called from within MassEnergyBalance routine at end of timestep!!!
///
/// \param i [in] Integer identifier of state variable
/// \param val [in] Double value of state variable to be set
//
void    CHydroUnit::SetStateVarValue      (const int i, const double &val)
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((i<0) || (i>=_pModel->GetNumStateVars()),
                   "CHydroUnit SetStateVarValue::improper index",BAD_DATA);
#endif
  _aStateVar[i]=val;
}

//////////////////////////////////////////////////////////////////
/// \brief Disables HRU
//
void CHydroUnit::Disable(){_Disabled=true;}

//////////////////////////////////////////////////////////////////
/// \brief Enables HRU
//
void CHydroUnit::Enable(){_Disabled=false;}

//////////////////////////////////////////////////////////////////
/// \brief links HRU to reservoir
//
void CHydroUnit::LinkToReservoir(const long long SBID){_res_linked=true;}


//////////////////////////////////////////////////////////////////
/// \brief Initializes HRU - converts coordinates of HRU centroid into UTM
/// \remark Called prior to simulation or any interpolation;
///
/// \param UTM_zone [out] Integer UTM zone
//
void    CHydroUnit::Initialize      (const int UTM_zone)
{
  LatLonToUTMXY(_Centroid.latitude,_Centroid.longitude,
                UTM_zone,
                _Centroid.UTM_x,   _Centroid.UTM_y);
}

//////////////////////////////////////////////////////////////////
/// \brief Recalculates derived parameters
/// \details Recalculates parameters that rely on time of year or forcing functons,
/// which have been updated since most recent calculation. Called from MassEnergyBalance() method
/// \note May wish to only recalculate daily or otherwise intermittently
///
/// \param &Options [in] Global model options information
/// \param &tt [in] current time
//
void CHydroUnit::RecalculateDerivedParams(const optStruct &Options,
                                          const time_struct &tt)
{
  CVegetationClass::RecalculateCanopyParams(_VegVar,this,_pModel,tt,Options);
  CVegetationClass::RecalculateRootParams  (_VegVar,this,_pModel,tt,Options);
}

//////////////////////////////////////////////////////////////////
/// \brief Updates forcing function
/// \note Called by model before each time step (Fnew generated by UpdateForcingFunctions routine)
///
/// \param &Fnew [in] New forcing functions to be updated into system
//
void CHydroUnit::UpdateForcingFunctions(const force_struct &Fnew)
{
  _Forcings=Fnew;
}

//////////////////////////////////////////////////////////////////
/// \brief Copies daily forcing values
/// \note Called by UpdateForcingFunctions for subdaily time steps to prevent repeated calculations
///
/// \param &F [out] New forcing functions being calculated
//
void CHydroUnit::CopyDailyForcings(force_struct &F)
{
  CopyDailyForcingItems(_Forcings,F);
}

//////////////////////////////////////////////////////////////////
/// \brief Sets the HRU-specific precipitation correction factor
///
//
void      CHydroUnit::SetPrecipMultiplier     (const double factor)
{
  _PrecipMult = factor;
}
//////////////////////////////////////////////////////////////////
/// \brief Sets the HRU-specific  specified gauge index to override interpolation weights
///
//
void     CHydroUnit::SetSpecifiedGaugeIndex(const int g)
{
  _SpecifiedGaugeIdx=g;
 }
//////////////////////////////////////////////////////////////////
/// \brief Changes the land use class mid-simulation
//
void CHydroUnit::ChangeLandUse(const CLandUseClass    *lult_class)
{
  _pSurface=lult_class->GetSurfaceStruct();
}
//////////////////////////////////////////////////////////////////
/// \brief Changes the vegetation class mid-simulation
//
void CHydroUnit::ChangeVegetation(const CVegetationClass *veg_class)
{
  _pVeg = veg_class->GetVegetationStruct();
}
//////////////////////////////////////////////////////////////////
/// \brief Changes the HRU Type mid-simulation
//
void CHydroUnit:: ChangeHRUType(const HRU_type typ)
{
  _HRUType=typ;
}
//////////////////////////////////////////////////////////////////
/// \brief Adjust HRU Forcing values mid-simulation
//
void CHydroUnit::AdjustHRUForcing(const forcing_type Ftyp,force_struct &F, const double& epsilon, const adjustment adj)
{
  if(Ftyp==F_PRECIP) {
    if      (adj==ADJ_MULTIPLICATIVE){F.precip*=epsilon;}
    else if (adj==ADJ_ADDITIVE      ){F.precip+=epsilon;}
    else if (adj==ADJ_REPLACE       ){F.precip =epsilon;}
    upperswap(F.precip,0.0);
  }
  else if(Ftyp==F_RAINFALL) {
    double sf=F.snow_frac;
    double Po=F.precip;
    if      (adj==ADJ_MULTIPLICATIVE){F.precip+=Po*(1.0-sf)*(epsilon-1.0);}
    else if (adj==ADJ_ADDITIVE      ){F.precip+=epsilon; }
    else if (adj==ADJ_REPLACE       ){F.precip+=epsilon-(1-sf)*Po;}
    upperswap(F.precip,0.0);
    if (F.precip==0){F.snow_frac=0;}
    else            {F.snow_frac=F.snow_frac*Po/F.precip;}
  }
  else if(Ftyp==F_SNOWFALL) {
    double sf=F.snow_frac;
    double Po=F.precip;
    if      (adj==ADJ_MULTIPLICATIVE){F.precip+=Po*(sf)*(epsilon-1.0);}
    else if (adj==ADJ_ADDITIVE      ){F.precip+=epsilon; }
    else if (adj==ADJ_REPLACE       ){F.precip+=epsilon-(sf)*Po;}
    upperswap(F.precip,0.0);
    if (F.precip==0){F.snow_frac=0;}
    else            {F.snow_frac=1.0-(1.0-F.snow_frac)*Po/F.precip;}
  }
  else if (Ftyp == F_TEMP_AVE) {
    //daily mean/min/max handled via AdjustDailyHRUForcings()
    if      (adj==ADJ_MULTIPLICATIVE){F.temp_ave*=epsilon; }
    else if (adj==ADJ_ADDITIVE      ){F.temp_ave+=epsilon; }
    else if (adj==ADJ_REPLACE       ){F.temp_ave =epsilon;}
  }
  else {
    ExitGracefully("CHydroUnit::AdjustHRUForcing: Only PRECIP, SNOWFALL, RAINFALL, and TEMP_AVE are supported for forcing perturbation. ",BAD_DATA);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Adjust HRU daily mean/min/max Forcing values mid-simulation (only used for temperature)
//
void CHydroUnit::AdjustDailyHRUForcings(const forcing_type Ftyp, force_struct &F,  const double* epsilon, const adjustment adj, const int nStepsPerDay)
{
  if (Ftyp == F_TEMP_AVE)
  {
    if      (adj==ADJ_MULTIPLICATIVE)
    {
      double adjust=0;
      for (int n=0;n<nStepsPerDay;n++){
        adjust+=1.0*epsilon[n]; //mean multiplicative perturbation over the day -exact for single time step
      }
      F.temp_daily_min*=adjust; //approximate - uses mean perturbation to adjust extremes
      F.temp_daily_max*=adjust; //approximate - uses mean perturbation to adjust extremes
      F.temp_daily_ave*=adjust; //approximate - uses mean perturbation to adjust extremes
    }
    else if (adj==ADJ_ADDITIVE      )
    {
      double adjust=0;
      for (int n=0;n<nStepsPerDay;n++){
        adjust+=epsilon[n];
      }
      F.temp_daily_ave+=adjust; //exact
      F.temp_daily_min+=adjust; //approximate - uses mean perturbation to adjust extremes
      F.temp_daily_max+=adjust; //approximate - uses mean perturbation to adjust extremes
    }
  }

  //doesn't make sense for other forcings or REPLACEMENT

}


//////////////////////////////////////////////////////////////////
/// \brief Set forcing values mid-simulation - used by BMI
//
void CHydroUnit::SetHRUForcing(const forcing_type Ftyp, const double &val)
{
  SetForcingFromType(Ftyp,_Forcings,val);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns maximum allowable value of state variable in model
///
/// \param i [in] Integer index of state variable type
/// \param *curr_state_var [in] Current array of state variables in HRU
/// \param &Options [in] Global model options information
/// \param ignorevar [in] true if variable values (such as those linked to vegetation)
///        should be disregarded (for applying initial conditions where these are not yet initialized)
/// \return Double representing maximum value of state variable
//
double        CHydroUnit::GetStateVarMax(const int      i,
                                         const double  *curr_state_var,
                                         const optStruct &Options,
                                         const bool ignorevar) const
{
  double max_var=ALMOST_INF;
  switch (_pModel->GetStateVarType(i))
  {
    case(SOIL):         {
      max_var=GetSoilCapacity(_pModel->GetStateVarLayer(i));  break;
    }
    case(CANOPY):       {
      if(!ignorevar) { max_var=_VegVar.capacity; }            break;
    }
    case(CANOPY_SNOW):  {
      if(!ignorevar) {max_var=_VegVar.snow_capacity;}         break;
    }
    case(GROUNDWATER):  {
      if(Options.modeltype != MODELTYPE_SURFACE){}            break; //GetAquiferHeadCapacity(0,HRUid);} //GWMIGRATE - HRUid k not accesible from here yet
    }
    case(DEPRESSION):   {
      if (_pSurface->dep_max>=0){max_var=_pSurface->dep_max;} break;
    }
    case(SNOW_COVER):   {
      max_var=1.0;                                            break;
    }
    case(SNOW_LIQ):
    {
      if(!ignorevar) {
        int iSNO  =_pModel->GetStateVarIndex(SNOW);
        double SWE = curr_state_var[iSNO];
        max_var = CalculateSnowLiquidCapacity(SWE, GetSnowDepth(), _pModel);
      }
      break;
    }
    default:{
      max_var=ALMOST_INF;break;
    }
  }/* end switch*/
  return max_var;
}
//////////////////////////////////////////////////////////////////
/// \brief returns snow albedo
/// \note uses default snow albedo if not tracked as state variable
/// \note lagged - information specific to start of time step only
///
/// \return current snow albedo in HRU [dimensionless]
//
double CHydroUnit::GetSnowAlbedo() const
{
  int    iSnAlb=_pModel->GetStateVarIndex(SNOW_ALBEDO);
  if     (iSnAlb==DOESNT_EXIST){return DEFAULT_SNOW_ALBEDO;}
  else                         {return GetStateVarValue(iSnAlb);}
}
/////////////////////////////////////////////////////////////////////
/// \brief Calculates snow depth
/// \return Snow depth [mm]
//
double CHydroUnit::GetSnowDepth() const
{
  int iSnow   =_pModel->GetStateVarIndex(SNOW);
  if (iSnow==DOESNT_EXIST){return 0.0; }

  int iSnowDep=_pModel->GetStateVarIndex(SNOW_DEPTH);
  if (iSnowDep!=DOESNT_EXIST){return GetStateVarValue(iSnowDep); }

  double SWE= GetStateVarValue(iSnow);
  return (SWE/TYPICAL_SNOW_DENS)*DENSITY_WATER; //[mm]
}

//////////////////////////////////////////////////////////////////
/// \brief returns snow surface temperature
/// \note uses default snow temperature if not tracked as state variable
/// \note lagged - information specific to start of time step only
///
/// \return current snow surface temperature in HRU [deg C]
//
double CHydroUnit::GetSnowTemperature() const
{
  int iSnTemp = _pModel->GetStateVarIndex(SNOW_TEMP);
  if (iSnTemp == DOESNT_EXIST){
    double sntmp = _pModel->GetGlobalParams()->GetParams()->snow_temperature;
    if (sntmp == NOT_NEEDED_AUTO){
      return 0.0;
    }
    return sntmp;
  }
  else {
    return this->GetStateVarValue(iSnTemp);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief returns snow storage as SWE, averaged over HRU
/// \note lagged - information specific to start of time step only
///
/// \return current snow storage in HRU, in mm
//
double  CHydroUnit::GetSnowSWE      () const
{
  int    iSnow=_pModel->GetStateVarIndex(SNOW);
  if (iSnow==DOESNT_EXIST){return 0.0;}
  else                    {return this->GetStateVarValue(iSnow);}
}
//////////////////////////////////////////////////////////////////
/// \brief returns snow cover fraction of SWE
/// \note lagged - information specific to start of time step only
///
/// \return current snow cover fraction, [0..1]
//
double  CHydroUnit::GetSnowCover () const
{
  const optStruct *Options=_pModel->GetOptStruct();

  int    iSnFrac=_pModel->GetStateVarIndex(SNOW_COVER);
  int    iSnow  =_pModel->GetStateVarIndex(SNOW);

  if (iSnow==DOESNT_EXIST){return 0.0;}
  double SWE=this->GetStateVarValue(iSnow);

  if(iSnFrac!=DOESNT_EXIST) { //snow cover explicitly simulated
    return this->GetStateVarValue(iSnFrac);
  }
  else  //snowcover depletion curves
  {
    if(Options->snow_depletion==SNOWCOV_NONE)
    {
      if(GetSnowSWE()<NEGLIGBLE_SNOW) { return 0.0; }
      else                            { return 1.0; }
    }
    else if(Options->snow_depletion==SNOWCOV_LINEAR)
    {
      //double snowthresh=GetSurfaceProps()->SDC_threshold;
      double snowthresh=200;
      return min(max(SWE,0.0)/snowthresh,1.0);
    }
  }
  return 0.0;
}
//////////////////////////////////////////////////////////////////
/// \brief returns surface temperature
/// \note Assumes snow is at snow surface temperature and ground is average temp. Surface temp weighted by snow cover.
/// \note lagged - information specific to start of time step only
///
/// \return current surface temperature in HRU [deg C]
//
double CHydroUnit::GetSurfaceTemperature() const
{
  double snow_cover  = GetSnowCover();
  double snow_temp   = GetSnowTemperature();
  double ground_temp = _Forcings.temp_ave; //temp - more proper methods

  if (_HRUType == HRU_GLACIER) {
    ground_temp = FREEZING_TEMP;
  }

  return (1 - snow_cover)*ground_temp + snow_cover*snow_temp;
}
//////////////////////////////////////////////////////////////////
/// \brief returns total land surface albedo
/// \note lagged - information specific to start of time step only
///
/// \return current total land surface albedo in HRU [dimensionless]
//
double  CHydroUnit::GetTotalAlbedo(const bool subcanopy) const
{
  double veg_albedo,land_albedo(0.0);

  double snow_albedo=GetSnowAlbedo();
  double snow_cover =GetSnowCover();
  double svf=_VegVar.skyview_fact;
  double Fc=_pSurface->forest_coverage;

  // if (Options.albedo_type==ALBEDO_DEFAULT)
  {
    if (_HRUType==HRU_STANDARD){
      int iTopSoil=_pModel->GetStateVarIndex(SOIL,0);
      double soil_sat=max(min(_aStateVar[iTopSoil]/GetSoilCapacity(0),1.0),0.0);
      land_albedo =(soil_sat    )*_pSoil[0]->albedo_wet    +(1.0-soil_sat    )*_pSoil[0]->albedo_dry;
    }
    else if (_HRUType==HRU_GLACIER){
      land_albedo=0.4;//GLACIER_ALBEDO;
    }
    else if (_HRUType==HRU_LAKE){
      land_albedo=0.1;//WATER_ALBEDO;
    }
    else if (_HRUType==HRU_WATER){
      land_albedo=0.1;//WATER_ALBEDO;
    }
    else if (_HRUType==HRU_ROCK){
      land_albedo=0.35;//ROCK_ALBEDO;
    }
    else if(_HRUType==HRU_WETLAND){
      land_albedo=0.15;//WETLAND_ALBEDO
    }
    land_albedo =(snow_cover  )*snow_albedo+(1.0-snow_cover  )*land_albedo;
    //correction for urban surfaces?

    if (subcanopy) {
       return land_albedo;
    }
    else // above canopy - handles
    {
      veg_albedo    =_pVeg->albedo; //correction for wetness?

      //JRC: checks put in just in case parameters not supplied (only cosmetic in Forcings.csv, since if not supplied, SW_RADIA_NET not used in calcs).
      if (veg_albedo<0 ){veg_albedo=0.14;}
      if (svf>1.0      ){svf=0.0;}
      if (land_albedo<0){land_albedo=0.3;}

      return (1.0-svf)*(Fc)*veg_albedo+((svf)*(Fc)+(1.0-Fc))*land_albedo;
    }
  }
  //else if (Options.albedo_type==ALBEDO_LANDUSE){
  //ground_albedo= _pSurface->albedo;
  //return (snow_cover  )*snow_albedo             +(1.0-snow_cover  )*ground_albedo;
  //}
}
