/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2025 the Raven Development Team
  ----------------------------------------------------------------*/
#ifndef MODEL_H
#define MODEL_H

#include "RavenInclude.h"
#include "ModelABC.h"
#include "StateVariables.h"
#include "HydroProcessABC.h"
#include "LateralExchangeABC.h"
#include "SubBasin.h"
#include "HydroUnits.h"
#include "TimeSeries.h"
#include "Gauge.h"
#include "CustomOutput.h"
#include "GlobalParams.h"
#include "TransientParam.h"
#include "SoilAndLandClasses.h"
#include "SoilProfile.h"
#include "ChannelXSect.h"
#include "CustomOutput.h"
#include "Reservoir.h"
#include "Transport.h"
#include "Diagnostics.h"
#include "ForcingGrid.h"
#include "ModelEnsemble.h"
#include "GroundwaterModel.h"
#include "GWSWProcesses.h"
#include "ChannelXSect.h"
#include "Convolution.h"
#include "DemandOptimization.h"

class CHydroProcessABC;
class CGauge;
class CCustomOutput;
class CCustomTable;
class CGroundwaterModel;
class CTransportModel;
class CEnsemble;
class CForcingGrid;
class CLandUseClass;  // defined in 'SoilAndLandClasses.h'
class CSubbasinGroup; // defined in 'SubBasin.h'
class CChannelXSect;  // defined in 'ChannelXSect.h'
class CSubBasin;      // defined in 'SubBasin.h'
struct class_change;
class CTransientParam;
class CDemandOptimizer;

////////////////////////////////////////////////////////////////////
/// \brief Data abstraction for water surface model
/// \details Stores and organizes HRUs and basins, provides access to all
/// details and functionality in the main routine and solver functions
//
class CModel: public CModelABC
{
private:/*------------------------------------------------------*/

  //Model attributes
  double         _WatershedArea;  ///< total area of all subbasins [km^2]

  int              _nHydroUnits;  ///< number of HRUs in model
  CHydroUnit     **_pHydroUnits;  ///< Array of pointers to HRUs

  int               _nHRUGroups;  ///< number of HRU groups in model
  CHRUGroup       **_pHRUGroups;  ///< Array of pointers to HRU groups
  int                _nSBGroups;  ///< number of Subbasin groups in model
  CSubbasinGroup   **_pSBGroups;  ///< Array of pointers to Subbasin groups

  int               _nSubBasins;  ///< number of subbasins
  CSubBasin       **_pSubBasins;  ///< array of pointers to subbasins [size:_nSubBasins]; each subbasin includes multiple HRUs/HydroUnits
  int          *_aSubBasinOrder;  ///< stores order of subbasin for routing [size:_nSubBasins] (may be relegated to local variable in InitializeRoutingNetwork)
  int         _maxSubBasinOrder;  ///< stores maximum subasin order for routing (may be relegated to local variable in InitializeRoutingNetwork)
  int           *_aOrderedSBind;  ///< stores list of subbasin indices ordered upstream to downstream [size:_nSubBasins]
  int         *_aDownstreamInds;  ///< stores list of downstream indices of basins (for speed) [size:_nSubBasins]

  int               _nStateVars;  ///< number of state variables: water and energy storage units, snow density, etc.
  sv_type       *_aStateVarType;  ///< type of state variable in unit i  [size:_nStateVars]
  int          *_aStateVarLayer;  ///< index of state variable for multilayer variables (e.g., SOIL); [size:_nStateVars] value=DOESNT_EXIST(-1) for unique variables (e.g. SURFACE_WATER)
  int         _aStateVarIndices[MAX_STATE_VAR_TYPES][MAX_SV_LAYERS]; ///< lookup table for state variable indices; the index of SOIL[1] in a state_var[] array may be returned by aStateVarIndices[(int)(SOIL)][1]

  int                _nSoilVars;  ///< number of soil layer storage units

  int               _nProcesses;  ///< number of hydrological processes that move water, mass, or energy from one storage unit to another
  CHydroProcessABC**_pProcesses;  ///< Array of pointers to hydrological processes
  bool   **_aShouldApplyProcess;  ///< array of flags for whether or not each process applies to each HRU [_nProcesses][_nHydroUnits]
  int           _nConvVariables;  ///< Number of convolution variables (a.k.a. processes) in model

  int                  _nGauges;  ///< number of precip/temp gauges for forcing interpolation
  CGauge             **_pGauges;  ///< array of pointers to gauges which store time series info [size:_nGauges]
  double       **_aGaugeWeights;  ///< array of weights for each gauge/HRU pair [_nHydroUnits][_nGauges]
  double        **_aGaugeWtTemp;
  double      **_aGaugeWtPrecip;

  int            _nForcingGrids;  ///< number of gridded forcing input data
  CForcingGrid **_pForcingGrids;  ///< gridded input data [size: _nForcingGrids]

  int                 _UTM_zone;  ///< model-wide UTM zone used for interpolation

  int                 _lake_sv;   ///< index of storage variable for lakes/wetlands (TMP?)

  CGlobalParams     *_pGlobalParams;  ///< pointer to global parameters  (used to be global, static)

  int                 _nTransParams;  ///< number of transient parameters
  CTransientParam   **_pTransParams;  ///< array of pointers to transient parameters with time series
  int                _nClassChanges;  ///< number of HRU Group class changes
  class_change     **_pClassChanges;  ///< array of pointers to class_changes
  int              _nParamOverrides;  ///< number of local parameter overrides
  param_override **_pParamOverrides;  ///< array of pointers to local parameter overrides

  CGroundwaterModel  *_pGWModel;  ///< pointer to corresponding groundwater model
  CTransportModel *_pTransModel;  ///< pointer to corresponding transport model
  CEnsemble         *_pEnsemble;  ///< pointer to model ensemble
  CDemandOptimizer        *_pDO;  ///< pointer to demand optimizer (NULL w/o optimization)

  //For Diagnostics Calculation
  CTimeSeriesABC **_pObservedTS;  ///< array of pointers of observation time series [size: _nObservedTS]
  CTimeSeries     **_pModeledTS;  ///< array of pointers of modeled time series corresponding to observations [size: _nObservedTS]
  int              _nObservedTS;  ///< number of observation time series
  int               *_aObsIndex;  ///< index of the next unprocessed observation

  CTimeSeriesABC**_pObsWeightTS;  ///< array of pointers of observation weight time series [size: _nObsWeightTS]
  int             _nObsWeightTS;  ///< number of observation weight time series

  CDiagnostic   **_pDiagnostics;  ///< array of pointers to diagnostics to be calculated [size: _nDiagnostics]
  int             _nDiagnostics;  ///< number of diagnostics to be calculated comparing obs. vs modeled
  CDiagPeriod   **_pDiagPeriods;  ///< array of pointers to diagnostic periods [size _nDiagPeriods]
  int             _nDiagPeriods;  ///< number of diagnostic periods
  agg_diag   **_pAggDiagnostics;  ///< array of pointers to aggregate diagnostic structures [size: _nAggDiagnostics]
  int          _nAggDiagnostics;  ///< number of aggregated diagnostics

  //Data Assimilation
  double             *_aDAscale; ///< array of data assimilation flow scaling parameters [size: _nSubBasins] (NULL w/o DA)
  double            *_aDAlength; ///< array of downstream distance to nearest DA observation [m] [size: _nSubBasins] (NULL w/o DA)
  double           *_aDAQadjust; ///< array of flow adjustments [m3/s] [size: _nSubBasins]
  double          *_aDADrainSum; ///< sum of assimilated drainage areas upstream of a subbasin outlet [km2] [size: _nSubBasins]
  double           *_aDADownSum; ///< the drainage area of the nearest downstream assimilated flow observation [km2] [size: _nSubBasins]
  double         *_aDAtimesince; ///< array of downstream time since most recent downstream DA observation [size: _nSubBasins] (NULL w/o DA)
  bool            *_aDAoverride; ///< array of booleans indicating if observation data is available for assimilation at basin p's outlet [size: _nSubBasins] (NULL w/o DA)
  double              *_aDAobsQ; ///< array of observed flow values in basins [size: _nSubBasins]  (NULL w/o DA)
  double        *_aDAscale_last; ///< array of scale factors from previous time step  [size: _nSubBasins]  (NULL w/o DA)

  force_perturb**_pPerturbations; ///< array of pointers to perturbation data; defines which forcing functions to perturb and how [size: _nPerturbations]
  int            _nPerturbations; ///< number of forcing functions to perturb

  //Water/Energy Balance information
  double      **_aCumulativeBal;  ///< cumulative amount of flowthrough [mm or MJ/m2 or mg/m2] for each process connection, each HRU [k][j*]
  double            **_aFlowBal;  ///< current time step flowthrough [mm or MJ/m2 or mg/m2] for each process connection, each HRU [k][j*]
  int        _nTotalConnections;  ///< total number of in-HRU connections in model
  double    *_aCumulativeLatBal;  ///< cumulative amount of flowthrough [mm-m2 or MJ or mg] for each lateral process connection [j**]
  double          *_aFlowLatBal;  ///< current time step flowthrough [mm-m2 or MJ or mg] for each lateral process connection [j**]
  int     _nTotalLatConnections;  ///< total number of between-HRU connections in model
  double            _CumulInput;  ///< cumulative water added to watershed (precipitation, basin inflows, etc.) [mm]
  double           _CumulOutput;  ///< cumulative outflow of water from system [mm]
  double             _initWater;  ///< initial water in system [mm]

  //Output
  CCustomOutput**_pCustomOutputs; ///< Array of pointers to custom output objects [size:_nCustomOutputs]
  int            _nCustomOutputs; ///< Nuber of custom output objects
  CCustomTable  **_pCustomTables; ///< Array of pointers to custom table objects [size:_nCustomTables]
  int             _nCustomTables; ///< Number of custom tables
  ofstream                _HYDRO; ///< output file stream for Hydrographs.csv
  ofstream              _STORAGE; ///< output file stream for WatershedStorage.csv
  ofstream             _FORCINGS; ///< output file stream for ForcingFunctions.csv
  ofstream             _RESSTAGE; ///< output file stream for ReservoirStages.csv
  ofstream              _DEMANDS; ///< output file stream for Demands.csv
  ofstream               _LEVELS; ///< output file stream for WaterLevels.csv
  int                _HYDRO_ncid; ///< output file ID for Hydrographs.nc
  int             _RESSTAGE_ncid; ///< output file ID for ReservoirStages.nc
  int              _STORAGE_ncid; ///< output file ID for WatershedStorage.nc
  int             _FORCINGS_ncid; ///< output file ID for ForcingFunctions.nc
  int                _RESMB_ncid; ///< output file ID for ReservoirMassBalance.nc

  double          *_aOutputTimes; ///< array of model major output times (LOCAL times at which full solution is written)
  int              _nOutputTimes; ///< size of array of model major output times
  int         _currOutputTimeInd; ///< index of current output time
  const CHRUGroup *_pOutputGroup; ///< pointer to HRU group for which storage histories are written to file (or NULL)

  const optStruct   *_pOptStruct; ///< pointer to model options information

  //Blended PET/potential melt members
  int              _PETBlends_N;       ///< Not the best place to store these rarely used members \todo : move to globals
  evap_method     *_PETBlends_type;
  double          *_PETBlends_wts;
  int              _PotMeltBlends_N;
  potmelt_method  *_PotMeltBlends_type;
  double          *_PotMeltBlends_wts;

  //Parameter classes
  CLandUseClass    **_pLandUseClasses;     ///< array of pointers to land use classes
  int                _nLandUseClasses;     ///< number of land use classes
  CSoilClass       **_pAllSoilClasses;     /// used to be static attribute of CSoilClass
  int                _nAllSoilClasses;     /// same of above
  CVegetationClass **_pAllVegClasses;      /// used to be static attribute of CVegetationClass
  int                _numVegClasses;       /// same of above
  CTerrainClass    **_pAllTerrainClasses;  ///< array of pointers to all terrain classes that have been created
  int                _nAllTerrainClasses;  ///< Number of terrain classes that have been created length of pAllTerrainClasses
  CSoilProfile     **_pAllSoilProfiles;    ///< Reference to array of all soil profiles in model
  int                _nAllSoilProfiles;    ///< Number of soil profiles in model (size of pAllSoilProfiles)
  CChannelXSect    **_pAllChannelXSects;
  int                _nAllChannelXSects;

  CStateVariable    *_pStateVar;           ///< pointer to state variable information object
  int                _nLatFlowProcesses;   ///< number of lateral flow processes

  //initialization subroutines:
  void           GenerateGaugeWeights (double **&aWts, const forcing_type forcing, const optStruct 	 &Options);
  void       InitializeRoutingNetwork ();
  void         InitializeObservations (const optStruct 	 &Options);
  void     InitializeDataAssimilation (const optStruct   &Options);

  void      WriteEnsimStandardHeaders (const optStruct 	 &Options);
  void     WriteNetcdfStandardHeaders (const optStruct 	 &Options);
  void          WriteEnsimMinorOutput (const optStruct 	 &Options,
                                       const time_struct &tt);
  void         WriteNetcdfMinorOutput (const optStruct   &Options,
                                       const time_struct &tt);
  void    InitializeParameterOverrides();

  //private routines used during simulation:
  force_struct      GetAverageForcings() const;
  double       GetTotalChannelStorage () const;
  double      GetTotalReservoirStorage() const;
  double       GetTotalRivuletStorage () const;

  void                     CorrectPET(const optStruct &Options,
                                      force_struct &F,
                                      const CHydroUnit *pHRU,
                                      const double elev,
                                      const double ref_elev_temp,
                                      const int k);
  void                  CorrectPrecip(const optStruct &Options,
                                      force_struct &F,
                                      const double elev,
                                      const double ref_elev_precip,
                                      const int    k,
                                      const time_struct  &tt) const;
  void                    CorrectTemp(const optStruct &Options,
                                      force_struct &F,
                                      const double elev,
                                      const double ref_elev_temp,
                                      const time_struct &tt);
  double                  EstimatePET(const force_struct &F,
                                      const CHydroUnit   *pHRU,
                                      const double       &wind_measurement_ht,
                                      const double       &ref_elevation,
                                      const evap_method   evap_type,
                                      const optStruct    &Options,
                                      const time_struct  &tt,
                                      const bool          open_water);
  double         EstimateWindVelocity(const optStruct    &Options,
                                      const CHydroUnit   *pHRU,
                                      const force_struct &F,
                                      const double       &wind_measurement_ht,
                                      const int           k);
  double            WindspeedAtHeight(const double &z,
                                      const optStruct    &Options,
                                      const CHydroUnit   *pHRU,
                                      const force_struct &F,
                                      const double       &wind_measurement_ht);
  double           EstimateCloudCover(const optStruct    &Options,
                                      const force_struct &F,
                                      const int           k);
  double  CalculateSubDailyCorrection(const force_struct &F,
                                      const optStruct    &Options,
                                      const double       &elev,
                                      const double       &ref_elev_temp,
                                      const time_struct  &tt,
                                      const int          k);
  double        EstimatePotentialMelt(const force_struct *F,
                                      const potmelt_method method,
                                      const optStruct    &Options,
                                      const CHydroUnit   *pHRU,
                                      const time_struct  &tt);
  double         EstimateSnowFraction(const rainsnow_method  method,
                                      const CHydroUnit* pHRU,
                                      const force_struct* F,
                                      const optStruct& Options);

  double       CalculateAggDiagnostic(const int ii, const int j,
                                      const double &starttime, const double &endtime,
                                      const comparison compare,const double &thresh,
                                      const optStruct &Options);

  //Routines for deriving missing data based on gridded data provided

  CForcingGrid *ForcingCopyCreate(const CForcingGrid *pGrid, const forcing_type typ, const double &interval, const int nVals, const optStruct &Options);


  void         GenerateAveSubdailyTempFromMinMax        (const optStruct &Options);
  void         GenerateMinMaxAveTempFromSubdaily        (const optStruct &Options);
  void         GenerateMinMaxSubdailyTempFromAve        (const optStruct &Options);
  void         GeneratePrecipFromSnowRain               (const optStruct &Options);
  void         GenerateRainFromPrecip                   (const optStruct &Options);
  void         GenerateZeroSnow                         (const optStruct &Options);

  bool         ForcingGridIsInput                       (const forcing_type &ftype) const;
  bool         ForcingGridIsAvailable                   (const forcing_type &ftype) const;
  double       GetAverageSnowFrac                       (const int idx, const double t, const int n) const;

  void         AddFromPETParamList                      (string *aP,class_type *aPC,int &nP,
                                                         const evap_method &evaporation, const netSWRad_method &SW_radia_net) const;
  void         AddFromPotMeltParamList                  (string *aP,class_type *aPC,int &nP,
                                                         const potmelt_method &pot_melt) const;

public:/*-------------------------------------------------------*/

  //Constructor/Destructor:
  CModel(const int nsoillayers, const optStruct &Options);
  ~CModel();

  //Inherited Accessor functions (from ModelABC.h)
  bool              StateVarExists     (sv_type type) const;

  CGlobalParams    *GetGlobalParams    () const;
  int               GetNumStateVars    () const;
  sv_type           GetStateVarType    (const int i) const;
  int               GetStateVarIndex   (sv_type type) const; //assumes layer=0
  int               GetStateVarIndex   (sv_type type, int layer) const;//overriden for multilayer variables
  int               GetStateVarLayer   (const int i) const; //for multilayer variables
  double            GetConcentration   (const int k, const int i) const;

  double            GetFlux            (const int k, const int js, const optStruct &Options) const;
  double            GetLatFlow         (const int js, const optStruct &Options) const;
  double            GetCumulativeFlux  (const int k, const int i, const bool to) const;
  double            GetCumulFluxBetween(const int k,const int iFrom,const int iTo) const;

  double            GetAvgStateVar     (const int i) const;
  double            GetAvgConcentration(const int i) const;
  double            GetAvgForcing      (const forcing_type &ftype) const;
  double            GetAvgCumulFlux    (const int i, const bool to) const;
  double            GetAvgCumulFluxBet (const int iFrom, const int iTo) const;

  int               GetNumSoilLayers   () const;
  int               GetLakeStorageIndex() const;

  /* below are functions that were static in the past */
  // CLandUseClass
  CLandUseClass          *StringToLUClass(const string s);
  int                     GetNumLanduseClasses () const;
  CLandUseClass          *GetLanduseClass(int);
  void                    SummarizeLUClassesToScreen();
  void                    DestroyAllLanduseClasses();
  const CLandUseClass    *StringToLUClass  (const string s) const;
  const CLandUseClass    *GetLanduseClass  (const int    c) const;
  const int               GetLandClassIndex(const string s) const;
  // CSoilClass
  CSoilClass             *StringToSoilClass(const string s);
  int                     GetNumSoilClasses();
  const CSoilClass       *GetSoilClass(int c);
  void                    AddSoilClass(CSoilClass *pSoilClass);
  void                    SummarizeSoilClassesToScreen();
  void                    DestroyAllSoilClasses();
  // CVegetationClass
  CVegetationClass       *StringToVegClass(const string s);
  int                     GetNumVegClasses();
  const CVegetationClass *GetVegClass(int c);
  void                    AddVegClass(CVegetationClass *pVegClass);
  void                    SummarizeVegClassesToScreen();
  void                    DestroyAllVegClasses();
  // CTerrainClass
  CTerrainClass          *StringToTerrainClass(const string s);
  int                     GetNumTerrainClasses();
  const CTerrainClass    *GetTerrainClass(int c);
  void                    AddTerrainClass(CTerrainClass *pTerrainClass);
  void                    SummarizeTerrainClassesToScreen();
  void                    DestroyAllTerrainClasses();
  // CSoilProfile
  CSoilProfile           *StringToSoilProfile(const string s);
  int                     GetNumSoilProfiles();
  const CSoilProfile     *GetSoilProfile(int c);
  void                    AddSoilProfile(CSoilProfile *pSoilProfile);
  void                    SummarizeSoilProfilesToScreen();
  void                    DestroyAllSoilProfiles();
  // ChannelXSects
  CChannelXSect         *StringToChannelXSect(const string s);
  int                    GetNumChannelXSects();
  void                   AddChannelXSect(CChannelXSect *pXSect);
  void                   SummarizeChannelXSectToScreen();
  void                   DestroyAllChannelXSections();
  void                   CheckForChannelXSectsDuplicates(const optStruct &Options);
  void                   WriteRatingCurves(const optStruct& Options) const;
  void                   WriteRatingCurves(const optStruct* Options) const;
  // Convolution variables
  int                    GetNumConvolutionVariables() const;
  void                   IncrementConvolutionCount();
  // StateVariable
  CStateVariable        *GetStateVarInfo() const;
  void                   SetStateVarInfo(CStateVariable *pStateVar);
  // LateralExchangeABC
  int                    GetNumLatFlowProcesses();
  void                   CountOneMoreLatFlowProcess();

  /*--below are only available to global routines--*/
  //Accessor functions
  int               GetNumHRUs                        () const;
  int               GetNumHRUGroups                   () const;
  int               GetNumSubBasins                   () const;
  int               GetNumSubBasinGroups              () const;
  CHydroUnit       *GetHydroUnit                      (const int k ) const;
  CHydroUnit       *GetHRUByID                        (const long long int HRUID) const;
  CHRUGroup        *GetHRUGroup                       (const int kk) const;
  CHRUGroup        *GetHRUGroup                       (const string name) const;
  CSubBasin        *GetSubBasin                       (const int p ) const;
  CSubBasin        *GetSubBasinByID                   (const long long ID) const;
  CSubbasinGroup   *GetSubBasinGroup                  (const int pp) const;
  CSubbasinGroup   *GetSubBasinGroup                  (const string name) const;
  CHydroProcessABC *GetProcess                        (const int j ) const;
  CGauge           *GetGauge                          (const int g) const;
  CForcingGrid     *GetForcingGrid                    (const forcing_type &ftype) const;
  int               GetNumGauges                      () const;
  int               GetNumForcingGrids                () const;
  int               GetNumProcesses                   () const;
  process_type      GetProcessType                    (const int j ) const;
  int               GetNumConnections                 (const int j ) const;
  int               GetNumForcingPerturbations        () const;
  double            GetAveragePrecip                  () const;
  double            GetAverageSnowfall                () const;
  int               GetOrderedSubBasinIndex           (const int pp) const;
  int               GetDownstreamBasin                (const int p ) const;
  int               GetSubBasinIndex                  (const long long SBID) const;
  int               GetGaugeIndexFromName             (const string name) const;
  int               GetForcingGridIndexFromType       (const forcing_type &ty) const;
  const CSubBasin **GetUpstreamSubbasins              (const long long SBID, int &nUpstream) const;
  bool              IsSubBasinUpstream                (const long long SBID,const long long SBIDdown) const;

  double            GetWatershedArea                  () const;
  bool              IsInHRUGroup                      (const int k,
                                                       const string HRUGroupName) const;
  bool              IsInSubBasinGroup                 (const long long SBID,
                                                       const string SBGroupName) const;

  int                   GetNumObservedTS              () const;
  const CTimeSeriesABC *GetObservedTS                 (const int i) const;
  const CTimeSeriesABC* GetSimulatedTS                (const int i) const;

  double               GetObjFuncVal                  (long long calib_SBID,diag_type calib_Obj, const string calib_period) const;

  const optStruct     *GetOptStruct                   () const;
  CTransportModel     *GetTransportModel              () const;
  CGroundwaterModel   *GetGroundwaterModel            () const;
  CEnsemble           *GetEnsemble                    () const;
  CDemandOptimizer    *GetManagementOptimizer         () const;

  void              GetParticipatingParamList         (string *aP,
                                                       class_type *aPC,
                                                       int &nP,
                                                       const optStruct &Options) const;
  class_type        ParamNameToParamClass             (const string param_str, const string class_name) const;

  //Manipulator Functions: called by Parser
  void    AddProcess                (        CHydroProcessABC  *pMov            );
  void    AddHRU                    (        CHydroUnit        *pHRU            );
  void    AddHRUGroup               (        CHRUGroup         *pHRUGrp         );
  void    AddSubBasin               (        CSubBasin         *pWS             );
  void    AddSubBasinGroup          (        CSubbasinGroup    *pSBGrp          );
  void    AddGauge                  (        CGauge            *pGage           );
  void    AddForcingGrid            (        CForcingGrid      *pGrid           , forcing_type typ);
  void    AddStateVariables         (const sv_type             *aSV,
                                     const int                 *aLev,
                                     const int                  nSV             );
  void    AddLandUseClass           (       CLandUseClass      *pLU             );
  void    AddCustomOutput           (      CCustomOutput       *pCO             );
  void    AddCustomTable            (      CCustomTable        *pTab            );
  void    AddTransientParameter     (      CTransientParam     *pTP             );
  void    AddParameterOverride      (      param_override      *pPO             );
  void    AddPropertyClassChange    (const string             HRUgroup,
                                     const class_type         tclass,
                                     const string             new_class,
                                     const time_struct       &tt,
                                     const optStruct         &Options           );
  void    AddObservedTimeSeries     (        CTimeSeriesABC    *pTS             );
  void    AddObservedWeightsTS      (        CTimeSeriesABC    *pTS             );
  void    AddDiagnostic             (        CDiagnostic       *pDiag           );
  void    AddDiagnosticPeriod       (        CDiagPeriod       *pDiagPer        );
  void    AddAggregateDiagnostic    (agg_stat stat, string datatype, int group_ind);
  void    AddForcingPerturbation    (forcing_type type, disttype distrib, double* distpars, int group_index, adjustment adj, int nStepsPerDay);

  void    AddModelOutputTime        (const time_struct       &tt_out,
                                     const optStruct         &Options           );
  void    SetLakeStorage            (const sv_type            SV,
                                     const int                lev               );//TMP?
  void    SetOutputGroup            (const CHRUGroup         *pOut              );
  void    SetNumSnowLayers          (const int                nLayers           );

  void    OverrideStreamflow        (const long long SBID);
  void    SetEnsembleMode           (CEnsemble *pEnsemble);
  void    AddDemandOptimization     (CDemandOptimizer *pDO);

  void    SetPETBlendValues         (const int N, const evap_method    *aET, const double *wts);
  void    SetPotMeltBlendValues     (const int N, const potmelt_method *aPM, const double *wts);

  int     GetBlendedForcingsNumWeights(const string label);
  void    DeleteCustomOutputs       ();

  /*--Other Functions: mostly called by Solver--*/
  //called only once prior to simulation:
  void        Initialize                 (const optStruct &Options);
  void        InitializeBasins           (const optStruct &Options,const bool re_init);
  void        InitializePostRVM          (const optStruct &Options);
  void        WriteOutputFileHeaders     (const optStruct &Options);
  void        GenerateGriddedPrecipVars  (const optStruct &Options);
  void        GenerateGriddedTempVars    (const optStruct &Options);
  void        ClearTimeSeriesData        (const optStruct &Options);

  //called during simulation:
  void        UpdateParameter            (const class_type  &ctype,
                                          const string      pname,
                                          const string      cname,
                                          const double      &value);
  void        ApplyLocalParamOverrrides  (const int         k,
                                          const bool        revert);

  //called during simulation:
  //critical simulation routines (called once during each timestep):
  void        UpdateTransientParams      (const optStruct   &Options,
                                          const time_struct &tt);
  void        UpdateHRUForcingFunctions  (const optStruct   &Options,
                                          const time_struct &tt); //declaration in UpdateForcings.cpp
  void        UpdateDiagnostics          (const optStruct   &Options,
                                          const time_struct &tt);
  void        RecalculateHRUDerivedParams(const optStruct   &Options,
                                          const time_struct &tt);
  bool        ApplyProcess               (const int          j,
                                          const double      *state_var,
                                          const CHydroUnit  *pHRU,
                                          const optStruct   &Options,
                                          const time_struct &tt,
                                                int         *iFrom,
                                                int         *iTo,
                                                int         &nConnections,
                                                double      *rates_of_change) const;
  bool        ApplyLateralProcess        (const int          j,
                                          const double* const* state_vars,
                                          const optStruct   &Options,
                                          const time_struct &tt,
                                                int         *kFrom,
                                                int         *kTo,
                                                int         *iFrom,
                                                int         *iTo,
                                                int         &nLatConnections,
                                                double      *exchange_rates) const;
  void         AssimilationOverride      (const int p,
                                          const optStruct &Options, const time_struct &tt);
  void         PrepareAssimilation       (const optStruct &Options, const time_struct &tt);
  void         AssimilationBackPropagate (const optStruct &Options, const time_struct &tt);

  void         PrepareForcingPerturbation(const optStruct &Options, const time_struct &tt);
  void         ApplyForcingPerturbation  (const forcing_type f, force_struct &F, const int k, const optStruct& Options, const time_struct& tt);

  //water/energy/mass balance routines
  void   CalculateInitialWaterStorage (const optStruct   &Options);
  void        IncrementBalance        (const int q_star,
                                       const int k,
                                       const double moved);//[mm] or [MJ/m2] or [mg]
  void        IncrementLatBalance     (const int j_star,
                                       const double moved);//[mm] or [MJ/m2] or [mg]
  void        IncrementCumulInput     (const optStruct &Options, const time_struct &tt);
  void        IncrementCumOutflow     (const optStruct &Options, const time_struct &tt);

  //output routines
  void        WriteMinorOutput        (const optStruct &Options, const time_struct &tt);
  void        WriteSimpleOutput       (const optStruct &Options, const time_struct &tt);
  void        WriteMajorOutput        (const time_struct &tt,string solfile,bool final) const;
  void        WriteMajorOutput        (const optStruct &Options, const time_struct &tt,string solfile,bool final) const;
  void        WriteProgressOutput     (const optStruct &Options, clock_t elapsed_time, int elapsed_steps, int total_steps);
  void        CloseOutputStreams      ();
  void        SummarizeToScreen       (const optStruct &Options) const;
  void        RunDiagnostics          (const optStruct &Options);
};

#endif
