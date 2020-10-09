/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2019 the Raven Development Team
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
#include "GroundwaterClass.h"

class CHydroProcessABC;
class CGauge;
class CCustomOutput;
class CGroundwaterModel;
class CTransportModel;
class CEnsemble;
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
  int         _aStateVarIndices[MAX_STATE_VARS][MAX_SV_LAYERS]; ///< lookup table for state variable indices; the index of SOIL[1] in a state_var[] array may be returned by aStateVarIndices[(int)(SOIL)][1]

  int                _nSoilVars;  ///< number of soil layer storage units
  int           _nAquiferLayers;  ///< number of aquifer layer storage units
  int              _nSnowLayers;  ///< number of snow layer storage units (default: 1)

  int               _nProcesses;  ///< number of hydrological processes that move water, mass, or energy from one storage unit to another
  CHydroProcessABC**_pProcesses;  ///< Array of pointers to hydrological processes

  bool   **_aShouldApplyProcess; ///< array of flags for whether or not each process applies to each HRU [_nProcesses][_nHydroUnits]

  int                  _nGauges;  ///< number of precip/temp gauges for forcing interpolation
  CGauge             **_pGauges;  ///< array of pointers to gauges which store time series info [size:_nGauges]
  double       **_aGaugeWeights;  ///< array of weights for each gauge/HRU pair [_nHydroUnits][_nGauges]
  double        **_aGaugeWtTemp;
  double      **_aGaugeWtPrecip;

  int            _nForcingGrids;  ///< number of gridded forcing input data
  CForcingGrid **_pForcingGrids;  ///< gridded input data [size: _nForcingGrids]

  int                 _UTM_zone;  ///< model-wide UTM zone used for interpolation

  int                 _lake_sv;   ///< index of storage variable for lakes/wetlands (TMP?)

  int             _nTransParams;  ///< number of transient parameters
  CTransientParam**_pTransParams; ///< array of pointers to transient parameters with time series
  int            _nClassChanges;  ///< number of HRU Group class changes
  class_change **_pClassChanges;  ///< array of pointers to class_changes

  CGroundwaterModel  *_pGWModel;  ///< pointer to corresponding groundwater model
  CTransportModel *_pTransModel;  ///< pointer to corresponding transport model
  CEnsemble         *_pEnsemble;  ///< pointer to model ensemble

  //For Diagnostics Calculation
  CTimeSeriesABC **_pObservedTS;  ///< array of pointers of observation time series [size: _nObservedTS]
  CTimeSeries     **_pModeledTS;  ///< array of pointers of modeled time series corresponding to observations [size: _nObservedTS]
  int              _nObservedTS;  ///< number of observation time series
  int               *_aObsIndex;  ///< index of the next unprocessed observation

  CTimeSeriesABC**_pObsWeightTS;  ///< array of pointers of observation weight time series [size: _nObsWeightTS]
  int             _nObsWeightTS;  ///< number of observation weight time series

  CDiagnostic   **_pDiagnostics;  ///< array of pointers to diagnostics to be calculated [size: _nDiagnostics]
  int             _nDiagnostics;  ///< number of diagnostics to be calculated comparing obs. vs modeled

  //Data Assimilation 
  double             *_aDAscale; ///< array of data assimilation flow scaling parameters [size: _nSubBasins] (NULL w/o DA)
  double            *_aDAlength; ///< array of downstream distance to nearest DA observation [m] [size: _nSubBasins] (NULL w/o DA)
  double         *_aDAtimesince; ///< array of downstream time since most recent downstream DA observation [size: _nSubBasins] (NULL w/o DA)
  
  //Water/Energy Balance information
  double      **_aCumulativeBal;  ///< cumulative amount of flowthrough [mm or MJ/m2 or mg/m2] for each process connection, each HRU [k][j*]
  double            **_aFlowBal;  ///< current time step flowthrough [mm or MJ/m2 or mg/m2] for each process connection, each HRU [k][j*]
  int        _nTotalConnections;  ///< total number of in-HRU connections in model
  double    *_aCumulativeLatBal;  ///< cumulative amount of flowthrough [mm-m2 or MJ or mg] for each lateral process connection [j**]
  double          *_aFlowLatBal;  ///< current time step flowthrough [mm-m2 or MJ or mg] for each lateral process connection [j**]
  int     _nTotalLatConnections;  ///< total number of between-HRU connections in model

  double            _CumulInput;  ///< cumulative water added to watershed (precipitation, basin inflows, etc.) [mm]
  double           _CumulOutput;  ///< cumulative outflow of water from system [mm]
  double         _CumEnergyGain;  ///< cumulative area-averaged average energy gain [MJ/m2]
  double         _CumEnergyLoss;  ///< cumulative area-averaged energy loss [MJ/m2]
  double             _initWater;  ///< initial water in system [mm]

  //Output
  CCustomOutput**_pCustomOutputs; ///< Array of pointers to custom output objects
  int            _nCustomOutputs; ///< Nuber of custom output objects
  ofstream                _HYDRO; ///< output file stream for Hydrographs.csv
  ofstream              _STORAGE; ///< output file stream for WatershedStorage.csv
  ofstream             _FORCINGS; ///< output file stream for ForcingFunctions.csv
  ofstream             _RESSTAGE; ///< output file stream for ReservoirStages.csv
  ofstream              _DEMANDS; ///< output file stream for Demands.csv
  int                _HYDRO_ncid; ///< output file ID  for Hydrographs.nc;     
  int              _STORAGE_ncid; ///< output file ID  for WatershedStorage.nc; 
  int             _FORCINGS_ncid; ///< output file ID  for ForcingFunctions.nc;   
  double          *_aOutputTimes; ///< array of model major output times (LOCAL times at which full solution is written)
  int              _nOutputTimes; ///< size of array of model major output times
  int         _currOutputTimeInd; ///< index of current output time
  const CHRUGroup *_pOutputGroup; ///< pointer to HRU group for which storage histories are written to file (or NULL)

  const optStruct   *_pOptStruct; ///< pointer to model options information

  //initialization subroutines:
  void           GenerateGaugeWeights (double **&aWts, const forcing_type forcing, const optStruct 	 &Options);
  void       InitializeRoutingNetwork ();
  void           InitializeBasinFlows (const optStruct 	 &Options);
  void         InitializeObservations (const optStruct 	 &Options);
  void     InitializeDataAssimilation (const optStruct   &Options);

  void      WriteEnsimStandardHeaders (const optStruct 	 &Options);
  void     WriteNetcdfStandardHeaders (const optStruct 	 &Options);
  void          WriteEnsimMinorOutput (const optStruct 	 &Options,
                                       const time_struct &tt);
  void         WriteNetcdfMinorOutput (const optStruct   &Options,
                                       const time_struct &tt);

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
  double         EstimateWindVelocity(const optStruct &Options,
                                      const force_struct &F,
                                      const double &wind_measurement_ht,
                                      const int k);
  double         EstimateCloudCover  (const optStruct &Options,
                                      const force_struct &F,
                                      const int k);
  double  CalculateSubDailyCorrection(const force_struct &F,
                                      const optStruct    &Options,
                                      const double       &elev,
                                      const double       &ref_elev_temp,
                                      const time_struct  &tt,
                                      const int          k);
  double        EstimatePotentialMelt(const force_struct *F,
                                      const optStruct &Options,
                                      const CHydroUnit *pHRU,
                                      const time_struct &tt);


  //Routines for deriving missing data based on gridded data provided

  CForcingGrid *ForcingCopyCreate(const CForcingGrid *pGrid, const forcing_type typ, const double &interval, const int nVals);


  void         GenerateAveSubdailyTempFromMinMax        (const optStruct &Options);
  void         GenerateMinMaxAveTempFromSubdaily        (const optStruct &Options);
  void         GenerateMinMaxSubdailyTempFromAve        (const optStruct &Options);
  void         GeneratePrecipFromSnowRain               (const optStruct &Options);
  void         GenerateRainFromPrecip                   (const optStruct &Options);
  void         GenerateZeroSnow                         (const optStruct &Options);

  bool         ForcingGridIsInput                       (const forcing_type &ftype) const;
  bool         ForcingGridIsAvailable                   (const forcing_type &ftype) const;
  double       GetAverageSnowFrac                       (const int idx, const double t, const int n) const;


public:/*-------------------------------------------------------*/

  //Constructor/Destructor:
  CModel(const soil_model SM, const int nsoillayers, const optStruct &Options);
  ~CModel();

  //Inherited Accessor functions (from ModelABC.h)
  bool              StateVarExists     (sv_type type) const;

  int               GetNumStateVars    () const;
  sv_type           GetStateVarType    (const int i) const;
  int               GetStateVarIndex   (sv_type type) const; //assumes layer=0
  int               GetStateVarIndex   (sv_type type, int layer) const;//overriden for multilayer variables
  int               GetStateVarLayer   (const int i) const; //for multilayer variables

  double            GetFlux            (const int k, const int js, const optStruct &Options) const;
  double            GetLatFlow         (const int js, const optStruct &Options) const;
  double            GetCumulativeFlux  (const int k, const int i, const bool to) const;
  double            GetCumulFluxBetween(const int k,const int iFrom,const int iTo) const;

  double            GetAvgStateVar     (const int i) const;
  double            GetAvgForcing      (const string &forcing_string) const;
  double            GetAvgCumulFlux    (const int i, const bool to) const;
  double            GetAvgCumulFluxBet (const int iFrom, const int iTo) const;

  int               GetNumSoilLayers   () const;
  int               GetNumAquiferLayers() const;
  int               GetLakeStorageIndex() const; //TMP?

  /*--below are only available to global routines--*/
  //Accessor functions
  int               GetNumHRUs                        () const;
  int               GetNumHRUGroups                   () const;
  int               GetNumSubBasins                   () const;
  int               GetNumSubBasinGroups              () const;
  CHydroUnit       *GetHydroUnit                      (const int k ) const;
  CHydroUnit       *GetHRUByID                        (const int HRUID) const;
  CHRUGroup        *GetHRUGroup                       (const int kk) const;
  CHRUGroup        *GetHRUGroup                       (const string name) const;
  CSubBasin        *GetSubBasin                       (const int p ) const;
  CSubBasin        *GetSubBasinByID                   (const long ID) const;
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
  double            GetAveragePrecip                  () const;
  double            GetAverageSnowfall                () const;
  int               GetOrderedSubBasinIndex           (const int pp) const;
  int               GetDownstreamBasin                (const int p ) const;
  int               GetSubBasinIndex                  (const long ID) const;
  int               GetGaugeIndexFromName             (const string name) const;
  int               GetForcingGridIndexFromType       (const forcing_type &ty) const;
  const CSubBasin **GetUpstreamSubbasins              (const int SBID, int &nUpstream) const;

  double            GetWatershedArea                  () const;
  bool              IsInHRUGroup                      (const int k,
                                                       const string HRUGroupName) const;
  bool              IsInSubBasinGroup                 (const long SBID,
                                                       const string SBGroupName) const;

  const optStruct  *GetOptStruct                      () const;
  CTransportModel  *GetTransportModel                 () const;
  CGroundwaterModel*GetGroundwaterModel               () const;
  CEnsemble        *GetEnsemble                       () const;

  void              GetParticipatingParamList         (string *aP,
                                                       class_type *aPC,
                                                       int &nP,
                                                       const optStruct &Options) const;

  //Manipulator Functions: called by Parser
  void    AddProcess                (        CHydroProcessABC  *pMov            );
  void    AddHRU                    (        CHydroUnit        *pHRU            );
  void    AddHRUGroup               (        CHRUGroup         *pHRUGrp         );
  void    AddSubBasin               (        CSubBasin         *pWS             );
  void    AddSubBasinGroup          (        CSubbasinGroup    *pSBGrp          );
  void    AddGauge                  (        CGauge            *pGage           );
  void    AddForcingGrid            (        CForcingGrid      *pGrid           , forcing_type typ);
  void    AddStateVariables         (const sv_type           *aSV,
                                     const int               *aLev,
                                     const int                nSV               );
  void    AddAquiferStateVars       (const int                nLayers           );
  void    AddCustomOutput           (        CCustomOutput     *pCO             );
  void    AddTransientParameter     (        CTransientParam   *pTP             );
  void    AddPropertyClassChange    (const string             HRUgroup,
                                     const class_type         tclass,
                                     const string             new_class,
                                     const time_struct       &tt,
                                     const optStruct         &Options           );
  void    AddObservedTimeSeries     (        CTimeSeriesABC    *pTS             );
  void    AddObservedWeightsTS      (        CTimeSeriesABC    *pTS             );
  void    AddDiagnostic             (        CDiagnostic       *pDiag           );
  void    AddModelOutputTime        (const time_struct       &tt_out,
                                     const optStruct         &Options           );
  void    SetLakeStorage            (const sv_type            SV,
                                     const int                lev               );//TMP?
  void    SetAggregatedVariable     (const sv_type            SV,
                                     const int                lev,
                                     const string             group_name        );
  void    SetOutputGroup            (const CHRUGroup         *pOut              );
  void    SetNumSnowLayers          (const int                nLayers           );

  void    OverrideStreamflow        (const long SBID);
  void    SetEnsembleMode           (CEnsemble *pEnsemble);

  /*--Other Functions: mostly called by Solver--*/
  //called only once prior to simulation:
  void        Initialize                 (const optStruct &Options);
  void        WriteOutputFileHeaders     (const optStruct &Options);
  void        GenerateGriddedPrecipVars  (const optStruct &Options);
  void        GenerateGriddedTempVars    (const optStruct &Options);

  //called during simulation:
  void        UpdateParameter            (const class_type  &ctype,
                                          const string      pname,
                                          const string      cname,
                                          const double      &value);
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

  void         AssimilateStreamflow      (const optStruct &Options,
                                          const time_struct &tt);

  //water/energy/mass balance routines
  void        IncrementBalance        (const int q_star,
                                       const int k,
                                       const double moved);//[mm] or [MJ/m2] or [mg]
  void        IncrementLatBalance     (const int j_star,  
                                       const double moved);//[mm] or [MJ/m2] or [mg]
  void        IncrementCumulInput     (const optStruct &Options, const time_struct &tt);
  void        IncrementCumOutflow     (const optStruct &Options);

  //output routines
  void        WriteMinorOutput        (const optStruct &Options, const time_struct &tt);
  void        WriteMajorOutput        (string solfile, const optStruct &Options, const time_struct &tt, bool final) const;
  void        WriteProgressOutput     (const optStruct &Options, clock_t elapsed_time, int elapsed_steps, int total_steps);
  void        CloseOutputStreams      ();
  void        SummarizeToScreen       (const optStruct &Options) const;
  void        RunDiagnostics          (const optStruct &Options);
};

#endif
