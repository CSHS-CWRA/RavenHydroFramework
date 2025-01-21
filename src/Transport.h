/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2022 the Raven Development Team
----------------------------------------------------------------*/

#ifndef TRANSPORTMODEL_H
#define TRANSPORTMODEL_H

#include "RavenInclude.h"
#include "Model.h"

enum constit_type {
  AQUEOUS,
  SORBED,
  ENTHALPY,
  ISOTOPE,
  TRACER
};

struct constit_source
{
  bool   dirichlet;       ///< =true for dirichlet, false for neumann
  int    constit_index;   ///< constituent index (c)
  int    i_stor;          ///< index of water storage compartment
  int    kk;              ///< index of HRU group to which source is applied (default is all for DOESNT_EXIST)
  double concentration;   ///< fixed concentration [mg/m2] or enthalpy [MJ/m2] (or DOESNT_EXIST=-1 if time series should be used)
  double concentration2;  ///< fixed concentration [mg/m2] or enthalpy [MJ/m2] for blend-type condition
  double flux;            ///< fixed mass flux [mg/m2/d] or heat flux [MJ/m2/d] (or DOESNT_EXIST=-1 if time series should be used)
  const  CTimeSeries *pTS; ///< time series of fixed concentration or mass/heat flux (or NULL if fixed should be used)
};

enum gparam_type
{
  PAR_DECAY_COEFF,        ///< linear decay coefficient  [1/d]
  PAR_MASS_LOSS_RATE,     ///< zero-order decay or loss [mg/m2/d]
  PAR_TRANSFORM_COEFF,    ///< linear transformation coeff [1/d for n=1, (mg/L)^-n /d for n!=1]
  PAR_STOICHIO_RATIO,     ///< ratio of mass reactant to mass product [mg/mg]
  PAR_UPTAKE_COEFF,       ///< uptake coefficeint [-]
  PAR_TRANSFORM_N,        ///< transformation exponent [-]
  PAR_EQFIXED_RATIO,      ///< fixed equilibrium coefficient [-] (B=cc*A)
  PAR_EQUIL_COEFF,        ///< linear equilibrium coefficient [1/d] k*(B-cc*A)
  PAR_SORPT_COEFF         ///< sorption coefficient [l/kg]
};
struct geochem_param
{
  gparam_type type;         ///< parameter type
  double      value;        ///< parameter value (units determined by type)
  int         constit_ind;  ///< constituent index
  int         constit_ind2; ///< 2nd constituent index (for reactions/transformations)
  int         process_ind;  ///< process index
  string      class_name;   ///< soil/veg/LU class name
  int         ii_water_stor; ///< water storage compartment index ii

  geochem_param(gparam_type ty) { //default constructor
    type         =ty;
    value        =0.0;
    constit_ind  =0;
    constit_ind2 =DOESNT_EXIST;
    process_ind  =DOESNT_EXIST;
    ii_water_stor=DOESNT_EXIST;
    class_name   ="";
  }
};

///////////////////////////////////////////////////////////////////
/// \brief Class for coordinating transport simulation
/// \details Implemented in Transport.cpp
//
class CModel;
class CConstituentModel;
class CEnthalpyModel;

class CTransportModel
{
private:/*------------------------------------------------------*/
  CModel *pModel;                    ///< pointer to model object

  int  _nAdvConnections;             ///< number of advective/dispersive transport connections between water storage units
  int *_iFromWater;                  ///< state variable indices of water compartment source      [size: nAdvConnections]
  int *_iToWater;                    ///< state variable indices of water compartment destination [size: nAdvConnections]
  int *_js_indices;                  ///< process index (j*) of connection                        [size: nAdvConnections]

  int  _nLatConnections;             ///< number of lateral advective transport connections between water storage units
  int  *_iLatFromWater;              ///< state variable indices of water compartment source      [size: _nLatConnections]
  int  *_iLatToWater;                ///< state variable indices of water compartment destination [size: _nLatConnections]
  int  *_iLatFromHRU;                ///< state variable indices of water compartment source      [size: _nLatConnections]
  int  *_iLatToHRU;                  ///< state variable indices of water compartment destination [size: _nLatConnections]
  int  *_latqss_indices;             ///< process index (q**) of connection                       [size: _nLatConnections]

  int  _nWaterCompartments;          ///< number of water storage compartments which may contain constituent
  int *_iWaterStorage;               ///< state variable indices of water storage compartments which may contain constituent [size: _nWaterCompartments]
  int *_aIndexMapping;               ///< lookup table to convert global state variable index i to local water storage index ii [size: pModel::_nStateVars prior to transport variables being included]
                                     ///< basically inverse of _iWaterStorage
  int  _nIndexMapping;               ///< size of index mapping array [=pModel::_nStateVars prior to transport variables being included]

  string         *_aProcessNames;    ///< list of recognized geochemical process names  [size: _nProcessNames]
  int             _nProcessNames;    ///< size of process name list

  geochem_param **_pGeochemParams;   ///< array of pointers to geochem parameter structures [size: _nTransParams]
  int             _nGeochemParams;   ///< number of geochemistry paramers

  int                 _nConstituents;  ///< number of transported constituents [c=0.._nConstituents-1]
  CConstituentModel **_pConstitModels; ///< array of pointers to constituent models [size:_nConstituents]
  CEnthalpyModel     *_pEnthalpyModel; ///< special pointer to enthalpy constituent, if present (or NULL if not)

  void m_to_cj(const int layerindex,int &c,int &j) const;

  string GetConstituentLongName_loc(const int layerindex) const; //e.g., "Nitrogen in Soil Water[2]"


public:/*-------------------------------------------------------*/
  CTransportModel(CModel *pMod);
  ~CTransportModel();

  //Accessors
  int          GetLayerIndexFromName2(const string name,const int comp_m) const;     //non-static version

  string GetConstituentTypeName(const int layerindex);//e.g., "Nitrogen" (m=layer index)
  string GetConstituentTypeName2(const int c) const;         //e.g., "Nitrogen" (c=constituent index)
  string GetConstituentLongName(const int layerindex) const;//e.g., "Nitrogen in Soil Water[2]"
  string GetConstituentShortName(const int layerindex) const; //e.g., "!Nitrogen_SOIL[2]"

  int    GetNumConstituents() const;
  int    GetConstituentIndex(const string name) const;

  CConstituentModel *GetConstituentModel(const int c)  const;
  CConstituentModel *GetConstituentModel2(const int c) const;

  int    GetNumWaterCompartments() const;
  int    GetNumAdvConnections() const;
  int    GetNumLatAdvConnections() const;

  int    GetFromIndex   (const int c,const int q) const;
  int    GetToIndex     (const int c,const int q) const;
  int    GetStorIndex   (const int c,const int ii) const;
  int    GetLatFromIndex(const int c,const int qq) const;
  int    GetLatToIndex  (const int c,const int qq) const;

  int    GetFromWaterIndex   (const int q) const;
  int    GetToWaterIndex     (const int q) const;
  int    GetJsIndex          (const int q) const;
  int    GetLatFromHRU       (const int qq) const;
  int    GetLatToHRU         (const int qq) const;
  int    GetLatFromWaterIndex(const int qq) const;
  int    GetLatToWaterIndex  (const int qq) const;
  int    GetLatqsIndex       (const int qq) const;

  int    GetStorWaterIndex(const int ii) const;
  int    GetWaterStorIndexFromLayer(const int m) const;
  int    GetWaterStorIndexFromSVIndex(const int i) const;

  int    GetLayerIndex(const int c,const int i_stor) const;

  double GetGeochemParam    (const gparam_type gtyp, const int c, const int ii, const int procind, const CHydroUnit *pHRU) const;
  double GetGeochemParam    (const gparam_type gtyp, const int c, const int c2, const int ii,const int procind,const CHydroUnit* pHRU) const;

  const CEnthalpyModel *GetEnthalpyModel() const;

  double GetConcentration           (const int k,const int sv_index) const;

  double GetAdvectionCorrection     (const int c,const CHydroUnit* pHRU,const int iFromWater,const int iToWater,const double& C) const;

  int    GetProcessIndex            (const string name) const;

  //Manipulators
  void   AddConstituent             (string name,constit_type type,bool is_passive);
  void   AddProcessName             (string name);
  void   AddGeochemParam            (const gparam_type typ,const string constit_name,const string constit_name2,
                                     const string procname,const int watcompart, const string pclass, const double &value);

  //Routines called during model execution
  void   Prepare                    (const optStruct &Options);
  void   Initialize                 (const optStruct& Options);
  void   InitializeParams           (const optStruct& Options); //called at end of .rvi file read
  void   CalculateLateralConnections();

  void   IncrementCumulInput        (const optStruct &Options,const time_struct &tt);
  void   IncrementCumulOutput       (const optStruct &Options);

  void   WriteOutputFileHeaders     (const optStruct &Options) const;
  void   WriteMinorOutput           (const optStruct &Options,const time_struct &tt) const;
  void   WriteMajorOutput           (ofstream& RVC) const;
  void   CloseOutputFiles           () const;
};
///////////////////////////////////////////////////////////////////
/// \brief Class for coordinating transport simulation for specific constituent
/// \details Implemented in ConstituentModel.cpp
//
class CConstituentModel
{
protected:
  const CModel               *_pModel;
  const CTransportModel *_pTransModel;

  string                     _name;  ///< constituent name (e.g., "Nitrogen")
  int               _constit_index;  ///< master constituent index, c,  of this constituent
  constit_type               _type;  ///< AQUEOUS [mg], ENTHALPY [MJ], ISOTOPE [mg/mg] or TRACER [-]
  bool                 _is_passive;  ///< doesn't transport via advection (default: false)

  bool         _lateral_via_convol;  ///< true if lateral mass/energy transfer handled within routing convolution (default: false)

  // Routing/state var storage
  int                   *_nMinHist;  ///< size of upstream loading history in each basin [size: nSubBasins]
  double               **_aMinHist;  ///< array used for storing routing upstream loading history [mg/d] or [MJ/d] [size: nSubBasins x _nMinHist[p]]
  int                  *_nMlatHist;  ///< size of lateral loading history in each basin [size: nSubBasins]
  double             ** _aMlatHist;  ///< array used for storing routing lateral loading history [mg/d] or [MJ/d] [size: nSubBasins  x _nMlatHist[p]]
  double                  **_aMout;  ///< array storing current mass flow at points along channel [mg/d] or [MJ/d] [size: nSubBasins x _nSegments(p)]
  double              *_aMout_last;  ///< array used for storing mass outflow from channel at start of timestep [mg/d] or [MJ/d] [size: nSubBasins ]
  double              *_aMlat_last;  ///< array storing mass/energy outflow from start of timestep [size: nSubBasins]
  double                 *_aMlocal;  ///< local contribution to current subbasin mass/energy outflow [mg/d] or [MJ/d] [size: nSubBasins] (just from in-catchment mass/energy routing)
  double               *_aMlocLast;  ///< last local contribution [mg/d] or [MJ/d] [size: nSubBasins]

  double                   *_aMres;  ///< array used for storing reservoir masses [mg] or enthalpy [MJ] [size: nSubBasins]
  double              *_aMres_last;  ///< array storing reservoir mass [mg] or enthalpy [MJ] at start of timestep [size: nSubBasins]
  double                   *_aMsed;  ///< array storing reservoir sediment mass [mg] [size: nSubBasins]
  double             * _aMsed_last;  ///< array storing reservoir sediment mass [mg] at start of timestep [size: nSubBasins]
  double               *_aMout_res;  ///< array storing reservoir mass outflow [mg/d] or heat loss [MJ/d] [size: nSubBasins]
  double          *_aMout_res_last;  ///< array storing reservoir mass outflow [mg/d] or enthalpy outflow  [MJ/d] at start of timestep  [size: nSubBasins]
  double               *_aMresRain;  ///< array storing reservoir rain inputs [mg/d] or enthalpy input [MJ/d] [size: nSubBasins]

  // Mass balance tracking variables
  double         *_channel_storage;  ///< array storing channel storage [mg] or [MJ] [size: nSubBasins]
  double         *_rivulet_storage;  ///< array storing rivulet storage [mg] or [MJ] [size: nSubBasins]
  double              _cumul_input;  ///< cumulative mass [mg] or enthalpy [MJ] added to system
  double             _cumul_output;  ///< cumulative mass [mg] or enthalpy [MJ] lost from system
  double             _initial_mass;  ///< initial mass [mg] or enthalpy [MJ] in system

  // Source information
  constit_source       **_pSources;  ///< array of pointers to constituent sources [size: nSources]
  int                    _nSources;  ///< number of constituent sources

  int            **_aSourceIndices;  ///< lookup table to convert (global) water storage index to corresponding source, if any [size: nStateVariables][size: nHRUs]

  int              _nSpecFlowConcs;  ///< number of specified flow concentration/temperature time series [mg/L]
  CTimeSeries    **_pSpecFlowConcs;  ///< array of pointers to time series of specified flow concentration/temperatures - TS tag corresponds to SBID
  int              _nMassLoadingTS;  ///< number of specified mass/energy loading time series [kg/d]
  CTimeSeries    **_pMassLoadingTS;  ///< array of pointers to time series of mass loadings [kg/d] - TS tag corresponds to SBID

  ofstream                 _OUTPUT;  ///< output stream for Concentrations.csv/Temperatures.csv
  ofstream                 _POLLUT;  ///< output stream for Pollutograph.csv/StreamTemperatures.csv
  ofstream                _LOADING;  ///< output stream for MassLoadings.csv
  int                   _CONC_ncid;  ///< NetCDF id for Concentrations.nc/Temperatures.nc
  int                 _POLLUT_ncid;  ///< NetCDF id for Pollutograph.nc/StreamTemperatures.nc
  int                _LOADING_ncid;  ///< NetCDF id for MassLoadings.nc

  // private member funcctions
  void   DeleteRoutingVars();

  // mass balance routines
  double GetTotalRivuletConstituentStorage() const;
  double GetTotalChannelConstituentStorage() const;
  double GetMassAddedFromInflowSources(const double &t,const double &tstep) const;

  virtual double GetNetReachLosses        (const int p);
  virtual double GetCatchmentTransitLosses(const int p) const;

public:/*-------------------------------------------------------*/
  // Constructor/destructor
  CConstituentModel(CModel *pMod,CTransportModel *pTMod,string name,constit_type type,bool is_passive,const int c);
  ~CConstituentModel();

  // Accessors
  constit_type            GetType() const;
  string                  GetName() const;
  bool                    IsPassive() const;

  virtual double GetOutflowConcentration(const int p) const;
  double GetIntegratedMassOutflow(const int p,const double &tstep) const;

  bool   IsDirichlet             (const int i_stor,const int k,const time_struct &tt,double &Cs, const double blend=1.0) const;
  double GetSpecifiedMassFlux    (const int i_stor,const int k,const time_struct &tt) const;

  virtual double GetAdvectionCorrection(const CHydroUnit* pHRU,const int iFromWater,const int iToWater,const double& C) const;

  virtual double CalculateReportingConcentration(const double &M,const double &V) const;
  virtual double ConvertConcentration  (const double &Cs) const; //converts T->MJ/mm/m2 or C->mg/mm/m2 or Ciso->mg/mm/m2

  // Manipulators
  void   AddDirichletCompartment (const int i_stor,const int kk,const double Cs, const double Cs2);
  void   AddDirichletTimeSeries  (const int i_stor,const int kk,const CTimeSeries *pTS);
  void   AddInfluxSource         (const int i_stor,const int kk,const double flux);
  void   AddInfluxTimeSeries     (const int i_stor,const int kk,const CTimeSeries *pTS);
  void   AddInflowConcTimeSeries (const CTimeSeries *pTS);
  void   AddMassLoadingTimeSeries(const CTimeSeries *pTS);

  void   SetChannelMass          (const int p,const double mass);
  void   SetRivuletMass          (const int p,const double mass);
  void   SetMoutArray            (const int p,const int nsegs,   const double *aMout,const double MoutLast);
  void   SetMlatHist             (const int p,const int histsize,const double *aMlat,const double MlatLast);
  void   SetMinHist              (const int p,const int histsize,const double *aMin);
  void   SetInitialReservoirMass (const int p,const double res_mass,const double res_mass_last);
  void   SetInitReservoirSedMass (const int p, const double res_mass, const double res_mass_last);
  void   SetReservoirPrecipLoad  (const int p,const double precip_load_rate);
  void   SetReservoirMassOutflow (const int p,const double Mout,    const double MoutLast);

  void   ClearTimeSeriesData     (const optStruct& Options);

          void   Prepare                    (const optStruct &Options);
  virtual void   Initialize                 (const optStruct& Options);
          void   InitializeRoutingVars      ();
          void   CalculateInitialMassStorage(const optStruct& Options);

          void   IncrementCumulInput        (const optStruct &Options,const time_struct &tt);
          void   IncrementCumulOutput       (const optStruct &Options);

  virtual double ApplyInCatchmentRouting  (const int p,const double *aUnitHydro, const double *aQlatHist,const double *aMlatHist,const int nMlatHist,const double &tstep) const;
  virtual void   ApplyConvolutionRouting  (const int p,const double *aRouteHydro,const double *aQinHist, const double *aMinHist, const int nSegments,const int nMinHist,const double &tstep,double *aMout_new) const;
          void   ApplySpecifiedMassInflows(const int p,const double t,double &Minnew);
          void   SetMassInflows           (const int p,const double Minnew);
          void   SetLateralInfluxes       (const int p,const double RoutedMass);

  virtual void   InCatchmentRoute         (const int p,double &Mlat_new, const optStruct &Options);
  virtual void   PrepareForInCatchmentRouting(const int p);
  virtual void   PrepareForRouting        (const int p);
          void   RouteMass                (const int p,      double *aMoutnew,double &Mlat_new,double &ResMass,double &ResSedMass, const optStruct &Options,const time_struct &tt) const;
  virtual void   RouteMassInReservoir     (const int p,const double *aMoutnew,                 double &ResMass,double &ResSedMass, const optStruct &Options,const time_struct &tt) const;
          void   UpdateMassOutflows       (const int p,const double *aMoutnew,const double &Mlat_new,const double &ResMass,const double &ResSedMass, double &MassOutflow,
                                           const optStruct &Options,const time_struct &tt,const bool initialize);

  virtual void   WriteOutputFileHeaders      (const optStruct &Options);
  virtual void   WriteMinorOutput            (const optStruct &Options,const time_struct &tt);
  virtual void   WriteEnsimOutputFileHeaders (const optStruct &Options);
  virtual void   WriteEnsimMinorOutput       (const optStruct &Options,const time_struct &tt);
  virtual void   WriteNetCDFOutputFileHeaders(const optStruct &Options);
  virtual void   WriteNetCDFMinorOutput      (const optStruct &Options,const time_struct& tt);
          void   WriteMajorOutput            (ofstream& RVC) const;
  virtual void   CloseOutputFiles            ();
};

#endif
