/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2021 the Raven Development Team
----------------------------------------------------------------*/

#ifndef TRANSPORTMODEL_H
#define TRANSPORTMODEL_H

#include "RavenInclude.h"
#include "Model.h"

enum constit_type {
  AQUEOUS,
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
  double flux;            ///< fixed mass flux [mg/m2/d] or heat flux [MJ/m2/d] (or DOESNT_EXIST=-1 if time series should be used)
  const  CTimeSeries *pTS; ///< time series of fixed concentration or mass/heat flux (or NULL if fixed should be used)
};

struct transport_params
{
  double decay_coeff;     ///< constituent base linear decay coefficient [1/d]
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
  int *_iFromWater;                  ///< state variable indices of water compartment source [size: nAdvConnections]
  int *_iToWater;                    ///< state variable indices of water compartment destination [size: nAdvConnections]
  int *_js_indices;                  ///< process index (j*) of connection [size: nAdvConnections]

  int  _nLatConnections;             ///< number of lateral advective transport connections between water storage units
  int  *_iLatFromWater;              ///< state variable indices of water compartment source [size: _nLatConnections]
  int  *_iLatToWater;                ///< state variable indices of water compartment destination [size: _nLatConnections]
  int  *_iLatFromHRU;                ///< state variable indices of water compartment source [size: _nLatConnections]
  int  *_iLatToHRU;                  ///< state variable indices of water compartment destination [size: _nLatConnections]
  int  *_latqss_indices;             ///< process index (q**) of connection [size: nLatConnections]

  int  _nWaterCompartments;          ///< number of water storage compartments which may contain constituent
  int *_iWaterStorage;               ///< state variable indices of water storage compartments which may contain constituent [size: _nWaterCompartments]
  int *_aIndexMapping;               ///< lookup table to convert state variable index i to local water storage index [size: pModel::_nStateVars prior to transport variables being included]
                                     ///< basically inverse of iWaterStorage
  int  _nIndexMapping;               ///< size of index mapping array [=pModel::_nStateVars prior to transport variables being included]

  int                 _nConstituents;  ///< number of transported constituents [c=0.._nConstituents-1]
  CConstituentModel **_pConstitModels; ///< array of pointers to constituent models [size:_nConstituents]
  CEnthalpyModel     *_pEnthalpyModel; ///< special pointer to enthalpy constituent, if present (or NULL if not)

  void m_to_cj(const int layerindex,int &c,int &j) const;

  string GetConstituentLongName_loc(const int layerindex) const; //e.g., "Nitrogen in Soil Water[2]"

  static CTransportModel *_pTransModel;    ///< pointer to self (needed for access to transport model through static functions)

public:/*-------------------------------------------------------*/
  CTransportModel(CModel *pMod);
  ~CTransportModel();

  //Accessors
  static int    GetLayerIndexFromName(const string name,const int comp_m);
  int          GetLayerIndexFromName2(const string name,const int comp_m) const;     //non-static version

  static string GetConstituentTypeName(const int layerindex);//e.g., "Nitrogen" (m=layer index)
  static string GetConstituentTypeName2(const int c);         //e.g., "Nitrogen" (c=constituent index)
  static string GetConstituentLongName(const int layerindex);//e.g., "Nitrogen in Soil Water[2]"
  string GetConstituentShortName(const int layerindex) const; //e.g., "!Nitrogen_SOIL[2]"

  int    GetNumConstituents() const;
  const transport_params *GetConstituentParams(const int c) const;
  int    GetConstituentIndex(const string name) const;

  CConstituentModel *GetConstituentModel(const int c);
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

  int    GetLayerIndex(const int c,const int i_stor) const;

  const CEnthalpyModel *GetEnthalpyModel() const;

  double GetConcentration           (const int k,const int sv_index) const;

  double GetAdvectionCorrection     (const int c,const CHydroUnit* pHRU,const int iFromWater,const int iToWater,const double& C) const;


  //Manipulators
  void   AddConstituent             (string name,constit_type type,bool is_passive);
  //
  void   Prepare                    (const optStruct &Options);
  void   Initialize                 (const optStruct& Options);
  void   CalculateLateralConnections();

  void   IncrementCumulInput        (const optStruct &Options,const time_struct &tt);
  void   IncrementCumulOutput       (const optStruct &Options);

  void   SetGlobalParameter(const string const_name,const string param_name,const double &value,bool noisy);

  void   WriteOutputFileHeaders     (const optStruct &Options) const;
  void   WriteMinorOutput           (const optStruct &Options,const time_struct &tt) const;
  void   WriteMajorOutput           (ofstream &RVC) const;
  void   WriteEnsimOutputFileHeaders(const optStruct &Options) const;
  void   WriteEnsimMinorOutput      (const optStruct &Options,const time_struct &tt) const;
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
  int               _constit_index;  ///< master constituent index of this constituent 
  constit_type               _type;  ///< AQUEOUS [mg], ENTHALPY [MJ], or TRACER [-]
  bool                 _is_passive;  ///< doesn't transport via advection (default: false)

  transport_params *_pConstitParams; ///< pointer to constituent parameters 

  // Routing/state var storage
  double               **_aMinHist;  ///< array used for storing routing upstream loading history [mg/d] or [MJ/d] [size: nSubBasins x nMinhist(p)]
  double              **_aMlatHist;  ///< array used for storing routing lateral loading history [mg/d] or [MJ/d] [size: nSubBasins  x nMlathist(p)]
  double                  **_aMout;  ///< array storing current mass flow at points along channel [mg/d] or [MJ/d] [size: nSubBasins x _nSegments(p)]
  double              *_aMout_last;  ///< array used for storing mass outflow from channel at start of timestep [mg/d] or [MJ/d] [size: nSubBasins ] 
  double              *_aMlat_last;  ///< array storing mass/energy outflow from start of timestep [size: nSubBasins]

  double                   *_aMres;  ///< array used for storing reservoir masses [mg] or enthalpy [MJ] [size: nSubBasins]
  double              *_aMres_last;  ///< array storing reservoir mass [mg] or enthalpy [MJ] at start of timestep [size: nSubBasins]
  double               *_aMout_res;  ///< array storing reservoir mass outflow [mg/d] or heat loss [MJ/d] [size: nSubBasins]
  double          *_aMout_res_last;  ///< array storing reservoir mass outflow [mg/d] or enthalpy outflow  [MJ/d] at start of timestep  [size: nSubBasins]

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

  // private member funcctions
  void   DeleteRoutingVars();
  void   InitializeConstitParams(transport_params *P);

  // mass balance routines
  double GetTotalRivuletConstituentStorage() const;
  double GetTotalChannelConstituentStorage() const;
  double GetMassAddedFromInflowSources(const double &t,const double &tstep) const;
  
  virtual double GetNetReachLosses(const int p) const;

public:/*-------------------------------------------------------*/
  // Constructor/destructor
  CConstituentModel(CModel *pMod,CTransportModel *pTMod,string name,constit_type type,bool is_passive,const int c);
  ~CConstituentModel();

  // Accessors
  const transport_params *GetConstituentParams() const;
  constit_type            GetType() const;
  string                  GetName() const;
  bool                    IsPassive() const;
  
  virtual double GetOutflowConcentration(const int p) const;
  double GetIntegratedMassOutflow(const int p,const double &tstep) const;

  bool   IsDirichlet         (const int i_stor,const int k,const time_struct &tt,double &Cs) const;
  double GetSpecifiedMassFlux(const int i_stor,const int k,const time_struct &tt) const;

  double GetDecayCoefficient (const CHydroUnit *pHRU,const int iStorWater) const;
  virtual double GetAdvectionCorrection(const CHydroUnit* pHRU,const int iFromWater,const int iToWater,const double& C) const;

  virtual double CalculateConcentration(const double &M,const double &V) const;

  // Manipulators
  void   AddDirichletCompartment (const int i_stor,const int kk,const double Cs);
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
  void   SetReservoirMassOutflow (const int p,const double Mout,    const double MoutLast);
  
          void   Prepare              (const optStruct &Options);
  virtual void   Initialize           (const optStruct& Options);
          void   InitializeRoutingVars();

          void   IncrementCumulInput  (const optStruct &Options,const time_struct &tt);
          void   IncrementCumulOutput (const optStruct &Options);

  virtual void   ApplyConvolutionRouting  (const int p,const double *aRouteHydro,const double *aQinHist,const double *aMinHist,const int nSegments,const int nMinHist,const double &tstep,double *aMout_new) const;
          void   ApplySpecifiedMassInflows(const int p,const double t,double &Minnew);
          void   SetMassInflows           (const int p,const double Minnew);
          void   SetLateralInfluxes       (const int p,const double RoutedMass);
          void   RouteMass                (const int p,double *aMoutnew,double &ResMass,const optStruct &Options,const time_struct &tt) const;
  virtual void   UpdateMassOutflows       (const int p,double *aMoutnew,double &ResMass,double &ResMassOutflow,const optStruct &Options,const time_struct &tt,bool initialize);

  virtual void   WriteOutputFileHeaders     (const optStruct &Options);
  virtual void   WriteMinorOutput           (const optStruct &Options,const time_struct &tt);
          void   WriteMajorOutput           (ofstream &RVC) const;
  virtual void   WriteEnsimOutputFileHeaders(const optStruct &Options);
  virtual void   WriteEnsimMinorOutput      (const optStruct &Options,const time_struct &tt);
  virtual void   CloseOutputFiles();
};

///////////////////////////////////////////////////////////////////
/// \brief Class for coordinating transport simulation for nutrients (child of ConstituentModel)
/// \details Implemented in NutrientTransport.cpp
//
class CNutrientModel:public CConstituentModel
{
public:/*-------------------------------------------------------*/
  //Accessors specific to Nutrient/Contaminant Transport
  double GetDecayCoefficient   (const CHydroUnit *pHRU,const int iStorWater) const;
  double GetAdvectionCorrection(const CHydroUnit *pHRU,const int iFromWater,const int iToWater, const double &C) const;

  //double GetTransformCoefficient(const int c,const int c2,const CHydroUnit *pHRU,const int iStorWater) const;
  //double GetStoichioCoefficient(const int c,const int c2,const CHydroUnit *pHRU,const int iStorWater) const;
};
#endif
