/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2020 the Raven Development Team
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
struct constituent
{
  string   name;          ///< constituent name (e.g., "Nitrogen")

  constit_type type;      ///< AQUEOUS [mg], ENTHALPY [MJ], or TRACER [-]
  bool     can_evaporate; ///< true if constituent can be transported through evaporation
  bool     is_passive;    ///< doesn't transport via advection

  double   decay_rate;    ///< independent linear decay rate of constituent (happens everywhere; not environmentally mediated) [1/d]

  double   cumul_input;   ///< cumulative mass [mg] or enthalpy [MJ] lost from system  
  double   cumul_output;  ///< cumulative mass [mg] or enthalpy [MJ] lost from system 
  double   initial_mass;  ///< initial mass [mg] or enthalpy [MJ] in system 

  ofstream OUTPUT;        ///< for concentrations.csv
  ofstream POLLUT;        ///< for pollutograph.csv

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
class CTransportModel
{
private:/*------------------------------------------------------*/
  CModel *pModel;                    ///< pointer to model object
  static CTransportModel *_pTransModel;    ///< pointer to only transport model (needed for access to transport model through static functions)

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

  int                _nConstituents;  ///< number of transported constituents [c=0.._nConstituents-1]
  constituent      **_pConstituents;  ///< array of pointers to constituent structures [size: _nConstituents]
  transport_params **_pConstitParams; ///< array of pointers to constituent parameters [size: _nConstituents]

  double ***_aMinHist;                ///< array used for storing routing upstream loading history [mg/d] or [MJ/d] [size: nSubBasins x _nConstituents x nMinhist(p)]
  double ***_aMlatHist;               ///< array used for storing routing lateral loading history [mg/d] or [MJ/d] [size: nSubBasins x _nConstituents x nMlathist(p)]
  double ***_aMout;                   ///< array storing current mass flow at points along channel [mg/d] or [MJ/d] [size: nSubBasins x _nConstituents x _nSegments(p)]
  double  **_aMout_last;              ///< array used for storing mass outflow from channel at start of timestep [mg/d] or [MJ/d] [size: nSubBasins x _nConstituents] 

  double  **_aMres;                   ///< array used for storing reservoir masses [mg] or enthalpy [MJ] [size: nSubBasins x _nConstituents]
  double  **_aMres_last;              ///< array storing reservoir mass [mg] or enthalpy [MJ] at start of timestep [size: nSubBasins x _nConstituents]
  double  **_aMout_res;               ///< array storing reservoir mass outflow [mg/d] or heat loss [MJ/d] [size: nSubBasins x _nConstituents]
  double  **_aMout_res_last;          ///< array storing reservoir mass outflow [mg/d] or enthalpy outflow  [MJ/d] at start of timestep  [size: nSubBasins x _nConstituents]

  double  **_aMlat_last;              ///< array storing mass outflow from start of timestep [size: nSubBasins x _nConstituents]

  double  **_channel_storage;         ///< array storing channel storage [mg] or [MJ] [size: nSubBasins x _nConstituents] 
  double  **_rivulet_storage;         ///< array storing rivulet storage [mg] or [MJ] [size: nSubBasins x _nConstituents] 
  
  constit_source **pSources;         ///< array of pointers to constituent sources [size: nSources]
  int              nSources;         ///< number of constituent sources
  int            **_aSourceIndices;  ///< lookup table to convert constitutent and (global) water storage index to corresponding source, if any [size: _nConstituents x nStateVariables]
  
  void m_to_cj(const int layerindex, int &c, int &j) const;
  void DeleteRoutingVars();
  void InitializeConstitParams(transport_params *P);

  double GetTotalRivuletConstituentStorage(const int c) const;
  double GetTotalChannelConstituentStorage(const int c) const;

  //Specific to Enthalpy Transport ------------------------------- 
  bool    _EnthalpyIsSimulated;      ///< boolean; true if ENTHALPY is one of the constituents in the model
  double **_aEnthalpySource;         ///< recent time history of reach source term [MJ/m3/d] [size: nSubBasins x nMinhist(p)]
  double  *_aEnthalpyBeta;           ///< array of beta terms for reach energy exchange [1/d] [size: nSubBasins] 

  double GetAvgLatentHeatFlux() const;

public:/*-------------------------------------------------------*/
  CTransportModel(CModel *pMod);
  ~CTransportModel();

  //Accessors
  static int    GetLayerIndexFromName(const string name,const int comp_m);
  int          GetLayerIndexFromName2(const string name,const int comp_m) const;     //non-static version

  static string GetConstituentTypeName (const int m);         //e.g., "Nitrogen" (m=layer index)
  static string GetConstituentTypeName2(const int c);         //e.g., !Nitrogen (c=constituent index)
  static string GetConstituentLongName (const int layerindex);//e.g., "Nitrogen in Soil Water[2]"
  string GetConstituentName     (const int layerindex) const; //e.g., "Nitrogen in Soil Water[2]"
  string GetConstituentShortName(const int layerindex) const; //e.g., "!Nitrogen_SOIL[2]"

  int    GetNumConstituents() const;
  const constituent      *GetConstituent(const int c) const;
  const transport_params *GetConstituentParams(const int c) const;
  int    GetConstituentIndex(const string name) const;

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
  
  int    GetStorWaterIndex         (const int ii) const;
  int    GetWaterStorIndexFromLayer(const int m) const;

  int    GetLayerIndex             (const int c, const int i_stor) const;

  double GetIceContent           (const double *state_vars, const int iWater) const;
  double GetWaterTemperature     (const double *state_vars, const int iWater) const;

  double GetOutflowConcentration (const int p, const int c) const;
  double GetIntegratedMassOutflow(const int p, const int c,const double &tstep) const;
  double GetOutflowIceFraction   (const int p, const int c) const;

  bool   ConstituentIsPassive   (const int c) const;

  double GetDecayCoefficient    (const int c,const CHydroUnit *pHRU, const int iStorWater) const;
  double GetRetardationFactor   (const int c,const CHydroUnit *pHRU, const int iFromWater,const int iToWater) const;
  double GetTransformCoefficient(const int c, const int c2, const CHydroUnit *pHRU, const int iStorWater) const;
  double GetStoichioCoefficient (const int c, const int c2, const CHydroUnit *pHRU, const int iStorWater) const;

  bool   IsDirichlet         (const int i_stor, const int c, const int k, const time_struct &tt, double &Cs) const;
  double GetSpecifiedMassFlux(const int i_stor, const int c, const int k, const time_struct &tt) const;

  //Manipulators
  void   AddConstituent         (string name, constit_type type, bool is_passive); 
  void   AddDirichletCompartment(const string const_name, const int i_stor, const int kk, const double Cs);
  void   AddDirichletTimeSeries (const string const_name, const int i_stor, const int kk, const CTimeSeries *pTS);
  void   AddInfluxSource        (const string const_name, const int i_stor, const int kk, const double flux);
  void   AddInfluxTimeSeries    (const string const_name, const int i_stor, const int kk, const CTimeSeries *pTS);
  //
  void   Prepare(const optStruct &Options);
  void   CalculateLateralConnections();
  void   Initialize();
  void   InitializeRoutingVars();

  void   IncrementCumulInput (const optStruct &Options, const time_struct &tt);
  void   IncrementCumulOutput(const optStruct &Options);

  void   SetGlobalParameter(const string const_name, const string param_name, const double &value, bool noisy);
  void   SetMassInflows    (const int p, const double *aMinnew);
  void   SetLateralInfluxes(const int p, const double *aRoutedMass);
  void   RouteMass         (const int p,       double **aMoutnew, double *aResMass, const optStruct &Options,const time_struct &tt) const;
  void   UpdateMassOutflows(const int p,       double **aMoutnew, double *aResMass, double *aResMassOutflow, const optStruct &Options,const time_struct &tt,bool initialize);
  
  void   UpdateReachEnergySourceTerms(const int p);

  void   WriteOutputFileHeaders     (const optStruct &Options) const;
  void   WriteMinorOutput           (const optStruct &Options, const time_struct &tt) const;
  void   WriteEnsimOutputFileHeaders(const optStruct &Options) const;
  void   WriteEnsimMinorOutput      (const optStruct &Options, const time_struct &tt) const;
  void   CloseOutputFiles           () const;
};
#endif
