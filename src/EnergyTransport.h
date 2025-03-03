/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2023 the Raven Development Team
----------------------------------------------------------------*/

#ifndef ENERGY_TRANSPORT_H
#define ENERGY_TRANSPORT_H

#include "RavenInclude.h"
#include "Transport.h"

///////////////////////////////////////////////////////////////////
/// \brief Class for coordinating transport simulation for enthalpy (child of ConstituentModel)
/// \details Implemented in EnergyTransport.cpp
//
class CEnthalpyModel :public CConstituentModel
{
  double **_aEnthalpySource;   ///< recent time history of reach source term [MJ/m3/d] [size: nSubBasins x nMinhist(p)]
  double    *_aEnthalpyBeta;   ///< array of beta terms for reach energy exchange [1/d] [size: nSubBasins]

  double **_aEnthalpySource2;  ///< recent time history of in-catchment source term [MJ/m3/d] [size: nSubBasins x nMlathist(p)]
  double        *_aInCatch_a;  ///< array of a terms for in-catchment energy exchange [1/d] [size: nSubBasins]
  double        *_aInCatch_b;  ///< array of b terms for in-catchment energy exchange [1/d] [size: nSubBasins]
  double             *_aKbed;  ///< bed heat transfer coefficient [MJ/m2/d/K] [size: nSubBasins]

  int            _mHyporheic;  ///< model soil layer corresponding to hyporheic mixing layer (default=2)
  double       *_aMinResTime;  ///< minimum residence time [d] (<1 dt) in each stream reach [size: nSubBasins]

  double          *_aBedTemp;  ///< array of riverbed temperatures (*state variable*) [C] [size: nSubBasins]

  double       *_aTave_reach;  ///< time-lagged average reach water temperature [C] [size: nSubBasins]
  double   *_aSS_temperature;  ///< array of steady state target temperatures in each basin [C] [size: nSubBasins]

  ofstream        _STREAMOUT;  ///< Output stream for StreamReachEnergyBalances.csv
  ofstream          _LAKEOUT;  ///< Output stream for LakeEnergyBalances.csv

  bool      _anyGaugedLakes;   ///< true if any gauged lakes/reservoirs exist for writing LakeEnergyBalances.csv

  double GetReachFrictionHeat(const double &Q,const double &slope,const double &top_width) const;

  void   UpdateReachEnergySourceTerms    (const int p);
  void   UpdateCatchmentEnergySourceTerms(const int p);
  double GetNetReachLosses               (const int p);
  double GetCatchmentTransitLosses       (const int p) const;

  double FunkyTemperatureIntegral(const int k,const int m,const double *S,
                                  const double &h1,const double &h2,const double &h3,const double &beta,const double &dt) const;

public:/*-------------------------------------------------------*/
  CEnthalpyModel(CModel *pMod,CTransportModel *pTMod,string name,const int c);
  ~CEnthalpyModel();

  double CalculateReportingConcentration(const double &M,const double &V) const;
  double GetOutflowConcentration        (const int p) const;
  double GetBedTemperature              (const int p) const;
  double ConvertConcentration           (const double &Cs) const;

  //Accessors specific to Enthalpy Transport
  double GetIceContent           (const double *state_vars, const int iWater) const;
  double GetWaterTemperature     (const double *state_vars, const int iWater) const;

  double GetEnergyLossesInTransit(const int p,double &Q_sens,double &Q_GW) const;
  double GetEnergyLossesFromReach(const int p,double &Q_sens,double &Q_cond,double &Q_lat,double &Q_GW,double &Q_rad_in,double &Q_lw_in, double &Q_lw_out,double &Q_lateral, double &Q_fric, double &Tave) const;
  double GetEnergyLossesFromLake (const int p,double &Q_sens,double &Q_cond,double &Q_lat,double &Q_rad_in,double &Q_lw_in, double &Q_lw_out, double &Q_rain) const;

  double GetOutflowIceFraction   (const int p) const;
  double GetAvgLatentHeatFlux    () const;
  double GetDirichletEnthalpy    (const CHydroUnit *pHRU, const double &Cs) const;

  //Manipulators
  //
  void   SetBedTemperature       (const int p, const double &T);
  void   SetHyporheicLayer       (const int m);

  //Manipulators (inherited from CConstitModel)
  void   Initialize              (const optStruct& Options);

  void   PrepareForInCatchmentRouting(const int p);
  void   PrepareForRouting       (const int p);
  double ApplyInCatchmentRouting (const int p,
                                  const double *aUnitHydro,
                                  const double *aQlatHist,
                                  const double *aMlatHist,
                                  const int     nMlatHist,
                                  const double &tstep) const;
  void   ApplyConvolutionRouting (const int     p,
                                  const double *aRouteHydro,
                                  const double *aQinHist,
                                  const double *aMinHist,
                                  const int     nSegments,
                                  const int     nMinHist,
                                  const double &tstep,double *aMout_new) const;
  void   RouteMassInReservoir    (const int    p ,
                                  const double* aMoutnew, double& ResMass, double& ResSedMass,
                                  const optStruct   &Options,
                                  const time_struct &tt) const;

  void   WriteOutputFileHeaders      (const optStruct& Options);
  void   WriteMinorOutput            (const optStruct& Options,const time_struct& tt);
  void   WriteEnsimOutputFileHeaders (const optStruct &Options);
  void   WriteEnsimMinorOutput       (const optStruct &Options,const time_struct &tt);
  void   CloseOutputFiles            ();
};
#endif
