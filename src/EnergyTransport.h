/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2021 the Raven Development Team
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

  ofstream       _STREAMOUT;   ///< Output stream for StreamReachEnergyBalances.csv

  double GetReachFrictionHeat(const double &Q,const double &slope,const double &perim) const;
  void   UpdateReachEnergySourceTerms(const int p);
  double GetNetReachLosses           (const int p) const;

public:/*-------------------------------------------------------*/
  CEnthalpyModel(CModel *pMod,CTransportModel *pTMod,string name,const int c);
  ~CEnthalpyModel();

  double CalculateConcentration(const double &M,const double &V) const;
  double GetOutflowConcentration(const int p) const;

  //Accessors specific to Enthalpy Transport
  double GetIceContent           (const double *state_vars, const int iWater) const;
  double GetWaterTemperature     (const double *state_vars, const int iWater) const;

  double GetEnergyLossesFromReach(const int p,double &Q_sens,double &Q_lat,double &Q_GW,double &Q_rad,double &Q_fric) const;
  double GetOutflowIceFraction   (const int p) const;
  double GetAvgLatentHeatFlux    () const;
  double GetDirichletEnthalpy    (const CHydroUnit *pHRU, const double &Cs) const;

  //Manipulators (inherited from CConstitModel)
  void   Initialize              (const optStruct& Options);
  void   ApplyConvolutionRouting (const int     p,
                                  const double *aRouteHydro,
                                  const double *aQinHist,
                                  const double *aMinHist,
                                  const int     nSegments,
                                  const int     nMinHist,
                                  const double &tstep,double *aMout_new) const;

  void   UpdateMassOutflows(const int          p,
                                  double      *aMoutnew,
                                  double      &ResMass,
                                  double      &MassOutflow,
                            const optStruct   &Options,
                            const time_struct &tt,
                                  bool         initialize);

  void   WriteOutputFileHeaders      (const optStruct& Options);
  void   WriteMinorOutput            (const optStruct& Options,const time_struct& tt);
  void   WriteEnsimOutputFileHeaders (const optStruct &Options);
  void   WriteEnsimMinorOutput       (const optStruct &Options,const time_struct &tt);
  void   CloseOutputFiles            ();
};
#endif