/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2021 the Raven Development Team
----------------------------------------------------------------*/

#ifndef ISOTOPE_TRANSPORT_H
#define ISOTOPE_TRANSPORT_H

#include "RavenInclude.h"
#include "Transport.h"

/// Isotope global constants
const double TO_PER_MILLE=1000.0;   //convert from raw ratio to per mille (per thousand)
const double RV_VMOW     =0.002228; // Vienna Mean Ocean Water Mass Ratio  (O18/O16)  

/// Isotope types
enum iso_type 
{
  ISO_O18, //Oxygen-18
  ISO_H2   //Hydrogen-2 (Deuterium)
};

///////////////////////////////////////////////////////////////////
/// \brief Class for coordinating transport simulation for isotopes (child of ConstituentModel)
/// \details Implemented in IsotopeTransport.cpp
//
class CIsotopeModel :public CConstituentModel
{
private:
  iso_type _isotope;

  double RvalToConcentration(const double &Rv  ) const;
  double ConcentrationToRval(const double &conc) const;
  double ConcToComposition  (const double &conc) const;
  double CompositionToConc  (const double &d   ) const;

public:/*-------------------------------------------------------*/
  CIsotopeModel(CModel *pMod,CTransportModel *pTMod,string name,const int c, const iso_type iso);
  ~CIsotopeModel();

  double GetConcentration(const double &mass, const double &vol); //returns [o/oo]

  double GetAdvectionCorrection(const CHydroUnit* pHRU,const int iFromWater,const int iToWater,const double &C) const;

  //Manipulators (inherited from CConstitModel)
  void   Initialize                 (const optStruct& Options);
  void   WriteEnsimOutputFileHeaders(const optStruct &Options);
  void   WriteEnsimMinorOutput      (const optStruct &Options,const time_struct &tt);
};
#endif