/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ----------------------------------------------------------------
  Class CSoilProfile
  Class CAquiferStack
  ----------------------------------------------------------------*/
#ifndef SOIL_PROFILE_H
#define SOIL_PROFILE_H

#include "SoilAndLandClasses.h"

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for a stack of soil horizons
/// \details Each horizon is defined
///   by a soil class (e.g., sand, clay, guelph loam) and thickness
///   horizon m=0 is top horizon, m=nHorizons-1 is bottom
//
class CSoilProfile
{
protected:/*-------------------------------------------------------*/
  string       tag;            ///< nickname for profile
  int          nHorizons;      ///< Number of soil horizons in stack
  CSoilClass **pSoilClasses;   ///< Array of soil horizones making up profile
  double      *thicknesses;    ///< Array of Thickness of horizons [m]

  static CSoilProfile **pAllSoilProfiles; ///< Reference to array of all soil profiles in model
  static int            NumSoilProfiles;  ///< Number of soil profiles in model (size of pAllSoilProfiles)

public:/*-------------------------------------------------------*/
  //Constructors:
  CSoilProfile(const string name);
  ~CSoilProfile();

  //Accessors
  string             GetTag        ()            const;
  double             GetThickness  (const int m) const;
  const soil_struct *GetSoilStruct (const int m) const;
  string             GetSoilTag    (const int m) const;
  double             GetNumHorizons()            const;

  void           AllocateSoilLayers(const int           nSoilLayers,
                                    const soil_struct **pSoil,
                                    double             *thickness) const;

  void               AddHorizon    (double            thickness, //[m]
                                    const CSoilClass *pHorizonClass);

  static int                 GetNumProfiles        ();
  static const CSoilProfile *StringToSoilProfile   (const string s);

  static void                DestroyAllSoilProfiles();

  static void                SummarizeToScreen     ();
};

#endif
