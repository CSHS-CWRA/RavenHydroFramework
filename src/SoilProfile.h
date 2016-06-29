/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright © 2008-2014 the Raven Development Team
------------------------------------------------------------------
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
    int          nHorizons;			 ///< Number of soil horizons in stack
    CSoilClass **pSoilClasses;	 ///< Array of soil horizones making up profile
    double      *thicknesses;    ///< Array of Thickness of horizons [m]

    static CSoilProfile **pAllSoilProfiles; ///< Reference to array of all soil profiles in model
    static int            NumSoilProfiles;	///< Number of soil profiles in model (size of pAllSoilProfiles)

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

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for a stack of aquifers/aquitards
/// \details Each layer is defined
///   by a soil class (e.g., sand, clay, guelph loam) and thickness
///   horizon m=0 is top horizon, m=nHorizons-1 is bottom. Each aquifer 
///   layer is underlain by a less conductive aquitard with specified thickness.
//

class CAquiferStack
{
  private:/*-------------------------------------------------------*/  
    string       tag; 
    int          nLayers;

    CSoilClass **pAquiferSoils;	  ///< array of pointers to aquifer soils
    CSoilClass **pAquitardSoils;	///< array of pointers to aquitard soils
    double      *aAquifer_thick;  ///< array of aquitfer thicknesses [m] 
    double      *aAquitard_thick; ///< array of aquitard thicknesses [m]

    static CAquiferStack **pAllAqStacks; ///< Reference to array of all aquifer stacks in model
    static int             NumAqStacks;	///< Number of aquifer stacks in model (size of pAllAquiferStacks)

  public:/*-------------------------------------------------------*/
    //Constructors:
    CAquiferStack(const string name);
    ~CAquiferStack();

    //Accessors
    string             GetTag               ()            const;
    double             GetNumLayers         ()            const;

    double             GetAquiferThickness  (const int m) const;
    const soil_struct *GetAquiferSoil       (const int m) const;
    string             GetAquiferSoilTag    (const int m) const;
    

    double             GetAquitardThickness (const int m) const;
    const soil_struct *GetAquitardSoil      (const int m) const;
    string             GetAquitardSoilTag   (const int m) const;
    
    void               AddLayer             (double            thickness, //[m]
                                             const CSoilClass *pLayerSoil,
                                             double            aquitard_thick,
                                             const CSoilClass *pAquitardSoil);
    void               AddLayer             (double            thickness, //[m]
                                             const CSoilClass *pLayerSoil);

    void               AllocateAqLayers     (const int           nAqLayers,
                                             const soil_struct **pSoil,
                                             double             *thickness) const; 

    static int                  GetNumAqStacks         ();
    static const CAquiferStack *StringToAqStack   (const string s);

    static void                 DestroyAllAqStacks     ();

};
#endif