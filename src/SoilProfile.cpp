/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2022 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Properties.h"
#include "SoilAndLandClasses.h"
#include "SoilProfile.h"
/*****************************************************************
   Constructor / Destructor
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the soil profile constructor
/// \param name [in] Soil profile identifier
CSoilProfile::CSoilProfile(const string name)
{
  tag=name;
  if (!DynArrayAppend((void**&)(pAllSoilProfiles),(void*)(this),NumSoilProfiles)){
    ExitGracefully("CSoilProfile::Constructor: creating NULL soil profile",BAD_DATA);};

  nHorizons   =0;
  pSoilClasses=NULL;
  thicknesses =NULL;
}

///////////////////////////////////////////////////////////////////
/// \brief Implentation of the destructor
//
CSoilProfile::~CSoilProfile()
{
  if (DESTRUCTOR_DEBUG){cout<<"  DELETING SOIL PROFILE "<<endl;}
  delete [] pSoilClasses;
  delete [] thicknesses;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns tag of soil profile (e.g. "ALL_SILT")
/// \return Tag of soil profile
//
string CSoilProfile::GetTag() const{return tag;}

///////////////////////////////////////////////////////////////////
/// \brief Returns pointer to soil horizon in profile
/// \param m [in] Index of soil horizon to be returned
/// \return Soil horizon corresponding to index m
//
const soil_struct *CSoilProfile::GetSoilStruct(const int m) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((m<0) || (m>=nHorizons),
                   "CSoilProfile::GetSoilStruct: bad soil horizon index",BAD_DATA);
#endif
  return pSoilClasses[m]->GetSoilStruct();
}

///////////////////////////////////////////////////////////////////
/// \brief Returns thickness of soil horizon in meters [m]
/// \param m [in] Index of a specific soil horizon
/// \return Thickness of soil horizon [m]
//
double CSoilProfile::GetThickness (const int m) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((m<0) || (m>=nHorizons),
                   "CSoilProfile::GetSoilStruct: bad soil horizon index",BAD_DATA);
#endif
  return thicknesses[m];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns tag of soil horizon specified by m
/// \param m [in] Soil horizon specifier
/// \return soil horizon specified by m
//
string CSoilProfile::GetSoilTag (const int m) const
{
  ExitGracefullyIf((m<0) || (m>=nHorizons),
                   "CSoilProfile::GetSoilStruct: bad soil horizon index",BAD_DATA);
  return pSoilClasses[m]->GetTag();
}

///////////////////////////////////////////////////////////////////
/// \brief Returns the number of horizons in soil profile
/// \return Number of horizons in soil profile
double CSoilProfile::GetNumHorizons () const
{
  return nHorizons;
}

//////////////////////////////////////////////////////////////////
/// \brief Adds horizon of specified soil class to bottom of profile
/// \param thickness [in] Thickness of horizon to be added
/// \param *pHorizonClass [in] pointer to horizon object to be added to collection
//
void CSoilProfile::AddHorizon(double thickness, //[m]
                              const CSoilClass *pHorizonClass)
{
  if (!DynArrayAppend((void**&)(pSoilClasses),(void*)(pHorizonClass),nHorizons)){
    ExitGracefully("CSoilProfile::Constructor: creating NULL soil profile",BAD_DATA);};

  double *tmpthick=new double [nHorizons];
  for (int m=0;m<nHorizons-1;m++){//copy previous array
    tmpthick[m]=thicknesses[m];
  }
  tmpthick[nHorizons-1]=thickness;//add new data
  delete [] thicknesses;//delete previous array
  thicknesses=tmpthick;
}

//////////////////////////////////////////////////////////////////
/// \brief Allocate soil layers
//
/// \details Discretizes soil profile into nSoilLayers layers for numerical
///     solution, e.g., takes a 7 meter soil profile with three horizons
///     and breaks it up into 12 layers.
///     The algorithm shoots for at least one layer per profile, then
///     aims for higher discretization in top layers, coarser
///     discretization at depth. Thicknesses are in metres
//
/// \param nSoilLayers [in] Number of soil layers in profile
/// \param **pSoil [out] Reference to soil properties
/// \param *thickness [in & out] Reference to array of layer thicknesses [m]
//
void CSoilProfile::AllocateSoilLayers(const int           nSoilLayers,
                                      const soil_struct  **pSoil,
                                      double              *thickness) const
{
  static soil_struct blank_soil_struct;

  CSoilClass::InitializeSoilProperties(blank_soil_struct, false,0);
  blank_soil_struct.porosity =0;
  blank_soil_struct.hydraul_cond=0;
  //A problem if nSoilLayers<nHorizons!!-average layers!?
  int m,divi;
  if (nHorizons==0){ //special case for lakes and glaciers
    for (m=0;m<nSoilLayers;m++){
      pSoil[m]=&blank_soil_struct;
      thickness[m]=REAL_SMALL;
    }
    return;
  }
  else if (nSoilLayers==nHorizons){//e.g., 3 soil layers and 3 horizons
    for (m=0;m<nHorizons;m++){
      pSoil    [m]=pSoilClasses[m]->GetSoilStruct();
      thickness[m]=thicknesses [m];
    }
    return;
  }
  else if ((nSoilLayers%nHorizons)==0)//e.g., 9 soil layers and 3 horizons
  {
    divi=(nSoilLayers/nHorizons);
    for (m=0;m<nHorizons;m++){
      for (int i=0;i<divi;i++){
        pSoil    [m*divi+i]=pSoilClasses[m]->GetSoilStruct();
        thickness[m*divi+i]=thicknesses [m]/divi;
      }
    }
    return;
  }
  else
  {
    cout <<"nSoilLayers:"<<nSoilLayers<<" nHorizons:"<<nHorizons<<endl;
    ExitGracefully("CSoilProfile::AllocateSoilLayers: can currently handle a number of soil profile layers which are an even multiple of the number of soil horizons (usually 1:1)",BAD_DATA);
  }
}
/*****************************************************************
   Static Initialization, Accessors, Destructors
*****************************************************************/
CSoilProfile **CSoilProfile::pAllSoilProfiles=NULL;
int            CSoilProfile::NumSoilProfiles=0;

//////////////////////////////////////////////////////////////////
/// \brief Returns number of soil profiles in model
/// \return Number of soil profiles in model
//
int CSoilProfile::GetNumProfiles()
{
  return NumSoilProfiles;
}

//////////////////////////////////////////////////////////////////
/// \brief Destroy all soil profiles in model
//
void CSoilProfile::DestroyAllSoilProfiles()
{
  if (DESTRUCTOR_DEBUG){cout <<"DESTROYING ALL SOIL PROFILES"<<endl;}
  for (int p=0; p<NumSoilProfiles;p++){
    delete pAllSoilProfiles[p];
  }
  delete [] pAllSoilProfiles;
}

///////////////////////////////////////////////////////////////////
/// \brief Summarize soil profile information to screen
//
void CSoilProfile::SummarizeToScreen()
{
  cout<<"==================="<<endl;
  cout<<"Soil Profile Summary:"<<NumSoilProfiles<<" soils in database"<<endl;
  for (int p=0; p<NumSoilProfiles;p++)
  {
    cout<<"-Soil profile \""<<pAllSoilProfiles[p]->GetTag()<<"\" "<<endl;
    cout<<"    #horizons:"<<pAllSoilProfiles[p]->GetNumHorizons()<<endl;
    for (int m=0; m<pAllSoilProfiles[p]->GetNumHorizons(); m++){
      cout <<"      -layer #"<<m+1<<": "<<pAllSoilProfiles[p]->GetSoilTag(m)<<" (thickness: "<< pAllSoilProfiles[p]->GetThickness(m)<<" m)"<<endl;
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Converts string to soil profile
/// \details Converts string (e.g., "ALL_SILT" in HRU file) to soil profile
///  can accept either soilprofile index or soilprofile tag
///  if string is invalid, returns NULL
/// \param s [in] String identifier of soil profile
/// \return Pointer to Soil profile to which passed string corresponds
//
const CSoilProfile *CSoilProfile::StringToSoilProfile(const string s)
{
  string sup=StringToUppercase(s);
  for (int p=0;p<NumSoilProfiles;p++)
  {
    if (!sup.compare(StringToUppercase(pAllSoilProfiles[p]->GetTag()))){return pAllSoilProfiles[p];}
    else if (s_to_i(s.c_str())==(p+1))            {return pAllSoilProfiles[p];}
  }
  return NULL;
}
