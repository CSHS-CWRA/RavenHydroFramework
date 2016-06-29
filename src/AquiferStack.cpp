/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright © 2008-2014 the Raven Development Team
----------------------------------------------------------------*/
#include "Properties.h"
#include "SoilAndLandClasses.h"
#include "SoilProfile.h"
/*****************************************************************
   Constructor / Destructor
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the aquifer stack constructor
/// \param name [in] Aquifer stack identifier
CAquiferStack::CAquiferStack(const string name)
{
  tag=name;
	if (!DynArrayAppend((void**&)(pAllAqStacks),(void*)(this),NumAqStacks)){
		 ExitGracefully("CAquiferStack::Constructor: creating NULL soil profile",BAD_DATA);};

  nLayers        =0;
	pAquiferSoils  =NULL;
	aAquifer_thick =NULL;
	pAquitardSoils =NULL;
	aAquitard_thick=NULL;
}

///////////////////////////////////////////////////////////////////
/// \brief Implentation of the destructor
//
CAquiferStack::~CAquiferStack()
{
	if (DESTRUCTOR_DEBUG){cout<<"  DELETING AQUIFER STACK "<<endl;}
	delete [] pAquiferSoils;
	delete [] aAquifer_thick;
	delete [] pAquitardSoils;
	delete [] aAquitard_thick;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns tag of soil profile (e.g. "ALL_SILT")
/// \return Tag of soil profile 
//
string CAquiferStack::GetTag() const{return tag;}

///////////////////////////////////////////////////////////////////
/// \brief Returns pointer to soil layer in aquifer stack
/// \param m [in] Index of aquifer soil layer to be returned
/// \return aquifer soil corresponding to layer m
//
const soil_struct *CAquiferStack::GetAquiferSoil(const int m) const
{
	ExitGracefullyIf((m<0) || (m>=nLayers),
		"CAquiferStack::GetAquiferSoil: bad aquifer layer index",BAD_DATA);
	return pAquiferSoils[m]->GetSoilStruct();
}

///////////////////////////////////////////////////////////////////
/// \brief Returns thickness of aquifer layer in meters [m]
/// \param m [in] Index of a specific aquifer layer
/// \return Thickness of aquifer layer [m]
//
double CAquiferStack::GetAquiferThickness(const int m) const
{
	ExitGracefullyIf((m<0) || (m>=nLayers),
		"CAquiferStack::GetAquiferThickness: bad aquifer layer index",BAD_DATA);
	return aAquifer_thick[m];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns tag of aquifer layer specified by m
/// \param m [in] Aquifer layer specifier
/// \return aquifer layer specified by m
//
string CAquiferStack::GetAquiferSoilTag (const int m) const
{
	ExitGracefullyIf((m<0) || (m>=nLayers),
		"CAquiferStack::GetAquiferSoilTag: bad aquifer layer index",BAD_DATA);
	return pAquiferSoils[m]->GetTag();
}

///////////////////////////////////////////////////////////////////
/// \brief Returns the number of horizons in soil profile
/// \return Number of horizons in soil profile
double CAquiferStack::GetNumLayers () const
{
	return nLayers;
}

//////////////////////////////////////////////////////////////////
/// \brief Adds aquifer  of specified soil class with specified thickness to bottom of aquifer stack
/// \details adds null aquitard pointer and aquitard thickness of zero
/// \param thickness [in] Thickness of layer to be added [m]
/// \param *pLayerSoil [in] pointer to soil class to be added to collection
//
void CAquiferStack::AddLayer(double thickness, //[m]
                             const CSoilClass *pLayerSoil)
{
	if (!DynArrayAppend((void**&)(pAquiferSoils),(void*)(pLayerSoil),nLayers)){
		 ExitGracefully("CAquiferStack::Constructor: creating NULL aquifer layer",BAD_DATA);};
	DynArrayAppend((void**&)(pAquitardSoils),NULL,nLayers); //Aquitard is null,

  double *tmpthick  =new double [nLayers];
  double *tmpaqthick=new double [nLayers];
	for (int m=0;m<nLayers-1;m++){//copy previous arrays
		tmpthick  [m]=aAquifer_thick[m];
    tmpaqthick[m]=aAquitard_thick[m];
	}
	tmpthick  [nLayers-1]=thickness;//add new data
  tmpaqthick[nLayers-1]=0.0;      //thickness of 0.0 for null aquifer
	delete [] aAquifer_thick;//delete previous arrays
  delete [] aAquitard_thick;
	aAquifer_thick=tmpthick;
  aAquitard_thick=tmpaqthick;
}
//////////////////////////////////////////////////////////////////
/// \brief Adds aquifer and aquitard of specified soil classes with specified thicknesses to bottom of aquifer stack
/// \param thickness [in] Thickness of aquifer layer to be added [m]
/// \param *pLayerSoil [in] pointer to soil class to be added to collection
/// \param aquitard_thick [in] Thickness of aquitard to be added [m]
/// \param *pAquitardSoil [in] pointer to aquitard soil class
//
void CAquiferStack::AddLayer(double            thickness, //[m]
                             const CSoilClass *pLayerSoil,
                             double            aquitard_thick,
                             const CSoilClass *pAquitardSoil)
{
	if (!DynArrayAppend((void**&)(pAquiferSoils),(void*)(pLayerSoil),nLayers)){
		 ExitGracefully("CAquiferStack::Constructor: creating NULL aquifer layer",BAD_DATA);};
	if (!DynArrayAppend((void**&)(pAquitardSoils),(void*)(pAquitardSoil),nLayers)){
		 ExitGracefully("CAquiferStack::Constructor: creating NULL aquitard layer",BAD_DATA);};

  double *tmpthick  =new double [nLayers];
  double *tmpaqthick=new double [nLayers];
	for (int m=0;m<nLayers-1;m++){//copy previous arrays
		tmpthick  [m]=aAquifer_thick[m];
    tmpaqthick[m]=aAquitard_thick[m];
	}
	tmpthick  [nLayers-1]=thickness;//add new data
  tmpaqthick[nLayers-1]=aquitard_thick;//add new data
	delete [] aAquifer_thick;//delete previous arrays
  delete [] aAquitard_thick;
	aAquifer_thick=tmpthick;
  aAquitard_thick=tmpaqthick;
}
//////////////////////////////////////////////////////////////////
/// \brief Allocate soil layers
//
/// \details Discretizes soil profile into nAqLayers layers for numerical 
///	solution, e.g., takes a 7 meter soil profile with three horizons
///	and breaks it up into 12 layers. 
///	The algorithm shoots for at least one layer per profile, then 
///	aims for higher discretization in top layers, coarser 
///	discretization at depth. Thicknesses are in metres
//
/// \param nAqLayers [in] Number of aquifer layers in model (may be more finely discretized than physical aquifer layers)
/// \param **pSoils [out] Array of pointers to aquifer soil properties
/// \param *thickness [out] Reference to array of layer aAquifer_thick [m]
//
void CAquiferStack::AllocateAqLayers     (const int           nAqLayers,
                                          const soil_struct **pSoils,
                                          double             *thickness) const 
{
  static soil_struct blank_soil_struct;
  
  CSoilClass::InitializeSoilProperties(blank_soil_struct, false);
  blank_soil_struct.porosity =0;
  blank_soil_struct.hydraul_cond=0;
	
	int m;
  if (nLayers==0){ //special case for lakes and glaciers
    for (m=0;m<nAqLayers;m++){
      pSoils   [m]=&blank_soil_struct;
      thickness[m]=REAL_SMALL;
    }
    return;
  }
	else if (nAqLayers>=nLayers)
  {
		for (m=0;m<nLayers;m++){
			pSoils   [m]=pAquiferSoils[m]->GetSoilStruct();
			thickness[m]=aAquifer_thick [m];
		}
		for (m=nLayers;m<nAqLayers;m++){
			pSoils   [m]=&blank_soil_struct;
			thickness[m]=0.0;
		}
		return;
  }
  else
  {
    cout <<"nAqLayers:"<<nAqLayers<<" nLayers:"<<nLayers<<endl;
	  ExitGracefully("CAquiferStack::AllocateAqLayers: cannot currently handle # of aquifer layers combination",BAD_DATA);
  }
}
/*****************************************************************
   Static Initialization, Accessors, Destructors
*****************************************************************/
CAquiferStack **CAquiferStack::pAllAqStacks=NULL;
int            CAquiferStack::NumAqStacks=0;

//////////////////////////////////////////////////////////////////
/// \brief Returns number of soil profiles in model
/// \return Number of soil profiles in model
//
int CAquiferStack::GetNumAqStacks()
{
  return NumAqStacks;
}

//////////////////////////////////////////////////////////////////
/// \brief Destroy all soil profiles in model
//
void CAquiferStack::DestroyAllAqStacks()
{
	if (DESTRUCTOR_DEBUG){cout <<"DESTROYING ALL AQUIFER STACKS"<<endl;}
  for (int p=0; p<NumAqStacks;p++){
    delete pAllAqStacks[p];
  }
  delete [] pAllAqStacks;
}

//////////////////////////////////////////////////////////////////
/// \brief Converts string to soil profile
/// \details Converts string (e.g., "ALL_SILT" in HRU file) to soil profile
///  can accept either soilprofile index or soilprofile tag
///  if string is invalid, returns NULL
/// \param s [in] String identifier of soil profile
/// \return Pointer to Soil profile to which passed string corresponds
//
const CAquiferStack *CAquiferStack::StringToAqStack(const string s)
{
  string name=StringToUppercase(s);
  for (int p=0;p<NumAqStacks;p++)
  {
    if (!name.compare(pAllAqStacks[p]->GetTag())){return pAllAqStacks[p];}
    else if (s_to_i(name.c_str())==(p+1))        {return pAllAqStacks[p];}
  }
  return NULL;
}