/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team, James R. Craig, Andy Snowdon
  ------------------------------------------------------------------
  Class CPropertyClassABC
  Class CVegetationClass
  Class CRiverClass
  Class CLandUseClass
  ----------------------------------------------------------------*/
#ifndef PROPERTY_CLASSES_H
#define PROPERTY_CLASSES_H

#include "RavenInclude.h"
#include "Properties.h"

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for property classification
/// \details Parent class fo soil, land use, river, etc. classes
//
class CPropertyClass
{
protected:/*----------------------------------------------------*/

  string                  tag;          ///< nickname for class, e.g., "BROADLEAF_FOREST"

  static CPropertyClass **pAllClasses;  ///< array of all classes that have been created
  static int              NumClasses;   ///< Number of classes

public:/*-------------------------------------------------------*/
  //Constructors:
  CPropertyClass(const string name);
  ~CPropertyClass();

  //Accessors
  string                        *GetTag() const;

  //routines
  static int                     GetNumClasses();
  static const CPropertyClass   *StringToClass(const string s);
  static virtual void            SummarizeToScreen()=0;
  static void                    DestroyAllClasses();
};
