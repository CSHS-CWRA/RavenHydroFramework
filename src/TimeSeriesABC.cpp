/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------*/
#include "TimeSeriesABC.h"


// Should be moved to RavenInclude??


/*****************************************************************
   Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/
///////////////////////////////////////////////////////////////////
/// \brief Implementation of the constructor of an abstract time serie
///
/// \param type      [in] type of time series being constructed (subclass)
/// \param Name      [in] name of time series
/// \param tag       [in] data tag
/// \param filename  [in] original source file
//
CTimeSeriesABC::CTimeSeriesABC(ts_type type,
                               string  Name,
                               string  tag,
                               string  filename)
{
  _type      =type;
  _name      =Name;
  _tag       =tag;
  _srcfile   =filename;
}

///////////////////////////////////////////////////////////////////
/// \brief Implementation of copy constructor (for which the address of a time series is passed)
/// \param &t [in] Address of a time series of which a "copy" is made
//
CTimeSeriesABC::CTimeSeriesABC(string Name,
                               const CTimeSeriesABC &t)
{
  _type=t.GetType();
  _name=Name;
  _tag       =t.GetTag();
  _srcfile = "";
}
///////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
//
CTimeSeriesABC::~CTimeSeriesABC()
{
}

/*****************************************************************
   Accessors
------------------------------------------------------------------
*****************************************************************/

///////////////////////////////////////////////////////////////////
/// \brief Returns type (subclass) of time series
/// \return type of time series
//
CTimeSeriesABC::ts_type CTimeSeriesABC::GetType() const{return _type;}

///////////////////////////////////////////////////////////////////
/// \brief Returns name of time series
/// \return name of time series
//
string CTimeSeriesABC::GetName()  const{return _name;}

//////////////////////////////////////////////////////////////////
/// \brief Returns data tag
/// \return data tag (e.g., HRU ID or Basin ID of observation data)
//
string   CTimeSeriesABC::GetTag()      const{return _tag;}

//////////////////////////////////////////////////////////////////
/// \brief Returns source input file
/// \return source file of time series (if one exists)
//
string   CTimeSeriesABC::GetSourceFile()      const{return _srcfile;}
