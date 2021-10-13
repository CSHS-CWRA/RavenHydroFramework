/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ----------------------------------------------------------------*/
#include "TimeSeriesABC.h"

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
                               long    loc_ID,
                               string  filename="")
{
  _type       =type;
  _name       =Name;
  _loc_ID     =loc_ID;
  _srcfile    =filename;
  _constit_ind=DOESNT_EXIST;
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
  _loc_ID=t.GetLocID();
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
/// \brief Returns location identifier
/// \return location identifier (e.g., HRU ID or Basin ID of observation data)
//
long   CTimeSeriesABC::GetLocID()      const{return _loc_ID;}

//////////////////////////////////////////////////////////////////
/// \brief Returns location identifier
/// \return constituent index (or DOESNT_EXIST, otherwise) 
//
int    CTimeSeriesABC::GetConstitInd() const { return _constit_ind; }

//////////////////////////////////////////////////////////////////
/// \brief Sets location identifier
/// \param ID - HRU ID or SubBasin ID linked to observations
//
void   CTimeSeriesABC::SetLocID(long ID) { _loc_ID=ID; }

//////////////////////////////////////////////////////////////////
/// \brief Sets constituent index
/// \param c - constituent index
//
void   CTimeSeriesABC::SetConstitInd(const int c) {_constit_ind=c; }

//////////////////////////////////////////////////////////////////
/// \brief Returns source input file
/// \return source file of time series (if one exists)
//
string   CTimeSeriesABC::GetSourceFile()      const{return _srcfile;}
