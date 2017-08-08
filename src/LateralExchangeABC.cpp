/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------*/
#include "LateralExchangeABC.h"
#include "HydroProcessABC.h"
/*****************************************************************
   Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/
CLateralExchangeProcessABC::CLateralExchangeProcessABC(const process_type ptype) : CHydroProcessABC(ptype)
{

  _nLatConnections=0;
  _kFrom=NULL;
  _kTo=NULL;
  _iFromLat=NULL;
  _iToLat=NULL;

  _LatFlowIndex=_nLatFlowProcesses;
  _nLatFlowProcesses++;
}
/*****************************************************************
   Static members
------------------------------------------------------------------
*****************************************************************/
int CLateralExchangeProcessABC::_nLatFlowProcesses=0;
const CModel *CLateralExchangeProcessABC::_pModel=NULL;
void CLateralExchangeProcessABC::SetModel(const CModel *pM){_pModel=pM;}
/*****************************************************************
   Private Member Functions
------------------------------------------------------------------
*****************************************************************/
void CLateralExchangeProcessABC::DynamicSpecifyLatConnections(const int nLatConnects)
{
  delete [] _kFrom; //if these have already been created
  delete [] _kTo;
  delete [] _iFromLat; 
  delete [] _iToLat;
  _iToLat=NULL;
  _nLatConnections=nLatConnects;
  _kFrom   =new int [_nLatConnections];
  _kTo     =new int [_nLatConnections];
  _iFromLat=new int [_nLatConnections];
  _iToLat  =new int [_nLatConnections];
  for (int q=0;q<_nLatConnections;q++)
  {
    _kFrom   [q]    =DOESNT_EXIST;
    _kTo     [q]    =DOESNT_EXIST;
    _iFromLat[q]    =DOESNT_EXIST;
    _iToLat  [q]    =DOESNT_EXIST;
  }
}
/*****************************************************************
   Accessors
------------------------------------------------------------------
*****************************************************************/
//////////////////////////////////////////////////////////////////
/// \brief Returns global lateral flow process index (ranges from 0 ot _nLatFlowProcesses)
///
/// \return global lateral flow process index
//
int CLateralExchangeProcessABC::GetLateralFlowIndex()
{
  return _LatFlowIndex;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of lateral connections between HRUs
///
/// \return (integer) number of lateral connections between HRUs
//
int CLateralExchangeProcessABC::GetNumLatConnections()
{
  return _nLatConnections;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns array of global indices (k's) of 'from' HRUs
///
/// \return (int*) array of global indices (k's) of 'from' HRUs
//
const int *CLateralExchangeProcessABC::GetFromHRUIndices()
{
  return _kFrom;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns array of global indices (k's) of 'to' HRUs
///
/// \return (int*) array of global indices (k's) of 'to' HRUs
//
const int *CLateralExchangeProcessABC::GetToHRUIndices()
{
  return _kTo;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns array of global indices (i's) of 'from' state variables
///
/// \return (int*) array of global indices (i's) of 'from' state variables
//
const int *CLateralExchangeProcessABC::GetLateralFromIndices()
{
  return _iFromLat;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns array of global indices (i's) of 'to' state variables
///
/// \return (int*) array of global indices (i's) of 'to' state variables
//
const int *CLateralExchangeProcessABC::GetLateralToIndices()
{
  return _iToLat;
}

/*

  void CLateralExchangeProcessABC::GetRatesOfChange(const double      *state_vars,
                              const CHydroUnit  *pHRU,
                              const optStruct   &Options,
                              const time_struct &tt,
                                    double      *rates) const=0;//still purely virtual

  void CLateralExchangeProcessABC::ApplyConstraints(const double      *state_vars,
                              const CHydroUnit  *pHRU,
                              const optStruct   &Options,
                              const time_struct &tt,
                                   double      *rates) const=0;//still purely virtual*/
