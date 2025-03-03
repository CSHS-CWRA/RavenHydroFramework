/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2023 the Raven Development Team
  ------------------------------------------------------------------
  Frozen lake ice
  ----------------------------------------------------------------*/
#include "HydroProcessABC.h"
#include "FrozenLake.h"
#include "EnergyTransport.h"

/*****************************************************************
   FrozenGround Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the standard constructor
/// \param lakefreeze_type [in] Model of lake freeze selected
/// \param pTransMod [in] pointer to transport model
/// \param pModel [in] - ponter to hydrologic model
//
CmvFrozenLake::CmvFrozenLake(lakefreeze_type       type,
                             const CTransportModel *pTransMod,
                             CModelABC             *pModel)
  :CHydroProcessABC(LAKE_FREEZING, pModel)
{
  _type =type;
  _pTransModel=pTransMod;

  CHydroProcessABC::DynamicSpecifyConnections(1);
  iFrom[0]=pModel->GetStateVarIndex(ICE_THICKNESS);
  iTo  [0]=pModel->GetStateVarIndex(ICE_THICKNESS);
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of default destructor
//
CmvFrozenLake::~CmvFrozenLake(){}

//////////////////////////////////////////////////////////////////
/// \brief does nothing
//
void CmvFrozenLake::Initialize()
{
  if (_type == LFREEZE_THERMAL) {
    const CEnthalpyModel *pEnthalpyModel=_pTransModel->GetEnthalpyModel();

    ExitGracefullyIf(pEnthalpyModel==NULL,"CmvFrozenLake:: Enthalpy is not simulated- must use :Transport TEMPERATURE",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for capillary rise algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by capillary rise  algorithm (size of aP[] and aPC[])
//
void CmvFrozenLake::GetParticipatingParamList(string  *aP , class_type *aPC , int &nP) const
{
  nP=0;
  if (_type==LFREEZE_BASIC)
  {
    nP=1;
    aP [0]="LAKESNOW_BUFFER_HT";     aPC [0]=CLASS_LANDUSE;
  }
  else if (_type == LFREEZE_THERMAL)
  {

  }
  else
  {
    ExitGracefully("CmvFrozenLake::GetParticipatingParamList: undefined frozen lake algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
/// \details User specifies from and to compartments, levels not known before construction
///
/// \param lakefreeze_type [in] algorithm used
/// \param *aSV  [out] Reference to state variable types needed by algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param  &nSV [out] Number of state variables required by algorithm (size of aSV[] and aLev[] arrays)
//
void CmvFrozenLake::GetParticipatingStateVarList(lakefreeze_type  btype,sv_type *aSV, int *aLev, int &nSV)
{
  nSV=0;
  aSV[0]=ICE_THICKNESS; aLev[0]=DOESNT_EXIST; nSV++;
  aSV[1]=SNOW;          aLev[1]=DOESNT_EXIST; nSV++;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rate of change of ice thickness [mm/d]
///
/// \param *state_vars [in] array of state variables
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of change of ice thickness [mm/day]
//
void   CmvFrozenLake::GetRatesOfChange( const double      *state_vars,
                                          const CHydroUnit  *pHRU,
                                          const optStruct   &Options,
                                          const time_struct &tt,
                                                double      *rates) const
{
  if ((pHRU->GetHRUType()!=HRU_LAKE) &&
      (pHRU->GetHRUType()!=HRU_WETLAND)){return;}//Lakes and wetlands only (includes reservoirs)

  double ice_thick=state_vars[iFrom[0]];

  //-----------------------------------------------------------------
  if (_type==LFREEZE_BASIC)
  {
    double pot_freeze=-pHRU->GetForcingFunctions()->potential_melt; //[mm/d]
    double buff_ht   = pHRU->GetSurfaceProps()->lakesnow_buffer_ht;

    int iSWE=pModel->GetStateVarIndex(SNOW);
    double SWE=state_vars[iSWE];

    double corr=1.0;
    if      (SWE>=buff_ht){corr=0.0;            }
    else if (SWE<=0.0    ){corr=1.0;            }
    else                  {corr=1.0-SWE/buff_ht;}

    rates[0]=corr*pot_freeze;
  }
  //-----------------------------------------------------------------
  else if (_type == LFREEZE_THERMAL)
  {
    //just queries thermal model and updates ice thickness variable
    // \todo[funct] need to allow snow buffering

    ice_thick=state_vars[iFrom[0]];

    const CEnthalpyModel *pEnthalpyModel=_pTransModel->GetEnthalpyModel();

    int    iSW     =pModel->GetStateVarIndex(SURFACE_WATER);
    double ice_frac=pEnthalpyModel->GetIceContent(state_vars,iSW);

    double new_ice_thick=ice_frac*state_vars[iSW];

    rates[0]=(new_ice_thick-ice_thick)/Options.timestep;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
///
/// \param *state_vars [in] state variable array for this HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of loss from soil to other soil layer [mm/day]
//
void   CmvFrozenLake::ApplyConstraints( const double     *state_vars,
                                           const CHydroUnit *pHRU,
                                           const optStruct  &Options,
                                           const time_struct &tt,
                                           double     *rates) const
{
  double ice_thick=state_vars[iFrom[0]];

  rates[0]=max(rates[0],-ice_thick/Options.timestep); //thickness cannot be less than zero
}
