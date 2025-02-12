/*------------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ------------------------------------------------------------------
  Recharge
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "SoilWaterMovers.h"
#include "GroundwaterModel.h"
#include "Model.h"

/*****************************************************************
   Recharge Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the standard constructor
/// \param cr_type [in] Model of capillary rise selected
/// \param In_index [in] Soil storage unit index from which water is lost
/// \param Out_index [in] Soil storage unit index to which water rises
//
CmvRecharge::CmvRecharge(recharge_type rech_type,
                         int           to_index,
                         int           junk,
                         CModelABC     *pModel)
  :CHydroProcessABC(RECHARGE, pModel)
{
  ExitGracefullyIf(to_index==DOESNT_EXIST,
                   "CmvRecharge Constructor: invalid 'to' compartment specified",BAD_DATA);

  CHydroProcessABC::DynamicSpecifyConnections(1);
  iFrom[0]=pModel->GetStateVarIndex(ATMOS_PRECIP);
  iTo  [0]=to_index;
  _type = rech_type;
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the recharge constructor
/// \param se_type [in] Model of recharge selected
//
CmvRecharge::CmvRecharge(recharge_type rech_type,
                         int           nConns,
                         CModelABC     *pModel)
  :CHydroProcessABC(RECHARGE, nConns, pModel)
{
  int iGW, i;
  _type = rech_type;
  if (_type==RECHARGE_CONSTANT)
  {
    CHydroProcessABC::DynamicSpecifyConnections(1);

    iGW  =pModel->GetStateVarIndex(GROUNDWATER);   //****fix to not always go to top layer of aquifer?
    iFrom[0]=pModel->GetStateVarIndex(PONDED_WATER);     iTo[0]=iGW;
  }
  else if(_type==RECHARGE_CONSTANT_OVERLAP)
  {
    CHydroProcessABC::DynamicSpecifyConnections(nConns);

    iGW  =pModel->GetStateVarIndex(GROUNDWATER);   //****fix to not always go to top layer of aquifer?

    for(i = 0; i < nConns; i++)
    {
      iFrom[i]=pModel->GetStateVarIndex(PONDED_WATER);
      iTo[i]=iGW;
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Implementation of default destructor
//
CmvRecharge::~CmvRecharge()
{
}
//////////////////////////////////////////////////////////////////
/// \brief Verifies iFrom - iTo connectivity
/// \details Ensures that water rises from soil layer/ groundwater
/// to another soil layer or groundwater unit
//
void CmvRecharge::Initialize()
{
  ExitGracefullyIf((pModel->GetStateVarType(iTo[0])!=SOIL) &&
                   (pModel->GetStateVarType(iTo[0])!=GROUNDWATER),
                   "CmvRecharge::Initialize:Recharge must go to soil or groundwater unit",BAD_DATA);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for capillary rise algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by capillary rise  algorithm (size of aP[] and aPC[])
//
void CmvRecharge::GetParticipatingParamList(string  *aP , class_type *aPC , int &nP) const
{
  if (_type == RECHARGE_CONSTANT)
  {
    nP = 0;
  }
  else if (_type == RECHARGE_FROMFILE)
  {
    nP=0;
  }
  else if(_type==RECHARGE_CONSTANT_OVERLAP)
  {
    nP=0;
  }
  else
  {
    ExitGracefully("CmvRecharge::GetParticipatingParamList: undefined recharge algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
/// \details User specifies from and to compartments, levels not known before construction
///
/// \param *aSV [out] Reference to state variable types needed by cr algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by cr algorithm (size of aSV[] and aLev[] arrays)
//
void CmvRecharge::GetParticipatingStateVarList(recharge_type	r_type,sv_type *aSV, int *aLev, int &nSV)
{
  if (r_type==RECHARGE_FROMFILE){
    nSV=0;
    //soil and atmos_precip will already be in model by default
  }
  else if (r_type==RECHARGE_CONSTANT)
  {
    nSV=2;
    aSV [0]=PONDED_WATER; aSV	 [1]=GROUNDWATER;
    aLev[0]=DOESNT_EXIST; aLev [1]=0;
  }
  else if (r_type==RECHARGE_CONSTANT_OVERLAP)
  {
    nSV=2;
    aSV [0]=PONDED_WATER; aSV	 [1]=GROUNDWATER;
    aLev[0]=DOESNT_EXIST; aLev [1]=0;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rate to aquifer due to Recharge[mm/day]
/// \details Note that the formatting issues below will be resolved once the references have been extracted.. \n \n
///
/// \param *storage [in] Reference to soil storage from which water rises
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of loss of water from soil to another soil layer [mm/day]
//
void   CmvRecharge::GetRatesOfChange( const double      *storage,
                                      const CHydroUnit  *pHRU,
                                      const optStruct   &Options,
                                      const time_struct &tt,
                                      double            *rates) const
{
  if (pHRU->GetHRUType()!=HRU_STANDARD){return;}//Lakes & glaciers

  if( _type==RECHARGE_FROMFILE){
    rates[0]=pHRU->GetForcingFunctions()->recharge;
  }
  //------------------------------------------------------------
  else if(_type==RECHARGE_CONSTANT)
  {
    double rech = 1.0;  // GWUPDATE How's that for constant?
    rates[0] = rech;     //[mm/d]

  }
  else if(_type==RECHARGE_CONSTANT_OVERLAP)
  {
    ExitGracefully("CmvRecharge::GetRatesOfChange: undefined recharge type:RECHARGE_CONSTANT_OVERLAP",STUB);
  }
  else
  {
    ExitGracefully("CmvRecharge::GetRatesOfChange: undefined recharge type",BAD_DATA);
  }//end recharge type select*/
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow cannot drain "from" compartment over timestep.
///
/// \param *storage [in] state variable array for this HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of loss from soil to other soil layer [mm/day]
//
void   CmvRecharge::ApplyConstraints( const double     *state_vars,
                                           const CHydroUnit *pHRU,
                                           const optStruct  &Options,
                                           const time_struct &tt,
                                           double     *rates) const
{
  if (pHRU->GetHRUType()!=HRU_STANDARD){return;}//Lakes & glaciers

  if (_type==RECHARGE_FROMFILE)
  {
    //exceedance of max "to" compartment
    //water flow simply slows (or stops) so that receptor will not overfill during tstep
    if (!Options.allow_soil_overfill){
      rates[0]=min(rates[0],(pHRU->GetStateVarMax(iTo[0],state_vars,Options)-state_vars[iTo[0]])/Options.timestep);
    }
  }
  else
  {
     //cant remove more than is there
  	 rates[0]=min(rates[0],state_vars[iFrom[0]]/Options.timestep);
  }
}
