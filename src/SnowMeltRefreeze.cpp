/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------
  Melting
  Freezing
  ----------------------------------------------------------------*/
#include "HydroProcessABC.h"
#include "SnowMovers.h"

/*****************************************************************
   Snowmelt Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the Snow melt constructor
/// \param melt_type [in] Model of snow melt selected
/// \param Out_index [in] Index of storage compartment to which meltwater enters
//
CmvSnowMelt::CmvSnowMelt(snowmelt_type melt_type,int Out_index):
  CHydroProcessABC(SNOWMELT)
{
  type=melt_type;
  CHydroProcessABC::DynamicSpecifyConnections(1);
  ExitGracefullyIf(Out_index==DOESNT_EXIST,
                   "CmvSnowMelt Constructor: invalid 'to' compartment specified",BAD_DATA);

  iFrom[0]=pModel->GetStateVarIndex(SNOW);
  iTo  [0]=Out_index;
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation fo the default destructor
//
CmvSnowMelt::~CmvSnowMelt(){}

//////////////////////////////////////////////////////////////////
/// \brief Initializes snow melt modelling object
//
void CmvSnowMelt::Initialize(){}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for snow melt algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by snow melt algorithm (size of aP[] and aPC[])
//
void CmvSnowMelt::GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const
{
  if (type==MELT_POTMELT)
  {
    nP=0;
  }
  else
  {
    ExitGracefully("CmvSnowMelt::GetParticipatingParamList: undefined snow melt algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
///
/// \param melt_type [in] Snow melt modelling type selected
/// \param *aSV [out] Array of state variable types needed by snow melt algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by sublimation algorithm (size of aSV[] and aLev[] arrays)
//
void CmvSnowMelt::GetParticipatingStateVarList(snowmelt_type melt_type,sv_type *aSV, int *aLev, int &nSV)
{
  nSV=1;
  aSV[0]=SNOW; aLev[0]=DOESNT_EXIST;
  //variable user-specified output
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rates of loss of water from snowpack due to snowmelt $[mm/d]
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Rate of loss from "from" compartment [mm/d]
//
void CmvSnowMelt::GetRatesOfChange( const double                 *state_vars,
                                    const CHydroUnit *pHRU,
                                    const optStruct      &Options,
                                    const time_struct &tt,
                                    double     *rates) const
{
  if (pHRU->GetHRUType()==HRU_LAKE){return;}


  //------------------------------------------------------------
  if (type==MELT_POTMELT)
  {
    rates[0]=threshPositive(pHRU->GetForcingFunctions()->potential_melt);
  }
};

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow cannot drain "from" compartment over timestep
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Rate of loss from "from" compartment [mm/d]
//
void   CmvSnowMelt::ApplyConstraints( const double               *state_vars,
                                      const CHydroUnit *pHRU,
                                      const optStruct    &Options,
                                      const time_struct &tt,
                                      double     *rates) const
{
  if (pHRU->GetHRUType()==HRU_LAKE){return;}

  if (state_vars[iFrom[0]]<=0){rates[0]=0.0;}
  if (rates[0]<0.0){rates[0]=0.0;}//positivity constraint

  //cant remove more than is there
  rates[0]=threshMin(rates[0],state_vars[iFrom[0]]/Options.timestep,0.0);

}

//************************************************************************************************
//************************************************************************************************
//************************************************************************************************


/*****************************************************************
   SnowSqueeze Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/
//////////////////////////////////////////////////////////////////
/// \brief Implementation of the snow squeeze constructor
/// \param Out_index [in] Index of the compartment to which water is lost
//
CmvSnowSqueeze::CmvSnowSqueeze(int Out_index):
  CHydroProcessABC(SNOWSQUEEZE)
{
  ExitGracefullyIf(Out_index==DOESNT_EXIST,
                   "CmvSnowSqueeze Constructor: invalid 'to' compartment specified",BAD_DATA);

  CHydroProcessABC::DynamicSpecifyConnections(1);

  iFrom[0]=pModel->GetStateVarIndex(SNOW_LIQ);
  iTo  [0]=Out_index;
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation fo the default destructor
//
CmvSnowSqueeze::~CmvSnowSqueeze(){}

//////////////////////////////////////////////////////////////////
/// \brief Initializes snow squeeze modelling object
//
void CmvSnowSqueeze::Initialize(){}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for snow squeeze algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by snow squeeze algorithm (size of aP[] and aPC[])
//
void CmvSnowSqueeze::GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const
{
  nP=1;
  aP[0]="SNOW_SWI";              aPC[0]=CLASS_GLOBAL;
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
///
/// \param *aSV [out] Array of state variable types needed by snow squeeze algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by snow squeeze algorithm (size of aSV[] and aLev[] arrays)
//
void CmvSnowSqueeze::GetParticipatingStateVarList(sv_type *aSV, int *aLev, int &nSV)
{
  nSV=1;
  aSV[0]=SNOW_LIQ; aLev[0]=DOESNT_EXIST;
  //variable user-specified output
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rates of loss of water due to reduced snow porosity storage [mm/d]
//
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Rate of loss from "from" compartment [mm/d]
//
void CmvSnowSqueeze::GetRatesOfChange( const double              *state_vars,
                                       const CHydroUnit *pHRU,
                                       const optStruct   &Options,
                                       const time_struct &tt,
                                       double     *rates) const
{
  if (pHRU->GetHRUType()==HRU_LAKE){return;}

  double liq_cap,SD(0.0),S,SL;

  S=state_vars[pModel->GetStateVarIndex(SNOW)];
  SL=state_vars[pModel->GetStateVarIndex(SNOW_LIQ)];
  SD=S/0.2; //20% snow density

  if (pModel->GetStateVarIndex(SNOW_DEPTH)!=DOESNT_EXIST){
    SD=state_vars[pModel->GetStateVarIndex(SNOW_DEPTH)];
  }
  liq_cap=CalculateSnowLiquidCapacity(S,SD,Options);

  rates[0]=max(SL-liq_cap,0.0)/Options.timestep;

};

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow cannot drain "from" compartment over timestep
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Rate of loss from "from" compartment [mm/d]
//
void   CmvSnowSqueeze::ApplyConstraints( const double            *state_vars,
                                         const CHydroUnit *pHRU,
                                         const optStruct         &Options,
                                         const time_struct &tt,
                                         double     *rates) const
{
  if (pHRU->GetHRUType()==HRU_LAKE){return;}

  if (state_vars[iFrom[0]]<=0){rates[0]=0.0;}
  if (rates[0]<0.0){rates[0]=0.0;}//positivity constraint
}

//************************************************************************************************
//************************************************************************************************
//************************************************************************************************

/*****************************************************************
   Refreezing Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the snow refreeze constructor
/// \param frz_type [in] Refreeze algorithm type
//
CmvSnowRefreeze::CmvSnowRefreeze(refreeze_type frz_type):
  CHydroProcessABC(REFREEZE)
{
  type =frz_type;
  if (type==FREEZE_DEGREE_DAY){
    CHydroProcessABC::DynamicSpecifyConnections(1);

    iFrom[0]=pModel->GetStateVarIndex(SNOW_LIQ);
    iTo  [0]=pModel->GetStateVarIndex(SNOW);
  }
  else{
    ExitGracefully("CmvSnowRefreeze::Constructor: undefined snow refreeze type",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation fo the default destructor
//
CmvSnowRefreeze::~CmvSnowRefreeze(){}

//////////////////////////////////////////////////////////////////
/// \brief Initializes refreeze modelling object
//
void CmvSnowRefreeze::Initialize(){}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for refreeze algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by refreeze algorithm (size of aP[] and aPC[])
//
void CmvSnowRefreeze::GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const
{
  if (type==FREEZE_DEGREE_DAY)
  {
    nP=1;
    aP[0]="REFREEZE_FACTOR";       aPC[0]=CLASS_LANDUSE;
  }
  else
  {
    ExitGracefully("CmvSnowRefreeze::GetParticipatingParamList: undefined refreeze algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
///
/// \param frz_type [in] Refreeze modelling type selected
/// \param *aSV [out] Array of state variable types needed by snow melt algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by refreeze algorithm (size of aSV[] and aLev[] arrays)
//
void CmvSnowRefreeze::GetParticipatingStateVarList(refreeze_type frz_type,
                                                   sv_type *aSV, int *aLev, int &nSV)
{
  nSV=2;
  aSV[0]=SNOW_LIQ;    aLev[0]=DOESNT_EXIST;
  aSV[1]=SNOW;        aLev[1]=DOESNT_EXIST;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rates of loss of water to refreeze [mm/d]
//
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Rate of loss from "from" compartment [mm/d]
//
void CmvSnowRefreeze::GetRatesOfChange( const double             *state_vars,
                                        const CHydroUnit *pHRU,
                                        const optStruct  &Options,
                                        const time_struct &tt,
                                        double     *rates) const
{
  if (pHRU->GetHRUType()==HRU_LAKE){return;}

  double Kf=pHRU->GetSurfaceProps()->refreeze_factor;//[mm/K/d]

  if (type==FREEZE_DEGREE_DAY)
  {
    double Ta=pHRU->GetForcingFunctions()->temp_daily_ave;

    rates[0]=threshPositive(Kf*(FREEZING_TEMP-Ta));
  }
};

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details For all methods, ensures that rate of refreeze cannot freeze more
/// liquid water than is available
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Rate of loss from "from" compartment [mm/d]
//
void  CmvSnowRefreeze::ApplyConstraints(const double             *state_vars,
                                        const CHydroUnit *pHRU,
                                        const optStruct  &Options,
                                        const time_struct &tt,
                                        double     *rates) const
{
  if (pHRU->GetHRUType()==HRU_LAKE){return;}

  if (state_vars[iFrom[0]]<=0){rates[0]=0.0;}//reality check

  //can't remove more than is there
  rates[0]=threshMin(rates[0],state_vars[iFrom[0]]/Options.timestep,0.0);

  if (rates[0]<0)          {rates[0]=0.0;}//positivity constraint
}
