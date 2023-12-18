/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2022 the Raven Development Team
  ----------------------------------------------------------------
  Canopy Evaporation
  Canopy Snow Sublimation
  Canopy Drip
  Canopy Snow Melt
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "VegetationMovers.h"
#include "Model.h"

/*****************************************************************
   Canopy Evaporation Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/
//////////////////////////////////////////////////////////////////
/// \brief Implementation of the canopy evaporation constructor
/// \param cetype [in] Model of canopy evaporation
//
CmvCanopyEvap::CmvCanopyEvap(canevap_type cetype,
                             CModelABC    *pModel)
  :CHydroProcessABC(CANOPY_EVAPORATION, pModel)
{
  type =cetype;

  CHydroProcessABC::DynamicSpecifyConnections(2);//nConnections=2
  iFrom[0]=pModel->GetStateVarIndex(CANOPY);
  iTo  [0]=pModel->GetStateVarIndex(ATMOSPHERE);
  iFrom[1]=pModel->GetStateVarIndex(AET);
  iTo  [1]=pModel->GetStateVarIndex(AET);
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the standard destructor
//
CmvCanopyEvap::~CmvCanopyEvap(){}

//////////////////////////////////////////////////////////////////
/// \brief Validate iTo/iFrom connectivity of the evaporation process
//
void CmvCanopyEvap::Initialize()
{
  ExitGracefullyIf(pModel->GetStateVarType(iFrom[0])!=CANOPY,
                   "CmvCanopyEvap::Initialize:Canopy evaporation must come from canopy unit",BAD_DATA);
  ExitGracefullyIf(pModel->GetStateVarType(iTo[0])!=ATMOSPHERE,
                   "CmvCanopyEvap::Initialize:Canopy evaporation must go to atmosphere",BAD_DATA);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for canopy evaporation algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by canopy evaporation algorithm (size of aP[] and aPC[])
//
void CmvCanopyEvap::GetParticipatingParamList(string  *aP , class_type *aPC , int &nP) const
{
  if (type==CANEVP_RUTTER)
  {
    nP=3;
    aP[0]="FOREST_COVERAGE"; aPC[0]=CLASS_LANDUSE;  //JRCFLAG
    aP[1]="MAX_CAPACITY";    aPC[1]=CLASS_VEGETATION;
    aP[2]="TRUNK_FRACTION";  aPC[2]=CLASS_VEGETATION;
  }
  else if (type==CANEVP_MAXIMUM)
  {
    nP=1;
    aP[0]="FOREST_COVERAGE"; aPC[0]=CLASS_LANDUSE; //JRCFLAG
  }
  else if (type==CANEVP_ALL)
  {
    nP=0;
  }
  else
  {
    ExitGracefully("CmvCanopyEvap::GetParticipatingParamList: undefined canopy evaporation algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
///
/// \param cetype [in] Model of canopy evaporation used
/// \param *aSV [out] Array of state variable types needed by canopy evaporation algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by canopy evaporation algorithm (size of aSV[] and aLev[] arrays)
//
void CmvCanopyEvap::GetParticipatingStateVarList(canevap_type cetype,sv_type *aSV, int *aLev, int &nSV)
{
  nSV=3;
  aSV[0]=CANOPY;     aLev[0]=DOESNT_EXIST;
  aSV[1]=ATMOSPHERE; aLev[1]=DOESNT_EXIST;
  aSV[2]=AET;        aLev[2]=DOESNT_EXIST;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rate of loss of water from canopy to atmosphere [mm/d]
/// \details  if type == CANEVP_RUTTER (rutter model), evaporation is proportional to canopy storage \n
/// elseif type == CANEVP_MAXIMUM, evaporation is at PET
///
/// \param *state_vars [in] Array of current state variables for HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this process takes place
/// \param *rates [out] Rate of water loss from canopy to atmosphere [mm/d]
//
void CmvCanopyEvap::GetRatesOfChange( const double      *state_vars,
                                      const CHydroUnit  *pHRU,
                                      const optStruct   &Options,
                                      const time_struct &tt,
                                      double     *rates) const
{
  if ((pHRU->GetHRUType()!=HRU_STANDARD) &&
      (pHRU->GetHRUType()!=HRU_WETLAND)){return;}

  double Fc=pHRU->GetSurfaceProps()->forest_coverage;
  double cap=pHRU->GetVegVarProps()->capacity;
  rates[0]=0.0;//default
  if (Fc==0){return;}

  double PET=max(pHRU->GetForcingFunctions()->PET,0.0) ;
  double stor=min(max(state_vars[iFrom[0]],0.0),cap*Fc); //correct for potentially invalid storage

  double PETused=0.0;//[mm/d]
  if(!Options.suppressCompetitiveET) {
    PET-=(state_vars[pModel->GetStateVarIndex(AET)]/Options.timestep);
    PET=max(PET,0.0);
  }

  if (type==CANEVP_RUTTER)//-------------------------------------
  {
    double Ft =pHRU->GetVegetationProps()->trunk_fraction;
    if (pModel->GetStateVarIndex(TRUNK)==DOESNT_EXIST){Ft=0.0;}//overrides if trunk not explicitly modeled

    rates[0]=(1.0-Ft)*Fc*PET*(stor/(cap*Fc));
    PETused=rates[0];
  }
  else if (type==CANEVP_MAXIMUM)//-------------------------------
  {
    rates[0]=Fc*PET;
    PETused=rates[0];
  }
  else if (type==CANEVP_ALL)//-----------------------------------
  {
    //all canopy mass evaporates 'instantaneously'
    rates[0]=state_vars[iFrom[0]]/Options.timestep;
    PETused=rates[0];
  }
  else//---------------------------------------------------------
  {
    ExitGracefully("CmvCanopyEvap: this process not coded yet",STUB);
  }

  rates[1]=PETused;
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow cannot drain "from" compartment over timestep
///
/// \param *state_vars [in] Array of current state variables for HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this process takes place
/// \param *rates [out] Rate of water loss from canopy to atmosphere [mm/d]

//
void CmvCanopyEvap::ApplyConstraints( const double      *state_vars,
                                      const CHydroUnit  *pHRU,
                                      const optStruct   &Options,
                                      const time_struct &tt,
                                      double      *rates) const
{
  if ((pHRU->GetHRUType()!=HRU_STANDARD) &&
      (pHRU->GetHRUType()!=HRU_WETLAND)){return;}

  double oldRates=rates[0];

  //must be positive
  rates[0]=max(rates[0],0.0);

  //cant remove more than is there
  rates[0]=min(rates[0],state_vars[iFrom[0]]/Options.timestep);

  //update AET
  rates[1]-=(oldRates-rates[0]);
}

//*****************************************************************************************
//*****************************************************************************************
//*****************************************************************************************

/*****************************************************************
   Canopy Snow Sublimation Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the snow evaporation constructor
/// \param cetype [in] Model of canopy snow evaporation
//
CmvCanopySublimation::CmvCanopySublimation(sublimation_type cetype,
                                           CModelABC       *pModel)
  :CHydroProcessABC(CANOPY_SNOW_EVAPORATION, pModel)
{
  type =cetype;

  CHydroProcessABC::DynamicSpecifyConnections(2);//nConnections=2
  iFrom[0]=pModel->GetStateVarIndex(CANOPY_SNOW);
  iTo  [0]=pModel->GetStateVarIndex(ATMOSPHERE);
  iFrom[1]=iTo[1]=pModel->GetStateVarIndex(AET);
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the standard destructor
//
CmvCanopySublimation::~CmvCanopySublimation(){}

//////////////////////////////////////////////////////////////////
/// \brief Validate iTo/iFrom connectivity of the snowpack evaporation
//
void CmvCanopySublimation::Initialize()
{
  ExitGracefullyIf(pModel->GetStateVarType(iFrom[0])!=CANOPY_SNOW,
                   "CmvCanopySublimation::Initialize:Canopy evaporation must come from canopy unit",BAD_DATA);
  ExitGracefullyIf(pModel->GetStateVarType(iTo[0])!=ATMOSPHERE,
                   "CmvCanopySublimation::Initialize:Canopy evaporation must go to atmosphere",BAD_DATA);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for canopy snowpack sublimation algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by canopy snowpack sublimation algorithm (size of aP[] and aPC[])
//
void CmvCanopySublimation::GetParticipatingParamList(string  *aP , class_type *aPC , int &nP) const
{
  if (type==SUBLIM_MAXIMUM)
  {
    nP=1;
    aP[0]="FOREST_COVERAGE"; aPC[0]=CLASS_LANDUSE; //JRCFLAG
  }
  else if(type==SUBLIM_SVERDRUP)
  {
    nP=1;
    aP[0]="SNOW_ROUGHNESS";   aPC[0]=CLASS_GLOBAL;
  }
  else //most don't have parameters
  {
    nP=0;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
///
/// \param cetype [in] Model of canopy snowpack sublimation used
/// \param *aSV [out] Array of state variable types needed by canopy snowpack sublimation algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by canopy snowpack sublimation algorithm (size of aSV[] and aLev[] arrays)
//
void CmvCanopySublimation::GetParticipatingStateVarList(sublimation_type cs_type,sv_type *aSV, int *aLev, int &nSV)
{
  nSV=3;
  aSV[0]=CANOPY_SNOW; aLev[0]=DOESNT_EXIST;
  aSV[1]=ATMOSPHERE;  aLev[1]=DOESNT_EXIST;
  aSV[2]=AET;         aLev[2]=DOESNT_EXIST;
}

double  SublimationRate(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const CModelABC   *pModel,
                        const time_struct &tt,
                        const double      &wind_vel,
                        sublimation_type   type); //defined in Sublimation.cpp

//////////////////////////////////////////////////////////////////
/// \brief Returns rate of loss of water from canopy to atmosphere [mm/d]
/// \details  if type == CANEVP_MAXIMUM, evaporation is at PET
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of loss from canopy SWE to atmosphere [mm/d]
//
void CmvCanopySublimation::GetRatesOfChange(  const double      *state_vars,
                                              const CHydroUnit  *pHRU,
                                              const optStruct   &Options,
                                              const time_struct &tt,
                                              double      *rates) const
{
  if ((pHRU->GetHRUType()!=HRU_STANDARD) &&
      (pHRU->GetHRUType()!=HRU_WETLAND)){return;}

  double Fc=pHRU->GetSurfaceProps()->forest_coverage;
  rates[0]=0.0;//default
  if (Fc==0){return;}

  double PET=max(pHRU->GetForcingFunctions()->PET,0.0);

  double PETused=0.0;//[mm/d]
  if(!Options.suppressCompetitiveET) {
    PET-=(state_vars[pModel->GetStateVarIndex(AET)]/Options.timestep);
    PET=max(PET,0.0);
  }

  if (type==SUBLIM_MAXIMUM)//----------------------------------
  {
    //all canopy mass sublimates 'instantaneously' (up to threshold) based upon PET
    rates[0]=Fc*PET;
    PETused=rates[0];
  }
  else if (type==SUBLIM_ALL)//----------------------------------
  {
    //all canopy mass sublimates 'instantaneously' to atmosphere
    rates[0]=state_vars[iFrom[0]]/Options.timestep;
    PETused=rates[0];
  }
  else//--------------------------------------------------------
  {
    double wind_vel=pHRU->GetForcingFunctions()->wind_vel;//pModel->WindspeedAtHeight(pHRU->GetVegVarProps()->height, Options,pHRU,F,2.0);

    ExitGracefully("SUBLIMATION CANOPY - must adjust wind velocity", STUB);
    rates[0] = Fc*SublimationRate(state_vars, pHRU, pModel, tt, wind_vel, type);
  }
  rates[1]=PETused;
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow cannot drain "from" compartment over timestep
///
/// \param *state_vars [in] Array of state variables for HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of loss from canopy to atmosphere [mm/d]
//
void CmvCanopySublimation::ApplyConstraints( const double      *state_vars,
                                          const CHydroUnit  *pHRU,
                                          const optStruct   &Options,
                                          const time_struct &tt,
                                          double      *rates) const
{
  if ((pHRU->GetHRUType()!=HRU_STANDARD) &&
      (pHRU->GetHRUType()!=HRU_WETLAND)){return;}

  double oldRates=rates[0];

  //cant remove more than is there
  rates[0]=min(rates[0],state_vars[iFrom[0]]/Options.timestep);

  //update AET
  rates[1]-=(oldRates-rates[0]);

}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the canopy drip constructor
/// \param cdtype [in] Model of canopy drip
/// \param to_index [in] Index of storage compartment to which water is lost
//
CmvCanopyDrip::CmvCanopyDrip(candrip_type cdtype,
                             int          to_index,
                             CModelABC    *pModel)
  :CHydroProcessABC(CANOPY_DRIP, pModel)
{
  int iCan;
  type =cdtype;
  ExitGracefullyIf(to_index==DOESNT_EXIST,
                   "CmvCanopyDrip Constructor: invalid to compartment specified",BAD_DATA);
  iCan     =pModel->GetStateVarIndex(CANOPY);
  CHydroProcessABC::DynamicSpecifyConnections(1);//nConnections=1
  iFrom[0]=iCan;      iTo[0]=to_index;
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the standard destructor
//
CmvCanopyDrip::~CmvCanopyDrip(){}

//////////////////////////////////////////////////////////////////
/// \brief Verifies that canopy drip comes from a canopy unit
//
void     CmvCanopyDrip::Initialize()
{
  ExitGracefullyIf(pModel->GetStateVarType(iFrom[0])!=CANOPY,
                   "CmvCanopyDrip::Initialize:Canopy drip must come from canopy unit",BAD_DATA);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for canopy drip algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by canopy drip algorithm (size of aP[] and aPC[])
//
void CmvCanopyDrip::GetParticipatingParamList(string  *aP , class_type *aPC , int &nP) const
{
  if (type==CANDRIP_RUTTER)
  {
    nP=3;
    aP[0]="FOREST_COVERAGE"; aPC[0]=CLASS_LANDUSE; //JRCFLAG
    aP[1]="MAX_CAPACITY";    aPC[1]=CLASS_VEGETATION;
    aP[2]="STEMFLOW_FRAC";   aPC[2]=CLASS_VEGETATION;
  }
  else if (type==CANDRIP_SLOWDRAIN)
  {
    nP=3;
    aP[0]="DRIP_PROPORTION"; aPC[0]=CLASS_VEGETATION;
    aP[1]="MAX_CAPACITY";    aPC[1]=CLASS_VEGETATION;
    aP[2]="FOREST_COVERAGE"; aPC[2]=CLASS_LANDUSE; //JRCFLAG
  }
  else
  {
    ExitGracefully("CmvCanopyDrip::GetParticipatingParamList: undefined canopy drip algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
/// \remark 'To' compartment is user specified
///
/// \param cdtype [in] Model of canopy drip used
/// \param *aSV [out] Array of state variable types needed by canopy drip algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by canopy drip algorithm (size of aSV[] and aLev[] arrays)
//
void CmvCanopyDrip::GetParticipatingStateVarList(candrip_type cdtype,sv_type *aSV, int *aLev, int &nSV)
{
  nSV=1;
  aSV[0]=CANOPY;  aLev[0]=DOESNT_EXIST;
  //'to' compartment is user-specified
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rate of loss of water from canopy to (typically) ponded water [mm/d]
/// \details   if type==CANDRIP_RUTTER (rutter model), drip rate is calculated from storage overflow, as in Brook90 \cite Federer2010
///
/// \param *state_vars [in] Array of state variables for HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of loss from canopy to land  surface [mm/day]
//
void CmvCanopyDrip::GetRatesOfChange( const double      *state_vars,
                                      const CHydroUnit  *pHRU,
                                      const optStruct   &Options,
                                      const time_struct &tt,
                                      double            *rates) const
{

  if ((pHRU->GetHRUType()!=HRU_STANDARD) &&
      (pHRU->GetHRUType()!=HRU_WETLAND)){return;}

  double Fc,p;
  rates[0]=0.0;//default
  Fc =pHRU->GetSurfaceProps()->forest_coverage;
  if (Fc==0){return;}

  double stor=state_vars[iFrom[0]];
  double cap =pHRU->GetVegVarProps()->capacity;

  if (type==CANDRIP_RUTTER)//-------------------------------
  {
    p  =pHRU->GetVegetationProps()->stemflow_frac;

    if (pModel->GetStateVarIndex(TRUNK)<0){p=0.0;}//overrides if trunk not modeled

    //if storage is greater than capacity, then overflow occurs at rate d(S-C)/dt
    //this means storage cannot be exceeded for a full timestep

    rates[0]=(1.0-p)*threshMax((stor-Fc*cap)/Options.timestep,0.0,0.0);
  }
  else if (type==CANDRIP_SLOWDRAIN)//----------------------------
  {
    double drip=pHRU->GetVegetationProps()->drip_proportion;
    rates[0]=threshPositive((stor-Fc*cap)/Options.timestep)+//overflow
      threshMin(drip*(stor/Fc),stor/Fc/Options.timestep,0.0);//slow drip //THRESHOLD BEHAVIOR
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow cannot drain "from" compartment over timestep
///
/// \param *state_vars [in] Array of state variables for HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of loss from canopy to land  surface [mm/day]
//
void CmvCanopyDrip::ApplyConstraints( const double      *state_vars,
                                      const CHydroUnit  *pHRU,
                                      const optStruct   &Options,
                                      const time_struct &tt,
                                      double      *rates) const
{
  if ((pHRU->GetHRUType()!=HRU_STANDARD) &&
      (pHRU->GetHRUType()!=HRU_WETLAND)){return;}

  //cant remove more than is there
  rates[0]=threshMin(rates[0],state_vars[iFrom[0]]/Options.timestep,0.0);
}
