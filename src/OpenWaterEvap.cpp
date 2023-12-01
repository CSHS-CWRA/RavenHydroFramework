/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team

  Open Water Evap
  Lake Evap
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "OpenWaterEvap.h"
#include "Model.h"

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the Open water evaporation constructor
/// \param owtype [in] Model of open water evaporation used
//
CmvOWEvaporation::CmvOWEvaporation(owevap_type owtype,
                                   const int   i_from,
                                   CModelABC   *pModel)
  :CHydroProcessABC(OPEN_WATER_EVAPORATION, pModel)
{
  _type =owtype;

  if(_type==OPEN_WATER_UWFS) {
    CHydroProcessABC::DynamicSpecifyConnections(3);
    iFrom[0]=i_from;
    iTo  [0]=pModel->GetStateVarIndex(ATMOSPHERE);                          //rates[0]: PONDED_WATER/DEPRESSION->ATMOSPHERE
    iFrom[1]=pModel->GetStateVarIndex(MIN_DEP_DEFICIT);  iTo[1]=iFrom[1];   //rates[1]:
    iFrom[2]=pModel->GetStateVarIndex(AET);              iTo[2]=iFrom[2];   //rates[2]: AET->AET
  }
  else { //Default
    CHydroProcessABC::DynamicSpecifyConnections(2);//nConnections=2
    iFrom[0]=i_from;                          iTo[0]=pModel->GetStateVarIndex(ATMOSPHERE);    //rates[0]: PONDED_WATER/DEPRESSION->ATMOSPHERE
    iFrom[1]=pModel->GetStateVarIndex(AET);   iTo[1]=iFrom[1];                                  //rates[1]: AET->AET
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvOWEvaporation::~CmvOWEvaporation(){}

//////////////////////////////////////////////////////////////////
/// \brief Verify that iFrom[] - iTo[] connection is DEPRESSION-ATMOSPHERE
//
void CmvOWEvaporation::Initialize()
{
  sv_type typ=pModel->GetStateVarType(iFrom[0]);
  ExitGracefullyIf((typ!=DEPRESSION) && (typ!=PONDED_WATER) && (typ!=SURFACE_WATER),
                   "CmvOWEvaporation::Initialize:Open Water evaporation must come from depression, ponded water, or surface water",BAD_DATA);
  ExitGracefullyIf((_type==OPEN_WATER_UWFS) && (typ!=DEPRESSION),
                   "CmvOWEvaporation::Initialize:OPEN_WATER_UWFS must come from depression storage only",BAD_DATA);
  ExitGracefullyIf(pModel->GetStateVarType(iTo[0])!=ATMOSPHERE,
                   "CmvOWEvaporation::Initialize:Open Water evaporation must go to atmosphere",BAD_DATA);
}
//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for baseflow algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by baseflow algorithm (size of aP[] and aPC[])
//
void CmvOWEvaporation::GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const
{
  nP=1;
  aP[0]="OW_PET_CORR";   aPC[0]=CLASS_LANDUSE;

  if(_type==OPEN_WATER_RIPARIAN) {
    nP=2;
    aP[1]="STREAM_FRACTION"; aPC[1]=CLASS_LANDUSE;
  }
  else if(_type==OPEN_WATER_UWFS) {
    nP=5;
    aP[0]="MAX_DEP_AREA_FRAC";      aPC[0]=CLASS_LANDUSE;
    aP[1]="DEP_MAX";                aPC[1]=CLASS_LANDUSE;
    aP[2]="UWFS_B";                 aPC[2]=CLASS_LANDUSE;
    aP[3]="UWFS_BETAMIN";           aPC[3]=CLASS_LANDUSE;
    aP[4] = "OW_PET_CORR";          aPC[4] = CLASS_LANDUSE;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Sets reference to state variable types needed by evaporation algorithm
///
/// \param owtype [in] Model of open water evaporation used
/// \param *aSV [out] Array of state variable types needed by ow evaporation algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by ow evaporation algorithm (size of aSV[] and aLev[] arrays)
//
void CmvOWEvaporation::GetParticipatingStateVarList(owevap_type owtype, sv_type *aSV, int *aLev, int &nSV)
{
  nSV=2;
  aSV[0]=ATMOSPHERE;  aLev[0]=DOESNT_EXIST;
  aSV[1]=AET;         aLev[1]=DOESNT_EXIST;
  //other state var is specified in constructor
  if(owtype==OPEN_WATER_UWFS)
  {
    nSV=3;
    aSV[0] = ATMOSPHERE;      aLev[0]=DOESNT_EXIST;
    aSV[1] = AET;             aLev[1]=DOESNT_EXIST;
    aSV[2]=MIN_DEP_DEFICIT;   aLev[2]=DOESNT_EXIST;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rate of loss of water from open water to atmosphere [mm/d]
///
/// \param *state_vars [in] Current array of state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time strucutre
/// \param *rates [out] rates[0] is rate of loss from open water to atmosphere [mm/day]
//
void CmvOWEvaporation::GetRatesOfChange(const double* state_vars,
                                        const CHydroUnit* pHRU,
                                        const optStruct& Options,
                                        const time_struct& tt,
                                        double* rates) const
{
    double OWPET;
    OWPET = pHRU->GetForcingFunctions()->OW_PET;            //open water PET rate [mm/d]

    if(pHRU->IsLinkedToReservoir()) { return; }//reservoir-linked HRUs handle ET via reservoir MB

    if (!Options.suppressCompetitiveET) {
        //competitive ET - reduce PET by AET
        OWPET -= (state_vars[pModel->GetStateVarIndex(AET)] / Options.timestep);
        OWPET = max(OWPET, 0.0);
    }

    //-------------------------------------------------------------------
    if (_type == OPEN_WATER_EVAP)
    {
        rates[0] = pHRU->GetSurfaceProps()->ow_PET_corr * OWPET;
    }
    //-------------------------------------------------------------------
    else if (_type == OPEN_WATER_RIPARIAN)
    {
        double Fstream = pHRU->GetSurfaceProps()->stream_fraction;
        rates[0] = pHRU->GetSurfaceProps()->ow_PET_corr * Fstream * OWPET;
    }
    //-------------------------------------------------------------------
    else if (_type == OPEN_WATER_UWFS)
    {
      double Dmin = state_vars[iFrom[1]];        //=Dmin [mm] if positive,=-Pf if negative

      rates[0] = pHRU->GetSurfaceProps()->max_dep_area_frac * pHRU->GetSurfaceProps()->ow_PET_corr * OWPET;//DEPRESSION->ATMOS
      if (state_vars[iFrom[1]] > 0.0) // if Dmin is positive changes in deficit will be sames as rate [0]
      {
        rates[1] = rates[0];
      }
      else  //for the case that Dmin is zero or negative Dminnew (after evaporation) will be equal to rate [0] and changes in dmin will be calculated as follows
      {
        double Dminnew;
        Dminnew = rates[0] * Options.timestep;
        rates[1] = (Dminnew - state_vars[iFrom[1]]) / Options.timestep;
      }
    }
    //-------------------------------------------------------------------

    rates[_nConnections-1]=rates[0]; //handles AET
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
///
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &t [in] Current model time structure
/// \param *rates [out] rates[0] is rate of loss from open water to atmosphere [mm/day]
//
void  CmvOWEvaporation::ApplyConstraints( const double            *state_vars,
                                          const CHydroUnit  *pHRU,
                                          const optStruct   &Options,
                                          const time_struct &t,
                                          double      *rates) const
{
  if (rates[0]<0)             {rates[0]=0.0;}//positivity constraint

  //can't remove more than is there (with exception of surface water in reach or lake HRU)
  if(iFrom[0]!=pModel->GetStateVarIndex(SURFACE_WATER)) {
    rates[0]=threshMin(rates[0],state_vars[iFrom[0]]/Options.timestep,0.0);
    if(state_vars[iFrom[0]]<=0) { rates[0]=0.0; }//reality check
  }
  rates[_nConnections-1]=rates[0]; //correct AET accordingly
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the lake evaporation constructor
/// \param lktype [in] Model of open water evaporation used
/// \param fromIndex [in] Index of lake from which water evaporates
//
CmvLakeEvaporation::CmvLakeEvaporation(lakeevap_type lktype,
                                       const int     fromIndex,
                                       CModelABC     *pModel)
  :CHydroProcessABC(LAKE_EVAPORATION, pModel)
{
  type =lktype;
  ExitGracefullyIf(fromIndex==DOESNT_EXIST,
                   "CmvLakeEvaporation Constructor: invalid 'from' compartment specified",BAD_DATA);

  CHydroProcessABC::DynamicSpecifyConnections(2);//nConnections=1
  iFrom[0]=fromIndex;
  iTo  [0]=pModel->GetStateVarIndex(ATMOSPHERE);;     //rates[0]: LAKE->ATMOSPHERE
  iFrom[1]=pModel->GetStateVarIndex(AET);
  iTo  [1]=pModel->GetStateVarIndex(AET);             //rates[0]: AET->AET
}
//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvLakeEvaporation::~CmvLakeEvaporation(){}

//////////////////////////////////////////////////////////////////
/// \brief Verify that process moves water to atmosphere
//
void CmvLakeEvaporation::Initialize()
{
  ExitGracefullyIf(pModel->GetStateVarType(iTo[0])!=ATMOSPHERE,
                   "CmvLakeEvaporation::Initialize:Open Water evaporation must go to atmosphere",BAD_DATA);
}
//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for baseflow algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by baseflow algorithm (size of aP[] and aPC[])
//
void CmvLakeEvaporation::GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const
{
  if(type==LAKE_EVAP_BASIC)//-------------------------------------
  {
    nP=1;
    aP[0]="LAKE_PET_CORR"; aPC[0]=CLASS_LANDUSE;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Sets reference to state variable types needed by evaporation algorithm
/// \note "From" compartment specified by :LakeStorage command
///
/// \param lktype [in] Model of open water evaporation used
/// \param *aSV [out] Array of state variable types needed by ow evaporation algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by lake evaporation algorithm (size of aSV[] and aLev[] arrays)
//
void CmvLakeEvaporation::GetParticipatingStateVarList(lakeevap_type lktype, sv_type *aSV, int *aLev, int &nSV)
{
  nSV=2;
  aSV[0]=ATMOSPHERE;  aLev[0]=DOESNT_EXIST;
  aSV[1]=AET;         aLev[1]=DOESNT_EXIST;
  //'from' compartment specified by :LakeStorage command
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rate of loss of water from lake to atmosphere [mm/d]
/// \details  if type==LAKE_EVAP_BASIC, evaporation is calculated using fraction of PET method
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time structure
/// \param *rates [out] rates[0] is rate of loss from lake to atmosphere [mm/d]
//
void CmvLakeEvaporation::GetRatesOfChange(const double                  *state_vars,
                                          const CHydroUnit      *pHRU,
                                          const optStruct         &Options,
                                          const time_struct &tt,
                                          double      *rates) const
{
  if ((!pHRU->IsLake()) && (pModel->GetLakeStorageIndex()!=iFrom[0])){
    return;
  } //only works for lakes OR special storage units designated using :LakeStorage (the latter is to support HBV)

  if(pHRU->IsLinkedToReservoir()) { return; }//reservoir-linked HRUs handle ET via reservoir MB

  double OWPET;
  OWPET = pHRU->GetForcingFunctions()->OW_PET;          //calls PET rate [mm/d]

  if (type==LAKE_EVAP_BASIC)//-------------------------------------
  {
    rates[0]= pHRU->GetSurfaceProps()->lake_PET_corr*OWPET;
  }
  rates[_nConnections-1]=rates[0]; //handles AET
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow ccannot drain lake over timestep
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time structure
/// \param *rates [out] rates[0] is rate of loss from lake to atmosphere [mm/d]
//
void  CmvLakeEvaporation::ApplyConstraints( const double                *state_vars,
                                            const CHydroUnit  *pHRU,
                                            const optStruct   &Options,
                                            const time_struct &t,
                                            double      *rates) const
{
  if (state_vars[iFrom[0]]<=0){rates[0]=0.0; }//reality check
  if (rates[0]<0)             {rates[0]=0.0; }//positivity constraint

  //can't remove more than is there
  rates[0]=threshMin(rates[0],state_vars[iFrom[0]]/Options.timestep,0.0);

  rates[_nConnections-1]=rates[0]; //correct AET accordingly
}
