/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ------------------------------------------------------------------
  Crop growth & evolution routines
  ----------------------------------------------------------------*/
#include "HydroProcessABC.h"
#include "CropGrowth.h"

//////////////////////////////////////////////////////////////////
/// \brief Implentation of the CropHeatUnitEvolve constructor
/// \param CHU_type [in] Method of modelling crop heat units selected
//
CmvCropHeatUnitEvolve::CmvCropHeatUnitEvolve(CHUevolve_type CHU_type):
  CHydroProcessABC(CROP_HEAT_UNIT_EVOLVE)
{
  type =CHU_type;

  int iCropHeat; //used by all crop heat routines
  iCropHeat   =pModel->GetStateVarIndex(CROP_HEAT_UNITS);
  ExitGracefullyIf(iCropHeat==DOESNT_EXIST,
                   "CROP HEAT EVOLVE: crop heat unit state variable required",BAD_DATA);

  if (type==CHU_ONTARIO)
  {
    CHydroProcessABC::DynamicSpecifyConnections(1);
    iFrom[0]=iCropHeat;       iTo[0]=iCropHeat;//rates[0]: CROP_HEAT_UNITS->CROP_HEAT_UNITS
  }
  else{
    ExitGracefully("CmvCropHeatUnitEvolve::Constructor: undefined crop heat evolution algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default CropHeatUnitEvolve destructor
//
CmvCropHeatUnitEvolve::~CmvCropHeatUnitEvolve(){}

//////////////////////////////////////////////////////////////////
/// \brief Initializes CropHeatUnitEvolve object
//
void CmvCropHeatUnitEvolve::Initialize(){}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for CHUevolve algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by CHUevolve algorithm (size of aP[] and aPC[])
//
void CmvCropHeatUnitEvolve::GetParticipatingParamList(string  *aP , class_type *aPC , int &nP) const
{
  if (type==CHU_ONTARIO)
  {
    nP=0;
  }
  else
  {
    ExitGracefully("CmvCropHeatUnitEvolve::GetParticipatingParamList: undefined CHUevolve algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
///
/// \param CHU_type [in] Method of modeling CHU evolution
/// \param *aSV [out] Array of state variable types needed by CHUevolve algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by CHUevolve algorithm (size of aSV[] and aLev[] arrays)
//
void CmvCropHeatUnitEvolve::GetParticipatingStateVarList(CHUevolve_type CHU_type,
                                                         sv_type *aSV, int *aLev, int &nSV)
{
  nSV=1;
  aSV[0]=CROP_HEAT_UNITS;   aLev[0]=DOESNT_EXIST;
}

//////////////////////////////////////////////////////////////////
/// \brief For modeling changes in crop heat units over time
/// \details standard Ontario CHU evolution approach (University of Guelph)
/// \ref UoG
///
/// \param *state_var [in] array of current state variables
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Rate of change in crop heat units and related state variables
//
void CmvCropHeatUnitEvolve::GetRatesOfChange (const double               *state_var,
                                              const CHydroUnit *pHRU,
                                              const optStruct      &Options,
                                              const time_struct &tt,
                                              double     *rates) const
{
  if (pHRU->GetHRUType()!=HRU_STANDARD){return;}//Lakes & Glaciers

  double CHU,old_CHU;
  CHU    =state_var[pModel->GetStateVarIndex(CROP_HEAT_UNITS)];
  old_CHU=CHU;

  if (type==CHU_ONTARIO)
  {

    double CHU_day, CHU_night;

    double temp_ave=pHRU->GetForcingFunctions()->temp_daily_ave;
    double temp_max=pHRU->GetForcingFunctions()->temp_daily_max;
    double temp_min=pHRU->GetForcingFunctions()->temp_daily_min;
    bool growing_season=(CHU>-0.5);

    if (!growing_season)
    {
      if (temp_ave>12.8) // \todo [fix hack] This is a stopgap, should use 3_day min temp (non-existent F.temp_3daymin??)
      {
        CHU=0.0; //CHU reset to zero at beginning of growing season
      }
    }
    else if (growing_season)
    {
      CHU_day  =max(3.33*(temp_max-10)-0.084*pow(temp_max-10.0,2.0),0.0);
      CHU_night=max(1.8*(temp_min-4.4),0.0);
      CHU+=0.5*(CHU_day+CHU_night);

      //end of growing season: CHU set to -1
      if (temp_ave<-2.0){CHU=-1.0;}
      if (tt.month>9   ){CHU=-1.0;}//>september 30th
    }
    rates[0]=(CHU-old_CHU)/Options.timestep;

  }//end ONTARIO evolve
  else{
    rates[0]=0.0;
  }
};
//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that  CHU >= 0
///
/// \param *state_vars [in] Reference to pertinent state variables
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &t [in] Current model time
/// \param *rates [out] Rate of change in crop heat units
//
void CmvCropHeatUnitEvolve::ApplyConstraints( const double      *state_vars,
                                              const CHydroUnit  *pHRU,
                                              const optStruct       &Options,
                                              const time_struct &t,
                                              double      *rates) const
{
  if (pHRU->GetHRUType()!=HRU_STANDARD){return;}//Lakes & Glaciers
  //no constraints, allow to be =-1 prior to start of growing season
}
