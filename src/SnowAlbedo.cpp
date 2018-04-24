/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2018 the Raven Development Team
  ----------------------------------------------------------------
  Snow albedo evolution routines
  ----------------------------------------------------------------*/
#include "HydroProcessABC.h"
#include "Albedo.h"

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the snow albedo evolution constructor
//
CmvSnowAlbedoEvolve::CmvSnowAlbedoEvolve(snowalb_type snalb_type):
  CHydroProcessABC(SNOW_ALBEDO_EVOLVE)
{
  type =snalb_type;

  int iSnowAlbedo; //used by all snow albedo routines
  iSnowAlbedo   =pModel->GetStateVarIndex(SNOW_ALBEDO);
  ExitGracefullyIf(iSnowAlbedo==DOESNT_EXIST,
                   "SNOW ALBEDO EVOLVE: snow albedo state variable required",BAD_DATA);

  if (type==SNOALB_UBCWM)
  {
    CHydroProcessABC::DynamicSpecifyConnections(1);
    iFrom[0]=iSnowAlbedo;       iTo[0]=iSnowAlbedo;//rates[0]: SNOW_ALBEDO->SNOW_ALBEDO
  }
  else{
    ExitGracefully("CmvSnowAlbedoEvolve::Constructor: undefined snow albedo algorithm",BAD_DATA);
  }
}

///////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvSnowAlbedoEvolve::~CmvSnowAlbedoEvolve(){}

///////////////////////////////////////////////////////////////////
/// \brief Function to initialize CmvSnowAlbedoEvolve objects
//
void CmvSnowAlbedoEvolve::Initialize(){}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for snow albedo algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by snow albedo algorithm (size of aP[] and aPC[])
//
void CmvSnowAlbedoEvolve::GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const
{
  if (type==SNOALB_UBCWM)
  {
    nP=6;
    aP[0]="MAX_SNOW_ALBEDO";     aPC[0]=CLASS_GLOBAL;
    aP[1]="MIN_SNOW_ALBEDO";     aPC[1]=CLASS_GLOBAL;
    aP[2]="UBC_ALBASE";          aPC[2]=CLASS_GLOBAL;
    aP[3]="UBC_ALBREC";          aPC[3]=CLASS_GLOBAL;
    aP[4]="UBC_MAX_CUM_MELT";    aPC[4]=CLASS_GLOBAL;
    aP[5]="UBC_ALBSNW";          aPC[5]=CLASS_GLOBAL;
  }
  else
  {
    ExitGracefully("CmvSnowAlbedoEvolve::GetParticipatingParamList: undefined snow albedo algorithm",BAD_DATA);
  }
}

///////////////////////////////////////////////////////////////////
/// \brief Gets participating state variable list
///
/// \param bal_type [in] Snow albedo algorithm type
/// \param *aSV [out] Array of state variable types modified by this process
/// \param *aLev [out] Array of indicies corresponding to state variable layers modified by this process. ( = DOESNT_EXIST for non-multiple layers)
/// \param &nSV [out] Number of state variables modified by this process (i.e. length of *aSV)
//
void CmvSnowAlbedoEvolve::GetParticipatingStateVarList(snowalb_type bal_type,
                                                       sv_type *aSV, int *aLev, int &nSV)
{
  nSV=1;
  aSV[0]=SNOW_ALBEDO;   aLev[0]=DOESNT_EXIST;

  if (bal_type==SNOALB_UBCWM)
  {
    nSV=3;
    aSV[1]=SNOW;            aLev[1]=DOESNT_EXIST;
    aSV[2]=CUM_SNOWMELT;    aLev[2]=DOESNT_EXIST;
  }
}

///////////////////////////////////////////////////////////////////
/// \brief Calculates rates of change in snow albedo over time
/// \remark SNOALB_UBC algorithm adapted from UBC Watershed Model
/// \copyright Michael Quick
///
/// \param *state_var [in] Array of state variables in specified HRU (size=CModel::_nStateVars)
/// \param *pHRU [in] Pointer to an instatiated HRU wherein the process is occuring
/// \param &Options [in] Global model options information
/// \param &tt [in] Reference to the instance in time to which the calculated rates of change pertain
/// \param *rates [out] Array (size=nConnects) which contains calculate rates of change of modified variables.
//
void CmvSnowAlbedoEvolve::GetRatesOfChange (const double                 *state_var,
                                            const CHydroUnit *pHRU,
                                            const optStruct      &Options,
                                            const time_struct &tt,
                                            double     *rates) const
{
  if (pHRU->GetHRUType()==HRU_LAKE){return;}

  if (type==SNOALB_UBCWM)
  {
    double snow,albedo,cum_melt,snowfall,old_albedo;

    UBC_snow_par PP=CGlobalParams::GetParams()->UBC_snow_params;
    double min_alb=CGlobalParams::GetParams()->min_snow_albedo;
    double max_alb=CGlobalParams::GetParams()->max_snow_albedo;

    snow    =state_var[pModel->GetStateVarIndex(SNOW)];
    albedo  =state_var[pModel->GetStateVarIndex(SNOW_ALBEDO)];
    cum_melt=state_var[pModel->GetStateVarIndex(CUM_SNOWMELT)];
    snowfall=pHRU->GetForcingFunctions()->precip*pHRU->GetForcingFunctions()->snow_frac;

    old_albedo=albedo;
    //snow albedo decay
    if (snow>REAL_SMALL || snowfall>REAL_SMALL)
    {
      lowerswap(albedo,max_alb);

      if (albedo>PP.ALBASE){
        static double albrec_factor = pow(PP.ALBREC,Options.timestep);
        albedo*= albrec_factor;
      }

      else if (albedo<PP.ALBASE){ //separate if in case the previous change goes below ALBASE
        albedo = PP.ALBASE*exp(-cum_melt/PP.MAX_CUM_MELT);
      }
      upperswap(albedo,min_alb);

      if (snowfall >0)//increase in snow albedo
      {
        albedo += min(snowfall*Options.timestep / PP.ALBSNW,1.0) * (max_alb - albedo);
      }
    }

    rates[0]=(albedo-old_albedo)/Options.timestep;

  }//end UBC evolve
  else{
    rates[0]=0.0;
  }
};

///////////////////////////////////////////////////////////////////
/// \brief Applies constraints - ensures that albedo is on interval [0..1]
///
/// \param *state_vars [in] Array of state variables in specified HRU (size=CModel::_nStateVars)
/// \param *pHRU [in] Pointer to an instantiated HRU
/// \param &Options [in] Global model options information
/// \param &t [in] Reference to a specific instance in time
/// \param *rates [out] Array (size=nConnects) which contains calculate rates of change of modified variables.
//
void CmvSnowAlbedoEvolve::ApplyConstraints( const double      *state_vars,
                                            const CHydroUnit  *pHRU,
                                            const optStruct       &Options,
                                            const time_struct &t,
                                            double      *rates) const
{
  if (pHRU->GetHRUType()==HRU_LAKE){return;}

  double albedo  =state_vars[pModel->GetStateVarIndex(SNOW_ALBEDO)];
  lowerswap(rates[0],(1.0-albedo)/Options.timestep);//albedo<-=1
  upperswap(rates[0],     -albedo/Options.timestep);//albedo>=0
}
