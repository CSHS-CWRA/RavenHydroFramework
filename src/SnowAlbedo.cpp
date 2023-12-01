/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2018 the Raven Development Team
  ----------------------------------------------------------------
  Snow albedo evolution routines
  ----------------------------------------------------------------*/
#include "HydroProcessABC.h"
#include "Albedo.h"
#include "Model.h"

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the snow albedo evolution constructor
//
CmvSnowAlbedoEvolve::CmvSnowAlbedoEvolve(snowalb_type snalb_type, CModelABC *pModel)
  :CHydroProcessABC(SNOW_ALBEDO_EVOLVE, pModel)
{
  type =snalb_type;

  int iSnowAlbedo; //used by all snow albedo routines
  iSnowAlbedo   =pModel->GetStateVarIndex(SNOW_ALBEDO);
  ExitGracefullyIf(iSnowAlbedo==DOESNT_EXIST,
                   "SNOW ALBEDO EVOLVE: snow albedo state variable required",BAD_DATA);

  if((type==SNOALB_UBCWM) || (type==SNOALB_CRHM_ESSERY))
  {
    CHydroProcessABC::DynamicSpecifyConnections(1);
    iFrom[0]=iSnowAlbedo;       iTo[0]=iSnowAlbedo;//rates[0]: SNOW_ALBEDO->SNOW_ALBEDO
  }
  else if(type==SNOALB_BAKER)
  {
    int iSnowAge=pModel->GetStateVarIndex(SNOW_AGE);
    CHydroProcessABC::DynamicSpecifyConnections(2);
    iFrom[0]=iSnowAlbedo;       iTo[0]=iSnowAlbedo;//rates[0]: SNOW_ALBEDO->SNOW_ALBEDO
    iFrom[1]=iSnowAge;          iTo[1]=iSnowAge;   //rates[1]: SNOW_AGE->SNOW_AGE
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
  else if(type==SNOALB_CRHM_ESSERY)
  {
    nP=6;
    aP[0]="ALB_DECAY_COLD";     aPC[0]=CLASS_GLOBAL;
    aP[1]="ALB_DECAY_MELT";     aPC[1]=CLASS_GLOBAL;
    aP[2]="MIN_SNOW_ALBEDO";    aPC[2]=CLASS_GLOBAL;
    aP[3]="MAX_SNOW_ALBEDO";    aPC[3]=CLASS_GLOBAL;
    aP[4]="BARE_GROUND_ALBEDO"; aPC[4]=CLASS_GLOBAL;
    aP[5]="SNOWFALL_ALBTHRESH"; aPC[5]=CLASS_GLOBAL;
  }
  else if(type==SNOALB_BAKER)
  {
    nP=2;
    aP[0]="BARE_GROUND_ALBEDO"; aPC[0]=CLASS_GLOBAL;
    aP[1]="SNOWFALL_ALBTHRESH"; aPC[1]=CLASS_GLOBAL;
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
  else if(bal_type==SNOALB_CRHM_ESSERY)
  {
    nSV=3;
    aSV[1]=SNOW;            aLev[1]=DOESNT_EXIST;
    aSV[2]=SNOW_TEMP;       aLev[2]=DOESNT_EXIST;
  }
  else if(bal_type==SNOALB_BAKER)
  {
    nSV=3;
    aSV[1]=SNOW;            aLev[1]=DOESNT_EXIST;
    aSV[2]=SNOW_AGE;        aLev[2]=DOESNT_EXIST;
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
void CmvSnowAlbedoEvolve::GetRatesOfChange (const double *state_vars,
                                            const CHydroUnit *pHRU,
                                            const optStruct &Options,
                                            const time_struct &tt,
                                            double *rates) const
{

  const force_struct *F=pHRU->GetForcingFunctions();

  if (type==SNOALB_UBCWM)//===================================================================
  {
    double snow,albedo,cum_melt,snowfall,old_albedo;

    UBC_snow_par PP = pModel->GetGlobalParams()->GetParams()->UBC_snow_params;
    double min_alb  = pModel->GetGlobalParams()->GetParams()->min_snow_albedo;
    double max_alb  = pModel->GetGlobalParams()->GetParams()->max_snow_albedo;

    snow    =state_vars[pModel->GetStateVarIndex(SNOW)];
    albedo  =state_vars[pModel->GetStateVarIndex(SNOW_ALBEDO)];
    cum_melt=state_vars[pModel->GetStateVarIndex(CUM_SNOWMELT)];
    snowfall=F->precip*F->snow_frac;

    old_albedo=albedo;
    //snow albedo decay
    if (snow>REAL_SMALL || snowfall>REAL_SMALL)
    {
      lowerswap(albedo,max_alb);

      if (albedo>PP.ALBASE){
        double albrec_factor = pow(PP.ALBREC,Options.timestep);
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
  else if(type==SNOALB_CRHM_ESSERY)//=========================================================
  { //ported from CRHM (Pomeroy, 2007) routine ClassalbedoRichard::run
    double tstep=Options.timestep;

    double a1          = pModel->GetGlobalParams()->GetParams()->alb_decay_cold;    //Albedo decay time constant for cold snow (~0.008/d)
    double a2          = pModel->GetGlobalParams()->GetParams()->alb_decay_melt;    //Albedo decay time constant for melting snow (~0.12/d)
    double albmin      = pModel->GetGlobalParams()->GetParams()->min_snow_albedo;
    double albmax      = pModel->GetGlobalParams()->GetParams()->max_snow_albedo;
    double albbare     = pModel->GetGlobalParams()->GetParams()->bare_ground_albedo;
    double snowfall_th = pModel->GetGlobalParams()->GetParams()->snowfall_albthresh;//Minimum snowfall to refresh snow albedo [mm/d] (~10 mm/d)

    double albedo     =state_vars[pModel->GetStateVarIndex(SNOW_ALBEDO)];
    double SWE        =state_vars[pModel->GetStateVarIndex(SNOW)];
    double snow_temp  =state_vars[pModel->GetStateVarIndex(SNOW_TEMP)];
    double old_albedo=albedo;

    if (SWE > 1.0) //1mm threshold
    {
      if(snow_temp < 0) {albedo-=a1*tstep; } // cold snow - 0th order decay
      else              {albedo = (albedo - albmin)*exp(-a2*tstep) + albmin; } // melting snow -1st order decay

      albedo += (albmax - albedo)*((F->snow_frac*F->precip)/snowfall_th);

      upperswap(albedo,albmin);
      lowerswap(albedo,albmax);
    }
    else {
      albedo = albbare;
    }
    rates[0]=(albedo-old_albedo)/tstep;
  }
  else if(type==SNOALB_BAKER)//===============================================================
  {
    //ported from CRHM (Pomeroy, 2007) -ClassalbedoBaker::run
    // albedo decay formulation after
    // Baker, D.G., Ruschy, D.L., Wall, D.B., 1990. The albedo decay of prairie snows. J. Appl. Meteor. 29 _2, 179-187
    double tstep=Options.timestep;

    double SWE      = state_vars[pModel->GetStateVarIndex(SNOW)];
    double albedo   = state_vars[pModel->GetStateVarIndex(SNOW_ALBEDO)];
    double snow_age = state_vars[pModel->GetStateVarIndex(SNOW_AGE)];

    double snowfall_th = pModel->GetGlobalParams()->GetParams()->snowfall_albthresh;//Minimum snowfall to refresh snow albedo [mm/d] (~10 mm/d)
    double albbare     = pModel->GetGlobalParams()->GetParams()->bare_ground_albedo;

    double old_albedo =albedo;
    double old_snowage=snow_age;//time [d] since albedo refresh

    if (SWE > 1.0)
    {
      if ((F->snow_frac*F->precip)>snowfall_th){snow_age=0.0;}
      else                                     {snow_age+=tstep; };
      albedo = 0.9 - 0.0473 * pow(snow_age,0.1);
      upperswap(albedo,albbare);
    }
    else
    {
      snow_age=0.0;
      albedo = albbare;
    }
    rates[0]=(albedo  -old_albedo )/tstep;
    rates[1]=(snow_age-old_snowage)/tstep;
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

  double albedo  =state_vars[pModel->GetStateVarIndex(SNOW_ALBEDO)];
  lowerswap(rates[0],(1.0-albedo)/Options.timestep);//albedo<-=1
  upperswap(rates[0],     -albedo/Options.timestep);//albedo>=0
}
