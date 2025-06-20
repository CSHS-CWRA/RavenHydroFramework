/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2023 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Properties.h"
#include "GlobalParams.h"

/*****************************************************************
   Constructor / Destructor
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the global parameter constructor
//
CGlobalParams::CGlobalParams()
{
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
//
CGlobalParams::~CGlobalParams()
{
  if (DESTRUCTOR_DEBUG){cout<<"  DELETING GLOBAL PARAMETER CLASS "<<endl;}
}

///////////////////////////////////////////////////////////////////
/// \brief Returns pointer to global parameters
/// \return pointer to global parameters
//
const global_struct *CGlobalParams::GetParams() {return &G;}


//////////////////////////////////////////////////////////////////
/// \brief Summarize global parameter information to screen
//
void CGlobalParams::SummarizeToScreen()
{
  //Nothing, for now
}

//////////////////////////////////////////////////////////////////
/// \brief Automatically calculates global parameters
/// \details Sets terrain properties based upon simple terrain parameters
///     Input [Ttmp] has been read from .rvp file - if parameter ==
///     AUTO_COMPLETE, then empirical relationships are used to estimate
///     parameters
///
/// \param &Gtmp [in] Input global class parameters (read from .rvp file)
/// \param &Gtemplate [in] Default global class parameters
//
void CGlobalParams::AutoCalculateGlobalParams(const global_struct &Gtmp, const global_struct &Gtemplate)
{
  bool autocalc=false;
  string warn;
  bool chatty=true;
  //Standard global parameters
  //----------------------------------------------------------------------------

  /*
  FUTURE GENERALIZED IMPLEMENTATION ::
  for(int i=0; i<_nGPParams; i++){
    if(_isAutocalculable[i]){
      autocalc=SetCalculableValue(G.max_snow_albedo,Gtmp.max_snow_albedo,Gtemplate.max_snow_albedo);
      if(autocalc){
        bool warning=true;
        if     (_aGPNames[i]=="MAX_SNOW_ALBEDO"){ _aGPvalues[i]=0.95; } //[C/km] (UBC default)
        else if(_aGPNames[i]=="SNOW_SWI"       ){ _aGPvalues[i]=0.05; }
        else {}//...


        warn="The required global parameter "+_aGPNames[i]+" was autogenerated with value "+to_string(_aGPvalue[i]);
        if(chatty && warning){ WriteAdvisory(warn,false); }
      }
      else{
        SetSpecifiedValue(_aGPvalues[i],aGPTempValues[i],aGPTemplateValues[i],false,_aGPNames);
      }
    }*/


  autocalc=SetCalculableValue(G.snow_SWI,Gtmp.snow_SWI,Gtemplate.snow_SWI);
  if (autocalc)
  {
    G.snow_SWI=0.05;
    warn="The required global parameter SNOW_SWI (:IrreducibleSnowSaturation) was autogenerated with value "+to_string(G.snow_SWI);
    if (chatty){WriteAdvisory(warn,false);}
  }
  autocalc=SetCalculableValue(G.snow_SWI_min,Gtmp.snow_SWI_min,Gtemplate.snow_SWI_min);
  if(autocalc)
  {
    G.snow_SWI_min=0.05;
    warn="The required global parameter SNOW_SWI_MIN was autogenerated with value "+to_string(G.snow_SWI_min);
    if(chatty) { WriteAdvisory(warn,false); }
  }
  autocalc=SetCalculableValue(G.snow_SWI_max,Gtmp.snow_SWI_max,Gtemplate.snow_SWI_max);
  if(autocalc)
  {
    G.snow_SWI_max=0.05;
    warn="The required global parameter SNOW_SWI_MAX was autogenerated with value "+to_string(G.snow_SWI_min);
    if(chatty) { WriteAdvisory(warn,false); }
  }

  autocalc=SetCalculableValue(G.snow_temperature,Gtmp.snow_temperature,Gtemplate.snow_temperature);
  if (autocalc)
  {
    G.snow_temperature=FREEZING_TEMP;//default=freezing
    warn="The required global parameter SNOW_TEMPERATURE was autogenerated with value "+to_string(G.snow_temperature);
    if (chatty){WriteAdvisory(warn,false);}
  }
  autocalc=SetCalculableValue(G.snow_roughness,Gtmp.snow_roughness,Gtemplate.snow_roughness);
  if (autocalc)
  {
    G.snow_roughness=0.0002;//m
    warn="The required global parameter SNOW_ROUGHNESS was autogenerated with value "+to_string(G.snow_roughness)+" [m]";
    if (chatty){WriteAdvisory(warn,false);}
  }
  autocalc=SetCalculableValue(G.rainsnow_temp,Gtmp.rainsnow_temp,Gtemplate.rainsnow_temp);
  if (autocalc)
  {
    G.rainsnow_temp=0.15;
    warn="The required global parameter RAINSNOW_TEMP was autogenerated with value "+to_string(G.rainsnow_temp);
    if (chatty){WriteAdvisory(warn,false);}
  }
  autocalc=SetCalculableValue(G.rainsnow_delta,Gtmp.rainsnow_delta,Gtemplate.rainsnow_delta);
  if (autocalc)
  {
    G.rainsnow_delta=2.0;
    warn="The required global parameter RAINSNOW_DELTA was autogenerated with value "+to_string(G.rainsnow_delta);
    if (chatty){WriteAdvisory(warn,false);}
  }

  autocalc=SetCalculableValue(G.max_SWE_surface,Gtmp.max_SWE_surface,Gtemplate.max_SWE_surface);
  if (autocalc)
  {
    G.max_SWE_surface=100.0;
    warn="The required global parameter MAX_SWE_SURFACE was autogenerated with value "+to_string(G.max_SWE_surface);
    if (chatty){WriteAdvisory(warn,false);}
  }
  autocalc=SetCalculableValue(G.TOC_multiplier,Gtmp.TOC_multiplier,Gtemplate.TOC_multiplier);
  if (autocalc)
  {
    G.TOC_multiplier=1.0;//default=nothing (no warning needed)
  }
  autocalc=SetCalculableValue(G.TIME_TO_PEAK_multiplier,Gtmp.TIME_TO_PEAK_multiplier,Gtemplate.TIME_TO_PEAK_multiplier);
  if (autocalc)
  {
    G.TIME_TO_PEAK_multiplier=1.0;//default=nothing (no warning needed)
  }
  autocalc=SetCalculableValue(G.GAMMA_SHAPE_multiplier,Gtmp.GAMMA_SHAPE_multiplier,Gtemplate.GAMMA_SHAPE_multiplier);
  if (autocalc)
  {
    G.GAMMA_SHAPE_multiplier=1.0;//default=nothing (no warning needed)
  }
  autocalc=SetCalculableValue(G.GAMMA_SCALE_multiplier,Gtmp.GAMMA_SCALE_multiplier,Gtemplate.GAMMA_SCALE_multiplier);
  if (autocalc)
  {
    G.GAMMA_SCALE_multiplier=1.0;//default=nothing (no warning needed)
  }
  autocalc=SetCalculableValue(G.reference_flow_mult,Gtmp.reference_flow_mult,Gtemplate.reference_flow_mult);
  if (autocalc)
  {
    G.reference_flow_mult=10;
    warn="The required global parameter REFERENCE_FLOW_MULT is set to the value "+to_string(G.reference_flow_mult);
    if (chatty){WriteAdvisory(warn,false);}
  }
  autocalc=SetCalculableValue(G.adiabatic_lapse,Gtmp.adiabatic_lapse,Gtemplate.adiabatic_lapse);
  if (autocalc)
  {
    G.adiabatic_lapse=6.4; //[C/km] (UBC default)
    warn="The required global parameter ADIABATIC_LAPSE was autogenerated with value "+to_string(G.adiabatic_lapse);
    if (chatty){WriteAdvisory(warn,false);}
  }

  autocalc=SetCalculableValue(G.wet_adiabatic_lapse,Gtmp.wet_adiabatic_lapse,Gtemplate.wet_adiabatic_lapse);
  if (autocalc)
  {
    G.wet_adiabatic_lapse=6.4; //[C/km] (UBC default)
    warn="The required global parameter WET_ADIABATIC_LAPSE was autogenerated with value "+to_string(G.wet_adiabatic_lapse);
    if (chatty){WriteAdvisory(warn,false);}
  }

  autocalc=SetCalculableValue(G.precip_lapse,Gtmp.precip_lapse,Gtemplate.precip_lapse);
  if (autocalc)
  {
    G.precip_lapse=0.0; //[mm/d/km] (default - no correction)
    //no warning - default is no correction
  }
  for (int i=0;i<12;i++){
    autocalc=SetCalculableValue(G.UBC_s_corr[i],Gtmp.UBC_s_corr[i],Gtemplate.UBC_s_corr[i]);
    if (autocalc)
    {
      G.UBC_s_corr[i]=1.0;//reasonable default- no correction
      //no warning - default is no correction
    }
    autocalc=SetCalculableValue(G.UBC_n_corr[i],Gtmp.UBC_n_corr[i],Gtemplate.UBC_n_corr[i]);
    if (autocalc)
    {
      G.UBC_n_corr[i]=1.0;//reasonable default- no correction
      //no warning - default is no correction
    }
  }

  autocalc=SetCalculableValue(G.max_snow_albedo,Gtmp.max_snow_albedo,Gtemplate.max_snow_albedo);
  if (autocalc)
  {
    G.max_snow_albedo=0.95; //[C/km] (UBC default)
    warn="The required global parameter MAX_SNOW_ALBEDO was autogenerated with value "+to_string(G.max_snow_albedo);
    if (chatty){WriteAdvisory(warn,false);}
  }

  autocalc=SetCalculableValue(G.min_snow_albedo,Gtmp.min_snow_albedo,Gtemplate.min_snow_albedo);
  if (autocalc)
  {
    G.min_snow_albedo=0.3; //[C/km] (UBC default)
    warn="The required global parameter MIN_SNOW_ALBEDO was autogenerated with value "+to_string(G.min_snow_albedo);
    if (chatty){WriteAdvisory(warn,false);}
  }

  autocalc=SetCalculableValue(G.alb_decay_cold,Gtmp.alb_decay_cold,Gtemplate.alb_decay_cold);
  if (autocalc)
  {
    G.alb_decay_cold=0.008; //[1/d] (CRHM default)
    warn="The required global parameter ALB_DECAY_COLD was autogenerated with value "+to_string(G.alb_decay_cold);
    if (chatty){WriteAdvisory(warn,false);}
  }
  autocalc=SetCalculableValue(G.alb_decay_melt,Gtmp.alb_decay_melt,Gtemplate.alb_decay_melt);
  if (autocalc)
  {
    G.alb_decay_melt=0.12; //[1/d] (CRHM default)
    warn="The required global parameter ALB_DECAY_MELT was autogenerated with value "+to_string(G.alb_decay_melt);
    if (chatty){WriteAdvisory(warn,false);}
  }
  autocalc=SetCalculableValue(G.bare_ground_albedo,Gtmp.bare_ground_albedo,Gtemplate.bare_ground_albedo);
  if (autocalc)
  {
    G.bare_ground_albedo=0.25; //( default)
    warn="The required global parameter BARE_GROUND_ALBEDO was autogenerated with value "+to_string(G.bare_ground_albedo);
    if (chatty){WriteAdvisory(warn,false);}
  }
  autocalc=SetCalculableValue(G.snowfall_albthresh,Gtmp.snowfall_albthresh,Gtemplate.snowfall_albthresh);
  if (autocalc)
  {
    G.snowfall_albthresh=10; //[mm/d] (CRHM default)
    warn="The required global parameter SNOWFALL_ALBTHRESH was autogenerated with value "+to_string(G.snowfall_albthresh);
    if (chatty){WriteAdvisory(warn,false);}
  }
  autocalc=SetCalculableValue(G.UBC_exposure_fact,Gtmp.UBC_exposure_fact,Gtemplate.UBC_exposure_fact);
  if (autocalc)
  {
    G.UBC_exposure_fact=0.01;  //(UBC default)
    warn="The required global parameter UBC_EXPOSURE_FACT was autogenerated with value "+to_string(G.UBC_exposure_fact);
    if (chatty){WriteAdvisory(warn,false);}
  }
  autocalc=SetCalculableValue(G.UBC_cloud_penet,Gtmp.UBC_cloud_penet,Gtemplate.UBC_cloud_penet);
  if (autocalc)
  {
    G.UBC_cloud_penet=0.25;  //(UBC default)
    warn="The required global parameter UBC_CLOUD_PENET was autogenerated with value "+to_string(G.UBC_cloud_penet);
    if (chatty){WriteAdvisory(warn,false);}
  }
  autocalc=SetCalculableValue(G.UBC_LW_forest_fact,Gtmp.UBC_LW_forest_fact,Gtemplate.UBC_LW_forest_fact);
  if (autocalc)
  {
    G.UBC_LW_forest_fact=0.76;  //(UBC default)
    warn="The required global parameter UBC_LW_FOREST_FACT was autogenerated with value "+to_string(G.UBC_LW_forest_fact);
    if (chatty){WriteAdvisory(warn,false);}
  }
  autocalc=SetCalculableValue(G.UBC_flash_ponding,Gtmp.UBC_flash_ponding,Gtemplate.UBC_flash_ponding);
  if (autocalc)
  {
    G.UBC_flash_ponding=36;  //[mm] (UBC default)
    warn="The required global parameter UBC_FLASH_PONDING was autogenerated with value "+to_string(G.UBC_flash_ponding);
    if (chatty){WriteAdvisory(warn,false);}
  }

  autocalc=SetCalculableValue(G.max_reach_seglength,Gtmp.max_reach_seglength,Gtemplate.max_reach_seglength);
  if (autocalc)
  {
    G.max_reach_seglength=DEFAULT_MAX_REACHLENGTH;  //[km] - defaults to one segment per reach
    warn="The required global parameter MAX_REACH_SEGLENGTH was autogenerated with value "+to_string(G.max_reach_seglength);
    if (chatty){WriteAdvisory(warn,false);}
  }

  autocalc=SetCalculableValue(G.reservoir_relax,Gtmp.reservoir_relax,Gtemplate.reservoir_relax);
  if(autocalc)
  {
    G.reservoir_relax=0.4;
  }

  autocalc=SetCalculableValue(G.assimilation_fact,Gtmp.assimilation_fact,Gtemplate.assimilation_fact);
  if(autocalc)
  {
    G.assimilation_fact=1.0;
  }

  autocalc=SetCalculableValue(G.assim_upstream_decay,Gtmp.assim_upstream_decay,Gtemplate.assim_upstream_decay);
  if(autocalc)
  {
    G.assim_upstream_decay=0.0; //[1/km]
  }
  autocalc=SetCalculableValue(G.assim_time_decay,Gtmp.assim_time_decay,Gtemplate.assim_time_decay);
  if(autocalc)
  {
    G.assim_time_decay=0.2; //[1/d]
  }
  autocalc=SetCalculableValue(G.reservoir_demand_mult,Gtmp.reservoir_demand_mult,Gtemplate.reservoir_demand_mult);
  if(autocalc)
  {
    G.reservoir_demand_mult=1.0; //default
  }
  autocalc=SetCalculableValue(G.HBVEC_lapse_rate,Gtmp.HBVEC_lapse_rate,Gtemplate.HBVEC_lapse_rate);
  if (autocalc)
  {
    G.HBVEC_lapse_rate=0.08; //Default [mm/d/km]
  }
  autocalc = SetCalculableValue(G.HBVEC_lapse_upper, Gtmp.HBVEC_lapse_upper, Gtemplate.HBVEC_lapse_upper);
  if (autocalc)
  {
    G.HBVEC_lapse_upper=0.0;//Default [mm/d/km]
  }
  autocalc=SetCalculableValue(G.HBVEC_lapse_elev,Gtmp.HBVEC_lapse_elev,Gtemplate.HBVEC_lapse_elev);
  if (autocalc)
  {
    G.HBVEC_lapse_elev=5000.0;//Default [m]
  }


  //Model-specific global parameters - cannot be autocomputed, must be specified by user
  //----------------------------------------------------------------------------
  SetSpecifiedValue(G.UBC_lapse_params.max_range_temp,Gtmp.UBC_lapse_params.max_range_temp,Gtemplate.UBC_lapse_params.max_range_temp,false,"UBC_MAX_RANGE_TEMP");
  SetSpecifiedValue(G.UBC_lapse_params.A0PPTP,Gtmp.UBC_lapse_params.A0PPTP,Gtemplate.UBC_lapse_params.A0PPTP,false,"UBC_A0PPTP");
  SetSpecifiedValue(G.UBC_lapse_params.A0STAB,Gtmp.UBC_lapse_params.A0STAB,Gtemplate.UBC_lapse_params.A0STAB,false,"UBC_A0STAB");
  SetSpecifiedValue(G.UBC_lapse_params.A0PELA,Gtmp.UBC_lapse_params.A0PELA,Gtemplate.UBC_lapse_params.A0PELA,false,"UBC_A0PELA");
  SetSpecifiedValue(G.UBC_lapse_params.A0TLXM,Gtmp.UBC_lapse_params.A0TLXM,Gtemplate.UBC_lapse_params.A0TLXM,false,"UBC_A0TLXM");
  SetSpecifiedValue(G.UBC_lapse_params.A0TLNH,Gtmp.UBC_lapse_params.A0TLNH,Gtemplate.UBC_lapse_params.A0TLNH,false,"UBC_A0TLNH");
  SetSpecifiedValue(G.UBC_lapse_params.A0TLNM,Gtmp.UBC_lapse_params.A0TLNM,Gtemplate.UBC_lapse_params.A0TLNM,false,"UBC_A0TLNM");
  SetSpecifiedValue(G.UBC_lapse_params.A0TLXH,Gtmp.UBC_lapse_params.A0TLXH,Gtemplate.UBC_lapse_params.A0TLXH,false,"UBC_A0TLXH");
  SetSpecifiedValue(G.UBC_lapse_params.E0LHI,Gtmp.UBC_lapse_params.E0LHI,Gtemplate.UBC_lapse_params.E0LHI,false,"UBC_E0LHI");
  SetSpecifiedValue(G.UBC_lapse_params.E0LLOW,Gtmp.UBC_lapse_params.E0LLOW,Gtemplate.UBC_lapse_params.E0LLOW,false,"UBC_E0LLOW");
  SetSpecifiedValue(G.UBC_lapse_params.E0LMID,Gtmp.UBC_lapse_params.E0LMID,Gtemplate.UBC_lapse_params.E0LMID,false,"UBC_E0LMID");
  SetSpecifiedValue(G.UBC_lapse_params.P0GRADL,Gtmp.UBC_lapse_params.P0GRADL,Gtemplate.UBC_lapse_params.P0GRADL,false,"UBC_P0GRADL");
  SetSpecifiedValue(G.UBC_lapse_params.P0GRADM,Gtmp.UBC_lapse_params.P0GRADM,Gtemplate.UBC_lapse_params.P0GRADM,false,"UBC_P0GRADM");
  SetSpecifiedValue(G.UBC_lapse_params.P0GRADU,Gtmp.UBC_lapse_params.P0GRADU,Gtemplate.UBC_lapse_params.P0GRADU,false,"UBC_P0GRADU");
  SetSpecifiedValue(G.UBC_lapse_params.P0TEDL,Gtmp.UBC_lapse_params.P0TEDL,Gtemplate.UBC_lapse_params.P0TEDL,false,"UBC_P0TEDL");
  SetSpecifiedValue(G.UBC_lapse_params.P0TEDU,Gtmp.UBC_lapse_params.P0TEDU,Gtemplate.UBC_lapse_params.P0TEDU,false,"UBC_P0TEDU");

  SetSpecifiedValue(G.UBC_snow_params.ALBASE,Gtmp.UBC_snow_params.ALBASE,Gtemplate.UBC_snow_params.ALBASE,false,"UBC_ALBASE");
  SetSpecifiedValue(G.UBC_snow_params.ALBREC,Gtmp.UBC_snow_params.ALBREC,Gtemplate.UBC_snow_params.ALBREC,false,"UBC_ALBREC");
  SetSpecifiedValue(G.UBC_snow_params.ALBSNW,Gtmp.UBC_snow_params.ALBSNW,Gtemplate.UBC_snow_params.ALBSNW,false,"UBC_ALBSNW");
  SetSpecifiedValue(G.UBC_snow_params.MAX_CUM_MELT,Gtmp.UBC_snow_params.MAX_CUM_MELT,Gtemplate.UBC_snow_params.MAX_CUM_MELT,false,"UBC_MAX_CUM_MELT");

  SetSpecifiedValue(G.UBC_GW_split,Gtmp.UBC_GW_split,Gtemplate.UBC_GW_split,false,"UBC_GW_SPLIT");
  SetSpecifiedValue(G.airsnow_coeff,Gtmp.airsnow_coeff,Gtemplate.airsnow_coeff,false,"AIRSNOW_COEFF");
  SetSpecifiedValue(G.avg_annual_snow,Gtmp.avg_annual_snow,Gtemplate.avg_annual_snow,false,"AVG_ANNUAL_SNOW");
  SetSpecifiedValue(G.avg_annual_runoff,Gtmp.avg_annual_runoff,Gtemplate.avg_annual_runoff,false,"AVG_ANNUAL_RUNOFF");
  SetSpecifiedValue(G.MOHYSE_PET_coeff,Gtmp.MOHYSE_PET_coeff,Gtemplate.MOHYSE_PET_coeff,false,"MOHYSE_PET_COEFF");
  SetSpecifiedValue(G.SWI_reduct_coeff,Gtmp.SWI_reduct_coeff,Gtemplate.SWI_reduct_coeff,false,"SWI_REDUCT_COEFF");
  SetSpecifiedValue(G.init_stream_temp,Gtmp.init_stream_temp,Gtemplate.init_stream_temp,false,"INIT_STREAM_TEMP");
  SetSpecifiedValue(G.windvel_icept,Gtmp.windvel_icept,Gtemplate.windvel_icept,false,"WINDVEL_ICEPT");
  SetSpecifiedValue(G.windvel_scale,Gtmp.windvel_scale,Gtemplate.windvel_scale,false,"WINDVEL_SCALE");

}

//////////////////////////////////////////////////////////////////
/// \brief Sets default global parameters properties
/// \param &T [out] Terrain properties class
/// \param is_template [in] True if the default value being set is for the template class
//
void CGlobalParams::InitializeGlobalParameters(global_struct &g, bool is_template)//static
{


  /*
  //FUTURE GENERALIZED IMPLEMENTATION
  for(int i=0; i<_nGPParams; i++){
    _aGPvalues[i]=DefaultParameterValue(is_template,true);
    if (_aGPNames[i]==MAX_REACH_SEGLENGTH){
      _aGPvalues[i]=DEFAULT_MAX_REACHLENGTH;
    }
  }
  */

  g.max_reach_seglength  =DEFAULT_MAX_REACHLENGTH;//defaults to one segment per reach
  g.reservoir_relax      =0.4;
  g.assim_upstream_decay =0.0;
  g.assim_time_decay     =0.2;
  g.assimilation_fact    =1.0;
  g.reservoir_demand_mult=1.0;
  g.reference_flow_mult  =10.0;

  //calculable /estimable parameters
  g.snow_SWI            =DefaultParameterValue(is_template,true);
  g.snow_SWI_min        =DefaultParameterValue(is_template,true);
  g.snow_SWI_max        =DefaultParameterValue(is_template,true);
  g.snow_temperature    =DefaultParameterValue(is_template,true);
  g.snow_roughness      =DefaultParameterValue(is_template,true);
  g.rainsnow_temp       =DefaultParameterValue(is_template,true);
  g.rainsnow_delta      =DefaultParameterValue(is_template,true);
  g.adiabatic_lapse     =DefaultParameterValue(is_template,true);
  g.wet_adiabatic_lapse =DefaultParameterValue(is_template,true);
  g.precip_lapse        =DefaultParameterValue(is_template,true);
  g.airsnow_coeff       =DefaultParameterValue(is_template,true);
  g.MOHYSE_PET_coeff    =DefaultParameterValue(is_template,true);
  g.max_snow_albedo     =DefaultParameterValue(is_template,true);
  g.min_snow_albedo     =DefaultParameterValue(is_template,true);
  g.alb_decay_cold      =DefaultParameterValue(is_template,true);
  g.alb_decay_melt      =DefaultParameterValue(is_template,true);
  g.bare_ground_albedo  =DefaultParameterValue(is_template,true);
  g.snowfall_albthresh  =DefaultParameterValue(is_template,true);

  g.HBVEC_lapse_rate    =DefaultParameterValue(is_template,true);
  g.HBVEC_lapse_upper   =DefaultParameterValue(is_template,true);
  g.HBVEC_lapse_elev    =DefaultParameterValue(is_template,true);

  g.max_SWE_surface         =DefaultParameterValue(is_template,true);
  g.TOC_multiplier          =DefaultParameterValue(is_template,true);
  g.TIME_TO_PEAK_multiplier =DefaultParameterValue(is_template,true);
  g.GAMMA_SHAPE_multiplier  =DefaultParameterValue(is_template,true);
  g.GAMMA_SCALE_multiplier  =DefaultParameterValue(is_template,true);

  //model-specific parameters
  g.avg_annual_snow     =DefaultParameterValue(is_template,false);
  g.avg_annual_runoff   =DefaultParameterValue(is_template,false);
  g.init_stream_temp    =DefaultParameterValue(is_template,false);
  g.windvel_icept       =DefaultParameterValue(is_template,false);
  g.windvel_scale       =DefaultParameterValue(is_template,false);

  for (int i=0;i<12;i++){
    g.UBC_s_corr[i]=DefaultParameterValue(is_template,true);
    g.UBC_n_corr[i]=DefaultParameterValue(is_template,true);
  }

  g.UBC_exposure_fact      =DefaultParameterValue(is_template,true);
  g.UBC_cloud_penet        =DefaultParameterValue(is_template,true);
  g.UBC_LW_forest_fact     =DefaultParameterValue(is_template,true);
  g.UBC_flash_ponding      =DefaultParameterValue(is_template,true);

  g.UBC_GW_split           =DefaultParameterValue(is_template,false);
  g.UBC_lapse_params.A0PELA=DefaultParameterValue(is_template,false);
  g.UBC_lapse_params.A0PPTP=DefaultParameterValue(is_template,false);
  g.UBC_lapse_params.A0STAB=DefaultParameterValue(is_template,false);
  g.UBC_lapse_params.A0TLXM =DefaultParameterValue(is_template,false);
  g.UBC_lapse_params.A0TLNH =DefaultParameterValue(is_template,false);
  g.UBC_lapse_params.A0TLNM =DefaultParameterValue(is_template,false);
  g.UBC_lapse_params.A0TLXH =DefaultParameterValue(is_template,false);
  g.UBC_lapse_params.E0LHI  =DefaultParameterValue(is_template,false);
  g.UBC_lapse_params.E0LLOW =DefaultParameterValue(is_template,false);
  g.UBC_lapse_params.E0LMID =DefaultParameterValue(is_template,false);
  g.UBC_lapse_params.P0GRADL=DefaultParameterValue(is_template,false);
  g.UBC_lapse_params.P0GRADM=DefaultParameterValue(is_template,false);
  g.UBC_lapse_params.P0GRADU=DefaultParameterValue(is_template,false);
  g.UBC_lapse_params.P0TEDL =DefaultParameterValue(is_template,false);
  g.UBC_lapse_params.P0TEDU =DefaultParameterValue(is_template,false);
  g.UBC_lapse_params.max_range_temp =DefaultParameterValue(is_template,false);

  g.UBC_snow_params.ALBASE  =DefaultParameterValue(is_template,false);
  g.UBC_snow_params.ALBREC  =DefaultParameterValue(is_template,false);
  g.UBC_snow_params.ALBSNW  =DefaultParameterValue(is_template,false);
  g.UBC_snow_params.MAX_CUM_MELT =DefaultParameterValue(is_template,false);

  g.SWI_reduct_coeff =DefaultParameterValue(is_template,false);

}
//////////////////////////////////////////////////////////////////
/// \brief Sets the value of the terrain property corresponding to param_name
/// \param &G [out] Global Property structure
/// \param param_name [in] Parameter identifier
/// \param value [in] Value of parameter to be set
/// \note not really needed, but can support future parse options
//
void  CGlobalParams::SetGlobalProperty (const string &param_name,
                                        const double &value)
{
  SetGlobalProperty(G,param_name,value);
}
//////////////////////////////////////////////////////////////////
/// \brief Sets the value of the terrain property corresponding to param_name
/// \param &G [out] Global Property structure
/// \param param_name [in] Parameter identifier
/// \param value [in] Value of parameter to be set
/// \note not really needed, but can support future parse options
//
void  CGlobalParams::SetGlobalProperty (global_struct &G,
                                        string       param_name,
                                        const double value)
{
  string name;
  name = StringToUppercase(param_name);

  /*for(int i=0; i<_nGPParams; i++){
    if(!name.compare(_aGPNames[i])){_aGPvalues[i]=value; }
  }*/

  if      (!name.compare("SNOW_SWI"            )){G.snow_SWI=value;}
  else if (!name.compare("SNOW_SWI_MIN"        )){G.snow_SWI_min=value; }
  else if (!name.compare("SNOW_SWI_MAX"        )){G.snow_SWI_max=value; }
  else if (!name.compare("SWI_REDUCT_COEFF"    )){G.SWI_reduct_coeff=value; }
  else if (!name.compare("SNOW_TEMPERATURE"    )){G.snow_temperature=value;}
  else if (!name.compare("SNOW_ROUGHNESS"      )){G.snow_roughness=value;}
  else if (!name.compare("RAINSNOW_TEMP"       )){G.rainsnow_temp=value;}
  else if (!name.compare("RAINSNOW_DELTA"      )){G.rainsnow_delta=value;}
  else if (!name.compare("ADIABATIC_LAPSE"     )){G.adiabatic_lapse=value;}
  else if (!name.compare("REFERENCE_FLOW_MULT" )){G.reference_flow_mult=value;}
  else if (!name.compare("WET_ADIABATIC_LAPSE" )){G.wet_adiabatic_lapse=value;}
  else if (!name.compare("PRECIP_LAPSE"        )){G.precip_lapse=value;}

  else if (!name.compare("TOC_MULTIPLIER"          )){G.TOC_multiplier=value;}
  else if (!name.compare("TIME_TO_PEAK_MULTIPLIER" )){G.TIME_TO_PEAK_multiplier=value;}
  else if (!name.compare("GAMMA_SHAPE_MULTIPLIER"  )){G.GAMMA_SHAPE_multiplier=value;}
  else if (!name.compare("GAMMA_SCALE_MULTIPLIER"  )){G.GAMMA_SCALE_multiplier=value;}

  //WARNING: this sets *all* 12 SW correction parameters to "value".
  else if (!name.compare("UBC_SW_N_CORR"       )){for (int i=0;i<12;i++){G.UBC_n_corr[i]=value;}}
  else if (!name.compare("UBC_SW_S_CORR"       )){for (int i=0;i<12;i++){G.UBC_s_corr[i]=value;}}

  else if (!name.compare("MAX_SNOW_ALBEDO"     )){G.max_snow_albedo=value;}
  else if (!name.compare("MIN_SNOW_ALBEDO"     )){G.min_snow_albedo=value;}
  else if (!name.compare("ALB_DECAY_COLD"      )){G.alb_decay_cold=value;}
  else if (!name.compare("ALB_DECAY_MELT"      )){G.alb_decay_melt=value;}
  else if (!name.compare("BARE_GROUND_ALBEDO"  )){G.bare_ground_albedo=value;}
  else if (!name.compare("SNOWFALL_ALBTHRESH"  )){G.snowfall_albthresh=value;}

  else if (!name.compare("UBC_ALBASE"          )){G.UBC_snow_params.ALBASE=value;}
  else if (!name.compare("UBC_ALBREC"          )){G.UBC_snow_params.ALBREC=value;}
  else if (!name.compare("UBC_ALBSNW"          )){G.UBC_snow_params.ALBSNW=value;}
  else if (!name.compare("UBC_MAX_CUM_MELT"    )){G.UBC_snow_params.MAX_CUM_MELT=value;}
  else if (!name.compare("UBC_GW_SPLIT"        )){G.UBC_GW_split=value;}
  else if (!name.compare("UBC_FLASH_PONDING"   )){G.UBC_flash_ponding=value;}
  else if (!name.compare("UBC_EXPOSURE_FACT"   )){G.UBC_exposure_fact=value;}
  else if (!name.compare("UBC_CLOUD_PENET"     )){G.UBC_cloud_penet=value;}
  else if (!name.compare("UBC_LW_FOREST_FACT"  )){G.UBC_LW_forest_fact=value;}

  else if (!name.compare("UBC_A0PELA"  )){G.UBC_lapse_params.A0PELA=value;}
  else if (!name.compare("UBC_A0PPTP"  )){G.UBC_lapse_params.A0PPTP=value;}
  else if (!name.compare("UBC_A0STAB"  )){G.UBC_lapse_params.A0STAB=value;}
  else if (!name.compare("UBC_A0TLXM"  )){G.UBC_lapse_params.A0TLXM=value;}
  else if (!name.compare("UBC_A0TLNH"  )){G.UBC_lapse_params.A0TLNH=value;}
  else if (!name.compare("UBC_A0TLNM"  )){G.UBC_lapse_params.A0TLNM=value;}
  else if (!name.compare("UBC_A0TLXH"  )){G.UBC_lapse_params.A0TLXH=value;}
  else if (!name.compare("UBC_E0LHI"   )){G.UBC_lapse_params.E0LHI=value;}
  else if (!name.compare("UBC_E0LLOW"  )){G.UBC_lapse_params.E0LLOW=value;}
  else if (!name.compare("UBC_E0LMID"  )){G.UBC_lapse_params.E0LMID=value;}
  else if (!name.compare("UBC_P0GRADL" )){G.UBC_lapse_params.P0GRADL=value;}
  else if (!name.compare("UBC_P0GRADM" )){G.UBC_lapse_params.P0GRADM=value;}
  else if (!name.compare("UBC_P0GRADU" )){G.UBC_lapse_params.P0GRADU=value;}
  else if (!name.compare("UBC_P0TEDL"  )){G.UBC_lapse_params.P0TEDL=value;}
  else if (!name.compare("UBC_P0TEDU"  )){G.UBC_lapse_params.P0TEDU=value;}
  else if (!name.compare("UBC_MAX_RANGE_TEMP"  )){G.UBC_lapse_params.max_range_temp=value;}

  else if (!name.compare("AIRSNOW_COEFF"       )){G.airsnow_coeff=value;}
  else if (!name.compare("AVG_ANNUAL_SNOW"     )){G.avg_annual_snow=value;}
  else if (!name.compare("AVG_ANNUAL_RUNOFF"   )){G.avg_annual_runoff=value;}
  else if (!name.compare("INIT_STREAM_TEMP"    )){G.init_stream_temp=value;}
  else if (!name.compare("MAX_SWE_SURFACE"     )){G.max_SWE_surface=value;}
  else if (!name.compare("MOHYSE_PET_COEFF"    )){G.MOHYSE_PET_coeff=value;}
  else if (!name.compare("MAX_REACH_SEGLENGTH" )){G.max_reach_seglength=value;}
  else if (!name.compare("RESERVOIR_RELAX"     )){G.reservoir_relax=value; }
  else if (!name.compare("ASSIMILATION_FACT"   )){G.assimilation_fact=value; }
  else if (!name.compare("ASSIM_UPSTREAM_DECAY")){G.assim_upstream_decay=value; }
  else if (!name.compare("ASSIM_TIME_DECAY"    )){G.assim_time_decay=value; }
  else if (!name.compare("RESERVOIR_DEMAND_MULT")){G.reservoir_demand_mult=value; }
  else if (!name.compare("WINDVEL_ICEPT"       )){G.windvel_icept=value; }
  else if (!name.compare("WINDVEL_SCALE"       )){G.windvel_scale=value; }
  else if (!name.compare("HBVEC_LAPSE_RATE"    )){G.HBVEC_lapse_rate=value; }
  else if (!name.compare("HBVEC_LAPSE_UPPER"   )){G.HBVEC_lapse_upper=value; }
  else if (!name.compare("HBVEC_LAPSE_ELEV"    )){G.HBVEC_lapse_elev=value; }
  else{
    WriteWarning("CGlobalParams::SetGlobalProperty: Unrecognized/invalid global parameter name ("+name+") in .rvp file",false);

  }
}

//////////////////////////////////////////////////////////////////
/// \brief gets global property corresponding to param_name
/// \param param_name [in] Parameter identifier
/// \returns value of parameter
//
double CGlobalParams::GetParameter(const string param_name)
{
  return CGlobalParams::GetGlobalProperty(G,param_name);
}

//////////////////////////////////////////////////////////////////
/// \brief gets global property corresponding to param_name
/// \param param_name [in] Parameter identifier
/// \returns value of parameter
//
double *CGlobalParams::GetAddress(const string name)
{
  if      (!name.compare("SNOW_SWI"            )){return &(G.snow_SWI);}
  else if (!name.compare("SNOW_TEMPERATURE"    )){return &(G.snow_temperature);}
  else if (!name.compare("RAINSNOW_TEMP"       )){return &(G.rainsnow_temp);}
  else if (!name.compare("RAINSNOW_DELTA"      )){return &(G.rainsnow_delta);}
  else if (!name.compare("ADIABATIC_LAPSE"     )){return &(G.adiabatic_lapse);}
  else if (!name.compare("REFERENCE_FLOW_MULT" )){return &(G.reference_flow_mult);}
  else if (!name.compare("WET_ADIABATIC_LAPSE" )){return &(G.wet_adiabatic_lapse);}
  else if (!name.compare("PRECIP_LAPSE"        )){return &(G.precip_lapse);}

  return NULL;
}

///////////////////////////////////////////////////////////////////////////
/// \brief Returns global parameter value corresponding to param_name from global variable structure provided
/// \param &G [in] global variable structure
/// \param param_name [in] Parameter name
/// \return parameter value corresponding to parameter name
//
double CGlobalParams::GetGlobalProperty(const global_struct &G, string  param_name, const bool strict)
{

  string name;
  name = StringToUppercase(param_name);

  /*for(int i=0; i<_nGPParams; i++){
    if(!name.compare(_GPNames[i])){ return _GPvals[i]; break;}
  }*/

  if      (!name.compare("RAINSNOW_TEMP"       )){return G.rainsnow_temp;}
  else if (!name.compare("SNOW_SWI"            )){return G.snow_SWI;}
  else if (!name.compare("SNOW_SWI_MIN"        )){return G.snow_SWI_min; }
  else if (!name.compare("SNOW_SWI_MAX"        )){return G.snow_SWI_max; }
  else if (!name.compare("SWI_REDUCT_COEFF"    )){return G.SWI_reduct_coeff; }
  else if (!name.compare("SNOW_TEMPERATURE"    )){return G.snow_temperature;}
  else if (!name.compare("SNOW_ROUGHNESS"      )){return G.snow_roughness;}
  else if (!name.compare("RAINSNOW_DELTA"      )){return G.rainsnow_delta;}
  else if (!name.compare("ADIABATIC_LAPSE"     )){return G.adiabatic_lapse;}
  else if (!name.compare("REFERENCE_FLOW_MULT" )){return G.reference_flow_mult;}
  else if (!name.compare("WET_ADIABATIC_LAPSE" )){return G.wet_adiabatic_lapse;}
  else if (!name.compare("PRECIP_LAPSE"        )){return G.precip_lapse;}

  else if (!name.compare("TOC_MULTIPLIER"          )){return G.TOC_multiplier;}
  else if (!name.compare("TIME_TO_PEAK_MULTIPLIER" )){return G.TIME_TO_PEAK_multiplier;}
  else if (!name.compare("GAMMA_SHAPE_MULTIPLIER"  )){return G.GAMMA_SHAPE_multiplier;}
  else if (!name.compare("GAMMA_SCALE_MULTIPLIER"  )){return G.GAMMA_SCALE_multiplier;}

  //WARNING: this get only returns the first (zero-index) element of the 12 SW correction parameters
  else if (!name.compare("UBC_SW_N_CORR"        )){return G.UBC_n_corr[0];}
  else if (!name.compare("UBC_SW_S_CORR"        )){return G.UBC_s_corr[0];}

  else if (!name.compare("MAX_SNOW_ALBEDO"     )){return G.max_snow_albedo;}
  else if (!name.compare("MIN_SNOW_ALBEDO"     )){return G.min_snow_albedo;}
  else if (!name.compare("ALB_DECAY_COLD"      )){return G.alb_decay_cold;}
  else if (!name.compare("ALB_DECAY_MELT"      )){return G.alb_decay_melt;}
  else if (!name.compare("BARE_GROUND_ALBEDO"  )){return G.bare_ground_albedo;}
  else if (!name.compare("SNOWFALL_ALBTHRESH"  )){return G.snowfall_albthresh;}
  else if (!name.compare("UBC_ALBASE"          )){return G.UBC_snow_params.ALBASE;}
  else if (!name.compare("UBC_ALBREC"          )){return G.UBC_snow_params.ALBREC;}
  else if (!name.compare("UBC_ALBSNW"          )){return G.UBC_snow_params.ALBSNW;}
  else if (!name.compare("UBC_MAX_CUM_MELT"    )){return G.UBC_snow_params.MAX_CUM_MELT;}
  else if (!name.compare("UBC_GW_SPLIT"        )){return G.UBC_GW_split;}
  else if (!name.compare("UBC_FLASH_PONDING"   )){return G.UBC_flash_ponding;}
  else if (!name.compare("UBC_EXPOSURE_FACT"   )){return G.UBC_exposure_fact;}
  else if (!name.compare("UBC_CLOUD_PENET"     )){return G.UBC_cloud_penet;}
  else if (!name.compare("UBC_LW_FOREST_FACT"  )){return G.UBC_LW_forest_fact;}

  else if (!name.compare("UBC_A0PELA"  )){return G.UBC_lapse_params.A0PELA;}
  else if (!name.compare("UBC_A0PPTP"  )){return G.UBC_lapse_params.A0PPTP;}
  else if (!name.compare("UBC_A0STAB"  )){return G.UBC_lapse_params.A0STAB;}
  else if (!name.compare("UBC_A0TLXM"  )){return G.UBC_lapse_params.A0TLXM;}
  else if (!name.compare("UBC_A0TLNH"  )){return G.UBC_lapse_params.A0TLNH;}
  else if (!name.compare("UBC_A0TLNM"  )){return G.UBC_lapse_params.A0TLNM;}
  else if (!name.compare("UBC_A0TLXH"  )){return G.UBC_lapse_params.A0TLXH;}
  else if (!name.compare("UBC_E0LHI"   )){return G.UBC_lapse_params.E0LHI;}
  else if (!name.compare("UBC_E0LLOW"  )){return G.UBC_lapse_params.E0LLOW;}
  else if (!name.compare("UBC_E0LMID"  )){return G.UBC_lapse_params.E0LMID;}
  else if (!name.compare("UBC_P0GRADL" )){return G.UBC_lapse_params.P0GRADL;}
  else if (!name.compare("UBC_P0GRADM" )){return G.UBC_lapse_params.P0GRADM;}
  else if (!name.compare("UBC_P0GRADU" )){return G.UBC_lapse_params.P0GRADU;}
  else if (!name.compare("UBC_P0TEDL"  )){return G.UBC_lapse_params.P0TEDL;}
  else if (!name.compare("UBC_P0TEDU"  )){return G.UBC_lapse_params.P0TEDU;}
  else if (!name.compare("UBC_MAX_RANGE_TEMP"  )){return G.UBC_lapse_params.max_range_temp;}

  else if (!name.compare("AIRSNOW_COEFF"        )){return G.airsnow_coeff;}
  else if (!name.compare("AVG_ANNUAL_SNOW"      )){return G.avg_annual_snow;}
  else if (!name.compare("AVG_ANNUAL_RUNOFF"    )){return G.avg_annual_runoff;}
  else if (!name.compare("INIT_STREAM_TEMP"     )){return G.init_stream_temp;}
  else if (!name.compare("MAX_SWE_SURFACE"      )){return G.max_SWE_surface;}
  else if (!name.compare("MOHYSE_PET_COEFF"     )){return G.MOHYSE_PET_coeff;}
  else if (!name.compare("MAX_REACH_SEGLENGTH"  )){return G.max_reach_seglength;}
  else if (!name.compare("RESERVOIR_RELAX"      )){return G.reservoir_relax; }
  else if (!name.compare("ASSIMILATION_FACT"    )){return G.assimilation_fact; }
  else if (!name.compare("ASSIM_UPSTREAM_DECAY" )){return G.assim_upstream_decay; }
  else if (!name.compare("ASSIM_TIME_DECAY"     )){return G.assim_time_decay; }
  else if (!name.compare("RESERVOIR_DEMAND_MULT")){return G.reservoir_demand_mult; }
  else if (!name.compare("WINDVEL_ICEPT"        )){return G.windvel_icept; }
  else if (!name.compare("WINDVEL_SCALE"        )){return G.windvel_scale; }
  else if (!name.compare("HBVEC_LAPSE_RATE"     )){return G.HBVEC_lapse_rate; }
  else if (!name.compare("HBVEC_LAPSE_UPPER"    )){return G.HBVEC_lapse_upper; }
  else if (!name.compare("HBVEC_LAPSE_ELEV"     )){return G.HBVEC_lapse_elev; }
  else{
    if (strict){
      string msg="CGlobalParams::GetParameter: Unrecognized/invalid global parameter name in .rvp file: "+name;
      ExitGracefully(msg.c_str(),BAD_DATA_WARN);
      return 0.0;
    }
    else {
      return INDEX_NOT_FOUND;
    }
  }

}
