/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Properties.h"
#include "SoilAndLandClasses.h"
/*****************************************************************
   Constructor / Destructor
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the LandUseClass constructor
/// \param name [in] String nickname for land use class
//
CLandUseClass::CLandUseClass(const string name)
{
  this->S.landuse_name=name;
  if (!DynArrayAppend((void**&)(pAllLUClasses),(void*)(this),NumLUClasses)){
    ExitGracefully("CLandUseClass::Constructor: creating NULL land use class",BAD_DATA);};
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CLandUseClass::~CLandUseClass()
{
  if (DESTRUCTOR_DEBUG){cout<<"  DELETING LAND USE CLASS "<<endl;}
}

///////////////////////////////////////////////////////////////////
/// \brief Returns reference to surface properties
/// \return Surface properties associated with land use class
//
const surface_struct *CLandUseClass::GetSurfaceStruct() const{return &S;}

///////////////////////////////////////////////////////////////////
/// \brief Return nick name identifier of land use class
/// \return nick name identifier of land use class
//
string                CLandUseClass::GetLanduseName  () const{return S.landuse_name;}

/*****************************************************************
   Static Initialization, Accessors, Destructors
*****************************************************************/
CLandUseClass **CLandUseClass::pAllLUClasses=NULL;
int             CLandUseClass::NumLUClasses=0;

//////////////////////////////////////////////////////////////////
/// \brief Return number of land use classes
/// \return Number of land use classes
//
int CLandUseClass::GetNumClasses(){
  return NumLUClasses;
}

//////////////////////////////////////////////////////////////////
/// \brief Summarize LU class information to screen
//
void CLandUseClass::SummarizeToScreen()
{
  cout<<"==================="<<endl;
  cout<<"Land Use Class Summary:"<<NumLUClasses<<" LU/LT classes in database"<<endl;
  for (int c=0; c<NumLUClasses;c++){
    cout<<"-LULT. class \""<<pAllLUClasses[c]->GetLanduseName()<<"\" "<<endl;
    cout<<"    impermeable: "<<pAllLUClasses[c]->GetSurfaceStruct()->impermeable_frac*100<<" %"<<endl;
    cout<<"       forested: "<<pAllLUClasses[c]->GetSurfaceStruct()->forest_coverage*100<<" %"<<endl;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Destroy all LU classes
//
void CLandUseClass::DestroyAllLUClasses()
{
  if (DESTRUCTOR_DEBUG){cout <<"DESTROYING ALL LULT CLASSES"<<endl;}
  for (int c=0; c<NumLUClasses;c++){
    delete pAllLUClasses[c];
  }
  delete [] pAllLUClasses;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns the LU class corresponding to passed string
/// \details Converts string (e.g., "AGRICULTURAL" in HRU file) to LU class
///  can accept either lultclass index or lultclass tag
///  if string is invalid, returns NULL
/// \param s [in] LU class identifier (tag or index)
/// \return Pointer to LU class corresponding to identifier string s
//
CLandUseClass *CLandUseClass::StringToLUClass(const string s)
{
  string sup=StringToUppercase(s);
  for (int c=0;c<NumLUClasses;c++)
  {
    if (!sup.compare(StringToUppercase(pAllLUClasses[c]->GetLanduseName()))){return pAllLUClasses[c];}
    else if (s_to_i(s.c_str())==(c+1))                                      {return pAllLUClasses[c];}
  }
  return NULL;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns the land use  class corresponding to the passed index
///  if index is invalid, returns NULL
/// \param c [in] Soil class index
/// \return Reference to land use class corresponding to index c
//
const CLandUseClass *CLandUseClass::GetLUClass(int c)
{
  if ((c<0) || (c>=NumLUClasses)){return NULL;}
  return pAllLUClasses[c];
}
//////////////////////////////////////////////////////////////////
/// \brief Automatically calculates surface propeties
/// \details  Sets surface properties based upon simple lu/lt parameters
///  Input [Stmp & Rtmp] has been read from .rvp file - if parameter ==
///  AUTO_COMPLETE, then empirical relationships are used to estimate
///  parameters
///
/// \param &Stmp [in] Input LU parameters (read from .rvp file)
/// \param &Sdefault [in] Default LU parameters
//
void CLandUseClass::AutoCalculateLandUseProps(const surface_struct &Stmp,
                                              const surface_struct &Sdefault)
//const surface_struct &needed_params
{
  bool autocalc;
  string warn;
  bool chatty=true;

  //these parameters are required
  S.landuse_name    =Stmp.landuse_name;
  S.impermeable_frac=Stmp.impermeable_frac;
  ExitGracefullyIf(S.impermeable_frac<0.0 || S.impermeable_frac>1.0,"Invalid parameter value for IMPERMEABLE_FRAC: must be between 0 and 1",BAD_DATA_WARN);

  //Forest coverage
  autocalc=SetCalculableValue(S.forest_coverage,Stmp.forest_coverage,Sdefault.forest_coverage);
  if (autocalc)
  {
    //RELATIONSHIP REQUIRED
    S.forest_coverage =0.0;
    //no warning -default is no forest cover
  }
  ExitGracefullyIf(S.forest_coverage>1.0,"Invalid parameter value for FOREST_COVERAGE: must be between 0 and 1",BAD_DATA_WARN);

  autocalc=SetCalculableValue(S.forest_sparseness ,Stmp.forest_sparseness,Sdefault.forest_sparseness);
  if (autocalc)
  {
    S.forest_sparseness =0.0;  //Default - not sparse
    //no warning -default is no sparseness
  }
  //Standard surface properties
  //----------------------------------------------------------------------------
  //Roughness
  autocalc=SetCalculableValue(S.roughness,Stmp.roughness,Sdefault.roughness);
  if (autocalc)
  {
    //RELATIONSHIP REQUIRED
    S.roughness =0.0;
    //no warning -default is no roughness
  }

  autocalc=SetCalculableValue(S.max_sat_area_frac,Stmp.max_sat_area_frac,Sdefault.max_sat_area_frac);
  if (autocalc)
  {
    S.max_sat_area_frac =1.0;
    //no warning -default is no max
  }

  //Snow properties
  autocalc=SetCalculableValue(S.melt_factor,Stmp.melt_factor,Sdefault.melt_factor);
  if (autocalc)
  {
    //RELATIONSHIP REQUIRED?
    S.melt_factor =5.04;//[mm/K/d]/// \ref GAWSER ??
    warn="The required parameter MELT_FACTOR for land use class "+S.landuse_name+" was autogenerated with value "+to_string(S.melt_factor);
    if (chatty){WriteAdvisory(warn,false);}
  }
  autocalc=SetCalculableValue(S.DD_melt_temp,Stmp.DD_melt_temp,Sdefault.DD_melt_temp);
  if (autocalc)
  {
    S.DD_melt_temp =FREEZING_TEMP;//no warning -default is zero
  }
  autocalc=SetCalculableValue(S.min_melt_factor,Stmp.min_melt_factor,Sdefault.min_melt_factor);
  if (autocalc)
  {
    S.min_melt_factor =S.melt_factor;//[mm/K/d]
    //no warning
  }
  autocalc=SetCalculableValue(S.refreeze_factor,Stmp.refreeze_factor,Sdefault.refreeze_factor);
  if (autocalc)
  {
    S.refreeze_factor =5.04;//[mm/K/d]// \ref GAWSER
    warn="The required parameter REFREEZE_FACTOR for land use class "+S.landuse_name+" was autogenerated with value "+to_string(S.refreeze_factor);
    if (chatty){WriteAdvisory(warn,false);}
  }
  autocalc=SetCalculableValue(S.DD_refreeze_temp,Stmp.DD_refreeze_temp,Sdefault.DD_refreeze_temp);
  if(autocalc)
  {
    S.DD_refreeze_temp =0.0;//[C]
    warn="The required parameter DD_REFREEZE_TEMP for land use class "+S.landuse_name+" was autogenerated with value "+to_string(S.DD_refreeze_temp);
    if(chatty) { WriteAdvisory(warn,false); }
  }
  autocalc=SetCalculableValue(S.refreeze_exp,Stmp.refreeze_exp,Sdefault.refreeze_exp);
  if(autocalc)
  {
    S.refreeze_exp =1.0;//[C]
    warn="The required parameter REFREEZE_EXP for land use class "+S.landuse_name+" was autogenerated with value "+to_string(S.refreeze_exp);
    if(chatty) { WriteAdvisory(warn,false); }
  }
  autocalc=SetCalculableValue(S.HBV_melt_for_corr,Stmp.HBV_melt_for_corr,Sdefault.HBV_melt_for_corr);
  if (autocalc)
  {
    S.HBV_melt_for_corr =1.0;  //Default - no forest correction
    //no warning
  }
  autocalc=SetCalculableValue(S.HBV_melt_asp_corr,Stmp.HBV_melt_asp_corr,Sdefault.HBV_melt_asp_corr);
  if (autocalc)
  {
    S.HBV_melt_asp_corr =0.0;  //Default - no aspect correction
    //no warning
  }
  autocalc=SetCalculableValue(S.HBV_melt_glacier_corr,Stmp.HBV_melt_glacier_corr,Sdefault.HBV_melt_glacier_corr);
  if (autocalc)
  {
    S.HBV_melt_glacier_corr =1.0;  //Default - no glacier correction
    //no warning
  }
  autocalc=SetCalculableValue(S.ow_PET_corr,Stmp.ow_PET_corr,Sdefault.ow_PET_corr);
  if (autocalc)
  {
    S.ow_PET_corr =1.0;  //Default - no correction
    //no warning
  }
  autocalc=SetCalculableValue(S.lake_PET_corr,Stmp.lake_PET_corr,Sdefault.lake_PET_corr);
  if (autocalc)
  {
    S.lake_PET_corr =1.0;  //Default - no correction
    //no warning
  }
  autocalc=SetCalculableValue(S.forest_PET_corr,Stmp.forest_PET_corr,Sdefault.forest_PET_corr);
  if (autocalc)
  {
    S.forest_PET_corr =1.0;  //Default - no correction
    //no warning
  }
  autocalc=SetCalculableValue(S.priestleytaylor_coeff,Stmp.priestleytaylor_coeff,Sdefault.priestleytaylor_coeff);
  if(autocalc)
  {
    S.priestleytaylor_coeff =1.28;  //Default - no correction
  }
  autocalc=SetCalculableValue(S.SCS_Ia_fraction,Stmp.SCS_Ia_fraction,Sdefault.SCS_Ia_fraction);
  if (autocalc)
  {
    S.SCS_Ia_fraction=0.2;  //Default SCS Ia=0.2*S
    warn="The required parameter SCS_IA_FRACTION for land use class "+S.landuse_name+" was autogenerated with value "+to_string(S.SCS_Ia_fraction);
    if (chatty){WriteAdvisory(warn,false);}
  }

  autocalc=SetCalculableValue(S.snow_patch_limit,Stmp.snow_patch_limit,Sdefault.snow_patch_limit);
  if (autocalc)
  {
    S.snow_patch_limit=0.0;  //Default snow patch limit = 0.0 mm
    //no warning - default is no patching
  }

  autocalc = SetCalculableValue(S.conv_melt_mult, Stmp.conv_melt_mult, Sdefault.conv_melt_mult);
  if(autocalc)
  {
	  S.conv_melt_mult = 0.113; //UBCWM default
  }

  autocalc = SetCalculableValue(S.cond_melt_mult, Stmp.cond_melt_mult, Sdefault.cond_melt_mult);
  if (autocalc)
  {
	  S.cond_melt_mult = 0.44; //UBCWM default
  }

  autocalc = SetCalculableValue(S.rain_melt_mult, Stmp.rain_melt_mult, Sdefault.rain_melt_mult);
  if (autocalc)
  {
	  S.rain_melt_mult = 1.0; //UBCWM default
  }

  autocalc = SetCalculableValue(S.UBC_icept_factor, Stmp.UBC_icept_factor, Sdefault.UBC_icept_factor);
  if (autocalc){
    S.UBC_icept_factor = 0.0; //Temporary for UBCWM until translator repaired - should not have default
    warn="The required parameter UBC_ICEPT_FACTOR for land use class "+S.landuse_name+" was autogenerated with value "+to_string(S.UBC_icept_factor);
    if(chatty) { WriteAdvisory(warn,false); }
  }

  autocalc = SetCalculableValue(S.sky_view_factor,Stmp.sky_view_factor,Sdefault.sky_view_factor);
  if(autocalc) {
    S.sky_view_factor = 0.8; //80% of hemisphere  is visible from location --> 1.0 for really flat terrain
    warn="The required parameter SKY_VIEW_FACTOR for land use class "+S.landuse_name+" was autogenerated with value "+to_string(S.sky_view_factor);
    if(chatty) { WriteAdvisory(warn,false); }
  }

  autocalc = SetCalculableValue(S.geothermal_grad,Stmp.geothermal_grad,Sdefault.geothermal_grad);
  if(autocalc) {
    S.geothermal_grad = 0.0; //default is that geothermal gradient is ignored
  }
  autocalc = SetCalculableValue(S.min_wind_speed,Stmp.min_wind_speed,Sdefault.min_wind_speed);
  if(autocalc) {
    S.min_wind_speed = 1.0; 
    warn="The required parameter MIN_WIND_SPEED for land use class "+S.landuse_name+" was autogenerated with value 1.0 m/s";
    if(chatty) { WriteAdvisory(warn,false); }
  }
  autocalc = SetCalculableValue(S.max_wind_speed,Stmp.max_wind_speed,Sdefault.max_wind_speed);
  if(autocalc) {
    S.max_wind_speed = 8.0; //m/s
    warn="The required parameter MAX_WIND_SPEED for land use class "+S.landuse_name+" was autogenerated with value 8.0 m/s";
    if(chatty) { WriteAdvisory(warn,false); }
  }
  autocalc = SetCalculableValue(S.stream_fraction,Stmp.stream_fraction,Sdefault.stream_fraction);
  if(autocalc) {
    S.stream_fraction = 0.0;
    warn="The required parameter STREAM_FRACTION for land use class "+S.landuse_name+" was autogenerated with value 0.0";
    if(chatty) { WriteAdvisory(warn,false); }
  }
  /*for(i=0;i<N_LU_PARAMETERS;i++) {
  if (!S.params[i].iscomputable){
    SetSpecifiedValue(S.params[i].value,Stmp.params[i].value,Sdefault.params[i].value,params[i].name);//(needed_params.partition_coeff>0.0)
  }
  }*/ // \todo[funct] - PARAMETEROVERHAUL (replaces below)


  //Model-specific LULT properties - cannot be autocomputed, must be specified by user
  //----------------------------------------------------------------------------
  bool needed=false;//TMP DEBUG (to be removed- this handled elsewhere in code)
  
  SetSpecifiedValue(S.partition_coeff,Stmp.partition_coeff,Sdefault.partition_coeff,needed,"PARTITION_COEFF");//(needed_params.partition_coeff>0.0)
  SetSpecifiedValue(S.SCS_CN,Stmp.SCS_CN,Sdefault.SCS_CN,needed,"SCS_CN");
  SetSpecifiedValue(S.dep_max,Stmp.dep_max,Sdefault.dep_max,needed,"DEP_MAX");
  SetSpecifiedValue(S.dep_max_flow,Stmp.dep_max_flow,Sdefault.dep_max_flow,needed,"DEP_MAX_FLOW");
  SetSpecifiedValue(S.dep_n,Stmp.dep_n,Sdefault.dep_n,needed,"DEP_N");
  SetSpecifiedValue(S.dep_threshhold,Stmp.dep_threshhold,Sdefault.dep_threshhold,needed,"DEP_THRESHHOLD");
  SetSpecifiedValue(S.dep_k,Stmp.dep_k,Sdefault.dep_k,needed,"DEP_K");
  SetSpecifiedValue(S.dep_seep_k,Stmp.dep_seep_k,Sdefault.dep_seep_k,needed,"DEP_SEEP_K");
  SetSpecifiedValue(S.dep_crestratio,Stmp.dep_crestratio,Sdefault.dep_crestratio,needed,"DEP_CRESTRATIO");
  SetSpecifiedValue(S.PDMROF_b,Stmp.PDMROF_b,Sdefault.PDMROF_b,needed,"PDMROF_B");
  SetSpecifiedValue(S.PDM_b,Stmp.PDM_b,Sdefault.PDM_b,needed,"PDM_B");
  SetSpecifiedValue(S.max_dep_area_frac,Stmp.max_dep_area_frac,Sdefault.max_dep_area_frac,needed,"MAX_DEP_AREA_FRAC");
  SetSpecifiedValue(S.ponded_exp,Stmp.ponded_exp,Sdefault.ponded_exp,needed,"PONDED_EXP");
  SetSpecifiedValue(S.uwfs_b,Stmp.uwfs_b,Sdefault.uwfs_b,needed,"UWFS_B");
  SetSpecifiedValue(S.uwfs_betamin,Stmp.uwfs_betamin,Sdefault.uwfs_betamin,needed,"UWFS_BETAMIN");
  SetSpecifiedValue(S.bf_loss_fraction,Stmp.bf_loss_fraction,Sdefault.bf_loss_fraction,needed,"BF_LOSS_FRACTION");
  SetSpecifiedValue(S.AWBM_areafrac1,Stmp.AWBM_areafrac1,Sdefault.AWBM_areafrac1,needed,"AWBM_AREAFRAC1");
  SetSpecifiedValue(S.AWBM_areafrac2,Stmp.AWBM_areafrac2,Sdefault.AWBM_areafrac2,needed,"AWBM_AREAFRAC2");
  SetSpecifiedValue(S.AWBM_bflow_index,Stmp.AWBM_bflow_index,Sdefault.AWBM_bflow_index,needed,"AWBM_BFLOW_INDEX");
  SetSpecifiedValue(S.lake_rel_coeff,Stmp.lake_rel_coeff,Sdefault.lake_rel_coeff,needed,"LAKE_REL_COEFF");
  SetSpecifiedValue(S.abst_percent,Stmp.abst_percent,Sdefault.abst_percent,needed,"ABST_PERCENT");
  SetSpecifiedValue(S.HBV_glacier_Kmin,Stmp.HBV_glacier_Kmin,Sdefault.HBV_glacier_Kmin,needed,"HBV_GLACIER_KMIN");
  SetSpecifiedValue(S.glac_storage_coeff,Stmp.glac_storage_coeff,Sdefault.glac_storage_coeff,needed,"GLAC_STORAGE_COEFF");
  SetSpecifiedValue(S.HBV_glacier_Ag,Stmp.HBV_glacier_Ag,Sdefault.HBV_glacier_Ag,needed,"HBV_GLACIER_AG");
  SetSpecifiedValue(S.CC_decay_coeff,Stmp.CC_decay_coeff,Sdefault.CC_decay_coeff,needed,"CC_DECAY_COEFF");
  SetSpecifiedValue(S.GR4J_x4,Stmp.GR4J_x4,Sdefault.GR4J_x4,needed,"GR4J_X4");
  SetSpecifiedValue(S.wind_exposure, Stmp.wind_exposure, Sdefault.wind_exposure, needed,"WIND_EXPOSURE");
  SetSpecifiedValue(S.fetch,Stmp.fetch,Sdefault.fetch,needed,"FETCH");
  SetSpecifiedValue(S.AET_coeff,Stmp.AET_coeff,Sdefault.AET_coeff,needed,"AET_COEFF");
  SetSpecifiedValue(S.max_melt_factor,Stmp.max_melt_factor,Sdefault.max_melt_factor,needed,"MAX_MELT_FACTOR");
  SetSpecifiedValue(S.DD_aggradation,Stmp.DD_aggradation,Sdefault.DD_aggradation,needed,"DD_AGGRADATION");
  SetSpecifiedValue(S.gamma_scale,Stmp.gamma_scale,Sdefault.gamma_scale,needed,"GAMMA_SCALE");
  SetSpecifiedValue(S.gamma_scale2,Stmp.gamma_scale2,Sdefault.gamma_scale2,needed,"GAMMA_SCALE2");
  SetSpecifiedValue(S.gamma_shape,Stmp.gamma_shape,Sdefault.gamma_shape,needed,"GAMMA_SHAPE");
  SetSpecifiedValue(S.gamma_shape2,Stmp.gamma_shape2,Sdefault.gamma_shape2,needed,"GAMMA_SHAPE2");
  SetSpecifiedValue(S.HMETS_runoff_coeff,Stmp.HMETS_runoff_coeff,Sdefault.HMETS_runoff_coeff,needed,"HMETS_RUNOFF_COEFF");
  SetSpecifiedValue(S.bsnow_distrib,Stmp.bsnow_distrib,Sdefault.bsnow_distrib,needed,"BSNOW_DISTRIB");
  SetSpecifiedValue(S.convection_coeff,Stmp.convection_coeff,Sdefault.convection_coeff,needed,"CONVECTION_COEFF");
  SetSpecifiedValue(S.pet_lin_coeff,Stmp.pet_lin_coeff,Sdefault.pet_lin_coeff,needed,"PET_LIN_COEFF");
}

//////////////////////////////////////////////////////////////////
/// \brief Sets default Surface properties
/// \details Initializes all surface properties to DEFAULT_VALUE
///  if is_template==true, initializes instead to NOT_SPECIFIED or AUTO_CALCULATE
/// \param &S [out] Surface properties class
/// \param is_template [in] True if the default value being set is for the template class
//
void CLandUseClass::InitializeSurfaceProperties(string name, surface_struct &S, bool is_template)
{
  
  S.landuse_name     =name;

  /*for(i=0;i<N_LU_PARAMETERS;i++) {
    S.params[i].value=DefaultParameterValue(is_template,S.params[i].iscomputable);
  }*/ // \todo[funct] - PARAMETEROVERHAUL (replaces below)

  //required parameters
  S.impermeable_frac =0.0;

  //Autocalculable parameters
  S.forest_coverage  =DefaultParameterValue(is_template,true);
  S.forest_sparseness=DefaultParameterValue(is_template,true);
  S.roughness        =DefaultParameterValue(is_template,true);
  S.melt_factor      =DefaultParameterValue(is_template,true);
  S.min_melt_factor  =DefaultParameterValue(is_template,true);
  S.DD_melt_temp     =DefaultParameterValue(is_template,true);
  S.refreeze_factor  =DefaultParameterValue(is_template,true);
  S.DD_refreeze_temp =DefaultParameterValue(is_template,true);
  S.refreeze_exp     =DefaultParameterValue(is_template,true);
  S.HBV_melt_asp_corr=DefaultParameterValue(is_template,true);
  S.HBV_melt_for_corr=DefaultParameterValue(is_template,true);
  S.HBV_melt_glacier_corr=DefaultParameterValue(is_template,true);//1.64
  S.max_sat_area_frac=DefaultParameterValue(is_template,true);//0.250;    //default [-]
  S.ow_PET_corr      =DefaultParameterValue(is_template,true);//1.0;      //[-]
  S.lake_PET_corr    =DefaultParameterValue(is_template,true);//1.0;      //[-]
  S.forest_PET_corr  =DefaultParameterValue(is_template,true);//1.0;      //[-]
  S.priestleytaylor_coeff  =DefaultParameterValue(is_template,true);//1.28;      //[-]
  S.SCS_Ia_fraction  =DefaultParameterValue(is_template,true);//0.2
  S.snow_patch_limit = DefaultParameterValue(is_template,true);//0.0
  S.conv_melt_mult = DefaultParameterValue(is_template, true);
  S.cond_melt_mult = DefaultParameterValue(is_template, true);
  S.rain_melt_mult = DefaultParameterValue(is_template, true);
  S.sky_view_factor   =DefaultParameterValue(is_template,true);
  S.geothermal_grad   =DefaultParameterValue(is_template,true);//~0.03
  S.min_wind_speed   =DefaultParameterValue(is_template,true);//~1 m/s
  S.max_wind_speed   =DefaultParameterValue(is_template,true);//~8 m/s
  S.stream_fraction  =DefaultParameterValue(is_template,true);//0

  //User-specified parameters
  S.partition_coeff   =DefaultParameterValue(is_template,false);//0.4;//needs reasonable defaults
  S.SCS_CN            =DefaultParameterValue(is_template,false);//50
  S.dep_max           =DefaultParameterValue(is_template,false);//6.29;     //[mm]
  S.dep_max_flow      =DefaultParameterValue(is_template,false);            //[mm/d]
  S.dep_n             =DefaultParameterValue(is_template,false);//1.0;      //[-]
  S.dep_crestratio    =DefaultParameterValue(is_template,false);//1.5;      //[mm]
  S.PDMROF_b          =DefaultParameterValue(is_template,false);//4;        //[-]
  S.PDM_b             =DefaultParameterValue(is_template,false);//4;        //[-]
  S.max_dep_area_frac =DefaultParameterValue(is_template,false);//0;        //[0..1]
  S.ponded_exp        =DefaultParameterValue(is_template,false);//2         //[-]
  S.uwfs_b            =DefaultParameterValue(is_template,false);//~1-10     //[-]
  S.uwfs_betamin      =DefaultParameterValue(is_template,false);//1+        //[-]
  S.bf_loss_fraction  =DefaultParameterValue(is_template,false);//0.0       //[0..1]
  S.AWBM_areafrac1    =DefaultParameterValue(is_template,false);//0.134     //[0..1]
  S.AWBM_areafrac2    =DefaultParameterValue(is_template,false);//0.433     //[0..1]
  S.AWBM_bflow_index  =DefaultParameterValue(is_template,false);//0.3       //[0..1]

  S.dep_threshhold    =DefaultParameterValue(is_template,false);//1.5;      //[mm]
  S.lake_rel_coeff    =DefaultParameterValue(is_template,false);//0.3;      //[1/d]
  S.dep_k             =DefaultParameterValue(is_template,false);//0.1;      //[1/d]
  S.dep_seep_k        =DefaultParameterValue(is_template,false);//0.1;      //[1/d]
  S.abst_percent      =DefaultParameterValue(is_template,false);//0.1;
  S.HBV_glacier_Kmin  =DefaultParameterValue(is_template,false);//0.05
  S.glac_storage_coeff=DefaultParameterValue(is_template,false);//0.10
  S.HBV_glacier_Ag    =DefaultParameterValue(is_template,false);//0.05 /mm
  S.CC_decay_coeff    =DefaultParameterValue(is_template,false);//0.0? [1/d]
  S.GR4J_x4           =DefaultParameterValue(is_template,false);//1.7 [d]
  S.UBC_icept_factor  =DefaultParameterValue(is_template,false);//0 [0..1]
  S.wind_exposure     =DefaultParameterValue(is_template,false);//
  S.fetch             =DefaultParameterValue(is_template,false);//300+ [m]
  S.AET_coeff         =DefaultParameterValue(is_template,false);//0.5 [0..1]
  S.max_melt_factor   =DefaultParameterValue(is_template,false);//~6
  S.DD_aggradation    =DefaultParameterValue(is_template,false);//~0.01
  S.gamma_scale       =DefaultParameterValue(is_template,false);//
  S.gamma_shape       =DefaultParameterValue(is_template,false);//
  S.gamma_scale2      =DefaultParameterValue(is_template,false);//
  S.gamma_shape2      =DefaultParameterValue(is_template,false);//
  S.HMETS_runoff_coeff=DefaultParameterValue(is_template,false);//0.4
  S.bsnow_distrib     =DefaultParameterValue(is_template,false);//0.4
  S.convection_coeff  =DefaultParameterValue(is_template,false);//~2
  S.pet_lin_coeff     =DefaultParameterValue(is_template,false);//~0.1-0.3
}
//////////////////////////////////////////////////////////////////
/// \brief Sets the value of the surface property corresponding to param_name
/// \param &S [out] Surface properties class
/// \param param_name [in] Parameter identifier
/// \param value [in] Value of parameter to be set
//
void  CLandUseClass::SetSurfaceProperty(const string &param_name,
                                        const double &value)
{
  SetSurfaceProperty(S,param_name,value);
}
//////////////////////////////////////////////////////////////////
/// \brief Sets the value of the surface property corresponding to param_name
/// \param &S [out] Surface properties class
/// \param param_name [in] Parameter identifier
/// \param value [in] Value of parameter to be set
//
void  CLandUseClass::SetSurfaceProperty(surface_struct &S,
                                        const string    param_name,
                                        const double value)
{
  string name;
  name = StringToUppercase(param_name);

  /*for(i=0;i<N_LU_PARAMETERS;i++) {
     if (!name.compare(S.params[i].name)){S.params[i].value=value;}
  }*/ // \todo[funct] - PARAMETEROVERHAUL (replaces below)

  if      (!name.compare("IMPERMEABLE_FRAC"       )){S.impermeable_frac=value;}
  else if (!name.compare("FOREST_COVERAGE"        )){S.forest_coverage=value;}
  else if (!name.compare("ROUGHNESS"              )){S.roughness=value;}

  else if (!name.compare("FOREST_SPARSENESS"      )){S.forest_sparseness=value;}
  else if (!name.compare("MELT_FACTOR"            )){S.melt_factor=value;}
  else if (!name.compare("MIN_MELT_FACTOR"        )){S.min_melt_factor=value;}
  else if (!name.compare("MAX_MELT_FACTOR"        )){S.max_melt_factor=value;}
  else if (!name.compare("DD_AGGRADATION"         )){S.DD_aggradation=value;}
  else if (!name.compare("DD_MELT_TEMP"           )){S.DD_melt_temp=value;}
  else if (!name.compare("REFREEZE_FACTOR"        )){S.refreeze_factor=value;}
  else if (!name.compare("DD_REFREEZE_TEMP"       )){S.DD_refreeze_temp=value; }
  else if (!name.compare("REFREEZE_EXP"           )){S.refreeze_exp=value; }
  else if (!name.compare("HBV_MELT_ASP_CORR"      )){S.HBV_melt_asp_corr=value;}
  else if (!name.compare("HBV_MELT_FOR_CORR"      )){S.HBV_melt_for_corr=value;}
  else if (!name.compare("MAX_SAT_AREA_FRAC"      )){S.max_sat_area_frac=value;}
  else if (!name.compare("HBV_MELT_GLACIER_CORR"  )){S.HBV_melt_glacier_corr=value;}
  else if (!name.compare("HBV_GLACIER_KMIN"       )){S.HBV_glacier_Kmin=value;}
  else if (!name.compare("GLAC_STORAGE_COEFF"     )){S.glac_storage_coeff=value;}
  else if (!name.compare("HBV_GLACIER_AG"         )){S.HBV_glacier_Ag=value;}
  else if (!name.compare("SNOW_PATCH_LIMIT"		  )){S.snow_patch_limit = value; }
  else if (!name.compare("CONV_MELT_MULT"		  )){S.conv_melt_mult = value; }
  else if (!name.compare("COND_MELT_MULT"		  )){S.cond_melt_mult = value; }
  else if (!name.compare("RAIN_MELT_MULT"		  )){S.rain_melt_mult = value; }
  else if (!name.compare("CC_DECAY_COEFF"         )){S.CC_decay_coeff=value;}
  else if (!name.compare("PARTITION_COEFF"        )){S.partition_coeff=value;}
  else if (!name.compare("SCS_CN"                 )){S.SCS_CN=value;}
  else if (!name.compare("SCS_IA_FRACTION"        )){S.SCS_Ia_fraction=value;}
  else if (!name.compare("DEP_MAX"                )){S.dep_max =value;}
  else if (!name.compare("DEP_MAX_FLOW"           )){S.dep_max_flow =value;}
  else if (!name.compare("DEP_N"                  )){S.dep_n =value;}
  else if (!name.compare("DEP_THRESHHOLD"         )){S.dep_threshhold =value;}
  else if (!name.compare("DEP_CRESTRATIO"         )){S.dep_crestratio =value;}
  else if (!name.compare("PDMROF_B"               )){S.PDMROF_b =value; }
  else if (!name.compare("PDM_B"                  )){S.PDM_b =value; }
  else if (!name.compare("MAX_DEP_AREA_FRAC"      )){S.max_dep_area_frac =value; }
  else if (!name.compare("PONDED_EXP"             )){S.ponded_exp =value; }
  else if (!name.compare("UWFS_B"                 )){S.uwfs_b =value; }
  else if (!name.compare("UWFS_BETAMIN"           )){S.uwfs_betamin =value; }
  else if (!name.compare("BF_LOSS_FRACTION"       )){S.bf_loss_fraction =value; }
  else if (!name.compare("AWBM_AREAFRAC1"         )){S.AWBM_areafrac1 =value; }
  else if (!name.compare("AWBM_AREAFRAC2"         )){S.AWBM_areafrac2 =value; }
  else if (!name.compare("AWBM_BFLOW_INDEX"       )){S.AWBM_bflow_index =value; }
  else if (!name.compare("LAKE_REL_COEFF"         )){S.lake_rel_coeff =value;}
  else if (!name.compare("DEP_K"                  )){S.dep_k =value;}
  else if (!name.compare("DEP_SEEP_K"             )){S.dep_seep_k =value;}
  else if (!name.compare("ABST_PERCENT"           )){S.abst_percent =value;}
  else if (!name.compare("OW_PET_CORR"            )){S.ow_PET_corr=value;}
  else if (!name.compare("LAKE_PET_CORR"          )){S.lake_PET_corr=value;}
  else if (!name.compare("FOREST_PET_CORR"        )){S.forest_PET_corr=value;}
  else if (!name.compare("PRIESTLEYTAYLOR_COEFF"  )){S.priestleytaylor_coeff=value;}
  else if (!name.compare("PET_LIN_COEFF"          )){S.pet_lin_coeff=value;}
  else if (!name.compare("GR4J_X4"                )){S.GR4J_x4=value;}
  else if (!name.compare("UBC_ICEPT_FACTOR"       )){S.UBC_icept_factor=value;}
  else if (!name.compare("WIND_EXPOSURE"          )){S.wind_exposure=value;}
  else if (!name.compare("FETCH"                  )){S.fetch=value;}
  else if (!name.compare("AET_COEFF"              )){S.AET_coeff=value;}
  else if (!name.compare("GAMMA_SCALE"            )){S.gamma_scale=value;}
  else if (!name.compare("GAMMA_SHAPE"            )){S.gamma_shape=value;}
  else if (!name.compare("GAMMA_SCALE2"           )){S.gamma_scale2=value;}
  else if (!name.compare("GAMMA_SHAPE2"           )){S.gamma_shape2=value;}
  else if (!name.compare("HMETS_RUNOFF_COEFF"     )){S.HMETS_runoff_coeff=value; }
  else if (!name.compare("BSNOW_DISTRIB"          )){S.bsnow_distrib=value;}
  else if (!name.compare("SKY_VIEW_FACTOR"        )){S.sky_view_factor=value;}
  else if (!name.compare("CONVECTION_COEFF"       )){S.convection_coeff=value;}
  else if (!name.compare("GEOTHERMAL_GRAD"        )){S.geothermal_grad=value;}
  else if (!name.compare("MIN_WIND_SPEED"         )){S.min_wind_speed=value;}
  else if (!name.compare("MAX_WIND_SPEED"         )){S.max_wind_speed=value;}
  else if (!name.compare("STREAM_FRACTION"        )){S.stream_fraction=value;}  
  else{
    WriteWarning("Trying to set value of unrecognized/invalid land use/land type parameter "+ name,false);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief gets surface property corresponding to param_name
/// \param param_name [in] Parameter identifier
/// \returns value of parameter
//
double CLandUseClass::GetSurfaceProperty(string param_name) const
{
  return GetSurfaceProperty(S,param_name);
}

///////////////////////////////////////////////////////////////////////////
/// \brief Returns land surface property value corresponding to param_name from structure provided
/// \param surface_struct [in] land surface structure
/// \param param_name [in] Parameter name
/// \return LULT property corresponding to parameter name
//
double CLandUseClass::GetSurfaceProperty(const surface_struct &S, string param_name, const bool strict)
{
  string name;
  name = StringToUppercase(param_name);

  /*for(i=0;i<N_LU_PARAMETERS;i++) {
  if (!name.compare(S.params[i].name)){return S.params[i].value;}
  }*/ // \todo[funct] - PARAMETEROVERHAUL (replaces below)

  if      (!name.compare("IMPERMEABLE_FRAC"       )){return S.impermeable_frac;}
  else if (!name.compare("FOREST_COVERAGE"        )){return S.forest_coverage;}
  else if (!name.compare("ROUGHNESS"              )){return S.roughness;}

  else if (!name.compare("FOREST_SPARSENESS"      )){return S.forest_sparseness;}
  else if (!name.compare("MELT_FACTOR"            )){return S.melt_factor;}
  else if (!name.compare("MIN_MELT_FACTOR"        )){return S.min_melt_factor;}
  else if (!name.compare("DD_MELT_TEMP"           )){return S.DD_melt_temp;}
  else if (!name.compare("MAX_MELT_FACTOR"        )){return S.max_melt_factor;}
  else if (!name.compare("DD_AGGRADATION"         )){return S.DD_aggradation;}
  else if (!name.compare("REFREEZE_FACTOR"        )){return S.refreeze_factor;}
  else if (!name.compare("DD_REFREEZE_TEMP"       )){return S.DD_refreeze_temp; }
  else if (!name.compare("REFREEZE_EXP"           )){return S.refreeze_exp; }
  else if (!name.compare("HBV_MELT_ASP_CORR"      )){return S.HBV_melt_asp_corr;}
  else if (!name.compare("HBV_MELT_FOR_CORR"      )){return S.HBV_melt_for_corr;}
  else if (!name.compare("MAX_SAT_AREA_FRAC"      )){return S.max_sat_area_frac;}
  else if (!name.compare("HBV_MELT_GLACIER_CORR"  )){return S.HBV_melt_glacier_corr;}
  else if (!name.compare("HBV_GLACIER_KMIN"       )){return S.HBV_glacier_Kmin;}
  else if (!name.compare("SNOW_PATCH_LIMIT"		    )){return S.snow_patch_limit; }
  else if (!name.compare("CONV_MELT_MULT"		      )){return S.conv_melt_mult; }
  else if (!name.compare("COND_MELT_MULT"		      )){return S.cond_melt_mult; }
  else if (!name.compare("RAIN_MELT_MULT"		      )){return S.rain_melt_mult; }
  else if (!name.compare("GLAC_STORAGE_COEFF"     )){return S.glac_storage_coeff;}
  else if (!name.compare("HBV_GLACIER_AG"         )){return S.HBV_glacier_Ag;}
  else if (!name.compare("CC_DECAY_COEFF"         )){return S.CC_decay_coeff;}
  else if (!name.compare("PARTITION_COEFF"        )){return S.partition_coeff;}
  else if (!name.compare("SCS_IA_FRACTION"        )){return S.SCS_Ia_fraction;}
  else if (!name.compare("SCS_CN"                 )){return S.SCS_CN;}
  else if (!name.compare("DEP_MAX"                )){return S.dep_max ;}
  else if (!name.compare("DEP_MAX_FLOW"           )){return S.dep_max_flow;}
  else if (!name.compare("DEP_N"                  )){return S.dep_n;}
  else if (!name.compare("DEP_THRESHHOLD"         )){return S.dep_threshhold;}
  else if (!name.compare("DEP_K"                  )){return S.dep_k;}
  else if (!name.compare("DEP_SEEP_K"             )){return S.dep_seep_k;}
  else if (!name.compare("DEP_CRESTRATIO"         )){return S.dep_crestratio;}
  else if (!name.compare("PDMROF_B"               )){return S.PDMROF_b; }
  else if (!name.compare("PDM_B"                  )){return S.PDM_b; }
  else if (!name.compare("MAX_DEP_AREA_FRAC"      )){return S.max_dep_area_frac; }
  else if (!name.compare("PONDED_EXP"             )){return S.ponded_exp; }
  else if (!name.compare("UWFS_B"                 )){return S.uwfs_b; }
  else if (!name.compare("UWFS_BETAMIN"           )){return S.uwfs_betamin; }
  else if (!name.compare("BF_LOSS_FRACTION"       )){return S.bf_loss_fraction; }
  else if (!name.compare("AWBM_AREAFRAC1"         )){return S.AWBM_areafrac1; }
  else if (!name.compare("AWBM_AREAFRAC2"         )){return S.AWBM_areafrac2; }
  else if (!name.compare("AWBM_BFLOW_INDEX"       )){return S.AWBM_bflow_index; }
  else if (!name.compare("LAKE_REL_COEFF"         )){return S.lake_rel_coeff;}
  else if (!name.compare("ABST_PERCENT"           )){return S.abst_percent;}
  else if (!name.compare("OW_PET_CORR"            )){return S.ow_PET_corr;}
  else if (!name.compare("LAKE_PET_CORR"          )){return S.lake_PET_corr;}
  else if (!name.compare("FOREST_PET_CORR"        )){return S.forest_PET_corr;}
  else if (!name.compare("PRIESTLEYTAYLOR_COEFF"  )){return S.priestleytaylor_coeff;}
  else if (!name.compare("PET_LIN_COEFF"          )){return S.pet_lin_coeff;}
  else if (!name.compare("GR4J_X4"                )){return S.GR4J_x4;}
  else if (!name.compare("UBC_ICEPT_FACTOR"       )){return S.UBC_icept_factor;}
  else if (!name.compare("WIND_EXPOSURE"          )){return S.wind_exposure;}
  else if (!name.compare("FETCH"                  )){return S.fetch;}
  else if (!name.compare("AET_COEFF"              )){return S.AET_coeff;}
  else if (!name.compare("GAMMA_SCALE"            )){return S.gamma_scale;}
  else if (!name.compare("GAMMA_SHAPE"            )){return S.gamma_shape;}
  else if (!name.compare("GAMMA_SCALE2"           )){return S.gamma_scale2;}
  else if (!name.compare("GAMMA_SHAPE2"           )){return S.gamma_shape2;}
  else if (!name.compare("HMETS_RUNOFF_COEFF"     )){return S.HMETS_runoff_coeff; }
  else if (!name.compare("BSNOW_DISTRIB"          )){return S.bsnow_distrib; }
  else if (!name.compare("SKY_VIEW_FACTOR"        )){return S.sky_view_factor; }
  else if (!name.compare("CONVECTION_COEFF"       )){return S.convection_coeff; }
  else if (!name.compare("GEOTHERMAL_GRAD"        )){return S.geothermal_grad; }
  else if (!name.compare("MIN_WIND_SPEED"         )){return S.min_wind_speed; }
  else if (!name.compare("MAX_WIND_SPEED"         )){return S.max_wind_speed; }
  else if (!name.compare("STREAM_FRACTION"        )){return S.stream_fraction; }
  else{
    if (strict){
      string msg="CLandUseClass::GetSurfaceProperty: Unrecognized/invalid LU/LT parameter name in .rvp file: "+name;
      ExitGracefully(msg.c_str(),BAD_DATA_WARN);
      return 0.0;
    }
    else {
      return INDEX_NOT_FOUND;
    }
  }
}
