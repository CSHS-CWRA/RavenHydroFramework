/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2018 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Properties.h"
#include "SoilAndLandClasses.h"

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the vegetation class constructor
/// \param name [in] String nickname for vegetation class
//
CVegetationClass::CVegetationClass(const string name)
{
  this->V.vegetation_name=name;
  if (!DynArrayAppend((void**&)(pAllVegClasses),(void*)(this),NumVegClasses)){
    ExitGracefully("CVegetationClass::Constructor: creating NULL vegetation class",BAD_DATA);};
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
//
CVegetationClass::~CVegetationClass()
{
  if (DESTRUCTOR_DEBUG){cout<<"  DELETING VEGETATION CLASS "<<endl;}
}
/*****************************************************************
   Get Canopy, Root Structures, Tag
*****************************************************************/
//////////////////////////////////////////////////////////////////
/// \brief Returns vegetation properties corresponding to vegetation class
/// \return vegetation properties corresponding to vegetation class
//
const veg_struct *CVegetationClass::GetVegetationStruct() const{return &V;}

///////////////////////////////////////////////////////////////////
/// \brief Return nick name identifier of vegetation class
/// \return nick name identifier of vegetation class
//
string               CVegetationClass::GetVegetationName         () const{return V.vegetation_name;}

/*****************************************************************
   Static Initialization, Accessors, Destructors
*****************************************************************/
CVegetationClass **CVegetationClass::pAllVegClasses=NULL;
int                CVegetationClass::NumVegClasses=0;

//////////////////////////////////////////////////////////////////
/// \brief Return number of vegetation classes
/// \return Number of vegetation classes
//
int CVegetationClass::GetNumClasses(){
  return NumVegClasses;
}

//////////////////////////////////////////////////////////////////
/// \brief Summarize vegetation class information to screen
//
void CVegetationClass::SummarizeToScreen()
{
  cout<<"==================="<<endl;
  cout<<"Vegetation Class Summary:"<<NumVegClasses<<" vegetation classes in database"<<endl;
  for (int c=0; c<NumVegClasses;c++){
    cout<<"-Veg. class \""<<pAllVegClasses[c]->GetVegetationName()<<"\" "<<endl;
    cout<<"    max. height: "<<pAllVegClasses[c]->GetVegetationStruct()->max_height<<" m"<<endl;
    cout<<"       max. LAI: "<<pAllVegClasses[c]->GetVegetationStruct()->max_LAI  <<endl;
    cout<<"   max. conduct: "<<pAllVegClasses[c]->GetVegetationStruct()->max_leaf_cond<<" mm/s"<<endl;
    cout<<"         albedo: "<<pAllVegClasses[c]->GetVegetationStruct()->albedo<<endl;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Destroy all vegetation classes
//
void CVegetationClass::DestroyAllVegClasses()
{
  if (DESTRUCTOR_DEBUG){cout <<"DESTROYING ALL VEGETATION CLASSES"<<endl;}
  for (int c=0; c<NumVegClasses;c++){
    delete pAllVegClasses[c];
  }
  delete [] pAllVegClasses;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns the vegetation class corresponding to passed string
/// \details Converts string (e.g., "BROADLEAF" in HRU file) to vegetation class
///  can accept either vegclass index or vegclass tag
///  if string is invalid, returns NULL
/// \param s [in] Vegetation class identifier (tag or index)
/// \return Reference to vegetation class corresponding to identifier s
//
CVegetationClass *CVegetationClass::StringToVegClass(const string s)
{
  string sup;
  sup=StringToUppercase(s);
  for (int c=0;c<NumVegClasses;c++)
  {
    if (!sup.compare(StringToUppercase(pAllVegClasses[c]->GetVegetationName()))){return pAllVegClasses[c];}
    else if (s_to_i(sup.c_str())==(c+1))                                        {return pAllVegClasses[c];}
  }
  return NULL;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns the vegetation class corresponding to the passed index
///  if index is invalid, returns NULL
/// \param c [in] Soil class index
/// \return Reference to vegetation class corresponding to index c
//
const CVegetationClass *CVegetationClass::GetVegClass(int c)
{
  if ((c<0) || (c>=NumVegClasses)){return NULL;}
  return pAllVegClasses[c];
}
/*****************************************************************
   Accessors
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Automatically calculates vegetation propeties
/// \details  Sets vegetation properties based upon simple vegetation parameters
///      Input [Ctmp & Rtmp] has been read from .rvp file - if parameter ==
///      AUTO_COMPLETE, then empirical relationships are used to estimate
///      physical parameters
///
/// \param &Ctmp [in] Input canopy properties
/// \param &Rtmp [in] Input root properties
/// \param &Cdefault [in] Default canopy properties
/// \param &Rdefault [in] Default root properties
//
string vegautolist[]={"SVF_EXTINCTION","ALBEDO","ALBEDO_WET","RAIN_ICEPT_FACT","SNOW_ICEPT_FACT"};
void CVegetationClass::AutoCalculateVegetationProps(const veg_struct &Vtmp, const veg_struct &Vdefault)
{
  int mon;
  bool autocalc;
  string warn;
  bool chatty(true);

  //these parameters are required
  V.vegetation_name=Vtmp.vegetation_name;
  V.max_LAI        =Vtmp.max_LAI;
  V.max_leaf_cond  =Vtmp.max_leaf_cond;
  V.max_height     =Vtmp.max_height;
  
  //Standard vegetation properties
  //----------------------------------------------------------------------------
  //Skyview extinction coefficient
  autocalc=SetCalculableValue(V.svf_extinction,Vtmp.svf_extinction,Vdefault.svf_extinction);
  if (autocalc)
  {
    //RELATIONSHIP REQUIRED
    V.svf_extinction=0.5;
    warn="The required parameter SVF_EXTINCTION for vegetation class "+V.vegetation_name+" was autogenerated with value "+to_string(V.svf_extinction);
    if(chatty){WriteAdvisory(warn,false);}
  }
  //Albedos
  autocalc=SetCalculableValue(V.albedo,Vtmp.albedo,Vdefault.albedo);
  if (autocalc)
  {
    //RELATIONSHIP REQUIRED
    V.albedo=0.14;
    warn="The required parameter ALBEDO for vegetation class "+V.vegetation_name+" was autogenerated with value "+to_string(V.albedo);
    if(chatty){WriteAdvisory(warn,false);}
  }
  autocalc=SetCalculableValue(V.albedo_wet,Vtmp.albedo_wet,Vdefault.albedo_wet);
  if (autocalc)
  {
    V.albedo_wet=V.albedo+0.04; ///< From WATCLASS CLASSA       \cite verseghy2008
    warn="The required parameter ALBEDO_WET for vegetation class "+V.vegetation_name+" was autogenerated with value "+to_string(V.albedo_wet);
    if(chatty){WriteAdvisory(warn,false);}
  }
  //Interception factors
  autocalc=SetCalculableValue(V.rain_icept_fact,Vtmp.rain_icept_fact,Vdefault.rain_icept_fact);
  if (autocalc)
  {
    V.rain_icept_fact=0.06;     ///< From Dingman       \cite Dingman2002
    warn="The required parameter RAIN_ICEPT_FACT for vegetation class "+V.vegetation_name+" was autogenerated with value "+to_string(V.rain_icept_fact);
    if(chatty){WriteAdvisory(warn,false);}
  }
  autocalc=SetCalculableValue(V.snow_icept_fact,Vtmp.snow_icept_fact,Vdefault.snow_icept_fact);
  if (autocalc)
  {
    V.snow_icept_fact=0.04;     ///< From Dingman       \cite Dingman2002
    warn="The required parameter SNOW_ICEPT_FACT for vegetation class "+V.vegetation_name+" was autogenerated with value "+to_string(V.snow_icept_fact);
    if(chatty){WriteAdvisory(warn,false);}
  }
  autocalc=SetCalculableValue(V.rain_icept_pct,Vtmp.rain_icept_pct,Vdefault.rain_icept_pct);
  if (autocalc)
  {
    V.rain_icept_pct=0.12;      //reasonable default
    warn="The required parameter RAIN_ICEPT_PCT for vegetation class "+V.vegetation_name+" was autogenerated with value "+to_string(V.rain_icept_pct);
    if(chatty){WriteAdvisory(warn,false);}
  }
  ExitGracefullyIf(V.rain_icept_pct>1.0,"Invalid value for RAIN_ICEPT_PCT (>1.0)",BAD_DATA_WARN);

  autocalc=SetCalculableValue(V.snow_icept_pct,Vtmp.snow_icept_pct,Vdefault.snow_icept_pct);
  if (autocalc)
  {
    V.snow_icept_pct=0.10;      //reasonable default
    warn="The required parameter SNOW_ICEPT_PCT for vegetation class "+V.vegetation_name+" was autogenerated with value "+to_string(V.snow_icept_pct);
    if(chatty){WriteAdvisory(warn,false);}
  }
  ExitGracefullyIf(V.snow_icept_pct>1.0,"Invalid value for SNOW_ICEPT_PCT (>1.0)",BAD_DATA_WARN);

  //Ratio of SAI to vegetation height
  autocalc=SetCalculableValue(V.SAI_ht_ratio,Vtmp.SAI_ht_ratio,Vdefault.SAI_ht_ratio);
  if (autocalc)
  {
    V.SAI_ht_ratio=0.0; //Default is that SAI is ignored
    //no warning - default is zero
  }
  //Trunk Fraction
  autocalc=SetCalculableValue(V.trunk_fraction,Vtmp.trunk_fraction,Vdefault.trunk_fraction);
  if (autocalc)
  {
    V.trunk_fraction=0.0;
    //no warning - default is zero
  }
  //Stemflow percentage
  autocalc=SetCalculableValue(V.stemflow_frac,Vtmp.stemflow_frac,Vdefault.stemflow_frac);
  if (autocalc)
  {
    V.stemflow_frac=0.0;//0.03;
    //no warning - default is zero
  }
  //PET vegetation correction
  autocalc=SetCalculableValue(V.PET_veg_corr,Vtmp.PET_veg_corr,Vdefault.PET_veg_corr);
  if (autocalc)
  {
    V.PET_veg_corr=1.0;
    //no warning - default is one
  }
    

  double max_LAI=V.max_LAI;
  double max_SAI=V.SAI_ht_ratio*V.max_height;

  //Canopy Capacity [mm]
  autocalc=SetCalculableValue(V.max_capacity,Vtmp.max_capacity,Vdefault.max_capacity);
  if (autocalc)
  {
    V.max_capacity=CAP_LAI_RATIO*(max_LAI+max_SAI);     ///< maximum storage capacity [mm] (Brook90 / Dingman 7-2) \cite Federer2010
    warn="The required parameter MAX_CAPACITY for vegetation class "+V.vegetation_name+" was autogenerated with value "+to_string(V.max_capacity);
    if (chatty){WriteAdvisory(warn,false);}
  }
  //Canopy Snow Capacity
  autocalc=SetCalculableValue(V.max_snow_capacity,Vtmp.max_snow_capacity,Vdefault.max_snow_capacity);
  if (autocalc)
  {
    V.max_snow_capacity=SCAP_LAI_RATIO*(max_LAI+max_SAI);       //Brook90/Dingman box 5.1/CLM
    ///< V.max_snow_capacity=6.0 *(0.27+46*0.119)*(V.max_LAI+maxSAI);//WATCLASS (Schmidt & Gluns, 1991)[units seem off] \cite schmidt1991CJFR
    warn="The required parameter MAX_SNOW_CAPACITY for vegetation class "+V.vegetation_name+" was autogenerated with value "+to_string(V.max_snow_capacity);
    if (chatty){WriteAdvisory(warn,false);}
  }
  //Canopy Maximum Snow Load
  autocalc=SetCalculableValue(V.max_snow_load,Vtmp.max_snow_load,Vdefault.max_snow_load);
  if(autocalc)
  {
    V.max_snow_load=6.0;       //Hedstrom and Pomeroy
    warn="The required parameter MAX_SNOW_LOAD for vegetation class "+V.vegetation_name+" was autogenerated with value "+to_string(V.max_snow_load);
    if(chatty) { WriteAdvisory(warn,false); }
  }

  //Relative canopy height
  for (mon=0;mon<12;mon++)
  {
    autocalc=SetCalculableValue(V.relative_ht[mon],Vtmp.relative_ht[mon],Vdefault.relative_ht[mon]);
    if (autocalc)
    {
      V.relative_ht[mon]=1.0;//reasonable (no warning - zero correction factor)
    }
  }


  //Relative LAI
  bool is_auto=false;
  for (mon=0;mon<12;mon++)
  {
    autocalc=SetCalculableValue(V.relative_LAI[mon],Vtmp.relative_LAI[mon],Vdefault.relative_LAI[mon]);
    if (autocalc)
    {
      V.relative_LAI[mon]=1.0;//0.5*(1.0+sin((mon+1.5)/12.0*2*PI));//peaks in mid-aug, lowest in mid-feb
      is_auto=true;
    }
  }
  if (is_auto){
    warn="The required array parameter RELATIVE_LAI[] for vegetation class "+V.vegetation_name+" was autogenerated with values of 1.0";
    WriteAdvisory(warn,false);
  }



  // Transport Parameters (default to zero)
  //----------------------------------------------------------------------------
  for(int c=0; c<MAX_CONSTITUENTS;c++){
    autocalc=SetCalculableValue(V.uptake_moderator[c],Vtmp.uptake_moderator[c],Vdefault.uptake_moderator[c]);
    if(autocalc)
    {
      V.uptake_moderator[c]=1.0; //reasonable default
    }
  }

  //Model-specific vegetation properties - cannot be autocomputed, must be specified by user
  //----------------------------------------------------------------------------
  bool needed=false;
  bool bad;
  bad=SetSpecifiedValue(V.drip_proportion,Vtmp.drip_proportion,Vdefault.drip_proportion,needed,"DRIP_PROPORTION");
  bad=SetSpecifiedValue(V.max_intercept_rate,Vtmp.max_intercept_rate,Vdefault.max_intercept_rate,needed,"MAX_INTERCEPT_RATE");
  bad=SetSpecifiedValue(V.CHU_maturity,Vtmp.CHU_maturity,Vdefault.CHU_maturity,needed,"CHU_MATURITY");
  bad=SetSpecifiedValue(V.root_extinct,Vtmp.root_extinct,Vdefault.root_extinct,needed,"ROOT_EXTINCT");
  bad=SetSpecifiedValue(V.max_root_length,Vtmp.max_root_length,Vdefault.max_root_length,needed,"MAX_ROOT_LENGTH");
  bad=SetSpecifiedValue(V.min_resistivity,Vtmp.min_resistivity,Vdefault.min_resistivity,needed,"MIN_RESISTIVITY");
  bad=SetSpecifiedValue(V.psi_critical,Vtmp.psi_critical,Vdefault.psi_critical,needed,"PSI_CRITICAL");
  bad=SetSpecifiedValue(V.rootradius,Vtmp.rootradius,Vdefault.rootradius,needed,"ROOTRADIUS");
  bad=SetSpecifiedValue(V.xylem_frac,Vtmp.xylem_frac,Vdefault.xylem_frac,needed,"XYLEM_FRAC");
  bad=SetSpecifiedValue(V.veg_dens,Vtmp.veg_dens,Vdefault.veg_dens,needed,"VEG_DENS");
  bad=SetSpecifiedValue(V.veg_diam,Vtmp.veg_diam,Vdefault.veg_diam,needed,"VEG_DIAM");
  bad=SetSpecifiedValue(V.veg_mBeta,Vtmp.veg_mBeta,Vdefault.veg_mBeta,needed,"VEG_MBETA");
  bad=SetSpecifiedValue(V.veg_conv_coeff,Vtmp.veg_conv_coeff,Vdefault.veg_conv_coeff,needed,"VEG_CONV_COEFF");
  if(bad) {
    //do something.
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets default vegetation properties
/// \details Sets vegetation properties based upon vegetation type
/// \param &V [out] Parameter values to be set
/// \param is_template [in] True if the default value being set is for the template class
//
void CVegetationClass::InitializeVegetationProps(string name, veg_struct &V, bool is_template)
{
  V.vegetation_name =name;

  //Required parameters
  V.max_height    =25;
  V.max_leaf_cond =5.3; //From broadleaf, Dingman
  V.max_LAI       =6.0;

  //Autocalculable parameters
  V.svf_extinction    =DefaultParameterValue(is_template,true);//0.5

  V.albedo            =DefaultParameterValue(is_template,true);//0.14;
  V.albedo_wet        =DefaultParameterValue(is_template,true);//0.18;

  V.rain_icept_fact    =DefaultParameterValue(is_template,true);
  V.snow_icept_fact    =DefaultParameterValue(is_template,true);
  V.rain_icept_pct     =DefaultParameterValue(is_template,true);
  V.snow_icept_pct     =DefaultParameterValue(is_template,true);

  for (int mon=0;mon<12;mon++)
  {
    V.relative_ht [mon]=DefaultParameterValue(is_template,true);
    V.relative_LAI[mon]=DefaultParameterValue(is_template,true);
  }

  V.trunk_fraction    =DefaultParameterValue(is_template,true);
  V.stemflow_frac     =DefaultParameterValue(is_template,true);
  V.SAI_ht_ratio      =DefaultParameterValue(is_template,true);

  V.max_capacity      =DefaultParameterValue(is_template,true);
  V.max_snow_capacity =DefaultParameterValue(is_template,true);
  V.max_snow_load     =DefaultParameterValue(is_template,true);

  V.PET_veg_corr=1.0;

  //No autocalculation here (OR SUFFICIENT WARNING IF PARAMETER NOT SPECIFIED) \todo [funct]
  for(int c=0;c<MAX_CONSTITUENTS;c++){
    V.uptake_moderator[c]                       =1.0; //1.0
  }

  //User-specified parameters
  V.drip_proportion   =DefaultParameterValue(is_template,false);//0.2;//[1/day]
  V.max_intercept_rate=DefaultParameterValue(is_template,false);//10.0;//[mm/day]
  V.CHU_maturity      =DefaultParameterValue(is_template,false);//2800;
  V.max_root_length   =DefaultParameterValue(is_template,false);//200;  //[mm/m2]
  V.min_resistivity   =DefaultParameterValue(is_template,false);//1.0;  //[d/mm]
  V.xylem_frac        =DefaultParameterValue(is_template,false);//0.01; //fraction of plant resistance in xylem
  V.rootradius        =DefaultParameterValue(is_template,false);//4;    //[mm]
  V.psi_critical      =DefaultParameterValue(is_template,false);//0.0;  //[-mm] minimum plant leaf water potential
  V.root_extinct      =DefaultParameterValue(is_template,false);        //[-]
  V.veg_dens          =DefaultParameterValue(is_template,false);        //[/m2]
  V.veg_diam          =DefaultParameterValue(is_template,false);        //[m]
  V.veg_mBeta         =DefaultParameterValue(is_template,false);        //[-]
  V.veg_conv_coeff    =DefaultParameterValue(is_template,false);//2.0   //[MJ/m2/d/K]
}
////////////////////////////////////////////////////////////////////
/// \brief Sets vegetation property
/// \param param_name [in] Identifier of parameter to be set
/// \param value [in] Value of parameter to be set
//
void CVegetationClass::SetVegetationProperty(const string &param_name, const double &value)
{
  SetVegetationProperty(V,param_name,value);
}

////////////////////////////////////////////////////////////////////
/// \brief Sets vegetation property
/// \details Sets a single parameter associated with vegetation class
/// \param &V [out] Reference to vegetation properties associated with vegetation class
/// \param param_name [in] Identifier of parameter to be set
/// \param value [in] Value of parameter to be set
//
void  CVegetationClass::SetVegetationProperty(veg_struct  &V,
                                              const string param_name,
                                              const double &value)
{
  string name;
  name = StringToUppercase(param_name);

  //Canopy params
  if      (!name.compare("MAX_HEIGHT"           )){V.max_height=value;}
  else if (!name.compare("MAX_LEAF_COND"        )){V.max_leaf_cond=value;}
  else if (!name.compare("MAX_LAI"              )){V.max_LAI=value;}

  else if (!name.compare("SVF_EXTINCTION"       )){V.svf_extinction=value;}
  else if (!name.compare("ALBEDO"               )){V.albedo=value;}
  else if (!name.compare("ALBEDO_WET"           )){V.albedo_wet=value;}
  else if (!name.compare("RAIN_ICEPT_FACT"      )){V.rain_icept_fact=value;}
  else if (!name.compare("SNOW_ICEPT_FACT"      )){V.snow_icept_fact=value;}
  else if (!name.compare("TRUNK_FRACTION"       )){V.trunk_fraction=value;}
  else if (!name.compare("STEMFLOW_FRAC"        )){V.stemflow_frac=value;}
  else if (!name.compare("SAI_HT_RATIO"         )){V.SAI_ht_ratio=value;}
  else if (!name.compare("MAX_CAPACITY"         )){V.max_capacity=value;}
  else if (!name.compare("MAX_SNOW_CAPACITY"    )){V.max_snow_capacity=value;}
  else if (!name.compare("MAX_SNOW_LOAD"        )){V.max_snow_load=value;}
  else if (!name.compare("TFRAIN"               )){V.rain_icept_pct=1.0-value;}
  else if (!name.compare("TFSNOW"               )){V.snow_icept_pct=1.0-value;}
  else if (!name.compare("RAIN_ICEPT_PCT"       )){V.rain_icept_pct=value;}
  else if (!name.compare("SNOW_ICEPT_PCT"       )){V.snow_icept_pct=value;}
  else if (!name.compare("DRIP_PROPORTION"      )){V.drip_proportion=value;}
  else if (!name.compare("MAX_INTERCEPT_RATE"   )){V.max_intercept_rate=value;}
  else if (!name.compare("CHU_MATURITY"         )){V.CHU_maturity=value;}
  else if (!name.compare("VEG_DIAM"             )){V.veg_diam=value;}
  else if (!name.compare("VEG_MBETA"            )){V.veg_mBeta=value;}
  else if (!name.compare("VEG_DENS"             )){V.veg_dens=value;}
  else if (!name.compare("PET_VEG_CORR"         )){V.PET_veg_corr=value;}
  else if (!name.compare("VEG_CONV_COEFF"       )){V.veg_conv_coeff=value;}

  else if (!name.compare("MAX_ROOT_LENGTH"      )){V.max_root_length=value;}
  else if (!name.compare("MIN_RESISTIVITY"      )){V.min_resistivity=value;}
  else if (!name.compare("XYLEM_FRAC"           )){V.xylem_frac=value;}
  else if (!name.compare("ROOTRADIUS"           )){V.rootradius=value;}
  else if (!name.compare("PSI_CRITICAL"         )){V.psi_critical=value;}
  else if (!name.compare("ROOT_EXTINCT"         )){V.root_extinct=value;}

  else if (!name.compare("RELATIVE_HT"          )){for (int mon=0;mon<12;mon++){V.relative_ht [mon]=value;}}//special case
  else if (!name.compare("RELATIVE_LAI"         )){for (int mon=0;mon<12;mon++){V.relative_LAI[mon]=value;}}//special case
  else{
    WriteWarning("Trying to set value of unrecognized/invalid vegetation parameter \""+name+"\"",true);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Sets the value of the vegetation property corresponding to param_name
/// \note This is declared as a static member because vegetation class
/// is not instantiated prior to read of .rvp file
/// \param constit_ind [in] consitutent index, c
/// \param constit_ind [in] constituent index of second constituent (optional - needed for some parameters)
/// \param &V [out] veg properties class
/// \param param_name [in] Parameter identifier
/// \param value [in] Value of parameter to be set
//
void  CVegetationClass::SetVegTransportProperty(int          constit_ind,
                                                int          constit_ind2,
                                                veg_struct  &V,
                                                string       param_name,
                                                const double value)
{
  string name;
  name = StringToUppercase(param_name);

  if      (!name.compare("UPTAKE_MODIFIER"         )){V.uptake_moderator   [constit_ind]=value;}
  else{
    WriteWarning("CVegetationClass::SetVegTransportProperty: Unrecognized/invalid soil parameter name ("+name+") in .rvp file",false);
  }
}
////////////////////////////////////////////////////////////////////
/// \brief Returns vegetation property specified by parameter
/// \param param_name [in] Identifier of parameter to be returned
/// \return Value corresponding to param_name
//
double CVegetationClass::GetVegetationProperty(string param_name) const
{
  return GetVegetationProperty(V,param_name);
}
///////////////////////////////////////////////////////////////////////////
/// \brief Returns veg transport property value corresponding to param_name from structure provided
/// \param constit_ind [in] constituent index c
/// \param veg_struct [in] Soil structure
/// \param param_name [in] Parameter name
/// \return veg property corresponding to parameter name
//
double CVegetationClass::GetVegTransportProperty(const int constit_ind, const veg_struct &V, string &param_name)
{
  string name;
  name = StringToUppercase(param_name);

  if      (!name.compare("UPTAKE_MODERATOR"         )){return V.uptake_moderator[constit_ind];}
  else{
    string msg="CVegetationClass::GetVegTransportProperty: Unrecognized/invalid soil parameter name  ("+name+") in .rvp file: "+name;
    ExitGracefully(msg.c_str(),BAD_DATA);
    return 0.0;
  }
}
///////////////////////////////////////////////////////////////////////////
/// \brief Returns vegetation property value corresponding to param_name from structure provided
/// \param veg_struct [in] Vegetation structure
/// \param param_name [in] Parameter name
/// \return Vegetation property corresponding to parameter name
//
double CVegetationClass::GetVegetationProperty(const veg_struct &V, string param_name, const bool strict)
{
  string name;
  name = StringToUppercase(param_name);

  //Canopy params
  if      (!name.compare("MAX_HEIGHT"           )){return V.max_height;}
  else if (!name.compare("MAX_LEAF_COND"        )){return V.max_leaf_cond;}
  else if (!name.compare("MAX_LAI"              )){return V.max_LAI;}

  else if (!name.compare("SVF_EXTINCTION"       )){return V.svf_extinction;}
  else if (!name.compare("ALBEDO"               )){return V.albedo;}
  else if (!name.compare("ALBEDO_WET"           )){return V.albedo_wet;}
  else if (!name.compare("RAIN_ICEPT_FACT"      )){return V.rain_icept_fact;}
  else if (!name.compare("SNOW_ICEPT_FACT"      )){return V.snow_icept_fact;}
  else if (!name.compare("RAIN_ICEPT_PCT"       )){return V.rain_icept_pct;}
  else if (!name.compare("SNOW_ICEPT_PCT"       )){return V.snow_icept_pct;}
  else if (!name.compare("TRUNK_FRACTION"       )){return V.trunk_fraction;}
  else if (!name.compare("STEMFLOW_FRAC"        )){return V.stemflow_frac;}
  else if (!name.compare("SAI_HT_RATIO"         )){return V.SAI_ht_ratio;}
  else if (!name.compare("MAX_CAPACITY"         )){return V.max_capacity;}
  else if (!name.compare("MAX_SNOW_CAPACITY"    )){return V.max_snow_capacity;}
  else if (!name.compare("MAX_SNOW_LOAD"        )){return V.max_snow_load;}

  else if (!name.compare("DRIP_PROPORTION"      )){return V.drip_proportion;}
  else if (!name.compare("MAX_INTERCEPT_RATE"   )){return V.max_intercept_rate;}
  else if (!name.compare("CHU_MATURITY"         )){return V.CHU_maturity;}
  else if (!name.compare("VEG_DENS"             )){return V.veg_dens;}
  else if (!name.compare("VEG_MBETA"            )){return V.veg_mBeta;}
  else if (!name.compare("VEG_DIAM"             )){return V.veg_diam;}
  else if (!name.compare("PET_VEG_CORR"         )){return V.PET_veg_corr;}
  else if (!name.compare("VEG_CONV_COEFF"       )){return V.veg_conv_coeff;}

  else if (!name.compare("MAX_ROOT_LENGTH"      )){return V.max_root_length;}
  else if (!name.compare("MIN_RESISTIVITY"      )){return V.min_resistivity;}
  else if (!name.compare("XYLEM_FRAC"           )){return V.xylem_frac;}
  else if (!name.compare("ROOTRADIUS"           )){return V.rootradius;}
  else if (!name.compare("PSI_CRITICAL"         )){return V.psi_critical;}
  else if (!name.compare("ROOT_EXTINCT"         )){return V.root_extinct;}

  else if (!name.compare("RELATIVE_HT"          )){return V.relative_ht[0];} //special case: only used for autocalc
  else if (!name.compare("RELATIVE_LAI"         )){return V.relative_LAI[0];} //special case: only used for autocalc

  else{
    if (strict){
      string msg="CVegetationClass::GetVegetationProperty: Unrecognized/invalid vegetation parameter name: "+name;
      ExitGracefully(msg.c_str(),BAD_DATA_WARN);
      return 0.0;
    }
    else {
      return INDEX_NOT_FOUND;
    }
  }
}


