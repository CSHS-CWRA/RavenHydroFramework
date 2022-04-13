/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2022 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Properties.h"
#include "SoilAndLandClasses.h"

/*****************************************************************
   Constructor / Destructor
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the soil class constructor
/// \param name [in] String nickname for soil (e.g., "SILTY_SAND")
//
CSoilClass::CSoilClass(const string name,const int nConstit)
{
  _tag=name;
  if (!DynArrayAppend((void**&)(_pAllSoilClasses),(void*)(this),_nAllSoilClasses)){
    ExitGracefully("CSoilClass::Constructor: creating NULL soil class",BAD_DATA);};
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
//
CSoilClass::~CSoilClass(){
  if (DESTRUCTOR_DEBUG){cout<<"  DELETING SOIL CLASS "<<endl;}
}

///////////////////////////////////////////////////////////////////
/// \brief Returns reference to soil properties
/// \return Soil properties associated with soil class
//
const soil_struct *CSoilClass::GetSoilStruct() const{return &_Soil;}

///////////////////////////////////////////////////////////////////
/// \brief Return nick name identifier of soil class
/// \return nick name identifier of soil class
//
string             CSoilClass::GetTag       () const{return _tag;}

/*****************************************************************
   Static Initialization, Accessors, Destructors
*****************************************************************/
CSoilClass **CSoilClass::_pAllSoilClasses=NULL;
int          CSoilClass::_nAllSoilClasses=0;


//////////////////////////////////////////////////////////////////
/// \brief Return number of soil classes
/// \return Number of soil classes
//
int CSoilClass::GetNumClasses(){
  return _nAllSoilClasses;
}


//////////////////////////////////////////////////////////////////
/// \brief Summarize soil class information to screen
//
void CSoilClass::SummarizeToScreen()
{
  cout<<"==================="<<endl;
  cout<<"Soil Class Summary:"<<_nAllSoilClasses<<" soils in database"<<endl;
  for (int c=0; c<_nAllSoilClasses;c++){
    cout<<"-Soil class \""<<_pAllSoilClasses[c]->GetTag()<<"\" "<<endl;
    cout<<"       %sand: "<<_pAllSoilClasses[c]->GetSoilStruct()->sand_con<<endl;
    cout<<"       %clay: "<<_pAllSoilClasses[c]->GetSoilStruct()->clay_con<<endl;
    cout<<"    %organic: "<<_pAllSoilClasses[c]->GetSoilStruct()->org_con<<endl;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Destroy all soil classes
//
void CSoilClass::DestroyAllSoilClasses()
{
  if (DESTRUCTOR_DEBUG){cout <<"DESTROYING ALL SOIL CLASSES"<<endl;}
  for (int c=0; c<_nAllSoilClasses;c++){
    delete _pAllSoilClasses[c];
  }
  delete [] _pAllSoilClasses;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns the soil class corresponding to passed string
/// \details Converts string (e.g., "SILT" in HRU file) to soilclass
///  can accept either soilclass index or soilclass _tag
///  if string is invalid, returns NULL
/// \param s [in] Soil class identifier (_tag or index)
/// \return Reference to soil class corresponding to identifier s
//
CSoilClass *CSoilClass::StringToSoilClass(const string s)
{
  string sup=StringToUppercase(s);
  for (int c=0;c<_nAllSoilClasses;c++)
  {
    if (!sup.compare(StringToUppercase(_pAllSoilClasses[c]->GetTag()))){return _pAllSoilClasses[c];}
    else if (s_to_i(s.c_str())==(c+1))                                 {return _pAllSoilClasses[c];}
  }
  return NULL;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns the soil class corresponding to the passed index
///  if index is invalid, returns NULL
/// \param c [in] Soil class index
/// \return Reference to soil class corresponding to index c
//
const CSoilClass *CSoilClass::GetSoilClass(int c)
{
  if ((c<0) || (c>=CSoilClass::_nAllSoilClasses)){return NULL;}
  return _pAllSoilClasses[c];
}
//////////////////////////////////////////////////////////////////
/// \brief Automatically calculates soil propeties
/// \details Sets soil properties based upon simple soil parameters (sand,clay
///  &org pct) using simple pedotransfer functions. Input [Stmp] has been
///  read from .rvp file - if parameter == AUTO_COMPLETE, then empirical
///  relationships are used to estimate parameter from sand, clay &
///  organic content
///  \todo [QA/QC] should add some reasonable checks for physical realism (e.g., positivity)
///
/// \param &Stmp [in] Input soil parameters (read from .rvp file)
/// \param &Sdefault [in] Default soil parameters
/// \param nConstit [in] number of constitutents 
//
void CSoilClass::AutoCalculateSoilProps(const soil_struct &Stmp,
                                        const soil_struct &Sdefault,
                                        const int nConstit)
{
  bool autocalc;
  string warn;

  //these parameters are required
  _Soil.sand_con=Stmp.sand_con;
  _Soil.clay_con=Stmp.clay_con;
  _Soil.org_con =Stmp.org_con;
  ExitGracefullyIf((_Soil.sand_con<0) || (_Soil.sand_con>1.0),"AutoCalculateSoilProps: SAND_CON must be between 0 and 1",BAD_DATA_WARN);
  ExitGracefullyIf((_Soil.clay_con<0) || (_Soil.clay_con>1.0),"AutoCalculateSoilProps: CLAY_CON must be between 0 and 1",BAD_DATA_WARN);
  ExitGracefullyIf((_Soil.org_con <0) || (_Soil.org_con >1.0), "AutoCalculateSoilProps: ORG_CON must be between 0 and 1",BAD_DATA_WARN);

  double Vsand=_Soil.sand_con/DENSITY_SAND;
  double Vorg =_Soil.org_con /DENSITY_OM  ;
  double Vclay=_Soil.clay_con/DENSITY_CLAY;
  double Vtot =Vsand+Vclay+Vorg;
  bool chatty(true);

  //empiricism at its finest:

  //Standard soil properties
  //----------------------------------------------------------------------------
  //Porosity
  autocalc=SetCalculableValue(_Soil.porosity,Stmp.porosity,Sdefault.porosity);
  if (autocalc)
  {
    _Soil.porosity=0.489-0.126*_Soil.sand_con;         ///< \ref CLM eqn 7.72, WATCLASS \cite Oleson2012 \cite verseghy2008
    if (_Soil.sand_con == 0.4111111){ _Soil.porosity = 1.0; }///special default for conceptual models without sand content specified
    warn="The required parameter POROSITY for soil class "+_tag+" was autogenerated with value "+to_string(_Soil.porosity);
    if (chatty){WriteAdvisory(warn,false);}
  }
  ExitGracefullyIf((_Soil.porosity<0) || (_Soil.porosity>1.0),"AutoCalculateSoilProps: POROSITY must be between 0 and 1",BAD_DATA);

  //Stone fraction
  autocalc=SetCalculableValue(_Soil.stone_frac,Stmp.stone_frac,Sdefault.stone_frac);
  if (autocalc)
  {
    _Soil.stone_frac=0.0;
    //no warning - default is no correction
  }
  ExitGracefullyIf((_Soil.stone_frac<0) || (_Soil.stone_frac>1.0),"AutoCalculateSoilProps: STONE_FRAC must be between 0 and 1",BAD_DATA);

  //Bulk density
  autocalc=SetCalculableValue(_Soil.bulk_density,Stmp.bulk_density,Sdefault.bulk_density);
  if (autocalc)
  {
    _Soil.bulk_density =DENSITY_SAND*(1-_Soil.porosity);
    warn="The required parameter BULK_DENSITY for soil class "+_tag+" was autogenerated with value "+to_string(_Soil.bulk_density);
    if (chatty){WriteAdvisory(warn,false);}
  }

  _Soil.cap_ratio = _Soil.porosity*(1.0-_Soil.stone_frac);    ///< \ref [-] cap ratio (from Brook90 SOILPAR routine) \cite Dingman2002

  double SandPoro=(1.0-_Soil.porosity )*(Vsand/Vtot);
  double ClayPoro=(1.0-_Soil.porosity )*(Vclay/Vtot);
  double OrgnPoro=(1.0-_Soil.porosity )*(Vorg /Vtot);

  //Thermal properties
  //----------------------------------------------------------------------------
  ///< \ref arithmetic means (e.g., CLM 6.61,6.68, also from WATCLASS) \cite Oleson2012
  autocalc=SetCalculableValue(_Soil.heat_capacity,Stmp.heat_capacity,Sdefault.heat_capacity);
  if (autocalc)
  {
    _Soil.heat_capacity =
      (HCP_SAND   *SandPoro+
       HCP_CLAY   *ClayPoro+
       HCP_ORGANIC*OrgnPoro)/(1.0-_Soil.porosity );
    warn="The required parameter HEAT_CAPACITY for soil class "+_tag+" was autogenerated with value "+to_string(_Soil.heat_capacity);
    if (chatty){WriteAdvisory(warn,false);}
  }

  autocalc=SetCalculableValue(_Soil.thermal_cond,Stmp.thermal_cond,Sdefault.thermal_cond);
  if (autocalc)
  {
    _Soil.thermal_cond  =
      (TC_SAND    *SandPoro+
       TC_CLAY    *ClayPoro+
       TC_ORGANIC *OrgnPoro)/(1.0-_Soil.porosity );
    //geometric mean might also be useful (JRC)
    warn="The required parameter THERMAL_COND for soil class "+_tag+" was autogenerated with value "+to_string(_Soil.thermal_cond);
    if (chatty){WriteAdvisory(warn,false);}
  }

  //Unsaturated flow properties
  //----------------------------------------------------------------------------
  //Hydraulic Conductivity
  autocalc=SetCalculableValue(_Soil.hydraul_cond,Stmp.hydraul_cond,Sdefault.hydraul_cond);
  if (autocalc)
  {
    //RELATIONSHIP REQUIRED
    _Soil.hydraul_cond=0.001;//[mm/d]
    //_Soil.hydraul_cond=pow(10.0,(1.53*_Soil.sand_con-0.884))* 7.0556E-6;//WATCLAS_Soil. Units unknown.
    warn="The required parameter HYDRAUL_COND for soil class "+_tag+" was autogenerated with value "+to_string(_Soil.hydraul_cond);
    if (chatty){WriteAdvisory(warn,false);}
  }

  //Saturation at field capacity
  autocalc=SetCalculableValue(_Soil.field_capacity,Stmp.field_capacity,Sdefault.field_capacity);
  if (autocalc)
  {
    ///< \ref From Brakensiek et al. (1984), as reported in HELP manual \cite Brakensiek1984APP
    _Soil.field_capacity = 0.1535 - 0.18*_Soil.sand_con+ 0.39*_Soil.clay_con+0.1943*_Soil.porosity;
    warn="The required parameter FIELD_CAPACITY for soil class "+_tag+" was autogenerated with value "+to_string(_Soil.field_capacity);
    if (chatty){WriteAdvisory(warn,false);}
  }

  //Wilting point saturation
  autocalc=SetCalculableValue(_Soil.sat_wilt,Stmp.sat_wilt,Sdefault.sat_wilt);
  if (autocalc)
  {
    ///< \ref From Brakensiek et al. (1984), as reported in HELP manual \cite Brakensiek1984APP
    _Soil.sat_wilt =0.037 -0.04*_Soil.sand_con + 0.44*_Soil.clay_con+0.0482*_Soil.porosity;
    warn="The required parameter SAT_WILT for soil class "+_tag+" was autogenerated with value "+to_string(_Soil.sat_wilt);
    if (chatty){WriteAdvisory(warn,false);}
    //_Soil.sat_wilt =0.0; //common for conceptual models
  }

  //Minimum saturation
  autocalc=SetCalculableValue(_Soil.sat_res,Stmp.sat_res,Sdefault.sat_res);
  if (autocalc)
  {
    _Soil.sat_res =0.0;
    //_Soil.sat_res =0.04; //WATCLASS

    ///< From Rawls et al 1982, as reported in HELP manual \cite Rawls1982TA
    //if (_Soil.sat_wilt> 0.04){_Soil.sat_res =0.000+0.60*_Soil.sat_wilt;}
    //else                 {_Soil.sat_res =0.014+0.25*_Soil.sat_wilt;}
    // no warning - default is no residual correction
  }

  //Air entry pressure
  autocalc=SetCalculableValue(_Soil.air_entry_pressure,Stmp.air_entry_pressure,Sdefault.air_entry_pressure);
  if (autocalc)
  {
    _Soil.air_entry_pressure=10*pow(10,1.88-1.31*_Soil.sand_con);//CLM eqn 7.75, WATCLASS
    warn="The required parameter AIR_ENTRY_PRESSURE for soil class "+_tag+" was autogenerated with value "+to_string(_Soil.air_entry_pressure);
    if (chatty){WriteAdvisory(warn,false);}
  }

  //Wilting pressure
  autocalc=SetCalculableValue(_Soil.wilting_pressure,Stmp.wilting_pressure,Sdefault.wilting_pressure);
  if (autocalc)
  {
    _Soil.wilting_pressure =40;
    //RELATIONSHIP REQUIRED
    //calculate from wilting point saturation?
    warn="The required parameter WILTING_PRESSURE for soil class "+_tag+" was autogenerated with value "+to_string(_Soil.wilting_pressure);
    if (chatty){WriteAdvisory(warn,false);}
  }

  ///< \ref Clapp-Hornberger parameters \cite Clapp1978WRR
  // for Clapp-Hornberger parabolic transition at high saturations
  // \math \f$ S = 0.5*(n+1)+sqrt(0.25*(n-1)^2 + psi/m) \f$
  // \math \f$ psi = m * (S - n) * (S - 1) \f$

  autocalc=SetCalculableValue(_Soil.clapp_b,Stmp.clapp_b,Sdefault.clapp_b);
  if (autocalc)
  {
    _Soil.clapp_b=2.91+15.9*_Soil.clay_con;///< CLM eqn 7.73, WATCLASS \cite Oleson2012
    warn="The required parameter CLAPP_B for soil class "+_tag+" was autogenerated with value "+to_string(_Soil.clapp_b);
    if (chatty){WriteAdvisory(warn,false);}
  }
  double psi_inf=_Soil.air_entry_pressure*pow(SAT_INF,-_Soil.clapp_b);
  autocalc=SetCalculableValue(_Soil.clapp_m,Stmp.clapp_m,Sdefault.clapp_m);
  if (autocalc)
  {
    _Soil.clapp_m=(psi_inf/ pow(1.0-SAT_INF,2))- _Soil.clapp_b*(psi_inf)/(SAT_INF*(1-SAT_INF));
    warn="The required parameter CLAPP_M for soil class "+_tag+" was autogenerated with value "+to_string(_Soil.clapp_m);
    if (chatty){WriteAdvisory(warn,false);}
  }
  //Wetting front Matric potential (for Green-Ampt)
  autocalc=SetCalculableValue(_Soil.wetting_front_psi,Stmp.wetting_front_psi,Sdefault.wetting_front_psi);
  if (autocalc)
  {
    ///< \ref SWAT 2005 documentation 2:1.2.5 (from Rawls and Brakensiek, 1985) \cite Rawls1985
    /*
      double mc    =_Soil.clay_con;
      double ms    =_Soil.sand_con;
      double coeff =-7.32561-4.9837*ms;
      coeff+=-3.479*(mc*mc)-7.99*(ms*ms);
      coeff+=(3.809479+16.08*(ms*ms)+16.02*(mc*mc))*_Soil.porosity;
      _Soil.wetting_front_psi=10.0*exp(6.5309+15.83*mc*mc+3.44*ms*mc-13.6*ms*ms*mc+coeff*_Soil.porosity);//[-mm]
    */

    ///< WATFLOOD Documentation (2007), supposedly from Philip, 1954 \cite philip1954SS
    //_Soil.wetting_front_psi=250.0*log(_Soil.hydraul_cond/SEC_PER_DAY)+100.0; //dumb

    //From Neuman integration of Clapp-Hornberger saturation relationship
    double bb  =_Soil.clapp_b;
    double psia=_Soil.air_entry_pressure;
    _Soil.wetting_front_psi=(2*bb+3)/(2*bb+6)*psia;
    warn="The required parameter WETTING_FRONT_PSI for soil class "+_tag+" was autogenerated with value "+to_string(_Soil.wetting_front_psi);
    if (chatty){WriteAdvisory(warn,false);}
    ExitGracefullyIf(_Soil.wetting_front_psi<0,"AutocalculateSoilProperties: wetting front suction must be positive and non-zero.",BAD_DATA_WARN);
  }

  autocalc=SetCalculableValue(_Soil.clapp_n,Stmp.clapp_n,Sdefault.clapp_n);
  if (autocalc)
  {
    _Soil.clapp_n=2*SAT_INF-1.0-(psi_inf*_Soil.clapp_b/(_Soil.clapp_m*SAT_INF));
    warn="The required parameter CLAPP_N for soil class "+_tag+" was autogenerated with value "+to_string(_Soil.clapp_n);
    if (chatty){WriteAdvisory(warn,false);}
  }

  //Lateral Heterogeneity in Hydraulic Conductivity
  autocalc=SetCalculableValue(_Soil.ksat_std_deviation,Stmp.ksat_std_deviation,Sdefault.ksat_std_deviation);
  if (autocalc)
  {
    _Soil.ksat_std_deviation=0.0;
    //no warning - default is zero variability
  }

  //Evaporation properties
  //----------------------------------------------------------------------------
  //Soil evaporation resistance at field capacity
  autocalc=SetCalculableValue(_Soil.evap_res_fc,Stmp.evap_res_fc,Sdefault.evap_res_fc);
  if (autocalc)
  {
    _Soil.evap_res_fc=1.0;//?
    //RELATIONSHIP REQUIRED
    warn="The required parameter EVAP_RES_FC for soil class "+_tag+" was autogenerated with value "+to_string(_Soil.evap_res_fc);
    if (chatty){WriteAdvisory(warn,false);}
  }
  //Shuttleworth b parameter
  autocalc=SetCalculableValue(_Soil.shuttleworth_b,Stmp.shuttleworth_b,Sdefault.shuttleworth_b);
  if (autocalc)
  {
    _Soil.shuttleworth_b=1.0;//?
    //RELATIONSHIP REQUIRED
    //warn="The required parameter SHUTTLEWORTH_B for soil class "+_tag+" was autogenerated with value "+to_string(_Soil.shuttleworth_b);
    //if (chatty){WriteAdvisory(warn,false);}
  }
  //PET Correction factor
  autocalc=SetCalculableValue(_Soil.PET_correction,Stmp.PET_correction,Sdefault.PET_correction);
  if (autocalc)
  {
    _Soil.PET_correction=1.0; //no warning - no correction
  }
  //unavailable fraction
  autocalc=SetCalculableValue(_Soil.unavail_frac,Stmp.unavail_frac,Sdefault.unavail_frac);
  if(autocalc)
  {
    _Soil.unavail_frac=0.0; //no warning
  }
  

  //Albedo
  //----------------------------------------------------------------------------
  autocalc=SetCalculableValue(_Soil.albedo_wet,Stmp.albedo_wet,Sdefault.albedo_wet);
  if (autocalc)
  {
    _Soil.albedo_wet=0.08+0.06*_Soil.sand_con; //WATCLASS CLASSB
    warn="The required parameter ALBEDO_WET for soil class "+_tag+" was autogenerated with value "+to_string(_Soil.albedo_wet);
    if (chatty){WriteAdvisory(warn,false);}
  }
  autocalc=SetCalculableValue(_Soil.albedo_dry,Stmp.albedo_dry,Sdefault.albedo_dry);
  if (autocalc)
  {
    _Soil.albedo_dry=0.14+0.24*_Soil.sand_con; //WATCLASS CLASSB
    warn="The required parameter ALBEDO_DRY for soil class "+_tag+" was autogenerated with value "+to_string(_Soil.albedo_dry);
    if (chatty){WriteAdvisory(warn,false);}
  }

  //Model-specific soil properties - cannot be autocomputed, must be specified by user
  //----------------------------------------------------------------------------
  bool needed=false;
  bool bad;
  bad=SetSpecifiedValue(_Soil.VIC_zmin,Stmp.VIC_zmin,Sdefault.VIC_zmin,needed,"VIC_ZMIN");
  bad=SetSpecifiedValue(_Soil.VIC_zmax,Stmp.VIC_zmax,Sdefault.VIC_zmax,needed,"VIC_ZMAX");
  bad=SetSpecifiedValue(_Soil.VIC_alpha,Stmp.VIC_alpha,Sdefault.VIC_alpha,needed,"VIC_ALPHA");
  bad=SetSpecifiedValue(_Soil.VIC_evap_gamma,Stmp.VIC_evap_gamma,Sdefault.VIC_evap_gamma,needed,"VIC_EVAP_GAMMA");
  bad=SetSpecifiedValue(_Soil.VIC_b_exp,Stmp.VIC_b_exp,Sdefault.VIC_b_exp,needed,"B_EXP");
  bad=SetSpecifiedValue(_Soil.max_perc_rate,Stmp.max_perc_rate,Sdefault.max_perc_rate,needed,"MAX_PERC_RATE");
  bad=SetSpecifiedValue(_Soil.perc_n,Stmp.perc_n,Sdefault.perc_n,needed,"PERC_N");
  bad=SetSpecifiedValue(_Soil.perc_coeff,Stmp.perc_coeff,Sdefault.perc_coeff,needed,"PERC_COEFF");
  bad=SetSpecifiedValue(_Soil.SAC_perc_alpha,Stmp.SAC_perc_alpha,Sdefault.SAC_perc_alpha,needed,"SAC_PERC_ALPHA");
  bad=SetSpecifiedValue(_Soil.SAC_perc_expon,Stmp.SAC_perc_expon,Sdefault.SAC_perc_expon,needed,"SAC_PERC_EXPON");
  bad=SetSpecifiedValue(_Soil.SAC_perc_pfree,Stmp.SAC_perc_pfree,Sdefault.SAC_perc_pfree,needed,"SAC_PERC_PFREE");
  bad=SetSpecifiedValue(_Soil.perc_aspen,Stmp.perc_aspen,Sdefault.perc_aspen,needed,"PERC_ASPEN");
  bad=SetSpecifiedValue(_Soil.max_interflow_rate,Stmp.max_interflow_rate,Sdefault.max_interflow_rate,needed,"MAX_INTERFLOW_RATE");
  bad=SetSpecifiedValue(_Soil.interflow_coeff,Stmp.interflow_coeff,Sdefault.interflow_coeff,needed,"INTERFLOW_COEFF");
  bad=SetSpecifiedValue(_Soil.max_cap_rise_rate,Stmp.max_cap_rise_rate,Sdefault.max_cap_rise_rate,needed,"MAX_CAP_RISE_RATE");
  bad=SetSpecifiedValue(_Soil.max_baseflow_rate,Stmp.max_baseflow_rate,Sdefault.max_baseflow_rate,needed,"MAX_BASEFLOW_RATE");
  bad=SetSpecifiedValue(_Soil.baseflow_n,Stmp.baseflow_n,Sdefault.baseflow_n,needed,"BASEFLOW_N");
  bad=SetSpecifiedValue(_Soil.baseflow_coeff,Stmp.baseflow_coeff,Sdefault.baseflow_coeff,needed,"BASEFLOW_COEFF");
  bad=SetSpecifiedValue(_Soil.HBV_beta,Stmp.HBV_beta,Sdefault.HBV_beta,needed,"HBV_BETA");
  bad=SetSpecifiedValue(_Soil.UBC_evap_soil_def,Stmp.UBC_evap_soil_def,Sdefault.UBC_evap_soil_def,needed,"UBC_EVAP_SOIL_DEF");
  bad=SetSpecifiedValue(_Soil.UBC_infil_soil_def,Stmp.UBC_infil_soil_def,Sdefault.UBC_infil_soil_def,needed,"UBC_INFIL_SOIL_DEF");
  bad=SetSpecifiedValue(_Soil.GR4J_x2,Stmp.GR4J_x2,Sdefault.GR4J_x2,needed,"GR4J_X2");
  bad=SetSpecifiedValue(_Soil.GR4J_x3,Stmp.GR4J_x3,Sdefault.GR4J_x3,needed,"GR4J_X3");
  bad=SetSpecifiedValue(_Soil.baseflow_thresh,Stmp.baseflow_thresh,Sdefault.baseflow_thresh,needed,"BASEFLOW_THRESH");
  bad=SetSpecifiedValue(_Soil.exchange_flow,Stmp.exchange_flow,Sdefault.exchange_flow,needed,"EXCHANGE_FLOW");
  bad=SetSpecifiedValue(_Soil.baseflow_coeff2,Stmp.baseflow_coeff2,Sdefault.baseflow_coeff2,needed,"BASEFLOW_COEFF2");
  bad=SetSpecifiedValue(_Soil.storage_threshold,Stmp.storage_threshold,Sdefault.storage_threshold,needed,"STORAGE_THRESHOLD");

  if(bad) {
    //do something...
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets default soil properties
/// \remark Assumes pure sand with zero organic content
/// \param &S [out] Soil properties class
/// \param is_template [in] True if the default value being set is for the template class
//
void CSoilClass::InitializeSoilProperties(soil_struct &S, bool is_template,int nConstituents)//static
{
  S.org_con =0.0;
  S.clay_con=0.0;
  S.sand_con=1.0;

  //-----------------------------------------------------
  S.porosity                               =DefaultParameterValue(is_template,true);;//AUTO_COMPUTE; //always needed
  S.bulk_density                           =DefaultParameterValue(is_template,true);
  S.stone_frac                             =AUTO_COMPUTE; //always needed
  S.heat_capacity                          =DefaultParameterValue(is_template,true);
  S.thermal_cond                           =DefaultParameterValue(is_template,true);

  S.hydraul_cond                           =DefaultParameterValue(is_template,true);
  S.clapp_b =S.clapp_m =S.clapp_n          =DefaultParameterValue(is_template,true);
  S.sat_res =S.sat_wilt=S.field_capacity   =DefaultParameterValue(is_template,true);

  S.air_entry_pressure =S.wilting_pressure =DefaultParameterValue(is_template,true);
  S.wetting_front_psi                      =DefaultParameterValue(is_template,true);
  S.ksat_std_deviation                     =DefaultParameterValue(is_template,true);
  S.unavail_frac                           =DefaultParameterValue(is_template,true);

  S.evap_res_fc =S.shuttleworth_b          =DefaultParameterValue(is_template,true);
  S.PET_correction                         =DefaultParameterValue(is_template,true);

  S.albedo_wet                             =DefaultParameterValue(is_template,true);
  S.albedo_dry                             =DefaultParameterValue(is_template,true);

  //Conceptual parameters
  //-----------------------------------------------------
  S.VIC_zmin          =DefaultParameterValue(is_template,false);
  S.VIC_zmax          =DefaultParameterValue(is_template,false);
  S.VIC_alpha         =DefaultParameterValue(is_template,false);
  S.VIC_evap_gamma    =DefaultParameterValue(is_template,false);//
  S.VIC_b_exp         =DefaultParameterValue(is_template,false);//
  S.max_perc_rate     =DefaultParameterValue(is_template,false);//
  S.perc_n            =DefaultParameterValue(is_template,false);//
  S.perc_coeff        =DefaultParameterValue(is_template,false);//
  S.SAC_perc_alpha    =DefaultParameterValue(is_template,false);//100;
  S.SAC_perc_expon    =DefaultParameterValue(is_template,false);//3.0;
  S.SAC_perc_pfree    =DefaultParameterValue(is_template,false);//0.06;
  S.perc_aspen        =DefaultParameterValue(is_template,false);//
  S.max_interflow_rate=DefaultParameterValue(is_template,false);//500;    //[mm/d]
  S.interflow_coeff   =DefaultParameterValue(is_template,false);//0.1;
  S.max_baseflow_rate =DefaultParameterValue(is_template,false);//5000;   //[mm/d]
  S.baseflow_n        =DefaultParameterValue(is_template,false);//5;      //[-]
  S.baseflow_coeff    =DefaultParameterValue(is_template,false);//0.100;  //[1/day]
  S.max_cap_rise_rate =DefaultParameterValue(is_template,false);//0.0;    //[mm/d]
  S.HBV_beta          =DefaultParameterValue(is_template,false);//1.0;
  S.UBC_evap_soil_def =DefaultParameterValue(is_template,false);//100.0;  //[mm]
  S.UBC_infil_soil_def=DefaultParameterValue(is_template,false);//100.0;//[mm]
  S.GR4J_x2           =DefaultParameterValue(is_template,false);//0.0 //[mm/d]
  S.GR4J_x3           =DefaultParameterValue(is_template,false);//90 //[mm]
  S.baseflow_thresh   =DefaultParameterValue(is_template,false);//0 //[-]
  S.exchange_flow     =DefaultParameterValue(is_template,false);//0 //[mm/d]
  S.storage_threshold =DefaultParameterValue(is_template,false);//0 //[-]
  S.baseflow_coeff2   =DefaultParameterValue(is_template,false);//0.100;  //[1/day]
}
//////////////////////////////////////////////////////////////////
/// \brief Sets the value of the soil property corresponding to param_name
/// \param param_name [in] Parameter identifier
/// \param value [in] Value of parameter to be set
//
void  CSoilClass::SetSoilProperty(string        param_name,
                                  const double &value)
{
  SetSoilProperty(_Soil,param_name,value);
}
//////////////////////////////////////////////////////////////////
/// \brief Sets the value of the soil property corresponding to param_name
/// \note This is declared as a static member because soil class
/// is not instantiated prior to read of .rvp file
/// \param &S [out] Soil properties class
/// \param param_name [in] Parameter identifier
/// \param value [in] Value of parameter to be set
//
void  CSoilClass::SetSoilProperty(soil_struct &S,
                                  string       param_name,
                                  const double value)
{
  string name;
  name = StringToUppercase(param_name);

  if      (!name.compare("ORG_CON"             )){S.org_con=value;}
  else if (!name.compare("CLAY_CON"            )){S.clay_con=value;}
  else if (!name.compare("SAND_CON"            )){S.sand_con=value;}
  else if (!name.compare("POROSITY"            )){S.porosity=value;}
  else if (!name.compare("STONE_FRAC"          )){S.stone_frac=value;}
  else if (!name.compare("BULK_DENSITY"        )){S.bulk_density=value;}
  else if (!name.compare("HEAT_CAPACITY"       )){S.heat_capacity=value;}
  else if (!name.compare("THERMAL_COND"        )){S.thermal_cond=value;}
  else if (!name.compare("HYDRAUL_COND"        )){S.hydraul_cond=value;}
  else if (!name.compare("CLAPP_B"             )){S.clapp_b=value;}
  else if (!name.compare("CLAPP_M"             )){S.clapp_m=value;}
  else if (!name.compare("CLAPP_N"             )){S.clapp_n=value;}
  else if (!name.compare("SAT_RES"             )){S.sat_res=value;}
  else if (!name.compare("SAT_WILT"            )){S.sat_wilt=value;}
  else if (!name.compare("FIELD_CAPACITY"      )){S.field_capacity=value;}
  else if (!name.compare("AIR_ENTRY_PRESSURE"  )){S.air_entry_pressure=value;}
  else if (!name.compare("WILTING_PRESSURE"    )){S.wilting_pressure=value;}
  else if (!name.compare("WETTING_FRONT_PSI"   )){S.wetting_front_psi=value;}
  else if (!name.compare("KSAT_STD_DEVIATION"  )){S.ksat_std_deviation=value;}
  else if (!name.compare("UNAVAIL_FRAC"        )){S.unavail_frac=value;}

  else if (!name.compare("EVAP_RES_FC"         )){S.evap_res_fc=value;}
  else if (!name.compare("SHUTTLEWORTH_B"      )){S.shuttleworth_b=value;}
  else if (!name.compare("PET_CORRECTION"      )){S.PET_correction=value;}
  else if (!name.compare("ALBEDO_WET"          )){S.albedo_wet=value;}
  else if (!name.compare("ALBEDO_DRY"          )){S.albedo_dry=value;}
  else if (!name.compare("VIC_ZMIN"            )){S.VIC_zmin=value;}
  else if (!name.compare("VIC_ZMAX"            )){S.VIC_zmax=value;}
  else if (!name.compare("VIC_ALPHA"           )){S.VIC_alpha=value;}
  else if (!name.compare("VIC_EVAP_GAMMA"      )){S.VIC_evap_gamma=value;}
  else if (!name.compare("B_EXP"               )){S.VIC_b_exp=value;} //backward compat
  else if (!name.compare("VIC_B_EXP"           )){S.VIC_b_exp=value;}
  else if (!name.compare("MAX_PERC_RATE"       )){S.max_perc_rate=value;}
  else if (!name.compare("PERC_N"              )){S.perc_n=value;}
  else if (!name.compare("PERC_COEFF"          )){S.perc_coeff=value;}
  else if (!name.compare("SAC_PERC_ALPHA"      )){S.SAC_perc_alpha=value;}
  else if (!name.compare("SAC_PERC_EXPON"      )){S.SAC_perc_expon=value;}
  else if (!name.compare("SAC_PERC_PFREE"      )){S.SAC_perc_pfree=value;}
  else if (!name.compare("PERC_ASPEN"          )){S.perc_aspen=value;}
  else if (!name.compare("MAX_INTERFLOW_RATE"  )){S.max_interflow_rate=value;}
  else if (!name.compare("INTERFLOW_COEFF"     )){S.interflow_coeff=value;}
  else if (!name.compare("MAX_BASEFLOW_RATE"   )){S.max_baseflow_rate=value;}
  else if (!name.compare("BASEFLOW_N"          )){S.baseflow_n=value;}
  else if (!name.compare("BASE_STOR_COEFF"     )){S.baseflow_coeff=value;}
  else if (!name.compare("BASEFLOW_COEFF"      )){S.baseflow_coeff=value;}
  else if (!name.compare("MAX_CAP_RISE_RATE"   )){S.max_cap_rise_rate=value;}
  else if (!name.compare("HBV_BETA"            )){S.HBV_beta=value;}
  else if (!name.compare("UBC_EVAP_SOIL_DEF"   )){S.UBC_evap_soil_def =value;}
  else if (!name.compare("UBC_INFIL_SOIL_DEF"  )){S.UBC_infil_soil_def=value;}
  else if (!name.compare("GR4J_X2"             )){S.GR4J_x2 =value;}
  else if (!name.compare("GR4J_X3"             )){S.GR4J_x3 =value;}
  else if (!name.compare("BASEFLOW_THRESH"     )){S.baseflow_thresh =value;}
  else if (!name.compare("EXCHANGE_FLOW"       )){S.exchange_flow=value;}
  else if (!name.compare("BASEFLOW_COEFF2"     )){S.baseflow_coeff2=value; }
  else if (!name.compare("STORAGE_THRESHOLD"   )){S.storage_threshold=value;}
  else{
    WriteWarning("CSoilClass::SetSoilProperty: Unrecognized/invalid soil parameter name ("+name+") in .rvp file",false);
  }
}
///////////////////////////////////////////////////////////////////////////
/// \brief Returns soil property value corresponding to param_name
/// \param param_name [in] Parameter name
/// \return Soil propertie corresponding to parameter name
//
double CSoilClass::GetSoilProperty(string &param_name) const
{
  return GetSoilProperty(_Soil,param_name);
}
///////////////////////////////////////////////////////////////////////////
/// \brief Returns soil property value corresponding to param_name from structure provided
/// \param soil_struct [in] Soil structure
/// \param param_name [in] Parameter name
/// \return Soil propertie corresponding to parameter name
//
double CSoilClass::GetSoilProperty(const soil_struct &S, string param_name,  const bool strict)
{
  string name;
  name = StringToUppercase(param_name);

  if      (!name.compare("ORG_CON"             )){return S.org_con;}
  else if (!name.compare("CLAY_CON"            )){return S.clay_con;}
  else if (!name.compare("SAND_CON"            )){return S.sand_con;}
  else if (!name.compare("POROSITY"            )){return S.porosity;}
  else if (!name.compare("STONE_FRAC"          )){return S.stone_frac;}
  else if (!name.compare("BULK_DENSITY"        )){return S.bulk_density;}
  else if (!name.compare("HEAT_CAPACITY"       )){return S.heat_capacity;}
  else if (!name.compare("THERMAL_COND"        )){return S.thermal_cond;}
  else if (!name.compare("HYDRAUL_COND"        )){return S.hydraul_cond;}
  else if (!name.compare("CLAPP_B"             )){return S.clapp_b;}
  else if (!name.compare("CLAPP_M"             )){return S.clapp_m;}
  else if (!name.compare("CLAPP_N"             )){return S.clapp_n;}
  else if (!name.compare("SAT_RES"             )){return S.sat_res;}
  else if (!name.compare("SAT_WILT"            )){return S.sat_wilt;}
  else if (!name.compare("FIELD_CAPACITY"      )){return S.field_capacity;}
  else if (!name.compare("AIR_ENTRY_PRESSURE"  )){return S.air_entry_pressure;}
  else if (!name.compare("WILTING_PRESSURE"    )){return S.wilting_pressure;}
  else if (!name.compare("WETTING_FRONT_PSI"   )){return S.wetting_front_psi;}
  else if (!name.compare("KSAT_STD_DEVIATION"  )){return S.ksat_std_deviation;}
  else if (!name.compare("UNAVAIL_FRAC"        )){return S.unavail_frac;}

  else if (!name.compare("EVAP_RES_FC"         )){return S.evap_res_fc;}
  else if (!name.compare("SHUTTLEWORTH_B"      )){return S.shuttleworth_b;}
  else if (!name.compare("PET_CORRECTION"      )){return S.PET_correction;}
  else if (!name.compare("ALBEDO_WET"          )){return S.albedo_wet;}
  else if (!name.compare("ALBEDO_DRY"          )){return S.albedo_dry;}
  else if (!name.compare("VIC_ZMIN"            )){return S.VIC_zmin;}
  else if (!name.compare("VIC_ZMAX"            )){return S.VIC_zmax;}
  else if (!name.compare("VIC_ALPHA"           )){return S.VIC_alpha;}
  else if (!name.compare("VIC_EVAP_GAMMA"      )){return S.VIC_evap_gamma;}
  else if (!name.compare("B_EXP"               )){return S.VIC_b_exp;}//backwards compat
  else if (!name.compare("VIC_B_EXP"           )){return S.VIC_b_exp;}
  else if (!name.compare("MAX_PERC_RATE"       )){return S.max_perc_rate;}
  else if (!name.compare("PERC_N"              )){return S.perc_n;}
  else if (!name.compare("PERC_COEFF"          )){return S.perc_coeff;}
  else if (!name.compare("SAC_PERC_ALPHA"      )){return S.SAC_perc_alpha;}
  else if (!name.compare("SAC_PERC_EXPON"      )){return S.SAC_perc_expon;}
  else if (!name.compare("SAC_PERC_PFREE"      )){return S.SAC_perc_pfree;}
  else if (!name.compare("PERC_ASPEN"          )){return S.perc_aspen;}
  else if (!name.compare("MAX_INTERFLOW_RATE"  )){return S.max_interflow_rate;}
  else if (!name.compare("INTERFLOW_COEFF"     )){return S.interflow_coeff;}
  else if (!name.compare("MAX_BASEFLOW_RATE"   )){return S.max_baseflow_rate;}
  else if (!name.compare("BASEFLOW_N"          )){return S.baseflow_n;}
  else if (!name.compare("BASE_STOR_COEFF"     )){return S.baseflow_coeff;}
  else if (!name.compare("BASEFLOW_COEFF"      )){return S.baseflow_coeff;}
  else if (!name.compare("MAX_CAP_RISE_RATE"   )){return S.max_cap_rise_rate;}
  else if (!name.compare("HBV_BETA"            )){return S.HBV_beta;}
  else if (!name.compare("UBC_EVAP_SOIL_DEF"   )){return S.UBC_evap_soil_def;}
  else if (!name.compare("UBC_INFIL_SOIL_DEF"  )){return S.UBC_infil_soil_def;}
  else if (!name.compare("GR4J_X2"             )){return S.GR4J_x2;}
  else if (!name.compare("GR4J_X3"             )){return S.GR4J_x3;}
  else if (!name.compare("BASEFLOW_THRESH"     )){return S.baseflow_thresh;}
  else if (!name.compare("EXCHANGE_FLOW"       )){return S.exchange_flow;}
  else if (!name.compare("BASEFLOW_COEFF2"     )){return S.baseflow_coeff2; }
  else if (!name.compare("STORAGE_THRESHOLD"   )){return S.storage_threshold;}

  //transport parameters NOT HANDLED HERE!
  else{
    if (strict){
      string msg="CSoilClass::GetSoilProperty: Unrecognized/invalid soil parameter name ("+name+") in .rvp file: "+name;
      ExitGracefully(msg.c_str(),BAD_DATA_WARN);
      return 0.0;
    }
    else {
      return INDEX_NOT_FOUND;
    }

  }
}
