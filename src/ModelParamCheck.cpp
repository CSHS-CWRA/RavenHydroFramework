/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2019 the Raven Development Team
----------------------------------------------------------------*/
#include "Model.h"
#include "Properties.h"

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list for all forcing estimation/correction algorithms
///
/// \param *aP [out] array of parameter names needed by all forcing algorithms
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required (size of aP[] and aPC[])
//
void CModel::GetParticipatingParamList(string *aP,class_type *aPC,int &nP,const optStruct &Options) const
{

  /// \todo [clean] - tear out this ugly routine and replace with file read of this information
  nP=0;

  //Just assume needed:
  aP[nP]="FOREST_COVERAGE";   aPC[nP]=CLASS_LANDUSE; nP++;
  aP[nP]="POROSITY";          aPC[nP]=CLASS_SOIL;    nP++;

  // Interpolation Method parameters
  //----------------------------------------------------------------------
  if((Options.interpolation==INTERP_NEAREST_NEIGHBOR) || (Options.interpolation==INTERP_AVERAGE_ALL))
  {
    // no parameter required
  }
  else if(Options.interpolation==INTERP_FROM_FILE)
  {
    // timeseries at gauge
  }
  else if(Options.interpolation==INTERP_INVERSE_DISTANCE)
  {
    // This method have not been tested yet
  }

  // Routing algorithms parameters
  //----------------------------------------------------------------------
  // Routing Method
  if((Options.routing==ROUTE_NONE) || (Options.routing==ROUTE_DIFFUSIVE_WAVE) || (Options.routing==ROUTE_PLUG_FLOW)
    || (Options.routing==ROUTE_STORAGECOEFF) || (Options.routing==ROUTE_MUSKINGUM) || (Options.routing==ROUTE_MUSKINGUM_CUNGE))
  {
    // Parameters are located in the RVH file?
    // Channel Geometry
    // Manning's n
  }

  // Catchment Route Method
  //----------------------------------------------------------------------
  if(Options.catchment_routing==ROUTE_DELAYED_FIRST_ORDER)
  {
    // This method have not been tested yet
    // Parameters are located in the RVH file, and there's no checking routine for that file yet.
    //aP[nP]="TIME_LAG"; aPC[nP]=CLASS_; nP++;
    //aP[nP]="RES_CONSTANT"; aPC[nP]=CLASS_; nP++;
  }
  //else if (Options.catchment_routing==ROUTE_LAG)
  //{
  //ROUTE_LAG method not implemented yet (See ParseInput.cpp)
  //aP[nP]="TIME_LAG"; aPC[nP]=CLASS_L; nP++;
  //}
  else if(Options.catchment_routing==ROUTE_GAMMA_CONVOLUTION)
  {
    // Parameters are located in the RVH file, and there's no checking routine for that file yet.
    //aP[nP]="TIME_CONC"; aPC[nP]=CLASS_SUBBASIN; nP++;
    //aP[nP]="GAMMA_SHAPE"; aPC[nP]=CLASS_SUBBASIN; nP++;
  }
  else if(Options.catchment_routing==ROUTE_TRI_CONVOLUTION)
  {
    // Parameters are located in the RVH file, and there's no checking routine for that file yet.
    //aP[nP]="TIME_TO_PEAK"; aPC[nP]=CLASS_; nP++;
    //aP[nP]="TIME_CONC"; aPC[nP]=CLASS_; nP++;
  }
  else if(Options.catchment_routing==ROUTE_RESERVOIR_SERIES)
  {
    // Parameters are located in the RVH file, and there's no checking routine for that file yet.
    //aP[nP]="RES_CONSTANT"; aPC[nP]=CLASS_; nP++;
    //aP[nP]="NUM_RESERVOIRS"; aPC[nP]=CLASS_; nP++;
  }
  else if(Options.catchment_routing==ROUTE_DUMP)
  {
    // no parameter required
  }

  // Evaporation algorithms parameters
  //----------------------------------------------------------------------

  // Evaporation Method
  if(Options.evaporation==PET_DATA)
  {
    // timeseries at gauge
  }
  else if(Options.evaporation==PET_FROMMONTHLY)
  {
    // Parameters are located in the RVT file, and there's no checking routine for that file yet.
    //aP[nP]=":MonthlyAveEvaporation"; aPC[nP]=CLASS_; nP++;
    //aP[nP]=":MonthlyAveTemperature"; aPC[nP]=CLASS_; nP++;
  }
  else if(Options.evaporation==PET_MONTHLY_FACTOR)
  {
    aP[nP]="FOREST_PET_CORR"; aPC[nP]=CLASS_LANDUSE; nP++;
    // Parameters are located in the RVT file, and there's no checking routine for that file yet.
    //aP[nP]=":MonthlyEvapFactor";   aPC[nP]=CLASS_; nP++;
  }
  else if(Options.evaporation==PET_PENMAN_MONTEITH)
  {
    aP[nP]="MAX_HEIGHT";    aPC[nP]=CLASS_VEGETATION; nP++;
    aP[nP]="RELATIVE_HT";   aPC[nP]=CLASS_VEGETATION; nP++;
    aP[nP]="MAX_LAI";       aPC[nP]=CLASS_VEGETATION; nP++; //JRCFLAG
    aP[nP]="RELATIVE_LAI";  aPC[nP]=CLASS_VEGETATION; nP++;
    aP[nP]="MAX_LEAF_COND"; aPC[nP]=CLASS_VEGETATION; nP++;
    aP[nP]="FOREST_SPARSENESS";    aPC[nP]=CLASS_LANDUSE; nP++;
    aP[nP]="ROUGHNESS";     aPC[nP]=CLASS_LANDUSE; nP++;
  }
  else if(Options.evaporation==PET_PENMAN_COMBINATION)
  {
    aP[nP]="MAX_HEIGHT";  aPC[nP]=CLASS_VEGETATION; nP++;
    aP[nP]="RELATIVE_HT"; aPC[nP]=CLASS_VEGETATION; nP++;
  }
  else if(Options.evaporation==PET_HARGREAVES)
  {
    // Need Max and Min Monthly Temp
    //aP[nP]="TEMP_MONTH_MAX";  aPC[nP]=CLASS_; nP++;
    //aP[nP]="TEMP_MONTH_MIN";  aPC[nP]=CLASS_; nP++;
  }
  else if(Options.evaporation==PET_MOHYSE)
  {
    aP[nP]="MOHYSE_PET_COEFF"; aPC[nP]=CLASS_GLOBAL; nP++;
  }
  else if((Options.evaporation==PET_CONSTANT) || (Options.evaporation==PET_HAMON) || (Options.evaporation==PET_HARGREAVES_1985)
    || (Options.evaporation==PET_TURC_1961) || (Options.evaporation==PET_MAKKINK_1957) || (Options.evaporation==PET_PRIESTLEY_TAYLOR) || (Options.evaporation==PET_SHUTTLEWORTH_WALLACE))
  {
    // no parameter required
  }

  //Anywhere Albedo needs to be calculated for SW_Radia_net (will later be moved to albedo options)
  if((Options.evaporation==PET_PRIESTLEY_TAYLOR) || (Options.evaporation==PET_SHUTTLEWORTH_WALLACE) ||
    (Options.evaporation==PET_PENMAN_MONTEITH)  || (Options.evaporation==PET_PENMAN_COMBINATION) ||
    (Options.evaporation==PET_JENSEN_HAISE) || (Options.evaporation==PET_GRANGERGRAY))
  {
    if(Options.SW_radia_net == NETSWRAD_CALC) {
      aP[nP]="ALBEDO";            aPC[nP]=CLASS_VEGETATION; nP++; //for SW_Radia_net
      aP[nP]="ALBEDO_WET";        aPC[nP]=CLASS_SOIL; nP++; //for SW_Radia_net
      aP[nP]="ALBEDO_DRY";        aPC[nP]=CLASS_SOIL; nP++; //for SW_Radia_net
      aP[nP]="SVF_EXTINCTION";    aPC[nP]=CLASS_VEGETATION; nP++; //for SW_Radia_net
    }
  }

  // OW_Evaporation Method
  //----------------------------------------------------------------------
  if(Options.ow_evaporation==PET_DATA)
  {
    // timeseries at gauge
  }
  else if(Options.ow_evaporation==PET_FROMMONTHLY)
  {
    // \todo [funct] Parameters are located in the RVT file, and there's no checking routine for that file yet.
    //aP[nP]=":MonthlyAveEvaporation"; aPC[nP]=CLASS_; nP++;
    //aP[nP]=":MonthlyAveTemperature"; aPC[nP]=CLASS_; nP++;
  }
  else if(Options.ow_evaporation==PET_MONTHLY_FACTOR)
  {
    aP[nP]="FOREST_PET_CORR"; aPC[nP]=CLASS_LANDUSE; nP++;
    // \todo [funct] Parameters are located in the RVT file, and there's no checking routine for that file yet.
    //aP[nP]=":MonthlyEvapFactor";   aPC[nP]=CLASS_; nP++;
  }
  else if(Options.ow_evaporation==PET_PENMAN_MONTEITH)
  {
    aP[nP]="MAX_HEIGHT";    aPC[nP]=CLASS_VEGETATION; nP++;
    aP[nP]="RELATIVE_HT";   aPC[nP]=CLASS_VEGETATION; nP++;
    aP[nP]="MAX_LAI";       aPC[nP]=CLASS_VEGETATION; nP++; //JRCFLAG
    aP[nP]="RELATIVE_LAI";  aPC[nP]=CLASS_VEGETATION; nP++;
    aP[nP]="MAX_LEAF_COND"; aPC[nP]=CLASS_VEGETATION; nP++;
    aP[nP]="FOREST_SPARSENESS";    aPC[nP]=CLASS_LANDUSE; nP++;
    aP[nP]="ROUGHNESS";     aPC[nP]=CLASS_LANDUSE; nP++;
  }
  else if(Options.ow_evaporation==PET_PENMAN_COMBINATION)
  {
    aP[nP]="MAX_HEIGHT";  aPC[nP]=CLASS_VEGETATION; nP++;
    aP[nP]="RELATIVE_HT"; aPC[nP]=CLASS_VEGETATION; nP++;
  }
  else if((Options.ow_evaporation==PET_CONSTANT) || (Options.ow_evaporation==PET_HAMON) ||
    (Options.ow_evaporation==PET_HARGREAVES) || (Options.ow_evaporation==PET_HARGREAVES_1985) ||
    (Options.ow_evaporation==PET_TURC_1961) || (Options.ow_evaporation==PET_MAKKINK_1957) ||
    (Options.ow_evaporation==PET_PRIESTLEY_TAYLOR) || (Options.ow_evaporation==PET_GRANGERGRAY))
  {
    // no parameter required/listed
  }
  // Orographic PET Correction Method
  //------------------------------------------------------------
  if(Options.orocorr_PET==OROCORR_UBCWM)
  {
    // hardcoded for now
  }
  else if((Options.orocorr_PET==OROCORR_NONE) || (Options.orocorr_PET==OROCORR_SIMPLELAPSE))
  {
    // no parameter required
  }

  // Radiation algorithms parameters
  //----------------------------------------------------------------------
  // SW Radiation Method
  if(Options.SW_radiation==SW_RAD_UBCWM)
  {
    // currently hardcoded
    // Parameters are located in the RVH file, and there's no checking routine for that file yet.
    //aP[nP]="HORIZON_CORR"; aPC[nP]=CLASS_; nP++;
    //aP[nP]="RES_CONSTANT"; aPC[nP]=CLASS_; nP++;
    aP[nP]="UBC_EXPOSURE_FACT"; aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="UBC_SW_S_CORR";     aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="UBC_SW_N_CORR";     aPC[nP]=CLASS_GLOBAL; nP++;
  }
  else if(Options.SW_radiation==SW_RAD_DEFAULT)
  {
    // SLOPE and ASPECT are requried, but a check is unnecessary
  }
  else if(Options.SW_radiation==SW_RAD_DATA)
  {
    // timeseries at gauge
  }
  else if(Options.SW_radiation==SW_RAD_NONE)
  {
  }

  // LW Radiation Method
  //----------------------------------------------------------------------
  if(Options.LW_radiation==LW_RAD_DEFAULT)
  {
    // FOREST_COVERAGE is required, but a check is unnecessary
  }
  else if(Options.LW_radiation==LW_RAD_UBCWM)
  {

    aP[nP]="UBC_LW_FOREST_FACT"; aPC[nP]=CLASS_GLOBAL; nP++;
  }
  else if(Options.LW_radiation==LW_RAD_DATA)
  {
    // timeseries at gauge
  }

  // Cloud Cover Method
  //----------------------------------------------------------------------
  if(Options.cloud_cover==CLOUDCOV_DATA)
  {
    // timeseries at gauge
  }
  else if(Options.cloud_cover==CLOUDCOV_NONE)
  {
    //no parameters needed
  }
  else if(Options.cloud_cover==CLOUDCOV_UBCWM)
  {
    //no parameters needed except Gauge->_cloud_min_temp and Gauge->_cloud_max_temp, which defaults to no cloud cover
  }
  //else if (Options.cloud_cover==CLOUDCOV_HARGREAVES)
  //{
  // Not Implemented yet
  //}

  //Canopy cover SW correction method
  //----------------------------------------------------------------------
  if(Options.SW_canopycorr==SW_CANOPY_CORR_STATIC) {
    aP[nP]="RELATIVE_LAI";     aPC[nP]=CLASS_VEGETATION; nP++;
    aP[nP]="FOREST_SPARSENESS";aPC[nP]=CLASS_LANDUSE; nP++;
    aP[nP]="SVF_EXTINCTION";  aPC[nP]=CLASS_VEGETATION; nP++;
  }
  else if(Options.SW_canopycorr ==SW_CANOPY_CORR_UBCWM) {
    aP[nP]="UBC_EXPOSURE_FACT";  aPC[nP]=CLASS_GLOBAL; nP++;
  }

  //Cloud cover SW correction method
  //----------------------------------------------------------------------
  if(Options.SW_cloudcovercorr ==SW_CLOUD_CORR_UBCWM) {
    aP[nP]="UBC_CLOUD_PENET";  aPC[nP]=CLASS_GLOBAL; nP++;
  }

  // ET Radiation Method
  //----------------------------------------------------------------------
  // Not Implemented yet

  // Precipitation algorithms parameters
  //----------------------------------------------------------------------
  // Rain Snow Fraction Method
  if(Options.rainsnow==RAINSNOW_DATA)
  {
    // timeseries at gauge
  }
  else if(Options.rainsnow==RAINSNOW_DINGMAN)
  {
    aP[nP]="RAINSNOW_TEMP"; aPC[nP]=CLASS_GLOBAL; nP++;
  }
  else if((Options.rainsnow==RAINSNOW_HBV) || (Options.rainsnow==RAINSNOW_UBCWM))
  {
    aP[nP]="RAINSNOW_TEMP";  aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="RAINSNOW_DELTA"; aPC[nP]=CLASS_GLOBAL; nP++;
  }
  else if(Options.rainsnow==RAINSNOW_HARDER) {

  }


  // Precipitation Interception Fraction Method
  //----------------------------------------------------------------------
  //handled by CmvPrecipitation::GetParticipatingParamList

  // Orographic Precipitation Correction Method
  //----------------------------------------------------------------------
  if(Options.orocorr_precip==OROCORR_HBV)
  {
    // Parameters are located in the RVT file, and there's no checking routine for that file yet.
    //aP[nP]=":RainCorrection"; aPC[nP]=CLASS_; nP++;
    //aP[nP]=":SnowCorrection"; aPC[nP]=CLASS_; nP++;
  }
  else if((Options.orocorr_precip==OROCORR_UBCWM) || (Options.orocorr_precip==OROCORR_UBCWM2))
  {
    aP[nP]="UBC_E0LHI"; aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="UBC_E0LLOW"; aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="UBC_E0LMID"; aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="UBC_P0GRADL"; aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="UBC_P0GRADM"; aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="UBC_P0GRADU"; aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="UBC_A0STAB"; aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="UBC_MAX_RANGE_TEMP"; aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="UBC_A0PPTP";       aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="MAX_INTERCEPT_RATE";aPC[nP]=CLASS_VEGETATION; nP++;
    aP[nP]="RAIN_ICEPT_PCT";    aPC[nP]=CLASS_VEGETATION; nP++;
    aP[nP]="FOREST_COVERAGE";   aPC[nP]=CLASS_LANDUSE; nP++; //JRCFLAG
    aP[nP]="UBC_ICEPT_FACTOR";   aPC[nP]=CLASS_LANDUSE; nP++;
  }
  else if(Options.orocorr_precip==OROCORR_SIMPLELAPSE)
  {
    aP[nP]="PRECIP_LAPSE"; aPC[nP]=CLASS_GLOBAL; nP++;
  }
  else if(Options.orocorr_precip==OROCORR_NONE)
  {
    // no parameter required
  }

  // Temperature algorithms parameters
  //----------------------------------------------------------------------
  // Orographic Temperature Correction Method
  if((Options.orocorr_temp==OROCORR_UBCWM)  || (Options.orocorr_temp==OROCORR_UBCWM2))
  {
    aP[nP]="UBC_A0TLXH"; aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="UBC_A0TLNM"; aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="UBC_A0TLXM"; aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="UBC_A0TLNH"; aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="UBC_P0TEDL"; aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="UBC_P0TEDU"; aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="UBC_MAX_RANGE_TEMP"; aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="ADIABATIC_LAPSE";     aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="WET_ADIABATIC_LAPSE"; aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="MAX_INTERCEPT_RATE"; aPC[nP]=CLASS_VEGETATION; nP++;
    aP[nP]="RAIN_ICEPT_PCT";     aPC[nP]=CLASS_VEGETATION; nP++;
  }
  else if((Options.orocorr_temp==OROCORR_HBV) || (Options.orocorr_temp==OROCORR_SIMPLELAPSE))
  {
    aP[nP]="ADIABATIC_LAPSE"; aPC[nP]=CLASS_GLOBAL; nP++;
  }
  else if(Options.orocorr_temp==OROCORR_NONE)
  {
    // no parameter required
  }

  // Energy algorithms parameters
  //----------------------------------------------------------------------
  // Potential Melt Method
  if(Options.pot_melt == POTMELT_DATA)
  {
    //none
  }
  else if(Options.pot_melt == POTMELT_DEGREE_DAY)
  {
    aP[nP]="MELT_FACTOR"; aPC[nP]=CLASS_LANDUSE; nP++;
    aP[nP]="DD_MELT_TEMP";aPC[nP]=CLASS_LANDUSE; nP++;
  }
  else if(Options.pot_melt==POTMELT_RESTRICTED)
  {
    aP[nP]="MELT_FACTOR"; aPC[nP]=CLASS_LANDUSE; nP++;
    aP[nP]="DD_MELT_TEMP";aPC[nP]=CLASS_LANDUSE; nP++;
  }
  else if(Options.pot_melt==POTMELT_HBV)
  {
    aP[nP]="MELT_FACTOR";       aPC[nP]=CLASS_LANDUSE; nP++;
    aP[nP]="DD_MELT_TEMP";      aPC[nP]=CLASS_LANDUSE; nP++;
    aP[nP]="MIN_MELT_FACTOR";   aPC[nP]=CLASS_LANDUSE; nP++;
    aP[nP]="HBV_MELT_ASP_CORR"; aPC[nP]=CLASS_LANDUSE; nP++;
    aP[nP]="HBV_MELT_FOR_CORR"; aPC[nP]=CLASS_LANDUSE; nP++;
  }
  else if(Options.pot_melt==POTMELT_HBV_ROS)
  {
    aP[nP]="MELT_FACTOR";       aPC[nP]=CLASS_LANDUSE; nP++;
    aP[nP]="DD_MELT_TEMP";      aPC[nP]=CLASS_LANDUSE; nP++;
    aP[nP]="MIN_MELT_FACTOR";   aPC[nP]=CLASS_LANDUSE; nP++;
    aP[nP]="HBV_MELT_ASP_CORR"; aPC[nP]=CLASS_LANDUSE; nP++;
    aP[nP]="HBV_MELT_FOR_CORR"; aPC[nP]=CLASS_LANDUSE; nP++;
    aP[nP]="RAIN_MELT_MULT";    aPC[nP]=CLASS_LANDUSE; nP++;
  }
  else if(Options.pot_melt==POTMELT_UBCWM)
  {
    aP[nP]="FOREST_COVERAGE";     aPC[nP]=CLASS_LANDUSE; nP++; //JRCFLAG
    aP[nP]="RAIN_MELT_MULT";      aPC[nP]=CLASS_LANDUSE; nP++;
    aP[nP]="CONV_MELT_MULT";      aPC[nP]=CLASS_LANDUSE; nP++;
    aP[nP]="COND_MELT_MULT";      aPC[nP]=CLASS_LANDUSE; nP++;

    aP[nP]="MIN_SNOW_ALBEDO";     aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="UBC_SW_S_CORR";       aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="UBC_SW_N_CORR";       aPC[nP]=CLASS_GLOBAL; nP++;
  }
  else if(Options.pot_melt==POTMELT_USACE)
  {
    aP[nP]="WIND_EXPOSURE";       aPC[nP]=CLASS_LANDUSE; nP++;
  }
  else if(Options.pot_melt==POTMELT_EB)
  {
    aP[nP]="MAX_HEIGHT";       aPC[nP]=CLASS_VEGETATION; nP++;
    aP[nP]="RELATIVE_HT";      aPC[nP]=CLASS_VEGETATION; nP++;
    aP[nP]="ROUGHNESS";        aPC[nP]=CLASS_LANDUSE; nP++;
    aP[nP]="SNOW_TEMPERATURE"; aPC[nP]=CLASS_GLOBAL; nP++;
  }
  else if(Options.pot_melt==POTMELT_CRHM_EBSM)
  {
  }
  else if(Options.pot_melt==POTMELT_HMETS)
  {
    aP[nP]="MIN_MELT_FACTOR";   aPC[nP]=CLASS_LANDUSE; nP++;
    aP[nP]="MAX_MELT_FACTOR";   aPC[nP]=CLASS_LANDUSE; nP++;
    aP[nP]="DD_MELT_TEMP";      aPC[nP]=CLASS_LANDUSE; nP++;
    aP[nP]="DD_AGGRADATION";    aPC[nP]=CLASS_LANDUSE; nP++;
  }
  // Sub Daily Method
  //----------------------------------------------------------------------
  if((Options.subdaily==SUBDAILY_NONE) || (Options.subdaily==SUBDAILY_SIMPLE) || (Options.subdaily==SUBDAILY_UBC))
  {
    // no parameter required
  }

  // Atmospheric Variables algorithms parameters
  //----------------------------------------------------------------------
  // Wind Speed Method
  if(Options.wind_velocity==WINDVEL_DATA)
  {
    // timeseries at gauge
  }
  else if(Options.wind_velocity==WINDVEL_UBCWM)
  {
    aP[nP]="FOREST_COVERAGE"; aPC[nP]=CLASS_LANDUSE; nP++; 
    aP[nP]="UBC_P0TEDL"; aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="UBC_P0TEDU"; aPC[nP]=CLASS_GLOBAL; nP++;
    aP[nP]="UBC_MAX_RANGE_TEMP"; aPC[nP]=CLASS_GLOBAL; nP++;
  }
  else if(Options.wind_velocity==WINDVEL_CONSTANT)
  {
    // no parameter required
  }

  // Relative Humidity Method parameters
  //----------------------------------------------------------------------
  if((Options.rel_humidity==RELHUM_CONSTANT) || (Options.rel_humidity==RELHUM_DATA) || (Options.rel_humidity==RELHUM_MINDEWPT))
  {
    // no parameter required
  }

  // Air Pressure Method parameters
  //----------------------------------------------------------------------
  if(Options.air_pressure==AIRPRESS_DATA)
  {
    // timeseries at gauge
  }
  else if((Options.air_pressure==AIRPRESS_BASIC) || (Options.air_pressure==AIRPRESS_UBC) || (Options.air_pressure==AIRPRESS_CONST))
  {
    // no parameter required
  }

  // Monthly Interpolation Method parameters
  //----------------------------------------------------------------------
  if((Options.month_interp==MONTHINT_UNIFORM) || (Options.month_interp==MONTHINT_LINEAR_FOM)
    || (Options.month_interp==MONTHINT_LINEAR_MID) || (Options.month_interp==MONTHINT_LINEAR_21))
  {
    // no parameter required
  }


  //...
}

