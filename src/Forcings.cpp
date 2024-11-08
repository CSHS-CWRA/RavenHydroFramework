/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2024 the Raven Development Team
  ----------------------------------------------------------------*/

#include "RavenInclude.h"

/*****************************************************************
   Forcing functions
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Sets the value of all entries in the input force structure to 0.0
///
/// \param &F [out] HRU forcing functions structure
//
void ZeroOutForcings(force_struct &F)
{
  F.precip=0.0;
  F.precip_daily_ave = 0.0;
  F.precip_5day=0.0;
  F.snow_frac=0.0;
  F.irrigation=0.0;

  F.temp_ave=0.0;
  F.temp_daily_min=0.0;
  F.temp_daily_max=0.0;
  F.temp_daily_ave=0.0;
  F.temp_month_max=0.0;
  F.temp_month_min=0.0;
  F.temp_month_ave=0.0;
  F.temp_ave_unc  =0.0;
  F.temp_max_unc  =0.0;
  F.temp_min_unc  =0.0;


  F.air_dens    =0.0;
  F.air_pres    =0.0;
  F.rel_humidity=0.0;

  F.cloud_cover=0.0;
  F.ET_radia=0.0;
  F.SW_radia=0.0;
  F.SW_radia_subcan=0.0;
  F.SW_subcan_net=0.0;
  F.LW_incoming=0.0;
  F.SW_radia_net=0.0;
  F.LW_radia_net=0.0;

  F.day_length=0.0;
  F.day_angle=0.0;

  F.wind_vel=0.0;

  F.PET=0.0;
  F.OW_PET=0.0;
  F.PET_month_ave=0.0;

  F.potential_melt=0.0;

  F.recharge=0.0;
  F.precip_temp=0.0;
  F.precip_conc=0.0;

  F.subdaily_corr=0.0;
}
//////////////////////////////////////////////////////////////////
/// \brief Copys only the forcings that are constant over the course of a day from stcuture Ffrom to structure Fto
///
//
void CopyDailyForcingItems(force_struct& Ffrom, force_struct& Fto)
{
  Fto.temp_daily_min = Ffrom.temp_daily_min;
  Fto.temp_daily_max = Ffrom.temp_daily_max;
  Fto.temp_daily_ave = Ffrom.temp_daily_ave;
  Fto.temp_month_max = Ffrom.temp_month_max;
  Fto.temp_month_min = Ffrom.temp_month_min;
  Fto.temp_month_ave = Ffrom.temp_month_ave;

  Fto.day_length     = Ffrom.day_length;
  Fto.day_angle      = Ffrom.day_angle;

  //Fto.SW_radia_unc   = Ffrom.SW_radia_unc;
  //Fto.ET_radia       = Ffrom.ET_radia;
  //Fto.ET_radia_flat  = Ffrom.ET_radia_flat;

  Fto.PET_month_ave  = Ffrom.PET_month_ave;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns enumerated forcing type from corresponding string
///
/// \param &forcing_string [in] String identifier of a forcing type
/// \returns enumerated forcing type, F_UNRECOGNIZED if string does not correspond to forcing
//
forcing_type GetForcingTypeFromString(const string &forcing_string)
{
  string f=StringToUppercase(forcing_string);
  if      (f=="PRECIP"           ){return F_PRECIP;}
  else if (f=="PRECIP_DAILY_AVE" ){return F_PRECIP_DAILY_AVE; }
  else if (f=="PRECIP_5DAY"      ){return F_PRECIP_5DAY;}
  else if (f=="SNOW_FRAC"        ){return F_SNOW_FRAC;}
  else if (f=="SNOWFALL"         ){return F_SNOWFALL;}
  else if (f=="RAINFALL"         ){return F_RAINFALL;}
  else if (f=="IRRIGATION"       ){return F_IRRIGATION; }
  else if (f=="TEMP_AVE"         ){return F_TEMP_AVE;}
  else if (f=="TEMP_MIN"         ){return F_TEMP_DAILY_MIN;}
  else if (f=="MIN_TEMPERATURE"  ){return F_TEMP_DAILY_MIN;}//alias
  else if (f=="TEMP_DAILY_MIN"   ){return F_TEMP_DAILY_MIN;}//alias
  else if (f=="TEMP_MAX"         ){return F_TEMP_DAILY_MAX;}
  else if (f=="MAX_TEMPERATURE"  ){return F_TEMP_DAILY_MAX;}
  else if (f=="TEMP_DAILY_MAX"   ){return F_TEMP_DAILY_MAX;}//alias
  else if (f=="TEMP_DAILY_AVE"   ){return F_TEMP_DAILY_AVE;}//alias
  else if (f=="TEMP_MONTH_MAX"   ){return F_TEMP_MONTH_MAX;}
  else if (f=="TEMP_MONTH_MIN"   ){return F_TEMP_MONTH_MIN;}
  else if (f=="TEMP_MONTH_AVE"   ){return F_TEMP_MONTH_AVE;}
  else if (f=="TEMP_AVE_UNC"     ){return F_TEMP_AVE_UNC;}
  else if (f=="TEMP_MAX_UNC"     ){return F_TEMP_MAX_UNC;}
  else if (f=="TEMP_MIN_UNC"     ){return F_TEMP_MIN_UNC;}

  else if (f=="AIR_DENS"         ){return F_AIR_DENS;}
  else if (f=="AIR_DENSITY"      ){return F_AIR_DENS;}//alias
  else if (f=="AIR_PRES"         ){return F_AIR_PRES;}
  else if (f=="AIR_PRESSURE"     ){return F_AIR_PRES;}//alias
  else if (f=="REL_HUMIDITY"     ){return F_REL_HUMIDITY;}

  else if (f=="ET_RADIA"         ){return F_ET_RADIA;}
  else if (f=="ET_RADIA_FLAT"    ){return F_ET_RADIA_FLAT;}
  else if (f=="SW_RADIA"         ){return F_SW_RADIA;}
  else if (f=="SHORTWAVE"        ){return F_SW_RADIA;}//alias
  else if (f=="SW_RADIA_NET"     ){return F_SW_RADIA_NET;}
  else if (f=="SW_RADIA_SUBCAN"  ){return F_SW_RADIA_SUBCAN;}
  else if (f=="SW_SUBCAN_NET"    ){return F_SW_SUBCAN_NET;}
  else if (f=="LW_INCOMING"      ){return F_LW_INCOMING;}
  else if (f=="LW_RADIA"         ){return F_LW_RADIA_NET;}//alias
  else if (f=="LONGWAVE"         ){return F_LW_RADIA_NET;}//alias
  else if (f=="LW_RADIA_NET"     ){return F_LW_RADIA_NET;}
  else if (f=="CLOUD_COVER"      ){return F_CLOUD_COVER;}

  else if (f=="DAY_LENGTH"       ){return F_DAY_LENGTH;}
  else if (f=="DAY_ANGLE"        ){return F_DAY_ANGLE;}

  else if (f=="WIND_VEL"         ){return F_WIND_VEL;}

  else if (f=="PET"              ){return F_PET;}
  else if (f=="OW_PET"           ){return F_OW_PET;}
  else if (f=="PET_MONTH_AVE"    ){return F_PET_MONTH_AVE;}

  else if (f=="POTENTIAL_MELT"   ){return F_POTENTIAL_MELT;}

  else if (f=="RECHARGE"         ){return F_RECHARGE;}
  else if (f=="PRECIP_TEMP"      ){return F_PRECIP_TEMP; }
  else if (f=="PRECIP_CONC"      ){return F_PRECIP_CONC; }
  else if (f=="SUBDAILY_CORR"    ){return F_SUBDAILY_CORR;}

  return F_UNRECOGNIZED;
}
/////////////////////////////////////////////////////////////////////
/// \brief Return double value of forcing function specified by passed string parameter of force structure f
///
/// \param &ftype [in] forcing function as enumerateed type
/// \param &f [out] HRU forcing functions structure
/// \return Double forcing function value corresponding to ftype
//
double GetForcingFromType(const forcing_type &ftype, const force_struct &f)
{
  if      (ftype==F_PRECIP          ){return f.precip;}
  else if (ftype==F_PRECIP_DAILY_AVE){return f.precip_daily_ave;}
  else if (ftype==F_PRECIP_5DAY     ){return f.precip_5day;}
  else if (ftype==F_SNOW_FRAC       ){return f.snow_frac;}
  else if (ftype==F_SNOWFALL        ){return (    f.snow_frac)*f.precip;}
  else if (ftype==F_RAINFALL        ){return (1.0-f.snow_frac)*f.precip;}
  else if (ftype==F_IRRIGATION      ){return f.irrigation;}

  else if (ftype==F_TEMP_AVE        ){return f.temp_ave;}
  else if (ftype==F_TEMP_DAILY_MIN  ){return f.temp_daily_min;}
  else if (ftype==F_TEMP_DAILY_MAX  ){return f.temp_daily_max;}
  else if (ftype==F_TEMP_DAILY_AVE  ){return f.temp_daily_ave;}
  else if (ftype==F_TEMP_MONTH_MAX  ){return f.temp_month_max;}
  else if (ftype==F_TEMP_MONTH_MIN  ){return f.temp_month_min;}
  else if (ftype==F_TEMP_MONTH_AVE  ){return f.temp_month_ave;}

  else if (ftype==F_TEMP_AVE_UNC    ){return f.temp_ave_unc ;}
  else if (ftype==F_TEMP_MAX_UNC    ){return f.temp_max_unc ;}
  else if (ftype==F_TEMP_MIN_UNC    ){return f.temp_min_unc ;}


  else if (ftype==F_AIR_DENS        ){return f.air_dens;}
  else if (ftype==F_AIR_PRES        ){return f.air_pres;}
  else if (ftype==F_REL_HUMIDITY    ){return f.rel_humidity;}

  else if (ftype==F_CLOUD_COVER     ){return f.cloud_cover;}
  else if (ftype==F_ET_RADIA        ){return f.ET_radia;}
  else if (ftype==F_ET_RADIA_FLAT   ){return f.ET_radia_flat; }
  else if (ftype==F_SW_RADIA        ){return f.SW_radia;}
  else if (ftype==F_SW_RADIA_UNC    ){return f.SW_radia_unc;}
  else if (ftype==F_SW_RADIA_SUBCAN ){return f.SW_radia_subcan;}
  else if (ftype==F_SW_SUBCAN_NET   ){return f.SW_subcan_net;}
  else if (ftype==F_SW_RADIA_NET    ){return f.SW_radia_net;}
  else if (ftype==F_LW_INCOMING     ){return f.LW_incoming;}
  else if (ftype==F_LW_RADIA_NET    ){return f.LW_radia_net;}

  else if (ftype==F_DAY_LENGTH      ){return f.day_length;}
  else if (ftype==F_DAY_ANGLE       ){return f.day_angle;}

  else if (ftype==F_WIND_VEL        ){return f.wind_vel;}

  else if (ftype==F_PET             ){return f.PET;}
  else if (ftype==F_OW_PET          ){return f.OW_PET;}
  else if (ftype==F_PET_MONTH_AVE   ){return f.PET_month_ave;}

  else if (ftype==F_POTENTIAL_MELT  ){return f.potential_melt;}

  else if (ftype==F_RECHARGE        ){return f.recharge;}
  else if (ftype==F_PRECIP_TEMP     ){return f.precip_temp;}
  else if (ftype==F_PRECIP_CONC     ){return f.precip_conc;}

  else if (ftype==F_SUBDAILY_CORR   ){return f.subdaily_corr;}

  //else if (ftype==F_UNRECOGNIZED   ){return 0;}

  else
  {
    ExitGracefully("GetForcingFromType: invalid forcing type",RUNTIME_ERR);
  }
  return 0.0;
}
/////////////////////////////////////////////////////////////////////
/// \brief Modify value of forcing function specified by forcing type within force structure f
///
/// \param &ftype [in] forcing function as enumerateed type
/// \param &f [out] HRU forcing functions structure
/// \param val - value to set
//
void  SetForcingFromType(const forcing_type &ftype, force_struct &f, const double &val)
{
  if      (ftype==F_PRECIP          ){f.precip=val;}
  else if (ftype==F_PRECIP_DAILY_AVE){f.precip_daily_ave=val;}
  else if (ftype==F_PRECIP_5DAY     ){f.precip_5day=val;}
  else if (ftype==F_SNOW_FRAC       ){f.snow_frac=val;}
  else if (ftype==F_SNOWFALL        ){} //derived- cant set directly
  else if (ftype==F_RAINFALL        ){} //derived- cant set directly
  else if (ftype==F_IRRIGATION      ){f.irrigation=val;}

  else if (ftype==F_TEMP_AVE        ){f.temp_ave=val;}
  else if (ftype==F_TEMP_DAILY_MIN  ){f.temp_daily_min=val;}
  else if (ftype==F_TEMP_DAILY_MAX  ){f.temp_daily_max=val;}
  else if (ftype==F_TEMP_DAILY_AVE  ){f.temp_daily_ave=val;}
  else if (ftype==F_TEMP_MONTH_MAX  ){f.temp_month_max=val;}
  else if (ftype==F_TEMP_MONTH_MIN  ){f.temp_month_min=val;}
  else if (ftype==F_TEMP_MONTH_AVE  ){f.temp_month_ave=val;}

  else if (ftype==F_TEMP_AVE_UNC    ){f.temp_ave_unc =val;}
  else if (ftype==F_TEMP_MAX_UNC    ){f.temp_max_unc =val;}
  else if (ftype==F_TEMP_MIN_UNC    ){f.temp_min_unc =val;}

  else if (ftype==F_AIR_DENS        ){f.air_dens=val;}
  else if (ftype==F_AIR_PRES        ){f.air_pres=val;}
  else if (ftype==F_REL_HUMIDITY    ){f.rel_humidity=val;}

  else if (ftype==F_CLOUD_COVER     ){f.cloud_cover=val;}
  else if (ftype==F_ET_RADIA        ){f.ET_radia=val;}
  else if (ftype==F_ET_RADIA_FLAT   ){f.ET_radia_flat=val; }
  else if (ftype==F_SW_RADIA        ){f.SW_radia=val;}
  else if (ftype==F_SW_RADIA_UNC    ){f.SW_radia_unc=val;}
  else if (ftype==F_SW_RADIA_SUBCAN ){f.SW_radia_subcan=val;}
  else if (ftype==F_SW_SUBCAN_NET   ){f.SW_subcan_net=val;}
  else if (ftype==F_SW_RADIA_NET    ){f.SW_radia_net=val;}
  else if (ftype==F_LW_INCOMING     ){f.LW_incoming=val;}
  else if (ftype==F_LW_RADIA_NET    ){f.LW_radia_net=val;}

  else if (ftype==F_DAY_LENGTH      ){f.day_length=val;}
  else if (ftype==F_DAY_ANGLE       ){f.day_angle=val;}

  else if (ftype==F_WIND_VEL        ){f.wind_vel=val;}

  else if (ftype==F_PET             ){f.PET=val;}
  else if (ftype==F_OW_PET          ){f.OW_PET=val;}
  else if (ftype==F_PET_MONTH_AVE   ){f.PET_month_ave=val;}

  else if (ftype==F_POTENTIAL_MELT  ){f.potential_melt=val;}

  else if (ftype==F_RECHARGE        ){f.recharge=val;}
  else if (ftype==F_PRECIP_TEMP     ){f.precip_temp=val;}
  else if (ftype==F_PRECIP_CONC     ){f.precip_conc=val;}

  else if (ftype==F_SUBDAILY_CORR   ){f.subdaily_corr=val;}

  //else if (ftype==F_UNRECOGNIZED   ){0=val;}

  else
  {
    ExitGracefully("GetForcingFromType: invalid forcing type",RUNTIME_ERR);
  }
}
/////////////////////////////////////////////////////////////////////
/// \brief Return double value of forcing function specified by passed string parameter of force structure f
///
/// \param &forcing_string [in] String value of forcing function
/// \param &f [out] HRU forcing functions structure
/// \return Double forcing function value corresponding to &forcing_string
//
double GetForcingFromString(const string &forcing_string, const force_struct &f)
{
  forcing_type ftype;
  ftype=GetForcingTypeFromString(forcing_string);
  return GetForcingFromType(ftype,f);
}

/////////////////////////////////////////////////////////////////////
/// \brief Return a string containing the units of the given forcing type
///
/// \param &ftype [in] Ientifier of the forcing function type
/// \return String containing forcing function units
//
string GetForcingTypeUnits(forcing_type ftype)
{
  string units = "none";
  switch(ftype)
  {
  case F_PRECIP:          {units="mm/d"; break;}
  case F_PRECIP_DAILY_AVE:{units="mm/d"; break; }
  case F_PRECIP_5DAY:     {units="mm";   break;}
  case F_SNOW_FRAC:       {units="0-1";  break;}
  case F_SNOWFALL:        {units="mm/d"; break;}
  case F_RAINFALL:        {units="mm/d"; break;}
  case F_IRRIGATION:      {units="mm/d"; break;}

  case F_TEMP_AVE:        {units="C"; break;}
  case F_TEMP_DAILY_MIN:  {units="C"; break;}
  case F_TEMP_DAILY_MAX:  {units="C"; break;}
  case F_TEMP_DAILY_AVE:  {units="C"; break;}
  case F_TEMP_MONTH_MAX:  {units="C"; break;}
  case F_TEMP_MONTH_MIN:  {units="C"; break;}
  case F_TEMP_MONTH_AVE:  {units="C"; break;}
  case F_TEMP_AVE_UNC:    {units="C"; break;}
  case F_TEMP_MIN_UNC:    {units="C"; break;}
  case F_TEMP_MAX_UNC:    {units="C"; break;}

  case F_AIR_DENS:        {units="kg/m3"; break;}
  case F_AIR_PRES:        {units="kPa"; break;}
  case F_REL_HUMIDITY:    {units="0-1"; break;}

  case F_CLOUD_COVER:     {units="0-1"; break;}
  case F_ET_RADIA:        {units="MJ/m2/d"; break;}
  case F_ET_RADIA_FLAT:   {units="MJ/m2/d"; break;}
  case F_SW_RADIA:        {units="MJ/m2/d"; break;}
  case F_SW_RADIA_UNC:    {units="MJ/m2/d"; break;}
  case F_SW_RADIA_SUBCAN: {units="MJ/m2/d"; break;}
  case F_SW_SUBCAN_NET:   {units="MJ/m2/d"; break;}
  case F_SW_RADIA_NET:    {units="MJ/m2/d"; break;}
  case F_LW_INCOMING:     {units="MJ/m2/d"; break;}
  case F_LW_RADIA_NET:    {units="MJ/m2/d"; break;}

  case F_DAY_LENGTH:      {units="d"; break;}
  case F_DAY_ANGLE:       {units="rad"; break;}

  case F_WIND_VEL:        {units="m/s"; break;}

  case F_PET:             {units="mm/d"; break;}
  case F_OW_PET:          {units="mm/d"; break;}
  case F_PET_MONTH_AVE:   {units="mm/d"; break;}

  case F_POTENTIAL_MELT:  {units="mm/d"; break;}

  case F_RECHARGE:        {units="mm/d"; break;}
  case F_PRECIP_TEMP:     {units="C";    break;}
  case F_PRECIP_CONC:     {units="mg/L"; break;}

  case F_SUBDAILY_CORR:   {units="none"; break;}
  default:
    //ExitGracefully("GetForcingTypeUnits: invalid forcing string",RUNTIME_ERR);
    break;
  }
  return units;
}
/////////////////////////////////////////////////////////////////////
/// \brief Return a string equivalent of the given forcing type
///
/// \param &ftype [in] Identifier of the forcing function type
/// \return String containing forcing function
//
string ForcingToString(const forcing_type ftype)
{
  string fstring = "none";
  switch(ftype)
  {
  case F_PRECIP:          {fstring="PRECIP"; break;}
  case F_PRECIP_DAILY_AVE:{fstring="PRECIP_DAILY_AVE"; break; }
  case F_PRECIP_5DAY:     {fstring="PRECIP_5DAY"; break;}
  case F_SNOW_FRAC:       {fstring="SNOW_FRAC"; break;}
  case F_SNOWFALL:        {fstring="SNOWFALL"; break;}
  case F_RAINFALL:        {fstring="RAINFALL"; break;}
  case F_IRRIGATION:      {fstring="IRRIGATION"; break;}

  case F_TEMP_AVE:        {fstring="TEMP_AVE"; break;}
  case F_TEMP_DAILY_MIN:  {fstring="TEMP_DAILY_MIN"; break;}
  case F_TEMP_DAILY_MAX:  {fstring="TEMP_DAILY_MAX"; break;}
  case F_TEMP_DAILY_AVE:  {fstring="TEMP_DAILY_AVE"; break;}
  case F_TEMP_MONTH_MAX:  {fstring="TEMP_MONTH_MAX"; break;}
  case F_TEMP_MONTH_MIN:  {fstring="TEMP_MONTH_MIN"; break;}
  case F_TEMP_MONTH_AVE:  {fstring="TEMP_MONTH_AVE"; break;}
  case F_TEMP_AVE_UNC:    {fstring="TEMP_AVE_UNC"; break;}
  case F_TEMP_MIN_UNC:    {fstring="TEMP_MIN_UNC"; break;}
  case F_TEMP_MAX_UNC:    {fstring="TEMP_MAX_UNC"; break;}

  case F_AIR_DENS:        {fstring="AIR_DENSITY"; break;}
  case F_AIR_PRES:        {fstring="AIR_PRES"; break;}
  case F_REL_HUMIDITY:    {fstring="REL_HUMIDITY"; break;}

  case F_CLOUD_COVER:     {fstring="CLOUD_COVER"; break;}
  case F_ET_RADIA:        {fstring="ET_RADIA"; break;}
  case F_ET_RADIA_FLAT:   {fstring="ET_RADIA_FLAT";break;}
  case F_SW_RADIA:        {fstring="SW_RADIA"; break;}
  case F_SW_RADIA_UNC:    {fstring="SW_RADIA_UNC"; break;}
  case F_SW_RADIA_SUBCAN: {fstring="SW_RADIA_SUBCAN"; break;}
  case F_SW_SUBCAN_NET:   {fstring="SW_SUBCAN_NET"; break;}
  case F_LW_INCOMING:     {fstring="LW_INCOMING"; break;}
  case F_SW_RADIA_NET:    {fstring="SW_RADIA_NET"; break;}
  case F_LW_RADIA_NET:    {fstring="LW_RADIA_NET"; break;}

  case F_DAY_LENGTH:      {fstring="DAY_LENGTH"; break;}
  case F_DAY_ANGLE:       {fstring="DAY_ANGLE"; break;}

  case F_WIND_VEL:        {fstring="WIND_VEL"; break;}

  case F_PET:             {fstring="PET"; break;}
  case F_OW_PET:          {fstring="OW_PET"; break;}
  case F_PET_MONTH_AVE:   {fstring="PET_MONTH_AVE"; break;}

  case F_POTENTIAL_MELT:  {fstring="POTENTIAL_MELT"; break;}

  case F_RECHARGE:        {fstring="RECHARGE"; break;}
  case F_PRECIP_TEMP:     {fstring="PRECIP_TEMP"; break; }
  case F_PRECIP_CONC:     {fstring="PRECIP_CONC"; break; }

  case F_SUBDAILY_CORR:   {fstring="SUBDAILY_CORR"; break;}
  default:
    break;
  }
  return fstring;
}
