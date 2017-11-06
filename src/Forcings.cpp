/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
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
  F.SW_radia_net=0.0;
  F.LW_radia=0.0;

  F.day_length=0.0;
  F.day_angle=0.0;

  F.wind_vel=0.0;

  F.PET=0.0;
  F.OW_PET=0.0;
  F.PET_month_ave=0.0;

  F.potential_melt=0.0;

  F.subdaily_corr=0.0;
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

  else if (f=="TEMP_AVE"         ){return F_TEMP_AVE;}
  else if (f=="TEMP_MIN"         ){return F_TEMP_DAILY_MIN;}
  else if (f=="MIN_TEMPERATURE"  ){return F_TEMP_DAILY_MIN;}
  else if (f=="TEMP_DAILY_MIN"   ){return F_TEMP_DAILY_MIN;}
  else if (f=="TEMP_MAX"         ){return F_TEMP_DAILY_MAX;}
  else if (f=="MAX_TEMPERATURE"  ){return F_TEMP_DAILY_MAX;}
  else if (f=="TEMP_DAILY_MAX"   ){return F_TEMP_DAILY_MAX;}
  else if (f=="TEMP_DAILY_AVE"   ){return F_TEMP_DAILY_AVE;}
  else if (f=="TEMP_MONTH_MAX"   ){return F_TEMP_MONTH_MAX;}
  else if (f=="TEMP_MONTH_MIN"   ){return F_TEMP_MONTH_MIN;}
  else if (f=="TEMP_MONTH_AVE"   ){return F_TEMP_MONTH_AVE;}
  else if (f=="TEMP_AVE_UNC"     ){return F_TEMP_AVE_UNC;}
  else if (f=="TEMP_MAX_UNC"     ){return F_TEMP_MAX_UNC;}
  else if (f=="TEMP_MIN_UNC"     ){return F_TEMP_MIN_UNC;}

  else if (f=="AIR_DENS"         ){return F_AIR_DENS;}
  else if (f=="AIR_PRES"         ){return F_AIR_PRES;}
  else if (f=="REL_HUMIDITY"     ){return F_REL_HUMIDITY;}

  else if (f=="ET_RADIA"         ){return F_ET_RADIA;}
  else if (f=="SW_RADIA"         ){return F_SW_RADIA;}
  else if (f=="SHORTWAVE"        ){return F_SW_RADIA;}
  else if (f=="SW_RADIA_NET"     ){return F_SW_RADIA_NET;}
  else if (f=="LW_RADIA"         ){return F_LW_RADIA;}
  else if (f=="LONGWAVE"         ){return F_LW_RADIA;}
  else if (f=="CLOUD_COVER"      ){return F_CLOUD_COVER;}

  else if (f=="DAY_LENGTH"       ){return F_DAY_LENGTH;}
  else if (f=="DAY_ANGLE"        ){return F_DAY_ANGLE;}

  else if (f=="WIND_VEL"         ){return F_WIND_VEL;}

  else if (f=="PET"              ){return F_PET;}
  else if (f=="OW_PET"           ){return F_OW_PET;}
  else if (f=="PET_MONTH_AVE"    ){return F_PET_MONTH_AVE;}

  else if (f=="POTENTIAL_MELT"   ){return F_POTENTIAL_MELT;}

  else if (f=="SUBDAILY_CORR"    ){return F_SUBDAILY_CORR;}

  else
  {
    return F_UNRECOGNIZED;
    //cout <<"Forcing string:|"<<f<<"|"<<endl;
    //ExitGracefully("GetForcingTypeFromString: invalid forcing string",RUNTIME_ERR);
  }
  return F_UNRECOGNIZED;
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

  if      (ftype==F_PRECIP          ){return f.precip;}
  else if (ftype==F_PRECIP_DAILY_AVE){return f.precip_daily_ave;}
  else if (ftype==F_PRECIP_5DAY     ){return f.precip_5day;}
  else if (ftype==F_SNOW_FRAC       ){return f.snow_frac;}
  else if (ftype==F_SNOWFALL        ){return (    f.snow_frac)*f.precip;}
  else if (ftype==F_RAINFALL        ){return (1.0-f.snow_frac)*f.precip;}

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
  else if (ftype==F_SW_RADIA        ){return f.SW_radia;}
  else if (ftype==F_SW_RADIA_UNC    ){return f.SW_radia_unc;}
  else if (ftype==F_SW_RADIA_NET    ){return f.SW_radia_net;}
  else if (ftype==F_LW_RADIA        ){return f.LW_radia;}

  else if (ftype==F_DAY_LENGTH      ){return f.day_length;}
  else if (ftype==F_DAY_ANGLE       ){return f.day_angle;}

  else if (ftype==F_WIND_VEL        ){return f.wind_vel;}

  else if (ftype==F_PET             ){return f.PET;}
  else if (ftype==F_OW_PET          ){return f.OW_PET;}
  else if (ftype==F_PET_MONTH_AVE   ){return f.PET_month_ave;}

  else if (ftype==F_POTENTIAL_MELT  ){return f.potential_melt;}

  else if (ftype==F_SUBDAILY_CORR   ){return f.subdaily_corr;}

  //else if (ftype==F_UNRECOGNIZED   ){return 0;}

  else
  {
    ExitGracefully("GetForcingFromString: invalid forcing string",RUNTIME_ERR);
  }
  return 0.0;
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
  case F_PRECIP_5DAY:     {units="mm"; break;}
  case F_SNOW_FRAC:       {units="0-1"; break;}
  case F_SNOWFALL:        {units="mm"; break;}
  case F_RAINFALL:        {units="mm"; break;}

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
  case F_SW_RADIA:        {units="MJ/m2/d"; break;}
  case F_SW_RADIA_NET:    {units="MJ/m2/d"; break;}
  case F_LW_RADIA:        {units="MJ/m2/d"; break;}

  case F_DAY_LENGTH:      {units="d"; break;}
  case F_DAY_ANGLE:       {units="rad"; break;}

  case F_WIND_VEL:        {units="m/s"; break;}

  case F_PET:             {units="mm/d"; break;}
  case F_OW_PET:          {units="mm/d"; break;}
  case F_PET_MONTH_AVE:   {units="mm/d"; break;}

  case F_POTENTIAL_MELT:  {units="mm/d"; break;}

  case F_SUBDAILY_CORR:   {units="none"; break;}
  default:
    //ExitGracefully("GetForcingFromString: invalid forcing string",RUNTIME_ERR);
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
  case F_SW_RADIA:        {fstring="SW_RADIA"; break;}
  case F_SW_RADIA_NET:    {fstring="SW_RADIA_NET"; break;}
  case F_LW_RADIA:        {fstring="LW_RADIA"; break;}

  case F_DAY_LENGTH:      {fstring="DAY_LENGTH"; break;}
  case F_DAY_ANGLE:       {fstring="DAY_ANGLE"; break;}

  case F_WIND_VEL:        {fstring="WIND_VEL"; break;}

  case F_PET:             {fstring="PET"; break;}
  case F_OW_PET:          {fstring="OW_PET"; break;}
  case F_PET_MONTH_AVE:   {fstring="PET_MONTH_AVE"; break;}

  case F_POTENTIAL_MELT:  {fstring="POTENTIAL_MELT"; break;}

  case F_SUBDAILY_CORR:   {fstring="SUBDAILY_CORR"; break;}
  default:
    break;
  }
  return fstring;
}
