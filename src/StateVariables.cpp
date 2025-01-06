/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2023 the Raven Development Team
  ----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "StateVariables.h"
#include "Transport.h"

//--Initialize Static Variables-----------------------------------
// int     CStateVariable::_nAliases        =0;
// string *CStateVariable::_aAliases        =NULL;
// string *CStateVariable::_aAliasReferences=NULL;

//////////////////////////////////////////////////////////////////
/// Constructor
//
CStateVariable::CStateVariable()
{
  this->Initialize();
  _pTransportModel=NULL;
}

//////////////////////////////////////////////////////////////////
/// \brief Initializes static arrays of CStateVariable class
//
void CStateVariable::Initialize()
{
  _nAliases = 0;
  _aAliases = NULL;
  _aAliasReferences = NULL;
}

//////////////////////////////////////////////////////////////////
/// Destructor
//
CStateVariable::~CStateVariable()
{
  delete [] _aAliases;
  delete [] _aAliasReferences;
}

//////////////////////////////////////////////////////////////////
/// \brief Delete references to static alias arrays
//
void CStateVariable::Destroy(){}

//////////////////////////////////////////////////////////////////
/// \brief Dynamically adds additional string, s, onto dynamic array of strings, pArr. Increments size of array by one
///
/// \param *&pArr [in & out] Array of strings to which a new string will be appended
/// \param s [in] String to be appended to **pArr
/// \param &size [in & out] Integer size of array of strings
/// \return Boolean variable to indicate success of appending
//
bool StrArrayAppend(string *& pArr, string s,int &size)
{
  string *tmp=NULL;
  if ((pArr==NULL) && (size>0)) {return false;}
  size=size+1;                                                //increment size
  tmp=new string [size+1];                                      //allocate memory
  if (tmp==NULL){ExitGracefully("StrArrayAppend::Out of memory",OUT_OF_MEMORY);}
  for (int i=0; i<(size-1); i++){                             //copy array
    tmp[i]=pArr[i];
  }
  tmp[size-1]=s;                                              //add new string
  if (size>1){delete [] pArr; pArr=NULL;}                     //delete old array of pointers
  pArr=tmp;                                                   //redirect pointer
  return true;
}

///////////////////////////////////////////////////////////////////
/// \brief Adds an alias to static array of aliases and updates array of references to aliases
/// \param s1 [in] Alias to be appended to the static array of aliases
/// \param s2 [in] String referenced by the new alias
//
void CStateVariable::AddAlias(const string s1, const string s2)
{
  StrArrayAppend(_aAliases,         s1, _nAliases);
  _nAliases--;
  StrArrayAppend(_aAliasReferences, s2, _nAliases);
}

//////////////////////////////////////////////////////////////////
/// \brief Sets pointer to transport model
/// \param pTransportModel [in] Pointer to transport model
//
void CStateVariable::SetTransportModel(CTransportModel *pTransportModel)
{
  this->_pTransportModel = pTransportModel;
}

//////////////////////////////////////////////////////////////////
/// \brief Checks if input string is in alias list
/// \details If the input string is in the alias list, then it returns the
///  string which it is an alias for, otherwise returns the unchanged input
///
/// \param s [in] Input string
/// \return Value for which input string is an alias for, or input string if this value doesn't exist
//
string CStateVariable::CheckAliasList(const string s)
{
  string sup=StringToUppercase(s);
  for (int i=0;i<_nAliases;i++){
    if (!sup.compare(StringToUppercase(_aAliases[i])))
    {
      return _aAliasReferences[i];
    }
  }
  return s;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns a string describing the state variable
/// \param typ [in] Type of state variable
/// \param layerindex [in] Layer index (when applicable)
/// \return String describing state variable
//
string CStateVariable::GetStateVarLongName(const sv_type typ, const int layerindex,
                                           const CTransportModel *pTransportModel)
{

  string name;
  switch(typ)
  {
    //Water Storage [mm]
  case(SURFACE_WATER):      {name="Surface Water";              break;}
  case(ATMOSPHERE):         {name="Cum. Losses to Atmosphere";  break;}
  case(ATMOS_PRECIP):       {name="Cum. Precipitation";         break;}
  case(PONDED_WATER):       {name="Ponded Water";               break;}

  case(SOIL):               {name="Soil Water";                 break;}
  case(CANOPY):             {name="Canopy";                     break;}
  case(CANOPY_SNOW):        {name="Canopy Snow";                break;}
  case(TRUNK):              {name="Trunks of trees";            break;}
  case(ROOT):               {name="Root Storage";               break;}
  case(GROUNDWATER):        {name="Deep Groundwater";           break;}
  case(DEPRESSION):         {name="Depression";                 break;}
  case(SNOW):               {name="Snow";                       break;}
  case(NEW_SNOW):           {name="New Snow";                   break;}
  case(SNOW_LIQ):           {name="Snow Melt (Liquid)";         break;}
  case(TOTAL_SWE):          {name="Total SWE";                  break;}
  case(WETLAND):            {name="Wetlands";                   break;}
  case(CUM_INFIL):          {name="Cumulative infiltration";    break;}
  case(GA_MOISTURE_INIT):   {name="Green Ampt initial soil Water"; break;}
  case(LATERAL_EXCHANGE):   {name="Lateral exchange storage";   break;}
  case(SNOW_DRIFT):         {name="Blowing Snow";               break;}
  case(LAKE_STORAGE):       {name="Net Lake Storage";           break;}
  case(MIN_DEP_DEFICIT):    {name="Minimum depression deficit"; break;}

  //Temperature/Energy storage
  case(ENERGY_LOSSES):      {name="Energy Losses";              break;}
  case(SURFACE_WATER_TEMP): {name="Surface Water Temperature";  break;}
  case(SNOW_TEMP):          {name="Temperature of snow";        break;}
  case(COLD_CONTENT):       {name="Cold Content";               break;}
  case(SOIL_TEMP):          {name="Temperature of soil";        break;}
  case(CANOPY_TEMP):        {name="Temperature of canopy";      break;}

    //Snow variables
  case(SNOW_DEPTH):         {name="Snow depth";                 break;}
  case(PERMAFROST_DEPTH):   {name="Permafrost depth";           break;}
  case(THAW_DEPTH):         {name="Ground thaw depth";          break;}
  case(SNOW_DEPTH_STDDEV):  {name="Snow depth sigma";           break;}
  case(SNOW_COVER):         {name="Fractional Snow Cover";      break;}
  case(CUM_SNOWMELT):       {name="Cumulative Snowmelt";        break;}
  case(SNOW_DEFICIT):       {name="Snow Deficit";               break;}
  case(ICE_THICKNESS):      {name="Ice Thickness";              break;}

  case(SNOW_AGE):           {name="Snow Age";                   break;}
  case(SNODRIFT_TEMP):      {name="Blowing Snow Temperature";   break;}

  case(GLACIER):            {name="Glacier Liquid Storage";     break;}
  case(GLACIER_ICE):        {name="Glacier Ice";                break;}
  case(GLACIER_CC):         {name="Glacier Cold Content";       break;}

  case(SNOW_ALBEDO):        {name="Snow Albedo";                break;}
  case(CROP_HEAT_UNITS):    {name="Crop Heat Units";            break;}

  case(CONVOLUTION):        {name="Convolution Storage";        break;}
  case(CONV_STOR):          {name="Convolution Substorage";     break;}

  case(AET):                {name="Actual Evapotranspiration";  break;}
  case(RUNOFF):             {name="Runoff";                     break;}
  case(STREAMFLOW):         {name="Streamflow";                 break;}
  case(RESERVOIR_STAGE):    {name="Reservoir Stage";            break;}

    //Transport variables
  case(CONSTITUENT):        {name="Constituent";                break;}  //overwritten below
  case(CONSTITUENT_SRC):    {name="Constituent Source";         break;}  //overwritten below
  case(CONSTITUENT_SINK):   {name="Constituent Sink";           break;}  //overwritten below
  case(CONSTITUENT_SW):     {name="Constituent in Surface Water"; break;}  //overwritten below
    //..
  default:
  {
    cout<<typ<<endl;
    name="Unknown State Variable Type";
    ExitGracefully("CStateVariable::GetStateVarLongName: Unknown State Variable Type",BAD_DATA);//STRICT
    break;
  }
  }
  if ((typ==SOIL) || (typ==SOIL_TEMP) ||
      (typ==CONVOLUTION) || (typ==CONV_STOR) || (typ==LATERAL_EXCHANGE))
  {
    name=name+"["+to_string(layerindex)+"]";
  }
  if (((typ==SNOW) || (typ==GROUNDWATER))  && (layerindex>=0)){  //GROUNDWATER Removed LS 2/25/2020
    name = name + "["+to_string(layerindex)+"]";
  }
  if (typ==CONSTITUENT){
    name = pTransportModel->GetConstituentLongName(layerindex);
  }
  else if (typ==CONSTITUENT_SRC){
    int c=layerindex;
    name = "Source of " + pTransportModel->GetConstituentTypeName2(c);
  }
  else if (typ==CONSTITUENT_SINK){
    int c=layerindex;
    name = "Sink of " + pTransportModel->GetConstituentTypeName2(c);
  }
  else if (typ==CONSTITUENT_SW){
    int c=layerindex;
    name = "SW Sink of " + pTransportModel->GetConstituentTypeName2(c);
  }
  return name;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns a string containing the state variable units
/// \param typ [in] Type of state variable
/// \return String containing state variable units
//
string CStateVariable::GetStateVarUnits(const sv_type typ)
{
  string units;
  switch(typ)
  {
    //Water Storage [mm]
  case(SURFACE_WATER):    {units="mm"; break;}
  case(ATMOSPHERE):       {units="mm"; break;}
  case(ATMOS_PRECIP):     {units="mm"; break;}
  case(PONDED_WATER):     {units="mm"; break;}

  case(SOIL):             {units="mm"; break;}
  case(CANOPY):           {units="mm"; break;}
  case(CANOPY_SNOW):      {units="mm"; break;}
  case(TRUNK):            {units="mm"; break;}
  case(ROOT):             {units="mm"; break;}
  case(GROUNDWATER):      {units="mm"; break;}
  case(DEPRESSION):       {units="mm"; break;}
  case(SNOW):             {units="mm"; break;}
  case(NEW_SNOW):         {units="mm"; break;}
  case(SNOW_LIQ):         {units="mm"; break;}
  case(TOTAL_SWE):        {units="mm"; break;}
  case(WETLAND):          {units="mm"; break;}
  case(CUM_INFIL):        {units="mm"; break;}
  case(GA_MOISTURE_INIT): {units="mm"; break;}
  case(LATERAL_EXCHANGE): {units="mm"; break;}
  case(SNOW_DRIFT):       {units="mm"; break;}
  case(LAKE_STORAGE):     {units="mm"; break;}
  case(MIN_DEP_DEFICIT):  {units="mm"; break;}

    //Temperature/Energy storage [C] or [MJ/m^2]
  case(ENERGY_LOSSES):    {units="MJ/m2"; break;}
  case(SURFACE_WATER_TEMP):{units="C";  break;}
  case(SNOW_TEMP):        {units="C";   break;}
  case(COLD_CONTENT):     {units="C";   break;}
  case(SOIL_TEMP):        {units="C";   break;}
  case(CANOPY_TEMP):      {units="C";   break;}

    //Snow variables
  case(SNOW_DEPTH):       {units="mm";   break;}
  case(PERMAFROST_DEPTH): {units="mm";   break;}
  case(THAW_DEPTH):       {units="mm";   break;}
  case(SNOW_DEPTH_STDDEV):{units="log(mm)"; break;}
  case(SNOW_COVER):       {units="0-1";  break;}
  case(CUM_SNOWMELT):     {units="mm";   break;}
  case(SNOW_DEFICIT):     {units="mm";   break;}
  case(SNOW_AGE):         {units="d";    break;}
  case(SNODRIFT_TEMP):    {units="C";    break;}
  case(ICE_THICKNESS):    {units="mm";   break;}

  case(GLACIER):          {units="mm";   break;}
  case(GLACIER_ICE):      {units="mm";   break;}
  case(GLACIER_CC):       {units="mm";   break;}

  case(SNOW_ALBEDO):      {units="none"; break;}

  case(CROP_HEAT_UNITS):  {units="none"; break;}

  case(CONVOLUTION):      {units="mm";   break;}
  case(CONV_STOR):        {units="mm";   break;}

  case(AET):              {units="mm";   break;}
  case(RUNOFF):           {units="mm";   break;}
  case(STREAMFLOW):       {units="m3/s"; break;}
  case(RESERVOIR_STAGE):  {units="m";    break;}

  case(CONSTITUENT):      {units="mg/m2"; break;}
  case(CONSTITUENT_SRC):  {units="mg/m2"; break;}
  case(CONSTITUENT_SINK): {units="mg/m2"; break;}
  case(CONSTITUENT_SW):   {units="mg/m2"; break;}
    //..
  default:
  {
    cout<<typ<<endl;
    units="Unknown State Variable units";
    ExitGracefully("CStateVariable::GetStateVarUnits: Unknown State Variable Type",BAD_DATA);//STRICT
    break;
  }
  }

  return units;
}


//////////////////////////////////////////////////////////////////
/// \brief Converts a string parameter to an enumerated state variable type
/// \details If the string is an array item (e.g. "SOIL[2]"), also returns the index.
/// Otherwise this index = DOESNT_EXIST
/// \note if the string begins with '!' it is assumed to be a constituent
///
/// \param s [in] String-formatted state variable type
/// \param &layer_index [out] Array index corresponding to a layer
/// \param strict [in] Flag for if state variable type must be recognized to continue operation of program
sv_type CStateVariable::StringToSVType(const string s, int &layer_index,bool strict)
{
  sv_type typ;
  string stmp;
  stmp = StringToUppercase(s);
  stmp = this->CheckAliasList(stmp);//checks for alias

  layer_index=0;//default
  string tmp=SVStringBreak(stmp,layer_index);

  if      (!tmp.compare("SURFACEWATER"    )){typ=SURFACE_WATER;}
  else if (!tmp.compare("SURFACE_WATER"   )){typ=SURFACE_WATER;}
  else if (!tmp.compare("PONDED_WATER"    )){typ=PONDED_WATER;}
  else if (!tmp.compare("ATMOSPHERE"      )){typ=ATMOSPHERE;}
  else if (!tmp.compare("ATMOS_PRECIP"    )){typ=ATMOS_PRECIP;}
  else if (!tmp.compare("SOIL"            )){typ=SOIL;}
  else if (!tmp.compare("GROUNDWATER"     )){typ=GROUNDWATER;}
  else if (!tmp.compare("GROUND_WATER"    )){typ=GROUNDWATER;}
  else if (!tmp.compare("GWMODEL"         )){typ=GROUNDWATER;}
  else if (!tmp.compare("AQUIFER"         )){typ=GROUNDWATER;}
  else if (!tmp.compare("CANOPY"          )){typ=CANOPY;}
  else if (!tmp.compare("CANOPYSNOW"      )){typ=CANOPY_SNOW;}
  else if (!tmp.compare("CANOPY_SNOW"     )){typ=CANOPY_SNOW;}
  else if (!tmp.compare("TRUNK"           )){typ=TRUNK;}
  else if (!tmp.compare("ROOT"            )){typ=ROOT;}
  else if (!tmp.compare("SNOW"            )){typ=SNOW;}
  else if (!tmp.compare("NEW_SNOW"        )){typ=NEW_SNOW;}
  else if (!tmp.compare("SNOWLIQ"         )){typ=SNOW_LIQ;}
  else if (!tmp.compare("SNOW_LIQ"        )){typ=SNOW_LIQ;}
  else if (!tmp.compare("TOTAL_SWE"       )){typ=TOTAL_SWE;}
  else if (!tmp.compare("SNOW_DEPTH"      )){typ=SNOW_DEPTH;}
  else if (!tmp.compare("SNOW_TEMP"       )){typ=SNOW_TEMP;}
  else if (!tmp.compare("COLD_CONTENT"    )){typ=COLD_CONTENT;}
  else if (!tmp.compare("WETLAND"         )){typ=WETLAND;}
  else if (!tmp.compare("DEPRESSION"      )){typ=DEPRESSION;}
  else if (!tmp.compare("LAKE_STORAGE"    )){typ=LAKE_STORAGE;}
  else if (!tmp.compare("MIN_DEP_DEFICIT" )){typ=MIN_DEP_DEFICIT;}
  else if (!tmp.compare("ENERGY_LOSSES"   )){typ=ENERGY_LOSSES;}
  else if (!tmp.compare("SNOW_COVER"      )){typ=SNOW_COVER;}
  else if (!tmp.compare("SNOW_DEFICIT"    )){typ=SNOW_DEFICIT;}
  else if (!tmp.compare("SNOW_AGE"        )){typ=SNOW_AGE;}
  else if (!tmp.compare("SNODRIFT_TEMP"   )){typ=SNODRIFT_TEMP;}
  else if (!tmp.compare("SNOW_DRIFT"      )){typ=SNOW_DRIFT;}
  else if (!tmp.compare("ICE_THICKNESS"   )){typ=ICE_THICKNESS;}
  else if (!tmp.compare("GLACIER"         )){typ=GLACIER;}
  else if (!tmp.compare("GLACIER_ICE"     )){typ=GLACIER_ICE;}
  else if (!tmp.compare("GLACIER_CC"      )){typ=GLACIER_CC;}
  else if (!tmp.compare("GLACIER_COLD_CONTENT")){typ=GLACIER_CC;}
  else if (!tmp.compare("CUM_SNOWMELT"    )){typ=CUM_SNOWMELT;}
  else if (!tmp.compare("SNOW_TEMP"       )){typ=SNOW_TEMP;}
  else if (!tmp.compare("SNOW_ALBEDO"     )){typ=SNOW_ALBEDO;}
  else if (!tmp.compare("CROP_HEAT_UNITS" )){typ=CROP_HEAT_UNITS;}
  else if (!tmp.compare("CUM_INFIL"       )){typ=CUM_INFIL;}
  else if (!tmp.compare("CONVOLUTION"     )){typ=CONVOLUTION;}
  else if (!tmp.compare("CONV_STOR"       )){typ=CONV_STOR;}
  else if (!tmp.compare("GA_MOISTURE_INIT")){typ=GA_MOISTURE_INIT;}
  else if (!tmp.compare("LATERAL_EXCHANGE")){typ=LATERAL_EXCHANGE;}

  else if (!tmp.compare("AET"             )){typ=AET; }
  else if (!tmp.compare("RUNOFF"          )){typ=RUNOFF; }
  else if (!tmp.compare("MULTIPLE"        )){typ=MULTIPLE_SVTYPE;}
  else if (!tmp.compare("STREAMFLOW"      )){typ=STREAMFLOW;}
  else if (!tmp.compare("RESERVOIR_STAGE" )){typ=RESERVOIR_STAGE;}

  else if (!tmp.compare("CONSTITUENT"     )){typ=CONSTITUENT;}
  else if (!tmp.compare("CONSTITUENT_SRC" )){typ=CONSTITUENT_SRC;}
  else if (!tmp.compare("CONSTITUENT_SINK")){typ=CONSTITUENT_SINK;}
  else if (!tmp.compare("CONSTITUENT_SW"  )){typ=CONSTITUENT_SW;}
  else if (tmp.c_str()[0]=='!'             ){
    if      (!tmp.substr(tmp.length()-4,4).compare("_SRC" )){typ=CONSTITUENT_SRC;}
    else if (!tmp.substr(tmp.length()-5,5).compare("_SINK")){typ=CONSTITUENT_SINK;}
    else if (!tmp.substr(tmp.length()-3,3).compare( "_SW" )){typ=CONSTITUENT_SW;}
    else                                                    {typ=CONSTITUENT;}
  }
  else                                      {typ=UNRECOGNIZED_SVTYPE;}

  if ((typ==CONSTITUENT) && ((int)(tmp.find_first_of("|"))!=-1)) //only used if e.g., !Nitrogen|SOIL[1] (rather than CONSTITUENT[32] or !Nitrogen[32]) is used
  {
    layer_index = _pTransportModel->GetLayerIndexFromName2(tmp, layer_index);
    if (layer_index==DOESNT_EXIST){typ=UNRECOGNIZED_SVTYPE;}
  }

  if ((strict) && (typ==UNRECOGNIZED_SVTYPE)){
    cout<<tmp<<endl;
    ExitGracefully("StringToSVtype: Unrecognized State Variable (2)",BAD_DATA);//strict
  }
  return typ;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns a string describing the state variable
/// \param typ [in] Type of state variable
/// \param layerindex [in] Layer index (when applicable)
/// \return String describing state variable
//
string CStateVariable::SVTypeToString(const sv_type typ, const int layerindex)
{

  string name;
  switch(typ)
  {
    //Water Storage [mm]
    case(SURFACE_WATER):      {name="SURFACE_WATER";            break;}
    case(ATMOSPHERE):         {name="ATMOSPHERE";               break;}
    case(ATMOS_PRECIP):       {name="ATMOS_PRECIP";             break;}
    case(PONDED_WATER):       {name="PONDED_WATER";             break;}
    case(SOIL):               {name="SOIL";                     break;}
    case(CANOPY):             {name="CANOPY";                   break;}
    case(CANOPY_SNOW):        {name="CANOPY_SNOW";              break;}
    case(TRUNK):              {name="TRUNK";                    break;}
    case(ROOT):               {name="ROOT";                     break;}
    case(GROUNDWATER):        {name="GROUNDWATER";              break;}
    case(DEPRESSION):         {name="DEPRESSION";               break;}
    case(SNOW):               {name="SNOW";                     break;}
    case(NEW_SNOW):           {name="NEW_SNOW";                 break;}
    case(SNOW_LIQ):           {name="SNOW_LIQ";                 break;}
    case(TOTAL_SWE):          {name="TOTAL_SWE";                break;}
    case(WETLAND):            {name="WETLAND";                  break;}
    case(CUM_INFIL):          {name="CUM_INFIL";                break;}
    case(GA_MOISTURE_INIT):   {name="GA_MOISTURE_INIT";         break;}
    case(SNOW_DRIFT):         {name="SNOW_DRIFT";               break;}
    case(LAKE_STORAGE):       {name="LAKE_STORAGE";             break;}
    case(MIN_DEP_DEFICIT):    {name="MIN_DEP_DEFICIT";          break;}

    //Temperature/Energy storage
    case(ENERGY_LOSSES):      {name="ENERGY_LOSSES";            break;}
    case(SURFACE_WATER_TEMP): {name="SURFACE_WATER_TEMP";       break;}
    case(SNOW_TEMP):          {name="SNOW_TEMP";                break;}
    case(COLD_CONTENT):       {name="COLD_CONTENT";             break;}
    case(SOIL_TEMP):          {name="SOIL_TEMP";                break;}
    case(CANOPY_TEMP):        {name="CANOPY_TEMP";              break;}

    //Snow variables
    case(SNOW_DEPTH):         {name="SNOW_DEPTH";               break;}
    case(THAW_DEPTH):         {name="THAW_DEPTH";               break;}
    case(SNOW_DEPTH_STDDEV):  {name="SNOW_DEPTH_STDDEV";        break;}
    case(SNOW_COVER):         {name="SNOW_COVER";               break;}
    case(CUM_SNOWMELT):       {name="CUM_SNOWMELT";             break;}
    case(SNOW_DEFICIT):       {name="SNOW_DEFICIT";             break;}
    case(SNOW_AGE):           {name="SNOW_AGE";                 break;}
    case(SNODRIFT_TEMP):      {name="SNODRIFT_TEMP";            break;}
    case(ICE_THICKNESS):      {name="ICE_THICKNESS";            break;}

    //Glacier variables
    case(GLACIER):            {name="GLACIER";                  break;}
    case(GLACIER_ICE):        {name="GLACIER_ICE";              break;}
    case(GLACIER_CC):         {name="GLACIER_CC";               break;}

    case(SNOW_ALBEDO):        {name="SNOW_ALBEDO";              break;}
    case(CROP_HEAT_UNITS):    {name="CROP_HEAT_UNITS";          break;}

    //Convolution Variables
    case(CONVOLUTION):        {name="CONVOLUTION";              break;}
    case(CONV_STOR):          {name="CONV_STOR";                break;}

    //Lateral exchange
    case(LATERAL_EXCHANGE):   {name="LATERAL_EXCHANGE";         break;}

    case(AET):                {name="AET";                      break;}
    case(RUNOFF):             {name="RUNOFF";                   break;}
    case(STREAMFLOW):         {name="STREAMFLOW";               break;}
    case(RESERVOIR_STAGE):    {name="RESERVOIR_STAGE";          break;}

    case(MULTIPLE_SVTYPE):    {name="MULTIPLE";                 break;}

    //Transport variables
    case(CONSTITUENT):    {
      name = "!" + this->_pTransportModel->GetConstituentTypeName(layerindex); //e.g., !Nitrogen
      break;
    }
    case(CONSTITUENT_SRC):    {
      int c=layerindex;
      name = "!" + this->_pTransportModel->GetConstituentTypeName2(c) + "_SRC"; //e.g., !Nitrogen_SRC
      break;
    }
    case(CONSTITUENT_SINK):    {
      int c=layerindex;
      name = "!" + this->_pTransportModel->GetConstituentTypeName2(c) + "_SINK"; //e.g., !Nitrogen_SINK
      break;
    }
    case(CONSTITUENT_SW):    {
      int c=layerindex;
      name = "!" + this->_pTransportModel->GetConstituentTypeName2(c) + "_SW"; //e.g., !Nitrogen_SW
      break;
    }
      //..
    default:
    {
      cout<<typ<<endl;
      name="Unknown State Variable Type";
      ExitGracefully("CStateVariable::SVTypeToString: Unknown State Variable Type",BAD_DATA);//STRICT
      break;
    }
  }
  //multilayer variables
  if ((typ==SOIL) || (typ==SOIL_TEMP) || (typ==CONSTITUENT)  ||
    (typ==CONVOLUTION) || (typ==CONV_STOR) || (typ==LATERAL_EXCHANGE))
  {
    name=name+"["+to_string(layerindex)+"]";
  }
  if (((typ==SNOW) || (typ==SNOW_LIQ) || (typ==COLD_CONTENT) || (typ==GROUNDWATER)) && (layerindex>0)){
    name=name+"["+to_string(layerindex)+"]";
  }
  return name;
}
//////////////////////////////////////////////////////////////////
/// \brief Utiltiy method for parsing strings of the form "NAME[12]"
/// \details Converts a string, s, of the form "NAME[12]" and returns the name
/// ("NAME") and number (num=12). If the string is of the form NAME (i.e.,
/// without a bracketed number), returns the original string "NAME"
/// with num=-1 (DOESNT_EXIST)
///
/// \param s [in] Input string of the form "NAME[12]" to be parsed
/// \param &num [out] Array index corresponding to a layer
/// \return String without array index (i.e., "NAME")
//
string CStateVariable::SVStringBreak(const string s, int &num)
{
  const char *ss=s.c_str();
  if (strrchr(ss,'[')==NULL){num=-1; return s;}
  const char *pch;
  const char *pch2;
  static char tmp[50];
  static char tmp2[50];
  memset(&tmp[0], 0, sizeof(tmp)); //clear arrays
  memset(&tmp2[0], 0, sizeof(tmp2));
  char key1[] = "[";
  char key2[] = "]";
  pch  = strpbrk (ss, key1); //pch is now pointer to '[' character
  pch2 = strpbrk (ss, key2); //pch is now pointer to ']' character
  if ((pch==NULL) || (pch2==NULL)){num=0;return s;}//one or both brackets missing - layer = 0
#if defined (__APPLE__) || defined (__CYGWIN__)
  strxfrm(tmp,ss,pch-ss+1);
  strxfrm(tmp2,pch+1,pch2-pch);
#else
  strxfrm(tmp,ss,pch-ss);          //Extract type (e.g., "SOIL" from "SOIL[13]")
  strxfrm(tmp2,pch+1,pch2-pch-1);  //Extract index (e.g., "13" from "SOIL[13]")
#endif
  //cout<<tmp2<<endl;
  for(int i = 0; i < (int)(strlen(tmp2)); ++i) {
    if(!isdigit(tmp2[i])) {
      string warn="SVStringBreak: non-integer ["+to_string(tmp2)+"] used in brackets of state variable in input file";
      ExitGracefully(warn.c_str(),BAD_DATA_WARN);
    }
  }
  num=s_to_i(tmp2);
  return string(tmp);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns true if the state variable is a water storage unit
/// \param typ [in] State variable identifier
/// \param conv_coverup [in] [default=true] - reports CONVOLUTION as water storage, but NOT CONV_STOR (for simplified reporting) - false for transport calcs
/// \return Boolean indicating if the passed state variable is a water storage unit
bool  CStateVariable::IsWaterStorage (sv_type      typ, bool conv_coverup)
{
  switch(typ)
  {
  case(SURFACE_WATER):   {return true;}
  case(PONDED_WATER):    {return true;}
  case(ATMOSPHERE):      {return true;}
  case(ATMOS_PRECIP):    {return true;}
  case(SNOW):            {return true;}
  case(GROUNDWATER):     {return true;}
  case(SOIL):            {return true;}
  case(CANOPY):          {return true;}
  case(CANOPY_SNOW):     {return true;}
  case(ROOT):            {return true;}
  case(DEPRESSION):      {return true;}
  case(SNOW_LIQ):        {return true;}
  case(WETLAND):         {return true;}
  case(GLACIER):         {return true;}
  case(GLACIER_ICE):     {return true;}
  case(CONVOLUTION):     {return conv_coverup;}
  case(CONV_STOR):       {return !conv_coverup;}
  case(NEW_SNOW):        {return true;}
  case(LATERAL_EXCHANGE):{return true;}
  case(SNOW_DRIFT):      {return true;}
  case(LAKE_STORAGE):    {return true;}
  //..
  default:
  {
    return false;
  }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns true if the state variable is an energy storage unit
/// \param typ [in] State variable identifier
/// \return Boolean indicating if the passed state variable is an energy storage unit
bool  CStateVariable::IsEnergyStorage (sv_type      typ)
{
  switch(typ)
  {
  case(COLD_CONTENT):   {return true;}
  case(ENERGY_LOSSES):  {return true;}
    //..
  default:
  {
    return false;
  }
  }
}
