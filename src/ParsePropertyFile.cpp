/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright © 2008-2014 the Raven Development Team
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Properties.h"
#include "Model.h"
#include "ParseLib.h"
#include "GlobalParams.h"
#include "SoilAndLandClasses.h"
#include <string>


bool ParsePropArray(CParser *p,int *indices,double **properties,
                    int &num_read,string *tags,const int line_length,const int max_classes);
void  RVPParameterWarning   (string  *aP, class_type *aPC, int &nP, const optStruct &Options);
void  AddToMasterParamList   (string        *&aPm, class_type       *&aPCm, int       &nPm, 
                              const string  *aP , const class_type *aPC , const int &nP); 
//////////////////////////////////////////////////////////////////
/// \brief This method parses the class properties .rvp file
///
/// \details Parses input \b model.rvp file that defines soil, land use, terrain, and vegetation classes and their properties \n
///   ":SoilClasses"  \n
///   {string tag, double pct_sand,double pct_clay,double pct_org}xNumSoilClasses  \n
///   :EndSoilClasses  \n
///   ":SoilProfiles"  \n
///   {string tag, int numhorizons, \n
///        soil1,thick1,soil2,thick2,...,soilN,thickN}xNumHorizons [where thickness is in meters]  \n
///   :EndSoilProfiles  \n
///   ":VegetationClasses"  \n
///   {string tag, max height [m], max LAI, max leaf cond.}xNumVegetationClasses  \n
///   :EndVegetationClasses  \n
///   ":LandUseClasses"  \n
///   {string tag, imperm. frac[-], SCS}xNumLUClasses  \n
///   :EndLandUseClasses  \n
///   ":TerrainClasses" \n
///
/// \param *&pModel [in] The input model object
/// \param &Options [in] Global model options information
/// \return Boolean variable indicating success of parsing
//
bool ParseClassPropertiesFile(CModel         *&pModel,
                              const optStruct &Options)
{

  global_struct     global_template;
  global_struct     parsed_globals;

  int               num_parsed_veg=0;
  CVegetationClass *pVegClasses[MAX_VEG_CLASSES];
  veg_struct        parsed_veg [MAX_VEG_CLASSES];
  string            vegtags    [MAX_VEG_CLASSES];

  int               num_parsed_lult=0;
  CLandUseClass    *pLUClasses  [MAX_LULT_CLASSES];
  surface_struct    parsed_surf [MAX_LULT_CLASSES];
  string            lulttags    [MAX_LULT_CLASSES];

  int               num_parsed_soils=0;
  CSoilClass       *pSoilClasses[MAX_SOIL_CLASSES];
  soil_struct       parsed_soils[MAX_SOIL_CLASSES];
  string            soiltags    [MAX_SOIL_CLASSES];

  int               num_parsed_terrs=0;
  CTerrainClass    *pTerrClasses[MAX_TERRAIN_CLASSES];
  terrain_struct    parsed_terrs[MAX_TERRAIN_CLASSES];
  string            terraintags [MAX_TERRAIN_CLASSES];

  int               num_parsed_profiles=0;
  CSoilProfile     *pProfiles   [MAX_SOIL_PROFILES];

  int               num_parsed_aqstacks=0;
  CAquiferStack    *pAqStacks   [MAX_AQUIFER_STACKS];

  const int         MAX_PROPERTIES_PER_LINE=20;

  int               MAX_NUM_IN_CLASS=max(max(MAX_SOIL_CLASSES,MAX_VEG_CLASSES),MAX_LULT_CLASSES);
  bool              invalid_index;
  int               num_read;
  int              *indices;
  double          **properties;
  string            aParamStrings[MAXINPUTITEMS];
  int               nParamStrings=0;


  

  CGlobalParams::InitializeGlobalParameters(global_template,true);
  CGlobalParams::InitializeGlobalParameters(parsed_globals,false); 
  
  CVegetationClass::InitializeVegetationProps (parsed_veg[0],true);//zero-index canopy is template
  vegtags    [0]="[DEFAULT]";
  num_parsed_veg++;

  CSoilClass::InitializeSoilProperties(parsed_soils[0],true);//zero-index soil is template - does not create class
  soiltags    [0]="[DEFAULT]";
  num_parsed_soils++;

  CLandUseClass::InitializeSurfaceProperties(parsed_surf[0],true);//zero-index LULT is template
  lulttags    [0]="[DEFAULT]";
  num_parsed_lult++;

  CTerrainClass::InitializeTerrainProperties(parsed_terrs[0],true);//zero-index terrain is template
  terraintags [0]="[DEFAULT]";
  num_parsed_terrs++;

  indices   =new int     [MAX_NUM_IN_CLASS];
  properties=new double *[MAX_NUM_IN_CLASS];
  for (int i=0;i<MAX_NUM_IN_CLASS;i++){
    properties[i]=new double [MAX_PROPERTIES_PER_LINE];
  }

  const int MAX_PARAMS=100;
  int         nP;
  string     *aP  =new string [MAX_PARAMS];
  class_type *aPC =new class_type [MAX_PARAMS];
  
  int         nPmaster(0);   //number of parameters in master list of required parameters
  string     *aPmaster=NULL;  //array of parameter names in master list of required parameters 
  class_type *aPCmaster=NULL; //array of parameter class types in master list of required parameters 
   
  if (!Options.silent){
    cout<<"Generating Master Parameter List..."<<endl;}

  // Check parameters for each estimation/correction algorithm:
  //--------------------------------------------------------------------------
  pModel->GetParticipatingParamList(aP,aPC,nP,Options);
  AddToMasterParamList(aPmaster,aPCmaster,nPmaster,aP,aPC,nP);

  // Check parameters for each hydrological process algorithm:
  //--------------------------------------------------------------------------
  for (int j=0;j<pModel->GetNumProcesses();j++)
  {
    pModel->GetProcess(j)->GetParticipatingParamList(aP,aPC,nP);
    AddToMasterParamList(aPmaster,aPCmaster,nPmaster,aP,aPC,nP);
  }
  // Add parameters needed to support calculation of maximum state variable value (in CHydroUnit::GetStateVarMax())
  //--------------------------------------------------------------------------
  if (pModel->GetStateVarIndex(SNOW_LIQ) != DOESNT_EXIST){
    aP [0]="SNOW_SWI"; aPC[0]=CLASS_GLOBAL;
    AddToMasterParamList(aPmaster, aPCmaster,nPmaster, aP, aPC, 1);
  }
  if (pModel->GetStateVarIndex(DEPRESSION) != DOESNT_EXIST){
    aP [0]="DEP_MAX"; aPC[0]=CLASS_LANDUSE;
    AddToMasterParamList(aPmaster, aPCmaster,nPmaster, aP, aPC, 1);
  }
  //canopy handled in CmvPrecipitation::GetParticipatingParamList

  //Prior to this, all parameters default to NOT_NEEDED or NOT_NEEDED_AUTO
  //this overrides default values of required variables to NOT_SPECIFIED and AUTO_COMPUTE
  double val;
  for (int p=0;p<nPmaster;p++)
  {
    if      (aPCmaster[p]==CLASS_SOIL      )
    {      
      val=CSoilClass::GetSoilProperty(parsed_soils[0],aPmaster[p]);
      if (val==NOT_NEEDED_AUTO){val=AUTO_COMPUTE;}
      else if (val==NOT_NEEDED){val=NOT_SPECIFIED;}
      CSoilClass::SetSoilProperty    (parsed_soils[0],aPmaster[p],val);
    }
    else if (aPCmaster[p]==CLASS_VEGETATION)
    {
      val=CVegetationClass::GetVegetationProperty(parsed_veg[0],aPmaster[p]);
      if (val==NOT_NEEDED_AUTO){val=AUTO_COMPUTE;}
      else if (val==NOT_NEEDED){val=NOT_SPECIFIED;}
      CVegetationClass::SetVegetationProperty    (parsed_veg[0],aPmaster[p],val);
    }
    else if (aPCmaster[p]==CLASS_LANDUSE   )
    {   
      val=CLandUseClass::GetSurfaceProperty(parsed_surf[0],aPmaster[p]);
      if (val==NOT_NEEDED_AUTO){val=AUTO_COMPUTE;}
      else if (val==NOT_NEEDED){val=NOT_SPECIFIED;}
      CLandUseClass::SetSurfaceProperty   (parsed_surf [0],aPmaster[p],val);
    }
    else if (aPCmaster[p]==CLASS_TERRAIN   )
    {   
      val=CTerrainClass::GetTerrainProperty(parsed_terrs[0],aPmaster[p]);
      if (val==NOT_NEEDED_AUTO){val=AUTO_COMPUTE;}
      else if (val==NOT_NEEDED){val=NOT_SPECIFIED;}
      CTerrainClass::SetTerrainProperty    (parsed_terrs[0],aPmaster[p],val);
    }
    else if (aPCmaster[p]==CLASS_GLOBAL    )
    {   
      val=CGlobalParams::GetGlobalProperty(global_template,aPmaster[p]);
      if (val==NOT_NEEDED_AUTO){val=AUTO_COMPUTE;}
      else if (val==NOT_NEEDED){val=NOT_SPECIFIED;}
      CGlobalParams::SetGlobalProperty    (global_template,aPmaster[p],val);
    }
  }

  if (Options.noisy){
    cout <<"======================================================"<<endl;
    cout << "Parsing Properties File " << Options.rvp_filename <<"..."<<endl;
    cout <<"======================================================"<<endl;
  }

  int   code;
  bool  done;
  bool  ended(false);
  int   Len,line(0);
  char *s[MAXINPUTITEMS];

  ifstream INPUT;
  INPUT.open(Options.rvp_filename.c_str());  
  if (INPUT.fail()){ 
    cout << "ERROR opening file: "<< Options.rvp_filename<<endl; return false;}
  
  ExitGracefullyIf(pModel==NULL,
    "ParseClassPropertiesFile: model has not yet been created.",BAD_DATA);

  CParser *p=new CParser(INPUT,Options.rvp_filename,line);

  ifstream INPUT2;           //For Secondary input
  CParser *pMainParser=NULL; //for storage of main parser while reading secondary files

  //--Sift through file-----------------------------------------------
  bool end_of_file=p->Tokenize(s,Len);
  while (!end_of_file) 
  {
    if (ended){break;}
    if (Options.noisy){ cout << "reading line " << p->GetLineNumber() << ": ";}
    /*assign code for switch statement
    ------------------------------------------------------------------
    <100         : ignored/special
    0   thru 100 : Options
    100 thru 200 : Soil Classes
    200 thru 300 : Land Use Classes
    300 thru 400 : Vegetation classes
    400 thru 500 : 
    500 thru 600 : Snow Classes
    ------------------------------------------------------------------
    */

    code=0;   
    //---------------------SPECIAL -----------------------------
    if       (Len==0)                                  {code=-1; }
    else if  (!strcmp(s[0],"*"                       )){code=-2; }//comment
    else if  (!strcmp(s[0],"%"                       )){code=-2; }//comment
    else if  (!strcmp(s[0],"#"                       )){code=-2; }//comment
    else if  (s[0][0]=='#')                            {code=-2; }//comment
    else if  (!strcmp(s[0],":RedirectToFile"         )){code=-3; }//redirect to secondary file
    else if  (!strcmp(s[0],":End"                    )){code=-4; }//stop reading
    //--------------------SOIL PARAMS --------------------------
    else if  (!strcmp(s[0],":SoilClasses"            )){code=1;  }//REQUIRED
    else if  (!strcmp(s[0],":SoilProperties"         )){code=2;  }
    else if  (!strcmp(s[0],":SoilProfiles"           )){code=6;  }//REQUIRED
    else if  (!strcmp(s[0],":SoilParameterList"      )){code=57; }
    //--------------------LU/LT PARAMS -------------------------
    else if  (!strcmp(s[0],":LandUseClasses"         )){code=100;}//REQUIRED
    else if  (!strcmp(s[0],":LandUseParameterList"   )){code=101;}
    else if  (!strcmp(s[0],":SnowParameterList"      )){code=101;}// \todo TO BECOME OBSOLETE
    else if  (!strcmp(s[0],":GlacierParameterList"   )){code=101;}// \todo TO BECOME OBSOLETE
    else if  (!strcmp(s[0],":LandUseChange"          )){code=102;}
    //--------------------VEGETATION PARAMS --------------------
    else if  (!strcmp(s[0],":VegetationClasses"      )){code=200;}//REQUIRED
    else if  (!strcmp(s[0],":SeasonalCanopyLAI"      )){code=201;}
    else if  (!strcmp(s[0],":SeasonalRelativeLAI"    )){code=201;}
    else if  (!strcmp(s[0],":SeasonalCanopyHeight"   )){code=202;}
    else if  (!strcmp(s[0],":SeasonalRelativeHeight" )){code=202;}
    else if  (!strcmp(s[0],":VegetationParameterList")){code=206;}
    else if  (!strcmp(s[0],":VegetationChange"       )){code=207;}
    //--------------------AQUIFER PARAMS -----------------------'
    else if  (!strcmp(s[0],":AquiferClasses"         )){code=300;}
    else if  (!strcmp(s[0],":AquiferProfiles"        )){code=301;}
    else if  (!strcmp(s[0],":AquiferStacks"          )){code=301;}
    //--------------------RIVER CHANNEL PARAMS -------------------
    else if  (!strcmp(s[0],":ChannelProfile"         )){code=500;}
    else if  (!strcmp(s[0],":ChannelRatingCurves"    )){code=501;}
    //--------------------TERRAIN PARAMS -----------------------
    else if  (!strcmp(s[0],":TerrainClasses"         )){code=600;}//REQUIRED
    else if  (!strcmp(s[0],":TerrainParameterList"   )){code=601;}
    //--------------------GLOBAL PARAMS ------------------------
    else if  (!strcmp(s[0],":AdiabaticLapseRate"     )){code=700;}
    else if  (!strcmp(s[0],":UBCTempLapseRates"      )){code=701;}
    else if  (!strcmp(s[0],":UBCPrecipLapseRates"    )){code=702;}
    else if  (!strcmp(s[0],":UBCEvapLapseRates"      )){code=703;}
    else if  (!strcmp(s[0],":WetAdiabaticLapse"      )){code=704;}   
    else if  (!strcmp(s[0],":ReferenceMaxTemperatureRange")){code=705;} 
    else if  (!strcmp(s[0],":UBCSnowParams"          )){code=706;}
    else if  (!strcmp(s[0],":RainSnowTransition"     )){code=708;}
    else if  (!strcmp(s[0],":IrreducibleSnowSaturation")){code=709;}
    else if  (!strcmp(s[0],":UBCGroundwaterSplit"    )){code=710;}
    else if  (!strcmp(s[0],":PrecipitationLapseRate" )){code=711;}
    else if  (!strcmp(s[0],":UBCExposureFactor"      )){code=712;}
    else if  (!strcmp(s[0],":UBCCloudPenetration"    )){code=713;}
    else if  (!strcmp(s[0],":UBCLWForestFactor"      )){code=714;}
    else if  (!strcmp(s[0],":UBCFlashPonding"        )){code=715;}
    else if  (!strcmp(s[0],":AirSnowCoeff"           )){code=716;}
    else if  (!strcmp(s[0],":AvgAnnualSnow"          )){code=717;}
    else if  (!strcmp(s[0],":AvgAnnualRunoff"        )){code=718;}
    else if  (!strcmp(s[0],":UBCNorthSWCorr"         )){code=750;}  //TEMP
    else if  (!strcmp(s[0],":UBCSouthSWCorr"         )){code=751;}  //TEMP
    else if  (!strcmp(s[0],":GlobalParameter"        )){code=720;}
    //--------------------TRANSPORT-------------------------------
    else if  (!strcmp(s[0],":DecayCoefficient"       )){code=900;}
    else if  (!strcmp(s[0],":SoilTransportParamList" )){code=901;}


    //--------------------OTHER ------- ------------------------
    else if  (!strcmp(s[0],":TransientParameter"     )){code=800;}

    switch(code)
    {
      case(-1):  //----------------------------------------------
      {/*Blank Line*/
        if (Options.noisy) {cout <<""<<endl;}break;
      }
      case(-2):  //----------------------------------------------
      {/*Comment*/
        if (Options.noisy) {cout <<"*"<<endl;} break;
      }
      case(-3):  //----------------------------------------------
      {/*:RedirectToFile*/
        string filename="";
        for (int i=1;i<Len;i++){filename+=s[i]; if(i<Len-1){filename+=' ';}}
        if (Options.noisy) {cout <<"Redirect to file: "<<filename<<endl;}
        
        string filedir = GetDirectoryName(Options.rvp_filename); //if a relative path name, e.g., "/path/model.rvp", only returns e.g., "/path"
        if (StringToUppercase(filename).find(StringToUppercase(filedir)) == string::npos){ //checks to see if absolute dir already included in redirect filename
          filename = filedir + "//" + filename;
        }

        INPUT2.open(filename.c_str()); 
        if (INPUT2.fail()){
          ostrstream FILENAME;
          FILENAME<<":RedirectToFile: Cannot find file "<<filename<<ends;
          ExitGracefully(FILENAME.str() ,BAD_DATA); 
        }
        else{
          pMainParser=p;    //save pointer to primary parser
          p=new CParser(INPUT2,filename,line);//open new parser
        }
        break;
      }
      case(-4):  //----------------------------------------------
      {/*:End*/
        if (Options.noisy) {cout <<"EOF"<<endl;} ended=true; break;
      }
      case(1):  //----------------------------------------------
      {/*:SoilClasses
        string "SoilClasses" 
        {string tag, double pct_sand,double pct_clay,double pct_silt, double pct_org}xNumSoilClasses
        :EndSoilClases*/
        if (Options.noisy) {cout <<"Soil Classes"<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        done=false;
        while (!done)
        { 
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Attributes")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
		      else if (Len >= 1)
		      {
			      if (num_parsed_soils >= MAX_SOIL_CLASSES - 1){
				      ExitGracefully("ParseClassPropertiesFile: exceeded maximum # of soil classes", BAD_DATA);
			      }

			      CSoilClass::InitializeSoilProperties(parsed_soils[num_parsed_soils], false);
			      pSoilClasses[num_parsed_soils - 1] = new CSoilClass(s[0]);
			      soiltags[num_parsed_soils] = s[0];

			      if (Len >= 5){
				      parsed_soils[num_parsed_soils].sand_con = s_to_d(s[1]);
				      parsed_soils[num_parsed_soils].clay_con = s_to_d(s[2]);
				      double silt_con = s_to_d(s[3]);
				      if (fabs(silt_con + s_to_d(s[1]) + s_to_d(s[2]) - 1.0) > REAL_SMALL){
					      WriteWarning("sand, clay, and silt content should add up to one for soil class " + to_string(s[0]) + ". Silt content will be recalculated.", Options.noisy);
				      }
				      parsed_soils[num_parsed_soils].org_con = s_to_d(s[4]); 
			      }
			      else{//default - common for conceptual models
              parsed_soils[num_parsed_soils].sand_con = 0.4111111; //special code for autogeneration of soil params
				      parsed_soils[num_parsed_soils].clay_con = 0.0;
				      parsed_soils[num_parsed_soils].org_con  = 0.0;
			      }
			    num_parsed_soils++;
          }
          else{
            WriteWarning("Incorrect number of terms in :SoilClasses entry",Options.noisy );
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndSoilClasses")){done=true;}
        }
        break;
      }
      case(2):  //----------------------------------------------
      {/*SoilProperties
        :SoilProperties 
        {string soil_tag, double porosity,double stone_frac,double rho_bulk}x<=NumSoilClasses
        :EndSoilProperties*/
        if (Options.noisy) {cout <<"Soil Properties"<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        invalid_index=ParsePropArray(p,indices,properties,num_read,soiltags,4,num_parsed_soils);
        ExitGracefullyIf(invalid_index,
              "ParseClassPropertiesFile: Invalid soiltype code in SoilProperties command",BAD_DATA);
        for (int i=0;i<num_read;i++)
        {
          parsed_soils[indices[i]].porosity    =properties[i][0]; 
          parsed_soils[indices[i]].stone_frac  =properties[i][1];
          parsed_soils[indices[i]].bulk_density=properties[i][2];
        } 
        break;
      }
      case(6):  //----------------------------------------------
      {/*SoilProfiles
        ":SoilProfiles" 
        {string tag, int numhorizons, soil1,thick1,soil2,thick2,...,soilN,thickN}xNumHorizons
        :EndSoilProfiles
        [[thicknesses in meters]]*/
        if (Options.noisy) {cout <<"Soil Profiles"<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        done=false;
        while (!done)
        { 
          if      (IsComment(s[0], Len)){}//comment line
          else if (Len>=2)
          {
            if (num_parsed_profiles>=MAX_SOIL_PROFILES-1){
              ExitGracefully("ParseClassPropertiesFile: exceeded maximum # of soil profiles",BAD_DATA);}
            
            pProfiles[num_parsed_profiles]=new CSoilProfile(s[0]);

            int nhoriz=s_to_i(s[1]);
            ExitGracefullyIf(nhoriz<0,
              "ParseClassPropertiesFile: invalid number of soil horizons in profile",BAD_DATA);
            //zero is a valid entry for lakes & glaciers (zero horizons)
            ExitGracefullyIf(Len!=(nhoriz*2+2),
              "ParseClassPropertiesFile:  :SoilProfiles invalid command length",BAD_DATA);   
            bool is_special=((!string(s[0]).compare("LAKE")) ||
                             (!string(s[0]).compare("GLACIER")) ||
                             (!string(s[0]).compare("ROCK")));
            ExitGracefullyIf((nhoriz==0) && (!is_special),
              "ParseClassPropertiesFile:  only special soil profiles (LAKE,GLACIER, or ROCK) can have zero horizons",BAD_DATA);   
              
            for (int m=0;m<nhoriz;m++)
            {
              const CSoilClass *pHorizon;
              double thisthick;
              thisthick=s_to_d(s[2*m+2+1]);
              pHorizon =CSoilClass::StringToSoilClass(string(s[2*m+2]));
              if (pHorizon==NULL){cout<<"Offending soilcode: "<<string(s[2*m+2])<<endl;}
              ExitGracefullyIf(pHorizon==NULL,
                "ParseClassPropertiesFile: bad soiltype code in soil profile",BAD_DATA);
              pProfiles[num_parsed_profiles]->AddHorizon(thisthick,pHorizon);
            }
            num_parsed_profiles++;
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndSoilProfiles")){done=true;}
        }

        break;
      }
      case(57)://----------------------------------------------
      {/*:SoilParameterList (General) [optional comments]
        string ":SoilParameterList" 
        :Parameters, NAME_1,NAME_2,...,NAME_N
        :Units,      unit_1,unit_2,...,unit_N
        {string soil_tag, double param_1, param_2, ... ,param_3}x<=NumSoilClasses
        :EndSoilParameterList*/
        if (Options.noisy) {cout <<"Soil Parameter List"<<endl;}
        done=false;
        nParamStrings=0;
        while (!done)
        {
          p->Tokenize(s,Len);
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Parameters")){
            for (int i=0;i<Len;i++){
              aParamStrings[i]=s[i];
            }
            nParamStrings=Len;
            //done=true;
          }
          else if (!strcmp(s[0],":Units")){
            //Do nothing with units for now
            done=true;
          }
          else {p->ImproperFormat(s); break;} 
        }
        invalid_index=ParsePropArray(p,indices,properties,num_read,soiltags,nParamStrings,num_parsed_soils);
        ExitGracefullyIf(invalid_index,
              "ParseClassPropertiesFile: Invalid soiltype code in SoilParameterList command",BAD_DATA);
        if (Options.noisy){
          for (int j=0;j<nParamStrings-1;j++){cout<<"  "<<aParamStrings[j+1]<<endl;}
        }
        for (int i=0;i<num_read;i++)
        {
          for (int j=0;j<nParamStrings-1;j++)
          {
            CSoilClass::SetSoilProperty(parsed_soils[indices[i]],aParamStrings[j+1],properties[i][j]); 
          }
        }
        break;
      }
      //===========================================================================================
      //===========================================================================================
      case(100):  //----------------------------------------------
      {/*:LandUseClasses
        {string LULT_tag, double frac_impermeable {opt forest_coverage}}
        :EndLandUseClasses*/
        if (Options.noisy) {cout <<"Land Use Classes"<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        done=false;
        while (!done)
        { 
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Attributes")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (Len>=2)
          {
            if (num_parsed_lult>=MAX_LULT_CLASSES-1){
              ExitGracefully("ParseClassPropertiesFile: exceeded maximum # of LU/LT classes",BAD_DATA);}
           
            CLandUseClass::InitializeSurfaceProperties(parsed_surf[num_parsed_lult],false);
            pLUClasses [num_parsed_lult-1]=new CLandUseClass(s[0]);
            lulttags   [num_parsed_lult]=s[0];
            parsed_surf[num_parsed_lult].impermeable_frac=s_to_d(s[1]); 
            if (Len>=3){
              parsed_surf[num_parsed_lult].forest_coverage=s_to_d(s[2]);
            }
            num_parsed_lult++;
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndLandUseClasses")){done=true;}
        }
        break;
      }
      case(101)://----------------------------------------------
      {/*:LandUseParameterList (General) [optional comments]
        string ":LandUseParameterList" 
        :Parameters, NAME_1,NAME_2,...,NAME_N
        {string lult_tag, double param_1, param_2, ... ,param_3}x<=NumLULTClasses
        :LandUseParameterList*/
        if (Options.noisy) {cout <<"Land Use / Land Type Parameter List"<<endl;}
        done=false;
        while (!done)
        {
          p->Tokenize(s,Len);
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Parameters")){
            for (int i=0;i<Len;i++){
              aParamStrings[i]=s[i];
            }
            nParamStrings=Len;
            //done=true;
          }
          else if (!strcmp(s[0],":Units")){
            //Do nothing with units for now
            done=true;
          }
          else {p->ImproperFormat(s); break;} 
        }
        invalid_index=ParsePropArray(p,indices,properties,num_read,lulttags,nParamStrings,num_parsed_lult);
        ExitGracefullyIf(invalid_index,
              "ParseClassPropertiesFile: Invalid land use code in LandUseParameterList command",BAD_DATA);
        if (Options.noisy){
          for (int j=0;j<nParamStrings-1;j++){cout<<"  "<<aParamStrings[j+1]<<endl;}
        }
        for (int i=0;i<num_read;i++)
        {
          for (int j=0;j<nParamStrings-1;j++){
            CLandUseClass::SetSurfaceProperty(parsed_surf[indices[i]],aParamStrings[j+1],properties[i][j]);
          }
        }
        break;
      }
      case(102):  //----------------------------------------------
      {/*:LandUseChange [HRU group] [new LULT tag] [YYYY-mm-dd] */
        if (Options.noisy) {cout <<"Change in Land Use Class"<<endl;}
        if (Len<4){p->ImproperFormat(s); break;} 
        time_struct tt;
        tt=DateStringToTimeStruct(string(s[3]),string("00:00:00"));
        pModel->AddPropertyClassChange(s[1],CLASS_LANDUSE,s[2], tt);
        break;
      }
      //==========================================================
      //==========================================================
      case(200):  //----------------------------------------------
      {/*Vegetation Classes
        :VegetationClasses
        {string tag, double max height [m], max LAI, max leaf cond}xNumVegetationClasses
        :EndVegetationClasses*/
        if (Options.noisy) {cout <<"Vegetation Classes"<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        done=false;
        while (!done)
        { 
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Attributes")){}//attributes are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (Len==4)
          {
            if (num_parsed_veg>=MAX_VEG_CLASSES-1){
              ExitGracefully("ParseClassPropertiesFile: exceeded maximum # of vegetation classes",BAD_DATA);}
           
            CVegetationClass::InitializeVegetationProps(parsed_veg [num_parsed_veg],false);
            pVegClasses[num_parsed_veg-1]=new CVegetationClass(s[0]);
            vegtags    [num_parsed_veg]=s[0];
            parsed_veg [num_parsed_veg].max_height   =s_to_d(s[1]); 
            parsed_veg [num_parsed_veg].max_LAI      =s_to_d(s[2]); 
            parsed_veg [num_parsed_veg].max_leaf_cond=s_to_d(s[3]);
            
            num_parsed_veg++;
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndVegetationClasses")){done=true;}
        }
        break;
      }
      case(201):  //----------------------------------------------
      {/*SeasonalCanopyLAI
        :SeasonalCanopyLAI 
        {string veg_tag, double LAI_Jan..LAI_dec}x<=NumVegClasses
        :EndSeasonalCanopyLAI*/
        if (Options.noisy) {cout <<"Seasonal Canopy LAI"<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        invalid_index=ParsePropArray(p,indices,properties,num_read,vegtags,13,num_parsed_veg);
        ExitGracefullyIf(invalid_index,
              "ParseClassPropertiesFile: Invalid vegetation code in SeasonalCanopyLAI command",BAD_DATA);
        for (int i=0;i<num_read;i++){
          for (int mon=0;mon<12;mon++){
            parsed_veg[indices[i]].relative_LAI[mon]=properties[i][mon]; 
          }
        }
        break;
      }
      case(202):  //----------------------------------------------
      {/*SeasonalCanopyHeight
        :SeasonalCanopyHeight
        {string veg_tag, double Ht_Jan..Ht_dec}x<=NumVegClasses
        :EndSeasonalCanopyHeight*/
        if (Options.noisy) {cout <<"Seasonal Canopy Height"<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        invalid_index=ParsePropArray(p,indices,properties,num_read,vegtags,13,num_parsed_veg);
        ExitGracefullyIf(invalid_index,
              "ParseClassPropertiesFile: Invalid vegetation code in SeasonalCanopyHeight command",BAD_DATA);
        for (int i=0;i<num_read;i++){
          for (int mon=0;mon<12;mon++){
            parsed_veg[indices[i]].relative_ht[mon]=properties[i][mon]; 
          }
        }
        break;
      }
      case(203):  //----------------------------------------------
      {
        break;
      }
      case(204):  //----------------------------------------------
      {
        break;
      }
      case(205):  //----------------------------------------------
      {
        break;
      }
      case(206)://----------------------------------------------
      {/*:VegetationParameterList (General) [optional comments]
       ":VegetationParameterList" 
        :Parameters, NAME_1,NAME_2,...,NAME_N
        :Units     , unit_1,unit_2,...,unit_N
        {string veg_tag, double param_1, param_2, ... ,param_N}x[<=NumVegClasses]
        :EndVegetationParameterList*/
        if (Options.noisy) {cout <<"Vegetation Parameter List"<<endl;}
        done=false;
        while (!done)
        {
          p->Tokenize(s,Len);
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Parameters")){
            for (int i=0;i<Len;i++){
              aParamStrings[i]=s[i];
            }
            nParamStrings=Len;
          }
          else if (!strcmp(s[0],":Units")){
            //Do nothing with units for now
            done=true;
          }
          else {p->ImproperFormat(s); break;} 
        }
        invalid_index=ParsePropArray(p,indices,properties,num_read,vegtags,nParamStrings,num_parsed_veg);
        ExitGracefullyIf(invalid_index,
              "ParseClassPropertiesFile: Invalid vegetation code in VegetationParameterList command",BAD_DATA);
        if (Options.noisy){
          for (int j=0;j<nParamStrings-1;j++){cout<<"  "<<aParamStrings[j+1]<<endl;}
        }
        for (int i=0;i<num_read;i++)
        {
          for (int j=0;j<nParamStrings-1;j++)
          {
            CVegetationClass::SetVegetationProperty(parsed_veg   [indices[i]],
                                                    aParamStrings[j+1],
                                                    properties   [i][j]);
          }
        }
        break;
      }
      case(207):  //----------------------------------------------
      {/*:VegetationChange [HRU group] [new Veg tag] [YYYY-mm-dd] */
        if (Options.noisy) {cout <<"Change in Vegetation"<<endl;}
        if (Len<4){p->ImproperFormat(s); break;} 
        time_struct tt;
        tt=DateStringToTimeStruct(string(s[3]),string("00:00:00"));
        pModel->AddPropertyClassChange(s[1],CLASS_VEGETATION,s[2], tt);
        break;
      }
      //==========================================================
      //==========================================================
      case(300):  //----------------------------------------------
      {/*AquiferClasses
        :AquiferClasses 
        {string tag, soil_type, thickness}xNumAquiferClasses
        :EndAquiferClasses*/
        if (Options.noisy) {cout <<"Aquifer Classes (OBSOLETE)"<<endl;}
        break;
      }
      case(301):  //----------------------------------------------
      {/*AquiferProfiles / AquiferStacks
        ":AquiferProfiles" 
        {string tag, int numhorizons, soil1,thick1,soil2,thick2,...,soilN,thickN}xNumHorizons
        :EndAquiferProfiles
        [[thicknesses in meters]]*/
        if (Options.noisy) {cout <<"Aquifer Profiles"<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        done=false;
        while (!done)
        { 
          if      (IsComment(s[0], Len)){}//comment line
          else if (Len>=2)
          {
            if (num_parsed_profiles>=MAX_AQUIFER_STACKS-1){
              ExitGracefully("ParseClassPropertiesFile: exceeded maximum # of aquifer stacks",BAD_DATA);}
            
            pAqStacks[num_parsed_aqstacks]=new CAquiferStack(s[0]);

            int nlayers=s_to_i(s[1]);
            ExitGracefullyIf(nlayers<0,
              "ParseClassPropertiesFile: invalid number of aquifer layers in stack",BAD_DATA);
            //zero is a valid entry for lakes & glaciers (zero horizons)
            ExitGracefullyIf(Len!=(nlayers*2+2),
              "ParseClassPropertiesFile:  :AquiferProfiles invalid command length",BAD_DATA);       
            for (int m=0;m<nlayers;m++)
            {
              const CSoilClass *pLayerSoil;
              double thisthick;
              thisthick=s_to_d(s[2*m+2+1]);
              pLayerSoil =CSoilClass::StringToSoilClass(string(s[2*m+2]));
              if (pLayerSoil==NULL){cout<<"Offending soilcode: "<<string(s[2*m+2])<<endl;}
              ExitGracefullyIf(pLayerSoil==NULL,
                "ParseClassPropertiesFile: bad soiltype code in aquifer stack",BAD_DATA);
              pAqStacks[num_parsed_aqstacks]->AddLayer(thisthick,pLayerSoil);
            }
            num_parsed_aqstacks++;
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndAquiferProfiles")){done=true;}
        }

        break;
      }
      case(500):  //----------------------------------------------
      {/*ChannelProfile
        :ChannelProfile [optional string name]
          :Bedslope {double slope}
          :SurveyPoints
            {double x double bed_elev}x num survey points
          :EndSurveyPoints
          :RoughnessZones
            {double x_z double mannings_n}xnum roughness zones
          :EndRoughnessZones
        :EndChannelProfile*/
        string tag;
        double slope=0.0;
        if (Options.noisy) {cout <<"Channel Profile"<<endl;}
        tag=s[1]; 
        double *x =new double [MAX_SURVEY_PTS];
        double *y =new double [MAX_SURVEY_PTS];
        double *xz=new double [MAX_SURVEY_PTS];
        double *n =new double [MAX_SURVEY_PTS];
        double *nn=new double [MAX_SURVEY_PTS];

        int countSP=0; //number of survey points
        int countRS=0; //number of roughness segments

        p->Tokenize(s,Len);
        done=false;
        while (!done)
        { 
          if      (IsComment(s[0], Len)){}//comment line
          else if ((Len==2) && (!strcmp(s[0],":Bedslope")))   
          {
            slope=s_to_d(s[1]);
          }
          else if ((Len==1) && (!strcmp(s[0],":SurveyPoints")))
          {
            p->Tokenize(s,Len);
            done=false;
            while (!done)
            { 
              if      (IsComment(s[0], Len)){}//comment line
              else if (Len==2)           {x[countSP]=s_to_d(s[0]); y[countSP]=s_to_d(s[1]);countSP++;}
              else                       {p->ImproperFormat(s); break;}
              p->Tokenize(s,Len);
              if (!strcmp(s[0],":EndSurveyPoints")){done=true;}
            }
            done=false;
          }
          else if ((Len==1) && (!strcmp(s[0],":RoughnessZones")))
          {
            p->Tokenize(s,Len);
            done=false;
            while (!done)
            { 
              if      (IsComment(s[0], Len)){}//comment line
              else if (Len==2           ){xz[countRS]=s_to_d(s[0]);n [countRS]=s_to_d(s[1]);countRS++;}
              else                       {p->ImproperFormat(s); break;}
              p->Tokenize(s,Len);
              if (!strcmp(s[0],":EndRoughnessZones")){done=true;}
            }
            done=false;
          }
          else if ((Len==1) && (!strcmp(s[0],":EndChannelProfile")))
          {
            done=true;
          }
          else                       {p->ImproperFormat(s); break;}
          if (done!=true){p->Tokenize(s,Len);}
        }

        //Create Channel
        CChannelXSect *pChannel=NULL;

        //map roughnesses to survey segments
        if (xz[0]>x[0]){
          WriteWarning(":ChannelProfile command: leftmost mannings zone to right of leftmost survey point",Options.noisy);
		    }
       
		    int j=0;
        for (int i=0; i<countSP;i++) //go through segments (segment i between x[i] and x[i+1])
        {
          if      (j==countRS-1){nn[i]=n[j];} //in rightmost zone
          else if (xz[j+1]>x[i+1]){nn[i]=n[j];}//next zone switch not in this segment, use most recent
          else { //zone switch (potentially more than one) in this segment
            double Li=(x[i+1]-x[i]);
            nn[i]=(xz[j+1]-x[i])/Li*n[j];
            while ((j<countRS-1) && (xz[j+1]<x[i+1]))
            {
              j++;
              if (j==countRS-1){nn[i]+=(x[i+1]-xz[j])/Li*n[j];}
              else            {nn[i]+=(min(x[i+1],xz[j+1])-xz[j])/Li*n[j];}
            }
          }
          //cout<<"n["<<i<<"]:"<<nn[i]<<endl;
        }

        pChannel=new CChannelXSect(tag,countSP,x,y,nn,slope); 

        delete [] x; delete [] n; delete [] xz; delete [] y;delete []nn;
        break;
      }
      case(501):  //----------------------------------------------
      {/*ChannelRatingCurves
        :ChannelRatingCurves [optional string name]
          :Bedslope {double slope}
          :SurveyPoints
            {double x double bed_elev}x num survey points
          :EndSurveyPoints
          :RoughnessZones
            {double x_z double mannings_n}xnum roughness zones
          :EndRoughnessZones
        :EndChannelProfile*/
        string tag;
        double slope=0.0;
        if (Options.noisy) {cout <<"Channel Profile"<<endl;}
        tag=s[1]; 
        double *st=new double [MAX_SURVEY_PTS];
        double *A =new double [MAX_SURVEY_PTS];
        double *W =new double [MAX_SURVEY_PTS];
        double *Q =new double [MAX_SURVEY_PTS];
        double *P =new double [MAX_SURVEY_PTS];
        int count=0;
        p->Tokenize(s,Len);
        done=false;
        while (!done)
        { 
          if      (IsComment(s[0], Len)){}//comment line
          else if ((Len>=2) && (!strcmp(s[0],":Bedslope")))   
          {
            slope=s_to_d(s[1]);
          }
          else if ((Len>=1) && (!strcmp(s[0],":StageRelationships")))
          {
            p->Tokenize(s,Len);
            done=false;
            while (!done)
            { 
              if      (IsComment(s[0], Len)){}//comment line
              else if (Len>=4)  {
                st[count]=s_to_d(s[0]); 
                A[count]=s_to_d(s[1]);
                W[count]=s_to_d(s[2]);
                Q[count]=s_to_d(s[3]);
                if (Len>=5){
                  P[count]=s_to_d(s[4]);
                }
                else{
                  P[count]=2*st[count]+W[count];
                }

                count++;
              }
              else{p->ImproperFormat(s); break;}
              p->Tokenize(s,Len);
              if (!strcmp(s[0],":EndStageRelationships")){done=true;}
            }
            done=false;
          }
          else if ((Len>=1) && (!strcmp(s[0],":EndChannelRatingCurves")))
          {
            done=true;
          }
          else                       {p->ImproperFormat(s); break;}
          if (done!=true){p->Tokenize(s,Len);}
        }

        //Create Channel
        CChannelXSect *pChannel=NULL;
        
        pChannel=new CChannelXSect(tag,count,Q,st,W,A,P,slope); 

        delete [] st; delete [] A; delete [] W; delete [] Q;
        break;
      }
      //==========================================================
      //==========================================================
      case(600):  //----------------------------------------------
      {/*TerrainClasses
        :TerrainClasses
        {string tag, hillslope_len, drainage_dens[,lambda]}xNumTerrainClasses
        :EndTerrainClasses*/
        if (Options.noisy) {cout <<"Terrain Classes"<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        done=false;
        while (!done)
        { 
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Attributes")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (Len>=3)
          {
            if ((num_parsed_terrs>=(MAX_TERRAIN_CLASSES-1)) || (num_parsed_terrs<0)){
              ExitGracefully("ParseClassPropertiesFile: exceeded maximum # of terrain classes",BAD_DATA);}
           
            CTerrainClass::InitializeTerrainProperties(parsed_terrs[num_parsed_terrs],false);//sets all to autocompute
            pTerrClasses[num_parsed_terrs-1]=new CTerrainClass(s[0]);
            terraintags [num_parsed_terrs]=s[0];
            parsed_terrs[num_parsed_terrs].hillslope_length =s_to_d(s[1]); 
            parsed_terrs[num_parsed_terrs].drainage_density =s_to_d(s[2]); 
						if (Len==4)
							parsed_terrs[num_parsed_terrs].lambda =s_to_d(s[3]);
						num_parsed_terrs++;
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndTerrainClasses")){done=true;}
        }
        break;
      }
      case(601)://----------------------------------------------
      {/*:TerrainParameterList (General) [optional comments]
        string ":TerrainParameterList" 
        :Parameters, NAME_1,NAME_2,...,NAME_N
        {string terrain_tag, double param_1, param_2, ... ,param_3}x<=NumTerrainClasses
        :EndTerrainParameterList*/
        if (Options.noisy) {cout <<"Terrain Parameter List"<<endl;}
        done=false;
        while (!done)
        {
          p->Tokenize(s,Len);
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Parameters")){
            for (int i=0;i<Len;i++){
              aParamStrings[i]=s[i];
            }
            nParamStrings=Len;
          }
          else if (!strcmp(s[0],":Units")){
            //Do nothing with units for now
            done=true;
          }
          else {p->ImproperFormat(s); break;} 
        }
        invalid_index=ParsePropArray(p,indices,properties,num_read,terraintags,nParamStrings,num_parsed_terrs);
        ExitGracefullyIf(invalid_index,
              "ParseClassPropertiesFile: Invalid terrain class code in TerrainParameterList command",BAD_DATA);
        if (Options.noisy){
          for (int j=0;j<nParamStrings-1;j++){cout<<"  "<<aParamStrings[j+1]<<endl;}
        }
        for (int i=0;i<num_read;i++)
        {
          for (int j=0;j<nParamStrings-1;j++){
            CTerrainClass::SetTerrainProperty(parsed_terrs[indices[i]],aParamStrings[j+1],properties[i][j]);
          }
        }
        break;
      }
      case(700):  //----------------------------------------------
      {/*:AdiabaticLapseRate
         string ":AdiabaticLapseRate" value*/
        if (Options.noisy){cout <<"Adiabatic Lapse Rate"<<endl;}
        if (Len<2){p->ImproperFormat(s); break;}
        parsed_globals.adiabatic_lapse=s_to_d(s[1]);
        break;
      }
      case(701):  //----------------------------------------------
      {/*:UBCTempLapseRates
         string ":UBCTempLapseRates" A0TLXM A0TLNM A0TLXH A0TLNH P0TEDL P0TEDU */
        if (Options.noisy){cout <<"UBC Temperature Lapse Rate Params"<<endl;}
        if (Len<7){p->ImproperFormat(s); break;}
        parsed_globals.UBC_lapse_params.A0TLXM=s_to_d(s[1]);
        parsed_globals.UBC_lapse_params.A0TLNM=s_to_d(s[2]);
        parsed_globals.UBC_lapse_params.A0TLXH=s_to_d(s[3]);
        parsed_globals.UBC_lapse_params.A0TLNH=s_to_d(s[4]);
        parsed_globals.UBC_lapse_params.P0TEDL=s_to_d(s[5]);
        parsed_globals.UBC_lapse_params.P0TEDU=s_to_d(s[6]);
        break;
      }
      case(702):  //----------------------------------------------
      {/*:UBCPrecipLapseRates
         string ":UBCPrecipLapseRates" E0LLOW E0LMID  E0LHI  P0GRADL P0GRADM P0GRADU A0STAB*/
        if (Options.noisy){cout <<"UBC Precipitation Lapse Rate Params"<<endl;}
        if (Len<8){p->ImproperFormat(s); break;}
        parsed_globals.UBC_lapse_params.E0LLOW =s_to_d(s[1]);
        parsed_globals.UBC_lapse_params.E0LMID =s_to_d(s[2]);
        parsed_globals.UBC_lapse_params.E0LHI  =s_to_d(s[3]);
        parsed_globals.UBC_lapse_params.P0GRADL=s_to_d(s[4]);
        parsed_globals.UBC_lapse_params.P0GRADM=s_to_d(s[5]);
        parsed_globals.UBC_lapse_params.P0GRADU=s_to_d(s[6]);
        parsed_globals.UBC_lapse_params.A0STAB =s_to_d(s[7]);
        break;
      }
      case(703):  //----------------------------------------------
      {/*:UBCEvapLapseRates
         string ":UBCEvapLapseRates" A0PELA */
        if (Options.noisy){cout <<"UBC Evaporation Lapse Rate Params"<<endl;}
        if (Len<2){p->ImproperFormat(s); break;}
        parsed_globals.UBC_lapse_params.A0PELA =s_to_d(s[1]);
        break;
      }
      case (704):  //----------------------------------------------
      {/*:WetAdiabaticLapse       
         string ":WetAdiabaticLapse" [wet adiabatic lapse rate]  [A0PPTP] */
        if (Options.noisy){cout <<"UBC Precipitation Wet Adiabatic Lapse Rate Params"<<endl;}
        if (Len<3){p->ImproperFormat(s); break;}
        parsed_globals.wet_adiabatic_lapse =s_to_d(s[1]);
        parsed_globals.UBC_lapse_params.A0PPTP =s_to_d(s[2]);
        break;
      }
      case (705):  //----------------------------------------------
      {/*:ReferenceMaxTemperatureRange    
         string ":ReferenceMaxTemperatureRange" [range] */
        if (Options.noisy){cout <<"UBC Precipitation Wet Adiabatic Lapse Rate Params"<<endl;}
        if (Len<2){p->ImproperFormat(s); break;}
        parsed_globals.UBC_lapse_params.max_range_temp =s_to_d(s[1]);
        break;
      }
      case (706):  //----------------------------------------------
      {/*:UBCSnowParams    
         string ":UBCSnowParams" P0ALBMIN P0ALBMAX P0ALBREC P0ALBASE P0ALBSNW P0ALBMLX */
        if (Options.noisy){cout <<"UBC Snow Parameters"<<endl;}
        if (Len<7){p->ImproperFormat(s); break;}
        parsed_globals.UBC_snow_params.MIN_SNOW_ALBEDO=s_to_d(s[1]);
        parsed_globals.UBC_snow_params.MAX_SNOW_ALBEDO=s_to_d(s[2]);
        parsed_globals.UBC_snow_params.ALBREC         =s_to_d(s[3]);
        parsed_globals.UBC_snow_params.ALBASE         =s_to_d(s[4]);
        parsed_globals.UBC_snow_params.ALBSNW         =s_to_d(s[5]);
        parsed_globals.UBC_snow_params.MAX_CUM_MELT   =s_to_d(s[6]);
        break;
      }
      case (707):  //----------------------------------------------
      {
         //Obsolete
        break;
      }
      case (708):  //----------------------------------------------
      {/*:RainSnowTransition    
         string ":RainSnowTransition"  rainsnow_temp rainsnow_delta*/
        if (Options.noisy){cout <<"Rain/Snow Transition Parameters"<<endl;}
        if (Len<3){p->ImproperFormat(s); break;}
        parsed_globals.rainsnow_temp =s_to_d(s[1]);
        parsed_globals.rainsnow_delta=s_to_d(s[2]);
        break;
      }
      case (709):  //----------------------------------------------
      {/*:IrreducibleSnowSaturation    
         string ":IrreducibleSnowSaturation"  snow_SWI*/
        if (Options.noisy){cout <<"Maximum liquid water content in snow"<<endl;}
        if (Len<2){p->ImproperFormat(s); break;}
        parsed_globals.snow_SWI =s_to_d(s[1]);
        break;
      }    
      case (710):  //----------------------------------------------
      {/*:UBCGroundwaterSplit    
         string ":UBCGroundwaterSplit"  split*/
        if (Options.noisy){cout <<"UBC Groundwater split fraction (PODZSH)"<<endl;}
        if (Len<2){p->ImproperFormat(s); break;}
        parsed_globals.UBC_GW_split =s_to_d(s[1]);
        break;
      }
      case(711):  //----------------------------------------------
      {/*:PrecipitationLapseRate
         string ":PrecipitationLapseRate" value*/
        if (Options.noisy){cout <<"Precipitation Lapse Rate"<<endl;}
        if (Len<2){p->ImproperFormat(s); break;}
        parsed_globals.precip_lapse=s_to_d(s[1]);
        break;
      }
      case (712):  //----------------------------------------------
      {/*:UBCExposureFactor    
         string ":UBCExposureFactor"  value*/
        if (Options.noisy){cout <<"UBC exposure factor (F0ERGY)"<<endl;}
        if (Len<2){p->ImproperFormat(s); break;}
        parsed_globals.UBC_exposure_fact =s_to_d(s[1]);
        break;
      }
      case (713):  //----------------------------------------------
      {/*:UBCCloudPenetration    
         string ":UBCCloudPenetration"  value*/
        if (Options.noisy){cout <<"UBC cloud penetration factor (P0CAST)"<<endl;}
        if (Len<2){p->ImproperFormat(s); break;}
        parsed_globals.UBC_cloud_penet =s_to_d(s[1]);
        break;
      }
      case (714):  //----------------------------------------------
      {/*:UBCLWForestFactor    
         string ":UBCLWForestFactor"  value*/
        if (Options.noisy){cout <<"UBC longwave forest estimation factor (P0BLUE*P0LWVF)"<<endl;}
        if (Len<2){p->ImproperFormat(s); break;}
        parsed_globals.UBC_LW_forest_fact =s_to_d(s[1]);
        break;
      }
      case (715):  //----------------------------------------------
      {/*:UBCFlashPonding   
         string ":UBCFlashPonding"  value*/
        if (Options.noisy){cout <<"UBC Flash Factor ponding threshold (V0FLAS)"<<endl;}
        if (Len<2){p->ImproperFormat(s); break;}
        parsed_globals.UBC_flash_ponding =s_to_d(s[1]);
        break;
      }
      case (716):  //----------------------------------------------
      {/*:AirSnowCoeff   
         string ":AirSnowCoeff"  value*/
        if (Options.noisy){cout <<"air/snow heat transfer coefficient"<<endl;}
        if (Len<2){p->ImproperFormat(s); break;}
        parsed_globals.airsnow_coeff=s_to_d(s[1]);
        break;
      }
      case (717):  //----------------------------------------------
      {/*:AvgAnnualSnow   
         string ":AvgAnnualSnow"  value*/
        if (Options.noisy){cout <<"average annual snowfall"<<endl;}
        if (Len<2){p->ImproperFormat(s); break;}
        parsed_globals.avg_annual_snow=s_to_d(s[1]);
        break;
      }
      case (718):  //----------------------------------------------
      {/*:AvgAnnualRunoff
         string ":AvgAnnualRunoff"  value*/
        if (Options.noisy){cout <<"average annual runoff"<<endl;}
        if (Len<2){p->ImproperFormat(s); break;}
        parsed_globals.avg_annual_runoff=s_to_d(s[1]);
        break;
      }
      case(750):  //----------------------------------------------
      {/*:UBCNorthSWCorr
         string ":UBCNorthSWCorr" J F M A M J J A S O N D*/
        if (Options.noisy){cout <<"UBC North-facing shortwave correction"<<endl;}
        if (Len<13){p->ImproperFormat(s); break;}
        for (int i=0;i<12;i++){
          parsed_globals.UBC_n_corr[i]=s_to_d(s[i+1]);
        }
        break;
      }
      case(751):  //----------------------------------------------
      {/*:UBCSouthSWCorr
         string ":UBCSouthSWCorr" J F M A M J J A S O N D*/
        if (Options.noisy){cout <<"UBC South-facing shortwave correction"<<endl;}
        if (Len<13){p->ImproperFormat(s); break;}
        for (int i=0;i<12;i++){
          parsed_globals.UBC_s_corr[i]=s_to_d(s[i+1]);
        }
        break;
      }
      case(720):  //----------------------------------------------
      {/*:GlobalParameter
         string ":GlobalParameter" {string PARAMETER_NAME} {double value}*/
        if (Options.noisy){cout <<"Global Parameter"<<endl;}
        if (Len<3){p->ImproperFormat(s); break;}
        CGlobalParams::SetGlobalProperty(parsed_globals,s[1],s_to_d(s[2]));
        break;
      }
      case(800):  //----------------------------------------------
      {/*:TransientParameter {PARAM_NAME} {Parameter_class} [ClassName]
           {yyyy-mm-dd hh:mm:ss double tstep int nMeasurements}
           {double value} x nMeasurements
         :EndTransientParameter*/
        if (Options.noisy){cout <<"Transient Parameter"<<endl;}
        if (Len<3){p->ImproperFormat(s); break;}
        string param_name,class_name;
        class_type  ptype(CLASS_GLOBAL);
        param_name=(string)(s[1]);
        if      (!strcmp(s[2],"SOIL"        )){ptype=CLASS_SOIL;}
        else if (!strcmp(s[2],"VEGETATION"  )){ptype=CLASS_VEGETATION;}
        else if (!strcmp(s[2],"LANDUSE"     )){ptype=CLASS_LANDUSE;}
        else if (!strcmp(s[2],"TERRAIN"     )){ptype=CLASS_TERRAIN;}
        else if (!strcmp(s[2],"GLOBALS"     )){ptype=CLASS_GLOBAL;}
        else{
          ExitGracefully("ParsePropertyFile: invalid parameter class in :TransientParameter command",BAD_DATA);
        }
        if ((ptype!=CLASS_GLOBAL) && (Len>=4)){class_name=(string)(s[3]);}

        CTimeSeries     *pTimeSer=NULL;
        pTimeSer=CTimeSeries::Parse(p,true,param_name+"_"+s[2]+"_"+s[3],"",Options);
        if (pTimeSer!=NULL)
        {
          CTransientParam *pTransParam=NULL;
          pTransParam = new CTransientParam(pTimeSer,param_name,ptype,class_name);
          pModel->AddTransientParameter(pTransParam); 
        }
        else{
          WriteWarning("ParseClassPropertiesFile: unable to read :TransientParameter time series",true);
        }
        break;
      }
      case(900):  //----------------------------------------------
      {/*:DecayCoefficient
         :DecayCoefficient [constit_name] [value]*/
        if (Options.noisy){cout <<"DecayCoefficient"<<endl;}
        if (Len<3){p->ImproperFormat(s); break;}
        if (pModel->GetTransportModel()!=NULL){
        //pModel->GetTransportModel()->SetParameter(s[1],"DECAY_COEFF",s_to_d(s[2]),0)
        }
        break;
      }
      case(901)://----------------------------------------------
      {/*:SoilTransportParameterList (General) [optional comments]
        string ":SoilTransportParameterList" [constit_name]
        :Parameters, NAME_1,NAME_2,...,NAME_N
        :Units,      unit_1,unit_2,...,unit_N
        {string soil_tag, double param_1, param_2, ... ,param_3}x<=NumSoilClasses
        :EndSoilTransportParameterList*/
        if (Options.noisy) {cout <<"Soil Transport Parameter List"<<endl;}
        done=false;
        if (Len<2){p->ImproperFormat(s); break;}
        int constit_ind=pModel->GetTransportModel()->GetConstituentIndex(s[1]);
        if (constit_ind==DOESNT_EXIST)
        {
          ExitGracefully("Parsing property file: invalid constituent name in :SoilTransportParameterList command",BAD_DATA_WARN);
          break;
        }
        nParamStrings=0;
        while (!done)
        {
          p->Tokenize(s,Len);
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Parameters")){
            for (int i=0;i<Len;i++){
              aParamStrings[i]=s[i];
            }
            nParamStrings=Len;
            //done=true;
          }
          else if (!strcmp(s[0],":Units")){
            //Do nothing with units for now
            done=true;
          }
          else {p->ImproperFormat(s); break;} 
        }
        invalid_index=ParsePropArray(p,indices,properties,num_read,soiltags,nParamStrings,num_parsed_soils);
        ExitGracefullyIf(invalid_index,
              "ParseClassPropertiesFile: Invalid soiltype code in SoilParameterList command",BAD_DATA);
        if (Options.noisy){
          for (int j=0;j<nParamStrings-1;j++){cout<<"  "<<aParamStrings[j+1]<<endl;}
        }
        for (int i=0;i<num_read;i++)
        {
          for (int j=0;j<nParamStrings-1;j++)
          {
            CSoilClass::SetSoilTransportProperty(constit_ind,parsed_soils[indices[i]],aParamStrings[j+1],properties[i][j]); 
          }
        }
        break;
      }
      default://----------------------------------------------
      {
        char firstChar = *(s[0]);
        switch(firstChar)
        {
        case ':':
          {
            if     (!strcmp(s[0],":FileType"))    {if (Options.noisy){cout<<"Filetype"<<endl;}}//do nothing
            else if(!strcmp(s[0],":Application")) {if (Options.noisy){cout<<"Application"<<endl;}}//do nothing
            else if(!strcmp(s[0],":Version"))     {if (Options.noisy){cout<<"Version"<<endl;}}//do nothing
            else if(!strcmp(s[0],":WrittenBy"))   {if (Options.noisy){cout<<"WrittenBy"<<endl;}}//do nothing
            else if(!strcmp(s[0],":CreationDate")){if (Options.noisy){cout<<"CreationDate"<<endl;}}//do nothing
            else if(!strcmp(s[0],":SourceFile"))  {if (Options.noisy){cout<<"SourceFile"<<endl;}}//do nothing
            else if(Options.noisy)
            {
							string warn;
							warn="Ignoring unrecognized command starting with "+string(s[0])+" in .rvp file";
							WriteWarning(warn,Options.noisy);
						}
          }
          break;
        default:
          {
			  string errString = "Unrecognized command in .rvp file:\n   " + string(s[0]);
            ExitGracefully(errString.c_str(),BAD_DATA);//STRICT 
          }
          break;
        }
      }
    }//end switch
    end_of_file=p->Tokenize(s,Len);

    //return after file redirect
    if ((end_of_file) && (pMainParser!=NULL))
    {
      if (Options.noisy) {cout <<"...Return to main file"<<endl;}
      INPUT2.clear();INPUT2.close();
      delete p;
      p=pMainParser;
      pMainParser=NULL;
      end_of_file=p->Tokenize(s,Len);
    }
  } //end while (!end_of_file)

  INPUT.close();
  
   
  if (!Options.silent){cout<<"Autocalculating Model Parameters..."<<endl;}

  CGlobalParams::AutoCalculateGlobalParams          (parsed_globals,global_template);
 
  for (int c=1;c<num_parsed_veg;c++){
    pVegClasses [c-1]->AutoCalculateVegetationProps (parsed_veg[c],parsed_veg[0]); 
  }
  for (int c=1;c<num_parsed_soils;c++){
    pSoilClasses[c-1]->AutoCalculateSoilProps       (parsed_soils[c],parsed_soils[0]); 
  }
  for (int c=1;c<num_parsed_lult;c++){
    pLUClasses  [c-1]->AutoCalculateLandUseProps    (parsed_surf[c],parsed_surf[0]); 
  }
  for (int c=1;c<num_parsed_terrs;c++){
    pTerrClasses[c-1]->AutoCalculateTerrainProps    (parsed_terrs[c],parsed_terrs[0]); 
  }

  if (!Options.silent){
    cout<<"...done Autocalculating."<<endl;
    cout<<"Checking for Required Model Parameters..."<<endl;
}

  RVPParameterWarning(aPmaster,aPCmaster,nPmaster,Options);

  //Check for existence of classes
  //--------------------------------------------------------------------------
  int terrain_param_count(0);
  int lult_param_count(0);
  int veg_param_count(0);
  int soil_param_count(0);
  for (int ii=0;ii<nP;ii++)
  {
    if      (aPC[ii]==CLASS_SOIL      ){soil_param_count++;}
    else if (aPC[ii]==CLASS_VEGETATION){ veg_param_count++;}
    else if (aPC[ii]==CLASS_LANDUSE   ){lult_param_count++;}
    else if (aPC[ii]==CLASS_TERRAIN   ){terrain_param_count++;}
  }

  ExitGracefullyIf((soil_param_count>0) && (num_parsed_soils==0),
    "ParsePropertyFile: Soil parameters needed for model, but no soil classes have been defined",BAD_DATA_WARN);
  ExitGracefullyIf((lult_param_count>0) && (num_parsed_lult==0),
    "ParsePropertyFile: Land Use parameters needed for model, but no terrain classes have been defined",BAD_DATA_WARN);
  ExitGracefullyIf((veg_param_count>0) && (num_parsed_veg==0),
    "ParsePropertyFile: Vegetation parameters needed for model, but no vegetation classes have been defined",BAD_DATA_WARN);
  ExitGracefullyIf((terrain_param_count>0) && (num_parsed_terrs==0),
    "ParsePropertyFile: Terrain parameters needed for model, but no terrain classes have been defined",BAD_DATA_WARN);

  delete [] aP;
  delete [] aPC;
  
  if (!Options.silent){cout<<"...Done Checking"<<endl;}
   
  //Check if there is at least one class for each major classification genre
  //Won't be needed if GetParticipatingParamList() populated for all processes
  ExitGracefullyIf(num_parsed_soils==0,
    "No soil classes specified in .rvp file. Cannot proceed.",BAD_DATA_WARN);
  ExitGracefullyIf(num_parsed_veg==0,
    "No vegetation classes specified in .rvp file. Cannot proceed.",BAD_DATA_WARN);
  ExitGracefullyIf(num_parsed_lult==0,
    "No land use classes specified in .rvp file. Cannot proceed.",BAD_DATA_WARN);
  ExitGracefullyIf(num_parsed_profiles==0,
    "No soil profiles specified in .rvp file. Cannot proceed.",BAD_DATA_WARN);

  //CVegetationClass::SummarizeToScreen();

  
  ExitGracefullyIf((CChannelXSect::GetNumChannelXSects()==0) && (Options.routing!=ROUTE_NONE),
  "No channel profiles specified in .rvp file. Cannot proceed.",BAD_DATA_WARN);
  delete [] indices;
  for (int i=0;i<MAX_NUM_IN_CLASS;i++){delete [] properties[i];}delete [] properties;

  return true;
}

///////////////////////////////////////////////////////////////////
/// \brief Parses array of class properties
/// \details Parses array of class properties in the form (starting at CLASS_TAG1): \n
/// :Properties\n
///     CLASS_TAG1, v1, v2, v3, v4,... \n
///     CLASS_TAG2, v1, v2, v3, v4,... \n
///     ...\n
///     CLASS_TAGM, v1, v2, v3, v4,... \n
///   :EndProperties
/// \param *p [in] File parsing object
/// \param *indices [out] Index number of class (with reference to tags[] array)
/// \param **properties [out] array of model properties
/// \param &num_read [out] Number of rows read
/// \param *tags [in] Array of tag names
/// \param line_length [in] Number of expected columns in array
/// \param max_classes [in] Maximum allowable rows in array
/// \return True if operation is successful
//
bool ParsePropArray(CParser   *p,           //parser
                    int       *indices,     //output: index number of class (with reference to tags[] array)
                    double   **properties,  //output: array of properties
                    int       &num_read,    //output: number of rows read
                    string    *tags,        //array of tag names
                    const int  line_length, //number of expected columns in array
                    const int  max_classes) //maximum allowable rows in array
{
  int Len;
  char *s[MAXINPUTITEMS];
  p->Tokenize(s,Len);
  bool done=false;
  num_read=0;

  while (!done)
  { 
    if      (Len==0){}
    else if (!strcmp(s[0],"*")){} 
    else if (!strcmp(s[0],"#")){}
    else if (Len>=line_length)
    {
      int index=DOESNT_EXIST;
      for (int i=0;i<max_classes;i++){
        //cout<<tags[i]<<"<-->"<<s[0]<<endl;
        if (!tags[i].compare(s[0])){index=i;break;}
      }
      if (index==DOESNT_EXIST){
        cout <<"bad tag:"<<s[0]<<endl;
        return true;
      }//invalid class index
      indices[num_read]=index;

      for (int j=1;j<line_length;j++){
        properties[num_read][j-1]=AutoOrDouble(s[j]);
      }
      num_read++;
    }
    else{
      p->ImproperFormat(s); break;
    }
    p->Tokenize(s,Len);
    if (!strncmp(s[0],":End",4)){done=true;}
  }
  return false;
}
///////////////////////////////////////////////////////////////////
/// \brief Given list of nP needed parameters aP[] of type aPC[], warns user if some are not found associated with classes
///
/// \param *aP [out] array of needed parameter names needed
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters in list (size of aP[] and aPC[])
/// \param &Options global options structure
//
void  RVPParameterWarning   (string  *aP, class_type *aPC, int &nP, const optStruct &Options) 
{
  for (int ii=0;ii<nP;ii++)
  {
    if (Options.noisy){cout<<"    checking availability of parameter "<<aP[ii]<<endl;}

    if (aPC[ii]==CLASS_SOIL){
      for (int c=0;c<CSoilClass::GetNumClasses();c++){
        if (CSoilClass::GetSoilClass(c)->GetSoilProperty(aP[ii])==NOT_SPECIFIED){
          string warning="ParsePropertyFile: required soil property "+aP[ii]+" not included in .rvp file for class "+CSoilClass::GetSoilClass(c)->GetTag();
          ExitGracefully(warning.c_str(),BAD_DATA_WARN);
        }
      }
    }
    else if (aPC[ii]==CLASS_VEGETATION){
      for (int c=0;c<CVegetationClass::GetNumClasses();c++){
        if (CVegetationClass::GetVegClass(c)->GetVegetationProperty(aP[ii])==NOT_SPECIFIED){
          string warning="ParsePropertyFile: required vegetation property "+aP[ii]+" not included in .rvp file for vegetation class "+CVegetationClass::GetVegClass(c)->GetTag();
          ExitGracefully(warning.c_str(),BAD_DATA_WARN);
        }
      }
    }
    else if (aPC[ii]==CLASS_LANDUSE){
      for (int c=0;c<CLandUseClass::GetNumClasses();c++){
        if (CLandUseClass::GetLUClass(c)->GetSurfaceProperty(aP[ii])==NOT_SPECIFIED){
            
          string warning="ParsePropertyFile: required land use/land type property "+aP[ii]+" not included in .rvp file for land use class "+CLandUseClass::GetLUClass(c)->GetTag();
          ExitGracefully(warning.c_str(),BAD_DATA_WARN);
        }
      }
    }
    else if (aPC[ii]==CLASS_TERRAIN){
      for (int c=0;c<CTerrainClass::GetNumClasses();c++){
        if (CTerrainClass::GetTerrainClass(c)->GetTerrainProperty(aP[ii])==NOT_SPECIFIED){ 
          string warning="ParsePropertyFile: required terrain property "+aP[ii]+" not included in .rvp file for terrain class "+CTerrainClass::GetTerrainClass(c)->GetTag();
          ExitGracefully(warning.c_str(),BAD_DATA_WARN);
        }
      }
    }
    else if (aPC[ii]==CLASS_GLOBAL){
      if (CGlobalParams::GetParameter(aP[ii])==NOT_SPECIFIED){
            
          string warning="ParsePropertyFile: required global parameter "+aP[ii]+" not included in .rvp file.";
          ExitGracefully(warning.c_str(),BAD_DATA_WARN);
        }
    }
  }
}
///////////////////////////////////////////////////////////////////
/// \brief Append list of nP needed parameters aP[] of type aPC[] to master list aPm[],aPCm[] of size aPm
/// \note may later want to only append non-repetitive parameters
///
/// \param *aPm [in/out] master array of needed parameter names 
/// \param *aPCm [in/out] Class type (soil, vegetation, landuse or terrain) corresponding to each needed parameter in master list
/// \param &nPm [in/out] Number of parameters in master list (initial size of aPm[] and aPCm[])
/// \param *aP [in] array of needed parameter names to be appended
/// \param *aPC [in] Class type (soil, vegetation, landuse or terrain) corresponding to each needed parameter in appended list
/// \param &nP [in] Number of parameters in list to be appended (size of aP[] and aPC[])
//
void  AddToMasterParamList   (string        *&aPm, class_type       *&aPCm, int       &nPm, 
                              const string  *aP , const class_type *aPC , const int &nP) 
{
  if (nP==0){return;}
  
  string      *aPm_new=NULL;
  class_type *aPCm_new=NULL;

  aPm_new=new string     [nPm+nP];
  aPCm_new=new class_type [nPm+nP];
  if (aPCm_new==NULL){ExitGracefully("DynArrayAppend::Out of memory",OUT_OF_MEMORY);}

  for (int i=0;i<nPm;i++){
    aPm_new [i]=aPm [i];
    aPCm_new[i]=aPCm[i];
  }
  for (int i=0;i<nP;i++){
    aPm_new [nPm+i]=aP[i];
    aPCm_new[nPm+i]=aPC[i];
  }
  if (aPm!=NULL){delete [] aPm; aPm=NULL;}
  if (aPm!=NULL){delete [] aPCm;aPCm=NULL;}
  aPm =aPm_new; 
  aPCm=aPCm_new;
  nPm=nPm+nP;
}

  