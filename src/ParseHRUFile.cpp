/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Model.h"
#include "StateVariables.h"
#include "HydroUnits.h"
#include "ParseLib.h"

CReservoir *ReservoirParse(CParser *p,string name,int &HRUID,const optStruct &Options);

//////////////////////////////////////////////////////////////////
/// \brief Parses HRU properties file
/// \details model.rvh: input file that defines HRU, subbasin properties, number of HRUs \n
/// Rules: \n
/// - SubBasins command before HRUs command
/// - SubBasins command before SubBasinProperties, etc.
/// - HRUs command before HRUInitialConditions
///
/// \param *&pModel [out] Reference to model object
/// \param &Options [out] Global model options information
/// \return True if operation is successful
//
bool ParseHRUPropsFile(CModel *&pModel, const optStruct &Options, bool terrain_required)
{
  int         i;                //counters
  long        SBID;             //subbasin ID
  CHydroUnit *pHRU;             //temp pointers
  CSubBasin  *pSB;
  bool        ended=false;
  bool        in_ifmode_statement=false;

  ifstream    HRU;
  HRU.open(Options.rvh_filename.c_str());
  if (HRU.fail()){
    cout << "ERROR opening file: "<<Options.rvh_filename<<endl; return false;}

  int   Len,line(0),code;
  char *s[MAXINPUTITEMS];
  CParser *pp=new CParser(HRU,Options.rvh_filename,line);

  string            aParamStrings[MAXINPUTITEMS];
  int               nParamStrings=0;

  ifstream INPUT2;           //For Secondary input
  CParser *pMainParser=NULL; //for storage of main parser while reading secondary files

  if (Options.noisy){
    cout <<"======================================================"<<endl;
    cout <<"Parsing HRU Input File "<< Options.rvh_filename <<"..."<<endl;
    cout <<"======================================================"<<endl;
  }

  //--Sift through file-----------------------------------------------
  bool end_of_file=pp->Tokenize(s,Len);
  while (!end_of_file)
  {
    if (ended){break;}
    if (Options.noisy){ cout << "reading line " << pp->GetLineNumber() << ": ";}

    /*assign code for switch statement
      ------------------------------------------------------------------
      <100         : ignored/special
      0   thru 100 : All other
      ------------------------------------------------------------------
    */

    code=0;
    //---------------------SPECIAL -----------------------------
    if       (Len==0)                                    {code=-1; }//blank line
    else if  (IsComment(s[0],Len))                       {code=-2; }//comment
    else if  (!strcmp(s[0],":End"                      )){code=-4; }//stop reading
    else if  (!strcmp(s[0],":IfModeEquals"             )){code=-5; }
    else if  (in_ifmode_statement)                       {code=-6; }
    else if  (!strcmp(s[0],":RedirectToFile"           )){code=-3; }//redirect to secondary file
    //--------------------MODEL OPTIONS ------------------------
    else if  (!strcmp(s[0],":SubBasins"                )){code=1;  }
    else if  (!strcmp(s[0],":HRUs"                     )){code=2;  }
    else if  (!strcmp(s[0],":Reservoir"                )){code=3;  }
    else if  (!strcmp(s[0],":SubBasinProperties"       )){code=7;  }
    else if  (!strcmp(s[0],":HRUGroup"                 )){code=8;  }
    else if  (!strcmp(s[0],":PopulateHRUGroup"         )){code=9;  }
    else if  (!strcmp(s[0],":SubBasinGroup"            )){code=10; }
    else if  (!strcmp(s[0],":SBGroupPropertyMultiplier")){code=11; }
    else if  (!strcmp(s[0],":SBGroupPropertyOverride"  )){code=12; }
    else if  (!strcmp(s[0],":DisableHRUGroup"          )){code=13; }
    else if  (!strcmp(s[0],":DisableSubBasinGroup"     )){code=14; }
    else if  (!strcmp(s[0],":PopulateSubBasinGroup"    )){code=15; } 
    
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

      filename=CorrectForRelativePath(filename,Options.rvt_filename);

      INPUT2.open(filename.c_str());
      if (INPUT2.fail()){
        string warn=":RedirectToFile: Cannot find file "+filename;
        ExitGracefully(warn.c_str(),BAD_DATA);
      }
      else{
        pMainParser=pp;   //save pointer to primary parser
        pp=new CParser(INPUT2,filename,line);//open new parser
      }
      break;
    }
    case(-4):  //----------------------------------------------
    {/*:End*/
      if (Options.noisy) {cout <<"EOF"<<endl;} ended=true; break;
    }
    case(-5):  //----------------------------------------------
    {/*:IfModeEquals*/
      if(Len>1) {
        if(Options.noisy) { cout <<"Mode statement start..."<<endl; }
        char testmode=s[1][0];
        if(testmode!=Options.run_mode) {
          in_ifmode_statement=true;
        }
      }
      break;
    }
    case(-6):  //----------------------------------------------
    {/*in_ifmode_statement*/
      if(Options.noisy) { cout <<"* skip"<<endl; }
      if(!strcmp(s[0],":EndIfModeEquals"))
      {
        if(Options.noisy) { cout <<"...Mode statement end"<<endl; }
        in_ifmode_statement=false;
      }
      break;
    }
    case(1):  //----------------------------------------------
    { /*
        ":SubBasins"
        {int ID , string (no spaces) name, int down_ID, string profile, double length [km], bool gauged,{optional double Qref}}x nSubBasins
        :EndSubBasins
        -down_ID=-1 for basins not draining into other modeled basins
        -profile can be 'NONE' for ROUTE_NONE option, but must be linked to actual profile otherwise
      */
      if (Options.noisy) {cout <<"Subbasin data..."<<endl;}
      if (Len!=1){pp->ImproperFormat(s);}
      else{
        while (((Len==0) || (strcmp(s[0],":EndSubBasins"))) && (!end_of_file))
        {
          end_of_file=pp->Tokenize(s,Len);
          if      (IsComment(s[0],Len))          {}//comment line
          else if (!strcmp(s[0],":Attributes"  )){}//ignored by Raven - needed for GUIs
          else if (!strcmp(s[0],":Units"       )){}//ignored by Raven - needed for GUIs
          else if (!strcmp(s[0],":EndSubBasins")){}//done
          else
          {
            if (Len<6){pp->ImproperFormat(s);}

            CChannelXSect const *pChan=NULL;
            pChan=CChannelXSect::StringToChannelXSect(string(s[3]));
            string error,error2;
            error="Parse HRU File: Unrecognized Channel profile code ("+string(s[3])+") in SubBasins command";
            error2="Parse HRU File: NONE cannot be used as channel code if routing method is anything other than ROUTE_NONE";
            ExitGracefullyIf((pChan==NULL) && (string(s[3])!="NONE") && (Options.routing!=ROUTE_NONE),error.c_str() ,BAD_DATA_WARN);
            ExitGracefullyIf((pChan==NULL) && (string(s[3])=="NONE") && (Options.routing!=ROUTE_NONE),error2.c_str(),BAD_DATA_WARN);

            if(!StringIsLong(s[0])){
              error="Parse HRU File: Subbasin ID \""+string(s[0])+"\" must be unique integer or long integer";
              ExitGracefully(error.c_str(),BAD_DATA_WARN);
            }
           
            double length;
            length=AutoOrDouble(s[4]);
            ExitGracefullyIf(length>2000,
              "ParseHRUPropsFile: length of river greater than 1000km in :SubBasin command block. Units issue? Reach length should be provided in km.",BAD_DATA_WARN);
            if (length!=AUTO_COMPUTE){length*=M_PER_KM;}//convert to m from km

            bool gaged;
            gaged=s_to_b(s[5]);

            double Qref;
            if (Len==7){Qref=AutoOrDouble(s[6]);}
            else       {Qref=AUTO_COMPUTE;}

            pSB=NULL;
            pSB=new CSubBasin(s_to_l(s[0]),s[1], pModel,s_to_l(s[2]),pChan,length,Qref,gaged);
            ExitGracefullyIf(pSB==NULL,"ParseHRUPropsFile",OUT_OF_MEMORY);
            pModel->AddSubBasin(pSB);
          }
        }
      }
      break;
    }
    case(2):  //----------------------------------------------
    { /*
        ":HRUs"
        {ID, area, avg_elevation, lat, long, basinID, ..
        LU class,veg class, soil profile, aquifer profile, terrain class, slope[deg], aspect[deg from N]}x nHRUs
        :EndHRUs
      */
      if (Options.noisy) {cout <<"HRU data..."<<endl;}
      while (((Len==0) || (strcmp(s[0],":EndHRUs"))) && (!end_of_file))
      {
        end_of_file=pp->Tokenize(s,Len);
        if      (IsComment(s[0],Len))          {}//comment line
        else if (!strcmp(s[0],":Attributes"  )){}//ignored by Raven - needed for GUIs
        else if (!strcmp(s[0],":Units"       )){}//ignored by Raven - needed for GUIs
        else if (!strcmp(s[0],":EndHRUs")){}//done
        else
        {
          if (Len<13){pp->ImproperFormat(s);}

          string error;
          SBID =s_to_l(s[5]);//index must exist (in file, from 1 to nSB)

          if(!StringIsLong(s[0])){
            error="Parse HRU File: HRU ID \""+string(s[0])+"\" in :HRUs table must be unique integer or long integer";
            ExitGracefully(error.c_str(),BAD_DATA_WARN);
          }
          if(!StringIsLong(s[5])){
            error="Parse HRU File: Subbasin ID \""+string(s[5])+"\" in :HRUs table must be unique integer or long integer ";
            ExitGracefully(error.c_str(),BAD_DATA_WARN);
          }
          if((fabs(s_to_d(s[3]))>90) || (fabs(s_to_d(s[4]))>180)){
            error="Parse HRU File: longitude and latitude in :HRUs table (HRU "+to_string(s[0])+") must be valid (-90<lat<90; -180<lon<180)";
            ExitGracefully(error.c_str(),BAD_DATA_WARN);
          }
          if((s_to_d(s[11])>90) || (s_to_d(s[11])<0)){
            error="Parse HRU File: slope in :HRUs table (HRU "+to_string(s[0])+") must be valid (0<slope<90)";
            ExitGracefully(error.c_str(),BAD_DATA_WARN);
          }
          if((s_to_d(s[12])>360) || (s_to_d(s[12])<0)){
            error="Parse HRU File: aspect in :HRUs table (HRU "+to_string(s[0])+") must be valid (0<slope<360)";
            ExitGracefully(error.c_str(),BAD_DATA_WARN);
          }
          CLandUseClass const *pLULT=NULL;
          pLULT=CLandUseClass::StringToLUClass(string(s[6]));
          if (pLULT==NULL){
            error="Parse HRU File: Unrecognized Land Use Code/index: \""+string(s[6])+"\"";
            ExitGracefully(error.c_str(),BAD_DATA);
          }

          CVegetationClass const *pVegetation=NULL;
          pVegetation=CVegetationClass::StringToVegClass(string(s[7]));
          if (pVegetation==NULL){
            error="Parse HRU File: Unrecognized Vegetation Code/index: \""+string(s[7])+"\"";
            ExitGracefully(error.c_str(),BAD_DATA);
          }

          CSoilProfile const *pSoilProfile=NULL;
          pSoilProfile=CSoilProfile::StringToSoilProfile(string(s[8]));
          if (pSoilProfile==NULL){
            error="Parse HRU File: Unrecognized Soil Profile Code/index: \""+string(s[8])+"\"";
            ExitGracefully(error.c_str(),BAD_DATA);
          }

          void *pAqStack=NULL; //GWMIGRATE
          /*CAquiferStack const *pAqStack=NULL;
            if (string(s[9])!="[NONE]") //Aquifer profile can be NULL
            {
            pAqStack=CAquiferStack::StringToAquiferStack(string(s[9]));
            if (pAqStack==NULL){
            error="Parse HRU File: Unrecognized Soil Profile Code/index: \""+string(s[9])+"\"";
            ExitGracefully(error.c_str(),BAD_DATA_WARN);
            }
          }*/

          CTerrainClass const *pTerrain=NULL;
          if (string(s[10])!="[NONE]") //Terrain class can be NULL
          {
            pTerrain=CTerrainClass::StringToTerrainClass(string(s[10]));
            if (pTerrain==NULL){
              error="Parse HRU File: Unrecognized Terrain Code/index: \""+string(s[10])+"\"";
              ExitGracefully(error.c_str(),BAD_DATA_WARN);
            }
          }

          HRU_type HRUtype=HRU_STANDARD;
          if (!pSoilProfile->GetTag().substr(0,4).compare("LAKE"   )){HRUtype=HRU_LAKE;   }
          if (!pSoilProfile->GetTag().substr(0,7).compare("GLACIER")){HRUtype=HRU_GLACIER;}
          if (!pSoilProfile->GetTag().substr(0,4).compare("ROCK"   )){HRUtype=HRU_ROCK;   }
          if (!pSoilProfile->GetTag().substr(0,7).compare("WETLAND")){HRUtype=HRU_WETLAND;}
          pHRU=new CHydroUnit( pModel,
                               s_to_i(s[0]),//ID
                               pModel->GetNumHRUs(),//k - global model index
                               s_to_d(s[1]),//area
                               pModel->GetSubBasinIndex(SBID),
                               s_to_d(s[2]),//elev
                               s_to_d(s[3]),//lat
                               s_to_d(s[4]),//long
                               s_to_d(s[11])*DEGREES_TO_RADIANS,//slope (deg->radians)
                               s_to_d(s[12])*DEGREES_TO_RADIANS,//aspect (deg->radians)
                               HRUtype,
                               pSoilProfile,
							                 //pAqStack, GWMIGRATE
                               pVegetation,
                               pAqStack,
                               pTerrain,
                               pLULT);
          ExitGracefullyIf(pHRU==NULL,"ParseHRUPropsFile",OUT_OF_MEMORY);
          pModel->AddHRU(pHRU);
          pSB=pModel->GetSubBasinByID(SBID);
          if (pSB!=NULL){pSB->AddHRU(pHRU);}
          else          {
            ExitGracefully("ParseHRUProps: Bad Sub-basin index in :HRUs command",BAD_DATA);
          }
        }
      }//end while
      break;
    }
    case(3):  //----------------------------------------------
    {
      /*
        :Reservoir [name]
          :SubBasin [SBID]
          :HRUID [HRUID]
          :VolumeHeightRelation
            [DATA BLOCK]
          :EndVolumeHeightRelation
          :OutflowHeightRelation
          :EndOutflowHeightRelation
          ...
        :EndReservoir
        (and other formats - see ParseReservoir code)
      */
      if (Options.noisy) {cout <<":Reservoir"<<endl;}
      CReservoir *pRes;
      int HRUID;
      pRes=ReservoirParse(pp,s[1],HRUID,Options);
      pSB=pModel->GetSubBasinByID(pRes->GetSubbasinID());
      pRes->SetHRU(pModel->GetHRUByID(HRUID));
      if (pSB!=NULL){pSB->AddReservoir(pRes);}
      else          {
        ExitGracefully("ParseHRUProps: Bad Sub-basin index in :Reservoir command",BAD_DATA);
      }
      break;
    }
    case(7):  //----------------------------------------------
    { /*
        ":SubBasinProperties"
          :Parameters, paramname1,paramname2,...,paramnameN
          :Units     ,  units1, units2, ..., unitsN
          {ID,param1, param2,...,paramN} x nSubBasins (or a subset of SBs)
        :EndSubBasinProperties
      */
      if (Options.noisy) {cout <<"   Reading Basin Properties..."<<endl;}
      bool done;
      done=false;
      while (!done)
      {
        pp->Tokenize(s,Len);
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
        else {pp->ImproperFormat(s); break;}
      }
      if (Options.noisy){for (i=1;i<nParamStrings;i++){cout<<"    "<<aParamStrings[i]<<endl;}}
      bool good_string;
      while ((Len==0) || (strcmp(s[0],":EndSubBasinProperties")))
      {
        pp->Tokenize(s,Len);
        if      (IsComment(s[0],Len))          {}//comment line
        else if (!strcmp(s[0],":RedirectToFile")){
          ExitGracefully("Parse HRU File: :RedirectToFile cannot be inside at :SubBasinProperties block.",BAD_DATA_WARN);
        }//done
        else if (!strcmp(s[0],":EndSubBasinProperties")){}//done
        else
        {
          ExitGracefullyIf(Len<nParamStrings,
                           "Parse HRU File: incorrect number of terms in SubBasin properties",BAD_DATA);

          SBID=s_to_l(s[0]);

          pSB=pModel->GetSubBasinByID(SBID);
          if (pSB!=NULL){
            for (i=1;i<nParamStrings;i++){
              double in=AutoOrDouble(s[i]);

              //special handling of some parameters
              if(!aParamStrings[i].compare("TIME_CONC")    && (in!=AUTO_COMPUTE) && (in!=USE_TEMPLATE_VALUE)){
                in*=CGlobalParams::GetParameter("TOC_MULTIPLIER");
              }
              if(!aParamStrings[i].compare("TIME_TO_PEAK") && (in!=AUTO_COMPUTE) && (in!=USE_TEMPLATE_VALUE)){
                in*=CGlobalParams::GetParameter("TOC_MULTIPLIER");
              }
              if(!aParamStrings[i].compare("REACH_HRU_ID")) {
                if(pModel->GetHRUByID((int)in)!=NULL) {
                  in=pModel->GetHRUByID((int)in)->GetGlobalIndex(); //Convert ID to index
                }
                else {
                  ExitGracefully("ParseHRUPropsFile::invalid REACH_HRU_ID in :SubBasinProperties command",BAD_DATA_WARN);
                }
              }
              //end special handling 

              good_string=pSB->SetBasinProperties(aParamStrings[i],in);
              if (!good_string)
              {
                string err;
                err="Unknown parameter \""+aParamStrings[i]+"\" in :SubBasinProperties command";
                ExitGracefully(err.c_str(),BAD_DATA_WARN);
              }
            }
          }
          else{
            string warn;
            warn="Subbasin "+to_string(SBID)+" not in model, cannot set subbasin properties";
            WriteWarning(warn,Options.noisy);
          }
        }
      }
      break;
    }
    case(8):  //----------------------------------------------
    { /*
        ":HRUGroup" {name}
         {ID1,ID2,ID3,...} x nHRUs in group
        :EndHRUGroup
      */
      if (Options.noisy) {cout <<"   HRU Group..."<<endl;}
      if (Len!=2){pp->ImproperFormat(s);}
      CHRUGroup *pHRUGrp=NULL;
      pHRUGrp=pModel->GetHRUGroup(s[1]);
      if (pHRUGrp==NULL){//group not yet defined
        WriteAdvisory("HRU groups should ideally be defined in .rvi file (using :DefineHRUGroup(s) commands) before being populated in .rvh file",Options.noisy);
        pHRUGrp=new CHRUGroup(s[1],pModel->GetNumHRUGroups());
        pModel->AddHRUGroup(pHRUGrp);
      }
      while ((Len==0) || (strcmp(s[0],":EndHRUGroup")))
      {
        pp->Tokenize(s,Len);
        if      (IsComment(s[0], Len)){}//comment line
        else if (!strcmp(s[0],":EndHRUGroup")){}//done
        else
        {
          int k;
          int nHRUs=pModel->GetNumHRUs();
          for (i=0;i<Len;i++)
          {
            int ind1,ind2;
            s_to_range(s[i],ind1,ind2);
            bool found=false;
            ExitGracefullyIf((ind2-ind1)>10000,"Parsing :HRUGroup command: invalid range of HRU indices",BAD_DATA);
            bool gaps(false);
            for (int ii=ind1;ii<=ind2;ii++)
            {
              found=false;
              for (k=0;k<nHRUs;k++)
              {
                if (pModel->GetHydroUnit(k)->GetID()==ii)
                {
                  pHRUGrp->AddHRU(pModel->GetHydroUnit(k));found=true; break;
                }
              }
              if ((!found) && (ind2-ind1>1)){gaps=true;}
              if((!found) && (ind2==ind2)) {
                string warn="HRU ID "+to_string(ii)+" specified in :HRUGroup command for group "+pHRUGrp->GetName()+ " does not exist";
                WriteWarning(warn.c_str(),Options.noisy);
              }
            }
            if (gaps){
              WriteWarning("Range specified in :HRUGroup command has gaps",Options.noisy);
            }
          }
        }
      }
      break;
    }
    case(9):  //----------------------------------------------
    { /*
        ":PopulateHRUGroup" {name} "With" {conditionbase} {condition} {conditiondata}
        e.g.,
        :PopulateHRUGroup NotRock With HRUS NOTWITHIN RockHRUGroup
        
        :PopulateHRUGroup LowBand With ELEVATION BETWEEN  0 500
        :PopulateHRUGroup CroplandHRUs With LANDUSE EQUALS CROPLAND
        :PopulateHRUGroup NonCroplandHRUs With LANDUSE NOTEQUALS CROPLAND
        :PopulateHRUGroup BroadleafHRUs With VEGETATION  EQUALS  BROADLEAF
        :PopulateHRUGroup UpperWigwamHRUs With HRUS WITHIN_SBGROUP UpperWigwamSBs
        [not in here yet] :PopulateHRUGroup Rocks With HRUTYPE EQUALS ROCK
      */
      if (Options.noisy) {cout <<"   Populate HRU Group..."<<endl;}
      if (Len<6){pp->ImproperFormat(s);}

      CHRUGroup *pHRUGrp=NULL;
      pHRUGrp=pModel->GetHRUGroup(s[1]);

      if (pHRUGrp==NULL){//group not yet defined
        WriteWarning("HRU groups should ideally be defined in .rvi file (using :DefineHRUGroup(s) commands) before being populated in .rvh file (2)",Options.noisy);
        pHRUGrp=new CHRUGroup(s[1],pModel->GetNumHRUGroups());
        pModel->AddHRUGroup(pHRUGrp);
      }
      int k;
      if(!strcmp(s[3],"HRUS"))
      {
        if(!strcmp(s[4],"NOTWITHIN")){
          CHRUGroup* pHRUGrp2=NULL;
          pHRUGrp2=pModel->GetHRUGroup(s[5]);
          if(pHRUGrp2==NULL) {
            ExitGracefully(":PopulateHRUGroup: invalid HRU group reference used in command",BAD_DATA_WARN);
          }
          else {
            for(k=0;k<pModel->GetNumHRUs();k++)
            {
              if(!pHRUGrp2->IsInGroup(k)) { pHRUGrp->AddHRU(pModel->GetHydroUnit(k)); }
            }
          }
        }
        else if (!strcmp(s[4],"WITHIN_SBGROUP")) 
        {
          CSubbasinGroup *pSBGrp=NULL;
          pSBGrp=pModel->GetSubBasinGroup(s[5]);
          if(pSBGrp==NULL) {
            ExitGracefully(":PopulateHRUGroup: invalid Subbasin group reference used in command",BAD_DATA_WARN);
          }
          else {
            for(int p=0;p<pSBGrp->GetNumSubbasins();p++)
            {
              for(k=0;k<pSBGrp->GetSubBasin(p)->GetNumHRUs();k++)
              {
                pHRU=pModel->GetHRUByID(pSBGrp->GetSubBasin(p)->GetHRU(k)->GetID());
                pHRUGrp->AddHRU(pHRU);
              }
            }
          }
        }
      }
      else if(!strcmp(s[3],"LANDUSE"))
      {
        if(!strcmp(s[4],"EQUALS")){
          for(k=0;k<pModel->GetNumHRUs();k++)
          {
            if(pModel->GetHydroUnit(k)->GetSurfaceProps()->landuse_name==to_string(s[5])){ pHRUGrp->AddHRU(pModel->GetHydroUnit(k)); }
          }
        }
        if(!strcmp(s[4],"DOESNTEQUAL")){
          for(k=0;k<pModel->GetNumHRUs();k++)
          {
            if(pModel->GetHydroUnit(k)->GetSurfaceProps()->landuse_name!=to_string(s[5])){ pHRUGrp->AddHRU(pModel->GetHydroUnit(k)); }
          }
        }
      }
      else if(!strcmp(s[3],"VEGETATION"))
      {
        if(!strcmp(s[4],"EQUALS")){
          for(k=0;k<pModel->GetNumHRUs();k++)
          {
            if(pModel->GetHydroUnit(k)->GetVegetationProps()->vegetation_name==to_string(s[5])){ pHRUGrp->AddHRU(pModel->GetHydroUnit(k)); }
          }
        }
        if(!strcmp(s[4],"DOESNTEQUAL")){
          for(k=0;k<pModel->GetNumHRUs();k++)
          {
            if(pModel->GetHydroUnit(k)->GetVegetationProps()->vegetation_name!=to_string(s[5])){ pHRUGrp->AddHRU(pModel->GetHydroUnit(k)); }
          }
        }
      }
      else if(!strcmp(s[3],"ELEVATION"))
      {
        if(!strcmp(s[4],"BETWEEN")){
          for(k=0;k<pModel->GetNumHRUs();k++)
          {
            if((pModel->GetHydroUnit(k)->GetElevation()>=s_to_d(s[5])) &&
               (pModel->GetHydroUnit(k)->GetElevation()<s_to_d(s[6]))){ pHRUGrp->AddHRU(pModel->GetHydroUnit(k)); }
          }
        }
      }
      else if(!strcmp(s[3],"HRUTYPE"))
      {
        if(!strcmp(s[4],"EQUALS")){
          for(k=0;k<pModel->GetNumHRUs();k++)
          {
            if((pModel->GetHydroUnit(k)->GetHRUType()==HRU_ROCK) && (!strcmp(s[4],"ROCK"))){ pHRUGrp->AddHRU(pModel->GetHydroUnit(k)); }
            if((pModel->GetHydroUnit(k)->GetHRUType()==HRU_GLACIER) && (!strcmp(s[4],"GLACIER"))){ pHRUGrp->AddHRU(pModel->GetHydroUnit(k)); }
            if((pModel->GetHydroUnit(k)->GetHRUType()==HRU_LAKE) && (!strcmp(s[4],"LAKE"))){ pHRUGrp->AddHRU(pModel->GetHydroUnit(k)); }
          }
        }
        if(!strcmp(s[4],"DOESNTEQUAL")){
          for(k=0;k<pModel->GetNumHRUs();k++)
          {
            if((pModel->GetHydroUnit(k)->GetHRUType()!=HRU_ROCK) && (!strcmp(s[4],"ROCK"))){ pHRUGrp->AddHRU(pModel->GetHydroUnit(k)); }
            if((pModel->GetHydroUnit(k)->GetHRUType()!=HRU_GLACIER) && (!strcmp(s[4],"GLACIER"))){ pHRUGrp->AddHRU(pModel->GetHydroUnit(k)); }
            if((pModel->GetHydroUnit(k)->GetHRUType()!=HRU_LAKE) && (!strcmp(s[4],"LAKE"))){ pHRUGrp->AddHRU(pModel->GetHydroUnit(k)); }
          }
        }
      }
      break;
    }
    case(10):  //----------------------------------------------
    { /*
        ":SubBasinGroup" {name}
         {SBID1,SBID2,SBID3,...} x nBasins in group
        :EndSubBasinGroup
      */
      if (Options.noisy) {cout <<"   SubBasin Group..."<<endl;}
      if (Len!=2){pp->ImproperFormat(s);}
      CSubbasinGroup *pSBGrp=NULL;
      pSBGrp=pModel->GetSubBasinGroup(s[1]);
      if (pSBGrp==NULL){//group not yet defined
      //  WriteAdvisory("HRU groups should ideally be defined in .rvi file (using :DefineHRUGroup(s) commands) before being populated in .rvh file",Options.noisy);
        pSBGrp=new CSubbasinGroup(s[1],pModel->GetNumSubBasinGroups());
        pModel->AddSubBasinGroup(pSBGrp);
      }
      bool eof=false;
      while ((Len==0) || (strcmp(s[0],":EndSubBasinGroup")))
      {
        eof=pp->Tokenize(s,Len);
        if(eof) { break; }
        if      (IsComment(s[0], Len)){}//comment line
        else if (!strcmp(s[0],":EndSubBasinGroup")){}//done
        else
        {
          int p;
          int nSBs=pModel->GetNumSubBasins();
          for(i=0;i<Len;i++)
          {
            long SBID=s_to_l(s[i]);
            for(p=0;p<nSBs;p++)
            {
              if(pModel->GetSubBasin(p)->GetID()==SBID)
              {
                pSBGrp->AddSubbasin(pModel->GetSubBasin(p));break;
              }
            }
          }
        }
      }
      break;
    }
    case(11):  //----------------------------------------------
    { /*
      :SBGroupPropertyMultiplier [SBGROUP] [PROPERTY] [multiplier] 
      */
      CSubbasinGroup *pSBGrp;
      if(Len>=4) {
        pSBGrp=pModel->GetSubBasinGroup(s[1]);
        if(pSBGrp==NULL) {
          WriteWarning(":SBGroupPropertyMultiplier: invalid subbasin group specified",Options.noisy);
          break;
        }
        for(int p=0;p<pSBGrp->GetNumSubbasins();p++) {
          double val=pSBGrp->GetSubBasin(p)->GetBasinProperties(s[2]);
          if(val!=AUTO_COMPUTE) {
            pSBGrp->GetSubBasin(p)->SetBasinProperties(s[2],s_to_d(s[3])*val);
          }
        }
      }
      else {
        WriteWarning(":SBGroupPropertyMultiplier: improper syntax",Options.noisy);
      }
      break;
    }
    case(12):  //----------------------------------------------
    { /*
      :SBGroupPropertyOverride [SBGROUP] [PROPERTY] [value]
      */
      CSubbasinGroup *pSBGrp;
      if(Len>=4) {
        pSBGrp=pModel->GetSubBasinGroup(s[1]);
        if(pSBGrp==NULL) {
          WriteWarning(":SBGroupPropertyOverride: invalid subbasin group specified",Options.noisy);
          break;
        }
        for(int p=0;p<pSBGrp->GetNumSubbasins();p++) {
          pSBGrp->GetSubBasin(p)->SetBasinProperties(s[2],s_to_d(s[3]));          
        }
      }
      else {
        WriteWarning(":SBGroupPropertyOverride: improper syntax",Options.noisy);
      }
      break;
    }
    case(13):  //----------------------------------------------
    { /*
      :DisableHRUGroup [HRUGROUP] 
      */
      if(Options.noisy) { cout <<"Disabling HRU Group"<<endl; }
      if(Len>=2) {
        CHRUGroup *pHRUGrp=NULL;
        pHRUGrp=pModel->GetHRUGroup(s[1]);
        if(pHRUGrp==NULL) {
          ExitGracefully("Invalid HRU Group name supplied in :DisableHRUGroup command in .rvh file. Group must be defined first.",BAD_DATA_WARN);
          break;
        }
        else {
          pHRUGrp->DisableGroup();
        }
      }
      else {
        WriteWarning(":DisableHRUGroup: improper syntax",Options.noisy);
      }
      break;
    }
    case(14):  //----------------------------------------------
    { /*
      :DisableSubBasinGroup [HRUGROUP]
      */
      if(Options.noisy) { cout <<"Disabling HRU Group"<<endl; }
      if(Len>=2) {
        CSubbasinGroup *pSBGrp=NULL;
        pSBGrp=pModel->GetSubBasinGroup(s[1]);
        if(pSBGrp==NULL) {
          ExitGracefully("Invalid Subbasin Group name supplied in :DisableSubBasinGroup command in .rvh file. Group must be defined first.",BAD_DATA_WARN);
          break;
        }
        else {
          pSBGrp->DisableGroup();
        }
      }
      else {
        WriteWarning(":DisableSubBasinGroup: improper syntax",Options.noisy);
      }
      break;
    }
    case(15):  //----------------------------------------------
    { /*
      ":PopulateSubBasinGroup" {name} "With" {conditionbase} {condition} {conditiondata}
      e.g.,
      :PopulateSubBasinGroup NotRock With SUBBASINS NOTWITHIN RockSBGroup
      :PopulateSubBasinGroup UpstreamOfBasin2 With SUBBASINS UPSTREAM_OF 2
      :PopulateSubBasinGroup DownstreamOfBasin35 With SUBBASINS DOWNSTREAM_OF 35
      */
      if(Options.noisy) { cout <<"   Populate SubBasin Group..."<<endl; }
      if(Len<6) { pp->ImproperFormat(s); }

      CSubbasinGroup *pSBGroup=NULL;
      pSBGroup=pModel->GetSubBasinGroup(s[1]);

      string advice="SubBasinGroup "+to_string(s[1])+" was populated with these basins :";

      if(pSBGroup==NULL) {//group not yet defined
        pSBGroup=new CSubbasinGroup(s[1],pModel->GetNumSubBasinGroups());
        pModel->AddSubBasinGroup(pSBGroup);
      }
      if(!strcmp(s[3],"SUBBASINS"))
      {
        if(!strcmp(s[4],"NOTWITHIN")) 
        {
          CSubbasinGroup *pSBGroup2=NULL;
          pSBGroup2=pModel->GetSubBasinGroup(s[5]);
          if(pSBGroup2==NULL) {
            ExitGracefully(":PopulateSubBasinGroup: invalid SB group reference used in command",BAD_DATA_WARN);
          }
          else {
            for(int p=0;p<pModel->GetNumSubBasins();p++)
            {
              if(!pSBGroup2->IsInGroup(pModel->GetSubBasin(p)->GetID())) { pSBGroup->AddSubbasin(pModel->GetSubBasin(p)); }
            }
          }
        }
        else if(!strcmp(s[4],"UPSTREAM_OF")) 
        {
          int SBID=s_to_l(s[5]);
          for(int p=0;p<pModel->GetNumSubBasins();p++)
          {
            if(pModel->IsSubBasinUpstream(pModel->GetSubBasin(p)->GetID(),SBID)) { 
              pSBGroup->AddSubbasin(pModel->GetSubBasin(p));
              advice=advice+to_string(pModel->GetSubBasin(p)->GetID())+" ";
            }
          }
        }
        else if(!strcmp(s[4],"DOWNSTREAM_OF")) 
        {
          CSubBasin *pBasin=pModel->GetSubBasinByID(s_to_l(s[5]));
          if(pBasin==NULL) {
            ExitGracefully(":PopulateSubBasinGroup : invalid subbasin ID specified",BAD_DATA_WARN);
          }
          else {
            int iter=0;
            while ((pBasin->GetDownstreamID()!=DOESNT_EXIST) && (iter<1000)){
              pBasin=pModel->GetSubBasinByID(pBasin->GetDownstreamID());
              pSBGroup->AddSubbasin(pBasin);
              advice=advice+to_string(pBasin->GetID())+" ";
              iter++;
            }
            if(iter==1000) {
              ExitGracefully(":PopulateSubBasinGroup: cyclical downstream references in :SubBasins list",BAD_DATA_WARN);
            }
          }
        }
      }
	  WriteAdvisory(advice,Options.noisy);
      break;
    }
    default://------------------------------------------------
    {
      char firstChar = *(s[0]);
      switch(firstChar)
      {
      case ':':
      {
        if     (!strcmp(s[0],":FileType"    )) {if (Options.noisy){cout<<"Filetype"<<endl;}}//do nothing
        else if(!strcmp(s[0],":Application" )) {if (Options.noisy){cout<<"Application"<<endl;}}//do nothing
        else if(!strcmp(s[0],":Version"     )) {if (Options.noisy){cout<<"Version"<<endl;}}//do nothing
        else if(!strcmp(s[0],":WrittenBy"   )) {if (Options.noisy){cout<<"WrittenBy"<<endl;}}//do nothing
        else if(!strcmp(s[0],":CreationDate")) {if (Options.noisy){cout<<"CreationDate"<<endl;}}//do nothing
        else if(!strcmp(s[0],":SourceFile"  )) {if (Options.noisy){cout<<"SourceFile"<<endl;}}//do nothing
        else if(Options.noisy)
        {
          string warn ="IGNORING unrecognized command: " + string(s[0])+ " in .rvh file";
          WriteWarning(warn,Options.noisy);
        }
      }
      break;
      default:
      {
        string errString = "Unrecognized command in .rvh file:\n   " + string(s[0]);
        ExitGracefully(errString.c_str(),BAD_DATA);//STRICT
      }
      break;
      }
    }
    }//end switch(code)

    end_of_file=pp->Tokenize(s,Len);

    //return after file redirect, if in secondary file
    if ((end_of_file) && (pMainParser!=NULL))
    {
      INPUT2.clear();
      INPUT2.close();
      delete pp;
      pp=pMainParser;
      pMainParser=NULL;
      end_of_file=pp->Tokenize(s,Len);
    }
  } //end while !end_of_file
  HRU.close();

  //QA/QC
  //--------------------------------------------------------------------------
  if ((pModel->GetNumSubBasins() > 1) && (Options.routing == ROUTE_NONE))
  {
    string warn;
    warn = "ParseHRUPropsFile: ROUTE_NONE is being used for a model with more than one basin. Typically this method is used for lumped models only.";
    WriteWarning(warn,Options.noisy);
  }

  // Add parameters needed for discharge initialization/reference flow calculation
  //--------------------------------------------------------------------------
  if ((pModel->GetNumSubBasins()>1) && (CGlobalParams::GetParameter("AVG_ANNUAL_RUNOFF")<0))
  {
    // \todo: reduce generalization- only really needed if routing method requires Q_REF
    ExitGracefully("ParseHRUPropsFile:: AVG_ANNUAL_RUNOFF should be supplied (using :AvgAnnualRunoff command in .rvp file) if more than one basin is included in model",BAD_DATA_WARN);
  }

  // Check if terrain class is needed but not specified
  //--------------------------------------------------------------------------
  if(terrain_required) { //this is determined at start of .rvp file read
    for(int k=0;k<pModel->GetNumHRUs(); k++) {
      if(pModel->GetHydroUnit(k)->GetTerrainProps()==NULL) {
        ExitGracefully("ParseHRUFile:: At least one model process requires terrain parameters. Cannot have NULL ([NONE]) Terrain class in :HRUs block",BAD_DATA_WARN);
      }
    }
  }

  delete pp;
  pp=NULL;
  return true;
}
//////////////////////////////////////////////////////////////////
/// \brief parses the :Reservoir Command
/// \note the first line (:Reservoir [name]) has already been read in
/// \param p [in] parser (likely points to .rvh file)
/// \param name [in] name of reservoir
/// \param Options [in] model options structure
//
CReservoir *ReservoirParse(CParser *p,string name,int &HRUID,const optStruct &Options)
{

  /*Examples:

  :Reservoir ExampleReservoir
    :SubBasinID 23
    :HRUID 234
    :Type RESROUTE_STANDARD
    :VolumeStageRelation POWER_LAW
      0.1 2.3
    :EndVolumeStageRelation
    :OutflowStageRelation POWER_LAW
      0.1 2.3
    :EndOutflowStageRelation
    :AreaStageRelation LINEAR
      0.33
    :EndAreaStageRelation
  :EndReservoir

  :Reservoir ExampleReservoir
    :SubBasinID 23
    :HRUID 234
    :Type RESROUTE_STANDARD
    :StageRelations
      21 # number of points
      0.09 0 0 0.0  (h [m], Q [m3/s], V [m3], A [m2])
      0.1 2 43 0.2
      ...
      3.0 20000 3500 200
    :EndStageRelations
    :SeepageParameters coeff[m3/s/d] GW_stage [m]
  :EndReservoir

  :Reservoir ExampleReservoir
    :SubBasinID 23
    :HRUID 234
    :Type RESROUTE_TIMEVARYING
    :VaryingStageRelations
      21 # number of points
      [jul day1] [jul day2] [jul day3] ...
      0.09 0 0 0.0 0.0 (h [m], Q [m3/s], A [m2], V [m3])
      0.1 2 43 0.2 0.3
      ...
      3.0 20000 3500 200 25000
    :EndVaryingStageRelations
    :MaxCapacity 25000
  :EndReservoir

  :Reservoir ExampleReservoir # 'lake type format'
    :SubBasinID 23
    :HRUID 234
    :Type RESROUTE_STANDARD
    :WeirCoefficient 0.6
    :CrestWidth 10
    :MaxDepth 5.5 # relative to minimum weir crest elevation
    :LakeArea 10.5
    :AbsoluteCrestHeight 365 # absolute minimum weir crest height m.a.s.l. (optional)
  :EndReservoir


  */
  char *s[MAXINPUTITEMS];
  int    Len;

  long SBID=DOESNT_EXIST;

  double a_V(1000.0),b_V(1.0);
  double a_Q(10.0),b_Q(1.0);
  double a_A(1000.0),b_A(0.0);
  double *aQ(NULL),*aQ_ht(NULL); int NQ(0);
  double *aV(NULL),*aV_ht(NULL); int NV(0);
  double *aA(NULL),*aA_ht(NULL); int NA(0);
  double *aQund(NULL);
  double **aQQ=NULL;
  int     nDates(0);
  int    *aDates(NULL);
  double cwidth(DOESNT_EXIST),max_depth(DOESNT_EXIST);
  double weircoeff(DOESNT_EXIST),lakearea(DOESNT_EXIST),crestht(0.0);
  double max_capacity=0;
  double seep_coeff  =RAV_BLANK_DATA;
  double GW_stage    =0;
  double Smi[12],Sni[12],Sci[12];
  double Qmi[12],Qni[12],Qci[12];
  double Qmax,Smax;
  bool dztr=false;
  bool minstageDom=false;
  double demand_mult=1.0;

  curve_function type;
  type=CURVE_POWERLAW;

  CReservoir *pRes=NULL;
  HRUID=DOESNT_EXIST;

  while(!p->Tokenize(s,Len))
  {
    if(Options.noisy) { cout << "-->reading line " << p->GetLineNumber() << ": "; }
    if(Len == 0) {}//Do nothing
    else if(IsComment(s[0],Len)) { if(Options.noisy) { cout << "#" << endl; } }
    else if(!strcmp(s[0],":SubBasinID"))
    {
      if(Options.noisy) { cout << ":SubBasinID" << endl; }
      SBID = s_to_l(s[1]);
    }
    else if(!strcmp(s[0],":HRUID"))
    {
      if(Options.noisy) { cout << ":HRUID" << endl; }
      HRUID = s_to_i(s[1]);
    }
    else if(!strcmp(s[0],":Type"))
    {
      if(Options.noisy) { cout << ":Type" << endl; }
      //Obsolete - now ignored
    }
    else if(!strcmp(s[0],":CrestWidth"))
    {
      if(Options.noisy) { cout << ":CrestWidth" << endl; }
      cwidth=s_to_d(s[1]);
      type=CURVE_LAKE;
    }
    else if(!strcmp(s[0],":MaxCapacity"))
    {
      if(Options.noisy) { cout << ":MaxCapacity" << endl; }
      max_capacity=s_to_d(s[1]);
    }
    else if(!strcmp(s[0],":WeirCoefficient"))
    {
      if(Options.noisy) { cout << ":WeirCoefficient" << endl; }
      weircoeff=s_to_d(s[1]);
      type=CURVE_LAKE;
    }
    else if(!strcmp(s[0],":MaxDepth"))
    {
      if(Options.noisy) { cout << ":MaxDepth" << endl; }
      max_depth=s_to_d(s[1]);
      type=CURVE_LAKE;
    }
    else if(!strcmp(s[0],":LakeArea"))
    {
      if(Options.noisy) { cout << ":LakeArea" << endl; }
      lakearea=s_to_d(s[1]);
      if (lakearea<100.0){
        string warn=":LakeArea for lake-like reservoir in subbasin "+to_string(SBID)+" seems small. Note LakeArea should be in units of m2, not km2";
        WriteWarning(warn,Options.noisy);
      }
      type=CURVE_LAKE;
    }
    else if(!strcmp(s[0],":AbsoluteCrestHeight"))
    {
      if(Options.noisy) { cout << ":AbsoluteCrestHeight" << endl; }
      crestht=s_to_d(s[1]);
      type=CURVE_LAKE;
    }
    else if(!strcmp(s[0],":SeepageParameters"))
    {
      if(Options.noisy) { cout << ":SeepageParameters" << endl; }
      seep_coeff=s_to_d(s[1]);
      GW_stage  =s_to_d(s[2]);
    }
    else if(!strcmp(s[0],":MinStageConstraintDominant"))
    {
      if(Options.noisy) { cout << ":MinStageConstraintDominant" << endl; }
      minstageDom=true;
    }
    else if(!strcmp(s[0],":DemandMultiplier"))
    {
      if(Options.noisy) { cout << ":DemandMultiplier" << endl; }
      demand_mult=s_to_d(s[1]);
    }
    //----------------------------------------------------------------------------------------------
    else if(!strcmp(s[0],":VolumeStageRelation"))
    {
      if(Options.noisy) { cout << ":VolumeStageRelation" << endl; }
      if(Len >= 2) {
        if(!strcmp(s[1],"POWER_LAW"))
        {
          type = CURVE_POWERLAW;
          p->Tokenize(s,Len);
          if(Len >= 2) { a_V = s_to_d(s[0]); b_V = s_to_d(s[1]); }
          p->Tokenize(s,Len); //:EndVolumeStageRelation
        }
        else if(!strcmp(s[1],"LINEAR"))
        {
          p->Tokenize(s,Len);
          if(Len >= 1) { a_V = s_to_d(s[0]); b_V = 1.0; }
          p->Tokenize(s,Len); //:EndVolumeStageRelation
        }
        else if(!strcmp(s[1],"LOOKUP_TABLE"))
        {
          if(type!=CURVE_LAKE) { type = CURVE_DATA; } //enables :VolumeStageRelation to be used with lake-type 
          p->Tokenize(s,Len);
          if(Len >= 1) { NV = s_to_i(s[0]); }
          aV    = new double[NV];
          aV_ht = new double[NV];
          for(int i = 0; i < NV; i++) {
            p->Tokenize(s,Len);
            if(IsComment(s[0],Len)) { i--; }
            else {
              aV_ht[i] = s_to_d(s[0]);
              aV[i] = s_to_d(s[1]);
            }
          }
          p->Tokenize(s,Len); //:EndVolumeStageRelation
        }
      }
      else {
        //write warning
      }
    }
    //----------------------------------------------------------------------------------------------
    else if(!strcmp(s[0],":AreaStageRelation"))
    {
      if(Options.noisy) { cout << ":AreaStageRelation" << endl; }
      if(Len >= 2) {
        if(!strcmp(s[1],"POWER_LAW"))
        {
          type = CURVE_POWERLAW;
          p->Tokenize(s,Len);
          if(Len >= 2) { a_A = s_to_d(s[0]); b_A = s_to_d(s[1]); }
          p->Tokenize(s,Len); //:EndAreaStageRelation
        }
        else if(!strcmp(s[1],"CONSTANT"))
        {
          type = CURVE_LINEAR;
          p->Tokenize(s,Len);
          if(Len >= 1) { a_A = s_to_d(s[0]); b_A = 0.0; }
          p->Tokenize(s,Len); //:EndAreaStageRelation
        }
        else if(!strcmp(s[1],"LINEAR"))
        {
          type = CURVE_LINEAR;
          p->Tokenize(s,Len);
          if(Len >= 1) { a_A = s_to_d(s[0]); b_A = 1.0; }
          p->Tokenize(s,Len); //:EndAreaStageRelation
        }
        else if(!strcmp(s[1],"LOOKUP_TABLE"))
        {
//          type = CURVE_DATA;
          p->Tokenize(s,Len);
          if(Len >= 1) { NA = s_to_i(s[0]); }
          aA = new double[NA];
          aA_ht = new double[NA];
          for(int i = 0; i < NA; i++) {
            p->Tokenize(s,Len);
            if(IsComment(s[0],Len)) { i--; }
            else {
              aA_ht[i] = s_to_d(s[0]);
              aA[i] = s_to_d(s[1]);
            }
          }
          p->Tokenize(s,Len); //:EndAreaStageRelation
        }
      }
      else {
        //write warning
      }
    }
    //----------------------------------------------------------------------------------------------
    else if(!strcmp(s[0],":OutflowStageRelation"))
    {
      if(Options.noisy) { cout << ":OutflowStageRelation" << endl; }
      if(Len >= 2) {
        if(!strcmp(s[1],"POWER_LAW"))
        {
          type = CURVE_POWERLAW;
          p->Tokenize(s,Len);
          if(Len >= 2) { a_Q = s_to_d(s[0]); b_Q = s_to_d(s[1]); }
          p->Tokenize(s,Len); //:EndOutflowStageRelation
        }
        else if(!strcmp(s[1],"LINEAR"))
        {
          type = CURVE_LINEAR;
          p->Tokenize(s,Len);
          if(Len >= 1) { a_Q = s_to_d(s[0]); b_Q = 1.0; }
          p->Tokenize(s,Len); //:EndOutflowStageRelation
        }
        else if(!strcmp(s[0],"LOOKUP_TABLE"))
        {
          type = CURVE_DATA;
          p->Tokenize(s,Len);
          if(Len >= 1) { NQ = s_to_i(s[0]); }
          aQ = new double[NQ];
          aQ_ht = new double[NQ];
          for(int i = 0; i < NQ; i++) {
            p->Tokenize(s,Len);
            if(IsComment(s[0],Len)) { i--; }
            else {
              aQ_ht[i] = s_to_d(s[0]);
              aQ[i] = s_to_d(s[1]);
            }
          }
          p->Tokenize(s,Len); //:EndOutflowStageRelation
        }
      }
      else {
        //write warning
      }
    }
    //----------------------------------------------------------------------------------------------
    else if(!strcmp(s[0],":StageRelations"))
    {
      if(Options.noisy) { cout << ":StageRelations" << endl; }
      type = CURVE_DATA;
      p->Tokenize(s,Len);
      if(Len >= 1) { NQ = s_to_i(s[0]); }

      aQ_ht = new double[NQ];
      aQ = new double[NQ];
      aV = new double[NQ];
      aA = new double[NQ];
      aQund=new double[NQ];
      for(int i = 0; i < NQ; i++) {
        p->Tokenize(s,Len);
        if(IsComment(s[0],Len)) { i--; }
        else {
          if(Len>=4) {
            aQ_ht[i] = s_to_d(s[0]);
            aQ[i] = s_to_d(s[1]);
            aV[i] = s_to_d(s[2]);
            aA[i] = s_to_d(s[3]);
            aQund[i]=0.0;
            if(Len>=5) {
              aQund[i]=s_to_d(s[4]);
            }
          }
          else {
            WriteWarning("Incorrect line length (<4) in :Reservoir :StageRelations command",Options.noisy);
          }
        }
      }
      p->Tokenize(s,Len); //:EndStageRelations
    }
    //----------------------------------------------------------------------------------------------
    else if(!strcmp(s[0],":EndStageRelations"))
    {
    }
    //----------------------------------------------------------------------------------------------
    else if(!strcmp(s[0],":VaryingStageRelations"))
    {
      bool hasQund=false;
      int shift=0;
      if(Options.noisy) { cout << ":VaryingStageRelations" << endl; }
      type = CURVE_VARYING;
      p->Tokenize(s,Len);
      if(Len >= 1) { NQ = s_to_i(s[0]); } //# of flows
      if(Len>=2) { hasQund=true;shift=1; } //TMP DEBUG - shoudld have flag for extra column
      p->Tokenize(s,Len);
      nDates = Len;

      aDates = new int[nDates];
      for(int v = 0; v < nDates; v++) { aDates[v] = s_to_i(s[v]); }

      aQ_ht = new double[NQ];
      aQQ = new double *[nDates];
      for(int v = 0; v < nDates; v++) {
        aQQ[v] = new double[NQ];
      }
      aV = new double[NQ];
      aA = new double[NQ];
      aQund=new double[NQ];
      for(int i = 0; i < NQ; i++) {
        p->Tokenize(s,Len);
        if(IsComment(s[0],Len)) { i--; }
        else {
          if(Len==1) {//presumably :EndVaryingStageRelations
            ExitGracefully("CReservoir::Parse: improper number of rows in :VaryingStageRelations command",BAD_DATA);
          }
          else if(Len < 2 + nDates) {
            ExitGracefully("CReservoir::Parse: improper number of columns in :VaryingStageRelations command",BAD_DATA);
          }
          aQ_ht[i] = s_to_d(s[0]); //stage
          aV[i] = s_to_d(s[1]); //volume
          aA[i] = s_to_d(s[2]); //area
          if(hasQund) { aQund[i] = s_to_d(s[3]); } //Qund
          else { aQund[i]=0.0; }
          for(int v = 0; v < nDates; v++) {
            aQQ[v][i] = s_to_d(s[3 + v +shift]); //flows for each date
          }
        }
      }
      p->Tokenize(s,Len); //:EndVaryingStageRelations
    }

    else if(!strcmp(s[0],":DZTRResservoirModel"))
    {/*:DZTRResservoirModel
       :MaximumStorage           Vmax(m3)
       :MaximumChannelDischarge  Qmax(m3/s)
       :MonthlyMaxStorage        J F M A M J J A S N D  {or single constant value} (m3)
       :MonthlyNormalStorage     J F M A M J J A S N D  {or constant value} (m3)
       :MonthlyCriticalStorage   J F M A M J J A S N D  {or constant value} (m3)
       :MonthlyMaxDischarge      J F M A M J J A S N D  {or constant value} (m3/s)
       :MonthlyNormalDischarge   J F M A M J J A S N D  {or constant value} (m3/s)
       :MonthlyCriticalDischarge J F M A M J J A S N D  {or constant value} (m3/s)
     :EndDZTRReservoirModel
     */
     //Dynamically zoned target release model
      if(Options.noisy) { cout <<"DZTR Reservoir Model"<<endl; }

      dztr=true;
      
      p->Tokenize(s,Len);
      Smax=s_to_d(s[1]);

      p->Tokenize(s,Len);
      Qmax=s_to_d(s[1]);

      p->Tokenize(s,Len);
      if    (Len>=13) { for(int i=0;i<12;i++) { Smi[i]=s_to_d(s[i+1]); } }
      else if(Len>=2) { for(int i=0;i<12;i++) { Smi[i]=s_to_d(s[1]); } }

      p->Tokenize(s,Len);
      if    (Len>=13) { for(int i=0;i<12;i++) { Sni[i]=s_to_d(s[i+1]); } }
      else if(Len>=2) { for(int i=0;i<12;i++) { Sni[i]=s_to_d(s[1]); } }

      p->Tokenize(s,Len);
      if    (Len>=13) { for(int i=0;i<12;i++) { Sci[i]=s_to_d(s[i+1]); } }
      else if(Len>=2) { for(int i=0;i<12;i++) { Sci[i]=s_to_d(s[1]); } }

      p->Tokenize(s,Len);
      if    (Len>=13) { for(int i=0;i<12;i++) { Qmi[i]=s_to_d(s[i+1]); } }
      else if(Len>=2) { for(int i=0;i<12;i++) { Qmi[i]=s_to_d(s[1]); } }

      p->Tokenize(s,Len);
      if    (Len>=13) { for(int i=0;i<12;i++) { Qni[i]=s_to_d(s[i+1]); } }
      else if(Len>=2) { for(int i=0;i<12;i++) { Qni[i]=s_to_d(s[1]); } }

      p->Tokenize(s,Len);
      if    (Len>=13) { for(int i=0;i<12;i++) { Qci[i]=s_to_d(s[i+1]); } }
      else if(Len>=2) { for(int i=0;i<12;i++) { Qci[i]=s_to_d(s[1]); } }
      
      p->Tokenize(s,Len); //:EndDZTR

      break;
    }

    //----------------------------------------------------------------------------------------------
    else if(!strcmp(s[0],":EndReservoir")) {
      if(Options.noisy) { cout << ":EndReservoir" << endl; }
      break;
    }
    else {
      WriteWarning("Reservoir::Parse: unrecognized command (" + to_string(s[0]) + ") in :Reservoir-:EndReservoir Block",Options.noisy);
    }
  }

  if((type==CURVE_POWERLAW) || (type==CURVE_LINEAR))
  {
    pRes=new CReservoir(name,SBID,a_V,b_V,a_Q,b_Q,a_A,b_A);
  }
  else if(type==CURVE_DATA)
  {
    pRes=new CReservoir(name,SBID,aQ_ht,aQ,aQund,aA,aV,NQ);//presumes aQ_ht=aV_ht=aA_ht; NA=NV=NQ
  }
  else if(type==CURVE_VARYING)
  {
    pRes=new CReservoir(name,SBID,nDates,aDates,aQ_ht,aQQ,aQund,aA,aV,NQ);//presumes aQ_ht=aV_ht=aA_ht; NA=NV=NQ
  }
  else if(type==CURVE_LAKE)
  {
    ExitGracefullyIf(weircoeff==DOESNT_EXIST,"CReservoir::Parse: :WeirCoefficient must be specified for lake-type reservoirs",BAD_DATA_WARN);
    ExitGracefullyIf(cwidth   ==DOESNT_EXIST,"CReservoir::Parse: :CrestWidth must be specified for lake-type reservoirs",BAD_DATA_WARN);
    ExitGracefullyIf(lakearea ==DOESNT_EXIST,"CReservoir::Parse: :LakeArea  must be specified for lake-type reservoirs",BAD_DATA_WARN);
    ExitGracefullyIf(max_depth==DOESNT_EXIST,"CReservoir::Parse: :LakeDepth must be specified for lake-type reservoirs",BAD_DATA_WARN);

    pRes=new CReservoir(name,SBID,weircoeff,cwidth,crestht,lakearea,max_depth);
  }
  else {
    ExitGracefully("CReservoir::Parse: only currently supporting linear, powerlaw, or data reservoir rules",STUB);
  }

  if(max_capacity!=0) 
  { 
    pRes->SetMaxCapacity(max_capacity); //overrides estimates (which will typically be too large)
  }
  if(seep_coeff  !=RAV_BLANK_DATA) 
  { 
    pRes->SetGWParameters(seep_coeff,GW_stage); 
  }
  if((type==CURVE_LAKE) && (aV!=NULL) && (aV_ht!=NULL)) 
  {
    pRes->SetVolumeStageCurve(aV_ht,aV,NV);//allows user to override prismatic lake assumption
  }
  if((type==CURVE_LAKE) && (aA!=NULL) && (aA_ht!=NULL))
  {
    pRes->SetAreaStageCurve(aA_ht,aA,NA);//allows user to override prismatic lake assumption
  }
  if(dztr)
  {
    pRes->SetDZTRModel(Qmax,Smax,Sci,Sni,Smi,Qci,Qni,Qmi);
  }
  if(minstageDom) 
  {
    pRes->SetMinStageDominant();
  }
  if(demand_mult!=1.0) {
    pRes->SetDemandMultiplier(demand_mult);
  }
  
  for(int i = 0; i < nDates; i++) { delete[] aQQ[i]; }delete[] aQQ;
  delete[] aQ;
  delete[] aQ_ht;
  delete[] aQund;
  delete[] aV;
  delete[] aV_ht;
  delete[] aA;
  delete[] aA_ht;
  return pRes;
}