/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright © 2008-2014 the Raven Development Team
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Model.h"
#include "StateVariables.h"
#include "HydroUnits.h"
#include "ParseLib.h"

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
bool ParseHRUPropsFile(CModel *&pModel, const optStruct &Options)
{
  int         i;                //counters
  long        SBID;             //subbasin ID
  CHydroUnit *pHRU;             //temp pointers
  CSubBasin  *pSB;
  bool        ended=false;

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
    else if  (!strcmp(s[0],":RedirectToFile"           )){code=-3; }//redirect to secondary file
    else if  (!strcmp(s[0],":End"                      )){code=-4; }//stop reading
    //--------------------MODEL OPTIONS ------------------------
    else if  (!strcmp(s[0],":SubBasins"                )){code=1;  }
    else if  (!strcmp(s[0],":HRUs"                     )){code=2;  }
    else if  (!strcmp(s[0],":Reservoir"                )){code=3;  }
    else if  (!strcmp(s[0],":SubBasinProperties"       )){code=7;  }
    else if  (!strcmp(s[0],":HRUGroup"                 )){code=8;  }
    
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
        
        string filedir = GetDirectoryName(Options.rvt_filename); //if a relative path name, e.g., "/path/model.rvt", only returns e.g., "/path"
        if (StringToUppercase(filename).find(StringToUppercase(filedir)) == string::npos){ //checks to see if absolute dir already included in redirect filename
          filename = filedir + "\\" + filename;
        }

        INPUT2.open(filename.c_str());   
        if (INPUT2.fail()){
          ostrstream FILENAME;
          FILENAME<<":RedirectToFile: Cannot find file "<<filename<<ends;
          ExitGracefully(FILENAME.str() ,BAD_DATA); 
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
              string error;
              error="Parse HRU File: Unrecognized Channel profile code ("+string(s[3])+") in SubBasins command";
              ExitGracefullyIf((pChan==NULL) && (Options.routing!=ROUTE_NONE),error.c_str(),BAD_DATA);
              
              double length;
              length=AutoOrDouble(s[4]);
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

            void *pAqStack=NULL;
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
            pHRU=new CHydroUnit( pModel,
                                 s_to_i(s[0]),//ID
                                 pModel->GetNumHRUs(),//k - global model index
                                 s_to_d(s[1]),//area
                                 pModel->GetSubBasinIndex(SBID),
                                 s_to_d(s[2]),//elev
                                 s_to_d(s[3]),//lat 
                                 s_to_d(s[4]),//long
                                 s_to_d(s[11])*PI/180.0,//slope (deg->radians)
                                 s_to_d(s[12])*PI/180.0,//aspect (deg->radians)
                                 HRUtype,
                                 pSoilProfile,
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
        */
        if (Options.noisy) {cout <<":Reservoir"<<endl;}  
        CReservoir *pRes;
        int HRUID;
        pRes=CReservoir::Parse(pp,s[1],HRUID,Options);
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
         {ID,param1, param2,...,paramN} x nSubBasins
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
          else if (!strcmp(s[0],":EndSubBasinProperties")){}//done
          else
          {
            ExitGracefullyIf(Len<nParamStrings,
              "Parse HRU File: incorrect number of terms in SubBasin properties",BAD_DATA);
            
            SBID=s_to_l(s[0]);
            
            pSB=pModel->GetSubBasinByID(SBID);
            if (pSB!=NULL){
              for (i=1;i<nParamStrings;i++){
                good_string=pSB->SetBasinProperties(aParamStrings[i],AutoOrDouble(s[i]));
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
          WriteWarning("HRU groups should ideally be defined in .rvi file (using :DefineHRUGroup(s) commands) before being populated in .rvh file",Options.noisy);
          pHRUGrp=new CHRUGroup(s[1],pModel->GetNumStateVars());
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
                for (k=0;k<pModel->GetNumHRUs();k++)
                {
                  if (pModel->GetHydroUnit(k)->GetID()==ii)
                  {
                    pHRUGrp->AddHRU(pModel->GetHydroUnit(k));found=true;
                  }
                }
                if (!found){gaps=true;}
              }
              if (gaps){
                  WriteWarning("Range specified in HRUGroup command has gaps",Options.noisy);
               }
            }
          }
        }
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
  if ((pModel->GetNumSubBasins() > 1) && (Options.routing == ROUTE_NONE))
  {
    string warn;
    warn = "ParseHRUPropsFile: ROUTE_NONE is being used for a model with more than one basin. Typically this method is used for lumped models only.";
    WriteWarning(warn,Options.noisy);
  }

  delete pp; 
  pp=NULL;
  return true;
}