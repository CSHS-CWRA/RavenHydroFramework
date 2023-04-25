/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2023 the Raven Development Team
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Model.h"
#include "StateVariables.h"
#include "HydroUnits.h"
#include "ParseLib.h"

//////////////////////////////////////////////////////////////////
/// \brief Parses Live Communications File
/// \details model.rvl: input file that is read every N time steps 
///
/// \param *&pModel [out] Reference to model object
/// \param &Options [out] Global model options information
/// \return True if operation is successful
//
void ParseLiveFile(CModel *&pModel,const optStruct &Options, const time_struct &tt)
{
  if (Options.rvl_read_frequency==0.0){return;} //does nothing

  //if not evenly divided by frequency, return
  if(fabs(ffmod(tt.model_time,Options.rvl_read_frequency)) > 0.5*Options.timestep){return;}

  bool        ended(false);
  CHydroUnit *pHRU;
  CSubBasin  *pSB;

  ifstream    RVL;
  RVL.open(Options.rvl_filename.c_str());
  if(RVL.fail()) {
    string warn="ERROR opening model live file: "+Options.rvl_filename;
    ExitGracefully(warn.c_str(), BAD_DATA);return;
  }

  int   Len,line(0),code;
  char *s[MAXINPUTITEMS];
  CParser *pp=new CParser(RVL,Options.rvl_filename,line);

  //--Sift through file-----------------------------------------------
  bool end_of_file=pp->Tokenize(s,Len);
  while(!end_of_file)
  {
    if(ended) { break; }

    code=0;
    //---------------------SPECIAL -----------------------------
    if     (Len==0)                                 { code=-1; }//blank line
    else if(IsComment(s[0],Len))                    { code=-2; }//comment

    else if(!strcmp(s[0],":Echo"                 )) { code=-3; }
   //-------------------MODEL ENSEMBLE PARAMETERS----------------
    else if(!strcmp(s[0],":VegetationChange"     )) { code=1;  }
    else if(!strcmp(s[0],":LandUseChange"        )) { code=2;  }
    else if(!strcmp(s[0],":GroupLandUseChange"   )) { code=3;  }
    else if(!strcmp(s[0],":GroupVegetationChange")) { code=4;  }
    else if(!strcmp(s[0],":HRUTypeChange"        )) { code=5;  }
    else if(!strcmp(s[0],":GroupHRUTypeChange"   )) { code=6;  }
    else if(!strcmp(s[0],":SetStreamflow"        )) { code=10; }
    else if(!strcmp(s[0],":SetReservoirFlow"     )) { code=11; }
    else if(!strcmp(s[0],":SetReservoirStage"    )) { code=12; }
//    else if(!strcmp(s[0],":SetExtractionRate"     )) { code=13; }//set streamflow extraction
//    else if(!strcmp(s[0],":DistributedIrrigation" )) { code=14; }//add irrigation water (mm/d) to landscape
//    else if(!strcmp(s[0],":FlowDiversion"         )) { code=15; }//diverts flow from one basin to another
//    else if(!strcmp(s[0],":IrrigationDemand"      )) { code=16; }//removes flow from basin
//    else if(!strcmp(s[0],":UpdateStateVariable"   )) { code=17; }
//    else if(!strcmp(s[0],":UpdateParameter"       )) { code=18; }
    else if(!strcmp(s[0],":RepopulateHRUGroup"  )) { code=20; }
    switch(code)
    {
    case(-1):  //----------------------------------------------
    {/*Blank Line*/
      break;
    }
    case(-2):  //----------------------------------------------
    {/*Comment*/
      break;
    }
    case(-3):  //----------------------------------------------
    {/*:Echo [statement]*/
      for(int i=1;i<Len;i++) {
        cout<<s[i]<<" ";
      }
      break;
    }
    case(1):  //----------------------------------------------
    {/*:VegetationChange [HRU ID] [new vegetation tag] */
      pHRU=pModel->GetHRUByID(s_to_i(s[1]));
      if(pHRU!=NULL) {
        CVegetationClass *veg_class= CVegetationClass::StringToVegClass(s[2]);
        if(veg_class!=NULL) {
          pHRU->ChangeVegetation(veg_class);
        }
        else {
          WriteWarning("ParseLiveFile: invalid vegetation tag provided in :LandUseChange command",Options.noisy);
        }
      }
      else {
        WriteWarning("ParseLiveFile: invalid HRU ID provided in :VegetationChange command",Options.noisy);
      }
      break;
    }
    case(2):  //----------------------------------------------
    {/*:LandUseChange [HRU ID] [new LULT tag]*/
      pHRU=pModel->GetHRUByID(s_to_i(s[1]));
      if(pHRU!=NULL) {
        CLandUseClass *lult_class= CLandUseClass::StringToLUClass(s[2]);
        if(lult_class!=NULL) {
          pHRU->ChangeLandUse(lult_class);
        }
        else {
          WriteWarning("ParseLiveFile: invalid Land use tag provided in :LandUseChange command",Options.noisy);
        }
      }
      else {
        WriteWarning("ParseLiveFile: invalid HRU ID provided in :LandUseChange command",Options.noisy);
      }
      break;
    }
    case(3):  //----------------------------------------------
    {/*:GroupLandUseChange [HRU group] [new LULT tag]*/
      CHRUGroup *pHRUGroup=pModel->GetHRUGroup(s[1]);
      if(pHRUGroup!=NULL) {
        CLandUseClass *lult_class= CLandUseClass::StringToLUClass(s[2]);
        if(lult_class!=NULL) {
          for(int k=0;k<pHRUGroup->GetNumHRUs();k++) {
            pHRUGroup->GetHRU(k)->ChangeLandUse(lult_class);
          }
        }
        else {
          WriteWarning("ParseLiveFile: invalid Land use tag provided in :GroupLandUseChange command",Options.noisy);
        }
      }
      else {
        WriteWarning("ParseLiveFile: invalid HRU Group name provided in :GroupLandUseChange command",Options.noisy);
      }
      break;
    }
    case(4):  //----------------------------------------------
    {/*:GroupVegetationChange [HRU group] [new veg tag]*/
      CHRUGroup *pHRUGroup=pModel->GetHRUGroup(s[1]);
      if(pHRUGroup!=NULL) {
        CVegetationClass *veg_class=CVegetationClass::StringToVegClass(s[2]);
        if(veg_class!=NULL) {
          for(int k=0;k<pHRUGroup->GetNumHRUs();k++) {
            pHRUGroup->GetHRU(k)->ChangeVegetation(veg_class);
          }
        }
        else {
          WriteWarning("ParseLiveFile: invalid Land use tag provided in :GroupLandUseChange command",Options.noisy);
        }
      }
      else {
        WriteWarning("ParseLiveFile: invalid HRU Group name provided in :GroupLandUseChange command",Options.noisy);
      }
      break;
    }
    case(5):  //----------------------------------------------
    {/*:HRUTypeChange [HRU ID] [new HRU Type]*/
      pHRU=pModel->GetHRUByID(s_to_i(s[1]));
      if(pHRU!=NULL) {
        HRU_type hru_type=StringToHRUType(s[2]);
        if(hru_type!=HRU_INVALID_TYPE) {
          pHRU->ChangeHRUType(hru_type);
        }
        else {
          WriteWarning("ParseLiveFile: invalid HRU type tag provided in :HRUTypeChange command",Options.noisy);
        }
      }
      else {
        WriteWarning("ParseLiveFile: invalid HRU ID provided in :HRUTypeChange command",Options.noisy);
      }
      break;
    }
    case(6):  //----------------------------------------------
    {/*:GroupHRUTypeChange [HRU Group] [new HRU Type]*/
      CHRUGroup *pHRUGroup=pModel->GetHRUGroup(s[1]);
      if(pHRUGroup!=NULL) {
        HRU_type hru_type=StringToHRUType(s[2]);
        if(hru_type!=HRU_INVALID_TYPE) {
          for(int k=0;k<pHRUGroup->GetNumHRUs();k++) {
            pHRUGroup->GetHRU(k)->ChangeHRUType(hru_type);
          }
        }
        else {
          WriteWarning("ParseLiveFile: invalid HRU type tag provided in :HRUTypeChange command",Options.noisy);
        }
      }
      else {
        WriteWarning("ParseLiveFile: invalid HRU ID provided in :HRUTypeChange command",Options.noisy);
      }
      break;
    }
    case(10):  //----------------------------------------------
    { /*:SetStreamflow [SBID] [value]*/
      pSB=pModel->GetSubBasinByID(s_to_l(s[1]));
      double Q=fast_s_to_d(s[2]);
      pSB->SetQout(Q);
      //or pSB->SetCurrentOutflow(Q);
      break;
    }
    case(11):  //----------------------------------------------
    { /*:SetReservoirStage [SBID] [value]*/
      pSB=pModel->GetSubBasinByID(s_to_l(s[1]));
      pSB->GetReservoir()->SetReservoirStage(s_to_d(s[2]),s_to_d(s[2]));
      break;
    }
    case(12):  //----------------------------------------------
    { /*:SetReservoirFlow [SBID] [value]*/
      pSB=pModel->GetSubBasinByID(s_to_l(s[1]));
      //pSB->GetReservoir()->OverrideFlow(s_to_d(s[2]));
      //JRC: add CReservoir member _overrideQ which takes priority over override Q time series 
      ExitGracefully("ParseLiveFile:SetReservoirFlow",STUB);
      break;
    }
    case(20):  //----------------------------------------------
    { /*:RepopulateHRUGroup [HRUGroup]
      1, 8, 19, ..., 34
      :EndRepopulateHRUGroup
      */
      if (Options.noisy) {cout <<"  Repopulate HRU Group..."<<endl;}
      if (Len!=2){pp->ImproperFormat(s);}
      CHRUGroup *pHRUGrp=NULL;
      pHRUGrp=pModel->GetHRUGroup(s[1]);
      if (pHRUGrp==NULL){//group not yet defined
        ExitGracefully("ParseLiveFile: :RepopulateHRUGroup command: HRU Groups cannot be repopulated if they haven't been defined",BAD_DATA);
      }
      pHRUGrp->EmptyGroup();
      while ((Len==0) || (strcmp(s[0],":EndRepopulateHRUGroup")))
      {
        pp->Tokenize(s,Len);
        if      (IsComment(s[0], Len)){}//comment line
        else if (!strcmp(s[0],":EndRepopulateHRUGroup")){}//done
        else
        {
          int k;
          int nHRUs=pModel->GetNumHRUs();
          for (int i=0;i<Len;i++)
          {
            int ind1,ind2;
            s_to_range(s[i],ind1,ind2);
            bool found=false;
            ExitGracefullyIf((ind2-ind1)>10000,"Parsing :RepopulateHRUGroup command: invalid range of HRU indices",BAD_DATA);
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
              if((!found) && (ind1==ind2)) {
                string warn="ParseLiveFile: HRU ID "+to_string(ii)+" specified in :RepopulateHRUGroup command for group "+pHRUGrp->GetName()+ " does not exist";
                ExitGracefully(warn.c_str(),BAD_DATA);
              }
            }
          }
        }
      }
      break;
    }
    default://------------------------------------------------
    {
      string warning="Unrecognized command in live file at line "+to_string(Len);
      WriteWarning(warning,Options.noisy);
      break;
    }
    }//end switch(code)

    end_of_file=pp->Tokenize(s,Len);

  } //end while !end_of_file
  RVL.close();

  delete pp;
  pp=NULL;
}