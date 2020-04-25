/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2020 the Raven Development Team
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
    if     (Len==0)                                { code=-1; }//blank line
    else if(IsComment(s[0],Len))                   { code=-2; }//comment

    else if(!strcmp(s[0],":Echo"                )) { code=-3; }
   //-------------------MODEL ENSEMBLE PARAMETERS----------------
    else if(!strcmp(s[0],":VegetationChange"    )) { code=1; }
    else if(!strcmp(s[0],":LandUseChange"       )) { code=2; }
    else if(!strcmp(s[0],":GroupLandUseChange"       )) { code=3; }
    else if(!strcmp(s[0],":SetStreamflow"       )) { code=10; }
    else if(!strcmp(s[0],":SetReservoirFlow"    )) { code=11; }
    else if(!strcmp(s[0],":SetReservoirStage"   )) { code=12; }
//    else if(!strcmp(s[0],":SetExtractionRate"   )) { code=13; }//set streamflow extraction
//    else if(!strcmp(s[0],":DistributedIrrigation"   )) { code=14; }//add irrigation water (mm/d) to landscape

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
    case(10):  //----------------------------------------------
    { /*:SetStreamflow [SBID] [value]*/
      pSB=pModel->GetSubBasinByID(s_to_l(s[1]));
      double Q=fast_s_to_d(s[2]);
      double *aQ=new double[1];
      aQ[0]=Q;
      pSB->SetQoutArray(DOESNT_EXIST,aQ,Q);
      delete [] aQ;
      //or pSB->SetCurrentOutflow(Q);
      break;
    }
    case(11):  //----------------------------------------------
    { /*:SetReservoirStage [SBID] [value]*/
      pSB=pModel->GetSubBasinByID(s_to_l(s[1]));
      pSB->SetInitialReservoirStage(s_to_d(s[2]),s_to_d(s[2]));
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