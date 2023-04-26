/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2023 the Raven Development Team
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Model.h"
#include "StateVariables.h"
#include "HydroUnits.h"
#include "ParseLib.h"
#include "ModelEnsemble.h"
#include "EnKF.h"
bool IsContinuousFlowObs2(const CTimeSeriesABC* pObs,long SBID);
//////////////////////////////////////////////////////////////////
/// \brief Parses Ensemble Model file
/// \details model.rve: input file that defines ensemble member details for MC, calibration, etc.
///
/// \param *&pModel [out] Reference to model object
/// \param &Options [out] Global model options information
/// \return True if operation is successful
//
bool ParseEnsembleFile(CModel *&pModel,const optStruct &Options)
{
  int         i;              //counters
  bool        ended(false);
  bool        in_ifmode_statement=false;

  ifstream    ENS;
  ENS.open(Options.rve_filename.c_str());
  if(ENS.fail()) {
    cout << "ERROR opening model ensemble file: "<<Options.rve_filename<<endl; return false;
  }

  int   Len,line(0),code;
  char *s[MAXINPUTITEMS];
  CParser *pp=new CParser(ENS,Options.rve_filename,line);

  ifstream INPUT2;           //For Secondary input
  CParser *pMainParser=NULL; //for storage of main parser while reading secondary files

  CEnsemble *pEnsemble=pModel->GetEnsemble();

  if(Options.noisy) {
    cout <<"======================================================"<<endl;
    cout <<"Parsing Model Ensemble File " << Options.rve_filename <<"..."<<endl;
    cout <<"======================================================"<<endl;
  }
  
  //--Sift through file-----------------------------------------------
  bool end_of_file=pp->Tokenize(s,Len);
  while(!end_of_file)
  {
    if(ended) { break; }
    if(Options.noisy) { cout << "reading line " << pp->GetLineNumber() << ": "; }

    /*assign code for switch statement
    ------------------------------------------------------------------
    <0           : ignored/special
    0   thru 100 : All other
    ------------------------------------------------------------------
    */

    code=0;
    //---------------------SPECIAL -----------------------------
    if     (Len==0)                                       { code=-1; }//blank line
    else if(IsComment(s[0],Len))                          { code=-2; }//comment
    else if(!strcmp(s[0],":RedirectToFile"))              { code=-3; }//redirect to secondary file
    else if(!strcmp(s[0],":End"))                         { code=-4; }//stop reading
    else if(!strcmp(s[0],":IfModeEquals"                )){ code=-5; }
    else if(in_ifmode_statement)                          { code=-6; }
    else if(!strcmp(s[0],":EndIfModeEquals"             )){ code=-2; }//treat as comment - unused mode 

    //-------------------MODEL ENSEMBLE PARAMETERS----------------
    else if(!strcmp(s[0],":OutputDirectoryFormat"))       { code=1;  }
    else if(!strcmp(s[0],":RunNameFormat"))               { code=2;  }
    else if(!strcmp(s[0],":EnsembleRVCFormat"))           { code=3;  }
    else if(!strcmp(s[0],":ParameterDistributions"))      { code=10; }
    else if(!strcmp(s[0],":ObjectiveFunction"))           { code=11; }
    else if(!strcmp(s[0],":ForcingPerturbation"))         { code=12; }
    else if(!strcmp(s[0],":AssimilatedState"))            { code=13; }
    else if(!strcmp(s[0],":SolutionRunName"))             { code=14; }
    else if(!strcmp(s[0],":WindowSize"))                  { code=15; } 
    else if(!strcmp(s[0],":ObservationErrorModel"))       { code=16; }
    else if(!strcmp(s[0],":EnKFMode"))                    { code=18; }
    else if(!strcmp(s[0],":ExtraRVTFilename"))            { code=19; }
    else if(!strcmp(s[0],":AssimilateStreamflow"))        { code=101;}

    switch(code)
    {
    case(-1):  //----------------------------------------------
    {/*Blank Line*/
      if(Options.noisy) { cout <<""<<endl; }break;
    }
    case(-2):  //----------------------------------------------
    {/*Comment*/
      if(Options.noisy) { cout <<"*"<<endl; } break;
    }
    case(-3):  //----------------------------------------------
    {/*:RedirectToFile*/
      string filename="";
      for(int i=1;i<Len;i++) { filename+=s[i]; if(i<Len-1) { filename+=' '; } }
      if(Options.noisy) { cout <<"Redirect to file: "<<filename<<endl; }

      filename=CorrectForRelativePath(filename,Options.rve_filename);

      INPUT2.open(filename.c_str());
      if(INPUT2.fail()) {
        string warn;
        warn=":RedirectToFile (from .rve): Cannot find file "+filename;
        ExitGracefully(warn.c_str(),BAD_DATA);
      }
      else {
        if (pMainParser != NULL) {
          ExitGracefully("ParseEnsembleFile::nested :RedirectToFile commands (in already redirected files) are not allowed.",BAD_DATA);
        }
        pMainParser=pp;   //save pointer to primary parser
        pp=new CParser(INPUT2,filename,line);//open new parser
      }

      break;
    }
    case(-4):  //----------------------------------------------
    {/*:End*/
      if(Options.noisy) { cout <<"EOF"<<endl; } ended=true; 
      break;
    }
    case(-5):  //----------------------------------------------
    {/*:IfModeEquals [mode] {mode2} {mode3}*/
      if(Len>1) {
        if(Options.noisy) { cout <<"Mode statement start..."<<endl; }
        bool mode_match=false;
        for(int i=1; i<Len; i++) {
          if(s[i][0]==Options.run_mode) { mode_match=true; }
        }
        if(!mode_match) { in_ifmode_statement=true; }
      }
      break;
    }
    case(-6):  //----------------------------------------------
    {/*in_ifmode_statement*/
      if(!strcmp(s[0],":EndIfModeEquals"))
      {
        if(Options.noisy) { cout <<"...Mode statement end"<<endl; }
        in_ifmode_statement=false;
      }
      break;
    }
    case(1):  //----------------------------------------------
    {/*:OutputDirectoryFormat [dirname with * for ensemble ID]*/
      if(Options.noisy) { cout <<"OutputDirectoryFormat"<<endl; } 
      string outdir="";
      for (i=1;i<Len;i++){outdir=outdir+to_string(s[i]); }
      //outdir=CorrectForRelativePath(outdir,Options.rvi_filename); (not for directory?)
      pEnsemble->SetOutputDirectory(outdir);
      break;
    }
    case(2):  //----------------------------------------------
    {/*:RunNameFormat [runn name, with * for ensemble number locale] */
      if(Options.noisy) { cout <<"RunNameFormat"<<endl; } 
      pEnsemble->SetRunNames(s[1]); //No spaces!
      break;
    }
    case(3):  //----------------------------------------------
    {/*:EnsembleRVCFormat [solution file names, with * for ensemble number locale] */
      if(Options.noisy) { cout <<"Ensemble RVC Format"<<endl; } 
      pEnsemble->SetSolutionFiles(s[1]); //No spaces!
      break;
    }
    case(10):  //----------------------------------------------
    { /*
      ":ParameterDistributions"
      {param_name class_type class name base_value dist_type dist_param1 dist_param2... } x nDistributions
      :EndParameterDistributions
      */
      if(Options.noisy) { cout <<"Parameter Distributions..."<<endl; }
      while((Len==0) || (strcmp(s[0],":EndParameterDistributions")))
      {
        pp->Tokenize(s,Len);
        if     (IsComment(s[0],Len))         {}//comment line
        else if(!strcmp(s[0],":Attributes")) {}//ignored by Raven - needed for GUIs
        else if(!strcmp(s[0],":Units"     )) {}//ignored by Raven - needed for GUIs
        else if(!strcmp(s[0],":EndParameterDistributions")) {}//done
        else
        {
          ExitGracefullyIf(Len<7,
            "Parse Ensemble File: incorrect number of terms in ParameterDistributions command.",BAD_DATA);
          param_dist *dist=new param_dist;
          
          dist->param_name=s[0];
          
          class_type ptype;
          if     (!strcmp(s[1],"SOIL"      )) { ptype=CLASS_SOIL; }
          else if(!strcmp(s[1],"VEGETATION")) { ptype=CLASS_VEGETATION; }
          else if(!strcmp(s[1],"LANDUSE"   )) { ptype=CLASS_LANDUSE; }
          else if(!strcmp(s[1],"TERRAIN"   )) { ptype=CLASS_TERRAIN; }
          else if(!strcmp(s[1],"GLOBALS"   )) { ptype=CLASS_GLOBAL; }
          else if(!strcmp(s[1],"SUBBASIN"  )) { ptype=CLASS_SUBBASIN; }
          else if(!strcmp(s[1],"GAUGE"     )) { ptype=CLASS_GAUGE; }
          else {
            ExitGracefully("ParseEnsembleFile: invalid parameter class in :ParameterDistributions command",BAD_DATA);
          }
          dist->param_class=ptype;
          dist->class_group=s[2];
          dist->default_val=s_to_d(s[3]);

          if     (!strcmp(s[4],"DIST_UNIFORM")) { dist->distribution=DIST_UNIFORM; }
          else if(!strcmp(s[4],"DIST_NORMAL" )) { dist->distribution=DIST_NORMAL; }
          else if(!strcmp(s[4],"DIST_GAMMA"  )) { dist->distribution=DIST_GAMMA; }
          else {
            ExitGracefully("ParseEnsembleFile: invalid distribution type in :ParameterDistributions command",BAD_DATA);
          }
          dist->distpar[0]=s_to_d(s[5]);
          dist->distpar[1]=s_to_d(s[6]);
          if(Len>7) { dist->distpar[2]=s_to_d(s[7]); }

          if (pEnsemble->GetType() == ENSEMBLE_MONTECARLO){
            ((CMonteCarloEnsemble*)(pEnsemble))->AddParamDist(dist);
          }
          else if(pEnsemble->GetType()==ENSEMBLE_DDS) {
            ((CDDSEnsemble*)(pEnsemble))->AddParamDist(dist);
          }
        }
      }
      break;
    }
    case(11):  //----------------------------------------------
    {/*:ObjectiveFunction [SBID] [Diagnostic] {Period}*/
      /*ObjectiveFunction READFROMFILE [filename] //LATER! */
      if(Options.noisy) { cout <<":ObjectiveFunction"<<endl; }

      int width=DOESNT_EXIST;
      string tmp = CStateVariable::SVStringBreak(s[2], width); //using other routine to grab width
      diag_type diag=StringToDiagnostic(tmp);
      
      if (diag==DIAG_UNRECOGNIZED) {
        ExitGracefully("ParseEnsembleFile::unknown diagnostic in :ObjectiveFunction command.",BAD_DATA_WARN);
      }
      string per_string="ALL";
      if (Len>=4){per_string=to_string(s[3]); }
      
      //if diagnostic doesn't exist,
      //CDiagnostic *pDiag=CDiagnostic(diag);
      //pModel->AddDiagnostic(pDiag);
      if(pEnsemble->GetType()==ENSEMBLE_DDS) {
        CDDSEnsemble *pDDS=((CDDSEnsemble*)(pEnsemble));
        pDDS->SetCalibrationTarget(s_to_l(s[1]),diag,per_string);
      }
      else {
        WriteWarning(":ObjectiveFunction command will be ignored; no calibration method specified.",Options.noisy);
      }
      break;
    }
    case(12):  //----------------------------------------------
    { //:ForcingPerturbation [forcingtype] [distrib] [dist param1] [dist param2] [adj_type] {HRU_Group}
      //:ForcingPerturbation RAINFALL DIST_GAMMA [dist param1] [dist param2] [adj_type] {HRU_Group}
      if(Options.noisy) { cout <<":ForcingPerturbation"<<endl; }
      if(pEnsemble->GetType()==ENSEMBLE_ENKF) {
        
        forcing_type ftyp=GetForcingTypeFromString(s[1]);
        disttype distrib=DIST_NORMAL;
        if     (!strcmp(s[2],"DIST_UNIFORM")) { distrib=DIST_UNIFORM; }
        else if(!strcmp(s[2],"DIST_NORMAL" )) { distrib=DIST_NORMAL; }
        else if(!strcmp(s[2],"DIST_GAMMA"  )) { distrib=DIST_GAMMA; }
        else {
          ExitGracefully("ParseEnsembleFile: invalid distribution type in :ForcingPerturbation command",BAD_DATA);
        }
        double distpars[3];
        adjustment adj=ADJ_ADDITIVE;
        distpars[0]=s_to_d(s[3]);
        distpars[1]=s_to_d(s[4]);
        distpars[2]=0.0;
        int kk=DOESNT_EXIST;
        if     (!strcmp(s[5],"ADDITIVE"       )) { adj=ADJ_ADDITIVE; }
        else if(!strcmp(s[5],"MULTIPLICATIVE" )) { adj=ADJ_MULTIPLICATIVE; }
        else {
          ExitGracefully("ParseEnsembleFile: invalid perturbation adjustment type in :ForcingPerturbation command",BAD_DATA);
        }
        if(Len>=7) {
          kk=pModel->GetHRUGroup(s[6])->GetGlobalIndex();
        }
        CEnKFEnsemble* pEnKF=((CEnKFEnsemble*)(pEnsemble));
        pEnKF->AddPerturbation(ftyp,distrib,distpars,kk, adj);
      }
      else {
        WriteWarning(":ForcingPerturbation command will be ignored; only valid for EnKF ensemble simulation.",Options.noisy);
      }
      break;
    }
    case(13):  //----------------------------------------------
    { //:AssimilatedState STREAMFLOW [SubBasinGroup] # only STREAMFLOW to be supported initially
      //:AssimilatedState SNOW       [HRUGroup]
      if(Options.noisy) { cout <<":ForcingPerturbation"<<endl; }
      if(pEnsemble->GetType()==ENSEMBLE_ENKF) {
        sv_type sv;
        int lay;
        sv=CStateVariable::StringToSVType(s[1],lay,true);
        int kk=DOESNT_EXIST;

        int i=pModel->GetStateVarIndex(sv,lay);
        if ((i == DOESNT_EXIST) && (sv!=STREAMFLOW) && (sv!=RESERVOIR_STAGE) ){
          string warn="State variable "+to_string(s[1])+" does not exist in this model and will be ignored in the :AssimilatedState command";
          WriteWarning(warn,Options.noisy);
          break;
        }
        if ((sv==STREAMFLOW) || (sv==RESERVOIR_STAGE)) {
          if(pModel->GetSubBasinGroup(s[2])==NULL) {
            ExitGracefully("ParseEnsembleFile: :AssimilatedState subbasin/HRU group does not exist",BAD_DATA_WARN);
            break;
          }
          kk=pModel->GetSubBasinGroup(s[2])->GetGlobalIndex();
        }

        else {
          if(pModel->GetHRUGroup(s[2])==NULL) {
            ExitGracefully("ParseEnsembleFile: :AssimilatedState subbasin/HRU group does not exist",BAD_DATA_WARN);
            break;
          }
          kk=pModel->GetHRUGroup(s[2])->GetGlobalIndex();
        }
        
        CEnKFEnsemble* pEnKF=((CEnKFEnsemble*)(pEnsemble));
        pEnKF->AddAssimilationState(sv,lay,kk);
      }
      else {
        WriteWarning(":AssimilatedState command will be ignored; only valid for EnKF ensemble simulation.",Options.noisy);
      }
      break;
    }
    case(14):  //----------------------------------------------
    {/*:SolutionRunName {runname}*/
      if(Options.noisy) { cout <<":SolutionRunName"<<endl; }
      if(pEnsemble->GetType()==ENSEMBLE_ENKF) {
        CEnKFEnsemble* pEnKF=((CEnKFEnsemble*)(pEnsemble));
        string runname=Options.run_name; //default
        if (Len>1){runname=to_string(s[1]); }
        if (Options.warm_ensemble_run!=""){ runname= Options.warm_ensemble_run;} //overwritten if handed in via command line
        pEnKF->SetWarmRunname(runname); 
      }
      else {
        WriteWarning(":SolutionRunName command will be ignored; only valid for EnKF ensemble simulation.",Options.noisy);
      }
      break;
    }
    case(15):  //----------------------------------------------
    {/*:WindowSize [# of timesteps][*/
      if(Options.noisy) { cout <<":WindowSize"<<endl; }
      if(pEnsemble->GetType()==ENSEMBLE_ENKF) {
        CEnKFEnsemble* pEnKF=((CEnKFEnsemble*)(pEnsemble));
        pEnKF->SetWindowSize(s_to_i(s[1])); 
      }
      else {
        WriteWarning(":WindowSize command will be ignored; only valid for EnKF ensemble simulation.",Options.noisy);
      }
      break;
    }
    case(16):  //----------------------------------------------
    { //:ObservationErrorModel [statetype] [distrib] [dist param1] [dist param2] [adj_type]
      //:ObservationErrorModel STREAMFLOW DIST_NORMAL 1 0.07 MULTIPLICATIVE
      if(Options.noisy) { cout <<":ObservationErrorModel"<<endl; }
      if(pEnsemble->GetType()==ENSEMBLE_ENKF) {
        
        sv_type sv;
        int lay;
        sv=CStateVariable::StringToSVType(s[1],lay,true);

        disttype distrib=DIST_NORMAL;
        if     (!strcmp(s[2],"DIST_UNIFORM")) { distrib=DIST_UNIFORM; }
        else if(!strcmp(s[2],"DIST_NORMAL" )) { distrib=DIST_NORMAL; }
        else if(!strcmp(s[2],"DIST_GAMMA"  )) { distrib=DIST_GAMMA; }
        else {
          ExitGracefully("ParseEnsembleFile: invalid distribution type in :ObservationErrorModel command",BAD_DATA);
        }

        double distpars[3];
        distpars[0]=s_to_d(s[3]);
        distpars[1]=s_to_d(s[4]);
        distpars[2]=0.0;

        adjustment adj=ADJ_ADDITIVE;
        if     (!strcmp(s[5],"ADDITIVE"       )) { adj=ADJ_ADDITIVE; }
        else if(!strcmp(s[5],"MULTIPLICATIVE" )) { adj=ADJ_MULTIPLICATIVE; }
        else {
          ExitGracefully("ParseEnsembleFile: invalid perturbation adjustment type in :ObservationErrorModel command",BAD_DATA);
        }

        CEnKFEnsemble* pEnKF=((CEnKFEnsemble*)(pEnsemble));
        pEnKF->AddObsPerturbation(sv,distrib,distpars,adj);
      }
      else {
        WriteWarning(":ObservationErrorModel command will be ignored; only valid for EnKF ensemble simulation.",Options.noisy);
      }
      break;
    }
    case(17):  //----------------------------------------------
    {
      break;
    }
    case(18):  //----------------------------------------------
    {/*:EnKFMode [mode]*/
      if(Options.noisy) { cout <<":EnKFMode"<<endl; }
      if(pEnsemble->GetType()==ENSEMBLE_ENKF) {
        CEnKFEnsemble* pEnKF=((CEnKFEnsemble*)(pEnsemble));
        if (pEnKF->GetEnKFMode()!=ENKF_UNSPECIFIED){break;} //set in run_info.nc, do not override

        EnKF_mode mode=ENKF_SPINUP;

        if      (!strcmp(s[1],"ENKF_SPINUP"       )){mode=ENKF_SPINUP;}
        else if (!strcmp(s[1],"ENKF_CLOSED_LOOP"  )){mode=ENKF_CLOSED_LOOP;}
        else if (!strcmp(s[1],"ENKF_OPEN_LOOP"    )){mode=ENKF_OPEN_LOOP;}
        else if (!strcmp(s[1],"ENKF_FORECAST"     )){mode=ENKF_FORECAST;}
        else if (!strcmp(s[1],"ENKF_OPEN_FORECAST")){mode=ENKF_OPEN_FORECAST;}
   
        else if (!strcmp(s[1],"ENKF_CLOSEDLOOP"   )){mode=ENKF_CLOSED_LOOP;}
        else if (!strcmp(s[1],"ENKF_OPENLOOP"     )){mode=ENKF_OPEN_LOOP;}
        else {
          cout<<s[1]<<endl;
          ExitGracefully("ParseEnsembleFile: EnKFMode-  invalid mode specified",BAD_DATA_WARN);
        }
        pEnKF->SetEnKFMode(mode); 
      }
      else {
        WriteWarning(":EnKFMode command will be ignored; only valid for EnKF ensemble simulation.",Options.noisy);
      }
      break;
    }
    case(19):  //----------------------------------------------
    {/*:ExtraRVTFilename [filename.rvt]*/
      if(Options.noisy) { cout <<":ExtraRVTFilename"<<endl; }
      if(pEnsemble->GetType()==ENSEMBLE_ENKF) {
        CEnKFEnsemble* pEnKF=((CEnKFEnsemble*)(pEnsemble));
        string file=CorrectForRelativePath(s[1],Options.rvi_filename);//with .rvt extension! (may include * wildcard operator for filename)
        pEnKF->SetExtraRVTFile(file);
      }
      else {
        WriteWarning(":ExtraRVTFilename command will be ignored; only valid for EnKF ensemble simulation.",Options.noisy);
      }
      break;
    }
    case(101)://----------------------------------------------
    {/*:AssimilateStreamflow  [SBID]*/
      if(Options.noisy) { cout <<"Assimilate streamflow"<<endl; }
      long SBID=s_to_l(s[1]);
      if(pModel->GetSubBasinByID(SBID)==NULL) {
        WriteWarning("ParseTimeSeries::Trying to assimilate streamflow at non-existent subbasin "+to_string(SBID),Options.noisy);
        break;
      }
      bool ObsExists=false;
      for(int i=0; i<pModel->GetNumObservedTS(); i++) {
        if(IsContinuousFlowObs2(pModel->GetObservedTS(i),SBID)) {
          ObsExists=true;
          break;
        }
      }
      if(ObsExists) {
        pModel->GetSubBasinByID(SBID)->IncludeInAssimilation();
      }
      else {
        string warn="ParseEnsembleSeriesFile::AssimilateStreamflow: no observation time series associated with subbasin "+to_string(SBID)+". Cannot assimilate flow in this subbasin.";
        WriteWarning(warn.c_str(),Options.noisy);
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
        if     (!strcmp(s[0],":FileType"    )) { if(Options.noisy) { cout<<"Filetype"   <<endl; } }//do nothing
        else if(!strcmp(s[0],":Application" )) { if(Options.noisy) { cout<<"Application"<<endl; } }//do nothing
        else if(!strcmp(s[0],":Version"     )) { if(Options.noisy) { cout<<"Version"    <<endl; } }//do nothing
        else if(!strcmp(s[0],":WrittenBy"   )) { if(Options.noisy) { cout<<"WrittenBy"  <<endl; } }//do nothing
        else if(!strcmp(s[0],":CreationDate")) { if(Options.noisy) { cout<<"CreationDate"<<endl; } }//do nothing
        else if(!strcmp(s[0],":SourceFile"  )) { if(Options.noisy) { cout<<"SourceFile"  <<endl; } }//do nothing
        else 
        {
          string warn ="IGNORING unrecognized command: " + string(s[0])+ " in .rve file";
          WriteWarning(warn,Options.noisy);
        }
      }
      break;
      default:
      {
        string errString = "Unrecognized command in .rve file:\n   " + string(s[0]);
        ExitGracefully(errString.c_str(),BAD_DATA_WARN);//STRICT
      }
      break;
      }
    }
    }//end switch(code)

    end_of_file=pp->Tokenize(s,Len);

    //return after file redirect, if in secondary file
    if((end_of_file) && (pMainParser!=NULL))
    {
      INPUT2.clear();
      INPUT2.close();
      delete pp;
      pp=pMainParser;
      pMainParser=NULL;
      end_of_file=pp->Tokenize(s,Len);
    }
  } //end while !end_of_file
  ENS.close();

  delete pp;
  pp=NULL;
  return true;
}
