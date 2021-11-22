/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2021 the Raven Development Team
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Model.h"
#include "StateVariables.h"
#include "HydroUnits.h"
#include "ParseLib.h"
#include "ModelEnsemble.h"

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
    <100         : ignored/special
    0   thru 100 : All other
    ------------------------------------------------------------------
    */

    code=0;
    //---------------------SPECIAL -----------------------------
    if     (Len==0)                                       { code=-1; }//blank line
    else if(IsComment(s[0],Len))                          { code=-2; }//comment
    else if(!strcmp(s[0],":RedirectToFile"))              { code=-3; }//redirect to secondary file
    else if(!strcmp(s[0],":End"))                         { code=-4; }//stop reading
    
    //-------------------MODEL ENSEMBLE PARAMETERS----------------
    else if(!strcmp(s[0],":OutputDirectoryFormat"))       { code=1; }
    else if(!strcmp(s[0],":RunNameFormat"))               { code=2; }
    else if(!strcmp(s[0],":ParameterDistributions"))      { code=10; }
    else if(!strcmp(s[0],":ObjectiveFunction"))           { code=11; }


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
    case(1):  //----------------------------------------------
    {/*:OutputDirectoryFormat*/
      if(Options.noisy) { cout <<"OutputDirectoryFormat"<<endl; } 
      string outdir="";
      for (i=1;i<Len;i++){outdir=outdir+to_string(s[i]); }
      //outdir=CorrectForRelativePath(outdir,Options.rvi_filename); (not for directory?)
      pEnsemble->SetOutputDirectory(outdir);
      break;
    }
    case(2):  //----------------------------------------------
    {/*:RunNameFormat*/
      if(Options.noisy) { cout <<"RunNameFormat"<<endl; } 
      pEnsemble->SetRunNames(s[1]); //No spaces!
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

      diag_type diag=DIAG_NASH_SUTCLIFFE;
      if     (!strcmp(s[2],"NASH_SUTCLIFFE"    )) { diag=DIAG_NASH_SUTCLIFFE; }
      else if(!strcmp(s[2],"RMSE"              )) { diag=(DIAG_RMSE); }
      else if(!strcmp(s[2],"PCT_BIAS"          )) { diag=(DIAG_PCT_BIAS); }
      else if(!strcmp(s[2],"ABSERR"            )) { diag=(DIAG_ABSERR); }
      else if(!strcmp(s[2],"ABSMAX"            )) { diag=(DIAG_ABSMAX); }
      else if(!strcmp(s[2],"PDIFF"             )) { diag=(DIAG_PDIFF); }
      else if(!strcmp(s[2],"TMVOL"             )) { diag=(DIAG_TMVOL); }
      else if(!strcmp(s[2],"RCOEF"             )) { diag=(DIAG_RCOEF); }
      else if(!strcmp(s[2],"NSC"               )) { diag=(DIAG_NSC); }
      else if(!strcmp(s[2],"RSR"               )) { diag=(DIAG_RSR); }
      else if(!strcmp(s[2],"R2"                )) { diag=(DIAG_R2); }
      else if(!strcmp(s[2],"CUMUL_FLOW"        )) { diag=(DIAG_CUMUL_FLOW); }
      else if(!strcmp(s[2],"LOG_NASH"          )) { diag=(DIAG_LOG_NASH); }
      else if(!strcmp(s[2],"KLING_GUPTA"       )) { diag=(DIAG_KLING_GUPTA); }
      else if(!strcmp(s[2],"KLING_GUPTA_DEVIATION")) { diag=(DIAG_KLING_GUPTA_DEVIATION); }
      else if(!strcmp(s[2],"NASH_SUTCLIFFE_DER")) { diag=(DIAG_NASH_SUTCLIFFE_DER); }
      else if(!strcmp(s[2],"RMSE_DER"          )) { diag=(DIAG_RMSE_DER); }
      else if(!strcmp(s[2],"KLING_GUPTA_DER"   )) { diag=(DIAG_KLING_GUPTA_DER); }
      else if(!strcmp(s[2],"MBF"               )) { diag=(DIAG_MBF); }
      else {
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
