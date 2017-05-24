/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Model.h"
#include "StateVariables.h"
#include "HydroUnits.h"
#include "ParseLib.h"

//////////////////////////////////////////////////////////////////
/// \brief Parses Initial conditions file
/// \details model.rvc: input file that defines HRU and Subbasin initial conditions\n
///
/// \param *&pModel [out] Reference to model object
/// \param &Options [out] Global model options information
/// \return True if operation is successful
//
bool ParseInitialConditionsFile(CModel *&pModel, const optStruct &Options)
{
  int         i,k;              //counters
  long        SBID;             //subbasin ID
  //CHydroUnit *pHRU;           //temporary pointers to HRUs, Subbasins
  CSubBasin  *pSB;
  bool        ended(false);

  ifstream    IC;
  IC.open(Options.rvc_filename.c_str());
  if (IC.fail()){
    cout << "ERROR opening initial conditions file: "<<Options.rvc_filename<<endl; return false;}

  int   Len,line(0),code;
  char *s[MAXINPUTITEMS];
  CParser *pp=new CParser(IC,Options.rvc_filename,line);

  ifstream INPUT2;           //For Secondary input
  CParser *pMainParser=NULL; //for storage of main parser while reading secondary files

  if (Options.noisy){
    cout <<"======================================================"<<endl;
    cout <<"Parsing Initial conditions File " << Options.rvc_filename <<"..."<<endl;
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
    //-------------------MODEL INITIAL CONDITIONS----------------
    else if  (!strcmp(s[0],":BasinInitialConditions"   )){code=1;  }
    else if  (!strcmp(s[0],":HRUInitialConditions"     )){code=2;  }
    else if  (!strcmp(s[0],":EndHRUInitialConditions"  )){code=-2; }
    else if  (!strcmp(s[0],":ALL"                      )){code=3;  }
    else if  (!strcmp(s[0],":UniformInitialConditions" )){code=3;  }
    else if  (!strcmp(s[0],":HRUStateVariableTable"    )){code=4;  }
    else if  (!strcmp(s[0],":EndHRUStateVariableTable" )){code=-2;  }
    else if  (!strcmp(s[0],":InitialConditions"        )){code=5;  }
    else if  (!strcmp(s[0],":EndInitialConditions"     )){code=-2; }
    else if  (!strcmp(s[0],":BasinStateVariables"      )){code=6;  }
    else if  (!strcmp(s[0],":EndBasinStateVariables"   )){code=-2; }
    else if  (!strcmp(s[0],":InitialReservoirFlow"     )){code=7; }
    else if  (!strcmp(s[0],":InitialReservoirStage"    )){code=8; }
    else if  (!strcmp(s[0],":TimeStamp"                )){code=10; }

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
        filename = filedir + "//" + filename;
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
        ":BasinInitialConditions"
        {SBID, {flow x num_segments}} x nSubBasins
        :EndBasinInitialConditions
      */
      if (Options.noisy) {cout <<"Basin Initial Conditions..."<<endl;}
      while ((Len==0) || (strcmp(s[0],":EndBasinInitialConditions")))
      {
        pp->Tokenize(s,Len);
        if      (IsComment(s[0], Len)){}//comment line
        else if (!strcmp(s[0],":Attributes"  )){}//ignored by Raven - needed for GUIs
        else if (!strcmp(s[0],":Units"       )){}//ignored by Raven - needed for GUIs
        else if (!strcmp(s[0],":EndBasinInitialConditions")){}//done
        else
        {
          ExitGracefullyIf(Len<2,
                           "Parse HRU File: incorrect number of terms in SubBasin initial conditions",BAD_DATA);
          SBID=s_to_l(s[0]);
          pSB=NULL;
          pSB=pModel->GetSubBasinByID(SBID);
          if (pSB!=NULL){
            double *aQout=new double [1];
            aQout[0]=s_to_d(s[1]);
            pSB->SetQoutArray(DOESNT_EXIST,aQout,aQout[0]);
            delete [] aQout;
          }
          else          {
            WriteWarning("Subbasin "+to_string(SBID)+" not in model, cannot set initial conditions",Options.noisy);
          }
        }
      }
      break;
    }
    case(2):  //----------------------------------------------
    { /* :HRUInitialConditions
       */
      if (Options.noisy) {cout <<"   Reading HRU Initial Conditions..."<<endl;}
      //does nothing, just an input file placeholder
      break;
    }
    case(3):  //----------------------------------------------
    { /* :UniformInitialConditions [svtype] [svval]
         -sets blanket IC values for a single state variable across all HRUs
         -only needed if IC not equal to zero
      */
      if (Options.noisy) {cout <<"Initial Conditions (Uniform)"<<endl;}
      if (Len==3)
      {
        int     layer_ind(-1);
        int     SVind;
        sv_type typ;

        typ=CStateVariable::StringToSVType(s[1],layer_ind,false);

        if (typ==UNRECOGNIZED_SVTYPE){
          WriteWarning(":UniformInitialConditions: unrecognized state variable type "+to_string(s[1]),Options.noisy);
          break;
        }

        SVind=pModel->GetStateVarIndex(typ,layer_ind);

        if (SVind!=DOESNT_EXIST){
          for (k=0;k<pModel->GetNumHRUs();k++){
            pModel->GetHydroUnit(k)->SetStateVarValue(SVind,s_to_d(s[2]));
          }
        }
        else{
          WriteWarning("Initial conditions specified for state variable not in model ("+to_string(s[1])+")",Options.noisy);
        }
      }
      else{
        pp->ImproperFormat(s);
      }
      break;
    }
    case(4):  //----------------------------------------------
    { /* :HRUStateVariableTable (formerly :IntialConditionsTable)
         ":Attributes" {list of state variables}
         ":Units" {list of units}
         {HRUID, SV1,SV2,SV3} x nHRUs
         :EndHRUStateVariableTable
         - sets unique IC values for different HRUs and one or multiple state variables
      */
      if (Options.noisy) {cout <<"   Reading HRU Initial Condition Table..."<<endl;}

      pp->Tokenize(s,Len);

      if (strcmp(s[0],":Attributes")){
        WriteWarning(":HRUStateVariableTable command: first line must begin with :Attributes",Options.noisy);break;}
      int nSV=Len-1;
      int *SVinds=new int [nSV];
      for (i=0;i<nSV;i++)
      {
        int     layer_ind(DOESNT_EXIST);
        sv_type typ;
        typ=CStateVariable::StringToSVType(s[i+1],layer_ind,false);

        if (typ==UNRECOGNIZED_SVTYPE){
          WriteWarning(":HRUStateVariableTable: unrecognized state variable type "+to_string(s[i+1]),Options.noisy);
          SVinds[i]=DOESNT_EXIST;

          //cout<<"!State Var:"<<s[i+1];
          //cout<<" State Var Ind:"<<SVinds[i]<<endl<<g_output_directory<<endl;
        }
        else{
          SVinds[i]=pModel->GetStateVarIndex(typ,layer_ind);

          //cout<<"State Var:"<<s[i+1];
          //cout<<" State Var Ind:"<<SVinds[i]<<endl;
          if (SVinds[i]==DOESNT_EXIST){
            WriteWarning("Initial conditions specified for state variable not in model ("+to_string(s[i+1])+")",Options.noisy);
          }
          if ((typ==ATMOS_PRECIP) || (typ==ATMOSPHERE) || (typ==GLACIER_ICE)){//initial conditions of cumulative precip, evap, and glacier loss ignored, left at zero
            SVinds[i]=DOESNT_EXIST;
          }
          /// \todo [funct] must also zero out concentrations linked to atmos_precip, atmosphere, and glacier ice

        }
      }
      k=0;
      while ((Len==0) || (strcmp(s[0],":EndHRUStateVariableTable")))
      {
        pp->Tokenize(s,Len);

        if      (IsComment(s[0], Len)){}//comment line
        else if (!strcmp(s[0],":Units")){}//ignored by raven
        else if (!strcmp(s[0],":EndHRUStateVariableTable")){}//end command
        else
        {
          ExitGracefullyIf(Len!=(nSV+1),
                           "Parse Initial Conditions File: incorrect number of columns in HRU State Variable Table ",BAD_DATA);
          ExitGracefullyIf(k>=pModel->GetNumHRUs(),
                           "Parse: HRUStateVariableTable: # of rows more than # of HRUs", BAD_DATA);
          for (i=0;i<nSV;i++){
            if (SVinds[i]!=DOESNT_EXIST){
              pModel->GetHydroUnit(k)->SetStateVarValue(SVinds[i],s_to_d(s[i+1]));
            }
          }
          k++;
        }
      }//end while //end for k=0 to nHRUs
      // concern: not properly referenced by HRU ID!
      delete [] SVinds;
      break;
    }
    case(5):  //----------------------------------------------
    { /*
        ":InitialConditions" {SV_NAME}
        {v1,v2,v3,v4} x nHRUs
        :EndInitialConditions
      */
      if (Len<2){pp->ImproperFormat(s);}
      else{
        if (Options.noisy) {cout <<"   Reading Initial Conditions for "<<s[1]<<endl;}
        string stateVariable= s[1];
        sv_type SVtype;
        int     SVlayer=0;
        SVtype=CStateVariable::StringToSVType(s[1],SVlayer,false);
        if (SVtype==UNRECOGNIZED_SVTYPE){
          WriteWarning("Unrecognized State Variable type " + string(s[1]) +" in :InitialConditions command",Options.noisy);
          break;
        }
        int nHRUs=pModel->GetNumHRUs();
        double *v=new double [nHRUs];
        int count =0;
        while (count<pModel->GetNumHRUs())
        {
          if(!strcmp(s[0],":EndInitialConditions"))
          { break; }
          else
          {
            pp->Tokenize(s,Len);
            if      (IsComment(s[0], Len)){}//comment line
            else if (!strcmp(s[0],":EndInitialConditions")){}//done
            else
            {
              for (int i=0;i<Len;i++)
              {
                if (count<nHRUs){
                  v[count]=s_to_d(s[i]);
                }
                count++;
              }
            }
          }

          int iSV=pModel->GetStateVarIndex(SVtype,SVlayer);
          if (iSV!=DOESNT_EXIST){
            if ((SVtype!=ATMOS_PRECIP) && (SVtype!=ATMOSPHERE)){//initial conditions of cumulative precip & evap ignored, left at zero
              for (int k=0;k<pModel->GetNumHRUs();k++){
                pModel->GetHydroUnit(k)->SetStateVarValue(iSV,v[k]);
              }
            }
          }
        }/*end while*/
        if (count!=nHRUs)
        {
          WriteWarning("Initial condition count is incorrect for state variable \""+(stateVariable)+"\"",Options.noisy);
          //cout<<"       NumHRUs="<<pModel->GetNumHRUs()<<" values found="<<count<<endl;
        }
        delete [] v;
        int iSV=pModel->GetStateVarIndex(SVtype,SVlayer);
        if (iSV==DOESNT_EXIST){
          WriteWarning("Unused state Variable " + string(s[1]) +" in :InitialConditions command will be ignored",Options.noisy);
        }
      }

      break;
    }
    case(6):  //----------------------------------------------
    {
      /*
        ":BasinStateVariables"
        :BasinIndex ID, name
        :ChannelStorage [val]
        :RivuletStorage [val]
        :Qout [nsegs] [aQout x nsegs] [aQoutLast]
        :Qlat [nQlatHist] [aQlatHist x nQlatHist] [QlatLast]
        :Qin  [nQinHist] [aQinHist x nQinHist]
        {:ResStage [stage] [last stage]}
        :BasinIndex 1D, name
        ...
        :EndBasinStateVariables
      */
      if (Options.noisy) {cout <<"   Reading Basin State Variables"<<endl;}
      bool done=false;
      CSubBasin *pBasin=NULL;
      int SBID=DOESNT_EXIST;
      do
      {
        bool eof=pp->Tokenize(s,Len);
        if (eof){done=true;}
        if      (Len==0){}//do nothing
        else if (IsComment(s[0],Len)){}//do nothing
        else if (!strcmp(s[0],":BasinIndex"))
        {
          SBID=s_to_l(s[1]);
          pBasin=pModel->GetSubBasinByID(SBID);
          ExitGracefullyIf(pBasin==NULL,
                           "ParseInitialConditionsFile: bad basin index in .rvc file",BAD_DATA);
          if (Options.noisy) {cout <<"     Reading Basin "<<pBasin->GetID()<<": "<<pBasin->GetName()<<endl;}
        }
        else if (!strcmp(s[0],":ChannelStorage"))
        {
          if (Len>2){
            pBasin->SetChannelStorage(s_to_d(s[1]));
          }
        }
        else if (!strcmp(s[0],":RivuletStorage"))
        {
          if (Len>2){
            pBasin->SetRivuletStorage(s_to_d(s[1]));
          }
        }
        else if (!strcmp(s[0],":Qout"))
        {
          if (Len>2){
            int nsegs=s_to_i(s[1]);
            if (Len>=nsegs+3){
              double *aQout=new double [nsegs+1];
              for (int i=0;i<=nsegs;i++){
                aQout[i]=s_to_d(s[i+2]);
              }
              pBasin->SetQoutArray(nsegs,aQout,aQout[nsegs]);
              delete [] aQout;
            }
          }
        }
        else if (!strcmp(s[0],":Qlat"))
        {
          if (Len>2){
            int histsize=s_to_i(s[1]);
            if (Len>=histsize+3){
              double *aQlat=new double [histsize+1];
              for (int i=0;i<=histsize;i++){
                aQlat[i]=s_to_d(s[i+2]);
              }
              pBasin->SetQlatHist(histsize,aQlat,aQlat[histsize]);
              delete [] aQlat;
            }
          }
        }
        else if (!strcmp(s[0],":Qin"))
        {
          if (Len>2){
            int histsize=s_to_i(s[1]);
            if (Len>=histsize+2){
              double *aQin=new double [histsize];
              for (int i=0;i<histsize;i++){
                aQin[i]=s_to_d(s[i+2]);
              }
              pBasin->SetQinHist(histsize,aQin);
              delete [] aQin;
            }
          }
        }
        else if (!strcmp(s[0],":ResStage"))
        {
          if (Len>=3){
            pBasin->SetReservoirStage(s_to_d(s[1]));
          }
        }
        else if (!strcmp(s[0],":EndBasinStateVariables"))
        {
          done=true;
        }
      } while (!done);
      break;
    }
    case(7):  //----------------------------------------------
    {
      //:InitialReservoirFlow [SBID] [flow in m3/s]
      int SBID=s_to_l(s[1]);
      CSubBasin *pBasin=pModel->GetSubBasinByID(SBID);
      ExitGracefullyIf(pBasin==NULL,
                       "ParseInitialConditionsFile: bad basin index in :InitialReservoirFlow command (.rvc file)",BAD_DATA);
      pBasin->SetReservoirFlow(AutoOrDouble(s[2]));
      break;
    }
    case(8):  //----------------------------------------------
    {
      //:InitialReservoirStage [SBID] [stage]
      int SBID=s_to_l(s[1]);
      CSubBasin *pBasin=pModel->GetSubBasinByID(SBID);
      ExitGracefullyIf(pBasin==NULL,
                       "ParseInitialConditionsFile: bad basin index in :InitialReservoirFlow command (.rvc file)",BAD_DATA);
      pBasin->SetReservoirStage(s_to_d(s[2]));
      break;
    }
    case(10):  //----------------------------------------------
    {
      /*
        ":TimeStamp" [date] [hr]
      */
      //purely QA/QC check
      if (Len<3){break;}
      time_struct tt;
      tt=DateStringToTimeStruct(s[1],s[2]);
      if ((Options.julian_start_day!=tt.julian_day) || (Options.julian_start_year!=tt.year)){
        WriteWarning(":Timestamp command in initial conditions (.rvc) file is not consistent with :StartDate command in model (.rvi) file",Options.noisy);
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
          string warn ="IGNORING unrecognized command: " + string(s[0])+ " in .rvc file";
          WriteWarning(warn,Options.noisy);
        }
      }
      break;
      default:
      {
        string errString = "Unrecognized command in .rvc file:\n   " + string(s[0]);
        ExitGracefully(errString.c_str(),BAD_DATA_WARN);//STRICT
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
  IC.close();

  // check quality of initial state variables
  //======================================================================
  CHydroUnit *pHRU;
  double *v=new double [pModel->GetNumStateVars()];
  for (int k=0;k<pModel->GetNumHRUs();k++)
  {
    pHRU=pModel->GetHydroUnit(k);
    for (int i=0; i<pModel->GetNumStateVars(); i++)
    {
      v[i]=pHRU->GetStateVarValue(i);
    }
    for (int i=0; i<pModel->GetNumStateVars(); i++)
    {
      double maxv=max(pHRU->GetStateVarMax(i,v,Options),0.0);
      if (v[i]>maxv)// check for capacity
      {
        string name=CStateVariable::GetStateVarLongName(pModel->GetStateVarType(i),pModel->GetStateVarLayer(i));
        string warn ="maximum state variable limit exceeded in initial conditions for " + name+ " (in HRU "+to_string(pHRU->GetID())+") in .rvc file";
        WriteWarning(warn,Options.noisy);

        if (!Options.keepUBCWMbugs){
          pHRU->SetStateVarValue(i,maxv);
        }
      }
    }
  }
  delete [] v;

  delete pp;
  pp=NULL;
  return true;
}
