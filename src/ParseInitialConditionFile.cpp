/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2022 the Raven Development Team
  ----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Model.h"
#include "StateVariables.h"
#include "HydroUnits.h"
#include "ParseLib.h"
#include "EnergyTransport.h"

void SetInitialStateVar(CModel *&pModel,const int SVind,const sv_type typ,const int m,const int k,const double &val);
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
  long long   SBID;             //subbasin ID
  //CHydroUnit *pHRU;           //temporary pointers to HRUs, Subbasins
  CSubBasin  *pSB;
  bool        ended(false);
  bool        in_ifmode_statement=false;

  time_struct tt;
  JulianConvert(0.0,Options.julian_start_day,Options.julian_start_year,Options.calendar,tt);

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


  //Initialize everything to zero (required for ensemble simulation)
  //--------------------------------------------------------------------------
  if (pModel->GetEnsemble()->GetNumMembers()>1)
  {
    // does not matter in EnKF/forecasting context, as it will get overwritten with .rvc state
    for (int i=0;i<pModel->GetNumStateVars();i++)
    {
      sv_type typ      =pModel->GetStateVarType(i);
      int     layer_ind=pModel->GetStateVarLayer(i);
      for(k=0;k<pModel->GetNumHRUs();k++) {
        SetInitialStateVar(pModel,i,typ,layer_ind,k,0.0);
      }
    }
    pModel->InitializeBasins(Options,true); //resets all flow rates to autocalc version
    for(int p=0;p<pModel->GetNumSubBasins();p++) {
      if(pModel->GetSubBasin(p)->GetReservoir()!=NULL) {
        pModel->GetSubBasin(p)->GetReservoir()->SetReservoirStage(0.0,0.0);
      }
    }
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
    string concname="UNKNOWN";
    code=0;
    //---------------------SPECIAL -----------------------------
    if       (Len==0)                                       {code=-1; }//blank line
    else if  (IsComment(s[0],Len))                          {code=-2; }//comment
    else if  (!strcmp(s[0],":End"                         )){code=-4; }//stop reading
    else if  (!strcmp(s[0],":IfModeEquals"                )){code=-5; }
    else if  (in_ifmode_statement)                          {code=-6; }
    else if  (!strcmp(s[0],":EndIfModeEquals"             )){code=-2; }//treat as comment - unused mode
    else if  (!strcmp(s[0],":RedirectToFile"              )){code=-3; }//redirect to secondary file
    //-------------------MODEL INITIAL CONDITIONS----------------
    else if  (!strcmp(s[0],":BasinInitialConditions"      )){code=1;  }

    else if  (!strcmp(s[0],":UniformInitialConditions"    )){code=2;  }
    else if  (!strcmp(s[0],":UniformInitialConcentration" )){code=3; concname=s[1];}
    else if  (!strcmp(s[0],":UniformInitialTemperature"   )){code=3; concname="TEMPERATURE";}

    else if  (!strcmp(s[0],":HRUStateVariableTable"       )){code=4; concname="";}
    else if  (!strcmp(s[0],":EndHRUStateVariableTable"    )){code=-2; }
    else if  (!strcmp(s[0],":InitialConcentrationTable"   )){code=4; concname=s[1]; }
    else if  (!strcmp(s[0],":EndInitialConcentrationTable")){code=-2; }
    else if  (!strcmp(s[0],":InitialTemperatureTable"     )){code=4; concname="TEMPERATURE";}
    else if  (!strcmp(s[0],":EndInitialTemperatureTable"  )){code=-2; }
    else if  (!strcmp(s[0],":InitialConditions"           )){code=5;  }//TO DEPRECATE
    else if  (!strcmp(s[0],":EndInitialConditions"        )){code=-2; }//TO DEPRECATE
    else if  (!strcmp(s[0],":BasinStateVariables"         )){code=6;  }
    else if  (!strcmp(s[0],":EndBasinStateVariables"      )){code=-2; }
    else if  (!strcmp(s[0],":BasinTransportVariables"     )){code=12; concname=s[1];}
    else if  (!strcmp(s[0],":EndBasinTransportVariables"  )){code=-2; }

    else if  (!strcmp(s[0],":InitialReservoirFlow"        )){code=7;  }
    else if  (!strcmp(s[0],":InitialReservoirStage"       )){code=8;  }

    //else if  (!strcmp(s[0],":InitialReservoirConcentration"       )){code=9;  concname=s[1];  }
    //else if  (!strcmp(s[0],":InitialReservoirTemperature"         )){code=9;  concname="TEMPERATURE"; }

    else if  (!strcmp(s[0],":TimeStamp"                   )){code=10; }
    else if  (!strcmp(s[0],":Nudge"                       )){code=11; }




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

      filename=CorrectForRelativePath(filename,Options.rvc_filename);

      INPUT2.open(filename.c_str());
      if (INPUT2.fail()){
        string warn=":RedirectToFile: Cannot find file "+filename;
        ExitGracefully(warn.c_str(),BAD_DATA);
      }
      else{
        if (pMainParser != NULL) {
          ExitGracefully("ParseInitialConditionsFile::nested :RedirectToFile commands (in already redirected files) are not allowed.",BAD_DATA);
        }
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
          SBID=s_to_ll(s[0]);
          pSB=NULL;
          pSB=pModel->GetSubBasinByID(SBID);
          if (pSB!=NULL){
            pSB->SetQout(s_to_d(s[1]));
          }
          else          {
            WriteWarning("Subbasin "+to_string(SBID)+" not in model, cannot set initial conditions",Options.noisy);
          }
        }
      }
      break;
    }
    case(2):  //----------------------------------------------
    { /* :UniformInitialConditions [svtype] [svval]
         -sets blanket IC values for a single state variable across all HRUs
         -only needed if IC not equal to zero
      */
      if (Options.noisy) {cout <<"Initial Conditions (Uniform)"<<endl;}
      if (Len>=3)
      {
        int     layer_ind(-1);
        int     SVind;
        sv_type typ;

        typ = pModel->GetStateVarInfo()->StringToSVType(s[1],layer_ind,false);

        if (typ==UNRECOGNIZED_SVTYPE){
          WriteWarning(":UniformInitialConditions: unrecognized state variable type "+to_string(s[1]),Options.noisy);
          break;
        }

        SVind=pModel->GetStateVarIndex(typ,layer_ind);

        if (SVind!=DOESNT_EXIST){
          for (k=0;k<pModel->GetNumHRUs();k++){
            SetInitialStateVar(pModel,SVind,typ,layer_ind,k,s_to_d(s[2]));
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
    case(3):  //----------------------------------------------
    { /*
      :UniformInitialConcentration [constit_name] [svtype] [value]
          or
      :UniformInitialTemperature [svtype] [value]
       specify uniform initial conditions in mg/L (or o/oo for isotopes) (or degrees C for temperature)
      */
      int shift=0;
      if(concname!="TEMPERATURE") { if(Options.noisy) { cout <<"Initial Concentrations (Uniform)"<<endl; } shift=1;}
      else                        { if(Options.noisy) { cout <<"Initial Temperatures (Uniform)"<<endl;   } }

      int c=DOESNT_EXIST;
      c=pModel->GetTransportModel()->GetConstituentIndex(concname);
      if(c==DOESNT_EXIST) {
        ExitGracefully("Invalid constituent name in :UniformInitialConcentration or :UniformInitialTemperature",BAD_DATA);
      }

      if(Len>=3)
      {
        int     layer_ind(-1);
        int     SVind;
        sv_type typ;

        typ = pModel->GetStateVarInfo()->StringToSVType(s[1+shift],layer_ind,false); //water state var
        if(typ==UNRECOGNIZED_SVTYPE) {
          WriteWarning(":UniformInitialConcentration: unrecognized state variable type "+to_string(s[1+shift]),Options.noisy);
          break;
        }

        SVind=pModel->GetStateVarIndex(typ,layer_ind); //water storage SVIND
        if(SVind==DOESNT_EXIST) {
          WriteWarning("Uniform initial conditions specified for water storage state variable not in model ("+to_string(s[1+shift])+")",Options.noisy);
          break;
        }

        int m=pModel->GetTransportModel()->GetLayerIndex(c,SVind);  //constituent layer
        SVind=pModel->GetStateVarIndex(CONSTITUENT,m); //constituent SVIND

        for(k=0;k<pModel->GetNumHRUs();k++) {
          SetInitialStateVar(pModel,SVind,CONSTITUENT,m,k,s_to_d(s[2+shift]));
        }
      }
      else {
        pp->ImproperFormat(s);
      }
      break;
    }
    case(4):  //----------------------------------------------
    { /* :HRUStateVariableTable (formerly :IntialConditionsTable)
           :Attributes {list of state variables}
           :Units {list of units}
           {HRUID, SV1,SV2,SV3,...} x (<=)nHRUs
           OR
           {HRUGroup, SV1,SV2,SV3,...}
         :EndHRUStateVariableTable
         - sets unique IC values for different HRUs and one or multiple state variables

         or

         :InitialTemperatureTable
           :Attributes {list of WATER state variables}
           :Units {list of units}
           {HRUID, T1,T2,T3,...} x (<=)nHRUs
           OR
           {HRUGroup, T1,T2,T3,...}
         :EndInitialTemperatureTable

         or
         :InitialConcentrationTable [constituent name]
           :Attributes {list of WATER state variables}
           :Units {list of units}
           {HRUID, C1,C2,C3,...} x (<=)nHRUs
            OR
           {HRUGroup, C1,C2,C3,...}
         :EndInitialConcentrationTable
      */
      if      ((concname==""           ) && (Options.noisy)) {cout <<"   Reading HRU Initial Condition Table..."<<endl;}
      else if ((concname=="TEMPERATURE") && (Options.noisy)) {cout <<"   Reading Temperature Initial Condition Table..."<<endl; }
      else if (                             (Options.noisy)) {cout <<"   Reading Concentration Initial Condition Table for "<<concname<<"..."<<endl; }

      pp->Tokenize(s,Len);
      int c=DOESNT_EXIST;
      if(concname!="") {
        c=pModel->GetTransportModel()->GetConstituentIndex(concname);
        if(c==DOESNT_EXIST) {
          ExitGracefully("Invalid constituent name in :InitialConcentrationTable or :InitialTemperatureTable",BAD_DATA);
        }
      }

      // Read :Attributes header-------------------------------------------------
      if (strcmp(s[0],":Attributes")){
        WriteWarning(":HRUStateVariableTable command: first line must begin with :Attributes",Options.noisy);break;}
      int nSV=Len-1;
      int *SVinds=new int [nSV];
      for (i=0;i<nSV;i++)
      {
        int     layer_ind(DOESNT_EXIST);
        sv_type typ;

        typ = pModel->GetStateVarInfo()->StringToSVType(s[i+1],layer_ind,false);

        if (typ==UNRECOGNIZED_SVTYPE){
          WriteWarning(":HRUStateVariableTable: unrecognized state variable type "+to_string(s[i+1]),Options.noisy);
          SVinds[i]=DOESNT_EXIST;
        }
        else{
          SVinds[i]=pModel->GetStateVarIndex(typ,layer_ind);

          if(c!=DOESNT_EXIST) { //convert water storage index to corresponding constituent index
            int m=pModel->GetTransportModel()->GetLayerIndex(c,SVinds[i]);
            SVinds[i]=DOESNT_EXIST;
            if(m!=DOESNT_EXIST) {
              SVinds[i]=pModel->GetStateVarIndex(CONSTITUENT,m);
            }
          }

          if (SVinds[i]==DOESNT_EXIST){
            WriteWarning("Initial conditions specified for state variable not in model ("+to_string(s[i+1])+")",Options.noisy);
          }

          if ((typ==ATMOS_PRECIP) || (typ==ATMOSPHERE) || (typ==GLACIER_ICE)){//initial conditions of cumulative precip, evap, and glacier loss ignored, left at zero
            SVinds[i]=DOESNT_EXIST;
          }

          ///  zero out concentrations/temperatures linked to atmos_precip, atmosphere, and glacier ice
          if(typ==CONSTITUENT) {
            int ii=pModel->GetTransportModel()->GetWaterStorIndexFromLayer(layer_ind);
            if (ii!=DOESNT_EXIST){
              sv_type wattyp=pModel->GetStateVarType(ii);
              if((wattyp==ATMOS_PRECIP) || (wattyp==ATMOSPHERE) || (wattyp==GLACIER_ICE)) {//initial concentrations conditions of cumulative precip, evap, and glacier loss ignored, left at zero
                SVinds[i]=DOESNT_EXIST;
              }
            }
          }
        }
      }

      // Read body of Table -----------------------------------------------------
      long long int HRUID;
      CHydroUnit *pHRU;
      CHRUGroup *pHRUGrp;
      while ( (Len==0) ||
              ((strcmp(s[0],":EndHRUStateVariableTable")) &&
               (strcmp(s[0],":EndInitialTemperatureTable")) &&
               (strcmp(s[0],":EndInitialConcentrationTable"))) )
      {
        pp->Tokenize(s,Len);

        if      (IsComment(s[0], Len)){}//comment line
        else if (!strcmp(s[0],":Units")){}//ignored by Raven
        else if (!strcmp(s[0],":EndHRUStateVariableTable")){}//end command
        else if (!strcmp(s[0],":EndInitialTemperatureTable")){}//end command
        else if (!strcmp(s[0],":EndInitialConcentrationTable")){}//end command
        else //row in SV table
        {
          ExitGracefullyIf(Len!=(nSV+1),
            "Parse :HRUStateVariableTable: incorrect number of columns in initial conditions table row (.rvc file)",BAD_DATA);

          pHRUGrp=pModel->GetHRUGroup(s[0]);

          HRUID=s_to_ll(s[0]);
          pHRU=pModel->GetHRUByID(HRUID);

          if((pHRU==NULL) && (pHRUGrp==NULL)){
            string warn="HRU ID or HRU Group ["+to_string(HRUID)+"] in .rvc initial conditions table not found in model";
            WriteWarning(warn,Options.noisy);
          }
          else{
            for(i=0;i<nSV;i++){
              if(SVinds[i]!=DOESNT_EXIST){
                if(c==DOESNT_EXIST) { //not concentration/temperature
                  if (pHRU!=NULL){
                    pHRU->SetStateVarValue(SVinds[i],s_to_d(s[i+1]));
                  }
                  else {
                    for (int kl = 0; kl < pHRUGrp->GetNumHRUs(); kl++) {
                      pHRUGrp->GetHRU(kl)->SetStateVarValue(SVinds[i],s_to_d(s[i+1]));
                    }
                  }
                }
                else{ //concentration/temperature
                  double val=s_to_d(s[i+1]);
                  int      m=pModel->GetStateVarLayer(SVinds[i]);
                  if (pHRU != NULL) {
                    SetInitialStateVar(pModel,SVinds[i],CONSTITUENT,m,pHRU->GetGlobalIndex(),val);
                  }
                  else{
                    for (int kl = 0; kl < pHRUGrp->GetNumHRUs(); kl++) {
                      SetInitialStateVar(pModel,SVinds[i],CONSTITUENT,m,pHRUGrp->GetHRU(kl)->GetGlobalIndex(),val);
                    }
                  }
                }
              }
            }
          }
        }//if Iscomment...
      }//end while

      delete [] SVinds;
      break;
    }
    case(5):  //----------------------------------------------
    { /* \todo[clean]: make this command obsolete - retained SOLELY for backward compatibility
        ":InitialConditions" {SV_NAME}
        {v1,v2,v3,v4} x nHRUs
        :EndInitialConditions
      */
	    WriteWarning(":InitialConditions command: this command will be deprecated. Please use :HRUStateVariableTable instead.",Options.noisy);
      if (Len<2){pp->ImproperFormat(s);}
      else{
        if (Options.noisy) {cout <<"   Reading Initial Conditions for "<<s[1]<<endl;}
        string stateVariable= s[1];
        sv_type SVtype;
        int     SVlayer=0;
        SVtype = pModel->GetStateVarInfo()->StringToSVType(s[1], SVlayer, false);
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
          //zero out rest
          int iSV=pModel->GetStateVarIndex(SVtype,SVlayer);
          if(iSV!=DOESNT_EXIST){
            if((SVtype!=ATMOS_PRECIP) && (SVtype!=ATMOSPHERE)){//initial conditions of cumulative precip & evap ignored, left at zero
              for(int k=count;k<pModel->GetNumHRUs();k++){
                pModel->GetHydroUnit(k)->SetStateVarValue(iSV,0);
              }
            }
          }
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
    { /* :BasinStateVariables
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
      long long SBID=DOESNT_EXIST;
      do
      {
        bool eof=pp->Tokenize(s,Len);
        if (eof){done=true;}
        if      (Len==0){}//do nothing
        else if (IsComment(s[0],Len)){}//do nothing
        else if (!strcmp(s[0],":BasinIndex"))
        {
          SBID=s_to_ll(s[1]);
          pBasin=pModel->GetSubBasinByID(SBID);
          ExitGracefullyIf(pBasin==NULL,
                           "ParseInitialConditionsFile: bad basin index in .rvc file",BAD_DATA);
          //if (Options.noisy) {cout <<"     Reading Basin "<<pBasin->GetID()<<": "<<pBasin->GetName()<<endl;}
        }
        else if (!strcmp(s[0],":ChannelStorage"))
        {
          if (Len>=2){
            pBasin->SetChannelStorage(s_to_d(s[1]));
          }
        }
        else if (!strcmp(s[0],":RivuletStorage"))
        {
          if (Len>=2){
            pBasin->SetRivuletStorage(s_to_d(s[1]));
          }
        }
        else if (!strcmp(s[0],":Qout"))
        {
          if (Len>2){
            int nsegs=s_to_i(s[1]);
            double *aQout=new double [nsegs+1];
            int i=2;
            int j=0;
            do {
              aQout[j]=s_to_d(s[i]);
              i++; j++;
              if ((i==Len) && (j<nsegs+1)) { end_of_file = pp->Tokenize(s, Len); i = 0; }
            } while ((j<=nsegs) && (!end_of_file));
            pBasin->SetQoutArray(nsegs,aQout,aQout[nsegs]);
            delete [] aQout;
          }
        }
        else if (!strcmp(s[0],":Qlat"))
        {
          if (Len>2){
            int histsize=s_to_i(s[1]);
            double *aQlat=new double [histsize+1];
            int i=2;
            int j=0;
            do {
              aQlat[j]=s_to_d(s[i]);
              i++; j++;
              if ((i==Len) && (j<histsize+1)){end_of_file=pp->Tokenize(s,Len); i=0;}
            } while ((j<=histsize) && (!end_of_file));
            pBasin->SetQlatHist(histsize,aQlat,aQlat[histsize]);
            delete [] aQlat;
          }
        }
        else if (!strcmp(s[0],":Qin"))
        {
          if (Len>2){
            int histsize=s_to_i(s[1]);
            double *aQin=new double [histsize];
            int i=2;
            int j=0;
            do {
              aQin[j]=s_to_d(s[i]);
              i++; j++;
              if ((i==Len) && (j<histsize)){end_of_file=pp->Tokenize(s,Len); i=0;}
            } while ((j<histsize) && (!end_of_file));

            pBasin->SetQinHist(histsize,aQin);
            delete [] aQin;
          }
        }
        else if (!strcmp(s[0],":ResStage"))
        {
          if (Len>=3){
            pBasin->GetReservoir()->SetReservoirStage(s_to_d(s[1]),s_to_d(s[2]));
          }
        }
        else if(!strcmp(s[0],":ResFlow"))
        {
          if(Len>=3) {
            pBasin->GetReservoir()->SetInitialFlow(s_to_d(s[1]),s_to_d(s[2]),tt,Options);
          }
        }
        else if (!strcmp(s[0], ":ControlFlow"))
        {
          if (Len >= 4) {
            pBasin->GetReservoir()->SetControlFlow(s_to_i(s[1]),s_to_d(s[2]),s_to_d(s[3]));
          }
        }
        else if(!strcmp(s[0],":ResDAscale"))
        {
          if(Len>=3) {
            pBasin->GetReservoir()->SetDataAssimFactors(s_to_d(s[1]),s_to_d(s[2]));
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
      SBID=s_to_ll(s[1]);
      CSubBasin *pBasin=pModel->GetSubBasinByID(SBID);
      if(pBasin==NULL) {
        ExitGracefully("ParseInitialConditionsFile: bad basin index in :InitialReservoirFlow command (.rvc file)",BAD_DATA_WARN);break;
      }
      time_struct tt;
      JulianConvert(0.0,Options.julian_start_day,Options.julian_start_year,Options.calendar,tt);
      pBasin->GetReservoir()->UpdateReservoir(tt,Options); //ensures correct discharge rating curve is used to calculate flow
      pBasin->GetReservoir()->SetInitialFlow(AutoOrDouble(s[2]),AutoOrDouble(s[2]),tt,Options);
      break;
    }
    case(8):  //----------------------------------------------
    {
      //:InitialReservoirStage [SBID] [stage]
      SBID=s_to_ll(s[1]);
      CSubBasin *pBasin=pModel->GetSubBasinByID(SBID);
      if(pBasin==NULL) {
        ExitGracefully("ParseInitialConditionsFile: bad basin index in :InitialReservoirFlow command (.rvc file)",BAD_DATA_WARN); break;
      }
      if(pBasin->GetReservoir()==NULL) {
        //warning to handle commented out reservoirs
        string warn="ParseInitialConditionsFile: bad basin index in :InitialReservoirStage command (.rvc file), no reservoir exists in basin " +to_string(SBID)+". Command will be ignored.";
        WriteWarning(warn.c_str(),Options.noisy);
      }
      else{
        pBasin->GetReservoir()->SetReservoirStage(s_to_d(s[2]),s_to_d(s[2]));
      }
      break;
    }
    case(10):  //----------------------------------------------
    {
      //:TimeStamp [yyyy-mm-dd] [00:00:00]
      //purely QA/QC check
      if (Len<3){break;}
      time_struct tt;
      tt=DateStringToTimeStruct(s[1],s[2],Options.calendar);
      if ((Options.julian_start_day!=tt.julian_day) || (Options.julian_start_year!=tt.year)){
        WriteWarning(":Timestamp command in initial conditions (.rvc) file is not consistent with :StartDate command in model (.rvi) file",Options.noisy);
      }
      break;
    }
    case(11): //------------------------------------------------
    {
      //:Nudge NUDGE_MULTIPLY [state_var] [mutiplier/add-on] [HRU Group]
      // e.g.,
      //:Nudge NUDGE_MULTIPLY SNOW 1.6 Band_1250
      //
      if(Len<5) {
        WriteWarning("Incorrect syntax for :Nudge command",Options.noisy);break;
      }
      if(pModel->GetHRUGroup(s[4])==NULL) {
        WriteWarning("Invalid HRU group in :Nudge command",Options.noisy);break;
      }
      int     SVlayer,iSV,nHRUs;
      sv_type SVtype;
      double  val;
      SVtype  = pModel->GetStateVarInfo()->StringToSVType(s[2], SVlayer, true);
      iSV     = pModel->GetStateVarIndex(SVtype,SVlayer);
      nHRUs   = pModel->GetNumHRUs();

      if(!strcmp(s[1],"NUDGE_MULTIPLY")) {//----------------------------------
        for(int k=0;k<nHRUs;k++) {
          if(pModel->IsInHRUGroup(k,s[4])) {
            val=pModel->GetHydroUnit(k)->GetStateVarValue(iSV);
            val*=s_to_d(s[3]);
            pModel->GetHydroUnit(k)->SetStateVarValue(iSV,val);
          }
        }
      }
      else if(!strcmp(s[1],"NUDGE_ADD")) {//----------------------------------
        for(int k=0;k<nHRUs;k++) {
          if(pModel->IsInHRUGroup(k,s[4])) {
            val=pModel->GetHydroUnit(k)->GetStateVarValue(iSV);
            val+=s_to_d(s[3]);
            pModel->GetHydroUnit(k)->SetStateVarValue(iSV,val);
          }
        }
      }
      break;
    }
    case(12):  //----------------------------------------------
    {
      /*
      :BasinTransportVariables [constit_name]
        :BasinIndex ID
          :ChannelMass [val]
          :RivuletMass [val]
          :Mout [nsegs]     [aMout     x nsegs    ] [aMoutLast]
          :Mlat [nQlatHist] [aMlatHist x nQlatHist] [MlatLast]
          :Min  [nQinHist]  [aMinHist  x nQinHist ]
          {:ResMassOut [Mout_res] [last_Mout_res]}
          {:ResMass [mass] [last mass]}
          {:ResSedMass [sed mass] [sed last mass]}
          {:BedTemperature [bed temp]}
        :BasinIndex ID
        ...
      :EndBasinTransportVariables
      */
      int c=pModel->GetTransportModel()->GetConstituentIndex(concname);
      if(c==DOESNT_EXIST) {
        WriteWarning("Unrecognized constituent entry in :BasinTransportVariables within .rvc file. Command was ignored.",Options.noisy);
        break;
      }
      CConstituentModel *pConstit=pModel->GetTransportModel()->GetConstituentModel(c);

      int p=-1;
      if(Options.noisy) { cout <<"   Reading Basin Transport Variables for constituent "<<pConstit->GetName()<<endl; }
      bool done=false;
      CSubBasin *pBasin=NULL;
      SBID=DOESNT_EXIST;
      do
      {
        bool eof=pp->Tokenize(s,Len);
        if(eof) { done=true; }
        if(Len==0) {}//do nothing
        else if(IsComment(s[0],Len)) {}//do nothing
        else if(!strcmp(s[0],":BasinIndex"))
        {
          SBID=s_to_ll(s[1]);
          pBasin=pModel->GetSubBasinByID(SBID);
          ExitGracefullyIf(pBasin==NULL,
            "ParseInitialConditionsFile: bad basin index in .rvc file",BAD_DATA);
          p++;
          if(Options.noisy) { cout <<"     Reading Transport Vars for Basin "<<pBasin->GetID()<<": "<<pBasin->GetName()<<endl; }
        }
        else if(!strcmp(s[0],":ChannelMass"))
        {
          if(Len>=2) {
            pConstit->SetChannelMass(p,s_to_d(s[1]));
          }
        }
        else if(!strcmp(s[0],":RivuletMass"))
        {
          if(Len>=2) {
            pConstit->SetRivuletMass(p,s_to_d(s[1]));
          }
        }
        else if(!strcmp(s[0],":Mout"))
        {
          if(Len>2) {
            int nsegs=s_to_i(s[1]);
            if(Len>=nsegs+3) {
              double *aMout=new double[nsegs+1];
              for(int i=0;i<=nsegs;i++) {
                aMout[i]=s_to_d(s[i+2]);
              }
              pConstit->SetMoutArray(p,nsegs,aMout,aMout[nsegs]);
              delete[] aMout;
            }
            else {
              WriteWarning("ParseInitialConditionsFile: incorrect number of terms in :Mout initial conditions item",Options.noisy);
            }
          }
        }
        else if(!strcmp(s[0],":Mlat"))
        {
          if(Len>2) {
            int histsize=s_to_i(s[1]);
            if(Len>=histsize+3) {
              double *aMlat=new double[histsize+1];
              for(int i=0;i<=histsize;i++) {
                aMlat[i]=s_to_d(s[i+2]);
              }
              pConstit->SetMlatHist(p,histsize,aMlat,aMlat[histsize]);
              delete[] aMlat;
            }
            else {
              WriteWarning("ParseInitialConditionsFile: incorrect number of terms in :Mlat initial conditions item",Options.noisy);
            }
          }
        }
        else if(!strcmp(s[0],":Min"))
        {
          if(Len>2) {
            int histsize=s_to_i(s[1]);
            if(Len>=histsize+2) {
              double *aMin=new double[histsize];
              for(int i=0;i<histsize;i++) {
                aMin[i]=s_to_d(s[i+2]);
              }
              pConstit->SetMinHist(p,histsize,aMin);
              delete[] aMin;
            }
            else {
              WriteWarning("ParseInitialConditionsFile: incorrect number of terms in :Min initial conditions item",Options.noisy);
            }
          }
        }
        else if(!strcmp(s[0],":ResMass"))
        {
          if(Len>=3) {
            pConstit->SetInitialReservoirMass(p,s_to_d(s[1]),s_to_d(s[2]));
          }
        }
        else if(!strcmp(s[0],":ResMassOut"))
        {
          if(Len>=3) {
            pConstit->SetReservoirMassOutflow(p,s_to_d(s[1]),s_to_d(s[2]));
          }
        }
        else if (!strcmp(s[0], ":ResSedMass"))
        {
          if (Len >= 3) {
            pConstit->SetInitReservoirSedMass(p, s_to_d(s[1]), s_to_d(s[2]));
          }
        }
        else if (!strcmp(s[0], ":BedTemperature"))
        {
          if (Len >= 2) {
            CEnthalpyModel *pEnth=(CEnthalpyModel*)(pConstit);
            pEnth->SetBedTemperature (p,s_to_d(s[1]));
          }
        }
        else if(!strcmp(s[0],":EndBasinTransportVariables"))
        {
          done=true;
        }
      } while(!done);
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
        else
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
      double maxv=max(pHRU->GetStateVarMax(i,v,Options,true),0.0); //ignores all variable maximum thresholds that are dependent upon model state
      if ((v[i]-maxv)>PRETTY_SMALL)// check for capacity
      {
        string name = CStateVariable::GetStateVarLongName(pModel->GetStateVarType(i),
                                                          pModel->GetStateVarLayer(i),
                                                          pModel->GetTransportModel());
        string warn ="maximum state variable limit exceeded in initial conditions for " + name+ " (in HRU "+to_string(pHRU->GetHRUID())+") in .rvc file";
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

void SetInitialStateVar(CModel *&pModel,const int SVind,const sv_type typ,const int m,const int k,const double &val)
{
  if(typ!=CONSTITUENT) //native units, no problem
  {
    pModel->GetHydroUnit(k)->SetStateVarValue(SVind,val);
  }
  else { // specifying mg/L or degrees C, must be converted to mg/m2 or MJ/m2
    string name =pModel->GetTransportModel()->GetConstituentTypeName(m);
    int       c=pModel->GetTransportModel()->GetConstituentIndex(name);

    if(c==DOESNT_EXIST) {
      string warn="Constituent '"+name+"' in .rvc file does not exist in the model.";
      WriteWarning(warn.c_str(), true);
      return;
    }

    int   iStor=pModel->GetTransportModel()->GetWaterStorIndexFromLayer(m);
    double  vol=pModel->GetHydroUnit(k)->GetStateVarValue(iStor); //[mm] Assumes this has already been initialized (this will be true for .rvc file)

    double conc=pModel->GetTransportModel()->GetConstituentModel(c)->ConvertConcentration(val); //(mg/L or o/oo or C) -> mg/mm-m2 or MJ/mm-m2
    double mass=conc*vol; //mg/mm-m2*mm ->mg/m2
    //specified in mg/L, convert to mg/m2
    // OR
    // specified in o/oo, convert to mg/m2
    // OR
    // specified in degC, convert to mJ/m2
    pModel->GetHydroUnit(k)->SetStateVarValue(SVind,mass);
  }
}
