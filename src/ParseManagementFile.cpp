/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2024 the Raven Development Team
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Model.h"
//#include "StateVariables.h"
#include "ParseLib.h"
#include "DemandOptimization.h"

//////////////////////////////////////////////////////////////////
/// \brief Parses Demand Management file
/// \details model.rvm: input file that defines management optimization problem and solution options 
///
/// \param *&pModel [out] Reference to model object
/// \param &Options [out] Global model options information
/// \return True if operation is successful
//
bool ParseManagementFile(CModel *&pModel,const optStruct &Options)
{
  int         i;              //counters
  bool        ended(false);
  bool        in_ifmode_statement=false;

  ifstream    RVM;
  RVM.open(Options.rvm_filename.c_str());
  if(RVM.fail()) {
    cout << "ERROR opening model management file: "<<Options.rvm_filename<<endl; return false;
  }

  int   Len,line(0),code;
  char *s[MAXINPUTITEMS];
  CParser *pp=new CParser(RVM,Options.rve_filename,line);

  ifstream INPUT2;           //For Secondary input
  CParser *pMainParser=NULL; //for storage of main parser while reading secondary files

  CDemandOptimizer *pDO=pDO=pModel->GetDemandOptimizer();

  pDO->Initialize(pModel,Options); //only requires rvh read

  if(Options.noisy) {
    cout <<"======================================================"<<endl;
    cout <<"Parsing Model Management File " << Options.rvm_filename <<"..."<<endl;
    cout <<"======================================================"<<endl;
  }
  string firstword;
  //--Sift through file-----------------------------------------------
  firstword=pp->Peek();
  if (firstword == ":DefineDecisionVariable") {
    pp->NextIsMathExp();
  }
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
    bool is_goal=false;
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
    else if(!strcmp(s[0],":LookbackDuration"))            { code=1;  }
    else if(!strcmp(s[0],":DemandGroup"))                 { code=11; }
    else if(!strcmp(s[0],":DemandMultiplier"))            { code=12; }
    else if(!strcmp(s[0],":DemandGroupMultiplier"))       { code=13; }
    else if(!strcmp(s[0],":DemandPenalty"))               { code=14; }
    else if(!strcmp(s[0],":DemandResetDate"))             { code=15; }
    else if(!strcmp(s[0],":DemandPriority"))              { code=16; }
    else if(!strcmp(s[0],":NamedConstant"))               { code=20; }
    else if(!strcmp(s[0],":DefineDecisionVariable"))      { code=21; }
    else if(!strcmp(s[0],":DecisionVariableBounds"))      { code=22; }
    else if(!strcmp(s[0],":ManagementConstraint"))        { code=23; is_goal=false; }
    else if(!strcmp(s[0],":ManagementGoal"))              { code=23; is_goal=true;}
    else if(!strcmp(s[0],":LookupTable"))                 { code=30; }

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
      for(i=1;i<Len;i++) { filename+=s[i]; if(i<Len-1) { filename+=' '; } }
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
        for(i=1; i<Len; i++) {
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
    {/*:LookbackDuration [duration, in days]*/
      if(Options.noisy) { cout <<"LookbackDuration"<<endl; }
      double duration=s_to_d(s[1]);
      int n=static_cast <int>(rvn_floor(duration/Options.timestep+REAL_SMALL));
      pDO->SetHistoryLength(n);
      break;
    }
    case(11):  //----------------------------------------------
    {/*:DemandGroup [groupname] 
         [demand1] [demand2] ... [demandN]
       :EndDemandGroup     
     */
      if(Options.noisy) { cout <<"Demand Group"<<endl; }
      //pDO->AddDemandGroup(s[1]);
      //...
      break;
    }
    case(12):  //----------------------------------------------
    {/*:DemandMultiplier [demand] [multiplier] */
      if(Options.noisy) { cout <<"Demand Multiplier"<<endl; }
      double mult = s_to_d(s[2]);
      //pDO->MultiplyDemand(s[1],mult);
      break;
    }
    case(13):  //----------------------------------------------
    {/*:DemandGroupMultiplier [groupname] [multiplier] */
      if(Options.noisy) { cout <<"Demand Group Multiplier"<<endl; }
      double mult = s_to_d(s[2]);
      //pDO->MultiplyGroupDemand(s[1],mult);
      break;
    }
    case(14):  //----------------------------------------------
    {/*:DemandPenalty [demand1] [penalty] */
      if(Options.noisy) { cout <<"Demand Penalty"<<endl; }
      //pDO->SetDemandPenalty(s[1],s_to_d(s[2]);
      break;
    }
    case(15):  //----------------------------------------------
    {/*:DemandResetDate [demand1] [julian day] */
      if(Options.noisy) { cout <<"Demand Reset Date"<<endl; }
      pDO->SetCumulativeDate(s_to_i(s[2]), s[1]);
      break;
    }
    case(16):  //----------------------------------------------
    {/*:DemandPriority [demand1] [priority] */
      //pDO->SetDemandPriority(s[1], s_to_d(s[2]));
      break;
    }
    case(20):  //----------------------------------------------
    {/*:NamedConstant [name] [value] */
      if(Options.noisy) { cout <<"Demand Penalty"<<endl; }
      pDO->AddUserConstant(s[1],s_to_d(s[2]));
      break;
    }
    case(21):  //----------------------------------------------
    { /*:DefineDecisionVariable [name] = [expressionRHS] */
      if(Options.noisy) { cout <<"Define Decision Variable"<<endl; }
      expressionStruct *pExp;      
      dv_constraint    *pConst=NULL;
      decision_var     *pDV = new decision_var(s[1],DOESNT_EXIST,DV_USER);
      
      pDO->AddDecisionVar(pDV);
      pExp=pDO->ParseExpression((const char**)(s),Len,pp->GetLineNumber(),pp->GetFilename());
      if (pExp!=NULL){
        pConst = pDO->AddConstraint(s[1], pExp, false);
      } 
      else {
        string warn ="Invalid expression in :DefineDecisionVariable command at line " + pp->GetLineNumber();
        WriteWarning(warn.c_str(),Options.noisy);
      }
      break;
    }
    case(22):  //----------------------------------------------
    {/*:DecisionVariableBounds [name] [low] [upp] */
      if(Options.noisy) { cout <<"Decision variable bounds"<<endl; }
      pDO->SetDecisionVarBounds(s[1], s_to_d(s[2]), s_to_d(s[3]));
      break;
    }
    case(23):  //----------------------------------------------
    {/*:ManagementConstraint [name] 
          :Expresion [expression]
          :Conditional [condition]
          :Conditional [condition]
       :EndManagementConstraint
       or
       :ManagementGoal [name] 
          :Expresion [expression]
          :Conditional [condition]
          :Conditional [condition]
          :Penalty [value] {value2} 
       :EndManagementGoal
     */
      if(Options.noisy) { cout <<"Management Constraint or Management Goal"<<endl; }
      string name=s[1];
      expressionStruct *pExp;      
      dv_constraint    *pConst=NULL;
      firstword=pp->Peek();
      if (firstword == ":Expression") {
        pp->NextIsMathExp();
      }

      while(!pp->Tokenize(s,Len))
      {
        if(Options.noisy) { cout << "-->reading line " << pp->GetLineNumber() << ": "; }
        if     (Len == 0)            { if(Options.noisy) { cout << "#" << endl; } }//Do nothing
        else if(IsComment(s[0],Len)) { if(Options.noisy) { cout << "#" << endl; } }
        else if(!strcmp(s[0], ":Expression")) 
        {
          pExp=pDO->ParseExpression((const char**)(s),Len,pp->GetLineNumber(),pp->GetFilename());
          if (pExp!=NULL){
            pConst=pDO->AddConstraint(name,pExp,is_goal);
          }
          else {
            string warn ="Invalid expression in :Expression command at line " + pp->GetLineNumber();
            WriteWarning(warn.c_str(),Options.noisy);
          }
        }
        else if (!strcmp(s[0], ":Condition")) 
        {
          //TODO: handle more complex, expression-based conditions 
          if (pConst!=NULL){
            dv_condition *pCond = new dv_condition();
            pCond->dv_name=s[1];
            if      (!strcmp(s[2],"IS_BETWEEN"     )){pCond->compare=COMPARE_BETWEEN;}
            else if (!strcmp(s[2],"IS_GREATER_THAN")){pCond->compare=COMPARE_GREATERTHAN;}
            else if (!strcmp(s[2],"IS_LESS_THAN"   )){pCond->compare=COMPARE_LESSTHAN;}
            pCond->value=s_to_d(s[3]);
            if (Len>=5){
              pCond->value2 = s_to_d(s[4]);
            }
            pConst->AddCondition(pCond);
          } 
          else{
            ExitGracefully("ParseManagementFile: :Condition statement must appear after valid :Expression in :ManagementConstraint command",BAD_DATA_WARN);
          }
        } 
        else if (!strcmp(s[0], ":Penalty")) {
          pConst->penalty_under = s_to_d(s[1]);
          pConst->penalty_over = s_to_d(s[1]);
          if (Len >= 3) {
            pConst->penalty_over = s_to_d(s[2]);
          }
        }
        else if (!strcmp(s[0], ":Priority")) {
          pConst->priority = s_to_i(s[1]); //for later 
        }
        else if (!strcmp(s[0], ":EndManagementGoal")) {
          break;
        }
        else if (!strcmp(s[0], ":EndManagementConstraint")) {
          break;
        } 
        else {
          WriteWarning("ParseManagementFile: Unrecognized command in :ManagementConstraint command block",Options.noisy);
        }
        firstword=pp->Peek();
        if (firstword == ":Expression") {
          pp->NextIsMathExp();
        }
      }
      break;
    }
    case(30):  //----------------------------------------------
    {/*:LookupTable [name]  
         N 
         x1 y1 
         x2 y2
         ...
         xN yN 
       :EndLookupTable
     */
      if(Options.noisy) { cout <<"LookupTable"<<endl; }
      
      string name = s[1];

      pp->Tokenize(s,Len);
      int N = s_to_i(s[0]);

      double *aX = new double[N];
      double *aY = new double[N];
      for(int i = 0; i < N; i++) 
      {
        pp->Tokenize(s,Len);
        if(IsComment(s[0],Len)) { i--; }
        else {
          if(Len>=2) {
            aX[i] = s_to_d(s[0]);  
            aY[i] = s_to_d(s[1]);
          }
          else {
            WriteWarning("Incorrect line length (<2) in :LookupTable command",Options.noisy);
          }
        }
      }
      pp->Tokenize(s,Len); //:EndLookupTable

      CLookupTable *pLUT = new CLookupTable(name,aX,aY,N);
      pDO->AddUserLookupTable(pLUT);
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
          string errString = "Unrecognized command in .rvm file:\n   " + string(s[0]);
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
  RVM.close();

  pDO->InitializePostRVMRead(pModel,Options);

  delete pp;
  pp=NULL;

  return true;
}