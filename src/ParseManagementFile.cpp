/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2025 the Raven Development Team
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Model.h"
#include "ParseLib.h"
#include "DemandOptimization.h"
#include "SubBasin.h"
#include "Demands.h"
#include "LookupTable.h"

void SummarizeExpression(const char **s, const int Len, expressionStruct* exp); //defined in DemandExpressionHandling.cpp

void swapWildcards(const char **s,const int Len, string **aWildcards,const int &nWildcards)
{
  //return; //
  if (nWildcards==0){return;}
  if (!strcmp(s[0],":LoopVector")){return;}
  if (!strcmp(s[0],":EndLoopThrough")){return;}

  //cout<<"PRE-SWAP:  ";for (int i=0;i<Len;i++){cout<<s[i]<<" ";} cout<<endl;
  string tmp,tmp2,wc;
  const char* tmp3;
  for (int i=0;i<Len;i++)
  {
    tmp = s[i];
    for (int j=0;j<nWildcards;j++)
    {
      wc=aWildcards[j][0];
      size_t ind = tmp.find(wc);
      if (ind!=string::npos)
      {
        tmp2 = tmp.substr(0,ind)+aWildcards[j][1]+tmp.substr(ind+wc.length(),string::npos)+'\0';
        tmp3 = tmp2.c_str();
        memcpy((void*)s[i],tmp3,tmp2.length());
      }
    }
  }
  //cout<<"POST-SWAP: ";for (int i=0;i<Len;i++){cout<<s[i]<<" ";} cout<<endl;
}

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

  CDemand    *pDemand=NULL;
  int         demand_ind=0;
  long long   demandSBID;
  long long   demand_ID;
  string      demand_name;

  ifstream    INPUT2;                //For Secondary input
  CParser    *pMainParser=NULL;      //for storage of main parser while reading secondary files
  ifstream    INPUT3;                //For tertiary input
  CParser    *pSecondaryParser=NULL; //for storage of secondary parser while reading tertiary files

  ifstream    RVM;
  RVM.open(Options.rvm_filename.c_str(),ios::binary);
  if(RVM.fail()) {
    cout << "ERROR opening model management file: "<<Options.rvm_filename<<endl; return false;
  }

  string  firstword;
  int     Len,line(0),code;
  char   *s[MAXINPUTITEMS];
  CParser *pp=new CParser(RVM,Options.rvm_filename,line);

  // - looping variables --------------------------------------------
  const int       MAX_WILDCARDS=10;

  int             loopCount=0;
  streampos       loopStartPos;
  int             loopStartLine;
  CSubbasinGroup *pLoopSBGroup=NULL;
  CDemandGroup   *pLoopDemandGroup=NULL;
  int             loopListSize=0;
  string **aWildcards  = new string *[MAX_WILDCARDS  ];
  string **aLoopVector = new string *[MAX_WILDCARDS  ];
  for (int i = 0; i < MAX_WILDCARDS; i++) {
    aWildcards [i]=new string [2];
    aWildcards [i][0] = "";
    aWildcards [i][1] = "";
    aLoopVector[i]=NULL;
  }
  int    nWildcards=0;

  CDemandOptimizer *pDO=pModel->GetManagementOptimizer();

  pDO->Initialize(pModel,Options); //only requires rvh,.rvt read

  if(Options.noisy) {
    cout <<"======================================================"<<endl;
    cout <<"Parsing Model Management File " << Options.rvm_filename <<"..."<<endl;
    cout <<"======================================================"<<endl;
  }

  //--Sift through file-----------------------------------------------
  firstword=pp->Peek();
  if ((firstword == ":DefineDecisionVariable") || (firstword == ":DemandExpression") || (firstword == ":ReturnExpression") || (firstword == ":DefineWorkflowVariable"))
  {
    pp->NextIsMathExp();
  }
  bool end_of_file=pp->Tokenize(s,Len);

  while(!end_of_file)
  {
    swapWildcards((const char**)(s),Len,aWildcards,nWildcards);

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

    //-------------------MODEL MANAGEMENT PARAMETERS------------
    else if(!strcmp(s[0],":LookbackDuration"))            { code=1;  }
    else if(!strcmp(s[0],":DebugLevel"))                  { code=2;  }
    else if(!strcmp(s[0],":DemandGroup"))                 { code=11; }
    else if(!strcmp(s[0],":DemandMultiplier"))            { code=12; }
    else if(!strcmp(s[0],":DemandGroupMultiplier"))       { code=13; }
    else if(!strcmp(s[0],":DemandPenalty"))               { code=14; }
    else if(!strcmp(s[0],":DemandResetDate"))             { code=15; }
    else if(!strcmp(s[0],":DemandPriority"))              { code=16; }
    else if(!strcmp(s[0],":DemandIsUnrestricted"))        { code=17; }

    else if(!strcmp(s[0],":NamedConstant"))               { code=20; }
    else if(!strcmp(s[0],":DefineDecisionVariable"))      { code=21; }
    else if(!strcmp(s[0],":DecisionVariableBounds"))      { code=22; }
    else if(!strcmp(s[0],":ManagementConstraint"))        { code=23; is_goal=false; }
    else if(!strcmp(s[0],":ManagementGoal"))              { code=23; is_goal=true;  }
    else if(!strcmp(s[0],":DeclareDecisionVariable"))     { code=24; }
    else if(!strcmp(s[0],":DefineWorkflowVariable"))      { code=25; }
    else if(!strcmp(s[0],":WorkflowVarDefinition"))       { code=26; }
    else if(!strcmp(s[0],":DeclareWorkflowVariable"))     { code=27; }
    else if(!strcmp(s[0],":LookupTable"))                 { code=30; }
    else if(!strcmp(s[0],":OverrideStageDischargeCurve")) { code=31; }
    else if(!strcmp(s[0],":LookupTableFromCSV"))          { code=32; }

    else if(!strcmp(s[0],":LoopThrough"))                 { code=40; }
    else if(!strcmp(s[0],":EndLoopThrough"))              { code=41; }
    else if(!strcmp(s[0],":LoopVector"))                  { code=42; }

    else if(!strcmp(s[0],":WaterDemand"))                 { code=50; }
    else if(!strcmp(s[0],":EndWaterDemand"))              { code=51; }
    else if(!strcmp(s[0],":DemandTimeSeries"))            { code=52; }
    else if(!strcmp(s[0],":ReturnTimeSeries"))            { code=53; }
    else if(!strcmp(s[0],":ReturnDestination"))           { code=54; }
    else if(!strcmp(s[0],":IrrigationDestination"))       { code=55; }
    else if(!strcmp(s[0],":ResetDate"))                   { code=56; }
    else if(!strcmp(s[0],":Penalty"))                     { code=57; }
    else if(!strcmp(s[0],":ReturnFraction"))              { code=58; }
    else if(!strcmp(s[0],":ReturnExpression"))            { code=59; }
    else if(!strcmp(s[0],":DemandFlowFraction"))          { code=60; }
    else if(!strcmp(s[0],":DemandLookupTable"))           { code=61; }
    else if(!strcmp(s[0],":DemandExpression"))            { code=62; }
    //else if(!strcmp(s[0],":AnnualLicense"))             { code=63; }
    else if(!strcmp(s[0],":ReservoirWaterDemand"))        { code=64; }
    else if(!strcmp(s[0],":EndReservoirWaterDemand"))     { code=65; }
    else if(!strcmp(s[0],":IsUnrestricted"))              { code=66; }

    else if(!strcmp(s[0],":UserTimeSeries"))              { code=70; }

    else if(!strcmp(s[0],":NonlinearVariable"))           { code=80; }
    else if(!strcmp(s[0],":MaxNLSolverIterations"))       { code=81; }
    else if(!strcmp(s[0],":NLSolverTolerance"))           { code=82; }
    else if(!strcmp(s[0],":NLRelaxationCoeff"))           { code=83; }

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

      filename=CorrectForRelativePath(filename,Options.rvm_filename);

      if (pSecondaryParser != NULL){
        ExitGracefully("ParseEnsembleFile::nested :RedirectToFile commands are not allowed to be nested more than two levels (e.g., rvm file to rvm file to rvm file to rvm file)",BAD_DATA);
      }
      if (pMainParser == NULL) {
        INPUT2.open(filename.c_str(),ios::binary); //binary enables tellg() to work correctly for Unix files in parseLib::Peek()
        if(INPUT2.fail()) {
          string warn;
          warn=":RedirectToFile (from .rvm): Cannot find file "+filename;
          ExitGracefully(warn.c_str(),BAD_DATA);
        }
        pMainParser=pp;
        pp=new CParser(INPUT2,filename,line);//open new parser
      } //from base .rvm file
      else {
        INPUT3.open(filename.c_str(),ios::binary); //binary enables tellg() to work correctly for Unix files in parseLib::Peek()
        if(INPUT3.fail()) {
          string warn;
          warn=":RedirectToFile (from .rvm): Cannot find file "+filename;
          ExitGracefully(warn.c_str(),BAD_DATA);
        }
        pSecondaryParser=pp;
        pp=new CParser(INPUT3,filename,line);//open new parser
      } //from already redirected .rvm file
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
    case(2):  //----------------------------------------------
    {/*:DebugLevel [int level]*/
      if(Options.noisy) { cout <<":DebugLevel"<<endl; }
      pDO->SetDebugLevel(s_to_i(s[1]));
      break;
    }
    case(11):  //----------------------------------------------
    {/*:DemandGroup [groupname]
         [demand1] [demand2] ... [demandN]
       :EndDemandGroup
     */
      if(Options.noisy) { cout <<"Demand Group"<<endl; }
      /*if (pDO->DemandsAreInitialized()) {
        ExitGracefully("ParseManagementFile: all :Demand commands must be specified before management goals/constraints/decision vars in the .rvm file",BAD_DATA_WARN);
      }*/
      if (Len!=2){pp->ImproperFormat(s);}

      pDO->AddDemandGroup(s[1]);
      CDemandGroup *pGrp;
      pGrp=pDO->GetDemandGroup(pDO->GetNumDemandGroups()-1);

      bool eof=false;
      while ( (Len==0) || (strcmp(s[0],":EndDemandGroup")) )
      {
        eof=pp->Tokenize(s,Len);
        if(eof) { break; }
        if      (IsComment(s[0], Len)){}//comment line
        else if (!strcmp(s[0],":EndDemandGroup")){}//done
        else if (s[0][0] == ':') {
          string warn="ParseManagementFile: Command found between :DemandGroup...:EndDemandGroup commands at line "+to_string(pp->GetLineNumber())+" of .rvm file. Only demand IDs should be between these two commands";
          ExitGracefully(warn.c_str(),BAD_DATA_WARN);
        }
        else
        {
          for(i=0;i<Len;i++)
          {
            int d=pDO->GetDemandIndexFromName(s[i]);
            if (d!=DOESNT_EXIST){
              pGrp->AddDemand(pDO->GetWaterDemand(d));
            }
            else {
              string warn="ParseManagementFile: invalid or disabled demand in :DemandGroup command. Will be ignored";
              WriteWarning(warn.c_str(),Options.noisy);
            }
          }
        }
      }
      break;
    }
    case(12):  //----------------------------------------------
    {/*:DemandMultiplier [demand] [multiplier] */
      if(Options.noisy) { cout <<"Demand Multiplier"<<endl; }
      double mult = s_to_d(s[2]);
      if (mult < 0) {
        ExitGracefully("ParseManagementFile: :DemandMultiplier cannot be negative",BAD_DATA_WARN);
      }
      int d=pDO->GetDemandIndexFromName(s[1]);
      if (d == DOESNT_EXIST) {
        WriteWarning("ParseManagementFile: invalid or disabled demand index in :DemandMultiplier command. Will be ignored",Options.noisy);
      }
      else{
        pDO->GetWaterDemand(d)->SetMultiplier(mult);
      }
      break;
    }
    case(13):  //----------------------------------------------
    {/*:DemandGroupMultiplier [groupname] [multiplier] */
      if(Options.noisy) { cout <<"Demand Group Multiplier"<<endl; }
      double mult = s_to_d(s[2]);
      if (mult < 0) {
        ExitGracefully("ParseManagementFile: :DemandGroupMultiplier cannot be negative",BAD_DATA_WARN);
      }
      CDemandGroup *pDG=pDO->GetDemandGroupFromName(s[1]);
      if (pDG == NULL) {
        WriteWarning("ParseManagementFile: invalid group name in :DemandGroupMultiplier command. Will be ignored",Options.noisy);
      }
      else {
        for (int ii = 0; ii < pDG->GetNumDemands(); ii++) {
          pDG->GetDemand(ii)->SetMultiplier(s_to_d(s[2]));
        }
      }
      break;
    }
    case(14):  //----------------------------------------------
    {/*:DemandPenalty [demand1] [penalty] */ //demand1 is 1234 or FarmerBob, not !D1234 or !D.FarmerBob
      if(Options.noisy) { cout <<"Demand Penalty"<<endl; }
      pDO->SetDemandPenalty(s[1],s_to_d(s[2]));
      break;
    }
    case(15):  //----------------------------------------------
    {/*:DemandResetDate [demand1] [julian day] */ //demand1 is 1234 or FarmerBob, not !D1234 or !D.FarmerBob
      if(Options.noisy) { cout <<"Demand Reset Date"<<endl; }
      pDO->SetCumulativeDate(s_to_i(s[2]), s[1]);
      break;
    }
    case(16):  //----------------------------------------------
    {/*:DemandPriority [demand1] [priority] */ //demand1 is 1234 or FarmerBob, not !D1234 or !D.FarmerBob
      if(Options.noisy) { cout <<"Demand Priority "<<endl; }
      //pDO->SetDemandPriority(s[1], s_to_d(s[2]));
      break;
    }
    case(17):  //----------------------------------------------
    {/*:DemandIsUnrestricted [demand1] */ //demand1 is 1234 or FarmerBob, not !D1234 or !D.FarmerBob
      if(Options.noisy) { cout <<"Unrestricted Demand"<<endl; }
      pDO->SetDemandAsUnrestricted(s[1]);
      break;
    }
    case(20):  //----------------------------------------------
    {/*:NamedConstant [name] [value] */
      if(Options.noisy) { cout <<"Named Constant "<<endl; }
      if (Len > 3) {
        WriteWarning("The :NamedConstant command should have only two arguments. Command will be ignored",Options.noisy);
      }
      else{
        pDO->AddUserConstant(s[1],s_to_d(s[2]));
      }
      break;
    }
    case(21):  //----------------------------------------------
    { /*:DefineDecisionVariable [name] = [expressionRHS] */
      if(Options.noisy) { cout <<"Define Decision Variable"<<endl; }
      expressionStruct *pExp;
      managementGoal    *pConst=NULL;
      decision_var     *pDV = new decision_var(s[1],DOESNT_EXIST,DV_USER,pDO->GetNumUserDVs());

      pDO->AddUserDecisionVar(pDV);
      pExp=pDO->ParseExpression((const char**)(s),Len,pp->GetLineNumber(),pp->GetFilename());

      if (pDO->GetDebugLevel()>=1){
        SummarizeExpression((const char**)(s),Len,pExp);
      }

      if (pExp!=NULL){
        pConst=new managementGoal();
        pConst->name=s[1];
        pConst->is_goal=false;
        pConst->AddExpression(pExp);
        pDO->AddGoalOrConstraint(pConst);
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
      pDO->SetUserDecisionVarBounds(s[1], s_to_d(s[2]), s_to_d(s[3]));
      break;
    }
    case(23):  //----------------------------------------------
    {/*:ManagementConstraint [name]
         :OperatingRegime A
           :Expression [expression]
           :Condition [condition]
           :Condition [condition]
         :EndOperatingRegime
         :OperatingRegime B
           :Expression [expression]
         :EndOperatingRegime
       :EndManagementConstraint
       or
       :ManagementGoal [name]
         :OperatingRegime A
          :Expression [expression]
          :Condition [condition]
          :Condition [condition]
          :Penalty [value] {value2} #specific to goal
         :EndOperatingRegime
         :Penalty [value] {value2} #default penalty if outside op block
       :EndManagementGoal
       or
       :ManagementGoal [name]
         :Expression [expression]
         :Penalty [value] {value2}
         :UseStageUnitsCorrection [SBID]
         :OverrideStageDischargeCurve [SBID]
       :EndManagementGoal
     */
      if(Options.noisy) { cout <<"Management Constraint or Management Goal"<<endl; }

      managementGoal *pGoal=new managementGoal();
      pGoal->name=s[1];
      pGoal->is_goal=is_goal;

      expressionStruct *pExp;
      bool is_first=true;
      bool in_op_block=false;
      firstword=pp->Peek();
      if (firstword == ":Expression") {pp->NextIsMathExp();}
      if (firstword == ":Condition")  {pp->NextIsMathExp();}

      while(!pp->Tokenize(s,Len))
      {
        swapWildcards((const char**)(s),Len,aWildcards,nWildcards);

        if(Options.noisy) { cout << "-->reading line " << pp->GetLineNumber() << ": "; }

        if     (Len == 0)             { if(Options.noisy) { cout << "#" << endl; } }//Do nothing
        else if (IsComment(s[0],Len)) { if(Options.noisy) { cout << "#" << endl; } }
        //----------------------------------------------
        else if (!strcmp(s[0], ":OperatingRegime"))
        {
          if (Options.noisy){cout<<" Operating regime "<<endl; }
          if (Len < 2) {
            ExitGracefully("ParseManagementFile: OperatingRegime name missing.",BAD_DATA_WARN);
          }
          op_regime *pOR=new op_regime(s[1]);//[OPUPDATE]
          pGoal->AddOperatingRegime(pOR,is_first);
          is_first=false;
          in_op_block=true;
        }
        //----------------------------------------------
        else if (!strcmp(s[0], ":EndOperatingRegime"))
        {
          if (Options.noisy){cout<<" End operating regime "<<endl; }
          in_op_block=false;
        }
        //----------------------------------------------
        else if(!strcmp(s[0], ":Expression"))
        {
          if (Options.noisy){cout<<" Expression "<<endl; }
          if (pGoal->GetCurrentExpression() != NULL) {
            ExitGracefully("ParseManagementFile: only one :Expression allowed in each :OperatingRegime command block (or only one if no :OperatingRegime blocks used).",BAD_DATA_WARN);
            break;
          }
          pExp=pDO->ParseExpression((const char**)(s),Len,pp->GetLineNumber(),pp->GetFilename());
          if (pExp!=NULL){
            pGoal->AddExpression(pExp);
          }
          else {
            string warn ="Invalid expression in :Expression command at line " + to_string(pp->GetLineNumber());
            WriteWarning(warn.c_str(),Options.noisy);
          }
          if (pDO->GetDebugLevel()>=1){
            SummarizeExpression((const char**)(s),Len,pExp);
          }
        }
        //----------------------------------------------
        else if (!strcmp(s[0], ":Condition"))
        {
          if (Options.noisy){cout<<" Condition "<<endl; }
          if (pGoal==NULL){
            ExitGracefully("ParseManagementFile: :Condition statement must appear after valid :Expression in :ManagementConstraint command",BAD_DATA_WARN);
          }
          exp_condition *pCond;
          pCond=pDO->ParseCondition((const char**)(s),Len,pp->GetLineNumber(),pp->GetFilename());
          pGoal->AddOpCondition(pCond);
        }
        //----------------------------------------------
        else if (!strcmp(s[0], ":Penalty"))
        {
          if (Options.noisy){cout<<" Penalty "<<endl; }
          if (!pGoal->is_goal) {
            ExitGracefully("ParseManagementFile: :Penalty found within :ManagementConstraint command. It will be ignored, as all constraints have infinite penalty.",BAD_DATA_WARN);
          }
          if (in_op_block) {
            int k=pGoal->nOperRegimes;
            pGoal->pOperRegimes[k]->penalty_under = s_to_d(s[1]);
            pGoal->pOperRegimes[k]->penalty_over  = s_to_d(s[1]);
            if (Len >= 3) {
              pGoal->pOperRegimes[k]->penalty_over = s_to_d(s[2]);
            }
          }
          else{
            pGoal->penalty_under = s_to_d(s[1]);
            pGoal->penalty_over  = s_to_d(s[1]);
            if (Len >= 3) {
              pGoal->penalty_over = s_to_d(s[2]);
            }
          }
        }
        //----------------------------------------------
        else if (!strcmp(s[0], ":Priority"))
        {
          if (Options.noisy){cout<<" Priority "<<endl; }
          pGoal->priority = s_to_i(s[1]); //for later
        }
        //----------------------------------------------
        else if (!strcmp(s[0], ":UseStageUnitsCorrection"))
        {
          if (Options.noisy){cout<<" Use Stage Units Correction "<<endl; }
          CSubBasin *pSB=pModel->GetSubBasinByID(s_to_ll(s[1]));

          ExitGracefullyIf(pSB->GetGlobalIndex()==DOESNT_EXIST,"ParseManagementFile: subbasin ID in :UseStageUnitsCorrection is invalid",BAD_DATA_WARN);

          if (pSB->GetReservoir()==NULL){
            string advice="ParseManagementFile:The reservoir in subbasin "+to_string(pSB->GetID()) + " doesnt exist and cannot be used to calculate a stage units correction.";
            ExitGracefully(advice.c_str(), BAD_DATA_WARN);
          }
          else{
            int k=pSB->GetReservoir()->GetHRUIndex();
            if (k==DOESNT_EXIST){
                string advice="ParseManagementFile:The reservoir in subbasin "+to_string(pSB->GetID()) + " doesnt have an :HRUID and cannot be used to calculate a stage units correction.";
                ExitGracefully(advice.c_str(), BAD_DATA_WARN);
            }
            else {
              pGoal->use_stage_units=true;
              pGoal->reservoir_index=pSB->GetGlobalIndex();
            }
          }
        }
        //----------------------------------------------
        else if (!strcmp(s[0], ":EndManagementGoal")) {
          if (Options.noisy){cout<<endl; }
          break;
        }
        //----------------------------------------------
        else if (!strcmp(s[0], ":EndManagementConstraint")) {
          if (Options.noisy){cout<<endl; }
          break;
        }
        else {
          string warn="ParseManagementFile: Unrecognized command "+to_string(s[0]) + " in :ManagementConstraint command block";
          WriteWarning(warn.c_str(), Options.noisy);
        }
        firstword=pp->Peek();
        if (firstword == ":Expression") {pp->NextIsMathExp();}
        if (firstword == ":Condition")  {pp->NextIsMathExp();}
      }
      //any invalid expressions have to shut down simulation
      bool badgoal=false;
      for (int k = 0; k < pGoal->nOperRegimes; k++) {
        if (pGoal->pOperRegimes[k]->pExpression == NULL) { badgoal=true;}
      }
      if (!badgoal) {
        pDO->AddGoalOrConstraint(pGoal);
      }

      break;
    }
    case(24):  //----------------------------------------------
    { /*:DeclareDecisionVariable [name]  */
      if(Options.noisy) { cout <<"Declare Decision Variable"<<endl; }
      decision_var     *pDV = new decision_var(s[1],DOESNT_EXIST,DV_USER,pDO->GetNumUserDVs());
      pDO->AddUserDecisionVar(pDV);
      break;
    }
    case(25):  //----------------------------------------------
    { /*:DefineWorkflowVariable [name] = [expressionRHS] */
      if(Options.noisy) { cout <<"Define Workflow Variable"<<endl; }
      expressionStruct *pExp;
      workflowVar *pWV;

      if (pDO->GetWorkflowVarStruct(to_string(s[1])) != NULL){//already declared
        pWV=pDO->GetWorkflowVarStruct(to_string(s[1]));
      }
      else{
        pWV= new workflowVar();
        pWV->name = to_string(s[1]);
        pDO->AddWorkflowVariable(pWV);
      }

      pExp=pDO->ParseExpression((const char**)(s),Len,pp->GetLineNumber(),pp->GetFilename());

      if (pExp!=NULL){
        pWV->AddExpression(pExp);
      }
      else{
        string warn ="Invalid expression in :DefineWorkflowVariable command at line " + pp->GetLineNumber();
        WriteWarning(warn.c_str(),Options.noisy);
        break;
      }

      if (pDO->GetDebugLevel()>=1){
        SummarizeExpression((const char**)(s),Len,pExp);
      }

      break;
    }
    case(26):  //----------------------------------------------
    {/*:WorkflowVarDefinition [name]
         :OperatingRegime A
           :Expression [name] = [expression]
           :Condition [condition]
           :Condition [condition]
         :EndOperatingRegime
         :OperatingRegime B
           :Expression [expression]
         :EndOperatingRegime
       :EndWorkflowVarDefinition
       or
       :WorkflowVarDefinition [name]
         :Expression [name] = [expression]
       :EndWorkflowVarDefinition
     */
      if (Options.noisy) { cout << "Workflow Variable Definition Statement" << endl; }

      workflowVar* pWV;
      if (pDO->GetWorkflowVarStruct(to_string(s[1])) != NULL){//already declared
        pWV=pDO->GetWorkflowVarStruct(to_string(s[1]));
      }
      else{
        pWV= new workflowVar();
        pWV->name = to_string(s[1]);
        pDO->AddWorkflowVariable(pWV);
      }

      expressionStruct *pExp;
      bool is_first=true;
      bool in_op_block=false;
      firstword=pp->Peek();
      if (firstword == ":Expression") {pp->NextIsMathExp();}
      if (firstword == ":Condition")  {pp->NextIsMathExp();}

      while(!pp->Tokenize(s,Len))
      {
        swapWildcards((const char**)(s),Len,aWildcards,nWildcards);

        if(Options.noisy) { cout << "-->reading line " << pp->GetLineNumber() << ": "; }

        if     (Len == 0)             { if(Options.noisy) { cout << "#" << endl; } }//Do nothing
        else if (IsComment(s[0],Len)) { if(Options.noisy) { cout << "#" << endl; } }
        //----------------------------------------------
        else if (!strcmp(s[0], ":OperatingRegime"))
        {
          if (Options.noisy){cout<<" Operating regime "<<endl; }
          if (Len < 2) {
            ExitGracefully("ParseManagementFile: OperatingRegime name missing.",BAD_DATA_WARN);
          }
          op_regime *pOR=new op_regime(s[1]);
          pWV->AddOperatingRegime(pOR,is_first);
          is_first=false;
          in_op_block=true;
        }
        //----------------------------------------------
        else if (!strcmp(s[0], ":EndOperatingRegime"))
        {
          if (Options.noisy){cout<<" End operating regime "<<endl; }
          in_op_block=false;
        }
        //----------------------------------------------
        else if(!strcmp(s[0], ":Expression"))
        {
          if (Options.noisy){cout<<" Expression "<<endl; }
          if (pWV->GetCurrentExpression() != NULL) {
            ExitGracefully("ParseManagementFile: only one :Expression allowed in each :OperatingRegime command block (or only one if no :OperatingRegime blocks used).",BAD_DATA_WARN);
            break;
          }
          pExp=pDO->ParseExpression((const char**)(s),Len,pp->GetLineNumber(),pp->GetFilename());
          if (pExp!=NULL){
            pWV->AddExpression(pExp);
          }
          else {
            string warn ="Invalid expression in :Expression command at line " + to_string(pp->GetLineNumber());
            WriteWarning(warn.c_str(),Options.noisy);
          }
          if (pDO->GetDebugLevel()>=1){
            SummarizeExpression((const char**)(s),Len,pExp);
          }
        }
        //----------------------------------------------
        else if (!strcmp(s[0], ":Condition"))
        {
          if (Options.noisy){cout<<" Condition "<<endl; }

          exp_condition *pCond;
          if (pWV!=NULL){
            pCond=pDO->ParseCondition((const char**)(s),Len,pp->GetLineNumber(),pp->GetFilename());
            if (pCond != NULL) {
              pWV->AddOpCondition(pCond);
            }
          }
          else{
            ExitGracefully("ParseManagementFile: :Condition statement must appear after valid :Expression in :WorkflowVarDefinition command",BAD_DATA_WARN);
          }
        }
        //----------------------------------------------
        else if (!strcmp(s[0], ":EndWorkflowVarDefinition")) {
          if (Options.noisy){cout<<endl; }
          break;
        }
        else {
          WriteWarning("ParseManagementFile: Unrecognized command in :WorkflowVarDefinition command block",Options.noisy);
        }
        firstword=pp->Peek();
        if (firstword == ":Expression") {pp->NextIsMathExp();}
        if (firstword == ":Condition")  {pp->NextIsMathExp();}
      }
      //any invalid expressions have to shut down simulation
      bool baddef=false;
      for (int k = 0; k < pWV->nOperRegimes; k++) {
        if (pWV->pOperRegimes[k]->pExpression == NULL) { baddef=true;}
      }
      if (baddef) {
        ExitGracefully("ParseManagementFile: Bad or missing :Expression in one or more operating regimes within workflow variable definition",BAD_DATA_WARN);
      }
      break;
    }
    case(27):  //----------------------------------------------
    { /*:DeclareWorkflowVariable [name]*/
       if(Options.noisy) { cout <<"Declare Workflow Variable"<<endl; }
       workflowVar *pWV=new workflowVar();
       pWV->name=to_string(s[1]);
       pDO->AddWorkflowVariable(pWV);
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
      delete [] aX;
      delete [] aY;
      break;
    }
    case(31):  //----------------------------------------------
    { /*:OverrideStageDischargeCurve [SBID]  */

      if (Options.noisy){cout<<"Override Stage Discharge Curve"<<endl; }
      CSubBasin *pSB=pModel->GetSubBasinByID(s_to_ll(s[1]));
      if (pSB == NULL) {
        ExitGracefullyIf(pSB==NULL,"ParseManagementFile: subbasin ID in :OverrideStageDischargeCurve is invalid",BAD_DATA_WARN);
        break;
      }
      
      if (pSB->GetReservoir()==NULL){
        string advice="ParseManagementFile:The reservoir in subbasin "+to_string(pSB->GetID()) + " doesnt exist and stage discharge curve cannot be overridden.";
        ExitGracefully(advice.c_str(), BAD_DATA_WARN);
        break;
      }
      else{
        pDO->OverrideSDCurve(pSB->GetGlobalIndex());
      }

      break;
    }
    case(32):  //----------------------------------------------
    {/*:LookupTableFromCSV [name] [filename.csv]
       format:
        HeaderX, HeaderY
        x1, y1
        x2, y2
        ...
        xN, yN
        (eof)
     */
      if(Options.noisy) { cout <<"LookupTableFromCSV"<<endl; }

      const int MAX_LOOKUP_TABLE=256;

      int     N;
      double *aX = new double[MAX_LOOKUP_TABLE];
      double *aY = new double[MAX_LOOKUP_TABLE];
      int     line2(0);
      bool    eof(false);
      string  name = s[1];
      string  filename="";
      for (int i=2;i<Len;i++){filename+=to_string(s[i]); }

      ifstream    LUCSV;
      LUCSV.open(filename.c_str());
      if (LUCSV.fail()){cout << "ERROR opening file: "<<filename<<endl; return false;}

      CParser *ppCSV=new CParser(LUCSV,filename,line2);

      ppCSV->Tokenize(s,Len); //headers

      N=0;
      do {
        eof=ppCSV->Tokenize(s,Len);
        if      ((eof) || (IsComment(s[0],Len))) { }
        else {
          if(Len>=2) {
            aX[N] = s_to_d(s[0]);
            aY[N] = s_to_d(s[1]);
            N++;
          }
          else {
            delete [] aX;delete [] aY;delete ppCSV;
            WriteWarning("Incorrect line length (<2) in :LookupTableFromCSV table",Options.noisy);
            break;
          }
        }
      } while (!eof);

      LUCSV.close();
      CLookupTable *pLUT = new CLookupTable(name,aX,aY,N);
      pDO->AddUserLookupTable(pLUT);
      delete [] aX;
      delete [] aY;
      delete ppCSV;
      break;
    }

    case(40):  //----------------------------------------------
    { /*:LoopThrough [SB_GROUP or DEMAND_GROUP] [group name]  */
      if(Options.noisy) { cout <<"Start Loop"<<endl; }
      //if (Len<2){ImproperFormatWarning(":NumericalMethod",p,Options.noisy); break;}

      if ((pLoopSBGroup != NULL) || (pLoopDemandGroup !=NULL) || (loopListSize!=0)){
        ExitGracefully("ParseManagementFile: Cannot have nested :LoopThrough statements",BAD_DATA);
      }

      loopStartPos=pp->GetPosition();
      loopStartLine=pp->GetLineNumber();

      if       (!strcmp(s[1],"SB_GROUP"    ))
      {
        pLoopSBGroup=pModel->GetSubBasinGroup(s[2]);
        if (pLoopSBGroup == NULL) {
          ExitGracefully("ParseManagementFile: bad subbasin group name in :LoopThrough command",BAD_DATA_WARN);
        }
        else {
          loopCount=0;
          aLoopVector[0]=new string [pLoopSBGroup->GetNumSubbasins()];
          aLoopVector[1]=new string [pLoopSBGroup->GetNumSubbasins()];
          for (int p = 0; p < pLoopSBGroup->GetNumSubbasins(); p++)
          {
            long long SBID = pLoopSBGroup->GetSubBasin(p)->GetID();
            string SBName = pLoopSBGroup->GetSubBasin(p)->GetName();
            aLoopVector[0][p]=to_string(SBID);
            aLoopVector[1][p]=SBName;
          }

          aWildcards[0][0] = "$ID$";   aWildcards[0][1] = aLoopVector[0][0];
          aWildcards[1][0] = "$NAME$"; aWildcards[1][1] = aLoopVector[1][0];
          nWildcards=2;
        }
      }
      else if  (!strcmp(s[1],"DEMAND_GROUP"))
      {
        pLoopDemandGroup=pDO->GetDemandGroupFromName(s[2]);
        if (pLoopDemandGroup == NULL) {
          ExitGracefully("ParseManagementFile: bad demand group name in :LoopThrough command",BAD_DATA_WARN);
        }
        else {
          loopCount=0;
          aLoopVector[0]=new string [pLoopDemandGroup->GetNumDemands()];
          aLoopVector[1]=new string [pLoopDemandGroup->GetNumDemands()];
          for (int p = 0; p < pLoopDemandGroup->GetNumDemands(); p++)
          {
            long long  ID = pLoopDemandGroup->GetDemand(loopCount)->GetDemandID();
            string  DName = pLoopDemandGroup->GetDemand(loopCount)->GetName();
            aLoopVector[0][p]=to_string(ID);
            aLoopVector[1][p]=DName;
          }

          aWildcards[0][0] = "$ID$";   aWildcards[0][1] = aLoopVector[0][0];
          aWildcards[1][0] = "$NAME$"; aWildcards[1][1] = aLoopVector[1][0];
          nWildcards=2;
        }
      }
      else if (!strcmp(s[1], "LIST"))
      {
        loopListSize=s_to_i(s[2]);
        loopCount=0;
        nWildcards=0;
      }
      else {
        ExitGracefully(":LoopThrough- invalid format",BAD_DATA_WARN);
      }
      break;
    }
    case(41):  //----------------------------------------------
    { /*:EndLoopThrough   */
      if(Options.noisy) { cout <<"End Loop"<<endl; }
      bool loopDone=false;

      loopCount++;

      if ((pLoopSBGroup!=NULL) && (loopCount == pLoopSBGroup->GetNumSubbasins()))//iterating  subbasin loop
      {
        pLoopSBGroup=NULL;
        loopDone=true;
      }
      else if ((pLoopDemandGroup != NULL) && (loopCount == pLoopDemandGroup->GetNumDemands()))
      {
        pLoopDemandGroup=NULL;
        loopDone=true;
      }
      else if ((loopListSize != 0) && (loopCount == loopListSize))
      {
        loopListSize=0;
        loopDone=true;
      }

      if (loopDone) //Finished. Reset loop variables
      {
        for (int j = 0; j < nWildcards; j++) {
          delete [] aLoopVector[j];
        }
        loopCount       =0;
        loopListSize    =0;
        pLoopSBGroup    =NULL;
        pLoopDemandGroup=NULL;
      }
      else //not finished yet - update all wildcards and jump back to start of loop
      {
        for (int j = 0; j < nWildcards; j++) {
          aWildcards[j][1]=aLoopVector[j][loopCount];
        }

        // jump back to start of loop
        pp->SetPosition(loopStartPos);
        pp->SetLineCounter(loopStartLine);
      }
      break;
    }
    case(42):  //----------------------------------------------
    { /*:LoopVector $wildcard$ v1 v2 v3 ... vN   */ /*vN are strings*/
      if(Options.noisy) { cout <<"Loop Vector Data"<<endl; }
      if (loopCount==0){//only read on first loop through
        int size=Len-2;
        if ((pLoopSBGroup != NULL) && (size != pLoopSBGroup->GetNumSubbasins())) {
          ExitGracefully(":LoopVector: invalid number of terms in vector. Should be equal to size of subbasin group",BAD_DATA_WARN);
        }
        if ((pLoopDemandGroup != NULL) && (size != pLoopDemandGroup->GetNumDemands())) {
          ExitGracefully(":LoopVector: invalid number of terms in vector. Should be equal to size of demand group",BAD_DATA_WARN);
        }
        if ((loopListSize != 0) && (size != loopListSize)) {
          ExitGracefully(":LoopVector: invalid number of terms in vector. Should be equal to size of List",BAD_DATA_WARN);
        }
        if ((pLoopSBGroup == NULL) && (pLoopDemandGroup == NULL) && (loopListSize==0)) {
          ExitGracefully(":LoopVector command must be within :LoopThrough-:EndLoopThrough command block",BAD_DATA_WARN);
        }
        aLoopVector[nWildcards]=new string [size];
        for (int i=0;i<size; i++){
          aLoopVector[nWildcards][i] = s[i];
        }
        aWildcards[nWildcards][0] = s[1];
        aWildcards[nWildcards][1] = aLoopVector[nWildcards][0];
        nWildcards++;
      }
      break;
    }
    case (50): //--------------------------------------------
    { /*
      :WaterDemand [SBID] [ID] [Name]
         ...
      :EndWaterDemand
      */
      if(Options.noisy) { cout <<"Water demand object"<<endl; }
      if (Len<4){
        ExitGracefully("Incorrect number of terms in :WaterDemand command header.",BAD_DATA_WARN);
      }
      else{
        demandSBID =s_to_ll(s[1]);
        demand_ID  =s_to_i(s[2]);
        demand_name=s[3];
        CSubBasin *pSB=pModel->GetSubBasinByID(demandSBID);
        if (pSB!=NULL) {
          demand_ind=pSB->GetNumWaterDemands();
          pDemand=new CDemand(demand_ID,demand_name,demandSBID,false,pModel);
          pDemand->SetLocalIndex(demand_ind);
          pSB->AddWaterDemand(pDemand);
          pModel->GetManagementOptimizer()->AddWaterDemand(pDemand);
        }
        else {
          string warn="Invalid subbasin ID ("+to_string(demandSBID)+") in :WaterDemand command for demand "+to_string(s[3]);
           ExitGracefully(warn.c_str(), BAD_DATA_WARN);
           break;
        }
      }
      break;
    }
    case (51): //--------------------------------------------
    {/*:EndWaterDemand*/
      if(Options.noisy) { cout <<"..end water demand object"<<endl; }
      pDemand=NULL;
      break;
    }
    case (52): //--------------------------------------------
    {/*:DemandTimeSeries
       {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
       {double Qdem} x nMeasurements [m3/s]
     :EndDemandTimeSeries
     */
      if(Options.noisy) { cout <<"Demand time series"<<endl; }
      if (pDemand==NULL){
        ExitGracefully(":DemandTimeSeries must be between :WaterDemand and :EndWaterDemand commands.",BAD_DATA_WARN);
      }
      else {
        CTimeSeries *pTS = CTimeSeries::Parse(pp, false, demand_name, demandSBID, "none", Options);
        pTS->SetDemandID(demand_ID);
        pDemand->SetDemandTimeSeries(pTS);
      }
      break;
    }
    case (53): //--------------------------------------------
    {/*:ReturnTimeSeries
       {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
       {double Qret} x nMeasurements [m3/s]
     :EndReturnTimeSeries
     */
      if(Options.noisy) { cout <<"Return time series"<<endl; }
      if (pDemand==NULL){
        ExitGracefully(":ReturnTimeSeries must be between :WaterDemand and :EndWaterDemand commands.",BAD_DATA_WARN);
      }
      else {
        CTimeSeries *pTS = CTimeSeries::Parse(pp, false, demand_name, demandSBID, "none", Options);
        pTS->SetDemandID(demand_ID);
        pDemand->SetReturnTimeSeries(pTS);
      }
      break;
    }
    case (54): //--------------------------------------------
    {/*:ReturnDestination [SBID] */
      if(Options.noisy) { cout <<"Return destination ID"<<endl; }
      if (pDemand==NULL){
        ExitGracefully(":ReturnDestination must be between :WaterDemand and :EndWaterDemand commands.",BAD_DATA_WARN);
      }
      else {
        if (pModel->GetSubBasinByID(s_to_ll(s[1])) == NULL) {
          ExitGracefully("ParseManagementFile: Bad subbasin ID supplied to :ReturnDestination command.",BAD_DATA_WARN);
        }
        pDemand->SetTargetSBID(s_to_ll(s[1]));
      }
      break;
    }
    case (55): //--------------------------------------------
    {/*:IrrigationDestination [HRUGroup]*/
      if(Options.noisy) { cout <<"Irrigation destination"<<endl; }
      if (pDemand==NULL){
        ExitGracefully(":IrrigationDestination must be between :WaterDemand and :EndWaterDemand commands.",BAD_DATA_WARN);
      }
      else {
        CHRUGroup *pGrp=pModel->GetHRUGroup(s[1]);
        if (pGrp == NULL) {
          ExitGracefully(":IrrigationDestination uses invalid HRU Group.",BAD_DATA_WARN);
        }
        else {
          pDemand->SetIrrigationGroup(pGrp->GetGlobalIndex());
        }
      }
      break;
    }
    case (56): //--------------------------------------------
    {/*:ResetDate [julian_day]*/
      if(Options.noisy) { cout <<"Demand reset date"<<endl; }
      if (pDemand==NULL){
        ExitGracefully(":ResetDate must be between :WaterDemand and :EndWaterDemand commands.",BAD_DATA_WARN);
      }
      else {

        pDO->SetCumulativeDate(s_to_i(s[1]), to_string(demand_ID));
      }
      break;
    }
    case (57): //--------------------------------------------
    {/*:Penalty [value]*/
      if(Options.noisy) { cout <<"Demand penalty"<<endl; }
      if (pDemand==NULL){
        ExitGracefully(":Penalty must be between :WaterDemand and :EndWaterDemand commands.",BAD_DATA_WARN);
      }
      else {
        pDO->SetDemandPenalty(to_string(demand_ID),s_to_d(s[1]));
      }
      break;
    }
    case (58): //--------------------------------------------
    {/*:ReturnFraction [value]*/
      if(Options.noisy) { cout <<"Return flow fraction"<<endl; }
      if (pDemand==NULL){
        ExitGracefully(":ReturnFraction must be between :WaterDemand and :EndWaterDemand commands.",BAD_DATA_WARN);
      }
      else {
        pDemand->SetReturnFraction(s_to_d(s[1]));
      }
      break;
    }
    case (59): //--------------------------------------------
    {/*:ReturnExpression [expression]*/
      if(Options.noisy) { cout <<"Return flow expression"<<endl; }
      if (pDemand==NULL){
        ExitGracefully(":ReturnExpression must be between :WaterDemand and :EndWaterDemand commands.",BAD_DATA_WARN);
      }
      else {
        expressionStruct *pExp;

        pExp=pDO->ParseExpression((const char**)(s),Len,pp->GetLineNumber(),pp->GetFilename());

        if (pDO->GetDebugLevel()>=1){
          SummarizeExpression((const char**)(s),Len,pExp);
        }

        if (pExp!=NULL){
          pDemand->SetReturnExpression(pExp);
        }
        else {
          string warn ="Invalid expression in :ReturnExpression command at line " + pp->GetLineNumber();
          WriteWarning(warn.c_str(),Options.noisy);
        }
      }
      break;
    }
    case (60): //--------------------------------------------
    {/*:DemandFlowFraction [value]*/
      if(Options.noisy) { cout <<"Demand flow fraction"<<endl; }
      if (pDemand==NULL){
        ExitGracefully(":DemandFlowFraction must be between :WaterDemand and :EndWaterDemand commands.",BAD_DATA_WARN);
      }
      else {
        pDemand->SetDemandFraction(s_to_d(s[1]));
      }
      break;
    }
    case (61): //--------------------------------------------
    {/*:DemandLookupTable
       N
       {Q_i D_i} x N
       :EndDemandLookupTable
     */
      if(Options.noisy) { cout <<"Demand flow fraction"<<endl; }
      ExitGracefully(":DemandLookupTable.",STUB);
      break;
    }
    case (62): //--------------------------------------------
    {/*:DemandExpression [expression]*/
      if(Options.noisy) { cout <<"Demand expression"<<endl; }
      if (pDemand==NULL){
        ExitGracefully(":DemandExpression must be between :WaterDemand and :EndWaterDemand commands.",BAD_DATA_WARN);
      }
      else {
        expressionStruct *pExp;

        pExp=pDO->ParseExpression((const char**)(s),Len,pp->GetLineNumber(),pp->GetFilename());

        if (pDO->GetDebugLevel()>=1){
          SummarizeExpression((const char**)(s),Len,pExp);
        }

        if (pExp!=NULL){
          pDemand->SetDemandExpression(pExp);
        }
        else {
          string warn ="Invalid expression in :DemandExpression command at line " + pp->GetLineNumber();
          WriteWarning(warn.c_str(),Options.noisy);
        }
      }
      break;
    }
    case (63): //--------------------------------------------
    { /*
      :ReservoirWaterDemand [SBID] [ID] [Name]
         ...
      :EndReservoirWaterDemand
      */
      if(Options.noisy) { cout <<"Reservoir Water demand object"<<endl; }
      if (Len<4){
        ExitGracefully("Incorrect number of terms in :ReservoirWaterDemand command header.",BAD_DATA_WARN);
      }
      else{
        demandSBID =s_to_ll(s[1]);
        demand_ID  =s_to_ll(s[2]);
        demand_name=s[3];
        CSubBasin *pSB=pModel->GetSubBasinByID(demandSBID);
        if (pSB!=NULL) {
          if (pSB->GetReservoir() != NULL) {
            demand_ind=pSB->GetReservoir()->GetNumWaterDemands();
            pDemand=new CDemand(demand_ID,demand_name,demandSBID,true,pModel);
            pDemand->SetLocalIndex(demand_ind);
            pSB->GetReservoir()->AddDemand(pDemand);
            pModel->GetManagementOptimizer()->AddWaterDemand(pDemand);
          }
          else
          {
            ExitGracefully("There is no reservoir in the subbasin specified within the :ReservoirWaterDemand command header.",BAD_DATA_WARN);
            break;
          }
        }
        else {
           ExitGracefully("Invalid subbasin ID in :ReservoirWaterDemand command header.",BAD_DATA_WARN);
           break;
        }
      }
      break;
    }
    case (64): //--------------------------------------------
    {/*:EndReservoirWaterDemand*/
      if(Options.noisy) { cout <<"..end reservoir water demand object"<<endl; }
      pDemand=NULL;
      break;
    }
    case (65): //--------------------------------------------
    {/*:IsUnrestricted*/
      if(Options.noisy) { cout <<"Set demand as Unrestricted"<<endl; }
      if (pDemand == NULL) {
        ExitGracefully(":IsUnrestricted command must be between :WaterDemand and :EndWaterDemand commands.",BAD_DATA_WARN);
      }
      else{
        pDO->SetDemandAsUnrestricted(pDemand->GetName());
      }
      break;
    }
    case (70): //---------------------------------------------
    {/*:UserTimeSeries {name} [units]
        {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
        {double val} x nValues
      :EndUserTimeSeries
     */
      CTimeSeries *pTimeSer;
      if(Options.noisy) { cout <<"User-specified Time Series"<<endl; }
      pTimeSer=CTimeSeries::Parse(pp,true,s[1], DOESNT_EXIST, "none", Options);
      pModel->GetManagementOptimizer()->AddUserTimeSeries(pTimeSer);
      break;
    }
    case (80): //---------------------------------------------
    {/*:NonlinearVariable [?name] [corresponding DV]
     e.g.,
     :NonlinearVariable ?Q130 Q130
     */
      if(Options.noisy) { cout <<"Non-linear variable declaration"<<endl; }
      pModel->GetManagementOptimizer()->AddNonLinVariable(s[1],s[2]);
      break;
    }
    case (81): //---------------------------------------------
    {/*:MaxNLSolverIterations [#]*/
      if(Options.noisy) { cout <<"Maximum solver iterations"<<endl; }
      pModel->GetManagementOptimizer()->SetMaxIterations(s_to_i(s[1]));
      break;
    }
    case (82): //---------------------------------------------
    {/*:NLSolverTolerance [#]*/
      if(Options.noisy) { cout <<"Non-linear solver tolerance"<<endl; }
      pModel->GetManagementOptimizer()->SetSolverTolerance(s_to_d(s[1]));
      break;
    }
    case (83): //---------------------------------------------
    {/*:NLRelaxationCoeff [#]*/
      if(Options.noisy) { cout <<"Non-linear solver relaxation coefficient"<<endl; }
      pModel->GetManagementOptimizer()->SetRelaxationCoeff(s_to_d(s[1]));
      break;
    }
    default://------------------------------------------------
    {
      /*cout << "UNREC LINE: (Len=" << Len << ")";
      for (int i = 0; i < Len; i++) {
        cout<<s[i]<<"|";
      }
      cout<<endl;*/
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
            string warn ="IGNORING unrecognized command: " + string(s[0])+ " in .rvm file";
            WriteWarning(warn,Options.noisy);
          }
          break;
        }
        default:
        {
          string errString = "Unrecognized command in .rvm file:\n   " + string(s[0]);
          ExitGracefully(errString.c_str(),BAD_DATA_WARN);//STRICT
          break;
        }
      }
    }
    }//end switch(code)

    firstword=pp->Peek();
    if ((firstword == ":DefineDecisionVariable") || (firstword == ":DemandExpression") || (firstword == ":ReturnExpression") || (firstword == ":DefineWorkflowVariable"))
    {
      pp->NextIsMathExp();
    }

    end_of_file=pp->Tokenize(s,Len);
    swapWildcards((const char**)(s),Len,aWildcards,nWildcards);

     //return after file redirect, if in tertiary file
    if ((end_of_file) && (pSecondaryParser != NULL))
    {
      INPUT3.clear();
      INPUT3.close();
      delete pp;
      pp=pSecondaryParser;
      pSecondaryParser=NULL;
      end_of_file=pp->Tokenize(s,Len);
      swapWildcards((const char**)(s),Len,aWildcards,nWildcards);
    }
    else{
      //return after file redirect, if in secondary file
      if((end_of_file) && (pMainParser!=NULL))
      {
        INPUT2.clear();
        INPUT2.close();
        delete pp;
        pp=pMainParser;
        pMainParser=NULL;
        end_of_file=pp->Tokenize(s,Len);
        swapWildcards((const char**)(s),Len,aWildcards,nWildcards);
      }
    }
  } //end while !end_of_file
  RVM.close();

  pDO->InitializeDemands    (pModel,Options);
  pDO->InitializePostRVMRead(pModel,Options);

  if (loopCount != 0) {
    ExitGracefully("Unclosed loop statement in .rvm file.", BAD_DATA_WARN);
  }

  delete [] aLoopVector;
  for (int i = 0; i < MAX_WILDCARDS; i++) {
    delete [] aWildcards[i];
  }
  delete [] aWildcards;

  delete pp;
  pp=NULL;

  return true;
}
