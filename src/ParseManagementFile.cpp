/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2024 the Raven Development Team
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

  ifstream    RVM;
  RVM.open(Options.rvm_filename.c_str());
  if(RVM.fail()) {
    cout << "ERROR opening model management file: "<<Options.rvm_filename<<endl; return false;
  }

  int   Len,line(0),code;
  char *s[MAXINPUTITEMS];
  CParser *pp=new CParser(RVM,Options.rvm_filename,line);

  ifstream INPUT2;           //For Secondary input
  CParser *pMainParser=NULL; //for storage of main parser while reading secondary files

  const int MAX_WILDCARDS=10;

  int loopCount=0;
  streampos loopStartPos;
  int       loopStartLine;
  CSubbasinGroup *pLoopSBGroup=NULL;
  CDemandGroup   *pLoopDemandGroup=NULL;
  string **aWildcards  = new string *[MAX_WILDCARDS  ];
  string **aLoopVector = new string *[MAX_WILDCARDS  ];
  for (int i = 0; i < MAX_WILDCARDS; i++) {
    aWildcards [i]=new string [2];
    aWildcards [i][0] = "";
    aWildcards [i][1] = "";
    aLoopVector[i]=NULL;
  }
  int    nWildcards=0;

  CDemandOptimizer *pDO=pModel->GetDemandOptimizer();

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

    //-------------------MODEL ENSEMBLE PARAMETERS----------------
    else if(!strcmp(s[0],":LookbackDuration"))            { code=1;  }
    else if(!strcmp(s[0],":DebugLevel"))                  { code=2; }
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
    else if(!strcmp(s[0],":DefineControlVariable"))       { code=25; }
    else if(!strcmp(s[0],":LookupTable"))                 { code=30; }

    else if(!strcmp(s[0],":LoopThrough"))                 { code=40; }
    else if(!strcmp(s[0],":EndLoopThrough"))              { code=41; }
    else if(!strcmp(s[0],":LoopVector"))                  { code=42; }

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

      INPUT2.open(filename.c_str());
      if(INPUT2.fail()) {
        string warn;
        warn=":RedirectToFile (from .rvm): Cannot find file "+filename;
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
      manConstraint    *pConst=NULL;
      decision_var     *pDV = new decision_var(s[1],DOESNT_EXIST,DV_USER,pDO->GetNumUserDVs());
      
      pDO->AddDecisionVar(pDV);
      pExp=pDO->ParseExpression((const char**)(s),Len,pp->GetLineNumber(),pp->GetFilename());

      if (pDO->GetDebugLevel()>=1){
        SummarizeExpression((const char**)(s),Len,pExp); 
      }

      if (pExp!=NULL){
        pConst = pDO->AddGoalOrConstraint(s[1], false);
        pConst->AddExpression(pExp); 
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
         :EndOperatingRegime
         :Penalty [value] {value2} 
       :EndManagementGoal
       :ManagementGoal [name] 
         :Expression [expression]
         :Penalty [value] {value2} 
       :EndManagementGoal
     */
      if(Options.noisy) { cout <<"Management Constraint or Management Goal"<<endl; }

      manConstraint *pConst=pDO->AddGoalOrConstraint(s[1], is_goal);

      expressionStruct *pExp;      
      bool is_first=true;
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
          op_regime *pOR=new op_regime(s[1]);//[OPUPDATE] 
          pConst->AddOperatingRegime(pOR,is_first);
          is_first=false;
        }
        //----------------------------------------------
        else if (!strcmp(s[0], ":EndOperatingRegime")) 
        {
          //does nothing
        }
        //----------------------------------------------
        else if(!strcmp(s[0], ":Expression")) 
        {
          if (pConst->GetCurrentExpression() != NULL) {
            ExitGracefully("ParseManagementFile: only one :Expression allowed in each :OperatingRegime command block (or only one if no :OperatingRegime blocks used).",BAD_DATA_WARN);
            break;
          }
          pExp=pDO->ParseExpression((const char**)(s),Len,pp->GetLineNumber(),pp->GetFilename());
          if (pExp!=NULL){
            pConst->AddExpression(pExp);
          }
          else {
            string warn ="Invalid expression in :Expression command at line " + pp->GetLineNumber();
            WriteWarning(warn.c_str(),Options.noisy);
          }
          if (pDO->GetDebugLevel()>=1){
            SummarizeExpression((const char**)(s),Len,pExp); 
          }
        }
        //----------------------------------------------
        else if (!strcmp(s[0], ":Condition")) 
        {
          //TODO: Would it be better to support @date(), @between, @day_of_year() in general expression??
          //:Condition !Q32[0] < 300 + @ts(myTs,0)
          //:Condition DATE IS_BETWEEN 1975-01-02 and 2010-01-02
          //:Condition DATE > @date(1975-01-02) //[NOT YET SUPPORTED] 
          //:Condition DATE < @date(2010-01-02) //[NOT YET SUPPORTED] 
          //:Condition MONTH = 2
          //:Condition DAY_OF_YEAR IS_BETWEEN 173 and 210
          //:Condition DAY_OF_YEAR > 174
          //:Condition DAY_OF_YEAR < 210
          //:Condition DAY_OF_YEAR IS_BETWEEN 300 20 //wraps around 
          //:Condition @is_between(DAY_OF_YEAR,300,20) = 1  //[NOT YET SUPPORTED] 
          if (pConst!=NULL){
            bool badcond=false;
            exp_condition *pCond = new exp_condition();
            pCond->dv_name=s[1];
            bool is_exp=false;
            for (int i = 0; i < Len; i++) {
              if ((s[i][0]=='+') || (s[i][0]=='-') || (s[i][0]=='*') || (s[i][0]=='/') || (s[i][0]=='=') || (s[i][0]=='<') || (s[i][0]=='>')){
                is_exp=true;
              }
            }
            if (is_exp) {
              pCond->pExp=pDO->ParseExpression((const char**)(s),Len,pp->GetLineNumber(),pp->GetFilename());
              pConst->AddOpCondition(pCond);
            }
            else{
              if      (!strcmp(s[2],"IS_BETWEEN"     )){pCond->compare=COMPARE_BETWEEN;}
              else if (!strcmp(s[2],"IS_GREATER_THAN")){pCond->compare=COMPARE_GREATERTHAN;}
              else if (!strcmp(s[2],"IS_LESS_THAN"   )){pCond->compare=COMPARE_LESSTHAN;}
              else if (!strcmp(s[2],"IS_EQUAL_TO"    )){pCond->compare=COMPARE_IS_EQUAL;}
              else if (!strcmp(s[2],"IS_NOT_EQUAL_TO")){pCond->compare=COMPARE_NOT_EQUAL;} 
              else {
                ExitGracefully("ParseManagementFile: unrecognized comparison operator in :Condition statement",BAD_DATA_WARN);
                break;
              }
              pCond->value=s_to_d(s[3]);
              if (Len>=5){
                pCond->value2 = s_to_d(s[4]);
              }
              if      (!strcmp(s[1],"DATE"     )){
                pCond->date_string=s[3];
                if (Len>=5){
                  pCond->date_string2 = s[4];
                }
              }

              if (pCond->dv_name[0] == '!') { //decision variable 
                char   tmp =pCond->dv_name[1];
                string tmp2=pCond->dv_name.substr(2);
                char code=pCond->dv_name[1];
                if ((code=='Q') || (code=='h') || (code=='I')){
                  long SBID=s_to_l(tmp2.c_str());
                  if (pModel->GetSubBasinByID(SBID) == NULL) {
                    ExitGracefully("ParseManagementFile: Subbasin ID in :Condition statement is invalid.",BAD_DATA_WARN);
                  }
                  else if (!pModel->GetSubBasinByID(SBID)->IsEnabled()) {
                    WriteWarning("ParseManagementFile: Subbasin in :Condition statement is disabled in this model configuration. Conditional will be assumed true.",Options.noisy);
                    badcond=true;
                  }
                  else if ((code == 'h') || (code == 'I')) {
                    if (pModel->GetSubBasinByID(SBID)->GetReservoir() == NULL) {
                      ExitGracefully("ParseManagementFile: !h or !I used in :Condition statement for subbasin without lake or reservoir",BAD_DATA_WARN);
                    }
                  }
                }
                else { //demand 
                  int d=pDO->GetDemandIndexFromName(tmp2);
                  if (d == DOESNT_EXIST) {
                    ExitGracefully("ParseManagementFile: !D or !C used in :Condition statement has invalid demand ID",BAD_DATA_WARN);
                  }
                }
              }
              if (!badcond)
              {
                pCond->p_index=pDO->GetIndexFromDVString(pCond->dv_name); 
                pConst->AddOpCondition(pCond);
              }
            }
          } 
          else{
            ExitGracefully("ParseManagementFile: :Condition statement must appear after valid :Expression in :ManagementConstraint command",BAD_DATA_WARN);
          }
        } 
        //----------------------------------------------
        else if (!strcmp(s[0], ":Penalty")) {
          if (!pConst->is_goal) {
            ExitGracefully("ParseManagementFile: :Penalty found within :ManagementConstraint command. It will be ignored, as all constraints have infinite penalty.",BAD_DATA_WARN);
          }
          pConst->penalty_under = s_to_d(s[1]);
          pConst->penalty_over  = s_to_d(s[1]);
          if (Len >= 3) {
            pConst->penalty_over = s_to_d(s[2]);
          }
        }
        //----------------------------------------------
        else if (!strcmp(s[0], ":Priority")) {
          pConst->priority = s_to_i(s[1]); //for later 
        }
        //----------------------------------------------
        else if (!strcmp(s[0], ":EndManagementGoal")) {
          break;
        }
        //----------------------------------------------
        else if (!strcmp(s[0], ":EndManagementConstraint")) {
          break;
        } 
        else {
          WriteWarning("ParseManagementFile: Unrecognized command in :ManagementConstraint command block",Options.noisy);
        }
        firstword=pp->Peek();
        if (firstword == ":Expression") {pp->NextIsMathExp();}
        if (firstword == ":Condition")  {pp->NextIsMathExp();}
      }
      break;
    }
    case(24):  //----------------------------------------------
    { /*:DeclareDecisionVariable [name]  */
      if(Options.noisy) { cout <<"Declare Decision Variable"<<endl; }      
      decision_var     *pDV = new decision_var(s[1],DOESNT_EXIST,DV_USER,pDO->GetNumUserDVs());
      pDO->AddDecisionVar(pDV);
      break;
    }
    case(25):  //----------------------------------------------
    { /*:DefineControlVariable [name] = [expressionRHS] */
      if(Options.noisy) { cout <<"Define Control Variable"<<endl; }
      expressionStruct *pExp;      
      pExp=pDO->ParseExpression((const char**)(s),Len,pp->GetLineNumber(),pp->GetFilename());

      if (pDO->GetDebugLevel()>=1){
        SummarizeExpression((const char**)(s),Len,pExp); 
      }

      if (pExp!=NULL){
        pDO->AddControlVariable(s[1],pExp);
      } 
      else {
        string warn ="Invalid expression in :DefineControlVariable command at line " + pp->GetLineNumber();
        WriteWarning(warn.c_str(),Options.noisy);
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
      delete [] aX;
      delete [] aY;
      break;
    }
    case(40):  //----------------------------------------------
    { /*:LoopThrough [SB_GROUP or DEMAND_GROUP] [group name]  */
      if(Options.noisy) { cout <<"Start Loop"<<endl; }      
      //if (Len<2){ImproperFormatWarning(":NumericalMethod",p,Options.noisy); break;}

      loopStartPos=pp->GetPosition();
      loopStartLine=pp->GetLineNumber();

      if       (!strcmp(s[1],"SB_GROUP"    )){
       
        if ((pLoopSBGroup != NULL) || (pLoopDemandGroup !=NULL)){
          ExitGracefully("ParseManagementFile: Cannot have nested :LoopThrough statements",BAD_DATA);
        }
        pLoopSBGroup=pModel->GetSubBasinGroup(s[2]);
        if (pLoopSBGroup == NULL) {
          ExitGracefully("ParseManagementFile: bad subbasin group name in :LoopThrough command",BAD_DATA_WARN);
        } 
        else {
          loopCount=0;
          long     SBID = pLoopSBGroup->GetSubBasin(loopCount)->GetID();
          string SBName = pLoopSBGroup->GetSubBasin(loopCount)->GetName();
          
          aWildcards[0][0] = "$ID$";   aWildcards[0][1] = to_string(SBID);
          aWildcards[1][0] = "$NAME$"; aWildcards[1][1] = SBName;

          nWildcards=2;
          //look for $SBID$ 
        }
      }
      else if  (!strcmp(s[1],"DEMAND_GROUP")){
        //\todo[funct]
        ExitGracefully("Loop through demand group not yet supported",STUB);
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

      if (pLoopSBGroup!=NULL)//iterating  subbasin loop 
      {
        if (loopCount == pLoopSBGroup->GetNumSubbasins())
        {
          pLoopSBGroup=NULL;
          loopDone=true;
        }
        else
        {
          long     SBID = pLoopSBGroup->GetSubBasin(loopCount)->GetID();
          string SBName = pLoopSBGroup->GetSubBasin(loopCount)->GetName();
          
          aWildcards[0][0] = "$ID$";   aWildcards[0][1] = to_string(SBID);
          aWildcards[1][0] = "$NAME$"; aWildcards[1][1] = SBName;
        }
      }
      else if (pLoopDemandGroup != NULL) {
        
        if (loopCount == pLoopDemandGroup->GetNumDemands()){
          pLoopDemandGroup=NULL;
          loopDone=true;
        }
        else {

        }
      }
      else {
        //shouldnt happen 
        loopDone=true;
      }

      if (loopDone) //Finished. Reset loop variables
      {
        for (int j = 0; j < nWildcards - 2; j++) {
          delete [] aLoopVector[j];
        }
        loopCount       =0;
        pLoopSBGroup    =NULL;
        pLoopDemandGroup=NULL;
      }
      else {
        if (nWildcards>=2){
          for (int j = 0; j < nWildcards - 2; j++) {
            aWildcards[j+2][1]=aLoopVector[j][loopCount];
          }
        }

        cout<<" LOOP COUNT: "<<loopCount<<" of "<< endl;
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
        if ((pLoopSBGroup == NULL) && (pLoopDemandGroup == NULL)) {
          ExitGracefully(":LoopVector command must be within :LoopThrough-:EndLoopThrough command block",BAD_DATA_WARN);
        }
        if (nWildcards>=2){
          aLoopVector[nWildcards-2]=new string [size];
          for (int i=0;i<size; i++){
            aLoopVector[nWildcards-2][i] = s[i + 2];
          }
          aWildcards[nWildcards][0] = s[1];
          aWildcards[nWildcards][1] = aLoopVector[nWildcards-2][0];
          nWildcards++;
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
    if (firstword == ":DefineDecisionVariable") {
      pp->NextIsMathExp();
    }

    end_of_file=pp->Tokenize(s,Len);
    swapWildcards((const char**)(s),Len,aWildcards,nWildcards);

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
  } //end while !end_of_file
  RVM.close();

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