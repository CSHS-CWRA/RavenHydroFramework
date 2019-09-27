/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright © 2008-2014 the Raven Development Team
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Properties.h"
#include "Model.h"
#include "ParseLib.h"
#include "GlobalParams.h"
#include "GroundwaterClass.h"
#include <string>

bool ParseGWPropArray(CParser *p,int *indices,double **properties,
                    int &num_read,string *tags,const int line_length,const int max_classes);
void  RVGParameterWarning   (string  *aP, class_type *aPC, int &nP, const optStruct &Options);
void  AddToMasterGWPropsParamList   (string        *&aPm, class_type       *&aPCm, int       &nPm, 
                                     const string  *aP , const class_type *aPC , const int &nP); 

//////////////////////////////////////////////////////////////////
/// \brief This method parses the groundwater stress period class properties .rvg file
///
/// \details Parses input \b model.rvg file that defines ET and drain properties for GW \n
///   ":GWSPProperties" string name \n
///   ":StressPeriod"  string or int "stress period" \n
///   ":ETProperties"  \n
///      int id, int numcell,double pet_max,double ext_depth  \n
///   :EndETProperties  \n
///   ":DrainProperties"  \n
///      int id, int numcell, double drain_elev, double conductance \n
///   :EndDrainProperties  \n
///   ":Recharge"  \n
///      int id, int numcell, double recharge rate \n
///   :EndRecharge  \n
///
/// \param *&pModel [in] The input model object
/// \param &Options [in] Global model options information
/// \return Boolean variable indicating success of parsing
//
bool ParseGWPropsFile(CModel         *&pModel,
                      const optStruct &Options)
{

  global_struct         global_template;
  global_struct         parsed_globals;

  int                   num_parsed_gw=0;
  CGWStressPeriodClass *pGWClasses  [MAX_SP_CLASSES];
  gw_struct             parsed_gw   [MAX_SP_CLASSES];
  string                gwtags      [MAX_SP_CLASSES];

  const int             MAX_PROPERTIES_PER_LINE=20;

  int                   MAX_NUM_IN_CLASS=MAX_GW_CELLS;
  int                  *indices;
  double              **properties;
  string                aParamStrings[MAXINPUTITEMS];
  int                   nParamStrings=0;

  string               *gwName;

  CGWStressPeriodClass *pGWSP=pModel->GetGroundwaterModel()->GetGWSPs(0);

  CGlobalParams::InitializeGlobalParameters(global_template,true);
  CGlobalParams::InitializeGlobalParameters(parsed_globals,false); 

  CGWStressPeriodClass::InitializeGWProperties(parsed_gw[0],true);//zero-index gw is template - does not create class
  gwtags    [0]="[DEFAULT]";
  num_parsed_gw++;

  indices   =new int     [MAX_NUM_IN_CLASS];
  properties=new double *[MAX_NUM_IN_CLASS];
  for (int i=0;i<MAX_NUM_IN_CLASS;i++){
    properties[i]=new double [MAX_PROPERTIES_PER_LINE];
  }

  const int MAX_PARAMS=100;
  string     *aP  =new string [MAX_PARAMS];
  class_type *aPC =new class_type [MAX_PARAMS];
  
  int         nPmaster(0);   //number of parameters in master list of required parameters
  string     *aPmaster=NULL;  //array of parameter names in master list of required parameters 
  class_type *aPCmaster=NULL; //array of parameter class types in master list of required parameters 
   
  cout<<"Generating Master Parameter List..."<<endl;

  // Check parameters for each estimation/correction algorithm:
  //--------------------------------------------------------------------------
 // pModel->GetParticipatingParamList(aP,aPC,nP,Options);
 // AddToMasterGWPropsParamList(aPmaster,aPCmaster,nPmaster,aP,aPC,nP);

  // Check parameters for each hydrological process algorithm:
  //--------------------------------------------------------------------------
 // for (int j=0;j<pModel->GetNumProcesses();j++)
 // {
  //  pModel->GetProcess(j)->GetParticipatingParamList(aP,aPC,nP);
 //   AddToMasterGWPropsParamList(aPmaster,aPCmaster,nPmaster,aP,aPC,nP);
 // }
  
  //Prior to this, all parameters default to NOT_NEEDED or NOT_NEEDED_AUTO
  //this overrides default values of required variables to NOT_SPECIFIED and AUTO_COMPUTE
 double val;
 for (int p=0;p<nPmaster;p++)
 {
 //   val=CGWStressPeriodClass::GetGWProperty(parsed_gw[0],aPmaster[p]);
   val=pGWSP->GetGWProperty(parsed_gw[num_parsed_gw],aPmaster[p],p);
 //   if (val==NOT_NEEDED_AUTO){val=AUTO_COMPUTE;}
 //   else if (val==NOT_NEEDED){val=NOT_SPECIFIED;}
 //   CGWStressPeriodClass::SetGWProperty    (parsed_gw[0],aPmaster[p],val,0);
   pGWSP->SetGWProperty(parsed_gw[num_parsed_gw],aPmaster[p],val,p);
 }

  cout <<"======================================================"<<endl;
  cout << "Parsing Groundwater Properties File " << Options.rvg_filename <<"..."<<endl;
  cout <<"======================================================"<<endl;

  int   code;
  bool  done;
  bool  ended(false);
  int   Len,line(0);
  char *s[MAXINPUTITEMS];

  ifstream INPUT;
  INPUT.open(Options.rvg_filename.c_str());  
  if (INPUT.fail()){ 
    cout << "ERROR opening file: "<< Options.rvg_filename<<endl; return false;}
  
  ExitGracefullyIf(pModel==NULL,
    "ParseGWClassPropertiesFile: model has not yet been created.",BAD_DATA);

  CParser *p=new CParser(INPUT,line);

  ifstream INPUT2;           //For Secondary input
  CParser *pMainParser=NULL; //for storage of main parser while reading secondary files

  //--Sift through file-----------------------------------------------
  bool end_of_file=p->Tokenize(s,Len);
  while (!end_of_file) 
  {
    if (ended){break;}
    if (Options.noisy){ cout << "reading line " << p->GetLineNumber() << ": ";}
    /*assign code for switch statement
    ------------------------------------------------------------------
    <100         : ignored/special
    0   thru 100 : Options
    100 thru 200 : Groundwater Classes
    ------------------------------------------------------------------
    */

    code=0;   
    //---------------------SPECIAL -----------------------------
    if       (Len==0)                                  {code=-1; }
    else if  (!strcmp(s[0],"*"                       )){code=-2; }//comment
    else if  (!strcmp(s[0],"%"                       )){code=-2; }//comment
    else if  (!strcmp(s[0],"#"                       )){code=-2; }//comment
    else if  (s[0][0]=='#')                            {code=-2; }//comment
    else if  (!strcmp(s[0],":RedirectToFile"         )){code=-3; }//redirect to secondary file
    else if  (!strcmp(s[0],":End"                    )){code=-4; }//stop reading
    //--------------------DECLARATION PARAMS --------------------------
    else if  (!strcmp(s[0],":GWSPProperties"         )){code=1;  }//REQUIRED
    else if  (!strcmp(s[0],":EndGWSPProperties"      )){code=2;  }//REQUIRED
    else if  (!strcmp(s[0],":StressPeriod"           )){code=3;  }//REQUIRED
    else if  (!strcmp(s[0],":EndStressPeriod"        )){code=4;  }//REQUIRED
    //--------------------GROUNDWATER PARAMS -------------------------
    else if  (!strcmp(s[0],":ETProperties"           )){code=100;}
    else if  (!strcmp(s[0],":DrainProperties"        )){code=101;}
    else if  (!strcmp(s[0],":Recharge"               )){code=102;}
    else if  (!strcmp(s[0],":HydraulicConductivity"  )){code=103;}
    else if  (!strcmp(s[0],":Storage"                )){code=104;}

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
        INPUT2.open(filename.c_str()); 
        if (INPUT2.fail()){
          ostrstream FILENAME;
          FILENAME<<":RedirectToFile: Cannot find file "<<filename<<ends;
          ExitGracefully(FILENAME.str() ,BAD_DATA); 
        }
        else{
          pMainParser=p;    //save pointer to primary parser
          p=new CParser(INPUT2,line);//open new parser
        }
        break;
      }
      case(-4):  //----------------------------------------------
      {/*:End*/
        if (Options.noisy) {cout <<"EOF"<<endl;} ended=true; break;
      }
      case(1):  //----------------------------------------------
      {/*:GWSPPRoperties [name]*/
        if (Options.noisy) {cout <<"Reading Groundwater Stress Periods..."<<endl;}
        char tmp[64];
        sprintf(tmp, "GroundwaterStressPeriod%i", CGWStressPeriodClass::GetNumGWSPClasses()+1);
        string gName = tmp;
        if (Len > 1){                 // A gw name may have been specified
          if ((char)*(s[1]) != '#'){  // first character a '#' comment character?
            gName   = s[1];             // Valid gw name, so we use it
            gwName  = &gName;
          }
        }
        if ((num_parsed_gw>=(MAX_SP_CLASSES-1)) || (num_parsed_gw<0)){
              ExitGracefully("ParseGWPropsFile: exceeded maximum # of SP classes",BAD_DATA);}

        //CGWStressPeriodClass::InitializeGWProperties(parsed_gw[num_parsed_gw],false);
        pGWSP->InitializeGWProperties(parsed_gw[num_parsed_gw],false);
        num_parsed_gw++;
        pGWClasses[num_parsed_gw-1]=new CGWStressPeriodClass(gName);
        pGWSP = pGWClasses[num_parsed_gw-1];
        gwtags    [num_parsed_gw-1]=s[1];

        //cout<<"xxx: "<<pGWClasses[num_parsed_gw-1]<<endl;
        pModel->GetGroundwaterModel()->AddGWSP(pGWClasses[num_parsed_gw-1]);
        cout << "    " << "GW PROPERTY GROUP " << gName << " CREATED" << endl; 
        break;
      }
      case(2):  //----------------------------------------------
      {/*":EndGWSPProperties"*/
        if (Options.noisy) {cout <<"End Groundwater Properties"<<endl;}
        if(num_parsed_gw==0)
        {
          ExitGracefully("ParseGWPropsFile:: :EndGWSPProperties encountered before :GWSPPRoperties",BAD_DATA);
        }
        ended=true;
        break;
      }
      case(3):  //----------------------------------------------
      {/*:StressPeriod*/
        if (Options.noisy) {cout <<"Reading Stress Period data..."<<endl;}
        if (Len!=2){p->ImproperFormat(s); break;} 
        //p->Tokenize(s,Len);
        string spe = "stressperiod";
        pGWSP->SetGWProperty(spe,s_to_d(s[1]),0);
        
        //pGWClasses[num_parsed_gw-1]->SetGWProperty(spe,s_to_d(s[1]),0);
        //parsed_gw [num_parsed_gw].stressperiod =s_to_d(s[1]);
        //ExitGracefullyIf(pGWClasses==NULL,
        //  "ParseGWPropsFile::Stress Period specified outside of a :Groundwater-:EndGroundwater statement",BAD_DATA);
        //string paramName = "STRESS_PERIOD";
        //pGWClasses[num_parsed_gw]->SetGWProperty(paramName,s_to_d(s[1]),0);
        cout << "    " << "STRESS PERIOD " << parsed_gw [num_parsed_gw].stressperiod << " CREATED" << endl;
        break;
      }
      case(4):  //----------------------------------------------
      {/*:EndStressPeriod*/
        if (Options.noisy) {cout <<"End Stress Period"<<endl;}
        break;
      }
      //===========================================================================================
      //===========================================================================================
      case(100):  //----------------------------------------------
      {/*:ETProperties 
        string ":NumberRecords" int amount
        :Parameters, NAME_1,NAME_2,...,NAME_N
        {int ID, double numCELL, double PET_MAX, double EXT_Depth
        :EndETProperties*/
        if (Options.noisy) {cout <<"ET Properties..."<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        done = false;
        int i = 0;
        while(!done)
        { 
          
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Parameters")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":NumberRecords")){
            string nrETe = "numrec_ET";
            pGWSP->SetGWProperty(nrETe,s_to_d(s[1]),0);
            pGWSP->AllocateArrays(nrETe,s_to_i(s[1]));
            //pGWClasses[num_parsed_gw-1]->SetGWProperty(nrETe,s_to_d(s[1]),0);
            //parsed_gw [num_parsed_gw].sp_props.numrec_ET = s_to_d(s[1]);
          }
          else if (Len==4)
          {
            if (num_parsed_gw>=MAX_SP_CLASSES-1){
              ExitGracefully("ParseGWStressPeriodClassFile: exceeded maximum # of SP classes",BAD_DATA);}
            
            string ncETe = "numcell_ET";
            string pmaxe = "pet_max";
            string edepe = "ext_depth";
            pGWSP->SetGWProperty(ncETe,s_to_d(s[1]),i);
            pGWSP->SetGWProperty(pmaxe,s_to_d(s[2]),i);
            pGWSP->SetGWProperty(edepe,s_to_d(s[3]),i);
            //pGWClasses[num_parsed_gw-1]->SetGWProperty(ncETe,s_to_d(s[1]),i);
            //pGWClasses[num_parsed_gw-1]->SetGWProperty(pmaxe,s_to_d(s[2]),i);
            //pGWClasses[num_parsed_gw-1]->SetGWProperty(edepe,s_to_d(s[3]),i);
            //parsed_gw [num_parsed_gw].sp_props.numcell_ET[i] =s_to_d(s[1]); 
            //parsed_gw [num_parsed_gw].sp_props.pet_max[i]    =s_to_d(s[2]); 
            //parsed_gw [num_parsed_gw].sp_props.ext_depth[i]  =s_to_d(s[3]);
            i++;
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndETProperties")){
            done=true;
            cout << "    "  << i << " ET PROPS Read..." << endl;
          }
        }
        
        break;
      }
      case(101):  //----------------------------------------------
      {/*:DrainProperties 
        string ":NumberRecords" int amount
        :Parameters, NAME_1,NAME_2,...,NAME_N
        {int ID, double numCELL, double DRAIN_ELEVATION, double CONDUCTANCE
        :EndETProperties*/
        if (Options.noisy) {cout <<"Drain Properties..."<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        done = false;
        int i = 0;
        while(!done)
        { 
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Parameters")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":NumberRecords")){
            string nrde = "numrec_D";
            pGWSP->SetGWProperty(nrde,s_to_d(s[1]),0);
            pGWSP->AllocateArrays(nrde,s_to_i(s[1]));
            //pGWClasses[num_parsed_gw-1]->SetGWProperty(nrde,s_to_d(s[1]),0);
            //parsed_gw [num_parsed_gw].sp_props.numrec_D = s_to_d(s[1]);
          }
          else if (Len==4)
          {
            if (num_parsed_gw>=MAX_SP_CLASSES-1){
              ExitGracefully("ParseGWStressPeriodClassFile: exceeded maximum # of SP classes",BAD_DATA);}
            
            string ncDe  = "numcell_Drains";
            string drele = "drain_elev";
            string conde = "conductance";

            pGWSP->SetGWProperty(ncDe,s_to_d(s[1]),i);
            pGWSP->SetGWProperty(drele,s_to_d(s[2]),i);
            pGWSP->SetGWProperty(conde,s_to_d(s[3]),i);
            //pGWClasses[num_parsed_gw-1]->SetGWProperty(ncDe,s_to_d(s[1]),i);
            //pGWClasses[num_parsed_gw-1]->SetGWProperty(drele,s_to_d(s[2]),i);
            //pGWClasses[num_parsed_gw-1]->SetGWProperty(conde,s_to_d(s[3]),i);
            //parsed_gw [num_parsed_gw].sp_props.numcell_Drains[i]=s_to_d(s[1]); 
            //parsed_gw [num_parsed_gw].sp_props.drain_elev[i]    =s_to_d(s[2]); 
            //parsed_gw [num_parsed_gw].sp_props.conductance[i]   =s_to_d(s[3]);
            i++;
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndDrainProperties")){
            done=true;
            cout << "    "  << i << " DRAIN PROPS Read..." << endl;
          }
        }
                
        break;
      }
      case(102):  //----------------------------------------------
      {/*:Recharge 
        string ":NumberRecords" int amount
        :Parameters, NAME_1,NAME_2,...,NAME_N
        {int ID, double numCELL, double RECHARGE
        :EndRecharge*/
        if (Options.noisy) {cout <<"Recharge Properties..."<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        done = false;
        int i = 0;
        while(!done)
        { 
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Parameters")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":NumberRecords")){
            string nrRe = "numrec_R";
            pGWSP->SetGWProperty(nrRe,s_to_d(s[1]),0);
            pGWSP->AllocateArrays(nrRe,s_to_i(s[1]));
            //pGWClasses[num_parsed_gw-1]->SetGWProperty(nrRe,s_to_d(s[1]),0);
            //parsed_gw [num_parsed_gw].sp_props.numrec_R = s_to_d(s[1]);
          }
          else if (Len==3)
          {
            if (num_parsed_gw>=MAX_SP_CLASSES-1){
              ExitGracefully("ParseGWStressPeriodClassFile: exceeded maximum # of SP classes",BAD_DATA);}
            
            string ncRe = "numcell_Rech";
            string reche= "recharge";

            pGWSP->SetGWProperty(ncRe,s_to_d(s[1]),i);
            pGWSP->SetGWProperty(reche,s_to_d(s[2]),i);
            //pGWClasses[num_parsed_gw-1]->SetGWProperty(ncRe,s_to_d(s[1]),i);
            //pGWClasses[num_parsed_gw-1]->SetGWProperty(reche,s_to_d(s[2]),i);
            //parsed_gw [num_parsed_gw].sp_props.numcell_Rech[i]=s_to_d(s[1]); 
            //parsed_gw [num_parsed_gw].sp_props.recharge[i]    =s_to_d(s[2]); 
            i++;
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndRecharge")){
            done=true;
            cout << "    " << i << " RECHARGE CONDITIONS Read..." << endl;
          }
        }
        break;
      }
      case(103):  //----------------------------------------------
      {/*:HydraulicConductivity 
        string ":NumberRecords" int amount
        string ":Anisotropic" or ":Isotropic"
        :Parameters, NAME_1,NAME_2,...,NAME_N
        {int ID, double numCELL, double CONDUCTIVITY_X, double CONDUCTIVITY_Y, double CONDUCTIVITY_Z,
        :EndHydraulicConductivity*/
        if (Options.noisy) {cout <<"Hydraulic Conductivity..."<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        done = false;
        int i = 0;
        while(!done)
        { 
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Parameters")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Anisotropic")){
            ExitGracefully("STUB: Anisotropy is not yet integrated",BAD_DATA);
          }
          else if (!strcmp(s[0],":Isotropic")){
            //NOT IMPLEMENTED YET
          }
          else if (!strcmp(s[0],":NumberRecords")){
            string nrK = "numrec_K";
            pGWSP->SetGWProperty(nrK,s_to_d(s[1]),0);
            pGWSP->AllocateArrays(nrK,s_to_i(s[1]));
            //pGWClasses[num_parsed_gw-1]->SetGWProperty(nrRe,s_to_d(s[1]),0);
            //parsed_gw [num_parsed_gw].sp_props.numrec_R = s_to_d(s[1]);
          }
          else if (Len==3)
          {
            if (num_parsed_gw>=MAX_SP_CLASSES-1){
              ExitGracefully("ParseGWStressPeriodClassFile: exceeded maximum # of SP classes",BAD_DATA);}
            string ncK = "numcell_K";
            string Kval= "hyd_cond";
            
            pGWSP->SetGWProperty(ncK,s_to_d(s[1]),i);
            pGWSP->SetGWProperty(Kval,s_to_d(s[2]),i);
            //pGWClasses[num_parsed_gw-1]->SetGWProperty(ncRe,s_to_d(s[1]),i);
            //pGWClasses[num_parsed_gw-1]->SetGWProperty(reche,s_to_d(s[2]),i);
            //parsed_gw [num_parsed_gw].sp_props.numcell_Rech[i]=s_to_d(s[1]); 
            //parsed_gw [num_parsed_gw].sp_props.recharge[i]    =s_to_d(s[2]); 
            i++;
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndHydraulicConductivity")){
            done=true;
            cout << "    " << i << " HYDRAULIC CONDUCTIVITY CONDITIONS Read..." << endl;
          }
        }
        break;
      }
      case(104):  //----------------------------------------------
      {/*:Storage 
        string ":NumberRecords" int amount
        :Parameters, NAME_1,NAME_2,...,NAME_N
        {int ID, double numCELL, double SS
        :EndStorage*/
        if (Options.noisy) {cout <<"Storage Properties..."<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        done = false;
        int i = 0;
        while(!done)
        { 
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Parameters")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":NumberRecords")){
            string nrS = "numrec_S";
            pGWSP->SetGWProperty(nrS,s_to_d(s[1]),0);
            pGWSP->AllocateArrays(nrS,s_to_i(s[1]));
            //pGWClasses[num_parsed_gw-1]->SetGWProperty(nrRe,s_to_d(s[1]),0);
            //parsed_gw [num_parsed_gw].sp_props.numrec_R = s_to_d(s[1]);
          }
          else if (Len==3)
          {
            if (num_parsed_gw>=MAX_SP_CLASSES-1){
              ExitGracefully("ParseGWStressPeriodClassFile: exceeded maximum # of SP classes",BAD_DATA);}
            
            string ncS = "numcell_S";
            string Sval= "ss";

            pGWSP->SetGWProperty(ncS,s_to_d(s[1]),i);
            pGWSP->SetGWProperty(Sval,s_to_d(s[2]),i);
            //pGWClasses[num_parsed_gw-1]->SetGWProperty(ncRe,s_to_d(s[1]),i);
            //pGWClasses[num_parsed_gw-1]->SetGWProperty(reche,s_to_d(s[2]),i);
            //parsed_gw [num_parsed_gw].sp_props.numcell_Rech[i]=s_to_d(s[1]); 
            //parsed_gw [num_parsed_gw].sp_props.recharge[i]    =s_to_d(s[2]); 
            i++;
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndStorage")){
            done=true;
            cout << "    " << i << " STORAGE CONDITIONS Read..." << endl;
          }
        }
        break;
      }
    }//end switch
    end_of_file=p->Tokenize(s,Len);

    //return after file redirect
    if ((end_of_file) && (pMainParser!=NULL))
    {
      if (Options.noisy) {cout <<"...Return to main file"<<endl;}
      INPUT2.clear();INPUT2.close();
      delete p;
      p=pMainParser;
      pMainParser=NULL;
      end_of_file=p->Tokenize(s,Len);
    }
  } //end while (!end_of_file)

  INPUT.close();
  
  cout<<"Checking for Required Model Parameters..."<<endl;

 // RVPParameterWarning(aPmaster,aPCmaster,nPmaster,Options);

  //Check for existence of classes
  //--------------------------------------------------------------------------
  /*int gw_param_count(0);
  for (int ii=0;ii<nP;ii++)
  {
    if      (aPC[ii]==CLASS_GWSTRESSPERIOD      ){gw_param_count++;}
  }

  ExitGracefullyIf((gw_param_count>0) && (num_parsed_gw==0),
    "ParsePropertyFile: Groundwater parameters needed for model, but no gw classes have been defined",BAD_DATA_WARN);
  
  delete [] aP;
  delete [] aPC;*/
  
  cout<<"..Done Checking"<<endl;
   
  //Check if there is at least one class for each major classification genre
  //Won't be needed if GetParticipatingParamList() populated for all processes
  ExitGracefullyIf(num_parsed_gw==0,
    "No groundwater classes specified in .rvg file. Cannot proceed.",BAD_DATA_WARN);

  delete [] indices;
  for (int i=0;i<MAX_NUM_IN_CLASS;i++){delete [] properties[i];}delete [] properties;

  return true;
}

///////////////////////////////////////////////////////////////////
/// \brief Parses array of class properties
/// \details Parses array of class properties in the form (starting at CLASS_TAG1): \n
/// :Properties\n
///     CLASS_TAG1, v1, v2, v3, v4,... \n
///     CLASS_TAG2, v1, v2, v3, v4,... \n
///     ...\n
///     CLASS_TAGM, v1, v2, v3, v4,... \n
///   :EndProperties
/// \param *p [in] File parsing object
/// \param *indices [out] Index number of class (with reference to tags[] array)
/// \param **properties [out] array of model properties
/// \param &num_read [out] Number of rows read
/// \param *tags [in] Array of tag names
/// \param line_length [in] Number of expected columns in array
/// \param max_classes [in] Maximum allowable rows in array
/// \return True if operation is successful
//
bool ParseGWPropArray(CParser   *p,           //parser
                    int       *indices,     //output: index number of class (with reference to tags[] array)
                    double   **properties,  //output: array of properties
                    int       &num_read,    //output: number of rows read
                    string    *tags,        //array of tag names
                    const int  line_length, //number of expected columns in array
                    const int  max_classes) //maximum allowable rows in array
{
  int Len;
  char *s[MAXINPUTITEMS];
  p->Tokenize(s,Len);
  bool done=false;
  num_read=0;

  while (!done)
  { 
    if      (Len==0){}
    else if (!strcmp(s[0],"*")){} 
    else if (!strcmp(s[0],"#")){}
    else if (Len==line_length)
    {
      int index=DOESNT_EXIST;
      for (int i=0;i<max_classes;i++){
        //cout<<tags[i]<<"<-->"<<s[0]<<endl;
        if (!tags[i].compare(s[0])){index=i;break;}
      }
      if (index==DOESNT_EXIST){
        cout <<"bad tag:"<<s[0]<<endl;
        return true;
      }//invalid class index
      indices[num_read]=index;

      for (int j=1;j<line_length;j++){
        properties[num_read][j-1]=AutoOrDouble(s[j]);
      }
      num_read++;
    }
    else{
      p->ImproperFormat(s); break;
    }
    p->Tokenize(s,Len);
    if (!strncmp(s[0],":End",4)){done=true;}
  }
  return false;
}
///////////////////////////////////////////////////////////////////
/// \brief Given list of nP needed parameters aP[] of type aPC[], warns user if some are not found associated with classes
///
/// \param *aP [out] array of needed parameter names needed
/// \param *aPC [out] Class type (groundwater) corresponding to each parameter
/// \param &nP [out] Number of parameters in list (size of aP[] and aPC[])
/// \param &Options global options structure
//
/*void  RVGParameterWarning   (string  *aP, class_type *aPC, int &nP, const optStruct &Options) 
{
  for (int ii=0;ii<nP;ii++)
  {
    if (Options.noisy){cout<<"    checking availability of parameter "<<aP[ii]<<endl;}

    if (aPC[ii]==CLASS_GWSTRESSPERIOD){
      for (int c=0;c<CGWStressPeriodClass::GetNumGWClasses();c++){
        if (CGWStressPeriodClass::GetGWClass(c)->GetGWProperty(aP[ii])==NOT_SPECIFIED){
          string warning="ParseGWPropertyFile: required groundwater property "+aP[ii]+" not included in .rvg file for class "+CGWStressPeriodClass::GetGWClass(c)->GetTag();
          ExitGracefully(warning.c_str(),BAD_DATA_WARN);
        }
      }
    }
    else if (aPC[ii]==CLASS_GLOBAL){
      if (CGlobalParams::GetParameter(aP[ii])==NOT_SPECIFIED){
            
          string warning="ParsePropertyFile: required global parameter "+aP[ii]+" not included in .rvp file.";
          ExitGracefully(warning.c_str(),BAD_DATA_WARN);
        }
    }
  }
}*/
///////////////////////////////////////////////////////////////////
/// \brief Append list of nP needed parameters aP[] of type aPC[] to master list aPm[],aPCm[] of size aPm
/// \note may later want to only append non-repetitive parameters
///
/// \param *aPm [in/out] master array of needed parameter names 
/// \param *aPCm [in/out] Class type (groundwater) corresponding to each needed parameter in master list
/// \param &nPm [in/out] Number of parameters in master list (initial size of aPm[] and aPCm[])
/// \param *aP [in] array of needed parameter names to be appended
/// \param *aPC [in] Class type (groundwater) corresponding to each needed parameter in appended list
/// \param &nP [in] Number of parameters in list to be appended (size of aP[] and aPC[])
//
void  AddToMasterGWPropsParamList   (string        *&aPm, class_type       *&aPCm, int       &nPm, 
                              const string  *aP , const class_type *aPC , const int &nP) 
{
  if (nP==0){return;}
  
  string      *aPm_new=NULL;
  class_type *aPCm_new=NULL;

  aPm_new=new string     [nPm+nP];
  aPCm_new=new class_type [nPm+nP];
  if (aPCm_new==NULL){ExitGracefully("DynArrayAppend::Out of memory",OUT_OF_MEMORY);}

  for (int i=0;i<nPm;i++){
    aPm_new [i]=aPm [i];
    aPCm_new[i]=aPCm[i];
  }
  for (int i=0;i<nP;i++){
    aPm_new [nPm+i]=aP[i];
    aPCm_new[nPm+i]=aPC[i];
  }
  if (aPm!=NULL){delete [] aPm; aPm=NULL;}
  if (aPm!=NULL){delete [] aPCm;aPCm=NULL;}
  aPm =aPm_new; 
  aPCm=aPCm_new;
  nPm=nPm+nP;
}

  