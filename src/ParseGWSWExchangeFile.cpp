/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright � 2008-2019 the Raven Development Team
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Properties.h"
#include "Model.h"
#include "ParseLib.h"
#include "GlobalParams.h"
#include "GroundwaterClass.h"
#include <string>

bool ParseExPropArray(CParser *p,int *indices,double **properties,
                    int &num_read,string *tags,const int line_length,const int max_classes);
void  RVEParameterWarning   (string  *aP, class_type *aPC, int &nP, const optStruct &Options);
void  AddToMasterExPropsParamList   (string        *&aPm, class_type       *&aPCm, int       &nPm, 
                                     const string  *aP , const class_type *aPC , const int &nP); 

//////////////////////////////////////////////////////////////////
/// \brief Parses Overlap/Exchange properties file
/// \details model.rvs: input file that defines exchange fluxes between SW and GW \n
///
/// \param *&pModel [out] Reference to model object
/// \param &Options [out] Global model options information
/// \return True if operation is successful
//
bool ParseGWSWExchangeFile(CModel *&pModel, 
           
                  const optStruct  &Options)
{
	//global_struct         global_template;
  //global_struct         parsed_globals;

  COverlapExchangeClass  *pOE=pModel->GetGroundwaterModel()->GetOEs(0);

  int                   num_parsed_oe=0;
  COverlapExchangeClass *pOEClasses  [MAX_OE_CLASSES];
  oe_struct             parsed_oe   [MAX_OE_CLASSES];
  string                oetags      [MAX_OE_CLASSES];

  const int             MAX_PROPERTIES_PER_LINE=20;

  int                   MAX_NUM_IN_CLASS=MAX_GW_CELLS;
  int                  *indices;
  double              **properties;
  string                aParamStrings[MAXINPUTITEMS];
  int                   nParamStrings=0;

//  string               *oeName;


 // CGlobalParams::InitializeGlobalParameters(global_template,true);
  //CGlobalParams::InitializeGlobalParameters(parsed_globals,false); 

  //COverlapExchangeClass::InitializeOEProperties(parsed_oe[0],true);//zero-index oe is template - does not create class
  //oetags    [0]="[DEFAULT]";
  num_parsed_oe++;

  indices   =new int     [MAX_NUM_IN_CLASS];
  properties=new double *[MAX_NUM_IN_CLASS];
  for (int i=0;i<MAX_NUM_IN_CLASS;i++){
    properties[i]=new double [MAX_PROPERTIES_PER_LINE];
  }

  
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
 for (int p=0;p<nPmaster;p++)
 {
 //   val=CGWStressPeriodClass::GetGWProperty(parsed_gw[0],aPmaster[p]);
   //*****val=pOE->GetOEProperty(parsed_oe[num_parsed_oe],aPmaster[p],p);
 //   if (val==NOT_NEEDED_AUTO){val=AUTO_COMPUTE;}
 //   else if (val==NOT_NEEDED){val=NOT_SPECIFIED;}
 //   CGWStressPeriodClass::SetGWProperty    (parsed_gw[0],aPmaster[p],val,0);
   //****** pOE->SetOEProperty(parsed_oe[num_parsed_oe],aPmaster[p],val,p,0);
 }

  cout <<"======================================================"<<endl;
  cout << "Parsing Exchange Properties File " << Options.rvs_filename <<"..."<<endl;
  cout <<"======================================================"<<endl;

  int   code;
  bool  done;
  bool  ended(false);
  int   Len,line(0);
  char *s[MAXINPUTITEMS];

  ifstream INPUT;
  INPUT.open(Options.rvs_filename.c_str());  
  if (INPUT.fail()){ 
    cout << "ERROR opening file: "<< Options.rvs_filename<<endl; return false;}
  
  ExitGracefullyIf(pModel==NULL,
    "ParseOEClassPropertiesFile: model has not yet been created.",BAD_DATA);

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
    0   thru 100 : Overlap Allocation
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
    else if  (!strcmp(s[0],":Cell"                   )){code=1;  }//REQUIRED
    else if  (!strcmp(s[0],":EndCell"                )){code=2;  }//REQUIRED
    else if  (!strcmp(s[0],":RelationshipCount"      )){code=3;  }//REQUIRED
    //--------------------EXCHANGE FLUX RELATIONSHIPS------------------
    else if  (!strcmp(s[0],":Evapotranspiration"     )){code=100;  }
    else if  (!strcmp(s[0],":Drain"                  )){code=101;  }
    else if  (!strcmp(s[0],":SaturatedArea"          )){code=102;  }
    


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
      {/*:Cell {String or number stating which cell this UFRs apply to} */
        if (Options.noisy) {cout <<"Exchange Relationships..."<<endl;}
        if (Len!=2){p->ImproperFormat(s); break;} 

        int cell_num;
        if (!strcmp(s[1],"ALL")){cell_num = -9999;}
        else                    {cell_num = s_to_i(s[1]);}

        string par_name = "CELL";               //check if ALL, if not get cell number
        pOE->SetOEProperty(par_name,cell_num);
        break;
      }
      case(2):  //----------------------------------------------
      {/*:EndCell*/
        if (Options.noisy) {cout <<"End Cell"<<endl;}
        break;
      }
      case(3):  //----------------------------------------------
      {/*:RelationshipCount {Number stating total number of UFRs in file} */
        if (Options.noisy) {cout <<"Relationship Count"<<endl;}
        if (Len!=2){p->ImproperFormat(s); break;} 

        pOE->SetOEProperty(to_string("RELATIONSHIPCOUNT"),s_to_i(s[1]));
        
        break;
      }
      case(100):  //----------------------------------------------
      {/*:Evapotranspiration 
        :Relationship_Type string type
        double coeff_a, double power_b  (if Power Law)*/
        if (Options.noisy) {cout <<" Exchange Flux Properties for Evapotranspiration..."<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        done = false;
        string par_cell = "CELL";
        while(!done)
        { 
          
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Parameters")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (Len==2)
          {
            string par_name = "ET_TYPE";
            if(!strcmp(s[0],":Relationship_Type"))
            {
              if     (!strcmp(s[1],"PowerLaw"      )){pOE->SetOEProperty(par_name,"POWER_LAW");}
              else{ExitGracefully("ParseGWSWExchangeFile: ET Relationship Type not defined",BAD_DATA);}
            }
            else
            {
              string test = pOE->GetOEProperty(par_name, par_name, 0);
              string par_name2 = "POWER_LAW";
              if(!strcmp(test.c_str(),par_name2.c_str()))
              {
                string ETa = "ET_A";
                string ETb = "ET_B";
                pOE->SetOEProperty(ETa,s_to_d(s[0]),(int)(pOE->GetOEProperty(par_cell))-1);     //****FIX - MASSIVE HACK JOB - FIX ASAP
                pOE->SetOEProperty(ETb,s_to_d(s[1]),(int)(pOE->GetOEProperty(par_cell))-1);
                //cout<<"a: "<<s[0]<<"   b: "<<s[1]<<endl;
              }
            }
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndEvapotranspiration")){
            done=true;
            //cout << "    Evapotranspiration Flux Relationship Read..." << endl;
          }
        }
        
        break;
      }
      case(101):  //----------------------------------------------
      {/*:Drain 
        :Relationship_Type string type
        double coeff_a, double power_b  (if Power Law)*/
        if (Options.noisy) {cout <<" Exchange Flux Properties for Drains..."<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        done = false;
        string par_cell = "CELL";
        while(!done)
        { 
          
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Parameters")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (Len==2)
          {
            string par_name = "D_TYPE";
            if(!strcmp(s[0],":Relationship_Type"))
            {
              if     (!strcmp(s[1],"PowerLaw"      )){pOE->SetOEProperty(par_name,"POWER_LAW");}
              else{ExitGracefully("ParseGWSWExchangeFile: Drain Relationship Type not defined",BAD_DATA);}
            }
            else
            {
              string test = pOE->GetOEProperty(par_name, par_name,0);
              string par_name2 = "POWER_LAW";
              //cout<<"name: "<<test<<endl;
              //if(par_name2.compare(test) == 0)
              if(!strcmp(test.c_str(),par_name2.c_str()))
              {
                string Da = "D_A";
                string Db = "D_B";
                pOE->SetOEProperty(Da,s_to_d(s[0]),(int)(pOE->GetOEProperty(par_cell))-1);
                pOE->SetOEProperty(Db,s_to_d(s[1]),(int)(pOE->GetOEProperty(par_cell))-1);
                //cout<<"a: "<<s[0]<<"   b: "<<s[1]<<endl;
              }
            }
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndDrain")){
            done=true;
            //cout << "    Drain Flux Relationship Read..." << endl;
          }
        }
        break;
      }
      case(102):  //----------------------------------------------
      {/*:SaturatedArea
        :Relationship_Type string type
        double coeff_a, double power_b  (if Power Law)*/
        if (Options.noisy) {cout <<" Exchange Flux Properties for Saturated Area..."<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        done = false;
        string par_cell = "CELL";
        while(!done)
        { 
          
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Parameters")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (Len==2)
          {
            string par_name = "SA_TYPE";
            if(!strcmp(s[0],":Relationship_Type"))
            {
              if     (!strcmp(s[1],"PowerLaw"      )){pOE->SetOEProperty(par_name,"POWER_LAW");}
              else{ExitGracefully("ParseGWSWExchangeFile: Saturated Area Relationship Type not defined",BAD_DATA);}
            }
            else
            {
              string test = pOE->GetOEProperty(par_name, par_name,0);
              string par_name2 = "POWER_LAW";
              if(!strcmp(test.c_str(),par_name2.c_str()))
              {
                string SAa = "SA_A";
                string SAb = "SA_B";
                pOE->SetOEProperty(SAa,s_to_d(s[0]),(int)(pOE->GetOEProperty(par_cell))-1);
                pOE->SetOEProperty(SAb,s_to_d(s[1]),(int)(pOE->GetOEProperty(par_cell))-1);
              }
            }
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndSaturatedArea")){
            done=true;
            //cout << "    Saturated Area Relationship Read..." << endl;
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
  ExitGracefullyIf(num_parsed_oe==0,
    "No overlap/exchange classes specified in .rvv file. Cannot proceed.",BAD_DATA_WARN);

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
bool ParseExPropArray(CParser   *p,           //parser
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
/*void  RVEParameterWarning   (string  *aP, class_type *aPC, int &nP, const optStruct &Options) 
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
void  AddToMasterExPropsParamList   (string        *&aPm, class_type       *&aPCm, int       &nPm, 
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
