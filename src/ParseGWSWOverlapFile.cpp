/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright © 2008-2019 the Raven Development Team
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Properties.h"
#include "Model.h"
#include "ParseLib.h"
#include "GlobalParams.h"
#include "GroundwaterClass.h"
#include <string>

bool ParseOEPropArray(CParser *p,int *indices,double **properties,
                    int &num_read,string *tags,const int line_length,const int max_classes);
void  RVVParameterWarning   (string  *aP, class_type *aPC, int &nP, const optStruct &Options);
void  AddToMasterOEPropsParamList   (string        *&aPm, class_type       *&aPCm, int       &nPm, 
                                     const string  *aP , const class_type *aPC , const int &nP); 

//////////////////////////////////////////////////////////////////
/// \brief Parses Overlap/Exchange properties file
/// \details model.rvv: input file that defines Subbasin to GW Cell connections \n
///
/// \param *&pModel [out] Reference to model object
/// \param &Options [out] Global model options information
/// \return True if operation is successful
//
bool ParseGWSWOverlapFile(CModel *&pModel, 
                  const optStruct &Options)
{
	//global_struct         global_template;
  //global_struct         parsed_globals;

  COverlapExchangeClass *pOE;

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

  //string               *oeName;


 // CGlobalParams::InitializeGlobalParameters(global_template,true);
  //CGlobalParams::InitializeGlobalParameters(parsed_globals,false); 

  //COverlapExchangeClass::InitializeOEProperties(parsed_oe[0],true);//zero-index oe is template - does not create class
  oetags    [0]="[DEFAULT]";
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
//******   val=pOE->GetOEProperty(parsed_oe[num_parsed_oe],aPmaster[p],p);
 //   if (val==NOT_NEEDED_AUTO){val=AUTO_COMPUTE;}
 //   else if (val==NOT_NEEDED){val=NOT_SPECIFIED;}
 //   CGWStressPeriodClass::SetGWProperty    (parsed_gw[0],aPmaster[p],val,0);
//******    pOE->SetOEProperty(parsed_oe[num_parsed_oe],aPmaster[p],val,p,0);
 }

  cout <<"======================================================"<<endl;
  cout << "Parsing Overlap Properties File " << Options.rvv_filename <<"..."<<endl;
  cout <<"======================================================"<<endl;

  int   code;
  bool  done;
  bool  ended(false);
  int   Len,line(0);
  char *s[MAXINPUTITEMS];

  ifstream INPUT;
  INPUT.open(Options.rvv_filename.c_str());  
  if (INPUT.fail()){ 
    cout << "ERROR opening file: "<< Options.rvv_filename<<endl; return false;}
  
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
    else if  (!strcmp(s[0],":HRUs"                   )){code=1;  }//REQUIRED
    else if  (!strcmp(s[0],":GWCells"                )){code=2;  }//REQUIRED
    else if  (!strcmp(s[0],":OverlapArea"            )){code=3;  }//REQUIRED


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
          string warn;
          warn=":RedirectToFile: Cannot find file "+filename;
          ExitGracefully(warn.c_str(),BAD_DATA);
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
      {/*:HRUs 
        {int number of HRUs/Subbasins connected to GW Cells} */
        if (Options.noisy) {cout <<"HRU Connections..."<<endl;}
        if (Len!=2){p->ImproperFormat(s); break;} 
        //p->Tokenize(s,Len);
        done = false;
        while(!done)
        { 
          
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Parameters")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (Len==2)
          {
            pOE->InitializeOEProperties(parsed_oe[num_parsed_oe-1],false); //GWMIGRATE - THIS SHOULD BE A STATIC ROUTINE!!
            pOEClasses[num_parsed_oe-1]=new COverlapExchangeClass();
            pOE = pOEClasses[num_parsed_oe-1];
            oetags    [num_parsed_oe-1]="DEFAULT";
            
            string par_name = "NHRU";
            pOE->SetOEProperty(par_name,s_to_i(s[1]));
          }
          else{
            p->ImproperFormat(s); break;
          }
          pModel->GetGroundwaterModel()->AddOE(pOE);
          done=true;
          cout << "    "  << s_to_i(s[1]) << " CONNECTED HRU/SUBBASINs Read..." << endl;
        }
        
        break;
      }
      case(2):  //----------------------------------------------
      {/*:GWCells 
        {int number of GW cells connected to SW system} */
        if (Options.noisy) {cout <<"Number of GW Cells in contact with SW System..."<<endl;}
        if (Len!=2){p->ImproperFormat(s); break;} 
        //p->Tokenize(s,Len);
        done = false;

        while(!done)
        { 
          
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Parameters")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (Len==2)
          {
            string par_name = "NGWCELLS";
            
            pOE->SetOEProperty(par_name,s_to_i(s[1]));
            pOE->AllocateOEArrays(s_to_i(s[1]));
          }
          else{
            p->ImproperFormat(s); break;
          }
          done=true;
          cout << "    "  << s_to_i(s[1]) << " CONNECTED GW CELLS Read..." << endl;
        }
        
        break;
      }
      case(3):  //----------------------------------------------
      {/*:OverlapArea 
        {HRU#, {% Overlap of GW Cell with HRU} x num GW cells}
        :EndOverlapArea*/
        if (Options.noisy) {cout <<"Overlap Properties..."<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        done = false;
        int i = 0;
        while(!done)
        { 
          string nHRU = "NHRU";
          string nGWC = "NGWCELLS";
          
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Parameters")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (Len==(pOE->GetOEProperty(nGWC))+1)
          {
            if (num_parsed_oe>MAX_OE_CLASSES){
              ExitGracefully("ParseGWSWOverlapFile: exceeded maximum # of OE classes",BAD_DATA);}
    
            for(int k = 0; k < (pOE->GetOEProperty(nGWC)); k++)
            {
              string par_name = "OVERLAP";
              pOE->SetOEProperty(par_name, s_to_d(s[k+1]), i, k);
              //cout<<"overlap: "<<s_to_d(s[k+1])<<"   i: "<<i<<"   k: "<<k<<endl;
            }
            i++;
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndOverlapArea")){
            done=true;
            cout << "    "  << i << " GW-SW OVERLAP PROPS Read..." << endl;
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
bool ParseOEPropArray(CParser   *p,           //parser
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
void  AddToMasterOEPropsParamList   (string        *&aPm, class_type       *&aPCm, int       &nPm, 
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
