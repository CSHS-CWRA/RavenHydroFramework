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

bool ParseGWGeoPropArray(CParser *p,int *indices,double **properties,
                    int &num_read,string *tags,const int line_length,const int max_classes);
void  RVDParameterWarning   (string  *aP, class_type *aPC, int &nP, const optStruct &Options);
void  AddToMasterGWGeoParamList   (string        *&aPm, class_type       *&aPCm, int       &nPm, 
                              const string  *aP , const class_type *aPC , const int &nP); 

//////////////////////////////////////////////////////////////////
/// \brief This method parses the groundwater geometry class properties .rvd file
///
/// \details Parses input \b model.rvd file that defines the geometry properties for GW \n
///   ":Overview" \n
///      int nodes, int layers, int horizontal conn, int vert dis, int matrix format
///   ":Basement"  \n
///      int confined/unconfined  \n
///   ":Layers"  \n
///      int nodes/layer \n
///   :TopElevations  \n
///      [n x m] array of elevations  \n
///   :BottomElevations  \n
///      [n x m] array of elevations  \n
///   :CellArea  \n
///      [n x m] array of vell areas  \n
///   :NumberConnections  \n
///      [n x m] array of number of connections  \n
///   :Connections  \n
///      int numCell, [array of cell numbers connected to numCell]  \n
///   :VerticalConnections  \n
///      int numCell, [array indicating vertical or horizontal connection to numCell]  \n
///   :CL12  \n
///      int numCell, [array of distances from cell centre to interface for numCell]  \n
///   :InterfaceArea  \n
///      int numCell, [array of interface areas for numCell]  \n
///
/// \param *&pModel [in] The input model object
/// \param &Options [in] Global model options information
/// \return Boolean variable indicating success of parsing
//
bool ParseGWGeometryFile(CModel         *&pModel,
                         const optStruct &Options)
{
  //global_struct         global_template;
  //global_struct         parsed_globals;

  int                   num_parsed_gwGeo=0;
  CGWGeometryClass     *pGWGeoClasses  [MAX_GW_CLASSES];
  gw_dis_struct         parsed_gwGeo   [MAX_GW_CLASSES];
  string                gwGeotags      [MAX_GW_CLASSES];

  const int             MAX_PROPERTIES_PER_LINE=20;

  int                   MAX_NUM_IN_CLASS=MAX_GW_CELLS;
  int                  *indices;
  double              **properties;
  string                aParamStrings[MAXINPUTITEMS];
  int                   nParamStrings=0;

  CGWGeometryClass *pGWGeomet=pModel->GetGroundwaterModel()->GetGWGeom();

  //CGlobalParams::InitializeGlobalParameters(global_template,true);
  //CGlobalParams::InitializeGlobalParameters(parsed_globals,false); 

  //CGWGeometryClass::InitializeGWGeoProperties(parsed_gwGeo[num_parsed_gwGeo],true);//zero-index gw is template - does not create class
  gwGeotags    [0]="[DEFAULT]";
  //num_parsed_gwGeo++;

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
 // AddToMasterGWGeoParamList(aPmaster,aPCmaster,nPmaster,aP,aPC,nP);

  // Check parameters for each hydrological process algorithm:
  //--------------------------------------------------------------------------
  for (int j=0;j<pModel->GetNumProcesses();j++)
  {
    //pModel->GetProcess(j)->GetParticipatingParamList(aP,aPC,nP);
    //AddToMasterGWGeoParamList(aPmaster,aPCmaster,nPmaster,aP,aPC,nP);
  }
  
  //Prior to this, all parameters default to NOT_NEEDED or NOT_NEEDED_AUTO
  //this overrides default values of required variables to NOT_SPECIFIED and AUTO_COMPUTE
  double val;
  for (int p=0;p<nPmaster;p++)
  {
    //val=CGWGeometryClass::GetGWGeoProperty(parsed_gwGeo[num_parsed_gwGeo],aPmaster[p]);
    val=pGWGeomet->GetGWGeoProperty(parsed_gwGeo[num_parsed_gwGeo],aPmaster[p]);
    //if (val==NOT_NEEDED_AUTO){val=AUTO_COMPUTE;}
    //else if (val==NOT_NEEDED){val=NOT_SPECIFIED;}
    //CGWGeometryClass::SetGWGeoProperty    (parsed_gwGeo[num_parsed_gwGeo],aPmaster[p],val);
    //CGWGeometryClass::SetGWGeoProperty    (parsed_gwGeo[num_parsed_gwGeo],aPmaster[p],val,0,0);
    pGWGeomet->SetGWGeoProperty    (parsed_gwGeo[num_parsed_gwGeo],aPmaster[p],val);
    pGWGeomet->SetGWGeoProperty    (parsed_gwGeo[num_parsed_gwGeo],aPmaster[p],val,0,0);
  }

  cout <<"======================================================"<<endl;
  cout << "Parsing Groundwater Geometry File " << Options.rvd_filename <<"..."<<endl;
  cout <<"======================================================"<<endl;

  int   code;
  bool  done;
  bool  ended(false);
  int   Len,line(0);
  char *s[MAXINPUTITEMS];

  ifstream INPUT;
  INPUT.open(Options.rvd_filename.c_str());  
  if (INPUT.fail()){ 
    cout << "ERROR opening file: "<< Options.rvd_filename<<endl; return false;}
  
  ExitGracefullyIf(pModel==NULL,
    "ParseClassGWGeometryFile: model has not yet been created.",BAD_DATA);

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
    100 thru 200 : Groundwater Geometry Properties
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
    //--------------------OPTION PARAMS --------------------------
    else if  (!strcmp(s[0],":Overview"               )){code=1;  }//REQUIRED

    //--------------------Geometry PARAMS -------------------------
    else if  (!strcmp(s[0],":Basement"              )){code=100;}
    else if  (!strcmp(s[0],":Layers"                )){code=101;}
    else if  (!strcmp(s[0],":TopElevations"         )){code=102;}
    else if  (!strcmp(s[0],":BottomElevations"      )){code=103;}
    else if  (!strcmp(s[0],":CellArea"              )){code=104;}
    else if  (!strcmp(s[0],":NumberConnections"     )){code=105;}
    else if  (!strcmp(s[0],":Connections"           )){code=106;}
    else if  (!strcmp(s[0],":VerticalConnections"   )){code=107;}
    else if  (!strcmp(s[0],":CL12"                  )){code=108;}
    else if  (!strcmp(s[0],":FaceWidth"             )){code=109;}
    else if  (!strcmp(s[0],":InterfaceArea"         )){code=110;}
    else if  (!strcmp(s[0],":TopographicMinimum"    )){code=111;}

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
        num_parsed_gwGeo++;
        if (Options.noisy) {/*cout <<"EOF"<<endl;*/} ended=true; break;
      }
      case(1):  //----------------------------------------------
      {/*:Overview */
        if (Options.noisy) {cout <<"Groundwater Geometry Class"<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        done = false;
        while(!done)
        {
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Attributes")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (Len==6)
          {
            if ((num_parsed_gwGeo>=(MAX_GW_CLASSES-1)) || (num_parsed_gwGeo<0)){
              ExitGracefully("ParseGWGeometryClassFile: exceeded maximum # of gw classes",BAD_DATA);}

            //CGWGeometryClass::InitializeGWGeoProperties(parsed_gwGeo[num_parsed_gwGeo],false);
            pGWGeomet->InitializeGWGeoProperties(parsed_gwGeo[num_parsed_gwGeo],false);
            pGWGeoClasses[num_parsed_gwGeo]=new CGWGeometryClass(s[0]);
            //pGWGeomet[num_parsed_gwGeo]=new CGWGeometryClass(s[0]);
            pGWGeomet = pGWGeoClasses[num_parsed_gwGeo];
            gwGeotags    [num_parsed_gwGeo]=s[0];

            string nodee  = "num_nodes";
            string layere = "num_layers";
            string nconne = "num_connections";
            string vdise  = "vert_dis";
            string mforme = "matrix_format";

            /*pGWGeoClasses[num_parsed_gwGeo]->SetGWGeoProperty(nodee,s_to_d(s[1]));
            pGWGeoClasses[num_parsed_gwGeo]->SetGWGeoProperty(layere,s_to_d(s[2]));
            pGWGeoClasses[num_parsed_gwGeo]->SetGWGeoProperty(nconne,s_to_d(s[3]));
            pGWGeoClasses[num_parsed_gwGeo]->SetGWGeoProperty(vdise,s_to_d(s[4]));
            pGWGeoClasses[num_parsed_gwGeo]->SetGWGeoProperty(mforme,s_to_d(s[5]));*/
            pGWGeomet->SetGWGeoProperty(nodee,s_to_d(s[1]));
            pGWGeomet->SetGWGeoProperty(layere,s_to_d(s[2]));
            pGWGeomet->SetGWGeoProperty(nconne,s_to_d(s[3]));
            pGWGeomet->SetGWGeoProperty(vdise,s_to_d(s[4]));
            pGWGeomet->SetGWGeoProperty(mforme,s_to_d(s[5]));

            pGWGeomet->AllocateArrays(s_to_i(s[1]));

            string head = "GW_HEAD";
            for(int i = 0; i < s_to_i(s[1]); i++)
            {
              pGWGeomet->SetGWGeoProperty(head, 0, i, 0);
            }
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndOverview")){
            //pModel->AddGWGeo(pGWGeoClasses[num_parsed_gwGeo]);
            pModel->GetGroundwaterModel()->SetGWGeometry(pGWGeomet);
            done=true;
            cout << "    " << "GW GEOMETRY " << gwGeotags[num_parsed_gwGeo] << " Read..." << endl;
          }
        }     
        break;
      }
      
      //===========================================================================================
      //===========================================================================================
      case(100):  //----------------------------------------------
      {/*:Basement */
        if (Options.noisy) {cout <<"Geometry Basement"<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        done = false;
        while(!done)
        {
          if      (IsComment(s[0], Len)){}//comment line
          else if (Len==1)
          {
            string basee = "basement";
            //pGWGeoClasses[num_parsed_gwGeo]->SetGWGeoProperty(basee,s_to_d(s[0]));
            pGWGeomet->SetGWGeoProperty(basee,s_to_d(s[0]));
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndBasement")){done=true;}
        }
        break;
      }
      case(101):  //----------------------------------------------
      {/*:Layers */
        if (Options.noisy) {cout <<"Geometry Nodes per Layer"<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        done = false;
        while(!done)
        {
          if      (IsComment(s[0], Len)){}//comment line
          else if (Len==1)
          {
            string nodlaye = "nodes_layer";
            //pGWGeoClasses[num_parsed_gwGeo]->SetGWGeoProperty(nodlaye,s_to_d(s[0]));
            pGWGeomet->SetGWGeoProperty(nodlaye,s_to_d(s[0]));
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndLayers")){done=true;}
        }
        break;
      }
      case(102):  //----------------------------------------------
      {/*:TopElevation */
        if (Options.noisy) {cout <<"Geometry Top Elevation"<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        int i = 0;
        done = false;
        while(!done)
        { 
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Attributes")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (Len==2)
          {
            string tope = "top_elev";
            //pGWGeoClasses[num_parsed_gwGeo]->SetGWGeoProperty(tope,s_to_d(s[1]),i,0);
            pGWGeomet->SetGWGeoProperty(tope,s_to_d(s[1]),i,0);
            i++;
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndTopElevations")){done=true;}
        }
        break;
      }
      case(103):  //----------------------------------------------
      {/*:BottomElevation */
        if (Options.noisy) {cout <<"Geometry Bottom Elevation"<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        int i = 0;
        done=false;
        while(!done)
        { 
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Attributes")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (Len==2)
          {
            string bote = "bot_elev";
            //pGWGeoClasses[num_parsed_gwGeo]->SetGWGeoProperty(bote,s_to_d(s[1]),i,0);
            pGWGeomet->SetGWGeoProperty(bote,s_to_d(s[1]),i,0);
            i++;
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndBottomElevations")){done=true;}
        }
        break;
      }
      case(104):  //----------------------------------------------
      {/*:CellArea */
        if (Options.noisy) {cout <<"Geometry Cell Area"<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        int i = 0;
        done = false;
        while(!done)
        { 
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Attributes")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (Len==2)
          {
            string careae = "cell_area";
            //pGWGeoClasses[num_parsed_gwGeo]->SetGWGeoProperty(careae,s_to_d(s[1]),i,0);
            pGWGeomet->SetGWGeoProperty(careae,s_to_d(s[1]),i,0);
            i++;
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndCellArea")){done=true;}
        }
        break;
      }
      case(105):  //----------------------------------------------
      {/*:NumberConnections */
        if (Options.noisy) {cout <<"Geometry Connections per Node"<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        int i = 0;
        done = false;
        while(!done)
        { 
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Attributes")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (Len==2)
          {
            string ncconne = "num_cell_connections";
            //pGWGeoClasses[num_parsed_gwGeo]->SetGWGeoProperty(ncconne,s_to_d(s[1]),i,0);
            pGWGeomet->SetGWGeoProperty(ncconne,s_to_d(s[1]),i,0);
            i++;
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndNumberConnections")){done=true;}
        }
        break;
      }
      case(106):  //----------------------------------------------
      {/*:Connections */
        if (Options.noisy) {cout <<"Geometry Cell Connections"<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        int i = 0;
        done = false;
        while(!done)
        { 
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Attributes")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (Len>1)
          {
            int j = 0;
            //double numConnMax = pGWGeoClasses[num_parsed_gwGeo]->GetGWGeoProperty("num_cell_connections",i,0)-1;
            double numConnMax = pGWGeomet->GetGWGeoProperty("num_cell_connections",i,0)-1;
            while(j < numConnMax)
            {
              string cconne = "cell_connections";
              //pGWGeoClasses[num_parsed_gwGeo]->SetGWGeoProperty(cconne,s_to_d(s[j+1]),i,j);
              pGWGeomet->SetGWGeoProperty(cconne,s_to_d(s[j+1]),i,j);
              j++;
            }
            i++;
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndConnections")){done=true;}
        }
        break;
      }
      case(107):  //----------------------------------------------
      {/*:VerticalConnections */
        if (Options.noisy) {cout <<"Geometry Vertical Connections"<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        int i = 0;
        done = false;
        while(!done)
        { 
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Attributes")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (Len>1)
          {
            int j = 0;
            //double numConnMax = pGWGeoClasses[num_parsed_gwGeo]->GetGWGeoProperty("num_cell_connections",i,0)-1;
            double numConnMax = pGWGeomet->GetGWGeoProperty("num_cell_connections",i,0)-1;
            while(j < numConnMax)
            {
              string vcconne = "vert_cell_connections";
              //pGWGeoClasses[num_parsed_gwGeo]->SetGWGeoProperty(vcconne,s_to_d(s[j+1]),i,j);
              pGWGeomet->SetGWGeoProperty(vcconne,s_to_d(s[j+1]),i,j);
              j++;
            }
            i++;
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndVerticalConnections")){done=true;}
        }
        break;
      }
      case(108):  //----------------------------------------------
      {/*:CL12 */
        if (Options.noisy) {cout <<"Geometry Cell Distance - Node to Interface"<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        int i = 0;
        done = false;
        while(!done)
        { 
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Attributes")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (Len>1)
          {
            int j = 0;
            //double numConnMax = pGWGeoClasses[num_parsed_gwGeo]->GetGWGeoProperty("num_cell_connections",i,0)-1;
            double numConnMax = pGWGeomet->GetGWGeoProperty("num_cell_connections",i,0)-1;
            while(j < numConnMax)
            {
              string nintle = "node_interface_length";
              //pGWGeoClasses[num_parsed_gwGeo]->SetGWGeoProperty(nintle,s_to_d(s[j+1]),i,j);
              pGWGeomet->SetGWGeoProperty(nintle,s_to_d(s[j+1]),i,j);
              j++;
            }
            i++;
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndCL12")){done=true;}
        }
        break;
      }
      case(109):  //----------------------------------------------
      {/*:FaceWidth */
        if (Options.noisy) {cout <<"Geometry Interface Width"<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        int i = 0;
        done = false;
        while(!done)
        { 
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Attributes")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (Len>1)
          {
            int j = 0;
            //double numConnMax = pGWGeoClasses[num_parsed_gwGeo]->GetGWGeoProperty("num_cell_connections",i,0)-1;
            double numConnMax = pGWGeomet->GetGWGeoProperty("num_cell_connections",i,0)-1;
            while(j < numConnMax)
            {
              string nintle = "face_width";
              //pGWGeoClasses[num_parsed_gwGeo]->SetGWGeoProperty(nintle,s_to_d(s[j+1]),i,j);
              pGWGeomet->SetGWGeoProperty(nintle,s_to_d(s[j+1]),i,j);
              j++;
            }
            i++;
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndFaceWidth")){done=true;}
        }
        break;
      }
      case(110):  //----------------------------------------------
      {/*:InterfaceArea */
        if (Options.noisy) {cout <<"Geometry Interface Area"<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        int i = 0;
        done = false;
        while(!done)
        { 
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Attributes")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (Len>1)
          {
            int j = 0;
            //double numConnMax = pGWGeoClasses[num_parsed_gwGeo]->GetGWGeoProperty("num_cell_connections",i,0)-1;
            double numConnMax = pGWGeomet->GetGWGeoProperty("num_cell_connections",i,0)-1;
            while(j < numConnMax)
            {
              string intare = "interface_area";
              //pGWGeoClasses[num_parsed_gwGeo]->SetGWGeoProperty(intare,s_to_d(s[j+1]),i,j);
              pGWGeomet->SetGWGeoProperty(intare,s_to_d(s[j+1]),i,j);
              j++;
            }
            i++;
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndInterfaceArea")){done=true;}
        }
        break;
      }
      case(111):  //----------------------------------------------
      {/*:TopographicMinimum */
        if (Options.noisy) {cout <<"Geometry Topographic Minimum"<<endl;}
        if (Len!=1){p->ImproperFormat(s); break;} 
        p->Tokenize(s,Len);
        int i = 0;
        done = false;
        while(!done)
        { 
          if      (IsComment(s[0], Len)){}//comment line
          else if (!strcmp(s[0],":Attributes")){}//columns are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (Len==2)
          {
            string tope = "min_topog";
            //pGWGeoClasses[num_parsed_gwGeo]->SetGWGeoProperty(tope,s_to_d(s[1]),i,0);
            pGWGeomet->SetGWGeoProperty(tope,s_to_d(s[1]),i,0);
            i++;
          }
          else{
            p->ImproperFormat(s); break;
          }
          p->Tokenize(s,Len);
          if (!strcmp(s[0],":EndTopographicMinimum")){done=true;}
        }
        break;
      }

//===========================================================================================
//===========================================================================================

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
  /*int gwGeo_param_count(0);
  for (int ii=0;ii<nP;ii++)
  {
    if      (aPC[ii]==CLASS_GWGEOMETRY      ){gwGeo_param_count++;}
  }

  ExitGracefullyIf((gwGeo_param_count>0) && (num_parsed_gwGeo==0),
    "ParsePropertyFile: Groundwater parameters needed for model, but no gw classes have been defined",BAD_DATA_WARN);
  
  delete [] aP;
  delete [] aPC;
  */
  cout<<"..Done Checking"<<endl;
   
  //Check if there is at least one class for each major classification genre
  //Won't be needed if GetParticipatingParamList() populated for all processes
  ExitGracefullyIf(num_parsed_gwGeo==0,
    "No groundwater classes specified in .rvd file. Cannot proceed.",BAD_DATA_WARN);

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
bool ParseGWGeoPropArray(CParser   *p,           //parser
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
void  RVDParameterWarning   (string  *aP, class_type *aPC, int &nP, const optStruct &Options) 
{
  for (int ii=0;ii<nP;ii++)
  {
    if (Options.noisy){cout<<"    checking availability of parameter "<<aP[ii]<<endl;}

    if (aPC[ii]==CLASS_GWGEOMETRY){
      for (int c=0;c<CGWGeometryClass::GetNumGWGeoClasses();c++){
        if (CGWGeometryClass::GetGWGeoClass(c)->GetGWGeoProperty(aP[ii])==NOT_SPECIFIED){
          string warning="ParseGWGeometryFile: required groundwater property "+aP[ii]+" not included in .rvg file for class "+CGWGeometryClass::GetGWGeoClass(c)->GetTag();
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
}
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
void  AddToMasterGWGeoParamList   (string        *&aPm, class_type       *&aPCm, int       &nPm, 
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

  