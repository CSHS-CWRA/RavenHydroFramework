/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2019 the Raven Development Team
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Properties.h"
#include "Model.h"
#include "ParseLib.h"
#include "GlobalParams.h"
#include "GroundwaterModel.h"
#include "GWSWProcesses.h"
#include "GWRiverConnection.h"
#include <string>

/*
bool ParseGWGeoPropArray(CParser* p, int* indices, double** properties,
                       int& num_read, string* tags, const int line_length, const int max_classes);
void  RVDParameterWarning(string* aP, class_type* aPC, int& nP, const optStruct& Options);
void  AddToMasterGWGeoParamList(string*& aPm, class_type*& aPCm, int& nPm,
const string* aP, const class_type* aPC, const int& nP);
*/

//////////////////////////////////////////////////////////////////
/// \brief Parse main groundwater model file, model.rvg
///
/// \param *&pModel [in] The input model object
/// \param &Options [in] Global model options information
/// \return Boolean variable indicating success of parsing
//////////////////////////////////////////////////////////////////
bool ParseGWFile(CModel*& pModel, const optStruct& Options)
{
  //-- No Model No Go
  ExitGracefullyIf(pModel == NULL,"ParseClassGWFile: model has not yet been created.", BAD_DATA);

  int  code, Len, line(0);
  CGroundwaterModel*    pGWModel = pModel->GetGroundwaterModel();
  char *s[MAXINPUTITEMS];
  bool ended(false);
  bool done;

  //Open file
  ifstream RVG;
  RVG.open(Options.rvg_filename.c_str());
  if (RVG.fail()) {
    cout << "ERROR opening *.rvg file:" << Options.rvg_filename << endl; return false;
  }

  CParser* p = new CParser(RVG, Options.rvg_filename, line);

  ifstream INPUT2;             //For Secondary input
  ifstream MFUSG;              // For testing modflow name file

  CParser* pMainParser = NULL; //for storage of main parser while reading secondary files

  if (Options.noisy) {
    cout << "======================================================" << endl;
    cout << "Parsing Groundwater File " << Options.rvg_filename << "..." << endl;
    cout << "======================================================" << endl;
  }

  //Default Values
  //--Sift through file-----------------------------------------------
  bool end_of_file = p->Tokenize(s, Len);
  while (!end_of_file)
  {
    if (ended) { break; }
    if (Options.noisy) { cout << "reading line " << p->GetLineNumber() << ": "; }

    /*assign code for switch statement
      --------------------------------------------------------------------------------------
       <0          : ignored/special
        1 thru  99 : Options
      100 thru 199 : Geometry Params (e.g., Elevations, Area, Connections)
      200 thru 299 : Cell Params (e.g., Hydraulic Conductivity, Storage)
      300 thru 399 : Stress Period Params (
      --------------------------------------------------------------------------------------
    */

    code = 0;
    //---------------------SPECIAL ----------------------------------
    if       (Len == 0)                                {code=-1; }
    else if  (!strcmp(s[0], "*"                      )){code=-2; }//comment
    else if  (!strcmp(s[0], "%"                      )){code=-2; }//comment
    else if  (!strcmp(s[0], "#"                      )){code=-2; }//comment
    else if  (s[0][0] == '#'                         ) {code=-2; }//comment
    else if  (!strcmp(s[0], ":RedirectToFile"        )){code=-3; }//redirect to secondary file
    else if  (!strcmp(s[0], ":End"                   )){code=-4; }//premature end of file
    //--------------------Modflow-USG -------------------------------
    else if (!strcmp(s[0], ":NameFile"               )){code=1;  }//REQUIRED reads in MF model
    else if (!strcmp(s[0], ":RavenCBCUnitNumber"     )){code=2;  }
    else if (!strcmp(s[0], ":RechargeOptionCode"     )){code=3;  }
    //--------------------Model Connection PARAMS -------------------
    else if  (!strcmp(s[0], ":OverlapWeights"        )){code=125; }
    //--------------------GW-SW Process PARAMS ----------------------
    else if  (!strcmp(s[0],":Recharge"               )){code=305; }
    else if  (!strcmp(s[0],":ETProperties"           )){code=306; }
    else if  (!strcmp(s[0],":Drains"                 )){code=307; }
    //--------------------River (PBJ) PARAMS ----------------------
    else if (!strcmp(s[0], ":GWRiverConnection"      )){code=400; }

    switch (code)
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

        filename = CorrectForRelativePath(filename ,Options.rvg_filename);

        INPUT2.open(filename.c_str());
        if (INPUT2.fail()){
          string warn;
          warn=":RedirectToFile: Cannot find file "+filename;
          ExitGracefully(warn.c_str() ,BAD_DATA);
        }
        else{
          if (pMainParser != NULL) {
            ExitGracefully("ParseGWFile::nested :RedirectToFile commands (in already redirected files) are not allowed.",BAD_DATA);
          }
          pMainParser=p;    //save pointer to primary parser
          p=new CParser(INPUT2,filename,line);//open new parser
        }

        break;
      }
      case(-4):  //----------------------------------------------
      {/*:End*/
        if (Options.noisy) {/*cout <<"EOF"<<endl;*/} ended=true;
        break;
      }
      //===========================================================================================
      //===========================================================================================
      case(1):    //----------------------------------------------
      {/*NameFile:*/
        string filename = "";
        if (Options.noisy) { cout << "Namefile filename: " << s[1] << endl; }

        filename = CorrectForRelativePath(s[1], Options.rvg_filename);

        //-- Ensure Name File exists (easier to deal with here than within Fortran)
        MFUSG.open(filename.c_str());
        if (MFUSG.fail()){
          string warn;
          warn= "ParseGWFile: Cannot find Modflow Name file " +filename;
          ExitGracefully(warn.c_str(), BAD_DATA);
        }
        MFUSG.clear();
        MFUSG.close();

        pGWModel->InitMFUSG(CorrectForRelativePath(s[1], Options.rvg_filename));
        break;
      }

      case(2):    //----------------------------------------------
      {/*:RavenCBCUnitNumber
       MODFLOW-USG Unit number to which cell-by-cell flows will be written
       (Should correspond to a data file in the MF-USG name file) */
        int unit_no;
        if (Options.noisy) { cout << "Raven Cell-by-Cell Flow Output Unit Number: " << s[1] << endl; }
        if (Len >= 2) {
          unit_no = s_to_i(s[1]);
          pGWModel->SetCBCUnitNo(unit_no);
        }
        else {
          p->ImproperFormat(s); break;
        }
      break;
      }

      case (3):    //----------------------------------------------
      {/*:RechargeOptionCode - See Groundwater Model for details*/
        int rech_opt;
        if (Options.noisy) { cout << "Raven Recharge Option Code: " << s[1] << endl; }
        if (Len >= 2) {
          rech_opt = s_to_i(s[1]);
          if (rech_opt == 1 || rech_opt == 3) {  // This perhaps should be a string: 1 = TOP, 3 = ACTIVE
            pGWModel->SetRCHOptCode(rech_opt);
          } else {
            ExitGracefully("Invalid Recharge Option Code - Must be 1 or 3", BAD_DATA);
          }
        }
        else {
          p->ImproperFormat(s); break;
        }
      break;
      }

    //===========================================================================================
    //===========================================================================================

      case(125):  //----------------------------------------------
      {/*OverlapWeights - based on Forcing Grid Weights
       # [HRU ID] [Cell #] [w_ki]
       #   w_ki =A_ki/A_i where A_ki is area of overlap between HRU k and cell i and A_i is area of cell i
       #  //GWUPDATE should cell weights be forced to add to 1??? Or HRU???
       # [HRUID] [CELL#] [w_kl]
          1       2       0.3
          1       3       0.5
          1       23      0.2
          2       2       0.4
          2       3       0.6
          3       5       1.0*/
        // \todo check for GW Model setup
        long long int HRUid;
        int           nodeid;
        int nNodesLayOne = pGWModel->GetNumNodesByLay(1);
        double  wnode;

        //-- Need to keep track of total weights assigned to nodes
        //   Could do it later (cleaner) but this way limits loops over nNodes
        double *cellweights;
        cellweights = new double[nNodesLayOne];
        for (int n = 0; n < nNodesLayOne; n++) { cellweights[n] = 0.0; } // initialize to 1

        if (Options.noisy) {cout <<"Overlap weights..."<<endl;}
        while (((Len==0) || (strcmp(s[0],":EndOverlapWeights"))) && (!(p->Tokenize(s,Len))))
        {
          if      (IsComment(s[0],Len))            {}//comment line
          else if (!strcmp(s[0],":Attributes"    )){}//ignored by Raven - needed for GUIs
          else if (!strcmp(s[0],":Units"         )){}//ignored by Raven - needed for GUIs
          else if (!strcmp(s[0],":EndOverlapWeights")){}//done
          else if (Len >= 3) {
            //-- Assume line of data - HRUID NODE WEIGHT
            HRUid = s_to_ll(s[0]);
            nodeid = s_to_i(s[1]);
            wnode = s_to_d(s[2]);

            // Don't bother with checking/adding if weight is zero
            if (wnode != 0.0) {
              //-- Check if HRU exists //--Could be slow, keep track of HRUs that have been checked?
              if (pModel->GetHRUByID(HRUid) == NULL) {
                cout << "Bad HRU id for cell " << nodeid << endl;
                ExitGracefully("ParseGWFile: :OverlapWeights Bad HRU id",BAD_DATA);
              }
              //-- Check if Cell ID exists
              if (nodeid > nNodesLayOne) {
                ExitGracefully("ParseGWFile: :OverlapWeights Bad node id (greater than Layer 1 node count)",BAD_DATA);
              }
              //-- Check if weight is reasonable
              if (wnode > 1.0 || wnode < 0.0) {
                cout << "Bad weight for cell " << nodeid << endl;
                ExitGracefully("ParseGWFile: :OverlapWeights cell weights must be between 0 and 1.0",BAD_DATA);
              }
              //-- Check if cell total weight is reasonable (within limits, should cutoff be user defined?)
              cellweights[nodeid - 1] += wnode;  // nodeid is 1-index
              if (cellweights[nodeid - 1] - 1.0 > 1e-8) {
                cout << "Bad total weight for cell " << nodeid << ". Total: " << cellweights[nodeid -1] << endl;
                ExitGracefully("ParseGWFile: :OverlapWeights cell weights add to above 1.0",BAD_DATA);
              }
              //-- All good? Set the weight in the sparse "matrix"
              pGWModel->SetOverlapWeight((int)(HRUid), nodeid, wnode, Options);
            }
          }
          else {
            p->ImproperFormat(s); break;
          }
        }
        //-- Warn about first layer cells with weights < 1 (within limits, like above, could be user-defined)
        for (int n = 0; n < nNodesLayOne; n++) {
          if (cellweights[n] - 1.0 > 1e-8) {
            string warn ="GW Model Overlap Weights for cell " + to_string(n+1) + " add to " + to_string(cellweights[n]);
            warn += " (Less than 1.0)";
            WriteWarning(warn,Options.noisy);
          }
        }

        //-- Delete pointer to cell weights
        delete [] cellweights; cellweights = NULL;
        break;
      }

      //===========================================================================================
      //===========================================================================================
      case(305):  //----------------------------------------------
      {/*:Recharge
       string ":NumberRecords" int amount
       :Parameters, NAME_1,NAME_2,...,NAME_N
       {int ID, double numCELL, double RECHARGE
       :EndRecharge*/
        int i=0, nodeid;

        gwrecharge_type rech_type=RECHARGE_FLUX;

        // Get Recharge obj from GWModel
        CGWRecharge *pGWRecharge = pGWModel->GetRechProcess();

        if (Options.noisy) { cout << "Recharge Properties..." << endl; }
        while (((Len == 0) || (strcmp(s[0], ":EndRecharge"))) && (!(p->Tokenize(s, Len))))
        {
          if (IsComment(s[0], Len)) {}               //comment line
          else if (!strcmp(s[0], ":Parameters")) {}  //ignored by Raven - needed for GUIs
          else if (!strcmp(s[0], ":Attributes")) {}//ignored by Raven - needed for GUIs
          else if (!strcmp(s[0], ":Units")) {}       //ignored by Raven - needed for GUIs
          else if (!strcmp(s[0], ":EndRecharge")) {} //done
          else if (!strcmp(s[0], ":NumNodes")) {
            if (Len >= 2) {
              pGWRecharge->InitializeRechargeClass(s_to_i(s[1]));
              rech_type = pGWRecharge->getRechargeType();
            }
          }
          else if (Len >= 2)
          {
            // Has the recharge process been initialized? IS THIS A GOOD WAY TO DO THIS?
            if (pGWRecharge->getNumNodes() <= 0) {
              ExitGracefully("ParseGWPropsFile:: Number of Recharge Nodes not set or set to zero", BAD_DATA);
            }
            // Read in properties
            nodeid = s_to_i(s[0]);

            // Check node inputs (LS - what about lower bound? e.g. n < 0)
            if (nodeid > pGWModel->GetNumNodesByLay(1)) {
              ExitGracefully("ParseGWFile: Recharge nodes must be in the first layer", BAD_DATA);
            }
            if (i >= pGWRecharge->getNumNodes()) {
              ExitGracefully("ParseGWFile: :Recharge Too many nodes", BAD_DATA);
            }

            // Inputs all good? New Recharge Node
            if (rech_type == RECHARGE_FLUX) {
              pGWRecharge->addNodeFlux(i, nodeid, s_to_d(s[1]));
            }
            i++;

          }
          else {
            if (pGWRecharge->getNumNodes() <= 0) {
              ExitGracefully("ParseGWPropsFile:: Number of Recharge Nodes not set or set to zero", BAD_DATA);
            }
            p->ImproperFormat(s); break;
          }
        }
        break;
      }
      case(306):  //----------------------------------------------
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
          else if (!strcmp(s[0], ":Attributes")) {}//ignored by Raven - needed for GUIs
          else if (!strcmp(s[0],":Units")){} //units are explicit within Raven - useful for GUIs
          else if (!strcmp(s[0],":NumberRecords")){
            //GWSPUPDATE
          }
          else if (Len==4)
          {


            //GWSPUPDATE Add ET values
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
      case(307):  //----------------------------------------------
      {/*Drains - Move water from GW model when specific cell elevation is reached,
                  limited by drain conductance value
       :Drains
          :Parameters, ELEVATION, CONDUCTIVITY
          {int ID, double drain_elev, double drain_cond}
       :EndDrains
       */
        int    i = 0;
        int    nodeid;
        double drain_elev;
        double drain_cond;

        // Get Drain obj from GWModel
        CGWDrain* pDrain = pGWModel->GetDrainProcess();

        if (Options.noisy) { cout << "GW Drains..." << endl; }
        while (((Len == 0) || (strcmp(s[0], ":EndDrains"))) && (!(p->Tokenize(s, Len))))
        {
          if (IsComment(s[0], Len)) {}//comment line
          else if (!strcmp(s[0], ":Parameters")) {}//ignored by Raven - needed for GUIs
          else if (!strcmp(s[0], ":Attributes")) {}//ignored by Raven - needed for GUIs
          else if (!strcmp(s[0], ":Units")) {}//ignored by Raven - needed for GUIs
          else if (!strcmp(s[0], ":EndDrains")) {}//done
          else if (!strcmp(s[0], ":NumDrains")) {
            if (Len >= 2) { pDrain->InitializeDrainClass(s_to_i(s[1])); }
          }
          else if (Len >= 3)
          {
            // Has the drain been initialized? IS THIS A GOOD WAY TO DO THIS?
            if (pDrain->getNumNodes() <= 0) {
              ExitGracefully("ParseGWPropsFile:: NumDrains not set or set to zero", BAD_DATA);
            }
            // Read in properties
            nodeid = s_to_i(s[0]); drain_elev = s_to_d(s[1]); drain_cond = s_to_d(s[2]);

            // Check inputs
            if (nodeid > pGWModel->GetNumNodesByLay(1)) {
              //LS - as implemented, GWSWProcesses must be in the first layer
              // (we don't have/require HRU-node overlap weights for lower layers)
              ExitGracefully("ParseGWFile: Drain nodes must be in the first layer", BAD_DATA);
            }
            if (i >= pDrain->getNumNodes()) {
              ExitGracefully("ParseGWFile: :Drains Too many drains",BAD_DATA);
            }

            // Inputs all good? Set drain
            pDrain->addDrain(i, nodeid, drain_elev, drain_cond);
            i++;

          } else {
            p->ImproperFormat(s); break;
          }
        }
        break;
      }
      case(400):  //----------------------------------------------
      {/*GWRiverConnection - Basin-level connection between specified GW cells and SW River
         River is specified by line segments within triangles that connect 3 GW cells.
         The connections are based upon segment start/end distances from the 3 cell centers,
         i.e. the barycentric coordinates of the triangles. Conductances are given at the
         start and end of the segments as well, either as a straight condutance
         (condtype: CONDUCTANCE) or per unit length (condtype: UNITCOND).
         The connection has two modes: Drain or River. River uses the basin river elevation
         in calculating the flow to/from the stream.
       :GWRiverConnection
          :Mode     RIVER
          :CondType UNITCOND
          :Attributes
          :Units
          {3x int Nodes, 6x double weights, double seg_length, 2x double elevation, 2x double conductivity, int basinID}
       :EndGWRiverConnection
       */
        int    i = 0, segid = 0, nNodes;
        int    nodeid[3];
        int    SBID, SBp;        // Potential poor variable naming - SBp is the SB index
        double node_weights[6];  // First 3 segment start (A), last 3 segment end (B)
        double elev[2];
        double cond[2];
        double length;

        CSubBasin  *pSB;
        nNodes = pGWModel->GetNumNodes();

        // Get GWRiverConnection object from GW Model
        CGWRiverConnection* pRiver = pGWModel->GetRiverConnection();

        if (Options.noisy) { cout << "GW River Connection..." << endl; }
        while (((Len == 0) || (strcmp(s[0], ":EndGWRiverConnection"))) && (!(p->Tokenize(s, Len))))
        {
          if (IsComment(s[0], Len))                        {}//comment line
          else if (!strcmp(s[0], ":Parameters"))           {}//ignored by Raven - needed for GUIs
          else if (!strcmp(s[0], ":Attributes"))           {}//ignored by Raven - needed for GUIs
          else if (!strcmp(s[0], ":Units"))                {}//ignored by Raven - needed for GUIs
          else if (!strcmp(s[0], ":EndGWRiverConnection")) {}//done
          else if (!strcmp(s[0], ":Segments"))  {
            if (Len >= 2) {
              pRiver->Initialize(s_to_i(s[1]), pModel->GetNumSubBasins()); }
          }
          else if (!strcmp(s[0], ":Mode")) {
            pRiver->SetMode(s[1]);
          }
          else if (!strcmp(s[0], ":CondType")) {
            pRiver->SetCondType(s[1]);
          }
          else if (Len >= 15)
          {
            // Read in (numerous) attributes
            for (i=0; i<3; i++) {
              nodeid      [i]   = s_to_i(s[i]);
              node_weights[i]   = s_to_d(s[3+i]);
              node_weights[i+3] = s_to_d(s[6+i]);
            }
            length = s_to_d(s[9]);
            for (i=0; i<2; i++) {
              elev[i] = s_to_d(s[10+i]);
              cond[i] = s_to_d(s[12+i]);
            }
            SBID = s_to_i(s[14]);

            // Check inputs
            for (i=0; i<3; i++) {
              if (nodeid[i] > nNodes) {
                ExitGracefully("ParseGWFile: :CGWRiverConnection Bad cell id (greater than node count)", BAD_DATA);
              }
            }
            // Check SubBasin ID
            pSB=pModel->GetSubBasinByID(SBID);
            if (pSB==NULL){
              ExitGracefully("ParseGWFile: Bad Sub-basin index in :CGWRiverConnection",BAD_DATA);
            }
            SBp = pModel->GetSubBasinIndex(SBID);

            // Inputs all good? Add Segment
            pRiver->AddSegment(segid, nodeid, node_weights, elev, cond, length, SBp);
            segid++;

          }
          else {
            p->ImproperFormat(s); break;
          }
        }
        break;
      }
      //===========================================================================================
      //===========================================================================================
      default://----------------------------------------------
      {
        char firstChar = *(s[0]);
        if (firstChar==':')
        {
          if     (!strcmp(s[0],":FileType"))    {if (Options.noisy){cout<<"Filetype"<<endl;}}//do nothing
          else if(!strcmp(s[0],":Application")) {if (Options.noisy){cout<<"Application"<<endl;}}//do nothing
          else if(!strcmp(s[0],":Version"))     {if (Options.noisy){cout<<"Version"<<endl;}}//do nothing
          else if(!strcmp(s[0],":WrittenBy"))   {if (Options.noisy){cout<<"WrittenBy"<<endl;}}//do nothing
          else if(!strcmp(s[0],":CreationDate")){if (Options.noisy){cout<<"CreationDate"<<endl;}}//do nothing
          else if(!strcmp(s[0],":SourceFile"))  {if (Options.noisy){cout<<"SourceFile"<<endl;}}//do nothing
          else if(!strcmp(s[0],":Name"))        {if (Options.noisy){cout<<"Name"<<endl;}}//do nothing
          else
          {
            string warn ="IGNORING unrecognized command: " + string(s[0])+ " in .rvg file";
            WriteWarning(warn,Options.noisy);
          }
        }
        else
        {
          string errString = "Unrecognized command in .rvg file:\n   " + string(s[0]);
          ExitGracefully(errString.c_str(),BAD_DATA_WARN);
        }
        break;
      }
      //===========================================================================================
      //===========================================================================================
    } // End of switch

    end_of_file=p->Tokenize(s,Len);

    //return after file redirect, if in secondary file
    if ((end_of_file) && (pMainParser!=NULL))
    {
      INPUT2.clear();
      INPUT2.close();
      delete p;
      p=pMainParser;
      pMainParser=NULL;
      end_of_file=p->Tokenize(s,Len);
    }
  } //end while (!end_of_file)

  delete p; p = NULL;

  return true;

}
