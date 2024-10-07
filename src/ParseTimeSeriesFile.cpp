/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2024 the Raven Development Team
  ----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Model.h"
#include "TimeSeries.h"
#include "IrregularTimeSeries.h"
#include "ParseLib.h"

void AllocateReservoirDemand(CModel *&pModel,const optStruct &Options,long SBID, long SBIDres,double pct_met,int jul_start,int jul_end);
bool IsContinuousFlowObs2(const CTimeSeriesABC* pObs,long SBID);
void GetNetCDFStationArray(const int ncid, const string filename,int &stat_dimid,int &stat_varid, long *&aStations, string *&aStat_string,int &nStations);
//////////////////////////////////////////////////////////////////
/// \brief Parse input time series file, model.rvt
///
/// \param *&pModel [out] Reference to the model object
/// \param Options [in] Global model options
/// \return True if operation was successful
//
bool ParseTimeSeriesFile(CModel *&pModel, const optStruct &Options)
{

  CTimeSeries *pTimeSer;
  CGauge *pGage=NULL;
  CForcingGrid  *pGrid=NULL;

  int   code;
  int   Len,line(0);
  char *s[MAXINPUTITEMS];
  double monthdata[12];
  string warn;
  bool   in_ifmode_statement=false;
  int    ncid=-9;
  string const_name;

  bool ended            = false;
  bool has_irrig        = false;
  bool grid_initialized = false;
  bool is_3D            = false;  // true if gridded forcing is 3D

  ifstream RVT;
  ifstream INPUT2;           //For Secondary input
  CParser* pMainParser=NULL; //for storage of main parser while reading secondary files

  if (Options.in_bmi_mode && (strcmp(Options.rvt_filename.c_str(), "") == 0)) {  // an RVT may not be specified for a BMI run
    return (true);
  }

  RVT.open(Options.rvt_filename.c_str());
  if (RVT.fail()){
    cout << "ERROR opening *.rvt file: "<<Options.rvt_filename<<endl; return false;}

  CParser *p=new CParser(RVT,Options.rvt_filename,line);

  if (Options.noisy)
  {
    cout <<"==========================================================="<<endl;
    cout << "Parsing Forcing Data File " << Options.rvt_filename <<"..."<<endl;
    cout <<"==========================================================="<<endl;
  }

  //--Sift through file-----------------------------------------------
  bool end_of_file=p->Tokenize(s,Len);
  while (!end_of_file)
  {
    if (ended){break;}
    if (Options.noisy){ cout << "reading line " << p->GetLineNumber() << ": ";}

    code=0;
    //---------------------SPECIAL -----------------------------
    if       (Len==0)                                       {code=-1; }
    else if  (IsComment(s[0],Len))                          {code=-2; }//comment
    else if  (!strcmp(s[0],":End"                         )){code=-4; }//premature end of file
    else if  (!strcmp(s[0],":IfModeEquals"                )){code=-5; }
    else if  (in_ifmode_statement)                          {code=-6; }
    else if  (!strcmp(s[0],":EndIfModeEquals"             )){code=-2; }//treat as comment - unused mode
    else if  (!strcmp(s[0],":RedirectToFile"              )){code=-3; }//redirect to secondary file
    //--------------------GAUGE BASIC DATA- --------------------
    else if  (!strcmp(s[0],":Gauge"                       )){code=1;  }
    else if  (!strcmp(s[0],":EndGauge"                    )){code=2;  }
    else if  (!strcmp(s[0],":Latitude"                    )){code=9;  }
    else if  (!strcmp(s[0],":Longitude"                   )){code=10; }
    else if  (!strcmp(s[0],":Elevation"                   )){code=11; }
    else if  (!strcmp(s[0],":MeasurementHeight"           )){code=16; }
    //--------------------FORCING FUNCTIONS --------------------
    else if  (!strcmp(s[0],":Data"                        )){code=12; }
    else if  (!strcmp(s[0],":MultiData"                   )){code=13; }
    else if  (!strcmp(s[0],":GaugeDataTable"              )){code=15; }
    else if  (!strcmp(s[0],":EnsimTimeSeries"             )){code=20; }
    //----------GAUGE-SPECIFIC CORRECTION TERMS / PARAMETERS----
    else if  (!strcmp(s[0],":RainCorrection"              )){code=30; }
    else if  (!strcmp(s[0],":SnowCorrection"              )){code=31; }
    else if  (!strcmp(s[0],":CloudTempRanges"             )){code=32; }
    else if  (!strcmp(s[0],":TemperatureCorrection"       )){code=33;}
    //-------------------OBSERVATIONS---------------------------
    else if  (!strcmp(s[0],":ObservationData"             )){code=40; }
    else if  (!strcmp(s[0],":IrregularObservations"       )){code=41; }
    else if  (!strcmp(s[0],":ObservationWeights"          )){code=42; }
    else if  (!strcmp(s[0],":IrregularWeights"            )){code=43; }
    //----------HYDROGRAPHS /RESERVOIR PARAMS -------------------
    else if  (!strcmp(s[0],":BasinInflowHydrograph"       )){code=50; }
    else if  (!strcmp(s[0],":BasinInflowHydrograph2"      )){code=51; }
    else if  (!strcmp(s[0],":ReservoirExtraction"         )){code=52; }
    else if  (!strcmp(s[0],":ReservoirDemand"             )){code=52; }
    else if  (!strcmp(s[0],":VariableWeirHeight"          )){code=53; }
    else if  (!strcmp(s[0],":ReservoirMaxStage"           )){code=54; }
    else if  (!strcmp(s[0],":ReservoirMinStage"           )){code=55; }
    else if  (!strcmp(s[0],":ReservoirMinStageFlow"       )){code=56; }
    else if  (!strcmp(s[0],":ReservoirTargetStage"        )){code=57; }
    else if  (!strcmp(s[0],":ReservoirMaxQDelta"          )){code=58; }
    else if  (!strcmp(s[0],":ReservoirMaxQDecrease"       )){code=59; }
    else if  (!strcmp(s[0],":ReservoirMinFlow"            )){code=60; }
    else if  (!strcmp(s[0],":ReservoirMaxFlow"            )){code=61; }
    else if  (!strcmp(s[0],":OverrideReservoirFlow"       )){code=62; }
    //----------DEMAND PARAMS ----------------------------------
    else if  (!strcmp(s[0],":IrrigationDemand"            )){code=63; }
    else if  (!strcmp(s[0],":WaterDemand"                 )){code=63; }
    else if  (!strcmp(s[0],":ReservoirDownstreamFlow"     )){code=64; }
    else if  (!strcmp(s[0],":ReservoirDownstreamDemand"   )){code=65; }
    else if  (!strcmp(s[0],":FlowDiversion"               )){code=66; }
    else if  (!strcmp(s[0],":FlowDiversionLookupTable"    )){code=67; }
    else if  (!strcmp(s[0],":EnvironmentalMinFlow"        )){code=68; }
    else if  (!strcmp(s[0],":UnusableFlowPercentage"      )){code=69; }
    //--------------------Other --------------------------------
    else if  (!strcmp(s[0],":MonthlyAveTemperature"       )){code=70; }
    else if  (!strcmp(s[0],":MonthlyAveEvaporation"       )){code=71; }
    else if  (!strcmp(s[0],":MonthlyEvapFactor"           )){code=71; }
    else if  (!strcmp(s[0],":MonthlyMinTemperature"       )){code=72; }
    else if  (!strcmp(s[0],":MonthlyMaxTemperature"       )){code=73; }
    else if  (!strcmp(s[0],":UserTimeSeries"              )){code=74; }
    //-----------------CONTROLS ---------------------------------
    else if  (!strcmp(s[0],":OverrideStreamflow"          )){code=100;}
    else if  (!strcmp(s[0],":AssimilateStreamflow"        )){code=101;}
    else if  (!strcmp(s[0],":SpecifyGaugeForHRU"          )){code=102;}
    //-----------------TRANSPORT--------------------------------
    else if  (!strcmp(s[0],":FixedConcentrationTimeSeries")){code=300; const_name=s[1];}
    else if  (!strcmp(s[0],":FixedTemperatureTimeSeries"  )){code=300; const_name="TEMPERATURE";}
    else if  (!strcmp(s[0],":MassFluxTimeSeries"          )){code=301;}
    else if  (!strcmp(s[0],":SpecifiedInflowConcentration")){code=302;}
    else if  (!strcmp(s[0],":SpecifiedInflowTemperature"  )){code=303;}
    else if  (!strcmp(s[0],":MassLoading"                 )){code=304;}
    //---------GRIDDED INPUT (lat,lon,time)---------------------
    else if  (!strcmp(s[0],":GriddedForcing"              )){code=400;}
    else if  (!strcmp(s[0],":ForcingType"                 )){code=401;}
    else if  (!strcmp(s[0],":FileNameNC"                  )){code=402;}
    else if  (!strcmp(s[0],":VarNameNC"                   )){code=403;}
    else if  (!strcmp(s[0],":DimNamesNC"                  )){code=404;}
    else if  (!strcmp(s[0],":GridWeights"                 )){code=405;}
    else if  (!strcmp(s[0],":EndGriddedForcing"           )){code=406;}
    else if  (!strcmp(s[0],":Deaccumulate"                )){code=407;}
    else if  (!strcmp(s[0],":TimeShift"                   )){code=408;}
    else if  (!strcmp(s[0],":LinearTransform"             )){code=409;}
    else if  (!strcmp(s[0],":PeriodEndingNC"              )){code=410;}
    else if  (!strcmp(s[0],":LatitudeVarNameNC"           )){code=411;}
    else if  (!strcmp(s[0],":LongitudeVarNameNC"          )){code=412;}
    else if  (!strcmp(s[0],":ElevationVarNameNC"          )){code=413;}
    else if  (!strcmp(s[0],":GridWeightsByAttribute"      )){code=414;}
    else if  (!strcmp(s[0],":StationElevationsByAttribute")){code=415;}
    else if  (!strcmp(s[0],":StationIDNameNC"             )){code=416;}
    else if  (!strcmp(s[0],":StationElevationsByIdx"      )){code=417;}
    else if  (!strcmp(s[0],":MapStationsTo"               )){code=418;}//Alternate to :GridWeights for :StationForcing command

    //---------STATION DATA INPUT AS NETCDF (stations,time)------
    //             code 401-405 & 407-414 are shared between
    //             GriddedForcing and StationForcing
    else if  (!strcmp(s[0],":StationForcing"              )){code=500;}
    else if  (!strcmp(s[0],":EndStationForcing"           )){code=506;}

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

      filename =CorrectForRelativePath(filename ,Options.rvt_filename);

      INPUT2.open(filename.c_str());
      if (INPUT2.fail()){
        warn=":RedirectToFile: Cannot find file "+filename;
        ExitGracefully(warn.c_str(),BAD_DATA);
      }
      else{
        if (pMainParser != NULL) {
          ExitGracefully("ParseTimeSeriesFile::nested :RedirectToFile commands (in already redirected files) are not allowed.",BAD_DATA);
        }
        pMainParser=p;    //save pointer to primary parser
        p=new CParser(INPUT2,filename,line);//open new parser
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
      if(Options.noisy) { cout <<"* skip"<<endl; }
      if(!strcmp(s[0],":EndIfModeEquals"))
      {
        if(Options.noisy) { cout <<"...Mode statement end"<<endl; }
        in_ifmode_statement=false;
      }
      break;
    }
    case(1):  //----------------------------------------------
    {/*:Gauge [optional name (no spaces)]*/
      if (Options.noisy) {cout <<"Gauge Information Begin"<<endl;}
      char tmp[64];
      sprintf(tmp, "Gauge%i", pModel->GetNumGauges()+1);
      string gName = tmp;
      if (Len > 1){                 // A gauge name may have been specified
        if ((char)*(s[1]) != '#'){  // first character a '#' comment character?
          gName = s[1];             // Valid gauge name, so we use it
        }
      }
      pGage=new CGauge(gName,NOT_SPECIFIED,NOT_SPECIFIED,0.0);
      pModel->AddGauge(pGage);
      break;
    }
    case(2):  //----------------------------------------------
    {/*:EndGauge*/
      if (Options.noisy) {cout <<"End Gauge"<<endl;}
      if(pGage == NULL)
      {
        ExitGracefully("ParseTimeSeriesFile:: :EndGauge encountered before :Gauge",BAD_DATA);
      }
      else
      {
        ExitGracefullyIf(pGage->GetLocation().latitude==NOT_SPECIFIED,
                         "ParseTimeSeriesFile::latitude not specified for gauge station",BAD_DATA);
        ExitGracefullyIf(pGage->GetLocation().longitude==NOT_SPECIFIED,
                         "ParseTimeSeriesFile::longitude not specified for gauge station",BAD_DATA);
      }
      break;
    }
    case(9):  //----------------------------------------------
    {/*:Latitude*/
      if (Options.noisy) {cout <<"Gauge latitude"<<endl;}
      ExitGracefullyIf(pGage==NULL,
                       "ParseTimeSeriesFile::Latitude specified outside of a :Gauge-:EndGauge statement",BAD_DATA);
      ExitGracefullyIf(Len<=1,"ParseTimeSeriesFile - no value provided after :Latitude command",BAD_DATA_WARN);
      pGage->SetLatitude(s_to_d(s[1]));
      break;
    }
    case(10):  //----------------------------------------------
    {/*:Longitude*/
      if (Options.noisy) {cout <<"Gauge longitude"<<endl;}
      ExitGracefullyIf(pGage==NULL,
                       "ParseTimeSeriesFile::Longitude specified outside of a :Gauge-:EndGauge statement",BAD_DATA);
      ExitGracefullyIf(Len<=1,"ParseTimeSeriesFile - no value provided after :Longitude command",BAD_DATA_WARN);
      pGage->SetLongitude(s_to_d(s[1]));
      break;
    }
    case(11):  //----------------------------------------------
    {/*:Elevation*/
      if (Options.noisy) {cout <<"Gauge elevation"<<endl;}
      ExitGracefullyIf(pGage==NULL,
                       "ParseTimeSeriesFile::Elevation specified outside of a :Gauge-:EndGauge statement",BAD_DATA);
      ExitGracefullyIf(Len<=1,"ParseTimeSeriesFile - no value provided after :Elevation command",BAD_DATA_WARN);
      pGage->SetElevation(s_to_d(s[1]));
      break;
    }
    case(12):  //----------------------------------------------
    {/*:Data [DATA TYPE] {units string}
       {yyyy-mm-dd hh:mm:ss double tstep int nMeasurements}
       {double values} x nMeasurements
       :EndData

       or

       :Data [DATA TYPE] {units string}
          :ReadFromNetCDF
             :FileNameNC     {NetCDF file name}
             :VarNameNC      {name of variable in NetCDF file; variable must be 2D (nstations x ntime) or (ntime x nstations)}
             :DimNamesNC     {dimension name of stations ; dimension name of time}
             #:RedirectToFile {weights of station contribution to each HRU; usually redirect to txt-file}  --> is coming from Interp_from_file (*.rvi)
          :EndReadFromNetCDF
       :EndData
     */
      if(Options.noisy) { cout <<"Time Series Data: "<<s[1]<<endl; }
      string name=to_string(s[1]);
      ExitGracefullyIf(pGage==NULL,
        "ParseTimeSeriesFile::Time Series Data added outside of a :Gage-:EndGauge statement",BAD_DATA);
      pTimeSer=CTimeSeries::Parse(p,true,name,DOESNT_EXIST,pGage->GetName(),Options);
      forcing_type ftype=GetForcingTypeFromString(name);
      if(ftype==F_UNRECOGNIZED) { ExitGracefully("Unrecognized forcing type string in :Data command",BAD_DATA_WARN); }
      if(pGage==NULL){ExitGracefully("ParseTimeSeriesFile:: :Data command must be within a :Gauge-:EndGauge block.",BAD_DATA_WARN); }
      else {
        pGage->AddTimeSeries(pTimeSer,ftype);
      }
      break;
    }
    case(13):  //----------------------------------------------
    {/*:MultiData
       {yyyy-mm-dd hh:mm:ss double tstep int nMeasurements}
       :Parameters, PARAM_TAG_1, PARAM_TAG_2, ...
       :Units, unit1, unit2,...
       {double p1, double p2, ... } x nMeasurements [mm/d]
       :EndMultiData
     */
      if (Options.noisy) {cout <<"Multiple Time Series"<<endl;}
      ExitGracefullyIf(pGage==NULL,
                       "ParseTimeSeriesFile:: :MultiData command added before specifying a gauge station and its properties",BAD_DATA);

      int nSeries=0;
      CTimeSeries **pTimeSerArray=NULL;
      forcing_type       aTypes[MAX_MULTIDATA];
      for (int i=0;i<MAX_MULTIDATA;i++){aTypes[i]=F_UNRECOGNIZED;}

      pTimeSerArray=CTimeSeries::ParseMultiple(p,nSeries,aTypes,true,Options);

      ExitGracefullyIf(nSeries<=0,"Parse :Multidata command - no series specified",BAD_DATA);
      if(pGage==NULL) { ExitGracefully("ParseTimeSeriesFile:: :MultiData command must be within a :Gauge-:EndGauge block.",BAD_DATA_WARN); }
      else {
        for(int i=0; i<nSeries; i++)
        {
          pGage->AddTimeSeries(pTimeSerArray[i],aTypes[i]);
        }
      }
      break;
    }
    case(16):  //----------------------------------------------
    {/*:MeasurementHeight [value] */
      if (Options.noisy) {cout <<"Gauge measurement height (above land surface)"<<endl;}
      ExitGracefullyIf(pGage==NULL,
                       "ParseTimeSeriesFile::MeasurementHeight specified outside of a :Gauge-:EndGauge statement",BAD_DATA);
      pGage->SetMeasurementHt(s_to_d(s[1]));
      break;
    }
    case(20):  //----------------------------------------------
    {/*:EnsimTimeSeries*/
      ExitGracefullyIf(pGage==NULL,
                       "ParseTimeSeriesFile:: ensim time series file specified before specifying a gauge station and its properties",BAD_DATA);

      int nSeries=0;
      CTimeSeries **pTimeSerArray=NULL;
      forcing_type       aTypes[MAX_MULTIDATA];
      for (int i=0;i<MAX_MULTIDATA;i++){aTypes[i]=F_UNRECOGNIZED;}

      string filename;
      for (int i=1;i<Len;i++){filename+=s[i]; if(i<Len-1){filename+=' ';}}
      if (Options.noisy) {cout <<"Ensim Time Series: "<<filename<<endl;}

      filename =CorrectForRelativePath(filename ,Options.rvt_filename);

      pTimeSerArray=CTimeSeries::ParseEnsimTb0(filename,nSeries,aTypes,Options);

      ExitGracefullyIf(nSeries<=0,"Parse :EnsimTimeSeries command - no series specified- likely incorrect .tb0 format",BAD_DATA);
      for (int i=0; i<nSeries; i++)
      {
        pGage->AddTimeSeries(pTimeSerArray[i],aTypes[i]);
      }
      break;
    }
    case(30):  //----------------------------------------------
    {/* :RainCorrection [double value] */
      if (Options.noisy) {cout <<"Rainfall Correction"<<endl;}
      ExitGracefullyIf(pGage==NULL && pGrid==NULL,
                       "ParseTimeSeriesFile::Precipitation correction added before specifying a gauge station/ gridded forcing and its properties",BAD_DATA);
      if (Len<2){p->ImproperFormat(s); break;}
      if (pGage != NULL) { pGage->SetGaugeProperty("RAINFALL_CORR",s_to_d(s[1])); }
      if (pGrid != NULL) { pGrid->SetRainfallCorr(s_to_d(s[1])); }
      break;
    }
    case(31):  //----------------------------------------------
    {/* :SnowCorrection [double value] */
      if (Options.noisy) {cout <<"Snowfall Correction"<<endl;}
      ExitGracefullyIf(pGage==NULL && pGrid==NULL,
                       "ParseTimeSeriesFile::Snowfall correction added before specifying a gauge station/ gridded forcing and its properties",BAD_DATA);
      if (Len<2){p->ImproperFormat(s); break;}
      if (pGage != NULL) { pGage->SetGaugeProperty("SNOWFALL_CORR",s_to_d(s[1])); }
      if (pGrid != NULL) { pGrid->SetSnowfallCorr(s_to_d(s[1])); }
      break;
    }
    case(32):  //----------------------------------------------
    {/* :CloudTempRanges [double cloud_temp_min C] [cloud_temp_max C] */
      if (Options.noisy) {cout <<"Cloud cover temperature ranges"<<endl;}
      ExitGracefullyIf(pGage==NULL,
                       "ParseTimeSeriesFile::Cloudcover temperature ranges added before specifying a gauge station and its properties",BAD_DATA);
      if (Len<3){p->ImproperFormat(s); break;}
      pGage->SetGaugeProperty("CLOUD_MIN_RANGE",s_to_d(s[1]));
      pGage->SetGaugeProperty("CLOUD_MAX_RANGE",s_to_d(s[2]));
      break;
    }
    case(33):  //----------------------------------------------
    {/* :TemperatureCorrection [double value] */
        if (Options.noisy) { cout << "Temperature Correction" << endl; }
        ExitGracefullyIf(pGage == NULL && pGrid == NULL,
            "ParseTimeSeriesFile::Temperature correction added before specifying a gauge station/ gridded forcing and its properties", BAD_DATA);
        if (Len < 2) { p->ImproperFormat(s); break; }
        if (pGage != NULL) { pGage->SetGaugeProperty("TEMP_CORR", s_to_d(s[1])); }
        if (pGrid != NULL) { pGrid->SetTemperatureCorr(s_to_d(s[1])); }
        break;
    }
    case (40): //---------------------------------------------
    {/*:ObservationData [data type] [long SBID or int HRUID] {constituent name if data type=STREAM_CONCENTRATION }
       {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
       {double value} x nMeasurements
       :EndObservationData
     */
      if (Options.noisy) {cout <<"Observation data"<<endl;}
      if (Len<3){p->ImproperFormat(s); break;}

      bool ishyd      =!strcmp(s[1], "HYDROGRAPH");
      bool isstage    =!strcmp(s[1], "RESERVOIR_STAGE");
      bool isinflow   =!strcmp(s[1], "RESERVOIR_INFLOW");
      bool isnetinflow=!strcmp(s[1], "RESERVOIR_NETINFLOW");
      bool isconc     =!strcmp(s[1], "STREAM_CONCENTRATION");
      bool istemp     =!strcmp(s[1], "STREAM_TEMPERATURE");
      bool islevel    =!strcmp(s[1], "WATER_LEVEL");
      bool islakearea =!strcmp(s[1], "LAKE_AREA");
      bool invalidSB=(pModel->GetSubBasinByID(s_to_l(s[2]))==NULL);

      bool period_ending =ishyd;
      //Hydrographs are internally stored as period-ending!

      pTimeSer=CTimeSeries::Parse(p,true,to_string(s[1]),s_to_l(s[2]),"none",Options,period_ending);

      if(isconc) {
        ExitGracefullyIf(Len<4,"ParseTimeSeriesFile: STREAM_CONCENTRATION observation must include constituent name",BAD_DATA_WARN);
        int c = pModel->GetTransportModel()->GetConstituentIndex(s[3]);
        if(c==DOESNT_EXIST) {
          warn="ParseTimeSeries:: Invalid/unused constituent name in observation stream concentration time series ["+pTimeSer->GetSourceFile()+"]. Will be ignored";
          WriteWarning(warn.c_str(),Options.noisy); break;
        }
        pTimeSer->SetConstitInd(c);
      }
      else if(istemp) {
        int c = pModel->GetTransportModel()->GetConstituentIndex("TEMPERATURE");
        pTimeSer->SetConstitInd(c);
      }

      if (ishyd && invalidSB){
        warn="ParseTimeSeries:: Invalid subbasin ID in observation hydrograph time series ["+pTimeSer->GetSourceFile()+"]. Will be ignored";
        WriteWarning(warn.c_str(),Options.noisy); break;
      }
      if(isstage && invalidSB){
        warn="ParseTimeSeries:: Invalid subbasin ID in observed reservoir stage time series ["+pTimeSer->GetSourceFile()+"]. Will be ignored";
        WriteWarning(warn.c_str(),Options.noisy); break;
      }
      if(isinflow && invalidSB){
        warn="ParseTimeSeries:: Invalid subbasin ID in observed reservoir inflow time series ["+pTimeSer->GetSourceFile()+"]. Will be ignored";
        WriteWarning(warn.c_str(),Options.noisy); break;
      }
      if(isnetinflow && invalidSB){
        warn="ParseTimeSeries:: Invalid subbasin ID in observed reservoir net inflow time series ["+pTimeSer->GetSourceFile()+"]. Will be ignored";
        WriteWarning(warn.c_str(),Options.noisy);  break;
      }
      if(islevel && invalidSB){
        warn="ParseTimeSeries:: Invalid subbasin ID in observed water level time series ["+pTimeSer->GetSourceFile()+"]. Will be ignored";
        WriteWarning(warn.c_str(),Options.noisy);  break;
      }
      if(islakearea && invalidSB){
        warn="ParseTimeSeries:: Invalid subbasin ID in observed lake area time series ["+pTimeSer->GetSourceFile()+"]. Will be ignored";
        WriteWarning(warn.c_str(),Options.noisy); break;
      }
      pModel->AddObservedTimeSeries(pTimeSer);
      break;
    }
    case (41): //---------------------------------------------
    {/*:IrregularObservations {data type} {long SBID or int HRUID} {int nMeasurements}
       {yyyy-mm-dd} {hh:mm:ss.0} {double value} x nMeasurements
       :EndIrregularObservations
     */
      if (Options.noisy) {cout <<"Irregular Observation"<<endl;}
      if (Len<4){p->ImproperFormat(s); break;}
      CTimeSeriesABC *pIrregTimeSers;
      pIrregTimeSers=CIrregularTimeSeries::Parse(p,to_string(s[1]),s_to_l(s[2]),s_to_i(s[3]));
      pModel->AddObservedTimeSeries(pIrregTimeSers);
      break;
    }
    case (42): //---------------------------------------------
    {/*:ObservationWeights {data type} {long SBID or int HRUID}  {constituent name if data type=STREAM_CONCENTRATION }
       {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
       {double value} x nMeasurements
       :EndObservationWeights
     */
      if (Options.noisy) {cout <<"Observation weights"<<endl;}
      if (Len<3){p->ImproperFormat(s); break;}

      bool ishyd      =!strcmp(s[1], "HYDROGRAPH");
      bool isstage    =!strcmp(s[1], "RESERVOIR_STAGE");
      bool isinflow   =!strcmp(s[1], "RESERVOIR_INFLOW");
      bool isnetinflow=!strcmp(s[1], "RESERVOIR_NETINFLOW");
      bool isconc     =!strcmp(s[1], "STREAM_CONCENTRATION");
      bool istemp     =!strcmp(s[1], "STREAM_TEMPERATURE");
      bool islevel    =!strcmp(s[1], "WATER_LEVEL");
      bool islakearea =!strcmp(s[1], "LAKE_AREA");
      bool invalidSB=(pModel->GetSubBasinByID(s_to_l(s[2]))==NULL);

      pTimeSer=CTimeSeries::Parse(p,true,to_string(s[1]),s_to_l(s[2]),"none",Options);

      if(isconc) {
        ExitGracefullyIf(Len<4,"ParseTimeSeriesFile: STREAM_CONCENTRATION observation must include constituent name",BAD_DATA_WARN);
        int c = pModel->GetTransportModel()->GetConstituentIndex(s[3]);
        if(c==DOESNT_EXIST) {
          warn="ParseTimeSeries:: Invalid/unused constituent name in observation stream concentration time series ["+pTimeSer->GetSourceFile()+"]. Will be ignored";
          WriteWarning(warn.c_str(),Options.noisy); break;
        }
        pTimeSer->SetConstitInd(c);
      }
      else if(istemp) {
        int c = pModel->GetTransportModel()->GetConstituentIndex("TEMPERATURE");
        pTimeSer->SetConstitInd(c);
      }

      if (ishyd && invalidSB){
        warn="ParseTimeSeries:: Invalid subbasin ID in observation hydrograph weights time series ["+pTimeSer->GetSourceFile()+"]. Will be ignored";
        WriteWarning(warn.c_str(),Options.noisy); break;
      }
      if(isstage && invalidSB){
        warn="ParseTimeSeries:: Invalid subbasin ID in observed reservoir stage weights time series ["+pTimeSer->GetSourceFile()+"]. Will be ignored";
        WriteWarning(warn.c_str(),Options.noisy); break;
      }
      if(isinflow && invalidSB){
        warn="ParseTimeSeries:: Invalid subbasin ID in observed reservoir inflow weights time series ["+pTimeSer->GetSourceFile()+"]. Will be ignored";
        WriteWarning(warn.c_str(),Options.noisy); break;
      }
      if(isnetinflow && invalidSB){
        warn="ParseTimeSeries:: Invalid subbasin ID in observed reservoir net inflow weights time series ["+pTimeSer->GetSourceFile()+"]. Will be ignored";
        WriteWarning(warn.c_str(),Options.noisy);  break;
      }
      if(islevel && invalidSB){
        warn="ParseTimeSeries:: Invalid subbasin ID in observed water level weights time series ["+pTimeSer->GetSourceFile()+"]. Will be ignored";
        WriteWarning(warn.c_str(),Options.noisy);  break;
      }
      if(islakearea && invalidSB){
        warn="ParseTimeSeries:: Invalid subbasin ID in observed lake area weights time series ["+pTimeSer->GetSourceFile()+"]. Will be ignored";
        WriteWarning(warn.c_str(),Options.noisy);  break;
      }
      pModel->AddObservedWeightsTS(pTimeSer);
      break;
    }
    case (43): //---------------------------------------------
    {/*:IrregularWeights {data type} {long SBID or int HRUID} {int nMeasurements}
       {yyyy-mm-dd} {hh:mm:ss.0} {double weight} x nMeasurements
       :EndIrregularWeights
     */
      if (Options.noisy) {cout <<"Irregular weights"<<endl;}
      if (Len<4){p->ImproperFormat(s); break;}
      CTimeSeriesABC *pIrregTimeSers;
      pIrregTimeSers=CIrregularTimeSeries::Parse(p,to_string(s[1]),s_to_l(s[2]),s_to_i(s[3]));
      pModel->AddObservedWeightsTS(pIrregTimeSers);
      break;
    }
    case (50): //---------------------------------------------
    {/*:BasinInflowHydrograph {long SBID}
       {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
       {double Qin} x nMeasurements [m3/s]
       :EndBasinInflowHydrograph
     */
      if (Options.noisy) {cout <<"Basin Inflow Hydrograph"<<endl;}
      long SBID=DOESNT_EXIST;
      CSubBasin *pSB;
      if (Len>=2){SBID=s_to_l(s[1]);}

      pSB=pModel->GetSubBasinByID(SBID);
      pTimeSer=CTimeSeries::Parse(p,false,"Inflow_Hydrograph_"+to_string(SBID),SBID,"none",Options);
      if (pSB!=NULL){
        pSB->AddInflowHydrograph(pTimeSer);
      }
      else
      {
        warn=":BasinInflowHydrograph Subbasin "+to_string(SBID)+" not in model, cannot set inflow hydrograph";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case (51): //---------------------------------------------
    {/*:BasinInflowHydrograph2 {long Basincode}
       {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
       {double Qin} x nMeasurements [m3/s]
       :EndBasinInflowHydrograph2
     */
      if (Options.noisy) {cout <<"Basin Inflow Hydrograph(2)"<<endl;}
      long SBID=DOESNT_EXIST;
      CSubBasin *pSB;
      if (Len>=2){SBID=s_to_l(s[1]);}
      pSB=pModel->GetSubBasinByID(SBID);
      pTimeSer=CTimeSeries::Parse(p,false,"Inflow_Hydrograph_"+to_string(SBID),SBID,"none",Options);
      if (pSB!=NULL){
        pSB->AddDownstreamInflow(pTimeSer);
      }
      else
      {
        warn=":BasinInflowHydrograph2 Subbasin "+to_string(SBID)+" not in model, cannot set inflow hydrograph";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case (52): //---------------------------------------------
    {/*:ReservoirDemand {long SBID} [int demandID] [string demandName] (or :ReservoirExtraction)
         {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
         {double Qout} x nMeasurements [m3/s]
       :EndReservoirDemand
     */
      if (Options.noisy) {cout <<"Reservoir Water Demand Time Series"<<endl;}
      long SBID=DOESNT_EXIST;
      int  demand_ID=DOESNT_EXIST;
      string demand_name = "";

      CSubBasin *pSB;
      if (Len>=2){SBID=s_to_l(s[1]);}
      if (Len>=3){demand_ID=s_to_i(s[2]);}
      if (Len>=4){demand_name=s[3];}

      pSB=pModel->GetSubBasinByID(SBID);

      if ((pSB!=NULL) && (pSB->GetReservoir()!=NULL))
      {
	    has_irrig=true;

        pTimeSer=CTimeSeries::Parse(p,true,demand_name,SBID,"none",Options);
        pTimeSer->SetIDTag(demand_ID);

        int ii=pSB->GetReservoir()->GetNumWaterDemands();

        CDemand *pDem=new CDemand(demand_ID,demand_name,SBID,true,pTimeSer);
        pDem->SetLocalIndex(ii);
        pSB->GetReservoir()->AddDemand(pDem);
        if (Options.management_optimization){
          pModel->GetDemandOptimizer()->AddWaterDemand(pDem);
        }
      }
      else
      {
        warn=":ReservoirDemand Subbasin "+to_string(SBID)+" not in model or doesn't have reservoir, cannot set Reservoir Demand";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case (53): //---------------------------------------------
    {/*:VariableWeirHeight {long SBID}
       {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
       {double delta_h} x nMeasurements [m]
       :EndVariableWeirHeight
     */
      if (Options.noisy) {cout <<"Weir Height Time Series"<<endl;}
      long SBID=DOESNT_EXIST;
      CSubBasin *pSB;
      if (Len>=2){SBID=s_to_l(s[1]);}
      pSB=pModel->GetSubBasinByID(SBID);
      pTimeSer=CTimeSeries::Parse(p,true,"WeirHeight_"+to_string(SBID),SBID,"none",Options);
      if ((pSB!=NULL) && (pSB->GetReservoir()!=NULL)){
        pSB->GetReservoir()->AddWeirHeightTS(pTimeSer);
      }
      else
      {
        warn=":VariableWeirHeight Subbasin "+to_string(SBID)+" not in model or doesn't have reservoir, cannot set Reservoir weir height";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case (54): //---------------------------------------------
    {/*:ReservoirMaxStage {long SBID}
       {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
       {double stage} x nMeasurements [m]
       :EndReservoirMaxStage
     */
      if (Options.noisy) {cout <<"Maximum stage time series "<<endl;}
      long SBID=DOESNT_EXIST;
      CSubBasin *pSB;
      if (Len>=2){SBID=s_to_l(s[1]);}
      pSB=pModel->GetSubBasinByID(SBID);
      pTimeSer=CTimeSeries::Parse(p,true,"_MaxStage_"+to_string(SBID),SBID,"none",Options);
      if ((pSB!=NULL) && (pSB->GetReservoir()!=NULL)){
        if (pModel->GetDemandOptimizer() != NULL) {
          pModel->GetDemandOptimizer()->AddUserTimeSeries(pTimeSer);
        }
        else {
          pSB->GetReservoir()->AddMaxStageTimeSeries(pTimeSer);
        }
      }
      else
      {
        warn=":ReservoirMaxStage Subbasin "+to_string(SBID)+" not in model or doesn't have reservoir, cannot set Reservoir maximum stage";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case (55): //---------------------------------------------
    {/*:ReservoirMinStage {long SBID}
       {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
       {double stage} x nMeasurements [m]
       :EndReservoirMinStage
     */
      if (Options.noisy) {cout <<"Reservoir Minimum Stage Time Series"<<endl;}
      long SBID=DOESNT_EXIST;
      CSubBasin *pSB;
      if (Len>=2){SBID=s_to_l(s[1]);}
      pSB=pModel->GetSubBasinByID(SBID);
      pTimeSer=CTimeSeries::Parse(p,true,"_MinStage_"+to_string(SBID),SBID,"none",Options);
      if ((pSB!=NULL) && (pSB->GetReservoir()!=NULL)){
        if (pModel->GetDemandOptimizer() != NULL) {
          pModel->GetDemandOptimizer()->AddUserTimeSeries(pTimeSer);
        }
        else {
          pSB->GetReservoir()->AddMinStageTimeSeries(pTimeSer);
        }
      }
      else
      {
        warn=":ReservoirMinStage Subbasin "+to_string(SBID)+" not in model or doesn't have reservoir, cannot set minimum stage";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case (56): //---------------------------------------------
    {/*:ReservoirMinStageFlow {long SBID}
       {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
       {double Q_min} x nMeasurements [m3/s]
       :EndReservoirMinStageFlow
     */
      if (Options.noisy) {cout <<"Reservoir Minimum Stage Discharge Time Series"<<endl;}
      long SBID=DOESNT_EXIST;
      CSubBasin *pSB;
      if (Len>=2){SBID=s_to_l(s[1]);}
      pSB=pModel->GetSubBasinByID(SBID);
      pTimeSer=CTimeSeries::Parse(p,true,"_MinStageFlow_"+to_string(SBID),SBID,"none",Options);
      if ((pSB!=NULL) && (pSB->GetReservoir()!=NULL)){
        if (pModel->GetDemandOptimizer() != NULL) {
          pModel->GetDemandOptimizer()->AddUserTimeSeries(pTimeSer);
        }
        else {
          pSB->GetReservoir()->AddMinStageFlowTimeSeries(pTimeSer);
        }
      }
      else
      {
        warn=":ReservoirMinStageFlow Subbasin "+to_string(SBID)+" not in model or doesn't have reservoir, cannot set minimum stage discharge";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case (57): //---------------------------------------------
    {/*:ReservoirTargetStage {long SBID}
       {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
       {double h_target} x nMeasurements [m]
       :EndReservoirTargetStage
     */
      if (Options.noisy) {cout <<"Reservoir Target Stage Time Series"<<endl;}
      long SBID=DOESNT_EXIST;
      CSubBasin *pSB;
      if (Len>=2){SBID=s_to_l(s[1]);}
      pSB=pModel->GetSubBasinByID(SBID);
      pTimeSer=CTimeSeries::Parse(p,true,"_TargetStage_"+to_string(SBID),SBID,"none",Options);
      if((pSB!=NULL) && (pSB->GetReservoir()!=NULL)) {
        if (pModel->GetDemandOptimizer() != NULL) {
          pModel->GetDemandOptimizer()->AddUserTimeSeries(pTimeSer);
        }
        else {
          pSB->GetReservoir()->AddTargetStageTimeSeries(pTimeSer);
        }
      }
      else
      {
        warn=":ReservoirTargetStage Subbasin "+to_string(SBID)+" not in model or doesnt have reservoir, cannot set target stage";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case (58): //---------------------------------------------
    {/*:ReservoirMaxQDelta {long SBID}
       {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
       {double Qdelta} x nMeasurements [m3/s/d]
       :EndReservoirMaxQDelta
     */
      if (Options.noisy) {cout <<"Reservoir Maximum Discharge Delta Time Series"<<endl;}
      long SBID=DOESNT_EXIST;
      CSubBasin *pSB;
      if (Len>=2){SBID=s_to_l(s[1]);}
      pSB=pModel->GetSubBasinByID(SBID);
      pTimeSer=CTimeSeries::Parse(p,true,"_MaxQDelta_"+to_string(SBID),SBID,"none",Options);
      if ((pSB!=NULL) && (pSB->GetReservoir()!=NULL)){
        if (pModel->GetDemandOptimizer() != NULL) {
          pModel->GetDemandOptimizer()->AddUserTimeSeries(pTimeSer);
        }
        else {
          pSB->GetReservoir()->AddMaxQIncreaseTimeSeries(pTimeSer);
        }
      }
      else
      {
        warn=":ReservoirMaxQDelta Subbasin "+to_string(SBID)+" not in model or doesn't have reservoir, cannot set maximum Qdelta";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case (59): //---------------------------------------------
    {/*:ReservoirMaxQDecrease {long SBID}
     {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
     {double Qdelta} x nMeasurements [m3/s/d]
     :EndReservoirMaxQDecrease
     */
      if(Options.noisy) { cout <<"Reservoir Maximum Discharge Decrease Time Series"<<endl; }
      long SBID=DOESNT_EXIST;
      CSubBasin *pSB;
      if(Len>=2) { SBID=s_to_l(s[1]); }
      pSB=pModel->GetSubBasinByID(SBID);
      if((pSB!=NULL) && (pSB->GetReservoir()!=NULL)) {
        pTimeSer=CTimeSeries::Parse(p,true,"_MaxQDecrease_"+to_string(SBID),SBID,"none",Options);
        if (pModel->GetDemandOptimizer() != NULL) {
          pModel->GetDemandOptimizer()->AddUserTimeSeries(pTimeSer);
        }
        else {
          pSB->GetReservoir()->AddMaxQDecreaseTimeSeries(pTimeSer);
        }
      }
      else
      {
        warn=":ReservoirMaxQDecrease Subbasin "+to_string(SBID)+" not in model or doesn't have reservoir, cannot set maximum Qdelta decrease";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case (60): //---------------------------------------------
    {/*:ReservoirMinFlow {long Basincode}
     {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
     {double Qmin} x nMeasurements [m3/s]
     :EndReservoirMinFlow
     */
      if(Options.noisy) { cout <<"Reservoir Minimum Flow time series"<<endl; }
      long SBID=DOESNT_EXIST;
      CSubBasin *pSB;
      if(Len>=2) { SBID=s_to_l(s[1]); }
      pSB=pModel->GetSubBasinByID(SBID);
      if((pSB!=NULL) && (pSB->GetReservoir()!=NULL)){
        pTimeSer=CTimeSeries::Parse(p,true,"_ResQmin_"+to_string(SBID),SBID,"none",Options);
        if (pModel->GetDemandOptimizer() != NULL) {
          pModel->GetDemandOptimizer()->AddUserTimeSeries(pTimeSer);
        }
        else {
          pSB->GetReservoir()->AddMinQTimeSeries(pTimeSer);
        }
      }
      else
      {
        warn=":ReservoirMinFlow Subbasin "+to_string(SBID)+" not in model or doesn't have reservoir, cannot set minimum flow";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case (61): //---------------------------------------------
    {/*:ReservoirMaxFlow {long Basincode}
     {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
     {double Qmin} x nMeasurements [m3/s]
     :EndReservoirMinFlow
     */
      if(Options.noisy) { cout <<"Reservoir Maximum Flow time series"<<endl; }
      long SBID=DOESNT_EXIST;
      CSubBasin *pSB;
      if(Len>=2) { SBID=s_to_l(s[1]); }
      pSB=pModel->GetSubBasinByID(SBID);
      if((pSB!=NULL) && (pSB->GetReservoir()!=NULL)) {
        pTimeSer=CTimeSeries::Parse(p,true,"_ResQmax_"+to_string(SBID),SBID,"none",Options);
        if (pModel->GetDemandOptimizer() != NULL) {
          pModel->GetDemandOptimizer()->AddUserTimeSeries(pTimeSer);
        }
        else{
          pSB->GetReservoir()->AddMaxQTimeSeries(pTimeSer);
        }
      }
      else
      {
        warn=":ReservoirMaxFlow Subbasin "+to_string(SBID)+" not in model or doesn't have reservoir, cannot set maximum flow";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case (62): //---------------------------------------------
    {/*:OverrideReservoirFlow {long SBID}
       {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
       {double Q} x nMeasurements [m3/s]
       :EndOverrideReservoirFlow
     */
      if (Options.noisy) {cout <<"Forced Reservoir Outflow Time Series"<<endl;}
      long SBID=DOESNT_EXIST;
      CSubBasin *pSB;
      if (Len>=2){SBID=s_to_l(s[1]);}
      pSB=pModel->GetSubBasinByID(SBID);
      pTimeSer=CTimeSeries::Parse(p,true,"_ResFlow_"+to_string(SBID),SBID,"none",Options);
      if ((pSB!=NULL) && (pSB->GetReservoir()!=NULL)){
        if (pModel->GetDemandOptimizer() != NULL) {
          pModel->GetDemandOptimizer()->AddUserTimeSeries(pTimeSer);
        }
        else{
          pSB->GetReservoir()->AddOverrideQTimeSeries(pTimeSer);
        }
      }
      else
      {
        warn=":OverrideReservoirFlow Subbasin "+to_string(SBID)+" not in model or doesn't have reservoir, cannot set overriden discharge";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case (63): //---------------------------------------------
    {/*:IrrigationDemand or :WaterDemand {long SBID} [int DemandID] [string demand_name/alias]
       {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
       {double Qin} x nMeasurements [m3/s]
     :EndIrrigationDemand or :EndWaterDemand
     */
      if(Options.noisy) { cout <<"Irrigation/Water use demand time series "<<endl; }
      long SBID     =DOESNT_EXIST;
      int  demand_ID=DOESNT_EXIST;
      int  ii       =0;
      string demand_name = "";

      if(Len>=2) {SBID       =s_to_l(s[1]); }
      if(Len>=3) {demand_ID  =s_to_i(s[2]);}
      if(Len>=4) {demand_name=s[3];}

      CSubBasin *pSB=pModel->GetSubBasinByID(SBID);
      if (pSB!=NULL) {ii=pSB->GetNumWaterDemands();}
      if (demand_ID   == DOESNT_EXIST) { demand_ID = SBID*10+ii;} //DEFAULT - SBIDx
      if (demand_name == "")           { demand_name = "D"+to_string(SBID)+"abcdefghijklmnopqrstuvwxyz"[ii];} //e.g., D120b; }

      pTimeSer = CTimeSeries::Parse(p, false, demand_name, SBID, "none", Options);
      pTimeSer->SetIDTag(demand_ID);

      if(pSB!=NULL) {
        has_irrig=true;

        CDemand *pDem=new CDemand(demand_ID,demand_name,SBID,false,pTimeSer);
        pDem->SetLocalIndex(ii);
        pSB->AddWaterDemand(pDem);
        if (Options.management_optimization){
          pModel->GetDemandOptimizer()->AddWaterDemand(pDem);
        }
      }
      else
      {
        warn=":IrrigationDemand: Subbasin "+to_string(SBID)+" not in model, cannot set irrigation/water demand time series";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case (64): //---------------------------------------------
    {/*:ReservoirDownstreamFlow {long Basincode} {long downstreamSBID} {double range}
     {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
     {double Qtarget} x nMeasurements [m3/s]
     :EndReservoirDownstreamFlow
     */
      if(Options.noisy) { cout <<"Reservoir Dowstream Flow  target time series"<<endl; }
      long SBID=DOESNT_EXIST;
      long SBID_down=DOESNT_EXIST;
      CSubBasin *pSB,*pSBdown;
      double range=0;
      if(Len>=4) {
        SBID=s_to_l(s[1]);
        SBID_down=s_to_l(s[2]);
        range=s_to_d(s[3]);
      }
      else {
        ExitGracefully(":ReservoirDownstreamFlow: incorrect number of terms in command",BAD_DATA);
      }
      pSB    =pModel->GetSubBasinByID(SBID);
      pSBdown=pModel->GetSubBasinByID(SBID_down);
      pTimeSer=CTimeSeries::Parse(p,true,"_ResDownQ_"+to_string(SBID),SBID,"none",Options);
      if((pSB!=NULL) && (pSBdown!=NULL) && (pSB->GetReservoir()!=NULL)){
        pSB->GetReservoir()->AddDownstreamTargetQ(pTimeSer,pSBdown,range);
      }
      else
      {
        warn=":ReservoirDownstreamFlow Subbasin "+to_string(SBID)+" or downstream subbasin "+to_string(SBID_down)+" not in model or doesn't have reservoir, cannot set downstream flow target";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case (65): //---------------------------------------------
    {/*:ReservoirDownstreamDemand {long downstream SBID} {long reservoir SBID} {double percent_met} [julian_start] [julian_end]*/
      if(Options.noisy) { cout <<"Reservoir downstream demand"<<endl; }
      long   SBID     =DOESNT_EXIST;
      long   SBIDres  =DOESNT_EXIST;
      int    jul_start=0;
      int    jul_end  =365;
      double pct_met  =0.0;
      if(Len>=4) {
        SBID=s_to_l(s[1]);
        double tmp=AutoOrDouble(s[2]);
        if(tmp==AUTO_COMPUTE) { SBIDres=AUTO_COMPUTE_LONG; }
        else                  { SBIDres=s_to_l(s[2]); }
        pct_met=s_to_d(s[3]);
      }
      if(Len>=6) { //optional date commands
        jul_start=s_to_i(s[4])-1;
        jul_end  =s_to_i(s[5])-1; //converts to internal Raven convention
      }

      AllocateReservoirDemand(pModel,Options,SBID,SBIDres,pct_met,jul_start,jul_end);

      break;
    }

    case (66): //---------------------------------------------
    {/*:FlowDiversion [fromSBID] [toSBID] [fract. diverted] [Qmin] {start_day} {end_day}*/
      if(Options.noisy) { cout <<"Flow diversion"<<endl; }
      long SBID=DOESNT_EXIST;
      CSubBasin *pSB;
      if(Len>=2) { SBID=s_to_l(s[1]); }
      pSB=pModel->GetSubBasinByID(SBID);

      if(pSB!=NULL) {
        int start=0; //default - all year
        int end  =366; //default - all year
        if(Len>=7) { start=s_to_i(s[5]); end=s_to_i(s[6]); }
        int target_p=pModel->GetSubBasinIndex(s_to_l(s[2]));
        if ((s_to_l(s[2])==-1) || (target_p!=DOESNT_EXIST)) {
          pSB->AddFlowDiversion(start,end,target_p,s_to_d(s[4]),s_to_d(s[3])); has_irrig=true;
        }
        else {
          warn=":FlowDiversion command: Target subbasin "+to_string(s_to_l(s[2]))+" not found in model, cannot add diversion";
          WriteWarning(warn,Options.noisy);
        }
      }
      else
      {
        warn=":FlowDiversion command: Subbasin "+to_string(SBID)+" not in model, cannot add diversion";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case (67): //---------------------------------------------
    {/* :FlowDiversionLookupTable [fromSBID] [toSBID] {start_day} {end_day}
     nPoints
     {Qsource_i Qdivert_i} x nPoints
     :EndFlowDiversionLookupTable
     */
      if(Options.noisy) { cout << ":FlowDiversionRatingCurve" << endl; }
      long SBID=DOESNT_EXIST;
      long SBID2=0;
      int target_p(DOESNT_EXIST),start, end;
      int NQ(0);
      CSubBasin *pSB;
      start=0;   //default - all year
      end  =366; //default - all year

      if(Len>=2) { SBID=s_to_l(s[1]); }
      pSB=pModel->GetSubBasinByID(SBID);

      if(pSB!=NULL) {
        if(Len>=5) { start=s_to_i(s[3]); end=s_to_i(s[4]); }
        target_p=pModel->GetSubBasinIndex(s_to_l(s[2]));
        SBID2=s_to_l(s[2]);
      }
      else if (SBID == DOESNT_EXIST) {
        warn=":FlowDiversionLookupTable command: no source subbasin ID provided, cannot add diversion";
        WriteWarning(warn,Options.noisy);
      }
      else
      {
        warn=":FlowDiversionLookupTable command: Subbasin "+to_string(SBID)+" not found in model, cannot add diversion";
        WriteWarning(warn,Options.noisy);
      }

      p->Tokenize(s,Len);
      if(Len >= 1) { NQ = s_to_i(s[0]); }

      double *aQ1 = new double[NQ];
      double *aQ2 = new double[NQ];
      for(int i = 0; i < NQ; i++) {
        p->Tokenize(s,Len);
        if(IsComment(s[0],Len)) { i--; }
        else {
          if(Len>=2) {
            aQ1[i] = s_to_d(s[0]);  aQ2[i] = s_to_d(s[1]);
          }
          else {
            WriteWarning("Incorrect line length (<2) in :FlowDiversionLookupTable command",Options.noisy);
          }
        }
      }
      p->Tokenize(s,Len); //:EndFlowDiversionRatingCurve

      if((SBID2==-1) || (target_p!=DOESNT_EXIST)) {
        pSB->AddFlowDiversion(start,end,target_p,aQ1,aQ2,NQ);
      }
      else {
        warn=":FlowDiversionLookupTable command: Target subbasin "+to_string(s_to_l(s[2]))+" not in model, cannot add diversion";
        WriteWarning(warn,Options.noisy);
      }
      delete [] aQ1; delete [] aQ2;
      break;
    }
    case (68): //---------------------------------------------
    {/*:EnvironmentalMinFlow {long SBID}
     {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
     {double Qmin} x nMeasurements [m3/s]
     :EndEnvironmentalMinFlow
     */
      if(Options.noisy) { cout <<"Environmental Minimum Flow"<<endl; }
      long SBID=DOESNT_EXIST;
      CSubBasin *pSB;
      if(Len>=2) { SBID=s_to_l(s[1]); }
      pSB=pModel->GetSubBasinByID(SBID);
      pTimeSer=CTimeSeries::Parse(p,false,"EnviroMinFlow_"+to_string(SBID),SBID,"none",Options);
      if(pSB!=NULL) {
        pSB->AddEnviroMinFlow(pTimeSer);
      }
      else
      {
        warn=":EnvironmentalMinFlow: Subbasin "+to_string(SBID)+" not in model, cannot set minimum flow time series";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case (69): //---------------------------------------------
    {/*:UnusableFlowPercentage [SBID] [fraction]*/
      if(Options.noisy) { cout <<"Unusable flow percentage"<<endl; }
      if(Len<3) { p->ImproperFormat(s); break; }
      CSubBasin *pBasin=pModel->GetSubBasinByID(s_to_l(s[1]));
      if(pBasin==NULL) {
        ExitGracefully("Invalid subbasin ID in :UnusableFlowPercentage command.",BAD_DATA_WARN);
      }
      else {
        pBasin->SetUnusableFlowPercentage(s_to_d(s[2]));
      }
      break;
    }
    case (70): //---------------------------------------------
    {/*:MonthlyAveTemperature {temp x 12}*/
      if (Options.noisy) {cout <<"Monthly Average Temperatures"<<endl;}
      ExitGracefullyIf(pGage==NULL,
                       "ParseTimeSeriesFile::Monthly Average Temperatures declared before specifying a gauge station and its properties",BAD_DATA);
      if (Len!=13){p->ImproperFormat(s); break;}
      for (int i=0;i<12;i++){monthdata[i]=s_to_d(s[i+1]);}
      pGage->SetMonthlyAveTemps(monthdata);
      break;
    }
    case (71): //---------------------------------------------
    {/*:MonthlyAveEvaporation {PET x 12}*/
      if (Options.noisy) {cout <<"Monthly Average Evaporation"<<endl;}
      ExitGracefullyIf(pGage==NULL,
                       "ParseTimeSeriesFile::Monthly Average Evaporation declared before specifying a gauge station and its properties",BAD_DATA);
      if (Len < 13){p->ImproperFormat(s); break;}
      for (int i=0;i<12;i++){monthdata[i]=s_to_d(s[i+1]);}
      pGage->SetMonthlyPET(monthdata);
      break;
    }
    case (72): //---------------------------------------------
    {/*:MonthlyMinTemperature {temp x 12}*/
      if (Options.noisy) {cout <<"Monthly Minimum Temperatures"<<endl;}
      ExitGracefullyIf(pGage==NULL,
                       "ParseTimeSeriesFile::Monthly Minimum Temperatures declared before specifying a gauge station and its properties",BAD_DATA);
      if (Len!=13){p->ImproperFormat(s); break;}
      for (int i=0;i<12;i++){monthdata[i]=s_to_d(s[i+1]);}
      pGage->SetMonthlyMinTemps(monthdata);
      break;
    }
    case (73): //---------------------------------------------
    {/*:MonthlyMaxTemperature {temp x 12}*/
      if (Options.noisy) {cout <<"Monthly Maximum Temperatures"<<endl;}
      ExitGracefullyIf(pGage==NULL,
                       "ParseTimeSeriesFile::Monthly Maximum Temperatures declared before specifying a gauge station and its properties",BAD_DATA);
      if (Len!=13){p->ImproperFormat(s); break;}
      for (int i=0;i<12;i++){monthdata[i]=s_to_d(s[i+1]);}
      pGage->SetMonthlyMaxTemps(monthdata);
      break;
    }
    case (74): //---------------------------------------------
    {/*:UserTimeSeries {name} [units]
        {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
        {double val} x nValues
      :EndUserTimeSeries
     */
      if(Options.noisy) { cout <<"User-specified Time Series"<<endl; }
      pTimeSer=CTimeSeries::Parse(p,false,s[1], DOESNT_EXIST, "none", Options);
      if(pModel->GetDemandOptimizer()!=NULL) {
        pModel->GetDemandOptimizer()->AddUserTimeSeries(pTimeSer);
      }
      else
      {
        warn=":User-specified time series found without being in management optimization mode. This command will be ignored.";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case(100)://----------------------------------------------
    {/*:OverrideStreamflow  [SBID]*/
      if (Options.noisy){cout <<"Override streamflow"<<endl;}
      long SBID=s_to_l(s[1]);
      if (pModel->GetSubBasinByID(SBID)==NULL){
        WriteWarning("ParseTimeSeries::Trying to override streamflow at non-existent subbasin "+to_string(SBID),Options.noisy);
        break;
      }
      pModel->OverrideStreamflow(SBID);
      break;
    }
    case(101)://----------------------------------------------
    {/*:AssimilateStreamflow  [SBID]*/
      if(Options.noisy) { cout <<"Assimilate streamflow"<<endl; }
      long SBID=s_to_l(s[1]);
      if(pModel->GetSubBasinByID(SBID)==NULL) {
        WriteWarning("ParseTimeSeries::Trying to assimilate streamflow at non-existent subbasin "+to_string(SBID),Options.noisy);
        break;
      }
      bool ObsExists=false;
      for(int i=0; i<pModel->GetNumObservedTS(); i++) {
        if(IsContinuousFlowObs2(pModel->GetObservedTS(i),SBID)) {
          ObsExists=true;
          break;
        }
      }
      if(ObsExists) {
        pModel->GetSubBasinByID(SBID)->IncludeInAssimilation();
      }
      else {
        warn="ParseTimeSeriesFile::AssimilateStreamflow: no observation time series associated with subbasin "+to_string(SBID)+". Cannot assimilate flow in this subbasin.";
        WriteWarning(warn.c_str(),Options.noisy);
      }
      break;
    }
    case(102)://----------------------------------------------
    {/*:SpecifyGaugeForHRU [HRUID] [GaugeName]*/
      if (Options.noisy){cout <<"Specify gauge for HRU"<<endl;}
      long long int HRUID=s_to_ll(s[1]);
      if (pModel->GetHRUByID(HRUID)==NULL){
        WriteWarning("ParseTimeSeries::invalid HRU ID in :SpecifyGaugeForHRU command: "+to_string(HRUID),Options.noisy);
        break;
      }
      int g=pModel->GetGaugeIndexFromName(s[2]);
      if (g==DOESNT_EXIST){
        WriteWarning("ParseTimeSeries: invalid gauge name in :SpecifyGaugeForHRU command ("+to_string(s[2])+")",Options.noisy);
        break;
      }
      pModel->GetHRUByID(HRUID)->SetSpecifiedGaugeIndex(g);
      break;
    }
    case (300)://----------------------------------------------
    {/*:FixedConcentrationTimeSeries
       :FixedConcentrationTimeSeries [string constit_name] [string state_var (storage compartment)] {optional HRU Group name}
         {yyyy-mm-dd hh:mm:ss double tstep int nMeasurements}
         {double concentration values} x nMeasurements
       :EndFixedConcentrationTimeSeries
       or
       :FixedTemperatureTimeSeries [string state_var (storage compartment)] {optional HRU Group name}
         {yyyy-mm-dd hh:mm:ss double tstep int nMeasurements}
         {double temperature values} x nMeasurements
       :EndFixedConcentrationTimeSeries
     */
      if (Options.noisy){cout <<"Fixed concentration/temperature time series"<<endl;}

      int layer_ind;
      int i_stor;
      int shift=0;
      if(const_name=="TEMPERATURE") { shift=-1; }
      sv_type typ = pModel->GetStateVarInfo()->StringToSVType(s[2+shift],layer_ind,false);
      if (typ==UNRECOGNIZED_SVTYPE){
        WriteWarning(":ConcentrationTimeSeries command: unrecognized storage variable name: "+to_string(s[2]),Options.noisy);
        break;
      }
      i_stor=pModel->GetStateVarIndex(typ,layer_ind);
      if (i_stor!=DOESNT_EXIST){
        int kk=DOESNT_EXIST;
        if (Len>3+shift){
          CHRUGroup *pSourceGrp;
          pSourceGrp=pModel->GetHRUGroup(s[3+shift]);
          if (pSourceGrp==NULL){
            ExitGracefully("Invalid HRU Group name supplied in :FixedConcentrationTimeSeries or :FixedTemperatureTimeSeries command in .rvt file",BAD_DATA_WARN);
            break;
          }
          else{
            kk=pSourceGrp->GetGlobalIndex();
          }
        }
        CTimeSeries *pTS;
        pTS=CTimeSeries::Parse(p,true,const_name+"_"+to_string(s[2+shift]),DOESNT_EXIST,"none",Options);//name=constitutent name+storage name

        int c=pModel->GetTransportModel()->GetConstituentIndex(const_name);
        if(c!=DOESNT_EXIST) {
          pModel->GetTransportModel()->GetConstituentModel(c)->AddDirichletTimeSeries(i_stor,kk,pTS);
          pTS->SetConstitInd(c);
        }
        else {
          ExitGracefully("ParseMainInputFile: invalid constiuent in :FixedConcentrationTimeSeries or :FixedTemperatureTimeSeries command in .rvt file",BAD_DATA_WARN);
        }
      }
      else{
        ExitGracefully(":FixedConcentrationTimeSeries or :FixedTemperatureTimeSeries command: invalid state variable name",BAD_DATA_WARN);
      }
      break;
    }
    case (301)://----------------------------------------------
    {/*:MassFluxTimeSeries
       :MassFluxTimeSeries [string constit_name] [string state_var (storage compartment)] {optional HRU Group name}
         {yyyy-mm-dd hh:mm:ss double tstep int nMeasurements}
         {double flux values, in mg/m2/d} x nMeasurements
       :EndMassFluxTimeSeries
     */
      if (Options.noisy){cout <<"Mass flux time series"<<endl;}

      int layer_ind;
      int i_stor;
      string const_name = to_string(s[1]);
      sv_type typ = pModel->GetStateVarInfo()->StringToSVType(s[2], layer_ind, false);
      if (typ==UNRECOGNIZED_SVTYPE){
        WriteWarning(":MassFluxTimeSeries command: unrecognized storage variable name: "+to_string(s[2]),Options.noisy);
        break;
      }
      i_stor=pModel->GetStateVarIndex(typ,layer_ind);
      if (i_stor!=DOESNT_EXIST){
        int kk=DOESNT_EXIST;
        if (Len>3)
        {
          CHRUGroup *pSourceGrp;
          pSourceGrp=pModel->GetHRUGroup(s[3]);
          if (pSourceGrp==NULL){
            ExitGracefully("Invalid HRU Group name supplied in :MassFluxTimeSeries command in .rvt file",BAD_DATA_WARN);
            break;
          }
          else{
            kk=pSourceGrp->GetGlobalIndex();
          }
        }
        int c=pModel->GetTransportModel()->GetConstituentIndex(s[1]);

        CTimeSeries *pTS;
        pTS=CTimeSeries::Parse(p,true,const_name+"_"+to_string(s[2]),DOESNT_EXIST,"none",Options);//name=constitutent name

        if(c!=DOESNT_EXIST) {
          pModel->GetTransportModel()->GetConstituentModel(c)->AddInfluxTimeSeries(i_stor,kk,pTS);
          pTS->SetConstitInd(c);
        }
        else {
          cout<<"Bad constituent: "<< s[1]<<endl;
          ExitGracefully("ParseTimeSeriesInputFile: invalid constiuent in :MassFluxTimeSeries command in .rvt file",BAD_DATA_WARN);
        }
      }
      else{
        ExitGracefully(":MassFluxTimeSeries command: invalid state variable name",BAD_DATA_WARN);
      }
      break;
    }
    case(302):
    {/*:SpecifiedInflowConcentration [constit_name] [SBID]
      {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
      {double value} x nMeasurements
       :EndSpecifiedInflowConcentration
        */
      if(Options.noisy) { cout <<"Specified stream concentration"<<endl; }

      int c=pModel->GetTransportModel()->GetConstituentIndex(s[1]);
      if(c==DOESNT_EXIST) {
        ExitGracefully("ParseTimeSeriesFile: :SpecifiedStreamConcentration: invalid constituent name. Command will be ignored.",BAD_DATA_WARN); break;
      }
      long SBID=s_to_l(s[2]);
      if(pModel->GetSubBasinByID(SBID)==NULL) {
        ExitGracefully("ParseTimeSeriesFile: :SpecifiedStreamConcentration: invalid subbasin ID. Command will be ignored.",BAD_DATA_WARN); break;
      }

      CTimeSeries *pTS;
      pTS=CTimeSeries::Parse(p,true,"Specified Conc "+pModel->GetTransportModel()->GetConstituentTypeName(c)+"_"+to_string(SBID),SBID,"none",Options);
      pTS->SetLocID(SBID);

      pModel->GetTransportModel()->GetConstituentModel(c)->AddInflowConcTimeSeries(pTS);

      break;

      /*
      :OverrideStreamConcentration [constit_ind] [SBID]
        {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
        {double value} x nMeasurements
      :EndOverrideStreamConcentration

      long SBID=s_to_l(s[2]);
      if(pModel->GetSubBasinByID(SBID)==NULL) {
        WriteWarning("ParseTimeSeries::Trying to override stream concentration at non-existent subbasin "+to_string(SBID),Options.noisy);
        break;
      }
      int c=pModel->GetTransportModel()->GetConstituentIndex(s[2]);
      if(c==DOESNT_EXIST) {
        WriteWarning("ParseTimeSeriesFile: :OverrideStreamConcentration: invalid constituent name. Command will be ignored.",Options.noisy);
      }
      else {
        pModel->GetTransportModel()->OverrideStreamConcentration(SBID,c);
      }
      break;*/
    }
    case(303):
    {/*:SpecifiedInflowTemperature [SBID]
      {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
      {double value} x nMeasurements
       :EndSpecifiedInflowTemperature
        */
      if(Options.noisy) { cout <<"Specified stream temperature"<<endl; }

      int c=pModel->GetTransportModel()->GetConstituentIndex("TEMPERATURE");
      if(c==DOESNT_EXIST) {
        ExitGracefully("ParseTimeSeriesFile: :SpecifiedStreamTemperature: temperature is not being simulated. ",BAD_DATA_WARN);
      }
      long SBID=s_to_l(s[1]);
      if(pModel->GetSubBasinByID(SBID)==NULL) {
        ExitGracefully("ParseTimeSeriesFile: :SpecifiedStreamConcentration: invalid subbasin ID. ",BAD_DATA_WARN);
      }

      CTimeSeries *pTS;
      pTS=CTimeSeries::Parse(p,true,"Specified Conc "+pModel->GetTransportModel()->GetConstituentTypeName(c)+"_"+to_string(SBID),SBID,"none",Options);
      pTS->SetLocID(SBID);

      pModel->GetTransportModel()->GetConstituentModel(c)->AddInflowConcTimeSeries(pTS);

      break;
    }
    case(304)://----------------------------------------------
    {/*:MassLoading [constit_name] [SBID]
         {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
         {double value [kg/d]} x nMeasurements
       :EndMassLoading
       */
      if(Options.noisy) { cout <<"Specified mass loading to stream"<<endl; }

      int c=pModel->GetTransportModel()->GetConstituentIndex(s[1]);
      if(c==DOESNT_EXIST) {
        ExitGracefully("ParseTimeSeriesFile: :MassLoading: invalid constituent name. Command will be ignored.",BAD_DATA_WARN); break;
      }
      long SBID=s_to_l(s[2]);
      if(pModel->GetSubBasinByID(SBID)==NULL) {
        ExitGracefully("ParseTimeSeriesFile: :MassLoading: invalid subbasin ID. Command will be ignored.",BAD_DATA_WARN); break;
      }

      CTimeSeries* pTS;
      pTS=CTimeSeries::Parse(p,true,"Mass Loading "+pModel->GetTransportModel()->GetConstituentTypeName(c)+"_"+to_string(SBID),SBID,"none",Options);
      pTS->SetLocID(SBID);

      pModel->GetTransportModel()->GetConstituentModel(c)->AddMassLoadingTimeSeries(pTS);

      break;
    }

    case (400)://----------------------------------------------
    {/*:GriddedForcing
         :ForcingType PRECIP
         :FileNameNC  /Users/mai/Desktop/Tolson_SVN/btolson/trunk/basins/york_lumped/York_daily.nc
         :VarNameNC   pre
         :DimNamesNC  nlon nlat ntime
         :GridWeights # GRIDWEIGHTS COMMAND ALWAYS LAST!
           [#HRUs] [#GRID CELLS NC]
           [CELL#] [HRUID] [w_lk]
           ....
         :EndGridWeights
       :EndGriddedForcing
     */
      if (Options.noisy){cout <<":GriddedForcing"<<endl;}

#ifndef _RVNETCDF_
      ExitGracefully("ParseTimeSeriesFile: :GriddedForcing blocks are only allowed when NetCDF library is available!",BAD_DATA);
#endif
      string tmp[3];
      tmp[0] = "NONE";
      tmp[1] = "NONE";
      tmp[2] = "NONE";
      is_3D  = true;
      pGrid=new CForcingGrid("NONE",     // ForcingType,
                             "NONE",     // FileNameNC
                             "NONE",     // VarNameNC,
                             tmp,        // DimNames[3]
                             is_3D);     // is_3D

      grid_initialized = false;

      break;
    }
    case (406)://----------------------------------------------
    {/*:EndGriddedForcing*/
      if(Options.noisy) { cout <<":EndGriddedForcing "<<endl; }
#ifndef _RVNETCDF_
      ExitGracefully("ParseTimeSeriesFile: :GriddedForcing blocks are only allowed when NetCDF library is available!",BAD_DATA);
#endif
      if(pGrid->GetNumberNonZeroGridCells()==0) { //Never called SetIdxNonZeroGridCells()
        ExitGracefully("ParseTimeSeriesFile: :GriddedForcing command is missing :GridWeights or :StationWeightsByAttribute entry",BAD_DATA_WARN);
      }
      else {
        if(!grid_initialized) { pGrid->ForcingGridInit(Options);grid_initialized = true; }
        pModel->AddForcingGrid(pGrid,pGrid->GetForcingType());
      }
      pGrid=NULL;
      break;
    }
    case (401)://----------------------------------------------
    {/*:ForcingType [forcing type, e.g., PRECIP] */
      if (Options.noisy){cout <<"   :ForcingType"<<endl;}
#ifndef _RVNETCDF_
      ExitGracefully("ParseTimeSeriesFile: :GriddedForcing and :StationForcing blocks are only allowed when NetCDF library is available!",BAD_DATA);
#endif
      ExitGracefullyIf(pGrid==NULL, "ParseTimeSeriesFile: :ForcingType command must be within a :GriddedForcing or :StationForcing block",BAD_DATA);
      ExitGracefullyIf(grid_initialized,"ParseTimeSeriesFile: :ForcingType command must be before :GaugeWeights command",BAD_DATA);

      pGrid->SetForcingType(GetForcingTypeFromString(s[1]));

      break;
    }
    case (402)://----------------------------------------------
    {/*:FileNameNC  [filename] */
      if (Options.noisy){cout <<"   :FileNameNC"<<endl;}
#ifndef _RVNETCDF_
      ExitGracefully("ParseTimeSeriesFile: :GriddedForcing and :StationForcing blocks are only allowed when NetCDF library is available!",BAD_DATA);
#else
      ExitGracefullyIf(pGrid==NULL,     "ParseTimeSeriesFile: :FileNameNC command must be within a :GriddedForcings or :StationForcing block",BAD_DATA);
      ExitGracefullyIf(grid_initialized,"ParseTimeSeriesFile: :FileNameNC command must be before :GaugeWeights command",BAD_DATA);

      string filename=s[1];
      filename =CorrectForRelativePath(filename ,Options.rvt_filename);

      // check for existence

      int    retval;                // error value for NetCDF routines
      retval = nc_open(filename.c_str(),NC_NOWRITE,&ncid);
      if (retval != 0){
	      string warn = "ParseTimeSeriesFile: :FileNameNC command: Cannot find gridded data file "+ filename;
	      ExitGracefully(warn.c_str(),BAD_DATA_WARN);
	      break;
      }
      HandleNetCDFErrors(retval);

      pGrid->SetFilename(filename);

#endif
      break;
    }
    case (403)://----------------------------------------------
    {/*:VarNameNC   [name of :ForcingType variable in NetCDF file] */
      if (Options.noisy){cout <<"   :VarNameNC"<<endl;}
#ifndef _RVNETCDF_
      ExitGracefully("ParseTimeSeriesFile: :GriddedForcing and :StationForcing blocks are only allowed when NetCDF library is available!",BAD_DATA);
#endif
      ExitGracefullyIf(pGrid==NULL,     "ParseTimeSeriesFile: :VarNameNC command must be within a :GriddedForcing or :StationForcing block",BAD_DATA);
      ExitGracefullyIf(grid_initialized,"ParseTimeSeriesFile: :VarNameNC command must be before :GaugeWeights command",BAD_DATA);
      pGrid->SetVarname(s[1]);
      break;
    }
    case (404)://----------------------------------------------
    {/*:DimNamesNC  [longitude netCDF alias] [latitude netCDF alias] [time netCDF alias]*/
      if (Options.noisy){cout <<"   :DimNamesNC"<<endl;}
#ifndef _RVNETCDF_
      ExitGracefully("ParseTimeSeriesFile: :GriddedForcing and :StationForcing blocks are only allowed when NetCDF library is available!",BAD_DATA);
#endif
      ExitGracefullyIf(pGrid==NULL,      "ParseTimeSeriesFile: :DimNamesNC command must be within a :GriddedForcing or :StationForcing block",BAD_DATA);
      ExitGracefullyIf(grid_initialized, "ParseTimeSeriesFile: :DimNamesNC command must be before :GaugeWeights command",BAD_DATA);
      string tmp[3];
      tmp[0] = s[1];
      tmp[1] = s[2];
      if (is_3D){ tmp[2] = s[3]; }

      pGrid->SetDimNames(tmp);

      break;
    }
    case (405)://----------------------------------------------
    {/*:GridWeights
       :NumberHRUs 3
       :NumberGridCells  22
       # [HRU ID] [Cell #] [w_kl]
       #     where w_kl is fraction of forcing for HRU k is from grid cell l=(i,j)
       #     and grid cell index l is derived by l = j * NC + i
       #     where i and j are the zero-indexed column and row of cell l respectively and
       #     NC is the total number of columns.
       #     minimum Cell # is 0, max is nrow*ncol-1
       #     Following contraint must be satisfied:
       #         sum(w_kl, {l=1,NC*NR}) = 1.0 for all HRUs k
       # [HRUID] [CELL#] [w_kl]
       1       2       0.3
       1       3       0.5
       1       23      0.2
       2       2       0.4
       2       3       0.6
       3       5       1.0
       :EndGridWeights*/
#ifndef _RVNETCDF_
      ExitGracefully("ParseTimeSeriesFile: :GriddedForcing and :StationForcing blocks are only allowed when NetCDF library is available!",BAD_DATA);
#endif

      ExitGracefullyIf(pGrid==NULL,
                       "ParseTimeSeriesFile: :GridWeights command must be within a :GriddedForcing or :StationForcing block",BAD_DATA);

      if (!grid_initialized) { //must initialize grid prior to adding grid weights
        grid_initialized = true;
        pGrid->ForcingGridInit(Options);
      }

      bool nHydroUnitsGiven = false;
      bool nGridCellsGiven  = false;
      int  nHydroUnits=0;
      int  nGridCells=0;
      CHydroUnit *pHRU=NULL;

      if (Options.noisy) {cout <<"GridWeights..."<<endl;}
      while (((Len==0) || (strcmp(s[0],":EndGridWeights"))) && (!(p->Tokenize(s,Len))))
      {

        if      (IsComment(s[0],Len))            {}//comment line
        else if (!strcmp(s[0],":Attributes"    )){}//ignored by Raven - needed for GUIs
        else if (!strcmp(s[0],":Units"         )){}//ignored by Raven - needed for GUIs
        else if (!strcmp(s[0],":NumberHRUs"    ))
        {
          nHydroUnits      = atoi(s[1]);
          nHydroUnitsGiven = true;

          ExitGracefullyIf(pModel->GetNumHRUs() < nHydroUnits,
                           "ParseTimeSeriesFile: :GridWeights :NumberHRUs exceeds number of HRUs in model",BAD_DATA);

          nHydroUnits=pModel->GetNumHRUs(); //no need to use this; should always use actual number of HRUs in model
          pGrid->SetnHydroUnits(nHydroUnits);
          if (nHydroUnitsGiven && nGridCellsGiven) {pGrid->AllocateWeightArray(nHydroUnits,nGridCells);}
        }
        else if (!strcmp(s[0],":NumberGridCells"    ) || !strcmp(s[0],":NumberStations"    ))
        {
          nGridCells      = atoi(s[1]);
          nGridCellsGiven = true;

          if (pGrid->GetCols() * pGrid->GetRows() != nGridCells) {
            printf(":NumberGridCells/:NumberStations   = %i\n",atoi(s[1]));
            printf("NetCDF cols * rows = %i\n",pGrid->GetCols() * pGrid->GetRows());
            ExitGracefully("ParseTimeSeriesFile: :NumberGridCells given does not agree with NetCDF file content",BAD_DATA);
          }

          if (nHydroUnitsGiven && nGridCellsGiven) {pGrid->AllocateWeightArray(nHydroUnits,nGridCells);}
        }
        else if (!strcmp(s[0],":EndGridWeights")){}//done
        else
        {
          if (nHydroUnitsGiven && nGridCellsGiven) {
            pHRU=NULL;
            pHRU = pModel->GetHRUByID(atoll(s[0]));
            if (pHRU == NULL) {
              printf("\n\n");
              printf("Wrong HRU ID in :GridWeights: HRU_ID = %s\n",s[0]);
              ExitGracefully("ParseTimeSeriesFile: HRU ID found in :GridWeights which does not exist in :HRUs!",BAD_DATA);
            }
            pGrid->SetWeightVal(pHRU->GetGlobalIndex(),atoi(s[1]),atof(s[2]));
          }
          else {
            ExitGracefully("ParseTimeSeriesFile: :NumberHRUs must be given in :GridWeights block",BAD_DATA);
          }
        } // end else
      } // end while

      // check that weightings sum up to one per HRU
      bool WeightArrayOK = pGrid->CheckWeightArray(nHydroUnits,nGridCells,pModel);
      ExitGracefullyIf(!WeightArrayOK,
                       "ParseTimeSeriesFile: Check of weights for gridded forcing failed. Sum of gridweights for ALL enabled HRUs must be 1.0.",BAD_DATA);

      // store (sorted) grid cell ids with non-zero weight in array
      pGrid->SetIdxNonZeroGridCells(nHydroUnits,nGridCells,Options);
      pGrid->CalculateChunkSize(Options);
      break;
    }

    case (407)://----------------------------------------------
    {/*:Deaccumulate */
      if (Options.noisy){cout <<"   :Deaccumulate"<<endl;}
      ExitGracefullyIf(pGrid==NULL, "ParseTimeSeriesFile: :Deaccumulate command must be within a :GriddedForcing or :StationForcing block",BAD_DATA);
      pGrid->SetToDeaccumulate();
      break;
    }
    case (408)://----------------------------------------------
    {/*:TimeShift [hh:mm:ss] or [days]*/
      if (Options.noisy){cout <<"   :TimeShift"<<endl;}
      ExitGracefullyIf(pGrid==NULL,     "ParseTimeSeriesFile: :TimeShift command must be within a :GriddedForcing or :StationForcing block",BAD_DATA);
      ExitGracefullyIf(Len<2,           "ParseTimeSeriesFile: :TimeShift expects at least one argument",BAD_DATA);
      ExitGracefullyIf(grid_initialized,"ParseTimeSeriesFile: :TimeShift argument in :GriddedForcing or :StationForcing block needs to be before :GridWeights block",BAD_DATA);

      string tString=s[1];
      double TimeShift=0.0;
      if((tString.length()>=2) && ((tString.substr(2,1)==":") || (tString.substr(1,1)==":"))) {//support for hh:mm:ss.00 format in timestep
        time_struct tt;
        tt=DateStringToTimeStruct("0000-01-01",tString,Options.calendar);
        TimeShift=(tt.julian_day);
      }
      else {
        TimeShift =(s_to_d(s[1]));//in days
      }

      pGrid->SetTimeShift(TimeShift);

      break;
    }
    case (409)://----------------------------------------------
    {/*:LinearTransform [a] [b] */
      if (Options.noisy){cout <<"   :LinearTransform"<<endl;}
      ExitGracefullyIf(pGrid==NULL,     "ParseTimeSeriesFile: :LinearTransform command must be within a :GriddedForcing or :StationForcing block",BAD_DATA);
      ExitGracefullyIf(Len!=3,          "ParseTimeSeriesFile: :LinearTransform expects exactly two arguments",BAD_DATA);
      ExitGracefullyIf(grid_initialized,"ParseTimeSeriesFile: :LinearTransform command must be before :GaugeWeights command",BAD_DATA);

      double LinTrans_a=atof(s[1]);
      double LinTrans_b=atof(s[2]);
      pGrid->SetLinearTransform(LinTrans_a,LinTrans_b);
      break;
    }
    case (410)://----------------------------------------------
    {/*:PeriodEndingNC  */
      if (Options.noisy){cout <<"   :PeriodEnding"<<endl;}
      ExitGracefullyIf(pGrid==NULL     ,"ParseTimeSeriesFile: :PeriodEnding command must be within a :GriddedForcing or :StationForcing block",BAD_DATA);
      ExitGracefullyIf(grid_initialized,"ParseTimeSeriesFile: :PeriodEndingNC command must be before :GaugeWeights command",BAD_DATA);
      pGrid->SetAsPeriodEnding();
      break;
    }
    case (411)://----------------------------------------------
    {/*:LatitudeVarNameNC  [variablename]*/
      if(Options.noisy) { cout <<"   :LatitudeVarNameNC"<<endl; }
      ExitGracefullyIf(pGrid==NULL,"ParseTimeSeriesFile: :LatitudeVarNameNC command must be within a :GriddedForcing or :StationForcing block",BAD_DATA);
      pGrid->SetAttributeVarName("Latitude",s[1]);
      break;
    }
    case (412)://----------------------------------------------
    {/*:LongitudeVarNameNC  [variablename]*/
      if(Options.noisy) { cout <<"   :LongitudeVarNameNC"<<endl; }
      ExitGracefullyIf(pGrid==NULL,"ParseTimeSeriesFile: :LongitudeVarNameNC command must be within a :GriddedForcing or :StationForcing block",BAD_DATA);
      pGrid->SetAttributeVarName("Longitude",s[1]);
      break;
    }
    case (413)://----------------------------------------------
    {/*:ElevationVarNameNC  [variablename]*/
      if(Options.noisy) { cout <<"   :ElevationVarNameNC"<<endl; }
      ExitGracefullyIf(pGrid==NULL,"ParseTimeSeriesFile: :ElevationVarNameNC command must be within a :GriddedForcing or :StationForcing block",BAD_DATA);
      pGrid->SetAttributeVarName("Elevation",s[1]);
      break;
    }
    case (416)://----------------------------------------------
    {/*:StationIDNameNC  [variablename]*/
      if(Options.noisy) { cout <<"   :StationIDNameNC"<<endl; }
      ExitGracefullyIf(pGrid==NULL,"ParseTimeSeriesFile: :StationIDNameNC command must be within a :StationForcing block",BAD_DATA);
      pGrid->SetAttributeVarName("StationIDs",s[1]);
      break;
    }
    case(414)://----------------------------------------------
    {/*:StationWeightsByAttribute
       :NumberHRUs 3
       :NumberStations  22
       # [HRU ID] [station Identifier] [w_kl]
       #     Following contraint must be satisfied:
       #         sum(w_kl, {l=1,NS}) = 1.0 for all HRUs k where NS=number of stations
       1       A3       0.3
       1       A34       0.5
       1       B2      0.2
       2       B0       0.4
       2       C2       0.6
       3       D15      1.0 # station identifier consistent with array in NetCDF consistent with StationIDNameNC
       :EndStationWeightsByAttribute*/
#ifndef _RVNETCDF_
      ExitGracefully("ParseTimeSeriesFile: :GriddedForcing and :StationForcing blocks are only allowed when NetCDF library is available!",BAD_DATA);
#endif
      ExitGracefullyIf(pGrid==NULL,
        "ParseTimeSeriesFile: :StationWeightsByAttribute command must be within a :StationForcing block",BAD_DATA);

      if(!grid_initialized) { //must initialize grid prior to adding grid weights
        grid_initialized = true;
        pGrid->ForcingGridInit(Options);
      }

      bool nHydroUnitsGiven = false;
      bool nGridCellsGiven  = false;
      bool attributeGiven=false;
      int  nHydroUnits=0;
      int  nGridCells=0;

      if(Options.noisy) { cout <<"StationWeightsByAttribute..."<<endl; }
      while(((Len==0) || (strcmp(s[0],":EndStationWeightsByAttribute"))) && (!(p->Tokenize(s,Len))))
      {

        if(IsComment(s[0],Len)) {}//comment line
        else if(!strcmp(s[0],":NumberHRUs"))
        {
          nHydroUnits      = atoi(s[1]);
          nHydroUnitsGiven = true;

          ExitGracefullyIf(pModel->GetNumHRUs() < nHydroUnits,
            "ParseTimeSeriesFile: :StationWeightsByAttribute :NumberHRUs exceeds number of HRUs in model",BAD_DATA);

          nHydroUnits=pModel->GetNumHRUs(); //no need to use this; should always use actual number of HRUs in model

          pGrid->SetnHydroUnits(nHydroUnits);
          if(nHydroUnitsGiven && nGridCellsGiven) { pGrid->AllocateWeightArray(nHydroUnits,nGridCells); }
        }
        else if( !strcmp(s[0],":NumberStations"))
        {
          nGridCells      = atoi(s[1]);
          nGridCellsGiven = true;

          if(pGrid->GetCols() * pGrid->GetRows() != nGridCells) {
            printf(":NumberStations   = %i\n",atoi(s[1]));
            printf("NetCDF cols * rows = %i\n",pGrid->GetCols() * pGrid->GetRows());
            ExitGracefully("ParseTimeSeriesFile: :NumberStations given does not agree with NetCDF file content",BAD_DATA);
          }

          if(nHydroUnitsGiven && nGridCellsGiven) { pGrid->AllocateWeightArray(nHydroUnits,nGridCells); }
        }
        else if(!strcmp(s[0],":EndGridWeights")) {}//done
        else
        {
          if(nHydroUnitsGiven && nGridCellsGiven && attributeGiven) {
            CHydroUnit *pHRU=NULL;
            pHRU = pModel->GetHRUByID(atoll(s[0]));
            if(pHRU == NULL) {
              printf("\n\n");
              printf("Wrong HRU ID in :StationWeightsByAttribute: HRU_ID = %s\n",s[0]);
              ExitGracefully("ParseTimeSeriesFile: HRU ID found in :StationWeightsByAttribute which does not exist in :HRUs!",BAD_DATA);
            }

            ExitGracefully("StationWeightsByAttribute",STUB);

            int ind=DOESNT_EXIST;//pGrid->GetStationIndex(s[1]);
            pGrid->SetWeightVal(pHRU->GetGlobalIndex(),ind,atof(s[2]));
          }
          else {
            ExitGracefully("ParseTimeSeriesFile: :NumberHRUs must be given in :StationWeightsByAttribute block",BAD_DATA);
          }
        } // end else
      } // end while

      // check that weightings sum up to one per HRU
      bool WeightArrayOK = pGrid->CheckWeightArray(nHydroUnits,nGridCells,pModel);
      ExitGracefullyIf(!WeightArrayOK,
        "ParseTimeSeriesFile: Check of weights for gridded forcing failed. Sum of gridweights for ALL enabled HRUs must be 1.0.",BAD_DATA);

      // store (sorted) grid cell ids with non-zero weight in array
      pGrid->SetIdxNonZeroGridCells(nHydroUnits,nGridCells,Options);
      pGrid->CalculateChunkSize(Options);
      break;

    }
    case (415)://----------------------------------------------
    {/*:StationElevations
     [int station ID1] [elev1]
     [station ID2] [elev1]
     ...
     [station IDN] [elevN]
     :EndStationElevations
     # where station IDs consistent with StationIDNameNC
     */
      ExitGracefullyIf(!grid_initialized,"ParseTimeSeriesFile: :StationElevations command must be after :GaugeWeights command",BAD_DATA);
      ExitGracefullyIf(pGrid==NULL,
        "ParseTimeSeriesFile: :StationElevations command must be within a :GriddedForcing or :StationForcing block",BAD_DATA);

      if(!grid_initialized) { //must initialize grid prior to adding grid cell/station elevations
        grid_initialized = true;
        pGrid->ForcingGridInit(Options);
      }


      if(Options.noisy) { cout <<"Station Elevations..."<<endl; }
      while(((Len==0) || (strcmp(s[0],":EndStationElevations"))) && (!(p->Tokenize(s,Len))))
      {
        if(IsComment(s[0],Len)) {}//comment line
        else if(!strcmp(s[0],":EndStationElevations")) {}//done
        else
        {
          int ind=DOESNT_EXIST;//pGrid->GetStationIndex(s[1]);
          ExitGracefully(":EndStationElevations",STUB);
          pGrid->SetStationElevation(ind,s_to_d(s[1]));
        } // end else
      } // end while
      break;
    }
    case (417)://----------------------------------------------
    {/*:StationElevationsByIdx
       [station index1] [elev1]
       [station index2] [elev1]
       ...
       [station indexN] [elevN]
     :EndStationElevationsByIdx
     # where station indices are same as nc file order
     */
      ExitGracefullyIf(!grid_initialized,"ParseTimeSeriesFile: :StationElevationsByIdx command must be after :GaugeWeights command",BAD_DATA);
      ExitGracefullyIf(pGrid==NULL,
        "ParseTimeSeriesFile: :StationElevations command must be within a :GriddedForcing or :StationForcing block",BAD_DATA);

      if(!grid_initialized) { //must initialize grid prior to adding grid cell/station elevations
        grid_initialized = true;
        pGrid->ForcingGridInit(Options);
      }

      if(Options.noisy) { cout <<"Station Elevations..."<<endl; }
      while(((Len==0) || (strcmp(s[0],":EndStationElevationsByIdx"))) && (!(p->Tokenize(s,Len))))
      {
        if(IsComment(s[0],Len)) {}//comment line
        else if(!strcmp(s[0],":EndStationElevationsByIdx")) {}//done
        else
        {
          int idx=s_to_i(s[0]);
          pGrid->SetStationElevation(idx,s_to_d(s[1]));
        } // end else
      } // end while
      break;
    }
    case (418)://----------------------------------------------
    {/*:MapStationsTo [SUBBASINS or HRUS]*/
      //Alternate to :GridWeights - used if station forcings map 1:1 with subbasins or HRUs
      //requires array attribute stations in same NetCDF file as forcing data

      if(Options.noisy) { cout <<"   :MapStationsTo command"<<endl; }
      ExitGracefullyIf(pGrid==NULL,"ParseTimeSeriesFile: :MapStationsTo command must be within a :StationForcing block",BAD_DATA);

      bool sb_command;
      if      (!strcmp(s[1],"SUBBASINS")){sb_command=true;}
      else if (!strcmp(s[1],"HRUS"     )){sb_command=false;}

      // Get vector of station IDs
      //====================================================================
      int   StatVecID;         // attribute ID of station vector
      int   StatDimID;         // nstations dimension ID
      long* StatIDs    =NULL;  // vector of Subbasin or HRU IDs (allocated in GetNetCDFStationArray)
                               //\todo[funct]: fix so StatIDs supports long long int
      string *junk=NULL;
      int   nStations=0;       // size of station Vector
      if (ncid == -9) {
        ExitGracefully(":FileNameNC command must precede :MapStationsTo command in :StationForcings block",BAD_DATA);
      }
      GetNetCDFStationArray(ncid, pGrid->GetFilename(),StatDimID,StatVecID, StatIDs,junk, nStations);
      if ((sb_command) && (nStations!=pModel->GetNumSubBasins())){
        ExitGracefully(":MapStationsTo command: Number of stations in NetCDF file not the same as the number of model subbasins",BAD_DATA);
      }
      if ((!sb_command) && (nStations!=pModel->GetNumHRUs())){
        ExitGracefully(":MapStationsTo command: Number of stations in NetCDF file not the same as the number of model HRUs",BAD_DATA);
      }
      /*cout<<":MapTo command: "<<endl;
      for (int i=0;i<nStations;i++){
        cout<<" "<<StatIDs[i]<<endl;
      }*/
      // grid weights preparation
      //====================================================================
      if (!grid_initialized) { //must initialize grid prior to adding grid weights
        grid_initialized = true;
        pGrid->ForcingGridInit(Options);
      }
      pGrid->SetnHydroUnits     (pModel->GetNumHRUs());
      pGrid->AllocateWeightArray(pModel->GetNumHRUs(),nStations);

      // generate grid weights
      //====================================================================
      int StationID=DOESNT_EXIST;
      for (int k=0;k<pModel->GetNumHRUs();k++)
      {
        if (sb_command){ //SUBBASINS
          StationID=DOESNT_EXIST;
          int p=pModel->GetHydroUnit(k)->GetSubBasinIndex();
          long SBID=pModel->GetSubBasin(p)->GetID();
          for (int i=0;i<nStations;i++){
            if (SBID==StatIDs[i]){
              StationID=i; break;
            }
          }
        }
        else { //HRUS
          StationID=DOESNT_EXIST;
          for (int i=0;i<nStations;i++){
            if (pModel->GetHydroUnit(k)->GetHRUID()==StatIDs[i]){
              StationID=i; break;
            }
          }
        }
        int p=pModel->GetHydroUnit(k)->GetSubBasinIndex();
        //cout<<"SETTING WEIGHT "<<k<<" "<<StationID<<" "<<pModel->GetSubBasin(p)->GetID()<<endl;
        pGrid->SetWeightVal(k,StationID,1.0);
      }
      // check that weightings sum up to one per HRU
      bool WeightArrayOK = pGrid->CheckWeightArray(pModel->GetNumHRUs(),nStations,pModel);
      ExitGracefullyIf(!WeightArrayOK,
                       "ParseTimeSeriesFile: Check of weights for gridded forcing in :MapToStations command failed. Sum of gridweights for ALL enabled HRUs must be 1.0.",BAD_DATA);

      // store (sorted) station ids with non-zero weight in array
      pGrid->SetIdxNonZeroGridCells(pModel->GetNumHRUs(),nStations,Options);
      pGrid->CalculateChunkSize(Options);

      delete [] StatIDs;
      delete [] junk;
      break;
    }
    case (500)://----------------------------------------------
    {/*:StationForcing
         :ForcingType PRECIP
         :FileNameNC  [filename.nc]
         :VarNameNC   pre
         :DimNamesNC  nstations ntime

         :GridWeights
           :NumberHRUs [#HRUs]
           :NumberStations [#STATIONS]
           [HRUID] [STATION#] [w_lk]
           ....
         :EndGridWeights
         #OR
         :MapStationsToHRUs
         #OR
         :MapStationsToSubBasins
       :EndStationForcing
     */
      if (Options.noisy){cout <<":StationForcing"<<endl;}

#ifndef _RVNETCDF_
      ExitGracefully("ParseTimeSeriesFile: :StationForcing blocks are only allowed when NetCDF library is available!",BAD_DATA);
#endif

      string tmp[3];
      tmp[0] = "NONE";
      tmp[1] = "NONE";
      tmp[2] = "NONE";
      is_3D  = false;
      pGrid=new CForcingGrid("NONE",     // ForcingType,
                             "NONE",     // FileNameNC
                             "NONE",     // VarNameNC,
                             tmp,        // DimNames[3],
                             is_3D);     // is_3D
      grid_initialized = false;

      break;
    }
    case (506)://----------------------------------------------
    {/*:EndStationForcing*/
      if (Options.noisy){cout <<":EndStationForcing"<<endl;}
#ifndef _RVNETCDF_
      ExitGracefully("ParseTimeSeriesFile: :StationForcing blocks are only allowed when NetCDF library is available!",BAD_DATA);
#endif
      if(!grid_initialized) { pGrid->ForcingGridInit(Options);grid_initialized = true; }
      pModel->AddForcingGrid(pGrid,pGrid->GetForcingType());
      pGrid=NULL;
      ncid=-9;
      break;
    }

    default: //----------------------------------------------
    {
      char firstChar = *(s[0]);
      switch(firstChar)
      {
      case ':':
      {
        if     (!strcmp(s[0],":FileType"))    {if (Options.noisy){cout<<"Filetype"<<endl;}}//do nothing
        else if(!strcmp(s[0],":Application")) {if (Options.noisy){cout<<"Application"<<endl;}}//do nothing
        else if(!strcmp(s[0],":Version"))     {if (Options.noisy){cout<<"Version"<<endl;}}//do nothing
        else if(!strcmp(s[0],":WrittenBy"))   {if (Options.noisy){cout<<"WrittenBy"<<endl;}}//do nothing
        else if(!strcmp(s[0],":CreationDate")){if (Options.noisy){cout<<"CreationDate"<<endl;}}//do nothing
        else if(!strcmp(s[0],":SourceFile"))  {if (Options.noisy){cout<<"SourceFile"<<endl;}}//do nothing
        else
        {
          warn ="IGNORING unrecognized command: |" + string(s[0])+ "| in .rvt file "+p->GetFilename()+" line: "+ to_string(p->GetLineNumber());
          WriteWarning(warn,Options.noisy);
        }
      }
      break;
      default:
      {
        string errString = "Unrecognized command in .rvt file: |" + string(s[0]) + "| file: "+p->GetFilename()+" line: "+ to_string(p->GetLineNumber());
        ExitGracefully(errString.c_str(),BAD_DATA);//STRICT
      }
      break;
      }
    }
    }//end switch(code)

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

  RVT.close();

  //QA/QC
  //--------------------------------
  if((has_irrig) && (pModel->GetTransportModel()->GetNumConstituents()>0)) {
    WriteWarning("ParseTimeSeriesFile: irrigation/diversions included with transport constituents. Since water demands are not currently simulated in the Raven transport module, transport results must be interpreted with care.",Options.noisy);
  }

  delete p; p=NULL;

  return true;
}

//////////////////////////////////////////////////////////////////
/// \brief handles allocation of reservoir demands from downstream
/// \param *&pModel [out] Reference to the model object
/// \param Options [in] Global model options
//
void AllocateReservoirDemand(CModel *&pModel, const optStruct &Options,long SBID,long SBIDres,double pct_met,int jul_start,int jul_end)
{
  double dmult;
  double mult = pModel->GetGlobalParams()->GetParameter("RESERVOIR_DEMAND_MULT");
  string warn;

  CSubBasin *pSB, *pSBres;
  pSB = pModel->GetSubBasinByID(SBID);
  if(pSB==NULL) {
    warn=":AllocateReservoirDemand: Subbasin "+to_string(SBID)+" not in model, cannot set reservoir downstream demand";
    WriteWarning(warn,Options.noisy);
    return;
  }
  if(SBIDres==AUTO_COMPUTE_LONG)
  {
    if(Options.res_demand_alloc==DEMANDBY_CONTRIB_AREA)//==================================================
    {
      int nUpstr=0;
      const CSubBasin **pUpstr = pModel->GetUpstreamSubbasins(SBID,nUpstr);
      double Atot=0;
      for(int p=0;p<nUpstr;p++) {
        if(pUpstr[p]->GetReservoir()!=NULL) {
          Atot+=pUpstr[p]->GetDrainageArea();
        }
      }
      double area;

      if(Atot>0.0) { //reservoirs exist upstream
        for(int p=0;p<nUpstr;p++) {
          if(pUpstr[p]->GetReservoir()!=NULL) {
            /*  for(int i=0;i<pUpstr[p]->GetNumDemands();i++) {
            if(pUpstr[p]->GetReservoir()->GetDemand(i).DownSB==pSB->GetID()) {
            //demand already assigned; this overrides AUTO for this subbasin
            }
            }*/
            dmult=pUpstr[p]->GetReservoir()->GetDemandMultiplier();
            area=pUpstr[p]->GetDrainageArea();
            pUpstr[p]->GetReservoir()->AddDownstreamDemand(pSB,area/Atot*mult*dmult,jul_start,jul_end);
          }
        }
      }
    }
    //else if(Options.res_demand_alloc==DEMANDBY_SURFACE_AREA) {//================================================

    //}
    else if(Options.res_demand_alloc==DEMANDBY_MAX_CAPACITY)//==================================================
    {
      int nUpstr=0;
      const CSubBasin **pUpstr = pModel->GetUpstreamSubbasins(SBID,nUpstr);
      double Vtot=0;
      for(int p=0;p<nUpstr;p++) {
        if(pUpstr[p]->GetReservoir()!=NULL) {
          Vtot+=pUpstr[p]->GetReservoir()->GetMaxCapacity();
        }
      }
      double volume;

      if(Vtot>0.0) { //reservoirs exist upstream
        for(int p=0;p<nUpstr;p++) {
          if(pUpstr[p]->GetReservoir()!=NULL) {
            dmult =pUpstr[p]->GetReservoir()->GetDemandMultiplier();
            volume=pUpstr[p]->GetReservoir()->GetMaxCapacity();
            pUpstr[p]->GetReservoir()->AddDownstreamDemand(pSB,(volume/Vtot)*dmult*mult,jul_start,jul_end);
          }
        }
      }
    }
  }
  else /* if SBIDres!=AUTO_COMPUTE_LONG */
  { //single connection //===============================================================================
    pSBres=pModel->GetSubBasinByID(SBIDres);
    if(pSBres==NULL) {
      warn=":AllocateReservoirDemand: Reservoir subbasin "+to_string(SBID)+" not in model, cannot set reservoir downstream demand";
      WriteWarning(warn,Options.noisy);
    }
    else if(pSBres->GetReservoir()==NULL) {
      warn=":AllocateReservoirDemand: subbasin "+to_string(SBID)+" does not have reservoir, cannot set reservoir downstream demand";
      WriteWarning(warn,Options.noisy);
    }
    else
    {
      dmult=pSBres->GetReservoir()->GetDemandMultiplier();
      pSBres->GetReservoir()->AddDownstreamDemand(pSB,pct_met*dmult*mult,jul_start,jul_end);
    }
  }
}
