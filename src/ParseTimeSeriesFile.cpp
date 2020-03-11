/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2018 the Raven Development Team
  ----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Model.h"
#include "TimeSeries.h"
#include "IrregularTimeSeries.h"
#include "ParseLib.h"

//////////////////////////////////////////////////////////////////
/// \brief Parse input time series file, model.rvt
/// \details
///model.rvt: input file that defines temperature, precip, and other external vars at a set of gauges \n
/// .rvt format: \n
///   ":Gauge" [optional name]
///     - :Latitude {lat}
///     - :Longitude {long}
///     - :Elevation {elev}
///     - :Rain
///        {int nMeasurements int julian day int julian year double tstep} // start of data collection
///        {double precip} x nMeasurements
///     - :EndRain
///     - :Snow
///        {int nMeasurements int julian day int julian year double tstep} // start of data collection
///        {double snow_precip} x nMeasurements
///     - :EndSnow
///     - :MinimumTemperature
///        {int nMeasurements int julian day int julian year double tstep}
///        {double min_temp} x nMeasurements
///     - :EndMinimumTemperature
///     - :MultiData
///        {int nMeasurements int julian day int julian year double timestep}
///      -   :Parameters, PARAM_1_TAG, PARAM_2_TAG, ...
///      -  :Units,      unit1,        unit2, ...
///     {double param_1, double param_2, ...} x nMeasurements
///     - ":EndMultiData" \n
///
///   ":EndGauge" \n
///   ...
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
  bool ended(false);

  // some tmp variables for gridded forcing input
  bool grid_initialized = false;
  bool is_3D            = false;  //

  ifstream RVT;
  RVT.open(Options.rvt_filename.c_str());
  if (RVT.fail()){
    cout << "ERROR opening *.rvt file: "<<Options.rvt_filename<<endl; return false;}

  CParser *p=new CParser(RVT,Options.rvt_filename,line);

  ifstream INPUT2;           //For Secondary input
  CParser *pMainParser=NULL; //for storage of main parser while reading secondary files

  if (Options.noisy){
    cout <<"======================================================"<<endl;
    cout << "Parsing Forcing Data File " << Options.rvt_filename <<"..."<<endl;
    cout <<"======================================================"<<endl;
  }

  //Default Values
  //--Sift through file-----------------------------------------------
  bool end_of_file=p->Tokenize(s,Len);
  while (!end_of_file)
  {
    if (ended){break;}
    if (Options.noisy){ cout << "reading line " << p->GetLineNumber() << ": ";}

    /*assign code for switch statement
      --------------------------------------------------------------------------------------
      <100         : ignored/special
      0   thru 100 : data
      --------------------------------------------------------------------------------------
    */

    code=0;
    //---------------------SPECIAL -----------------------------
    if       (Len==0)                                 {code=-1; }
    else if  (!strcmp(s[0],"*"                      )){code=-2; }//comment
    else if  (!strcmp(s[0],"%"                      )){code=-2; }//comment
    else if  (!strcmp(s[0],"#"                      )){code=-2; }//comment
    else if  (s[0][0]=='#')                           {code=-2; }//comment
    else if  (!strcmp(s[0],":RedirectToFile"        )){code=-3; }//redirect to secondary file
    else if  (!strcmp(s[0],":End"                   )){code=-4; }//premature end of file
    //--------------------GAUGE BASIC DATA- --------------------
    else if  (!strcmp(s[0],":Gauge"                 )){code=1;  }
    else if  (!strcmp(s[0],":EndGauge"              )){code=2;  }
    else if  (!strcmp(s[0],":Latitude"              )){code=9;  }
    else if  (!strcmp(s[0],":Longitude"             )){code=10; }
    else if  (!strcmp(s[0],":Elevation"             )){code=11; }
    else if  (!strcmp(s[0],":GaugeList"             )){code=14; }
    else if  (!strcmp(s[0],":MeasurementHeight"     )){code=16; }
    //--------------------FORCING FUNCTIONS --------------------
    else if  (!strcmp(s[0],":Rain"                  )){code=3;  } //should make obsolete
    else if  (!strcmp(s[0],":Snow"                  )){code=4;  } //should make obsolete
    else if  (!strcmp(s[0],":MinTemperature"        )){code=5;  } //should make obsolete
    else if  (!strcmp(s[0],":MaxTemperature"        )){code=6;  } //should make obsolete
    else if  (!strcmp(s[0],":TotalPrecip"           )){code=7;  } //should make obsolete
    else if  (!strcmp(s[0],":Precipitation"         )){code=7;  } //should make obsolete
    else if  (!strcmp(s[0],":AveTemperature"        )){code=12; } //should make obsolete

    else if  (!strcmp(s[0],":Data"                  )){code=8;  }
    else if  (!strcmp(s[0],":MultiData"             )){code=13; }
    else if  (!strcmp(s[0],":GaugeDataTable"        )){code=15; }
    else if  (!strcmp(s[0],":EnsimTimeSeries"       )){code=20; }
    //----------GAUGE-SPECIFIC CORRECTION TERMS / PARAMETERS----
    else if  (!strcmp(s[0],":RainCorrection"        )){code=30; }
    else if  (!strcmp(s[0],":SnowCorrection"        )){code=31; }
    else if  (!strcmp(s[0],":CloudTempRanges"       )){code=32; }
    //-------------------OBSERVATIONS---------------------------
    else if  (!strcmp(s[0],":ObservationData"       )){code=40; }
    else if  (!strcmp(s[0],":IrregularObservations" )){code=41; }
    else if  (!strcmp(s[0],":ObservationWeights"    )){code=42; }
    else if  (!strcmp(s[0],":IrregularWeights"      )){code=43; }
    //----------HYDROGRAPHS /RESERVOIR PARAMS -------------------
    else if  (!strcmp(s[0],":BasinInflowHydrograph" )){code=50; }
    else if  (!strcmp(s[0],":ReservoirExtraction"   )){code=51; }
    else if  (!strcmp(s[0],":VariableWeirHeight"    )){code=52; }
    else if  (!strcmp(s[0],":ReservoirMaxStage"     )){code=53; }
    else if  (!strcmp(s[0],":OverrideReservoirFlow" )){code=54; }
    else if  (!strcmp(s[0],":ReservoirMinStage"     )){code=55; }
    else if  (!strcmp(s[0],":ReservoirMinStageFlow" )){code=56; }
    else if  (!strcmp(s[0],":ReservoirTargetStage"  )){code=57; }
    else if  (!strcmp(s[0],":ReservoirMaxQDelta"    )){code=58; }
    else if  (!strcmp(s[0],":BasinInflowHydrograph2")){code=59; }
    else if  (!strcmp(s[0],":ReservoirMinFlow"      )){code=60; }
    else if  (!strcmp(s[0],":ReservoirDownstreamFlow")){code=61;}
    else if  (!strcmp(s[0],":ReservoirMaxQDecrease" )){code=62; }
    else if  (!strcmp(s[0],":IrrigationDemand"      )){code=63; }
    else if  (!strcmp(s[0],":ReservoirDownstreamDemand")){code=64; }
    else if  (!strcmp(s[0],":ReservoirMaxFlow"      )){code=65;}
    else if  (!strcmp(s[0],":FlowDiversion"         )){code=66;}
    else if  (!strcmp(s[0],":FlowDiversionLookupTable")){code=67;}    
    else if  (!strcmp(s[0],":DZTRReservoirModel"    )){code=68;}
    //--------------------Other --------------------------------
    else if  (!strcmp(s[0],":MonthlyAveTemperature" )){code=70; }
    else if  (!strcmp(s[0],":MonthlyAveEvaporation" )){code=71; }
    else if  (!strcmp(s[0],":MonthlyEvapFactor"     )){code=71; }
    else if  (!strcmp(s[0],":MonthlyMinTemperature" )){code=72; }
    else if  (!strcmp(s[0],":MonthlyMaxTemperature" )){code=73; }
    else if  (!strcmp(s[0],":MonthlyEvapFactors"    )){code=74; }
    else if  (!strcmp(s[0],":MonthlyAveEvaporations")){code=74; }
    else if  (!strcmp(s[0],":MonthlyMaxTemperatures")){code=75; }
    else if  (!strcmp(s[0],":MonthlyMinTemperatures")){code=76; }
    else if  (!strcmp(s[0],":MonthlyAveTemperatures")){code=77; }
    //-----------------CONTROLS ---------------------------------
    else if  (!strcmp(s[0],":OverrideStreamflow"    )){code=100; }
    //-----------------TRANSPORT--------------------------------
    else if  (!strcmp(s[0],":ConcentrationTimeSeries")){code=300; }
    else if  (!strcmp(s[0],":MassFluxTimeSeries"    )){code=300; }
    //---------GRIDDED INPUT (lat,lon,time)---------------------
    else if  (!strcmp(s[0],":GriddedForcing"          )){code=400;}
    else if  (!strcmp(s[0],":ForcingType"             )){code=401;}
    else if  (!strcmp(s[0],":FileNameNC"              )){code=402;}
    else if  (!strcmp(s[0],":VarNameNC"               )){code=403;}
    else if  (!strcmp(s[0],":DimNamesNC"              )){code=404;}
    else if  (!strcmp(s[0],":GridWeights"             )){code=405;}
    else if  (!strcmp(s[0],":EndGriddedForcing"       )){code=406;}
    else if  (!strcmp(s[0],":Deaccumulate"            )){code=407;}
    else if  (!strcmp(s[0],":TimeShift"               )){code=408;}
    else if  (!strcmp(s[0],":LinearTransform"         )){code=409;}
    else if  (!strcmp(s[0],":PeriodEndingNC"          )){code=410;}
    else if  (!strcmp(s[0],":LatitudeVarNameNC"       )){code=411;}
    else if  (!strcmp(s[0],":LongitudeVarNameNC"      )){code=412;}
    else if  (!strcmp(s[0],":ElevationVarNameNC"      )){code=413;}
    else if  (!strcmp(s[0],":GridWeightsByAttribute"  )){code=414;}
    else if  (!strcmp(s[0],":StationElevationsByAttribute")){code=415;}
    else if  (!strcmp(s[0],":StationIDNameNC"         )){code=416;}
    else if  (!strcmp(s[0],":StationElevationsByIdx"  )){code=417;}
    //---------STATION DATA INPUT AS NETCDF (stations,time)------
    //             code 401-405 & 407-409 are shared between
    //             GriddedForcing and StationForcing
    else if  (!strcmp(s[0],":StationForcing"          )){code=500;}
    else if  (!strcmp(s[0],":EndStationForcing"       )){code=506;}

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
        ostrstream FILENAME;
        FILENAME<<":RedirectToFile: Cannot find file "<<filename<<ends;
        ExitGracefully(FILENAME.str() ,BAD_DATA);
      }
      else{
        pMainParser=p;    //save pointer to primary parser
        p=new CParser(INPUT2,filename,line);//open new parser
      }
      break;
    }
    case(-4):  //----------------------------------------------
    {/*:End*/
      if (Options.noisy) {cout <<"EOF"<<endl;} ended=true; break;
    }
    case(1):  //----------------------------------------------
    {/*:Gauge [optional name]*/
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
    {/*":EndGauge"*/
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
    case(3):  //----------------------------------------------
    {/*":Rain"
       {int nMeasurements double start_day int start_year double tstep}
       {double rain} x nMeasurements [mm/d]
       :EndRain
     */
      if (Options.noisy) {cout <<"Rainfall data"<<endl;}
      ExitGracefullyIf(pGage==NULL,
                       "ParseTimeSeriesFile::Rainfall data added before specifying a gauge station and its properties",BAD_DATA);
      pTimeSer=CTimeSeries::Parse(p,true,"RAINFALL","",Options);
      pGage->AddTimeSeries(pTimeSer,F_RAINFALL);
      break;
    }
    case(4):  //----------------------------------------------
    {/*:Snow
       {int nMeasurements double start_day int start_year double tstep}
       {double snow} x nMeasurements [mm/d SWE]
       :EndSnow
     */
      if (Options.noisy) {cout <<"Snow data"<<endl;}
      ExitGracefullyIf(pGage==NULL,
                       "ParseTimeSeriesFile::Snow data added before specifying a gauge station and its properties",BAD_DATA);
      pTimeSer=CTimeSeries::Parse(p,true,"SNOWFALL","",Options);
      pGage->AddTimeSeries(pTimeSer,F_SNOWFALL);
      break;
    }
    case(5):  //----------------------------------------------
    {/*:MinTemperature
       {int nMeasurements double start_day int start_year double tstep}
       {double snow} x nMeasurements [C]
       :EndMinTemperature
     */
      if (Options.noisy) {cout <<"Minumum Temperature data"<<endl;}
      ExitGracefullyIf(pGage==NULL,
                       "ParseTimeSeriesFile::Temperature data added before specifying a gauge station and its properties",BAD_DATA);
      pTimeSer=CTimeSeries::Parse(p,true,"TEMP_DAILY_MIN","",Options);
      pGage->AddTimeSeries(pTimeSer,F_TEMP_DAILY_MIN);
      break;
    }
    case(6):  //----------------------------------------------
    {/*:MaxTemperature
       {int nMeasurements double start_day int start_year double tstep}
       {double snow} x nMeasurements [C]
       :EndMaxTemperature
     */
      if (Options.noisy) {cout <<"Maximum Temperature data"<<endl;}
      ExitGracefullyIf(pGage==NULL,
                       "ParseTimeSeriesFile::Temperature data added before specifying a gauge station and its properties",BAD_DATA);
      pTimeSer=CTimeSeries::Parse(p,true,"TEMP_DAILY_MAX","",Options);
      pGage->AddTimeSeries(pTimeSer,F_TEMP_DAILY_MAX);
      break;
    }
    case(7):  //----------------------------------------------
    {/*:TotalPrecip
       {int nMeasurements double start_day int start_year double tstep}
       {double rain} x nMeasurements [mm/d]
       :EndTotalPrecip
     */
      if (Options.noisy) {cout <<"Total Precipitation data"<<endl;}
      ExitGracefullyIf(pGage==NULL,
                       "ParseTimeSeriesFile::Precipitation data added outside of a :Gage-:EndGauge statement",BAD_DATA);
      pTimeSer=CTimeSeries::Parse(p,true,"PRECIP","",Options);
      pGage->AddTimeSeries(pTimeSer,F_PRECIP);
      //pGage->DeleteSnowData();//??
      break;
    }
    case(8):  //----------------------------------------------
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
      if (Options.noisy) {cout <<"Time Series Data: "<<s[1]<<endl;}
      string name=to_string(s[1]);
      ExitGracefullyIf(pGage==NULL,
                       "ParseTimeSeriesFile::Time Series Data added outside of a :Gage-:EndGauge statement",BAD_DATA);
      pTimeSer=CTimeSeries::Parse(p,true,name,"",Options);
      forcing_type ftype=GetForcingTypeFromString(name);
      if(ftype==F_UNRECOGNIZED){ ExitGracefully("Unrecognized forcing type string in :Data command",BAD_DATA_WARN); }
      pGage->AddTimeSeries(pTimeSer,ftype);
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
    {/*":AveTemperature"
       {int nMeasurements double start_day int start_year double tstep}
       {double temp} x nMeasurements [mm/d]
       :EndAveTemperature
     */
      if (Options.noisy) {cout <<"Ave. Temperature data"<<endl;}
      ExitGracefullyIf(pGage==NULL,
                       "ParseTimeSeriesFile::Temperature data added before specifying a gauge station and its properties",BAD_DATA);
      pTimeSer=CTimeSeries::Parse(p,true,"TEMP_DAILY_AVE","",Options);
      pGage->AddTimeSeries(pTimeSer,F_TEMP_DAILY_AVE);
      break;
    }
    case(13):  //----------------------------------------------
    {/*:MultiData
       {int nMeasurements double start_day int start_year double tstep} or
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
      for (int i=0; i<nSeries; i++)
      {
        pGage->AddTimeSeries(pTimeSerArray[i],aTypes[i]);
      }
      break;
    }
    case(14):  //----------------------------------------------
    {/*":GaugeList"
       :Attributes,  LATITUDE, LONGITUDE, ELEVATION, RAINCORRECTION,...
       :Units,   dec.deg, dec.deg, masl,none,...
       [name1, lat1, long1,elev1,...]
       ...
       [nameN, latN, longN,elevN,...]
       :EndGaugeList
     */
      if (Options.noisy) {cout <<"Gauge List"<<endl;}
      bool done=false;
      string aAttStrings[MAXINPUTITEMS];
      int nAttStrings=0;
      while (!done)
      {
        p->Tokenize(s,Len);
        if      (IsComment(s[0], Len)){}//comment line
        else if (!strcmp(s[0],":Attributes")){
          for (int i=1;i<Len;i++){
            aAttStrings[i-1]=s[i];
          }
          nAttStrings=Len-1;
        }
        else if (!strcmp(s[0],":Units")){
          //Do nothing with units for now
          done=true;
        }
        else {WriteWarning("Improper format of :GaugeList command",Options.noisy); break;}
      }

      double *values[MAX_GAUGES_IN_LIST];
      string gaugenames[MAX_GAUGES_IN_LIST];
      for (int g=0;g<MAX_GAUGES_IN_LIST;g++){values[g]=new double [nAttStrings];}
      p->Tokenize(s,Len);
      done=false;
      int nGauges=0;
      while (!done)
      {
        if      (IsComment(s[0],Len)){}
        else if (Len==nAttStrings+1)
        {
          gaugenames[nGauges]=s[0];
          for (int j=1;j<nAttStrings+1;j++){
            values[nGauges][j-1]=s_to_d(s[j]);
          }
          nGauges++;
        }
        else{
          p->ImproperFormat(s); break;
        }
        p->Tokenize(s,Len);
        if (!strcmp(s[0],":EndGaugeList")){done=true;}
      }
      //Create gauges, add to model
      for (int g=0;g<nGauges;g++)
      {
        CGauge *pGage=new CGauge(gaugenames[g],NOT_SPECIFIED,NOT_SPECIFIED,0);
        for (int j=0;j<nAttStrings;j++){
          pGage->SetProperty(aAttStrings[j],values[g][j]);
        }
        pModel->AddGauge(pGage);
      }
      for (int g=0;g<MAX_GAUGES_IN_LIST;g++){delete [] values[g];}
      break;
    }
    case(15):  //----------------------------------------------
    {/*":GaugeDataTable"
       :DataType PRECIP
       :Units         mm/d
       :StartTime 01-01-2012 00:00:00.0
       :TimeIncrement   01:00:00.0
       :NumMeasurements 730
       :Gauge, Cell_11,Cell_12,Cell_13, ... Cell_240360
       1,  0.0,0.0,0.0,                         ...,0.2
       2,  0.0,0.0,0.0,                         ...,0.1
       ...
       :EndGaugeDataTable
     */
      if (Options.noisy) {cout <<"Gauge Data Table"<<endl;}
      forcing_type datatype=F_UNRECOGNIZED;
      time_struct  tt;
      double   start_day=0;
      int      start_yr=1863;
      double   interval=0;
      int      numvalues;

      long     nMeasurements=0;
      bool     done=false;
      int      NG;
      double **values;
      int     *gauge_ind;

      while (!done)
      {
        p->Tokenize(s,Len);
        if      (IsComment(s[0],Len))           {/* do nothing */}
        else if (!strcmp(s[0],":Units"        )){/* do nothing */}
        else if (!strcmp(s[0],":DataType"     )){
          datatype=GetForcingTypeFromString(s[1]);
          ExitGracefullyIf(datatype==F_UNRECOGNIZED,
                           "ParseTimeSeriesFile:GaugeDataTable: invalid data type specifier",BAD_DATA);
        }
        else if (!strcmp(s[0],":NumMeasurements")){nMeasurements=s_to_l(s[1]);}
        else if (!strcmp(s[0],":TimeIncrement")){interval=s_to_d(s[1]);}
        else if (!strcmp(s[0],":StartTime"    )){
          tt=DateStringToTimeStruct(string(s[1]),string(s[2]),Options.calendar);
          start_day=tt.julian_day;
          start_yr =tt.year;
        }
        else if (!strcmp(s[0],":Gauge"))
        {
          ExitGracefullyIf(nMeasurements<=0,
                           "ParseTimeSeriesFile:GaugeDataTable: :NumMeasurements not specified or invalid",BAD_DATA);
          NG=Len-1;
          values=new double *[NG];
          gauge_ind=new int [NG];
          for (int g=0;g<NG;g++){
            values[g]=new double [nMeasurements];
          }
          numvalues=0;
          for (int g=0;g<NG;g++){
            gauge_ind[g]=pModel->GetGaugeIndexFromName(s[g+1]);
            ExitGracefullyIf(gauge_ind[g]==DOESNT_EXIST,
                             "ParseTimeSeriesFile:GaugeDataTable: unrecognized gauge name",BAD_DATA);
          }
          while (!done)
          {
            p->Tokenize(s,Len);
            if (IsComment(s[0],Len)){}
            else if (!strcmp(s[0],":EndGaugeDataTable")){done=true;}
            else if (Len==NG+1){
              ExitGracefullyIf(numvalues>=nMeasurements,
                               "ParseTimeSeriesFile:GaugeDataTable: more measurement entries than specified by :NumMeasurements",BAD_DATA);
              for (int g=0;g<NG;g++){
                values[g][numvalues]=s_to_d(s[g+1]);
              }
              numvalues++;
            }
            else{
              ExitGracefully("ParseTimeSeriesFile:GaugeDataTable: incorrect number of columns in table",BAD_DATA);
            }
          }
          ExitGracefullyIf(interval<=0,"ParseTimeSeriesFile:GaugeDataTable: invalid interval",BAD_DATA);
          ExitGracefullyIf(datatype==F_UNRECOGNIZED,"ParseTimeSeriesFile:GaugeDataTable: unrecognized datatype",BAD_DATA);
          //done reading table, now populate
          for (int g=0;g<NG;g++)
          {
            CTimeSeries *pTS=new CTimeSeries(ForcingToString(datatype),"",p->GetFilename(),start_day,start_yr,interval,values[g],numvalues,true);
            pModel->GetGauge(gauge_ind[g])->AddTimeSeries(pTS,datatype);
            delete [] values[g];
          }
          delete [] values;
          delete [] gauge_ind;
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
    {/*RainCorrection
       {":RainCorrection" double corr}
     */
      if (Options.noisy) {cout <<"Rainfall Correction"<<endl;}
      ExitGracefullyIf(pGage==NULL,
                       "ParseTimeSeriesFile::Precipitation correction added before specifying a gauge station and its properties",BAD_DATA);
      if (Len<2){p->ImproperFormat(s); break;}
      pGage->SetProperty("RAINFALL_CORR",s_to_d(s[1]));
      break;
    }
    case(31):  //----------------------------------------------
    {/*"SnowCorrection"
       {":SnowCorrection" double corr}
     */
      if (Options.noisy) {cout <<"Snowfall Correction"<<endl;}
      ExitGracefullyIf(pGage==NULL,
                       "ParseTimeSeriesFile::Precipitation correction added before specifying a gauge station and its properties",BAD_DATA);
      if (Len<2){p->ImproperFormat(s); break;}
      pGage->SetProperty("SNOWFALL_CORR",s_to_d(s[1]));
      break;
    }
    case(32):  //----------------------------------------------
    {/*":CloudTempRanges"
       {":CloudTempRanges" double cloud_temp_min [C] cloud_temp_max [C]}
     */
      if (Options.noisy) {cout <<"Cloud cover temperature ranges"<<endl;}
      ExitGracefullyIf(pGage==NULL,
                       "ParseTimeSeriesFile::Cloudcover temperature ranges added before specifying a gauge station and its properties",BAD_DATA);
      if (Len<3){p->ImproperFormat(s); break;}
      pGage->SetProperty("CLOUD_MIN_RANGE",s_to_d(s[1]));
      pGage->SetProperty("CLOUD_MAX_RANGE",s_to_d(s[2]));
      break;
    }
    case (40): //---------------------------------------------
    {/*:ObservationData {data type} {long SBID or int HRUID} {(optional) units}
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
      bool invalidSB=(pModel->GetSubBasinByID(s_to_l(s[2]))==NULL);
        
      bool period_ending =(!strcmp(s[1],"HYDROGRAPH")); 
      //Hydrographs are internally stored as period-ending!

      pTimeSer=CTimeSeries::Parse(p,true,to_string(s[1]),to_string(s[2]),Options,period_ending);

      if (ishyd && invalidSB){
        string warn="ParseTimeSeries:: Invalid subbasin ID in observation hydrograph time series ["+pTimeSer->GetSourceFile()+"]. Will be ignored";
        WriteWarning(warn.c_str(),Options.noisy); break;
      }
      if(isstage && invalidSB){
        string warn="ParseTimeSeries:: Invalid subbasin ID in observed reservoir stage time series ["+pTimeSer->GetSourceFile()+"]. Will be ignored";
        WriteWarning(warn.c_str(),Options.noisy); break;
      }
      if(isinflow && invalidSB){
        string warn="ParseTimeSeries:: Invalid subbasin ID in observed reservoir inflow time series ["+pTimeSer->GetSourceFile()+"]. Will be ignored";
        WriteWarning(warn.c_str(),Options.noisy); break;
      }
      if(isnetinflow && invalidSB){
        string warn="ParseTimeSeries:: Invalid subbasin ID in observed reservoir net inflow time series ["+pTimeSer->GetSourceFile()+"]. Will be ignored";
        WriteWarning(warn.c_str(),Options.noisy);  break;
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
      pIrregTimeSers=CIrregularTimeSeries::Parse(p,to_string(s[1]),to_string(s[2]),s_to_i(s[3]));
      pModel->AddObservedTimeSeries(pIrregTimeSers);
      break;
    }
    case (42): //---------------------------------------------
    {/*:ObservationWeights {data type} {long SBID or int HRUID}
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
      bool invalidSB=(pModel->GetSubBasinByID(s_to_l(s[2]))==NULL);

      pTimeSer=CTimeSeries::Parse(p,true,to_string(s[1]),to_string(s[2]),Options);

      if (ishyd && invalidSB){
        string warn="ParseTimeSeries:: Invalid subbasin ID in observation hydrograph weights time series ["+pTimeSer->GetSourceFile()+"]. Will be ignored";
        WriteWarning(warn.c_str(),Options.noisy); break;
      }
      if(isstage && invalidSB){
        string warn="ParseTimeSeries:: Invalid subbasin ID in observed reservoir stage weights time series ["+pTimeSer->GetSourceFile()+"]. Will be ignored";
        WriteWarning(warn.c_str(),Options.noisy); break;
      }
      if(isinflow && invalidSB){
        string warn="ParseTimeSeries:: Invalid subbasin ID in observed reservoir inflow weights time series ["+pTimeSer->GetSourceFile()+"]. Will be ignored";
        WriteWarning(warn.c_str(),Options.noisy); break;
      }
      if(isnetinflow && invalidSB){
        string warn="ParseTimeSeries:: Invalid subbasin ID in observed reservoir net inflow weights time series ["+pTimeSer->GetSourceFile()+"]. Will be ignored";
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
      pIrregTimeSers=CIrregularTimeSeries::Parse(p,to_string(s[1]),to_string(s[2]),s_to_i(s[3]));
      pModel->AddObservedWeightsTS(pIrregTimeSers);
      break;
    }
    case (50): //---------------------------------------------
    {/*:BasinInflowHydrograph {long Basincode}
       {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
       {double Qin} x nMeasurements [m3/s]
       :EndBasinInflowHydrograph
     */
      if (Options.noisy) {cout <<"Basin Inflow Hydrograph"<<endl;}
      long SBID=DOESNT_EXIST;
      CSubBasin *pSB;
      if (Len>=2){SBID=s_to_l(s[1]);}
      pSB=pModel->GetSubBasinByID(SBID);
      pTimeSer=CTimeSeries::Parse(p,false,"Inflow_Hydrograph_"+to_string(SBID),to_string(SBID),Options);
      if (pSB!=NULL){
        pSB->AddInflowHydrograph(pTimeSer);
      }
      else
      {
        string warn;
        warn=":BasinInflowHydrograph Subbasin "+to_string(SBID)+" not in model, cannot set inflow hydrograph";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }              
    case (51): //---------------------------------------------
    {/*:ReservoirExtraction {long SBID}
       {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
       {double Qout} x nMeasurements [m3/s]
       :EndReservoirExtraction
     */
      if (Options.noisy) {cout <<"Reservoir Extraction Time Series"<<endl;}
      long SBID=DOESNT_EXIST;
      CSubBasin *pSB;
      if (Len>=2){SBID=s_to_l(s[1]);}
      pSB=pModel->GetSubBasinByID(SBID);
      pTimeSer=CTimeSeries::Parse(p,true,"Extraction_"+to_string(SBID),to_string(SBID),Options);
      if ((pSB!=NULL) && (pSB->GetReservoir()!=NULL)){
        pSB->GetReservoir()->AddExtractionTimeSeries(pTimeSer);
      }
      else
      {
        string warn;
        warn=":ReservoirExtraction Subbasin "+to_string(SBID)+" not in model or doesn't have reservoir, cannot set Reservoir Extraction";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case (52): //---------------------------------------------
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
      pTimeSer=CTimeSeries::Parse(p,true,"WeirHeight_"+to_string(SBID),to_string(SBID),Options);
      if ((pSB!=NULL) && (pSB->GetReservoir()!=NULL)){
        pSB->GetReservoir()->AddWeirHeightTS(pTimeSer);
      }
      else
      {
        string warn;
        warn=":VariableWeirHeight Subbasin "+to_string(SBID)+" not in model or doesn't have reservoir, cannot set Reservoir weir height";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case (53): //---------------------------------------------
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
      pTimeSer=CTimeSeries::Parse(p,true,"MaxStage_"+to_string(SBID),to_string(SBID),Options);
      if ((pSB!=NULL) && (pSB->GetReservoir()!=NULL)){
        pSB->GetReservoir()->AddMaxStageTimeSeries(pTimeSer);
      }
      else
      {
        string warn;
        warn=":ReservoirMaxStage Subbasin "+to_string(SBID)+" not in model or doesn't have reservoir, cannot set Reservoir maximum stage";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case (54): //---------------------------------------------
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
      pTimeSer=CTimeSeries::Parse(p,true,"ResFlow_"+to_string(SBID),to_string(SBID),Options);
      if ((pSB!=NULL) && (pSB->GetReservoir()!=NULL)){
        pSB->GetReservoir()->AddOverrideQTimeSeries(pTimeSer);
      }
      else
      {
        string warn;
        warn=":OverrideReservoirFlow Subbasin "+to_string(SBID)+" not in model or doesn't have reservoir, cannot set overriden discharge";
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
      pTimeSer=CTimeSeries::Parse(p,true,"MinStage_"+to_string(SBID),to_string(SBID),Options);
      if ((pSB!=NULL) && (pSB->GetReservoir()!=NULL)){
        pSB->GetReservoir()->AddMinStageTimeSeries(pTimeSer);
      }
      else
      {
        string warn;
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
      pTimeSer=CTimeSeries::Parse(p,true,"MinStageFlow_"+to_string(SBID),to_string(SBID),Options);
      if ((pSB!=NULL) && (pSB->GetReservoir()!=NULL)){
        pSB->GetReservoir()->AddMinStageFlowTimeSeries(pTimeSer);
      }
      else
      {
        string warn;
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
      pTimeSer=CTimeSeries::Parse(p,true,"TargetStage_"+to_string(SBID),to_string(SBID),Options);
      if((pSB!=NULL) && (pSB->GetReservoir()!=NULL)) {
        pSB->GetReservoir()->AddTargetStageTimeSeries(pTimeSer);
      }
      else
      {
        string warn;
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
      pTimeSer=CTimeSeries::Parse(p,true,"QDelta_"+to_string(SBID),to_string(SBID),Options);
      if ((pSB!=NULL) && (pSB->GetReservoir()!=NULL)){
        pSB->GetReservoir()->AddMaxQIncreaseTimeSeries(pTimeSer);
      }
      else
      {
        string warn;
        warn=":ReservoirMaxQDelta Subbasin "+to_string(SBID)+" not in model or doesn't have reservoir, cannot set maximum Qdelta";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }  
    case (59): //---------------------------------------------
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
      pTimeSer=CTimeSeries::Parse(p,false,"Inflow_Hydrograph_"+to_string(SBID),to_string(SBID),Options);
      if (pSB!=NULL){
        pSB->AddDownstreamInflow(pTimeSer);
      }
      else
      {
        string warn;
        warn=":BasinInflowHydrograph2 Subbasin "+to_string(SBID)+" not in model, cannot set inflow hydrograph";
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
        pTimeSer=CTimeSeries::Parse(p,true,"Qmin_"+to_string(SBID),to_string(SBID),Options);
        pSB->GetReservoir()->AddMinQTimeSeries(pTimeSer);
      }
      else
      {
        string warn;
        warn=":ReservoirMinFlow Subbasin "+to_string(SBID)+" not in model or doesn't have reservoir, cannot set minimum flow";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case (61): //---------------------------------------------
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
      pTimeSer=CTimeSeries::Parse(p,true,"Qmin_"+to_string(SBID),to_string(SBID),Options);
      if((pSB!=NULL) && (pSBdown!=NULL) && (pSB->GetReservoir()!=NULL)){
        pSB->GetReservoir()->AddDownstreamTargetQ(pTimeSer,pSBdown,range);
      }
      else
      {
        string warn;
        warn=":ReservoirDownstreamFlow Subbasin "+to_string(SBID)+" or downstream subbasin "+to_string(SBID_down)+" not in model or doesn't have reservoir, cannot set downstream flow target";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case (62): //---------------------------------------------
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
        pTimeSer=CTimeSeries::Parse(p,true,"QDeltaDec_"+to_string(SBID),to_string(SBID),Options);
        pSB->GetReservoir()->AddMaxQDecreaseTimeSeries(pTimeSer);
      }
      else
      {
        string warn;
        warn=":ReservoirMaxQDecrease Subbasin "+to_string(SBID)+" not in model or doesn't have reservoir, cannot set maximum Qdelta decrease";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case (63): //---------------------------------------------
    {/*:IrrigationDemand {long Basincode}
     {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
     {double Qin} x nMeasurements [m3/s]
     :EndBasinInflowHydrograph2
     */
      if(Options.noisy) { cout <<"Irrigation/Water use demand"<<endl; }
      long SBID=DOESNT_EXIST;
      CSubBasin *pSB;
      if(Len>=2) { SBID=s_to_l(s[1]); }
      pSB=pModel->GetSubBasinByID(SBID);
      pTimeSer=CTimeSeries::Parse(p,false,"IrrigationDemand_"+to_string(SBID),to_string(SBID),Options);
      if(pSB!=NULL) {
        pSB->AddIrrigationDemand(pTimeSer);
      }
      else
      {
        string warn;
        warn=":IrrigationDemand: Subbasin "+to_string(SBID)+" not in model, cannot set irrigation demand time series";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case (64): //---------------------------------------------
    {/*:ReservoirDownstreamDemand {long downstream SBID} {long reservoir SBID} {double percent_met}*/
      if(Options.noisy) { cout <<"Reservoir downstream demand"<<endl; }
      long SBID=DOESNT_EXIST;
      long SBIDres=DOESNT_EXIST;
      CSubBasin *pSB,*pSBres;
      double mult=CGlobalParams::GetParameter("RESERVOIR_DEMAND_MULT");
      if(Len>=4) { 
        SBID=s_to_l(s[1]); 
        double tmp=AutoOrDouble(s[2]);
        if(tmp==AUTO_COMPUTE) { SBIDres=AUTO_COMPUTE_LONG; }
        else                  { SBIDres=s_to_l(s[2]); }
      }
      pSB=pModel->GetSubBasinByID(SBID);
      if(pSB==NULL) {
        string warn;
        warn=":ReservoirDownstreamDemand: Subbasin "+to_string(SBID)+" not in model, cannot set reservoir downstream demand";
        WriteWarning(warn,Options.noisy);
      }
      else 
      {
        if(SBIDres==AUTO_COMPUTE_LONG)
        {
          if (Options.res_demand_alloc==DEMANDBY_CONTRIB_AREA)
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
                  area=pUpstr[p]->GetDrainageArea();
                  pSBres=pModel->GetSubBasinByID(pUpstr[p]->GetID());
                  pSBres->AddReservoirDownstrDemand(pSB,area/Atot*mult);
                }
              }
            }
          }
          //else if(Options.res_demand_alloc==DEMANDBY_SURFACE_AREA) {

          //}
          else if(Options.res_demand_alloc==DEMANDBY_MAX_CAPACITY) {
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
                  volume=pUpstr[p]->GetReservoir()->GetMaxCapacity();
                  pSBres=pModel->GetSubBasinByID(pUpstr[p]->GetID());
                  pSBres->AddReservoirDownstrDemand(pSB,volume/Vtot*mult);
                }
              }
            }
          }
        }
        else 
        { //single connection
          pSBres=pModel->GetSubBasinByID(SBIDres);
          if(pSBres!=NULL) {
            pSBres->AddReservoirDownstrDemand(pSB,s_to_d(s[3])*mult);
          }
          else
          {
            string warn;
            warn=":ReservoirDownstreamDemand: Reservoir subbasin "+to_string(SBID)+" not in model, cannot set reservoir downstream demand";
            WriteWarning(warn,Options.noisy);
          }
        }
      }
      break;
    }
    case (65): //---------------------------------------------
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
        pTimeSer=CTimeSeries::Parse(p,true,"Qmax_"+to_string(SBID),to_string(SBID),Options);
        pSB->GetReservoir()->AddMaxQTimeSeries(pTimeSer);
      }
      else
      {
        string warn;
        warn=":ReservoirMaxFlow Subbasin "+to_string(SBID)+" not in model or doesn't have reservoir, cannot set maximum flow";
        WriteWarning(warn,Options.noisy);
      }
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
          pSB->AddFlowDiversion(start,end,target_p,s_to_d(s[4]),s_to_d(s[3]));
        }
        else {
          string warn;
          warn=":FlowDiversion command: Target subbasin "+to_string(s_to_l(s[2]))+" not found in model, cannot add diversion";
          WriteWarning(warn,Options.noisy);
        }
      }
      else
      {
        string warn;
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
      int target_p,start, end;
      int NQ;
      CSubBasin *pSB;
      start=0;   //default - all year
      end  =366; //default - all year

      if(Len>=2) { SBID=s_to_l(s[1]); }
      pSB=pModel->GetSubBasinByID(SBID);

      if(pSB!=NULL) {
        if(Len>=5) { start=s_to_i(s[3]); end=s_to_i(s[4]); }
        target_p=pModel->GetSubBasinIndex(s_to_l(s[2]));
      }
      else
      {
        string warn;
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

      if((s_to_l(s[2])==-1) || (target_p!=DOESNT_EXIST)) {
        pSB->AddFlowDiversion(start,end,target_p,aQ1,aQ2,NQ);
      }
      else {
        string warn;
        warn=":FlowDiversionLookupTable command: Target subbasin "+to_string(s_to_l(s[2]))+" not in model, cannot add diversion";
        WriteWarning(warn,Options.noisy);
      }
      delete [] aQ1; delete [] aQ2;
      break;
    }
    case (68): //---------------------------------------------
    {/*:DZTRResservoirModel [SBID]
        :MaximumStorage           Vmax(m3)  
        :MaximumChannelCapacity   Qmax(m3/s)
        :MonthlyMaxStorage        J F M A M J J A S N D  {or single constant value} (m3)
        :MonthlyNormalStorage     J F M A M J J A S N D  {or constant value} (m3)
        :MonthlyCriticalStorage   J F M A M J J A S N D  {or constant value} (m3)
        :MonthlyMaxDischarge      J F M A M J J A S N D  {or constant value} (m3/s)
        :MonthlyNormalDischarge   J F M A M J J A S N D  {or constant value} (m3/s)
        :MonthlyCriticalDischarge J F M A M J J A S N D  {or constant value} (m3/s)
      :EndDZTRReservoirModel
     */
      //Dynamically zoned target release model
      if(Options.noisy) { cout <<"DZTR Reservoir Model"<<endl; }
      double Smi[12],Sni[12],Sci[12];
      double Qmi[12],Qni[12],Qci[12];
      double Qmax,Smax;
      CSubBasin *pSB;
      long SBID;

      if(Len>=2) { SBID=s_to_l(s[1]); }
      pSB=pModel->GetSubBasinByID(SBID);
      ExitGracefullyIf(pSB==NULL,":DZTRResservoirModel: invalid subbasin ID",BAD_DATA);

      p->Tokenize(s,Len);
      Smax=s_to_d(s[1]);

      p->Tokenize(s,Len);
      Qmax=s_to_d(s[1]);

      p->Tokenize(s,Len);
      if     (Len>=13){ for(int i=0;i<12;i++){ Smi[i]=s_to_d(s[i+1]); }}
      else if (Len>=2){ for(int i=0;i<12;i++){ Smi[i]=s_to_d(s[1]  ); }}

      p->Tokenize(s,Len);
      if     (Len>=13){ for(int i=0;i<12;i++){ Sni[i]=s_to_d(s[i+1]); }}
      else if (Len>=2){ for(int i=0;i<12;i++){ Sni[i]=s_to_d(s[1]  ); }}

      p->Tokenize(s,Len);
      if     (Len>=13){ for(int i=0;i<12;i++){ Sci[i]=s_to_d(s[i+1]); }}
      else if (Len>=2){ for(int i=0;i<12;i++){ Sci[i]=s_to_d(s[1]  ); }}

      p->Tokenize(s,Len);
      if     (Len>=13){ for(int i=0;i<12;i++){ Qmi[i]=s_to_d(s[i+1]); }}
      else if (Len>=2){ for(int i=0;i<12;i++){ Qmi[i]=s_to_d(s[1]  ); }}

      p->Tokenize(s,Len);
      if     (Len>=13){ for(int i=0;i<12;i++){ Qni[i]=s_to_d(s[i+1]); }}
      else if (Len>=2){ for(int i=0;i<12;i++){ Qni[i]=s_to_d(s[1]  ); }}

      p->Tokenize(s,Len);
      if     (Len>=13){ for(int i=0;i<12;i++){ Qci[i]=s_to_d(s[i+1]); }}
      else if (Len>=2){ for(int i=0;i<12;i++){ Qci[i]=s_to_d(s[1]  ); }}

      if(pSB->GetReservoir()!=NULL) {
        pSB->GetReservoir()->SetDZTRModel(Qmax,Smax,Sni,Sci,Smi,Qni,Qci,Qmi);
      }
      else {
        ExitGracefully("ParseTimeSeriesFile: subbasin in :DZTRResservoirModel command does not have reservoir",BAD_DATA_WARN);
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
    {/*:MonthlyAveEvaporations
       gaugename1, J,F,M,A...O,N,D
       gaugename2, J,F,M,A...O,N,D
       ...
       :EndMonthlyAveEvaporations
     */
      if (Options.noisy) {cout <<"Monthly Average Evaporation"<<endl;}
      bool done=false;
      int g;
      while (!done)
      {
        p->Tokenize(s,Len);
        if (IsComment(s[0],Len)){}
        else if (!strncmp(s[0],":End",4)){done=true;}
        else
        {
          g=pModel->GetGaugeIndexFromName(s[0]);
          ExitGracefullyIf(g==DOESNT_EXIST,
                           "ParseTimeSeriesFile:: :MonthlyAveEvaporations: invalid gauge name",BAD_DATA);
          if (Len < 13){WriteWarning(":MonthlyAveEvaporations: bad number of entries",Options.noisy); break;}
          for (int i=0;i<12;i++){monthdata[i]=s_to_d(s[i+1]);}
          pModel->GetGauge(g)->SetMonthlyPET(monthdata);
        }
      }
      break;
    }
    case (75): //---------------------------------------------
    {/*:MonthlyMaxTemperatures
       gaugename1, J,F,M,A...O,N,D
       gaugename2, J,F,M,A...O,N,D
       ...
       :EndMonthlyMaxTemperatures
     */
      if (Options.noisy) {cout <<"Monthly Maximum temperatures"<<endl;}
      bool done=false;
      int g;
      while (!done)
      {
        p->Tokenize(s,Len);
        if (IsComment(s[0],Len)){}
        else if (!strncmp(s[0],":End",4)){done=true;}
        else
        {
          g=pModel->GetGaugeIndexFromName(s[0]);
          ExitGracefullyIf(g==DOESNT_EXIST,
                           "ParseTimeSeriesFile:: :MonthlyMaxTemperatures: invalid gauge name",BAD_DATA);
          if (Len < 13){WriteWarning(":MonthlyMaxTemperatures: bad number of entries",Options.noisy); break;}
          for (int i=0;i<12;i++){monthdata[i]=s_to_d(s[i+1]);}
          pModel->GetGauge(g)->SetMonthlyMaxTemps(monthdata);
        }
      }
      break;
    }
    case (76): //---------------------------------------------
    {/*:MonthlyMinTemperatures
       gaugename1, J,F,M,A...O,N,D
       gaugename2, J,F,M,A...O,N,D
       ...
       :EndMonthlyMinTemperatures
     */
      if (Options.noisy) {cout <<"Monthly Minimum temperatures"<<endl;}
      bool done=false;
      int g;
      while (!done)
      {
        p->Tokenize(s,Len);
        if (IsComment(s[0],Len)){}
        else if (!strncmp(s[0],":End",4)){done=true;}
        else
        {
          g=pModel->GetGaugeIndexFromName(s[0]);
          ExitGracefullyIf(g==DOESNT_EXIST,
                           "ParseTimeSeriesFile:: :MonthlyMinTemperatures: invalid gauge name",BAD_DATA);
          if (Len < 13){WriteWarning(":MonthlyMinTemperatures: bad number of entries",Options.noisy); break;}
          for (int i=0;i<12;i++){monthdata[i]=s_to_d(s[i+1]);}
          pModel->GetGauge(g)->SetMonthlyMinTemps(monthdata);
        }
      }
      break;
    }
    case (77): //---------------------------------------------
    {/*:MonthlyAveTemperatures
       gaugename1, J,F,M,A...O,N,D
       gaugename2, J,F,M,A...O,N,D
       ...
       :EndMonthlyAveTemperatures
     */
      if (Options.noisy) {cout <<"Monthly Average temperatures"<<endl;}
      bool done=false;
      int g;
      while (!done)
      {
        p->Tokenize(s,Len);
        if (IsComment(s[0],Len)){}
        else if (!strncmp(s[0],":End",4)){done=true;}
        else
        {
          g=pModel->GetGaugeIndexFromName(s[0]);
          ExitGracefullyIf(g==DOESNT_EXIST,
                           "ParseTimeSeriesFile:: :MonthlyAveTemperatures: invalid gauge name",BAD_DATA);
          if (Len < 13){WriteWarning(":MonthlyAveTemperatures: bad number of entries",Options.noisy); break;}
          for (int i=0;i<12;i++){monthdata[i]=s_to_d(s[i+1]);}
          pModel->GetGauge(g)->SetMonthlyAveTemps(monthdata);
        }
      }
      break;
    }
    case(100):
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
    case (300)://----------------------------------------------
    {/*:ConcentrationTimeSeries
       :ConcentrationTimeSeries [string constit_name] [string state_var (storage compartment)] {optional HRU Group name}
       {yyyy-mm-dd hh:mm:ss double tstep int nMeasurements}
       {double concentration values} x nMeasurements
       :EndConcentrationTimeSeries
     */
      if (Options.noisy){cout <<"Fixed concentration time series"<<endl;}

      int layer_ind;
      int i_stor;
      string const_name = to_string(s[1]);
      sv_type typ=CStateVariable::StringToSVType(s[2],layer_ind,false);
      if (typ==UNRECOGNIZED_SVTYPE){
        WriteWarning(":ConcentrationTimeSeries command: unrecognized storage variable name: "+to_string(s[2]),Options.noisy);
        break;
      }
      i_stor=pModel->GetStateVarIndex(typ,layer_ind);
      if (i_stor!=DOESNT_EXIST){
        int kk=DOESNT_EXIST;
        if (Len>3){
          CHRUGroup *pSourceGrp;
          pSourceGrp=pModel->GetHRUGroup(s[3]);
          if (pSourceGrp==NULL){
            ExitGracefully("Invalid HRU Group name supplied in :ConcentrationTimeSeries command in .rvt file",BAD_DATA_WARN);
            break;
          }
          else{
            kk=pSourceGrp->GetGlobalIndex();
          }
        }
        CTimeSeries *pTS;
        pTS=CTimeSeries::Parse(p,true,const_name+"_"+to_string(s[2]),"",Options);//name=constitutent name
        pModel->GetTransportModel()->AddDirichletTimeSeries(const_name, i_stor, kk,pTS);
      }
      else{
        ExitGracefully(":ConcentrationTimeSeries command: invalid state variable name",BAD_DATA_WARN);
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
      sv_type typ=CStateVariable::StringToSVType(s[2],layer_ind,false);
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
        CTimeSeries *pTS;
        pTS=CTimeSeries::Parse(p,true,const_name+"_"+to_string(s[2]),"",Options);//name=constitutent name
        pModel->GetTransportModel()->AddInfluxTimeSeries(const_name, i_stor, kk,pTS);
      }
      else{
        ExitGracefully(":MassFluxTimeSeries command: invalid state variable name",BAD_DATA_WARN);
      }
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
      if(!grid_initialized) { pGrid->ForcingGridInit(Options);grid_initialized = true; }
      pModel->AddForcingGrid(pGrid,pGrid->GetForcingType());
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
      
      //check for file existence
      ifstream TESTNETCDF;
      TESTNETCDF.open(filename.c_str());
      if(TESTNETCDF.fail()){
        string warn = "ParseTimeSeriesFile: :FileNameNC command: Cannot find gridded data file "+ filename; 
        ExitGracefully(warn.c_str(),BAD_DATA_WARN);
        break;
      }
      TESTNETCDF.close();

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
       #     and grid cell index l is derived by l = (i-1) * NC + j
       #     where i and j are the row and column of cell l respectively and
       #     NC is the total number of columns.
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

          ExitGracefullyIf(pModel->GetNumHRUs() != nHydroUnits,
                           "ParseTimeSeriesFile: :NumberHRUs given does not agree with HRUs given in *.rvh file",BAD_DATA);

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
            pHRU = pModel->GetHRUByID(atoi(s[0]));
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
      bool WeightArrayOK = pGrid->CheckWeightArray(nHydroUnits,nGridCells);
      ExitGracefullyIf(!WeightArrayOK,
                       "ParseTimeSeriesFile: Check of weights for gridded forcing failed. Sum per HRUID must be 1.0.",BAD_DATA);

      // store (sorted) grid cell ids with non-zero weight in array
      pGrid->SetIdxNonZeroGridCells(nHydroUnits,nGridCells);
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
      ExitGracefullyIf(pGrid==NULL,"ParseTimeSeriesFile: :StationIDNameNC command must be within a StationForcing block",BAD_DATA);
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

          ExitGracefullyIf(pModel->GetNumHRUs() != nHydroUnits,
            "ParseTimeSeriesFile: :NumberHRUs given does not agree with HRUs given in *.rvh file",BAD_DATA);

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
            pHRU = pModel->GetHRUByID(atoi(s[0]));
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
      bool WeightArrayOK = pGrid->CheckWeightArray(nHydroUnits,nGridCells);
      ExitGracefullyIf(!WeightArrayOK,
        "ParseTimeSeriesFile: Check of weights for gridded forcing failed. Sum per HRUID must be 1.0.",BAD_DATA);

      // store (sorted) grid cell ids with non-zero weight in array
      pGrid->SetIdxNonZeroGridCells(nHydroUnits,nGridCells);
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
          string warn ="IGNORING unrecognized command: |" + string(s[0])+ "| in .rvt file "+p->GetFilename()+" line: "+ to_string(p->GetLineNumber());
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

  delete p; p=NULL;

  return true;
}


