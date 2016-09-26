/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright © 2008-2014 the Raven Development Team
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

  int   code;
  int   Len,line(0);
  char *s[MAXINPUTITEMS];
  double monthdata[12];
  bool ended(false);

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
    //--------------------HYDROGRAPHS --------------------------
    else if  (!strcmp(s[0],":BasinInflowHydrograph" )){code=50; }
    else if  (!strcmp(s[0],":ReservoirExtraction"   )){code=51; }
    //--------------------Other --------------------------------
    else if  (!strcmp(s[0],":MonthlyAveTemperature" )){code=70; }
    else if  (!strcmp(s[0],":MonthlyAveEvaporation" )){code=71; }
    else if  (!strcmp(s[0],":MonthlyEvapFactor"     )){code=71; }
    else if  (!strcmp(s[0],":MonthlyMinTemperature" )){code=72; }
    else if  (!strcmp(s[0],":MonthlyMaxTemperature" )){code=73; }
    else if  (!strcmp(s[0],":MonthlyEvapFactors"    )){code=74; }
    else if  (!strcmp(s[0],":MonthlyAveEvaporations")){code=74; }
    else if  (!strcmp(s[0],":MonthlyMaxTemperatures" )){code=75; }
    else if  (!strcmp(s[0],":MonthlyMinTemperatures" )){code=76; }
    else if  (!strcmp(s[0],":MonthlyAveTemperatures" )){code=77; }
    //-----------------TRANSPORT--------------------------------
    else if  (!strcmp(s[0],":ConcentrationTimeSeries")){code=300; }
    else if  (!strcmp(s[0],":MassFluxTimeSeries     ")){code=300; }
    
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
        if (Options.noisy) {cout <<"Redirect to file: "<<filename<<endl;} /// \todo [funct]: should look in the same directory as .rvt file if incomplete path provided

        string filedir = GetDirectoryName(Options.rvt_filename); //if a relative path name, e.g., "/path/model.rvt", only returns e.g., "/path"
        if (StringToUppercase(filename).find(StringToUppercase(filedir)) == string::npos){ //checks to see if absolute dir already included in redirect filename
          filename = filedir + "//" + filename;
        }

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
        //if (Options.noisy) {cout <<"  Creating Gauge: "<< gName <<endl;}
        pGage=new CGauge(gName,NOT_SPECIFIED,NOT_SPECIFIED,0.0);
        
        pModel->AddGauge(pGage);

		  //todo [QA/QC]: require check here (or somewhere in the parsing process) to determine if the gauge count is greter than the specified global limit.
		  //		if number is exceeded we need to exit gracefully (WJ).

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
       */
        if (Options.noisy) {cout <<"Time Series Data: "<<s[1]<<endl;}
        ExitGracefullyIf(pGage==NULL,
          "ParseTimeSeriesFile::Time Series Data added outside of a :Gage-:EndGauge statement",BAD_DATA);
        pTimeSer=CTimeSeries::Parse(p,true,s[1],"",Options);
        pGage->AddTimeSeries(pTimeSer,GetForcingTypeFromString(s[1])); 
        break;
      }
      case(9):  //----------------------------------------------
      {/*:Latitude*/
        if (Options.noisy) {cout <<"Gauge latitude"<<endl;}
        ExitGracefullyIf(pGage==NULL,
          "ParseTimeSeriesFile::Latitude specified outside of a :Gauge-:EndGauge statement",BAD_DATA);
        pGage->SetLatitude(s_to_d(s[1]));
        break;
      }
      case(10):  //----------------------------------------------
      {/*:Longitude*/
        if (Options.noisy) {cout <<"Gauge longitude"<<endl;}
        ExitGracefullyIf(pGage==NULL,
          "ParseTimeSeriesFile::Longitude specified outside of a :Gauge-:EndGauge statement",BAD_DATA);
        // \todo [QA/QC] - should check if there is an s[1]! (also in lat/elev/etc in this file
        pGage->SetLongitude(s_to_d(s[1]));
        break;
      }
      case(11):  //----------------------------------------------
      {/*:Elevation*/
        if (Options.noisy) {cout <<"Gauge elevation"<<endl;}
        ExitGracefullyIf(pGage==NULL,
          "ParseTimeSeriesFile::Elevation specified outside of a :Gauge-:EndGauge statement",BAD_DATA);
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

        pTimeSerArray=CTimeSeries::ParseMultiple(p,nSeries,aTypes,true);

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
        
        double *values[MAX_GAUGES];
        string gaugenames[MAX_GAUGES];
        for (int g=0;g<MAX_GAUGES;g++){values[g]=new double [nAttStrings];}
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
        for (int g=0;g<MAX_GAUGES;g++){delete [] values[g];}
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
            tt=DateStringToTimeStruct(string(s[1]),string(s[2]));
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

        string filedir = GetDirectoryName(Options.rvt_filename); //if a relative path name, e.g., "/path/model.rvt", only returns e.g., "/path"
        if (StringToUppercase(filename).find(StringToUppercase(filedir)) == string::npos){ //checks to see if absolute dir already included in redirect filename
          filename = filedir + "//" + filename;
        }

        pTimeSerArray=CTimeSeries::ParseEnsimTb0(filename,nSeries,aTypes);
        
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
          pTimeSer=CTimeSeries::Parse(p,true,to_string(s[1]),to_string(s[2]),Options,(!strcmp(s[1], "HYDROGRAPH"))); // \todo[funct] should likley not be "is pulse" format. May need alternate time series class
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
          pTimeSer=CTimeSeries::Parse(p,true,to_string(s[1]),to_string(s[2]),Options);
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
          long SBID=-1; 
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
            warn="Subbasin "+to_string(SBID)+" not in model, cannot set inflow hydrograph";
            WriteWarning(warn,Options.noisy);
          }
          break;
        }
      case (51): //---------------------------------------------
        {/*:ReservoirExtraction {long Basincode}
          {yyyy-mm-dd} {hh:mm:ss.0} {double timestep} {int nMeasurements}
          {double Qout} x nMeasurements [m3/s]
          :EndReservoirExtraction
          */
          if (Options.noisy) {cout <<"Reservoir Extraction Time History"<<endl;}
          long SBID=-1; 
          CSubBasin *pSB;
          if (Len>=2){SBID=s_to_l(s[1]);}
          pSB=pModel->GetSubBasinByID(SBID);
          pTimeSer=CTimeSeries::Parse(p,true,"Extraction_"+to_string(SBID),to_string(SBID),Options);
          if (pSB!=NULL){
            pSB->AddReservoirExtract(pTimeSer);
          }
          else 
          {
            string warn;
            warn="Subbasin "+to_string(SBID)+" not in model, cannot set Reservoir Extraction";
            WriteWarning(warn,Options.noisy);
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
      default: 
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
            else if(Options.noisy)
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