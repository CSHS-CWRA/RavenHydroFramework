/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2020 the Raven Development Team
  ----------------------------------------------------------------*/

#include "RavenInclude.h"
#include "Properties.h"
#include "Model.h"
#include "StateVariables.h"
#include "Forcings.h"
#include "Precipitation.h"
#include "SoilWaterMovers.h"
#include "SnowMovers.h"
#include "VegetationMovers.h"
#include "GlacierProcesses.h"
#include "Albedo.h"
#include "CropGrowth.h"
#include "DepressionProcesses.h"
#include "OpenWaterEvap.h"
#include "ParseLib.h"
#include "Advection.h"
#include "Convolution.h"
#include "Diagnostics.h"
#include "CustomOutput.h"
#include "Decay.h"
#include "LatAdvection.h"
#include "PrairieSnow.h"
#include "ProcessGroup.h"
#include "HeatConduction.h"

bool ParseMainInputFile        (CModel *&pModel, optStruct &Options);
bool ParseClassPropertiesFile  (CModel *&pModel, const optStruct &Options);
bool ParseHRUPropsFile         (CModel *&pModel, const optStruct &Options);
bool ParseTimeSeriesFile       (CModel *&pModel, const optStruct &Options);
bool ParseInitialConditionsFile(CModel *&pModel, const optStruct &Options);
bool ParseEnsembleFile         (CModel *&pModel, const optStruct &Options);
bool ParseGWPropsFile          (CModel *&pModel, const optStruct &Options);
bool ParseGWGeometryFile       (CModel *&pModel, const optStruct &Options);
bool ParseGWSWExchangeFile     (CModel *&pModel, const optStruct &Options);
bool ParseGWSWOverlapFile      (CModel *&pModel, const optStruct &Options);

int  ParseSVTypeIndex          (string s,  CModel *&pModel);
int *ParseSVTypeArray          (char *string,  CModel *&pModel, int size);
void ImproperFormatWarning     (string command, CParser *p, bool noisy);
void AddProcess                (CModel *pModel, CHydroProcessABC* pMover, CProcessGroup *pProcGroup);
void AddNetCDFAttribute        (optStruct &Options,const string att,const string &val);

//////////////////////////////////////////////////////////////////
/// \brief This method is the primary Raven input routine that parses input files, called by main
///
/// \details This method provides an interface by which the main method can parse input files\n
///   Files used include:\n
///   - \b [modelname].rvi: input file that determines general model settings, nHRUs, timesteps,etc
///   - \b [modelname].rvp: default properties for LULT & soil type 
///   - \b [modelname].rvh: HRU/basin property files
///   - \b [modelname].rvt: time series precip/temp input
///   - \b [modelname].rvc: initial conditions file
///   - \b [modelname].rvg: groundwater properties file
///
/// \param *&pModel [in] The input model object
/// \param &Options [in] Global model options information
/// \return Boolean variable indicating success of parsing
//
bool ParseInputFiles (CModel      *&pModel,
                      optStruct    &Options)
{
  string filename;

  //Main input file (.rvi)
  if (!ParseMainInputFile        (pModel,Options)){
    if(Options.rvi_filename.compare("nomodel.rvi")==0){
      ExitGracefully("A model input file name must be supplied as an argument to the Raven executable.",BAD_DATA);return false;
    }
    ExitGracefully("Cannot find or read .rvi file",BAD_DATA);return false;}

  //Class Property file (.rvp)
  if (!ParseClassPropertiesFile  (pModel,Options)){
    ExitGracefully("Cannot find or read .rvp file",BAD_DATA);return false;}

  //HRU Property file (.rvh)
  if (!ParseHRUPropsFile         (pModel,Options)){
    ExitGracefully("Cannot find or read .rvh file",BAD_DATA);return false;}

  if((Options.modeltype == MODELTYPE_GROUNDWATER) || (Options.modeltype == MODELTYPE_COUPLED))
  {
    //GWMIGRATE - move all GW information to one GW file
    //Groundwater discretization input file
    if(!ParseGWGeometryFile(pModel,Options)) { return false; }
    if(!ParseGWPropsFile(pModel,Options)) { return false; }
  
    if(Options.modeltype == MODELTYPE_COUPLED)
    {
      //Groundwater Surface water interface input files
      if(!ParseGWSWOverlapFile(pModel,Options)) { return false; }

      if(Options.flux_exchange == FLUX_EX_UFR)
      {
        if(!ParseGWSWExchangeFile(pModel,Options)) { return false; }
      }
    }
  }
  
  for (int pp=0;pp<pModel->GetNumSubBasins(); pp++){
    if (pModel->GetSubBasin(pp)->GetReservoir()!=NULL){Options.write_reservoir=true;}
  }

  //Initial Conditions input file (.rvc)
  if (!ParseInitialConditionsFile(pModel,Options)){
    ExitGracefully("Cannot find or read .rvc file",BAD_DATA);return false;}

  CGroundwaterModel *pGW=pModel->GetGroundwaterModel();
  if (pGW!=NULL){pGW->ConvertInitialConditions(pModel);}//GWMIGRATE ****FIX - Shouldn't be here, need to place it somewhere else
    
  //Time series input file (.rvt)
  if (!ParseTimeSeriesFile       (pModel,Options)){
    ExitGracefully("Cannot find or read .rvt file",BAD_DATA);return false;}

  //Ensemble file (.rve)
  if(pModel->GetEnsemble()->GetType()!=ENSEMBLE_NONE) {
    if(!ParseEnsembleFile(pModel,Options)) {
      ExitGracefully("Cannot find or read .rve file",BAD_DATA);return false;
    }
  }
  
  if (!Options.silent){
    cout <<"...model input successfully parsed"             <<endl;
    cout <<endl;
  }

  return true;
}

///////////////////////////////////////////////////////////////////
/// \brief This local method (called by ParseInputFiles) reads an input .rvi file and generates a new model with all options and processes created.
///
/// \remark Custom output must be AFTER processes are specified.
///
/// \param *filename [in] The fully-qualified file name of .rvi input file
/// \param *&pModel [in] Input object that determines general model settings
/// \param &Options [in] Global model options information
/// \return Boolean   value indicating success of parsing
//
bool ParseMainInputFile (CModel     *&pModel,
                         optStruct   &Options)
{
  int i;
  CHydroProcessABC *pMover;
  CmvPrecipitation *pPrecip=NULL;
  CProcessGroup    *pProcGroup=NULL;
  CProcessGroup    *pProcGroupOuter=NULL; //for nested processgroups
  bool              transprepared(false);
  bool              runname_overridden(false);
  int               num_ensemble_members=1;
  ifstream INPUT;

  int      tmpN;
  sv_type *tmpS;
  int     *tmpLev;

  tmpS  =new sv_type[MAX_STATE_VARS];
  tmpLev=new int    [MAX_STATE_VARS];

  if (Options.noisy){
    cout <<"======================================================"<<endl;
    cout << "Parsing Input File " << Options.rvi_filename <<"..."<<endl;
    cout <<"======================================================"<<endl;
  }
 
  int   code;
  bool  ended(false);
  int   Len,line(0);
  char *s[MAXINPUTITEMS];
  INPUT.open(Options.rvi_filename.c_str());
  if (INPUT.fail()){cout << "Cannot find file "<<Options.rvi_filename <<endl; return false;}

  CParser *p=new CParser(INPUT,Options.rvi_filename,line);

  //Default Values---------------------------------------------------
  if(Options.run_name!=""){runname_overridden=true;}
  Options.julian_start_day        =0;//Jan 1
  Options.julian_start_year       =1666;
  Options.duration                =365;
  Options.calendar                =CALENDAR_PROLEPTIC_GREGORIAN; // Default calendar
  Options.timestep                =1;
  Options.output_interval         =1;
  Options.sol_method              =ORDERED_SERIES;
  Options.convergence_crit        =0.01;
  Options.max_iterations          =30;
  Options.ensemble                =ENSEMBLE_NONE; 
  Options.external_script         ="";

  Options.routing                 =ROUTE_STORAGECOEFF;
  Options.catchment_routing       =ROUTE_DUMP;
  Options.distrib_lat_inflow      =false;
  Options.res_demand_alloc        =DEMANDBY_CONTRIB_AREA;

  Options.interpolation           =INTERP_NEAREST_NEIGHBOR;
  Options.interp_file             ="";

  Options.soil_modeltype          =SOIL_ONE_LAYER;
  Options.num_soillayers          =-1;//used to check if SoilModel command is used
  Options.soil_representation     =BROOKS_COREY;

  //Forcing function estimation options:
  Options.evaporation             =PET_HARGREAVES_1985;
  Options.ow_evaporation          =PET_HARGREAVES_1985;
  Options.evap_infill             =PET_HARGREAVES_1985;
  Options.ow_evap_infill          =PET_HARGREAVES_1985;
  Options.orocorr_PET             =OROCORR_NONE;
  Options.orocorr_precip          =OROCORR_NONE;
  Options.orocorr_temp            =OROCORR_NONE;
  Options.LW_radiation            =LW_RAD_DEFAULT;
  Options.LW_incoming             =LW_INC_DEFAULT;
  Options.SW_radiation            =SW_RAD_DEFAULT;
  Options.cloud_cover             =CLOUDCOV_NONE;
  Options.SW_canopycorr           =SW_CANOPY_CORR_NONE;
  Options.SW_cloudcovercorr       =SW_CLOUD_CORR_NONE;
  Options.SW_radia_net            =NETSWRAD_CALC;
  Options.wind_velocity           =WINDVEL_CONSTANT;
  Options.rel_humidity            =RELHUM_CONSTANT;//RELHUM_MINDEWPT; [preferred]
  Options.air_pressure            =AIRPRESS_BASIC;
  Options.rainsnow                =RAINSNOW_DINGMAN;
  Options.month_interp            =MONTHINT_LINEAR_MID;
  Options.pot_melt                =POTMELT_DEGREE_DAY;
  Options.subdaily                =SUBDAILY_NONE;
  Options.interception_factor     =PRECIP_ICEPT_USER;
  Options.recharge                =RECHARGE_NONE;
  Options.direct_evap             =false;
  Options.keepUBCWMbugs           =false;
  Options.suppressCompetitiveET   =false;
  Options.snow_suppressPET        =false;
  Options.pavics                  =false;
  Options.deltaresFEWS            =false;
  Options.res_overflowmode        =OVERFLOW_ALL;
  
  //Groundwater model options
  Options.modeltype           =MODELTYPE_SURFACE;
  Options.gw_solver_outer     =GWSOL_NEWTONRAPHSON;
  Options.gw_solver_inner     =GWSOL_BICGSTAB;
  Options.gw_grid             =GRID_UNSTRUCTURED;
  Options.exchange_freq       =EXFREQ_EVERY_TIMESTEP;
  Options.flux_exchange       =FLUX_EX_STANDARD;
  Options.overlap_type        =DIRECT;
  
  //Output options:
  if (Options.silent!=true){ //if this wasn't overridden in flag to executable
    Options.noisy                 =false;
    Options.silent                =false;
  }
  Options.output_format       =OUTPUT_STANDARD;
  Options.write_energy        =false;
  Options.write_forcings      =false;
  Options.write_mass_bal      =false;
  Options.write_exhaustiveMB  =false;
  Options.write_channels      =false;
  Options.write_watershed_storage =true;
  Options.benchmarking        =false;
  Options.pause               =false;
  Options.debug_mode          =false;
  Options.ave_hydrograph      =true;
  Options.write_reservoir     =false;
  Options.write_reservoirMB   =false;
  Options.write_basinfile     =false;
  Options.write_interp_wts    =false;
  Options.write_demandfile    =false;
  Options.suppressICs         =false;
  Options.period_ending       =false;
  Options.period_starting     =false;//true;
  Options.write_group_mb      =DOESNT_EXIST;
  Options.diag_start_time     =-ALMOST_INF;
  Options.diag_end_time       = ALMOST_INF;
  Options.wateryr_mo          =10; //October
  Options.create_rvp_template =false;
  Options.write_constitmass   =false;
  Options.nNetCDFattribs      =0;
  Options.aNetCDFattribs      =NULL;
  Options.assimilation_on     =false;
  Options.assimilation_start  =0;
  Options.time_zone           =0;
  Options.rvl_read_frequency  =0.0; //do not read at all
  Options.custom_interval     =1.0; //daily
  Options.use_stopfile        =false;

  pModel=NULL;
  pMover=NULL;

  //--Sift through file-----------------------------------------------
  while (!(p->Tokenize(s,Len)))
  {
    if (ended){break;}
    if (Options.noisy){ cout << "reading line " << p->GetLineNumber() << ": ";}
    
    /*assign code for switch statement
      ------------------------------------------------------------------
      <100         : ignored/special
      0   thru 100 : Options
      100 thru 200 : System/Model Properties
      200 thru 300 : Hydrological Processes/Water movers
      300 thru 400 : Source/sinks
      ------------------------------------------------------------------
    */
    
    code=0;
    //---------------------SPECIAL -----------------------------
    if       (Len==0)                                     {code=-1; }
    else if  (!strcmp(s[0],"*"                          )){code=-2; }//comment
    else if  (!strcmp(s[0],"%"                          )){code=-2; }//comment
    else if  (!strcmp(s[0],"#"                          )){code=-2; }//comment
    else if  (s[0][0]=='#')                               {code=-2; }//comment
    else if  (!strcmp(s[0],":End"                       )){code=-3; }//premature end of file
    //--------------------MODEL OPTIONS ------------------------
    else if  (!strcmp(s[0],"?? "                        )){code=1;  }
    else if  (!strcmp(s[0],":JulianStartDay"            )){code=2;  }
    else if  (!strcmp(s[0],":JulianStartYear"           )){code=3;  }
    else if  (!strcmp(s[0],":Duration"                  )){code=4;  }

    else if  (!strcmp(s[0],":Method"                    )){code=5;  }
    else if  (!strcmp(s[0],":NumericalMethod"           )){code=5;  }
    else if  (!strcmp(s[0],":TimeStep"                  )){code=6;  }

    else if  (!strcmp(s[0],":Routing"                   )){code=8;  }
    else if  (!strcmp(s[0],":Evaporation"               )){code=9;  }
    else if  (!strcmp(s[0],":InterpolationMethod"       )){code=10; }
    else if  (!strcmp(s[0],":Interpolation"             )){code=10; }
    else if  (!strcmp(s[0],":MetGaugeInterpolation"     )){code=10; }
    else if  (!strcmp(s[0],":OW_Evaporation"            )){code=12; }
    else if  (!strcmp(s[0],":EndPause"                  )){code=13; }
    else if  (!strcmp(s[0],":SoilModel"                 )){code=14; }//REQUIRED- CREATES MODEL!
    else if  (!strcmp(s[0],":CatchmentRouting"          )){code=16; }
    else if  (!strcmp(s[0],":CatchmentRoute"            )){code=16; }
    else if  (!strcmp(s[0],":LakeStorage"               )){code=17; }//AFTER SoilModel Commmand
    else if  (!strcmp(s[0],":OroPETCorrect"             )){code=18; }
    else if  (!strcmp(s[0],":RainSnowMethod"            )){code=21; }
    else if  (!strcmp(s[0],":RainSnowFraction"          )){code=21; }
    else if  (!strcmp(s[0],":CloudCoverMethod"          )){code=22; }
    else if  (!strcmp(s[0],":LWRadiationMethod"         )){code=23; }
    else if  (!strcmp(s[0],":SWRadiationMethod"         )){code=24; }
    else if  (!strcmp(s[0],":MonthlyInterpolationMethod")){code=25; }
    else if  (!strcmp(s[0],":RelativeHumidityMethod"    )){code=26; }
    else if  (!strcmp(s[0],":StartDate"                 )){code=27; }
    else if  (!strcmp(s[0],":WindspeedMethod"           )){code=28; }
    else if  (!strcmp(s[0],":AirPressureMethod"         )){code=29; }
    else if  (!strcmp(s[0],":PrecipIceptFract"          )){code=30; }
    else if  (!strcmp(s[0],":OroTempCorrect"            )){code=31; }
    else if  (!strcmp(s[0],":OroPrecipCorrect"          )){code=32; }
    else if  (!strcmp(s[0],":AquiferLayers"             )){code=33; }//AFTER :SoilModel Commmand
    else if  (!strcmp(s[0],":PotentialMeltMethod"       )){code=34; }
    else if  (!strcmp(s[0],":SubdailyMethod"            )){code=35; }
    else if  (!strcmp(s[0],":SWCanopyCorrect"           )){code=36; }
    else if  (!strcmp(s[0],":SWCloudCorrect"            )){code=37; }
    else if  (!strcmp(s[0],":MultilayerSnow"            )){code=39; }//AFTER :SoilModel Commmand
    else if  (!strcmp(s[0],":RetainUBCWMBugs"           )){code=40; }
    else if  (!strcmp(s[0],":EndDate"                   )){code=41; }
    else if  (!strcmp(s[0],":RechargeMethod"            )){code=42; }
    else if  (!strcmp(s[0],":NetSWRadMethod"            )){code=43; }
    else if  (!strcmp(s[0],":DirectEvaporation"         )){code=44; }
    else if  (!strcmp(s[0],":Calendar"                  )){code=45; }
    else if  (!strcmp(s[0],":SnowCoverDepletion"        )){code=46; }
    else if  (!strcmp(s[0],":EnsembleMode"              )){code=47; }
    else if  (!strcmp(s[0],":SuppressCompetitiveET"     )){code=48; }
    else if  (!strcmp(s[0],":SnowSuppressesPET"         )){code=49; }
	//---I/O------------------------------------------------------
    else if  (!strcmp(s[0],":DebugMode"                 )){code=50; }
    else if  (!strcmp(s[0],":WriteMassBalanceFile"      )){code=51; }
    else if  (!strcmp(s[0],":WriteForcingFunctions"     )){code=52; }
    else if  (!strcmp(s[0],":WriteEnergyStorage"        )){code=53; }
    else if  (!strcmp(s[0],":WriteParametersFile"       )){code=54; }
    else if  (!strcmp(s[0],":WriteEnsimFormat"          )){code=55; }
    else if  (!strcmp(s[0],":WriteNetcdfFormat"         )){code=78; }
    else if  (!strcmp(s[0],":WriteNetCDFFormat"         )){code=78; }
    else if  (!strcmp(s[0],":RunName"                   )){code=56; }
    else if  (!strcmp(s[0],":NoisyMode"                 )){code=57; }
    else if  (!strcmp(s[0],":SilentMode"                )){code=58; }
    else if  (!strcmp(s[0],":PavicsMode"                )){code=88; }//some special options only for PAVICS
    else if  (!strcmp(s[0],":rvh_Filename"              )){code=59; }
    else if  (!strcmp(s[0],":rvp_Filename"              )){code=60; }
    else if  (!strcmp(s[0],":rvt_Filename"              )){code=61; }
    else if  (!strcmp(s[0],":rvc_Filename"              )){code=62; }
    else if  (!strcmp(s[0],":rve_Filename"              )){code=91; }
    else if  (!strcmp(s[0],":OutputDirectory"           )){code=63; }//Ideally, called before everything else
    else if  (!strcmp(s[0],":OutputDump"                )){code=64; }
    else if  (!strcmp(s[0],":MajorOutputInterval"       )){code=65; }
    else if  (!strcmp(s[0],":SnapshotHydrograph"        )){code=66; }
    else if  (!strcmp(s[0],":HRUStorageOutput"          )){code=67; }//After corresponding DefineHRUGroup(s) command
    else if  (!strcmp(s[0],":SuppressWarnings"          )){code=68; }
    else if  (!strcmp(s[0],":QuietMode"                 )){code=69; }
    else if  (!strcmp(s[0],":WriteExhaustiveMB"         )){code=70; }
    else if  (!strcmp(s[0],":EvaluationMetrics"         )){code=71; }
    else if  (!strcmp(s[0],":SuppressOutputICs"         )){code=72; }
    else if  (!strcmp(s[0],":SuppressOutput"            )){code=73; }
    else if  (!strcmp(s[0],":WriteHRUGroupMBFile"       )){code=74; }
    else if  (!strcmp(s[0],":OutputInterval"            )){code=15; }
    else if  (!strcmp(s[0],":EvaluationTime"            )){code=75; }//After StartDate or JulianStartDay and JulianStartYear commands
    else if  (!strcmp(s[0],":WaterYearStartMonth"       )){code=76; }
    else if  (!strcmp(s[0],":CreateRVPTemplate"         )){code=77; } 
    else if  (!strcmp(s[0],":WriteChannelInfo"          )){code=79; } 
    else if  (!strcmp(s[0],":BenchmarkingMode"          )){code=85; } 
    else if  (!strcmp(s[0],":WriteReservoirMBFile"      )){code=86; }
    else if  (!strcmp(s[0],":PeriodStartingFormatOff"   )){code=87; }
    else if  (!strcmp(s[0],":OutputConstituentMass"     )){code=89; }
    else if  (!strcmp(s[0],":NetCDFAttribute"           )){code=90; }
    else if  (!strcmp(s[0],":AssimilationStartTime"     )){code=92; }
    else if  (!strcmp(s[0],":AssimilateStreamflow"      )){code=93; }
    else if  (!strcmp(s[0],":DeltaresFEWSMode"          )){code=94; }
    else if  (!strcmp(s[0],":WriteSubbasinFile"         )){code=95; }
    else if  (!strcmp(s[0],":DontWriteWatershedStorage" )){code=96; }//avoid writing WatershedStorage.csv
    else if  (!strcmp(s[0],":TimeZone"                  )){code=97; }
    else if  (!strcmp(s[0],":WriteInterpolationWeights" )){code=101;}
    else if  (!strcmp(s[0],":CallExternalScript"        )){code=102;}
    else if  (!strcmp(s[0],":rvl_Filename"              )){code=103;}
    else if  (!strcmp(s[0],":ReadLiveFile"              )){code=104;}
    else if  (!strcmp(s[0],":CustomOutputInterval"      )){code=105;}
    else if  (!strcmp(s[0],":UseStopFile"               )){code=106;}
    else if  (!strcmp(s[0],":WriteDemandFile"           )){code=107;}

    else if  (!strcmp(s[0],":WriteGroundwaterHeads"     )){code=510;}
    else if  (!strcmp(s[0],":WriteGroundwaterFlows"     )){code=511;}
    else if  (!strcmp(s[0],":rvg_Filename"              )){code=512;}
    //-----------------------------------------------------------
    else if  (!strcmp(s[0],":DefineHRUGroup"            )){code=80; }//After :SoilModel command
    else if  (!strcmp(s[0],":DefineHRUGroups"           )){code=81; }//After :SoilModel command
    else if  (!strcmp(s[0],":DisableHRUGroup"           )){code=82; } 
    //-----------------------------------------------------------
    else if  (!strcmp(s[0],":Alias"                     )){code=98; }
    else if  (!strcmp(s[0],":CustomOutput"              )){code=99; }
    //--------------------SYSTEM OPTIONS -----------------------
    else if  (!strcmp(s[0],":StorageVars"               )){code=100;}//OBSOLETE
    else if  (!strcmp(s[0],":StateVariables"            )){code=100;}//OBSOLETE
    else if  (!strcmp(s[0],":AggregatedVariable"        )){code=150;}//After corresponding DefineHRUGroup(s) command
    //--------------------HYDROLOGICAL PROCESSES ---------------
    if       (!strcmp(s[0],":ProcessBegin"              )){code=200;}//REQUIRED
    else if  (!strcmp(s[0],":HydrologicProcesses"       )){code=200;}//REQUIRED
    else if  (!strcmp(s[0],":HydrologicalProcesses"     )){code=200;}//REQUIRED
    else if  (!strcmp(s[0],":Baseflow"                  )){code=201;}
    else if  (!strcmp(s[0],":CanopyEvaporation"         )){code=202;}
    else if  (!strcmp(s[0],":CanopyDrip"                )){code=203;}
    else if  (!strcmp(s[0],":Infiltration"              )){code=204;}//GENERALLY REQUIRED
    else if  (!strcmp(s[0],":Percolation"               )){code=205;}
    else if  (!strcmp(s[0],":Snowmelt"                  )){code=206;}
    else if  (!strcmp(s[0],":SnowMelt"                  )){code=206;}
    else if  (!strcmp(s[0],":SoilEvaporation"           )){code=208;}
    else if  (!strcmp(s[0],":SnowBalance"               )){code=209;}
    else if  (!strcmp(s[0],":Sublimation"               )){code=210;}
    else if  (!strcmp(s[0],":OpenWaterEvaporation"      )){code=211;}
    else if  (!strcmp(s[0],":Precipitation"             )){code=212;}//GENERALLY REQUIRED
    else if  (!strcmp(s[0],":Interflow"                 )){code=213;}
    else if  (!strcmp(s[0],":SnowRefreeze"              )){code=214;}
    else if  (!strcmp(s[0],":Refreeze"                  )){code=214;}
    else if  (!strcmp(s[0],":Flush"                     )){code=215;}
    else if  (!strcmp(s[0],":CapillaryRise"             )){code=216;}
    else if  (!strcmp(s[0],":LakeEvaporation"           )){code=217;}
    else if  (!strcmp(s[0],":SnowSqueeze"               )){code=218;}
    else if  (!strcmp(s[0],":GlacialMelt"               )){code=219;}
    else if  (!strcmp(s[0],":GlacierMelt"               )){code=219;}
    else if  (!strcmp(s[0],":GlacierRelease"            )){code=220;}
    else if  (!strcmp(s[0],":CanopySnowEvaporation"     )){code=221;}
    else if  (!strcmp(s[0],":CanopySnowEvap"            )){code=221;}
    else if  (!strcmp(s[0],":Overflow"                  )){code=222;}
    else if  (!strcmp(s[0],":-->Overflow"               )){code=222;}
    else if  (!strcmp(s[0],":SnowAlbedoEvolve"          )){code=223;}
    else if  (!strcmp(s[0],":CropHeatUnitEvolve"        )){code=224;}
    else if  (!strcmp(s[0],":Abstraction"               )){code=225;}
    else if  (!strcmp(s[0],":GlacierInfiltration"       )){code=226;}
    else if  (!strcmp(s[0],":Split"                     )){code=227;}
    else if  (!strcmp(s[0],":Convolve"                  )){code=228;}
    else if  (!strcmp(s[0],":SnowTempEvolve"            )){code=229;}
    else if  (!strcmp(s[0],":DepressionOverflow"        )){code=230;}
    else if  (!strcmp(s[0],":ExchangeFlow"              )){code=231;}
    else if  (!strcmp(s[0],":LateralFlush"              )){code=232;}
    else if  (!strcmp(s[0],":Seepage"                   )){code=233;}
    else if  (!strcmp(s[0],":Recharge"                  )){code=234;}
    else if  (!strcmp(s[0],":BlowingSnow"               )){code=235;}
    else if  (!strcmp(s[0],":LakeRelease"               )){code=236;}
    else if  (!strcmp(s[0],":-->RedirectFlow"           )){code=237;}
    else if  (!strcmp(s[0],":Drain"                     )){code=252;}
  
    //...
    else if  (!strcmp(s[0],":ProcessGroup"              )){code=295;}
    else if  (!strcmp(s[0],":EndProcessGroup"           )){code=296;}
    else if  (!strcmp(s[0],":-->Conditional"            )){code=297;}
    else if  (!strcmp(s[0],":EndHydrologicProcesses"    )){code=298;}
    else if  (!strcmp(s[0],":-->Cascade"                )){code=299;}
    //...
    //--------------------TRANSPORT PROCESSES ---------------
    if       (!strcmp(s[0],":Transport"                 )){code=300;}
    else if  (!strcmp(s[0],":FixedConcentration"        )){code=301;}//After corresponding DefineHRUGroup(s) command, if used
    else if  (!strcmp(s[0],":MassInflux"                )){code=302;}//After corresponding DefineHRUGroup(s) command, if used
    else if  (!strcmp(s[0],":GeochemicalProcesses"      )){code=303;}
    else if  (!strcmp(s[0],":EndGeochemicalProcesses"   )){code=304;}
    else if  (!strcmp(s[0],":Decay"                     )){code=305;}
    else if  (!strcmp(s[0],":Advection"                 )){code=306;}
    else if  (!strcmp(s[0],":Transformation"            )){code=307;}
    else if  (!strcmp(s[0],":Mineralization"            )){code=308;}
    else if  (!strcmp(s[0],":HeatConduction"            )){code=309;}
    //...
    //-------------------WATER MANAGEMENT---------------------
    else if  (!strcmp(s[0],":ReservoirDemandAllocation" )){code=400; }
    else if  (!strcmp(s[0],":ReservoirOverflowMode"     )){code=401; }
    //
    //-------------------GROUNDWATER -------------------------
    else if  (!strcmp(s[0],":ModelType"             )){code=500; }//AFTER SoilModel Commmand
    else if  (!strcmp(s[0],":GWMethodOuter"         )){code=501; }
    else if  (!strcmp(s[0],":GWMethodInner"         )){code=502; }
    else if  (!strcmp(s[0],":GWDiscretization"      )){code=503; }
    else if  (!strcmp(s[0],":GWSWExchangeFreq"      )){code=504; }
    else if  (!strcmp(s[0],":FluxExchange"          )){code=505; }
    else if  (!strcmp(s[0],":OverlapType"           )){code=506; }	

    ExitGracefullyIf((code>200) && (code<300) && (pModel==NULL),
                     "ParseMainInputFile: :HydrologicProcesses AND :SoilModel commands must be called before hydrologic processes are specified",BAD_DATA);
    ExitGracefullyIf((code>300) && (code<400) && (pModel->GetTransportModel()==NULL),
                     "ParseMainInputFile: :Transport command must be called before geochemical processes are specified",BAD_DATA);

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
    {/*:End*/
      if (Options.noisy) {cout <<"EOF"<<endl;} ended=true; break;
    }
    case(2):  //----------------------------------------------
    {/*Julian Start Day
       :JulianStartDay [double day] */
      if (Options.noisy) {cout <<"Julian Start Day"<<endl;}
      if (Len<2){ImproperFormatWarning(":JulianStartDay",p,Options.noisy); break;}
      Options.julian_start_day =s_to_d(s[1]);
      break;
    }
    case(3):  //----------------------------------------------
    {/*Julian Start Year
       :JulianStartYear [int year] */
      if (Options.noisy) {cout <<"Julian Start Year"<<endl;}
      if (Len<2){ImproperFormatWarning(":JulianStartYear",p,Options.noisy); break;}
      Options.julian_start_year =s_to_i(s[1]);
      break;
    }
    case(4):  //----------------------------------------------
    {/*Simulation Duration
       :Duration [double time]  /or/
       :Duration [double time] [hh:mm:ss.00]       
       */
      if (Options.noisy) {cout <<"Simulation duration"<<endl;}
      if (Len<2){ImproperFormatWarning(":Duration",p,Options.noisy);  break;}
      double partday=0.0;
      if ((Len>=3) &&  (*(s[2])!='#') && (*(s[2])!='*')){
        string tString=s[2];
        if((tString.length()>=2) && ((tString.substr(2,1)==":") || (tString.substr(1,1)==":")))//support for hh:mm:ss.00 format
        {
          partday=DateStringToTimeStruct("0000-01-01",tString,Options.calendar).julian_day;
        }
      }
      Options.duration =s_to_d(s[1])+partday;
      break;
    }
    case(5):  //----------------------------------------------
    {/*Numerical Simulation Method
       :NumericalMethod [string method] {optional more terms}*/
      if (Options.noisy) {cout <<"Numerical Simulation Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":NumericalMethod",p,Options.noisy); break;}
      if       (!strcmp(s[1],"EULER"            )){Options.sol_method =EULER;}
      else if  (!strcmp(s[1],"STANDARD"         )){Options.sol_method =ORDERED_SERIES;}
      else if  (!strcmp(s[1],"ORDERED_SERIES"   )){Options.sol_method =ORDERED_SERIES;}
      else if  (!strcmp(s[1],"ITER_HEUN")){
        Options.sol_method =ITERATED_HEUN;
        Options.convergence_crit = s_to_d(s[2]);  // take in convergence criteria
        Options.max_iterations   = s_to_d(s[3]);  // maximum number of iterations allowed
      }
      //else if  (!strcmp(s[1],"RUNGE_KUTTA"      )){Options.sol_method =RUNGE_KUTTA_4;}

      //...
      break;
    }
    case(6):  //----------------------------------------------
    {/*Simulation Method time Step
       string ":TimeStep", double tstep /or/
       string ":TimeStep", string hh:mm:ss.00
     */
      if (Options.noisy) {cout <<"Simulation Time Step"<<endl;}
      if (Len<2){ImproperFormatWarning(":TimeStep",p,Options.noisy); break;}
      string tString=s[1];
      if ((tString.length()>=2) && ((tString.substr(2,1)==":") || (tString.substr(1,1)==":"))){//support for hh:mm:ss.00 format
        time_struct tt;
        tt=DateStringToTimeStruct("0000-01-01",tString,Options.calendar);
        Options.timestep=FixTimestep(tt.julian_day);
      }
      else{
        Options.timestep =FixTimestep(s_to_d(s[1]));
      }
      break;
    }
    case(7):  //----------------------------------------------
    {
      break;
    }
    case(8):  //----------------------------------------------
    {/*Routing Method
       string ":Routing" string method */
      if (Options.noisy) {cout <<"Routing Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":Routing",p,Options.noisy); break;}
      if      (!strcmp(s[1],"MUSKINGUM"              )){Options.routing =ROUTE_MUSKINGUM;}
      else if (!strcmp(s[1],"ROUTE_MUSKINGUM"        )){Options.routing =ROUTE_MUSKINGUM;}
      else if (!strcmp(s[1],"ROUTE_MUSKINGUM_CUNGE"  )){Options.routing =ROUTE_MUSKINGUM_CUNGE;}
      else if (!strcmp(s[1],"ROUTE_STORAGECOEFF"     )){Options.routing =ROUTE_STORAGECOEFF;}
      else if (!strcmp(s[1],"ROUTE_PLUG_FLOW"        )){Options.routing =ROUTE_PLUG_FLOW;}
      else if (!strcmp(s[1],"ROUTE_DIFFUSIVE_WAVE"   )){Options.routing =ROUTE_DIFFUSIVE_WAVE;}
      else if (!strcmp(s[1],"ROUTE_HYDROLOGIC"       )){Options.routing =ROUTE_HYDROLOGIC;}
      else if (!strcmp(s[1],"ROUTE_NONE"             )){Options.routing =ROUTE_NONE;}
      else if (!strcmp(s[1],"ROUTE_TVD"              )){Options.routing =ROUTE_TVD;}
      else if (!strcmp(s[1],"NONE"                   )){Options.routing =ROUTE_NONE;}
      else{
        ExitGracefully("ParseMainInputFile: Unrecognized routing method",BAD_DATA_WARN);
      }
      break;
    }
    case(9):  //----------------------------------------------
    {/*Potential Evapotranspiration Method
       string ":Evaporation", string method */
      if (Options.noisy) {cout <<"Potential Evapotranspiration Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":Evaporation",p,Options.noisy); break;}
      if      (!strcmp(s[1],"CONSTANT"              )){Options.evaporation =PET_CONSTANT;}
      else if (!strcmp(s[1],"DATA"                  )){Options.evaporation =PET_DATA;}
      else if (!strcmp(s[1],"FROM_MONTHLY"          )){Options.evaporation =PET_FROMMONTHLY;}
      else if (!strcmp(s[1],"PENMAN_MONTEITH"       )){Options.evaporation =PET_PENMAN_MONTEITH;}
      else if (!strcmp(s[1],"PENMAN_COMBINATION"    )){Options.evaporation =PET_PENMAN_COMBINATION;}
      else if (!strcmp(s[1],"HAMON"                 )){Options.evaporation =PET_HAMON;}
      else if (!strcmp(s[1],"HARGREAVES"            )){Options.evaporation =PET_HARGREAVES;}
      else if (!strcmp(s[1],"HARGREAVES_1985"       )){Options.evaporation =PET_HARGREAVES_1985;}
      else if (!strcmp(s[1],"TURC_1961"             )){Options.evaporation =PET_TURC_1961;}
      else if (!strcmp(s[1],"MAKKINK_1957"          )){Options.evaporation =PET_MAKKINK_1957;}
      else if (!strcmp(s[1],"PRIESTLEY_TAYLOR"      )){Options.evaporation =PET_PRIESTLEY_TAYLOR;}
      else if (!strcmp(s[1],"MONTHLY_FACTOR"        )){Options.evaporation =PET_MONTHLY_FACTOR;}
      else if (!strcmp(s[1],"PET_CONSTANT"          )){Options.evaporation =PET_CONSTANT;}
      else if (!strcmp(s[1],"PET_DATA"              )){Options.evaporation =PET_DATA;}
      else if (!strcmp(s[1],"PET_FROMMONTHLY"       )){Options.evaporation =PET_FROMMONTHLY;}
      else if (!strcmp(s[1],"PET_PENMAN_MONTEITH"   )){Options.evaporation =PET_PENMAN_MONTEITH;}
      else if (!strcmp(s[1],"PET_PENMAN_COMBINATION")){Options.evaporation =PET_PENMAN_COMBINATION;}
      else if (!strcmp(s[1],"PET_HAMON"             )){Options.evaporation =PET_HAMON;}
      else if (!strcmp(s[1],"PET_HARGREAVES"        )){Options.evaporation =PET_HARGREAVES;}
      else if (!strcmp(s[1],"PET_HARGREAVES_1985"   )){Options.evaporation =PET_HARGREAVES_1985;}
      else if (!strcmp(s[1],"PET_TURC_1961"         )){Options.evaporation =PET_TURC_1961;}
      else if (!strcmp(s[1],"PET_MAKKINK_1957"      )){Options.evaporation =PET_MAKKINK_1957;}
      else if (!strcmp(s[1],"PET_PRIESTLEY_TAYLOR"  )){Options.evaporation =PET_PRIESTLEY_TAYLOR;}
      else if (!strcmp(s[1],"PET_MONTHLY_FACTOR"    )){Options.evaporation =PET_MONTHLY_FACTOR;}
      else if (!strcmp(s[1],"PET_PENMAN_SIMPLE33"   )){Options.evaporation =PET_PENMAN_SIMPLE33;}
      else if (!strcmp(s[1],"PET_PENMAN_SIMPLE39"   )){Options.evaporation =PET_PENMAN_SIMPLE39;}
      else if (!strcmp(s[1],"PET_GRANGERGRAY"       )){Options.evaporation =PET_GRANGERGRAY;}
      else if (!strcmp(s[1],"PET_MOHYSE"            )){Options.evaporation =PET_MOHYSE;}
      else if (!strcmp(s[1],"PET_OUDIN"             )){Options.evaporation =PET_OUDIN;}
      else if (!strcmp(s[1],"PET_LINACRE"           )){Options.evaporation =PET_LINACRE; }

      else{
        ExitGracefully("ParseMainInputFile: Unrecognized PET calculation method",BAD_DATA_WARN);
      }
      bool notspecified=true;
      if((Len>=3) && (Options.evaporation==PET_DATA))
      {
        notspecified=false;
        if      (!strcmp(s[2],"PET_CONSTANT"          )) { Options.evap_infill =PET_CONSTANT; }
        else if (!strcmp(s[2],"PET_FROMMONTHLY"       )) { Options.evap_infill =PET_FROMMONTHLY; }
        else if (!strcmp(s[2],"PET_PENMAN_MONTEITH"   )) { Options.evap_infill =PET_PENMAN_MONTEITH; }
        else if (!strcmp(s[2],"PET_PENMAN_COMBINATION")) { Options.evap_infill =PET_PENMAN_COMBINATION; }
        else if (!strcmp(s[2],"PET_HAMON"             )) { Options.evap_infill =PET_HAMON; }
        else if (!strcmp(s[2],"PET_HARGREAVES"        )) { Options.evap_infill =PET_HARGREAVES; }
        else if (!strcmp(s[2],"PET_HARGREAVES_1985"   )) { Options.evap_infill =PET_HARGREAVES_1985; }
        else if (!strcmp(s[2],"PET_TURC_1961"         )) { Options.evap_infill =PET_TURC_1961; }
        else if (!strcmp(s[2],"PET_MAKKINK_1957"      )) { Options.evap_infill =PET_MAKKINK_1957; }
        else if (!strcmp(s[2],"PET_PRIESTLEY_TAYLOR"  )) { Options.evap_infill =PET_PRIESTLEY_TAYLOR; }
        else if (!strcmp(s[2],"PET_MONTHLY_FACTOR"    )) { Options.evap_infill =PET_MONTHLY_FACTOR; }
        else if (!strcmp(s[2],"PET_PENMAN_SIMPLE33"   )) { Options.evap_infill =PET_PENMAN_SIMPLE33;}
        else if (!strcmp(s[2],"PET_PENMAN_SIMPLE39"   )) { Options.evap_infill =PET_PENMAN_SIMPLE39;}
        else if (!strcmp(s[2],"PET_GRANGERGRAY"       )) { Options.evap_infill =PET_GRANGERGRAY;}
        else if (!strcmp(s[2],"PET_MOHYSE"            )) { Options.evap_infill =PET_MOHYSE;}
        else if (!strcmp(s[2],"PET_OUDIN"             )) { Options.evap_infill =PET_OUDIN;}
        else if (!strcmp(s[2],"PET_LINACRE"           )) { Options.evap_infill =PET_LINACRE; }
        else { notspecified=true; }
      }
      if((Options.evaporation==PET_DATA) && (notspecified)) {
        string warning="PET_DATA was chosen, but no infilling method was specified. Blank PET data will be infilled using PET_HARGREAVES_1985";
        WriteAdvisory(warning,Options.noisy);
      }
      break;
    }
    case(10):  //----------------------------------------------
    {/*Interpolation Method
       string ":Interpolation" string method */
      if (Options.noisy) {cout <<"Interpolation Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":Interpolation",p,Options.noisy); break;}
      if      (!strcmp(s[1],"NEAREST_NEIGHBOR"           )){Options.interpolation=INTERP_NEAREST_NEIGHBOR;}
      else if (!strcmp(s[1],"AVERAGE_ALL"                )){Options.interpolation=INTERP_AVERAGE_ALL;}
      else if (!strcmp(s[1],"INVERSE_DISTANCE"           )){Options.interpolation=INTERP_INVERSE_DISTANCE;}
      else if (!strcmp(s[1],"FROM_FILE"                  )){Options.interpolation=INTERP_FROM_FILE;}
      else if (!strcmp(s[1],"INTERP_NEAREST_NEIGHBOR"    )){Options.interpolation=INTERP_NEAREST_NEIGHBOR;}
      else if (!strcmp(s[1],"INTERP_AVERAGE_ALL"         )){Options.interpolation=INTERP_AVERAGE_ALL;}
      else if (!strcmp(s[1],"INTERP_INVERSE_DISTANCE"    )){Options.interpolation=INTERP_INVERSE_DISTANCE;}
      else if (!strcmp(s[1],"INTERP_INVERSE_DISTANCE_ELEVATION")){ Options.interpolation=INTERP_INVERSE_DISTANCE_ELEVATION; }
      else if (!strcmp(s[1],"INTERP_FROM_FILE"           )){Options.interpolation=INTERP_FROM_FILE;}
      else{
        ExitGracefully("ParseMainInputFile: Unrecognized interpolation method",BAD_DATA_WARN);
      }
      if ((Options.interpolation==INTERP_FROM_FILE) && (Len>2))
      {
        Options.interp_file="";
        for (int i=2; i<Len-1;i++){
          Options.interp_file+=s[i];
          Options.interp_file+=" ";
        }
        Options.interp_file+=s[Len-1];

        Options.interp_file =CorrectForRelativePath(Options.interp_file ,Options.rvi_filename);
      }
      break;
    }
    case(12): //----------------------------------------------
    {/*Open Water Potential Evapotranspiration Method
       string ":OW_Evaporation", string method  {optional infill method}*/
      if (Options.noisy) {cout <<"Open Water Potential Evapotranspiration Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":OW_Evaporation",p,Options.noisy); break;}
      if      (!strcmp(s[1],"CONSTANT"              )){Options.ow_evaporation =PET_CONSTANT;}
      else if (!strcmp(s[1],"FROM_FILE"             )){Options.ow_evaporation =PET_DATA;}
      else if (!strcmp(s[1],"FROM_MONTHLY"          )){Options.ow_evaporation =PET_FROMMONTHLY;}
      else if (!strcmp(s[1],"PENMAN_MONTEITH"       )){Options.ow_evaporation =PET_PENMAN_MONTEITH;}
      else if (!strcmp(s[1],"PENMAN_COMBINATION"    )){Options.ow_evaporation =PET_PENMAN_COMBINATION;}
      else if (!strcmp(s[1],"PRIESTLEY_TAYLOR"      )){Options.ow_evaporation =PET_PRIESTLEY_TAYLOR;}
      else if (!strcmp(s[1],"HARGREAVES"            )){Options.ow_evaporation =PET_HARGREAVES;}
      else if (!strcmp(s[1],"HARGREAVES_1985"       )){Options.ow_evaporation =PET_HARGREAVES_1985;}
      else if (!strcmp(s[1],"MONTHLY_FACTOR"        )){Options.ow_evaporation =PET_MONTHLY_FACTOR;}
      else if (!strcmp(s[1],"PET_CONSTANT"          )){Options.ow_evaporation =PET_CONSTANT;}
      else if (!strcmp(s[1],"PET_DATA"              )){Options.ow_evaporation =PET_DATA;}
      else if (!strcmp(s[1],"PET_FROMMONTHLY"       )){Options.ow_evaporation =PET_FROMMONTHLY;}
      else if (!strcmp(s[1],"PET_PENMAN_MONTEITH"   )){Options.ow_evaporation =PET_PENMAN_MONTEITH;}
      else if (!strcmp(s[1],"PET_PENMAN_COMBINATION")){Options.ow_evaporation =PET_PENMAN_COMBINATION;}
      else if (!strcmp(s[1],"PET_HAMON"             )){Options.ow_evaporation =PET_HAMON;}
      else if (!strcmp(s[1],"PET_HARGREAVES"        )){Options.ow_evaporation =PET_HARGREAVES;}
      else if (!strcmp(s[1],"PET_HARGREAVES_1985"   )){Options.ow_evaporation =PET_HARGREAVES_1985;}
      else if (!strcmp(s[1],"PET_TURC_1961"         )){Options.ow_evaporation =PET_TURC_1961;}
      else if (!strcmp(s[1],"PET_MAKKINK_1957"      )){Options.ow_evaporation =PET_MAKKINK_1957;}
      else if (!strcmp(s[1],"PET_PRIESTLEY_TAYLOR"  )){Options.ow_evaporation =PET_PRIESTLEY_TAYLOR;}
      else if (!strcmp(s[1],"PET_MONTHLY_FACTOR"    )){Options.ow_evaporation =PET_MONTHLY_FACTOR;}
      else if (!strcmp(s[1],"PET_PENMAN_SIMPLE33"   )){Options.ow_evaporation =PET_PENMAN_SIMPLE33;}
      else if (!strcmp(s[1],"PET_PENMAN_SIMPLE39"   )){Options.ow_evaporation =PET_PENMAN_SIMPLE39;}
      else if (!strcmp(s[1],"PET_GRANGERGRAY"       )){Options.ow_evaporation =PET_GRANGERGRAY;}
      else if (!strcmp(s[1],"PET_MOHYSE"            )){Options.ow_evaporation =PET_MOHYSE;}
      else if (!strcmp(s[1],"PET_OUDIN"             )){Options.ow_evaporation =PET_OUDIN;}
      else if (!strcmp(s[1],"PET_LINACRE"           )){Options.ow_evaporation =PET_LINACRE; }
      else{
        ExitGracefully("ParseMainInputFile: Unrecognized PET calculation method",BAD_DATA_WARN);
      }
      bool notspecified=true;
      if((Len>=3) && (Options.ow_evaporation==PET_DATA))
      {
        notspecified=false;
        if      (!strcmp(s[2],"PET_CONSTANT"          )) { Options.ow_evap_infill =PET_CONSTANT; }
        else if (!strcmp(s[2],"PET_FROMMONTHLY"       )) { Options.ow_evap_infill =PET_FROMMONTHLY; }
        else if (!strcmp(s[2],"PET_PENMAN_MONTEITH"   )) { Options.ow_evap_infill =PET_PENMAN_MONTEITH; }
        else if (!strcmp(s[2],"PET_PENMAN_COMBINATION")) { Options.ow_evap_infill =PET_PENMAN_COMBINATION; }
        else if (!strcmp(s[2],"PET_HAMON"             )) { Options.ow_evap_infill =PET_HAMON; }
        else if (!strcmp(s[2],"PET_HARGREAVES"        )) { Options.ow_evap_infill =PET_HARGREAVES; }
        else if (!strcmp(s[2],"PET_HARGREAVES_1985"   )) { Options.ow_evap_infill =PET_HARGREAVES_1985; }
        else if (!strcmp(s[2],"PET_TURC_1961"         )) { Options.ow_evap_infill =PET_TURC_1961; }
        else if (!strcmp(s[2],"PET_MAKKINK_1957"      )) { Options.ow_evap_infill =PET_MAKKINK_1957; }
        else if (!strcmp(s[2],"PET_PRIESTLEY_TAYLOR"  )) { Options.ow_evap_infill =PET_PRIESTLEY_TAYLOR; }
        else if (!strcmp(s[2],"PET_MONTHLY_FACTOR"    )) { Options.ow_evap_infill =PET_MONTHLY_FACTOR; }
        else if (!strcmp(s[2],"PET_PENMAN_SIMPLE33"   )) { Options.ow_evap_infill =PET_PENMAN_SIMPLE33;}
        else if (!strcmp(s[2],"PET_PENMAN_SIMPLE39"   )) { Options.ow_evap_infill =PET_PENMAN_SIMPLE39;}
        else if (!strcmp(s[2],"PET_GRANGERGRAY"       )) { Options.ow_evap_infill =PET_GRANGERGRAY;}
        else if (!strcmp(s[2],"PET_MOHYSE"            )) { Options.ow_evap_infill =PET_MOHYSE;}
        else if (!strcmp(s[2],"PET_OUDIN"             )) { Options.ow_evap_infill =PET_OUDIN;}
        else if (!strcmp(s[1],"PET_LINACRE"           )) { Options.ow_evap_infill =PET_LINACRE; }
        else{notspecified=true;}
      }
      if ((Options.ow_evaporation==PET_DATA) && (notspecified)) {
        string warning="PET_DATA was chosen, but no infilling method was specified. Blank OW PET data will be infilled using PET_HARGREAVES_1985";
        WriteAdvisory(warning,Options.noisy);
      }
      break;
    }
    case(13): //----------------------------------------------
    {/*Pause at the end of the simulation
       string ":EndPause"*/
      if (Options.noisy) { cout <<"End of Simulation Pause "<<endl; }
      Options.pause = true;
      break;
    }
    case(14): //----------------------------------------------
    {/*Soil Model Approach
       string ":SoilModel", string method {optional vars}*/
      if (Options.noisy) {cout <<"Soil Model"<<endl;}
      if (Len<2){ImproperFormatWarning(":SoilModel",p,Options.noisy);  break;}
      Options.num_soillayers =1;
      if       (!strcmp(s[1],"SOIL_ONE_LAYER"   )){Options.soil_modeltype =SOIL_ONE_LAYER;}
      else if  (!strcmp(s[1],"SOIL_TWO_LAYER"   )){
        Options.soil_modeltype =SOIL_TWO_LAYER;
        Options.num_soillayers =2;
      }
      else if  (!strcmp(s[1],"SOIL_MULTILAYER"  )){
        Options.soil_modeltype = SOIL_MULTILAYER;
        Options.num_soillayers =s_to_i(s[2]);
      }
      else{
        ExitGracefully("ParseMainInputFile: Unrecognized Soil model",BAD_DATA);
      }
      //Builds model
      pModel=new CModel(Options.soil_modeltype,Options.num_soillayers,Options);
      break;
    }
    case(15): //----------------------------------------------
    {/*Write to output files every x timesteps
       string ":OutputInterval", double interval */
      if (Options.noisy) {cout <<"Output File Interval"<<endl;}
      if (Len<2){ImproperFormatWarning(":OutputInterval",p,Options.noisy);  break;}
      Options.output_interval = s_to_d(s[1]);
      break;
    }
    case(16): //----------------------------------------------
    {/*Catchment Routing Method
       string ":CatchmentRoute" string method */
      if (Options.noisy) {cout <<"Catchment Routing Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":CatchmentRoute",p,Options.noisy); break;}
      if      (!strcmp(s[1],"DUMP"                      )){Options.catchment_routing=ROUTE_DUMP;}
      else if (!strcmp(s[1],"DELAYED"                   )){Options.catchment_routing=ROUTE_DELAYED_FIRST_ORDER;}
      else if (!strcmp(s[1],"DELAYED_FIRST_ORDER"       )){Options.catchment_routing=ROUTE_DELAYED_FIRST_ORDER;}
      else if (!strcmp(s[1],"GAMMA"                     )){Options.catchment_routing=ROUTE_GAMMA_CONVOLUTION;}
      else if (!strcmp(s[1],"GAMMA_UH"                  )){Options.catchment_routing=ROUTE_GAMMA_CONVOLUTION;}
      else if (!strcmp(s[1],"GAMMA_CONVOLUTION"         )){Options.catchment_routing=ROUTE_GAMMA_CONVOLUTION;}
      else if (!strcmp(s[1],"TRIANGULAR_UH"             )){Options.catchment_routing=ROUTE_TRI_CONVOLUTION;}
      else if (!strcmp(s[1],"TRI_CONVOLUTION"           )){Options.catchment_routing=ROUTE_TRI_CONVOLUTION;}
      else if (!strcmp(s[1],"EXPONENTIAL_UH"            )){Options.catchment_routing=ROUTE_RESERVOIR_SERIES;}
      else if (!strcmp(s[1],"SERIES_LINEAR_RESERVOIRS"  )){Options.catchment_routing=ROUTE_RESERVOIR_SERIES;}
      else if (!strcmp(s[1],"RESERVOIRS_SERIES"         )){Options.catchment_routing=ROUTE_RESERVOIR_SERIES;}
      else if (!strcmp(s[1],"ROUTE_DUMP"                )){Options.catchment_routing=ROUTE_DUMP;}
      else if (!strcmp(s[1],"ROUTE_DELAYED_FIRST_ORDER" )){Options.catchment_routing=ROUTE_DELAYED_FIRST_ORDER;}
      else if (!strcmp(s[1],"ROUTE_GAMMA_CONVOLUTION"   )){Options.catchment_routing=ROUTE_GAMMA_CONVOLUTION;}
      else if (!strcmp(s[1],"ROUTE_TRI_CONVOLUTION"     )){Options.catchment_routing=ROUTE_TRI_CONVOLUTION;}
      else if (!strcmp(s[1],"ROUTE_RESERVOIR_SERIES"    )){Options.catchment_routing=ROUTE_RESERVOIR_SERIES;}
      else{
        ExitGracefully("ParseMainInputFile: Unrecognized catchment routing method",BAD_DATA_WARN);
      }
      break;
    }
    case(17): //----------------------------------------------
    {/*Determines where to send lake precipitation - must be after :SoilModel command
       ":LakeStorage" [sv_type lake_storage] */
      if (Options.noisy) {cout <<"Lake Storage"<<endl;}
      if (Len<2){ImproperFormatWarning(":LakeStorage",p,Options.noisy); break;}
      if (pModel==NULL){
        ExitGracefully(":LakeStorage command must be after :SoilModel command in .rvi file.",BAD_DATA_WARN); break;
      }
      tmpS[0]=CStateVariable::StringToSVType(s[1],tmpLev[0],true);
      pModel->AddStateVariables(tmpS,tmpLev,1);
      pModel->SetLakeStorage(tmpS[0],tmpLev[0]);
      break;
    }
    case(18): //----------------------------------------------
    {/*:OroPETCorrect" string method */
      if (Options.noisy) {cout <<"Orographic PET Correction Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":OroPETCorrect",p,Options.noisy); break;}
      if      (!strcmp(s[1],"HBV"             )){Options.orocorr_PET=OROCORR_HBV;}
      else if (!strcmp(s[1],"PRMS"            )){Options.orocorr_PET=OROCORR_PRMS;}
      else if (!strcmp(s[1],"NONE"            )){Options.orocorr_PET=OROCORR_NONE;}
      else if (!strcmp(s[1],"OROCORR_HBV"     )){Options.orocorr_PET=OROCORR_HBV;}
      else if (!strcmp(s[1],"OROCORR_PRMS"    )){Options.orocorr_PET=OROCORR_PRMS;}
      else if (!strcmp(s[1],"OROCORR_UBCWM"   )){Options.orocorr_PET=OROCORR_UBCWM;}
      else if (!strcmp(s[1],"OROCORR_NONE"    )){Options.orocorr_PET=OROCORR_NONE;}
      else{
        ExitGracefully("ParseMainInputFile: Unrecognized Orographic PET Correction Method",BAD_DATA_WARN);
      }
      break;
    }
    case(21): //----------------------------------------------
    {/*:RainSnowFraction" string method */
      if (Options.noisy) {cout <<"Snow - Rain mix calculation Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":RainSnowFraction",p,Options.noisy); break;}
      if      (!strcmp(s[1],"USE_DATA"           )){Options.rainsnow=RAINSNOW_DATA;}
      else if (!strcmp(s[1],"DINGMAN"            )){Options.rainsnow=RAINSNOW_DINGMAN;}
      else if (!strcmp(s[1],"HBV"                )){Options.rainsnow=RAINSNOW_HBV;}
      else if (!strcmp(s[1],"UBC"                )){Options.rainsnow=RAINSNOW_UBCWM;}
      else if (!strcmp(s[1],"HSPF"               )){Options.rainsnow=RAINSNOW_HSPF;}
      else if (!strcmp(s[1],"RAINSNOW_DATA"      )){Options.rainsnow=RAINSNOW_DATA;}
      else if (!strcmp(s[1],"RAINSNOW_DINGMAN"   )){Options.rainsnow=RAINSNOW_DINGMAN;}
      else if (!strcmp(s[1],"RAINSNOW_HBV"       )){Options.rainsnow=RAINSNOW_HBV;}
      else if (!strcmp(s[1],"RAINSNOW_UBCWM"     )){Options.rainsnow=RAINSNOW_UBCWM;}
      else if (!strcmp(s[1],"RAINSNOW_HSPF"      )){Options.rainsnow=RAINSNOW_HSPF;}
      else if (!strcmp(s[1],"RAINSNOW_HARDER"    )){Options.rainsnow=RAINSNOW_HARDER;}
      else {ExitGracefully("ParseInput:RainSnowMethod: Unrecognized method",BAD_DATA_WARN);}
      break;
    }
    case(22): //----------------------------------------------
    {/*:CloudCoverMethod" string method */
      if (Options.noisy) {cout <<"Cloud Cover Estimation Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":CloudCoverMethod",p,Options.noisy); break;}
      if      (!strcmp(s[1],"USE_DATA"           )){Options.cloud_cover=CLOUDCOV_DATA;}
      else if (!strcmp(s[1],"UBC"                )){Options.cloud_cover=CLOUDCOV_UBCWM;}
      else if (!strcmp(s[1],"NONE"               )){Options.cloud_cover=CLOUDCOV_NONE;}
      else if (!strcmp(s[1],"CLOUDCOV_DATA"      )){Options.cloud_cover=CLOUDCOV_DATA;}
      else if (!strcmp(s[1],"CLOUDCOV_UBCWM"     )){Options.cloud_cover=CLOUDCOV_UBCWM;}
      else if (!strcmp(s[1],"CLOUDCOV_NONE"      )){Options.cloud_cover=CLOUDCOV_NONE;}
      else {ExitGracefully("ParseInput:CloudCoverMethod: Unrecognized method",BAD_DATA_WARN);}
      break;
    }
    case(23): //----------------------------------------------
    {/*:LWRadiationMethod" string method */
      if (Options.noisy) {cout <<"Net Longwave Radiation Estimation Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":LWRadiationMethod",p,Options.noisy); break;}
      if      (!strcmp(s[1],"USE_DATA"         )){Options.LW_radiation=LW_RAD_DATA;}
      else if (!strcmp(s[1],"DEFAULT"          )){Options.LW_radiation=LW_RAD_DEFAULT;}
      else if (!strcmp(s[1],"UBC"              )){Options.LW_radiation=LW_RAD_UBCWM;}
      else if (!strcmp(s[1],"HSPF"             )){Options.LW_radiation=LW_RAD_HSPF;}
      else if (!strcmp(s[1],"LW_RAD_DATA"      )){Options.LW_radiation=LW_RAD_DATA;}
      else if (!strcmp(s[1],"LW_RAD_DEFAULT"   )){Options.LW_radiation=LW_RAD_DEFAULT;}
      else if (!strcmp(s[1],"LW_RAD_UBCWM"     )){Options.LW_radiation=LW_RAD_UBCWM;}
      else if (!strcmp(s[1],"LW_RAD_HSPF"      )){Options.LW_radiation=LW_RAD_HSPF;}
      else if (!strcmp(s[1],"LW_RAD_VALIANTZAS")){Options.LW_radiation=LW_RAD_VALIANTZAS;}
      else {ExitGracefully("ParseInput:LWRadiationMethod: Unrecognized method ",BAD_DATA_WARN);}
      break;
    }
    case(24): //----------------------------------------------
    {/*:SWRadiationMethod" string method */
      if (Options.noisy) {cout <<"Shortwave Radiation Estimation Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":SWRadiationMethod",p,Options.noisy); break;}
      if      (!strcmp(s[1],"SW_RAD_DATA"      )){Options.SW_radiation=SW_RAD_DATA;}
      else if (!strcmp(s[1],"SW_RAD_DEFAULT"   )){Options.SW_radiation=SW_RAD_DEFAULT;}
      else if (!strcmp(s[1],"SW_RAD_UBCWM"     )){Options.SW_radiation=SW_RAD_UBCWM;}
      else if (!strcmp(s[1],"UBC"              )){Options.SW_radiation=SW_RAD_UBCWM;}
      else if (!strcmp(s[1],"SW_RAD_VALIANTZAS")){Options.SW_radiation=SW_RAD_VALIANTZAS;}
      else if (!strcmp(s[1],"SW_RAD_NONE"      )){Options.SW_radiation=SW_RAD_NONE; }
      else {ExitGracefully("ParseInput:SWRadiationMethod: Unrecognized method",BAD_DATA_WARN);}

      if (Options.SW_radiation == SW_RAD_DATA){
        Options.SW_cloudcovercorr   =SW_CLOUD_CORR_NONE;
        WriteWarning("Cloud cover corrections have been set to 'NONE', since the shortwave radiation method is SW_RAD_DATA",Options.noisy);
      } //if data provided, then cloudcover corrections not needed
      if(Options.SW_radiation==SW_RAD_NONE) {
        WriteAdvisory("The shortwave radiation calculation method is SW_RAD_NONE. This may impact some snowmelt and PET algorithms which require radiation.",Options.noisy);
      }
      break;
    }
    case(25):  //--------------------------------------------
    {/*:MonthlyInterpolationMethod */
      if (Options.noisy) {cout <<"Monthly Interpolation"<<endl;}
      if (Len<2){ImproperFormatWarning(":MonthlyInterpolationMethod",p,Options.noisy); break;}
      if      (!strcmp(s[1],"UNIFORM"                )){Options.month_interp=MONTHINT_UNIFORM;}
      else if (!strcmp(s[1],"LINEAR_START_OF_MONTH"  )){Options.month_interp=MONTHINT_LINEAR_FOM;}
      else if (!strcmp(s[1],"LINEAR_FOM"             )){Options.month_interp=MONTHINT_LINEAR_FOM;}
      else if (!strcmp(s[1],"LINEAR_21"              )){Options.month_interp=MONTHINT_LINEAR_21;}
      else if (!strcmp(s[1],"LINEAR_MID"             )){Options.month_interp=MONTHINT_LINEAR_MID;}
      else if (!strcmp(s[1],"MONTHINT_UNIFORM"       )){Options.month_interp=MONTHINT_UNIFORM;}
      else if (!strcmp(s[1],"MONTHINT_LINEAR_FOM"    )){Options.month_interp=MONTHINT_LINEAR_FOM;}
      else if (!strcmp(s[1],"MONTHINT_LINEAR_21"     )){Options.month_interp=MONTHINT_LINEAR_21;}
      else if (!strcmp(s[1],"MONTHINT_LINEAR_MID"    )){Options.month_interp=MONTHINT_LINEAR_MID;}
      else {ExitGracefully("ParseInput:MonthlyInterpolationMethod: Unrecognized method",BAD_DATA_WARN);}
      break;
    }
    case(26): //----------------------------------------------
    {/*:RelativeHumidityMethod" string method */
      if (Options.noisy) {cout <<"Relative Humidity Estimation Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":RelativeHumidityMethod",p,Options.noisy); break;}
      if      (!strcmp(s[1],"CONSTANT"         )){Options.rel_humidity=RELHUM_CONSTANT;}
      else if (!strcmp(s[1],"MINDEWPT"         )){Options.rel_humidity=RELHUM_MINDEWPT;}
      else if (!strcmp(s[1],"RELHUM_CONSTANT"  )){Options.rel_humidity=RELHUM_CONSTANT;}
      else if (!strcmp(s[1],"RELHUM_MINDEWPT"  )){Options.rel_humidity=RELHUM_MINDEWPT;}
      else if (!strcmp(s[1],"RELHUM_DATA"      )){Options.rel_humidity=RELHUM_DATA;}
      else {ExitGracefully("ParseInput:RelativeHumidityMethod: Unrecognized method",BAD_DATA_WARN);}
      break;
    }
    case(27): //----------------------------------------------
    {/*:StartDate" string yyyy-mm-dd hh:mm:ss.00*/
      if (Options.noisy) {cout <<"Simulation Start Date"<<endl;}
      if (Len<3){ImproperFormatWarning(":StartDate",p,Options.noisy); break;}
      time_struct tt;
      tt=DateStringToTimeStruct(s[1],s[2],Options.calendar);
      Options.julian_start_day =tt.julian_day;
      Options.julian_start_year=tt.year;
      break;
    }
    case(28): //----------------------------------------------
    {/*:WindspeedMethod" string method */
      if (Options.noisy) {cout <<"Windspeed estimation Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":WindspeedMethod",p,Options.noisy); break;}

      if      (!strcmp(s[1],"CONSTANT"           )){Options.wind_velocity=WINDVEL_CONSTANT;}
      else if (!strcmp(s[1],"USE_DATA"           )){Options.wind_velocity=WINDVEL_DATA;}
      else if (!strcmp(s[1],"UBC"                )){Options.wind_velocity=WINDVEL_UBCWM;}
      else if (!strcmp(s[1],"WINDVEL_CONSTANT"   )){Options.wind_velocity=WINDVEL_CONSTANT;}
      else if (!strcmp(s[1],"WINDVEL_DATA"       )){Options.wind_velocity=WINDVEL_DATA;}
      else if (!strcmp(s[1],"WINDVEL_UBCWM"      )){Options.wind_velocity=WINDVEL_UBCWM;}
      else {ExitGracefully("ParseInput:WindspeedMethod: Unrecognized method",BAD_DATA_WARN);}
      break;
    }
    case(29): //----------------------------------------------
    {/*:AirPressureMethod" string method */
      if (Options.noisy) {cout <<"Air Pressure Estimation Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":AirPressureMethod",p,Options.noisy); break;}
      if      (!strcmp(s[1],"BASIC"         )){Options.air_pressure=AIRPRESS_BASIC;} //to make obsolete
      else if (!strcmp(s[1],"UBC"           )){Options.air_pressure=AIRPRESS_UBC;} //to make obsolete
      else if (!strcmp(s[1],"USE_DATA"      )){Options.air_pressure=AIRPRESS_DATA;} //to make obsolete
      else if (!strcmp(s[1],"CONSTANT"      )){Options.air_pressure=AIRPRESS_CONST;} //to make obsolete
      else if (!strcmp(s[1],"AIRPRESS_BASIC")){Options.air_pressure=AIRPRESS_BASIC;}
      else if (!strcmp(s[1],"AIRPRESS_UBC"  )){Options.air_pressure=AIRPRESS_UBC;}
      else if (!strcmp(s[1],"AIRPRESS_DATA" )){Options.air_pressure=AIRPRESS_DATA;}
      else if (!strcmp(s[1],"AIRPRESS_CONST")){Options.air_pressure=AIRPRESS_CONST;}
      else {ExitGracefully("ParseInput:AirPressureMethod: Unrecognized method",BAD_DATA_WARN);}
      break;
    }
    case(30): //----------------------------------------------
    {/*:PrecipIceptFract" string method */
      if (Options.noisy) {cout <<"Precipitation interception factor calculation Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":PrecipIceptFract",p,Options.noisy); break;}
      if      (!strcmp(s[1],"USER_SPECIFIED"       )){Options.interception_factor=PRECIP_ICEPT_USER;}
      else if (!strcmp(s[1],"ICEPT_USER"           )){Options.interception_factor=PRECIP_ICEPT_USER;}
      else if (!strcmp(s[1],"LAI_LINEAR"           )){Options.interception_factor=PRECIP_ICEPT_LAI;}
      else if (!strcmp(s[1],"ICEPT_LAI"            )){Options.interception_factor=PRECIP_ICEPT_LAI;}
      else if (!strcmp(s[1],"LAI_EXPONENTIAL"      )){Options.interception_factor=PRECIP_ICEPT_EXPLAI;}
      else if (!strcmp(s[1],"ICEPT_EXPLAI"         )){Options.interception_factor=PRECIP_ICEPT_EXPLAI;}
      else if (!strcmp(s[1],"PRECIP_ICEPT_USER"    )){Options.interception_factor=PRECIP_ICEPT_USER;}
      else if (!strcmp(s[1],"PRECIP_ICEPT_LAI"     )){Options.interception_factor=PRECIP_ICEPT_LAI;}
      else if (!strcmp(s[1],"PRECIP_ICEPT_EXPLAI"  )){Options.interception_factor=PRECIP_ICEPT_EXPLAI;}
      else if (!strcmp(s[1],"PRECIP_ICEPT_HEDSTROM")){Options.interception_factor=PRECIP_ICEPT_HEDSTROM;}
      else {ExitGracefully("ParseInput:PrecipIceptFract: Unrecognized method",BAD_DATA_WARN);}
      break;
    }
    case(31): //----------------------------------------------
    {/*:OroTempCorrect" string method */
      if (Options.noisy) {cout <<"Orographic Temperature Correction Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":OroTempCorrect",p,Options.noisy); break;}
      if      (!strcmp(s[1],"HBV"                )){Options.orocorr_temp=OROCORR_HBV;}
      else if (!strcmp(s[1],"UBC"                )){Options.orocorr_temp=OROCORR_UBCWM;}
      else if (!strcmp(s[1],"UBC_2"              )){Options.orocorr_temp=OROCORR_UBCWM2;}
      else if (!strcmp(s[1],"SIMPLE"             )){Options.orocorr_temp=OROCORR_SIMPLELAPSE;}
      else if (!strcmp(s[1],"NONE"               )){Options.orocorr_temp=OROCORR_NONE;}
      else if (!strcmp(s[1],"OROCORR_HBV"        )){Options.orocorr_temp=OROCORR_HBV;}
      else if (!strcmp(s[1],"OROCORR_UBCWM"      )){Options.orocorr_temp=OROCORR_UBCWM;}
      else if (!strcmp(s[1],"OROCORR_UBCWM2"     )){Options.orocorr_temp=OROCORR_UBCWM2;}
      else if (!strcmp(s[1],"OROCORR_SIMPLELAPSE")){Options.orocorr_temp=OROCORR_SIMPLELAPSE;}
      else if (!strcmp(s[1],"OROCORR_NONE"       )){Options.orocorr_temp=OROCORR_NONE;}
      else {ExitGracefully("ParseInput:OroTempCorrect: Unrecognized method",BAD_DATA_WARN);}
      break;
    }
    case(32): //----------------------------------------------
    {/*:OroPrecipCorrect" string method */
      if (Options.noisy) {cout <<"Orographic Precipitation Correction Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":OroPrecipCorrect",p,Options.noisy); break;}
      if      (!strcmp(s[1],"HBV"                )){Options.orocorr_precip=OROCORR_HBV;}
      else if (!strcmp(s[1],"UBC"                )){Options.orocorr_precip=OROCORR_UBCWM;}
      else if (!strcmp(s[1],"UBC_2"              )){Options.orocorr_precip=OROCORR_UBCWM2;}
      else if (!strcmp(s[1],"SIMPLE"             )){Options.orocorr_precip=OROCORR_SIMPLELAPSE;}
      else if (!strcmp(s[1],"NONE"               )){Options.orocorr_precip=OROCORR_NONE;}
      else if (!strcmp(s[1],"OROCORR_HBV"        )){Options.orocorr_precip=OROCORR_HBV;}
      else if (!strcmp(s[1],"OROCORR_UBCWM"      )){Options.orocorr_precip=OROCORR_UBCWM;}
      else if (!strcmp(s[1],"OROCORR_UBCWM2"     )){Options.orocorr_precip=OROCORR_UBCWM2;}
      else if (!strcmp(s[1],"OROCORR_SIMPLELAPSE")){Options.orocorr_precip=OROCORR_SIMPLELAPSE;}
      else if (!strcmp(s[1],"OROCORR_NONE"       )){Options.orocorr_precip=OROCORR_NONE;}
      else {ExitGracefully("ParseInput:OroPrecipCorrect: Unrecognized method",BAD_DATA_WARN);}
      break;
    }
    case(33): //----------------------------------------------
    {/*:AquiferLayers" int number_of_layers */
      if (Options.noisy) {cout <<"Number of Aquifer Layers"<<endl;}
      if (Len<2){ImproperFormatWarning(":AquiferLayers",p,Options.noisy); break;}
      pModel->AddAquiferStateVars(s_to_i(s[1]));
      break;
    }
    case(34): //----------------------------------------------
    {/*:PotentialMeltMethod" string method */
      if (Options.noisy) {cout <<"Potential Melt Calculation Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":PotentialMeltMethod",p,Options.noisy); break;}
      if      (!strcmp(s[1],"POTMELT_DEGREE_DAY")){Options.pot_melt=POTMELT_DEGREE_DAY;}
      else if (!strcmp(s[1],"POTMELT_EB"        )){Options.pot_melt=POTMELT_EB;}
      else if (!strcmp(s[1],"POTMELT_RESTRICTED")){Options.pot_melt=POTMELT_RESTRICTED;}
      else if (!strcmp(s[1],"POTMELT_DD_RAIN"   )){Options.pot_melt=POTMELT_DD_RAIN;}
      else if (!strcmp(s[1],"POTMELT_UBCWM"     )){Options.pot_melt=POTMELT_UBCWM;}
      else if (!strcmp(s[1],"POTMELT_HBV"       )){Options.pot_melt=POTMELT_HBV;}
      else if (!strcmp(s[1],"POTMELT_HBV_ROS"   )){Options.pot_melt=POTMELT_HBV_ROS; }
      else if (!strcmp(s[1],"POTMELT_DATA"      )){Options.pot_melt=POTMELT_DATA;}
      else if (!strcmp(s[1],"UBC"               )){Options.pot_melt=POTMELT_UBCWM;}
      else if (!strcmp(s[1],"HBV"               )){Options.pot_melt=POTMELT_HBV;}
      else if (!strcmp(s[1],"POTMELT_USACE"     )){Options.pot_melt=POTMELT_USACE;}
      else if (!strcmp(s[1],"POTMELT_CRHM_EBSM" )){Options.pot_melt=POTMELT_CRHM_EBSM; }
      else if (!strcmp(s[1],"POTMELT_HMETS"     )){Options.pot_melt=POTMELT_HMETS; }
      else {ExitGracefully("ParseInput:PotentialMelt: Unrecognized method",BAD_DATA_WARN);}
      break;
    }
    case(35): //----------------------------------------------
    {/*:SubdailyMethod" string method */
      if (Options.noisy) {cout <<"Subdaily Downscaling Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":SubdailyMethod",p,Options.noisy); break;}
      if      (!strcmp(s[1],"NONE"                 )){Options.subdaily=SUBDAILY_NONE;}
      else if (!strcmp(s[1],"SIMPLE"               )){Options.subdaily=SUBDAILY_SIMPLE;}
      else if (!strcmp(s[1],"UBC"                  )){Options.subdaily=SUBDAILY_UBC;}
      else if (!strcmp(s[1],"SUBDAILY_NONE"        )){Options.subdaily=SUBDAILY_NONE;}
      else if (!strcmp(s[1],"SUBDAILY_SIMPLE"      )){Options.subdaily=SUBDAILY_SIMPLE;}
      else if (!strcmp(s[1],"SUBDAILY_UBC"         )){Options.subdaily=SUBDAILY_UBC;}
      else {ExitGracefully("ParseInput:SubdailyMethod: Unrecognized method",BAD_DATA_WARN);}
      break;
    }
    case(36): //----------------------------------------------
    {/*:SWCanopyCorrect" string method */
      if (Options.noisy) {cout <<"Shortwave Canopy Transmittance Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":SWCanopyCorrect",p,Options.noisy); break;}
      if      (!strcmp(s[1],"NONE"                       )){Options.SW_canopycorr=SW_CANOPY_CORR_NONE;}
      else if (!strcmp(s[1],"STATIC"                     )){Options.SW_canopycorr=SW_CANOPY_CORR_STATIC;}
      else if (!strcmp(s[1],"DYNAMIC"                    )){Options.SW_canopycorr=SW_CANOPY_CORR_DYNAMIC;}
      else if (!strcmp(s[1],"UBC"                        )){Options.SW_canopycorr=SW_CANOPY_CORR_UBCWM;}
      else if (!strcmp(s[1],"SW_CANOPY_CORR_NONE"        )){Options.SW_canopycorr=SW_CANOPY_CORR_NONE;}
      else if (!strcmp(s[1],"SW_CANOPY_CORR_STATIC"      )){Options.SW_canopycorr=SW_CANOPY_CORR_STATIC;}
      else if (!strcmp(s[1],"SW_CANOPY_CORR_DYNAMIC"     )){Options.SW_canopycorr=SW_CANOPY_CORR_DYNAMIC;}
      else if (!strcmp(s[1],"SW_CANOPY_CORR_UBCWM"       )){Options.SW_canopycorr=SW_CANOPY_CORR_UBCWM;}
      else {ExitGracefully("ParseInput:SWCanopyCorrect: Unrecognized method",BAD_DATA_WARN);}
      break;
    }
    case(37): //----------------------------------------------
    {/*:SWCloudCorrect" string method */
      if (Options.noisy) {cout <<"Shortwave Cloud Cover correction Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":SWCloudCorrect",p,Options.noisy); break;}
      if      (!strcmp(s[1],"SW_CLOUD_CORR_NONE"     )){Options.SW_cloudcovercorr=SW_CLOUD_CORR_NONE;}
      else if (!strcmp(s[1],"SW_CLOUD_CORR_UBCWM"    )){Options.SW_cloudcovercorr=SW_CLOUD_CORR_UBCWM;}
      else if (!strcmp(s[1],"SW_CLOUD_CORR_DINGMAN"  )){Options.SW_cloudcovercorr=SW_CLOUD_CORR_DINGMAN;}
      else if (!strcmp(s[1],"SW_CLOUD_CORR_ANNANDALE")){Options.SW_cloudcovercorr=SW_CLOUD_CORR_ANNANDALE;}
      else {ExitGracefully("ParseInput:SWCloudCorrect: Unrecognized method",BAD_DATA_WARN);}
      break;
    }
    case(38): //----------------------------------------------
    {
      break;
    }
    case(39):  //--------------------------------------------
    {/* :MultilayerSnow [number of layers]*/
      if (Len<2){ImproperFormatWarning(":MultilayerSnow",p,Options.noisy); break;}
      if (pModel==NULL){
        ExitGracefully(":MultilayerSnow command must be after :SoilModel command in .rvi file.",BAD_DATA_WARN); break;
      }
      pModel->SetNumSnowLayers(s_to_i(s[1]));
      break;
    }
    case(40):  //--------------------------------------------
    {/* :RetainUBCWMBugs*/
      if (Options.noisy){cout<<":RetainUBCWMBugs"<<endl;}
      Options.keepUBCWMbugs=true;
      break;
    }
    case(41):  //--------------------------------------------
    {/* :EndDate string yyyy-mm-dd hh:mm:ss.00 */
      if (Options.noisy){cout<<":EndDate"<<endl;}
      ExitGracefullyIf(Options.julian_start_year==1666,":EndDate command must be after :StartDate command in .rvi file.",BAD_DATA_WARN);
      if (Len<3){ImproperFormatWarning(":EndDate",p,Options.noisy); break;}
      time_struct tt;
      tt=DateStringToTimeStruct(s[1],s[2],Options.calendar);
      Options.duration=TimeDifference(Options.julian_start_day,Options.julian_start_year,tt.julian_day,tt.year,Options.calendar);
      ExitGracefullyIf(Options.duration<=0, "ParseInput: :EndDate must be later than :StartDate.",BAD_DATA_WARN);
      break;
    }
    case(42)://----------------------------------------------
    {/*:RechargeMethod"  string method */
      if (Options.noisy) {cout <<"Recharge Calculation Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":RechargeMethod",p,Options.noisy); break;}
      if      (!strcmp(s[1],"RECHARGE_NONE"     )){Options.recharge=RECHARGE_NONE;}
      else if (!strcmp(s[1],"RECHARGE_DATA"     )){Options.recharge=RECHARGE_DATA;}
      else {ExitGracefully("ParseInput:RechargeMethod: Unrecognized method",BAD_DATA_WARN);}
      break;
    }
    case(43)://----------------------------------------------
    {/*:NetSWRadMethod"  string method */
      if (Options.noisy) {cout <<"Net Shortwave Radiation calculation Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":NetSWRadMethod",p,Options.noisy); break;}
      if      (!strcmp(s[1],"NETSWRAD_CALC"     )){Options.SW_radia_net=NETSWRAD_CALC;}
      else if (!strcmp(s[1],"NETSWRAD_DATA"     )){Options.SW_radia_net=NETSWRAD_DATA;}
      else {ExitGracefully("ParseInput :NetSWRadMethod: Unrecognized method",BAD_DATA_WARN);}
      break;
    }
    case(44)://----------------------------------------------
    {/*:DirectEvaporation"  */
      if(Options.noisy) { cout <<"Use Direct Evaporation"<<endl; }
      Options.direct_evap=true;
      break;
    }
    case(45)://----------------------------------------------
    {/*":Calendar" string calendar  */
      if(Options.noisy) { cout <<"Change model calendar"<<endl; }
      Options.calendar=StringToCalendar(s[1]);
      break;
    }
    case(46)://----------------------------------------------
    {/*:SnowCoverDepletion"  string method */
      if (Options.noisy) {cout <<"Snow cover depletion curve method"<<endl;}
      if (Len<2){ImproperFormatWarning(":SnowCoverDepletion",p,Options.noisy); break;}
      if      (!strcmp(s[1],"SNOWCOV_NONE"       )){Options.snow_depletion=SNOWCOV_NONE;}
      else if (!strcmp(s[1],"SNOWCOV_LINEAR"     )){Options.snow_depletion=SNOWCOV_LINEAR;}
      else {ExitGracefully("ParseInput :SnowCoverDepletion: Unrecognized method",BAD_DATA_WARN);}
      break;
    }
    case(47): //----------------------------------------------
    {/*:EnsembleMode" [string method] [double #members] */
      if(Options.noisy) { cout <<"Ensemble Mode"<<endl; }
      if(Len<3) { ImproperFormatWarning(":EnsembleMode",p,Options.noisy); break; }
      if      (!strcmp(s[1],"ENSEMBLE_NONE"      )) { Options.ensemble=ENSEMBLE_NONE; }
      else if (!strcmp(s[1],"ENSEMBLE_DDS"       )) { Options.ensemble=ENSEMBLE_DDS; }
      else if (!strcmp(s[1],"ENSEMBLE_MONTECARLO")) { Options.ensemble=ENSEMBLE_MONTECARLO; }
      else { ExitGracefully("ParseInput:EnsembleMode: Unrecognized ensemble simulation mode",BAD_DATA_WARN); }
      num_ensemble_members=s_to_i(s[2]);
      break;
    }
    case(48): //---------------------------------------------
    {/*:SuppressCompetitiveET */
      if(Options.noisy) { cout <<"Suppressing competitive ET"<<endl; }
      Options.suppressCompetitiveET =true;
      break;
    }
    case(49): //---------------------------------------------
    {/*:SnowSuppressesPET */
      if(Options.noisy) { cout <<"Suppressing PET if snow on ground"<<endl; }
      Options.snow_suppressPET =true;
      break;
    }
    case(50):  //--------------------------------------------
    {/*:DebugMode */
      if (Options.noisy){cout <<"Debug Mode ON"<<endl;}
      Options.debug_mode      =true;
      Options.write_mass_bal  =true;
      Options.write_energy    =true;
      Options.write_forcings  =true;
      Options.write_channels  =true;
      Options.write_gwhead    =true;
      Options.write_gwflow    =true;
    
      WriteWarning("Debug mode is ON: this will significantly slow down model progress. Note that ':DebugMode no' command is deprecated",Options.noisy);
      break;
    }
    case(51):  //--------------------------------------------
    {/*:WriteMassBalanceFile */
      if (Options.noisy) {cout <<"Write Mass Balance File ON"<<endl;}
      Options.write_mass_bal =true;
      break;
    }
    case(52):  //--------------------------------------------
    {/*:WriteForcingFunctions */
      if (Options.noisy) {cout <<"Write Forcing Functions File ON"<<endl;}
      Options.write_forcings =true;
      break;
    }
    case(53):  //--------------------------------------------
    {/*:WriteEnergyStorage */
      if (Options.noisy) {cout <<"Write Energy Storage File ON"<<endl;}
      Options.write_energy =true;
      break;
    }
    case(54):  //--------------------------------------------
    {/*:WriteParameters */
      WriteWarning(":WriteParameters command no longer supported by Raven after v2.8.1",Options.noisy);
      break;
    }
    case(55):  //--------------------------------------------
    {/*:WriteEnsimFormat */
      bool bVal = true;
      if(Len>1)
      {
        string sVal = StringToUppercase(string(s[1]));
        if ((sVal == "NO") || (sVal == "OFF") || (sVal == "FALSE")) { bVal = false; }
        if (sVal == "PERIODENDING"){Options.period_ending=true;}
      }

      if (bVal)
      {
        if (Options.noisy){cout <<"Write Ensim Format ON"<<endl;}
        Options.output_format=OUTPUT_ENSIM;
      }
      break;
    }
    case(56):  //--------------------------------------------
    {/*:RunName */
      if (Options.noisy) {cout <<"Using Run Name: "<<s[1]<<endl;}
      if(!runname_overridden){
        Options.run_name=s[1];
      }
      else{
        WriteWarning("ParseInputFile: when run_name is specified from command line, it cannot be overridden in the .rvi file",Options.noisy);
      }
      break;
    }
    case(57):  //--------------------------------------------
    {/*:NoisyMode */
      if (Options.noisy) {cout <<"Noisy Mode!!!!"<<endl;}
      Options.noisy=true;
      break;
    }
    case(58):  //--------------------------------------------
    {/*:SilentMode */
      if (Options.noisy) {cout<<endl;}
      Options.noisy=false;Options.silent=true;
      break;
    }
    case(59):  //--------------------------------------------
    {/*:rvh_Filename */
      if (Options.noisy) {cout <<"rvh filename: "<<s[1]<<endl;}
      Options.rvh_filename=CorrectForRelativePath(s[1] ,Options.rvi_filename);
      break;
    }
    case(60):  //--------------------------------------------
    {/*:rvp_Filename */
      if (Options.noisy) {cout <<"rvp filename: "<<s[1]<<endl;}
      Options.rvp_filename=CorrectForRelativePath(s[1] ,Options.rvi_filename);//with .rvp extension!
      break; 
    }
    case(61):  //--------------------------------------------
    {/*:rvt_Filename */
      if (Options.noisy) {cout <<"rvt filename: "<<s[1]<<endl;}
      Options.rvt_filename=CorrectForRelativePath(s[1] ,Options.rvi_filename);//with .rvt extension!
      break;
    }
    case(62):  //--------------------------------------------
    {/*:rvc_Filename */
      if (Options.noisy) {cout <<"rvc filename: "<<s[1]<<endl;}
      Options.rvc_filename=CorrectForRelativePath(s[1] ,Options.rvi_filename);//with .rvc extension!
      break;
    }
    case(91):  //--------------------------------------------
    {/*:rve_Filename */
      if(Options.noisy) { cout <<"rve filename: "<<s[1]<<endl; }
      Options.rve_filename=CorrectForRelativePath(s[1],Options.rvi_filename);//with .rve extension!
      break;
    }
    case(103):  //--------------------------------------------
    {/*:rvl_Filename */
      if(Options.noisy) { cout <<"rvl filename: "<<s[1]<<endl; }
      Options.rvl_filename=CorrectForRelativePath(s[1],Options.rvi_filename);//with .rve extension!
      break;
    }
    case(63):  //--------------------------------------------
    {/*:OutputDirectory */
      if (Options.noisy) {cout <<"Output directory: "<<s[1]<<"/"<<endl;}
      Options.output_dir="";
      for (int i=1;i<Len-1;i++){
        Options.output_dir+=to_string(s[i])+"/ ";   // append backslash to make sure it's a folder
      }
      Options.output_dir+=to_string(s[Len-1])+"/";  // append backslash to make sure it's a folder
      Options.main_output_dir=Options.output_dir;
      PrepareOutputdirectory(Options);

      ofstream WARNINGS((Options.main_output_dir+"Raven_errors.txt").c_str());
      WARNINGS.close();
      break;
    }
    case(64):  //--------------------------------------------
    {/*:OutputDump YYYY-MM-DD hh:mm:ss */
      if (Options.noisy) {cout <<"Output dump @ "<<s[1]<<endl;}
      if (IsValidDateString(s[1]))
      {
        time_struct tt_out=DateStringToTimeStruct(string(s[1]),string(s[2]),Options.calendar);
        pModel->AddModelOutputTime(tt_out,Options);
      }
      else{
        ExitGracefully("ParseMainInputFile: bad timestamp format in :OutputDump command",BAD_DATA_WARN);
      }
      break;
    }
    case(65):  //--------------------------------------------
    {/*:MajorOutputInterval step[d] */
      if (Options.noisy) {cout <<"Major model output interval of "<<s[1]<<" days"<<endl;}

      time_struct tt_out;
      double tstep=s_to_d(s[1]);
      for (double t=tstep;t<Options.duration;t+=tstep)
      {
        JulianConvert(t,Options.julian_start_day,Options.julian_start_year,Options.calendar,tt_out);
        pModel->AddModelOutputTime(tt_out,Options);
      }
      break;
    }
    case(66):  //--------------------------------------------
    {/*:SnapshotHydrograph */
      if (Options.noisy) {cout <<"Snapshot hydrographs"<<endl;}
      Options.ave_hydrograph=false;
      break;
    }
    case(67):  //--------------------------------------------
    {/*:HRUStorageOutput
       ":HRUStorageOutput"  HRU_Group
     */
      if (Options.noisy) {cout <<"HRU Storage Output"<<endl;}
      if (Len<2){ImproperFormatWarning(":HRUStorageOutput",p,Options.noisy); break;}
      for (int kk=0;kk<pModel->GetNumHRUGroups();kk++)
      {
        if (!pModel->GetHRUGroup(kk)->GetName().compare(s[1])){
          pModel->SetOutputGroup(pModel->GetHRUGroup(kk));
        }
      }
      break;
    }
    case(68):  //--------------------------------------------
    {/*:SuppressWarnings */
      if (Options.noisy) {cout <<"Suppressing Warnings"<<endl;}
      g_suppress_warnings=true;
      break;
    }
    case(69):  //--------------------------------------------
    {/*:QuietMode */ //(default)
      if (Options.noisy) {cout<<endl;}
      Options.noisy=false;Options.silent=false;
      break;
    }
    case(70):  //--------------------------------------------
    {/*:WriteExhaustiveMB */ //(default)
      if (Options.noisy) {cout<<"Exhaustive Mass Balance ON"<<endl;}
      Options.write_exhaustiveMB=true;
      break;
    }
    case(71):  //--------------------------------------------
    {/*:EvaluationMetrics [diag1] {diag2}...{diagN}*/
      if (Options.noisy) {cout<<"Evaluation Metrics"<<endl;}
      CDiagnostic *pDiag=NULL;
      bool invalid;
      for (int i=1; i<Len; i++)
      {
        invalid=false;pDiag=NULL;
        int width = DOESNT_EXIST;
        string tmp = CStateVariable::SVStringBreak(s[i], width);
        if      (!strcmp(s[i],"NASH_SUTCLIFFE"     )){pDiag=new CDiagnostic(DIAG_NASH_SUTCLIFFE);}
        else if (!strcmp(s[i],"RMSE"               )){pDiag=new CDiagnostic(DIAG_RMSE);}
        else if (!strcmp(s[i],"PCT_BIAS"           )){pDiag=new CDiagnostic(DIAG_PCT_BIAS);}
        else if (!strcmp(s[i],"ABSERR"             )){pDiag=new CDiagnostic(DIAG_ABSERR);}
        else if (!strcmp(s[i],"ABSMAX"             )){pDiag=new CDiagnostic(DIAG_ABSMAX);}
        else if (!strcmp(s[i],"PDIFF"              )){pDiag=new CDiagnostic(DIAG_PDIFF);}
        else if (!strcmp(s[i],"TMVOL"              )){pDiag=new CDiagnostic(DIAG_TMVOL);}
        else if (!strcmp(s[i],"RCOEF"              )){pDiag=new CDiagnostic(DIAG_RCOEF);}
        else if (!strcmp(s[i],"NSC"                )){pDiag=new CDiagnostic(DIAG_NSC);}
        else if (!strcmp(s[i],"RSR"                )){pDiag=new CDiagnostic(DIAG_RSR);}
        else if (!strcmp(s[i],"R2"                 )){pDiag=new CDiagnostic(DIAG_R2);}
        else if (!strcmp(s[i],"CUMUL_FLOW"         )){pDiag=new CDiagnostic(DIAG_CUMUL_FLOW);}
        else if (!strcmp(s[i],"LOG_NASH"           )){pDiag=new CDiagnostic(DIAG_LOG_NASH);}
        else if (!strcmp(s[i],"KLING_GUPTA"        )){pDiag=new CDiagnostic(DIAG_KLING_GUPTA);}
        else if (!strcmp(s[i],"NASH_SUTCLIFFE_DER" )){pDiag=new CDiagnostic(DIAG_NASH_SUTCLIFFE_DER);}
        else if (!strcmp(s[i],"RMSE_DER"           )){pDiag=new CDiagnostic(DIAG_RMSE_DER);}
        else if (!strcmp(s[i],"KLING_GUPTA_DER"    )){pDiag=new CDiagnostic(DIAG_KLING_GUPTA_DER);}
        else if (!strcmp(s[i],"MBF"                )){pDiag=new CDiagnostic(DIAG_MBF); }
        else if (!tmp.compare("NASH_SUTCLIFFE_RUN")) {pDiag=new CDiagnostic(DIAG_NASH_SUTCLIFFE_RUN, width); }
        else   {invalid=true;}
        if (!invalid){
          pModel->AddDiagnostic(pDiag);
        }
        else{
          WriteWarning("Invalid diagnostic ("+to_string(s[i])+") found in .rvi file",Options.noisy);
        }
      }
      break;
    }
    case (72): //--------------------------------------------
    {/*:SuppressOutputICs */
      Options.suppressICs=true;
      break;
    }
    case(73):  //--------------------------------------------
    {/*:SuppressOutput */
      if (Options.noisy){cout <<"Suppressing output"<<endl;}
      Options.output_format=OUTPUT_NONE;
      break;
    }
    case(74):  //--------------------------------------------
    {/*:WriteHRUGroupMBFile [HRU Group name]*/
      if (Options.noisy) {cout <<"Write HRU Group Mass Balance File ON"<<endl;}
      CHRUGroup *pHRUGroup;
      pHRUGroup = pModel->GetHRUGroup(s[1]);
      if (pHRUGroup != NULL){
        Options.write_group_mb = pModel->GetHRUGroup(s[1])->GetGlobalIndex();
      }
      else{
        WriteWarning("ParseMainInput: invalid HRU group specified in :WriteHRUGroupMBFile command. Please define groups using :DefineHRUGroups command prior to calling this command.",Options.noisy);
      }
      break;
    }          
    case(75):  //--------------------------------------------
    {/*:EvaluationTime [yyyy-mm-dd] [00:00:00] {yyyy-mm-dd} {00:00:00}*/ //AFTER StartDate or JulianStartDay and JulianStartYear commands
      if (Options.noisy) { cout << "Evaluation Time" << endl; }
      if (Len<3) { ImproperFormatWarning(":EvaluationTime", p, Options.noisy); break; }

      time_struct tt;
      tt = DateStringToTimeStruct(s[1], s[2], Options.calendar);
      Options.diag_start_time = TimeDifference(Options.julian_start_day,Options.julian_start_year,tt.julian_day, tt.year, Options.calendar);
      if (Len >= 5) // optional diagnostic end time
      {
        tt = DateStringToTimeStruct(s[3], s[4], Options.calendar);
        Options.diag_end_time = TimeDifference(Options.julian_start_day,Options.julian_start_year,tt.julian_day, tt.year, Options.calendar);
      }
      break;
    }
    case(76):  //--------------------------------------------
    {/*:WaterYearStartMonth [int month]*/
      if (Options.noisy) {cout <<"Water year starting month"<<endl;}
      int mo=s_to_i(s[1]);
      if((mo>0) && (mo<=12)){Options.wateryr_mo=mo;}
      else { WriteWarning("Invalid water year starting month in :WaterYearStartMonth command. Should be integer month between 1- and 12",Options.noisy); }
      break;
    }
    case(77):  //--------------------------------------------
    {/*:CreateRVPTemplate*/
      if (Options.noisy) {cout <<"Create RVP Template File"<<endl;}
      Options.create_rvp_template=true;
      break;
    }
    case(78):  //--------------------------------------------
    {/*:WriteNetCDFFormat */
      if (Options.noisy){cout <<"Write NetCDF Format ON"<<endl;}
      Options.output_format=OUTPUT_NETCDF;
      break;
    }
    case(79):  //--------------------------------------------
    {/*:WriteChannelInfo */
      if (Options.noisy) {cout <<"Write Channel Info file ON"<<endl;}
      Options.write_channels =true;
      break;
    }
    case(80):  //--------------------------------------------
    {/*:DefineHRUGroup */ //AFTER SoilModel Command and HydroProcesses commands
      if (Options.noisy) {cout <<"Defining HRU Group"<<endl;}
      if (Len<2){ImproperFormatWarning(":DefineHRUGroup",p,Options.noisy); break;}
      if(pModel==NULL) {
        ExitGracefully(":DefineHRUGroup command must be after :SoilModel command in .rvi file.",BAD_DATA_WARN); break;
      }
      CHRUGroup *pHRUGrp=NULL;
      pHRUGrp=new CHRUGroup(s[1],pModel->GetNumHRUGroups());
      pModel->AddHRUGroup(pHRUGrp);
      break;
    }
    case(81):  //--------------------------------------------
    {/*:DefineHRUGroups */ //AFTER SoilModel Command and HydroProcesses commands
      if (Options.noisy) {cout <<"Defining HRU Groups"<<endl;}
      if (Len<2){ImproperFormatWarning(":DefineHRUGroups",p,Options.noisy); break;}
      if(pModel==NULL) {
        ExitGracefully(":DefineHRUGroups command must be after :SoilModel command in .rvi file.",BAD_DATA_WARN); break;
      }
      CHRUGroup *pHRUGrp=NULL;
      for (int i=1;i<Len;i++){
        pHRUGrp=new CHRUGroup(s[i],pModel->GetNumHRUGroups());
        pModel->AddHRUGroup(pHRUGrp);
      }
      break;
    }
    case(82):  //--------------------------------------------
    {/*:DisableHRUGroup */ //AFTER DefineHRUGroup(s) commands 
      if (Options.noisy) {cout <<"Disabling HRU Group"<<endl;}
      if (Len<2){ImproperFormatWarning(":DisableHRUGroup",p,Options.noisy); break;}
      CHRUGroup *pHRUGrp=NULL;

      pHRUGrp=pModel->GetHRUGroup(s[1]);
      if (pHRUGrp==NULL){
        ExitGracefully("Invalid HRU Group name supplied in :DisableHRUGroup command in .rvi file",BAD_DATA_WARN);
        break;
      }
      else{
        pHRUGrp->DisableGroup();
      }
      break;
    }
    case(85):  //--------------------------------------------
    {/*:BenchmarkingMode */
      if (Options.noisy) {cout <<"Benchmarking Mode"<<endl;}
      Options.benchmarking=true;
      g_suppress_zeros=true;
      break;
    }
    case(86):  //--------------------------------------------
    {/*:WriteReservoirMBFile*/
      if (Options.noisy) {cout <<"Write Reservoir Mass Balance File ON"<<endl;}
      Options.write_reservoirMB=true;
      break;
    }
    case(87):  //--------------------------------------------
    {/*:PeriodStartingFormatOff*/
      if (Options.noisy) {cout <<"Backward compatible to version 2.7 output"<<endl;}
      Options.period_starting=false;
      break;
    }
    case(88):  //--------------------------------------------
    {/*:PavicsMode */
      if(Options.noisy) { cout<<endl; }
      Options.pavics=true;
      ofstream PROGRESS;
      PROGRESS.open((Options.main_output_dir+"Raven_progress.txt").c_str());
      if(PROGRESS.fail()) {
        ExitGracefully("ParseInput:: Unable to open Raven_progress.txt. Bad output directory specified?",RUNTIME_ERR);
      }
      PROGRESS.close();
      break;
    }
    case(89):  //--------------------------------------------
    {/*:OutputConstituentMass*/
      if (Options.noisy) {cout <<"Write constituent mass instead of concentrations"<<endl;}
      Options.write_constitmass=true;
      break;
    }
    case(90):  //--------------------------------------------
    {/*:NetCDFAttribute [attrib_name] [value (commas not allowed)]*/
      if (Options.noisy) {cout <<"Add NetCDF attribute"<<endl;}
      if (Len<3){ImproperFormatWarning(":NetCDFAttribute",p,Options.noisy); break;}
      string tmpstring=s[2];
      for(i=3;i<Len;i++){tmpstring=tmpstring+" "+s[i];}
      AddNetCDFAttribute(Options,s[1],tmpstring);
      break;
    }          
    case(92):  //--------------------------------------------
    {/*:AssimilationStartTime [yyyy-mm-dd] [00:00:00]*///AFTER StartDate or JulianStartDay and JulianStartYear commands
      if(Options.noisy) { cout << "Assimilation Start Time" << endl; }
      if(Len<3) { ImproperFormatWarning(":AssimilationStartTime",p,Options.noisy); break; }

      time_struct tt;
      tt = DateStringToTimeStruct(s[1],s[2],Options.calendar);
      Options.assimilation_start=TimeDifference(Options.julian_start_day,Options.julian_start_year,tt.julian_day,tt.year,Options.calendar);

      if(Options.assimilation_start>Options.duration) {
        string warn="Data assimilation start date is after model simulation end. No data assimilation will be applied.";
        WriteWarning(warn,Options.noisy);
      }
      break;
    }
    case(93):  //--------------------------------------------
    {/*:AssimilateStreamflow*/
      if(Options.noisy) { cout << "Assimilate streamflow on" << endl; }
      Options.assimilation_on=true;
      break;
    }
    case(94):  //--------------------------------------------
    {/*:DeltaresFEWSMode*/
      if(Options.noisy) { cout << "Deltares FEWS input ingestion mode" << endl; }
      Options.deltaresFEWS=true;
      break;
    }
    case(95):  //--------------------------------------------
    {/*:WriteSubbasinFile*/
      if(Options.noisy) { cout << "Write Subbasin file" << endl; }
      Options.write_basinfile=true;
      break;
    }
    case(96):  //--------------------------------------------
    {/*:DontWriteWatershedStorage */
      if(Options.noisy) { cout <<"Write WatershedStorage OFF"<<endl; }
      Options.write_watershed_storage=false;
      break;
    }
    case(97):  //--------------------------------------------
    {/*:TimeZone*/
      if(Options.noisy) { cout << "Set Time Zone" << endl; }
      Options.time_zone=s_to_i(s[1]);
      break;
    }
    case(98):  //--------------------------------------------
    {/*:Alias */
      if (Options.noisy) {cout <<"Alias"<<endl;}
      if (Len<3){ImproperFormatWarning(":Alias",p,Options.noisy); break;}
      CStateVariable::AddAlias(s[1],s[2]);
      break;
    }
    case(99):  //----------------------------------------------
    {/*:CustomOutput
       :CustomOutput [time_aggregation] [statistic] [parameter] [space_aggregation] {ONLY HRUGroup} {[hist_min] [hist_max] [#bins]} {filename} (optional)
     */

      if (Options.noisy) {cout <<"Custom Output"<<endl;}
      if (Len<5){ImproperFormatWarning(":CustomOutput",p,Options.noisy); break;}
      time_agg    ta;
      spatial_agg sa;
      agg_stat    stat;
      string force_str="";

      if      (!strcmp(s[1],"DAILY"          )){ta=DAILY;}
      else if (!strcmp(s[1],"MONTHLY"        )){ta=MONTHLY;}
      else if (!strcmp(s[1],"YEARLY"         )){ta=YEARLY;}
      else if (!strcmp(s[1],"ANNUAL"         )){ta=YEARLY;}
      else if (!strcmp(s[1],"WATER_YEARLY"   )){ta=WATER_YEARLY;}
      else if (!strcmp(s[1],"EVERY_NDAYS"    )){ta=EVERY_NDAYS; }
      else if (!strcmp(s[1],"CONTINUOUS"     )){ta=EVERY_TSTEP;}
      else{
        ta=DAILY;
        ExitGracefully("ParseMainInputFile: Unrecognized custom output temporal aggregation method",BAD_DATA);
      }
      if((ta==EVERY_NDAYS) && (Options.custom_interval==1.0)) {
        WriteWarning(":CustomOutput - EVERY_NDAYS option should only be used with :CustomInterval > 1, otherwise use of DAILY command preferred",Options.noisy);
      }

      //these statistics are always in time
      if      (!strcmp(s[2],"AVERAGE"         )){stat=AGG_AVERAGE;}
      else if (!strcmp(s[2],"MAXIMUM"         )){stat=AGG_MAXIMUM;}
      else if (!strcmp(s[2],"MINIMUM"         )){stat=AGG_MINIMUM;}
      else if (!strcmp(s[2],"MEDIAN"          )){stat=AGG_MEDIAN;}
      else if (!strcmp(s[2],"RANGE"           )){stat=AGG_RANGE;}
      else if (!strcmp(s[2],"HISTOGRAM"       )){stat=AGG_HISTOGRAM;}
      else if (!strcmp(s[2],"QUARTILES"       )){stat=AGG_QUARTILES;}
      else if (!strcmp(s[2],"95CI"            )){stat=AGG_95CI;}
      //
      else{
        stat=AGG_AVERAGE;
        ExitGracefully("ParseMainInputFile: Unrecognized custom output processing method",BAD_DATA);
      }

      //read in parameter information
      diagnostic diag;
      diag=VAR_STATE_VAR ;//For now, default
      int SV_ind,SV_ind2;
      sv_type sv_typ;
      SV_ind= ParseSVTypeIndex(s[3], pModel);
      SV_ind2=DOESNT_EXIST;

      //Special treatment of To:, From: and Between:1.And.2 fluxes
      string tmp = s[3];
      string right;
      if (!strcmp((tmp.substr(0, 3)).c_str(), "To:")){
        right = tmp.substr(3,string::npos);

        SV_ind=ParseSVTypeIndex(right, pModel);
        if (SV_ind == DOESNT_EXIST){
          WriteWarning("Custom output Flux variable " + right + " is unrecognized. No output will be written.", Options.noisy);
          break;
        }
        else
        {
          diag=VAR_TO_FLUX;
        }
      }
      if (!strcmp((tmp.substr(0, 5)).c_str(), "From:")){
        right = tmp.substr(5,string::npos);
        SV_ind=ParseSVTypeIndex(right, pModel);
        if (SV_ind == DOESNT_EXIST){
          WriteWarning("Custom output Flux variable " + right + " is unrecognized. No output will be written.", Options.noisy);
          break;
        }
        else
        {
          diag=VAR_FROM_FLUX;
        }
      }
      //e.g., :Between SOIL[0].And.ATMOSPHERE for AET
      if (!strcmp((tmp.substr(0, 8)).c_str(), "Between:")){
        right = tmp.substr(8,string::npos);
        string firstSV = right.substr(0,right.find(".And."));
        string lastSV  = right.substr(right.find(".And.")+5,string::npos);

        SV_ind =ParseSVTypeIndex(firstSV, pModel);
        SV_ind2=ParseSVTypeIndex(lastSV,  pModel);
        if (SV_ind == DOESNT_EXIST){
          WriteWarning("Custom output Flux variable " + firstSV + " is unrecognized. No output will be written.", Options.noisy);
          break;
        }
        else if (SV_ind2 == DOESNT_EXIST){
          WriteWarning("Custom output Flux variable " + lastSV + " is unrecognized. No output will be written.", Options.noisy);
          break;
        }
        else
        {
          diag=VAR_BETWEEN_FLUX;
        }
      }
      //not a state variable or a flux to/from state var - try forcing function
      if (SV_ind==DOESNT_EXIST){
        diag=VAR_FORCING_FUNCTION;
        sv_typ=UNRECOGNIZED_SVTYPE;
        force_str=s[3];
        if (GetForcingTypeFromString(force_str) == F_UNRECOGNIZED){
          WriteWarning("Custom output variable " + force_str + " is unrecognized. No output will be written.", Options.noisy);
          break;
        }
      }
      else{
        sv_typ=pModel->GetStateVarType(SV_ind);
      }


      if      (!strcmp(s[4],"BY_HRU"          )){sa=BY_HRU;}
      else if (!strcmp(s[4],"BY_BASIN"        )){sa=BY_BASIN;}
      else if (!strcmp(s[4],"BY_SUBBASIN"     )){sa=BY_BASIN;}
      else if (!strcmp(s[4],"ENTIRE_WATERSHED")){sa=BY_WSHED;}
      else if (!strcmp(s[4],"BY_WATERSHED"    )){sa=BY_WSHED;}
      else if (!strcmp(s[4],"BY_GROUP"        )){sa=BY_HRU_GROUP;}
      else if (!strcmp(s[4],"BY_HRU_GROUP"    )){sa=BY_HRU_GROUP;}
      else{
        sa=BY_HRU;
        ExitGracefully("ParseMainInputFile: Unrecognized custom output spatial aggregation method",BAD_DATA);
      }
      int kk_only=DOESNT_EXIST;
      string HRU_Group="";
      if ((sa==BY_HRU) && (Len>=7) && (string(s[5])=="ONLY")){
        sa=BY_SELECT_HRUS;
        HRU_Group=s[6];
        for (int kk=0;kk<pModel->GetNumHRUGroups();kk++){
          if (pModel->GetHRUGroup(kk)->GetName()==HRU_Group){
            kk_only=kk;
          }
        }
      }

      // get custom filename, if specified
      int start=5;
      if ((sa==BY_SELECT_HRUS) && (Len>=7) && (string(s[5])=="ONLY")){start=7;}
      else if ((Len>=8) && (stat==AGG_HISTOGRAM))                    {start=8;}
      string filename="";
      if (Len>start){for (int i=start;i<Len;i++){filename=filename+to_string(s[i]);}}

      CCustomOutput *pCustom;
      pCustom=new CCustomOutput(diag,sv_typ,SV_ind,SV_ind2,force_str,stat,ta,sa,filename,kk_only,pModel,Options);
      pModel->AddCustomOutput(pCustom);

      if ((Len>=8) && (stat==AGG_HISTOGRAM)){
        pCustom->SetHistogramParams(s_to_d(s[5]),s_to_d(s[6]), s_to_i(s[7]));
      }

      break;
    }
    case(100):  //--------------------------------------------
    {/*List of Included State Variables
       string ":StateVariables", int NumStateVars
       {string VarType} x NumStateVars
       &
     */
      if (Options.noisy) {cout <<"State Variables (OBSOLETE)"<<endl;}
      break;
    }
    case(101):  //--------------------------------------------
    {/*:WriteInterpolationWeights*/
      if(Options.noisy) { cout << "Write Interpolation Weights file" << endl; }
      Options.write_interp_wts=true;
      break;
    }
    case(102):  //--------------------------------------------
    {/*:CallExternalScript*/
      if(Options.noisy) { cout << "Call External Script" << endl; }
      for(int i=1;i<Len;i++) {
        Options.external_script+=s[i];
        if(i!=Len-1) { Options.external_script+=" "; } //assumes space delimiter - commas will cause issues
      }
      if(!system(NULL)) {
        WriteWarning("Unable to call external script on this machine - :CallExternalScript command will be ignored.",Options.noisy);
        Options.external_script="";
      }
      break;
    }
    case(104): //----------------------------------------------
    {/*:ReadLiveFile {optional frequency, d or hh:mm:ss}*/
      if(Options.noisy) { cout <<"Read Live File "<<endl; }
      Options.rvl_read_frequency=Options.timestep; //default if not overrun - read every timestep
      if(Len>=2) {
        string tString=s[1];
        if((tString.length()>=2) && 
          ((tString.substr(2,1)==":") || 
           (tString.substr(1,1)==":"))) //support for hh:mm:ss.00 format
        {
          time_struct tt=DateStringToTimeStruct("0000-01-01",tString,Options.calendar);
          Options.rvl_read_frequency=FixTimestep(tt.julian_day);
        }
        else {
          Options.rvl_read_frequency = s_to_d(s[1]);
        }
      }
      break;
    }
    case(105):
    {/*:CustomOutputInterval [interval, in days]*/
      if(Options.noisy) { cout <<"Custom output interval "<<endl; }
      if(Len<2) { ImproperFormatWarning(":CustomOutputInterval",p,Options.noisy); break; }
      Options.custom_interval=s_to_d(s[1]);
      if(Options.custom_interval==1.0) {
        WriteWarning("CustomOutputInterval: using interval of 1 day is equivalent to use the DAILY aggregation of custom output, which is preferred",Options.noisy);
      }
      if((Options.custom_interval<1.0) || (floor(Options.custom_interval)!=Options.custom_interval)) {
        ExitGracefully(":CustomOutputInterval: improper interval specified. Must be >0 and integer",BAD_DATA_WARN);
      }
      break;
    }
    case(106):
    {/*:UseStopFile*/
      if(Options.noisy) { cout <<"Use stopfile "<<endl; }
      Options.use_stopfile=true;
      break;
    }
    case(107):  //--------------------------------------------
    {/*:WriteDemandFile*/
      if(Options.noisy) { cout << "Write Irrigation Demand file" << endl; }
      Options.write_demandfile=true;
      break;
    }
    case(150):  //--------------------------------------------
    {/*:AggregatedVariable
       ":AggregatedVariable" [SV_TAG] {optional HRU_Group}
     */
      if (Options.noisy) {cout <<"Aggregated Variable"<<endl;}
      if (Len<2){ImproperFormatWarning(":AggregatedVariable",p,Options.noisy); break;}
      tmpS[0]=CStateVariable::StringToSVType(s[1],tmpLev[0],true);
      string group_name="ALL";
      if (Len==3){group_name=s[2];}
      pModel->SetAggregatedVariable(tmpS[0],tmpLev[0],group_name);
      break;
    }
    case(200):  //----------------------------------------------
    {/*HydrologicProcesses
       string ":HydrologicProcesses" */
      if (Options.noisy){cout <<"Begin Hydrologic Process List"<<endl;}
      break;
    }
    case(201):  //----------------------------------------------
    {/*Baseflow
       :Baseflow [string method] [int from_index] [SURFACE_WATER] */
      if (Options.noisy){cout <<"Baseflow Process"<<endl;}
      baseflow_type btype=BASE_LINEAR;
      if (Len<4){ImproperFormatWarning(":Baseflow",p,Options.noisy); break;}

      if      (!strcmp(s[1],"BASE_VIC"             )){btype=BASE_VIC;}
      else if (!strcmp(s[1],"BASE_TOPMODEL"        )){btype=BASE_TOPMODEL;}
      else if (!strcmp(s[1],"BASE_LINEAR"          )){btype=BASE_LINEAR;}
      else if (!strcmp(s[1],"BASE_LINEAR_CONSTRAIN")){btype=BASE_LINEAR_CONSTRAIN;}
      else if (!strcmp(s[1],"BASE_LINEAR_ANALYTIC" )){btype=BASE_LINEAR_ANALYTIC;}
      else if (!strcmp(s[1],"BASE_POWER_LAW"       )){btype=BASE_POWER_LAW;}
      else if (!strcmp(s[1],"BASE_CONSTANT"        )){btype=BASE_CONSTANT;}
      else if (!strcmp(s[1],"BASE_THRESH_POWER"    )){btype=BASE_THRESH_POWER;}
      else if (!strcmp(s[1],"BASE_THRESH_STOR"     )){btype=BASE_THRESH_STOR;}
      else if (!strcmp(s[1],"BASE_GR4J"            )){btype=BASE_GR4J;}
      else {
        ExitGracefully("ParseMainInputFile: Unrecognized baseflow process representation",BAD_DATA_WARN); break;
      }
      CmvBaseflow::GetParticipatingStateVarList(btype,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);
      tmpS[0]=CStateVariable::StringToSVType(s[2],tmpLev[0],true);
      pModel->AddStateVariables(tmpS,tmpLev,1);

      pMover=new CmvBaseflow(btype,ParseSVTypeIndex(s[2],pModel));
      AddProcess(pModel,pMover,pProcGroup);

      break;
    }
    case(202):  //----------------------------------------------
    {/*Canopy Evaporation
       :CanopyEvaporation [string method] CANOPY ATMOSPHERE*/
      if (Options.noisy){cout <<"Canopy Evaporation Process"<<endl;}
      canevap_type ce_type=CANEVP_RUTTER;
      if (Len<4){ImproperFormatWarning(":CanopyEvaporation",p,Options.noisy); break;}
      if      (!strcmp(s[1],"CANEVP_RUTTER"  )){ce_type=CANEVP_RUTTER;}
      else if (!strcmp(s[1],"CANEVP_MAXIMUM" )){ce_type=CANEVP_MAXIMUM;}
      else if (!strcmp(s[1],"CANEVP_ALL"     )){ce_type=CANEVP_ALL;}
      else
      {
        ExitGracefully("ParseMainInputFile: Unrecognized canopy evaporation process representation",BAD_DATA_WARN); break;
      }
      CmvCanopyEvap::GetParticipatingStateVarList(ce_type,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);

      pMover=new CmvCanopyEvap(ce_type);
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(203):  //----------------------------------------------
    {/*Canopy Drip
       :CanopyDrip [string method] [CANOPY] [int to_index]*/
      if (Options.noisy){cout <<"Canopy Drip Process"<<endl;}
      if (Len<4){ImproperFormatWarning(":CanopyDrip",p,Options.noisy); break;}
      candrip_type ctype=CANDRIP_RUTTER;
      if      (!strcmp(s[1],"CANDRIP_RUTTER"   )){ctype=CANDRIP_RUTTER;}
      else if (!strcmp(s[1],"CANDRIP_SLOWDRAIN")){ctype=CANDRIP_SLOWDRAIN;}
      else {
        ExitGracefully("ParseMainInputFile: Unrecognized canopy drip process representation",BAD_DATA_WARN); break;
      }
      CmvCanopyDrip::GetParticipatingStateVarList(ctype,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);
      tmpS[0]=CStateVariable::StringToSVType(s[3],tmpLev[0],true);
      pModel->AddStateVariables(tmpS,tmpLev,1);

      pMover=new CmvCanopyDrip(ctype,ParseSVTypeIndex(s[3],pModel));
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(204):  //----------------------------------------------
    {/*:Infiltration
       :Infiltration [string method] PONDED_WATER MULTIPLE*/
      if (Options.noisy){cout <<"Infiltration Process"<<endl;}
      infil_type itype=INF_RATIONAL;
      if (Len<4){ImproperFormatWarning(":Infiltration",p,Options.noisy); break;}
      if      (!strcmp(s[1],"INF_GREEN_AMPT"  )){itype=INF_GREEN_AMPT; }
      else if (!strcmp(s[1],"INF_GA_SIMPLE"   )){itype=INF_GA_SIMPLE; }
      else if (!strcmp(s[1],"INF_VIC_ARNO"    )){itype=INF_VIC_ARNO;  }
      else if (!strcmp(s[1],"INF_VIC"         )){itype=INF_VIC;       }
      else if (!strcmp(s[1],"INF_RATIONAL"    )){itype=INF_RATIONAL;  }
      else if (!strcmp(s[1],"INF_PRMS"        )){itype=INF_PRMS;      }
      else if (!strcmp(s[1],"INF_HBV"         )){itype=INF_HBV;       }
      else if (!strcmp(s[1],"INF_UBC"         )){itype=INF_UBC;       }
      else if (!strcmp(s[1],"INF_PARTITION"   )){itype=INF_RATIONAL;  }
      else if (!strcmp(s[1],"INF_GR4J"        )){itype=INF_GR4J;      }
      else if (!strcmp(s[1],"INF_SCS"         )){itype=INF_SCS;       }
      else if (!strcmp(s[1],"INF_SCS_NOABSTRACTION")){itype=INF_SCS_NOABSTRACTION;  }
      else if (!strcmp(s[1],"INF_HMETS"       )){itype=INF_HMETS;     }
      else if (!strcmp(s[1],"INF_ALL_INFILTRATES")) { itype=INF_ALL_INFILTRATES; }
      else if (!strcmp(s[1],"INF_XINANXIANG"  )) { itype=INF_XINANXIANG; }
      else {
        ExitGracefully("ParseMainInputFile: Unrecognized infiltration process representation",BAD_DATA_WARN); break;
      }
      CmvInfiltration::GetParticipatingStateVarList(itype,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);

      pMover=new CmvInfiltration(itype);
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(205):  //----------------------------------------------
    {/*:Percolation
       :Percolation [string method] [int from_index] [int to_index]*/
      if (Options.noisy){cout <<"Percolation Process"<<endl;}
      perc_type p_type=PERC_CONSTANT;
      if (Len<4){ImproperFormatWarning(":Percolation",p,Options.noisy); break;}
      if      (!strcmp(s[1],"PERC_POWER_LAW"  )){p_type=PERC_POWER_LAW;}
      else if (!strcmp(s[1],"POWER_LAW"       )){p_type=PERC_POWER_LAW;}
      else if (!strcmp(s[1],"PERC_GAWSER"     )){p_type=PERC_GAWSER;}
      else if (!strcmp(s[1],"PERC_GAWSER_CONSTRAIN")){p_type=PERC_GAWSER_CONSTRAIN;}
      else if (!strcmp(s[1],"PERC_PRMS"       )){p_type=PERC_PRMS;}
      else if (!strcmp(s[1],"PERC_SACRAMENTO" )){p_type=PERC_SACRAMENTO;}
      else if (!strcmp(s[1],"PERC_CONSTANT"   )){p_type=PERC_CONSTANT;}
      else if (!strcmp(s[1],"PERC_LINEAR"     )){p_type=PERC_LINEAR;}
      else if (!strcmp(s[1],"PERC_LINEAR_ANALYTIC")){p_type=PERC_LINEAR_ANALYTIC;}
      else if (!strcmp(s[1],"PERC_GR4J"       )){p_type=PERC_GR4J;}
      else if (!strcmp(s[1],"PERC_GR4JEXCH"   )){p_type=PERC_GR4JEXCH;}
      else if (!strcmp(s[1],"PERC_GR4JEXCH2"  )){p_type=PERC_GR4JEXCH2;}
      else if (!strcmp(s[1],"PERC_ASPEN"      )){p_type=PERC_ASPEN;}
      else {
        ExitGracefully("ParseMainInputFile: Unrecognized percolation process representation",BAD_DATA_WARN); break;
      }
      CmvPercolation::GetParticipatingStateVarList(p_type,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);
      tmpS[0]=CStateVariable::StringToSVType(s[2],tmpLev[0],true);
      tmpS[1]=CStateVariable::StringToSVType(s[3],tmpLev[1],true);
      pModel->AddStateVariables(tmpS,tmpLev,2);

      pMover=new CmvPercolation(p_type,
                                ParseSVTypeIndex(s[2],pModel),
                                ParseSVTypeIndex(s[3],pModel));
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(206):  //----------------------------------------------
    {/*SnowMelt [NOW OBSOLETE]
       :SnowMelt [string method] [SNOW] [int to_index]*/
      if (Options.noisy){cout <<"Snow Melt Process"<<endl;}
      snowmelt_type stype=MELT_POTMELT;
      if (Len<4){ImproperFormatWarning(":SnowMelt",p,Options.noisy); break;}
      if      (!strcmp(s[1],"MELT_POTMELT"           )){stype=MELT_POTMELT;}
      else {
        ExitGracefully("ParseMainInputFile: Unrecognized snowmelt process representation",BAD_DATA_WARN); break;
      }
      CmvSnowMelt::GetParticipatingStateVarList(stype,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);
      tmpS[0]=CStateVariable::StringToSVType(s[3],tmpLev[0],true);
      pModel->AddStateVariables(tmpS,tmpLev,1);

      pMover=new CmvSnowMelt(stype,
                             ParseSVTypeIndex(s[3],pModel));
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(207):  //----------------------------------------------
    {
      break;
    }
    case(208):  //----------------------------------------------
    {/*Soil Evaporation
       :SoilEvaporation [string method] MULTIPLE ATMOSPHERE
     */
      if (Options.noisy){cout <<"Soil Evaporation Process"<<endl;}
      soilevap_type se_type=SOILEVAP_VIC;
      if (Len<4){ImproperFormatWarning(":SoilEvaporation",p,Options.noisy); break;}
      if      (!strcmp(s[1],"SOILEVAP_VIC"          )){se_type=SOILEVAP_VIC;}
      else if (!strcmp(s[1],"SOILEVAP_TOPMODEL"     )){se_type=SOILEVAP_TOPMODEL;}
      else if (!strcmp(s[1],"SOILEVAP_SEQUEN"       )){se_type=SOILEVAP_SEQUEN;}
      else if (!strcmp(s[1],"SOILEVAP_ROOT"         )){se_type=SOILEVAP_ROOT;}
      else if (!strcmp(s[1],"SOILEVAP_ROOT_CONSTRAIN")){se_type=SOILEVAP_ROOT_CONSTRAIN;}
      else if (!strcmp(s[1],"SOILEVAP_HBV"          )){se_type=SOILEVAP_HBV;}
      else if (!strcmp(s[1],"SOILEVAP_HYPR"         )){se_type=SOILEVAP_HYPR;}
      else if (!strcmp(s[1],"SOILEVAP_UBC"          )){se_type=SOILEVAP_UBC;}
      else if (!strcmp(s[1],"SOILEVAP_CHU"          )){se_type=SOILEVAP_CHU;}
      else if (!strcmp(s[1],"SOILEVAP_GR4J"         )){se_type=SOILEVAP_GR4J;}
      else if (!strcmp(s[1],"SOILEVAP_LINEAR"       )){se_type=SOILEVAP_LINEAR;}
      else if (!strcmp(s[1],"SOILEVAP_ALL"          )){se_type=SOILEVAP_ALL;}
      else {
        ExitGracefully("ParseMainInputFile: Unrecognized soil evaporation process representation",BAD_DATA_WARN); break;
      }
      CmvSoilEvap::GetParticipatingStateVarList(se_type,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);

      pMover=new CmvSoilEvap(se_type);
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(209):  //----------------------------------------------
    {/*Snow Balance
       :SnowBalance [string method]    MULTIPLE MULTIPLE
       :SnowBalance SNOBAL_SIMPLE_MELT SNOW   PONDED_WATER/SNOW_LIQ
     */
      if (Options.noisy){cout <<"Snow Balance (Melt/Refreeze) Processes"<<endl;}
      snowbal_type sbtype=SNOBAL_COLD_CONTENT;
      if (Len<4){ImproperFormatWarning(":SnowBalance",p,Options.noisy); break;}

      if      (!strcmp(s[1],"SNOBAL_COLD_CONTENT" )){sbtype=SNOBAL_COLD_CONTENT;}
      else if (!strcmp(s[1],"SNOBAL_SIMPLE_MELT"  )){sbtype=SNOBAL_SIMPLE_MELT;}
      else if (!strcmp(s[1],"UBC"                 )){sbtype=SNOBAL_UBCWM;}      // backward compatible
      else if (!strcmp(s[1],"SNOBAL_UBCWM"        )){sbtype=SNOBAL_UBCWM;}
      else if (!strcmp(s[1],"SNOBAL_HBV"          )){sbtype=SNOBAL_HBV;}
      else if (!strcmp(s[1],"SNOBAL_CEMA_NIEGE"   )){sbtype=SNOBAL_CEMA_NEIGE;} // backward compatible
      else if (!strcmp(s[1],"SNOBAL_CEMA_NEIGE"   )){sbtype=SNOBAL_CEMA_NEIGE;}
      else if (!strcmp(s[1],"SNOBAL_TWO_LAYER"    )){sbtype=SNOBAL_TWO_LAYER;}
      else if (!strcmp(s[1],"SNOBAL_GAWSER"       )){sbtype=SNOBAL_GAWSER;}
      else if (!strcmp(s[1],"SNOBAL_CRHM_EBSM"    )){sbtype=SNOBAL_CRHM_EBSM;}
      else if (!strcmp(s[1],"SNOBAL_HMETS"        )){sbtype=SNOBAL_HMETS;}
      else {
        ExitGracefully("ParseMainInputFile: Unrecognized snow balance process representation",BAD_DATA_WARN); break;
      }
      CmvSnowBalance::GetParticipatingStateVarList(sbtype,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);

      if (sbtype == SNOBAL_SIMPLE_MELT) {
        tmpS[0]=CStateVariable::StringToSVType(s[3],tmpLev[0],true);
        pModel->AddStateVariables(tmpS,tmpLev,1);
        pMover=new CmvSnowBalance(sbtype, ParseSVTypeIndex(s[3],pModel)); 
      }
      else{
        pMover=new CmvSnowBalance(sbtype);
      }
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(210):  //----------------------------------------------
    {/*Sublimation
       :Sublimation [string method] SNOW ATMOSPHERE*/
      if (Options.noisy){cout <<"Sublimation Process"<<endl;}
      sublimation_type sub_type=SUBLIM_SVERDRUP;
      if (Len<4){ImproperFormatWarning(":Sublimation",p,Options.noisy); break;}
      if      (!strcmp(s[1],"SUBLIM_SVERDRUP"        )){sub_type=SUBLIM_SVERDRUP;}
      else if (!strcmp(s[1],"SUBLIM_KUZMIN"          )){sub_type=SUBLIM_KUZMIN;}
      else if (!strcmp(s[1],"SUBLIM_CENTRAL_SIERRA"  )){sub_type=SUBLIM_CENTRAL_SIERRA;}
      else if (!strcmp(s[1],"SUBLIM_PBSM"            )){sub_type=SUBLIM_PBSM;}
      else if (!strcmp(s[1],"SUBLIM_WILLIAMS"        )){sub_type=SUBLIM_WILLIAMS;}
      else if (!strcmp(s[1],"SUBLIM_CRHM"            )){sub_type=SUBLIM_CRHM; }
      else {
        ExitGracefully("ParseMainInputFile: Unrecognized sublimation process representation",BAD_DATA_WARN); break;
      }
      CmvSublimation::GetParticipatingStateVarList(sub_type,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);

      pMover=new CmvSublimation(sub_type);
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(211):  //----------------------------------------------
    {/*OpenWaterEvaporation
       :OpenWaterEvaporation [string method] [PONDED_WATER/DEPRESSION] ATMOSPHERE */
      if (Options.noisy){cout <<"Open Water Evaporation Process"<<endl;}
      owevap_type ow_type=OPEN_WATER_EVAP;
      if (Len<4){ImproperFormatWarning(":OpenWaterEvaporation",p,Options.noisy); break;}
      if      (!strcmp(s[1],"OPEN_WATER_EVAP"     )){ow_type=OPEN_WATER_EVAP;}
      else {
        ExitGracefully("ParseMainInputFile: Unrecognized Open Water Evaporation process representation",BAD_DATA_WARN); break;
      }
      CmvOWEvaporation::GetParticipatingStateVarList(ow_type,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);

      tmpS[0]=CStateVariable::StringToSVType(s[2],tmpLev[0],true);
      pModel->AddStateVariables(tmpS,tmpLev,1);

      pMover=new CmvOWEvaporation(ow_type,ParseSVTypeIndex(s[2],pModel));
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(212):  //----------------------------------------------
    {/*Precipitation
       :Precipitation PRECIP_RAVEN ATMOS_PRECIP MULTIPLE */
      if (Options.noisy){cout <<"Precipitation Process"<<endl;}
      if (Len<4){ImproperFormatWarning(":Precipitation",p,Options.noisy); break;}
      CmvPrecipitation::GetParticipatingStateVarList(tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);

      pMover=pPrecip=new CmvPrecipitation();
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(213):  //----------------------------------------------
    {/*Interflow
       :Interflow [string method] [int from_index] SURFACE_WATER */
      if (Options.noisy){cout <<"Interflow Process"<<endl;}
      interflow_type inttype=INTERFLOW_PRMS;
      if (Len<4){ImproperFormatWarning(":Interflow",p,Options.noisy); break;}
      if        (!strcmp(s[1],"INTERFLOW_PRMS"      )){inttype=INTERFLOW_PRMS;}
      else {
        ExitGracefully("ParseMainInputFile: Unrecognized interflow process representation",BAD_DATA_WARN); break;
      }
      CmvInterflow::GetParticipatingStateVarList(inttype,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);
      tmpS[0]=CStateVariable::StringToSVType(s[2],tmpLev[0],true);
      pModel->AddStateVariables(tmpS,tmpLev,1);

      pMover=new CmvInterflow(inttype,ParseSVTypeIndex(s[2],pModel));
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(214):  //----------------------------------------------
    {/*Snow Refreeze
       :SnowRefreeze [string method] SNOW_LIQ SNOW*/
      if (Options.noisy){cout <<"Snow Refreeze Process"<<endl;}
      refreeze_type rtype=FREEZE_DEGREE_DAY;
      if (Len<4){ImproperFormatWarning(":SnowRefreeze",p,Options.noisy); break;}
      if      (!strcmp(s[1],"FREEZE_DEGREE_DAY"  )){rtype=FREEZE_DEGREE_DAY;}
      else {
        ExitGracefully("ParseMainInputFile: Unrecognized snow refreeze process representation",BAD_DATA_WARN); break;
      }
      CmvSnowRefreeze::GetParticipatingStateVarList(rtype,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);

      pMover=new CmvSnowRefreeze(rtype);
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(215):  //----------------------------------------------
    {/*Flush
       :Flush RAVEN_DEFAULT [state_var from] [state_var to]*/
      if (Options.noisy){cout <<"Flushing Process"<<endl;}

      if (Len<4){ImproperFormatWarning(":Flush",p,Options.noisy); break;}
      tmpS[0]=CStateVariable::StringToSVType(s[2],tmpLev[0],true);
      tmpS[1]=CStateVariable::StringToSVType(s[3],tmpLev[1],true);
      pModel->AddStateVariables(tmpS,tmpLev,2);

      pMover=new CmvFlush(ParseSVTypeIndex(s[2],pModel),
                          ParseSVTypeIndex(s[3],pModel));
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(216):  //----------------------------------------------
    {/*:Capillary Rise
       :CapillaryRise [string method] [state_var from] [state_var to]*/
      if (Options.noisy){cout <<"Capillary Rise Process"<<endl;}
      if (Len<4){ImproperFormatWarning(":CapillaryRise",p,Options.noisy); break;}
      crise_type ctype=CRISE_HBV;
      if      (!strcmp(s[1],"RISE_HBV"   )){ctype=CRISE_HBV;}
      else if (!strcmp(s[1],"CRISE_HBV"  )){ctype=CRISE_HBV;}
      else {
        ExitGracefully("ParseMainInputFile: Unrecognized capillary rise process representation",BAD_DATA_WARN); break;
      }
      CmvCapillaryRise::GetParticipatingStateVarList(ctype,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);
      tmpS[0]=CStateVariable::StringToSVType(s[2],tmpLev[0],true);
      tmpS[1]=CStateVariable::StringToSVType(s[3],tmpLev[1],true);
      pModel->AddStateVariables(tmpS,tmpLev,2);

      pMover=new CmvCapillaryRise(ctype,
                                  ParseSVTypeIndex(s[2],pModel),
                                  ParseSVTypeIndex(s[3],pModel));
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(217):  //----------------------------------------------
    {/*LakeEvaporation
       :LakeEvaporation [string method] [state_var from] ATMOSPHERE*/
      if (Options.noisy){cout <<"Lake Evaporation Process"<<endl;}
      lakeevap_type lk_type=LAKE_EVAP_BASIC;
      if (Len<4){ImproperFormatWarning(":LakeEvaporation",p,Options.noisy); break;}
      if      (!strcmp(s[1],"BASIC"           )){lk_type=LAKE_EVAP_BASIC;}
      else if (!strcmp(s[1],"LAKE_EVAP_BASIC" )){lk_type=LAKE_EVAP_BASIC;}
      else {
        ExitGracefully("ParseMainInputFile: Unrecognized Lake Evaporation process representation",BAD_DATA_WARN); break;
      }
      CmvLakeEvaporation::GetParticipatingStateVarList(lk_type,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);
      int lake_ind;
      if (Len==3){
        tmpS[0]=CStateVariable::StringToSVType(s[2],tmpLev[0],true);
        pModel->AddStateVariables(tmpS,tmpLev,1);
        lake_ind=ParseSVTypeIndex(s[2],pModel);
      }
      else{
        lake_ind=pModel->GetLakeStorageIndex();
      }
      pMover=new CmvLakeEvaporation(lk_type,lake_ind);
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(218):  //----------------------------------------------
    {/*SnowSqueeze
       :SnowSqueeze SQUEEZE_RAVEN SNOW_LIQ [state_var to_index]*/
      if (Options.noisy){cout <<"Liquid Snow Release Process"<<endl;}
      if (Len<4){ImproperFormatWarning(":SnowSqueeze",p,Options.noisy); break;}

      CmvSnowSqueeze::GetParticipatingStateVarList(tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);
      tmpS[0]=CStateVariable::StringToSVType(s[3],tmpLev[0],true);
      pModel->AddStateVariables(tmpS,tmpLev,1);

      pMover=new CmvSnowSqueeze(ParseSVTypeIndex(s[3],pModel));
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(219):  //----------------------------------------------
    {/*GlacierMelt
       :GlacierMelt [string method] GLACIER_ICE GLACIER */
      if (Options.noisy){cout <<"Glacier Melt Process"<<endl;}
      if (Len<4){ImproperFormatWarning(":GlacierMelt",p,Options.noisy); break;}
      glacial_melt_type gm_type=GMELT_SIMPLE_MELT;
      if      (!strcmp(s[1],"HBV"                  )){gm_type=GMELT_HBV;}
      else if (!strcmp(s[1],"GMELT_HBV"            )){gm_type=GMELT_HBV;}
      else if (!strcmp(s[1],"UBC"                  )){gm_type=GMELT_UBC;}
      else if (!strcmp(s[1],"GMELT_UBC"            )){gm_type=GMELT_UBC;}
      else if (!strcmp(s[1],"SIMPLE"               )){gm_type=GMELT_SIMPLE_MELT;}
      else if (!strcmp(s[1],"GMELT_SIMPLE_MELT"    )){gm_type=GMELT_SIMPLE_MELT;}
      else {
        ExitGracefully("ParseMainInputFile: Unrecognized Glacier Melt process representation",BAD_DATA_WARN); break;
      }
      CmvGlacierMelt::GetParticipatingStateVarList(gm_type,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);

      pMover=new CmvGlacierMelt(gm_type);
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(220):  //----------------------------------------------
    {/*GlacialRelease
       :GlacierRelease [string method] GLACIER SURFACE_WATER*/
      if (Options.noisy){cout <<"Glacial Release Process"<<endl;}
      if (Len<4){ImproperFormatWarning(":GlacierRelease",p,Options.noisy); break;}
      glacial_release_type gm_type=GRELEASE_HBV_EC;
      if      (!strcmp(s[1],"HBV_EC"                  )){gm_type=GRELEASE_HBV_EC;}
      else if (!strcmp(s[1],"GRELEASE_HBV_EC"         )){gm_type=GRELEASE_HBV_EC;}
      else if (!strcmp(s[1],"LINEAR_STORAGE"          )){gm_type=GRELEASE_LINEAR;}
      else if (!strcmp(s[1],"GRELEASE_LINEAR"         )){gm_type=GRELEASE_LINEAR;}
      else if (!strcmp(s[1],"GRELEASE_LINEAR_ANALYTIC")){gm_type = GRELEASE_LINEAR_ANALYTIC;}
      else {
        ExitGracefully("ParseMainInputFile: Unrecognized Glacier Release process representation",BAD_DATA_WARN); break;
      }
      CmvGlacierRelease::GetParticipatingStateVarList(gm_type,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);

      pMover=new CmvGlacierRelease(gm_type);
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(221):  //----------------------------------------------
    {/*Canopy Snow Evaporation
       :CanopySnowEvaporation [string method] CANOPY_SNOW ATMOSPHERE
     */
      if (Options.noisy){cout <<"Canopy Snow Evaporation Process"<<endl;}
      canevap_type ce_type=CANEVP_ALL;
      if (Len<4){ImproperFormatWarning(":CanopySnowEvaporation",p,Options.noisy); break;}
      if      (!strcmp(s[1],"MAXIMUM" )){ce_type=CANEVP_MAXIMUM;}
      else if (!strcmp(s[1],"HBV"     )){ce_type=CANEVP_ALL;}
      else if (!strcmp(s[1],"ALL"     )){ce_type=CANEVP_ALL;}
      else if (!strcmp(s[1],"CANEVP_RUTTER"  )){ce_type=CANEVP_RUTTER;}
      else if (!strcmp(s[1],"CANEVP_MAXIMUM" )){ce_type=CANEVP_MAXIMUM;}
      else if (!strcmp(s[1],"CANEVP_ALL"     )){ce_type=CANEVP_ALL;}
      else
      {
        ExitGracefully("ParseMainInputFile: Unrecognized canopy snow evaporation process representation",BAD_DATA_WARN); break;
      }
      CmvCanopySnowEvap::GetParticipatingStateVarList(ce_type,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);

      pMover=new CmvCanopySnowEvap(ce_type);
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(222):  //----------------------------------------------
    {/*Overflow
       :Overflow OVERFLOW_RAVEN [sv from] [sv to]*/
      if (Options.noisy){cout <<"Overflow Process"<<endl;}
      if (Len<4){ImproperFormatWarning(":Overflow",p,Options.noisy); break;}

      tmpS[0]=CStateVariable::StringToSVType(s[2],tmpLev[0],true);
      tmpS[1]=CStateVariable::StringToSVType(s[3],tmpLev[1],true);
      pModel->AddStateVariables(tmpS,tmpLev,2);

      pMover=new CmvOverflow(ParseSVTypeIndex(s[2],pModel),
                             ParseSVTypeIndex(s[3],pModel));

      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(223):  //----------------------------------------------
    {/*Snow Albedo Evolution
       :SnowAlbedoEvolve [string method]*/
      if (Options.noisy){cout <<"Snow Albedo Evolution Process"<<endl;}
      snowalb_type snalb_type=SNOALB_UBCWM;
      if (Len<2){ImproperFormatWarning(":SnowAlbedoEvolve",p,Options.noisy); break;}
      if      (!strcmp(s[1],"UBC"                )){snalb_type=SNOALB_UBCWM;}
      else if (!strcmp(s[1],"SNOALB_UBCWM"       )){snalb_type=SNOALB_UBCWM;}
      else if (!strcmp(s[1],"SNOALB_CRHM_ESSERY" )){snalb_type=SNOALB_CRHM_ESSERY; }
      else if (!strcmp(s[1],"SNOALB_BAKER"       )){snalb_type=SNOALB_BAKER; }
      else
      {
        string message="ParseMainInputFile: Unrecognized snow albedo algorithm "+string(s[1]);
        ExitGracefully(message.c_str(),BAD_DATA_WARN); break;
      }
      CmvSnowAlbedoEvolve::GetParticipatingStateVarList(snalb_type,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);

      pMover=new CmvSnowAlbedoEvolve(snalb_type);
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(224):  //----------------------------------------------
    {/*Crop Heat Unit Evolution
       :CropHeatUnitEvolve [string method]*/
      if (Options.noisy){cout <<"Crop Heat Unit Evolution Process"<<endl;}
      CHUevolve_type CHU_type=CHU_ONTARIO;
      if (Len<2){ImproperFormatWarning(":CropHeatUnitEvolve",p,Options.noisy); break;}
      if      (!strcmp(s[1],"CHU_ONTARIO"  )){CHU_type=CHU_ONTARIO;}
      else
      {
        ExitGracefully("ParseMainInputFile: Unrecognized crop heat evolution algorithm",BAD_DATA_WARN); break;
      }
      CmvCropHeatUnitEvolve::GetParticipatingStateVarList(CHU_type,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);

      pMover=new CmvCropHeatUnitEvolve(CHU_type);
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }

    case(225):  //----------------------------------------------
    {/*Abstraction of rainfall/snowmelt
       :Abstraction [string method] PONDED_WATER DEPRESSION/MULTIPLE*/
      if (Options.noisy){cout <<"Rainfall/Snowmelt Abstraction Process"<<endl;}
      abstraction_type abst_type=ABST_FILL;
      if (Len<4){ImproperFormatWarning(":Abstraction",p,Options.noisy); break;}
      if      (!strcmp(s[1],"ABST_SCS"         )){abst_type=ABST_SCS;}
      else if (!strcmp(s[1],"ABST_PERCENTAGE"  )){abst_type=ABST_PERCENTAGE;}
      else if (!strcmp(s[1],"ABST_FILL"        )){abst_type=ABST_FILL;}
      else if (!strcmp(s[1],"ABST_PDMROF"      )){abst_type=ABST_PDMROF; }
      else
      {
        ExitGracefully("ParseMainInputFile: Unrecognized abstraction algorithm",BAD_DATA_WARN); break;
      }
      CmvAbstraction::GetParticipatingStateVarList(abst_type,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);

      pMover=new CmvAbstraction(abst_type);
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }

    case(226):  //----------------------------------------------
    {/*GlacierInfiltration
       :GlacierInfiltration [string method] PONDED_WATER MULTIPLE*/
      if (Options.noisy){cout <<"Glacial Infiltration Process"<<endl;}
      if (Len<4){ImproperFormatWarning(":GlacierInfiltration",p,Options.noisy); break;}
      glacial_infil_type gi_type=GINFIL_UBCWM;
      if      (!strcmp(s[1],"GINFIL_UBCWM"    )){gi_type=GINFIL_UBCWM;}
      else {
        ExitGracefully("ParseMainInputFile: Unrecognized Glacier Infiltration process representation",BAD_DATA_WARN); break;
      }
      CmvGlacierInfil::GetParticipatingStateVarList(gi_type,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);

      pMover=new CmvGlacierInfil(gi_type);
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }

    case(227):  //----------------------------------------------
    {/*Split
       :Split RAVEN_DEFAULT [from] [to1] [to2] [split pct]*/
      if (Options.noisy){cout <<"Split Process"<<endl;}
      if (Len<6){ImproperFormatWarning(":Split",p,Options.noisy); break;}

      tmpS[0]=CStateVariable::StringToSVType(s[2],tmpLev[0],true);
      tmpS[1]=CStateVariable::StringToSVType(s[3],tmpLev[1],true);
      tmpS[2]=CStateVariable::StringToSVType(s[4],tmpLev[2],true);

      pModel->AddStateVariables(tmpS,tmpLev,3);

      pMover=new CmvSplit(ParseSVTypeIndex(s[2],pModel),
                          ParseSVTypeIndex(s[3],pModel),
                          ParseSVTypeIndex(s[4],pModel),
                          s_to_d(s[5]));
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(228):  //----------------------------------------------
    {/*Convolve
       :Convolve [string method] CONVOLUTION[i] [TARGET]*/
      if (Options.noisy){cout <<"Convolution Process"<<endl;}
      convolution_type c_type=CONVOL_GR4J_1;
      if (Len<4){ImproperFormatWarning(":Convolve",p,Options.noisy); break;}

      if      (!strcmp(s[1],"CONVOL_GR4J_1"  )){c_type=CONVOL_GR4J_1;}
      else if (!strcmp(s[1],"CONVOL_GR4J_2"  )){c_type=CONVOL_GR4J_2;}
      else if (!strcmp(s[1],"CONVOL_GAMMA"   )){c_type=CONVOL_GAMMA;}
      else if (!strcmp(s[1],"CONVOL_GAMMA_2"  )){c_type=CONVOL_GAMMA_2;}
      else
      {
        ExitGracefully("ParseMainInputFile: Unrecognized convolution process representation",BAD_DATA_WARN); break;
      }
      CmvConvolution::GetParticipatingStateVarList(c_type,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);

      pMover=new CmvConvolution(c_type,ParseSVTypeIndex(s[3],pModel));
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(229):  //----------------------------------------------
    {/*SnowTempEvolve
       :SnowTempEvolve [string method] SNOW_TEMP*/
      if (Options.noisy){cout <<"Snow temperature evolution Process"<<endl;}
      snowtemp_evolve_type ste_type=SNOTEMP_NEWTONS;
      if (Len<3){ImproperFormatWarning(":SnowTempEvolve",p,Options.noisy); break;}
      if      (!strcmp(s[1],"SNOTEMP_NEWTONS"  )){ste_type=SNOTEMP_NEWTONS;}
      else
      {
        ExitGracefully("ParseMainInputFile: Unrecognized snow temp evolution process representation",BAD_DATA_WARN); break;
      }
      CmvSnowTempEvolve::GetParticipatingStateVarList(ste_type,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);

      pMover=new CmvSnowTempEvolve(ste_type);
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(230):  //----------------------------------------------
    {/*Overflow of depression/wetland storage
       :DepressionOverflow [string method] DEPRESSION SURFACE_WATER*/
      if (Options.noisy){cout <<"Overflow of depression/wetland storage process"<<endl;}
      depflow_type d_type=DFLOW_THRESHPOW;
      if (Len<4){ImproperFormatWarning(":DepressionOverflow",p,Options.noisy); break;}
      if      (!strcmp(s[1],"DFLOW_THRESHPOW"   )){d_type=DFLOW_THRESHPOW;}
      else if (!strcmp(s[1],"DFLOW_LINEAR"      )){d_type=DFLOW_LINEAR;}
      else if (!strcmp(s[1],"DFLOW_WEIR"        )){d_type=DFLOW_WEIR;}
      else
      {
        ExitGracefully("ParseMainInputFile: Unrecognized depression overflow algorithm",BAD_DATA_WARN); break;
      }
      CmvDepressionOverflow::GetParticipatingStateVarList(d_type,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);

      pMover=new CmvDepressionOverflow(d_type);
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(231):  //----------------------------------------------
    {/*Exchange flow with mixing zone
       :ExchangeFlow RAVEN_DEFAULT [state_var from] [state_var mixing_zone]*/
      if (Options.noisy){cout <<"Exchange flow with mixing zone Process"<<endl;}

      if (Len<4){ImproperFormatWarning(":ExchangeFlow",p,Options.noisy); break;}
      tmpS[0]=CStateVariable::StringToSVType(s[2],tmpLev[0],true);
      tmpS[1]=CStateVariable::StringToSVType(s[3],tmpLev[1],true);
      pModel->AddStateVariables(tmpS,tmpLev,2);

      pMover=new CmvExchangeFlow(ParseSVTypeIndex(s[2],pModel),
                                 ParseSVTypeIndex(s[3],pModel));
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(232):  //----------------------------------------------
    {/*Lateral Flush
       :LateralFlush RAVEN_DEFAULT [FROM_HRUGroup] [FROM SV] To [TO_HRUGROUP] [TO_SV] {'INTERBASIN'} */
      if(Options.noisy){ cout <<"Lateral Flush Process"<<endl; }
      bool interbasin=false;
      if(Len<7){ ImproperFormatWarning(":LateralFlush",p,Options.noisy); break; }

      tmpS[0]=CStateVariable::StringToSVType(s[3],tmpLev[0],true);
      tmpS[1]=CStateVariable::StringToSVType(s[6],tmpLev[1],true);
      pModel->AddStateVariables(tmpS,tmpLev,2);

      if ((Len>=8) && (!strcmp(s[7],"INTERBASIN"))){interbasin=true; }
 
      if((pModel->GetHRUGroup(s[2])==NULL) || (pModel->GetHRUGroup(s[2])==NULL)){
        ExitGracefully("ParseInput: Lateral Flush - invalid 'to' or 'from' HRU Group used. Must define using :DefineHRUGroups command.",BAD_DATA_WARN);
      }
      else{
        pMover=new CmvLatFlush( pModel->GetStateVarIndex(tmpS[0],tmpLev[0]),//from SV index
                                pModel->GetStateVarIndex(tmpS[1],tmpLev[1]),//to SV index
                                pModel->GetHRUGroup(s[2])->GetGlobalIndex(),
                                pModel->GetHRUGroup(s[5])->GetGlobalIndex(),
                                !interbasin);
        AddProcess(pModel,pMover,pProcGroup);
      }
      break;
    }
    case(233):  //----------------------------------------------
    {/*Seepage from depression/wetland storage to soil
       :Seepage [string method] DEPRESSION {SOIL[?]}*/
      if (Options.noisy){cout <<"Seepage from depression/wetland process"<<endl;}
      seepage_type s_type=SEEP_LINEAR;
      if (Len<4){ImproperFormatWarning(":Seepage",p,Options.noisy); break;}
      if      (!strcmp(s[1],"SEEP_LINEAR"      )){s_type=SEEP_LINEAR;}
      else
      {
        ExitGracefully("ParseMainInputFile: Unrecognized seepage algorithm",BAD_DATA_WARN); break;
      }
      CmvSeepage::GetParticipatingStateVarList(s_type,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);

      pMover=new CmvSeepage(s_type,ParseSVTypeIndex(s[3],pModel));
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(234):  //----------------------------------------------
    {/*recharge from other model or from deep groundwater to lower soil storage
       :Recharge RECHARGE_FROMFILE ATMOS_PRECIP SOIL[?]*/
      if(Options.noisy) { cout <<"Recharge process"<<endl; }
      if(Len<4) { ImproperFormatWarning(":Recharge",p,Options.noisy); break; }
      recharge_type rech_type=RECHARGE_FROMFILE;

      if     (!strcmp(s[1],"RECHARGE_CONSTANT"        )) { rech_type=RECHARGE_CONSTANT; }
      else if(!strcmp(s[1],"RECHARGE_FROMFILE"        )) { rech_type=RECHARGE_FROMFILE; }
      else if(!strcmp(s[1],"RAVEN_DEFAULT"            )) { rech_type=RECHARGE_FROMFILE; }
      else if(!strcmp(s[1],"RECHARGE_CONSTANT_OVERLAP")) { rech_type=RECHARGE_CONSTANT_OVERLAP; }
      else {
        ExitGracefully("ParseMainInputFile: Unrecognized recharge process representation",BAD_DATA);
      }
      CmvRecharge::GetParticipatingStateVarList(rech_type,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);

      tmpS[0]=CStateVariable::StringToSVType(s[3],tmpLev[0],true);
      pModel->AddStateVariables(tmpS,tmpLev,1);

      if(rech_type==RECHARGE_FROMFILE){
        pMover=new CmvRecharge(ParseSVTypeIndex(s[3],pModel));
      }
      else {
        int Conns=1;
        if     (Len == 2) { Conns = 1; }
        else if(Len == 3) { Conns = s_to_i(s[2]); } //GWMIGRATE - not sure what is happening here.
        pMover=new CmvRecharge(rech_type,Conns);
      }
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }

    case(235):  //----------------------------------------------
    {/*Blowing snow redistribution & Sublimation
       :BlowingSnow PBSM MULTIPLE MULTIPLE*/
      if (Options.noisy){cout <<"Blowing Snow process"<<endl;}
      if (Len<4){ImproperFormatWarning(":BlowingSnow",p,Options.noisy); break;}

      CmvPrairieBlowingSnow::GetParticipatingStateVarList(PBSM_FULL,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);

      pMover=new CmvPrairieBlowingSnow(PBSM_FULL);
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(236):  //----------------------------------------------
    {/*Lake Release 
       :LakeRelease LAKEREL_LINEAR LAKE SURFACE_WATER*/
      if (Options.noisy){cout <<"Lake Release process"<<endl;}
      if (Len<4){ImproperFormatWarning(":LakeRelease",p,Options.noisy); break;}
      lakerel_type l_type=LAKEREL_LINEAR;
      if      (!strcmp(s[1],"LAKEREL_LINEAR"      )){l_type=LAKEREL_LINEAR;}
      else
      {
        ExitGracefully("ParseMainInputFile: Unrecognized lake release algorithm",BAD_DATA_WARN); break;
      }
      CmvLakeRelease::GetParticipatingStateVarList(l_type,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);

      pMover=new CmvLakeRelease(l_type);
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(237):  //----------------------------------------------
    {/*:-->RedirectFlow
     :-->RedirectFlow  [old 'To' SV] [new 'To' SV]*/
      if(Options.noisy) { cout <<"Redirecting Flow"<<endl; }
      if(Len<3) { ImproperFormatWarning(":-->RedirectFlow",p,Options.noisy); break; }

      tmpS[0]=CStateVariable::StringToSVType(s[2],tmpLev[0],true);
      pModel->AddStateVariables(tmpS,tmpLev,1);

      // \todo [QA\QC]: check if this is a process that can support this (infiltration can support runoff redirects)

      pMover->Redirect(ParseSVTypeIndex(s[1],pModel),ParseSVTypeIndex(s[2],pModel));

      break;
    }
    case(252):  //----------------------------------------------
    {/*:Drain
      string ":Drain" string method */
      if (Options.noisy){cout <<"Drain Process"<<endl;}
        drain_type d_type=DRAIN_CONDUCTANCE;
      if (Len<2){ImproperFormatWarning(":Drain",p,Options.noisy); break;}
      if      (!strcmp(s[1],"CONDUCTANCE"        )){d_type=DRAIN_CONDUCTANCE;}
      else if (!strcmp(s[1],"UFR"                )){d_type=DRAIN_UFR;}
      else {
        ExitGracefully("ParseMainInputFile: Unrecognized drain process representation",BAD_DATA);
      }
      CmvDrain::GetParticipatingStateVarList(d_type,tmpS,tmpLev,tmpN);
      pModel->AddStateVariables(tmpS,tmpLev,tmpN);
        
      pMover=new CmvDrain(d_type);
      pModel->AddProcess(pMover);
      break;
    }
    case (295)://----------------------------------------------
    {/*ProcessGroup
       string ":ProcessGroup" [name]*/
      if(Options.noisy) { cout <<"Process Group Start"<<endl; }
      if(Len<1) { ImproperFormatWarning(":ProcessGroup",p,Options.noisy); break; }
      else {
        pProcGroupOuter=pProcGroup;
        pProcGroup=new CProcessGroup(s[1]);
      }
      break;
    }
    case (296)://----------------------------------------------
    {/*End ProcessGroup
       string ":EndProcessGroup {wt1 wt2 wt3 ... wtN}"
       or
       string ":EndProcessGroup CALCULATE_WTS {r1 r2 r3 ... rN-1}"*/
      if (Options.noisy){cout <<"Process Group End"<<endl;}
      
      int N=pProcGroup->GetGroupSize();
      double *aWts=new double[N];
      if(Len==N+1) {
        if(!strcmp(s[1],"CALCULATE_WTS")) { 
          for(i=0;i<N-1;i++) { aWts[i]=s_to_d(s[i+2]); }
          pProcGroup->CalcWeightsFromUniformNums(aWts,N-1);
        }
        else {    
          for(i=0;i<N;i++) { aWts[i]=s_to_d(s[i+1]); }
          pProcGroup->SetWeights(aWts,N);
        }
        delete aWts;
      }
      else if(Len>1) {
        WriteWarning("ParseMainInputFile: incorrect number of weights in :EndProcessGroup command. Contents ignored.",Options.noisy);
      }
      pModel->AddProcess(pProcGroup);
      pProcGroup=NULL;
      if(pProcGroupOuter!=NULL) {
        pProcGroup=pProcGroupOuter;
        pProcGroupOuter=NULL;
      }
      break;
    }
    case (297)://----------------------------------------------
    {/*->Conditional
       string ":-->Conditional" [basis_of_condition] [condition] [criterion]*/
      if (Options.noisy){cout <<"Hydro Process Conditional Statement"<<endl;}
      if (Len<4){ImproperFormatWarning(":-->Conditional",p,Options.noisy); break;}
      condition_basis basis;
      comparison       cond;
      if      (!strcmp(s[1],"HRU_TYPE"     )){basis=BASIS_HRU_TYPE;  }
      else if (!strcmp(s[1],"HRU_GROUP"    )){basis=BASIS_HRU_GROUP; }
      else if (!strcmp(s[1],"LAND_CLASS"   )){basis=BASIS_LANDCLASS; }
      else if (!strcmp(s[1],"VEGETATION"   )){basis=BASIS_VEGETATION; }
      else                                   {
        ExitGracefully("ParseMainInputFile: Conditional statement has invalid basis",BAD_DATA_WARN);break;

      }
      if      (!strcmp(s[2],"IS"           )){cond=COMPARE_IS_EQUAL;  }
      else if (!strcmp(s[2],"IS_NOT"       )){cond=COMPARE_NOT_EQUAL; }
      else{
        ExitGracefully("ParseMainInputFile: Conditional statement has invalid condition",BAD_DATA_WARN);break;
      }

      if (pMover!=NULL){
        pMover->AddCondition(basis,cond,s[3]);
      }
      else{
        ExitGracefully("ParseMainInputFile: Cannot add hydrological process condition without corresponding process",BAD_DATA_WARN);
      }
      break;
    }

    case (298)://----------------------------------------------
    {/* :EndHydrologicProcesses */
      if (Options.noisy){cout <<"End Hydrologic Process List"<<endl;}

      //Dynamically generate connections used in precipitation - requires specification of all other processes (and corresponding water SVs)
      if (pPrecip!=NULL){pPrecip->Initialize();}

      //must be done after processes are initialized (partition precip, in particular), but before constituents are added
      if (Options.noisy){ cout<<"Preparing Transport Model..."<<endl;}
      pModel->GetTransportModel()->Prepare(Options);
      transprepared=true;
      if (Options.noisy){ cout<<"  ...prepared"<<endl;}
      break;
    }

    case (299)://----------------------------------------------
    {/*->Cascade
       :-->Cascade [int start_variable] [next_variable] [next_variable]...*/
      if (Options.noisy){cout <<"Cascade"<<endl;}
      if (Len<=1){ImproperFormatWarning(":-->Cascade",p,Options.noisy); break;}
      int CascInds[MAX_STATE_VARS];

      for (i=1; i<Len;i++){//.add variables to model (if they arent there)
        tmpS    [i-1]=CStateVariable::StringToSVType  (s[i],tmpLev[i-1],true);
      }
      pModel->AddStateVariables(tmpS,tmpLev,Len-1);
      if (pMover!=NULL){
        for (i=1; i<Len;i++){
          CascInds[i-1]=ParseSVTypeIndex(s[i],pModel);
        }
        pMover->AddCascade  (CascInds,   Len-1);
      }
      else{
        ExitGracefully("ParseMainInputFile: Cannot add cascade without corresponding process",BAD_DATA_WARN);
      }
      break;
    }

    case (300)://----------------------------------------------
    {/*:Transport
       :Transport [string constituent_name] {PASSIVE} {TRACER}*/
      if (Options.noisy){cout <<"Transport constituent"<<endl;}

      ExitGracefullyIf(!transprepared,
                       ":Transport command must be after :EndHydrologicProcesses command in .rvi file", BAD_DATA);
      bool is_passive(false),is_tracer(false);
      if(Len>=3) {
        if(!strcmp(s[2],"PASSIVE")) { is_passive=true; }
        if(!strcmp(s[2],"TRACER" )) { is_tracer =true; }
      }
      if(Len>=4) {
        if(!strcmp(s[3],"PASSIVE")) { is_passive=true; }
        if(!strcmp(s[3],"TRACER" )) { is_tracer =true; }
      }
      constit_type ctype;
      if (is_tracer){ctype=TRACER;}
      else          {ctype=AQUEOUS;}
      if(!strcmp(s[2],"TEMPERATURE")) { ctype=ENTHALPY; }
      pModel->GetTransportModel()->AddConstituent(s[1],ctype,is_passive);

      pMover=new CmvAdvection(s[1],pModel->GetTransportModel());
      AddProcess(pModel,pMover,pProcGroup);

      pMover=new CmvLatAdvection(s[1],pModel->GetTransportModel());
      AddProcess(pModel,pMover,pProcGroup);

      if(ctype==ENTHALPY) {//add precipitation source condition, by default - Tprecip=Tair
        int iAtmPrecip=pModel->GetStateVarIndex(ATMOS_PRECIP);
        pModel->GetTransportModel()->AddDirichletCompartment(s[1],iAtmPrecip,DOESNT_EXIST,DIRICHLET_AIR_TEMP);
      }
      break;
    }

    case (301)://----------------------------------------------
    {/*:FixedConcentration
       :FixedConcentration [string constit_name] [string state_var (storage compartment)] [double concentration (mg/l) or .rvt file name] {optional HRU Group name}*/
      if (Options.noisy){cout <<"Fixed concentration transport constituent"<<endl;}

      if (!transprepared){
        ExitGracefully(":FixedConcentration command must be after :Transport command in .rvi file", BAD_DATA_WARN);
        break;
      }
      int layer_ind;
      int i_stor;
      sv_type typ=CStateVariable::StringToSVType(s[2],layer_ind,false);
      if (typ==UNRECOGNIZED_SVTYPE){
        WriteWarning(":FixedConcentration command: unrecognized storage variable name: "+to_string(s[2]),Options.noisy);
        break;
      }
      i_stor=pModel->GetStateVarIndex(typ,layer_ind);
      if (i_stor != DOESNT_EXIST){
        int kk = DOESNT_EXIST;
        if (Len > 4){
          CHRUGroup *pSourceGrp;
          pSourceGrp = pModel->GetHRUGroup(s[4]);
          if (pSourceGrp == NULL){
            ExitGracefully("Invalid HRU Group name supplied in :FixedConcentration command in .rvi file", BAD_DATA_WARN);
            break;
          }
          else{
            kk = pSourceGrp->GetGlobalIndex();
          }
        }

        pModel->GetTransportModel()->AddDirichletCompartment(s[1], i_stor, kk, s_to_d(s[3]));
        //if time series is specified, s_to_d(time series file) returns zero
      }
      else{
        string warn=":FixedConcentration command: invalid state variable name"+to_string(s[2]);
        ExitGracefully(warn.c_str(),BAD_DATA_WARN);
      }
      break;
    }

    case (302)://----------------------------------------------
    {/*:MassInflux
       :MassInflux [string constit_name] [string water_state_var] [double flux (mg/m2/d) or time series] {optional HRU Group name}*/
      if (Options.noisy){cout <<"Mass Influx of constituent"<<endl;}

      if (!transprepared){
        ExitGracefully(":MassInflux command must be after :Transport command in .rvi file", BAD_DATA_WARN);
        break;
      }
      int layer_ind;
      int i_stor;
      sv_type typ=CStateVariable::StringToSVType(s[2],layer_ind,false);
      if (typ==UNRECOGNIZED_SVTYPE){
        WriteWarning(":MassInflux command: unrecognized storage variable name: "+to_string(s[2]),Options.noisy);
        break;
      }
      i_stor=pModel->GetStateVarIndex(typ,layer_ind);
      if (i_stor!=DOESNT_EXIST){
        int kk=DOESNT_EXIST;
        if (Len>4){
          CHRUGroup *pSourceGrp;
          pSourceGrp=pModel->GetHRUGroup(s[4]);
          if (pSourceGrp==NULL){
            ExitGracefully("Invalid HRU Group name supplied in :MassInflux command in .rvi file",BAD_DATA_WARN);
            break;
          }
          else{
            kk=pSourceGrp->GetGlobalIndex();
          }
        }
        pModel->GetTransportModel()->AddInfluxSource(s[1],i_stor,kk,s_to_d(s[3]));
        //if time series is specified, s_to_d(time series file) returns zero
      }
      else{
        ExitGracefully(":MassInflux command: invalid state variable name",BAD_DATA_WARN);
      }
      break;
    }
    case (303)://----------------------------------------------
    {/*:GeochemicalProcesses*/
      if (Options.noisy){cout <<"Start Geochemical Processes"<<endl;}
      break;//Does nothing for now; placeholder
    }
    case (304)://----------------------------------------------
    {/*:EndGeochemicalProcesses*/
      if (Options.noisy){cout <<"End Geochemical Processes"<<endl;}
      break;//Does nothing for now; placeholder
    }
    case (305)://----------------------------------------------
    {/*:Decay
       :Decay [ALGORITHM] [Constituent]
     */
      if (Options.noisy){cout <<"Decay Process"<<endl;}
      decay_type dec_type=DECAY_BASIC;
      if (Len<3){ImproperFormatWarning(":Decay",p,Options.noisy); break;}
      if      (!strcmp(s[1],"DECAY_BASIC"    )){dec_type=DECAY_BASIC;}
      else if (!strcmp(s[1],"DECAY_ANALYTIC" )){dec_type=DECAY_ANALYTIC;}
      else
      {
        ExitGracefully("ParseMainInputFile: Unrecognized decay process representation",BAD_DATA_WARN); break;
      }
      pMover=new CmvDecay(s[2],dec_type,pModel->GetTransportModel());
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case (306)://----------------------------------------------
    {/*:Advection
       :Advection RAVEN_DEFAULT [constituent]
     */
      if (Options.noisy){cout <<"Advection process"<<endl;}
      break;//Does nothing for now; placeholder- advection automatically handled when transport is specified. Should always be first process
    }
    case (307)://----------------------------------------------
    {/*:Transformation
       :Transformation [algorithm] [constituent1] [constituent2]
     */
      if (Options.noisy){cout <<"Transformation process"<<endl;}
      ExitGracefully(":Transformation input command",STUB);
      
      transformation_type t_type=TRANSFORM_LINEAR_ANALYTIC;
      if (Len<3){ImproperFormatWarning(":Decay",p,Options.noisy); break;}
      if      (!strcmp(s[1],"TRANSFORM_LINEAR"          )){t_type=TRANSFORM_LINEAR;}
      else if (!strcmp(s[1],"TRANSFORM_LINEAR_ANALYTIC" )){t_type=TRANSFORM_LINEAR_ANALYTIC;}
      else if (!strcmp(s[1],"TRANSFORM_MINERALIZATION"  )){t_type=TRANSFORM_MINERALIZATION; }
      else
      {
        ExitGracefully("ParseMainInputFile: Unrecognized transformation process representation",BAD_DATA_WARN); break;
      }
      pMover=new CmvTransformation(s[2],s[3],t_type,pModel->GetTransportModel());
      AddProcess(pModel,pMover,pProcGroup);
      break;
    }
    case(309):  //----------------------------------------------
    {/*Heat Conduction
     :HeatConduction RAVEN_DEFAULT MULTIPLE MULTIPLE */
      if(Options.noisy) { cout <<"Heat Conduction Process"<<endl; }
      if(Len<4) { ImproperFormatWarning(":HeatConduction",p,Options.noisy); break; }

      if(!strcmp(s[1],"RAVEN_DEFAULT")) {  }
      else {
        ExitGracefully("ParseMainInputFile: Unrecognized heat conduction process representation",BAD_DATA_WARN); break;
      }

      pMover=new CmvHeatConduction(pModel->GetTransportModel());
      AddProcess(pModel,pMover,pProcGroup);

      break;
    }
    case(400):  //----------------------------------------------
    {/*Reservoir demand allocation method
     :ReservoirDemandAllocation [string method]*/
      if(Options.noisy) { cout <<"Reservoir Demand Allocation Method"<<endl; }
      if(Len<2) { ImproperFormatWarning(":ReservoirDemandAllocation",p,Options.noisy); break; }
      if     (!strcmp(s[1],"DEMANDBY_CONTRIB_AREA"    )) { Options.res_demand_alloc =DEMANDBY_CONTRIB_AREA; }
      else if(!strcmp(s[1],"DEMANDBY_MAX_CAPACITY"    )) { Options.res_demand_alloc =DEMANDBY_MAX_CAPACITY; }
      //else if(!strcmp(s[1],"DEMANDBY_STOR_DEFICIT"  )) { Options.res_demand_alloc =DEMANDBY_STOR_DEFICIT; }
      else
      {
        ExitGracefully("ParseMainInputFile: Unrecognized Reservoir Demand Allocation Method",BAD_DATA_WARN); break;
      }
      break;
    }
    case(401):  //----------------------------------------------
    {/*Reservoir overflow handling method
     :ReservoirOverflowMode [string method]*/
      if(Options.noisy) { cout <<"Reservoir Overflow handling Method"<<endl; }
      if(Len<2) { ImproperFormatWarning(":ReservoirOverflowMode",p,Options.noisy); break; }
      if     (!strcmp(s[1],"OVERFLOW_NATURAL"       )) { Options.res_overflowmode =OVERFLOW_NATURAL; }
      else if(!strcmp(s[1],"OVERFLOW_STAGEDISCHARGE")) { Options.res_overflowmode =OVERFLOW_NATURAL; }
      else if(!strcmp(s[1],"OVERFLOW_ALL"           )) { Options.res_overflowmode =OVERFLOW_ALL;     }
      else
      {
        ExitGracefully("ParseMainInputFile: Unrecognized Reservoir Overflow handling Method",BAD_DATA_WARN); break;
      }
      break;
    }
    case(500): //----------------------------------------------
    {/*:ModelType" string type */
      if (Options.noisy) {cout <<"Model Type"<<endl;}
      if (Len<2){ImproperFormatWarning(":ModelType",p,Options.noisy); break;}
      if      (!strcmp(s[1],"SURFACE"     )){Options.modeltype=MODELTYPE_SURFACE;}
      else if (!strcmp(s[1],"GROUNDWATER" )){
          Options.modeltype=MODELTYPE_GROUNDWATER;
          pModel->AddAquiferStateVars(1);}
      else if (!strcmp(s[1],"COUPLED"     )){
          Options.modeltype=MODELTYPE_COUPLED;
          pModel->AddAquiferStateVars(1);}
      else {ExitGracefully("ParseInput:ModelType: Unrecognized method",BAD_DATA_WARN);}
      break;
    }
    case(501): //----------------------------------------------
    {/*:GWMethodOuter" string method */
      if (Options.noisy) {cout <<"Groundwater Nonlinear Solver Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":GWMethodOuter",p,Options.noisy); break;}
        if      (!strcmp(s[1],"NEWTON_RAPHSON")){Options.gw_solver_outer=GWSOL_NEWTONRAPHSON;}
      else if (!strcmp(s[1],"PICARD"        )){Options.gw_solver_outer=GWSOL_PICARD;}
      else {ExitGracefully("ParseInput:GWMethodOuter: Unrecognized method",BAD_DATA_WARN);}
      break;
    }
    case(502): //----------------------------------------------
    {/*:GWMethodInner" string method */
      if (Options.noisy) {cout <<"Groundwater Linear Solver Method"<<endl;}
      if (Len<2){ImproperFormatWarning(":GWMethodInner",p,Options.noisy); break;}
      if      (!strcmp(s[1],"PCGU"        )){Options.gw_solver_inner=GWSOL_PCGU;}
        else if (!strcmp(s[1],"BICGStab"    )){Options.gw_solver_inner=GWSOL_BICGSTAB;}
      else {ExitGracefully("ParseInput:GWMethodInner: Unrecognized method",BAD_DATA_WARN);}
      break;
    }
    case(503): //----------------------------------------------
    {/*:GWDiscretization" string gridtype */
      if (Options.noisy) {cout <<"Groundwater Discretization"<<endl;}
      if (Len<2){ImproperFormatWarning(":GWDiscretization",p,Options.noisy); break;}
      if      (!strcmp(s[1],"UNSTRUCT"    )){Options.gw_grid=GRID_UNSTRUCTURED;}
      else {ExitGracefully("ParseInput:GWDiscretization: Unrecognized method",BAD_DATA_WARN);}
      break;
    }
    case(504): //----------------------------------------------
    {/*:GWSWExchangeFreq" string timestep */
      if (Options.noisy) {cout <<"Groundwater-Surface Water Exchange Frequency"<<endl;}
      if (Len<2){ImproperFormatWarning(":GWSWExchangeFreq",p,Options.noisy); break;}
      if      (!strcmp(s[1],"EVERY_TIMESTEP")){Options.exchange_freq=EXFREQ_EVERY_TIMESTEP;}
      else {ExitGracefully("ParseInput:GWSWExchangeFreq: Unrecognized method",BAD_DATA_WARN);}
      break;
    }
    case(505): //----------------------------------------------
    {/*:FluxExchange" string type */
      if (Options.noisy) {cout <<"Groundwater-Surface Water Exchange Type"<<endl;}
      if (Len<2){ImproperFormatWarning(":FluxExchange",p,Options.noisy); break;}
      if      (!strcmp(s[1],"STANDARD"      )){Options.flux_exchange=FLUX_EX_STANDARD;}
      else if (!strcmp(s[1],"UFR"           )){Options.flux_exchange=FLUX_EX_UFR;}
      else {ExitGracefully("ParseInput:FluxExchange: Unrecognized method",BAD_DATA_WARN);}
      break;
    }
    case(506): //----------------------------------------------
    {/*:OverlapType" string type */
      if (Options.noisy) {cout <<"HRU-Aquifer Overlap Type"<<endl;}
      if (Len<2){ImproperFormatWarning(":OverlapType",p,Options.noisy); break;}
      if      (!strcmp(s[1],"DIRECT"          )){Options.overlap_type=DIRECT;}
      else if (!strcmp(s[1],"AREA_WEIGHTED"   )){Options.overlap_type=AREA_WEIGHTED;}
      else {ExitGracefully("ParseInput:OverlapType: Unrecognized method",BAD_DATA_WARN);}
      break;
    }
    case(510):  //--------------------------------------------
    {/*:WriteGroundwaterHeads */
      if (Options.noisy) {cout <<"Write Groundwater Head File ON"<<endl;}
      Options.write_gwhead =true;
      break;
    }
      case(511):  //--------------------------------------------
    {/*:WriteGroundwaterFlows */
      if (Options.noisy) {cout <<"Write Groundwater Flows File ON"<<endl;}
      Options.write_gwflow =true;
      break;
    }
    case(512):  //--------------------------------------------
    {/*:rvg_Filename */
      if (Options.noisy) {cout <<"rvg filename: "<<s[1]<<endl;}
      Options.rvg_filename=s[1];
      break;
    }
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
          string warn ="IGNORING unrecognized command: " + string(s[0])+ " in .rvi file";
          WriteWarning(warn,Options.noisy);
        }
      }
      else
      {
        string errString = "Unrecognized command in .rvi file:\n   " + string(s[0]);
        ExitGracefully(errString.c_str(),BAD_DATA_WARN);
      }
      break;
    }
    }//switch
  } //end while (!(p->Tokenize(s,Len)))
  INPUT.close();

  //Check input quality
  //===============================================================================================
  ExitGracefullyIf(Options.timestep<=0,
                   "ParseMainInputFile::Must have a postitive time step",BAD_DATA);
  ExitGracefullyIf((pModel->GetStateVarIndex(CONVOLUTION,0)!=DOESNT_EXIST) && (pModel->GetTransportModel()->GetNumConstituents()>0),
                   "ParseMainInputFile: cannot currently perform transport with convolution processes",BAD_DATA);
  if((Options.nNetCDFattribs>0) && (Options.output_format!=OUTPUT_NETCDF)){
    WriteAdvisory("ParseMainInputFile: NetCDF attributes were specified but output format is not NetCDF.",Options.noisy);
  }
  for(int i=0; i<pModel->GetNumStateVars();i++) {
    if((pModel->GetStateVarType(i)==SOIL) && ((pModel->GetStateVarLayer(i))>(Options.num_soillayers-1))) {
      string warn="A soil variable with an index ("+to_string(pModel->GetStateVarLayer(i))+") greater than that allowed by the limiting number of layers indicated in the :SoilModel command ("+to_string(Options.num_soillayers)+") was included in the .rvi file";
      ExitGracefully(warn.c_str(),BAD_DATA);
    }
  }


  //Add Ensemble configuration to Model
  //===============================================================================================
  CEnsemble *pEnsemble=NULL;
  if     (Options.ensemble==ENSEMBLE_NONE      ) {pEnsemble=new CEnsemble(1,Options);}
  else if(Options.ensemble==ENSEMBLE_MONTECARLO) {pEnsemble=new CMonteCarloEnsemble(num_ensemble_members,Options);}
  pModel->SetEnsembleMode(pEnsemble);

  delete p; p=NULL;
  delete [] tmpS;
  delete [] tmpLev;
  return true;
}

///////////////////////////////////////////////////////////////////
/// \brief This local method checks that the passes string includes either integer of state variable of variable string
///
/// \param *s [in] String state variable name
/// \param *&pModel [in] Input model object
/// \return Integer index of state variable ins tate variable arrays, or DOESNT_EXIST (-1) if is is invalid
int  ParseSVTypeIndex(string s,  CModel *&pModel)
{
  int ind;
  int layer_ind(-1);
  sv_type typ=CStateVariable::StringToSVType(s,layer_ind,false);
  ExitGracefullyIf(pModel==NULL,"ParseSVTypeIndex: NULL model!?",RUNTIME_ERR);

  if (typ==UNRECOGNIZED_SVTYPE)
  {
    //if IsNumeric(s){
    return DOESNT_EXIST;
  }

  ind = pModel->GetStateVarIndex(typ,layer_ind);

  if (ind!=DOESNT_EXIST){
    return ind;
  }
  else{
    //cout<<"offending state variable: "<<s<<endl;
    //ExitGracefully("ParseSVTypeIndex: storage variable used in process description does not exist",BAD_DATA);
    return DOESNT_EXIST;
  }
}
///////////////////////////////////////////////////////////////////
/// \brief writes warning to Raven errors file and screen for improper formatting
///
/// \param command [in] string command name (e.g., ":CanopyDrip")
/// \param *p [in]  parser object
/// \param noisy [in] boolean: if true, warning written to screen
//
void ImproperFormatWarning(string command, CParser *p, bool noisy)
{
  string warn;
  warn=command+" command: improper line length at line "+to_string(p->GetLineNumber());
  WriteWarning(warn,noisy);
}
///////////////////////////////////////////////////////////////////
/// \brief add process to either model or process group
//
void AddProcess(CModel *pModel,CHydroProcessABC* pMover,CProcessGroup *pProcGroup)
{
  if(pProcGroup==NULL){ pModel->AddProcess(pMover); }
  else                { pProcGroup->AddProcess(pMover); }
}
///////////////////////////////////////////////////////////////////
/// \brief adds to list of NetCDF attributes in Options structure
//
void AddNetCDFAttribute(optStruct &Options,const string att,const string &val)
{
  netcdfatt *aTmp=new netcdfatt[Options.nNetCDFattribs+1];
  for(int i=0;i<Options.nNetCDFattribs;i++){
    aTmp[i].attribute=Options.aNetCDFattribs[i].attribute;
    aTmp[i].value    =Options.aNetCDFattribs[i].value;
  }
  aTmp[Options.nNetCDFattribs].attribute=att;
  aTmp[Options.nNetCDFattribs].value=val;
  delete[]Options.aNetCDFattribs;
  Options.aNetCDFattribs=&aTmp[0];
  Options.nNetCDFattribs++;
}
