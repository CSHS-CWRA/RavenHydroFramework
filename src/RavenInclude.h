/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2024 the Raven Development Team

  Includes declaration of global constants, enumerated types, and
  shared common & hydrological functions
  - visible to all classes
  ----------------------------------------------------------------*/
#ifndef RAVENINCLUDE_H
#define RAVENINCLUDE_H

#define _CRT_SECURE_NO_DEPRECATE 1
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS 1
#endif

//#define _MODFLOW_USG_ // uncomment if compiling MODFLOW-USG coupled version of Raven
//#define _STRICTCHECK_ // uncomment if strict checking should be enabled (slows down model)
#ifndef _LPSOLVE_
//#define _LPSOLVE_       // uncomment if compiling lpsolve Demand Optimization version of Raven
#endif
#define STANDALONE
#ifdef netcdf
#define _RVNETCDF_      // if Makefile is used this will be automatically be uncommented if netCDF library is available
#endif
#ifdef _RVNETCDF_
#include <netcdf.h>
#endif

#include <stdlib.h>
#include <cstring>
#include <algorithm>
#include <math.h>
#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <strstream>
#include <sstream>

using namespace std;

/*****************************************************************
 NAMING AND PROGRAMMING CONVENTIONS IN RAVEN LIBRARY
------------------------------------------------------------------

all classes are named beginning with a capital C
all water/energy mover process classes are named beginning with Cmv
all constants/enumerated types are ALL CAPITALS

pointers to OBJECTS (not basic data types) or arrays of pointers to OBJECTS start with p
array sizes start with n followed by capital letters (e.g, nStateVars)
dynamic arrays start with a followed by capital letters (e.g., aStateVars)

counters- i    loops through state variables
          j    loops through Hydrological/Energy Processes
               loops through local constituent indices (CTransport)
          k    loops through HRUs
          kk   loops through land surface elements, HRU groups
          m    loops through soil layers, soil horizons and multilayer variables
          q    loops through connections between storage units/state vars
          c    loops through transported constituents
               loops through soil, lult, terrain, & vegetation classes
          n    loops through time series pulses
               loops through time arrays [in subbasin history]
          p,pp loops through subbasins, soil profiles
          g    loops through gauges
UNITS: all time units in days
       all vertical flux rates in mm/day
       all routing flow rates in m3/s
       all energy fluxes in MJ/m2/d
       all energy storage in terms of temp, C or MJ/m2
       all water storage units in mm
       all constituent storage units in mg/m2
       property units are those most commonly used
    -units outside these conventions should be used LOCALLY
    -if a constant is used in units conversion, the named variable
     (e.g., MM_PER_INCH) is ALWAYS to be used. That way, someone
     can look at the code and figure out what the mystery number means!
-------------------------------------------------------------------
Other comments:
  -processes that rely upon system memory (i.e., something is
   triggered 6 days after snow starts melting) are STRONGLY discouraged:
   the model should be able to stop and be restarted with final conditions
   as initial conditions to the next run, as is consistent with a
   mechanistic model
******************************************************************/
//*****************************************************************
// Global Variables (necessary, but minimized, evils)
//*****************************************************************
extern string g_output_directory; ///< Had to be here to avoid passing Options structure around willy-nilly
extern double g_debug_vars[10];   ///< can store any variables used during debugging; written to raven_debug.csv if debug_mode is on
extern bool   g_suppress_warnings;///< Had to be here to avoid passing Options structure around willy-nilly
extern bool   g_suppress_zeros;   ///< converts all output numbers less than REAL_SMALL to zero
extern bool   g_disable_freezing; ///< disables freezing impacts in thermal wrapper code
extern double g_min_storage;      ///< minimum soil storage
extern int    g_current_e;        ///< current ensemble member index

// Model version
const std::string __RAVEN_VERSION__   ="3.8.1";
//*****************************************************************
// Global Constants
//*****************************************************************
const double  ALMOST_INF              =1e99;                                    ///< Largest possible double value to be used
const double  PRETTY_SMALL            =1e-8;                                    ///< useful small tolerance for zero tests
const double  REAL_SMALL              =1e-12;                                   ///< Smallest possible double value to be used
const double  PI                      =3.1415926535898;                         ///< Double approximation of pi

const double  GRAVITY                 =9.80616;                                 ///< Standard acceleration due to gravity [m/s2]
const double  ZERO_CELSIUS            =273.16;                                  ///< Zero degrees in Kelvin

const double  UNIV_GAS_CONST          =8314.47;                                 ///< Universal gas constant [J/K/kmol]
const double  DRY_GAS_CONST           =287.042;                                 ///< Specific gas constant for dry air [J/K/kg]
const double  VAP_GAS_CONST           =461.504;                                 ///< Specific gas constant for water vapour [J/K/kmol]
const double  STEFAN_BOLTZ            =4.90e-9;                                 ///< Stephan-Boltzmann constant [MJ/m2/d/K4]
const double  VON_KARMAN              =0.42;                                    ///< Von Karmann constant

//units conversion constants
const double  MM_PER_INCH             =25.4;                                    ///< [in] to [mm]
const double  MM_PER_METER            =1000;                                    ///< [m] to [mm]
const double  MM_PER_CM               =10;                                      ///< [cm] to [mm]
const double  CM_PER_METER            =100;                                     ///< [m] to [cm]
const double  CM3_PER_METER3          =1e6;                                     ///< [m3] to [ccm]
const double  METER_PER_CM            =0.01;                                    ///< [cm] to [m]
const double  M2_PER_KM2              =1e6;                                     ///< [km2] to [m2]
const double  MM2_PER_M2              =1e6;                                     ///< [m2] to [mm2]
const double  M_PER_KM                =1000;                                    ///< [km] to [m]
const double  GRAMS_PER_KG            =1000;                                    ///< [kg] to [g]
const double  MG_PER_KG               =1000000;                                 ///< [kg] to [mg]
const double  LITER_PER_M3            =1000;                                    ///< [m3] to [l]
const double  GPCM3_PER_KGPM3         =0.001;                                   ///< [kg/m3] to [g/ccm]
const double  MJ_PER_J                =1e-6;                                    ///< [J] to [MJ]
const double  KPA_PER_MPA             =1000;                                    ///< [MPa] to [KPa]
const double  PA_PER_KPA              =1000;                                    ///< [kPa] to [Pa]
const double  HPA_PER_KPA             =10;                                      ///< [kPa] to [hPa]
const double  MB_PER_KPA              =10;                                      ///< [KPa] to [millibars]
const double  KPA_PER_ATM             =101.325;                                 ///< [atm] to [KPa]
const double  SEC_PER_DAY             =86400;                                   ///< days to seconds
const double  MIN_PER_DAY             =1440;                                    ///< days to minutes
const double  SEC_PER_HR              =3600;                                    ///< hours to seconds
const double  DAYS_PER_YEAR           =365.25;                                  ///< years to days
const double  HR_PER_DAY              =24;                                      ///< days to hours
const double  MJ_PER_D_TO_WATT        =11.574;                                  ///< [MJ/d] to [W]
const double  WATT_TO_MJ_PER_D        =0.0864;                                  ///< [W] to [MJ/d]
const double  MJ_PER_M2_LANGLEY       =0.04184;                                 ///< Langley to [MJ/m2]
const double  INCH_PER_METER          =39.37;                                   ///< [m] to [in]
const double  FEET_PER_METER          =3.28;                                    ///< [m] to [ft]
const double  ACREFTD_PER_CMS         =70.0456;                                 ///< [acre-ft/d] to [m3/s]
const double  CDM_PER_DAY_PER_CMS     =86.4;                                    ///< [cdm/day] to [m3/s]
const double  MPH_PER_KPH             =1.609;                                   ///< [kph] to [mph]
const double  MPH_PER_MPS             =2.237;                                   ///< [m/s] to [mph]
const double  RADIANS_TO_DEGREES      =57.29578951;                             ///< [rad] to [deg]
const double  DEGREES_TO_RADIANS      =0.017453293;                             ///< [deg] to [rad]

/// \details 0.278*(24hr/d)*(1000^2m^2/km)*(0.001m/mm)*(1/86400s/day)
/// runoff=RATIONAL_CONV*C_R*rainfall intensity
/// [mm/d]=                  [mm/d]
const double  RATIONAL_CONV           =0.7722;                                  ///< Allows rational method to be applied
const double  DAYS_PER_MONTH[12]      ={31,28,31,30,31,30,31,31,30,31,30,31};   ///< Array of doubles containing the number of days in each month

const double  FREEZING_TEMP           =0.0;                                     ///< [C] Freezing temperature of water

const double  DENSITY_AIR             =1.2466;                                  ///< [kg/m3] Ambient Air Density (@ 10 C)
const double  DENSITY_WATER           =1.000e3;                                 ///< [kg/m3] Water Density
const double  DENSITY_ICE             =0.917e3;                                 ///< [kg/m3] Ice Density
const double  DENSITY_SAND            =2.650e3;                                 ///< [kg/m3] Sand Density
const double  DENSITY_CLAY            =2.650e3;                                 ///< [kg/m3] Clay Density
const double  DENSITY_OM              =1.300e3;                                 ///< [kg/m3] Organic Matter Density
const double  MAX_SNOW_DENS           =0.350e3;                                 ///< [kg/m3] maximum dry density of snowpack (GAWSER)
const double  FRESH_SNOW_DENS         =0.119e3;                                 ///< [kg/m3] fresh snow density @ 0 deg. C
const double  TYPICAL_SNOW_DENS       =0.250e3;                                 ///< [kg/m3] "typical" snow density

const double  TC_WATER                =0.0492;                                  ///< [MJ/m/d/K] Thermal conductivity of Water (0.57 W/m/k)
const double  TC_ICE                  =0.194;                                   ///< [MJ/m/d/K] Thermal conductivity of Ice (2.24 W/m/K)
const double  TC_SAND                 =0.691;                                   ///< [MJ/m/d/K] Thermal conductivity of sand (8.0 W/m/K vs 8.8 in CLM)
const double  TC_CLAY                 =0.216;                                   ///< [MJ/m/d/K] Thermal conductivity of clay (2.5 W/m/K vs 2.92 in CLM)
const double  TC_ORGANIC              =0.0216;                                  ///< [MJ/m/d/K] Thermal conductivity of organic matter (0.25 W/m/K)
const double  TC_AIR                  =0.00199;                                 ///< [MJ/m/d/K] Thermal conductivity of air (0.023 W/m/K)

const double  COM_WATER		            =4.58e-10;                                ///< [1/Pa] Compressiblity of Water
const double  COM_ICE                 =4.58e-10;                                ///< [1/Pa] Compressiblity of Ice
const double  HCP_WATER               =4.186;                                   ///< [MJ/m3/K] Volumetric Heat Capacity of Water
const double  HCP_ICE                 =1.938;                                   ///< [MJ/m3/K] Volumetric Heat Capacity of Ice
const double  HCP_CLAY                =2.380;                                   ///< [MJ/m3/K] Volumetric Heat Capacity of Clay
const double  HCP_SAND                =2.130;                                   ///< [MJ/m3/K] Volumetric Heat Capacity of Sand
const double  HCP_ORGANIC             =2.500;                                   ///< [MJ/m3/K] Volumetric Heat Capacity of Organic Matter
const double  HCP_AIR                 =0.00124;                                 ///< [MJ/m3/K] Volumetric Heat capacity of air

const double  SPH_ICE                 =2.100e-3;                                ///< [MJ/kg/K] Specific heat capacity of ice
const double  SPH_WATER               =4.186e-3;                                ///< [MJ/kg/K] Specific heat capacity of water
const double  SPH_SAND                =0.835e-3;                                ///< [MJ/kg/K] Specific heat capacity of water
const double  SPH_VEGETATION          =2.700e-3;                                ///< [MJ/kg/K] Specific heat capacity of vegetation
const double  SPH_AIR                 =1.012e-3;                                ///< [MJ/kg/K] Specific heat capacity of air

const double  LH_FUSION               =0.334;                                   ///< [MJ/kg]  Latent heat of fusion
const double  LH_VAPOR                =2.501;                                   ///< [MJ/kg]  Latent heat of vaporization
const double  LH_SUBLIM               =2.845;                                   ///< [MJ/kg]  Latent heat of sublimation

const double  EMISS_ATM               =0.985;                                   ///< [-] emissivity of the atmosphere and snowpack
const double  EMISS_CANOPY            =0.96;                                    ///< [-] emissivity of the canopy
const double  EMISS_WATER             =0.97;                                    ///< [-] emissivity of open water

const double  AIR_H20_MW_RAT          =0.622;                                   ///< ratio of molecular weight of air to that of water
const double  MOLWT_WATER             =18.01528;                                ///< [g/mol] molecular weight of H2O
const double  AMBIENT_AIR_PRESSURE    =101.3;                                   ///< [kPa]

const double  EARTH_RADIUS            =6.3712e6;                                ///< [m]
const double  EARTH_ANG_VEL           =6.283185;                                ///< Earth's angular velocity [rad/d]
const double  SOLAR_NOON              =0.5;                                     ///< The fraction of a day at which solar noon occurs (12:00 PM)
const double  SOLAR_CONSTANT          =118.1;                                   ///< [MJ/m2/d]
const double  GLOBAL_ALBEDO           =0.3;                                     ///< Global Albedo used to calculate the backscattered radiation according to Dingman (Dingman 2008, E-20)
const double  PEAK_TEMP_HR            =3;                                       ///< if =3, 3:00 PM is time of max temperature, 3:00 AM is time of min temp
const double  WINTER_SOLSTICE_ANG     =6.111043;                                ///< Dec 21 as day angle

//Hard-coded Empirical parameters
const double  SAT_INF                 =0.92;                                    ///< [0..1] cutoff saturation for parabolic calculation of phi in clapp-hornberger soil characteristics
const double  MIN_PRESSURE            =1e10;                                    ///< [-mm] minimum matric potential in soil

const double  DEFAULT_SNOW_ALBEDO     =0.8;                                     ///< [0..1]  default snow albedo used for all snow, if not tracked explicitly
const double  NEGLIGBLE_SNOW          =0.1;                                     ///< [mm]    SWE at which snow cover no longer impacts albedo/evap calculations (i.e., snow cover~0)
const double  CZS                     =0.13;                                    ///< [0..1]  ratio of roughness to height for smooth closed canopies (from Brook90)
const double  CZS_HEIGHT              =1.0;                                     ///< [m]     height below which CZS applies
const double  CZR                     =0.05;                                    ///< [0..1]  ratio of roughness to height for rough closed canopies
const double  CZR_HEIGHT              =10.0;                                    ///< [m]     height above which czr applies
const double  CLOSED_LAI              =4.0;                                     ///< [m2/m2] minimum LAI defining a closed canopy (Shuttleworth and Wallace (1985))
const double  Z_REF_ADJUST            =2.0;                                     ///< [m]     reference height for weather data above canopy height
const double  WIND_EXTINCT            =2.5;                                     ///< [-]     wind/diffusivity extinction coefficient

/// \note This is always 2 for broadleaves, and ranges from 2 for flat needles to pi for cylindrical needles.
const double  LEAF_PROJ_RAT           =2.2;                                     ///< [m2/m2] ratio of total leaf area to projected area.

const double HBV_REFERENCE_ELEV       =1000;                                    ///< [masl]  (zref in HBV)
const double HBV_PET_ELEV_CORR        =0.0005;                                  ///< [mm/m-d] PET lapse rate (ECALT in HBV-EC)
const double HBV_PET_TEMP_CORR        =0.5;                                     ///< [-]     (ETF in HBV-EC)
const double GLOBAL_PET_CORR          =1.0;                                     ///< [-]     (ECORR in HBV-EC)

//Flag variables
const int     DOESNT_EXIST            =-1;                                      ///< return value for nonexistent index
const int     INDEX_NOT_FOUND         =-2;                                      ///< return value for index not found

const double  AUTO_COMPUTE            =-11111.1;                                ///< arbitrary value indicating that a parameter is to be autocalculated
const long    AUTO_COMPUTE_LONG       =-11111;                                  ///< arbitrary value indicating that a parameter is to be autocalculated
const double  NOT_SPECIFIED           =-33333.3;                                ///< arbitrary value indicating that a parameter has not been specified
const double  USE_TEMPLATE_VALUE      =-55555.5;                                ///< arbitrary value indicating that a parameter should be set to the template value
const double  NOT_NEEDED              =-66666.6;                                ///< arbitrary value indicating that a non-auto parameter is not needed for the current model configuration
const double  NOT_NEEDED_AUTO         =-77777.7;                                ///< arbitrary value indicating that a autogeneratable parameter is not needed for the current model configuration
const double  NETCDF_BLANK_VALUE      =-9999.0;                                 ///< NetCDF flag for blank value
const double  RAV_BLANK_DATA          =-1.2345;                                 ///< double corresponding to blank/void data item (also used in input files)
const double  DIRICHLET_TEMP          =-9999.0;                                 ///< dirichlet concentration flag corresponding to air temperature
const int     FROM_STATION_VAR        =-55;                                     ///< special flag indicating that NetCDF indices should be looked up from station attribute table

//Decision constants
const double  HUGE_RESIST             =1e20;                                    ///< [d/mm]   essentially infinite resistance
const double  SMALL_ROOT_DENS         =0.00001;                                 ///< [mm/mm3] essentially negligible root density
const double  SMALL_ROOT_LENGTH       =0.1;                                     ///< [m]      essentially negligible root length
const double  SMALL_FLOWRATE          =0.0001;                                  ///< [m3/s]   essentially negligible flow rate
const double  TIME_CORRECTION         =0.0001;                                  ///< [d]      offset for time series min/max functions
const double  DEFAULT_MAX_REACHLENGTH =10000.0;                                 ///< [km]     very large maximum reach length (defaults to single segment per reach)

//Special symbols
const char  DEG_SYMBOL                ='o';                                     ///< degree symbol, (or \0xB0)

const string ravenASCII[]={
 "    ::                                                                          .  ::         \n",
 "     -%:  .                                                                     +* .##.   .=  \n",
 "  -#.:%%:.*:                                                                   :%%.=%#.  :#=  \n",
 " .%%+:@@:-@=                                                                 .*%*:%%= .*@*.   \n",
 "-. :%@%=*@+-%#:  :                                                           .+**=##-:*%%-..=%\n",
 "#%+:.=%@%*%%*#@*.-%:                                                        :*****+=%%*::*%%*.\n",
 "-**%*--#%@%%%#*%*=@*.                                                   :*=**#%%%%%%%%@@%-::  \n",
 " :=#%%####%%%%%%%%%%*:..                                      ..::..=%@%%%%#*##%%%%%%%%%%+.   \n",
 "    ....-+*%%%@@%%%%*%@@@#*%=:                               :=-+%%%%%%@@%********###%%%%#:   \n",
 "    :%@@%%%%%%%%%+=+#@@%%##%@@%=+:.           .::.      .:=*%%*#**%%@%%@@@%%%%%%%@@@@@@%*     \n",
 "    :*%%%%%%%%%#######%%@%%@@@@%@%%@@*=-.   .+*##+::=*##+==%%%%%%##%%%##***#%%%%%%##*=-.      \n",
 "     .+@@@@@@%%%%%%#%%%%@@%%@@@@@@%%@%*+%%%%@O@%%O##@@@%*#%%%%%%%%#%@%##%###%%%@@@@%+:        \n",
 "        :@%%%%%%%%%%%%%%#%*%@@@%%@%@@@@@@@@@@@|*|%@@@@@@@*%%%@@@@##%%%%###%%%%%%%%%*-...      \n",
 "          :+******#%%%#***%@%#%%@%@@@@@@@@|@@@%v#%%%@@|@@%%##+#@@@@##%%#*#%%%%%%%%%:          \n",
 "              :#@@%%%%%%%%%%@@@%+%%@%%@%%#%@@@%%@%%@@@%*@%%%%#**#%%%###%%#**#%@%*.            \n",
 "             .=**%@@@%%%%@@@@%%#%%%%%@@-.*#%@@%@@@@@@=..::#%@%*###%%%%%%%%%%#-                \n",
 "                 :=+%@%#@@@@%%%%@%+%@@%= .-#%%@@@@@@%.    =@@%#%+#%%%%%#+=.                   \n",
 "                      .%%@@@:%@#*-*=   .+%%%%%@@@@%%@#:   :#%:= -:..                          \n",
 "                        .-.  ..  ..  .=%#%@%%@@@@@@%@%%*:                                     \n",
 "                        .          .-%**%%#+%%@@@@%#%%%%@=.                                   \n",
 "                                   :+=#%%%%+%%#%@%@#@@@%%%:..                                 \n",
 "                                    .#*=**.*#%%%%*%=--..-==:.                                 \n",
 "                                              :-*:                                            \n"};

//*****************************************************************
//Exit Strategies
//*****************************************************************

///////////////////////////////////////////////////////////////////
/// \brief Series of codes containing possible reasons for exiting
//
enum exitcode
{
  BAD_DATA,       ///< For bad input provided by user (requires immediate exit from program)
  BAD_DATA_WARN,  ///< For bad input provided by user (requires shutdown prior to simulation)
  RUNTIME_ERR,    ///< For runtime error (bad programming)
  FILE_OPEN_ERR,  ///< For bad file open (requires immediate exit)
  RAVEN_OPEN_ERR, ///< for bad RavenErrors.txt file open
  STUB,           ///< Stub function
  OUT_OF_MEMORY,  ///< When out of memory
  SIMULATION_DONE ///< Upon completion of the simulation
};

void FinalizeGracefully(const char *statement, exitcode code);  //defined in RavenMain.cpp
void ExitGracefully(const char *statement, exitcode code);      //defined in RavenMain.cpp

/////////////////////////////////////////////////////////////////
/// \brief In-line function that calls ExitGracefully function in the case of condition
///
/// \param condition [in] Boolean indicating if program should exit gracefully
/// \param statement [in] String to print to user upon exit
/// \param code [in] Code to determine why the system is exiting
//
inline void ExitGracefullyIf(bool condition, const char *statement, exitcode code)
{
  if (condition){ExitGracefully(statement,code);}
}

//*****************************************************************
//  Global Constants
//*****************************************************************
const bool    DESTRUCTOR_DEBUG    =false;       ///< if true, screen output is generated when destructor is called
const int     MAX_SV_LAYERS       =160;         ///< Max number of layers per state variable (greater than MAX_SOILLAYERS)
const int     MAX_SOILLAYERS      =50;          ///< Max number of soil layers in profile
const int     MAX_STATE_VAR_TYPES =100;         ///< Max number of *types* of state variables in model
const int     MAX_STATE_VARS      =200;         ///< Max number of simulated state variables manipulable by one process (CAdvection worst offender)
const int     MAX_CONNECTIONS     =200;         ///< Max number of to/from connections in any single process (CAdvection worst offender)
const int     MAX_LAT_CONNECTIONS =4000;        ///< Max number of lateral HRU flow connections
const int     MAX_SOIL_PROFILES   =200;         ///< Max number of soil profiles
const int     MAX_VEG_CLASSES     =200;         ///< Max number of vegetation classes
const int     MAX_LULT_CLASSES    =200;         ///< Max number of lult classes
const int     MAX_TERRAIN_CLASSES =50;          ///< Max number of terrain classes
const int     MAX_SURVEY_PTS      =50;          ///< Max number of survey points
const int     MAX_CONSTITUENTS    =10;          ///< Max number of transported constituents
const int     MAX_RIVER_SEGS      =50;          ///< Max number of river segments
const int     MAX_FILENAME_LENGTH =256;         ///< Max filename length
const int     MAX_MULTIDATA       =10;          ///< Max multidata length
/******************************************************************
Enumerated Types
   found in optStruct - the structure of global model options
******************************************************************/

///////////////////////////////////////////////////////////////////
/// \brief The numerical methods with which a solution to differential equations can be found
//
enum numerical_method
{
  EULER,              ///< Euler's method
  ORDERED_SERIES,     ///< Conventional WB model method - processes in series
  ITERATED_HEUN       ///< 2nd Order Convergence Method
};

///////////////////////////////////////////////////////////////////
/// \brief The type of model that is being simulated
//
enum model_type
{
  MODELTYPE_SURFACE,          ///< Surface water simulation ONLY
  MODELTYPE_COUPLED,          ///< Coupled groundwater-surface water simulation
};

///////////////////////////////////////////////////////////////////
/// \brief The frequency of the exchange fluxes between SW and GW
//
enum flux_frequency
{
  EXFREQ_EVERY_TIMESTEP,      ///< Feedback between SW and GW occurs once every timestep
  NONE,                       ///< No feedback between SW and GW
};

///////////////////////////////////////////////////////////////////
/// \brief The method of determining the amount of flux exchange crossing between SW and GW
//
enum flux_method
{
  FLUX_EX_STANDARD,           ///< Regular process calculations
  FLUX_EX_UFR,      ///< All fluxes are determined using rating curves
};

////////////////////////////////////////////////////////////////////
/// \brief The methods which can be used to route water along channel between basins
//
enum routing_method
{
  ROUTE_NONE,            ///< No routing to be simulated
  ROUTE_PLUG_FLOW,       ///< Plug flow routing - no dissipation of wave travelling at reach celerity
  ROUTE_MUSKINGUM,       ///< standard Muskingum algorithm
  ROUTE_MUSKINGUM_CUNGE, ///< Muskingum-Cunge algorithm
  ROUTE_STORAGECOEFF,    ///< Storage coefficient approach
  ROUTE_DIFFUSIVE_WAVE,  ///< diffusive wave approximation
  ROUTE_DIFFUSIVE_VARY,  ///< variable celerity/diffusivity diffusive wave approach
  ROUTE_HYDROLOGIC,      ///< simple iterative mass balance approach dS/dt=I-O
  ROUTE_TVD,             ///< Total variation diminishing approach of Schwanenberg and Montero, 2016
  ROUTE_EXTERNAL         ///< routing handled via external script/live file support
};

////////////////////////////////////////////////////////////////////
/// \brief Methods for routing water (lateral flow) to catchment outlet
//
enum catchment_route//methods used for routing water (lateral flow) to catchment outlet
{
  ROUTE_DUMP,                ///< dumps all surface water directly in channel
  ROUTE_DELAYED_FIRST_ORDER, ///< Linear Unit Hydrograph with lag time
  ROUTE_GAMMA_CONVOLUTION,   ///< Gamma Unit Hydrograph
  ROUTE_TRI_CONVOLUTION,     ///< Triangular Unit Hydrograph
  ROUTE_RESERVOIR_SERIES     ///< Series of linear reservoirs (Nash Hydrograph)
};

////////////////////////////////////////////////////////////////////
/// \brief Methods for interpolating Met Station/Gauge data to HRUs
//
enum interp_method
{
  INTERP_AVERAGE_ALL,                ///< Interpolation by taking the average of all values
  INTERP_NEAREST_NEIGHBOR,           ///< Interpolation by assuming the value of the nearest-neighbour
  INTERP_INVERSE_DISTANCE,           ///< Interpolates by a average of all values, weighted by the inverse distance to the interpolated point
  INTERP_INVERSE_DISTANCE_ELEVATION, ///< Interpolates by a average of all values, weighted by the inverse elevation distance to the interpolated point
  INTERP_FROM_FILE                   ///< User-specified file used to specify interpolation weights for all HRUs
};

////////////////////////////////////////////////////////////////////
/// \brief Methods used for calculating potential evapotranspiration (PET)
/// \docminor Describe the PET methods
//
enum evap_method
{
  PET_CONSTANT,                 ///< constant uniform PET
  PET_DATA,                     ///< read PET from time series file
  PET_NONE,                     ///< all PET calculations disabled - PET==0 (for Routing only mode)
  PET_LINEAR_TEMP,              ///< linearly proportional to temperature above zero
  PET_FROMMONTHLY,              ///< PET estimated from specified monthly averages and daily temperature
  PET_MONTHLY_FACTOR,           ///< PET estimated from specified monthly averages and daily temperature (UBCWM-style)
  PET_PENMAN_MONTEITH,          ///< Penman-Monteith equation
  PET_PENMAN_COMBINATION,       ///< Penman Combination approach
  PET_PRIESTLEY_TAYLOR,         ///< Priestley Taylor
  PET_HARGREAVES,               ///< Hargreaves method
  PET_HARGREAVES_1985,          ///< Hargreaves (1985) method
  PET_HAMON,                    ///< Hamon (19??) method
  PET_JENSEN_HAISE,             ///<
  PET_TURC_1961,                ///<
  PET_MAKKINK_1957,             ///<
  PET_SHUTTLEWORTH_WALLACE,     ///<
  PET_PENMAN_SIMPLE33,          ///< Simplified Penman equation from eqn 33 of Valiantzas (2006)
  PET_PENMAN_SIMPLE39,          ///< Simplified Penman equation from eqn 39 of Valiantzas (2006)
  PET_GRANGERGRAY,              ///< Granger  Gray PET from CRHM (Granger and Gray, 1989)
  PET_MOHYSE,                   ///< MOHYSE algorithm (https://docplayer.fr/69668879-Le-modele-hydrologique-mohyse.html)
  PET_OUDIN,                    ///< Simple PET from Oudin et. al., 2005
  PET_LINACRE,                  ///< From Linacre, Agricultural Meteorology, 1977
  PET_VAPDEFICIT,               ///< linear function of vapour deficit, c*(e_sat-e_a) from Seitz and Moore, 2020
  PET_BLENDED,                  ///< a blended combination of 2 or more of the methods above
  PET_UNKNOWN                   ///< special PET type for unrecognized commands
};

bool IsDailyPETmethod(evap_method method); //DEfined in Evaporation.cpp

////////////////////////////////////////////////////////////////////
/// \brief Methods used for correcting forcing functions with elevation
//
enum orographic_corr
{
  OROCORR_NONE,          ///< no orographic corrections
  OROCORR_SIMPLELAPSE,   ///< simple linear adiabatic lapse
  OROCORR_HBV,           ///< HBV-EC style orographic corrections
  OROCORR_UBCWM,         ///< UBCWM-style orographic corrections
  OROCORR_UBCWM2,        ///< UBCWM-style orographic corrections (simpler)
  OROCORR_PRMS           ///< PRMS-style orographic corrections
};

////////////////////////////////////////////////////////////////////
/// \brief Methods used for converting total precipitation to rain/snow
//
enum rainsnow_method
{
  RAINSNOW_DATA,         ///< Use data instead of calculating rain/snow partition
  RAINSNOW_DINGMAN,      ///< from Dingman - based upon min & max daily temperatures
  RAINSNOW_HBV,          ///< Linear variation between two temperatures - corrects only rain portion
  RAINSNOW_HSPF,         ///< HSPF approach - variable transition temperature
  RAINSNOW_UBCWM,        ///< Linear variation between two temperatures
  RAINSNOW_HARDER,       ///< Harder & Pomeroy (2013) method ported over from CRHM (Pomeroy et al 2007)
  RAINSNOW_WANG,         ///< Wang et al. (2019) sigmoid function
  RAINSNOW_SNTHERM89,    ///< Jordan et al (1991) fixed function used in SNTHERM.89 model and Noah-MP3.6
  RAINSNOW_THRESHOLD
};

////////////////////////////////////////////////////////////////////
/// \brief Methods used for estimating cloud cover
//
enum cloudcov_method
{
  CLOUDCOV_NONE,         ///< Cloud cover corrections not used
  CLOUDCOV_DATA,         ///< Uses gauge data for cloud cover
  CLOUDCOV_UBCWM         ///< Temperature range-based approach from UBCWM
};

////////////////////////////////////////////////////////////////////
/// \brief Methods correting shortwave radiation for cloud cover
//
enum SW_cloudcover_corr
{
  SW_CLOUD_CORR_NONE,     ///< Cloud cover corrections not used (e.g, when shortwave is measured by radiometer)
  SW_CLOUD_CORR_UBCWM,    ///< Based on UBCWM apporach, which is identical to Dingman (2008) Eq. 5-31
  SW_CLOUD_CORR_DINGMAN,  ///< Dingman (2008) Eq. 5-30     (does not require a parameter)
  SW_CLOUD_CORR_ANNANDALE ///< Annandale et al., (2002) temperature based estimate of cloud cover influence
};

////////////////////////////////////////////////////////////////////
/// \brief Methods used for estimating changes to snow albedo
//
enum snalbedo_method
{
  SNOW_ALBEDO_UBC           ///< snow albedo evolution according to UBCWM
};

////////////////////////////////////////////////////////////////////
/// \brief Methods used for estimating canopy transmittance of shortwave radiation
//
enum SW_canopy_corr
{
  SW_CANOPY_CORR_NONE,      ///< Applies no canopy correction for shortwave radiation
  SW_CANOPY_CORR_STATIC,    ///< Bulk shortwave canopy transmittance (no hourly or seasonal variation)
  SW_CANOPY_CORR_DYNAMIC,   ///< Dynamic (varies hourly and seasonal) shortwave canopy transmittance
  SW_CANOPY_CORR_UBCWM      ///< shortwave canopy transmittance according to UBCWM  (simple factor)
};

////////////////////////////////////////////////////////////////////
/// \brief Methods used for estimating net longwave radiation
//
enum LW_method
{
  LW_RAD_DATA,              ///< Longwave radiation specified in time series files
  LW_RAD_NONE,              ///< Longwave radiation not simulated (equal to zero)
  LW_RAD_DEFAULT,           ///< from Dingman text: uses Kustas (1994) approach for effective emissivity \cite Moran1994WRR
  LW_RAD_UBCWM,             ///< UBCWM approach
  LW_RAD_HSPF,              ///< HSPF approach (U.S. Corps of Engineers, 1956)
  LW_RAD_VALIANTZAS         ///< From Valiantzas, 2006 via Doorenbos and Pruit (1977) and Shuttleworth and Wallace (1952)
};
////////////////////////////////////////////////////////////////////
/// \brief Methods used for estimating incoming longwave radiation
//
enum LWinc_method
{
  LW_INC_DEFAULT,      ///< default implementation - reverts to Dingman 2014 if LW_RAD_DEFAULT is used (for backward compatibility)
  LW_INC_DATA,         ///< specified in time series files
  LW_INC_SICART,       ///< From Sicart et al. (2005) as ported over from CRHM
  LW_INC_SKYVIEW,      ///< simple sky view factor approach for stream temperature modelling, Using Prata 1996 clear sky emissivity
  LW_INC_DINGMAN       ///< from Dingman (2014), using Brutsaert 1975 clear sky emissivity
};

////////////////////////////////////////////////////////////////////
/// \brief Methods used for estimating clear sky and extraterrestrial shortwave radiation
//
enum SW_method
{
  SW_RAD_NONE,
  SW_RAD_DATA,              ///< Shortwave radiation specified in time series files
  SW_RAD_DEFAULT,           ///< from Dingman text
  SW_RAD_UBCWM,             ///< UBCWM approach
  SW_RAD_VALIANTZAS         ///< From Valiantzas, 2006
};
////////////////////////////////////////////////////////////////////
/// \brief Methods used for subdaily temporal downscaling of daily average PET and snowmelt
//
enum subdaily_method
{
  SUBDAILY_NONE,            ///< no correction for daily average values used
  SUBDAILY_SIMPLE,          ///< Use half-sine wave pulse from dawn to dusk
  SUBDAILY_UBC,             ///< from UBCWM - based upon cumulative temperature hours above zero Celsius
};
////////////////////////////////////////////////////////////////////
/// \brief Representations of soil characteristic curves
//
enum soil_charact
{
  BROOKS_COREY,             ///< Brooks-Corey model
  CLAPP_HORNBERGER,         ///< Clapp-Hornberger model
  MUALEM_VANGENUCHTEN       ///< Mualem-Van Genuchten model
};

////////////////////////////////////////////////////////////////////
/// \brief Types of HRUs
/// \remark Lake, Rock, and Glacier types are special because there is no soil storage
//
enum HRU_type
{
  HRU_STANDARD,             ///< Standard HRU
  HRU_LAKE,                 ///< Lake HRU
  HRU_WATER,                ///< Non-lake Waterbody HRU (streams)
  HRU_GLACIER,              ///< Glacier HRU
  HRU_WETLAND,              ///< Wetland HRU
  HRU_ROCK,                 ///< Open Rock or Pavement (impermeable) HRUs
  HRU_INVALID_TYPE          ///< returned if type is invalid
};

////////////////////////////////////////////////////////////////////
/// \brief Methods for estimating relative humidity
//
enum relhum_method
{
  RELHUM_CONSTANT,          ///< naive: constant relative humidity of 0.5
  RELHUM_MINDEWPT,          ///< uses minimum daily temperature as estimate of dew point
  RELHUM_DATA               ///< relative humidity specfied as time series at gauge
};

////////////////////////////////////////////////////////////////////
/// \brief Methods for estimating air pressure
//
enum airpress_method
{
  AIRPRESS_CONST,           ///< standard atm pressure at 20C
  AIRPRESS_DATA,            ///< air pressure specified as time series at gauge
  AIRPRESS_BASIC,           ///< power law correction for elevation (source unknown)
  AIRPRESS_UBC              ///< from UBC Watershed model - simple elevation correction
};

////////////////////////////////////////////////////////////////////
/// \brief Methods of calculating wind velocity
//
enum windvel_method
{
  WINDVEL_CONSTANT,         ///< naive: constant wind velocity of 3 m/s
  WINDVEL_DATA,             ///< wind velocity specfied as time series at gauge
  WINDVEL_UBCWM,            ///< from UBC Watershed model: daily temperature range-based
  WINDVEL_UBC_MOD,          ///< simplified version of UBCWM algorithm (linear relationship with T_max-T_min)
  WINDVEL_SQRT,             ///< linear relationshipw with sqrt(T_max-T_min)
  WINDVEL_LOG,              ///< linear relationship with ln(T_max-T_min)
};
////////////////////////////////////////////////////////////////////
/// \brief Methods of approximating vertical wind profile
//
enum wind_prof_meth
{
  WINDPROF_UNIFORM,          ///< naive: wind velocity doesnt change with height
  WINDPROF_LOGARITHMIC,      ///< basic von Karman logarithmic profile
  WINDPROF_MAHAT             ///< informed by presence of vegetation from Mahat (????)
};

////////////////////////////////////////////////////////////////////
/// \brief Methods of calculating precipitation interception fraction
//
enum precip_icept_method
{
  PRECIP_ICEPT_NONE,        ///< interception is not represented (for routing mostly)
  PRECIP_ICEPT_USER,        ///< pct of precip captured by canopy is user specified (TFRAIN,TFSNOW)
  PRECIP_ICEPT_LAI,         ///< pct of precip captured by canopy is linearly proportional to LAI (Dingman)
  PRECIP_ICEPT_EXPLAI,      ///< pct of precip captured by canopy is proportional to exp(LAI) (CLM)
  PRECIP_ICEPT_HEDSTROM,    ///< pct of snow captured by canopy is proportional to LAI & snowfall rate (Hedstrom & Pomeroy 1998)
  PRECIP_ICEPT_STICKY       ///< pct of snow captured by canopy is function of temperature (SUMMA-like)
};

////////////////////////////////////////////////////////////////////
/// \brief Methods of estimating potential melt
//
enum potmelt_method
{
  POTMELT_DEGREE_DAY,       ///< simple degree day method
  POTMELT_EB,               ///< energy balance approach
  POTMELT_RESTRICTED,       ///< restricted degree-day method
  POTMELT_DD_RAIN,          ///< degree day with rain-on-snow
  POTMELT_UBCWM,            ///< UBC watershed model approach
  POTMELT_HBV,              ///< custom degree day model used in HBV-EC
  POTMELT_HBV_ROS,          ///< HBV-EC degree day model with rain-on-snow
  POTMELT_DATA,             ///< user-specified potential melt forcing
  POTMELT_USACE,            ///< US Army Corps of Engineers Snow Melt
  POTMELT_CRHM_EBSM,        ///< Energy balance snow model from the Cold Regions Hydrology Model (CRHM) (Pomeroy, 2007)
  POTMELT_HMETS,            ///< From HMETS model (Martel et al., 2017)
  POTMELT_RILEY,            ///< From Riley et al., 1972, as reported in HYDROTEL 2.1 manual
  POTMELT_BLENDED,          ///< weighted average of multiple methods
  POTMELT_NONE,             ///< Potential melt not calculated
  POTMELT_UNKNOWN           ///< special case - can't recognize melt method in input
};

////////////////////////////////////////////////////////////////////
/// \brief Methods of estimating/generating recharge
//
enum recharge_method
{
  RECHARGE_NONE,            ///< assumes recharge=0
  RECHARGE_MODEL,           ///< Recharge from GW SV Storage
  RECHARGE_DATA             ///< recharge from (usually gridded) data (e.g., from other model)
};

////////////////////////////////////////////////////////////////////
/// \brief Methods of estimating/generating net shortwave radiation
//
enum netSWRad_method
{
   NETSWRAD_DATA,          ///< data supplied
   NETSWRAD_CALC           ///< determined via calculations
};

////////////////////////////////////////////////////////////////////
/// \brief Methods of handling snow cover depletion
//
enum snowcov_method
{
   SNOWCOV_NONE,           ///< no SDC
   SNOWCOV_LINEAR          ///< linear SDC characterized by threshold mean snow depth
};
////////////////////////////////////////////////////////////////////
/// \brief Methods of handling upstream allocation of downstream irrigation demand
//
enum demand_alloc
{
  DEMANDBY_MAX_CAPACITY,   ///< downstream demand support weighted by maximum storage capacity of reservoir
  DEMANDBY_CONTRIB_AREA,   ///< downstream demand support weighted by contributing area to reservoir
  DEMANDBY_STOR_DEFICIT    ///< downstream demand support weighted inversely to storage deficit
};

////////////////////////////////////////////////////////////////////
/// \brief Methods of performing monthly interpolations
//
enum monthly_interp
{
  MONTHINT_UNIFORM,       ///< specified value is constant for month
  MONTHINT_LINEAR_FOM,    ///< linear interpolation between specified data on first of the month
  MONTHINT_LINEAR_21,     ///< linear interpolation between specified data on twenty-first of the month
  MONTHINT_LINEAR_MID     ///< linear interpolation between specified data on median day of the month (14th,15th,or 16th)
};

////////////////////////////////////////////////////////////////////
/// \brief Basis on which a condition for applying a method/process is built
/// \docminor Review description of this enumerated type and its methods
//
enum condition_basis
{
  BASIS_HRU_TYPE,         ///< condition is based upon HRU type (e.g., if is a lake...)
  BASIS_HRU_GROUP,        ///< condition is based upon HRU group
  BASIS_LANDCLASS,        ///< condition is based upon land use/ land type class (e.g., if urban...)
  BASIS_VEGETATION        ///< condition is based upon vegetation class (e.g., if broadleaf...)
};

////////////////////////////////////////////////////////////////////
/// \brief Statistics of an aggregation
//
enum agg_stat
{
  AGG_AVERAGE,   ///< Average of data set
  AGG_MAXIMUM,   ///< Maximum of data set
  AGG_MINIMUM,   ///< Minimum of data set
  AGG_MEDIAN,    ///< Median of data set
  AGG_CUMULSUM,  ///< Cumulative sum of data set (of rates)
  AGG_RANGE,     ///< Range of data set
  AGG_95CI,      ///< 5% and 95% quantiles of data set
  AGG_QUARTILES, ///< Quartiles of data set
  AGG_HISTOGRAM  ///< Full histogram of data set
  //AGG_INT_AVG_DAY///< Average of julian day values  (not currently implemented)
};
////////////////////////////////////////////////////////////////////
/// \brief Ensemble mode type
//
enum ensemble_type
{
  ENSEMBLE_NONE,         ///< standard single model run
  ENSEMBLE_MONTECARLO,   ///< basic Monte Carlo simulation
  ENSEMBLE_DDS,          ///< DDS optimization run
  ENSEMBLE_ENKF          ///< Ensemble Kalman Filter data assimilation run
};

////////////////////////////////////////////////////////////////////
/// \brief Possible comparison results
//
enum comparison
{
  COMPARE_IS_EQUAL,      ///< Compared entities are equal
  COMPARE_NOT_EQUAL,     ///< Compared entities are not equal
  COMPARE_GREATERTHAN,   ///< Entity 1 greater than entity 2
  COMPARE_LESSTHAN,      ///< Entity 1 less than entity 2
  COMPARE_BETWEEN        ///< Entity 1 within range
};

////////////////////////////////////////////////////////////////////
/// \brief Abstration of a condition for applying a method/process
/// \remark stores information about a conditional such as "land use is equal to urban" or "HRU type is not LAKE"
//
struct condition
{
  condition_basis basis;          ///< basis for condition (e.g., land use)
  comparison      compare_method; ///< type of condition (equal or not equal)
  string          data;           ///< conditional data (thing to which basis is compared)
};
////////////////////////////////////////////////////////////////////
/// \brief perturbation adjustment - adding or multiplying
//
enum adjustment {
  ADJ_ADDITIVE,
  ADJ_MULTIPLICATIVE,
  ADJ_REPLACE
};
////////////////////////////////////////////////////////////////////
/// \brief Desired output format
//
enum out_format
{
  OUTPUT_STANDARD,       ///< Output in default Raven format (.csv files)
  OUTPUT_ENSIM,          ///< Output in Ensim format (.tb0 files)
  OUTPUT_NETCDF,         ///< Output in NetCDF format (.nc files)
  OUTPUT_NONE            ///< Output is suppressed
};

////////////////////////////////////////////////////////////////////
/// \brief reservoir constraint conditions
//
enum res_constraint
{
  RC_MAX_STAGE,
  RC_MIN_STAGE,
  RC_NATURAL,
  RC_TARGET,
  RC_DOWNSTREAM_FLOW,
  RC_MAX_FLOW_INCREASE,
  RC_MAX_FLOW_DECREASE,
  RC_MIN_FLOW,
  RC_MAX_FLOW,
  RC_OVERRIDE_FLOW,
  RC_DRY_RESERVOIR,
  RC_MANAGEMENT,    //from management optimization
  RC_DZTR
};

////////////////////////////////////////////////////////////////////
/// \brief reservoir overflow handling options - how discharge is estimated once max reservoir stage is met
//
enum overflowmode
{
  OVERFLOW_ALL,     ///< calculates Q required to fix stage at max value
  OVERFLOW_NATURAL  ///< uses stage discharge curve to calculate Q
};
////////////////////////////////////////////////////////////////////
/// \brief Types of state variable
/// \note If an additional state variable type is added, the following routines must be revised:
/// - CStateVariables::GetStateVarLongName (defined in StateVariables.cpp)
/// - CStateVariables::SVTypeToString  (defined in StateVariables.cpp)
/// - CStateVariables::StringToSVType  (defined in StateVariables.cpp)
/// - CStateVariables::IsWaterStorage  (defined in StateVariables.cpp)
/// - CStateVariables::IsEnergyStorage (defined in StateVariables.cpp)
/// - CHydroUnit::GetStateVarMax
/// - CModel::PartitionPrecip (if it is a water storage unit recieving precipitation)
///
//
enum sv_type
{
  //Water Storage
  SURFACE_WATER,           ///< [mm] Streams & rivers: see surface_struct (REQUIRED)
  ATMOSPHERE,              ///< [mm] atmosphere : recieves water only!! (REQUIRED)
  ATMOS_PRECIP,            ///< [mm] atmosphere : provides water only!! (REQUIRED)
  PONDED_WATER,            ///< [mm] water (melt & precip) waiting to infiltrate/runoff (REQUIRED)

  SOIL,                    ///< [mm] Shallow subsurface/vadose zone
  CANOPY,                  ///< [mm] Trees & vegetation
  CANOPY_SNOW,             ///< [mm] snow in canopy
  TRUNK,                   ///< [mm] water stored in trunks of trees
  ROOT,                    ///< [mm] water stored in roots
  GROUNDWATER,             ///< [mm] Deep groundwater
  DEPRESSION,              ///< [mm] depression/surface storage
  SNOW,                    ///< [mm] frozen snow depth (mm SWE : snow water equivalent)
  NEW_SNOW,                ///< [mm] new snowfall waiting to be handled by snow balance (as SWE)
  SNOW_LIQ,                ///< [mm] liquid water content of snowpack
  WETLAND,                 ///< [mm] deep wetland depression storage
  GLACIER,                 ///< [mm] Glacier melt/reservoir storage
  GLACIER_ICE,             ///< [mm] Glacier ice - typically assumed to be infinite reservoir.
  LAKE_STORAGE,            ///< [mm] Net lake storage - relative to equilibrium datum - can go negative

  CONVOLUTION,             ///< [mm] Convolution storage - for conceptual models with intermediate convolution steps
  CONV_STOR,               ///< [mm] Convolution sub-storage - tracks internal water mass for convolution

  // Distribution tracking variables
  MIN_DEP_DEFICIT,         ///< [mm or -1..0] Minimum depression deficit (describes deficit distribution), =-percent full if negative

  // Memory variables
  CUM_INFIL,               ///< [mm] Cumulative infiltration to topsoil
  GA_MOISTURE_INIT,        ///< [mm] Initial topsoil moisture content for Green Ampt infiltration
  CUM_SNOWMELT,            ///< [mm] Cumulative snowmelt as SWE
  AET,                     ///< [mm] PET used up in given time step (diagnostic variable)
  RUNOFF,                  ///< [mm] net release of water to surface water in given times step (diagnostic variable)

  //Temperature/Energy storage [C] or [MJ/m^2]
  ENERGY_LOSSES,           ///< [MJ/m2] general energy losses

  SURFACE_WATER_TEMP,      ///< [C] Temperature of surface water
  SNOW_TEMP,               ///< [C] Temperature of snow
  COLD_CONTENT,            ///< [C] Cold content of snowpack
  SOIL_TEMP,               ///< [C] Temperature of soil
  CANOPY_TEMP,             ///< [C] Temperature fo canopy

  //Snow/Glacier variables
  SNOW_DEPTH,              ///< [mm] Snow depth - surrogate for density
  PERMAFROST_DEPTH,        ///< [mm] depth of permafrost
  THAW_DEPTH,              ///< [mm] depth of thaw
  SNOW_DEPTH_STDDEV,       ///< log([mm]) Snow depth standard deviation
  SNOW_COVER,              ///< [0..1] fractional snow cover
  GLACIER_CC,              ///< [mm] cold content of glacier
  SNOW_DEFICIT,            ///< [mm] remaining holding capacity of snowpack (surrogate for SNOW_LIQ)
  SNOW_AGE,                ///< [d] snow age, in days
  SNODRIFT_TEMP,           ///< [C] temperature of drifting snow
  SNOW_DRIFT,              ///< [mm as SWE] drifting snow storage
  ICE_THICKNESS,           ///< [mm as SWE] lake ice thickness

  SNOW_ALBEDO,             ///< [-] Snow Surface albedo

  //Crop variables
  CROP_HEAT_UNITS,         ///< [-] cumulative crop heat units

  //Transport variables
  CONSTITUENT,             ///< chemical species [mg/m2], enthalpy [MJ/m2], or tracer [-]
  CONSTITUENT_SRC,         ///< chemical species [mg/m2], enthalpy [MJ/m2], or tracer [-] cumulative source
  CONSTITUENT_SW,          ///< chemical species [mg/m2], enthalpy [MJ/m2], or tracer [-] dumped to surface water
  CONSTITUENT_SINK,        ///< chemical species [mg/m2], enthalpy [MJ/m2], or tracer [-] cumulative sink (e.g., decay)

  //Energy variables
  LATENT_HEAT,             ///< [MJ/m2] cumulative energy lost to evaporative phase change

  //Lateral exchange
  LATERAL_EXCHANGE,        ///< [mm] water storage in transit from HRU awaiting lateral transfer to other HRUs

  //Special - internal flags
  STREAMFLOW,              ///< only used for referencing in data assimilation
  RESERVOIR_STAGE,         ///< only used for referencing in data assimilation
  UNRECOGNIZED_SVTYPE,     ///< Unrecognized type of state variable
  MULTIPLE_SVTYPE,         ///< Multiple State variables (for parsing)
  USERSPEC_SVTYPE          ///< User specified state variable type (for parsing)
};

////////////////////////////////////////////////////////////////////
/// \brief Types of hydrological processes
/// \note If an additional process type is added, the following routines must be revised:
/// - An additional HydroProcess class must be created (e.g., CmvInterflow)
/// - ParseInput.cpp: ParseMainInputFile: case(200+) (An additional parse statement is neccesary)
/// - GetProcessName (in CommonFunctions.cpp)
//
enum process_type
{
  //In RichardsEquation.h
  RICHARDS,

  //In Precipitation.h:
  PRECIPITATION,

  //In Infiltration.h
  INFILTRATION,

  //in SoilWaterMovers.h:
  BASEFLOW,SOIL_EVAPORATION,INTERFLOW,PERCOLATION,CAPILLARY_RISE,RECHARGE,SOIL_BALANCE,

  //in VegetationMovers.h:
  CANOPY_EVAPORATION, CANOPY_SNOW_EVAPORATION, CANOPY_DRIP,
  OPEN_WATER_EVAPORATION, LAKE_EVAPORATION,

  //in SnowMovers.h
  SNOWMELT,REFREEZE,SUBLIMATION,SNOW_BALANCE,SNOWSQUEEZE,SNOWTEMP_EVOLVE,

  //in GlacerProcesses.h
  GLACIER_MELT,GLACIER_RELEASE,GLACIER_INFIL,

  //in HydroProcessABC.h
  FLUSH, SPLIT, OVERFLOW_PROC,CONVOLVE,EXCHANGE_FLOW,

  //in LateralExchangeABC.h
  LAT_FLUSH, LAT_EQUIL,

  //in Albedo.h
  SNOW_ALBEDO_EVOLVE,

  //In PrairieSnow.h
  BLOWING_SNOW,

  //in CropGrowth.h
  CROP_HEAT_UNIT_EVOLVE,

  //in DepressionProcesses.h
  ABSTRACTION, DEPRESSION_OVERFLOW, SEEPAGE, LAKE_RELEASE,

  //in Advection.h, MassLoading.h
  ADVECTION, LAT_ADVECTION, MASS_LOADING,

  //in Decay.h
  DECAY, TRANSFORMATION,

  //in SurfaceEnergyExchange.h
  PARTITION_ENERGY,

  //in HeatConduction.h
  HEATCONDUCTION,

  //In FrozenGround.h
  GROUND_FREEZING,

  //In FrozenLake.h
  LAKE_FREEZING,

  //in ProcessGroup.h
  PROCESS_GROUP,

  //in GWSWProcesses.h
  DRAIN, GWRECHARGE,

  //in CustomHydProcess.h
  RF_PROCESS,

  //..
  NULL_PROCESS_TYPE
};


/******************************************************************
  Global Structures
******************************************************************/
////////////////////////////////////////////////////////////////////
/// \brief NetCDF attribute structure
//
struct netcdfatt
{
  string attribute;
  string value;
};

////////////////////////////////////////////////////////////////////
/// \brief Stores all global model and solution method options
//
struct optStruct
{
  string           version;                   ///< Raven version - written to output file headers

  double           julian_start_day;          ///< julian day corresponding to t=0, simulation time
                                              ///< e.g., 45.25 corresponds to 6AM on Feb. 14
  int              julian_start_year;         ///< year corresponding to t=0
  double           duration;                  ///< simulation duration
  int              calendar;                  ///< enum int that specifies calendar of all dates used,
  //                                          ///< e.g., 2 = "CALENDAR_PROLEPTIC_GREGORIAN"
  int              time_zone;                 ///< int that specifies time zone relative to GMT for writing to NetCDF output
  double           forecast_shift;            ///< amount for start date/assimilation date/end date to be adjusted by from command line

  numerical_method sol_method;                ///< numerical solution method
  double           convergence_crit;          ///< convergence criteria
  double           max_iterations;            ///< maximum number of iterations for iterative solver method
  double           timestep;                  ///< numerical method timestep (in days)
  double           output_interval;           ///< write to output file every x number of timesteps
  ensemble_type    ensemble;                  ///< ensemble type (or ENSEMBLE_NONE if single model)
  string           external_script;           ///< call to external script/.exe once per timestep (or "" if none)
  double           rvl_read_frequency;        ///< frequency to read rvl file (in d, or 0.0 if not to be read)
  bool             use_stopfile;              ///< true if Raven should look for stopfile

  model_type       modeltype;                 ///< type of model being simulated
  flux_frequency   exchange_freq;             ///< frequency of exchange fluxes
  flux_method      flux_exchange;             ///< method of calculating exchange fluxes

  interp_method    interpolation;             ///< Method for interpolating Met Station/Gauge data to HRUs
  string           interp_file;               ///< name of file (in working directory) which stores interpolation weights

  string           run_name;                  ///< prefix to be used for all output files
  char             run_mode;                  ///< run mode - single character used to enable multiple model configs with if statements (default==' ')
  string           rvi_filename;              ///< fully qualified filename of rvi (main input) file
  string           rvh_filename;              ///< fully qualified filename of rvh (HRU-basin) file
  string           rvp_filename;              ///< fully qualified filename of rvp (parameters) file
  string           rvt_filename;              ///< fully qualified filename of rvt (time series) file
  string           rvc_filename;              ///< fully qualified filename of rvc (initial conditions) file
  string           rve_filename;              ///< fully qualified filename of rve (ensemble) file
  string           rvl_filename;              ///< fully qualified filename of rvl (live communications) file
  string           rvg_filename;              ///< fully qualified filename of rvg (groundwater properties) file
  string           rvm_filename;              ///< fully qualified filename of rvm (management) file
  string           runinfo_filename;          ///< fully qualified filename of runinfo.nc file from FEWS
  string           stateinfo_filename;        ///< fully qualified filename of state_mods.nc file from FEWS
  string           flowinfo_filename;         ///< fully qualified filename of flowstate_mods.nc file from FEWS
  string           paraminfo_filename;        ///< fully qualified filename of param_mods.nc file from FEWS
  string           warm_ensemble_run;         ///< run name prefix in warm ensemble solution files

  string           main_output_dir;           ///< primary output directory (RavenErrors.txt, =output_dir for non-ensemble)
  string           output_dir;                ///< output directory (can change during ensemble run)

  orographic_corr  orocorr_temp;              ///< method for correcting interpolated temperatures for elevation
  orographic_corr  orocorr_precip;            ///< method for correcting interpolated precipitation for elevation
  orographic_corr  orocorr_PET;               ///< method for correcting interpolated PET for elevation

  // Algorithm Choices
  rainsnow_method    rainsnow;                ///< method for converting total precip to rain/snow
  cloudcov_method    cloud_cover;             ///< cloud cover estimation method
  snalbedo_method    snow_albedo;             ///< method for estimating snow albedo
  LW_method          LW_radiation;            ///< net longwave radiation estimation method
  LWinc_method       LW_incoming;             ///< incoming longwave radiation estimation method
  SW_method          SW_radiation;            ///< shortwave radiation estimation method
  SW_cloudcover_corr SW_cloudcovercorr;       ///< method for cloudcover correction of shortwave radiation
  SW_canopy_corr     SW_canopycorr;           ///< method for estimating canopy transmittance of shortwave radiation
  netSWRad_method    SW_radia_net;            ///< method for calculating net shortwave radiation (calculated or data)
  evap_method        evaporation;             ///< PET estimation method
  evap_method        ow_evaporation;          ///< Open Water PET estimation method
  evap_method        evap_infill;             ///< PET estimation method when infilling blank PET_DATA
  evap_method        ow_evap_infill;          ///< Open Water PET estimation method when infilling blank PET_DATA
  relhum_method      rel_humidity;            ///< Relative humidity estimation method
  airpress_method    air_pressure;            ///< Air pressure estimation method
  windvel_method     wind_velocity;           ///< Wind velocity estimation mehtod
  wind_prof_meth     wind_profile;            ///< Wind vertical profile approximation method
  potmelt_method     pot_melt;                ///< Potential melt estimation method
  subdaily_method    subdaily;                ///< Subdaily PET/Snowmelt temporal downscaling correction
  recharge_method    recharge;                ///< aquifer/soil recharge method
  bool               direct_evap;             ///< true if PET is used to directly reduce precipitation
  snowcov_method     snow_depletion;          ///< method for handling snowcover depletion curve

  precip_icept_method interception_factor;    ///< method for calculating canopy interception factor

  routing_method     routing;                 ///< channel routing method
  catchment_route    catchment_routing;       ///< catchment routing method
  demand_alloc       res_demand_alloc;        ///< method used for allocating upstream reservoir support to meet downstream irrigation demand
  overflowmode       res_overflowmode;        ///< method used for handling outflow estimates when max stage exceeded in reservoir
  monthly_interp     month_interp;            ///< means of interpolating monthly data

  bool               keepUBCWMbugs;           ///< true if peculiar UBCWM bugs are retained (only really for BC Hydro use)
  bool               suppressCompetitiveET;   ///< true if competitive ET should be suppressed (for backward compatibility)
  bool               snow_suppressPET;        ///< true if presence of snow should set PET to zero

  // Soil model information
  int                num_soillayers;          ///< number of soil layers
  soil_charact       soil_representation;     ///< characteristic curves for unsaturated flow

  // Output Options
  bool             debug_mode;                ///< true if debugging mode is on
  bool             noisy;                     ///< true if parsing information written to screen
  bool             silent;                    ///< true if nothing should be written to screen (overrides noisy)
  bool             pavics;                    ///< true if specific setings for PAVICS system are applied
  //                                          ///< (such as simulation status JSON file) (default: FALSE)
  bool             deltaresFEWS;              ///< true if input is generated by Deltares FEWS
  out_format       output_format;             ///< output format (default: OUTPUT_STANDARD)
  double           custom_interval;           ///< custom output interval (i.e., for generating 10-day interval outputs)
  bool             write_forcings;            ///< true if ForcingFunctions.csv is written
  bool             write_mass_bal;            ///< true if WatershedMassEnergyBalance.csv is written
  bool             write_gwhead;      		    ///< true if GWHead.csv is written
  bool             write_gwflow;      		    ///< true if GWFlow.csv is written
  bool             write_reservoir;           ///< true if ReservoirStages.csv is written
  bool             write_reservoirMB;         ///< true if ReservoirMassBalance.csv is written
  bool             write_waterlevels;         ///< true if WaterLevels.csv is written
  bool             ave_hydrograph;            ///< true if average flows over timestep are reported in hydrograph output
  bool             write_exhaustiveMB;        ///< true if exhaustive mass balance diagnostics are written
  int              write_group_mb;            ///< index (kk) of HRU Group for MB writing, DOESNT_EXIST if not to be written
  bool             write_channels;            ///< true if writing channel rating curve information
  bool             write_watershed_storage;   ///< true if WatershedStorage.csv/tb0/nc is written (default: true)
  bool             write_constitmass;         ///< true if constituent mass [mg/m2] is written instead of concentration [mg/L] in output files
  bool             write_basinfile;           ///< true if subbasins params are written to SubbasinParams.csv
  bool             write_interp_wts;          ///< true if interpolation weights are written to InterpolationWeights.csv
  bool             write_demandfile;          ///< true if demands.csv file is to be written
  bool             write_simpleout;           ///< true if simple_out.csv file is to be written (for scripting)
  bool             write_massloading;         ///< true if MassLoadings.csv file is to be written
  bool             write_localflow;           ///< true if local flows are written to Hydrographs file (csv or nc)
  bool             write_netresinflow;        ///< true if reservoir net inflows are written to Hydrographs file (csv or nc)
  bool             benchmarking;              ///< true if benchmarking output - removes version/timestamps in output
  bool             suppressICs;               ///< true if initial conditions are suppressed when writing output time series
  bool             period_ending;             ///< true if period ending convention should be used for reading/writing Ensim files
  bool             period_starting;           ///< true if all timestep-averaged output is reported using starttime of timestep
  bool             pause;                     ///< determines whether the simulation pauses at the end of the model run
  string           working_dir;               ///< working directory
  int              wateryr_mo;                ///< starting month of water year (typically 10=October)
  bool             create_rvp_template;       ///< create an rvp template file after reading the .rvi

  // Diagnostic options
  double           diag_start_time;           ///< Model time to start diagnostics
  double           diag_end_time;             ///< Model time to start diagnostics

  // Other
  bool             assimilate_flow;           ///< turn on streamflow assimilation
  bool             assimilate_stage;          ///< turn on lake stage assimilation
  double           assimilation_start;        ///< assimilation start time (in model time [d])
  bool             management_optimization;   ///< apply water management optimization (default: false)
  netcdfatt       *aNetCDFattribs;            ///< array of NetCDF attrributes {attribute/value pair}
  int              nNetCDFattribs;            ///< size of array of NetCDF attributes
  int              NetCDF_chunk_mem;          ///< [MB] size of memory chunk for each forcing grid
  bool             in_bmi_mode;               ///< true if in BMI mode (no rvt files, no end time)
};

///////////////////////////////////////////////////////////////////
/// \brief Stores physical location for interpolation/radiation calculations
//
struct location
{
  double latitude;  ///< Latitude
  double longitude; ///< Longitude
  double UTM_x;     ///< x-coordinate in Universal Transverse Mercator coordinate system
  double UTM_y;     ///< y-coordinate in Universal Transverse Mercator coordinate system
};

////////////////////////////////////////////////////////////////////
/// \brief Stores information describing a specific instance in time
//
struct time_struct
{
  string date_string;  ///< String date
  double model_time;   ///< [d] time elapsed since model start time
  double julian_day;   ///< [d] Julian-format decimal date (time in days since 0:00 Jan 1 of current year)
  int    day_of_month; ///< Day of month
  int    month;        ///< [1..12] month of year
  int    year;         ///< year
  bool   leap_yr;      ///< Boolean flag that indicates leap year
  bool   day_changed;  ///< Boolean flag indicating change of day for subdaily time steps
};

///////////////////////////////////////////////////////////////////
/// \brief Possible calendars
//
enum calendars
{
  CALENDAR_GREGORIAN, //same as STANDARD
  CALENDAR_PROLEPTIC_GREGORIAN,
  CALENDAR_365_DAY, //=NO_LEAP
  CALENDAR_360_DAY,
  CALENDAR_JULIAN,
  CALENDAR_366_DAY //=ALL_LEAP
};

///////////////////////////////////////////////////////////////////
/// \brief Stores external forcing functions for single HRU over current time step
/// \note if an additional forcing function is added, the following routines must be revised: \n
/// - CModel::UpdateHRUForcingFunctions \n
/// - CModel::GetAverageForcings() \n
/// - GetForcingTypeUnits,ZeroOutForcings,GetForcingFromType, ForcingToString in Forcings.cpp\n
//
const int MAX_FORCING_TYPES=50;
enum forcing_type
{
  F_PRECIP,         F_PRECIP_DAILY_AVE, F_PRECIP_5DAY,    F_SNOW_FRAC,
  F_RAINFALL,       F_SNOWFALL,         F_IRRIGATION,
  F_TEMP_AVE,
  F_TEMP_DAILY_MIN, F_TEMP_DAILY_MAX,   F_TEMP_DAILY_AVE,
  F_TEMP_MONTH_MAX, F_TEMP_MONTH_MIN,   F_TEMP_MONTH_AVE,
  F_TEMP_AVE_UNC,   F_TEMP_MIN_UNC,     F_TEMP_MAX_UNC,
  F_AIR_PRES,       F_AIR_DENS,         F_REL_HUMIDITY,
  F_CLOUD_COVER,    F_SW_RADIA,         F_LW_RADIA_NET,   F_ET_RADIA,        F_ET_RADIA_FLAT,
  F_SW_RADIA_NET,   F_SW_RADIA_UNC,     F_LW_INCOMING,    F_SW_RADIA_SUBCAN, F_SW_SUBCAN_NET,
  F_DAY_LENGTH,     F_DAY_ANGLE,        F_WIND_VEL,
  F_PET,F_OW_PET,   F_PET_MONTH_AVE,
  F_SUBDAILY_CORR,  F_POTENTIAL_MELT,
  F_RECHARGE,
  F_PRECIP_CONC,
  F_PRECIP_TEMP,
  F_UNRECOGNIZED
};
////////////////////////////////////////////////////////////////////
/// \brief Contains information on forcing functions
//
struct force_struct
{
  double precip;              ///< precipitaiton rate over time step [mm/d]
  double precip_daily_ave;    ///< average precipitaiton over day (0:00-24:00) [mm/d]
  double precip_5day;         ///< 5-day precipitation total [mm] (needed for SCS)
  double snow_frac;           ///< fraction of precip that is snow [0..1]
  double irrigation;          ///< irrigation rate over time step [mm/d]
  double precip_temp;         ///< precipitation temperature [C]
  double precip_conc;         ///< precipitation concentration [C] (\todo[funct]: should make vector)

  double temp_ave;            ///< average air temp over time step [C]
  double temp_daily_min;      ///< minimum air temperature over day (0:00-24:00)[C]
  double temp_daily_max;      ///< maximum air temperature over day (0:00-24:00)[C]
  double temp_daily_ave;      ///< average air temp over day (0:00-24:00) [C]
  double temp_month_max;      ///< maximum air temp during month [C]
  double temp_month_min;      ///< minimum air temp during month [C]
  double temp_month_ave;      ///< average air temp during month [C]
  double temp_ave_unc;        ///< uncorrected daily average air temp  [C]
  double temp_min_unc;        ///< uncorrected daily min air temp  [C]
  double temp_max_unc;        ///< uncorrected daily max air temp  [C]

  double air_dens;            ///< Air density [kg/m3]
  double air_pres;            ///< Air pressure [kPa]
  double rel_humidity;        ///< relative humidity [0..1]

  double cloud_cover;         ///< Cloud cover [0..1]
  double ET_radia;            ///< uncorrected extraterrestrial shortwave radiation, after slope correction [MJ/m2/d]
  double ET_radia_flat;       ///< uncorrected extraterrestrial shortwave radiation on horizontal surface [MJ/m2/d]
  double SW_radia_unc;        ///< uncorrected shortwave radiation (before cloud and canopy corrections, after slope)[MJ/m2/d]
  double SW_radia;            ///< Incoming shortwave radiation (slope/air mass/horizon/cloud cover/canopy corrections applied, uncorrected for albedo) [MJ/m2/d]
  double SW_radia_net;        ///< Net shortwave radiation at canopy top (albedo corrected) [MJ/m2/d]
  double SW_radia_subcan;     ///< Shortwave incoming radiation at ground beneath canopy [MJ/m2/d]
  double SW_subcan_net;       ///< Net shortwave incoming radiation on ground beneath canopy [MJ/m2/d]
  double LW_radia_net;        ///< Net longwave radiation [MJ/m2/d]
  double LW_incoming;         ///< Incoming longwave radiation [MJ/m2/d]
  double day_length;          ///< day length [d]  (e.g., ~0.5 for equinox @ 45lat, 0.0 for the north pole during winter)
  double day_angle;           ///< day angle [0..2PI] =0 for Jan 1, 2pi for Dec 31

  double potential_melt;      ///< potential snowmelt rate [mm/d]

  double wind_vel;            ///< Wind velocity [m/s]

  double PET;                 ///< Potential Evapotranspiration [mm/d]
  double OW_PET;              ///< Open Water Potential Evapotranspiration [mm/d]
  double PET_month_ave;       ///< average PET during month [mm/d]

  double recharge;            ///< recharge to groundwater (typically from external model) [mm/d]

  double subdaily_corr;       ///< a subdaily correction factor to re-distribute daily average PET or snowmelt [-]
};
////////////////////////////////////////////////////////////////////
/// \brief probability distributions
//
enum disttype
{
  DIST_UNIFORM,     ///< Uniform distribution
  DIST_NORMAL,      ///< Gaussian (normal) distribution
  DIST_LOGNORMAL,   ///< Lognormal distribution
  DIST_GAMMA        ///< Gamma distribution
};
////////////////////////////////////////////////////////////////////
/// \brief Contains information on forcing function random perturbations
//
struct force_perturb
{
  forcing_type forcing;
  disttype     distribution;
  double       distpar[3];   ///< distribution parameters - conditional on distribution
  int          kk;           ///< HRU group index (or DOESNT_EXIST if perturbation should apply everywhere)
  adjustment   adj_type;     ///< additive or multiplicative adjustment
  double      *eps;          ///< stored adjustment factors for day - allows perturbations to be pre-calculated for day so that daily mean/min/max can be generated
};

/******************************************************************
  Other Functions (defined in CommonFunctions.cpp)
******************************************************************/
string GetProcessName(process_type ptype);

//Array Functions--------------------------------------------------
bool   DynArrayAppend(void **& pArr,void *xptr,int &size);
int    SmartIntervalSearch(const double &x, const double *ax, const int N,const int ilast);
int    NearSearchIndex(const int i, int guess_p, const int size);

//Threshold Correction Functions-----------------------------------
double threshPositive(const double &val);
double threshMin     (const double &val1, const double &val2, const double &smooth_coeff);
double threshMax     (const double &val1, const double &val2, const double &smooth_coeff);

//Time/Date Functions----------------------------------------------
bool        IsLeapYear(            const int         year,
                                   const int         calendar);
void        JulianConvert(               double      model_time,
                                   const double      start_date,
                                   const int         start_year,
                                   const int         calendar,
                                                     time_struct &tt);
string      DecDaysToHours(        const double      dec_date,
                                   const bool        truncate=false);
double      InterpolateMo(         const double      aVal[12],
                                   const time_struct &tt,
                                   const optStruct   &Options);
time_struct DateStringToTimeStruct(const string      sDate,
                                         string      sTime,
                                   const int         calendar);

string      TimeZoneToString      (const int tz);
double      TimeDifference        (const double      jul_day1,
                                   const int         year1,
                                   const double      jul_day2,
                                   const int         year2,
                                   const int         calendar);
void        AddTime               (const double      jul_day1,
                                   const int         year1,
                                   const double      &daysadded,
                                   const int         calendar,
                                         double      &jul_day_out,
                                         int         &year_out);
int         StringToCalendar      (      string      cal_chars);
string      GetCurrentMachineTime ();
double      FixTimestep           (      double      tstep);
bool        IsValidDateString     (const string      sDate);
double      RoundToNearestMinute  (const double& t);
bool        IsInDateRange         (const double &julian_day,
                                   const int    &julian_start,
                                   const int    &julian_end);
int      GetJulianDayFromMonthYear(const string &date_str, const int calendar);

bool        IsValidNetCDFTimeString   (const string unit_t_str);
time_struct TimeStructFromNetCDFString(const string unit_t_str,
                                       const string timestr,
                                       const int calendar,
                                       double &time_zone);
//Conversion Functions-------------------------------------------
double CelsiusToFarenheit       (const double &T);

//*****************************************************************
//Common Functions (inline)
//*****************************************************************

///////////////////////////////////////////////////////////////////
/// \brief Returns the maximum of two doubles
/// \param r [in] First double for comparison
/// \param l [in] Second double for comparison
/// \return Maximum of the two doubles passed as parameters
//
inline double        max    (double r, double l)  {return ((r>=l) ? r : l);}

///////////////////////////////////////////////////////////////////
/// \brief Returns the minimum of two doubles
/// \param r [in] First double for comparison
/// \param l [in] Second double for comparison
/// \return Minimum of the two doubles passed as parameters
//
inline double        min    (double r, double l)  {return ((r<=l) ? r : l);}

///////////////////////////////////////////////////////////////////
/// \brief Returns the maximum of two integers
/// \param r [in] First integer for comparison
/// \param l [in] Second integer for comparison
/// \return Maximum of the two integers passed as parameters
//
inline int           max    (int    r, int    l)  {return ((r>=l) ? r : l);}

///////////////////////////////////////////////////////////////////
/// \brief Returns the minimum of two integers
/// \param r [in] First integer for comparison
/// \param l [in] Second integer for comparison
/// \return Minimum of the two integers passed as parameters
//
inline int           min    (int    r, int    l)  {return ((r<=l) ? r : l);}

///////////////////////////////////////////////////////////////////
/// \brief Returns a boolean indicating if the passed doubles have opposite signs
/// \param r [in] First double for comparison
/// \param l [in] Second double for comparison
/// \return Boolean indicating if the passed doubles have opposite signs
//
inline bool          oppsign(double r, double l)  {return ((r*l) < 0.0);   }

///////////////////////////////////////////////////////////////////
/// \brief Raises int parameter l to the power of int param r and returns the result
/// \param r [in] Integer exponent
/// \param l [in] Integer base
/// \return Boolean indicating if the passed doubles have opposite signs
//
inline int           ipow   (int    l, int    r)  {return (int)(pow((double)(l),r)); }

///////////////////////////////////////////////////////////////////
/// \brief Returns the absolute value of an integer
/// \param i1 [in] Integer value whose absolute value will be determined
/// \return Absolute value of the integer parameter i1
//
inline int           iabs   (int i1)              {return (int)(fabs((double)(i1))); }

///////////////////////////////////////////////////////////////////
/// \brief Returns the square root of an integer
/// \param i1 [in] Integer value whose root will be determined
/// \return Square root (double) of the integer parameter i1
//
inline double        sqrt   (int i1)              {return sqrt((double)(i1));}

///////////////////////////////////////////////////////////////////
/// \brief Returns the exponential value of the integer parameter i1
/// \param i1 [in] Integer value whose exponential value will be determined
/// \return Exponential value (double) of the integer parameter i1
//
inline double        exp    (int i1)              {return exp((double)(i1));}

///////////////////////////////////////////////////////////////////
/// \brief Converts string parameter to integer type
/// \param *s1 [in] String to be converted to integer
/// \return Integer format of passed string
//
inline int           s_to_i (const char *s1)            {return (int)atof(s1);   }

///////////////////////////////////////////////////////////////////
/// \brief Converts string parameter to long long integer type
/// \param *s1 [in] String to be converted to long long integer
/// \return Integer format of passed string
//
inline long long int s_to_ll (const char *s1)            {return atoll(s1);   }

///////////////////////////////////////////////////////////////////
/// \brief Converts string parameter to double type
/// \param *s1 [in] String to be converted to double
/// \return Double format of passed string
//
inline double        s_to_d (const char *s1)            {return atof(s1);        }

///////////////////////////////////////////////////////////////////
/// \brief returns true if character string is long integer
/// \return true if character string is long integer
//
inline bool StringIsLong(const char *s1)
{
  char *p;
  strtoll(s1,&p,10);
  return !(*p);
}
///////////////////////////////////////////////////////////////////
/// \brief Converts string parameter (which represents an integer) to boolean type
/// \note This function takes string parameters which represent integers
/// \remark Assumes any non-zero value to be TRUE
/// \param *s1 [in] String to be converted to boolean
/// \return Boolean format of passed string
//
inline bool          s_to_b (const char *s1)            {return ((int)atof(s1)!=0);   }


///////////////////////////////////////////////////////////////////
/// \brief Converts string parameter (which represents an integer) to long int type
/// \note This function takes string parameters which represent long
/// \return long int format of passed string
//
inline long          s_to_l (const char *s1)            {return (long)atof(s1);   }

///////////////////////////////////////////////////////////////////
/// \brief Extracts long integer start and endpoints from string range of integers (i.e., string of the form "4-10")
/// \param *s1 [in] String range to be parsed
/// \param &v1 [out] Start point of range
/// \param &v2 [out] End point of range
//
inline void      s_to_range (const char *s1, long long int &v1, long long int &v2)
{
  string s=s1;
  int p=(int)(s.find_last_of("-"));
  if (p==-1){v1=v2=s_to_ll(s1);return;}
  else{
    v1=atoll(s.substr(0,(unsigned int)(p)).c_str());
    v2=atoll(s.substr((unsigned int)(p)+1,1000).c_str());
  };
}
/* //_finite causes all sorts of issues on multi-platform compilation
#ifndef _WIN32
#define _finite(v) finite(v)
#endif
*/

///////////////////////////////////////////////////////////////////
/// \brief converts any input to string
/// \param t [in] thing to be converted to a string
/// \return string equivalent of input t
//
template <class T>
inline std::string to_string (const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}

///////////////////////////////////////////////////////////////////
/// \brief Converts boolean to string
/// \param b [in] boolean to be converted to string
/// \return "true" or "false"
//
inline string          b_to_s (const bool b)            {return  (b ? "true" : "false");   }

///////////////////////////////////////////////////////////////////
/// \brief returns the modulo of parameter n with respect to d
/// \param &n [in] dividend
/// \param &d [in] divisor
/// \return modulo of parameter n with respect to d
//
inline double        ffmod  (const double &n,const double &d){
  double a=(double)((int)(n/d+0.5));//rounds n/d down to nearest int
  return (n/d-a)*d;
}

///////////////////////////////////////////////////////////////////
/// \brief Converts string parameter (which represents an integer) to boolean type
/// \note This function takes string parameters of the form "True" or "False" (not integers)
/// \remark Assumes any "non-true" string to be false
/// \param *s1 [in] String to be converted to boolean
/// \return Boolean format of passed string
//
inline bool          s_to_bool (const char *s1)   {
  if ((!strcmp(s1,"True")) || (!strcmp(s1,"TRUE")) || (!strcmp(s1,"true"))){return true;}
  return false;
}
///////////////////////////////////////////////////////////////////
/// \brief returns true if string is numeric
/// \param *s1 [in] String to be checked
/// \return true if string is numeric
/// \note possible (rare) exception: 00BAQ, or other string starting with 0
//
double fast_s_to_d             (const char *s);
inline bool          is_numeric(const char *s1){
  if ((fast_s_to_d(s1) == 0.0) && (s1[0] != '0')){return false;} //likely numeric
  return true;
}


///////////////////////////////////////////////////////////////////
/// \brief Assigns to double parameter u the value of v if v is lower than &u
/// \param &u [in & out] Double value to (potentially) be changed
/// \param v [in] Double value for comparison
//
inline void        lowerswap(double &u,const double v){if (v<u){u=v;}}

///////////////////////////////////////////////////////////////////
/// \brief Assigns to double parameter u the value of v if v is greater than &u
/// \param &u [in & out] Double value to (potentially) be changed
/// \param v [in] Double value for comparison
//
inline void        upperswap(double &u,const double v){if (v>u){u=v;}}

///////////////////////////////////////////////////////////////////
/// \brief Assigns to integer parameter &u the value of v if v is lower than &u
/// \param &u [in & out] Integer value to (potentially) be changed
/// \param v [in] Integer value for comparison
//
inline void        lowerswap(int    &u,const int    v){if (v<u){u=v;}}

///////////////////////////////////////////////////////////////////
/// \brief Assigns to integer parameter &u the value of v if v is greater than &u
/// \param &u [in & out] Integer value to (potentially) be changed
/// \param v [in] Integer value for comparison
//
inline void        upperswap(int    &u,const int    v){if (v>u){u=v;}}
///////////////////////////////////////////////////////////////////
/// \brief rounds v to 3 digits after decimal point
/// \param &u [in & out] Integer value to (potentially) be changed
/// \param v [in] Integer value for comparison
//
inline double roundTo3dig(const double &v) {
  return (double)((int)(v*1000.0+0.5)/1000.0);
}

//Parsing Functions-------------------------------------------
//defined in CommonFunctions.cpp
double   AutoOrDouble           (const string s);
string   StringToUppercase      (const string &s);
bool     IsComment              (const char *s, const int Len);
void     WriteWarning           (const string warn, bool noisy);
void     WriteAdvisory          (const string warn, bool noisy);
HRU_type StringToHRUType        (const string s);
double   fast_s_to_d            (const char *s);
double   FormatDouble           (const double &d);
void     SubstringReplace       (string& str,const string& from,const string& to);


//defined in NetCDFReading.cpp
int  GetCalendarFromNetCDF      (const int ncid,int varid_t,const string filename,const optStruct &Options);
void GetTimeVectorFromNetCDF    (const int ncid,const int varid_t,const int ntime,double *my_time);
void GetTimeInfoFromNetCDF      (const char *unit_t,int calendar,const double *time,const int ntime,const string filename,
                                 double &tstep,double &start_day,int &start_yr,double &time_zone);
void GetJulianDateFromNetCDFTime(const string unit_t,const int calendar,const double &time,
                                 double &start_day,int &start_yr);

//I/O Functions-----------------------------------------------
//defined in StandardOutput.cpp
void   PrepareOutputdirectory    (const optStruct &Options);
string GetDirectoryName          (const string &fname);
void   HandleNetCDFErrors        (int error_code);        ///< NetCDF error handling
string CorrectForRelativePath    (const string filename, const string relfile);
string GetFileExtension          (string filename);
string FilenamePrepare           (string filebase, const optStruct &Options);

#ifdef _WIN32
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

//--Autocompute Functions-------------------------------------------
//defined in CommonFunctions.cpp
bool   SetCalculableValue      (double &val, double set_val, double template_val);
bool   SetSpecifiedValue       (double &val, double set_val, double template_val, bool needed, string name);
double DefaultParameterValue   (bool is_template, bool is_computable);

//Mathematical Functions-------------------------------------------
//defined in CommonFunctions.cpp
double rvn_erfc         (const double &x);
double rvn_erf          (const double &x);
double rvn_trunc        (const double &x);
double rvn_floor        (const double &x);
double rvn_round        (const double &x);
bool   rvn_isnan        (const double &x);
double log_pdf          (const double &x, const double &mu, const double &sig);
double LambertN         (const double &x, const int N);
double GammaCumDist     (const double &t, const double &a, const double &b);
double IncompleteGamma  (const double &betat, const double &a);
double TriCumDist       (const double &t, const double &tc, const double &tp);
double NashCumDist      (const double &t, const double &k, const int &NR);
double ADRCumDist       (const double &t, const double &L, const double &v, const double &D);
double TimeVaryingADRCumDist(const double &t,const double &L,const double *v,int nv,const double &D,const double &dt);
void   CalcWeightsFromUniformNums(const double* aVals, double* aWeights, const int N);

//Array processing Functions-------------------------------------------------
//defined in CommonFunctions.cpp
void   quickSort        (double arr[], int left, int right) ;
double InterpolateCurve (const double x,const double *xx,const double *y,int N,bool extrapbottom);
void   getRanks         (const double *arr, const int N, int *ranks);
void   pushIntoIntArray (int*&a, const int &v, int &n);

//Geographic Conversion Functions-----------------------------------
//defined in UTM_to_LatLong.cpp
void LatLonToUTMXY (const double lat, //latitude, in decimal degrees
                    const double lon, //longitude, in decimal degrees
                    const int    zone,//UTM zone
                    double &x,
                    double &y);

//Hydrological Functions-------------------------------------------
//defined in CommonFunctions.cpp
double GetSaturatedVaporPressure(const double &T);
double GetSatVapSlope           (const double &T, const double &satvap);
double GetLatentHeatVaporization(const double &T);
double GetPsychrometricConstant (const double &P, const double &LH_vapor);
double GetWetBulbTemperature    (const double &P,const double &T,const double &rel_hum);
double GetAirDensity            (const double &T, const double &P);
double GetSpecificHumidity      (const double &rel_hum,const double &air_press,const double &T);
double GetVerticalTransportEfficiency     (const double &P,   //[kPa]
                                           const double &meas_ht,        //[m]
                                           const double &zero_pl,       //[m]
                                           const double &rough);        //[m]
double CalcAtmosphericConductance(const double &wind_vel,     //[m/d]
                                  const double &meas_ht,      //[m]
                                  const double &zero_pl,      //[m]
                                  const double &rough_ht,     //[m]
                                  const double &vap_rough_ht); //[m]
double GetDewPointTemp          (const double &Ta,const double &rel_hum);
double GetDewPointTemp          (const double &E);
double ConvertVolumetricEnthalpyToTemperature(const double &hv);
double ConvertTemperatureToVolumetricEnthalpy(const double &T, const double &pctfroz);
double ConvertVolumetricEnthalpyToIceContent (const double &hv);
double TemperatureEnthalpyDerivative         (const double &hv);

//Snow Functions---------------------------------------------------
//defined in SnowParams.cpp and PotentialMelt.cpp
class CModelABC;  // defined in ModelABC.h
double CalcFreshSnowDensity       (const double &air_temp);
double GetSnowThermCond           (const double &snow_dens);
double GetSensibleHeatSnow        (const double &air_temp,const double &surf_temp,const double &V, const double &ref_ht, const double &rough);
double GetLatentHeatSnow          (const double &P,const double &air_temp,const double &surf_temp,const double &rel_humid,const double &V,const double &ref_ht,const double &rough);
double GetRainHeatInput           (const double &surf_temp, const double &air_temp,const double &rain_rate,const double &rel_humid);
double GetSnowDensity             (const double &snowSWE,const double &snow_depth);
double CalculateSnowLiquidCapacity(const double &SWE, const double &snow_depth, const CModelABC* pModel);

//Crop Functions---------------------------------------------------
bool   IsGrowingSeason            (const time_struct &tt, const double &CHU);



#endif
