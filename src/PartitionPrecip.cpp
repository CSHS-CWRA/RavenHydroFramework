/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2020 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Model.h"
#include "Precipitation.h"

/*****************************************************************
   Precipitation Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default precipitation constructor
//
CmvPrecipitation::CmvPrecipitation():
									CHydroProcessABC(PRECIPITATION)
{}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvPrecipitation::~CmvPrecipitation(){}

//////////////////////////////////////////////////////////////////
/// \brief Initialize the connections within the precipitation model
//
void CmvPrecipitation::Initialize()
{
  
  int iAtmos=pModel->GetStateVarIndex(ATMOS_PRECIP);
  
  int N=12;

  //if (pModel->GetOptStruct()->simulate_lakesnow){N+=3;}

  CHydroProcessABC::DynamicSpecifyConnections(N);
  
  iFrom[0]=iAtmos; iTo[0]=pModel->GetStateVarIndex(PONDED_WATER);
  iFrom[1]=iAtmos; iTo[1]=pModel->GetStateVarIndex(CANOPY);
  iFrom[2]=iAtmos; iTo[2]=pModel->GetStateVarIndex(CANOPY_SNOW);
  iFrom[3]=iAtmos; iTo[3]=pModel->GetStateVarIndex(DEPRESSION);
  iFrom[4]=iAtmos; iTo[4]=pModel->GetStateVarIndex(ATMOSPHERE);
  iFrom[5]=iAtmos; iTo[5]=pModel->GetStateVarIndex(SNOW);
  iFrom[6]=iAtmos; iTo[6]=pModel->GetStateVarIndex(SNOW_LIQ);
  iFrom[7]=iAtmos; iTo[7]=pModel->GetStateVarIndex(SNOW_DEFICIT);
  iFrom[8]=iAtmos; iTo[8]=pModel->GetStateVarIndex(NEW_SNOW);
  iFrom[9]=iAtmos; iTo[9]=pModel->GetLakeStorageIndex(); //Lake handling
  
  iFrom[10]=iTo[10]=pModel->GetStateVarIndex(COLD_CONTENT);
  iFrom[11]=iTo[11]=pModel->GetStateVarIndex(SNOW_DEPTH);

  for(int q=0;q<12;q++) {
    if (iTo  [q]==DOESNT_EXIST){iTo  [q]=iAtmos; } //add dummy slot for missing connections
    if (iFrom[q]==DOESNT_EXIST){iFrom[q]=iAtmos; } //since this routine reacts to presence or absence of individual state variables
  }

  //lake snow handling
  /*if (pModel->GetOptStruct()->simulate_lakesnow)
  {
    iFrom[12]=pModel->GetStateVarIndex(SNOW);         iTo[12]=pModel->GetLakeStorageIndex();
    iFrom[13]=pModel->GetStateVarIndex(SNOW_LIQ);     iTo[13]=pModel->GetLakeStorageIndex();
    iFrom[14]=pModel->GetStateVarIndex(PONDED_WATER); iTo[14]=pModel->GetLakeStorageIndex(); 
  }*/
}
//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for baseflow algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by baseflow algorithm (size of aP[] and aPC[])
//
void CmvPrecipitation::GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const
{
  int iCan     = pModel->GetStateVarIndex(CANOPY     );
  int iSnowCan = pModel->GetStateVarIndex(CANOPY_SNOW);
  bool canopy_exists  = (iCan     != DOESNT_EXIST);
  bool cansnow_exists = (iSnowCan != DOESNT_EXIST);

  nP=0;
  aP[nP]="SAI_HT_RATIO";      aPC[nP]=CLASS_VEGETATION; nP++; //reasonable defaults
  aP[nP]="FOREST_SPARSENESS"; aPC[nP]=CLASS_LANDUSE;    nP++; //reasonable defaults
  aP[nP]="MAX_LAI";           aPC[nP]=CLASS_VEGETATION; nP++;

  //Interception factor parameters identified in CModel::GetParticipatingParamList()
  if (canopy_exists)
  {
    aP[nP] = "MAX_CAPACITY";      aPC[nP] = CLASS_VEGETATION; nP++;
    aP[nP] = "RELATIVE_LAI";      aPC[nP] = CLASS_VEGETATION; nP++;
  }
  if (cansnow_exists)
  {
    aP[nP] = "MAX_SNOW_CAPACITY"; aPC[nP] = CLASS_VEGETATION; nP++;
    aP[nP] = "RELATIVE_LAI";      aPC[nP] = CLASS_VEGETATION; nP++;
  }

}
//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
///
/// \param *aSV [out] Array of state variable types needed by algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by algorithm (size of aSV[] and aLev[] arrays)
//
void CmvPrecipitation::GetParticipatingStateVarList(sv_type *aSV, int *aLev, int &nSV)
{
  nSV=1;
  aSV[0]=ATMOS_PRECIP; aLev[0]=DOESNT_EXIST;
  //Fully dynamic-Precip expects any and all variables.
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rate of loss of water from precipitation to surface storage [mm/d]
/// \details Partitions Precip over the timestep to (mostly surficial) storage.
///     Partitioning depends upon which storage units exist.
/// - precip is always treated constant over a given time step,
/// - though partioning may vary based upon current system storage \n
///     returns precipitation (in mm/day) moving to each storage unit
///     over time step (stored in rates[]).
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time structure
/// \param *rates [out] rates of precip to various storage compartments [mm/d]
//
void CmvPrecipitation::GetRatesOfChange(const double             *state_vars,
                                        const CHydroUnit *pHRU,
                                        const optStruct  &Options,
                                        const time_struct &tt,
                                        double     *rates) const
{

  int    iCan,iSCan,iSnLiq,iLake;
  double total_precip,rainfall,snowfall;      //all [mm/day]
  double rain_int,snow_int,rainthru,snowthru; //all [mm/day]
  double capacity,snow_capacity;              //[mm]
  double rain_pct,snow_pct;                   //[0..1]
  double Fcan;                                //[0..1] fraction canopy covered
  double Fs;                                  //[0..1] fraction snow covered
  double Fsnow;                               //[0..1] percent precip as snow
  double air_temp;                            //[C]
  double tstep = Options.timestep;

  iCan    =pModel->GetStateVarIndex(CANOPY);
  iSCan   =pModel->GetStateVarIndex(CANOPY_SNOW);
  iSnLiq  =pModel->GetStateVarIndex(SNOW_LIQ);
  iLake   =pModel->GetLakeStorageIndex();

  //connection indices (hard-coded in Initialize())
  int qPond=0;
  int qCan=1;
  int qCanSnow=2;
  int qDep=3;
  int qAtmos=4;
  int qSnow=5;
  int qSnowLiq=6;
  int qSnowDef=7;
  int qNewSnow=8;
  int qLake=9;
  int qColdContent=10;
  int qSnowDepth=11;
    
  //get precipitation information from gauges
  total_precip=pHRU->GetForcingFunctions()->precip;//[mm/day]
  Fsnow       =pHRU->GetForcingFunctions()->snow_frac;
  air_temp    =pHRU->GetForcingFunctions()->temp_ave;

  //identify snow/rain fraction
  if (!pModel->StateVarExists(SNOW)){Fsnow=0;}
  snowfall=(    Fsnow)*total_precip;
  rainfall=(1.0-Fsnow)*total_precip;//[mm/day]

  double SWE=0.0;
  if(pModel->StateVarExists(SNOW)) {
    SWE=state_vars[pModel->GetStateVarIndex(SNOW)];
  }
  if (!pHRU->IsLake()) 
  {
    // Calculate Snow Cover
    //-----------------------------------------------------------------
    Fs = pHRU->GetSnowCover();

    // Throughfall & interception for rain & snow
    //-----------------------------------------------------------------
    //calculate rain interception-----------------------------------
    snow_int = rain_int = 0.0;
    Fcan = pHRU->GetSurfaceProps()->forest_coverage;

    //total rain interception rate over canopy covered portion of HRU
    rain_pct = pHRU->GetVegVarProps()->rain_icept_pct;
    rain_int = rainfall*rain_pct;
      
    if (pModel->StateVarExists(CANOPY)) 
    {
      capacity = pHRU->GetVegVarProps()->capacity;
      rain_int = threshMin(rain_int, max((capacity - state_vars[iCan]) / tstep, 0.0), 0.0);
      rates[qCan] = Fcan*rain_int;
    }
    else{//if no canopy SV is present, rain_pct moves to ATMOSPHERE (via "evaporation")
      rates[qAtmos] += Fcan*rain_int;
    }
    rainthru = rainfall - Fcan*rain_int;//[mm/day]
    
    //rate of snow interception over canopy portion of HRU
    snow_pct = pHRU->GetVegVarProps()->snow_icept_pct;
    snow_int = snowfall*snow_pct;

    //calculate snow interception-----------------------------------
    if (pModel->StateVarExists(CANOPY_SNOW))
    {
      snow_capacity = pHRU->GetVegVarProps()->snow_capacity;
      snow_int = max(threshMin(snow_int, max((snow_capacity - state_vars[iSCan]) / tstep, 0.0), 0.0), 0.0);
      rates[qCanSnow] = Fcan*snow_int;
    }
    else//if no canopy SV is present, snow_pct moves to ATMOSPHERE (via "sublimation")
    {
      rates[qAtmos] += Fcan*snow_int;
    }
    snowthru = snowfall - Fcan*snow_int;
    
    //Handle all snow variables
    //--------------------------------------------------------------------
    if(!pModel->StateVarExists(NEW_SNOW))
    {
      if (pModel->StateVarExists(SNOW) )
      {
        rates[qSnow]=snowthru;

        if (pModel->StateVarExists(SNOW_DEPTH)){
          rates[qSnowDepth] = snowthru*DENSITY_ICE / CalcFreshSnowDensity(air_temp);
        }
        if (pModel->StateVarExists(COLD_CONTENT)){
          rates[qColdContent] = snowthru*threshPositive(FREEZING_TEMP - air_temp)*HCP_ICE/ MM_PER_METER;//[mm/d]*[K]*[MJ/m3/K]*[m/mm]=MJ/m2/d
        }

        //move rain through to snowmelt--------------------------------
        if ((pModel->StateVarExists(SNOW_LIQ)) && (SWE > 0.0) && (!pModel->StateVarExists(SNOW_DEFICIT))) //Snow liquid in model
        {
          double melt_cap = pHRU->GetStateVarMax(iSnLiq, state_vars, Options);
          double fill_rate;
          fill_rate = Fs*threshMin(max(melt_cap - state_vars[iSnLiq], 0.0) / Options.timestep, rainthru, 0.0);
          rates[qSnowLiq] = fill_rate;
          rainthru -= fill_rate;
        }
        if ( (pModel->StateVarExists(SNOW_DEFICIT)) && (SWE > 0.0))  //Snow deficit in model (UBCWM)
        {
          double SWI = CGlobalParams::GetParams()->snow_SWI;
          rates[qSnowDef] = SWI*snowthru; //snowfall to snowpack
        }
      }
    }
    else{
      rates[qNewSnow]=snowthru;
    }

    //everything remaining goes to ponded water variable, usually handled by infiltration & abstraction algorithms
    if(pHRU->GetHRUType()!=HRU_WETLAND) 
    {
      rates[qPond]=rainthru;
    } //ponded water goes to depression in wetland 
    else 
    {
      rates[qDep] =rainthru;
      /*if(!wetland_frozen){ // \todo [funct] handle snow on frozen wetlands
        rates[qSnow]-=snowthru;
        rates[qDep] +=snowthru;
      }*/
    }
  }//End if isn't lake
  else /* in a lake HRU */ 
  {
    if (false){//(pModel->GetOptStruct()->simulate_lakesnow) {
    /*
      bool lake_frozen=false;
      if (lake_frozen)
      { // handle snow on frozen lakes
        if(!pModel->StateVarExists(NEW_SNOW))
        {
          if (pModel->StateVarExists(SNOW) )
          {
            rates[qSnow]=snowthru;

            if (pModel->StateVarExists(SNOW_DEPTH)){
              rates[qSnowDepth] = snowthru*DENSITY_ICE / CalcFreshSnowDensity(air_temp);
            }
            if (pModel->StateVarExists(COLD_CONTENT)){
              rates[qColdContent] = snowthru*threshPositive(FREEZING_TEMP - air_temp)*HCP_ICE/ MM_PER_METER;//[mm/d]*[K]*[MJ/m3/K]*[m/mm]=MJ/m2/d
            }

            //move rain through to snowmelt--------------------------------
            if ((pModel->StateVarExists(SNOW_LIQ)) && (SWE > 0.0) && (!pModel->StateVarExists(SNOW_DEFICIT))) //Snow liquid in model
            {
              double melt_cap = pHRU->GetStateVarMax(iSnLiq, state_vars, Options);
              double fill_rate;
              fill_rate = Fs*threshMin(max(melt_cap - state_vars[iSnLiq], 0.0) / Options.timestep, rainthru, 0.0);
              rates[qSnowLiq] = fill_rate;
              rainthru -= fill_rate;
            }
            if ( (pModel->StateVarExists(SNOW_DEFICIT)) && (SWE > 0.0))  //Snow deficit in model (UBCWM)
            {
              double SWI = CGlobalParams::GetParams()->snow_SWI;
              rates[qSnowDef] = SWI*snowthru; //snowfall to snowpack
            }
          }
        }
        else{
          rates[qNewSnow]=snowthru;
        }
      }
      else{ //recently unfrozen
        rates[qSnowToLake ]+=state_vars[iSnow]/Options.timestep;
        rates[qSnLiqToLake]+=state_vars[iSnLiq]/Options.timestep;

        //should also support multilayer snowpack, but that's an additional complication...
      }
    */
    //rates[qLake]+=rainthru;
    //rates[qPondToLake ]+=state_vars[iPond]/Options.timetep; //handles on-lake snowmelt 
    }
    else{
      rates[qLake]=snowfall+rainfall; // in lake, all water just added
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Applies constraints to baseflow
/// \details For all methods, ensures that rate of flow cannot drain "from" compartment over timestep.
///
/// \param *state_vars [in] array of state variable values for this HRU (of size CModel::_nStateVars)
/// \param *pHRU [in] Pointer to HRU object
/// \param &Options [in] Global model options information
/// \param &tt [in] Current input time structure
/// \param *rates [out] Rates of change of state variables (of size MAX_CONNECTIONS)
//
void   CmvPrecipitation::ApplyConstraints(const double     *state_vars,
                                          const CHydroUnit *pHRU,
                                          const optStruct  &Options,
                                          const time_struct &tt,
                                          double     *rates) const
{
  //handled in RatesOfChange routine
}
