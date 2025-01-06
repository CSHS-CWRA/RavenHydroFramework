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
CmvPrecipitation::CmvPrecipitation(CModelABC *pModel)
  :CHydroProcessABC(PRECIPITATION, pModel)
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

  int N=13;

  if (pModel->StateVarExists(ICE_THICKNESS)){
    N+=3;
    if (pModel->StateVarExists(DEPRESSION)){N+=3;}
  }
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
  iFrom[10]=iAtmos; iTo[10]=pModel->GetStateVarIndex(SURFACE_WATER); //reservoir handling

  iFrom[11]=iTo[11]=pModel->GetStateVarIndex(COLD_CONTENT);
  iFrom[12]=iTo[12]=pModel->GetStateVarIndex(SNOW_DEPTH);

  for(int q=0;q<13;q++) {
    if (iTo  [q]==DOESNT_EXIST){iTo  [q]=iAtmos; } //add dummy slot for missing connections
    if (iFrom[q]==DOESNT_EXIST){iFrom[q]=iAtmos; } //since this routine reacts to presence or absence of individual state variables
  }

  //frozen lake & wetland snow handling
  if (pModel->StateVarExists(ICE_THICKNESS))
  {
    iFrom[13]=pModel->GetStateVarIndex(SNOW);         iTo[13]=pModel->GetLakeStorageIndex();
    iFrom[14]=pModel->GetStateVarIndex(SNOW_LIQ);     iTo[14]=pModel->GetLakeStorageIndex();
    iFrom[15]=pModel->GetStateVarIndex(PONDED_WATER); iTo[15]=pModel->GetLakeStorageIndex();
    if (pModel->StateVarExists(DEPRESSION)){
      iFrom[16]=pModel->GetStateVarIndex(SNOW);         iTo[16]=pModel->GetStateVarIndex(DEPRESSION);
      iFrom[17]=pModel->GetStateVarIndex(SNOW_LIQ);     iTo[17]=pModel->GetStateVarIndex(DEPRESSION);
      iFrom[18]=pModel->GetStateVarIndex(PONDED_WATER); iTo[18]=pModel->GetStateVarIndex(DEPRESSION);
    }
  }
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

  int    iCan,iSCan,iSnLiq,iLake,iSWE,iPond;
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
  iSWE    =pModel->GetStateVarIndex(SNOW);
  iPond   =pModel->GetStateVarIndex(PONDED_WATER);
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
  int qSW=10;
  int qColdContent=11;
  int qSnowDepth=12;
  int qSnowToLake=13;
  int qSnLiqToLake=14;
  int qPondToLake=15;
  int qSnowToDep=16;
  int qSnLiqToDep=17;
  int qPondToDep=18;

  //get precipitation information from gauges
  total_precip=pHRU->GetForcingFunctions()->precip;//[mm/day]
  Fsnow       =pHRU->GetForcingFunctions()->snow_frac;
  air_temp    =pHRU->GetForcingFunctions()->temp_ave;

  //identify snow/rain fraction
  if (!pModel->StateVarExists(SNOW)){Fsnow=0;}
  snowfall=(    Fsnow)*total_precip;
  rainfall=(1.0-Fsnow)*total_precip;//[mm/day]

  rainfall+=pHRU->GetForcingFunctions()->irrigation;//[mm/day]

  double SWE=0.0;
  if(pModel->StateVarExists(SNOW)) {
    SWE=state_vars[pModel->GetStateVarIndex(SNOW)];
  }
  // Calculate Snow Cover
  //-----------------------------------------------------------------
  Fs = pHRU->GetSnowCover();

  if (!pHRU->IsLake())
  {
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
  }
  else //lake - no vegetation interception
  {
    snowthru = snowfall;
    rainthru = rainfall;
  }

  bool lake_frozen   =false;
  bool wetland_frozen=false;
  if (pModel->StateVarExists(ICE_THICKNESS))
  {
    int iIceThick=pModel->GetStateVarIndex(ICE_THICKNESS);
    if (state_vars[iIceThick]>10){ //10mm of ice can support snow?
      lake_frozen   =true;
      wetland_frozen=true;
    }
  }

  if ((pHRU->IsLake()) && (!lake_frozen))
  {
    if (pHRU->IsLinkedToReservoir()) {
      rates[qSW]=snowthru+rainthru; //handles case when both reservoirs and lake storage used
    }
    else{
      rates[qLake]=snowthru+rainthru; // in lake, all water just added
    }

    if (pModel->StateVarExists(ICE_THICKNESS))//all remaining snow and ponded water dumped to lake \todo[funct]: handle multilayer snowpack models
    {
       rates[qSnowToLake ]+=state_vars[iSWE  ]/Options.timestep;
       rates[qSnLiqToLake]+=state_vars[iSnLiq]/Options.timestep;
    }
  }
  else if ((pHRU->GetHRUType()==HRU_WETLAND) && (!wetland_frozen))
  {
    rates[qDep] =snowthru + rainthru; //all water goes to depression in wetland

    if (pModel->StateVarExists(ICE_THICKNESS))//all remaining snow and ponded water dumped to wetland/depression
    {
       rates[qSnowToDep ]+=state_vars[iSWE  ]/Options.timestep;
       rates[qSnLiqToDep]+=state_vars[iSnLiq]/Options.timestep;
    }
  }
  else if (pHRU->GetHRUType() == HRU_WATER)
  {
    rates[qSW] =snowthru + rainthru; //all water goes to surface water
  }
  else //regular land, frozen lakes or frozen wetlands - treat snow
  {
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
        if ((pModel->StateVarExists(SNOW_LIQ))  && (!pModel->StateVarExists(SNOW_DEFICIT)) && (SWE > 0.0)) //Snow liquid in model
        {
          double melt_cap = pHRU->GetStateVarMax(iSnLiq, state_vars, Options);
          double fill_rate;
          fill_rate = Fs*threshMin(max(melt_cap - state_vars[iSnLiq], 0.0) / Options.timestep, rainthru, 0.0);
          rates[qSnowLiq] = fill_rate;
          rainthru -= fill_rate;
        }
        if ( (pModel->StateVarExists(SNOW_DEFICIT)) && (SWE > 0.0))  //Snow deficit in model (UBCWM)
        {
          double SWI = pModel->GetGlobalParams()->GetParams()->snow_SWI;
          rates[qSnowDef] = SWI*snowthru; //snowfall to snowpack
        }
      }
    }
    else{
      rates[qNewSnow]=snowthru;
    }

    if (pHRU->IsLake())
    {
      if (pHRU->IsLinkedToReservoir()) {
        rates[qSW]=rainthru; //handles case when both reservoirs and lake storage used
      }
      else{
        rates[qLake]=rainthru; // in lake, all water just added
      }

      if (pModel->StateVarExists(ICE_THICKNESS))
      {
         rates[qPondToLake ]+=state_vars[iPond ]/Options.timestep; //handles on-lake snowmelt whether ice is present or not
      }
    }
    else if (pHRU->GetHRUType() == HRU_WETLAND)
    {
      rates[qDep] =rainthru; //ponded water goes to depression in wetland

      if (pModel->StateVarExists(ICE_THICKNESS))
      {
         rates[qPondToDep ]+=state_vars[iPond ]/Options.timestep; //handles on-wetland snowmelt whether ice is present or not
      }
    }
    else{
      //everything remaining goes to ponded water variable, usually handled by infiltration and/or abstraction algorithms
      rates[qPond]=rainthru;
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \param *state_vars [in] array of state variable values for this HRU (of size CModel::_nStateVars)
/// \param *pHRU [in] Pointer to HRU object
/// \param &Options [in] Global model options information
/// \param &tt [in] Current input time structure
/// \param *rates [out] Rates of change of state variables
//
void   CmvPrecipitation::ApplyConstraints(const double     *state_vars,
                                          const CHydroUnit *pHRU,
                                          const optStruct  &Options,
                                          const time_struct &tt,
                                          double     *rates) const
{
  //handled in RatesOfChange routine
}
