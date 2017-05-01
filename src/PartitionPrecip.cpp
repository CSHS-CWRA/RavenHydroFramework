
/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright © 2008-2014 the Raven Development Team
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
  CHydroProcessABC::DynamicSpecifyConnections(pModel->GetNumStateVars());
	//NOTE: nConnections==nStateVar 
	//(all water storage variables are potential receptacles for precip)
  /// \todo [funct] this should be changed somehow to all water vars
  for (int q=0;q<nConnections;q++)//all state var
  {
		iFrom[q]    =pModel->GetStateVarIndex(ATMOS_PRECIP);
		iTo  [q]    =q;
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
  int iCan=pModel->GetStateVarIndex(CANOPY);
  bool canopy_exists = (iCan != DOESNT_EXIST);
  int iSnowCan = pModel->GetStateVarIndex(CANOPY);
  bool cansnow_exists = (iSnowCan != DOESNT_EXIST);

  nP=0;
  aP[nP]="SAI_HT_RATIO";      aPC[nP]=CLASS_VEGETATION; nP++; //reasonable defaults
  aP[nP]="FOREST_SPARSENESS"; aPC[nP]=CLASS_LANDUSE; nP++; //reasonable defaults
  aP[nP]="MAX_LAI";           aPC[nP]=CLASS_VEGETATION; nP++;
  if (canopy_exists)
  {
    aP[nP] = "MAX_CAPACITY";  aPC[nP] = CLASS_VEGETATION; nP++;
    aP[nP] = "RELATIVE_LAI";  aPC[nP] = CLASS_VEGETATION; nP++; //reasonable defaults  
    if (pModel->GetOptStruct()->interception_factor == PRECIP_ICEPT_LAI)
    {
      aP[nP] = "RAIN_ICEPT_FACT";   aPC[nP] = CLASS_VEGETATION; nP++;
    }
    else if (pModel->GetOptStruct()->interception_factor == PRECIP_ICEPT_USER)
    {
      aP[nP] = "RAIN_ICEPT_PCT";    aPC[nP] = CLASS_VEGETATION; nP++;
    }
    else if (pModel->GetOptStruct()->interception_factor == PRECIP_ICEPT_EXPLAI)
    {
      // no parameter required
    }
    else if (pModel->GetOptStruct()->interception_factor == PRECIP_ICEPT_HEDSTROM)
    {
      // no parameter required
    }
    
  }
  if (cansnow_exists)
  {
    aP[nP] = "MAX_SNOW_CAPACITY"; aPC[nP] = CLASS_VEGETATION; nP++;

    if (pModel->GetOptStruct()->interception_factor == PRECIP_ICEPT_LAI)
    {
      aP[nP] = "SNOW_ICEPT_FACT";   aPC[nP] = CLASS_VEGETATION; nP++;
    }
    else if (pModel->GetOptStruct()->interception_factor == PRECIP_ICEPT_USER)
    {
      aP[nP] = "SNOW_ICEPT_PCT";    aPC[nP] = CLASS_VEGETATION; nP++;
    }
    else if (pModel->GetOptStruct()->interception_factor==PRECIP_ICEPT_EXPLAI)
    {
	  // no parameter required
    }
    else if (pModel->GetOptStruct()->interception_factor==PRECIP_ICEPT_HEDSTROM)
    {
	  // no parameter required
    }
    
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
///	Partitioning depends upon which storage units exist.
/// - precip is always treated constant over a given time step, 
/// - though partioning may vary based upon current system storage \n
///	returns precipitation (in mm/day) moving to each storage unit 
///	over time step (stored in rates[]).
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time structure
/// \param *rates [out] rates of precip to various storage compartments [mm/d]
//
void CmvPrecipitation::GetRatesOfChange(const double		 *state_vars, 
																				const CHydroUnit *pHRU, 
																				const optStruct	 &Options,
																				const time_struct &tt,
																				      double     *rates) const
{

  int    iPond, inewSnow;
	int    iCan,iSCan,iMelt,iLake;
  int    iSnow,iSnoD,iSnoT,iCC,iSnFrac;
  int    iSnDef,iAtm;
	double total_precip,rainfall,snowfall;			//all [mm/day]
	double rain_int,snow_int,rainthru,snowthru;	//all [mm/day]
	double capacity,snow_capacity;							//[mm]
	double rain_pct,snow_pct;										//[0..1]
	double Fcan;											          //[0..1] fraction canopy covered
  double Fs;                                  //[0..1] fraction snow covered
  double Fsnow;                               //[0..1] percent precip as snow
	double Temp;																//[C]
	double tstep = Options.timestep;

	iPond   =pModel->GetStateVarIndex(PONDED_WATER);
  inewSnow=pModel->GetStateVarIndex(NEW_SNOW);
  iCan	  =pModel->GetStateVarIndex(CANOPY);
	iSCan	  =pModel->GetStateVarIndex(CANOPY_SNOW);
	iSnow   =pModel->GetStateVarIndex(SNOW);
  iSnFrac =pModel->GetStateVarIndex(SNOW_COVER);
	iMelt	  =pModel->GetStateVarIndex(SNOW_LIQ);
	iSnoD   =pModel->GetStateVarIndex(SNOW_DEPTH);
  iCC     =pModel->GetStateVarIndex(COLD_CONTENT);
  iSnoT   =pModel->GetStateVarIndex(SNOW_TEMP);
  iSnDef  =pModel->GetStateVarIndex(SNOW_DEFICIT);
  iAtm    =pModel->GetStateVarIndex(ATMOSPHERE);

  iLake   =pModel->GetLakeStorageIndex();

  //get precipitation information from gauges
	total_precip=pHRU->GetForcingFunctions()->precip;//[mm/day]
	Fsnow       =pHRU->GetForcingFunctions()->snow_frac;
	Temp        =pHRU->GetForcingFunctions()->temp_ave; 

  //identify snow/rain fraction
	if (!pModel->StateVarExists(SNOW)){Fsnow=0;}
	snowfall=(    Fsnow)*total_precip;
	rainfall=(1.0-Fsnow)*total_precip;//[mm/day]

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

    if (pModel->StateVarExists(CANOPY))
    {
      //total rain interception rate over canopy covered portion of HRU
      rain_pct = pHRU->GetVegVarProps()->rain_icept_pct;
      rain_int = rainfall*rain_pct;
      
      capacity = pHRU->GetVegVarProps()->capacity;
      rain_int = threshMin(rain_int, max((capacity - state_vars[iCan]) / tstep, 0.0), 0.0);

      rates[iCan] = Fcan*rain_int;
    }
    //calculate snow interception-----------------------------------
    if (pModel->StateVarExists(CANOPY_SNOW))
    {
      //rate of snow interception over canopy portion of HRU
      snow_pct = pHRU->GetVegVarProps()->snow_icept_pct;
      snow_int = snowfall*snow_pct;

      snow_capacity = pHRU->GetVegVarProps()->snow_capacity;
      snow_int = max(threshMin(snow_int, max((snow_capacity - state_vars[iSCan]) / tstep, 0.0), 0.0), 0.0);
      rates[iSCan] = Fcan*snow_int;
    }
    /// \todo[funct] - if no canopy SV is present, snow_pct/rain_pct should move to ATMOSPHERE?
    snowthru = snowfall - Fcan*snow_int;
    rainthru = rainfall - Fcan*rain_int;//[mm/day]

    if(!pModel->StateVarExists(NEW_SNOW))
    {
      if (pModel->StateVarExists(SNOW))
			{
				rates[iSnow]=snowthru;
				if (pModel->StateVarExists(SNOW_DEPTH)){
				  rates[iSnoD] = snowthru*DENSITY_ICE / CalcFreshSnowDensity(Temp);
				}
				if (iCC != DOESNT_EXIST){
				  rates[iCC] = snowthru*threshPositive(FREEZING_TEMP - Temp)*HCP_ICE*MJ_PER_J / DENSITY_WATER;//[mm/d]*[K]*[J/m3/K]*[MJ/J]*[m3/kg]=MJ*mm/kg/d
				}
				
				//move rain through to snowmelt--------------------------------
				if ( (iMelt != DOESNT_EXIST) && (state_vars[iSnow] > 0.0) && (iSnDef == DOESNT_EXIST)) //Snow liquid in model
				{
					double melt_cap = pHRU->GetStateVarMax(iMelt, state_vars, Options);
					double fill_rate;
					fill_rate = Fs*threshMin(max(melt_cap - state_vars[iMelt], 0.0) / Options.timestep, rainthru, 0.0);
					//fill_rate=Fs*rainthru; //HBV implementation
					rates[iMelt] = fill_rate;
					rainthru -= fill_rate;
				}
				if ( (iSnDef != DOESNT_EXIST) && (state_vars[iSnow] > 0.0))  //Snow deficit in model
				{
					double SWI = CGlobalParams::GetParams()->snow_SWI;
					double rain_to_snow = 0;//min(rainthru, SWI*state_vars[iSnow] / Options.timestep); //adds rainfall to snowpack

					rates[iSnDef] = SWI*snowthru; //snowfall to snowpack
					rates[iSnow ] += rain_to_snow;
					rates[iSnDef] -= rain_to_snow;
					rainthru      -= rain_to_snow;
				} 
			}  
    }
    else{
      rates[inewSnow]=snowthru;
    }
    //everything remaining goes to ponded water variable, usually handled by infiltration & abstraction algorithms
    rates[iPond]=rainthru;

  }//End if is lake
  else // in lake, all water just added
  {
    /// \todo should check if iLake==DOESNT_EXIST 
    rates[iLake]=snowfall+rainfall;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Applies constraints to baseflow
/// \details For all methods, ensures that rate of flow cannot drain "from" compartment over timestep.
/// \note Presumes overfilling of "to" compartment is handled using cascade.
///
/// \param *state_vars [in] array of state variable values for this HRU (of size CModel::_nStateVars)
/// \param *pHRU [in] Pointer to HRU object
/// \param &Options [in] Global model options information
/// \param &tt [in] Current input time structure
/// \param *rates [out] Rates of change of state variables (of size MAX_STATE_VARS)
//
void   CmvPrecipitation::ApplyConstraints(const double     *state_vars, 
                                          const CHydroUnit *pHRU, 
                                          const optStruct  &Options,
                                          const time_struct &tt,
                                                double     *rates) const
{
  //handled in RatesOfChange routine
}