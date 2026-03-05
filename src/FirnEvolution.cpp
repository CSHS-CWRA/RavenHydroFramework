/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2026 the Raven Development Team
  ------------------------------------------------------------------
  Firn Evolution
  ----------------------------------------------------------------*/
#include "HydroProcessABC.h"
#include "GlacierProcesses.h"
#include "Model.h"

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the glacier ice melt constructor
/// \param melt_type [in] Model of firn evolution selected
//
CmvFirnEvolution::CmvFirnEvolution(firn_evolution_type fe_type,
                                   CModel             *pModel):
  CHydroProcessABC(FIRN_EVOLUTION, pModel)
{
  type =fe_type;

  CHydroProcessABC::DynamicSpecifyConnections(3);

  iFrom[0]=pModel->GetStateVarIndex(SNOW);
  iTo  [0]=pModel->GetStateVarIndex(FIRN);
  iFrom[1]=pModel->GetStateVarIndex(FIRN);
  iTo  [1]=pModel->GetStateVarIndex(GLACIER_ICE);
  iFrom[2]=pModel->GetStateVarIndex(FIRN_GRAVITY);
  iTo  [2]=pModel->GetStateVarIndex(FIRN_GRAVITY);
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvFirnEvolution::~CmvFirnEvolution(){}

//////////////////////////////////////////////////////////////////
/// \brief Initializes class
//
void CmvFirnEvolution::Initialize(){}


//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by algorithm (size of aP[] and aPC[])
//
void CmvFirnEvolution::GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const
{
  if (type==FIRNEVOL_SIMPLE)
  {
    nP=1;
    aP[0]="FIRN_COMPACTION_RATE"; aPC[0]=CLASS_LANDUSE;
  }
  else
  {
    ExitGracefully("CmvFirnEvolution::GetParticipatingParamList: undefined algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief returns list of state variables which are used by firn evolution algorithm
///
/// \param melt_type [in] Algorithm for firn evolution used
/// \param *aSV [out] Array of state variable types needed by firn evolution algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by firn evolution algorithm (size of aSV[] and aLev[] arrays)
//
void CmvFirnEvolution::GetParticipatingStateVarList(firn_evolution_type type,sv_type *aSV, int *aLev, int &nSV)
{
  if (type==FIRNEVOL_SIMPLE)
  {
    nSV=4;
    aSV[0]=GLACIER_ICE;  aLev[0]=DOESNT_EXIST;
    aSV[1]=FIRN;         aLev[1]=DOESNT_EXIST;
    aSV[2]=SNOW;         aLev[2]=DOESNT_EXIST; //\todo [funct]: multilayer snowpacks
    aSV[3]=FIRN_GRAVITY; aLev[3]=DOESNT_EXIST;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rates of melting of glacial ice [mm/day]
/// \details  if type==GMELT_DEGREE_DAY
///             <ul> <li> melt is proportional to temperature </ul>
/// if type==GMELT_UBC
///             <ul> <li> melt is proportional to storage; constant of proportionality depends upon surface snow/ice </ul>
///
/// \param *state_var [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which rates of change are calculated
/// \param *rates [out] Array of rates of change of modified state variables
///
//
void CmvFirnEvolution::GetRatesOfChange( const double             *state_var,
                                       const CHydroUnit  *pHRU,
                                       const optStruct   &Options,
                                       const time_struct &tt,
                                       double      *rates) const
{

  //----------------------------------------------------------------------
  if (type==FIRNEVOL_SIMPLE)
  {
    double snow_to_firn=0.0;    //[mm/d]
    double firn_to_glacier=0.0; //[mm/d]
    double firn_grav_change=0.0;//[1/d]

    double firnPct=1.0;               //percentage of snow remaining in summer converted to firn \todo[funct] - make this a parameter
    const double PACKED_GRAVITY=0.5;  //500 kg/m3 - density of snow converted to firn
    const double THRESH_GRAVITY=0.83; //830 kg/m3 - density at which firn converts to ice (pore close off density)
    const double DENS_GRADIENT =0.000005;//density gradient of firn [1/mm] (assumes ~300 kg/m3 over ~60m)

    double SWE =state_var[iFrom[0]];
    double firn=state_var[iFrom[1]];
    double grav=state_var[iFrom[2]];

    double firn_compaction_rate=pHRU->GetSurfaceProps()->firn_compaction_rate/DAYS_PER_YEAR;
    //double firn_compaction_rate=0.0000362; //[1/d] ~0.33 m/yr for 25m depth (~1.3%/yr)

double firn_old=firn;

    //if snow still present at end of summer, convert percentage to firn
    double grav_new;
    bool isMidnight=(tt.model_time-rvn_floor(tt.model_time+TIME_CORRECTION)<REAL_SMALL);
    if (((pHRU->GetLatRad()>=0.0) && (tt.julian_day==213) && (isMidnight)) || //Aug 1 in northern hemisphere
        ((pHRU->GetLatRad()< 0.0) && (tt.julian_day==32 ) && (isMidnight)))   //Feb 1 in southern hemisphere
    {
      snow_to_firn=(firnPct*SWE)/Options.timestep; //mm/d
      if ((firn+firnPct*SWE)>0.0){
        grav_new=(grav*(firn)+PACKED_GRAVITY*(firnPct*SWE))/(firn+firnPct*SWE);
      }
      else{
        grav_new=PACKED_GRAVITY;
      }
      firn+=snow_to_firn*Options.timestep;
      firn_grav_change=(grav_new-grav)/Options.timestep; //assume deep snow is at 500 kg/m3
      grav=grav_new;
    }

    double bottom_grav=(grav+0.5*DENS_GRADIENT*firn);
    double top_grav   =(grav-0.5*DENS_GRADIENT*firn);
    //all firn > density threshold converted to glacier ice
    if (THRESH_GRAVITY<bottom_grav)
    {
      firn_to_glacier=firn*(1.0-THRESH_GRAVITY/bottom_grav)/Options.timestep;
      firn-=firn_to_glacier*Options.timestep;
      grav_new=0.5*(THRESH_GRAVITY+top_grav);
      firn_grav_change=(grav_new-grav)/Options.timestep; //assume deep snow is at 500 kg/m3
      grav=grav_new;
    }

    //compact firn
    firn_grav_change+=firn_compaction_rate; //linear compaction rate [1/d]

    //cout<<" "<<pHRU->GetHRUID()<<" "<<snow_to_firn<<" "<<firn_to_glacier<<" "<<firn_grav_change<<" firn new/old: "<<firn<<" "<<firn_old<<" grav: "<<grav<<" "<<bottom_grav<<" "<<top_grav<<endl;
    rates[0]=snow_to_firn;     //SNOW->FIRN
    rates[1]=firn_to_glacier;  //FIRN->GLACIER_ICE
    rates[2]=firn_grav_change; //FIRN_GRAVITY change
  }
};

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow cannot drain "from" compartment over timestep
///
/// \param *state_var [in] Reference to set of participating state variables
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which rates of change are calculated
/// \param *rates [out] Corrected rates of state variable change
//
void   CmvFirnEvolution::ApplyConstraints( const double           *state_var,
                                         const CHydroUnit  *pHRU,
                                         const optStruct         &Options,
                                         const time_struct &tt,
                                         double      *rates) const
{
}
