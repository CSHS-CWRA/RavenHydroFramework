/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
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

  iFrom[0]=pModel->GetStateVarIndex(SNOW);
  iTo  [0]=pModel->GetStateVarIndex(FIRN);
  iFrom[1]=pModel->GetStateVarIndex(FIRN);
  iTo  [1]=pModel->GetStateVarIndex(GLACIER_ICE);
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
    nP=0;
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
    nSV=3;
    aSV[0]=GLACIER_ICE;  aLev[0]=DOESNT_EXIST;
    aSV[1]=FIRN;         aLev[1]=DOESNT_EXIST;
    aSV[2]=SNOW;         aLev[2]=DOESNT_EXIST; //\todo [funct]: multilayer snowpacks
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
  if (pHRU->GetHRUType()!=HRU_GLACIER){return;}

  //----------------------------------------------------------------------
  if (type==FIRNEVOL_SIMPLE)
  {
    double SWE=state_var[iFrom[0]];
    if (SWE>1000){}
    rates[0]=0.001*(SWE-1000);
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
  if (pHRU->GetHRUType()!=HRU_GLACIER){return;}

  if (rates[0]<0.0){rates[0]=0.0;}//positivity constraint on melt

  //cant remove more than is there
  //rates[0]=threshMin(rates[0],state_var[iFrom[0]]/Options.timestep,0.0);
}
