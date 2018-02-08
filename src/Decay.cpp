/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2018 the Raven Development Team
  ------------------------------------------------------------------
  Decay of soluble contaminant/tracer/nutrient
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "Decay.h"
#include "Model.h"
#include "StateVariables.h"

//////////////////////////////////////////////////////////////////
/// \brief Implentation of the Advection constructor
/// \param constituent [in] name of contaminant beign tracked
/// \param pFlow [in] flow process which drives advection (this acts as a wrapper for said process)
/// \param pModel [in] Model object
//
CmvDecay::CmvDecay(string           constit_name,
                   decay_type       dtyp,
                   CTransportModel *pTransportModel)
  :CHydroProcessABC(DECAY)
{
  _dtype=dtyp;
  _pTransModel=pTransportModel;
  _constit_ind=_pTransModel->GetConstituentIndex(constit_name);
  ExitGracefullyIf(_constit_ind==-1,
                   "CmvDecay constructor: invalid constituent name in :Decay command",BAD_DATA_WARN);

  int nWaterCompartments = _pTransModel->GetNumWaterCompartments();

  CHydroProcessABC::DynamicSpecifyConnections(nWaterCompartments);

  //decay occurs in all water storage compartments
  for (int ii=0;ii<nWaterCompartments;ii++)
  {
    iFrom[ii]=_pTransModel->GetStorIndex(_constit_ind,ii); //mass in water compartment
    iTo  [ii]=pModel->GetStateVarIndex(CONSTITUENT_SINK,_constit_ind); //'loss/sink' storage (for MB accounting)
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvDecay::~CmvDecay(){}

//////////////////////////////////////////////////////////////////
/// \brief Initializes Advection object
//
void   CmvDecay::Initialize(){}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for baseflow algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by advection algorithm (size of aP[] and aPC[])
//
void CmvDecay::GetParticipatingParamList(string  *aP, class_type *aPC, int &nP) const
{
  nP=0;
}

//////////////////////////////////////////////////////////////////
/// \brief For modeling changes contaminant mass/concentration linked to a flow process
///
/// \param *state_vars [in] array of current state variables
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] rates of change in each water storage compartment due to decay
//
void   CmvDecay::GetRatesOfChange(const double      *state_vars,
                                  const CHydroUnit  *pHRU,
                                  const optStruct   &Options,
                                  const time_struct &tt,
                                  double            *rates) const
{
  int    k=pHRU->GetGlobalIndex();
  double junk,mass;
  int iStor,iConstit;
  int nWaterCompartments = _pTransModel->GetNumWaterCompartments();
  double decay_coeff;

  for (int ii = 0; ii < nWaterCompartments; ii++)
  {
    iStor   =_pTransModel->GetStorWaterIndex(ii);
    iConstit=_pTransModel->GetStorIndex     (_constit_ind,ii);   //global state variable index of this constituent in this water storage

    mass=state_vars[iConstit];
    decay_coeff = _pTransModel->GetDecayCoefficient(_constit_ind,pHRU,iStor);

    if(_pTransModel->IsDirichlet(iStor,_constit_ind,k,tt,junk)){} //don't modify dirichlet source zones
    else {
      if (_dtype==DECAY_BASIC)
      {
        rates[ii]= -decay_coeff*mass; //[mg/m2/d]=[1/d]*[mg/m2]
      }
      else if (_dtype==DECAY_ANALYTIC)//analytical approach - definitely preferred - solution to dm/dt=-km integrated from t to t+dt
      {
        rates[ii] = mass * (1 - exp(-decay_coeff*Options.timestep))/Options.timestep; 
      }
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
///
/// \param *state_vars [in] array of current state variables
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] rates of change due to both associated flow process and advective transport
//
void   CmvDecay::ApplyConstraints(const double           *state_vars,
                                  const CHydroUnit *pHRU,
                                  const optStruct      &Options,
                                  const time_struct &tt,
                                  double     *rates) const
{
  int iConstit;
  int nWaterCompartments = _pTransModel->GetNumWaterCompartments();
  for (int ii = 0; ii < nWaterCompartments; ii++)
  {
    iConstit=_pTransModel->GetStorIndex(_constit_ind,ii);
    rates[ii] = min(rates[ii],state_vars[iConstit]/Options.timestep);//cannot remove more mass than is there
  }
}

