/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2022 the Raven Development Team
  ------------------------------------------------------------------
  Mass loading of soluble constituent
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "MassLoading.h"
#include "Model.h"

//////////////////////////////////////////////////////////////////
/// \brief Implentation of the MassLoading constructor
/// \param constit_name [in] name of constituent being tracked
/// \param pTransportModel [in] Transport Model object
//
CmvMassLoading::CmvMassLoading(string constit_name,
                           CTransportModel *pTransportModel)
  :CHydroProcessABC(MASS_LOADING)
{
  pTransModel=pTransportModel;
  _constit_ind=pTransModel->GetConstituentIndex(constit_name);

  int nWaterCompartments =pTransModel->GetNumWaterCompartments();
  CHydroProcessABC::DynamicSpecifyConnections(nWaterCompartments);


  for (int ii = 0; ii < nWaterCompartments; ii++) { //Neumann sources
    iFrom[ii]=pModel->GetStateVarIndex(CONSTITUENT_SRC,_constit_ind);
    int m=pTransModel->GetLayerIndex(_constit_ind,pTransModel->GetStorWaterIndex(ii));
    iTo  [ii]=pModel->GetStateVarIndex(CONSTITUENT,m);
    
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvMassLoading::~CmvMassLoading(){}

//////////////////////////////////////////////////////////////////
/// \brief Initializes mass loading object
//
void   CmvMassLoading::Initialize()
{
}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for baseflow algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by advection algorithm (size of aP[] and aPC[])
//
void CmvMassLoading::GetParticipatingParamList(string  *aP, class_type *aPC, int &nP) const
{
  nP=0; //All parameters are externally determined
}

//////////////////////////////////////////////////////////////////
/// \brief For modeling changes contaminant mass/concentration linked to a flow process
///
/// \param *state_vars [in] array of current state variables
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] rates of change due to both associated flow process and advective transport
//
void   CmvMassLoading::GetRatesOfChange(const double      *state_vars,
                                        const CHydroUnit  *pHRU,
                                        const optStruct   &Options,
                                        const time_struct &tt,
                                        double            *rates) const
{
  int    k=pHRU->GetGlobalIndex();
  
  //Handle Neumann influx conditions, if present
  //-------------------------------------------------------
  int iTo;
  for (int ii = 0; ii < pTransModel->GetNumWaterCompartments(); ii++)
  {
    iTo=pTransModel->GetStorWaterIndex(ii);
    rates[ii] = pTransModel->GetConstituentModel2(_constit_ind)->GetSpecifiedMassFlux(iTo, k, tt); //[mg/m2/d] or [MJ/m2/d]
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
void   CmvMassLoading::ApplyConstraints(const double               *state_vars,
                                        const CHydroUnit *pHRU,
                                        const optStruct      &Options,
                                        const time_struct &tt,
                                        double     *rates) const
{
  //does nothing
}
