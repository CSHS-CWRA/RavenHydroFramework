/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2018 the Raven Development Team
  ------------------------------------------------------------------
  Transformation of soluble contaminant/tracer/nutrient into another constituent
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "Decay.h"
#include "Model.h"
#include "StateVariables.h"

//////////////////////////////////////////////////////////////////
/// \brief Implentation of the Transfomation constructor
/// \param constit_name [in] name of reactant constituent
/// \param constit_name [in] name of product constituent
/// \param ttyp [in] transformation process type
/// \param pTransportModel [in] transport Model object
//
CmvTransformation::CmvTransformation(string           constit_name,
                                     string           constit_name2,
                                     transformation_type   ttyp,
                                     CTransportModel *pTransportModel)
  :CHydroProcessABC(TRANSFORMATION)
{
  _ttype=ttyp;
  _pTransModel=pTransportModel;
  _constit_ind1 =_pTransModel->GetConstituentIndex(constit_name);
  _constit_ind2=_pTransModel->GetConstituentIndex(constit_name2);
  ExitGracefullyIf(_constit_ind1==DOESNT_EXIST,
                   "CmvTransformation constructor: invalid constituent name in :Transformation command",BAD_DATA_WARN);
  ExitGracefullyIf(_constit_ind2==DOESNT_EXIST,
                   "CmvTransformation constructor: invalid second constituent name in :Transformation command",BAD_DATA_WARN);

  int nWaterCompartments = _pTransModel->GetNumWaterCompartments();
  CHydroProcessABC::DynamicSpecifyConnections(2*nWaterCompartments);

  //transformation occurs in all water storage compartments
  for (int ii=0;ii<nWaterCompartments;ii++)
  {
    iFrom[ii                   ]=_pTransModel->GetStorIndex(_constit_ind1,ii); //mass in water compartment
    iTo  [ii                   ]=_pTransModel->GetStorIndex(_constit_ind2,ii); //mass in water compartment
    iFrom[ii+nWaterCompartments]=_pTransModel->GetStorIndex(_constit_ind1,ii); //mass in water compartment
    iTo  [ii+nWaterCompartments]=pModel->GetStateVarIndex(CONSTITUENT_SINK,_constit_ind1); //'loss/sink' storage (for MB accounting)
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvTransformation::~CmvTransformation(){}

//////////////////////////////////////////////////////////////////
/// \brief Initializes Advection object
//
void   CmvTransformation::Initialize(){}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for baseflow algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by advection algorithm (size of aP[] and aPC[])
//
void CmvTransformation::GetParticipatingParamList(string  *aP, class_type *aPC, int &nP) const
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
void   CmvTransformation::GetRatesOfChange( const double      *state_vars,
                                            const CHydroUnit  *pHRU,
                                            const optStruct   &Options,
                                            const time_struct &tt,
                                            double            *rates) const
{
  int     k=pHRU->GetGlobalIndex();
  double  junk,mass1,mass2,vol1;
  int     iStor,iConstit1,iConstit2;
  
  double  transf_coeff;
  double  stoich_coeff; //stoichiometric coefficient, i.e., 1*A->alpha*B
  double  n;
  
  int     nWaterCompartments = _pTransModel->GetNumWaterCompartments();
  for (int ii = 0; ii < nWaterCompartments; ii++)
  {
    iStor    =_pTransModel->GetStorWaterIndex(ii);
    iConstit1=_pTransModel->GetStorIndex     (_constit_ind1,ii);   //global state variable index of reactant constituent in this water storage
    iConstit2=_pTransModel->GetStorIndex     (_constit_ind2,ii);   //global state variable index of product constituent in this water storage

    mass1=state_vars[iConstit1];
    mass2=state_vars[iConstit2];
    vol1 =state_vars[ii];

    transf_coeff = _pTransModel->GetTransformCoefficient(_constit_ind1,_constit_ind2,pHRU,iStor);
    stoich_coeff = _pTransModel->GetStoichioCoefficient (_constit_ind1,_constit_ind2,pHRU,iStor);
    n            = 1.01;// _pTransModel->GetTransformPower      (_constit_ind1,_constit_ind2,pHRU,iStor);

    if(_pTransModel->IsDirichlet(iStor,_constit_ind1,k,tt,junk)){} //don't reduce dirichlet source zones
    else {
      if (_ttype==TRANSFORM_LINEAR)
      {
        rates[ii]= -transf_coeff*stoich_coeff*mass1; 
        rates[ii+nWaterCompartments]=rates[ii]*(1-stoich_coeff)/stoich_coeff; //mass transformed to something else (sink)

        // dA/dt = - k * A
        // dB/dt = + k * A * s
        // dsink/dt = + k * A * (1-s) = dB/dt *(1-s)/s
      }
      else if (_ttype==TRANSFORM_LINEAR_ANALYTIC)//analytical approach - definitely preferred - solution to dm/dt=-km integrated from t to t+dt
      {
        rates[ii] = mass1*stoich_coeff * (1 - exp(-transf_coeff*Options.timestep))/Options.timestep; 
        rates[ii+nWaterCompartments]=rates[ii]*(1-stoich_coeff)/stoich_coeff; //mass transformed to something else (sink)
      }
      else if (_ttype==TRANSFORM_NONLINEAR)
      {
        rates[ii] = -transf_coeff*stoich_coeff *pow(mass1/vol1,n)*vol1; 
        rates[ii+nWaterCompartments]=rates[ii]*(1-stoich_coeff)/stoich_coeff; //mass transformed to something else (sink)
        ExitGracefully("TRANSFORM_NONLINEAR - need way of storing n",STUB);
      }
      else if (_ttype==TRANSFORM_NONLIN_ANALYTIC)//analytical approach - definitely preferred - solution to dm/dt=-km integrated from t to t+dt
      {
        double C0=mass1/vol1;
        rates[ii] = stoich_coeff *(C0- pow(pow(C0,1.0-n)+(n-1)*transf_coeff*Options.timestep,1.0/(1.0-n)))/Options.timestep*vol1; 
        rates[ii+nWaterCompartments]=rates[ii]*(1-stoich_coeff)/stoich_coeff; //mass transformed to something else (sink)
        ExitGracefully("TRANSFORM_NONLINEAR - need way of storing n",STUB);
      }
      else if(_ttype==TRANSFORM_MINERALIZATION)
      {
        double temp=pHRU->GetForcingFunctions()->temp_ave;
        double c1=1.0;
        double sat=1.0;
        double k;
        if      ((temp<=5 ) || (temp>=50)) { c1=0.0; }
        else if ((temp>=10) && (temp<=30)) { c1=1.0; }
        else if (temp< 10)                 { c1=(temp-5.0)/(10.0-5.0);}
        else if (temp> 30)                 { c1=(30-temp)/(30-50); }
        double vol_max=pHRU->GetStateVarMax(ii,state_vars,Options);

        sat=min(vol1/vol_max,1.0);
        double rn_min =10;// = _pTransModel->GetMineralizationCoefficient(_constit_ind1,_constit_ind2,pHRU,iStor);
        k = rn_min * c1 * sat;
        ExitGracefully("TRANSFORM_MINERALIZATION - need way of storing mineralization coeff",STUB);
        rates[ii] = mass1*(1 - exp(-k*Options.timestep))/Options.timestep;
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
void   CmvTransformation::ApplyConstraints( const double      *state_vars,
                                            const CHydroUnit  *pHRU,
                                            const optStruct   &Options,
                                            const time_struct &tt,
                                            double            *rates) const
{
  int iConstit1;
  int nWaterCompartments = _pTransModel->GetNumWaterCompartments();
  for (int ii = 0; ii < nWaterCompartments; ii++)
  {
    iConstit1=_pTransModel->GetStorIndex(_constit_ind1,ii);
    rates[ii] = min(rates[ii],state_vars[iConstit1]/Options.timestep);//cannot remove more mass than is there
  }
}

