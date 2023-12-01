/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2022 the Raven Development Team
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
/// \param iWatStor [in] state variable index of water storage compartment
/// \param pTransportModel [in] transport Model object
//
CmvTransformation::CmvTransformation(string           constit_name,
                                     string           constit_name2,
                                     transformation_type   ttyp,
                                     int              proc_ind,
                                     int              iWatStor,
                                     CTransportModel *pTransportModel,
                                     CModelABC       *pModel)
  :CHydroProcessABC(TRANSFORMATION, pModel)
{
  _ttype        =ttyp;
  _pTransModel  =pTransportModel;
  _process_ind  =proc_ind;
  _iWaterStore  =iWatStor;
  _constit_ind1 =_pTransModel->GetConstituentIndex(constit_name);
  _constit_ind2 =_pTransModel->GetConstituentIndex(constit_name2);
  ExitGracefullyIf(_constit_ind1==DOESNT_EXIST,
                   "CmvTransformation constructor: invalid constituent name in :Transformation command",BAD_DATA_WARN);
  ExitGracefullyIf(_constit_ind2==DOESNT_EXIST,
                   "CmvTransformation constructor: invalid second constituent name in :Transformation command",BAD_DATA_WARN);

  int m,m2,iWat;
  if(_iWaterStore==DOESNT_EXIST) { //occurs everywhere
    int nWaterCompartments = _pTransModel->GetNumWaterCompartments();
    CHydroProcessABC::DynamicSpecifyConnections(2*nWaterCompartments);

    //transformation occurs in all water storage compartments
    for(int ii=0;ii<nWaterCompartments;ii++)
    {
      iWat=_pTransModel->GetStorWaterIndex(ii);
      m =_pTransModel->GetLayerIndex(_constit_ind1,iWat);
      m2=_pTransModel->GetLayerIndex(_constit_ind2,iWat);

      iFrom[ii                   ]=pModel->GetStateVarIndex(CONSTITUENT,m);//mass in water compartment
      iTo  [ii                   ]=pModel->GetStateVarIndex(CONSTITUENT_SINK,_constit_ind1); //'loss/sink' storage (for MB accounting)
      iFrom[ii+nWaterCompartments]=pModel->GetStateVarIndex(CONSTITUENT_SRC ,_constit_ind1); //'loss/sink' storage (for MB accounting)
      iTo  [ii+nWaterCompartments]=pModel->GetStateVarIndex(CONSTITUENT,m2);//mass in water compartment
    }
  }
  else { //only occurs in one water storage unit
    CHydroProcessABC::DynamicSpecifyConnections(2);

    m =_pTransModel->GetLayerIndex(_constit_ind1,_iWaterStore);
    m2=_pTransModel->GetLayerIndex(_constit_ind2,_iWaterStore);
    iFrom[0]=pModel->GetStateVarIndex(CONSTITUENT,m);                 //mass in water compartment
    iTo  [0]=pModel->GetStateVarIndex(CONSTITUENT_SRC,_constit_ind1); //'loss/sink' storage (for MB accounting)
    iFrom[1]=pModel->GetStateVarIndex(CONSTITUENT_SINK,_constit_ind1);//'loss/sink' storage (for MB accounting)
    iTo  [1]=pModel->GetStateVarIndex(CONSTITUENT,m2);                //mass in water compartment
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

  if(_ttype==TRANS_LINEAR)
  {
    aP[0]="TRANSFORM_COEFF"; aPC[0]=CLASS_TRANSPORT; nP++;
    aP[1]="STOICHIO_RATIO";  aPC[1]=CLASS_TRANSPORT; nP++;
  }
  else if(_ttype==TRANS_NONLINEAR)
  {
    aP[0]="TRANSFORM_COEFF"; aPC[0]=CLASS_TRANSPORT; nP++;
    aP[1]="STOICHIO_RATIO";  aPC[1]=CLASS_TRANSPORT; nP++;
    aP[2]="TRANSFORM_N";     aPC[2]=CLASS_TRANSPORT; nP++;
  }

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
  double  junk,mass1,vol1;
  int     iStor;

  double  transf_coeff;
  double  stoich_ratio; //stoichiometric coefficient, i.e., 1*A->alpha*B

  int     nWaterCompartments = _pTransModel->GetNumWaterCompartments();
  int     ii_active          = _pTransModel->GetWaterStorIndexFromSVIndex(_iWaterStore);

  int     shift =nWaterCompartments;
  if(_iWaterStore!=DOESNT_EXIST) { shift=1; }

  int     q=0; //connection index
  for(int ii = 0; ii < nWaterCompartments; ii++)
  {
    iStor    =_pTransModel->GetStorWaterIndex(ii);

    if((_iWaterStore!=DOESNT_EXIST) && (ii!=ii_active)) { continue; } //only apply to one water compartment

    if(_pTransModel->GetConstituentModel2(_constit_ind1)->IsDirichlet(iStor,k,tt,junk)) { continue; } //don't modify dirichlet source zones

    mass1=state_vars[iFrom[q]]; //mg/m2
    vol1 =state_vars[iStor   ]; //mm

    transf_coeff = _pTransModel->GetGeochemParam(PAR_TRANSFORM_COEFF,_constit_ind1,_constit_ind2,ii,_process_ind,pHRU);
    if(transf_coeff==NOT_SPECIFIED) { continue; }

    stoich_ratio = _pTransModel->GetGeochemParam(PAR_STOICHIO_RATIO,_constit_ind1,_constit_ind2,ii,_process_ind,pHRU);
    if(stoich_ratio==NOT_SPECIFIED) { continue; }

    // for each reaction,
    // dA/dt      = - k * A
    // dsink/dt   = + k * A
    // dB/dt      = + k * A * s
    // dsource/dt = - k * A * s

    //-------------------------------------------------------------------------------
    if (_ttype==TRANS_LINEAR)//analytical approach -  solution to dm/dt=-km integrated from t to t+dt
    {                        //numerically robust equiv to rates[q]= transf_coeff*mass1;

      rates[q] = mass1 * (1 - exp(-transf_coeff*Options.timestep))/Options.timestep;

    }
    //-------------------------------------------------------------------------------
    else if (_ttype==TRANS_NONLINEAR)//analytical approach -  solution to dm/dt=-km integrated from t to t+dt
    {
      double C0=mass1/vol1; //mg/m2/mm
      C0=C0*MM_PER_METER/LITER_PER_M3; //mg/m2/mm-->mg/L

      double n = _pTransModel->GetGeochemParam(PAR_TRANSFORM_N,_constit_ind1,_constit_ind2,ii,_process_ind,pHRU);
      if(n==NOT_SPECIFIED) { n=1.001; }

      rates[q] = (C0- pow(pow(C0,1.0-n)+(n-1)*transf_coeff*Options.timestep,1.0/(1.0-n)))/Options.timestep; //mg/L/d
      rates[q] *=vol1*LITER_PER_M3/MM_PER_METER; //mg/L/d->mg/m2/d

      //numerically robust equiv to
      //rates[q] = transf_coeff *pow(mass1/vol1,n)*vol1;
    }

    rates[q+shift]=stoich_ratio*rates[q]; //handles proper accounting

    q++;
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
  int ii_active          = _pTransModel->GetWaterStorIndexFromSVIndex(_iWaterStore);

  int shift =nWaterCompartments;
  if(_iWaterStore!=DOESNT_EXIST) {shift=1;}

  int q=0;
  for (int ii = 0; ii < nWaterCompartments; ii++)
  {
    if(_iWaterStore!=DOESNT_EXIST) { //only apply to one water compartment
      if(ii!=ii_active) { break; }
    }
    double stoich_ratio = _pTransModel->GetGeochemParam(PAR_STOICHIO_RATIO,_constit_ind1,_constit_ind2,ii,_process_ind,pHRU);

    iConstit1=_pTransModel->GetStorIndex(_constit_ind1,ii);
    rates[q] = min(rates[q],state_vars[iConstit1]/Options.timestep);//cannot remove more mass than is there
    rates[q+shift]=-stoich_ratio*rates[q]; //handles proper accounting
    q++;
  }
}
