/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2022 the Raven Development Team
  ------------------------------------------------------------------
  Decay of soluble contaminant/tracer/nutrient
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "Decay.h"
#include "Model.h"
#include "StateVariables.h"

//////////////////////////////////////////////////////////////////
/// \brief Implentation of the Decay constructor
/// \param constit_name [in] name of decaying constituent
/// \param dtyp [in] decay process type
/// \param pTransportModel [in] transport Model object
//
CmvDecay::CmvDecay(string           constit_name,
                   decay_type       dtyp,
                   int              proc_ind,
                   int              iWatStor,
                   CTransportModel *pTransportModel)
  :CHydroProcessABC(DECAY)
{
  _dtype=dtyp;
  _pTransModel=pTransportModel;
  _constit_ind=_pTransModel->GetConstituentIndex(constit_name);
  _process_ind=proc_ind;
  ExitGracefullyIf(_constit_ind==-1,
                   "CmvDecay constructor: invalid constituent name in :Decay command",BAD_DATA_WARN);

  int nWaterCompartments = _pTransModel->GetNumWaterCompartments();
  int m;
  _iWaterStore=iWatStor;
  if(_iWaterStore==DOESNT_EXIST) { //decay occurs in all water storage compartments

    CHydroProcessABC::DynamicSpecifyConnections(nWaterCompartments);
    for(int ii=0;ii<nWaterCompartments;ii++)
    {
      int iWat=_pTransModel->GetStorWaterIndex(ii);
      m =_pTransModel->GetLayerIndex(_constit_ind,iWat);
      iFrom[ii]=pModel->GetStateVarIndex(CONSTITUENT,m); //mass in water compartment
      iTo  [ii]=pModel->GetStateVarIndex(CONSTITUENT_SINK,_constit_ind); //'loss/sink' storage (for MB accounting)
    }
  }
  else {
    CHydroProcessABC::DynamicSpecifyConnections(1);
    m =_pTransModel->GetLayerIndex(_constit_ind,_iWaterStore);
    iFrom[0]=pModel->GetStateVarIndex(CONSTITUENT,m); //mass in water compartment
    iTo  [0]=pModel->GetStateVarIndex(CONSTITUENT_SINK,_constit_ind); //'loss/sink' storage (for MB accounting)
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvDecay::~CmvDecay(){}

//////////////////////////////////////////////////////////////////
/// \brief Initializes Decay object
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
  aP[0]="DECAY_COEFF"; aPC[0]=CLASS_TRANSPORT; nP++;
}

//////////////////////////////////////////////////////////////////
/// \brief For modeling changes in contaminant mass/concentration linked to a decay process
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
  int     k=pHRU->GetGlobalIndex();
  double  junk,mass;
  int     iStor;
  int     q=0;
  int     nWaterCompartments = _pTransModel->GetNumWaterCompartments();
  int     ii_active          = _pTransModel->GetWaterStorIndexFromSVIndex(_iWaterStore);

  for (int ii = 0; ii < nWaterCompartments; ii++)
  {
    iStor   =_pTransModel->GetStorWaterIndex(ii);
    
    if((_iWaterStore!=DOESNT_EXIST) && (ii!=ii_active)) { continue; } //only apply to one water compartment 
    if(_pTransModel->GetConstituentModel2(_constit_ind)->IsDirichlet(iStor,k,tt,junk)) { continue; } //don't modify dirichlet source zones
    
    mass=state_vars[iFrom[q]];

    //------------------------------------------------------------------
    if      (_dtype==DECAY_BASIC)
    {
      //dm/dt=-km
      double decay_coeff = _pTransModel->GetGeochemParam(PAR_DECAY_COEFF,_constit_ind,ii,_process_ind,pHRU);
      if(decay_coeff==NOT_SPECIFIED) { continue; }

      rates[q]= decay_coeff*mass; 
    }
    //------------------------------------------------------------------
    else if (_dtype==DECAY_ZEROORDER)
    {
      //dm/dt=-J
      double loss_rate = _pTransModel->GetGeochemParam(PAR_MASS_LOSS_RATE,_constit_ind,ii,_process_ind,pHRU);
      if(loss_rate==NOT_SPECIFIED) { continue; }

      rates[q]= min(loss_rate,mass/Options.timestep);
    }
    //------------------------------------------------------------------
    else if (_dtype==DECAY_LINEAR) //analytical approach - preferred - solution to dm/dt=-km integrated from t to t+dt
    {
      //dm/dt=-km
      double decay_coeff = _pTransModel->GetGeochemParam(PAR_DECAY_COEFF,_constit_ind,ii,_process_ind,pHRU);
      if(decay_coeff==NOT_SPECIFIED) { continue; }

      rates[q] = mass * (1.0 - exp(-decay_coeff*Options.timestep))/Options.timestep; 
      //cout<<" DECAY: "<<mass<<" "<<rates[q]<<" "<<decay_coeff*mass<<" "<<
      //CStateVariable::GetStateVarLongName(pModel->GetStateVarType(iStor),pModel->GetStateVarLayer(iStor))<<endl;
    }
    //------------------------------------------------------------------
    else if(_dtype==DECAY_DENITRIF) 
    {
      double decay_coeff = _pTransModel->GetGeochemParam(PAR_DECAY_COEFF,_constit_ind,ii,_process_ind,pHRU);
      if(decay_coeff==NOT_SPECIFIED) { continue; }

      double temp=pHRU->GetForcingFunctions()->temp_ave;
      double c1     =1.0;
      double stor    =state_vars[iStor];
      double stor_max=pHRU->GetStateVarMax(iStor,state_vars,Options);

      if      ((temp<=5 ) || (temp>=50)) { c1=0.0; }
      else if ((temp>=10) && (temp<=30)) { c1=1.0; }
      else if (temp< 10)                 { c1=(temp-5.0 )/(10.0-5.0 );}
      else if (temp> 30)                 { c1=(30.0-temp)/(30.0-50.0);}
      
      decay_coeff *= c1 * min(stor/stor_max,1.0);

      rates[q] = mass * (1.0 - exp(-decay_coeff*Options.timestep))/Options.timestep;
    }
    //------------------------------------------------------------------
    /*else if(_dtype==DECAY_O2EXCHANGE)
    {
      //exchange of constituent (e.g., O2) with atmosphere
      //depression storage/surface water only
      double temp=pHRU->GetForcingFunctions()->temp_ave;
      sv_type typ=pModel->GetStateVarType(iStor);

      if((typ==SURFACE_WATER) || (typ==DEPRESSION) || (typ==PONDED_WATER)) {
        double C_sat = 10; //mg/l (should be function of temperature)
        C_sat*=LITER_PER_M3; //mg/m3

        double exch_coeff = _pTransModel->GetGeochemParam(PAR_O2_EXCHANGE,_constit_ind,ii,_process_ind,pHRU);
        if(exch_coeff==NOT_SPECIFIED) { continue; } //should be function of depth

        double vol=state_vars[iStor]/MM_PER_METER;
        rates[q] = (mass-C_sat*vol) * (1.0 - exp(-exch_coeff*Options.timestep))/Options.timestep;
      }
    }*/
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

