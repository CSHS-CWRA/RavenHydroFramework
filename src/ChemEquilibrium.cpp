/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2022 the Raven Development Team
  ------------------------------------------------------------------
  Chemical Equilibrium between two constituents
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "ChemEquilibrium.h"
#include "Model.h"
#include "StateVariables.h"

  //////////////////////////////////////////////////////////////////
  /// \brief Implentation of the chemical equilibrium constructor
  /// \param constit_name1 [in] name of constituent 1 
  /// \param constit_name2 [in] name of constituent 2
  /// \param etyp [in] decay process type
  /// \param proc_ind [in] transport process index
  /// \param iWatStor [in] global SV index of water storage (or DOESNT_EXIST, if process should be applied in all water compartments)
  /// \param pTransportModel [in] transport Model object
  //
CmvChemEquil::CmvChemEquil(string           constit_name1,
                           string           constit_name2,
                           chem_equil_type  etyp,
                           int              proc_ind,
                           int              iWatStor,
                           CTransportModel* pTransportModel)
  :CHydroProcessABC(DECAY)
{
  _eq_type     =etyp;
  _pTransModel =pTransportModel;
  _constit_ind1=_pTransModel->GetConstituentIndex(constit_name1);
  _constit_ind2=_pTransModel->GetConstituentIndex(constit_name2);
  _process_ind =proc_ind;
  _iWaterStore =iWatStor;

  ExitGracefullyIf(_constit_ind1==DOESNT_EXIST,
    "CmvChemEquil constructor: invalid constituent name in :Equilibrium command",BAD_DATA_WARN);
  ExitGracefullyIf(_constit_ind2==DOESNT_EXIST,
    "CmvChemEquil constructor: invalid constituent name in :Equilibrium command",BAD_DATA_WARN);

  int nWaterCompartments = _pTransModel->GetNumWaterCompartments();
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
CmvChemEquil::~CmvChemEquil() {}

//////////////////////////////////////////////////////////////////
/// \brief Initializes equilibrium object
//
void   CmvChemEquil::Initialize() {}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for baseflow algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by advection algorithm (size of aP[] and aPC[])
//
void CmvChemEquil::GetParticipatingParamList(string* aP,class_type* aPC,int& nP) const
{
  nP=0;
  if(_eq_type==EQUIL_FIXED_RATIO)
  {
    aP[0]="EQFIXED_RATIO"; aPC[0]=CLASS_TRANSPORT; nP++;
  }
  else if(_eq_type==EQUIL_LINEAR_SORPTION)
  {
    aP[0]="SORPT_COEFF";   aPC[0]=CLASS_TRANSPORT; nP++;
    aP[1]="POROSITY";      aPC[1]=CLASS_SOIL; nP++;
    aP[2]="BULK_DENSITY";  aPC[2]=CLASS_SOIL; nP++;
  }
  else if(_eq_type==EQUIL_LINEAR)
  {
    aP[0]="EQFIXED_RATIO";   aPC[0]=CLASS_TRANSPORT; nP++;
    aP[1]="EQUIL_COEFF";     aPC[1]=CLASS_TRANSPORT; nP++;
    aP[2]="STOICHIO_RATIO";  aPC[2]=CLASS_TRANSPORT; nP++;
  }
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
void   CmvChemEquil::GetRatesOfChange(const double     * state_vars,
                                      const CHydroUnit * pHRU,
                                      const optStruct  & Options,
                                      const time_struct& tt,
                                            double     * rates) const
{
  int     k=pHRU->GetGlobalIndex();
  double  junk,mass1,mass2,vol,cc;
  int     iStor;


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
    mass2=state_vars[iFrom[q]]; //mg/m2
    vol  =state_vars[iStor];    //mm

    //-------------------------------------------------------------------------------
    if (_eq_type==EQUIL_FIXED_RATIO)
    {
      cc = _pTransModel->GetGeochemParam(PAR_EQFIXED_RATIO,_constit_ind1,_constit_ind2,ii,_process_ind,pHRU);
      if(cc==NOT_SPECIFIED) { continue; }

      // C_1 = cc*C_2 (cc positive)
      // to enforce, 
      // dC_1=-(C_1^0-cc*C_2^0)/(1+cc) (equal to zero if C_1^0=cc*C_2^0)
      // dC_2=-dC_1

      rates[q] = -(mass1-cc*mass2)/(1+cc)/Options.timestep;
      rates[q+shift]=-rates[q]; //does stoichiometry not play a role here???
    }
    //-------------------------------------------------------------------------------
    else if(_eq_type==EQUIL_LINEAR_SORPTION)
    {
      if(pModel->GetStateVarType(iStor)==SOIL) 
      {
        //Cs = Kd*C_w
        //assumes constit 1 is aqueous, constit 2 is sorbed
        
        double Kd = _pTransModel->GetGeochemParam(PAR_SORPT_COEFF,_constit_ind1,_constit_ind2,ii,_process_ind,pHRU);
        if(Kd==NOT_SPECIFIED) { continue; }
        Kd/=LITER_PER_M3; //[m3/kg]

        int    m    =pModel->GetStateVarLayer(iStor);

        double thick=pHRU->GetSoilThickness(m); //[mm]
        double poro =pHRU->GetSoilProps    (m)->porosity;
        double rho_b=pHRU->GetSoilProps    (m)->bulk_density; //[kg/m3]
        double Vsoil=thick*(1-poro); //[mm]

        //m_s = cc*m_w
        cc=Kd*rho_b*Vsoil/vol;

        rates[q] = (mass2-cc*mass1)/(1+cc)/Options.timestep;
        rates[q+shift]=-rates[q]; 
      }
    }
    //-------------------------------------------------------------------------------
    else if(_eq_type==EQUIL_LINEAR)
    {
      // dC_1/dt = -k  *(C_1-cc*C_2)       (cc positive)
      // dC_2/dt = +k*s*(C_1-cc*C_2)
      double eq_coeff;
      double stoich_ratio; //stoichiometric ratio, i.e., A->s*B [mg]

      cc           = _pTransModel->GetGeochemParam(PAR_EQFIXED_RATIO ,_constit_ind1,_constit_ind2,ii,_process_ind,pHRU);
      eq_coeff     = _pTransModel->GetGeochemParam(PAR_EQUIL_COEFF   ,_constit_ind1,_constit_ind2,ii,_process_ind,pHRU);
      stoich_ratio = _pTransModel->GetGeochemParam(PAR_STOICHIO_RATIO,_constit_ind1,_constit_ind2,ii,_process_ind,pHRU);
      if(cc          ==NOT_SPECIFIED) { continue; }
      if(eq_coeff    ==NOT_SPECIFIED) { continue; }
      if(stoich_ratio==NOT_SPECIFIED) { continue; }

      rates[q] =  -(mass1-cc*mass2) * (1.0 - exp(-eq_coeff*Options.timestep))/Options.timestep;
      rates[q+shift]=-stoich_ratio*rates[q];
    }

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
void   CmvChemEquil::ApplyConstraints(const double     * state_vars,
                                      const CHydroUnit * pHRU,
                                      const optStruct  & Options,
                                      const time_struct& tt,
                                            double     * rates) const
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
    double stoich_ratio    = _pTransModel->GetGeochemParam(PAR_STOICHIO_RATIO,_constit_ind1,_constit_ind2,ii,_process_ind,pHRU);

    iConstit1=_pTransModel->GetStorIndex(_constit_ind1,ii);
    rates[q] = min(rates[q],state_vars[iConstit1]/Options.timestep);//cannot remove more mass than is there
    rates[q+shift]=-stoich_ratio*rates[q]; //handles proper accounting 
    q++;
  }
}

