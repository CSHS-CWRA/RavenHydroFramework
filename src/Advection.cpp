/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2022 the Raven Development Team
  ------------------------------------------------------------------
  Advection of soluble contaminant/tracer/nutrient
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "Advection.h"
#include "Model.h"
#include "StateVariables.h"
#include "EnergyTransport.h"

//////////////////////////////////////////////////////////////////
/// \brief Implentation of the Advection constructor
/// \param constit_name [in] name of constituent being tracked
/// \param pModel [in] Model object
//
CmvAdvection::CmvAdvection(string constit_name,
                           CTransportModel *pTransportModel)
  :CHydroProcessABC(ADVECTION)
{
  pTransModel=pTransportModel;
  _constit_ind=pTransModel->GetConstituentIndex(constit_name);

  int nAdvConnections    =pTransModel->GetNumAdvConnections();
  CHydroProcessABC::DynamicSpecifyConnections(3*nAdvConnections+1);//+2 once GW added

  for (int q=0;q<nAdvConnections;q++) //for regular advection
  {
    iFrom[q]=pTransModel->GetFromIndex(_constit_ind,q);
    iTo  [q]=pTransModel->GetToIndex  (_constit_ind,q);
  }
  for (int q=0;q<nAdvConnections;q++)//for Dirichlet source correction (from)
  {
    iFrom[nAdvConnections+q]=pModel->GetStateVarIndex(CONSTITUENT_SRC,_constit_ind);
    iTo  [nAdvConnections+q]=iFrom[q];
  }
  for (int q=0;q<nAdvConnections;q++)//for Dirichlet source correction (to)
  {
    iFrom[2*nAdvConnections+q]=pModel->GetStateVarIndex(CONSTITUENT_SRC,_constit_ind);
    iTo  [2*nAdvConnections+q]=iTo[q];
  }
  //advection into surface water
  int iSW=pModel->GetStateVarIndex(SURFACE_WATER);
  int   m=pTransModel->GetLayerIndex(_constit_ind,iSW);
  iFrom[3*nAdvConnections]=pModel->GetStateVarIndex(CONSTITUENT,m);
  iTo  [3*nAdvConnections]=pModel->GetStateVarIndex(CONSTITUENT_SW,0);

  //advection into groundwater
  //int iGW=pModel->GetStateVarIndex(GROUNDWATER);
  //int   m=pTransModel->GetLayerIndex(constit_ind,iGW);
  //iFrom[3*nAdvConnections+1]=pModel->GetStateVarIndex(CONSTITUENT,m);
  //iTo  [3*nAdvConnections+1]=pModel->GetStateVarIndex(CONSTITUENT_GW,0);

}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvAdvection::~CmvAdvection(){}

//////////////////////////////////////////////////////////////////
/// \brief Initializes Advection object
//
void   CmvAdvection::Initialize()
{
  ExitGracefullyIf(pModel->GetOptStruct()->sol_method!=ORDERED_SERIES,
    "CmvAdvection:Initalize: Advection only works with ordered series solution approach",BAD_DATA);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for baseflow algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by advection algorithm (size of aP[] and aPC[])
//
void CmvAdvection::GetParticipatingParamList(string  *aP, class_type *aPC, int &nP) const
{
  nP=0; //All parameters are externally determined
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
///
/// \param pModel [in] Model
/// \param *aSV [out] Array of state variable types needed by advection algorithm
/// \param *aLev [out] Array of layer of multilayer state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by CmvAdvection algorithm (size of aSV[] and aLev[] arrays)
//
/*void CmvAdvection::GetParticipatingStateVarList(sv_type *aSV, int *aLev, int &nSV)
  {
  //only returns added state vars, not those associated with flow model

  //nSV=pTransModel->GetNumWaterCompartments();
  //for (int i=0;i<nSV;i++){
  //  aSV[i  ]=CONSTITUENT; aLev[i  ]=i+constit_ind*nSV;
  //}
  }*/

//////////////////////////////////////////////////////////////////
/// \brief For modeling changes contaminant mass/concentration linked to a flow process
///
/// \param *state_vars [in] array of current state variables
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] rates of change due to both associated flow process and advective transport
//
void   CmvAdvection::GetRatesOfChange(const double      *state_vars,
                                      const CHydroUnit  *pHRU,
                                      const optStruct   &Options,
                                      const time_struct &tt,
                                      double            *rates) const
{
  int    q,iFromWater,iToWater,js;
  double mass,vol,Cs;
  double corr;               // advective correction (e.g.,1/retardation factor)
  
  double tstep=Options.timestep;
  int    nAdvConnections     =pTransModel->GetNumAdvConnections();
  CConstituentModel *pConstit=pTransModel->GetConstituentModel2(_constit_ind);

  bool   isEnthalpy=(pConstit->GetType()==ENTHALPY);
  bool   isIsotope =(pConstit->GetType()==ISOTOPE);
  int    k=pHRU->GetGlobalIndex();

  static double *Q =NULL; 
  static double *sv=NULL;
  if(Q==NULL){
    Q =new double [nAdvConnections]; // only done once at start of simulation for speed 
    sv=new double [pModel->GetNumStateVars()];
  }

  // copy all state variables into array
  memcpy(sv/*dest*/,state_vars/*src*/,sizeof(double)*(pModel->GetNumStateVars()));

  //get water fluxes, regenerate system state at start of timestep
  //-------------------------------------------------------------------------
  for (q=0;q<nAdvConnections;q++)
  {
    iFromWater=pTransModel->GetFromWaterIndex(q);
    iToWater  =pTransModel->GetToWaterIndex  (q);
    js        =pTransModel->GetJsIndex       (q);

    Q[q]=pModel->GetFlux(k,js,Options); //[mm/d]

    sv[iFromWater]+=Q[q]*tstep; //reverse determination of state variable history (i.e., rewinding flow calculations)
    sv[iToWater]  -=Q[q]*tstep;
  }

  //Calculate advective mass fluxes
  //-------------------------------------------------------------------------
  int iF,iT;
  for (q=0;q<nAdvConnections;q++)
  {
    rates[q]=0.0;
    iFromWater=pTransModel->GetFromWaterIndex(q);
    iToWater  =pTransModel->GetToWaterIndex  (q);
    //Advection rate calculation rates[q]=dm/dt=Q*C or dE/dt=Q*h
    mass=0;vol=1;
    iF=iFromWater;
    iT=iToWater;
    if      (Q[q]>0)
    {
      mass=sv[iFrom[q]];   //[mg/m2] or [MJ/m2]
      vol =sv[iFromWater]; //[mm]
    }
    else if (Q[q]<0)
    {
      mass=sv[iTo[q]];
      vol =sv[iToWater];
      iF=iToWater;iT=iFromWater;
    }
    if (vol>1e-6)//note: otherwise Q should generally be constrained to be <vol/tstep & 0.0<rates[q]<(m/tstep*corr)
    {
      corr=pTransModel->GetAdvectionCorrection(_constit_ind,pHRU,iF,iT,mass/vol*MM_PER_METER/LITER_PER_M3);//mg/m2/mm->mg/L

      rates[q]=Q[q]*mass/vol*corr; //[mg/m2/d] or [MJ/m2/d]
      if(fabs(rates[q])>mass/tstep) { rates[q]=(Q[q]/fabs(Q[q]))*mass/tstep; }//emptying out compartment

      if ((!isEnthalpy) && (!isIsotope)){ //negative enthalpy/temperature/isotopic composition allowed
        if(mass<-1e-9) { ExitGracefully("CmvAdvection - negative mass",RUNTIME_ERR); }
      }
    }
    
    //special consideration - atmospheric precip can have negative storage but still specified concentration or temperature
    //----------------------------------------------------------------------
    if((pModel->GetStateVarType(iF)==ATMOS_PRECIP) &&
      (pConstit->IsDirichlet(iF,k,tt,Cs,1.0-pHRU->GetForcingFunctions()->snow_frac)))
    {
      if(!isEnthalpy) {
        Cs=pConstit->ConvertConcentration(Cs);
      }
      else {
        Cs=pTransModel->GetEnthalpyModel()->GetDirichletEnthalpy(pHRU,Cs);
      }

      rates[q]=Q[q]*Cs; //[mg/m2/d] or [MJ/m2/d]
    }

    //update local mass and volume history
    sv[iFromWater]-=Q[q]*tstep;
    sv[iToWater  ]+=Q[q]*tstep;
    sv[iFrom[q]  ]-=rates[q]*tstep;
    sv[iTo  [q]  ]+=rates[q]*tstep;

    //Correct for Dirichlet conditions - update/override source mass and update MB tracking [CONSTITUENT_SRC]
    // two cases: contributor (iFromWater) or recipient (iToWater) water compartment is Dirichlet condition
    //----------------------------------------------------------------------
    double dirichlet_mass;
    if(pConstit->IsDirichlet(iFromWater,k,tt,Cs))
    {
      if(!isEnthalpy) {
        Cs=pConstit->ConvertConcentration(Cs);//[mg/L]->[mg/mm-m2]
      }
      else{
        Cs=pTransModel->GetEnthalpyModel()->GetDirichletEnthalpy(pHRU,Cs);
        if(pModel->GetStateVarType(iFromWater)==ATMOS_PRECIP) { Cs=0.0; }//don't explicitly try to track enthalpy content of atmospheric precip
      }
      mass          =sv[iFrom[q]]; 
      dirichlet_mass=Cs*sv[iFromWater]; //Cs*V 
      rates[nAdvConnections+q]+=(dirichlet_mass-mass)/Options.timestep; //From CONSTITUENT_SRC
      sv[iFrom[q]]=dirichlet_mass;      //override previous mass/enthalpy
    }

    if(pConstit->IsDirichlet(iToWater,k,tt,Cs))
    {
      if(!isEnthalpy) {
        Cs=pConstit->ConvertConcentration(Cs);
      }
      else{
        Cs=pTransModel->GetEnthalpyModel()->GetDirichletEnthalpy(pHRU,Cs);
      }
      mass          =sv[iTo[q]];
      dirichlet_mass=Cs*sv[iToWater]; //C*V
      rates[2*nAdvConnections+q]+=(dirichlet_mass-mass)/Options.timestep; //From CONSTITUENT_SRC
      sv[iTo[q]]=dirichlet_mass;      //override previous mass/enthalpy
    }

  } //ends "for (q=0;q<nAdvConnections;q++).."

  //delete static memory during last timestep
  if(tt.model_time>=Options.duration-Options.timestep/2)
  {
    delete [] Q;
    delete [] sv;
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
void   CmvAdvection::ApplyConstraints(const double               *state_vars,
                                      const CHydroUnit *pHRU,
                                      const optStruct      &Options,
                                      const time_struct &tt,
                                      double     *rates) const
{

}

/*if((vol>1e-6) && (mass/vol>1.0)) {
string tmpname;
if     (Q[q]>0){tmpname=CStateVariable::GetStateVarLongName(pModel->GetStateVarType(iFromWater),0);}
else if(Q[q]<0){tmpname=CStateVariable::GetStateVarLongName(pModel->GetStateVarType(iToWater),0);}
cout<<"BADCONC: "<<tmpname<<" "<<mass<<" "<<vol<<" "<<mass/vol<<" "<<tt.date_string<<endl;
}*/
/*if(vol>1e-6)//note: otherwise Q should generally be constrained to be <vol/tstep & 0.0<rates[q]<(m/tstep/Rf)
{
  rates[q]=Q[q]*mass/vol/Rf; //[mg/m2/d] or [MJ/m2/d]

  if(!isEnthalpy) { //negative enthalpy/temperature allowed (but shouldnt happen)
    if(mass<-1e-9) { ExitGracefully("CmvAdvection - negative mass",RUNTIME_ERR); }
  }
  if(fabs(rates[q])>mass/tstep) { rates[q]=(Q[q]/fabs(Q[q]))*mass/tstep; }//emptying out compartment
}
else if(mass>0) {
  string tmpname,tmp2;
  if(Q[q]>0)
  {
    tmpname=CStateVariable::GetStateVarLongName(pModel->GetStateVarType(iFromWater),0);
    tmp2=CStateVariable::GetStateVarLongName(pModel->GetStateVarType(iToWater),0);
  }
  else if(Q[q]<0)
  {
    tmpname=CStateVariable::GetStateVarLongName(pModel->GetStateVarType(iToWater),0);
    tmp2=CStateVariable::GetStateVarLongName(pModel->GetStateVarType(iFromWater),0);
  }
  if(vol<-0.01) {
    cout<<"TINYVOL: "<<Q[q]<<" "<<vol<<" "<<tmpname<<" to "<<tmp2<<endl;
  }
  //rates[q]=mass/tstep/Rf;
}*/