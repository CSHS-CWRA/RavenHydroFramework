/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2020 the Raven Development Team
  ------------------------------------------------------------------
  Advection of soluble contaminant/tracer/nutrient
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "Advection.h"
#include "Model.h"
#include "StateVariables.h"

//////////////////////////////////////////////////////////////////
/// \brief Implentation of the Advection constructor
/// \param constituent [in] name of contaminant beign tracked
/// \param pFlow [in] flow process which drives advection (this acts as a wrapper for said process)
/// \param pModel [in] Model object
//
CmvAdvection::CmvAdvection(string constit_name,
                           CTransportModel *pTransportModel)
  :CHydroProcessABC(ADVECTION)
{
  pTransModel=pTransportModel;
  _constit_ind=pTransModel->GetConstituentIndex(constit_name);

  int nAdvConnections=pTransModel->GetNumAdvConnections();

  CHydroProcessABC::DynamicSpecifyConnections(3*nAdvConnections+1);//+2 once GW added

  for (int q=0;q<nAdvConnections;q++) //for regular advection
  {
    iFrom[q]=pTransModel->GetFromIndex(_constit_ind,q);
    iTo  [q]=pTransModel->GetToIndex  (_constit_ind,q);
  }
  for (int q=0;q<nAdvConnections;q++)//for Dirichlet/Neumann source correction (from)
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
  double Rf;                 //retardation factor
  double sv[MAX_STATE_VARS]; //state variable history
  
  double tstep=Options.timestep;
  int    nAdvConnections=pTransModel->GetNumAdvConnections();
  bool   isEnthalpy=(pTransModel->GetConstituent(_constit_ind)->type==ENTHALPY);
  int    k=pHRU->GetGlobalIndex();

  static double    *Q=NULL; 
  if(Q==NULL){
    Q=new double [nAdvConnections]; // only done once at start of simulation for speed 
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
  for (q=0;q<nAdvConnections;q++)
  {
    rates[q]=0.0;
    iFromWater=pTransModel->GetFromWaterIndex(q);
    iToWater  =pTransModel->GetToWaterIndex  (q);

    Rf=1.0;
    //Rf=pTransModel->GetRetardationFactor(constit_ind,pHRU,iFromWater,iToWater);

    //Advection rate calculation rates[q]=dm/dt=Q*C or dE/dt=Q*h
    mass=0;vol=1;
    if      (Q[q]>0)
    {
      mass=sv[iFrom[q]];   //[mg/m2] or [MJ/m2]
      vol =sv[iFromWater]; //[mm]
    }
    else if (Q[q]<0)
    {
      mass=sv[iTo[q]];
      vol =sv[iToWater];
    }
    if (vol>1e-6)//note: otherwise Q should generally be constrained to be <vol/tstep & 0.0<rates[q]<(m/tstep/Rf)
    {
      rates[q]=Q[q]*mass/vol/Rf; //[mg/m2/d] or [MJ/m2/d]

      if(!isEnthalpy) { //negative enthalpy/temperature allowed (but shouldnt happen)
        if(mass<-1e-9) { ExitGracefully("CmvAdvection - negative mass",RUNTIME_ERR); }
      }
      if(fabs(rates[q])>mass/tstep) { rates[q]=(Q[q]/fabs(Q[q]))*mass/tstep; }//emptying out compartment
    }

    //special consideration - atmospheric precip can have negative storage but still specified concentration or temperature
    //----------------------------------------------------------------------
    if ((pModel->GetStateVarType(iFromWater)==ATMOS_PRECIP) &&
        (pTransModel->IsDirichlet(iFromWater,_constit_ind,k,tt,Cs)))
    {
      if(isEnthalpy) {
        if (Cs==DIRICHLET_AIR_TEMP) { //Special precip temperature condition
          double Tair    =pHRU->GetForcingFunctions()->temp_ave;
          double snowfrac=pHRU->GetForcingFunctions()->snow_frac;
          
          Cs=ConvertTemperatureToVolumetricEnthalpy(Tair,snowfrac)/MM_PER_METER; //[C]->[MJ/m3]->[MJ/mm-m2]
          Cs=max(Cs,0.0); //precip @ 0 degrees
        }
        else {
          double pctFroz=0.0; //assumes liquid water for flows
          Cs=ConvertTemperatureToVolumetricEnthalpy(Cs,pctFroz)/MM_PER_METER;    //[C]->[MJ/m3]->[MJ/mm-m2]
        }
      }
      else {
        Cs*=LITER_PER_M3/MM_PER_METER; //[mg/L]->[mg/mm-m2]
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
    if (pTransModel->IsDirichlet(iFromWater,_constit_ind,k,tt,Cs))
    {
      if(isEnthalpy) {
        if(Cs==DIRICHLET_AIR_TEMP) { //Special precip temperature condition
          Cs=0; //don't explicitly try to track enthalpy content of atmosphere
        }
        else {
          double pctFroz=0.0; //assumes liquid water for flows (reasonable)
          Cs=ConvertTemperatureToVolumetricEnthalpy(Cs,pctFroz)/MM_PER_METER;    //[C]->[MJ/m3]->[MJ/mm-m2]
        }
      }
      else {
        Cs*=LITER_PER_M3/MM_PER_METER; //[mg/L]->[mg/mm-m2]
      }
      mass          =sv[iFrom[q]]; 
      dirichlet_mass=Cs*sv[iFromWater]; //Cs*V 
      rates[nAdvConnections+q]+=(dirichlet_mass-mass)/Options.timestep; //From CONSTITUENT_SRC
      sv[iFrom[q]]=dirichlet_mass;      //override previous mass/enthalpy
    }

    if (pTransModel->IsDirichlet(iToWater,_constit_ind,k,tt,Cs))
    {
      if(isEnthalpy) {
        double pctFroz=0.0; //assumes liquid water for flows (reasonable)
        Cs=ConvertTemperatureToVolumetricEnthalpy(Cs,pctFroz)/MM_PER_METER;    //[C]->[MJ/m3]->[MJ/mm-m2]
      }
      else {
        Cs*=LITER_PER_M3/MM_PER_METER; //[mg/L]->[mg/mm-m2]
      }
      mass          =sv[iTo[q]];
      dirichlet_mass=Cs*sv[iToWater]; //C*V
      rates[2*nAdvConnections+q]+=(dirichlet_mass-mass)/Options.timestep; //From CONSTITUENT_SRC
      sv[iTo[q]]=dirichlet_mass;      //override previous mass/enthalpy
    }

    // \todo [funct]: handle dumping into surface water
    //if ROUTE_NONE, dump surface water compartment to CONSTIT_SW
    /*if (Options.catchment_routing==ROUTE_NONE){
      mass=sv[iFrom[3*nAdvConnections]];
      rates[3*nAdvConnections]=mass/Options.timestep;
      }*/
  } //ends "for (q=0;q<nAdvConnections;q++).."

  //Handle Neumann influx conditions, if present
  //-------------------------------------------------------
  /*for (q = 0; q < nAdvConnections; q++)
  {
    rates[q] = 0.0;
    int iFromWater = pTransModel->GetFromWaterIndex(q);
    rates[nAdvConnections + q] += pTransModel->GetSpecifiedMassFlux(iToWater, constit_ind, k, tt); //[mg/m2/d]
  }*/

  //delete static memory during last timestep
  if(tt.model_time>=Options.duration-Options.timestep/2)
  {
    delete [] Q;
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

