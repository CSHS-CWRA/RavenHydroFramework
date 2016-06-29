/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright © 2008-2014 the Raven Development Team
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
  constit_ind=pTransModel->GetConstituentIndex(constit_name);

  int nAdvConnections=pTransModel->GetNumAdvConnections();

  CHydroProcessABC::DynamicSpecifyConnections(3*nAdvConnections+1);

  for (int q=0;q<nAdvConnections;q++) //for regular advection
  {
    iFrom[q]=pTransModel->GetFromIndex(constit_ind,q);
    iTo  [q]=pTransModel->GetToIndex  (constit_ind,q);
  }
  for (int q=0;q<nAdvConnections;q++)//for Dirichlet source correction (from)
  {
    iFrom[nAdvConnections+q]=pModel->GetStateVarIndex(CONSTITUENT_SRC,constit_ind);
    iTo  [nAdvConnections+q]=iFrom[q];
  }
  for (int q=0;q<nAdvConnections;q++)//for Dirichlet source correction (to)
  {
    iFrom[2*nAdvConnections+q]=pModel->GetStateVarIndex(CONSTITUENT_SRC,constit_ind);
    iTo  [2*nAdvConnections+q]=iTo[q];
  }
  //advection into surface water
  int iSW=pModel->GetStateVarIndex(SURFACE_WATER,0);
  int   m=pTransModel->GetLayerIndex(constit_ind,iSW);
  iFrom[3*nAdvConnections]=pModel->GetStateVarIndex(CONSTITUENT,m);
  iTo  [3*nAdvConnections]=pModel->GetStateVarIndex(CONSTITUENT_SW,0);
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvAdvection::~CmvAdvection(){}

//////////////////////////////////////////////////////////////////
/// \brief Initializes Advection object
//
void   CmvAdvection::Initialize(){}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for baseflow algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by advection algorithm (size of aP[] and aPC[])
//
void CmvAdvection::GetParticipatingParamList(string  *aP, class_type *aPC, int &nP) const
{
  nP=0;
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
void   CmvAdvection::GetRatesOfChange(const double			*state_vars, 
																			const CHydroUnit	*pHRU, 
																			const optStruct	  &Options,
																			const time_struct &tt,
                                      double            *rates) const
{
  int    q,iFromWater,iToWater,js;
  double mass,vol,Cs;
  double Rf;                 //retardation factor
  double sv[MAX_STATE_VARS]; //state variable history

  double tstep=Options.timestep;
  int    nAdvConnections=pTransModel->GetNumAdvConnections();
  int    k=pHRU->GetGlobalIndex();
  double *Q=new double [nAdvConnections]; // \todo [optimize]: preallocate in initialize, save as member

  ExitGracefullyIf(Options.sol_method!=ORDERED_SERIES,
    "CmvAdvection: Advection only works with ordered series solution approach",BAD_DATA);// \todo [re-org] Should go in initialize

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
    
    sv[iFromWater]+=Q[q]*tstep;//reverse determination of state variable history (i.e., rewinding flow calculations)
    sv[iToWater]  -=Q[q]*tstep;
  }
  //if (false){//
  /*if (k==0) {
    cout <<"start ";
    for (i=1;i<pModel->GetNumStateVars();i++){
      if (CStateVariable::IsWaterStorage(pModel->GetStateVarType(i))){cout <<sv[i]<<" ";}
    }
    cout<<endl<<"-"<<endl;
    cout <<"end   ";
    for (i=1;i<pModel->GetNumStateVars();i++){
      if (CStateVariable::IsWaterStorage(pModel->GetStateVarType(i))){cout <<state_vars[i]<<" ";}
    }
    cout<<endl;
  }*/
  //Calculate advective mass fluxes
  //-------------------------------------------------------------------------
  for (q=0;q<nAdvConnections;q++)
  {
    rates[q]=0.0;
    int iFromWater=pTransModel->GetFromWaterIndex(q);
    int iToWater  =pTransModel->GetToWaterIndex  (q);

    Rf=1.0;
    //Rf=pTransModel->GetRetardationFactor(constit_ind,pHRU,iFromWater,iToWater);
    
    //Advection rate calculation rates[q]=dm/dt=Q*C
    mass=0;vol=1;
    if      (Q[q]>0)
    {
      mass=sv[iFrom[q]];   //[mg/m2]
      vol =sv[iFromWater]; //[mm]
    }
    else if (Q[q]<0)
    {
      mass=sv[iTo[q]];
      vol =sv[iToWater];
    }
    if (vol>1e-6){//note: otherwise Q should generally be constrained to be <vol/tstep & 0.0<rates[q]<(m/tstep/Rf)
      rates[q]=Q[q]*mass/vol/Rf; //[mg/m2/d] 
      if (mass<-1e-9){ ExitGracefully("CmvAdvection - negative mass",RUNTIME_ERR); }
      if (fabs(rates[q])>mass/tstep){rates[q]=(Q[q]/fabs(Q[q]))*mass/tstep;}//emptying out compartment
      /*if ((mass / vol)/LITER_PER_M3*MM_PER_METER > 1.0){
        cout << " conc > 1 !"<<pHRU->GetID()<<" "<<vol<<" "<<Q[q]<<" C:"<<(mass / vol)/LITER_PER_M3*MM_PER_METER<<endl;
      }*/
    }
    
    //special consideration - atmospheric precip can have negative storage but still specified concentration
    if ((pModel->GetStateVarType(iFromWater)==ATMOS_PRECIP) &&
         (pTransModel->IsDirichlet(iFromWater,constit_ind,k,tt,Cs))){ 
      Cs*=LITER_PER_M3/MM_PER_METER; //[mg/L]->[mg/mm-m2] 
      rates[q]=Q[q]*Cs;
    }

    //double Ctmp1 = sv[iTo  [q]  ]/sv[iToWater  ]/LITER_PER_M3*MM_PER_METER ;
    //double mtmp1=sv[iTo  [q]  ]/LITER_PER_M3*MM_PER_METER ; double vtmp1=sv[iToWater  ];

    //update local mass and volume history
    sv[iFromWater]-=Q[q]*tstep;
    sv[iToWater  ]+=Q[q]*tstep;
    sv[iFrom[q]  ]-=rates[q]*tstep;
    sv[iTo  [q]  ]+=rates[q]*tstep;

    //double Ctmp = sv[iTo  [q]  ]/sv[iToWater  ]/LITER_PER_M3*MM_PER_METER ;
    //double mtmp2=sv[iTo  [q]  ]/LITER_PER_M3*MM_PER_METER ; double vtmp2=sv[iToWater  ];
    //if ((Ctmp>1.0001) && (sv[iToWater  ]!=0)){cout<<"violation: "<<Ctmp1<<" ("<<mtmp1<<"/"<<vtmp1<<") "<<Ctmp<<" ("<<mtmp2<<"/"<<vtmp2<<") "<<rates[q]*tstep<<" "<< Q[q]*tstep<<" "<<CStateVariable::SVTypeToString(pModel->GetStateVarType(iFromWater),-1)<<" "<<CStateVariable::SVTypeToString(pModel->GetStateVarType(iToWater),-1)<<endl;}

    //Correct for Dirichlet conditions
    //----------------------------------------------------------------------
    if (pTransModel->IsDirichlet(iFromWater,constit_ind,k,tt,Cs))
    {
      Cs*=LITER_PER_M3/MM_PER_METER; //[mg/L]->[mg/mm-m2] 
      mass=sv[iFrom[q]];
      vol =sv[iFromWater];
      rates[nAdvConnections+q]+=(Cs*vol-mass)/Options.timestep;
      sv[iFrom[q]]=Cs*vol;
    }

    if (pTransModel->IsDirichlet(iToWater,constit_ind,k,tt,Cs))
    {
      Cs*=LITER_PER_M3/MM_PER_METER; //[mg/L]->[mg/mm-m2] 
      mass=sv[iTo[q]];
      vol =sv[iToWater];
      rates[2*nAdvConnections+q]+=(Cs*vol-mass)/Options.timestep;
      sv[iTo[q]]=Cs*vol;
    }

    // \todo [funct]: handle dumping into surface water
    //if ROUTE_NONE, dump surface water compartment to CONSTIT_SW
    /*if (Options.catchment_routing==ROUTE_NONE){
      mass=sv[iFrom[3*nAdvConnections]];
      rates[3*nAdvConnections]=mass/Options.timestep;
    }*/
  }
  delete [] Q;
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
void   CmvAdvection::ApplyConstraints(const double		 *state_vars, 
						                            const CHydroUnit *pHRU, 
						                            const optStruct	 &Options,
						                            const time_struct &tt,
                                              double     *rates) const
{
  
}

