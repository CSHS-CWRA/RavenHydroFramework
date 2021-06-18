/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2021 the Raven Development Team
------------------------------------------------------------------
Advection of soluble contaminant/tracer/nutrient
----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "LatAdvection.h"
#include "Model.h"
#include "StateVariables.h"

//////////////////////////////////////////////////////////////////
/// \brief Implentation of the LatAdvection constructor
/// \param constituent [in] name of contaminant beign tracked
/// \param pFlow [in] flow process which drives advection (this acts as a wrapper for said process)
/// \param pModel [in] Model object
//
CmvLatAdvection::CmvLatAdvection(string constit_name,
                           CTransportModel *pTransportModel)
  :CLateralExchangeProcessABC(LAT_ADVECTION)
{
  pTransModel=pTransportModel;
  _constit_ind=pTransModel->GetConstituentIndex(constit_name);
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvLatAdvection::~CmvLatAdvection() {}

//////////////////////////////////////////////////////////////////
/// \brief Initializes Advection object
//
void   CmvLatAdvection::Initialize() 
{
  int nAdvConnections=pTransModel->GetNumLatAdvConnections();

  if(nAdvConnections==0){return;}//first time through

  DynamicSpecifyLatConnections(3*nAdvConnections+1);
  for(int q=0;q<nAdvConnections;q++) //for regular advection
  {
    _iFromLat[q]=pTransModel->GetLatFromIndex(_constit_ind,q);
    _iToLat  [q]=pTransModel->GetLatToIndex  (_constit_ind,q);
    _kFrom   [q]=pTransModel->GetLatFromHRU(q);
    _kTo     [q]=pTransModel->GetLatToHRU(q);
  }
  for(int q=0;q<nAdvConnections;q++)//for Dirichlet/Neumann source correction (from)
  {
    _iFromLat[nAdvConnections+q]=pModel->GetStateVarIndex(CONSTITUENT_SRC,_constit_ind);
    _iToLat  [nAdvConnections+q]=_iFromLat[q];
    _kFrom   [nAdvConnections+q]=pTransModel->GetLatFromHRU(q);
    _kTo     [nAdvConnections+q]=pTransModel->GetLatFromHRU(q);
  }
  for(int q=0;q<nAdvConnections;q++)//for Dirichlet source correction (to)
  {
    _iFromLat[2*nAdvConnections+q]=pModel->GetStateVarIndex(CONSTITUENT_SRC,_constit_ind);
    _iToLat  [2*nAdvConnections+q]=_iToLat[q];
    _kFrom   [2*nAdvConnections+q]=pTransModel->GetLatToHRU(q);
    _kTo     [2*nAdvConnections+q]=pTransModel->GetLatToHRU(q);
  }
  //advection into surface water
  int iSW=pModel->GetStateVarIndex(SURFACE_WATER,0);
  int   m=pTransModel->GetLayerIndex(_constit_ind,iSW);
  _iFromLat[3*nAdvConnections]=pModel->GetStateVarIndex(CONSTITUENT,m);
  _iToLat  [3*nAdvConnections]=pModel->GetStateVarIndex(CONSTITUENT_SW,0);
  _kFrom   [3*nAdvConnections]=pTransModel->GetLatToHRU(0); // TMP DEBUG ???
  _kTo     [3*nAdvConnections]=pTransModel->GetLatToHRU(0);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for baseflow algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by advection algorithm (size of aP[] and aPC[])
//
void CmvLatAdvection::GetParticipatingParamList(string  *aP,class_type *aPC,int &nP) const
{
  nP=0;
}

//////////////////////////////////////////////////////////////////
/// \brief For modeling changes contaminant mass/concentration linked to a flow process
///
/// \param **state_vars [in] 2D array of current state variables [nHRUs][nSVs]
/// \param **pHRUs [in] array of pointers to HRUs
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *exchange_rates [out] Rate of loss from "from" compartment [mg/day] 
//
void   CmvLatAdvection::GetLateralExchange(const double * const *state_vars, //array of all SVs for all HRUs, [k][i]
                                            const CHydroUnit * const *pHRUs,
                                            const optStruct   &Options,
                                            const time_struct &tt,
                                            double      *exchange_rates) const
{

  int    q,iFromWater,iToWater,qs,kFrom,kTo;
  double mass,vol,Cs,Asource;
  double corr;                 //advective correction factor (e.g., retardation)

  double **sv,*Q;
  int    nHRUs=_pModel->GetNumHRUs();
  int    nSVs =_pModel->GetNumStateVars();
  int    nLatConnections=pTransModel->GetNumLatAdvConnections();
  double tstep=Options.timestep;

  sv=new double *[nHRUs];// \todo [optimize]: make a static function
  for(int k=0;k<nHRUs;k++) {
    sv[k]=new double [nSVs];
    for (int i=0;i<nSVs; i++){sv[k][i]=state_vars[k][i]; }
  }

  Q=new double[nLatConnections]; // \todo [optimize]: preallocate in initialize, save as member
  double Afrom, Ato;
  
  //get water fluxes, regenerate system state at start of timestep
  //-------------------------------------------------------------------------
  for(q=0;q<nLatConnections;q++)
  {
    kFrom     =pTransModel->GetLatFromHRU(q);
    kTo       =pTransModel->GetLatToHRU(q);
    iFromWater=pTransModel->GetLatFromWaterIndex(q);
    iToWater  =pTransModel->GetLatToWaterIndex(q);
    qs        =pTransModel->GetLatqsIndex(q);
    Afrom     =pHRUs[kFrom]->GetArea();
    Ato       =pHRUs[kTo  ]->GetArea();

    Q[q]=pModel->GetLatFlow(qs,Options); //volumetric flow [mm-m2/d]

    sv[kFrom][iFromWater]+=Q[q]*tstep/Afrom;//reverse determination of state variable history (i.e., rewinding flow calculations)
    sv[kTo  ][iToWater  ]-=Q[q]*tstep/Ato;

  }

  //Calculate advective mass fluxes
  //-------------------------------------------------------------------------
  for(q=0;q<nLatConnections;q++)
  {
    
    kFrom     =pTransModel->GetLatFromHRU(q);
    kTo       =pTransModel->GetLatToHRU(q);
    iFromWater=pTransModel->GetLatFromWaterIndex(q);
    iToWater  =pTransModel->GetLatToWaterIndex(q);
    Afrom     =pHRUs[kFrom]->GetArea();
    Ato       =pHRUs[kTo  ]->GetArea();

    //cout<<"lat con (iF,iT): ("<<iFromWater<<","<<iToWater<<") kF,kT:("<<kFrom<<","<<kTo<<") icF,icT:("<<_iFromLat[q]<<","<<_iToLat[q]<<")"<<endl;

    //Advection rate calculation rates[q]=dm/dt=Q*C
    mass=0;vol=1;Asource=1.0;
    if(Q[q]>=0)
    {
      mass=sv[kFrom][_iFromLat[q]]; //[mg/m2]
      vol =sv[kFrom][iFromWater]; //[mm]
      Asource=Afrom;
    }
    else if(Q[q]<0)
    {
      mass=sv[kTo][_iToLat[q]  ];
      vol =sv[kTo][iToWater];
      Asource=Ato; 
    }
    corr=1.0;
    corr=pTransModel->GetAdvectionCorrection(_constit_ind,pHRUs[kFrom],iFromWater,iToWater,mass/vol*MM_PER_METER/LITER_PER_M3);//mg/m2/mm->mg/L
    
    exchange_rates[q]=0.0;
    if(vol>1e-6) //note: otherwise Q should generally be constrained to be <vol/tstep & 0.0<rates[q]<(m/tstep/Rf)
    {
      exchange_rates[q]=(Q[q])*mass/vol*corr; //[mm-m2/d]*[mg/m2]/[mm]=[mg/d]

      ExitGracefullyIf(mass<-1e-9,"CmvLatAdvection - negative mass",RUNTIME_ERR);
      if(fabs(exchange_rates[q])>Asource*mass/tstep) { exchange_rates[q]=(Q[q]/fabs(Q[q]))*Asource*mass/tstep; }//emptying out compartment
    }

    //update local mass and volume history
    sv[kFrom][iFromWater  ]-=             Q[q]*tstep/Afrom;
    sv[kTo]  [iToWater    ]+=             Q[q]*tstep/Ato;
    sv[kFrom][_iFromLat[q]]-=exchange_rates[q]*tstep/Afrom;
    sv[kTo]  [_iToLat  [q]]+=exchange_rates[q]*tstep/Ato;

    //Correct for Dirichlet conditions
    //----------------------------------------------------------------------
    if(pTransModel->GetConstituentModel2(_constit_ind)->IsDirichlet(iFromWater,kFrom,tt,Cs))
    {
      Cs*=LITER_PER_M3/MM_PER_METER; //[mg/L]->[mg/mm-m2]
      mass=sv[kFrom][_iFromLat[q]];
      vol =sv[kFrom][iFromWater];
      exchange_rates[nLatConnections+q]+=(Cs*vol-mass)/Options.timestep*Afrom; //[mg/d]
      sv[kFrom][_iFromLat[q]]=Cs*vol;
    }
    if(pTransModel->GetConstituentModel2(_constit_ind)->IsDirichlet(iToWater,kTo,tt,Cs))
    {
      Cs*=LITER_PER_M3/MM_PER_METER; //[mg/L]->[mg/mm-m2]
      mass=sv[kTo][_iToLat[q]];
      vol =sv[kTo][iToWater];
      exchange_rates[2*nLatConnections+q]+=(Cs*vol-mass)/Options.timestep*Ato; //[mg/d]
      sv[kTo][_iToLat[q]]=Cs*vol;
    }

    // \todo [funct]: handle dumping into surface water
    //if ROUTE_NONE, dump surface water compartment to CONSTIT_SW
    /*if (Options.catchment_routing==ROUTE_NONE){
      Afrom=??;
      mass=sv[_iFromLat[3*nAdvConnections]];
      exchange_rates[3*nAdvConnections]=mass/Options.timestep*Afrom;
    }*/
  }

  delete[] Q;
  for (int k=0;k<nHRUs;k++){delete [] sv[k]; }  delete [] sv;
}
