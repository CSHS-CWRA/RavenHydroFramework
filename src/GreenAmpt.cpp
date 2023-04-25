/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Model.h"
#include "Infiltration.h"

double CalcX(const double &tstar);
double epsilon (const double &t,
                const double &alpha,
                const double &Ks,
                const double &w);

///////////////////////////////////////////////////////////////
/// \brief Calculates Green Ampt cumulative infiltration over time step
/// \docpriority Parameter descriptions should be reviewed for this entire file
/// \details Calculates GA Cumulative Infiltration, F(t), using high-order
/// approximation of Barry et al, given the time t, alpha,
/// saturated hydraulic conductivity k, and rainfall rate w
/// \ref Barry et al. \cite Barry2005Aiwr
/// \param &t [in] Time from start of constant rainfall (or time step) [d]
/// \param &alpha [in] GA alpha parameter [mm]
/// \param &Ks [in] Saturated conductivity [mm/d]
/// \param &w [in] Rainfall rate [mm/d]
//
double CmvInfiltration::GreenAmptCumInf(const double &t,     //[d], time from start of rainfall (or time step)
                                        const double &alpha, //[mm]
                                        const double &Ks,    //[mm/d]
                                        const double &w)      //rainfall rate, [mm/d]
{
  //explicit solution to green ampt formulae
  if (Ks==0){return 0.0;}
  if (w==0){return 0.0;}
  if (alpha==0){return 0.0;}

  double tp =alpha*Ks/w/(w-Ks); //ponding time
  if ((t<tp) || (w<Ks))
  {
    return w*t; //cumulative infiltration [mm]
  }
  else
  {
    double tstar,x;
    double Fp=w*tp;
    tstar=Ks/alpha*(t-tp);
    x=-(1.0+Fp/alpha)*exp(-(1.0+tstar+Fp/alpha));


    return alpha*(-1.0-LambertN(x,3)); //cumulative infiltration
  }
}

/////////////////////////////////////////////////////////////////
/// \brief Calclates Green Ampt Runoff
/// \details Calculates runoff [mm/d] from rainfall (rainthru) using GA infiltration method
/// \ref ??
/// \param *state_vars [in] Array of current model state variables for HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Calculated runoff  and infiltration from rainfall [mm/d]
/// \param &rainthru [in] Rate of rainfall [mm/day]
/// \return Runoff [mm/d] from rainfall
//
void CmvInfiltration::GetGreenAmptRunoff (const double            *state_vars,
                                          const CHydroUnit  *pHRU,
                                          const optStruct      &Options,
                                          const time_struct &tt,
                                          double      *rates,
                                          const double      &rainthru) const
{
  double finf; //[mm/day]
  double cumInf=0; //[mm] - cumulative infiltration at time t
  double Keff; //[mm/day] -effective hydraulic conductivity
  double alpha;//|psi_f|*(saturation deficit)
  double deficit;
  double runoff;

  const  soil_struct *S = pHRU->GetSoilProps(0);
  double Ksat                     =S->hydraul_cond;//[mm/d]
  double psi_wf   =S->wetting_front_psi;//[-mm]
  double poro                     =S->porosity;
  double stor     =state_vars[iTo[0]];
  double max_stor =pHRU->GetSoilCapacity(0);
  double Fimp     =pHRU->GetSurfaceProps()->impermeable_frac;
  double initStor =0;

  Keff=Ksat;

  if(type==INF_GA_SIMPLE)
  {
    cumInf=state_vars[pModel->GetStateVarIndex(CUM_INFIL)];//[mm]
    initStor = state_vars[pModel->GetStateVarIndex(GA_MOISTURE_INIT)];//[mm] Storage at start of rainfall

    //redistribute during periods of low rainfall
    if(rainthru <= Keff)
    {
      if(cumInf<=0 /*|| initStor>stor*/) //reset
      {
        cumInf=0;
        initStor=0;
      }
      else //redistribute
      {
        double distribution_rate;

        distribution_rate = max(cumInf*2/3,Keff); // 2/3 is an arbitrary number, and doesn't have a huge impact. Should be a parameter?

        cumInf -= distribution_rate*Options.timestep;
        initStor += distribution_rate*Options.timestep;
        initStor = min(initStor,stor);

        if(cumInf<=0)
        {
          cumInf = 0;
          initStor=0;
        }
      }
      rates[2]= (cumInf-state_vars[pModel->GetStateVarIndex(CUM_INFIL)])/Options.timestep;
      rates[3]= (initStor-state_vars[pModel->GetStateVarIndex(GA_MOISTURE_INIT)])/Options.timestep;
    }
  }

  //no rain, no infiltration
  if (rainthru <= REAL_SMALL){
    rates[0]=0;
    rates[1]=0;
    return;
  }

  deficit=(1.0-stor/max_stor)*poro;
  alpha=psi_wf*deficit;

  if (type==INF_GREEN_AMPT)
  {
    //high-order approximation
    cumInf=GreenAmptCumInf(Options.timestep,alpha,Keff,rainthru);//[mm]
  }
  else if (type==INF_GA_SIMPLE)
  {
    if(cumInf<=0) //First time step of rainfall event
    {
      cumInf=GreenAmptCumInf(Options.timestep,alpha,Keff,rainthru);//[mm]
      initStor=stor; //set initial storage
    }
    else
    {
      deficit=(1.0-initStor/max_stor)*poro;
      alpha=psi_wf*deficit;
    }
  }

  if(cumInf>0)
  {
    finf=Keff*(1.0+alpha/cumInf);//[mm/day]
  }
  else
  {
    finf = 0.0;
  }

  runoff = threshPositive(rainthru-finf); //THRESHOLD BEHAVIOUR
  runoff=(Fimp)*rainthru+(1-Fimp)*runoff; //correct for impermeable surfaces

  if (type==INF_GA_SIMPLE)
  {
    cumInf += (rainthru-runoff)*Options.timestep;
    rates[2]= (cumInf-state_vars[pModel->GetStateVarIndex(CUM_INFIL)])/Options.timestep;
    rates[3]= (initStor-state_vars[pModel->GetStateVarIndex(GA_MOISTURE_INIT)])/Options.timestep;
  }

  rates[0]=rainthru-runoff;
  rates[1]=runoff;
}

///////////////////////////////////////////////////////////////
/// \brief Returns dimensionless time, X
/// \details returns dimensionless time, X, from dimensionless time t*
/// \note used in upscaled GA algorithm (Craig et al 2010)
/// \math \f$ {t}^* ={wt}/{\alpha} \f$
/// \param &tstar Dimensionless time t*
/// \return dimensionless time X
//
double CalcX(const double &tstar)
{
  if (tstar==0.0){return 0.0;}
  return 1.0/(1.0/tstar+1.0);
}

///////////////////////////////////////////////////////////////
/// \brief Returns value of epsilon function used in upscaled GA
/// \ref (Craig et al 2010) \cite Craig2010HP
/// \param &t [in] time since onset of constant rainfall
/// \param &alpha [in] GA Alpha value [mm]
/// \param &Ks [in] Saturated conductivity [mm/d]
/// \param &w [in] rainfall rate [mm/d]
/// \return Epsilon approximation at t
//
double epsilon (const double &t,
                const double &alpha,
                const double &Ks,
                const double &w)
{
  double X=CalcX(w*t/alpha);
  if (X==0){return 0.0;}
  double kp=Ks/w/X;
  if ((kp<=0) || (kp>=1.0)){return 0.0;}
  return 0.135*pow(1-X,0.484)*2.69*pow(1.0-kp,1.74)*pow(kp,0.38);
}

///////////////////////////////////////////////////////////////
/// \brief Upscaled Green Ampt infiltration
/// \details analytical / semi-analytical upscaling of GA infiltration with
/// spatially variable conductivity\n
/// if eps_in:\n=0, then only linear approximation (eqn 10 from Craig et al)\n \cite Craig2010HP
///         =1, then use dirac approx of epsilon integral\n
///         =2, then use 2pt gauss integration of epsilon integral\n
///         =3, then numerically integrate epsilon integral
/// \ref Craig et al., 2010 \cite Craig2010HP
/// \param &t [in] time since onset of rainfall
/// \param &alpha [in] Alpha value
/// \param &mu_Y [in] Mean of log-conductivity
/// \param &sig_Y [in] Standard Deviation of log-conductivity
/// \param &w [in] rainfall rate
/// \param eps_in [in] Input epsilon value
/// \return Upscaled GA infiltration rate[mm/d]
//
double Smooth_GA_k( const double &t,  //time since onset of rainfall
                    const double &alpha,
                    const double &mu_Y,
                    const double &sig_Y,
                    const double &w,
                    int eps_in)
{
  /// \ref eqn 10 from text of Craig et al., 2010
  if (t    ==0){return w;}

  double sum =0.0;
  double X   =CalcX(t*w/alpha);
  double kbar=exp(mu_Y+sig_Y*sig_Y/2);  //average conductivity
  double A   =(log(w*X)-mu_Y)/sig_Y/sqrt(2.0);//upper limit (wX) in normalized Y space
  if (sig_Y==0){A=1e99;}

  sum+=0.5*w*rvn_erfc(A);
  sum+=0.5/X*kbar*rvn_erfc(sig_Y/sqrt(2.0)-A);

  if (eps_in==1)
  { //approximate:
    double P=0.5*(1.0+rvn_erf(A));
    double ktest=0.5*kbar*rvn_erfc(sig_Y/sqrt(2.0)-A)/P;

    sum+=w*P*epsilon(t,alpha,ktest,w);
  }
  else if (eps_in==2)
  { //2-point gauss integration
    double km,k1,k2,halfrange;
    halfrange=0.5*(w*X-exp(mu_Y-3.0*sig_Y));
    km=w*X-halfrange;
    k1=km-0.577*halfrange;
    k2=km+0.577*halfrange;
    if (halfrange>0){
      sum+=halfrange*w*epsilon(t,alpha,k1,w)*log_pdf(k1,mu_Y,sig_Y);
      sum+=halfrange*w*epsilon(t,alpha,k2,w)*log_pdf(k2,mu_Y,sig_Y);
    }
  }
  else if (eps_in==3)
  {
    //numerically integrated (working):
    double maxk=exp(mu_Y+4*sig_Y);
    //lowerswap(maxk,X);
    double mink=exp(mu_Y-4*sig_Y);
    double dk=(maxk-mink)/25000;
    double integ(0.0),ki;
    for (double k=mink;k<maxk;k+=dk)
    {
      ki=(k+0.5*dk);
      //just epsilon:
      integ+=(w*epsilon(t,alpha,ki,w))*log_pdf(ki,mu_Y,sig_Y)*dk;//valid and working!
    }
    sum+=integ;
  }

  return sum;
}


//////////////////////////////////////////////////////////////////
/// \brief Calculates heterogeneous Green Ampt runoff
/// \details Calculates runoff (mm/day) from rainfall (rainthru) using upscaled Green ampt infiltration method
/// \ref Craig et al., 2010 Craig, J.R., G.Liu, and E.D. Soulis, Runoff-infiltration partitioning using an upscaled Green-Ampt solution, Hydrologic Processes, 24(16), p2328-2334, 2010 \cite Craig2010HP
/// \param &rainthru [in] Rainthru [mm/d]
/// \param *S [in] Pointer to soil properties
/// \param thickness [in] Soil layer thickness [mm]
/// \param &soil_water_content [in] Soil water content [mm]
/// \param &tt [in] Current simulation time
/// \param Options [in] Global model options information
/// \return Heterogenous GA runoff [mm/d]
//
double CmvInfiltration::GetHeterogeneousGreenAmptRunoff( const double             &rainthru,//[mm/d]
                                                         const soil_struct *S,
                                                         const double      thickness,//[mm]
                                                         const double      &soil_water_content,//[mm]
                                                         const time_struct &tt,
                                                         const optStruct    &Options) const
{
  double finf;            //[mm/day]
  double alpha;           //|psi_f|*(saturation deficit)
  double mu_Y;

  double poro                     =S->porosity;
  double Ksat                     =S->hydraul_cond;//[mm/d]

  double psi_wf   =S->wetting_front_psi;//[-mm]
  double sig_Y            =S->ksat_std_deviation;//[-]

  alpha=psi_wf*poro*(1.0-soil_water_content/(poro*thickness));

  mu_Y=log(Ksat)-0.5*sig_Y*sig_Y;

  finf=Smooth_GA_k(Options.timestep*0.5,alpha,mu_Y,sig_Y,rainthru,1);//1 point gauss integration @ midpoint

  return threshPositive(rainthru-finf); //THRESHOLD BEHAVIOUR
}
