/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ------------------------------------------------------------------
  Abstraction (partitioning of ponded water/snowmelt to depression
  storage)
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "DepressionProcesses.h"
#include "Model.h"

/*****************************************************************
   Abstraction Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of abstraction constructor
/// \param absttype [in] Selected model of abstraction
//
CmvAbstraction::CmvAbstraction(abstraction_type absttype, CModelABC *pModel)
  :CHydroProcessABC(ABSTRACTION, pModel)
{
  _type=absttype;

  if((_type==ABST_PERCENTAGE) || (_type==ABST_FILL) || (_type==ABST_SCS)){
    CHydroProcessABC::DynamicSpecifyConnections(1);
    //abstraction (ponded-->depression)
    iFrom[0]=pModel->GetStateVarIndex(PONDED_WATER);
    iTo  [0]=pModel->GetStateVarIndex(DEPRESSION);
  }
  else if (_type==ABST_PDMROF) {
    CHydroProcessABC::DynamicSpecifyConnections(2);
    //abstraction (ponded-->depression)
    iFrom[0]=pModel->GetStateVarIndex(PONDED_WATER);
    iTo  [0]=pModel->GetStateVarIndex(DEPRESSION);
    //runoff (ponded-->surface water)
    iFrom[1]=pModel->GetStateVarIndex(PONDED_WATER);
    iTo  [1]=pModel->GetStateVarIndex(SURFACE_WATER);
  }
  else if (_type==ABST_UWFS) {
    CHydroProcessABC::DynamicSpecifyConnections(3);
    //abstraction (ponded-->depression)
    iFrom[0]=pModel->GetStateVarIndex(PONDED_WATER);
    iTo  [0]=pModel->GetStateVarIndex(DEPRESSION);
    //runoff (ponded-->surface water)
    iFrom[1]=pModel->GetStateVarIndex(PONDED_WATER);
    iTo  [1]=pModel->GetStateVarIndex(SURFACE_WATER);
    //adjust minimum deficit
    iFrom[2]=pModel->GetStateVarIndex(MIN_DEP_DEFICIT);
    iTo  [2]=pModel->GetStateVarIndex(MIN_DEP_DEFICIT);
  } 
  else if (_type == ABST_HGDM) {
    CHydroProcessABC::DynamicSpecifyConnections(4);
    //abstraction (ponded-->large or small depressions or outflow)
    iFrom[0]=pModel->GetStateVarIndex(PONDED_WATER);
    iTo  [0]=pModel->GetStateVarIndex(DEPRESSION,0);
    iFrom[1]=pModel->GetStateVarIndex(PONDED_WATER);
    iTo  [1]=pModel->GetStateVarIndex(DEPRESSION,1);
    iFrom[2]=pModel->GetStateVarIndex(PONDED_WATER);
    iTo  [2]=pModel->GetStateVarIndex(SURFACE_WATER);
    /*iFrom[3]=pModel->GetStateVarIndex(CONTRIB_FRAC);
    iTo  [3]=pModel->GetStateVarIndex(CONTRIB_FRAC);*/
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvAbstraction::~CmvAbstraction(){}

//////////////////////////////////////////////////////////////////
/// \brief Initializes abstraction object
//
void   CmvAbstraction::Initialize()
{
}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for abstraction algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by abstraction algorithm (size of aP[] and aPC[])
//
void CmvAbstraction::GetParticipatingParamList(string  *aP , class_type *aPC , int &nP) const
{
  if (_type==ABST_PERCENTAGE)
  {
    nP=1;
    aP[0]="ABST_PERCENT";           aPC[0]=CLASS_LANDUSE;
  }
  else if (_type==ABST_FILL)
  {
    nP=1;
    aP[0]="DEP_MAX";                aPC[0]=CLASS_LANDUSE;
  }
  else if (_type==ABST_SCS)
  {
    nP=2;
    aP[0]="SCS_CN";                 aPC[0]=CLASS_LANDUSE;
    aP[1]="SCS_IA_FRACTION";        aPC[1]=CLASS_LANDUSE;
  }
  else if(_type==ABST_PDMROF)
  {
    nP=2;
    aP[0]="PDMROF_B";               aPC[0]=CLASS_LANDUSE;
    aP[1]="DEP_MAX";                aPC[1]=CLASS_LANDUSE;
  }
  else if(_type==ABST_UWFS)
  {
    nP=4;
    aP[0]="MAX_DEP_AREA_FRAC";      aPC[0]=CLASS_LANDUSE;
    aP[1]="DEP_MAX";                aPC[1]=CLASS_LANDUSE;
    aP[2]="UWFS_B";                 aPC[2]=CLASS_LANDUSE;
    aP[3]="UWFS_BETAMIN";           aPC[3]=CLASS_LANDUSE;
  }
  else if(_type==ABST_HGDM)
  {
    nP=5;
    aP[0]="UPLAND_TO_SMALL";        aPC[0]=CLASS_LANDUSE;
    aP[1]="UPLAND_TO_LARGE";        aPC[1]=CLASS_LANDUSE;
    aP[2]="SMALL_TO_LARGE";         aPC[2]=CLASS_LANDUSE;
    aP[3]="HGDM_P";                 aPC[3]=CLASS_LANDUSE;
    aP[4]="HGDM_LP";                aPC[4]=CLASS_LANDUSE;
    aP[5]="DEP_MAX";                aPC[5]=CLASS_LANDUSE;
  }
  else
  {
    ExitGracefully("CmvAbstraction::GetParticipatingParamList: undefined abstraction algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variable list
///
/// \param absttype [in] Selected abstraction model
/// \param *aSV [out] Reference to array of state variables needed by abstraction algorithm
/// \param *aLev [out] Array of levels of multilevel state variables (or DOESNT_EXIST of single level)
/// \param &nSV [out] Number of participating state variables (length of aSV and aLev arrays)
//
void CmvAbstraction::GetParticipatingStateVarList(abstraction_type absttype, sv_type *aSV, int *aLev, int &nSV)
{
  nSV=2;
  aSV[0]=PONDED_WATER;        aLev[0]=DOESNT_EXIST;
  aSV[1]=DEPRESSION;          aLev[1]=DOESNT_EXIST;
  if(absttype==ABST_PDMROF) {
    nSV=3;
    aSV[2]=SURFACE_WATER;     aLev[2]=DOESNT_EXIST;
  }
  else if(absttype==ABST_UWFS) {
    nSV=4;
    aSV[2]=SURFACE_WATER;     aLev[2]=DOESNT_EXIST;
    aSV[3]=MIN_DEP_DEFICIT;   aLev[3]=DOESNT_EXIST;
  }
  else if(absttype==ABST_HGDM) {
    nSV=4;
    aSV[1]=DEPRESSION;        aLev[1]=0;
    aSV[2]=DEPRESSION;        aLev[2]=1;
    aSV[3]=SURFACE_WATER;     aLev[3]=DOESNT_EXIST;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Returns contributing fraction of depression filled landscape
//
/// \param current_storage [in] depression storage in mm at start of time step
/// \param delta_storage [in] water added/removed to depressions [mm]
/// \param max_storage [in] maximum depresion storage [mm]
/// \param current_contrib_frac [in] current contributing fraction [0..1]
/// \param threshold [in] delta threshold for triggering zero contributing area [mm]
/// from R routine linear_hysteresis_CF.R
//
double HGDMcontrib_fraction(const double &current_storage, 
                            const double &delta_storage,   
                            const double &max_storage,     
                            const double &current_contrib_frac, 
                            const double &threshold = -0.01) 
{  
  double vf1 = current_storage / max_storage; //[0..1]
  double cf1 = current_contrib_frac;          //[0..1] 
  double vf2,cf2;

  if (delta_storage == 0.0) 
  {
    return current_contrib_frac;
  } 
  else 
  {
    if (delta_storage < threshold) 
    {
      //vf2 = (current_storage + delta_storage) / max_storage; //never used
      return 0.0;
    } 
    else 
    {
      vf2 = vf1 + (1.0 - cf1 ) *  (delta_storage / max_storage);
      if (vf1 < 0.999) {
        cf2 = (((1.0 - cf1) * (vf2 - vf1)) / (1.0 - vf1)) + cf1;
        cf2 = min(max(min(cf2, vf2), 0.0), 1.0);
      } 
      else {
        cf2 = 1.0;
      }
      return cf2;
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rates of water loss to abstraction
//
/// \param *state_vars [in] Array of state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time at which abstraction is to be calculated
/// \param *rates [out] rates[0]= rate of abstraction [mm/d]
//
void   CmvAbstraction::GetRatesOfChange( const double        *state_vars,
                                         const CHydroUnit    *pHRU,
                                         const optStruct     &Options,
                                         const time_struct   &tt,
                                               double        *rates) const
{
  double ponded    =state_vars[iFrom[0]]; //[mm]
  double depression=state_vars[iTo  [0]]; //[mm]

  //----------------------------------------------------------------------------
  if      (_type==ABST_PERCENTAGE)
  {
    rates[0]=pHRU->GetSurfaceProps()->abst_percent*ponded/Options.timestep;
  }
  //----------------------------------------------------------------------------
  else if (_type==ABST_FILL)
  { //fills up storage, then stops
    double dep_space;           //[mm] available space left in depression storage

    dep_space = max(pHRU->GetSurfaceProps()->dep_max-depression,0.0);

    rates[0]=min(dep_space,max(ponded,0.0))/Options.timestep;
  }
  //----------------------------------------------------------------------------
  else if (_type==ABST_SCS)
  {
    double S,CN,TR,Ia;

    int condition=2;

    //S should really be calculated as auxilliary land surface param
    TR=pHRU->GetForcingFunctions()->precip_5day/MM_PER_INCH;
    CN=pHRU->GetSurfaceProps()->SCS_CN;

    //correct curve number for antecedent moisture conditions
    if ((tt.month>4) && (tt.month<9)){//growing season?? (northern hemisphere)
      //if (pHRU->GetForcingFunctions()->is_growing_season){
      if      (TR<1.4){condition=1;}
      else if (TR>2.1){condition=3;}
    }
    else{
      if      (TR<0.5){condition=1;}
      else if (TR>1.1){condition=3;}
    }
    if      (condition==1){CN = 5E-05 *pow(CN,3) + 0.0008*pow(CN,2) + 0.4431*CN;}//0.999R^2 with tabulated values (JRC)
    else if (condition==3){CN = 0.0003*pow(CN,3) - 0.0185*pow(CN,2) + 2.1586*CN;}//0.999R^2 with tabulated values (JRC)

    //calculate amount of runoff
    S   =MM_PER_INCH*(1000/CN-10);
    Ia  =(pHRU->GetSurfaceProps()->SCS_Ia_fraction)*S;

    rates[0]=max(Ia,ponded)/Options.timestep;
  }
  //----------------------------------------------------------------------------
  else if(_type==ABST_PDMROF)
  { //uses PDMROF algorithm from Mekonnen et al., Towards an improved land surface scheme for prairie landscapes, Journal of Hydrology 511, 2014
    double b;          //[-] pareto distribution parameter
    double dep_max;    //[mm] maximum depression storage across HRU
    double c_max;      //[mm] maximum local storage capacity in HRU
    double c_star;     //[mm] critical capacity, c*
    double abstracted; //[mm] total amount abstracted

    dep_max=pHRU->GetSurfaceProps()->dep_max;
    b      =pHRU->GetSurfaceProps()->PDMROF_b;

    c_max=(b+1)*dep_max;

    c_star=c_max*(1.0-pow(1.0-(depression/dep_max),1.0/(b+1.0)));

    //analytical evaluation of P*dt-equation 3 of Mekonnen; min() handles case where entire landscape sheds
    abstracted=dep_max*(pow(1.0-(c_star/c_max),b+1.0)-pow(1.0-min(c_star+ponded,c_max)/c_max,b+1.0));
    abstracted=min(abstracted,ponded);

    rates[0]=(       abstracted)/Options.timestep; //abstraction rate [PONDED->DEPRESSION]
    rates[1]=(ponded-abstracted)/Options.timestep; //runoff rate [PONDED->SURFACE_WATER]
  }
  //----------------------------------------------------------------------------
  else if(_type==ABST_UWFS)
  {
    double runoff;     //[mm] total runoff amount, <O_1> (averaged over HRU)
    double P=ponded;   //[mm] total "preciptiation" (averaged over HRU)
    double deltaDmin=0;//[mm] change in Dmin (in wetland; not averaged over HRU)

    double alpha = 1 ;//runoff ratio - here assume infiltration already applied, alpha=1.0 - all runs off
    double b       =pHRU->GetSurfaceProps()->uwfs_b;
    double beta_min=pHRU->GetSurfaceProps()->uwfs_betamin;
    double depfrac =pHRU->GetSurfaceProps()->max_dep_area_frac; //percentage of landscape covered in depressions
    double depmax  =pHRU->GetSurfaceProps()->dep_max;           //[mm] maximum total depression storage on landscape (averaged over HRU)
    if (depfrac==0){runoff=ponded;} //no abstraction, all will run off
    else
    {
      double Davg=(depmax-depression)/depfrac; //mm in wetland (not averaged over HRU as depression storage is)
      double Dmin=state_vars[iFrom[2]];        //=Dmin [mm] if positive,=-Pf if negative
      double Pf=0.0;                           //Fraction full
      double d;                                //[-] deficit distribution parameter
      double beta_ave=beta_min+1/b;            //<beta>
      double a    = b/(P+0.0001);
      double mult = alpha / (beta_ave - 1.0 + alpha);

      if (Dmin < 0.0) { //partially full
         Pf =-Dmin;
      }

      double Pstar=beta_min*P-Dmin;

      //backcalculate d from Dmin (or Pf) and Davg

      d = min((1.0-Pf) / (Davg - Dmin + 0.00001)+ 0.0001, 1000.0);
      //Calculate net runoff (in mm averaged over basin)

      if(Pstar>0.0) {
        runoff =mult*a*d/(a+d)*((1-Pf)*((d*Pstar-1+exp(-d*Pstar))/d/d)+(1-Pf)*(a*Pstar+1)/a/a+ (Pf) * (Pstar + 1.0 / a));
      }
      else {
        runoff=mult*a*d/(a+d)*(exp(a*Pstar)/a/a);
      }

      if (runoff > 0.001)
      {
        if (Pstar > 0.0)
        {
          Pf = 1 - (a / (a +d) * exp(-d*Pstar)) * (1 - Pf);
        }
        else
        {
          Pf = 1 - (a / (a + d) * (1 - Pf) + (d / (a + d)*(1 - Pf) + Pf) * (1 - exp(a * Pstar)));
        }
        deltaDmin = (-Pf) -Dmin;
        Dmin = 0;// Dmin after runoff will be zero means that at least some of the depressions are full
      }
      else
      {
        deltaDmin = -min(beta_ave * P,Dmin); //THis is wrong except in case where runoff~0.0 //Mahkameh -figure out how to update Dmin (done!)
                           //deltaDmin=Dmin (after runoff)- Dmin (before runoff)
                           // note that Dmin =-Pf if Pf at end is non-zero}
        Pf = 0.0; // with positive Dmin Of should be zero
      }
    }

    rates[0]=(ponded- runoff)/Options.timestep; //abstraction rate [PONDED->DEPRESSION]
    rates[1]=(        runoff)/Options.timestep; //runoff rate      [PONDED->SURFACE_WATER]
    rates[2]=(deltaDmin     )/Options.timestep; //change in Dmin   [MIN_DEP_DEFICIT->MIN_DEP_DEFICIT]]
  }
  //----------------------------------------------------------------------------
  else if(_type==ABST_HGDM)
  {
    //runoff from uplands and direct precip/melt on wetlands handled here; Evap handled in OWEVAP_HGDM
    /*int iPond=pModel->GetStateVarIndex(PONDED_WATER);
    int iDepS=pModel->GetStateVarIndex(DEPRESSION,0);
    int iDepL=pModel->GetStateVarIndex(DEPRESSION,1);
    int iSW  =pModel->GetStateVarIndex(SURFACE_WATER);
    int iCFac=pModel->GetStateVarIndex(CONTRIB_FRAC);

    double tstep=Options.timestep;

    double Vs_max  =pHRU->GetSurfaceProps()->dep_max;
    double Vl_max  =pHRU->GetSurfaceProps()->dep_max_large;     //[0..1] maximum depression storage of large gatekeeper
    double f_s_to_l=pHRU->GetSurfaceProps()->small_to_large;    //[0..1] fraction of depression landscape draining to large gatekeeper
    double f_s_to_o=(1-f_s_to_l);                               //[0..1] fraction of depression landscape draining to outlet 
    double fsu     =pHRU->GetSurfaceProps()->HGDM_frac_sm_dep;  //[0..1] fraction of HRU covered in small depressions + contrib areas of depressions
    double flu     =pHRU->GetSurfaceProps()->HGDM_frac_lr_dep;  //[0..1] fraction of HRU covered in large depression + contrib areas of depression
    double fl      =pHRU->GetSurfaceProps()->HGDM_frac_lr_are;  //[0..1] fraction of HRU covered by maximum area of large depression
    double p       =pHRU->GetSurfaceProps()->HGDM_P;
    double p_large =pHRU->GetSurfaceProps()->HGDM_P_LARGE;
    
    //rules: 
    // flu+fsu+fo=1.0
    // fl < flu

    double Atot  =pHRU->GetArea()*M2_PER_KM2;
    double Al_max=fl *Atot; //[km2] - maximum area of gatekeeper
    double Asu   =fsu*Atot; //[km2] - total contributing area + surface area of small depressions
    double Alu   =flu*Atot; //[km2] - total contributing area to large depression + surface area of large
    double As;              //[km2] - total water surface area of small and large depressions (As=A_w in Shook paper) 
    double Al;              //[km2] - surface area of large gatekeeper depression (<Al_max)

    double ponded  =state_vars[iPond];//[mm]
    double Vsmall  =state_vars[iDepS];//[mm]
    double Vlarge  =state_vars[iDepL];//[mm]
    double fsc_last=state_vars[iCFac];//[0..1] contributing fraction from previous time step

    double overflow, to_large(0.0), to_outlet(0.0), to_dep_s(0.0), to_dep_l,excess;
    
    //Small depressions ----------------------------------------------
    //As=pow(Vsmall/Vs_max,2/(p+2)); //from volfrac2areafrac_Clark
    //double dVs = runoff * (Asu - As) / Atot + ponded * (As/Atot);
    double dVs=ponded*(Asu/Atot);   //[mm]

    //As=pow(Vsmall/Vs_max,2/(p+2));
    //double AET=(As/Atot)*PET;
    //AET=min(AET,Vsmall);
    //Vsmall-=AET;
    //double fsc=HGDMcontrib_fraction(Vsmall,dVs,Vs_max,fsc_last); //fsc is only used in next timestep
    
    double fsc=HGDMcontrib_fraction(Vsmall,dVs,Vs_max,fsc_last); //fsc is only used in next timestep
    if (dVs>0)
    {
      Vsmall+=dVs*(1.0-fsc_last);
      overflow=max(Vsmall-Vs_max,0.0);
      Vsmall  =min(Vsmall,Vs_max);
      to_dep_s  =dVs*(1.0-fsc_last)-overflow;        //total ponded->small dep
      to_large  =(overflow+dVs*(fsc_last))*f_s_to_l; //total ponded->large
      to_outlet =(overflow+dVs*(fsc_last))*f_s_to_o; //total ponded->outflow
    }
    
    //Large depression-----------------------------------------------
    if      (Vlarge <= 0.0   ){Al=0.0;}
    else if (Vlarge >= Vl_max){Al=Al_max;}
    else                      {Al=Al_max*pow((Vlarge / Vl_max),2.0/(p_large + 2.0));}
    
    //double dVl=ponded*Al/Atot+to_large+runoff*(Alu-Al)/Atot+to_large;
    double dVl=ponded*(Alu/Atot)+to_large; //any water in local contributing area, which includes Al_max (i.e., Alu>=Al_max, always)
    Vlarge+=dVl;
    excess=0.0;
    if (Vlarge>Vl_max){ 
      excess= (Vlarge-Vl_max);
      to_outlet+=excess;    
      Vlarge=Vl_max;
    } 
    to_dep_l=dVl-excess;

    //No depression--------------------------------------------------
    to_outlet += ponded *(Atot-Asu-Alu)/Atot;

    rates[0]=to_dep_s/tstep;       //PONDED->SMALL DEPRESSIONS
    rates[1]=to_dep_l/tstep;       //PONDED->LARGE DEPRESSIONS
    rates[2]=to_outlet/tstep;      //PONDED->SURFACE_WATER
    rates[3]=(fsc-fsc_last)/tstep; //CONTRIB_FRAC

    //rates[4]=fs*(As+Au*(1.0-f_s_to_o)); //contributing area of depressons
    //rates[5]=Al_max+fsc_last*Asmax*f_s_to_l;        //if outflow from large

    //NOTES: to get runoff right, should correct infiltration algorithm to reduce infil by fs*(Aus/Atot)+(Al/Atot)
*/
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow cannot drain "from" compartment over timestep.
///
/// \param *state_vars [in] Array of state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time at which abstraction is to be calculated
/// \param *rates [out] rates[0]= rate of abstraction [mm/d]
//
void   CmvAbstraction::ApplyConstraints(const double             *state_vars,
                                        const CHydroUnit *pHRU,
                                        const optStruct      &Options,
                                        const time_struct &tt,
                                        double     *rates) const
{

  //cant remove more than is there (should never be an option)
  double pond=max(state_vars[iFrom[0]],0.0);
  rates[0]=min(rates[0],pond/Options.timestep);

  //reaching maximum depression storage
  double stor=max(state_vars[iTo[0]],0.0);
  double max_stor=pHRU->GetStateVarMax(iTo[0],state_vars,Options);
  double deficit = max(max_stor-stor,0.0);
  double abst=min(rates[0],deficit/Options.timestep);
  rates[0]=abst;

}
