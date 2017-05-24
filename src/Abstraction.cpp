/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ------------------------------------------------------------------
  Abstraction (partitioning of ponded water/snowmelt to depression
  storage)
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "DepressionProcesses.h"

/*****************************************************************
   Abstraction Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of abstraction constructor
/// \param absttype [in] Selected model of abstraction
//
CmvAbstraction::CmvAbstraction(abstraction_type absttype)
  :CHydroProcessABC(ABSTRACTION)
{
  type=absttype;

  CHydroProcessABC::DynamicSpecifyConnections(1);
  //abstraction (ponded-->depression)
  iFrom[0]=pModel->GetStateVarIndex(PONDED_WATER);
  iTo  [0]=pModel->GetStateVarIndex(DEPRESSION);

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
  if (type==ABST_PERCENTAGE)
  {
    nP=1;
    aP[0]="ABST_PERCENT";                       aPC[0]=CLASS_LANDUSE;
  }
  else if (type==ABST_FILL)
  {
    nP=1;
    aP[0]="DEP_MAX";                              aPC[0]=CLASS_LANDUSE;
  }
  else if (type==ABST_SCS)
  {
    nP=2;
    aP[0]="SCS_CN";                                             aPC[0]=CLASS_LANDUSE;
    aP[1]="SCS_IA_FRACTION";        aPC[1]=CLASS_LANDUSE;
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
  aSV[0]=PONDED_WATER;  aLev[0]=DOESNT_EXIST;
  aSV[1]=DEPRESSION;    aLev[1]=DOESNT_EXIST;
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
void   CmvAbstraction::GetRatesOfChange( const double                   *state_vars,
                                         const CHydroUnit       *pHRU,
                                         const optStruct        &Options,
                                         const time_struct &tt,
                                         double     *rates) const
{
  double ponded=state_vars[iFrom[0]];
  double depression=state_vars[iTo[0]];

  //----------------------------------------------------------------------------
  if (type==ABST_PERCENTAGE)
  {
    rates[0]=pHRU->GetSurfaceProps()->abst_percent*ponded/Options.timestep;
  }
  //----------------------------------------------------------------------------
  else if (type==ABST_FILL)
  { //fills up storage, then stops
    double dep_space;           //[mm] available space left in depression storage

    dep_space = max(pHRU->GetSurfaceProps()->dep_max-depression,0.0);

    rates[0]=min(dep_space,max(ponded,0.0))/Options.timestep;
  }
  //----------------------------------------------------------------------------
  else if (type==ABST_SCS)
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
  rates[0]=threshMin(rates[0],max(state_vars[iFrom[0]]/Options.timestep,0.0),0.0);

  //reaching maximum depression storage
  double max_stor=pHRU->GetStateVarMax(iTo[0],state_vars,Options);
  double abst=threshMin(rates[0],
                        max((max_stor-max(state_vars[iTo[0]],0.0))/Options.timestep,0.0),0.0);
  rates[0]=abst;

}
