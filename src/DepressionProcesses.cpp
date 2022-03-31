/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2019 the Raven Development Team
  ------------------------------------------------------------------
  DepressionOverflow
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "DepressionProcesses.h"

/*****************************************************************
   DepressionOverflow Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of depression overflow constructor
/// \param dtype [in] Selected model of depression overflow
//
CmvDepressionOverflow::CmvDepressionOverflow(depflow_type dtype)
  :CHydroProcessABC(DEPRESSION_OVERFLOW)
{
  _type=dtype;

  CHydroProcessABC::DynamicSpecifyConnections(1);
  //abstraction (ponded-->depression)
  iFrom[0]=pModel->GetStateVarIndex(DEPRESSION);
  iTo  [0]=pModel->GetStateVarIndex(SURFACE_WATER);

}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvDepressionOverflow::~CmvDepressionOverflow(){}

//////////////////////////////////////////////////////////////////
/// \brief Initializes abstraction object
//
void   CmvDepressionOverflow::Initialize()
{
}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for abstraction algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by abstraction algorithm (size of aP[] and aPC[])
//
void CmvDepressionOverflow::GetParticipatingParamList(string  *aP , class_type *aPC , int &nP) const
{
  if (_type==DFLOW_THRESHPOW)
  {
    nP=4;
    aP[0]="DEP_THRESHHOLD";     aPC[0]=CLASS_LANDUSE;
    aP[1]="DEP_N";              aPC[1]=CLASS_LANDUSE;
    aP[2]="DEP_MAX_FLOW";       aPC[2]=CLASS_LANDUSE;
    aP[3]="DEP_MAX";            aPC[3]=CLASS_LANDUSE;

  }
  else if (_type==DFLOW_LINEAR)
  {
    nP=2;
    aP[0]="DEP_THRESHHOLD";     aPC[0]=CLASS_LANDUSE;
    aP[1]="DEP_K";              aPC[1]=CLASS_LANDUSE;
  }
  else if (_type==DFLOW_WEIR)
  {
    nP=2;
    aP[0]="DEP_THRESHHOLD";     aPC[0]=CLASS_LANDUSE;
    aP[1]="DEP_CRESTRATIO";     aPC[1]=CLASS_LANDUSE;
  }
  else
  {
    ExitGracefully("CmvDepressionOverflow::GetParticipatingParamList: undefined depression overflow algorithm",BAD_DATA);
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
void CmvDepressionOverflow::GetParticipatingStateVarList(depflow_type dtype, sv_type *aSV, int *aLev, int &nSV)
{
  nSV=2;
  aSV[0]=DEPRESSION;     aLev[0]=DOESNT_EXIST;
  aSV[1]=SURFACE_WATER;  aLev[1]=DOESNT_EXIST;
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
void   CmvDepressionOverflow::GetRatesOfChange( const double      *state_vars,
                                                const CHydroUnit  *pHRU,
                                                const optStruct   &Options,
                                                const time_struct &tt,
                                                       double     *rates) const
{

  double stor=state_vars[iFrom[0]];

  //----------------------------------------------------------------------------
  if (_type==DFLOW_THRESHPOW)
  {
    double max_flow   =pHRU->GetSurfaceProps()->dep_max_flow;
    double n          =pHRU->GetSurfaceProps()->dep_n;
    double thresh_stor=pHRU->GetSurfaceProps()->dep_threshhold;
    double max_stor   =pHRU->GetSurfaceProps()->dep_max;

    if (max_stor<thresh_stor){
      rates[0] = max_flow;
    }
    else {
      rates[0] = max_flow*pow(max((stor - thresh_stor) / (max_stor - thresh_stor), 0.0), n);
    }
  }
  //----------------------------------------------------------------------------
  else if(_type==DFLOW_LINEAR)
  {
    double thresh_stor=pHRU->GetSurfaceProps()->dep_threshhold;
    double dep_k      =pHRU->GetSurfaceProps()->dep_k;
    
    //rates[0] = dep_k* max(stor-thresh_stor,0.0); //discrete
    rates[0]= max(stor-thresh_stor,0.0)*(1-exp(-dep_k*Options.timestep))/Options.timestep; //analytical formulation
  }
  //----------------------------------------------------------------------------
  else if(_type==DFLOW_WEIR)
  {
    double thresh_stor=pHRU->GetSurfaceProps()->dep_threshhold;
    double dep_rat    =pHRU->GetSurfaceProps()->dep_crestratio;
    double area = pHRU->GetArea();
    double c=dep_rat*sqrt(area); //crest length * weir coeff (C_d*L)
    
    rates[0]=0.666*c/area*sqrt(2*GRAVITY)*pow(max(stor-thresh_stor,0.0),1.5);
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
void   CmvDepressionOverflow::ApplyConstraints(const double              *state_vars,
                                               const CHydroUnit *pHRU,
                                               const optStruct      &Options,
                                               const time_struct &tt,
                                               double     *rates) const
{
  //cant remove more than is there
  rates[0]=threshMin(rates[0],max(state_vars[iFrom[0]]/Options.timestep,0.0),0.0);

}



/*****************************************************************
   Seepage Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of Seepage constructor
/// \param stype [in] Selected model of seepage
//
CmvSeepage::CmvSeepage(seepage_type stype, int iToSoil)
  :CHydroProcessABC(SEEPAGE)
{
  _type=stype;

  CHydroProcessABC::DynamicSpecifyConnections(1);
  //abstraction (ponded-->depression)
  iFrom[0]=pModel->GetStateVarIndex(DEPRESSION);
  iTo  [0]=iToSoil;
  ExitGracefullyIf(iToSoil==DOESNT_EXIST,
                   "CmvSeepage Constructor: invalid to compartment specified",BAD_DATA);
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvSeepage::~CmvSeepage(){}

//////////////////////////////////////////////////////////////////
/// \brief Initializes abstraction object
//
void   CmvSeepage::Initialize()
{
  ExitGracefullyIf((pModel->GetStateVarType(iTo[0])!=SOIL) &&
                  (pModel->GetStateVarType(iTo[0])!=GROUNDWATER),
                  "CmvSeepage::Initialize: seepage must move water to soil or groundwater unit",BAD_DATA);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for abstraction algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by abstraction algorithm (size of aP[] and aPC[])
//
void CmvSeepage::GetParticipatingParamList(string  *aP , class_type *aPC , int &nP) const
{
  if (_type==SEEP_LINEAR)
  {
    nP=1;
    aP[0]="DEP_SEEP_K";     aPC[0]=CLASS_LANDUSE;
  }
  else
  {
    ExitGracefully("CmvSeepage::GetParticipatingParamList: undefined seepage algorithm",BAD_DATA);
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
void CmvSeepage::GetParticipatingStateVarList(seepage_type dtype, sv_type *aSV, int *aLev, int &nSV)
{
  nSV=1;
  aSV[0]=DEPRESSION;     aLev[0]=DOESNT_EXIST;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rates of water loss to seepage
//
/// \param *state_vars [in] Array of state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time at which abstraction is to be calculated
/// \param *rates [out] rates[0]= rate of abstraction [mm/d]
//
void   CmvSeepage::GetRatesOfChange(const double      *state_vars,
                                    const CHydroUnit  *pHRU,
                                    const optStruct   &Options,
                                    const time_struct &tt,
                                           double     *rates) const
{

  double stor=state_vars[iFrom[0]];

  //----------------------------------------------------------------------------
   if(_type==SEEP_LINEAR)
  {
    double dep_k      =pHRU->GetSurfaceProps()->dep_seep_k;
    
    //rates[0] = dep_k* max(stor,0.0); //discrete
    rates[0]= max(stor,0.0)*(1-exp(-dep_k*Options.timestep))/Options.timestep; //analytical formulation
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
void   CmvSeepage::ApplyConstraints(const double      *state_vars,
                                    const CHydroUnit  *pHRU,
                                    const optStruct   &Options,
                                    const time_struct &tt,
                                          double      *rates) const
{
  //cant remove more than is there
  rates[0]=threshMin(rates[0],max(state_vars[iFrom[0]]/Options.timestep,0.0),0.0);

  //can't overfill target storage
  //water flow simply slows (or stops) so that receptor will not overfill during tstep
  double room;
  room=threshMax(pHRU->GetStateVarMax(iTo[0],state_vars,Options)-state_vars[iTo[0]],0.0,0.0);
  rates[0]=threshMin(rates[0],room/Options.timestep,0.0);
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the lake release constructor
/// \param lktype [in] Model of lake release used
/// \param fromIndex [in] Index of lake from which water is released
//
CmvLakeRelease::CmvLakeRelease(lakerel_type lktype)
  :CHydroProcessABC(LAKE_RELEASE)
{
  _type =lktype;
  
  CHydroProcessABC::DynamicSpecifyConnections(1);//nConnections=1
  iFrom[0]=pModel->GetLakeStorageIndex();
  iTo  [0]=pModel->GetStateVarIndex(SURFACE_WATER);     //rates[0]: LAKE->SURFACE_WATER

  ExitGracefullyIf(iFrom[0]==DOESNT_EXIST,
                  "CmvLakeRelease Constructor: invalid lake storage compartment specified",BAD_DATA);
  if(iFrom[0]==SURFACE_WATER){
    string warn="CmvLakeRelease Constructor: shouldn't use lake release if surface water is lake reservoir - it does nothing";
    WriteWarning(warn,true);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvLakeRelease::~CmvLakeRelease(){}

//////////////////////////////////////////////////////////////////
/// \brief Verify that process moves water to surface water from lake
//
void CmvLakeRelease::Initialize()
{
}
//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by baseflow algorithm (size of aP[] and aPC[])
//
void CmvLakeRelease::GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const
{
  if(_type==LAKEREL_LINEAR)//-------------------------------------
  {
    nP=1;
    aP[0]="LAKE_REL_COEFF"; aPC[0]=CLASS_LANDUSE;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Sets reference to state variable types needed by evaporation algorithm
/// \note "From" compartment specified by :LakeStorage command
///
/// \param lktype [in] Model of lake release used
/// \param *aSV [out] Array of state variable types needed 
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by algorithm (size of aSV[] and aLev[] arrays)
//
void CmvLakeRelease::GetParticipatingStateVarList(lakerel_type  lr_type,sv_type *aSV, int *aLev, int &nSV)
{
  nSV=1;
  aSV[0]=SURFACE_WATER;  aLev[0]=DOESNT_EXIST;
  //'from' compartment specified by :LakeStorage command
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rate of loss of water from lake to surface water [mm/d]
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time structure
/// \param *rates [out] rates[0] is rate of loss from lake to surface water [mm/d]
//
void CmvLakeRelease::GetRatesOfChange(const double      *state_vars,
                                      const CHydroUnit  *pHRU,
                                      const optStruct   &Options,
                                      const time_struct &tt,
                                      double            *rates) const
{
  if (!pHRU->IsLake()){return;}

  if(pHRU->IsLinkedToReservoir()){rates[0]= state_vars[iFrom[0]]/Options.timestep;}//reservoir-linked HRUs release directly to surface water

  if (_type==LAKEREL_LINEAR)//-------------------------------------
  {
    double K=pHRU->GetSurfaceProps()->lake_rel_coeff;
    rates[0]=state_vars[iFrom[0]]*(1-exp(-K*Options.timestep))/Options.timestep; // analytical formulation
    //rate can be positive or negative
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow ccannot drain lake over timestep
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time structure
/// \param *rates [out] rates[0] is rate of loss from lake to atmosphere [mm/d]
//
void  CmvLakeRelease::ApplyConstraints( const double      *state_vars,
                                        const CHydroUnit  *pHRU,
                                        const optStruct   &Options,
                                        const time_struct &t,
                                              double      *rates) const
{
  //can't create negative surface water storage
  rates[0]=threshMax(rates[0],-state_vars[iTo[0]]/Options.timestep,0.0);
}
