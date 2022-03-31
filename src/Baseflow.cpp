/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ------------------------------------------------------------------
  Baseflow
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "SoilWaterMovers.h"

/*****************************************************************
   Baseflow Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the bucket model constructor
///
/// \param btype [in] Baseflow algorithm type
/// \param In_index [in] Index of state variable providing baseflow to surface water
//
CmvBaseflow::CmvBaseflow(baseflow_type  btype,
                         int            In_index)
  :CHydroProcessABC(BASEFLOW)
{
  type =btype;
  ExitGracefullyIf(In_index==DOESNT_EXIST,
                   "CmvBaseflow Constructor: invalid from compartment specified",BAD_DATA);
  CHydroProcessABC::DynamicSpecifyConnections(1);
  iFrom[0]=In_index; iTo[0]=pModel->GetStateVarIndex(SURFACE_WATER);
}
//----------------------------------------------------------------

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
//
CmvBaseflow::~CmvBaseflow(){}


//////////////////////////////////////////////////////////////////
/// \brief Initializes baseflow
/// \details Initializes baseflow by checking for invalid "from" compartment types
/// (must be soil or GW) and unspecified parameters for all HRUs
///
//
void CmvBaseflow::Initialize()
{
  ExitGracefullyIf((pModel->GetStateVarType(iFrom[0])!=SOIL) &&
                   (pModel->GetStateVarType(iFrom[0])!=GROUNDWATER),
                   "CmvBaseflow::Initialize:Baseflow must come from soil or groundwater unit",BAD_DATA);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for baseflow algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by baseflow algorithm (size of aP[] and aPC[])
//
void CmvBaseflow::GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const
{
  if ((type == BASE_LINEAR) || (type == BASE_LINEAR_ANALYTIC) || (type == BASE_LINEAR_CONSTRAIN))
  {
    nP=1;
    aP[0]="BASEFLOW_COEFF"; aPC[0]=CLASS_SOIL;
  }
  else if (type==BASE_CONSTANT)
  {
    nP=1;
    aP[0]="MAX_BASEFLOW_RATE"; aPC[0]=CLASS_SOIL;
  }
  else if (type==BASE_POWER_LAW)
  {
    nP=2;
    aP[0]="BASEFLOW_COEFF";  aPC[0]=CLASS_SOIL;
    aP[1]="BASEFLOW_N";      aPC[1]=CLASS_SOIL;
  }
  else if (type==BASE_VIC)
  {
    nP=2;
    aP[0]="MAX_BASEFLOW_RATE"; aPC[0]=CLASS_SOIL;
    aP[1]="BASEFLOW_N";        aPC[1]=CLASS_SOIL;
  }
  else if (type==BASE_TOPMODEL)
  {
    nP=3;
    aP[0]="MAX_BASEFLOW_RATE"; aPC[0]=CLASS_SOIL;
    aP[1]="BASEFLOW_N";        aPC[1]=CLASS_SOIL;
    aP[2]="TOPMODEL_LAMBDA";   aPC[2]=CLASS_TERRAIN;
  }
  else if (type==BASE_GR4J)
  {
    nP=1;
    aP[0]="GR4J_X3"; aPC[0]=CLASS_SOIL;
  }
  else if (type==BASE_THRESH_POWER)
  {
    nP=3;
    aP[0]="MAX_BASEFLOW_RATE"; aPC[0]=CLASS_SOIL;
    aP[1]="BASEFLOW_N";        aPC[1]=CLASS_SOIL;
    aP[2]="BASEFLOW_THRESH";   aPC[2]=CLASS_SOIL;
  }
  else if(type==BASE_THRESH_STOR)
  {
    nP=2;
    aP[0]="BASEFLOW_COEFF2";   aPC[0]=CLASS_SOIL;
    aP[1]="STORAGE_THRESHOLD"; aPC[1]=CLASS_SOIL;
  }
  else{
    ExitGracefully("CmvBaseflow::GetParticipatingParamList: undefined baseflow algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating state variable list
///
/// \param btype [in] Baseflow algorithm type
/// \param *aSV [out] Array of state variable types needed by baseflow algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by baseflow algorithm (size of aSV[] and aLev[] arrays)
//
void CmvBaseflow::GetParticipatingStateVarList(baseflow_type  btype,sv_type *aSV, int *aLev, int &nSV)
{
  nSV=1;
  aSV[0]=SURFACE_WATER; aLev[0]=DOESNT_EXIST;
  //user specified 'from' compartment (Aquifer layer, lumped landform, or soil layer)
}

//////////////////////////////////////////////////////////////////
/// \brief Finds baseflow rate of change
/// \details returns  of loss from soil or aquifer to surface water [mm/d] 
/// if type==BASE_LINEAR (bucket model) or BASE_LINEAR_ANALYTIC (analytical solution over tstep for bucket model)
///     <ul> <li> flow is linearly  proportional to soil/lumped landform storage [mm] used in PRMS/HBV models, amongst others \cite Hamon1961 </ul>
/// if type==BASE_POWER_LAW
///     <ul> <li> flow is proportional to soil storage [mm] to a power </ul>
/// if type==BASE_CONSTANT
///             <ul> <li> flow is max_baseflow_rate (subject to availability) </ul>
/// if type==BASE_VIC \cite Gao2009
///             <ul> <li> flow is proportional to saturation to a power </ul>
/// if type==BASE_TOPMODEL
///             <ul> <li> flow is calculated using TOPMODEL baseflow method \cite Beven [mm/d] </ul>
///
/// \param *storage [in] Array of state variable values for this HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of loss from baseflow [mm/d]
///
//
void   CmvBaseflow::GetRatesOfChange( const double      *storage,
                                      const CHydroUnit  *pHRU,
                                      const optStruct &Options,
                                      const time_struct &tt,
                                      double     *rates) const
{
  const soil_struct *pSoil=NULL;

  double stor,max_stor=0;
  int    m(0);

  //Move below to CmvBaseflow::GetUsefulParams(pSoil,m,stor,max_stor)?
  //--Obtain critical parameters for calculation----------------------
  stor = storage[iFrom[0]];   //soil layer/groundwater water content [mm]

  sv_type fromType  = pModel->GetStateVarType(iFrom[0]);

  pSoil=NULL;
  max_stor=0.0;
  if      (fromType==SOIL)       {
    m        = pModel->GetStateVarLayer(iFrom[0]); //which soil layer
    pSoil    = pHRU->GetSoilProps(m);
    max_stor = pHRU->GetSoilCapacity(m);  //maximum storage of soil layer [mm]
  }

  stor=min(max(stor,0.0),max_stor); //correct for potentially invalid storage

  //--Calculate Rate of water loss from soil/GW reservoir------------
  //-----------------------------------------------------------------
  if (type==BASE_LINEAR)
  { //BUCKET MODEL, PRMS, HBV MODEL (slow reservoir)
    double K;
    K = pSoil->baseflow_coeff;  //baseflow rate [1/d]

    rates[0]= K * stor;
  }
  //-----------------------------------------------------------------
  else if (type == BASE_LINEAR_CONSTRAIN)
  {
    double K   = pSoil->baseflow_coeff;  //baseflow rate [1/d]
    double sat = stor / max_stor;

    if (sat > pHRU->GetSoilProps(m)->field_capacity)
    {
      rates[0] = K * stor;
    }
  }
  //-----------------------------------------------------------------
  else if (type==BASE_LINEAR_ANALYTIC)
  {
    double K;
    K = pSoil->baseflow_coeff;  //baseflow rate [1/d]

    rates[0]= stor*(1-exp(-K*Options.timestep))/Options.timestep; //Alternate analytical formulation
  }
  //-----------------------------------------------------------------
  else if (type==BASE_CONSTANT)
  { //Constant
    rates[0]= pSoil->max_baseflow_rate;
  }
  //-----------------------------------------------------------------
  else if (type==BASE_POWER_LAW)
  { // HBV MODEL (fast reservoir)
    double K,n;
    K = pSoil->baseflow_coeff;  //[1/d*mm^(-1/n)]
    n = pSoil->baseflow_n;

    rates[0]= K * pow(stor,n);
  }
  //-----------------------------------------------------------------
  else if (type==BASE_VIC)
  { // VIC Model
    double max_rate,n;
    max_rate = pSoil->max_baseflow_rate;
    n        = pSoil->baseflow_n;

    rates[0]=max_rate * pow(stor/max_stor,n);
  }
  //-----------------------------------------------------------------
  else if (type==BASE_TOPMODEL)
  { // TOPMODEL
    double K,n,lambda;
    K      = pSoil->max_baseflow_rate; //baseflow rate constant [mm/d]
    n      = pSoil->baseflow_n;        //baseflow exponent [unitless]
    lambda = pHRU ->GetTerrainProps()->topmodel_lambda;         //mean of the power-transformed topographic index [m]

    rates[0] = K * (max_stor / n * pow(lambda,-n)) * pow(stor/max_stor,n);    // baseflow rate [mm/d]
  }
  //-----------------------------------------------------------------
  else if (type==BASE_GR4J)
  {//GR4J
    double x3= pSoil->GR4J_x3; //GR4J reference storage amount [mm]
    rates[0]=stor*(1.0-pow(1.0+pow(max(stor/x3,0.0),4),-0.25))/Options.timestep;
  }
  //-----------------------------------------------------------------
  else if (type==BASE_THRESH_POWER)
  {

    double K,n,tsat;
    K    = pSoil->max_baseflow_rate;  //baseflow rate [mm/d]
    n    = pSoil->baseflow_n;
    tsat = pSoil->baseflow_thresh; //threshold saturation [-]

    double sat=stor / max_stor;

    rates[0]=0.0;
    if (sat>tsat){ 
      rates[0] = K*pow((sat - tsat) / (1.0-tsat),n); 
      rates[0] = min(rates[0],max_stor*(sat-tsat)/Options.timestep);
    }

  }
  //-----------------------------------------------------------------
  else if(type==BASE_THRESH_STOR)
  {

    double K,tstor;
    K     = pSoil->baseflow_coeff2;   //linear baseflow coeff [1/d]
    tstor = pSoil->storage_threshold; //threshold storage [mm]

    rates[0]=0.0;
    if(stor>tstor) { rates[0] = K*(stor-tstor); }

  }
  //-----------------------------------------------------------------
  else
  {
    ExitGracefully("CmvBaseflow::GetRatesOfChange: undefined baseflow type",BAD_DATA);
  }

}

//////////////////////////////////////////////////////////////////
/// \brief Applies constraints to baseflow
/// \details For all methods, ensures that rate of flow cannot drain "from" compartment over timestep.
///
/// \param *storage [in] array of state variable values for this HRU (of size CModel::_nStateVars)
/// \param *pHRU [in] Pointer to HRU object
/// \param &Options [in] Global model options information
/// \param &tt [in] Current input time structure
/// \param *rates [out] Rates of change of state variables 
//
void   CmvBaseflow::ApplyConstraints( const double     *storage,
                                      const CHydroUnit *pHRU,
                                      const optStruct  &Options,
                                      const time_struct &tt,
                                      double     *rates) const
{
  double min_stor=g_min_storage;

  //cant remove more than is there
  rates[0]=threshMin(rates[0],max(storage[iFrom[0]],min_stor)/Options.timestep,0.0);

  //cannot pull water from river
  rates[0]=max(rates[0],0.0);
}
