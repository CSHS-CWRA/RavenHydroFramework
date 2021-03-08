/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright © 2008-2012 the Raven Development Team
------------------------------------------------------------------
  Recharge
----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "SoilWaterMovers.h"
#include "GroundwaterClass.h"

/*****************************************************************
   Drain Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the drain constructor
/// \param se_type [in] Model of drain selected
//
CmvDrain::CmvDrain(drain_type d_type)
                         :CHydroProcessABC(DRAIN)
{
  int iSW;
  iSW  =pModel->GetStateVarIndex(SURFACE_WATER); 
  type = d_type;
  if (type==DRAIN_CONDUCTANCE)
  {  
    CHydroProcessABC::DynamicSpecifyConnections(1);

    iFrom[0]=pModel->GetStateVarIndex(GROUNDWATER); iTo[0]=iSW;
    //iFrom[1]=pModel->GetStateVarIndex(SOIL,0);     iTo[1]=iSW;
    //iFrom[2]=pModel->GetStateVarIndex(SOIL,1);     iTo[2]=iSW;
  }
  else if (type==DRAIN_UFR)
  {  
    CHydroProcessABC::DynamicSpecifyConnections(1);

    iFrom[0]=pModel->GetStateVarIndex(GROUNDWATER);   iTo[0]=iSW;
    //iFrom[1]=pModel->GetStateVarIndex(SOIL,0);        iTo[1]=iSW;
    //iFrom[2]=pModel->GetStateVarIndex(SOIL,1);        iTo[2]=iSW;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
//
CmvDrain::~CmvDrain()
{
}

//////////////////////////////////////////////////////////////////
/// \brief Initializes soil evaporation
//
void CmvDrain::Initialize(){}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for drain algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by drain algorithm (size of aP[] and aPC[])
//
void CmvDrain::GetParticipatingParamList(string  *aP , class_type *aPC , int &nP) const
{
  if (type==DRAIN_CONDUCTANCE)
  {
    nP=2;
    aP[0]="DRAIN_ELEV";       aPC[0]=CLASS_GWSTRESSPERIOD;
    aP[1]="CONDUCTANCE";      aPC[1]=CLASS_GWSTRESSPERIOD;
  }
  else if (type==DRAIN_UFR)
  {
    nP=2;
    aP[0]="D_A";       aPC[0]=CLASS_OVERLAPEXCHANGE;
    aP[1]="D_B";       aPC[1]=CLASS_OVERLAPEXCHANGE;
  }
  else
  {
    ExitGracefully("CmvDrain::GetParticipatingParamList: undefined drain algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
///
/// \param se_type [in] Model of drain used
/// \param *aSV [out] Array of state variable types needed by Drain algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by soil evaporation algorithm (size of aSV[] and aLev[] arrays)
//
void CmvDrain::GetParticipatingStateVarList(drain_type d_type,sv_type *aSV, int *aLev, int &nSV) 
{
  if (d_type==DRAIN_CONDUCTANCE)
  {  
    nSV=2;
    aSV [0]=GROUNDWATER;  aSV [1]=SURFACE_WATER;//aSV	[1]=SOIL;  aSV [2]=SOIL;   aSV [3]=SURFACE_WATER;  
    aLev[0]=0;            aLev[1]=DOESNT_EXIST;//aLev[1]=0;	   aLev[2]=1;      aLev[3]=DOESNT_EXIST;
  }
  else if (d_type==DRAIN_UFR)
  {  
    nSV=2;
    aSV [0]=GROUNDWATER;  aSV [1]=SURFACE_WATER;//aSV	[1]=SOIL;  aSV [2]=SOIL;   aSV [3]=SURFACE_WATER;  
    aLev[0]=0;            aLev[1]=DOESNT_EXIST;//aLev[1]=0;	   aLev[2]=1;      aLev[3]=DOESNT_EXIST;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rates of loss from set of soil layers & aquifer to SW due to Drains [mm/day]
/// \details Note that the formatting issues below will be resolved once the references have been extracted.. \n \n 
/// if type=DRAIN_CONDUCTANCE
///		<ul> <li> Drains work as in MODFLOW with a conductance layer at the topographic surface </ul>
///
/// if type=DISCHARGE_UFR
///		<ul> <li> Drain rate of change determined from Upscaled Flux Relationship </ul>
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of loss from "from" compartment [mm/day]
//
void CmvDrain::GetRatesOfChange (const double		 *state_vars, 
                                    const CHydroUnit *pHRU, 
                                    const optStruct	 &Options,
                                    const time_struct &tt,
                                          double     *rates) const
{
  if (pHRU->GetHRUType()!=HRU_STANDARD){return;}//Lake/Glacier case

  CGWStressPeriodClass *pGWSP = NULL;
  CGWGeometryClass     *pGWG = NULL;


  // MUST FIX!!================================
  ExitGracefully("CmvDrain::GetRatesOFChange",STUB);
  CGroundwaterModel *pGWModel=NULL; //\todo[funct] - must have mechanism for accessing groundwater class to get stress period info
  // MUST FIX!!================================

  pGWSP     = pGWModel->GetGWSPs(0);        // \todo[funct] - should be able to change Stress Periods
  pGWG      = pGWModel->GetGWGeom();

  //------------------------------------------------------------
  if(type==DRAIN_CONDUCTANCE)
  {
    //drainage is any amount of head over topographic surface controlled by conductance layer
    double conductance, drain_elev, head, poro, drain_rate, bot, area;
    int HRUid;
    const soil_struct *S = pHRU->GetAquiferProps(0);
    
    poro        = S->porosity;
    HRUid       = pHRU->GetID();
    conductance = pGWSP->GetGWStruct()->sp_props.conductance[HRUid-1];      
    drain_elev  = pGWSP->GetGWStruct()->sp_props.drain_elev[HRUid-1];
    bot         = pGWG->GetGWGeoProperty("BOT_ELEV",HRUid-1,0);
    area        = pGWG->GetGWGeoProperty("CELL_AREA",HRUid-1,0);

    head = ((state_vars[iFrom[0]] / poro) / MM_PER_METER) + bot;        //convert storage to head

    if(head > drain_elev)
    {
      drain_rate = ((conductance * head) - (conductance * drain_elev)) / area;         //modflow drainage  [m/d]
      rates[0]   = drain_rate * MM_PER_METER * poro;                                 //convert back to storage [mm]
    }
    else
    {
      rates[0] = 0;
    }
  }
  else if(type==DRAIN_UFR)
  {
    double poro;
    int HRUid;
    double head, head_norm, drain_norm, drain_rate;
    double topog_elev, bot, topog_min;
    double ext_depth, hyd_cond;
    double powerA, powerB;
    const soil_struct *S;
    COverlapExchangeClass *pOEP = NULL;
    pOEP      = pGWModel->GetOEs(0);          //****FIX - needs to get more than one overlap exchange class???

    HRUid     = pHRU->GetID();
    S         = pHRU->GetAquiferProps(0);
    bot       = pGWG->GetGWGeoProperty("BOT_ELEV",HRUid-1,0);

    topog_elev = pGWG->GetGWGeoProperty("TOP_ELEV",HRUid-1,0);
    topog_min  = pGWG->GetGWGeoProperty("MIN_TOPOG",HRUid-1,0);
    //cout<<"tmin: "<<topog_min<<endl;
    ext_depth  = pGWSP->GetGWStruct()->sp_props.ext_depth[HRUid-1];
    poro       = S->porosity;
    head       = ((state_vars[iFrom[0]] / poro) / MM_PER_METER) + bot;        //convert storage to head
    //cout<<"t: "<<topog_elev<<"   ed: "<<ext_depth<<"   p: "<<poro<<"   h: "<<head<<endl;
    if(head < (topog_elev - ext_depth))
    {
      rates[0] = 0;
    }
    else
    {
      head_norm  = (head - (topog_min - ext_depth))/(topog_elev - (topog_min - ext_depth));   //normalied head [0..1]
      powerA     = pOEP->GetOEProperty("D_A",HRUid-1);
      powerB     = pOEP->GetOEProperty("D_B",HRUid-1);
      hyd_cond   = pGWSP->GetGWProperty("HYD_COND",HRUid-1);

      //cout<<"hn: "<<head_norm<<"   A: "<<powerA<<"   B: "<<powerB<<"   K: "<<hyd_cond<<endl;
      //normalize head and use powerlaw to get drainage
      drain_norm = powerA * pow(head_norm, powerB);
      drain_rate = drain_norm * hyd_cond;
      rates[0]   = drain_rate * MM_PER_METER * poro;                     //convert back to storage
    }
    //cout<<"r1: "<<rates[0]<<endl;
  }
  else 
  {
    ExitGracefully("CmvDrain::GetRatesOfChange: undefined drain type",BAD_DATA);
  }//end drain type select
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function 
/// \details Ensures that the rate of flow cannot drain "from" compartment over timestep
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of loss from "from" compartment [mm/day]

//
void   CmvDrain::ApplyConstraints( const double		 *state_vars, 
                                      const CHydroUnit *pHRU, 
                                      const optStruct	 &Options,
                                      const time_struct &tt,
                                            double     *rates) const
{
  if (pHRU->GetHRUType()!=HRU_STANDARD){return;}//Lake/Glacier case

  //cant remove more than is there
  rates[0]=threshMin(rates[0],state_vars[iFrom[0]]/Options.timestep,0.0);
  //cout<<"dr: "<<rates[0]<<endl;
}
