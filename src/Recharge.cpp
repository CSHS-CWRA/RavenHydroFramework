/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2020 the Raven Development Team
  ------------------------------------------------------------------
  Recharge
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "SoilWaterMovers.h"
#include "GroundwaterClass.h"

/*****************************************************************
   Recharge Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the standard constructor
/// \param cr_type [in] Model of capillary rise selected
/// \param In_index [in] Soil storage unit index from which water is lost
/// \param Out_index [in] Soil storage unit index to which water rises
//
CmvRecharge::CmvRecharge(recharge_type rech_type,int to_index, int junk)
  :CHydroProcessABC(RECHARGE)
{
  ExitGracefullyIf(to_index==DOESNT_EXIST,
                   "CmvRecharge Constructor: invalid 'to' compartment specified",BAD_DATA);

  CHydroProcessABC::DynamicSpecifyConnections(1);
  iFrom[0]=pModel->GetStateVarIndex(ATMOS_PRECIP); 
  iTo[0]=to_index;
  _type = rech_type;
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the recharge constructor
/// \param se_type [in] Model of recharge selected
//
CmvRecharge::CmvRecharge(recharge_type rech_type, int nConns)
                         :CHydroProcessABC(RECHARGE, nConns)
{
  int iGW, i;
  _type = rech_type;
  if (_type==RECHARGE_CONSTANT)
  {  
    CHydroProcessABC::DynamicSpecifyConnections(1);

    iGW  =pModel->GetStateVarIndex(GROUNDWATER);   //****fix to not always go to top layer of aquifer? 
    iFrom[0]=pModel->GetStateVarIndex(PONDED_WATER);     iTo[0]=iGW;
  }
  else if(_type==RECHARGE_CONSTANT_OVERLAP)
  {
    CHydroProcessABC::DynamicSpecifyConnections(nConns);

    iGW  =pModel->GetStateVarIndex(GROUNDWATER);   //****fix to not always go to top layer of aquifer? 

    for(i = 0; i < nConns; i++)
    {
      iFrom[i]=pModel->GetStateVarIndex(PONDED_WATER);
      iTo[i]=iGW;
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Implementation of default destructor
//
CmvRecharge::~CmvRecharge()
{
}
//////////////////////////////////////////////////////////////////
/// \brief Verifies iFrom - iTo connectivity
/// \details Ensures that water rises from soil layer/ groundwater
/// to another soil layer or groundwater unit
//
void CmvRecharge::Initialize()
{
  ExitGracefullyIf((pModel->GetStateVarType(iTo[0])!=SOIL) &&
                   (pModel->GetStateVarType(iTo[0])!=GROUNDWATER),
                   "CmvRecharge::Initialize:Recharge must go to soil or groundwater unit",BAD_DATA);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for capillary rise algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by capillary rise  algorithm (size of aP[] and aPC[])
//
void CmvRecharge::GetParticipatingParamList(string  *aP , class_type *aPC , int &nP) const
{
  if ((_type==RECHARGE_CONSTANT) || (_type==RECHARGE_FROMFILE))
  {
    nP=0;
  }
  else if(_type==RECHARGE_CONSTANT_OVERLAP)
  {
    nP=0;
  }
  else
  {
    ExitGracefully("CmvRecharge::GetParticipatingParamList: undefined recharge algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
/// \details User specifies from and to compartments, levels not known before construction
///
/// \param *aSV [out] Reference to state variable types needed by cr algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by cr algorithm (size of aSV[] and aLev[] arrays)
//
void CmvRecharge::GetParticipatingStateVarList(recharge_type	r_type,sv_type *aSV, int *aLev, int &nSV)
{
  if (r_type==RECHARGE_FROMFILE){
    nSV=0;
    //soil and atmos_precip will already be in model by default
  }
  else if (r_type==RECHARGE_CONSTANT)
  {  
    nSV=2;
    aSV [0]=PONDED_WATER; aSV	 [1]=GROUNDWATER;  
    aLev[0]=DOESNT_EXIST; aLev [1]=0;	   
  }
  else if (r_type==RECHARGE_CONSTANT_OVERLAP)
  {  
    nSV=2;
    aSV [0]=PONDED_WATER; aSV	 [1]=GROUNDWATER;  
    aLev[0]=DOESNT_EXIST; aLev [1]=0;	   
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rate to aquifer due to Recharge[mm/day]
/// \details Note that the formatting issues below will be resolved once the references have been extracted.. \n \n 
///
/// \param *storage [in] Reference to soil storage from which water rises
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of loss of water from soil to another soil layer [mm/day]
//
void   CmvRecharge::GetRatesOfChange( const double      *storage,
                                      const CHydroUnit  *pHRU,
                                      const optStruct   &Options,
                                      const time_struct &tt,
                                      double            *rates) const
{
  if (pHRU->GetHRUType()!=HRU_STANDARD){return;}//Lakes & glaciers

  if( _type==RECHARGE_FROMFILE){
    rates[0]=pHRU->GetForcingFunctions()->recharge;
  }
  //------------------------------------------------------------
  else if(_type==RECHARGE_CONSTANT)
  {
    // MUST FIX!!================================
    ExitGracefully("CmvRecharge::ApplyConstraints",STUB); return;
    CGroundwaterModel *pGWModel=NULL; //\todo[funct] - must have mechanism for accessing groundwater class to get stress period info
    // MUST FIX!!================================

    CGWStressPeriodClass *pGWSP = NULL;
    pGWSP     = pGWModel->GetGWSPs(0);        //****FIX - should be able to change Stress Periods
    double Rech,poro, HRU_area, spec_stor;
    int HRUid            = pHRU->GetID();
    const soil_struct *S = pHRU->GetAquiferProps(0);
    poro                 = S->porosity;
    HRU_area             = pHRU->GetArea() * M2_PER_KM2 * MM2_PER_M2;        //[mm2]
    spec_stor            = pGWSP->GetGWProperty("SS",HRUid-1); //GWMIGRATE - issue! HRUID!= HRU index k

    Rech = pGWSP->GetGWStruct()->sp_props.recharge[HRUid-1];    //[mm/d]    
    //rates[0] = Rech * HRU_area; //[mm3/d]
    rates[0] = Rech * poro;     //[mm]
    //rates[0] = Rech * HRU_area * spec_stor;
    //cout<<"R: "<<rates[0]<<endl;
  }
  else if(_type==RECHARGE_CONSTANT_OVERLAP)
  {
    ExitGracefully("CmvRecharge::GetRatesOfChange: undefined recharge type:RECHARGE_CONSTANT_OVERLAP",STUB);
  }
  else 
  {
    ExitGracefully("CmvRecharge::GetRatesOfChange: undefined recharge type",BAD_DATA);
  }//end recharge type select*/
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow cannot drain "from" compartment over timestep.
/// Presumes overfilling of "to" compartment is handled using cascade
///
/// \param *storage [in] state variable array for this HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of loss from soil to other soil layer [mm/day]
//
void   CmvRecharge::ApplyConstraints( const double     *state_vars,
                                           const CHydroUnit *pHRU,
                                           const optStruct  &Options,
                                           const time_struct &tt,
                                           double     *rates) const
{
  if (pHRU->GetHRUType()!=HRU_STANDARD){return;}//Lakes & glaciers

  if (_type==RECHARGE_FROMFILE)
  {
    //exceedance of max "to" compartment
    //water flow simply slows (or stops) so that receptor will not overfill during tstep
    rates[0]=threshMin(rates[0],
                     (pHRU->GetStateVarMax(iTo[0],state_vars,Options)-state_vars[iTo[0]])/Options.timestep,0.0);
  }
  else 
  { //GWMIGRATE - CLEAN UP 
	  int HRUid;
	  double room, space, head, headMax, poro, bot;
    // MUST FIX!!================================
    ExitGracefully("CmvRecharge::ApplyConstraints",STUB);
    CGroundwaterModel *pGWModel=NULL; //\todo[funct] - must have mechanism for accessing groundwater class to get stress period info
    // MUST FIX!!================================

	  CGWGeometryClass *pGWGeo;
	  pGWGeo = pGWModel->GetGWGeom();

	  //cant remove more than is there
	  rates[0]=threshMin(rates[0],state_vars[iFrom[0]]/Options.timestep,0.0);

	  //exceedance of max "to" compartment
		//water flow simply slows (or stops) so that receptor will not overfill during tstep
	  if(iTo[0] == pModel->GetStateVarIndex(GROUNDWATER,0))
	  {
	    HRUid  = pHRU->GetID();
	    poro   = pHRU->GetAquiferProps(0)->porosity;
	    bot    = pGWGeo->GetGWGeoProperty("BOT_ELEV",HRUid-1,0); //GWMIGRATE - BIG ISSUE HERE - HRU ID != HRU index - can't just subtract 1
      
	    head     = ((state_vars[iTo[0]] / poro) / MM_PER_METER) + bot;        //convert storage to head
	    headMax  = pHRU->GetStateVarMax(iTo[0],state_vars,Options,HRUid-1);
	    space = threshMax(headMax - head,0.0,0.0);
	    room = (space * MM_PER_METER) * poro;                         //convert head to storage
	    //cout<<"ID: "<<HRUid<<"   poro: "<<poro<<"   head: "<<head<<"   headMax: "<<headMax<<"   space: "<<space<<"   room: "<<room<<endl;
	  }
	  else
	  {
	    room=threshMax(pHRU->GetStateVarMax(iTo[0],state_vars,Options)-state_vars[iTo[0]],0.0,0.0);
	  }
    rates[0]=threshMin(rates[0],room/Options.timestep,0.0);
  }
}
