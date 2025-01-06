/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2024 the Raven Development Team
  ----------------------------------------------------------------
  Fully coupled snow balance routines
  ----------------------------------------------------------------*/
#include "HydroProcessABC.h"
#include "SnowMovers.h"
#include "Model.h"

//////////////////////////////////////////////////////////////////
/// \brief   Implementation of the constructor
/// \details For combined modelling of melt, refreeze, and energy content
/// \param   bal_type Method of balancing energy selected
//
CmvSnowBalance::CmvSnowBalance(snowbal_type bal_type, CModelABC *pModel)
  :CHydroProcessABC(SNOW_BALANCE, pModel)
{
  type =bal_type;

  int iSnow; //used by all snow balance routines
  iSnow   =pModel->GetStateVarIndex(SNOW);

  if (type==SNOBAL_SIMPLE_MELT)
  {
    int iPond=pModel->GetStateVarIndex(PONDED_WATER);

    CHydroProcessABC::DynamicSpecifyConnections(1);

    iFrom[0]=iSnow;  iTo[0]=iPond;
  }
  else if (type==SNOBAL_COLD_CONTENT)
  {
    int iSnowLiq,iCC,iAtmosEn,iSW;
    iSnowLiq=pModel->GetStateVarIndex(SNOW_LIQ);
    iCC     =pModel->GetStateVarIndex(COLD_CONTENT);
    iAtmosEn=pModel->GetStateVarIndex(ENERGY_LOSSES);
    iSW     =pModel->GetStateVarIndex(SURFACE_WATER);

    CHydroProcessABC::DynamicSpecifyConnections(5);//nConnections=5

    iFrom[0]=iCC;       iTo[0]=iAtmosEn;//rates[0]: COLD_CONTENT->ENERGY_LOSSES
    iFrom[1]=iSnowLiq;  iTo[1]=iSnow;   //rates[1]: SNOW_LIQ->SNOW
    iFrom[2]=iAtmosEn;  iTo[2]=iCC;     //rates[2]: ENERGY_LOSSES->COLD_CONTENT
    iFrom[3]=iSnow;     iTo[3]=iSW;     //rates[3]: SNOW ->SURFACE_WATER
    iFrom[4]=iSnowLiq;  iTo[4]=iSW;     //rates[4]: SNOW_LIQ->SURFACE_WATER
  }
  else if (type==SNOBAL_HBV)
  {
    int iSnowLiq,iSoil,iPond;
    iSnowLiq  =pModel->GetStateVarIndex(SNOW_LIQ);
    iSoil     =pModel->GetStateVarIndex(SOIL,0);
    iPond     =pModel->GetStateVarIndex(PONDED_WATER);

    CHydroProcessABC::DynamicSpecifyConnections(5);

    iFrom[0]=iSnow;       iTo[0]=iSnowLiq;       //rates[0]: SNOW->SNOW_LIQ
    iFrom[1]=iSnow;       iTo[1]=iSoil;          //rates[1]: SNOW->SOIL
    iFrom[2]=iSnowLiq;    iTo[2]=iSoil;          //rates[2]: SNOW_LIQ->SOIL
    iFrom[3]=iSnow;       iTo[3]=iPond;          //rates[3]: SNOW->PONDED_WATER
    iFrom[4]=iSnowLiq;    iTo[4]=iPond;          //rates[4]: SNOW_LIQ->PONDED_WATER
  }
  else if (type==SNOBAL_UBCWM)
  {
    int iSnowLiq,iCumMelt,iPonded,iSnowCov,iColdCont;
    iSnowLiq  =pModel->GetStateVarIndex(SNOW_LIQ);
    iPonded   =pModel->GetStateVarIndex(PONDED_WATER);
    iSnowCov  =pModel->GetStateVarIndex(SNOW_COVER);
    iColdCont =pModel->GetStateVarIndex(COLD_CONTENT);
    iCumMelt  =pModel->GetStateVarIndex(CUM_SNOWMELT);

    int iSnowDef;
    iSnowDef  =pModel->GetStateVarIndex(SNOW_DEFICIT);

    CHydroProcessABC::DynamicSpecifyConnections(7); //nConnections=7

    iFrom[0]=iSnow;       iTo[0]=iSnowLiq;       //rates[0]: SNOW->SNOW_LIQ
    iFrom[1]=iSnowLiq;    iTo[1]=iPonded;        //rates[1]: SNOW_LIQ->PONDED
    iFrom[2]=iSnow;       iTo[2]=iPonded;        //rates[2]: SNOW->PONDED
    iFrom[3]=iColdCont;   iTo[3]=iColdCont;      //rates[3]: CC modification
    iFrom[4]=iSnowCov;    iTo[4]=iSnowCov;       //rates[4]: snow cover modification
    iFrom[5]=iCumMelt;    iTo[5]=iCumMelt;       //rates[5]: cumulative melt modification
    iFrom[6]=iSnowDef;    iTo[6]=iSnowDef;       //rates[6]: snow deficit modification
  }
  else if (type==SNOBAL_CEMA_NEIGE)
  {
    int iSnowCov,iPonded;
    iSnowCov  =pModel->GetStateVarIndex(SNOW_COVER);
    iPonded   =pModel->GetStateVarIndex(PONDED_WATER);

    CHydroProcessABC::DynamicSpecifyConnections(2); //nConnections=2

    iFrom[0]=iSnow;       iTo[0]=iPonded;            //rates[0]: SNOW->PONDED
    iFrom[1]=iSnowCov;    iTo[1]=iSnowCov;           //rates[1]: snow cover modification
  }
  else if(type==SNOBAL_TWO_LAYER)
  {
    int iSnowfall,iPonded,iSLSurf,iSLPack,iCCSurf,iCCPack,iSnowTemp,iCumMelt;
    iSnowfall =pModel->GetStateVarIndex(NEW_SNOW);
    iPonded   =pModel->GetStateVarIndex(PONDED_WATER);
    iSLSurf   =pModel->GetStateVarIndex(SNOW_LIQ,0);
    iSLPack   =pModel->GetStateVarIndex(SNOW_LIQ,1);
    iCCSurf   =pModel->GetStateVarIndex(COLD_CONTENT,0);
    iCCPack   =pModel->GetStateVarIndex(COLD_CONTENT,1);
    iSnowTemp =pModel->GetStateVarIndex(SNOW_TEMP);
    iCumMelt = pModel->GetStateVarIndex(CUM_SNOWMELT);

    CHydroProcessABC::DynamicSpecifyConnections(10); //nConnections=10

    iFrom[0]=iSnowfall;   iTo[0]=iSnow;      //rates[0]: SNOWFALL         -> SNOW
    iFrom[1]=iPonded;     iTo[1]=iSLSurf;    //rates[1]: Rain             -> SNOW_LIQ surface
    iFrom[2]=iSnow;       iTo[2]=iSLSurf;    //rates[2]: SNOW             -> SNOW_LIQ surface
    iFrom[3]=iSLSurf;     iTo[3]=iSLPack;    //rates[3]: SNOW_LIQ surface -> pack
    iFrom[4]=iSLPack;     iTo[4]=iSnow;      //rates[4]: SNOW_LIQ pack    -> SNOW
    iFrom[5]=iSLPack;     iTo[5]=iPonded;    //rates[5]: SNOW_LIQ pack    -> ponded water
    iFrom[6]=iCCSurf;     iTo[6]=iCCSurf;    //rates[6]: CC surface layer
    iFrom[7]=iCCPack;     iTo[7]=iCCPack;    //rates[7]: CC pack layer
    iFrom[8]=iSnowTemp;   iTo[8]=iSnowTemp;  //rates[8]: Snow Temp
    iFrom[9]=iCumMelt;    iTo[9]=iCumMelt;   //rates[9]: Cumulative Melt
  }
  else if(type==SNOBAL_GAWSER)
  {
    int iSnowfall, iPonded, iSWC, iSDEP, iLWC, iAtm;
    iSnowfall =pModel->GetStateVarIndex(NEW_SNOW);
    iPonded   =pModel->GetStateVarIndex(PONDED_WATER);
    iSWC      =pModel->GetStateVarIndex(SNOW);
    iSDEP     =pModel->GetStateVarIndex(SNOW_DEPTH);
    iLWC      =pModel->GetStateVarIndex(SNOW_LIQ);
    iAtm      =pModel->GetStateVarIndex(ATMOSPHERE);

    CHydroProcessABC::DynamicSpecifyConnections(10); //nConnections=10

    iFrom[0]=iSWC;        iTo[0]=iLWC;       //rates[0]: SNOW             -> SNOW_LIQ       (MELT_I in GAWSER)
    iFrom[1]=iPonded;     iTo[1]=iLWC;       //rates[1]: RAIN              -> LWC
    iFrom[2]=iLWC;        iTo[2]=iSWC;       //rates[2]: SNOW_LIQ         -> SNOW           (REFREZ in GAWSER)
    iFrom[3]=iSDEP;       iTo[3]=iSDEP;      //rates[3]: SNOW_DEPTH       -> SNOW_DEPTH     (COMPACTION in GAWSER)
    iFrom[4]=iSDEP;       iTo[4]=iSDEP;      //rates[4]: SNOW_DEPTH       -> SNOW_DEPTH     (SNOW_MELT and SUBLIMATION in GAWSER)
    iFrom[5]=iSDEP;       iTo[5]=iSDEP;      //rates[5]: SNOW_DEPTH       -> SNOW_DEPTH     (SNOWFALL in GAWSER)
    iFrom[6]=iSnowfall;   iTo[6]=iSWC;       //rates[6]: SNOWFALL         -> SNOW           (SNOWFALL in GAWSER)
    iFrom[7]=iSWC;        iTo[7]=iAtm;       //rates[7]: SNOW             -> ATMOSPHERE     (sublimation)
    iFrom[8]=iLWC;        iTo[8]=iPonded;    //rates[8]: LWC              -> PONDED         (EXCESS_LWC in GAWSER)
    iFrom[9]=iSWC;        iTo[9]=iPonded;    //rates[9]: SNOW             -> PONDED         (MELT_2 in GAWSER)

  }
  else if(type==SNOBAL_CRHM_EBSM)
  {
    int iSnowLiq,iPonded,iColdCont;
    iSnowLiq  =pModel->GetStateVarIndex(SNOW_LIQ);
    iPonded   =pModel->GetStateVarIndex(PONDED_WATER);
    iColdCont =pModel->GetStateVarIndex(COLD_CONTENT);

    CHydroProcessABC::DynamicSpecifyConnections(5); //nConnectons=5

    iFrom[0]=iSnow;       iTo[0]=iSnowLiq;       //rates[0]: SNOW->SNOW_LIQ
    iFrom[1]=iSnow;       iTo[1]=iPonded;        //rates[1]: SNOW->PONDED
    iFrom[2]=iSnowLiq;    iTo[2]=iPonded;        //rates[2]: SNOW_LIQ->PONDED
    iFrom[3]=iSnowLiq;    iTo[3]=iSnow;          //rates[3]: SNOW_LIQ->SNOW [refreeze]
    iFrom[4]=iColdCont;   iTo[4]=iColdCont;      //rates[4]: CC modification
  }
  else if(type==SNOBAL_HMETS)
  {
    int iSnowLiq,iCumMelt,iPonded;
    iSnowLiq  =pModel->GetStateVarIndex(SNOW_LIQ);
    iPonded   =pModel->GetStateVarIndex(PONDED_WATER);
    iCumMelt  =pModel->GetStateVarIndex(CUM_SNOWMELT);

    CHydroProcessABC::DynamicSpecifyConnections(5); //nConnectons=5

    iFrom[0]=iSnow;       iTo[0]=iSnowLiq;       //rates[0]: SNOW->SNOW_LIQ
    iFrom[1]=iSnow;       iTo[1]=iPonded;        //rates[1]: SNOW->PONDED
    iFrom[2]=iSnowLiq;    iTo[2]=iPonded;        //rates[2]: SNOW_LIQ->PONDED
    iFrom[3]=iSnowLiq;    iTo[3]=iSnow;          //rates[3]: SNOW_LIQ->SNOW [refreeze]
    iFrom[4]=iCumMelt;    iTo[4]=iCumMelt;       //rates[4]: modification of cumulative snow
  }
  else{
    ExitGracefully("CmvSnowBalance::Constructor: undefined snow balance type",BAD_DATA);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief   Implementation of the constructor
/// \details For combined modelling of melt, refreeze, and energy content
/// \param   bal_type Method of balancing energy selected
/// \param   iSnowTo index of 'To' state variable (should be ponded water or snow_liq)
//
CmvSnowBalance::CmvSnowBalance(snowbal_type bal_type, int iSnowTo, CModelABC *pModel) :
  CHydroProcessABC(SNOW_BALANCE, pModel)
{
  type = bal_type;

  int iSnow; //used by all snow balance routines
  iSnow = pModel->GetStateVarIndex(SNOW);

  if (type == SNOBAL_SIMPLE_MELT)
  {
    sv_type typ=pModel->GetStateVarType(iSnowTo);
    if ((typ != PONDED_WATER) && (typ != SNOW_LIQ)){
      ExitGracefully("CmvSnowBalance Constructor: SNOBAL_SIMPLE_MELT target should be either PONDED_WATER or SNOW_LIQ",BAD_DATA_WARN);
    }
    CHydroProcessABC::DynamicSpecifyConnections(1);
    iFrom[0] = iSnow;  iTo[0] = iSnowTo;
  }
  else{
    ExitGracefully("CmvSnowBalance::Constructor: incorrect constructor for this type.",RUNTIME_ERR);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvSnowBalance::~CmvSnowBalance(){}

//////////////////////////////////////////////////////////////////
/// \brief Initializes snow balance modelling object
//
void CmvSnowBalance::Initialize(){}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for snow balance algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by snow balance algorithm (size of aP[] and aPC[])
//
void CmvSnowBalance::GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const
{
  if (type==SNOBAL_SIMPLE_MELT)
  {
    nP=0;
  }
  else if (type==SNOBAL_COLD_CONTENT)
  {
    nP=1;
    aP[0]="SNOW_SWI";        aPC[0]=CLASS_GLOBAL;
  }
  else if (type==SNOBAL_HBV)
  {
    nP=2;
    aP[0]="REFREEZE_FACTOR"; aPC[0]=CLASS_LANDUSE;
    aP[1]="SNOW_SWI";        aPC[1]=CLASS_GLOBAL;
  }
  else if (type==SNOBAL_UBCWM)
  {
    nP=3;
    aP[0]="CC_DECAY_COEFF";   aPC[0]=CLASS_LANDUSE;
    aP[1]="SNOW_SWI";         aPC[1]=CLASS_GLOBAL;
    aP[2]="SNOW_PATCH_LIMIT"; aPC[2]=CLASS_LANDUSE;
  }
  else if (type==SNOBAL_CEMA_NEIGE)
  {
    nP=1;
    aP[0]="AVG_ANNUAL_SNOW";  aPC[0]=CLASS_GLOBAL;
  }
  else if (type==SNOBAL_TWO_LAYER)
  {
    nP=2;
    aP[0]="MAX_SWE_SURFACE";    aPC[0]=CLASS_GLOBAL;
    aP[1]="SNOW_SWI";           aPC[1]=CLASS_GLOBAL;
  }
  else if (type==SNOBAL_GAWSER)
  {
    nP=3;
    aP[0]="SNOW_SWI";         aPC[0]=CLASS_GLOBAL;
    aP[1]="REFREEZE_FACTOR";  aPC[1]=CLASS_LANDUSE;
    aP[2]="MELT_FACTOR";      aPC[2]=CLASS_LANDUSE;
  }
  else if(type==SNOBAL_CRHM_EBSM)
  {
    nP=1;
    aP[0]="SNOW_SWI";         aPC[0]=CLASS_GLOBAL;
  }
  else if(type==SNOBAL_HMETS)
  {
    nP=6;
    aP[0]="REFREEZE_FACTOR";     aPC[0]=CLASS_LANDUSE;
    aP[1]="REFREEZE_EXP";        aPC[1]=CLASS_LANDUSE;
    aP[2]="DD_REFREEZE_TEMP";    aPC[2]=CLASS_LANDUSE;
    aP[3]="SNOW_SWI_MIN";        aPC[3]=CLASS_GLOBAL;
    aP[4]="SNOW_SWI_MAX";        aPC[4]=CLASS_GLOBAL;
    aP[5]="SWI_REDUCT_COEFF";    aPC[5]=CLASS_GLOBAL;
  }
  else
  {
    ExitGracefully("CmvSnowBalance::GetParticipatingParamList: undefined snow balance algorithm",BAD_DATA);
  }
}


//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
///
/// \param bal_type [in] Snow energy/mass balance model selected
/// \param *aSV [out] Array of state variable types needed by energy balance algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by energy balance algorithm (size of aSV[] and aLev[] arrays)
//
void CmvSnowBalance::GetParticipatingStateVarList(snowbal_type bal_type,
                                                  sv_type *aSV, int *aLev, int &nSV)
{
  if (bal_type==SNOBAL_SIMPLE_MELT)
  {
    nSV=1;
    aSV[0]=SNOW;          aLev[0]=DOESNT_EXIST;
    ////user specified 'to' compartment (ponded water or snow_liq)
  }
  else if (bal_type==SNOBAL_COLD_CONTENT)
  {
    nSV=5;
    aSV[0]=SNOW_LIQ;      aLev[0]=DOESNT_EXIST;
    aSV[1]=COLD_CONTENT;  aLev[1]=DOESNT_EXIST;
    aSV[2]=ENERGY_LOSSES; aLev[2]=DOESNT_EXIST;
    aSV[3]=SURFACE_WATER; aLev[3]=DOESNT_EXIST;
    aSV[4]=SNOW;          aLev[4]=DOESNT_EXIST;
  }
  else if (bal_type==SNOBAL_HBV)
  {
    nSV=3;
    aSV[0]=SNOW;          aLev[0]=DOESNT_EXIST;
    aSV[1]=SNOW_LIQ;      aLev[1]=DOESNT_EXIST;
    aSV[2]=SNOW_DEPTH;    aLev[2]=DOESNT_EXIST;
    //plus top soil, which is always present
  }
  else if (bal_type==SNOBAL_UBCWM)
  {
    nSV=8;
    aSV[0]=SNOW;          aLev[0]=DOESNT_EXIST;
    aSV[1]=SNOW_LIQ;      aLev[1]=DOESNT_EXIST;
    aSV[2]=PONDED_WATER;  aLev[2]=DOESNT_EXIST;
    aSV[3]=COLD_CONTENT;  aLev[3]=DOESNT_EXIST;
    aSV[4]=SNOW_COVER;    aLev[4]=DOESNT_EXIST;
    aSV[5]=CUM_SNOWMELT;  aLev[5]=DOESNT_EXIST;
    aSV[6]=SNOW_ALBEDO;   aLev[6]=DOESNT_EXIST;//Used, but not modified
    aSV[7]=SNOW_DEFICIT;  aLev[7]=DOESNT_EXIST;
  }
  else if (bal_type==SNOBAL_CEMA_NEIGE)
  {
    nSV=3;
    aSV[0]=SNOW;          aLev[0]=DOESNT_EXIST;
    aSV[1]=PONDED_WATER;  aLev[1]=DOESNT_EXIST;
    aSV[2]=SNOW_COVER;    aLev[2]=DOESNT_EXIST;
  }
  else if (bal_type==SNOBAL_TWO_LAYER)
  {
    nSV=9;
    aSV[0]=NEW_SNOW;      aLev[0]=DOESNT_EXIST;
    aSV[1]=SNOW;          aLev[1]=DOESNT_EXIST;
    aSV[2]=SNOW_LIQ;      aLev[2]=0;
    aSV[3]=SNOW_LIQ;      aLev[3]=1;
    aSV[4]=COLD_CONTENT;  aLev[4]=0;
    aSV[5]=COLD_CONTENT;  aLev[5]=1;
    aSV[6]=PONDED_WATER;  aLev[6]=DOESNT_EXIST;
    aSV[7]=SNOW_TEMP;     aLev[7]=DOESNT_EXIST;
    aSV[8]=CUM_SNOWMELT;  aLev[8] = DOESNT_EXIST;
  }
  else if (bal_type==SNOBAL_GAWSER)
  {
    nSV=5;
    aSV[0]=NEW_SNOW;      aLev[0]=DOESNT_EXIST;
    aSV[1]=SNOW;          aLev[1]=DOESNT_EXIST;
    aSV[2]=SNOW_LIQ;      aLev[2]=DOESNT_EXIST;
    aSV[3]=SNOW_DEPTH;    aLev[3]=DOESNT_EXIST;
    aSV[4]=PONDED_WATER;  aLev[4]=DOESNT_EXIST;
  }
  else if(bal_type==SNOBAL_CRHM_EBSM)
  {
    nSV=4;
    aSV[0]=SNOW;          aLev[0]=DOESNT_EXIST;
    aSV[1]=SNOW_LIQ;      aLev[1]=DOESNT_EXIST;
    aSV[2]=PONDED_WATER;  aLev[2]=DOESNT_EXIST;
    aSV[3]=COLD_CONTENT;  aLev[3]=DOESNT_EXIST;
  }
  else if(bal_type==SNOBAL_HMETS)
  {
    nSV=4;
    aSV[0]=SNOW;          aLev[0]=DOESNT_EXIST;
    aSV[1]=SNOW_LIQ;      aLev[1]=DOESNT_EXIST;
    aSV[2]=PONDED_WATER;  aLev[2]=DOESNT_EXIST;
    aSV[3]=CUM_SNOWMELT;  aLev[3]=DOESNT_EXIST;
  }
  else{
    nSV=0;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rates of change in all state variables modeled during snow balance calculations
/// \details for combined modeling of melt, refreeze, and energy content \n
/// if type==SNOBAL_SIMPLE_MELT:
///             <ul> <li> simple application of potential melt to SWE; no other variables tracked </ul>
/// if type==SNOBAL_COLD_CONTENT:
///             <ul> <li> adapted from Brook90 SNOENRGY and SNOPACK routines </ul> \cite Federer2010
/// if type==SNOBAL_HBV-EC:
///             <ul> <li> adapted from HBV-Environment Canada </ul>
/// if type==SNOWBAL_UBC:
///             <ul> <li> adapted from the UBC Watershed Model (c) Michael Quick </ul>
///
/// \param *state_var [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Rate of loss of water from snow melt/refreeze combined [mm/d]
//
void CmvSnowBalance::GetRatesOfChange(const double               *state_var,
                                      const CHydroUnit *pHRU,
                                      const optStruct    &Options,
                                      const time_struct &tt,
                                      double     *rates) const
{
  //------------------------------------------------------------
  if (type==SNOBAL_SIMPLE_MELT)
  {
    rates[0]=max(pHRU->GetForcingFunctions()->potential_melt,0.0);
  }
  //------------------------------------------------------------
  else if (type==SNOBAL_CEMA_NEIGE)
  {
    double pot_melt,snow_cov,SWE;
    double avg_annual_snow = pModel->GetGlobalParams()->GetParams()->avg_annual_snow;
    double snotemp=pHRU->GetSnowTemperature();
    pot_melt=0.0;
    if (snotemp==FREEZING_TEMP){
      pot_melt=max(pHRU->GetForcingFunctions()->potential_melt,0.0);//positive
    }
    SWE=state_var[iFrom[0]];

    snow_cov=min(SWE/avg_annual_snow,1.0); //way too small for Canadian basins!

    rates[0]=(0.9*snow_cov+0.1)*min(SWE/Options.timestep,pot_melt);                       //melt
    rates[1]=(snow_cov-state_var[iFrom[1]])/Options.timestep;   //change in snow cover
  }
  //------------------------------------------------------------
  else if (type==SNOBAL_COLD_CONTENT)
  {
    ColdContentBalance(state_var,pHRU,Options,tt,rates);
  }
  //-----------------------------------------------------------
  else if(type==SNOBAL_TWO_LAYER)
  {
    TwoLayerBalance(state_var,pHRU,Options,tt,rates);
  }
  //-------------------------------------------------------------
  else if(type==SNOBAL_GAWSER)
  {
    GawserBalance(state_var,pHRU,Options,tt,rates);
  }
  //------------------------------------------------------------
  else if (type==SNOBAL_HBV)
  {
    double Ka    =pHRU->GetSurfaceProps()->refreeze_factor;//[mm/d/K]
    double Ta    =pHRU->GetForcingFunctions()->temp_daily_ave;
    double tstep =Options.timestep;

    double melt,refreeze, liq_cap,to_liq,overflow;

    double SWE  =state_var[iFrom[0]];//snow, SWE mm
    double SL   =state_var[iTo  [0]];//liquid snow, mm
    double SD   =state_var[pModel->GetStateVarIndex(SNOW_DEPTH)];//snow depth, mm

    melt=max(pHRU->GetForcingFunctions()->potential_melt,0.0);
    melt=min(melt,SWE/tstep); //positive, constrained by available snow
    SWE-=melt*tstep;

    refreeze=Ka*max(FREEZING_TEMP-Ta,0.0);//positively valued
    refreeze=max(min(SL/tstep,refreeze),0.0);
    SL-=refreeze*tstep;

    liq_cap=CalculateSnowLiquidCapacity(SWE, SD, pModel);
    to_liq=min(melt,max(liq_cap-SL,0.0)/tstep);
    SL+=to_liq*tstep;

    overflow=max((SL-liq_cap)/tstep,0.0);

    rates[0] =-refreeze;            //SNOW_LIQ-->SNOW
    rates[0]+=to_liq;               //SNOW->SNOW_LIQ
    rates[1] =max(melt-to_liq,0.0); //SNOW->SOIL[0]
    rates[2] =overflow;             //SNOW_LIQ->SOIL[0]

    if (pHRU->GetHRUType() != HRU_STANDARD) {//lakes & wetlands - don't send to soil
      rates[3]=rates[1]; rates[1]=0.0;
      rates[4]=rates[2]; rates[2]=0.0;
    }
  }
  //------------------------------------------------------------
  else if(type==SNOBAL_CRHM_EBSM)
  {
    CRHMSnowBalance(state_var,pHRU,Options,tt,rates);
  }
  //------------------------------------------------------------
  else if(type==SNOBAL_HMETS)
  {
    double SWE     =state_var[iFrom[0]];//snow, SWE mm
    double SL      =state_var[iFrom[3]];//liquid snow, mm
    double cum_melt=state_var[iFrom[4]];//cumulative melt mm

    double Kf     = pHRU->GetSurfaceProps()->refreeze_factor;
    double Tbf    = pHRU->GetSurfaceProps()->DD_refreeze_temp;
    double exp_fe = pHRU->GetSurfaceProps()->refreeze_exp;
    double fcmin  = pModel->GetGlobalParams()->GetParams()->snow_SWI_min;
    double fcmax  = pModel->GetGlobalParams()->GetParams()->snow_SWI_max;
    double Ccum   = pModel->GetGlobalParams()->GetParams()->SWI_reduct_coeff;

    double potmelt=max(pHRU->GetForcingFunctions()->potential_melt,0.0);
    double Tdiurnal=0.5*(pHRU->GetForcingFunctions()->temp_daily_ave+pHRU->GetForcingFunctions()->temp_daily_min);

    double pot_freeze,freeze,melt,SWI,deficit,to_liq;
    double tstep=Options.timestep;

    SL=max(SL,0.0);
    pot_freeze = Kf*pow(max(Tbf-Tdiurnal,0.0),exp_fe);   //refreeze rate, mm/d

    freeze = min(pot_freeze*tstep,SL);                   // Cannot freeze more water than available
    SL -=freeze;                                         // Adjust free water in snowpack
    SWE+=freeze;                                         // Water frozen becomes part of the snowpack

    melt = min(potmelt*tstep,SWE);
    SWE     -=melt;
    cum_melt+=melt;

    if(SWE<=0.0){ cum_melt=0.0; }

    SWI = max(fcmin,fcmax*(1-Ccum*cum_melt));        // Calculation of the water retention capacity of the snowpack.

    deficit=max(SWI*SWE-SL,0.0);
    to_liq=min(melt,deficit);                        //melt first fills deficit
    melt-=to_liq;
    SL  +=to_liq;

    //what if SL>SWI*SWE?
    double ripe=max(SL-SWI*SWE,0.0);
    SL-=ripe;

    rates[0]=to_liq/tstep;   //SNOW->SNOW_LIQ
    rates[1]=melt  /tstep;   //SNOW->PONDED_WATER (runoff)
    rates[2]=ripe;           //SNOW_LIQ->PONDED (not explicit here)
    rates[3]=freeze/tstep;   //SNOW_LIQ->SNOW (refreeze)
    rates[4]=(melt+to_liq)/tstep;   //CUM_SNOWMELT->CUM_SNOWMELT

    if(SWE<=0.0){ rates[4]=-state_var[iFrom[4]]/tstep; } //reset to zero as needed
  }
  //------------------------------------------------------------
  else if (type == SNOBAL_UBCWM)
  {
    double SWE = state_var[iFrom[0]]; //snow as SWE [mm]
    double Sliq = state_var[iFrom[1]]; //liquid snow [mm]
    double CC = state_var[iFrom[3]]; //cold content [mm]
    double coverage = state_var[iFrom[4]]; //snow coverage
    double cum_melt = state_var[iFrom[5]]; //cumulative melt [mm]

    double pot_melt;                     //potential melt [mm]
    double snowmelt = 0;                   //total snowmelt [mm]
    double tstep = Options.timestep;

    //reset cumulative melt to zero in October
    if ((pHRU->GetCentroid().latitude < 0) && (tt.month == 4) && (tt.day_of_month == 1)){
      cum_melt = 0.0;
    }
    else if ((tt.month == 10) && (tt.day_of_month == 1)){
      cum_melt = 0.0;
    }

    pot_melt = pHRU->GetForcingFunctions()->potential_melt*Options.timestep;

    //Reduce/increase cold content (always leads to pot_melt>=0, CC>=0)
    if (CC > pot_melt){
      CC -= pot_melt;
      pot_melt = 0;
    }
    else{
      pot_melt -= CC;
      CC = 0.0;
    }

    double snowpatch_limit = pHRU->GetSurfaceProps()->snow_patch_limit;

    if (snowpatch_limit > 0.0)
    {
      double snowmeltpatch = 0;
      double AF1, AF2, WSN;

      //if melting will expose ground, melt away snow up to the point of patching and reduce the potential melt
      if ((SWE > snowpatch_limit) && (SWE - pot_melt <= snowpatch_limit)){
        snowmeltpatch = SWE - snowpatch_limit; //calculate snowmelt above snow patch limit (snowmeltpatch)
        SWE = snowpatch_limit;
        pot_melt = pot_melt - snowmeltpatch; //reduce the potential melt amt of snow melted above the patch limit.
      }

      //if ground is already exposed (or on the cusp of being exposed)
      if (SWE <= snowpatch_limit)
      {
        //redistribute snow to a wedge
        WSN = 2.0 * sqrt(SWE * snowpatch_limit);
        lowerswap(pot_melt, WSN);

        if (SWE > 0)
        {
          AF1 = WSN / (2 * snowpatch_limit); //AF1 is area fraction of snowcover
          AF2 = (WSN - pot_melt) / (2 * snowpatch_limit);

          //so the snowmelt already allows for snow fraction
          snowmelt = (WSN * AF1) / 2 - (WSN - pot_melt) / 2 * AF2 + snowmeltpatch;
          WSN -= pot_melt;
          upperswap(WSN, 0.0);

          //reaccumulation of snow
          if (WSN < 2 * snowpatch_limit)
          {
            SWE = (WSN * WSN * 0.5) / (2 * snowpatch_limit);
          }
        }
        else{
          snowmelt = 0;
        }
      }
      else
      {
        snowmelt = max(min(pot_melt, SWE), 0.0); //assumes daily ts for now
        SWE -= snowmelt;
      }
    }
    else //snowpatch_limit = 0
    {
      snowmelt = max(min(pot_melt, SWE), 0.0);//assumes daily ts for now
      SWE -= snowmelt;
    }

    double loss = 0.0;
    double transfer = 0;//amount transferred from frozen to liquid snow
    double orig_SWE = state_var[iFrom[0]];
    bool   oldversion = false;

    double snow_d=state_var[iFrom[6]];//JRC : new

    if (oldversion)
    {
      if ((snowmelt > 0) && (SWE > REAL_SMALL)){ loss = snowmelt*(Sliq / orig_SWE); }
      Sliq -= loss;//loss of water holding capacity

      if (SWE < REAL_SMALL){
        loss += Sliq;
        Sliq = 0;
        CC = 0;
        transfer = 0.0;
      }
      else{
        double snDef = max(CalculateSnowLiquidCapacity(SWE, 0.0, pModel) - Sliq, 0.0);

        if (snowmelt > snDef)
        {
          Sliq += snDef;
          snowmelt -= snDef;
          transfer = snDef;
          snDef = 0.0;
        }
        else{
          snDef -= (snowmelt);
          Sliq += (snowmelt);
          transfer = (snowmelt);
          snowmelt = 0.0;
          //loss    =0.0;

        }
      }
    }
    else //new version
    {
      if (SWE < REAL_SMALL){
        snow_d=0.0;
        CC = 0;
      }
      else
      {
        if (snowmelt > snow_d){
          snowmelt-=snow_d;
          SWE     +=snow_d;
          snow_d=0;
        }
        else{
          snow_d-=snowmelt;
          SWE   +=snowmelt;
          snowmelt=0;
        }
      }
    }

    cum_melt+=-(SWE-state_var[iFrom[0]]);
    //Calculate snow coverage (this should be compared to AF1 which is the area fraction of snow cover, and non-linear)
    if (snowpatch_limit!=0.0){
      coverage=SWE/snowpatch_limit;
      upperswap(coverage,0.0);
      lowerswap(coverage,1.0);
    }
    else{
      if (SWE>0.0){coverage=1.0;}
      else        {coverage=0.0;}
    }

    //decay of cold content
    CC *=exp(-pHRU->GetSurfaceProps()->CC_decay_coeff*tstep);

    rates[0]=transfer/tstep;                       //rates[0]: SNOW->SNOW_LIQ
    rates[1]=loss/tstep;                           //rates[1]: SNOW_LIQ->PONDED  (net melt)
    rates[2]=-(SWE     -state_var[iFrom[2]]+transfer)/tstep; //rates[2]: SNOW->PONDED  (net melt)
    rates[3]=(CC      -state_var[iFrom[3]])/tstep; //rates[3]: CC modification
    rates[4]=(coverage-state_var[iFrom[4]])/tstep; //rates[4]: snow cover modification
    rates[5]=(cum_melt-state_var[iFrom[5]])/tstep; //rates[5]: cumulative melt modification

    if (!oldversion){
      rates[0] = 0.0;                                    //don't care about liquid content
      rates[1] = 0.0;
      rates[2] = snowmelt / tstep;                       //rates[2]: SNOW->PONDED
      rates[6] = (snow_d - state_var[iFrom[6]]) / tstep; //rates[6]: snow deficit modification
    }


  }//end UBC Melt

};
/*****************************************************************
   Cold Content Balance :Get Rates of Change
------------------------------------------------------------------
adapted from Brook90 SNOENRGY and SNOPACK routines
*****************************************************************/
//////////////////////////////////////////////////////////////////
/// \docminor This function needs better documentation
/// \brief Balances cold content from snow melt/refreeze
/// \ref Adapted from Brook90 SNOENRGY and SNOPACK routines \cite Federer2010
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Rate of change in state variables due to snow balance calculations
//
void CmvSnowBalance::ColdContentBalance(const double             *state_vars,
                                        const CHydroUnit *pHRU,
                                        const optStruct  &Options,
                                        const time_struct &tt,
                                        double     *rates) const
{
  const double CCFAC    = 0.3; //[MJ/m2-d-K]
  const double MELT_FAC = 1.5; //[MJ/m2-d-K],
  const double LAIMLT   = 0.2;
  const double SAIMLT   = 0.5;
  const double SWI = pModel->GetGlobalParams()->GetParams()->snow_SWI;

  double Ta         = pHRU->GetForcingFunctions()->temp_daily_ave;
  double day_length = pHRU->GetForcingFunctions()->day_length;

  double CC        =state_vars[iFrom[0]];//cold content MJ-mm/kg
  double SL        =state_vars[iFrom[1]];//liquid snow, mm
  double S         =state_vars[iFrom[3]];//snow, SWE mm

  double LAI       =pHRU->GetVegVarProps()->LAI;
  double SAI       =pHRU->GetVegVarProps()->SAI;
  double tstep     =Options.timestep;

  double incoming_snow_en;//incoming energy flux
  double pot_melt;        //incoming energy flux, in terms of water flux rate [mm/d]
  double CC_air;          //air temp, expressed in terms of cold content [MJ-mm/kg]
  double liq_snow_cap;    //maximum liquid snow capacity [mm]
  double rainthru=0.0;    //TMP DEBUG
  double Tsnow;           //[C]
  double instant_refreeze;
  double refreeze,cooling,warming,melt_to_liq,melt_to_SW,shrinking_poro,eq_en_avail;

  if (S==0){return;}

  //calculate instant equilibrium of CC and liquid water
  if ((CC>0) && (SL>0))
  {
    instant_refreeze=threshMin(SL/tstep,CC/LH_FUSION/tstep,0.0);//mm/day
    rates[1]+=instant_refreeze;            //SnowLiq->Snow
    rates[0]-=instant_refreeze*LH_FUSION;  //ColdCon->LostEn
    SL-=rates[1]*tstep;
  }

  Tsnow=FREEZING_TEMP-CC/S/SPH_ICE; //JRC: problematic if S is small
  //[C]=              [MJ-mm/kg]*[kg-K/MJ]*[1/mm]

  if (Ta <= FREEZING_TEMP)
  { //  snow warms or cools proportional to snow-air temperature difference
    incoming_snow_en=CCFAC*2.0*day_length*(Ta-Tsnow);
    //[MJ/m2/d]     =[MJ/m2-d-K]         *[K]
  }
  else
  { //energy input proportional to TA, modified by cover and daylength
    incoming_snow_en=MELT_FAC*2.0*day_length*(Ta-FREEZING_TEMP)*exp(-SAIMLT*SAI)*exp(-LAIMLT*LAI);
    //[MJ/m2/d]     =[MJ/m2-d-K]            *[K]
  }
  //incoming_snow_en*=30; //why is this necc to be consistent with hybrid??

  incoming_snow_en+=rainthru*threshPositive(Ta-FREEZING_TEMP)*SPH_WATER*DENSITY_WATER/MM_PER_METER;
  //[MJ/m2/d]      =[mm/d]  *[K]                             *[MJ/kg/K] *[kg/m3]      *[m/mm]

  pot_melt    =incoming_snow_en*(1.0/DENSITY_WATER/LH_FUSION*MM_PER_METER);
  //[mm/d]    =[MJ/m2/d]       *[m3/kg]      *[kg/MJ]   *[mm/m]

  CC_air      =(FREEZING_TEMP-Ta)*SPH_ICE*S;
  //MJ-mm/kg  =[K]               *[MJ/kg/K]*[mm]

  liq_snow_cap=CalculateSnowLiquidCapacity(S, 0.0, pModel);

  if (pot_melt<=0) //negative energy balance - snowpack cooling
  {
    refreeze=threshMin(SL/tstep,-pot_melt,0.0);             //cold used to refreeze liquid
    rates[1]+=refreeze;//SnowLiq->Snow [mm/d]
    pot_melt+=refreeze;

    cooling=threshMin(-pot_melt,0.0 ,0.0)*LH_FUSION;        //if any cold remaining, increases cold content
    cooling=threshMin(cooling,(CC_air-CC)/tstep,0.0);       //dont exceed maximum cold content (minimum temp Ta)
    rates[2]+=cooling;//Atm->CC [MJ-mm/kg-d]
  }
  else         //positive energy balance - snowpack warming
  {
    if ((pot_melt*LH_FUSION<CC/tstep) || (Ta<FREEZING_TEMP))
    {
      //energy merely reduces cold content, but not below equiv. air temp
      warming=pot_melt*LH_FUSION;
      warming=threshMin(warming,-(CC-CC_air)/tstep,0.0);//is negative sign correct here?

      rates[0]-=warming;//ColdCon<-LostEn
    }
    else
    {
      eq_en_avail=pot_melt;

      warming=CC/tstep;
      rates[0]+=warming; //ColdCon->LostEn//it's warm: all cold content lost
      eq_en_avail-=warming/LH_FUSION;

      melt_to_liq=threshMin(eq_en_avail,(liq_snow_cap-SL)/tstep,0.0);
      rates[1]-=melt_to_liq;//SnowLiq<-Snow (thus negative)
      eq_en_avail-=melt_to_liq;

      melt_to_SW=threshMin(eq_en_avail,S/tstep-melt_to_liq,0.0);
      rates[3]+=melt_to_SW;//snow->surface water

      shrinking_poro=melt_to_SW*SWI;
      rates[4]+=shrinking_poro;//snow_liq->surface water
    }
  }
  /*cout<<"iFrom[0]: "<<iFrom[0]<<endl;
    cout.precision(3);
    cout<<"Snow: "<<S<<endl;
    cout<<"Snow Liq:"<<SL<<endl;
    cout<<"Cold Content:"<<CC<<endl;
    cout<<"Day Length: "<<day_length<<endl;
    cout<<"Air Temp: "<<Ta<<endl;
    cout<<"Snow Temp: "<<Tsnow<<endl;
    cout<<"liq_snow_cap: "<<liq_snow_cap<<endl;
    cout<<"incoming_snow_en: "<<incoming_snow_en<<"[MJ/m2/d]"<<endl;
    cout<<"pot_melt: "<<pot_melt<<" mm/d"<<endl;
    if (pot_melt<0.0){cout<<"COOLING"<<endl;}
    else             {cout<<"WARMING"<<endl;}
    cout<<"CC_air: "<<CC_air<<endl;
    cout <<"rates[0]:"<<rates[0]<<endl;
    cout <<"rates[1]:"<<rates[1]<<endl;
    cout <<"rates[2]:"<<rates[2]<<endl;
    cout <<"rates[3]:"<<rates[3]<<endl;
    cout <<"rates[4]:"<<rates[4]<<endl<<endl;
    ExitGracefullyIf(S-(rates[4]-rates[1])*tstep<0,
    "Snowmelt Bad",RUNTIME_ERR);*/
}

//////////////////////////////////////////////////////////////////
/// \brief Balances cold content from snow melt/refreeze
///
/// \note Precipitation should come before snow balance for this to work as intended.
///    written by N. Sgro, 2016
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Rate of change in state variables due to snow balance calculations
//
void CmvSnowBalance::TwoLayerBalance(const double   *state_vars,
                                     const CHydroUnit *pHRU,
                                     const optStruct     &Options,
                                     const time_struct &tt,
                                     double     *rates) const
{
  //two layer snow mass balance based upon that of Wigmosta et al.,A distributed hydrology-vegetation model for complex terrain,  1994

  //initialize storage varaibles
  //------------------------------------------------------------------------
  double SweSurf; // Water equivalent of snow surface layer [mm]
  double SwePack; // Water equivalent of snow pack layer [mm]
  double newSnow  = state_vars[iFrom[0]]; //NEW SNOW [mm]
  double rainthru = state_vars[iFrom[1]]; //PONDED_WATER [mm]
  double Swe      = state_vars[iFrom[2]]; //total snow water equivalent [mm]
  double SlwcSurf = state_vars[iFrom[3]]; //liquid snow in snow surface layer [mm]
  double SlwcPack = state_vars[iFrom[4]]; //liquid snow in snow pack layer [mm]
  double CcSurf   = state_vars[iFrom[6]]; //cold content of snow surface layer [MJ/m2]
  double CcPack   = state_vars[iFrom[7]]; //cold content of snow pack layer [MJ/m2]
  double snowT    = state_vars[iFrom[8]]; //snow surface temperature [C]
  double cum_melt = state_vars[iFrom[9]]; //cumulative melt [mm]

  double gross_melt=0;

  if (Swe <= REAL_SMALL && newSnow <= REAL_SMALL) { return; } //no snow to balance

  //parameters
  //------------------------------------------------------------------------
  double MAXLIQ     = pModel->GetGlobalParams()->GetParams()->snow_SWI;         // maximum liquid water fraction of snow, dimensionless
  double MAXSWESURF = pModel->GetGlobalParams()->GetParams()->max_SWE_surface;  // maximum swe of surface layer of snowpack

  // forcings
  //------------------------------------------------------------------------
  double Mf = pHRU->GetForcingFunctions()->potential_melt*Options.timestep; //[mm]
  double Ta = pHRU->GetForcingFunctions()->temp_ave;

  //Rates
  //------------------------------------------------------------------------
  double meltSurf = 0.0;    // melt rate from snow surface layer
  double surfToPack = 0.0;  // liquid percolating from surface to pack
  double freezePack = 0.0;  // refreeze in the pack layer [mm]
  double SnowOutflow = 0.0; // outflow from the snow

  //INTERNAL
  //------------------------------------------------------------------------
  double ccSnowFall;

  //reset cumulative melt to zero in October
  //------------------------------------------------------------------------
  /*if ((pHRU->GetCentroid().latitude < 0) && (tt.month == 4) && (tt.day_of_month == 1)) {//southern hemisphere (April)
    cum_melt = 0.0;
  }
  else if ((tt.month == 10) && (tt.day_of_month == 1)) {//northern hemisphere
    cum_melt = 0.0;
  }*/

  //reconstruct two layer snowpack
  //------------------------------------------------------------------------
  if ( Swe > (MAXSWESURF-SlwcSurf))
  {
    SweSurf = MAXSWESURF-SlwcSurf;
  }
  else
  {
    SweSurf = Swe;
  }
  //SweSurf=min(MAXSWESURF-SlwcSurf,Swe);

  SwePack = Swe - SweSurf;

  //add snow throughfall to snowPack and compute its cold content
  //this also initializes the cold content with the first snow fall
  //------------------------------------------------------------------------
  ccSnowFall = HCP_ICE / MM_PER_METER * newSnow * max(-Ta, 0.0);

  //[MJ/m2] = [MJ/m3/K]  * [m/mm]      *  [mm]  *  [K]

  //distribute fresh snowfall
  //------------------------------------------------------------------------
  if ((newSnow + SweSurf+ SlwcSurf) < MAXSWESURF)// suface layer is below max
  {
    SweSurf += newSnow;
    CcSurf += ccSnowFall;
  }
  else
  {//snow surface layer is at max
    double deltaswePack = (SweSurf + SlwcSurf + newSnow) - MAXSWESURF; //JRC: If this is delta, shouldn't we subtract SwePack?
    double deltaccPack=0;
    if (SweSurf > 0) {
      deltaccPack = (deltaswePack/SweSurf) * CcSurf;//give pack the proportional coldconent of surface layer
    }
    SwePack += deltaswePack;
    SweSurf = MAXSWESURF - SlwcSurf;
    CcPack += deltaccPack;
    CcSurf += ccSnowFall - deltaccPack;
  }

  //add rainfall to surface snow liquid
  //------------------------------------------------------------------------
  SlwcSurf += rainthru;

  //add CC to Mf to get total energy available for melt/freeze
  // converts cold content to equivalent melt amount
  //------------------------------------------------------------------------
  Mf -= CcSurf / LH_FUSION / DENSITY_ICE*MM_PER_METER; // [MJ/m2] / [MJ/kg] /[kg/m3]*[mm/m]= [mm]

  //snowpack cooling or warming
  //engergy exchange occurs only between snow surface layer and atmosphere
  //snow surface layer and snow pack layer only exchange liquid water
  //------------------------------------------------------------------------
  //snowpack cooling or warming
  //------------------------------------------------------------------------
  if (Mf <= 0.0)  //   snow cooling
  {
    double posMf = -Mf;
    if (posMf < SlwcSurf)//refreeze part of liquid water
    {
      meltSurf += Mf;   // note that Mf is negative
      CcSurf = 0.0;
    }
    else// refreeze all liquid water and increase cold content
    {
      posMf    -= SlwcSurf;
      meltSurf -= SlwcSurf;
      CcSurf = posMf * LH_FUSION  * DENSITY_WATER / MM_PER_METER; // leftover energy is new CC ; this has to be water density since posMF is mmSWE

      if (SweSurf < 50.0) //dont cool below gruTa when Swe is below 50 mm (need this for stable run, otherwise CC and Tsnow start to fluctuate...)
      {
        CcSurf = min(CcSurf, -Ta * (SweSurf/MM_PER_METER) * HCP_ICE);
        CcSurf = max(CcSurf, 0.0);
      }
    }
  }
  else // snow warming (_Mf > 0)
  {
    CcSurf = 0.0; // can't have CC when snow is melting

    if (SweSurf < Mf) // All snow melts
    {
      meltSurf += SweSurf;
      gross_melt+=SweSurf;
    }
    else
    {
      meltSurf  += Mf;
      gross_melt +=Mf;
    }
  }
  SweSurf -= meltSurf; //new SWE after melt/freeze
  SlwcSurf+= meltSurf; //new Slwc after melt/freeze

  //Calculate how much liquid stays in surface layer, and pass the rest to the pack layer
  //------------------------------------------------------------------------
  if (SlwcSurf >= MAXLIQ * SweSurf)
  {
    surfToPack = SlwcSurf - MAXLIQ * SweSurf;
  }
  else
  {
    surfToPack = 0;
  }

  //snow pack layer
  //------------------------------------------------------------------------
  SlwcPack += surfToPack;//add the surface layer outflow to pack liquid water

  //re-freeze liquid water in snowpack layer
  //------------------------------------------------------------------------
  double ccLwcPack_eq=SlwcPack * LH_FUSION* DENSITY_WATER / MM_PER_METER; //[MJ/m2]
  if (CcPack > ccLwcPack_eq) //snow pack layer cold content not fully satisfied -freeze all water
  {
    CcPack     -= ccLwcPack_eq;
    freezePack += SlwcPack;
  }
  else //freeze only part of water
  {
    freezePack += CcPack / LH_FUSION / DENSITY_WATER * MM_PER_METER;
    CcPack = 0.0;
  }

  //update snow pack layer liquid water content
  //------------------------------------------------------------------------
  SwePack += freezePack; //new SWE after melt/freeze
  SlwcPack-= freezePack; //new Slwc after melt/freeze
  if (SlwcPack > SwePack * MAXLIQ)
  {
    SnowOutflow = SlwcPack - SwePack * MAXLIQ;
  }
  else
  {
    SnowOutflow = 0;
  }

  //Surface snow Temp - assumes an isothermal snow surface layer
  snowT = -CcSurf / (HCP_ICE* max(1.0, SweSurf)/MM_PER_METER);
  snowT *= 0.2; // NS: There's no physical basis for this, but snowT was always way too low.
  //            //     This correction puts snowT close to the air temperature.
  //            // \todo [bug] - identify source of this issue.

  cum_melt += gross_melt;

  rates[0] = newSnow     / Options.timestep;  // SNOWFALL -> SNOW [mm/d]
  rates[1] = rainthru    / Options.timestep;  // PONDEDWATER -> SNOW_LIQ surface [mm/d]
  rates[2] = meltSurf    / Options.timestep;  // SNOW -> SNOW_LIQ surface [mm/d]
  rates[3] = surfToPack  / Options.timestep;  // SNOW_LIQ surface -> SNOW_LIQ pack [mm/d]
  rates[4] = freezePack  / Options.timestep;  // SNOW_LIQ pack -> SNOW [mm/d]
  rates[5] = SnowOutflow / Options.timestep;  // SNOW_LIQ pack -> ponded water [mm/d]
  rates[6] = (CcSurf   - state_vars[iFrom[6]]) / Options.timestep; //CC pack layer  [mJ/m2/d]
  rates[7] = (CcPack   - state_vars[iFrom[7]]) / Options.timestep; //CC pack layer  [mJ/m2/d]
  rates[8] = (snowT    - state_vars[iFrom[8]]) / Options.timestep; //SNOW_TEMP   [C/d]
  rates[9] = (cum_melt - state_vars[iFrom[9]]) / Options.timestep; //Cumulative Melt [mm/d]
}

//////////////////////////////////////////////////////////////////
/// \brief GAWSER Snowmelt model
/// from object-GAWSER manual; coded by Konhee Lee, 2016
/// \note Precipitation should come before snow balance for this to work as intended.
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Rate of change in state variables due to snow balance calculations
//
void CmvSnowBalance::GawserBalance( const double      *state_vars,
                                    const CHydroUnit  *pHRU,
                                    const optStruct       &Options,
                                    const time_struct &tt,
                                    double            *rates) const
{
  double tstep = Options.timestep;         // Timestep [d]

  //initialize storage varaibles
  //------------------------------------------------------------------------
  double SWC        = state_vars[iFrom[0]];  // Solid Water Content
  double rainthru   = state_vars[iFrom[1]];  // PONDEDWATER [mm]
  double snow_liq   = state_vars[iFrom[2]];  // Liquid Water Content
  double snow_depth = state_vars[iFrom[3]];  // Snow Depth [mm]
  double newSnow    = state_vars[iFrom[6]];  // New Snow [mm]

  // Constants
  double KF =  pHRU->GetSurfaceProps()->refreeze_factor;  // refreeze factor [mm/d-degC]
  double KM =  pHRU->GetSurfaceProps()->melt_factor; // melt factor [mm/d-degC] //~5.04
  double SWI = pModel->GetGlobalParams()->GetParams()->snow_SWI; // Maximum fraction of pore space in snowpack for liquid snow

  double RHOICE = DENSITY_ICE/DENSITY_WATER;  // Relative density of ice
  double MRHO   = 0.35;                       // Maximum dry density for snowpack /// \todo [funct] - enable support of user-specified MRHO,
  double a = 0.1;                             // coefficient [1/degC]
  double b = 96;                              // coefficient [h]
  double NEWDEN = 0.1;                        // New snow density

  double TEMPs = pHRU->GetForcingFunctions()->temp_ave;
  double TBAS = 0;                         // Temperature of snowpack
  // Unclear in how its calculated
  // May want to use more intricate method after testing is completed

  // Calculate potential refreeze and melt
  double MELTP = 0;                        // Potential Melt [mm/h]
  double REFRZP = 0;                       // Potential Refreeze [mm/h]

  if (TEMPs < TBAS) {REFRZP = tstep*KF*(TBAS-TEMPs);}
  if (TEMPs > TBAS) {MELTP  = tstep*KM*(TEMPs-TBAS);}

  // Density Calculations
  double RHO = 0;                            // Dry density of snowpack
  if (SWC > 0) { RHO = min(MRHO, SWC / snow_depth); }

  double POR = 1- RHO / RHOICE;              // porosity of snowpack

  double snow_liqAP = POR * snow_depth * SWI;           // liquid water holding capactiy of snowpack [mm]

  double RAIN_IN_1 = 0;                      // Rain to snow_liq
  double MELT_1 = 0;                         // Melt to snow_liq
  double MELT_2 = 0;                         // Melt to Runoff
  double EXCESS_snow_liq = 0;                     // Excess snow_liq to runoff

  // Perform refreeze
  double REFRZ = 0;
  if (TBAS > TEMPs) { REFRZ = min(REFRZP, snow_liq); }

  // Update state variables after refreeze
  double SWC2 = SWC + REFRZ;
  double RHO2 = min(MRHO, SWC2 / snow_depth);
  double snow_liq2 = snow_liq - REFRZ;

  // Calculate Rainfall, Melt, and EXCESS snow_liq based on snow_liqAP
  if (snow_liq2 >= snow_liqAP)
  {
    EXCESS_snow_liq = snow_liq2 - snow_liqAP;
    MELT_2 = min(MELTP, SWC);
  }
  else {
    MELT_1 = min(MELTP,SWC);
    RAIN_IN_1 = rainthru;
  }

  // GAWSER's density based compaction calculation and snow_depth calculations
  double KC = b * exp(-1*a*TEMPs);                                           // Compaction time constant (h)
  double TRHO = (RHO2 * MRHO) / (RHO2 + (MRHO - RHO2)*exp(-1*tstep*24/KC));  // Delayed density

  double compaction = 0;                                                     // Compaction [mm]
  if (SWC2 > 0) { compaction = snow_depth - SWC2 / TRHO; }

  double snowmelt = 0;                                                       // Decrease in snow_depth from MELT [mm]
  if (SWC2 > 0) { snowmelt = (MELT_1 + MELT_2) / TRHO; }

  double SUBLIM = 0;                                                         // Sublimation [mm]

  if(TEMPs < 0){SUBLIM = 2.4 * tstep;} //2.4 mm/d
  if (SWC2 - MELT_1 - MELT_2 < SUBLIM) { SUBLIM = max(min(SUBLIM, SWC2 - MELT_1 - MELT_2),0.0); }

  double SUBLIM_snow_depth = 0;                                              // Decrease in snow_depth from Sublimation [mm]
  if (SWC2 > 0) { SUBLIM_snow_depth = SUBLIM / TRHO; }


  rates[0] = MELT_1 / tstep;                      // SNOW -> snow_liq
  rates[1] = RAIN_IN_1 / tstep;                   // RAIN -> snow_liq (IF snow_liq > snow_liqAP, RAIN_2 will remain in Ponded)
  rates[2] = REFRZ / tstep;                       // snow_liq  -> SNOW
  rates[3] = -compaction / tstep;                 // change in snowdepth from compaction
  rates[4] = -(snowmelt + SUBLIM_snow_depth) / tstep;   // change in snowdepth from snowmelt and sublimation
  rates[5] = newSnow / NEWDEN / tstep;            // change in snowdepth from snowfall
  rates[6] = newSnow / tstep;                     // SNOWFALL -> SNOW
  rates[7] = SUBLIM / tstep;                      // change in SNOW due to sublimation (physically makes no sense in GAWSER)
  rates[8] = EXCESS_snow_liq / tstep;             // snow_liq overflow to ponded
  rates[9] = MELT_2 / tstep;                      // Snowmelt going to ponded as snow_liq > snow_liqAP

}
//////////////////////////////////////////////////////////////////
/// \brief CRHM EBSM Snowmelt model
/// From ClassEBSM in CRHM
/// The model was first presented by Marks (1988), described conceptually by
/// Marks,et. al(1992) and Marks and Dozier(1992),and then described in great detail by Marks,et. al
/// (1997).
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Rate of change in state variables due to snow balance calculations
//
void CmvSnowBalance::CRHMSnowBalance(const double *state_vars,
                                     const CHydroUnit *pHRU,
                                     const optStruct &Options,
                                     const time_struct &tt,
                                     double *rates) const
{

  int iSnowLiq=pModel->GetStateVarIndex(SNOW_LIQ);
  int iSWE    =pModel->GetStateVarIndex(SNOW);
  int iSnowEn =pModel->GetStateVarIndex(COLD_CONTENT); //assumes [MJ/m2]

  double snow_liq    = state_vars[iSnowLiq];
  double SWE         = state_vars[iSWE];
  double snow_energy =-state_vars[iSnowEn]; //cold content [MJ/m2] is positive
  double snow_energy_old=snow_energy;

  if(SWE <= 0.0) { return; }

  double pot_melt =pHRU->GetForcingFunctions()->potential_melt;
  double T_min    =pHRU->GetForcingFunctions()->temp_daily_min;

  pot_melt*=LH_FUSION*DENSITY_WATER/MM_PER_METER*Options.timestep;  //[mm/d]->[MJ/m2]
  double snoliq_max = CalculateSnowLiquidCapacity(SWE, SWE/MAX_SNOW_DENS, pModel);
  double t_minus    = min(T_min,0.0);
  double Umin       = SWE*(2.115+0.00779*t_minus)*t_minus;//1000.0; //[MJ/m2] minimum snow energy (negative), documentation?? JRC- 1000 factor from CRHM unknown and leading to unreasonably small values of Umin.
  double refreeze   = 0.0;;
  double snowmelt0(0.0),snowmelt1(0.0),snowmelt2(0.0);


  //double Tmin=Umin/HCP_WATER/DENSITY_WATER/(SWE/MM_PER_METER); //Min cold content -> min snowpack mean temperature
  //cout<<Umin<<" "<<Tmin<<" "<<T_min<<endl;;
  snow_energy += pot_melt;
  snow_energy = max(Umin,snow_energy);
  if(snow_energy > 0.0)
  {
    double B=0.95; //thermal quality
    double melt = snow_energy/LH_FUSION/DENSITY_WATER*MM_PER_METER/B;   // [mm] = [MJ/m2]/[MJ/kg]/[kg/m3]*[mm/m]
    if(melt > (snoliq_max-snow_liq)) //ripening
    {
      if(SWE > melt) { //only snow melts, fills pores and overflows
        snowmelt0=snoliq_max-snow_liq;        //snow->snowliq
        snowmelt1=melt -(snoliq_max-snow_liq);//snow->surface w
        snowmelt2=0.0;                        //snow_liq->surface w
      }
      else { //everything melts and overflows
        snowmelt0=0;
        snowmelt1=SWE; //snow->sw
        snowmelt2=snow_liq;//snow_liq->sw
      }
    }
    else { // melt retained in snowpack
      snowmelt0=melt;
      snowmelt1=0;
      snowmelt2=0;
    }
    snow_energy = 0.0; //all energy used
  }
  else { // no melt - refreeze - reduce snow_energy accordingly
    snowmelt0=snowmelt1=snowmelt2=0;

    refreeze = -snow_energy/LH_FUSION/DENSITY_WATER*MM_PER_METER;
    if(snow_liq > refreeze) { //partial refreeze
      snow_energy   = 0.0;
      //refreeze=refreeze;
    }
    else { //total refreeze /cooling
      snow_energy  +=snow_liq*LH_FUSION*DENSITY_WATER/MM_PER_METER;
      refreeze=snow_liq;
    }
  } // no melt

  rates[0]=snowmelt0; //SNOW->LIQ_SNOW [mm]
  rates[1]=snowmelt1; //SNOW->PONDED_WATER [mm]
  rates[2]=snowmelt2; //SNOW_LIQ->PONDED_WATER [mm]
  rates[3]=refreeze;  //LIQ_SNOW->SNOW [mm]
  rates[4]=-(snow_energy-snow_energy_old)/Options.timestep; // COLD_CONTENT->COLD_CONTENT [MJ/m2]

  return;
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \note All constraints contained in routine itself- this routine is blank
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Rate of change in state variables due to snow balance calculations
//
void  CmvSnowBalance::ApplyConstraints( const double             *state_vars,
                                        const CHydroUnit *pHRU,
                                        const optStruct  &Options,
                                        const time_struct &tt,
                                        double     *rates) const
{

  if (type==SNOBAL_SIMPLE_MELT)
  {
    if (rates[0]<0.0){rates[0]=0.0;}//positivity constraint

    //cant remove more than is there
    rates[0]=threshMin(rates[0],max(state_vars[iFrom[0]]/Options.timestep,0.0),0.0);
  }
  else if (type == SNOBAL_UBCWM)
  {
    //snow cover should never be negative
    // subdaily correction sometimes puts it out of this range
    rates[4] = threshMax(rates[4], -state_vars[iFrom[4]] / Options.timestep, 0.0);
    rates[4] = threshMin(rates[4], (1-state_vars[iFrom[4]]) / Options.timestep, 0.0);
  }

  return;//most constraints contained in routine itself
}
