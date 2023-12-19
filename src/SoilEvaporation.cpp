/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2022 the Raven Development Team
  ----------------------------------------------------------------
  Soil Evaporation
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "SoilWaterMovers.h"
#include "Model.h"

/*****************************************************************
   Soil Evaporation Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the soil evaporation constructor
/// \param se_type [in] Model of soil evaporation selected
//
CmvSoilEvap::CmvSoilEvap(soilevap_type se_type,
                         CModelABC     *pModel)
  :CHydroProcessABC(SOIL_EVAPORATION, pModel)
{
  int iAtmos;
  iAtmos  =pModel->GetStateVarIndex(ATMOSPHERE);
  type = se_type;

  if(type==SOILEVAP_GAWSER)
  {
    CHydroProcessABC::DynamicSpecifyConnections(4);
    ExitGracefullyIf(pModel->GetNumSoilLayers()<2,
      "SOILEVAP_GAWSER algorithm requires at least 2 soil layers to operate. Please use a different :SoilModel or replace this evaporation algorithm.",BAD_DATA);

    iFrom[0]=pModel->GetStateVarIndex(SOIL,0);     iTo[0]=iAtmos;
    iFrom[1]=pModel->GetStateVarIndex(SOIL,1);     iTo[1]=iAtmos;
    iFrom[2]=pModel->GetStateVarIndex(DEPRESSION); iTo[2]=iAtmos;
    iFrom[3]=pModel->GetStateVarIndex(AET);        iTo[3]=iFrom[3];
  }
  else if((type==SOILEVAP_TOPMODEL) ||
          (type==SOILEVAP_VIC)    ||
          (type==SOILEVAP_HBV)    ||
          (type==SOILEVAP_UBC)    ||
          (type==SOILEVAP_PDM)    ||
          (type==SOILEVAP_CHU)    ||
          (type==SOILEVAP_GR4J)   ||
          (type==SOILEVAP_LINEAR) ||
          (type==SOILEVAP_HYMOD2) ||
          (type==SOILEVAP_ALL))
  {
    CHydroProcessABC::DynamicSpecifyConnections(2);

    iFrom[0]=pModel->GetStateVarIndex(SOIL,0);     iTo[0]=iAtmos;
    iFrom[1]=pModel->GetStateVarIndex(AET);        iTo[1]=iFrom[1];
  }
  else if (type==SOILEVAP_HYPR)
  {
    CHydroProcessABC::DynamicSpecifyConnections(3);

    iFrom[0]=pModel->GetStateVarIndex(SOIL,0);     iTo[0]=iAtmos;
    iFrom[1]=pModel->GetStateVarIndex(DEPRESSION); iTo[1]=iAtmos;
    iFrom[2]=pModel->GetStateVarIndex(AET);        iTo[2]=iFrom[2];
  }
  else if((type==SOILEVAP_SEQUEN) ||
    (type == SOILEVAP_ROOT) || (type == SOILEVAP_ROOT_CONSTRAIN))
  {
    CHydroProcessABC::DynamicSpecifyConnections(3);
    ExitGracefullyIf(pModel->GetNumSoilLayers()<2,
      "This soil infiltration algorithm requires at least 2 soil layers to operate. Please use a different :SoilModel or replace this evaporation algorithm.",BAD_DATA);

    iFrom[0]=pModel->GetStateVarIndex(SOIL,0);     iTo[0]=iAtmos;
    iFrom[1]=pModel->GetStateVarIndex(SOIL,1);     iTo[1]=iAtmos;
    iFrom[2]=pModel->GetStateVarIndex(AET);        iTo[2]=iFrom[2];
  }
  else if((type==SOILEVAP_ROOTFRAC) ||
          (type==SOILEVAP_FEDERER))
  {
    CHydroProcessABC::DynamicSpecifyConnections(nSoilLayers+1);
    ExitGracefully("CmvSoilEvap::Constructor:SOILEVAP_FEDERER",STUB);
    for(int m=0;m<nSoilLayers;m++) {
      iFrom[m]=pModel->GetStateVarIndex(SOIL,m);     iTo[m]=iAtmos;
    }
    iFrom[nSoilLayers]=pModel->GetStateVarIndex(AET);   iTo[nSoilLayers]=iFrom[nSoilLayers];
  }
  else if(type==SOILEVAP_SACSMA)
  {
    CHydroProcessABC::DynamicSpecifyConnections(7);
    iFrom[0]=pModel->GetStateVarIndex(SOIL,0);     iTo[0]=iAtmos;
    iFrom[1]=pModel->GetStateVarIndex(SOIL,1);     iTo[1]=iAtmos;
    iFrom[2]=pModel->GetStateVarIndex(SOIL,1);     iTo[2]=pModel->GetStateVarIndex(SOIL,0);
    iFrom[3]=pModel->GetStateVarIndex(SOIL,2);     iTo[3]=iAtmos;
    iFrom[4]=pModel->GetStateVarIndex(SOIL,4);     iTo[4]=pModel->GetStateVarIndex(SOIL,2);
    iFrom[5]=pModel->GetStateVarIndex(SOIL,5);     iTo[5]=iAtmos;

    iFrom[6]=pModel->GetStateVarIndex(AET);   iTo[6]=iFrom[6];
  }
  else if(type==SOILEVAP_AWBM) {
    CHydroProcessABC::DynamicSpecifyConnections(4);
    iFrom[0]=pModel->GetStateVarIndex(SOIL,0);     iTo[0]=iAtmos;
    iFrom[1]=pModel->GetStateVarIndex(SOIL,1);     iTo[1]=iAtmos;
    iFrom[2]=pModel->GetStateVarIndex(SOIL,2);     iTo[2]=iAtmos;
    iFrom[3]=pModel->GetStateVarIndex(AET);   iTo[3]=iFrom[3];
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
//
CmvSoilEvap::~CmvSoilEvap()
{
}

//////////////////////////////////////////////////////////////////
/// \brief Initializes soil evaporation
//
void CmvSoilEvap::Initialize(){

}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for soil evaporation algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by soil evaporation algorithm (size of aP[] and aPC[])
//
void CmvSoilEvap::GetParticipatingParamList(string  *aP , class_type *aPC , int &nP) const
{
  if (type==SOILEVAP_VIC)
  {
    nP=5;
    aP[0]="VIC_ALPHA";      aPC[0]=CLASS_SOIL;
    aP[1]="VIC_ZMAX";       aPC[1]=CLASS_SOIL;
    aP[2]="VIC_ZMIN";       aPC[2]=CLASS_SOIL;
    aP[3]="VIC_EVAP_GAMMA"; aPC[3]=CLASS_SOIL;
    aP[4]="POROSITY";       aPC[4]=CLASS_SOIL;
  }
  else if (type==SOILEVAP_GAWSER)
  {
    nP=1;
    aP[0]="POROSITY";       aPC[0]=CLASS_SOIL;
  }
  else if (type==SOILEVAP_FEDERER)
  {
    nP=0;
  }
  else if (type==SOILEVAP_ROOTFRAC)
  {
    nP=2;
    aP[0]="POROSITY";       aPC[0]=CLASS_SOIL;
    aP[1]="REL_ROOTDEN";    aPC[1]=CLASS_VEGETATION;
  }
  else if ((type==SOILEVAP_TOPMODEL) || (type==SOILEVAP_SEQUEN) || (type==SOILEVAP_ROOT) || (type == SOILEVAP_ROOT_CONSTRAIN))
  {
    nP=3;
    aP[0]="POROSITY";          aPC[0]=CLASS_SOIL;
    aP[1]="FIELD_CAPACITY";    aPC[1]=CLASS_SOIL;
    aP[2]="SAT_WILT";          aPC[2]=CLASS_SOIL;
  }
  else if (type==SOILEVAP_HBV)
  {
    nP=4;
    aP[0]="POROSITY";          aPC[0]=CLASS_SOIL;
    aP[1]="FIELD_CAPACITY";    aPC[1]=CLASS_SOIL;
    aP[2]="SAT_WILT";          aPC[2]=CLASS_SOIL;
    aP[3]="FOREST_COVERAGE";   aPC[3]=CLASS_LANDUSE; //JRCFLAG
  }
  else if(type==SOILEVAP_HYPR)
  {
    nP=8;
    aP[0]="POROSITY";             aPC[0]=CLASS_SOIL;
    aP[1]="FIELD_CAPACITY";       aPC[1]=CLASS_SOIL;
    aP[2]="SAT_WILT";             aPC[2]=CLASS_SOIL;
    aP[3]="FOREST_COVERAGE";      aPC[3]=CLASS_LANDUSE;
    aP[4]="MAX_DEP_AREA_FRAC";    aPC[4]=CLASS_LANDUSE;
    aP[5]="PONDED_EXP";           aPC[5]=CLASS_LANDUSE;
    aP[6]="DEP_MAX";              aPC[6]=CLASS_LANDUSE;
    aP[7]="PDMROF_B";             aPC[7]=CLASS_LANDUSE;
  }
  else if (type==SOILEVAP_UBC)
  {
    nP=4;
    aP[0]="IMPERMEABLE_FRAC";   aPC[0]=CLASS_LANDUSE;
    aP[1]="UBC_EVAP_SOIL_DEF";  aPC[1]=CLASS_SOIL;
    aP[2]="UBC_INFIL_SOIL_DEF"; aPC[2]=CLASS_SOIL;
    aP[3]="POROSITY";           aPC[3]=CLASS_SOIL;
  }
  else if (type==SOILEVAP_CHU)
  {
    nP=1;
    aP[0]="CHU_MATURITY";      aPC[0]=CLASS_VEGETATION;
  }
  else if(type==SOILEVAP_PDM)
  {
    nP=2;
    aP[0]="PDM_B";             aPC[0]=CLASS_LANDUSE;
    aP[1]="POROSITY";          aPC[1]=CLASS_SOIL;
  }
  else if(type==SOILEVAP_HYMOD2)
  {
    nP=5;
    aP[0]="PDM_B";             aPC[0]=CLASS_LANDUSE;
    aP[1]="POROSITY";          aPC[1]=CLASS_SOIL;
    aP[2]="HYMOD2_G";          aPC[2]=CLASS_LANDUSE;
    aP[3]="HYMOD2_KMAX";       aPC[3]=CLASS_LANDUSE;
    aP[4]="HYMOD2_EXP";        aPC[4]=CLASS_LANDUSE;
  }
  else if (type==SOILEVAP_LINEAR)
  {
    nP=1;
    aP[0]="AET_COEFF";         aPC[0]=CLASS_LANDUSE;
  }
  else if((type==SOILEVAP_GR4J) || (type==SOILEVAP_ALL))
  {
    nP=0;
  }
  else if(type==SOILEVAP_SACSMA)
  {
    nP=3;
    aP[0]="MAX_SAT_AREA_FRAC";       aPC[0]=CLASS_LANDUSE;
    aP[1]="IMPERMEABLE_FRAC";        aPC[1]=CLASS_LANDUSE;
    aP[2]="UNAVAIL_FRAC";            aPC[2]=CLASS_SOIL;
  }
  else if(type==SOILEVAP_AWBM)
  {
    nP=2;
    aP[0]="AWBM_AREAFRAC1";        aPC[0]=CLASS_LANDUSE;
    aP[1]="AWBM_AREAFRAC2";        aPC[1]=CLASS_LANDUSE;
  }
  else
  {
    ExitGracefully("CmvSoilEvap::GetParticipatingParamList: undefined soil evaporation algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
///
/// \param se_type [in] Model of soil evaporation used
/// \param *aSV [out] Array of state variable types needed by soil evaporation algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by soil evaporation algorithm (size of aSV[] and aLev[] arrays)
//
void CmvSoilEvap::GetParticipatingStateVarList(soilevap_type se_type,sv_type *aSV, int *aLev, int &nSV)
{
  if (se_type==SOILEVAP_GAWSER)
  {
    nSV=4;
    aSV [0]=SOIL;  aSV  [1]=SOIL;  aSV [2]=ATMOSPHERE;    aSV [3]=DEPRESSION;
    aLev[0]=0;     aLev [1]=1;     aLev[2]=DOESNT_EXIST;  aLev[3]=DOESNT_EXIST;
  }
  else if ((se_type==SOILEVAP_TOPMODEL) || (se_type==SOILEVAP_VIC) || (se_type==SOILEVAP_HBV) || (se_type==SOILEVAP_LINEAR) || (se_type==SOILEVAP_ALL) || (se_type==SOILEVAP_PDM) || (se_type==SOILEVAP_HYMOD2))
  {
    nSV=2;
    aSV [0]=SOIL;  aSV [1]=ATMOSPHERE;
    aLev[0]=0;     aLev[1]=DOESNT_EXIST;
  }
  else if(se_type==SOILEVAP_HYPR)
  {
    nSV=3;
    aSV[0]=SOIL;  aSV[1]=ATMOSPHERE;      aSV[2]=DEPRESSION;
    aLev[0]=0;    aLev[1]=DOESNT_EXIST;  aLev[2]=DOESNT_EXIST;
  }
  else if (se_type==SOILEVAP_UBC)
  {
    nSV=3;
    aSV [0]=SOIL;  aSV [1]=ATMOSPHERE;   aSV [2]=SNOW;
    aLev[0]=0;     aLev[1]=DOESNT_EXIST; aLev[2]=DOESNT_EXIST;

  }
  else if (se_type==SOILEVAP_CHU)
  {
    nSV=3;
    aSV [0]=SOIL;  aSV [1]=ATMOSPHERE;     aSV[2]=CROP_HEAT_UNITS;
    aLev[0]=0;     aLev[1]=DOESNT_EXIST;  aLev[2]=DOESNT_EXIST;
  }
  else if ((se_type==SOILEVAP_SEQUEN) || (se_type==SOILEVAP_ROOT) || (se_type==SOILEVAP_ROOT_CONSTRAIN))
  {
    nSV=3;
    aSV [0]=SOIL;       aLev[0]=0;
    aSV [1]=SOIL;       aLev[1]=1;
    aSV [2]=ATMOSPHERE; aLev[2]=DOESNT_EXIST;
  }
  else if ((se_type==SOILEVAP_ROOTFRAC) || (se_type==SOILEVAP_FEDERER))
  {
    nSV=0; //multilayer/user-specified
  }
  else if (se_type==SOILEVAP_GR4J)
  {
    nSV=2;
    aSV [0]=SOIL;  aSV [1]=ATMOSPHERE;
    aLev[0]=0;     aLev[1]=DOESNT_EXIST;
  }
  else if(se_type==SOILEVAP_SACSMA)
  {
    nSV=7;
    for(int m=0;m<=5;m++) {
      aSV[m]=SOIL; aLev[m]=m;
    }
    aSV[6]=ATMOSPHERE; aLev[6]=DOESNT_EXIST;
  }
  else if(se_type==SOILEVAP_AWBM)
  {
    nSV=4;
    for(int m=0;m<3;m++) {
      aSV[m]=SOIL; aLev[m]=m;
    }
    aSV[3]=ATMOSPHERE; aLev[3]=DOESNT_EXIST;
  }
  nSV++;
  aSV[nSV-1]=AET;
  aLev[nSV-1]=DOESNT_EXIST;

}

//////////////////////////////////////////////////////////////////
/// \brief Returns rates of loss from set of soil layers to atmosphere due to evapotranspiration/transpiration[mm/day]
///
/// \param *state_vars [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of loss from "from" compartment [mm/day]
//
void CmvSoilEvap::GetRatesOfChange (const double      *state_vars,
                                    const CHydroUnit  *pHRU,
                                    const optStruct   &Options,
                                    const time_struct &tt,
                                    double            *rates) const
{

  if (pHRU->GetHRUType()!=HRU_STANDARD){return;}//Lake/Glacier case

  double PET,PETused(0.0);
  const soil_struct *pSoil;

  PET=pHRU->GetForcingFunctions()->PET;

  if (!Options.suppressCompetitiveET){
    //competitive ET - reduce PET by AET
    PET-=(state_vars[pModel->GetStateVarIndex(AET)]/Options.timestep);
    PET=max(PET,0.0);
  }

  //------------------------------------------------------------
  if (type==SOILEVAP_ROOTFRAC)
  {
    ///  from Desborough, 1997 \cite Desborough1997MWR
    double      root_frac[MAX_SOILLAYERS];
    double      cap      [MAX_SOILLAYERS];
    int         m,q;

    double      rootsum=0.0;

    for (m=0;m<nSoilLayers;m++)
    {
      cap      [m]=pHRU->GetSoilCapacity(pModel->GetStateVarIndex(SOIL,m));
      root_frac[m]=pHRU->GetVegVarProps()->rel_rootden;

      rootsum+=root_frac[m];
    }
    for (q=0;q<_nConnections-1;q++)
    {
      m=q;
      rates[q]=PET*(root_frac[m]/rootsum)*threshMin(1.0,state_vars[iFrom[q]]/cap[m],0.0);
      PETused+=rates[q];
    }
  }
  //------------------------------------------------------------
  else if (type==SOILEVAP_LINEAR) //linear function of saturation
  {
    double stor    = state_vars[iFrom[0]];//[mm]
    double alpha   =pHRU->GetSurfaceProps()->AET_coeff;

    rates[0]  = min(alpha*stor,PET);  //evaporation rate [mm/d]
    PETused=rates[0];
  }
  //------------------------------------------------------------
  else if (type==SOILEVAP_ALL)
  {
    rates[0]  = PET;  //evaporation rate [mm/d]
    PETused=rates[0];
  }
  //------------------------------------------------------------
  else if ((type==SOILEVAP_TOPMODEL) || (type==SOILEVAP_HBV)  || (type==SOILEVAP_HYPR))
  {
    //From HBV Model (Bergstrom,1995)
    double stor,tens_stor; //[mm]

    stor      = state_vars[iFrom[0]];
    tens_stor = pHRU->GetSoilTensionStorageCapacity(0);

    rates[0]  = PET * min(stor/tens_stor,1.0);  //evaporation rate [mm/d]

    //correction for snow in non-forested areas (not in HYPR)
    if (type==SOILEVAP_HBV)
    {
      int iSnow=pModel->GetStateVarIndex(SNOW);
      double Fc=pHRU->GetSurfaceProps()->forest_coverage;
      if ((iSnow!=DOESNT_EXIST) && (state_vars[iSnow]>REAL_SMALL)) {rates[0]=(Fc)*rates[0];}//+(1.0-Fc)*0.0; (implied)
    }
    PETused=rates[0];

    //SOILEVAP_HYPR from
    //Ahmed et al., Toward Simple Modeling Practices in the Complex Canadian Prairie Watersheds,
    //Journal of Hydrologic Engineering 25(6), 04020024, doi:10.1061/(ASCE)HE.1943-5584.0001922, 2020
    if (type==SOILEVAP_HYPR)
    {
      int iDep=pModel->GetStateVarIndex(DEPRESSION);
      double maxPondedAreaFrac=pHRU->GetSurfaceProps()->max_dep_area_frac;
      double n                =pHRU->GetSurfaceProps()->ponded_exp;
      double dep_max          =pHRU->GetSurfaceProps()->dep_max;

      double area_ponded=maxPondedAreaFrac*pow(min(state_vars[iDep]/dep_max,1.0),n);

      rates[0]=(1.0-area_ponded)*rates[0]; //correct AET              //SOIL->ATMOS
      rates[1]=(    area_ponded)*pHRU->GetForcingFunctions()->OW_PET; //DEPRESSION->ATMOS

      PETused=(rates[0]+rates[1]);
    }

  }
  //------------------------------------------------------------
  else if(type==SOILEVAP_PDM)
  {
    double stor    =state_vars[iFrom[0]];
    double max_stor=pHRU->GetSoilCapacity(0);

    double b=pHRU->GetSurfaceProps()->PDM_b;

    double c_max =(b+1)*max_stor;
    double c_star=c_max*(1.0-pow(1.0-(stor/max_stor),1.0/(b+1.0)));

    rates[0] = min(c_star,(c_star/c_max)*PET);     //AET: SOIL->ATMOS

    PETused=rates[0];
  }
  //------------------------------------------------------------
  else if(type==SOILEVAP_HYMOD2)
  {
    //from Roy et al.  (2017), Using satellite-based evapotranspiration estimates to improve the structure of a simple conceptual rainfall-runoff model, HESS, 21(2), 879�896, doi:10.5194/hess-21-879-2017

    double stor    =state_vars[iFrom[0]];
    double max_stor=pHRU->GetSoilCapacity(0);

    double b   =pHRU->GetSurfaceProps()->PDM_b;
    double gp  =pHRU->GetSurfaceProps()->HYMOD2_G;    // [0..1] ET lower resistance parameter
    double Kmax=pHRU->GetSurfaceProps()->HYMOD2_Kmax; // [0..1] ET resistance parameter
    double ce  =pHRU->GetSurfaceProps()->HYMOD2_exp;  // [-] ET exponent parameter

    double c_max =(b+1)*max_stor;
    double c_star=c_max*(1.0-pow(1.0-(stor/max_stor),1.0/(b+1.0)));

    double K=Kmax*(gp+(1-gp)*pow(c_star/c_max,ce));
    rates[0] = K*min(c_star,PET);     //AET: SOIL->ATMOS

    PETused=rates[0];
  }
  //------------------------------------------------------------
  else if (type==SOILEVAP_CHU)
  { //Ontario Heat Crop Method: evaporation rated calculated using ratio of crop heat units to CHU maturity
    double CHU;

    CHU       = max(state_vars[pModel->GetStateVarIndex(CROP_HEAT_UNITS)],0.0);

    rates[0]  = PET*min(CHU/pHRU->GetVegetationProps()->CHU_maturity,1.0);  //evaporation rate [mm/d]
    PETused=rates[0];
  }
  //------------------------------------------------------------
  else if (type==SOILEVAP_UBC)
  { //From UBCWM Watershed Model (Quick, 1995)
    double P0AGEN=pHRU->GetSoilProps(0)->UBC_infil_soil_def; // [mm] - soil deficit at which effective impermeable fraction depletes to 0.1
    double P0EGEN=pHRU->GetSoilProps(0)->UBC_evap_soil_def;  // [mm] - soil deficit at which AET depletes to =0.1*PET
    double soil_deficit=max(pHRU->GetSoilCapacity(0)-state_vars[iFrom[0]],0.0);
    double Fimp        =pHRU->GetSurfaceProps()->impermeable_frac;

    //RFS Emulation - Calculate estimted soil deficit based on precip and AET from previous day
    double AET             =(PET)*pow(10.0,-soil_deficit/P0EGEN);
    double total_precip    =state_vars[pModel->GetStateVarIndex(PONDED_WATER,0)];
    double soil_deficit_est=max(soil_deficit - total_precip + AET*Options.timestep,0.0);

    //Calc relative impermeable fraction
    double b1=0.0;
    if (Fimp<1.0)
    {
      b1 = Fimp*pow(10.0, -soil_deficit_est / P0AGEN);
      // b1=0.0; //RFS: change in soil deficit not calculated using b1
    }
    //rates[0]=PET*pow(10.0,-soil_deficit/P0EGEN)*(1.0-b1); //actual ET
    rates[0]=(PET)*pow(10.0,-soil_deficit_est/P0EGEN)*(1.0-b1); //actual ET (RFS EMulation)
    PETused=rates[0];
  }
  //------------------------------------------------------------
  else if (type==SOILEVAP_VIC)
  {/// from (Woods et al 1992)
    double alpha,zmax,zmin,gamma2,Smax,Sat;
    double stor=state_vars[iFrom[0]];
    double stor_max;

    stor_max=pHRU->GetSoilCapacity(0);
    pSoil   =pHRU->GetSoilProps(0);

    alpha =pSoil->VIC_alpha;
    zmax  =pSoil->VIC_zmax;
    zmin  =pSoil->VIC_zmin;
    gamma2=pSoil->VIC_evap_gamma;

    Smax  =1.0/(alpha+1.0)*(alpha*zmax+zmin);
    Sat   =stor/stor_max;

    rates[0]=PET*(1.0-pow(1.0-Sat/Smax,gamma2));
    PETused=rates[0];
  }
  //------------------------------------------------------------------
  else if (type==SOILEVAP_GAWSER)
  { //from GAWSER Manual,1996 (Hinckley et al., 1996)
    double stor,stor2,dep_stor,max_stor1;
    double PETremain;

    stor     =state_vars[iFrom[0]];
    stor2    =state_vars[iFrom[1]];
    dep_stor =state_vars[iFrom[3]];
    max_stor1=pHRU->GetSoilCapacity(0);
    PETremain=PET;

    //first remove from depression storage
    rates[0]=min(PET,dep_stor/Options.timestep);
    PETremain-=rates[0];

    //then remove from top layer if storage>0.5 capacity
    rates[1]=min(PETremain,max(stor-0.5*max_stor1,0.0)/Options.timestep);
    PETremain-=rates[1];

    //then remove from both compartments equally
    double from_top=min(0.5*PETremain,stor/Options.timestep);
    rates[1]+=from_top;
    rates[2]=min(0.5*PETremain,stor2/Options.timestep);
    PETremain-=(from_top+rates[2]);

    //then just from bottom storage
    rates[2]=min(PETremain,stor2/Options.timestep);
    PETused=rates[0]+rates[1]+rates[2];
  }
  //------------------------------------------------------------
  else if (type==SOILEVAP_FEDERER)
  {
  ///  Adapted from Brook90 routine TBYLAYER based on model of Federer 1979, [A soil-plant-atmosphere model for transpiration and availability of soil water. Water Resour Res 15:555-562.]

    ExitGracefully("FedererSoilEvap::Not tested!",STUB);
    //FedererSoilEvap(PET,state_vars,pHRU,Options,tt,rates);
  }
  //------------------------------------------------------------------
  else if (type==SOILEVAP_SEQUEN)
  {
    double stor_u,stor_l;          //upper/lower soil layer storage
    double tens_stor_u,tens_stor_l;//maximum tension storage in soil layers [mm]

    stor_u     = max(state_vars[iFrom[0]],0.0);
    stor_l     = max(state_vars[iFrom[1]],0.0);
    tens_stor_u= pHRU->GetSoilTensionStorageCapacity(0);
    tens_stor_l= pHRU->GetSoilTensionStorageCapacity(1);

    rates[0]   = (PET           )*min(stor_u/tens_stor_u,1.0);  //upper layer evaporation rate [mm/d]
    rates[1]   = (PET - rates[0])*min(stor_l/tens_stor_l,1.0);  //lower layer evaporation rate [mm/d]
    //   if (rates[0]<-REAL_SMALL){ cout << stor_u << " " << tens_stor_u << " "<<PET<<endl; }
    PETused=rates[0]+rates[1];
  }
  //------------------------------------------------------------------
  else if (type==SOILEVAP_ROOT)
  {
    double stor_u,stor_l;                       //soil layer storage [mm]
    double tens_stor_u,tens_stor_l;             //maximum layer tension storage [mm]
    double rootfrac_u,rootfrac_l;               //relative root fraction soil layers [unitless]

    rootfrac_u  = 0.7;//pRootVar->rootfrac;
    rootfrac_l  = 1.0 - rootfrac_u;

    stor_u      = state_vars[iFrom[0]];
    tens_stor_u = pHRU->GetSoilTensionStorageCapacity(0);
    tens_stor_u=max(0.0001,tens_stor_u);

    rates[0] = PET * rootfrac_u * min(stor_u/tens_stor_u,1.0);  //upper layer evaporation rate [mm/d]

    stor_l      = state_vars[iFrom[1]];
    tens_stor_l = pHRU->GetSoilTensionStorageCapacity(1);
    tens_stor_l=max(0.0001,tens_stor_l);

    rates[1] = PET * rootfrac_l * min(stor_l/tens_stor_l,1.0);  //upper layer evaporation rate [mm/d]

    //cout<<rates[0]<<" "<<rates[1]<<" "<<PET<<" "<<stor_u<<" "<<" "<<tens_stor_u<<" "<<stor_l<<" "<<" "<<tens_stor_l<<endl;
    PETused=rates[0]+rates[1];
    //fix calculation of lower and upper - ITS WRONG SOMEWHERE (either ridiculously high or ridiculously low)!!!!
    //cout<<"PET: "<<PET<<"  root_frac: "<<rel_rootfrac_upper<<"  tension_stor: "<<tension_stor_upper<<"  max_ten: "<<max_tension_stor_upper<<endl;
    //cout<<" s_evap rate 0: "<<rates[0]<<"  s_evap rate 1: "<<rates[1]<<endl;
  }
  //------------------------------------------------------------------
  else if (type == SOILEVAP_ROOT_CONSTRAIN)
  {
    double stor_u, stor_l;              //soil layer storage [mm]
    double tens_stor_u, tens_stor_l;    //maximum layer tension storage [mm]
    double rootfrac_u, rootfrac_l;      //relative root fraction soil layers [unitless]
    double stor_wilt;

    rootfrac_u  = 0.7;//pRootVar->rootfrac;
    stor_u      = state_vars[iFrom[0]];
    tens_stor_u = pHRU->GetSoilTensionStorageCapacity(0);
    stor_wilt   = pHRU->GetSoilProps(0)->sat_wilt*pHRU->GetSoilCapacity(0);

    rates[0] = PET * rootfrac_u * min(stor_u / tens_stor_u, 1.0);  //upper layer evaporation rate [mm/d]

    if (rates[0] > max(stor_u - stor_wilt, 0.0))
    {
      rates[0] = max(stor_u - stor_wilt, 0.0);
    }

    rootfrac_l = 1.0 - rootfrac_u;
    stor_l = state_vars[iFrom[1]];
    tens_stor_l = pHRU->GetSoilTensionStorageCapacity(1);

    rates[1] = PET * rootfrac_l * 1;// min(stor_l / tens_stor_l, 1.0);  //upper layer evaporation rate [mm/d]
    PETused=rates[0]+rates[1];
  }
  //------------------------------------------------------------------
  /*else if (type==SOILEVAP_POWERLAW)
    { //similar to Watflood (n=1/2), except in watflood max_stor here could be between field capacity and max_stor
    double factor;
    double n;
    double stor,max_stor,wilting_pt;
    pSoil               =pHRU->GetSoilProps(0);

    n         =0.5;//pSoil->soilevap_n;
    max_stor  =pHRU->GetSoilCapacity(0);
    wilting_pt=(pSoil->sat_wilt*max_stor);

    factor=max(min((stor-wilting_pt)/(max_stor-wilting_pt),1.0),0.0);
    rates[0]=PET*pow(factor,n);
    PETused=rates[0];
    }*/
  //------------------------------------------------------------------
  else if (type==SOILEVAP_GR4J)
  { //from GR4J model (Perrin et al., 2003)
    double max_stor=pHRU->GetSoilCapacity(0);
    double stor    =state_vars[iFrom[0]];
    double sat     =stor/max_stor;

    double tmp=tanh(max(PET,0.0)/max_stor);

    rates[0]=stor*(2.0-sat)*tmp/(1.0+(1.0-sat)*tmp);
    PETused=rates[0];
  }
  //------------------------------------------------------------------
  else if(type==SOILEVAP_SACSMA)
  { //From Sacramento Soil Accounting Model
    double red;            // [mm] residual evap demand
    double e1,e2,e3,e5;    // [mm] different ETs from different regions

    double Adimp       =pHRU->GetSurfaceProps()->max_sat_area_frac; //ADIMP
    double Aimp        =pHRU->GetSurfaceProps()->impermeable_frac; //AIMP
    double rserv       =pHRU->GetSoilProps(3)->unavail_frac; //RSERV

    double Aperv = 1.0 - Adimp - Aimp;

    double uzt_stor_max =pHRU->GetSoilCapacity(0); //UZTM
    double uzf_stor_max =pHRU->GetSoilCapacity(1); //UZFM
    double lzt_stor_max =pHRU->GetSoilCapacity(2); //LZTM
    double lzfp_stor_max=pHRU->GetSoilCapacity(3); //LZFPM
    double lzfs_stor_max=pHRU->GetSoilCapacity(4); //LZFSM

    double uzt_stor  =state_vars[pModel->GetStateVarIndex(SOIL,0)]/Aperv;
    double uzf_stor  =state_vars[pModel->GetStateVarIndex(SOIL,1)]/Aperv;
    double lzt_stor  =state_vars[pModel->GetStateVarIndex(SOIL,2)]/Aperv;
    double lzfp_stor =state_vars[pModel->GetStateVarIndex(SOIL,3)]/Aperv;
    double lzfs_stor =state_vars[pModel->GetStateVarIndex(SOIL,4)]/Aperv;
    double adimc_stor=state_vars[pModel->GetStateVarIndex(SOIL,5)]/Adimp;
    if(Adimp==0) { adimc_stor=0; }

    e1 = min(PET * (uzt_stor/uzt_stor_max),uzt_stor);
    red = PET - e1;
    uzt_stor -= e1;
    rates[0]=e1*Aperv/Options.timestep; //UZT SOIL[0]->ATMOSPHERE

    e2=min(red,uzf_stor);
    red-=e2;
    uzf_stor-=e2;
    rates[1]=e2*Aperv/Options.timestep; //UZF SOIL[1]->ATMOSPHERE

    // equilibrate storage ratios in upper zone
    if((uzt_stor/uzt_stor_max) < (uzf_stor/uzf_stor_max))
    {
      double delta= uzt_stor_max*(uzt_stor + uzf_stor) / (uzt_stor_max + uzf_stor_max)-uzt_stor;
      //cout<<"uzf_stor "<<uzt_stor/uzt_stor_max<<" "<<uzf_stor/uzf_stor_max<<" "<<uzf_stor<<" "<<delta<<endl;
      uzt_stor += delta;
      uzf_stor -= delta;
      rates[2] = delta*Aperv/Options.timestep; //UZF SOIL[1]-> UZT SOIL[0]
    }

    e3 = min(red * (lzt_stor/(uzt_stor_max + lzt_stor_max)),lzt_stor);
    lzt_stor -= e3;
    rates[3]=e3*Aperv/Options.timestep;       //LZT SOIL[2]->ATMOSPHERE

    double ratlzt = lzt_stor/lzt_stor_max;
    double saved = rserv * (lzfp_stor_max + lzfs_stor_max);
    double ratlz = (lzt_stor+lzfp_stor+lzfs_stor-saved)/(lzt_stor_max+lzfp_stor_max+lzfs_stor_max-saved);
    if(ratlzt < ratlz)
    { // resupply lower zone tension water from lower zone free water if more water available there
      double del = min((ratlz - ratlzt) * lzt_stor_max,lzfs_stor);
      lzt_stor  += del;
      lzfs_stor -= del;
      //rates[4]=del*Aperv/Options.timestep;   //LZFS SOIL[4]->LZT SOIL[2]
    }

    // adjust adimc,additional impervious area storage, for evaporation
    e5 = min(e1 + (red+e2) * (adimc_stor-uzt_stor-e1) / (uzt_stor_max+lzt_stor_max),adimc_stor);
    adimc_stor -= e5;
    rates[5]=e5*Adimp/Options.timestep;     //ADIMC (SOIL[5])->ATMOS

    PETused=(e1+e2+e3)*Aperv+e5*Adimp;
  }
  //------------------------------------------------------------
  else if(type==SOILEVAP_AWBM)
  {//Boughton2004
    double a1 =pHRU->GetSurfaceProps()->AWBM_areafrac1;
    double a2 =pHRU->GetSurfaceProps()->AWBM_areafrac2;
    double a3=1.0-a1-a2;

    rates[0] =min(a1*PET,state_vars[iFrom[0]]/Options.timestep);
    rates[1] =min(a2*PET,state_vars[iFrom[1]]/Options.timestep);
    rates[2] =min(a3*PET,state_vars[iFrom[2]]/Options.timestep);

    PETused=rates[0]+rates[1]+rates[2];
  }
  else
  {
    ExitGracefully("CmvSoilEvaporation::GetRatesOfChange: undefined soil evaporation type",BAD_DATA);
  }//end soil_evap type select

  //updated used PET
  rates[_nConnections-1]=PETused;

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
void   CmvSoilEvap::ApplyConstraints( const double               *state_vars,
                                      const CHydroUnit *pHRU,
                                      const optStruct    &Options,
                                      const time_struct &tt,
                                      double     *rates) const
{
  if (pHRU->GetHRUType()!=HRU_STANDARD){return;}//Lake/Glacier case

  double corr=0,oldrate=0.0;
  for (int q=0;q<_nConnections-1;q++){
    oldrate=rates[q];
    //cant remove more than is there
    rates[q]=threshMin(rates[q],state_vars[iFrom[q]]/Options.timestep,0.0); //presumes these are all water storage
    corr+=oldrate-rates[q];
  }
  rates[_nConnections-1]-=corr;
}
