/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Model.h"
#include "Infiltration.h"



//////////////////////////////////////////////////////////////////
/// \brief Implementation of the infiltration process constructor
/// \param itype [in] Model of infiltration selected
//
CmvInfiltration::CmvInfiltration(infil_type itype,
                                 CModelABC *pModel)
  :CHydroProcessABC(INFILTRATION, pModel)
{
  type =itype;
  CHydroProcessABC::DynamicSpecifyConnections(2);
  //infiltration (ponded-->soil)
  iFrom[0]=pModel->GetStateVarIndex(PONDED_WATER);      iTo  [0]=pModel->GetStateVarIndex(SOIL,0);
  //runoff/remainder (ponded->surface water)
  iFrom[1]=pModel->GetStateVarIndex(PONDED_WATER);      iTo  [1]=pModel->GetStateVarIndex(SURFACE_WATER);

  if (itype==INF_GA_SIMPLE)
  {
    CHydroProcessABC::DynamicSpecifyConnections(4);

    iFrom[0]=pModel->GetStateVarIndex(PONDED_WATER);     iTo  [0]=pModel->GetStateVarIndex(SOIL,0);
    iFrom[1]=pModel->GetStateVarIndex(PONDED_WATER);     iTo  [1]=pModel->GetStateVarIndex(SURFACE_WATER);
    iFrom[2]=pModel->GetStateVarIndex(CUM_INFIL);        iTo  [2]=pModel->GetStateVarIndex(CUM_INFIL);
    iFrom[3]=pModel->GetStateVarIndex(GA_MOISTURE_INIT); iTo  [3]=pModel->GetStateVarIndex(GA_MOISTURE_INIT);
  }
  else if (type==INF_UBC){
    CHydroProcessABC::DynamicSpecifyConnections(5);
    ExitGracefullyIf(pModel->GetNumSoilLayers()<4,
                     "INF_UBCWM infiltration algorithm requires at least 4 soil layers to operate. Please use a different :SoilModel or replace this infiltration algorithm.",BAD_DATA);
    iFrom[0]=pModel->GetStateVarIndex(PONDED_WATER);    iTo  [0]=pModel->GetStateVarIndex(SOIL,0);
    iFrom[1]=pModel->GetStateVarIndex(PONDED_WATER);    iTo  [1]=pModel->GetStateVarIndex(SURFACE_WATER);
    iFrom[2]=pModel->GetStateVarIndex(PONDED_WATER);    iTo  [2]=pModel->GetStateVarIndex(SOIL,1);
    iFrom[3]=pModel->GetStateVarIndex(PONDED_WATER);    iTo  [3]=pModel->GetStateVarIndex(SOIL,2);
    iFrom[4]=pModel->GetStateVarIndex(PONDED_WATER);    iTo  [4]=pModel->GetStateVarIndex(SOIL,3);
  }
  else if(type==INF_HMETS) {
    CHydroProcessABC::DynamicSpecifyConnections(4);
    iFrom[0]=pModel->GetStateVarIndex(PONDED_WATER);    iTo[0]=pModel->GetStateVarIndex(SOIL,0);
    iFrom[1]=pModel->GetStateVarIndex(PONDED_WATER);    iTo[1]=pModel->GetStateVarIndex(SURFACE_WATER);
    iFrom[2]=pModel->GetStateVarIndex(PONDED_WATER);    iTo[2]=pModel->GetStateVarIndex(CONVOLUTION,0);
    iFrom[3]=pModel->GetStateVarIndex(PONDED_WATER);    iTo[3]=pModel->GetStateVarIndex(CONVOLUTION,1);
  }
  else if(type==INF_XINANXIANG) {
    CHydroProcessABC::DynamicSpecifyConnections(5);
    iFrom[0]=pModel->GetStateVarIndex(PONDED_WATER);    iTo[0]=pModel->GetStateVarIndex(SOIL,0);
    iFrom[1]=pModel->GetStateVarIndex(PONDED_WATER);    iTo[1]=pModel->GetStateVarIndex(SURFACE_WATER);
    iFrom[2]=pModel->GetStateVarIndex(PONDED_WATER);    iTo[2]=pModel->GetStateVarIndex(SOIL,1);
    iFrom[3]=pModel->GetStateVarIndex(SOIL,0);          iTo[3]=pModel->GetStateVarIndex(SOIL,1);
    iFrom[4]=pModel->GetStateVarIndex(SOIL,1);          iTo[4]=pModel->GetStateVarIndex(SURFACE_WATER);
  }
  else if(type==INF_AWBM) {
    CHydroProcessABC::DynamicSpecifyConnections(5);
    iFrom[0]=pModel->GetStateVarIndex(PONDED_WATER);    iTo[0]=pModel->GetStateVarIndex(SOIL,0);
    iFrom[1]=pModel->GetStateVarIndex(PONDED_WATER);    iTo[1]=pModel->GetStateVarIndex(SOIL,1);
    iFrom[2]=pModel->GetStateVarIndex(PONDED_WATER);    iTo[2]=pModel->GetStateVarIndex(SOIL,2);
    iFrom[3]=pModel->GetStateVarIndex(PONDED_WATER);    iTo[3]=pModel->GetStateVarIndex(SURFACE_WATER);
    iFrom[4]=pModel->GetStateVarIndex(PONDED_WATER);    iTo[4]=pModel->GetStateVarIndex(SOIL,3);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvInfiltration::~CmvInfiltration(){}

//////////////////////////////////////////////////////////////////
/// \brief Checks that process moves water to soil or lumped landform
//
void CmvInfiltration::Initialize()
{
  ExitGracefullyIf(pModel->GetStateVarType(iTo[0])!=SOIL,
                   "CmvInfiltration::Initialize:Infiltration must go to soil ",BAD_DATA);

  //Lumped landform only valid for SCS, partition coefficient

}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for infiltration algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by infiltration algorithm (size of aP[] and aPC[])
//
void CmvInfiltration::GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const
{

  if (type==INF_SCS)
  {
    nP=2;
    aP[0]="SCS_CN";           aPC[0]=CLASS_LANDUSE;
    aP[1]="SCS_IA_FRACTION";  aPC[1]=CLASS_LANDUSE;
  }
  else if (type==INF_SCS_NOABSTRACTION)
  {
    nP=1;
    aP[0]="SCS_CN";           aPC[0]=CLASS_LANDUSE;
  }
  else if (type==INF_RATIONAL)
  {
    nP=1;
    aP[0]="PARTITION_COEFF"; aPC[0]=CLASS_LANDUSE;
  }
  else if (type==INF_ALL_INFILTRATES)
  {
    nP=1;
    aP[0]="IMPERMEABLE_FRAC"; aPC[0]=CLASS_LANDUSE;
  }
  else if ((type==INF_GREEN_AMPT) || (type==INF_GA_SIMPLE))
  {
    nP=4;
    aP[0]="HYDRAUL_COND";      aPC[0]=CLASS_SOIL;
    aP[1]="WETTING_FRONT_PSI"; aPC[1]=CLASS_SOIL;
    aP[2]="POROSITY";          aPC[2]=CLASS_SOIL;
    aP[3]="IMPERMEABLE_FRAC";  aPC[3]=CLASS_LANDUSE;
  }
  else if (type==INF_UPSCALED_GREEN_AMPT)
  {
    nP=5;
    aP[0]="HYDRAUL_COND";      aPC[0]=CLASS_SOIL;
    aP[1]="WETTING_FRONT_PSI"; aPC[1]=CLASS_SOIL;
    aP[2]="KSAT_STD_DEVIATION";aPC[2]=CLASS_SOIL;
    aP[3]="POROSITY";          aPC[3]=CLASS_SOIL;
    aP[4]="IMPERMEABLE_FRAC";  aPC[4]=CLASS_LANDUSE;
  }
  else if (type==INF_PHILIP_1957)
  {
    nP=0;
  }
  else if(type==INF_HMETS)
  {
    nP=2;
    aP[0]="HMETS_RUNOFF_COEFF"; aPC[0]=CLASS_LANDUSE;
    aP[1]="IMPERMEABLE_FRAC";   aPC[1]=CLASS_LANDUSE;
  }
  else if (type==INF_VIC)
  {
    nP=5;
    aP[0]="VIC_ZMAX";         aPC[0]=CLASS_SOIL;
    aP[1]="VIC_ZMIN";         aPC[1]=CLASS_SOIL;
    aP[2]="VIC_ALPHA";        aPC[2]=CLASS_SOIL;
    aP[3]="POROSITY";         aPC[3]=CLASS_SOIL;
    aP[4]="IMPERMEABLE_FRAC"; aPC[4]=CLASS_LANDUSE;
  }
  else if (type==INF_VIC_ARNO)
  {
    nP=3;
    aP[0]="VIC_B_EXP";        aPC[0]=CLASS_SOIL;
    aP[1]="POROSITY";         aPC[1]=CLASS_SOIL;
    aP[2]="IMPERMEABLE_FRAC"; aPC[2]=CLASS_LANDUSE;
  }
  else if (type==INF_TOPMODEL)
  {
    nP=0;
    // algorithm not completed
  }
  else if (type==INF_PRMS)
  {
    nP=5;
    aP[0]="MAX_SAT_AREA_FRAC"; aPC[0]=CLASS_LANDUSE;
    aP[1]="POROSITY";          aPC[1]=CLASS_SOIL;
    aP[2]="FIELD_CAPACITY";    aPC[2]=CLASS_SOIL;
    aP[3]="SAT_WILT";          aPC[3]=CLASS_SOIL;
    aP[4]="IMPERMEABLE_FRAC";  aPC[4]=CLASS_LANDUSE;
  }
  else if (type==INF_HBV)
  {
    nP=3;
    aP[0]="POROSITY";         aPC[0]=CLASS_SOIL;
    aP[1]="HBV_BETA";         aPC[1]=CLASS_SOIL;
    aP[2]="IMPERMEABLE_FRAC"; aPC[2]=CLASS_LANDUSE;
  }
  else if (type==INF_UBC)
  {
    nP=6;
    aP[0]="IMPERMEABLE_FRAC";      aPC[0]=CLASS_LANDUSE;
    aP[1]="MAX_PERC_RATE";         aPC[1]=CLASS_SOIL;
    aP[2]="UBC_INFIL_SOIL_DEF";    aPC[2]=CLASS_SOIL;
    aP[3]="POROSITY";              aPC[3]=CLASS_SOIL;
    aP[4]="UBC_GW_SPLIT";          aPC[4]=CLASS_GLOBAL;
    aP[5]="UBC_FLASH_PONDING";     aPC[5]=CLASS_GLOBAL;
  }
  else if (type==INF_GR4J)
  {
    nP=2;
    aP[0]="POROSITY";              aPC[0]=CLASS_SOIL;
    aP[1]="IMPERMEABLE_FRAC";      aPC[1]=CLASS_LANDUSE;
  }
  else if(type==INF_XINANXIANG)
  {
    nP=6;
    aP[0]="IMPERMEABLE_FRAC";      aPC[0]=CLASS_LANDUSE;
    aP[1]="XINANXIANG_A";          aPC[1]=CLASS_SOIL;
    aP[2]="XINANXIANG_B";          aPC[2]=CLASS_SOIL;
    aP[3]="SAT_WILT";              aPC[3]=CLASS_SOIL;
    aP[4]="POROSITY";              aPC[4]=CLASS_SOIL;
    aP[5]="FIELD_CAPACITY";        aPC[5]=CLASS_SOIL;
    aP[6]="XINANXIANG_SHP";        aPC[6]=CLASS_SOIL;
  }
  else if(type==INF_PDM)
  {
    nP=3;
    aP[0]="IMPERMEABLE_FRAC";      aPC[0]=CLASS_LANDUSE;
    aP[1]="PDM_B";                 aPC[1]=CLASS_LANDUSE;
    aP[2]="POROSITY";              aPC[2]=CLASS_SOIL;
  }
  else if(type==INF_AWBM)
  {
    nP=3;
    aP[0]="AWBM_AREAFRAC1";        aPC[0]=CLASS_LANDUSE;
    aP[1]="AWBM_AREAFRAC2";        aPC[1]=CLASS_LANDUSE;
    aP[2]="AWBM_BFLOW_INDEX";      aPC[2]=CLASS_LANDUSE;
  }
  else
  {
    ExitGracefully("CmvInfiltration::GetParticipatingParamList: undefined infiltration algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
///
/// \param itype [in] User-specified infiltration model
/// \param *aSV [out] Array of state variable types needed by infiltration algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by infiltration algorithm (size of aSV[] and aLev[] arrays)
//
void CmvInfiltration::GetParticipatingStateVarList(infil_type  itype,sv_type *aSV, int *aLev, int &nSV)
{
  nSV=2;
  aSV[0]=PONDED_WATER;  aLev[0]=DOESNT_EXIST;
  aSV[1]=SURFACE_WATER; aLev[1]=DOESNT_EXIST;
  //user specified 'to' compartment (lumped landform or soil layer)

  //special cases - modify more than just ponded water, soil and surface water storage
  if (itype==INF_GA_SIMPLE)
  { //Tracks and updates cumulative infiltration
    nSV=4;
    aSV[2]=CUM_INFIL;        aLev[2]=DOESNT_EXIST;
    aSV[3]=GA_MOISTURE_INIT; aLev[3]=DOESNT_EXIST;
  }
  else if (itype==INF_UBC)
  { //requires 4 soil layers (topsoil, interflow, GW_U, GW_L)
    nSV=5;
    aSV[2]=SOIL;    aLev[2]=1;
    aSV[3]=SOIL;    aLev[3]=2;
    aSV[4]=SOIL;    aLev[4]=3;
  }
  else if (itype==INF_HMETS)
  { //requires 4 soil layers (topsoil, interflow, GW_U, GW_L)
    nSV=4;
    aSV[2]=CONVOLUTION;    aLev[2]=0;
    aSV[3]=CONVOLUTION;    aLev[3]=1;
  }
  //...
}

//////////////////////////////////////////////////////////////////
/// \brief Partitions ponded water between infiltration and runoff [mm/day]
/// \note  most routines completely deplete ponded water
/// \details
/// if type=INF_RATIONAL
///             <ul> <li> simple rational method - some percentage of ponded </ul>
/// if type==INF_SCS
///             <ul> <li> SCS curve number used, abstraction estimated as Ia ~0.2*S </ul>
/// if type==INF_SCS
///             <ul> <li> SCS curve number used, abstraction is modeled independently </ul>
/// if type==INF_GREEN_AMPT
///     <ul> <li> Green Ampt model (1911) - implicit formulation </ul> \cite green1911JoAS
/// if type==INF_GA_SIMPLE
///     <ul> <li> Green Ampt model -simple formulation </ul>
/// if type==INF_UBC
///     <ul> <li> Simple infiltration as performed in the UBC watershed model </ul>
/// \ref c) Michael Quick
/// if type==INF_PHILIP_1957
///     <ul> <li> Simple infiltration as calculated using the Philip 1957 solution </ul> \cite Philip1957SS
/// if type==INF_VIC
///     <ul> <li> Calculates infiltration using the VIC (Variable Infiltration Capacity) model </ul>
///    \ref of Wood et al. 1992. This is based upon the SPM formalism developed by Kavetski et al., 2003 \cite Wood1992JoGR \cite Kavetski2003WRR
/// if type==INF_VIC_ARNO
///     <ul> <li> Calculates infiltration using VIC/ARNO (Variable Infiltration Capacity) </ul>
///     if type==INF_PRMS
///             <ul> <li> Salculates infiltration using PRMS model as defined in FUSE </ul>
///  In all cases:
///   - infiltration=rates[0]
///   -  runoff      =rates[1]
/// most methods will completely deplete ponded water
/// \ref Clark et al. 2008) \cite Clark2008WRR
///
/// \param *state_vars [in] Current array of state variables
/// \param *pHRU [in] Pointer to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Array of rates of change of state variables [mm/day]
//
void CmvInfiltration::GetRatesOfChange (const double              *state_vars,
                                        const CHydroUnit  *pHRU,
                                        const optStruct   &Options,
                                        const time_struct &tt,
                                        double      *rates) const
{

  if ((pHRU->GetHRUType()!=HRU_STANDARD) && (pHRU->GetHRUType()!=HRU_ROCK)){return;}//Lakes & glaciers

  double runoff;
  double rainthru;
  double Fimp;
  double ponded_water;

  Fimp=pHRU->GetSurfaceProps()->impermeable_frac;

  ponded_water=max(state_vars[pModel->GetStateVarIndex(PONDED_WATER)],0.0);

  rainthru=(ponded_water/Options.timestep);//potential infiltration rate, mm/d

  if(pHRU->GetHRUType()==HRU_ROCK){ rates[1]=rainthru; return; } //if rock, nothing infiltrates, everything runs off

  int iTopSoil  =pModel->GetStateVarIndex(SOIL,0);

  //-----------------------------------------------------------------
  if (type==INF_RATIONAL)
  {
    runoff=pHRU->GetSurfaceProps()->partition_coeff*rainthru;
    rates[0]=rainthru-runoff;
    rates[1]=runoff;
  }
  //-----------------------------------------------------------------
  else if ((type==INF_SCS) || (type==INF_SCS_NOABSTRACTION))
  {
    runoff=GetSCSRunoff(pHRU,Options,tt,rainthru);
    rates[0]=rainthru-runoff;
    rates[1]=runoff;
  }
  //-----------------------------------------------------------------
  else if (type==INF_ALL_INFILTRATES)
  {
    runoff=0.0;
    runoff=(Fimp)*rainthru+(1-Fimp)*runoff; //correct for impermeable surfaces
    rates[0]=rainthru-runoff;
    rates[1]=runoff;
  }
  //-----------------------------------------------------------------
  else if (type==INF_PHILIP_1957)
  {
    //< \todo [add funct] Add Philip 1957 infiltration alg.
    ExitGracefully("CmvInfiltration::INF_PHILIP_1957",STUB);
  }
  //-----------------------------------------------------------------
  else if (type==INF_VIC)
  {
    double gamma,K1,sat,Smax;
    double zmax,zmin,alpha;
    double stor,max_stor;

    stor    =state_vars[iTopSoil];
    max_stor=pHRU->GetSoilCapacity(0);
    zmax    =pHRU->GetSoilProps(0)->VIC_zmax;
    zmin    =pHRU->GetSoilProps(0)->VIC_zmin;
    alpha   =pHRU->GetSoilProps(0)->VIC_alpha;

    sat     =min(stor/max_stor,1.0);
    gamma   =1.0/(alpha+1.0);
    K1      =pow((zmax-zmin)*alpha*gamma,-gamma);
    Smax    =gamma*(alpha*zmax+zmin);

    runoff=rainthru*threshMax(1.0-K1*pow(Smax-sat,gamma),1.0,0.0);

    runoff=(Fimp)*rainthru+(1-Fimp)*runoff; //correct for impermeable surfaces

    rates[0]=rainthru-runoff;
    rates[1]=runoff;
  }
  //-----------------------------------------------------------------
  else if (type==INF_HBV)
  {
    double beta,stor,max_stor,sat;
    stor    =state_vars[iTopSoil];
    max_stor=pHRU->GetSoilCapacity(0);
    beta    =pHRU->GetSoilProps(0)->HBV_beta;
    sat     =max(min(stor/max_stor,1.0),0.0);

    runoff=pow(sat,beta)*rainthru;

    runoff=(Fimp)*rainthru+(1-Fimp)*runoff; //correct for impermeable surfaces

    rates[0]=rainthru-runoff;
    rates[1]=runoff;
  }
  //-----------------------------------------------------------------
  else if (type==INF_VIC_ARNO)
  {
    double sat_area,stor,max_stor,b,sat;

    stor    =state_vars[iTopSoil];
    max_stor=pHRU->GetSoilCapacity(0);        //maximum storage of top soil layer [mm]
    b       =pHRU->GetSoilProps(0)->VIC_b_exp;//ARNO/VIC b exponent for runoff [-]
    sat     =min(stor/max_stor,1.0);          //soil saturation
    sat_area=1.0 - pow(1.0-sat,b);            //saturated area [-] fraction

    runoff  =sat_area*rainthru;

    runoff=(Fimp)*rainthru+(1-Fimp)*runoff; //correct for impermeable surfaces

    rates[0]=rainthru-runoff;
    rates[1]=runoff;
  }
  //-----------------------------------------------------------------
  else if (type==INF_PRMS)
  {
    double stor,tens_stor,sat_frac;
    sat_frac =pHRU->GetSurfaceProps()->max_sat_area_frac;
    tens_stor=pHRU->GetSoilTensionStorageCapacity(0);
    stor     =state_vars[iTopSoil];

    runoff =  sat_frac* min(stor/tens_stor,1.0) * rainthru;

    runoff=(Fimp)*rainthru+(1-Fimp)*runoff; //correct for impermeable surfaces

    rates[0]=rainthru-runoff;
    rates[1]=runoff;
  }
  //-----------------------------------------------------------------
  else if (type==INF_TOPMODEL)
  {
    double sat_area;                                                                                 //saturated area [-] fraction
    //double lambda    = pHRU->GetTerrainProps()->topmodel_lambda ;      //mean of the log-transformed topographical index [m]
    // double stor     = state_vars[iTopSoil];
    // double max_stor = pHRU->GetSoilCapacity(0);                //maximum storage of soil layer [mm]

    //double xi = lambda * pow(stor/max_stor,-1);

    sat_area = 0;///< \todo [add funct] calculate saturated area for TOPMODEL routine

    ExitGracefully("GetTOPMODELRUNOFF",STUB);

    runoff =  sat_area * rainthru;              //runoff rate [mm/d]

    runoff=(Fimp)*rainthru+(1-Fimp)*runoff; //correct for impermeable surfaces

    rates[0]=rainthru-runoff;
    rates[1]=runoff;
  }
  //-----------------------------------------------------------------
  else if ((type==INF_GREEN_AMPT) || (type==INF_GA_SIMPLE))
  {
    GetGreenAmptRunoff(state_vars,pHRU,Options,tt,rates,rainthru);
  }
  //-----------------------------------------------------------------
  else if (type==INF_UPSCALED_GREEN_AMPT)
  {
    runoff=GetHeterogeneousGreenAmptRunoff(rainthru,
                                           pHRU->GetSoilProps(0),
                                           pHRU->GetSoilThickness(0),
                                           state_vars[iTopSoil],
                                           tt,Options);
    runoff=(Fimp)*rainthru+(1-Fimp)*runoff; //correct for impermeable surfaces

    rates[0]=rainthru-runoff;
    rates[1]=runoff;
  }
  //-----------------------------------------------------------------
  else if (type==INF_UBC)
  {
    GetUBCWMRunoff(state_vars,pHRU,Options,tt,rates,rainthru);
  }
  //-----------------------------------------------------------------
  else if (type==INF_GR4J)
  {
    double x1=pHRU->GetSoilCapacity(0);
    double infil;
    double tmp;
    double stor=state_vars[iTopSoil];
    double sat=stor/x1;
    tmp=tanh(ponded_water/x1);
    infil=x1*(1.0-(sat*sat))*tmp/(1.0+sat*tmp);
    infil=(1.0-Fimp)*infil; //correct for impermeable surfaces
    rates[0]=infil;
    rates[1]=rainthru-infil;
  }
  //-----------------------------------------------------------------
  else if(type==INF_HMETS)
  {
    double infil,delayed,direct;

    double stor       =state_vars[iTopSoil];
    double max_stor   =pHRU->GetSoilCapacity(0);
    double coef_runoff=pHRU->GetSurfaceProps()->HMETS_runoff_coeff; //[-]

    double sat = min(max(stor/max_stor,0.0),1.0);

    direct=Fimp*rainthru;

    runoff=coef_runoff*(sat)*(1.0-Fimp)*rainthru; //[mm/d] 'horizontal transfer'

    infil=(1.0-Fimp)*rainthru-runoff;

    delayed=coef_runoff*pow(sat,2.0)*infil; //[mm/d]
    infil-=delayed;

    rates[0]=infil;     //PONDED->SOIL[0]
    rates[1]=direct;    //PONDED->SW
    rates[2]=runoff;    //PONDED->CONVOL[0]
    rates[3]=delayed;   //PONDED->CONVOL[1]
  }
  //-----------------------------------------------------------------
  else if(type==INF_XINANXIANG)
  {
    // from Jayawardena, AW and MC Zhou A modified spatial soil moisture storage capacity distribution curve for the Xinanjiang model.
    // Journal of Hydrology, 227(1-4), p93-113, 2000

    double b=1.0;//pHRU->GetSoilProps(0)->Xinanxiang_b;
    double c=1.0;//pHRU->GetSoilProps(0)->Xinanxiang_shp; [-0.5<c<0.5]
    double stor       =state_vars[iTopSoil];
    double tens_stor  =pHRU->GetSoilTensionStorageCapacity(0);
    double max_stor   =pHRU->GetSoilCapacity(0);
    double infil,sat_excess,direct;

    direct=(1.0-Fimp)*rainthru;

    double sat =max(0.0,stor/tens_stor);
    double sat1=max(0.0,stor/(max_stor-tens_stor));
    if(sat<=0.5-c) { infil=(    pow(0.5-c,1-b)*pow(  sat,b))*direct; }
    else           { infil=(1.0-pow(0.5+c,1-b)*pow(1-sat,b))*direct; }

    sat_excess=1.0;//(1.0-pow(1-sat1,n))*direct;

    ExitGracefully("INF_XINANXIANG",STUB);

    infil=(1.0-Fimp)*infil; //correct for impermeable surfaces (?)
    rates[0]=infil;          //PONDED->SOIL
    rates[1]=rainthru-infil; //PONDED->SW
  }
  //----------------------------------------------------------------------------
  else if(type==INF_PDM)
  { //uses PDM algorithm from Moore (1985)
    double b;          //[-] pareto distribution parameter
    double c_max;      //[mm] maximum local storage capacity in HRU
    double c_star;     //[mm] critical capacity, c*
    double infil;      //[mm] total amount abstracted

    b = pHRU->GetSurfaceProps()->PDM_b;

    double stor    =state_vars[iTopSoil];
    double max_stor=pHRU->GetSoilCapacity(0);

    c_max=(b+1)*max_stor;

    c_star=c_max*(1.0-pow(1.0-(stor/max_stor),1.0/(b+1.0)));

    //analytical evaluation of P*dt-equation 3 of Mekonnen et al. (2014); min() handles case where entire landscape sheds
    infil=max_stor*(pow(1.0-(c_star/c_max),b+1.0)-pow(1.0-min(c_star+ponded_water,c_max)/c_max,b+1.0));
    infil=min(infil,ponded_water);

    infil=(1.0-Fimp)*infil/Options.timestep; //correct for impermeable surfaces

    rates[0]=infil;          //PONDED->SOIL[0]
    rates[1]=rainthru-infil; //PONDED->SW
  }
  //----------------------------------------------------------------------------
  else if(type==INF_AWBM)
  { //uses Australia Water Balance algorithm from Boughton (2004)
    double infil1,infil2,infil3,runoff1,runoff2,runoff3,excess,to_GW;
    double def1,def2,def3; //[mm]
    double a1 =pHRU->GetSurfaceProps()->AWBM_areafrac1;
    double a2 =pHRU->GetSurfaceProps()->AWBM_areafrac1;
    double a3=1.0-a1-a2;
    double BFI=pHRU->GetSurfaceProps()->AWBM_bflow_index;

    def1=max(a1*pHRU->GetSoilCapacity(0)-state_vars[iFrom[0]],0.0); // deficit averaged over domain
    def2=max(a2*pHRU->GetSoilCapacity(1)-state_vars[iFrom[1]],0.0);
    def3=max(a3*pHRU->GetSoilCapacity(2)-state_vars[iFrom[2]],0.0);

    runoff1=max(a1*ponded_water-def1,0.0);
    runoff2=max(a2*ponded_water-def2,0.0);
    runoff3=max(a3*ponded_water-def3,0.0);
    infil1=a1*ponded_water-runoff1;
    infil2=a2*ponded_water-runoff2;
    infil3=a3*ponded_water-runoff3;

    excess=(1.0-BFI)*(runoff1+runoff2+runoff3);
    to_GW =(    BFI)*(runoff1+runoff2+runoff3);

    rates[0]=   infil1/Options.timestep;   //PONDED->SOIL[0]
    rates[1]=   infil2/Options.timestep;   //PONDED->SOIL[1]
    rates[2]=   infil3/Options.timestep;   //PONDED->SOIL[2]
    rates[3]=   excess/Options.timestep;   //PONDED->SURFACE_WATER
    rates[4]=   to_GW/Options.timestep;    //PONDED->SOIL[3]
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow cannot drain "from" compartment over timestep
///
/// \param *storage [in] Reference to compartments from which water infiltrates
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time
/// \param *rates [out] Corrected rate of infiltration [mm/day]
//
void   CmvInfiltration::ApplyConstraints( const double     *storage,
                                          const CHydroUnit *pHRU,
                                          const optStruct  &Options,
                                          const time_struct &tt,
                                          double     *rates) const
{
  if (pHRU->GetHRUType()!=HRU_STANDARD){return;}//Lakes & glaciers

  //cant remove more than is there (should never be an option)
  rates[0]=min(rates[0],storage[iFrom[0]]/Options.timestep);

  //reaching soil saturation level
  double inf=rates[0];
  double max_stor;
  if (!Options.allow_soil_overfill){
    max_stor=pHRU->GetStateVarMax(iTo[0],storage,Options);
    inf=min(rates[0],max(max_stor-storage[iTo[0]],0.0)/Options.timestep);
  }

  rates[1]+=(rates[0]-inf);
  rates[0]=inf;
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates runoff [mm/d] from rainfall using SCS curve number method
/// \note totalrain5days is total rainfall in the last five days [mm]
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current simulation time
/// \param &rainthru [in] Double rainfall/snowmelt rate [mm/d]
/// \return Runoff from rainfall [mm/d] according to SCS curve number method
//
double CmvInfiltration::GetSCSRunoff(const CHydroUnit *pHRU,
                                      const optStruct         &Options,
                                      const time_struct &tt,
                                      const double     &rainthru) const
{
  int    condition=2;     //antecedent moisture condition
  double S,CN,W,Weff;     //retention parameter [mm], CNII, rainfall [mm], runoff [mm]
  double TR;              //total rain over past 5 days [in]
  double Ia;              //initial abstraction, [mm]

  TR=pHRU->GetForcingFunctions()->precip_5day/MM_PER_INCH;
  CN=pHRU->GetSurfaceProps()->SCS_CN;

  //correct curve number for antecedent moisture conditions
  if((tt.month>4) && (tt.month<9)) {//growing season?? (northern hemisphere)
    //if (pHRU->GetForcingFunctions()->is_growing_season){
    //if (pHRU->GetStateVarValue(pModel->GetStateVarIndex(CROP_HEAT_UNITS))>-0.5){
    if     (TR<1.4) { condition=1; }
    else if(TR>2.1) { condition=3; }
  }
  else {
    if     (TR<0.5) { condition=1; }
    else if(TR>1.1) { condition=3; }
  }
  if     (condition==1) { CN = 5E-05 *pow(CN,3) + 0.0008*pow(CN,2) + 0.4431*CN; }//0.999R^2 with tabulated values (JRC)
  else if(condition==3) { CN = 7E-05 *pow(CN,3) - 0.0185*pow(CN,2) + 2.1586*CN; }//0.999R^2 with tabulated values (JRC)
  CN=min(CN,100.0);

  //calculate amount of runoff
  S   =MM_PER_INCH*(1000.0/CN-10.0);
  W   =rainthru*Options.timestep;//[mm]
  if(type!=INF_SCS_NOABSTRACTION)
  {
    Ia  =(pHRU->GetSurfaceProps()->SCS_Ia_fraction)*S;//SCS_Ia_fraction is =0.2 for standard SCS implementation
  }
  else
  {
    Ia=0.0;//abstraction handled by another routine
    //Note: when tied to an SCS Abstraction algorithm which explicitly tracks depression storage
    //"rainthru" is what remains after abstraction, i.e., W=W-Ia, or SCS_Ia_fraction=0.0
  }
  Weff=pow(threshPositive(W-Ia),2)/(W+(S-Ia));//[mm]
  if (W+(S-Ia)==0.0){Weff=0.0;}

  return Weff/Options.timestep;

}

//////////////////////////////////////////////////////////////////
/// \brief Calculates runoff [mm/d] from rainfall using UBCWM method
/// \copyright Michael Quick
///
/// \param *state_vars [in] Array of current state variables
/// \param *pHRU [in] Pointer to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current simulation time
/// \param *rates [out] Calculated runoff  and infiltration from rainfall [mm/d]
/// \param &rainthru [in] Double rainfall/snowmelt rate applied to ground surface [mm/d]
//
void CmvInfiltration::GetUBCWMRunoff    (const double             *state_vars,
                                         const CHydroUnit  *pHRU,
                                         const optStruct      &Options,
                                         const time_struct &tt,
                                         double      *rates,
                                         const double      &rainthru) const

{
  const double V0FLAX=1800.0; //[mm]

  double infil,to_GW,to_interflow,runoff;
  double flash_factor,b1,b2;

  double soil_deficit =max(0.0,pHRU->GetSoilCapacity(0)-state_vars[iTo[0]]);

  double Fimp          = pHRU->GetSurfaceProps()->impermeable_frac;
  double max_perc_rate = pHRU->GetSoilProps(1)->max_perc_rate;
  double P0AGEN        = pHRU->GetSoilProps(0)->UBC_infil_soil_def;
  double P0DSH         = pModel->GetGlobalParams()->GetParams()->UBC_GW_split;
  double V0FLAS        = pModel->GetGlobalParams()->GetParams()->UBC_flash_ponding;  //[mm]

  //calculate b1 parameter (effective permeable area without flash factor)
  b1=1.0;
  if (Fimp<1.0) {
    b1 = Fimp*pow(10.0, -soil_deficit / P0AGEN);
    //b1 = Fimp-(1.0-Fimp)*pow(10, -soil_deficit / P0AGEN);  //[GJ] better
  }

  //calculate flash factor (for high intensity storms)
  flash_factor=0.0;
  if (rainthru>V0FLAS){
    flash_factor = (1.0 + log(rainthru/V0FLAX)/log(V0FLAX/V0FLAS));
  }
  lowerswap(flash_factor,1.0);
  upperswap(flash_factor,0.0);

  //calculate b2 (effective permeable area pct)
  b2=(b1)*(1.0)+(1.0-b1)*flash_factor;

  g_debug_vars[0]=b2;//part of RFS Cheat to send b2 to glacial HRUs for infiltration calcs

  //distribute ponded water to surface water and soil stores
  infil       =min(soil_deficit/Options.timestep,rainthru      );//fills soil deficit
  to_GW       =min(max_perc_rate,                rainthru-infil); //overflows to GW
  to_interflow=(rainthru-infil-to_GW);//overflows to interflow

  infil*=1-b2;
  to_GW*=1-b2;
  to_interflow*=1-b2;

  runoff      =rainthru-infil-to_GW-to_interflow; //equivalently, b2*rainthru?

  //NS: This seems like a much more intuitive way to do this, but maybe breaks emulation?
  //runoff      = b2*rainthru;
  //infil       =min(soil_deficit/Options.timestep,rainthru-runoff      );//fills soil deficit
  //to_GW       =min(max_perc_rate,                rainthru-runoff-infil);//overflows to GW
  //to_interflow=(rainthru-runoff-infil-to_GW);                           //overflows to interflow

  // cout<<state_vars[iTo[0]]<<" "<<pHRU->GetSoilCapacity(0)<<" "<<pHRU->GetSoilThickness(0)<<" "<<soil_deficit<<" "<<rainthru<<" "<<infil<<" "<<runoff<<endl;
  rates[0]=infil;              // ponded water->soil deficit
  rates[1]=runoff;             // ponded water->surface water
  rates[2]=to_interflow;       // ponded water->interflow
  rates[3]=(1.0-P0DSH)*to_GW;  // ponded water->upper groundwater
  rates[4]=(    P0DSH)*to_GW;  // ponded water->lower groundwater

}
