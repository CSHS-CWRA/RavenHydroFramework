/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2022 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Properties.h"
#include "SoilAndLandClasses.h"
#include "HydroUnits.h"
#include "ModelABC.h"
//////////////////////////////////////////////////////////////////
/// \brief Calculates leaf conductance [mm/s]
/// \ref using Stewart (1988) model \cite Stewart1988AaFM
/// \ref Adapted from Dingman pg. 298-299 \cite Dingman1994
/// \param &max_leaf_cond Maximum leaf conductance [mm/s]
/// \param &min_leaf_cond Minimum leaf conductance [mm/s]
/// \param *F [in] Forcing functions for a specific HRU
/// \param &soil_moist_deficit [in] Soil moisture deficit [mm]
/// \return Current leaf conductance [mm/s]
//
inline double CalculateLeafConductance( const double       &max_leaf_cond,     //[mm/s]
                                        const double       &min_leaf_cond,     //[mm/s]
                                        const force_struct *F,
                                        const double       &soil_moist_deficit)//[mm]
{
  double f1,f2,f3,f4;
  double hum_deficit;

  if (F==NULL){return 0.5*max_leaf_cond;}//arbitrary

  //light effects
  f1=12.78*(F->SW_radia)/(11.57*(F->SW_radia)+104.4);

  ///< \ref Brook90: \cite Federer2010
  /// \todo [fix] - Brook90 corrections for radiation
  /*f1=0.0;
    double radiation=(F->LW_radia_net+F->SW_radia_net)*MJ_PER_D_TO_WATT;
    const double RM=1000;//[W/m2]
    const double R5=75;//[W/m2]
    double R0= RM * R5 / (RM - 2.0 * R5);
    double FS = 1.0 + 0.5*(VV->SAI)/(VV->LAI);
    f1=((RM+R0)/(RM*CAN_EXTINCT*FS))*log((R0+CAN_EXTINCT*radiation)/(R0+CAN_EXTINCT*radiation*exp(-CAN_EXTINCT*FS*LAI)));
  */

  //vapor-pressure deficit effects
  //double sat_vap=GetSaturatedVaporPressure(F->temp_ave);
  hum_deficit=(1.0-(F->rel_humidity))*AIR_H20_MW_RAT*(F->air_dens) /(F->air_pres);//eqn D8C [kg/m3]
  f2=1.0-66.6*hum_deficit;

  ///< \ref vapor deficit limitation (Brook90)- Lohammar et al. (1980) and Stannard (1993) \cite Lohammar1980EB \cite Stannard1993WRR
  //const double CVDP=2.0;
  //f2 = 1.0 / (1.0 + VPD / CVPD);//VPD - vapor deficit

  //leaf temperature
  f3=(F->temp_ave)*pow(max(40.0-(F->temp_ave),0.0),1.18)/691.0;

  ///< \ref temperature limitation (Brook90) \cite Federer2010
  /*double Ta=F->temp_daily;
    const double TL(0.0),T1(10.0),T2(30.0),TH(40.0);
    if      ((Ta<=TL)            ) {f3 = 0.0;}
    else if ((Ta> TL) && (Ta< T1)) {f3 = 1.0 - pow((T1 - Ta) / (T1 - TL),2);}
    else if ((Ta>=T1) && (Ta<=T2)) {f3 = 1.0;}
    else if ((Ta> T2) && (Ta< TH)) {f3 = 1.0 - pow((Ta - T2) / (TH - T2),2);}
    else                           {f3 = 0.0;}*/

  //leaf water content
  f4=1.0-0.00119*exp(8.1*soil_moist_deficit);   //needs to be fixed, doesn't work properly
  f4=1.0;

  //cout<<"f1-f4: "<<f1<<" "<<f2<<" "<<f3<<" "<<f4<<endl;
  //cout<<"soil water: "<<soil_moist_deficit<<" "<<max_leaf_cond<<" "<< min_leaf_cond<<endl;

  return min(max(f1*f2*f3*f4,0.0),1.0)*(max_leaf_cond-min_leaf_cond)+min_leaf_cond;

}

//////////////////////////////////////////////////////////////////
/// \brief Sets canopy properties based upon meteorological conditions and time of year
/// \remark Called at the start of each time step for each HRU
///
/// \param &VV [out] Vegetation properties strcture
/// \param *pHRU [in] HRU class object
/// \param *pModel [in] pointer to model object
/// \param &tt [in] time structure
/// \param &Options [in] Options structure
//
void CVegetationClass::RecalculateCanopyParams (      veg_var_struct    &VV,
                                                      const CHydroUnit       *pHRU,
                                                      const CModelABC        *pModel,
                                                      const time_struct      &tt,
                                                      const optStruct        &Options)
{
  VV.shelter_factor=0.5;// should be variable

  double max_height  =pHRU->GetVegetationProps()->max_height;
  double sparseness  =pHRU->GetSurfaceProps()->forest_sparseness;
  double max_LAI     =pHRU->GetVegetationProps()->max_LAI;
  double SAI_ht_ratio=pHRU->GetVegetationProps()->SAI_ht_ratio;
  double extinction  =pHRU->GetVegetationProps()->svf_extinction;
  double max_SAI     =max_height*SAI_ht_ratio;

  //Vegetation Height
  //------------------------------------------------------------
  VV.height=max_height*InterpolateMo(pHRU->GetVegetationProps()->relative_ht,tt,Options);

  //LAI/SAI
  //------------------------------------------------------------
  VV.LAI=(1.0-sparseness)*max_LAI*InterpolateMo(pHRU->GetVegetationProps()->relative_LAI,tt,Options);// \ref From Brook90 CANOPY Routine \cite Federer2010
  VV.SAI=(1.0-sparseness)*(VV.height*SAI_ht_ratio);/// \ref From Brook90 CANOPY Routine

  //LAI /SAI snow corrections From Brook90 CANOPY Routine
  /*
      double snowdepth        =0.0;
      int iSnowD=pModel->GetStateVarIndex(SNOW_DEPTH);
      if (iSnowD!=DOESNT_EXIST){snowdepth=pHRU->GetStateVarValue(iSnowD);}

      double height_above_snow;//[m]
      height_above_snow=threshPositive(VV.height - snowdepth/MM_PER_METER);//
      VV.LAI*=(height_above_snow/VV.height);
      VV.SAI*=(height_above_snow/VV.height); */

  //Rain/Snow storage capacity
  //------------------------------------------------------------
  double max_capacity    =pHRU->GetVegetationProps()->max_capacity;
  double max_sno_capacity=pHRU->GetVegetationProps()->max_snow_capacity;

  if ((max_LAI+max_SAI)>0){
    VV.capacity      =((VV.LAI+VV.SAI)/(max_LAI+max_SAI))*max_capacity;
    VV.snow_capacity =((VV.LAI+VV.SAI)/(max_LAI+max_SAI))*max_sno_capacity;
  }
  else{
    VV.capacity     =0.0;
    VV.snow_capacity=0.0;
  }

  //Skyview factor
  //------------------------------------------------------------
  VV.skyview_fact=exp(-extinction*(VV.LAI+VV.SAI)); //sky view factor (percentage of ground that recieves sunlight)

  //riparian skyview factor (from Moore, Leach, and Knudson 2014 with symmetric sides)
  //double D=pHRU->GetSurfaceProps()->streamwidth;
  //double H=VV.height;
  //VV.skyview_fact*=1/D*(sqrt(H*H+D*D)-sqrt(H*H));

  //Rain/Snow interception factors
  //------------------------------------------------------------
  if (Options.interception_factor==PRECIP_ICEPT_LAI)
  {///< Proportional to LAI \cite Dingman2002
    VV.rain_icept_pct=min(pHRU->GetVegetationProps()->rain_icept_fact*(VV.LAI+VV.SAI),1.0);///< \ref from Brook90 Dingman 7.2 (C.rain_icept_fact~0.06) \cite Federer2010
    VV.snow_icept_pct=min(pHRU->GetVegetationProps()->snow_icept_fact*(VV.LAI+VV.SAI),1.0);///< \ref from Brook90 Dingman 5.1 pg 217 (C.rain_icept_fact~0.04)
  }
  else if (Options.interception_factor==PRECIP_ICEPT_USER)
  {//User specified (UBC,HBV)
    VV.rain_icept_pct =pHRU->GetVegetationProps()->rain_icept_pct;
    VV.snow_icept_pct =pHRU->GetVegetationProps()->snow_icept_pct;
  }
  else if (Options.interception_factor==PRECIP_ICEPT_EXPLAI)
  {///< \ref from CLM eqn 7.2 \cite Oleson2012
    VV.rain_icept_pct=(1.0-exp(-0.5*(VV.LAI+VV.SAI)));
    VV.snow_icept_pct=(1.0-exp(-0.5*(VV.LAI+VV.SAI)));
  }
  else if(Options.interception_factor==PRECIP_ICEPT_NONE) {
    VV.rain_icept_pct=VV.snow_icept_pct=0.0;
  }

  int iCanSnow = pModel->GetStateVarIndex(CANOPY_SNOW);

  if (Options.interception_factor == PRECIP_ICEPT_HEDSTROM)
  {
    VV.rain_icept_pct=(1.0-exp(-0.5*(VV.LAI+VV.SAI)));
    if(iCanSnow != DOESNT_EXIST)
    { ///< \ref from Hedstrom & Pomeroy, 1998
      double max_snow_load= pHRU->GetVegetationProps()->max_snow_load; //[kg/m2] ~5.9-6.6
      double max_cap      = pHRU->GetVegetationProps()->max_snow_capacity;
      double rho_s        = CalcFreshSnowDensity(pHRU->GetForcingFunctions()->temp_ave);
      
      VV.snow_capacity = VV.LAI * max_snow_load *(MM_PER_METER/ DENSITY_WATER) * (0.27 + 46 / rho_s); //kg/m2*m3/kg->mm 

      if (max_cap > 0.0) {
       lowerswap(VV.snow_capacity, max_cap);//storage capacity override for floppy trees [M. Chernos]
      }

      double P       =pHRU->GetForcingFunctions()->snow_frac*pHRU->GetForcingFunctions()->precip* Options.timestep; //[mm]
      double stor    =pHRU->GetStateVarValue(iCanSnow);//[mm]
      double max_stor=VV.snow_capacity;//[mm]

      VV.snow_icept_pct = max(min(((max_stor-stor)/P) * (1 - exp(-(1.0-sparseness)* (P/max_stor))),1.0),0.0);
    }
    else { VV.snow_icept_pct = 0.0; }
  }
  else if(Options.interception_factor == PRECIP_ICEPT_STICKY) 
  {
    VV.rain_icept_pct=(1.0-exp(-0.5*(VV.LAI+VV.SAI)));

    if(iCanSnow != DOESNT_EXIST) 
    {///< \ref from SUMMA
      double gamma;
      double max_snow_load=pHRU->GetVegetationProps()->max_snow_load; //kg/m2 leaf
      double stor         =pHRU->GetStateVarValue(iCanSnow);//[mm]
      double Tave         =pHRU->GetForcingFunctions()->temp_ave;

      if      (Tave>-1.0){gamma=4.0;}
      else if (Tave>-3.0){gamma=1.5*Tave+5.5;}
      else               {gamma=1.0;}

      VV.snow_capacity =VV.LAI*max_snow_load*(MM_PER_METER/ DENSITY_WATER) *gamma;
      VV.snow_icept_pct=1.0-min(stor/VV.snow_capacity,1.0);
    }
    else{ VV.snow_icept_pct = 0.0; }
  }

  if ((Options.orocorr_precip==OROCORR_UBCWM) || (Options.orocorr_precip==OROCORR_UBCWM2))
  {
    //interception is (inelegantly, but necessarily) handled in the precipitation correction routine in UBCWM
    VV.rain_icept_pct=0;
    VV.snow_icept_pct=0;
  }

  //Leaf/Canopy Conductances
  //------------------------------------------------------------
  double soilmoist_deficit=0.5;//soil moisture deficit [mm]
  int iTopSoil=pModel->GetStateVarIndex(SOIL,0);
  if (iTopSoil!=DOESNT_EXIST){soilmoist_deficit=pHRU->GetSoilCapacity(0)-pHRU->GetStateVarValue(iTopSoil);}

  double max_leaf_cond=pHRU->GetVegetationProps()->max_leaf_cond;
  VV.leaf_cond      =CalculateLeafConductance(max_leaf_cond,0.0,pHRU->GetForcingFunctions(),soilmoist_deficit);

  VV.canopy_conductance = VV.shelter_factor*VV.leaf_cond*VV.LAI;///< \ref Dingman 7.54  \cite Dingman1994
  //VV.canopy_conductance = (5 * LAI); /// \ref Howell and Evett
  // \todo [funct]: may wish to attend to conductance in non-vegetated lands (conductance should go to infinity?)

  /*
  // Canopy Roughness parameters ( \ref from Howell, T.A and Evett, S.R., USDA-ARS)
  //------------------------------------------------------------
  VV.zero_pln_disp = 0.666 * VV.height;     ///< zero plane displacement height [m] (\ref  0.7 for Dingman eqn. 7-49) \cite Dingman1994
  VV.roughness     = 0.123 * VV.height;     /// momentum roughness length [m] ( \ref 0.1 for Dingman eqn. 7-49)
  double vap_rough_ht = 0.1 * VV.roughness; /// vapor roughness length [m] (\ref 1.0 for Dingman eqn. 7-49)
  */

  // Canopy Roughness parameters (\ref from Brook90 ROUGH routine)
  //------------------------------------------------------------
  double ratio,xx;

  //find figures for closed canopy
  double closed_roughness;//[m]       roughness length for closed canopy
  double closed_zerodisp; //[m]       zero-plane displacement for closed canopy
  closed_roughness=CVegetationClass::CalcClosedRoughness    (VV.height);
  closed_zerodisp =CVegetationClass::CalcClosedZeroPlaneDisp(VV.height,closed_roughness);
  upperswap(closed_roughness,pHRU->GetSurfaceProps()->roughness);

  ratio=(VV.LAI+VV.SAI)/(CLOSED_LAI+SAI_ht_ratio*VV.height); //(LAI + SAI) / (LAI + SAI)_closed canopy

  if (ratio>=1.0)
  {//closed canopy
    VV.zero_pln_disp =closed_zerodisp;
    VV.roughness     =closed_roughness;
  }
  else
  {///< Sparse canopy - \ref use Shuttleworth & Gurney, 1990 \cite shuttleworth1990QURMS
    xx=0.0;
    if (VV.height>REAL_SMALL){
      xx=ratio*pow(-1.0+exp(0.909-3.03*closed_zerodisp/VV.height),4);
    }
    VV.zero_pln_disp =1.1*VV.height*log(1.0+pow(xx,0.25));
    VV.roughness     =min(0.3*(VV.height-VV.zero_pln_disp),pHRU->GetSurfaceProps()->roughness+0.3*VV.height*sqrt(xx));
  }
  VV.reference_height = VV.height + Z_REF_ADJUST;
}

//////////////////////////////////////////////////////////////
/// \brief Calculates roughness length for closed canopy
/// \param &height [in] Canopy height [m]
/// \return Roughness length for closed canopy [m]
//
double CVegetationClass::CalcClosedRoughness(const double &height)
{
  double closed_roughness;//[m]       roughness length for closed canopy

  if      (height>=CZR_HEIGHT){closed_roughness=CZR*height;}
  else if (height<=CZS_HEIGHT){closed_roughness=CZS*height;}
  else                        {
    closed_roughness=CZS*CZS_HEIGHT;
    closed_roughness+=(CZR*CZR_HEIGHT-CZS*CZS_HEIGHT)*(height-CZS_HEIGHT)/(CZR_HEIGHT-CZS_HEIGHT);
  }
  return closed_roughness;
}

//////////////////////////////////////////////////////////////
/// \brief Calculates closed zero plane displacement
/// \param &height [in] Canopy height [m]
/// \param &closed_rough [in] closed roughness length [m]
/// \return Closed zero plane displacement [m]
//
double CVegetationClass::CalcClosedZeroPlaneDisp(const double &height, const double &closed_rough)
{
  return height-closed_rough/0.3;///< \ref (Monteith 1973; Shuttleworth and Gurney 1990). \cite monteith1973L \cite shuttleworth1990QURMS
}

//////////////////////////////////////////////////////////////
/// \brief Write canopy properties to screen
/// \param &V [in] Vegetation structure
/// \param &VV [in] Vegetation auxiliary properties
//
void CVegetationClass::WriteVegetationProps(const veg_struct &V, const veg_var_struct &VV)
{
  cout<<endl;
  cout<<"max_capacity:....."<<V.max_capacity<<endl;
  cout<<"max_snow_capacity."<<V.max_snow_capacity<<endl;
  cout<<"max_leaf_cond....."<<V.max_leaf_cond<<endl;
  cout<<"max_LAI..........."<<V.max_LAI<<endl;
  cout<<"max_height........"<<V.max_height<<endl;
  cout<<"SAI_ht_ratio......"<<V.SAI_ht_ratio<<endl;
  cout<<"trunk_fraction...."<<V.trunk_fraction<<endl;
  cout<<"---------"<<endl;
  cout<<"SAI..............."<<VV.SAI<<endl;               //Stem Area Index [m2/m2]
  cout<<"LAI..............."<<VV.LAI<<endl;               //Leaf Area Index [m2/m2]
  cout<<"height............"<<VV.height<<endl;            //vegetation height [m]
  cout<<"capacity.........."<<VV.capacity<<endl;          //rain storage capacity [mm]
  cout<<"snow_capacity....."<<VV.snow_capacity<<endl;     //snow storage capacity [mm SWE]
  cout<<"shelter_factor...."<<VV.shelter_factor<<endl;    //accounts for sheltered leaves : about 0->5-1->0
  cout<<"leaf_cond........."<<VV.leaf_cond<<endl;         //[mm/s]
  cout<<"canopy_conductance"<<VV.canopy_conductance<<endl;
  cout<<"rain_icept_pct...."<<VV.rain_icept_pct<<endl;    //percentage of rain intercepted [0->->1] (only on canopy portion)
  cout<<"snow_icept_pct...."<<VV.snow_icept_pct<<endl;    //percentage of snow intercepted [0..1] (only on canopy portion)
//cout<<"pot_intercept_rate; //potential interception rate [mm/d]
}

//////////////////////////////////////////////////////////////////
/// \brief Recalculates root parameters
/// \details calculates root length, plant, root, and xylem resistance,
/// cowan alpha (should be for each layer?)
/// \ref Adapted from Brook90 routine PLNTRES and CANOPY \cite Federer2010
///
/// \param &VV [out] Vegetation properties strcture
/// \param *pHRU [in] HRU class object
/// \param *pModel [in] pointer to model object
/// \param &tt [in] time structure
/// \param &Options [in] Options structure
//
void CVegetationClass::RecalculateRootParams(       veg_var_struct   &VV,
                                                    const CHydroUnit       *pHRU,
                                                    const CModelABC        *pModel,
                                                    const time_struct      &tt,
                                                    const optStruct        &Options)
{
  double delT;              //[0..1] root volume fraction in the layer
  double sum=0.0;           //[mm]   total relative length
  double DD[MAX_SOILLAYERS];//[0..1] stonefree layer thicknesses [mm]
  int m;

  double max_height     =pHRU->GetVegetationProps()->max_height;
  double sparseness     =pHRU->GetSurfaceProps()->forest_sparseness;
  double max_root_length=pHRU->GetVegetationProps()->max_root_length;
  double min_resistivity=pHRU->GetVegetationProps()->min_resistivity;
  double rootradius     =pHRU->GetVegetationProps()->rootradius;
  double xylem_frac     =pHRU->GetVegetationProps()->xylem_frac;

  /// root length and resistivity calculations (\ref from Brook90 CANOPY & PLNTRES routines)
  double frac=(1.0-sparseness)*(VV.height /max_height);
  VV.root_length   =frac         *max_root_length;
  VV.resistivity   =1.0/frac     *min_resistivity;


  /// Cowan alpha and root resistance calculations
  for (m=0;m<pModel->GetNumSoilLayers();m++)
  {
    DD[m]=(pHRU->GetSoilThickness(m)/MM_PER_METER)*(1.0-pHRU->GetSoilProps(m)->stone_frac);  //stonefree thickness
    sum+=DD[m]*VV.rel_rootden;
  }
  for (m=0;m<pModel->GetNumSoilLayers();m++)
  {
    //cout<<"  height: "<<VV.height<<"  max_h: "<<C.max_height<<"  max_r: "<<R.max_root_length<<"  root den: "<<VV.rel_rootden<<"  length: "<<VV.root_length<<endl;
    if ((VV.rel_rootden<SMALL_ROOT_DENS) || (VV.root_length<SMALL_ROOT_LENGTH))
    {
      VV.cowan_alpha    =1e20;
      VV.root_resistance=HUGE_RESIST;
    }
    else
    {
      double rootfrac_m;        //[0..1] fraction of total root length in layer m
      double rootdens_m;        //[mm/mm3] root density for layer,

      rootfrac_m = VV.rel_rootden * (DD[m]/sum);//this is a calculation for each layer!
      rootdens_m = rootfrac_m/MM_PER_METER/MM_PER_METER*VV.root_length / DD[m]; //(mm/m^2)*(m/mm)*(m/mm)/(mm)=(mm/mm^3)


      ///< \ref Cowan, I.R. 1965. Transport of water in the soil-plant-atmosphere system. J Appl Ecol 2:221-239. \cite Cowan1965JoAE
      delT=PI*pow(rootradius,2)*rootdens_m;//(mm^2*mm/mm^3)=[-]
      double alpha=(1.0/(8.0*PI*rootdens_m))*(delT-3.0-2*(log(delT))/ (1.0-delT));//dimensionless

      VV.cowan_alpha    =alpha/KPA_PER_MPA*DENSITY_WATER*GRAVITY/DD[m];
      VV.root_resistance=(1-xylem_frac)*VV.resistivity/rootfrac_m;

      //RV[m].cowan_alpha    =alpha/KPA_PER_MPA*DENSITY_WATER*GRAVITY/DD[m];
      //RV[m].root_resistance=(1-R.xylem_frac)*VV.resistivity/rootfrac_m;

    }
  }
}


