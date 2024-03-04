/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2021 the Raven Development Team
----------------------------------------------------------------*/
#include "SurfaceEnergyExchange.h"

/*****************************************************************
Atmospheric Convection Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of convection constructor
/// \param pTransMod [in] pointer to transport model
/// \param iStorTo [in] state variable index of water storage corresponding to surface energy store
//
CmvPartitionEnergy::CmvPartitionEnergy(const CTransportModel *pTransMod, CModelABC *pModel)
  :CHydroProcessABC(PARTITION_ENERGY, pModel)
{
  _pTransModel=pTransMod;
  int cTemp=_pTransModel->GetConstituentIndex("TEMPERATURE");
  ExitGracefullyIf(cTemp==DOESNT_EXIST,
    "CmvPartitionEnergy Constructor - TEMPERATURE transport must be turned on for atmospheric convection algorithm to be used",BAD_DATA_WARN);

  int N=5; //5 storage compartment which may exchange with atmosphere (for now)

  CHydroProcessABC::DynamicSpecifyConnections(2*N);

  int layer_To;
  int layer_Atmos=_pTransModel->GetLayerIndex(cTemp,pModel->GetStateVarIndex(ATMOSPHERE)); //layer index of atmospheric enthalpy
  int layer_LH   =pModel->GetStateVarIndex(LATENT_HEAT);

  layer_To=_pTransModel->GetLayerIndex(cTemp,pModel->GetStateVarIndex(SOIL,0)); //layer index of soil[0]
  iFrom[0  ]=pModel->GetStateVarIndex(CONSTITUENT,layer_Atmos); //global index of atmospheric enthalpy
  iTo  [0  ]=pModel->GetStateVarIndex(CONSTITUENT,layer_To   ); //global index of topsoil enthalpy
  iFrom[N+0]=pModel->GetStateVarIndex(CONSTITUENT,layer_LH   ); //global index of latent heat losses
  iTo  [N+0]=pModel->GetStateVarIndex(CONSTITUENT,layer_To   ); //global index of topsoil enthalpy

  if(pModel->GetStateVarIndex(SNOW)!=DOESNT_EXIST) {
    layer_To=_pTransModel->GetLayerIndex(cTemp,pModel->GetStateVarIndex(SNOW)); //layer index of snow

    iFrom[1  ]=pModel->GetStateVarIndex(CONSTITUENT,layer_Atmos); //global index of atmospheric enthalpy
    iTo  [1  ]=pModel->GetStateVarIndex(CONSTITUENT,layer_To   ); //global index of snow enthalpy
    iFrom[N+1]=pModel->GetStateVarIndex(CONSTITUENT,layer_LH   ); //global index of latent heat losses
    iTo  [N+1]=pModel->GetStateVarIndex(CONSTITUENT,layer_To   ); //global index of snow enthalpy
  }
  if(pModel->GetStateVarIndex(DEPRESSION)!=DOESNT_EXIST) {
    layer_To=_pTransModel->GetLayerIndex(cTemp,pModel->GetStateVarIndex(DEPRESSION)); //layer index of depression storage

    iFrom[2  ]=pModel->GetStateVarIndex(CONSTITUENT,layer_Atmos); //global index of atmospheric enthalpy
    iTo  [2  ]=pModel->GetStateVarIndex(CONSTITUENT,layer_To   ); //global index of depression enthalpy
    iFrom[N+2]=pModel->GetStateVarIndex(CONSTITUENT,layer_LH   ); //global index of latent heat losses
    iTo  [N+2]=pModel->GetStateVarIndex(CONSTITUENT,layer_To   ); //global index of depression enthalpy
  }
  if(pModel->GetStateVarIndex(CANOPY)!=DOESNT_EXIST) {
    layer_To=_pTransModel->GetLayerIndex(cTemp,pModel->GetStateVarIndex(CANOPY)); //layer index of canopy

    iFrom[3  ]=pModel->GetStateVarIndex(CONSTITUENT,layer_Atmos); //global index of atmospheric enthalpy
    iTo  [3  ]=pModel->GetStateVarIndex(CONSTITUENT,layer_To);    //global index of canopy enthalpy
    iFrom[N+3]=pModel->GetStateVarIndex(CONSTITUENT,layer_LH   ); //global index of latent heat losses
    iTo  [N+3]=pModel->GetStateVarIndex(CONSTITUENT,layer_To   ); //global index of canopy enthalpy
  }
  if(pModel->GetStateVarIndex(CANOPY_SNOW)!=DOESNT_EXIST) {
    layer_To=_pTransModel->GetLayerIndex(cTemp,pModel->GetStateVarIndex(CANOPY)); //layer index of canopy snow

    iFrom[4  ]=pModel->GetStateVarIndex(CONSTITUENT,layer_Atmos); //global index of atmospheric enthalpy
    iTo  [4  ]=pModel->GetStateVarIndex(CONSTITUENT,layer_To);    //global index of canopy snow enthalpy
    iFrom[N+4]=pModel->GetStateVarIndex(CONSTITUENT,layer_LH   ); //global index of latent heat losses
    iTo  [N+4]=pModel->GetStateVarIndex(CONSTITUENT,layer_To   ); //global index of canopy snow enthalpy
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvPartitionEnergy::~CmvPartitionEnergy(){} //do nothing

//////////////////////////////////////////////////////////////////
/// \brief Initializes convection object
//
void CmvPartitionEnergy::Initialize() {} //do nothing


//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for abstraction algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by abstraction algorithm (size of aP[] and aPC[])
//
void        CmvPartitionEnergy::GetParticipatingParamList(string  *aP,class_type *aPC,int &nP) const
{
  nP=0;
  aP[nP]="CONVECTION_COEFF";  aPC[nP]=CLASS_LANDUSE;    nP++;
  aP[nP]="VEG_CONV_COEFF";    aPC[nP]=CLASS_VEGETATION; nP++;
  aP[nP]="MAX_DEP_AREA_FRAC"; aPC[nP]=CLASS_LANDUSE;    nP++; //TMP debug - may change with variable depression area support
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variable list
///
/// \param absttype [in] Selected abstraction model
/// \param *aSV [out] Reference to array of state variables needed by abstraction algorithm
/// \param *aLev [out] Array of levels of multilevel state variables (or DOESNT_EXIST of single level)
/// \param &nSV [out] Number of participating state variables (length of aSV and aLev arrays)
//
void CmvPartitionEnergy::GetParticipatingStateVarList(sv_type *aSV,int *aLev,int &nSV)
{
  //do nothing
};

//////////////////////////////////////////////////////////////////
/// \brief local routine to retrieve storage temperatures given enthalpy and water content
// assumes that very thin/nonexistent stores are equilibrated with the atmosphere, i.e., T=Tair
//
double CmvPartitionEnergy::GetStorTemp(int iEnth,int iWatStor,const double &Tair,const double *state_vars) const{
  double Hv,stor,hv;
 // _pTransModel->GetConcentration(k,iWatStor)
  if(iWatStor!=DOESNT_EXIST) {
    Hv   =state_vars[iEnth]; //enthalpy, [MJ/m2]
    stor =state_vars[iWatStor]; //water storage [mm]
    if(stor>PRETTY_SMALL) {
      hv   =Hv/(stor/MM_PER_METER); //volumetric specific enthalpy [MJ/m3]
    }
    return ConvertVolumetricEnthalpyToTemperature(hv);
  }
  return Tair;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns rates of energy transfer from the atmosphere to land surface
//
/// \param *state_vars [in] Array of state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time at which abstraction is to be calculated
/// \param *rates [out] rates of energy transfer from atmosphere to land and land to latent heat transfer [MJ/m2/d]
//
void CmvPartitionEnergy::GetRatesOfChange(const double      *state_vars,
                                          const CHydroUnit  *pHRU,
                                          const optStruct   &Options,
                                          const time_struct &tt,
                                                double      *rates) const
{
  //double Hv,hv,stor;
  double Tsoil,Tsnow,Tdep,Tcan,TcanS;

  double Tair=pHRU->GetForcingFunctions()->temp_ave;
  double Rnet=pHRU->GetForcingFunctions()->LW_radia_net+
              pHRU->GetForcingFunctions()->SW_subcan_net;

  double LAI        =pHRU->GetVegVarProps()->LAI;
  double conv_coeff =pHRU->GetSurfaceProps()->convection_coeff; // [MJ/m2/d/K]
  double vconv_coeff=pHRU->GetVegetationProps()->veg_conv_coeff*LAI; // [MJ/m2/d/K]

  int iSoil=pModel->GetStateVarIndex(SOIL,0);
  int iSnow=pModel->GetStateVarIndex(SNOW);
  int iDep =pModel->GetStateVarIndex(DEPRESSION);
  int iCan =pModel->GetStateVarIndex(CANOPY);
  int iCanS=pModel->GetStateVarIndex(CANOPY_SNOW);

  const CEnthalpyModel *pEnthalpyModel=_pTransModel->GetEnthalpyModel();
  Tsoil=pEnthalpyModel->GetWaterTemperature(state_vars,iSoil);
  Tsnow=pEnthalpyModel->GetWaterTemperature(state_vars,iSnow);
  Tdep =pEnthalpyModel->GetWaterTemperature(state_vars,iDep );
  Tcan =pEnthalpyModel->GetWaterTemperature(state_vars,iCan );
  TcanS=pEnthalpyModel->GetWaterTemperature(state_vars,iCanS);
  Tsnow=min(Tsnow,0.0);
  TcanS=min(TcanS,0.0);

  double pct_cover=1.0;

  //evaluate partial coverage of HRU
  double dep_cover   =pHRU->GetSurfaceProps()->max_dep_area_frac;  //pHRU->GetDepressionAreaFraction(); //TMP DEBUG - may have to vary with depression storage
  double snow_cover  =pHRU->GetSnowCover();
  double forest_cover=pHRU->GetSurfaceProps()->forest_coverage;
  //double sparseness  =pHRU->GetSurfaceProps()->forest_sparseness;

  pct_cover=max(1.0-snow_cover-dep_cover,0.0);

  double sv_fact=(1.0-forest_cover)*1.0; //TMP DEBUG- should include some skyview factor impacts - these are already implicitly in Rnet??

  double K =(1.0-exp(- conv_coeff*Options.timestep))/Options.timestep; //analytical correction factors ensure proper asymptotic approach to Tair (replace just conv_coeff)
  double Kv=(1.0-exp(-vconv_coeff*Options.timestep))/Options.timestep;

  rates[0]=pct_cover   *( K*(Tair-Tsoil) + Rnet*sv_fact);   // ATMOSPHERE -> SOIL
  rates[1]=snow_cover  *( K*(Tair-Tsnow) + Rnet*sv_fact);   // ATMOSPHERE -> SNOW
  rates[2]=dep_cover   *( K*(Tair-Tdep ) + Rnet*sv_fact);   // ATMOSPHERE -> DEPRESSION
  rates[3]=forest_cover*(Kv*(Tair-Tcan ) + Rnet        );   // ATMOSPHERE -> CANOPY
  rates[4]=forest_cover*(Kv*(Tair-TcanS) + Rnet        );   // ATMOSPHERE -> CANOPY_SNOW

  //  dHsdHw=-stor*rho_soil*c_soil*dT/dH_w;
  //  rates[10]-=dHsdHw*rates[0]; // SOIL -> SOIL_EQUILIBRIUM
  //  rates[0]*=(1.0-dHsdHw);

  //Handle evaporative energy losses
  //----------------------------------------------------------------------
  int    nAdvConnections=_pTransModel->GetNumAdvConnections();
  int    k              =pHRU->GetGlobalIndex();
  double ET,LH_loss,LH_loss_sub;
  int N=5; //# of compartments
  for(int q=0;q<nAdvConnections;q++)
  {
    int iToWater  =_pTransModel->GetToWaterIndex(q);
    if(pModel->GetStateVarType(iToWater)==ATMOSPHERE)
    {
      int iFromWater=_pTransModel->GetFromWaterIndex(q);
      int js        =_pTransModel->GetJsIndex(q);

      ET=pModel->GetFlux(k,js,Options)/MM_PER_METER; //[m/d]

      LH_loss    =ET*DENSITY_WATER*LH_VAPOR; //[MJ/m2/d]
      LH_loss_sub=ET*DENSITY_WATER*LH_SUBLIM;

      if      (iFromWater==iSoil) { rates[N+0]-=LH_loss;     } // SOIL ->LATENT_HEAT (soil evaporation)
      else if (iFromWater==iSnow) { rates[N+1]-=LH_loss_sub; } // SNOW ->LATENT_HEAT (sublimation)
      else if (iFromWater==iDep ) { rates[N+2]-=LH_loss;     } // DEPRESSION ->LATENT_HEAT (OW evaporation)
      else if (iFromWater==iCan ) { rates[N+3]-=LH_loss;     } // CANOPY ->LATENT_HEAT (Canopy evap)
      else if (iFromWater==iCanS) { rates[N+4]-=LH_loss_sub; } // CANOPY_SNOW ->LATENT_HEAT (canopy sublimation)

      //implicitly assumes that no energy comes from cooling down evaporated/sublimated water
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
///
/// \param *state_vars [in] Array of state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Current model time at which abstraction is to be calculated
/// \param *rates [out] rates[0]= rate of abstraction [mm/d]
//
void CmvPartitionEnergy::ApplyConstraints(const double      *state_vars,
                                          const CHydroUnit  *pHRU,
                                          const optStruct   &Options,
                                          const time_struct &tt,
                                                double      *rates) const
{
  //do nothing

  //JRC: may have to handle very thin amounts of water getting too much incoming energy
  // dry soil ?
}
