/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ------------------------------------------------------------------
  Glacier Ice Melting
  Glacier Routing
  ----------------------------------------------------------------*/
#include "HydroProcessABC.h"
#include "GlacierProcesses.h"

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the glacier ice melt constructor
/// \param melt_type [in] Model of glacial melting selected
//
CmvGlacierMelt::CmvGlacierMelt(glacial_melt_type melt_type)://should just be out index?
  CHydroProcessABC(GLACIER_MELT)
{
  type =melt_type;

  if ((type==GMELT_HBV) || (type==GMELT_SIMPLE_MELT))
  {
    CHydroProcessABC::DynamicSpecifyConnections(1);

    iFrom[0]=pModel->GetStateVarIndex(GLACIER_ICE);
    iTo  [0]=pModel->GetStateVarIndex(GLACIER);
  }
  else if (type==GMELT_UBC)
  {
    CHydroProcessABC::DynamicSpecifyConnections(2);
    iFrom[0]=pModel->GetStateVarIndex(GLACIER_ICE);
    iTo  [0]=pModel->GetStateVarIndex(PONDED_WATER);
    iFrom[1]=pModel->GetStateVarIndex(GLACIER_CC);
    iTo  [1]=pModel->GetStateVarIndex(GLACIER_CC);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the glacier ice melt default destructor
//
CmvGlacierMelt::~CmvGlacierMelt(){}

//////////////////////////////////////////////////////////////////
/// \brief Initializes glacial melt class
//
void CmvGlacierMelt::Initialize(){}


//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for glacial melt algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by glacial melt algorithm (size of aP[] and aPC[])
//
void CmvGlacierMelt::GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const
{
  if (type==GMELT_SIMPLE_MELT)
  {
    nP=0;
  }
  else if (type==GMELT_HBV)
  {
    nP=1;
    aP[0]="HBV_MELT_GLACIER_CORR"; aPC[0]=CLASS_LANDUSE;
  }
  else if (type==GMELT_UBC)
  {
    nP=0;
  }
  else
  {
    ExitGracefully("CmvGlacierMelt::GetParticipatingParamList: undefined glacial melt algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief returns list of state variables which are used by glacial melt algorithm
///
/// \param melt_type [in] Algorithm for glacial melt used
/// \param *aSV [out] Array of state variable types needed by glacial melt algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by glacial melt algorithm (size of aSV[] and aLev[] arrays)
//
void CmvGlacierMelt::GetParticipatingStateVarList(glacial_melt_type melt_type,sv_type *aSV, int *aLev, int &nSV)
{
  if ((melt_type==GMELT_SIMPLE_MELT) || (melt_type==GMELT_HBV))
  {
    nSV=2;
    aSV[0]=GLACIER_ICE;  aLev[0]=DOESNT_EXIST;
    aSV[1]=GLACIER;      aLev[1]=DOESNT_EXIST;
  }
  else if (melt_type==GMELT_UBC)
  {
    nSV=3;
    aSV[0]=GLACIER_ICE;  aLev[0]=DOESNT_EXIST;
    aSV[1]=GLACIER;      aLev[1]=DOESNT_EXIST;
    aSV[2]=GLACIER_CC;   aLev[2]=DOESNT_EXIST;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rates of melting of glacial ice [mm/day]
/// \details  if type==GMELT_DEGREE_DAY
///             <ul> <li> melt is proportional to temperature </ul>
/// if type==GMELT_UBC
///             <ul> <li> melt is proportional to storage; constant of proportionality depends upon surface snow/ice </ul>
///
/// \param *state_var [in] Array of current state variables in HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which rates of change are calculated
/// \param *rates [out] Array of rates of change of modified state variables
///
//
void CmvGlacierMelt::GetRatesOfChange( const double             *state_var,
                                       const CHydroUnit  *pHRU,
                                       const optStruct   &Options,
                                       const time_struct &tt,
                                       double      *rates) const
{
  if (pHRU->GetHRUType()!=HRU_GLACIER){return;}

  //----------------------------------------------------------------------
  if (type==GMELT_SIMPLE_MELT)
  {
    rates[0]=pHRU->GetForcingFunctions()->potential_melt;
  }
  //----------------------------------------------------------------------
  else if (type==GMELT_HBV)
  {
    if (state_var[pModel->GetStateVarIndex(SNOW)]>REAL_SMALL){return;}//melt doesnt proceed until snow is gone

    rates[0]=pHRU->GetSurfaceProps()->HBV_melt_glacier_corr*threshPositive(pHRU->GetForcingFunctions()->potential_melt);
  }
  //----------------------------------------------------------------------
  else if (type==GMELT_UBC)
  {
    double CC           =state_var[iFrom[1]]; //glacier cold content [mm]
    double snow_coverage= pHRU->GetSnowCover(); //snow cover
    double pot_melt;                          //potential melt [mm/d]

    pot_melt =pHRU->GetForcingFunctions()->potential_melt;
    pot_melt*=pHRU->GetForcingFunctions()->subdaily_corr;

    if (pot_melt>0.0) {pot_melt*=(1.0-snow_coverage);}//glacier doesnt melt if snow is on top
    else              {pot_melt=0.0;} //issue: cold content for glaciers will never increase! (needed for emulation)

    //Reduce/increase cold content (always leads to pot_melt>=0, CC>=0)
    if (CC>pot_melt*Options.timestep){
      CC      -=pot_melt*Options.timestep;
      pot_melt=0;
    }
    else{
      pot_melt-=CC/Options.timestep;
      CC      =0.0;
    }

    //decay of cold content
    double P0CTKG=0.0;//basically sets CC=0 for this problem
    CC *= (1 - (1 / (1 + P0CTKG)));

    rates[0]=pot_melt;                                  //rates[0]: glacier ->ponded water
    rates[1]=(CC-state_var[iFrom[1]])/Options.timestep; //rates[1]: CC modification
  }
};

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow cannot drain "from" compartment over timestep
///
/// \param *state_var [in] Reference to set of participating state variables
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which rates of change are calculated
/// \param *rates [out] Corrected rates of state variable change
//
void   CmvGlacierMelt::ApplyConstraints( const double           *state_var,
                                         const CHydroUnit  *pHRU,
                                         const optStruct         &Options,
                                         const time_struct &tt,
                                         double      *rates) const
{
  if (pHRU->GetHRUType()!=HRU_GLACIER){return;}

  if (rates[0]<0.0){rates[0]=0.0;}//positivity constraint on melt

  //cant remove more than is there
  //rates[0]=threshMin(rates[0],state_var[iFrom[0]]/Options.timestep,0.0);
}



/*****************************************************************
   Glacier Release Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the glacier release constructor
/// \param mtype [in] Model of glacial release selected
//
CmvGlacierRelease::CmvGlacierRelease(glacial_release_type r_type):
  CHydroProcessABC(GLACIER_RELEASE)
{
  type =r_type;

  if ((type == GRELEASE_HBV_EC) || (type == GRELEASE_LINEAR) || (type == GRELEASE_LINEAR_ANALYTIC))
  {
    CHydroProcessABC::DynamicSpecifyConnections(1);

    iFrom[0]=pModel->GetStateVarIndex(GLACIER);
    iTo  [0]=pModel->GetStateVarIndex(SURFACE_WATER);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvGlacierRelease::~CmvGlacierRelease(){}

//////////////////////////////////////////////////////////////////
/// \brief Initializes glacial melt class
//
void CmvGlacierRelease::Initialize(){}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for glacial release algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by glacial release algorithm (size of aP[] and aPC[])
//
void CmvGlacierRelease::GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const
{
  if (type==GRELEASE_HBV_EC)
  {
    nP=3;
    aP[0]="HBV_GLACIER_KMIN";   aPC[0]=CLASS_LANDUSE;
    aP[1]="GLAC_STORAGE_COEFF"; aPC[1]=CLASS_LANDUSE;
    aP[2]="HBV_GLACIER_AG";     aPC[2]=CLASS_LANDUSE;
  }
  //-----------------------------------------------------------------
  else if (type == GRELEASE_LINEAR || type == GRELEASE_LINEAR_ANALYTIC)
  {
    nP=1;
    aP[0]="GLAC_STORAGE_COEFF"; aPC[0]=CLASS_LANDUSE;
  }
  else
  {
    ExitGracefully("CmvGlacierRelease::GetParticipatingParamList: undefined glacial release algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to state variables which are used by glacial release algorithm
///
/// \param melt_type [in] Model of glacial release used
/// \param *aSV [out] Array of state variable types needed by glacial release algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by glacial release algorithm (size of aSV[] and aLev[] arrays)
//
void CmvGlacierRelease::GetParticipatingStateVarList(glacial_release_type r_type,sv_type *aSV, int *aLev, int &nSV)
{
  nSV=1;
  aSV[0]=GLACIER;       aLev[0]=DOESNT_EXIST;

  if (r_type==GRELEASE_HBV_EC)
  {
    nSV=3;
    aSV[1]=SNOW;          aLev[1]=DOESNT_EXIST;
    aSV[2]=SNOW_LIQ;      aLev[2]=DOESNT_EXIST;
    aSV[3]=PONDED_WATER;  aLev[3]=DOESNT_EXIST;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rates of release of glacial water [mm/d] and other modified state variables
/// \details if type==GRELEASE_HBV_EC,
/// <ul> <li> melt is proportional to storage; constant of proportionality depends upon surface snow/ice </ul>
///
/// \param *storage [in] Reference to storage from which ice is released
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which rates of release are calculated
/// \param *rates [out] Array of rates of change in modified state variables
///
//
void CmvGlacierRelease::GetRatesOfChange( const double            *storage,
                                          const CHydroUnit  *pHRU,
                                          const optStruct     &Options,
                                          const time_struct &tt,
                                          double            *rates) const
{
  if (pHRU->GetHRUType()!=HRU_GLACIER){return;}

  double glacier_stor=storage[pModel->GetStateVarIndex(GLACIER)];

  //-----------------------------------------------------------------
  if (type==GRELEASE_HBV_EC)
  {
    double snow        =storage[pModel->GetStateVarIndex(SNOW)];
    double liq_snow    =storage[pModel->GetStateVarIndex(SNOW_LIQ)];

    double Kmin=pHRU->GetSurfaceProps()->HBV_glacier_Kmin;
    double Kmax=pHRU->GetSurfaceProps()->glac_storage_coeff;
    double Ag  =pHRU->GetSurfaceProps()->HBV_glacier_Ag;

    double K=Kmin+(Kmax-Kmin)*exp(-Ag*(snow+liq_snow));

    rates[0]=K*glacier_stor;
  }
  //-----------------------------------------------------------------
  else if (type == GRELEASE_LINEAR)
  {
    double K=pHRU->GetSurfaceProps()->glac_storage_coeff;

    rates[0]=K*glacier_stor;
  }
  else if (type == GRELEASE_LINEAR_ANALYTIC)
  {
    double K = pHRU->GetSurfaceProps()->glac_storage_coeff;
    rates[0] = glacier_stor*(1 - exp(-K*Options.timestep)) / Options.timestep; //Alternate analytical formulation
  }
};

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow cannot drain glacier over timestep
///
/// \param *storage [in] Reference to storage from which ice is released
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which rates of change are calculated
/// \param *rates [out] Corrected rates of change of state variables
//
void   CmvGlacierRelease::ApplyConstraints( const double                  *storage,
                                            const CHydroUnit  *pHRU,
                                            const optStruct   &Options,
                                            const time_struct &tt,
                                            double      *rates) const
{
  if (pHRU->GetHRUType()!=HRU_GLACIER){return;}

  if (rates[0]<0.0){rates[0]=0.0;}//positivity constraint
}





/*****************************************************************
   Glacier Infiltration Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the glacier infiltration constructor
/// \param mtype [in] Model of glacial infiltration selected
//
CmvGlacierInfil::CmvGlacierInfil(glacial_infil_type i_type):
  CHydroProcessABC(GLACIER_INFIL)
{
  type =i_type;

  if (type==GINFIL_UBCWM)
  {
    CHydroProcessABC::DynamicSpecifyConnections(3);

    iFrom[0]=pModel->GetStateVarIndex(PONDED_WATER);    iTo  [0]=pModel->GetStateVarIndex(GLACIER);
    iFrom[1]=pModel->GetStateVarIndex(PONDED_WATER);    iTo  [1]=pModel->GetStateVarIndex(SOIL,2);
    iFrom[2]=pModel->GetStateVarIndex(PONDED_WATER);    iTo  [2]=pModel->GetStateVarIndex(SOIL,3);

  }
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvGlacierInfil::~CmvGlacierInfil(){}

//////////////////////////////////////////////////////////////////
/// \brief Initializes glacial infiltration class
//
void CmvGlacierInfil::Initialize(){}

//////////////////////////////////////////////////////////////////
/// \brief Returns participating parameter list
///
/// \param *aP [out] array of parameter names needed for glacial release algorithm
/// \param *aPC [out] Class type (soil, vegetation, landuse or terrain) corresponding to each parameter
/// \param &nP [out] Number of parameters required by glacial release algorithm (size of aP[] and aPC[])
//
void CmvGlacierInfil::GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const
{
  if (type==GINFIL_UBCWM)
  {
    nP=6;
    aP[0]="IMPERMEABLE_FRAC";      aPC[0]=CLASS_LANDUSE;
    aP[1]="MAX_PERC_RATE";         aPC[1]=CLASS_SOIL;
    aP[2]="UBC_INFIL_SOIL_DEF";    aPC[2]=CLASS_SOIL;
    aP[3]="POROSITY";              aPC[3]=CLASS_SOIL;
    aP[4]="UBC_GW_SPLIT";          aPC[4]=CLASS_GLOBAL;
    aP[5]="UBC_FLASH_PONDING";     aPC[5]=CLASS_GLOBAL;
  }
  else
  {
    ExitGracefully("CmvGlacierInfil::GetParticipatingParamList: undefined glacial release algorithm",BAD_DATA);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to state variables which are used by glacial infiltration algorithm
///
/// \param melt_type [in] Model of glacial release used
/// \param *aSV [out] Array of state variable types needed by glacial release algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by glacial release algorithm (size of aSV[] and aLev[] arrays)
//
void CmvGlacierInfil::GetParticipatingStateVarList(glacial_infil_type r_type,sv_type *aSV, int *aLev, int &nSV)
{
  nSV=1;


  if (r_type==GINFIL_UBCWM)
  {
    nSV=2;
    aSV[0]=PONDED_WATER;   aLev[0]=DOESNT_EXIST;
    aSV[1]=SOIL;           aLev[1]=2;
    aSV[2]=SOIL;           aLev[2]=3;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief moves some portion of ponded water on a glacier to groundwater storage
///
/// \param *storage [in] Reference to storage from which ice is released
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which rates of release are calculated
/// \param *rates [out] Array of rates of change in modified state variables
///
//
void CmvGlacierInfil::GetRatesOfChange ( const double               *state_vars,
                                         const CHydroUnit  *pHRU,
                                         const optStruct     &Options,
                                         const time_struct &tt,
                                         double            *rates) const
{
  if (pHRU->GetHRUType()!=HRU_GLACIER){return;}

  //-----------------------------------------------------------------
  if (type==GINFIL_UBCWM)
  {
    double to_GW,runoff;
    double b2;

    double ponded       =state_vars[iFrom[0]];
    double rain_and_melt=ponded/Options.timestep;

    double max_perc_rate=pHRU->GetSoilProps(1)->max_perc_rate;
    double P0DSH        =CGlobalParams::GetParams()->UBC_GW_split;

    //.NET
    //to_GW=0;
    //runoff=rain_and_melt;

    //RFS:
    b2=g_debug_vars[0];//RFS Cheat - uses b2 calculated for adjacent soil
    to_GW       =min(max_perc_rate,rain_and_melt)*(1.0-b2); //overflows to GW
    runoff      =rain_and_melt*b2+max(rain_and_melt*(1-b2)-to_GW,0.0);

    rates[0]=runoff;              //ponded_water->glacier
    rates[1]=(1.0-P0DSH)*to_GW;   //ponded_water->upper grounwater
    rates[2]=(    P0DSH)*to_GW;   //ponded_water->lower groundwater
  }
};

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow cannot drain glacier over timestep
///
/// \param *storage [in] Reference to storage from which ice is released
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which rates of change are calculated
/// \param *rates [out] Corrected rates of change of state variables
//
void   CmvGlacierInfil::ApplyConstraints( const double            *storage,
                                          const CHydroUnit  *pHRU,
                                          const optStruct   &Options,
                                          const time_struct &tt,
                                          double      *rates) const
{
  if (pHRU->GetHRUType()!=HRU_GLACIER){return;}

  if (rates[0]<0.0){rates[0]=0.0;}//positivity constraint on runoff

  //cant remove more than is there
  rates[0]=threshMin(rates[0],storage[iFrom[0]]/Options.timestep,0.0);
}
