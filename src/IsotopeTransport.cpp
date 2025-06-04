/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2025 the Raven Development Team
----------------------------------------------------------------
CIsotopeModel routines related to isotope transport

Raven transports is in terms of mass concentration, mg/L or mg/mm/m2 (as is is for all species)
this is assumed to be equivalent to mass concentrations e.g., mg18O/(mgO); units conversions are IMPLICIT
and is internally converted to composition for advection corrections used to represent enrichment
Raven reports these concentrations as compositions in the standard output
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "IsotopeTransport.h"

//////////////////////////////////////////////////////////////////
/// \brief isotope model constructor
//
CIsotopeModel::CIsotopeModel(CModel* pMod,CTransportModel* pTMod,string name,const int c)
  :CConstituentModel(pMod,pTMod,name,ISOTOPE,false,c)
{
  _isotope=ISO_O18;
  if      (name=="18O"){_isotope=ISO_O18;}
  else if (name=="2H" ){_isotope=ISO_H2; }
  else {
    ExitGracefully("CIsotopeModel constructor: Invalid isotope name",BAD_DATA_WARN);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief isotope model destructor
//
CIsotopeModel::~CIsotopeModel(){}

//////////////////////////////////////////////////////////////////
/// \brief initializes isotope transport model
//
void CIsotopeModel::Initialize(const optStruct& Options)
{
  // Initialize base class members
  //--------------------------------------------------------------------
  CConstituentModel::Initialize(Options);
}
//////////////////////////////////////////////////////////////////
/// \brief initializes isotope transport model
//
void CIsotopeModel::UpdateInitialConditions(const optStruct& Options)
{
  // Set surface water initial conditions
  //--------------------------------------------------------------------
  double Cinit= ConvertConcentration(_initSWconc);
  const CSubBasin *pBasin;
  int nSB=_pModel->GetNumSubBasins();
  for(int p=0;p<nSB;p++)
  {
    pBasin=_pModel->GetSubBasin(p);

    int nSegments=pBasin->GetNumSegments();
    
    for(int i=0; i<pBasin->GetNumSegments(); i++)
    {
      _aMout[p][i]=pBasin->GetOutflowArray()[nSegments-1] * SEC_PER_DAY * Cinit; // assumes nSegs=1
    }
    _aMout_last[p]=_aMout[p][0];

    const double *aQin=_pModel->GetSubBasin(p)->GetInflowHistory();
    for(int i=0; i<_nMinHist[p]; i++) { 
      _aMinHist [p][i]=aQin[i]*SEC_PER_DAY*Cinit;
    }
    const double *aQlat=_pModel->GetSubBasin(p)->GetLatHistory();
    for(int i=0; i<_nMlatHist[p];i++) {
      _aMlatHist[p][i]=aQlat[i]*SEC_PER_DAY*Cinit;
    }
    _aMlat_last[p]=_aMlatHist[p][0];
    _aMlocal   [p]=_aMlat_last[p];
    _aMlocLast [p]=_aMlat_last[p];
    
    //_channel_storage[p]=0.0;
    //_rivulet_storage[p]=0.0;

    if (pBasin->GetReservoir()!=NULL)
    {
      _aMres     [p]=Cinit*pBasin->GetReservoir()->GetStorage();
      _aMres_last[p]=_aMres[p];

      _aMsed     [p]=Cinit*_pModel->GetHydroUnit(pBasin->GetReservoir()->GetHRUIndex())->GetArea()*M2_PER_KM2*pBasin->GetReservoir()->GetLakebedThickness();
      _aMsed_last[p]=_aMsed[p];

      _aMout_res      [p]=Cinit*pBasin->GetReservoir()->GetOutflowRate()*SEC_PER_DAY;
      _aMout_res_last [p]=Cinit*pBasin->GetReservoir()->GetOldOutflowRate()*SEC_PER_DAY;
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief set initial surface water concentration
//
void CIsotopeModel::SetSurfaceWaterConc(const double& delta) 
{
  _initSWconc=delta;
}
//////////////////////////////////////////////////////////////////
/// \brief returns COMPOSITION of isotope in [o/oo] (rather than concentration in mg/mg)
/// \param M [in] mass [mg/m2]
/// \param V [in] volume [mm]
//
double CIsotopeModel::CalculateReportingConcentration(const double& mass,const double& vol) const
{
  if(fabs(vol)<1e-6) {return 0;}
  double C=CConstituentModel::CalculateReportingConcentration(mass,vol); //mg/L
  return ConcToComposition(C); //mg/L-> o/oo
}
//////////////////////////////////////////////////////////////////
/// \brief returns outflow composition of isotope in subbasin p at current point in time
/// \notes only used for reporting; calculations exclusively in terms of mass/energy
/// \param p [in] subbasin index
//
double CIsotopeModel::GetOutflowConcentration(const int p) const
{
  double C=CConstituentModel::GetOutflowConcentration(p); //mg/L
  if (C==0.0){return 0.0;}
  return ConcToComposition(C); //mg/L-> o/oo
}
//////////////////////////////////////////////////////////////////
/// \brief Converts basic composition units [o/oo] to [mg/mm/m2]
/// \param delta [in] composition, in o/oo
/// \returns concentration, in mg/mm/m2
//
double CIsotopeModel::ConvertConcentration(const double &delta) const
{
  return CompositionToConc(delta)*LITER_PER_M3/MM_PER_METER; //[o/oo]->[mg/mm-m2]
}

//////////////////////////////////////////////////////////////////
/// \brief returns advection correction factor for isotope
/// handles enrichment factor for evaporating waters
/// \param pHRU       [in] pointer to HRU of interest
/// \param iFromWater [in] index of "from" water storage state variable
/// \param iToWater   [in] index of "to" water storage state variable
/// \param mass       [in] source mass prior to advection [mg]
/// \param vol        [in] source volume prior to advection [mm] or [mm-m2]
/// \param Q          [in] flow between source and recipient [mm/d] or [mm-m2] 
///
double CIsotopeModel::GetAdvectionCorrection(const CHydroUnit* pHRU,const int iFromWater,const int iToWater,const double& mass, const double &vol, const double &Q) const
{
  sv_type fromType,toType;

  fromType=_pModel->GetStateVarType(iFromWater);
  toType  =_pModel->GetStateVarType(iToWater);

  if (iFromWater==DOESNT_EXIST){ //special flag for reservoir->atmosphere exchange

  }

  if (toType==ATMOSPHERE) // handles all enrichment due to evaporation
  {
    //From Stadnyk-Falcone, PhD Thesis, University of Waterloo 2008
    double dE;     // evaporative composition    [o/oo]
    double dL;     // lake/wetland composition   [o/oo]
    double dP;     // meteoric water composition [o/oo]
    double dA;     // atmospheric composition    [o/oo]
    double dStar;  // limiting steady state composition of source water (not used here) [o/oo]
    double alpha_star(1), ep_star, tmp(0.0);
    double ekin; //kinetic fractionation fraction
    double CK0(28.6);
    double eta=1.0;    //turbulence parameter = 1 for soil, 0.5 for open water, 0.6 for wetlands
    double hprime=1.0; //1.0 for small water bodies, 0.88 for large lakes/reservoirs

    //Forcing functions
    double h=pHRU->GetForcingFunctions()->rel_humidity;
    double T=pHRU->GetForcingFunctions()->temp_ave+ZERO_CELSIUS; //[K]

    //precipitation composition
    dP=_pTransModel->GetConcentration(pHRU->GetGlobalIndex(),_constit_index,ATMOS_PRECIP,0);

    //source composition
    double C_L=mass/vol*MM_PER_METER/LITER_PER_M3;//mg/m2/mm-->mg/L
    dL=ConcToComposition(C_L);//[mg/mg]->[o/oo]

    //from Horita and Wesolowski,  Liquid-vapor fractionation of oxygen and hydrogen isotopes of water from the freezing to the critical temperature,
    //Geochimica et Cosmochimica Acta, 58(16), 1994,
    if     (_isotope==ISO_O18) {
      alpha_star=exp(0.001*(-7.685+6.7123*(1e3/T)-1.6664*(1e6/T/T)+0.35401*(1e9/T/T/T))); //eqn 6 of H&W
      CK0=28.6/TO_PER_MILLE;   //[-] Vogt, 1976
    }
    else if(_isotope==ISO_H2) {
      alpha_star=exp(0.001*(1158.8*(T*T*T/1e9)-1620.1*(T*T/1e6)+794.84*(T/1e3)+2.9992*(1e9/T/T/T)-161.04));  //eqn 5 of H&W
      CK0=25.0/TO_PER_MILLE;
    }
    ep_star   =alpha_star-1.0; //[-]
    
    if      (fromType==DEPRESSION   ) {eta=0.6;}
    else if (fromType==SOIL         ) {eta=1.0;}
    else if (fromType==SURFACE_WATER) {eta=0.5;}

    ekin      =eta*hprime*CK0*(1.0-h);  //[-]
    ekin*=0.6; //to represent % of transpiration

    dA   =(( dP/TO_PER_MILLE-ep_star     )/alpha_star)*TO_PER_MILLE;//[o/oo]

    // from Gat & Levy (1978), Gat (1981)
    dStar=(h*dA/TO_PER_MILLE+ekin+ep_star/alpha_star)/(h - ekin - ep_star/alpha_star)*TO_PER_MILLE; //[o/oo]

    dE   =((dL/TO_PER_MILLE-ep_star)/alpha_star - h*dA/TO_PER_MILLE - ekin) / (1.0-h+ekin)*TO_PER_MILLE; //[o/oo]
    
    if (dE>dL){dE=dL;}

    double C_E  =CompositionToConc(dE);
    double C_lim=CompositionToConc(dStar);

    double tstep=_pModel->GetOptStruct()->timestep;
    double Cmax = (mass-C_lim*(vol-fabs(Q)*tstep))/fabs(Q)/tstep; //implicit conversion from mg/L to mg/mm-m2
    if (fabs(Q)          <1e-6){Cmax=ALMOST_INF;}
    if (fabs(vol-Q*tstep)<1e-6){return 1.0;}//prevents residual mass in empty volume

    //cout<<" CE,L,lim,max: "<<dE<<" "<<dL<<" "<<dStar<<" "<<ConcToComposition(Cmax)<<" "<<(Cmax>C_E)<<endl;

    if      (C_E<  Cmax){return min(Cmax/C_L,1.0); } //reduced enrichment so as not to exceed C_lim at end of time step
    else if (C_E>=C_lim){return 1.0;               } //no additional enrichment
    else                {return min(C_E /C_L,1.0); } //enrichment
  }

  else if ( ((toType==SNOW_LIQ    ) || (fromType==SNOW)) ||
            ((toType==PONDED_WATER) || (fromType==SNOW)) ) //snow melt - enriches snow
  {
    // \todo[funct]  minus refreeze offset
  }
  else if ( ((fromType==SNOW_LIQ) || (toType==SNOW)) ) //refreeze - dilutes snow
  {
    // \todo[funct] 
  }

  return 1.0;
}
//////////////////////////////////////////////////////////////////
/// \brief Write transport output file headers in .tb0 format
/// \details Called prior to simulation (but after initialization) from CModel::Initialize()
/// \param &Options [in] Global model options information
//
void CIsotopeModel::WriteEnsimOutputFileHeaders(const optStruct& Options)
{
  this->WriteOutputFileHeaders(Options); //Ensim format not supported by isotope model
}

//////////////////////////////////////////////////////////////////
/// \brief Writes minor transport output to ensim formatted file at the end of each timestep (or multiple thereof)
/// \note only thing this modifies should be output streams; called from CModel::WriteMinorOutput()
/// \param &Options [in] Global model options information
/// \param &tt [in] Local (model) time at the end of the pertinent time step
//
void CIsotopeModel::WriteEnsimMinorOutput(const optStruct& Options,const time_struct& tt)
{
  this->WriteMinorOutput(Options,tt);//Ensim format not supported by isotope model
}

//concentrations are 18O/(16O+18O) =%18O by mass [mg/mg] (<1)
//compositions   are delta values = δ=(R/RV-1)*10^3

//////////////////////////////////////////////////////////////////
/// \brief converts from concentration [mg/mg] to delta value (composition) [o/oo]
/// \param conc [in] [mg/mg]
//
double CIsotopeModel::ConcToComposition(const double &conc) const
{
  double RV;

  if (conc >= 1.0) {
    ExitGracefully("CIsotopeModel::ConcToComposition: invalid concentration",RUNTIME_ERR);
  }
  if(_isotope==ISO_O18) { RV=RV_VMOW_O18; }
  else                  { RV=RV_VMOW_H2;  }
  return (conc/(1-conc)/RV-1.0)*TO_PER_MILLE; //δ

}
//////////////////////////////////////////////////////////////////
/// \brief converts from delta value (composition) [o/oo]  to concentration [mg/mg] to
/// \param d [in] δ-value (composition) [o/oo]
//
double CIsotopeModel::CompositionToConc(const double& delta) const
{
  double RV;
  if(_isotope==ISO_O18) { RV=RV_VMOW_O18; }
  else                  { RV=RV_VMOW_H2;  }

  double R=((delta/TO_PER_MILLE)+1.0)*RV;

  return R/(R+1);
}
