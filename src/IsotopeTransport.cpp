/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2022 the Raven Development Team
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
  else if (name=="2H" ){_isotope=ISO_H2;}
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
/// \param C [in] composition, in o/oo
/// \returns concentration, in mg/mm/m2
//
double CIsotopeModel::ConvertConcentration(const double &C) const
{
  return CompositionToConc(C)*LITER_PER_M3/MM_PER_METER; //[o/oo]->[mg/mm-m2]
}

//////////////////////////////////////////////////////////////////
/// \brief returns advection correction factor for isotope
/// handles enrichment factor for evaporating waters
/// \param pHRU       [in] pointer to HRU of interest
/// \param iFromWater [in] index of "from" water storage state variable
/// \param iToWater   [in] index of "to" water storage state variable
/// \param Cs         [in] isotopic concentration of source water [mg/mg]
///
double CIsotopeModel::GetAdvectionCorrection(const CHydroUnit* pHRU,const int iFromWater,const int iToWater,const double& Cs) const
{
  sv_type fromType,toType;

  fromType=_pModel->GetStateVarType(iFromWater);
  toType  =_pModel->GetStateVarType(iToWater);

  if (iFromWater==DOESNT_EXIST){ //special flag for reservoir->atmosphere exchange

  }

  if(toType==ATMOSPHERE) // handles all enrichment due to evaporation
  {
    //From Stadnyk-Falcone, PhD Thesis, University of Waterloo 2008
    double dE;     // evaporative composition    [o/oo]
    double dL;     // lake/wetland composition   [o/oo]
    double dP;     // meteoric water composition [o/oo]
    double dA;     // atmospheric composition    [o/oo]
    double dStar;  // limiting steady state composition of source water (not used here) [o/oo]

    double h=pHRU->GetForcingFunctions()->rel_humidity;
    double T=pHRU->GetForcingFunctions()->temp_ave+ZERO_CELSIUS; //[K]

    //TMP DEBUG - this bit should not be this complicated. Should have a _pTransModel->GetConcentration(k,c,ATMOS_PRECIP,m=0);
    //ATMOS_PRECIP should be a dirichlet condition!
    int iAtmPrecip=_pModel->GetStateVarIndex(ATMOS_PRECIP);
    int m         =_pTransModel->GetLayerIndex(_constit_index,iAtmPrecip);
    int i         =_pModel->GetStateVarIndex(CONSTITUENT,m);
    dP=_pTransModel->GetConcentration(pHRU->GetGlobalIndex(),i); //[o/oo]

    dL=ConcToComposition(Cs);//[mg/mg]->[o/oo] // is this necessary? I feel this should be specified in [o/oo]?

    double eta=1.0;        //turbulence parameter = 1 for soil, 0.5 for open water, 0.6 for wetlands
    double hprime=1.0;     //1.0 for small water bodies, 0.88 for large lakes/reservoirs
    double theta;          //advection correction
    if      (fromType==DEPRESSION   ) {eta=0.6;}
    else if (fromType==SOIL         ) {eta=1.0;}
    else if (fromType==SURFACE_WATER) {eta=0.5;}
    //else if (fromType==SNOW         )  //sublimation?
    //else if (fromType==CANOPY_SNOW  )  //sublimation?

    theta=(1-hprime)/(1-h);

    //from Horita and Wesolowski,  Liquid-vapor fractionation of oxygen and hydrogen isotopes of water from the freezing to the critical temperature,
    //Geochimica et Cosmochimica Acta, 58(16), 1994,
    double alpha_star,ep_star,tmp(0.0);
    double beta;
    double CK0(0.0);
    if     (_isotope==ISO_O18) {
      tmp=-7.685+6.7123*(1e3/T)-1.6664*(1e6/T/T)+0.35401*(1e9/T/T/T); //eqn 6 of H&W
      CK0=28.6/TO_PER_MILLE;   //[-] Vogt, 1976
    }
    else if(_isotope==ISO_H2) {
      tmp=1158.8*(T*T*T/1e9)-1620.1*(T*T/1e6)+794.84*(T/1e3)+2.9992*(1e9/T/T/T)-161.04;  //eqn 5 of H&W
      CK0=25.0/TO_PER_MILLE;
    }
    alpha_star=exp(0.001*tmp); //[-]
    ep_star   =alpha_star-1.0; //[-]
    beta      =eta*theta*CK0;  //[-]

    dA   =((  dP/TO_PER_MILLE-ep_star)/alpha_star)*TO_PER_MILLE;//[o/oo]

    //conversion - all compositions in o/oo must be converted to concentrations
    // from Gat & Levy (1978), Gat (1981)
    dStar=(h*dA/TO_PER_MILLE+beta*(1-h)+ep_star/alpha_star)/(h - beta*(1-h) - ep_star/alpha_star)*TO_PER_MILLE; //[o/oo]

    dE   =((dL/TO_PER_MILLE-ep_star)/alpha_star - h*dA/TO_PER_MILLE - beta*(1-h)) / ((1-h)*(1+beta))*TO_PER_MILLE; //[o/oo]

    double C_E=CompositionToConc(dE);
    double C_L=CompositionToConc(dL);

    //cout <<"C_E, C_L"<< C_E<<" "<<C_L<<endl;
    //return 1.0; //TMP DEBUG
    if (C_L==0.0){return 1.0;}
    return min(C_E/C_L,1.0);
  }

  else if ( ((toType==SNOW_LIQ    ) || (fromType==SNOW)) ||
            ((toType==PONDED_WATER) || (fromType==SNOW)) ) //snow melt - enriches snow
  {
    // minus refreeze offset

  }
  else if ( ((fromType==SNOW_LIQ) || (toType==SNOW)) ) //refreeze - dilutes snow
  {
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

//R-values       are O18/O16 mass ratio [mg/mg]
//concentrations are 18O/(16O+18O) =%18O by mass [mg/mg] (<1)
//compositions   are delta values = δ=(R/RV-1)*10^3

//////////////////////////////////////////////////////////////////
/// \brief converts from  R val (ratio of, e.g., 18O/16O) [mg/mg] to concentration [mg/mg]
/// \param R [in] R-value [mg/mg]
//
double CIsotopeModel::RvalToConcentration(const double &Rv) const
{
  return Rv/(Rv+1);
}
//////////////////////////////////////////////////////////////////
/// \brief converts from concentration [mg/mg] to R val (ratio of, e.g., 18O/16O) [mg/mg]
/// \param conc [in] [mg/mg]
//
double CIsotopeModel::ConcentrationToRval(const double &conc) const
{
  return conc/(1-conc);
}
//////////////////////////////////////////////////////////////////
/// \brief converts from concentration [mg/mg] to delta value (composition) [o/oo]
/// \param conc [in] [mg/mg]
//
double CIsotopeModel::ConcToComposition(const double &conc) const
{
  double RV;
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
  return RvalToConcentration(R);
}
