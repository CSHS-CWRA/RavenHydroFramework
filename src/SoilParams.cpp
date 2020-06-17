/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2020 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Properties.h"
#include "SoilAndLandClasses.h"
//////////////////////////////////////////////////////////////////
/// \brief Returns Shuttleworth Wallace soil resistance to evaporation [d/mm]
/// \ref From Brook90 FRSS subroutine \cite Federer2010
/// \param &psi [in] Matrix potential [mm]
/// \param &S [in] Soil properties structure
//
double CSoilClass::CalcSoilResistance(const double                      &psi,
                                      const soil_struct &S)
{
  double psi_fc=CalcSoilPotential(S.field_capacity,0.0,S);
  return S.evap_res_fc*pow(psi/psi_fc,S.shuttleworth_b);
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates soil hydraulic conductivity
/// \details Calculates wetted soil hydraulic conductivity [mm/d] based upon saturation
/// \remark Uses Brooks-Corey method
///
/// \param &sat_liq [in] Saturation fraction of liquid H20 [0..1]
/// \param &sat_ice [in] Saturated fraction of ice [0..1]
/// \param &S [in] Soil properties structure
/// \return Wetted soil hydraulic conductivity [mm/d]
//
double CSoilClass::CalcSoilHydCond(const double                 &sat_liq,
                                   const double                 &sat_ice,
                                   const soil_struct &S)
{
  return S.hydraul_cond*pow(sat_liq,2.0*S.clapp_b+3.0);//Brooks-Corey
}

///////////////////////////////////////////////////////////////////
/// \brief Calculates wetted soil saturation based upon matric potential, psi
/// \param &psi [in] Soil matric potential [-mm]
/// \param &S [in] Soil properties structure
/// return Corresponding soil saturation [-]
//
double CSoilClass::CalcSoilSaturation(const double                      &psi,
                                      const soil_struct &S)
{
  return max(1.0,pow(psi/S.air_entry_pressure,-1.0/S.clapp_b));  //Brooks-Corey
}

/////////////////////////////////////////////////////////////////////
/// \brief Calculates soil matric potential [-mm] based on saturation \cite Clapp1978WRR
/// \param &sat_liq [in] Saturation fraction of liquid H20 [0..1]
/// \param &sat_ice [in] Saturated fraction of ice [0..1]
/// \param &S [in] Soil properties structure
/// \return Soil matric potential [-mm]
//
double CSoilClass::CalcSoilPotential( const double &sat_liq,
                                      const double &sat_ice,
                                      const soil_struct &S)
{
  static double psi;

  psi=S.air_entry_pressure*pow(sat_liq,-S.clapp_b);//Brooks-Corey

  //for alternate parabolic transition  near saturation (Brook90- Clapp-Hornberger 1978)
  //if (sat_liq>SAT_INF){
  //    psi=pS->clapp_m*MM_PER_CM*(sat_liq-pS->clapp_n)*(sat_liq-1.0);
  //}

  return psi;
}


