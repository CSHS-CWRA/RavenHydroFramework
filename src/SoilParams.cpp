/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright © 2008-2014 the Raven Development Team
----------------------------------------------------------------*/
#include "Properties.h"
#include "SoilAndLandClasses.h"
//////////////////////////////////////////////////////////////////
/// \brief Returns Shuttleworth Wallace soil resistance to evaporation [d/mm]
/// \ref From Brook90 FRSS subroutine \cite Federer2010
/// \param &psi [in] Matrix potential [mm]
/// \param &S [in] Soil properties structure
//
double CSoilClass::CalcSoilResistance(const double			&psi,
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
double CSoilClass::CalcSoilHydCond(const double			&sat_liq, 
										 	 const double			&sat_ice, 
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
double CSoilClass::CalcSoilSaturation(const double			&psi,  
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
	//	psi=pS->clapp_m*MM_PER_CM*(sat_liq-pS->clapp_n)*(sat_liq-1.0);
	//}

	return psi;
}

/////////////////////////////////////////////////////////////////////
/// \brief Calculates wetted soil thermal conductivity based upon saturation & temperature
/// \ref adapted from Community Land Model Manual 3.0, Oleson et al., 2004 \cite Oleson2012 
///
/// \param &sat_liq [in] Saturation fraction of liquid H20 [0..1]
/// \param &sat_ice [in] Saturated fraction of ice [0..1]
/// \param &T [in] Temperature of soil layer [C]
/// \param &S [in] Soil properties structure
/// \return Soil thermal conductivity [W/m/K]
//
double CSoilClass::CalcSoilThermCond(const double			&sat_liq, 
												             const double			&sat_ice, 
												             const double			&T, 
												             const soil_struct &S)
{
	double therm_cond;
	double sat_cond;
	double dry_cond;
	double K;//kersten number
	double p	 =S.porosity;
	double rhob=S.bulk_density;

	//dry soil conductivity
	dry_cond=(0.135*rhob+64.7)/(2700-0.947*rhob); //CLM manual eqn. 6.62

	//saturated conductivity of water, ice & soil
	sat_cond=pow(S.thermal_cond,1.0-p)*pow(TC_WATER,p); //CLM manual eqn 6.60
	if (T<FREEZING_TEMP){
		sat_cond*=pow(TC_ICE,(1.0-sat_liq)*p);
	}

	//Kersten number
	K=sat_liq+sat_ice;
	if (T>FREEZING_TEMP){
		K = max(log(sat_liq+sat_ice)+1.0,0.0);
	}

	//final thermal conductivity
	therm_cond=dry_cond;
	if (sat_liq+sat_ice>1e-7){
		therm_cond=(K)*sat_cond+(1.0-K)*dry_cond;//CLM eqn 6.58
	}
	return therm_cond;
}

/////////////////////////////////////////////////////////////////////
/// \brief Calculates *volumetric* soil heat capacity [J/m^3/K]
/// \details This is just a simple weighted average of soil, ice and water capacities
/// 
/// \param &sat_liq [in] Saturation fraction of liquid H_{2}0 [0..1]
/// \param &sat_ice [in] Saturated fraction of ice [0..1]
/// \param *pS [in] Soil properties structure
/// \return Volumetric soil heat capacity [J/m^3/K]
//
double CSoilClass::CalcSoilHeatCap(const double			&sat_liq, 
											             const double			&sat_ice, 
											             const soil_struct *pS)
{
	return (1.0-pS->porosity)*pS->heat_capacity+
				 (sat_liq				  )*HCP_WATER+
				 (sat_ice				  )*HCP_ICE;
}
