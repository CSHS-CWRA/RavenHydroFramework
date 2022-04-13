/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2022 the Raven Development Team
  ----------------------------------------------------------------
  class definitions:
  CmvBaseflow
  CmvSoilEvap
  CmvInterflow
  CmvPercolation
  CmvCapillaryRise
  CmvDrain
  CmvRecharge
  CmvSoilBalance
  ----------------------------------------------------------------*/

#ifndef SOILWATERMOVERS_H
#define SOILWATERMOVERS_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"
/******************************************************************
   HYDROLOGICAL PROCESSES : CODING CONVENTIONS
-------------------------------------------------------------------
Each hydrological process should store all constants related only
to its (global) functioning
All units should be in mm, MJ/m2, mg/m2, and days except *within* RateOfChange
routines, where we can locally use other units
******************************************************************/

///////////////////////////////////////////////////////////////////
/// \brief Methods for modelling baseflow
enum baseflow_type
{
  BASE_CONSTANT,        ///< Constant baseflow method
  BASE_LINEAR,          ///< simple bucket model (HBV,PRMS,UBCWM,...)
  BASE_LINEAR_CONSTRAIN,///< simple bucket model but limited down to FC
  BASE_LINEAR_ANALYTIC, ///< simple bucket model, analytical sol'n over timestep
  BASE_VIC,             ///< VIC baseflow method
  BASE_TOPMODEL,        ///< TOPMODEL Baseflow method
  BASE_POWER_LAW,       ///< Power Law saturation
  BASE_GR4J,            ///< GR4J Baseflow method
  BASE_THRESH_POWER,    ///< power law saturation above threshold
  BASE_THRESH_STOR      ///< linear storage above threshold (HBV-Lite)
};

////////////////////////////////////////////////////////////////////
/// \brief Data abstraction for loss of water from soil/groundwater to surface water
//
class CmvBaseflow: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  baseflow_type  type; ///< Model of baseflow selected

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvBaseflow(baseflow_type btype,
              int           from_index);
  ~CmvBaseflow();

  //inherited functions
  void Initialize();
  void GetRatesOfChange(const double              *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double      *rates) const;
  void ApplyConstraints(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double      *rates) const;

  void        GetParticipatingParamList   (string  *aP , class_type *aPC , int &nP) const;
  static void GetParticipatingStateVarList(baseflow_type btype,
                                           sv_type *aSV, int *aLev, int &nSV);

};

///////////////////////////////////////////////////////////////////
/// \brief Methods of modelling evaporation from multi-layered soil to atmosphere
//
enum soilevap_type
{
  SOILEVAP_GAWSER,        ///< uses GAWSER approach
  SOILEVAP_FEDERER,       ///< uses Federer 1979 resistance calculations \ref Federer 1979 \cite federer1979WRR
  SOILEVAP_ROOTFRAC,      ///< linear relation between ET and tension storage, distributed by root fraction
  SOILEVAP_VIC,           ///< Variable Infiltration Capacity model
  SOILEVAP_TOPMODEL,      ///< linear relation between ET and tension storage
  SOILEVAP_SEQUEN,        ///< Sequential soil evaporation method for FUSE emulation - VIC ONLY
  SOILEVAP_ROOT,          ///< Root weighting soil evaporation method for FUSE emulation - VIC ONLY
  SOILEVAP_ROOT_CONSTRAIN,///< same as ROOT, but top layer constrained to be above sat_wilt
  SOILEVAP_HBV,           ///< Simple HBV model (Bergstrom, 1996) -linear relation between ET and tension storage, with snow correction
  SOILEVAP_HYPR,          ///< HYPR model with ponded area correction (Ahmed et al, 2020)
  SOILEVAP_UBC,           ///< UBCWM Model (Quick, 1996)
  SOILEVAP_CHU,           ///< Ontario Crop Heat Unit method
  SOILEVAP_PDM,           ///< From Probabilty Distributed Model (Moore, 1985)
  SOILEVAP_GR4J,          ///< GR4J model approach (Perrin et al., 2003)
  SOILEVAP_LINEAR,        ///< AET a linear function of soil moisture
  SOILEVAP_SACSMA,        ///< Sacramento Soil Moisture Accounting algorithm (should only be used with SOILBAL_SACSMA)
  SOILEVAP_AWBM,          ///< Australia Water Balance Model of Boughton (1993)
  SOILEVAP_ALL            ///< AET==PET
};
////////////////////////////////////////////////////////////////////
/// \brief Data abstraction of loss of water from multiple soil layers to atmosphere
/// \details Uses the root-weighting method
//
// defined in SoilEvaporation.cpp
class CmvSoilEvap: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  soilevap_type   type;         ///< Model of soil evaporation selected
  int            *soil_ind;     ///< array of soil indices
  int             nSoilLayers;  ///< number of soil layers subject to evaporation

  void FedererSoilEvap           (const double      &PET,
                                  const double            *storage,
                                  const CHydroUnit  *pHRU,
                                  const optStruct         &Options,
                                  const time_struct &tt,
                                  double      *rates) const;

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvSoilEvap(soilevap_type se_type);
  ~CmvSoilEvap();

  //inherited functions
  void Initialize();
  void GetRatesOfChange(const double              *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double      *rates) const;
  void ApplyConstraints(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double      *rates) const;

  void        GetParticipatingParamList   (string  *aP , class_type *aPC , int &nP) const;
  static void GetParticipatingStateVarList(soilevap_type se_type,
                                           sv_type *aSV,        int *aLev, int &nSV);

};

///////////////////////////////////////////////////////////////////
/// \brief Methods of modeling interflow
//
enum interflow_type
{
  INTERFLOW_PRMS                                 ///< PRMS inteflow model \ref defined in Clark et al 2007 \cite Clark2008WRR
};

////////////////////////////////////////////////////////////////////
/// \brief Data abstraction for lateral loss of water from soil to surface wayer
// Defined in Interflow.cpp
//
class CmvInterflow: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  interflow_type  type; ///< Model of interflow selected

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvInterflow( interflow_type  itype,
                int             from_index);
  ~CmvInterflow();

  //inherited functions
  void Initialize();
  void GetRatesOfChange(const double              *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double      *rates) const;
  void ApplyConstraints(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double      *rates) const;

  void        GetParticipatingParamList   (string  *aP , class_type *aPC , int &nP) const;
  static void GetParticipatingStateVarList(interflow_type       itype,
                                           sv_type *aSV, int *aLev, int &nSV);
};



///////////////////////////////////////////////////////////////////
/// \brief Method of modeling percolation between soil layers
//
enum perc_type
{
  PERC_GAWSER,    ///< percolation method used in GAWSER (Schroeter, 19..) \cite Schroeter1988
  PERC_GAWSER_CONSTRAIN,///< percolation method used in GAWSER limited down to FC
  PERC_POWER_LAW,       ///< percolation method used in VIC (clark et al., 2007), HBV (soil to fast res) \cite Clark2008WRR
  PERC_PRMS,                    ///< percolation methods used in PRMS (clark et al., 2007)
  PERC_SACRAMENTO,///< percolation method in Sacremento  (Clark et al., 2007)
  PERC_LINEAR,    ///< Linear storage approach
  PERC_LINEAR_ANALYTIC, ///< Linear storage approach, analytical sol'n over timestep
  PERC_CONSTANT,  ///< constant percolation rate (e.g., HBV)
  PERC_GR4J,      ///< percolation method from GR4J model
  PERC_GR4JEXCH,  ///< groundwater exchange from GR4J model
  PERC_GR4JEXCH2, ///< groundwater exchange from GR4J model
  PERC_ASPEN      ///< percolation to aspen trees (S. Grass, 2018)
};

////////////////////////////////////////////////////////////////////
/// \brief Data abstration of loss of water from one soil layer to a lower soil layer
//
class CmvPercolation: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  perc_type       type;        ///< Model of percolation selected
  int             nSoilLayers; ///< number of soil layers subject to percolation

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvPercolation(perc_type      p_type,
                 int            In_indices,                     //soil water storage
                 int            Out_index);
  ~CmvPercolation();

  //inherited functions
  void Initialize();
  void GetRatesOfChange(const double              *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double      *rates) const;
  void ApplyConstraints(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double      *rates) const;

  void        GetParticipatingParamList   (string  *aP , class_type *aPC , int &nP) const;
  static void GetParticipatingStateVarList(perc_type    p_type,
                                           sv_type *aSV, int *aLev, int &nSV);
};

///////////////////////////////////////////////////////////////////
/// \brief Methods of modeling capillary rise
//
enum crise_type
{
  CRISE_HBV   ///< HBV cappillary rise
};

////////////////////////////////////////////////////////////////////
/// \brief Data abstraction of capillary rise
/// \details Calculates loss of water from soil layers to upper soil layers
//
class CmvCapillaryRise: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  crise_type      _type;        ///< Model of capillary rise selected

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvCapillaryRise(crise_type   cr_type,
                   int                          In_index,                       //soil water storage
                   int                          Out_index);
  ~CmvCapillaryRise();

  //inherited functions
  void Initialize();
  void GetRatesOfChange(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                              double      *rates) const;
  void ApplyConstraints(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double            *rates) const;

  void        GetParticipatingParamList   (string  *aP , class_type *aPC , int &nP) const;
  static void GetParticipatingStateVarList(crise_type   cr_type,
                                           sv_type *aSV, int *aLev, int &nSV);
};

///////////////////////////////////////////////////////////////////
/// \brief Method of modeling drainage soil layers
//
enum drain_type
{
  DRAIN_CONDUCTANCE,    ///< drain method based upon MODFLOW approach using a river bed conductance layer
  DRAIN_UFR             ///< drain method using upscaled flux relationships
};

////////////////////////////////////////////////////////////////////
/// \brief Data abstration of loss of water from one soil layer to a lower soil layer
//
class CmvDrain: public CHydroProcessABC
{  
  private:/*------------------------------------------------------*/
		drain_type					type; ///< Model of drainage selected
		int					 nSoilLayers; ///< number of soil layers subject to drainage

  public:/*-------------------------------------------------------*/
		//Constructors/destructors:
		CmvDrain(drain_type	d_type);			
		~CmvDrain();

		//inherited functions
    void Initialize();
    void GetRatesOfChange(const double		  *state_vars, 
								          const CHydroUnit  *pHRU, 
								          const optStruct	  &Options,
								          const time_struct &tt,
                                double      *rates) const;
    void ApplyConstraints(const double      *state_vars,
											    const CHydroUnit  *pHRU, 
								          const optStruct	  &Options,
								          const time_struct &tt,
                                double      *rates) const;

    void        GetParticipatingParamList   (string  *aP , class_type *aPC , int &nP) const;
    static void GetParticipatingStateVarList(drain_type	d_type, sv_type *aSV, int *aLev, int &nSV);
};
///////////////////////////////////////////////////////////////////
/// \brief Method of modeling recharge to aquifers
//
enum recharge_type
{
  RECHARGE_FROMFILE,             ///< uses recharge from data 
  RECHARGE_CONSTANT,         ///< constant recharge method applied to aquifers
  RECHARGE_CONSTANT_OVERLAP, ///<constant recharge method applied to aquifers with area weighted separation to connected gw cells 
};
////////////////////////////////////////////////////////////////////
/// \brief Data abstraction of capillary rise
/// \details Calculates loss of water from soil layers to upper soil layers
//
class CmvRecharge: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  recharge_type	      _type;        ///< Model of recharge 

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvRecharge(recharge_type	rech_type, int to_index, int junk); //junk just to distinguish constructors
  CmvRecharge(recharge_type	rech_type, int nConns);			
  ~CmvRecharge();

  //inherited functions
  void Initialize();
  void GetRatesOfChange(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double            *rates) const;
  void ApplyConstraints(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double            *rates) const;

  void        GetParticipatingParamList   (string  *aP , class_type *aPC , int &nP) const;
  static void GetParticipatingStateVarList(recharge_type	r_type,sv_type *aSV, int *aLev, int &nSV);
};

///////////////////////////////////////////////////////////////////
/// \brief Methods for modelling soil balance
enum soilbal_type
{
  SOILBAL_SACSMA       ///< Sacramento Soil Moisture Accounting Model
};

////////////////////////////////////////////////////////////////////
/// \brief Data abstraction for loss of water from soil/groundwater to surface water
//
class CmvSoilBalance : public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  soilbal_type  _type; ///< Model of soil balance selected

public:/*-------------------------------------------------------*/
       //Constructors/destructors:
  CmvSoilBalance(soilbal_type sb_type);
  ~CmvSoilBalance();

  //inherited functions
  void Initialize();
  void GetRatesOfChange(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double            *rates) const;
  void ApplyConstraints(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double            *rates) const;

  void        GetParticipatingParamList   (string  *aP,class_type *aPC,int &nP) const;
  static void GetParticipatingStateVarList(soilbal_type sb_type,
                                           sv_type *aSV,int *aLev,int &nSV);
};
#endif
