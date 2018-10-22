/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------
  class definitions:
  CmvConvolution
  ----------------------------------------------------------------*/

#ifndef CONVOLUTION_H
#define CONVOLUTION_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"

const int MAX_CONVOL_STORES=50;

///////////////////////////////////////////////////////////////////
/// \brief Methods for modelling a convolution
enum convolution_type
{
  CONVOL_TRIANGLE,      ///< Triangular UH
  CONVOL_GR4J_1,        ///< GR4J UH type 1
  CONVOL_GR4J_2,        ///< GR4J UH type 2
  CONVOL_GAMMA,         ///< Gamma distribution unit hydrograph - uses params alpha_Gamma, beta_Gamma 
  CONVOL_GAMMA_2        ///< Gamma distribution unit hydrograph 2 - uses params alpha_Gamma2, beta_Gamma2
};

////////////////////////////////////////////////////////////////////
/// \brief Data abstraction for loss of water from soil/groundwater to surface water
//
class CmvConvolution: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  convolution_type  _type;        ///< Model of baseflow selected

  int               _iTarget;     ///< state variable index of outflow target


  static int        _nConv;       /// # of CONVOLUTION variables (a.k.a. processes) in model

  void GenerateUnitHydrograph(const CHydroUnit *pHRU, const optStruct &Options, double *aUnitHydro, int &N) const;

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvConvolution(convolution_type type,
                 const int        to_index);
  ~CmvConvolution();

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
  static void GetParticipatingStateVarList(convolution_type btype,
                                           sv_type *aSV, int *aLev, int &nSV);

};
#endif
