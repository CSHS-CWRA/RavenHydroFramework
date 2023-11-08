/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2023 the Raven Development Team
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
  static bool       _smartmode;   ///< True if smart algorithm is used to reduce stores
  static int        _nStores;     ///< Number of convolution stores used in smart mode

  convolution_type  _type;        ///< Model of convolution selected

  int               _iTarget;     ///< state variable index of outflow target

  double LocalCumulDist(const double &t, const CHydroUnit *pHRU) const;

  void GenerateUnitHydrograph(const CHydroUnit *pHRU, const optStruct &Options, double *aUnitHydro, int *aIntervals, int &N) const;

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvConvolution(convolution_type type,
                 const int        to_index,
                 CModelABC* pModel,
                 int conv_index);
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

  void GetParticipatingParamList   (string  *aP , class_type *aPC , int &nP) const;
  static void GetParticipatingStateVarList(convolution_type btype,
                                           sv_type *aSV, int *aLev,
                                           int &nSV, int conv_index);

};
#endif
