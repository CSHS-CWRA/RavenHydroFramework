/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2026 the Raven Development Team
  ----------------------------------------------------------------
  class definitions:
  CResidualErrorModel
  ----------------------------------------------------------------*/

#ifndef REM_H
#define REM_H

#include "RavenInclude.h"
#include "Diagnostics.h"

enum meantype{ //mean type used in AR regression
  MEAN_LINEAR,
  MEAN_CONSTANT,
  MEAN_ZERO
};

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for the redistribution of snow based on slope and snow SWE

class CResidualErrorModel
{
private:/*------------------------------------------------------*/
  
  meantype        _meanType; 

  CTimeSeriesABC *_pQ95;
  CTimeSeriesABC *_pQ05; 
  CTimeSeriesABC *_pError;

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CResidualErrorModel();
  ~CResidualErrorModel();

  void CalculateREM         (CTimeSeriesABC  *pTSmod,
                             CTimeSeriesABC  *pTSObs,
                             CTimeSeriesABC  *pTSWeights,
                             const double    &starttime,
                             const double    &endtime,
                             const optStruct &Options) const;

  void WriteOutput          (const optStruct &Options) const;
};

#endif