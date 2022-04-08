/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2022 the Raven Development Team
  ------------------------------------------------------------------
  Flush (abstract move of all water from one compartment to another)
  Split (move water to two compartments)
  Overflow
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "SoilWaterMovers.h"

/*****************************************************************
   Flush Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the constructor
/// \param In_index [in] Index of the storage compartment from which water is flushed
/// \param Out_index [in] Index of the storage compartment to which water is flushed
//
CmvFlush::CmvFlush(int In_index,                       //soil water storage
                   int Out_index, const double &pct)
  :CHydroProcessABC(FLUSH,In_index,Out_index)
{
  _percentage=pct;
  ExitGracefullyIf(In_index==DOESNT_EXIST,
                   "CmvFlush Constructor: invalid 'from' compartment specified",BAD_DATA);
  ExitGracefullyIf(Out_index==DOESNT_EXIST,
                   "CmvFlush Constructor: invalid 'to' compartment specified",BAD_DATA);
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvFlush::~CmvFlush(){}

//////////////////////////////////////////////////////////////////
/// \brief Flush initialization
//
void   CmvFlush::Initialize()
{
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
/// \details In this case, user specifies 'from' and 'to' compartments - levels are not known before construction
///
/// \param *aSV [out] Reference to state variable types needed by flush algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by flush algorithm (size of aSV[] and aLev[] arrays)
//
void CmvFlush::GetParticipatingStateVarList(sv_type *aSV, int *aLev, int &nSV)
{
  nSV=0;
  //user specified 'from' & 'to' compartment, Levels - not known before construction
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rate of transfer of water during flushing [mm/d]
///
/// \param *storage [in] Array of HRU State variables
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of exchange between compartments [mm/day]

//
void   CmvFlush::GetRatesOfChange( const double      *storage,
                                   const CHydroUnit  *pHRU,
                                   const optStruct   &Options,
                                   const time_struct &tt,
                                         double     *rates) const
{
  rates[0]=_percentage*storage[iFrom[0]]/Options.timestep;
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow cannot drain "from" compartment over timestep
///
/// \param *storage [in] Array of all state variables for HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of exchange between compartments [mm/day]
//
void   CmvFlush::ApplyConstraints(const double           *storage,
                                  const CHydroUnit *pHRU,
                                  const optStruct      &Options,
                                  const time_struct &tt,
                                  double     *rates) const
{
  //cant remove more than is there (already built in)
  //exceedance of max "to" compartment presumed impossible
}


/*****************************************************************
   Split Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the split constructor
/// \param In_index [in] Index of the storage compartment from which water is flushed
/// \param Out_index1 [in] Index of the first storage compartment to which water is flushed
/// \param Out_index2 [in] Index of the second storage compartment to which water is flushed
//
CmvSplit::CmvSplit(int         In_index,                       //soil water storage
                   int         Out_index1,
                   int         Out_index2,
                   double      split_amt)
  :CHydroProcessABC(SPLIT)
{
  ExitGracefullyIf(In_index==DOESNT_EXIST,
                   "CmvFlush Constructor: invalid 'from' compartment specified",BAD_DATA);
  ExitGracefullyIf(Out_index1==DOESNT_EXIST,
                   "CmvFlush Constructor: invalid 'to' compartment specified",BAD_DATA);
  ExitGracefullyIf(Out_index2==DOESNT_EXIST,
                   "CmvFlush Constructor: invalid 'to' compartment specified",BAD_DATA);
  _split_pct=split_amt;

  DynamicSpecifyConnections(2);

  iFrom[0]=In_index;  iTo  [0]=Out_index1;
  iFrom[1]=In_index;  iTo  [1]=Out_index2;
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvSplit::~CmvSplit(){}

//////////////////////////////////////////////////////////////////
/// \brief Split initialization
//
void   CmvSplit::Initialize()
{
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
/// \details In this case, user specifies 'from' and 'to' compartments - levels are not known before construction
///
/// \param *aSV [out] Reference to state variable types needed by flush algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by flush algorithm (size of aSV[] and aLev[] arrays)
//
void CmvSplit::GetParticipatingStateVarList(sv_type *aSV, int *aLev, int &nSV)
{
  nSV=0;
  //user specified 'from' & 'to' compartments, Levels - not known before construction
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rate of transfer of water during flushing [mm/d]
///
/// \param *storage [in] Array of HRU State variables
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of exchange between compartments [mm/day]

//
void   CmvSplit::GetRatesOfChange( const double                 *storage,
                                   const CHydroUnit       *pHRU,
                                   const optStruct        &Options,
                                   const time_struct &tt,
                                   double     *rates) const
{
  rates[0]=(    _split_pct)*storage[iFrom[0]]/Options.timestep;
  rates[1]=(1.0-_split_pct)*storage[iFrom[0]]/Options.timestep;
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow cannot drain "from" compartment over timestep
///
/// \param *storage [in] Array of all state variables for HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of exchange between compartments [mm/day]
//
void   CmvSplit::ApplyConstraints(const double           *storage,
                                  const CHydroUnit *pHRU,
                                  const optStruct      &Options,
                                  const time_struct &tt,
                                  double     *rates) const
{
  //cant remove more than is there (already built in)
  //exceedance of max "to" compartment presumed impossible
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the constructor
/// \param In_index [in] Index of the storage compartment from which water overflows
/// \param Out_index [in] Index of the storage compartment to which water overflows
//
CmvOverflow::CmvOverflow(int                                    In_index,                       //soil water storage
                         int                                       Out_index)
  :CHydroProcessABC(OVERFLOW_PROC,In_index,Out_index)
{
  ExitGracefullyIf(In_index==DOESNT_EXIST,
                   "CmvOverflow Constructor: invalid 'from' compartment specified",BAD_DATA);
  ExitGracefullyIf(Out_index==DOESNT_EXIST,
                   "CmvOverflow Constructor: invalid 'to' compartment specified",BAD_DATA);
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvOverflow::~CmvOverflow(){}

//////////////////////////////////////////////////////////////////
/// \brief Flush initialization
//
void   CmvOverflow::Initialize()
{
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
/// \details In this case, user specifies 'from' and 'to' compartments - levels are not known before construction
///
/// \param *aSV [out] Reference to state variable types needed by overflow algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by overflow algorithm (size of aSV[] and aLev[] arrays)
//
void CmvOverflow::GetParticipatingStateVarList(sv_type *aSV, int *aLev, int &nSV)
{
  nSV=0;
  //user specified 'from' & 'to' compartment, Levels - not known before construction
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rate of loss of water during overflow [mm/d]
///
/// \param *state_var [in] Array of state variables
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of exchange between compartments [mm/day]
//
void   CmvOverflow::GetRatesOfChange( const double                      *state_var,
                                      const CHydroUnit       *pHRU,
                                      const optStruct        &Options,
                                      const time_struct &tt,
                                      double     *rates) const
{

  double max_storage=max(pHRU->GetStateVarMax(iFrom[0],state_var,Options),0.0);
  //cout << "Overflow "<<max_storage<<endl;
  rates[0]=max(state_var[iFrom[0]]-max_storage,0.0)/Options.timestep;
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow cannot drain "from" compartment over timestep
///
/// \param *state_var [in] Array of state variables
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of exchange between compartments [mm/day]
//
void   CmvOverflow::ApplyConstraints(const double                *state_var,
                                     const CHydroUnit *pHRU,
                                     const optStruct      &Options,
                                     const time_struct &tt,
                                     double     *rates) const
{
  //cant remove more than is there (already built in)
  //exceedance of max "to" compartment handled using other overflows
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the constructor
/// \param In_index [in] Index of the main storage compartment
/// \param Out_index [in] Index of the mixing storage compartment
//
CmvExchangeFlow::CmvExchangeFlow(int   In_index,                       //soil water storage
                                 int   Out_index)//soil water storage
  :CHydroProcessABC(EXCHANGE_FLOW,In_index,Out_index)
{
  ExitGracefullyIf(In_index==DOESNT_EXIST,
                   "CmvOverflow Constructor: invalid 'from' compartment specified",BAD_DATA);
  ExitGracefullyIf(Out_index==DOESNT_EXIST,
                   "CmvOverflow Constructor: invalid 'to' compartment specified",BAD_DATA);
  CHydroProcessABC::DynamicSpecifyConnections(2);
  iFrom[0]=In_index;  iTo[0]=Out_index;
  iFrom[1]=Out_index; iTo[1]=In_index;
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CmvExchangeFlow::~CmvExchangeFlow(){}

//////////////////////////////////////////////////////////////////
/// \brief Exchange Flow initialization
//
void   CmvExchangeFlow::Initialize()
{
  sv_type fromType      = pModel->GetStateVarType(iFrom[0]);
  ExitGracefullyIf(fromType!=SOIL,
                   "CmvExchangeFlow::Initialize: exchange flow must be from a soil unit",BAD_DATA);
  sv_type toType        = pModel->GetStateVarType(iTo[0]);
  ExitGracefullyIf((toType!=SOIL),
                   "CmvExchangeFlow::Initialize: exchange flow  must be between two soil units",BAD_DATA);
}

//////////////////////////////////////////////////////////////////
/// \brief Sets reference to participating state variables
/// \details In this case, user specifies 'from' and 'to' compartments - levels are not known before construction
///
/// \param *aSV [out] Reference to state variable types needed by overflow algorithm
/// \param *aLev [out] Array of level of multilevel state variables (or DOESNT_EXIST, if single level)
/// \param &nSV [out] Number of state variables required by overflow algorithm (size of aSV[] and aLev[] arrays)
//
void CmvExchangeFlow::GetParticipatingStateVarList(sv_type *aSV, int *aLev, int &nSV)
{
  nSV=0;
  //user specified 'from' & 'to' compartment, Levels - not known before construction
}

void CmvExchangeFlow::GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const
{
  nP=1;
  aPC[0]=CLASS_SOIL;
  aP [0]="EXCHANGE_FLOW";
}

//////////////////////////////////////////////////////////////////
/// \brief Returns rate of loss of water during overflow [mm/d]
///
/// \param *state_var [in] Array of state variables
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of exchange between compartments [mm/day]
//
void   CmvExchangeFlow::GetRatesOfChange( const double                  *state_var,
                                          const CHydroUnit       *pHRU,
                                          const optStruct        &Options,
                                          const time_struct &tt,
                                          double     *rates) const
{
  int m   = pModel->GetStateVarLayer(iFrom[0]); //uses mixing layer property
  rates[0] = pHRU->GetSoilProps(m)->exchange_flow;
  rates[1] = pHRU->GetSoilProps(m)->exchange_flow;
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Ensures that the rate of flow cannot drain "from" compartment over timestep
///
/// \param *state_var [in] Array of state variables
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Specified point at time at which this accessing takes place
/// \param *rates [out] Rate of exchange between compartments [mm/day]
//
void   CmvExchangeFlow::ApplyConstraints(const double            *state_var,
                                         const CHydroUnit *pHRU,
                                         const optStruct     &Options,
                                         const time_struct &tt,
                                         double     *rates) const
{
  //cant remove more than is there (already built in)
  //exceedance of max "to" compartment handled using other overflows
}
