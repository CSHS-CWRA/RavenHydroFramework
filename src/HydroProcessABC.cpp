/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2021 the Raven Development Team
  ----------------------------------------------------------------*/
#include "HydroProcessABC.h"
#include "StateVariables.h"

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the constructor for an abstract hydrological process
/// \param ptype [in] Type of hydrological process
//
CHydroProcessABC::CHydroProcessABC(const process_type ptype)
{
  _process=ptype;
  _nConnections=0;
  iFrom=NULL;
  iTo  =NULL;

  _pConditions=NULL;
  _nConditions=0;

  ExitGracefullyIf(pModel==NULL,
                   "CHydroProcessABC::Constructor:no model associated with hydrologic process",BAD_DATA);
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the constructor for an abstract hydrological process
/// \param ptype [in] Type of hydrological process
//
CHydroProcessABC::CHydroProcessABC(const process_type ptype, const int nConns)
{
  _process  	  =ptype;
  _nConnections=nConns;
  iFrom=NULL;
  iTo  =NULL;

  _pConditions=NULL;
  _nConditions=0;

  ExitGracefullyIf(pModel==NULL,
    "CHydroProcessABC::Constructor:no model associated with hydrologic process",BAD_DATA);
}


//////////////////////////////////////////////////////////////////
/// \brief Implementation of the constructor for an abstract hydrological process with a single connection
/// \param ptype [in] Type of hydrological process
/// \param From_index [in] Index of state variable/ storage unit which (typically) loses water/energy/mass
/// \param To_index [in] Index of state variable/ storage unit which (typically) gain water/energy/mass
//
CHydroProcessABC::CHydroProcessABC(const process_type ptype,
                                   const int          From_index,
                                   const int          To_index)
{
  _process    =ptype;
  _nConnections=1;
  iFrom=new int [_nConnections];
  iTo  =new int [_nConnections];

  iFrom[0]    =From_index;
  iTo  [0]    =To_index;

  _pConditions=NULL;
  _nConditions=0;

  if ((iFrom[0]<0) || (iTo[0]<0)){
    ExitGracefully("CHydroProcessABC::Constructor:negative state variable index",BAD_DATA);}
  if ((iFrom[0]>=pModel->GetNumStateVars()) || (iTo[0]>=pModel->GetNumStateVars())){
    ExitGracefully("CHydroProcessABC::Constructor:invalid state variable index",BAD_DATA);}
  ExitGracefullyIf(pModel==NULL,
                   "CHydroProcessABC::Constructor:no model associated with hydrologic process",BAD_DATA);
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the constructor for an abstract hydrological process with multiple connections
///
/// \param ptype [in] Process type
/// \param *From_indices [in] Array of indices of state variable/ storage unit which (typically) loses water/energy/mass
/// \param *To_indices [in] Array of indices of state variable/ storage unit which (typically) gain water/energy/mass
/// \param nConnect [in] Number of connections (size of From_indices, To_indices)
//
CHydroProcessABC::CHydroProcessABC(const process_type ptype,
                                   const int         *From_indices,
                                   const int         *To_indices,
                                   const int          nConnect)
{
  _process    =ptype;
  _nConnections=nConnect;
  iFrom=new int [_nConnections];
  iTo  =new int [_nConnections];
  if ((From_indices!=NULL) && (To_indices!=NULL)){
    for (int q=0;q<_nConnections;q++)
    {
      iFrom[q]=From_indices[q];
      iTo  [q]=To_indices  [q];

      if ((iFrom[q]<0) || (iTo[q]<0)){
        ExitGracefully("CHydroProcessABC::Constructor:negative state variable index",BAD_DATA);}
      if ((iFrom[q]>=pModel->GetNumStateVars()) || (iTo[q]>=pModel->GetNumStateVars())){
        ExitGracefully("CHydroProcessABC::Constructor:invalid state variable index",BAD_DATA);}
    }
  }
  ExitGracefullyIf(pModel==NULL,
                   "CHydroProcessABC::Constructor:no model associated with hydrologic process",BAD_DATA);

  _pConditions=NULL;
  _nConditions=0;
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the destructor
//
CHydroProcessABC::~CHydroProcessABC()
{
  if (DESTRUCTOR_DEBUG){cout<<"  DELETING HYDROLOGIC PROCESS"<<endl;}
  delete [] iFrom;       iFrom      =NULL;
  delete [] iTo;         iTo        =NULL;
  for (int i=0;i<_nConditions;i++){delete _pConditions[i];} delete [] _pConditions; _pConditions=NULL;
}
/*****************************************************************
   Static Members
*****************************************************************/
CModelABC *CHydroProcessABC::pModel=NULL;

//////////////////////////////////////////////////////////////////
/// \brief Specify reference to surface water model pModel
/// \note used only in construction of model
/// \param *pM [in] Reference to model
//
void CHydroProcessABC::SetModel(CModelABC *pM)
{
  ExitGracefullyIf(pM    ==NULL,
                   "CHydroProcessABC::SetModel: NULL Model!",BAD_DATA);
  ExitGracefullyIf(pModel!=NULL,
                   "CHydroProcessABC::SetModel: Cannot specify more than one model in Raven",BAD_DATA);
  pModel=pM;
}

//////////////////////////////////////////////////////////////////
/// \brief Validate reference to model and iTo/iFrom connectivity of process
/// \remark Called before solution
//
void CHydroProcessABC::Initialize()
{
  ExitGracefullyIf(pModel==NULL,
                   "CHydroProcessABC::Initialize: NULL model!",BAD_DATA);
  for (int q=0;q<_nConnections;q++)
  {
    ExitGracefullyIf(iFrom[q]==DOESNT_EXIST,
                     "CHydroProcessABC::Initialize: From compartment doesnt exist",RUNTIME_ERR);
    ExitGracefullyIf(iTo  [q]==DOESNT_EXIST,
                     "CHydroProcessABC::Initialize: To compartment doesnt exist",RUNTIME_ERR);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Should return parameters that are not automatically generated by model
/// \remark Called prior to running model; does nothing for abstract base class, but is inherited by child classes
///
/// \param *aP Array of parameter names
/// \param *aPC Array of parameter types
/// \param &nP Number of parameters (size of aP, aPC)
//
void CHydroProcessABC::GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const
{
  nP=0;
}

/*****************************************************************
   Accessors
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Returns enumerated hydrological process type
/// \return Hydrological process type
//
process_type CHydroProcessABC::GetProcessType()     const{return _process;}

//////////////////////////////////////////////////////////////////
/// \brief Reserves memory (Creates variables) for dynamically generated to/from index arrays
/// \param nConnects [in] Number of connections
/// \remarks typically called from constructors of child classes
//
void CHydroProcessABC::DynamicSpecifyConnections(const int nConnects)
{
  delete [] iFrom; //if these have already been created
  delete [] iTo;
  _nConnections=nConnects;
  iFrom=new int [_nConnections];
  iTo  =new int [_nConnections];
  for (int q=0;q<_nConnections;q++)
  {
    iFrom[q]    =DOESNT_EXIST;
    iTo  [q]    =DOESNT_EXIST;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Returns the rates of change in state variables over timestep  (often, (iFrom[] storage) to another (iTo[] storage))
/// \note these rates are the rate of LOSS of "iFrom" Storage Units/state variables. This is equivalent to the gain of the "iTo" storage units/state variables
/// \note this is implemented for child classes, and is just a placeholder for this abstract base class
///
/// \param *state_vars [in] Array of state variables stored in current HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model option information
/// \param &tt [in] Current model time
/// \param *rates [out] Rates of water/energy moved between storage locations / rates of change of modified state variables (size: _nConnections)
//
void CHydroProcessABC::GetRatesOfChange(const double             *state_vars,
                                        const CHydroUnit *pHRU,
                                        const optStruct    &Options,
                                        const time_struct &tt,
                                        double     *rates) const
{
  for (int q=0;q<_nConnections;q++)
  {
    rates[q]=0.0;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Corrects rates of change (*rates) returned from RatesOfChange function
/// \details Performed to maintain constraints on state_variables.
/// \note this is implemented for child classes, and is just a placeholder for this abstract base class
///
/// \param *state_vars [in] Array of state variables stored in current HRU
/// \param *pHRU [in] Reference to pertinent HRU
/// \param &Options [in] Global model option information
/// \param &tt [in] Current model time
/// \param *rates [out] Rates of water/energy moved between storage locations / rates of change of modified state variables (size: _nConnections)
//
void CHydroProcessABC::ApplyConstraints(const double *state_vars,
                                        const CHydroUnit *pHRU,
                                        const optStruct  &Options,
                                        const time_struct &tt,
                                        double     *rates) const
{
  //default - no constraints
  return;
}

//////////////////////////////////////////////////////////////////
/// \brief Adds conditional statement required for hydrological process to be applied
///
/// \param basis [in] Condition basis
/// \param compare_method [in] Method of comparison
/// \param data [in] Data for comparison
//
void CHydroProcessABC::AddCondition( condition_basis basis,
                                     comparison      compare_method,
                                     string          data)
{
  condition *pCO=new condition;
  pCO->basis         =basis;
  pCO->compare_method=compare_method;
  pCO->data          =data;
  if (!DynArrayAppend((void**&)(_pConditions),(void*)(pCO),_nConditions)){
    ExitGracefully("CHydroProcessABC::AddCondition: adding NULL condition",BAD_DATA);}
}
//////////////////////////////////////////////////////////////////
/// \brief redefines state variable index of recieving compartment 
///
/// \param toSVindex [in] state variable index of original storage compartment 
/// \param newToSVindex [in] state variable index of new storage compartment (assumed valid)
//
void CHydroProcessABC::Redirect(const int toSVindex,const int newToSVindex) 
{
  bool found=false;
  for(int q=0;q<_nConnections;q++)
  {
   if (iTo[q]    ==toSVindex){iTo[q]=newToSVindex;found=true;}
  }
  if(!found) {
    WriteWarning("CHydroProcessABC::Redirect: :-->RedirectFlow command refers to state variable not used by process",false);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Tests conditional statements to determine whether this specific process should be applied in this particular HRU
/// \note based upon HRU type or other diagnostics
/// \param *pHRU [in] Reference to pertinent HRU
/// \return True if comparison is successful
//
bool CHydroProcessABC::ShouldApply(const CHydroUnit *pHRU) const
{

  if (!pHRU->IsEnabled()) { return false; }
  if (_nConditions==0   ) { return true;  }
  for (int i=0;i<_nConditions;i++)
  {
    if (_pConditions[i]->basis==BASIS_HRU_TYPE)
    {
      HRU_type typ=pHRU->GetHRUType();
      HRU_type typ2=HRU_STANDARD;

      typ2=StringToHRUType(_pConditions[i]->data);

      if      ((_pConditions[i]->compare_method==COMPARE_IS_EQUAL ) && (typ!=typ2)){return false;}
      else if ((_pConditions[i]->compare_method==COMPARE_NOT_EQUAL) && (typ==typ2)){return false;}
    }
    else if (_pConditions[i]->basis==BASIS_HRU_GROUP)
    {
      int global_k=pHRU->GetGlobalIndex();
      bool is_in_group=pModel->IsInHRUGroup(global_k,_pConditions[i]->data);

      if ((!is_in_group) && (_pConditions[i]->compare_method==COMPARE_IS_EQUAL )){return false;}
      if (( is_in_group) && (_pConditions[i]->compare_method==COMPARE_NOT_EQUAL )){return false;}

    }
    else if (_pConditions[i]->basis==BASIS_LANDCLASS){
      //
      if ((_pConditions[i]->data!=pHRU->GetSurfaceProps()->landuse_name) && (_pConditions[i]->compare_method==COMPARE_IS_EQUAL )){return false;}
      if ((_pConditions[i]->data==pHRU->GetSurfaceProps()->landuse_name) && (_pConditions[i]->compare_method==COMPARE_NOT_EQUAL)){return false;}
    }
    else if (_pConditions[i]->basis==BASIS_VEGETATION){
      //
      if ((_pConditions[i]->data!=pHRU->GetVegetationProps()->vegetation_name) && (_pConditions[i]->compare_method==COMPARE_IS_EQUAL )){return false;}
      if ((_pConditions[i]->data==pHRU->GetVegetationProps()->vegetation_name) && (_pConditions[i]->compare_method==COMPARE_NOT_EQUAL)){return false;}
    }
  }
  return true;
}
