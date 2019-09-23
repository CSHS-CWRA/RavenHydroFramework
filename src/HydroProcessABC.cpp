/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2019 the Raven Development Team
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

  _cascading  =false;
  _nCascades  =0;
  _iCascade   =NULL;

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

  _cascading  =false;
  _nCascades  =0;
  _iCascade   =NULL;
 
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

  _cascading  =false;
  _nCascades  =0;
  _iCascade   =NULL;

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

  _cascading  =false;
  _nCascades  =0;
  _iCascade   =NULL;

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
  delete [] _iCascade;  _iCascade   =NULL;
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
/// \brief Returns reference to array of 'from' indices
/// \return Reference to array of 'from' indices
//
const int*   CHydroProcessABC::GetFromIndices()     const{return &iFrom[0];}

//////////////////////////////////////////////////////////////////
/// \brief Returns reference to array of 'to' indices
/// \return Reference to array of 'to' indices
//
const int*   CHydroProcessABC::GetToIndices()       const{return &iTo[0];}

//////////////////////////////////////////////////////////////////
/// \brief Returns total number of connections between state variables manipulated by process
/// \return total number of connections between state variables manipulated by process
//
int          CHydroProcessABC::GetNumConnections()  const{return _nConnections+_nCascades;}

//////////////////////////////////////////////////////////////////
/// \brief Returns true if outflow cascades
/// \return True if outflow cascades
//
bool         CHydroProcessABC::HasCascade()         const{return _cascading;}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of cascade connections
/// \return Number of cascade connections
//
int          CHydroProcessABC::GetNumCascades()     const{return _nCascades;}

//////////////////////////////////////////////////////////////////
/// \brief Returns reference to list of indices of consecutive storage units from which cascade connections occur
/// \return Reference to list of indices of consecutive storage units from which cascade connections occur
//
int          CHydroProcessABC::GetCascadeFromIndex()const{return _iCascade[0];}//iFrom[0];}//default-works for single from/To- may wish to change

//////////////////////////////////////////////////////////////////
/// \brief Returns reference to list of indicies of consecutive storage units to which cascade connections occur
/// \return Reference to list of indicies of consecutive storage units to which cascade connections occur
//
const int*   CHydroProcessABC::GetCascadeToIndices()const{return &_iCascade[1];}

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
/// \param *rates [out] Rates of water/energy moved between storage locations / rates of change of modified state variables (size: _nConnections+_nCascades)
//
void CHydroProcessABC::GetRatesOfChange(const double             *state_vars,
                                        const CHydroUnit *pHRU,
                                        const optStruct    &Options,
                                        const time_struct &tt,
                                        double     *rates) const
{
  for (int q=0;q<_nConnections+_nCascades;q++)
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
/// \param *rates [out] Rates of water/energy moved between storage locations / rates of change of modified state variables (size: _nConnections+_nCascades)
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
/// \brief Tests conditional statements to determine whether this specific process should be applied in this particular HRU
/// \note based upon HRU type or other diagnostics
/// \param *pHRU [in] Reference to pertinent HRU
/// \return True if comparison is successful
//
bool CHydroProcessABC::ShouldApply(const CHydroUnit *pHRU) const
{

  if (_nConditions==0){return true;}
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

//////////////////////////////////////////////////////////////////
/// \brief Adds cascade to process
/// \param *indices [in] Reference to indicies of _cascading storage units
/// \param nIndices [in] Number of indicies in cascade
//
void CHydroProcessABC::AddCascade(const int *indices, const int nIndices)
{
  ExitGracefullyIf(_nConnections>1,
                   "CHydroProcessABC::AddCascade: cannot currently support _cascading for multi-connection processes",BAD_DATA);

  if (_iCascade!=NULL){
    WriteWarning("Cascade is being overwritten!",true);
    delete [] _iCascade; _iCascade=NULL;
  }
  _nCascades=nIndices-1;
  _iCascade = new int [nIndices];
  for (int i=0; i<nIndices; i++){
    _iCascade[i]=indices[i]; //first entry in _iCascade is original 'from' entry
    ExitGracefullyIf(_iCascade[i]==DOESNT_EXIST,
                     "CHydroProcessABC::AddCascade: invalid cascade compartment specified",BAD_DATA);

  }
  _cascading=true;

  int *iFromNew=new int [_nConnections+_nCascades];
  int *iToNew  =new int [_nConnections+_nCascades];
  for (int q=0;q<_nConnections;q++)
  {
    iFromNew[q]=iFrom[q];
    iToNew  [q]=iTo[q];
  }
  for (int q=0;q<_nCascades;q++)
  {
    iFromNew[_nConnections+q]=_iCascade[q];
    iToNew  [_nConnections+q]=_iCascade[q+1];
  }
  delete [] iFrom; iFrom=NULL; //a problem for some reason!
  delete [] iTo;   iTo=NULL;
  iFrom=iFromNew;
  iTo  =iToNew;
}

//////////////////////////////////////////////////////////////////
/// \brief Portions flow of material (typically water) through cascade of storage units until a unit with sufficient storage is reached
/// \details given a flow rate q1, into a _cascading, bucket-type system, portions this flow
/// downward through cascade until a receptacle with sufficient or infinite storage is reached \n
/// Example: \n
/// assume maxstorage =1mm for 3 buckets except the last, with infinite
/// storage. Each bucket starts out with 0.5 mm of water and a inflow rate of 5mm/day is
/// applied to the first bucket for tstep=1 day.
/// The distributed rates to each of the 4 buckets over the time step should be r[0]=0.5,
/// r[1]=0.5,r[2]=0.5;and r[3]=3.5 [mm/d]
///
/// \param *rates [out] Rates of change in state variable values modified by this process (size: _nConnections)
/// \param *state_vars [in] Array of state variables for this HRU
/// \param *maxstorage [in] maximum value of state variables (max storage in this case)
/// \param &tstep [in] Simulation time step
//
void CHydroProcessABC::Cascade(      double *rates,
                                     const double *state_vars,
                                     const double *maxstorage,
                                     const double &tstep)
{
  if (!_cascading){return;}

  double storage_remaining;
  int    q;
  q=_nConnections-1;//assumes last process connection is the one with overflow -may wish to correct
  for (int i=0;i<_nCascades;i++)
  {
    storage_remaining=maxstorage[_iCascade[i]]-state_vars[_iCascade[i]];
    if (rates[q]*tstep>storage_remaining)
    {
      rates[q+1]=rates[q]-storage_remaining/tstep;
      rates[q  ]=         storage_remaining/tstep;
    }
    q++;
  }
}
