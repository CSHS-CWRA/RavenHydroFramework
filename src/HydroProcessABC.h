/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2019 the Raven Development Team
  ----------------------------------------------------------------*/

#ifndef HYDROPROCESS_H
#define HYDROPROCESS_H

#include "RavenInclude.h"
#include "ModelABC.h"
#include "HydroUnits.h"
#include "GlobalParams.h"

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for physical hydrological processes (abstract base class)
/// \details Data Abstraction for physical processes that move water or energy
///   from one set of storage units (e.g., soil, indexed with the iFrom storage
///   variable) to another set (e.g., surface water, indexed with the iTo
///   storage variable) within the same HRU.
///
///   Most processes have only one "from" unit and one "to" unit. However,
///   others (such as soil water flow) require solution of a PDE
///   and cannot truly be uncoupled. Thus, processes such as infiltration
///   require simultaneous (e.g., finite difference) solution for the
///   change in water storage in all units of the soil column
///
///   Processes may also change auxiliary state variables, such as snow depth.
///   These auxiliary state variables are not part of the mass or energy
///   balance calculations. In this case, water mass and/or energy is not
///   directly changed, so the "from" and "to" units are the same
///
///   The Hydro Processes are unaware of the overall model (CModel), but can
///   obtain properties from individual HRUs. A process cannot
///   modify the system: it merely reports to the solver the appropriate
///   amount of water or energy to move, i.e., the vector \f$ d/dt{\phi} \f$,
///   where \f$ \phi \f$ is the vector of state variables
///
/// \note If an additional hydrological process is added, the following routines must be
///   revised:
///- ParseMainInputFile: an additional parse case is required
//
class CHydroProcessABC
{
protected:/*----------------------------------------------------*/

  process_type     _process;  ///< Hydrological proccess, e.g., sublimation, infiltration, etc.
  int                *iFrom;  ///< indices of state variables/storage units that (typically) lose water/energy/mass (size: nConnections)
  int                  *iTo;  ///< indices of state variables/storage units that (typically) gain water/energy/mass (size: nConnections) (==iFrom if no mass/energy is exchanged in change of state variable)
  int         _nConnections;  ///< usually 1, number of transfer routes between state variables/storage units

  condition  **_pConditions;  ///< array of conditions for process to occur (ALL must be satisfied)
  int          _nConditions;  ///< number of conditions

  static  CModelABC *pModel;  ///< pointer to model

  void DynamicSpecifyConnections(const int nConnects);

public:/*-------------------------------------------------------*/
  //Constructors:
  CHydroProcessABC(const process_type ptype);     //multiple connection dynamic constructor
  CHydroProcessABC(const process_type ptype,
                   const int          nConns);     //multiple connection dynamic constructor
  CHydroProcessABC(const process_type ptype,      //single connection constructor
                   const int          In_index,
                   const int          Out_index);
  CHydroProcessABC(const process_type ptype,      //multiple connection constructor
                   const int         *In_indices,
                   const int         *Out_indices,
                   const int          nConnect);
  virtual ~CHydroProcessABC();

  //accessors (some inlined for speed)
  inline const int*    GetFromIndices()       const { return &iFrom[0]; } ///< Returns pointer to array of 'from' indices
  inline const int*    GetToIndices()         const { return &iTo[0];   } ///< Returns pointer to array of 'to' indices
  inline int           GetNumConnections()    const { return _nConnections; } ///< Returns total number of connections between state variables manipulated by process
  process_type         GetProcessType()       const;

  virtual int          GetNumLatConnections() const { return 0; }

  bool                 ShouldApply(const CHydroUnit*pHRU) const;
  //functions
  static void          SetModel    (CModelABC *pM);/// \todo [reorg]: should really not be accessible

  void                 AddCondition(condition_basis basis,
                                    comparison      compare_method,
                                    string          data);

  void                 Redirect(const int toSVindex, const int newToSVindex);

  virtual void Initialize();
  virtual void GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const=0;

  //calculates and returns rates of water/energy LOSS of "iFrom" Storage Units/state variables (e.g., [mm/d] or [MJ/m2/d])
  //This is the gain of the "iTo" storage units [mm/d] or [MJ/m2/d] . For changes to state variables
  //that do not involve the movement of water/energy/mass between units (i.e., density changes), the iFrom and
  //iTo units are the same.
  virtual void GetRatesOfChange(const double      *state_vars,
                                const CHydroUnit  *pHRU,
                                const optStruct   &Options,
                                const time_struct &tt,
                                      double      *rates) const=0;

  virtual void ApplyConstraints(const double      *state_vars,
                                const CHydroUnit  *pHRU,
                                const optStruct   &Options,
                                const time_struct &tt,
                                      double      *rates) const=0;
};

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for the complete dump of water from one storage compartment to another
//
class CmvFlush: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  double _percentage;

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvFlush(int   from_index,
           int     to_index, const double &pct);
  ~CmvFlush();

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

  static void GetParticipatingStateVarList(sv_type *aSV, int *aLev, int &nSV);
  void GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const{}
};

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for the complete dump of water from one storage compartment to two others
//
class CmvSplit: public CHydroProcessABC
{
private:/*------------------------------------------------------*/

  double _split_pct;  ///<percentage of water going to target compartment #1

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvSplit(int   from_index,
           int     to_index1,
           int     to_index2,
           double  split_amt);
  ~CmvSplit();

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

  static void GetParticipatingStateVarList(sv_type *aSV, int *aLev, int &nSV);
  void GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const{}
};


///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for the overflow of water from one storage compartment to another
//
class CmvOverflow: public CHydroProcessABC
{
private:/*------------------------------------------------------*/

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvOverflow(int        from_index,
              int     to_index);
  ~CmvOverflow();

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

  static void GetParticipatingStateVarList(sv_type *aSV, int *aLev, int &nSV);
  void GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const{}
};

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for exchange flow with a mixing zone
//
class CmvExchangeFlow: public CHydroProcessABC
{
private:/*------------------------------------------------------*/

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvExchangeFlow(int    from_index,
                  int   mixingzone_index);
  ~CmvExchangeFlow();

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

  static void GetParticipatingStateVarList(sv_type *aSV, int *aLev, int &nSV);
  void GetParticipatingParamList(string *aP, class_type *aPC, int &nP) const;
};
#endif
