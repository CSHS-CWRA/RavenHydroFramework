/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------*/
#ifndef LATERALEXCHANGE_H
#define LATERALEXCHANGE_H

#include "RavenInclude.h"
#include "ModelABC.h"
#include "Model.h"
#include "HydroUnits.h"
#include "HydroProcessABC.h"
class CModel;

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for lateral exchange process. 
/// Inherits from CHydroProcessABC
//
class CLateralExchangeProcessABC: public CHydroProcessABC
{
protected:/*------------------------------------------------------*/

  static int _nLatFlowProcesses;

  int      _nLatConnections; //< number of HRU lateral connections
  int     *_kFrom;           //< array of HRU from indices [size: nLatConnections] (JRC: Usually 1?)
  int     *_kTo;             //< array of HRU to indices [size: nLatConnections]
  int     *_iFromLat;        //< array of 'From' state variable indices [size: nLatConnections] 
  int     *_iToLat;          //< array of 'To' state variable indices [size: nLatConnections] 

  int      _LatFlowIndex;    //< global index of lateral flow process 

  static const CModel *_pModel;

  void DynamicSpecifyLatConnections(const int nLatConnects);

public:/*-------------------------------------------------------*/

  static void SetModel(const CModel *pM);

  CLateralExchangeProcessABC(const process_type ptype);     //multiple connection dynamic constructor
  ~CLateralExchangeProcessABC();
  
  int GetLateralFlowIndex() const;

  int GetNumLatConnections() const;

  const int *GetFromHRUIndices() const;
  const int *GetToHRUIndices() const;
  const int *GetLateralFromIndices() const;
  const int *GetLateralToIndices() const;

  virtual void GetParticipatingParamList(string *aP,class_type *aPC,int &nP) const{ nP=0;return; }

  void GetRatesOfChange(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                              double      *rates) const{return;}//default- no vertical transfer

  void ApplyConstraints(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                              double      *rates) const{return;};//default - no vertical transfer

  //
  virtual void GetLateralExchange(const double * const *state_vars, //array of all SVs for all HRUs, [k][i]
                                  const CHydroUnit * const *pHRUs,    
                                  const optStruct   &Options,
                                  const time_struct &tt,
                                        double      *exchange_rates) const=0;//purely virtual - required

};

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for the complete dump of water from one storage 
/// compartment in one HRU to another in another HRU
//
class CmvLatFlush: public CLateralExchangeProcessABC
{
private:/*------------------------------------------------------*/
  int _iFlushFrom; //< global state variable index of source state var
  int _iFlushTo;   //< global state variable index of target state var
  int _kk_from;    //< HRU group index of source HRUs
  int _kk_to;      //< HRU group index of target HRUs
  
  bool _constrain_to_SBs; // all transfer is within one sub-basin; otherwise, requires only one recipient HRU in model


public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvLatFlush(int   from_sv_ind,
              int     to_sv_ind,
              int   from_HRU_grp,
              int     to_HRU_grp,
              bool  constrain_to_SBs);
  ~CmvLatFlush();

  //inherited functions
  void Initialize();
 
  static void GetParticipatingStateVarList(sv_type *aSV, int *aLev, int &nSV);

  void           GetParticipatingParamList(string *aP,class_type *aPC,int &nP) const;

  void GetLateralExchange(const double * const *state_vars, //array of all SVs for all HRUs, [k][i]
                          const CHydroUnit * const *pHRUs,    
                          const optStruct   &Options,
                          const time_struct &tt,
                                double      *exchange_rates) const;//purely virtual

};

///////////////////////////////////////////////////////////////////
/// \brief Data abstraction for the complete dump of water from one storage 
/// compartment in one HRU to another in another HRU
//
class CmvLatEquilibrate : public CLateralExchangeProcessABC
{
private:/*------------------------------------------------------*/
  int                 _iSV;   //< global state variable index of state var
  int                  _kk;   //< HRU group index (or DOESNT_EXIST, if all HRUs should be equlibrated)
  double            *_Asum;   //< sum of areas of HRUs in HRU group within subbasin  [size: nSubBasins]
  int          * _halfconn;   //< 1/2 number of connections in basin b (=nHRUs in kk within basin-1) [size: nSubBasins]
  bool   _constrain_to_SBs;   //< all transfer is within one sub-basin; otherwise, requires only one recipient HRU in model
  double       _mixing_rate;  //< [0..1] % of water mixed per day

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvLatEquilibrate(int   sv_ind,
                    int   HRU_grp,
                    double mixing_rate,
                    bool  constrain_to_SBs);
  ~CmvLatEquilibrate();

  //inherited functions
  void Initialize();

  static void GetParticipatingStateVarList(sv_type* aSV, int* aLev, int& nSV);

  void           GetParticipatingParamList(string* aP, class_type* aPC, int& nP) const;

  void GetLateralExchange(const double* const* state_vars, //array of all SVs for all HRUs, [k][i]
                          const CHydroUnit* const* pHRUs,
                          const optStruct& Options,
                          const time_struct& tt,
                          double* exchange_rates) const;//purely virtual

};
#endif