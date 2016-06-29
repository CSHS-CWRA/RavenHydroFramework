/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright © 2008-2014 the Raven Development Team
------------------------------------------------------------------
class definitions:
  CmvDecay
----------------------------------------------------------------*/

#ifndef DECAY_H
#define DECAY_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"
#include "Model.h"
enum decay_type
{
  DECAY_BASIC,    /// < basic decay, calculated as -k*C
  DECAY_ANALYTIC  /// < analytic treatment of decay over finite time step
};

////////////////////////////////////////////////////////////////////
/// \brief Calculates the decay of a substance 
//
class CmvDecay: public CHydroProcessABC
{  
  private:/*------------------------------------------------------*/
    const CTransportModel* _pTransModel;
    
    decay_type _dtype;              ///< decay algorithm type
    int _constit_ind;              ///< index of constituent which is decaying
 
  public:/*-------------------------------------------------------*/
		//Constructors/destructors:
		CmvDecay(string constit_name, decay_type dtyp, CTransportModel *pTransportModel);
		~CmvDecay(); 

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
								          const time_struct &t,
                                double      *rates) const;
    
    void        GetParticipatingParamList   (string  *aP, class_type *aPC, int &nP) const;

};
#endif