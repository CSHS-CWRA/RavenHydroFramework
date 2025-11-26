/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2025 the Raven Development Team
  ----------------------------------------------------------------*/
#ifndef LATCONNECT_H
#define LATCONNECT_H

#include "RavenInclude.h"
#include "ModelABC.h"
#include "Model.h"
#include "HydroUnits.h"
#include "HydroProcessABC.h"

class CModel;

///////////////////////////////////////////////////////////////////
/// \brief Class to store lateral connection information
//
class CLatConnect
{
public:
  long long int _source_hru_id;    //< HRU_ID of the source HRU
  long long int _recipient_hru_id; //< HRU_ID of the recipient HRU
  double        _weight;           //< Weight of the connection

  // Constructor
  CLatConnect(long long int hru, long long int connected_hru, double w);

  // Destructor
  ~CLatConnect();

  // Accessors
  int GetNumLatConnections() const;
  long long int GetSourceHRUID() const;
  long long int GetRecipientHRUID() const;
  double GetWeight() const;

  // Static methods to interact with the collection
  static int          GetNumLatConnections(const CModel *pModel);
  static CLatConnect* GetLatConnection(const CModel *pModel, int index);
  static bool CheckConnectionWeights (const CLatConnect* const connections[],
                                      const int nConnections,
                                      const CModel* pModel,
                                      const double tolerance = 0.05);
};


#endif // LATCONNECT_H
