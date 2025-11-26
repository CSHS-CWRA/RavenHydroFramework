/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2025 the Raven Development Team
  ------------------------------------------------------------------
  lateral connection class
  implemented by J. Berg, University of Bern
  ----------------------------------------------------------------*/

#include "HydroProcessABC.h"
#include "LatConnect.h"


///////////////////////////////////////////////////////////////////
/// \brief Constructor for CLatConnect
//
CLatConnect::CLatConnect(long long int hru, long long int connected_hru, double w)
{
  _source_hru_id=hru;
  _recipient_hru_id=connected_hru;
  _weight=w;

  if (_source_hru_id < 0) {
    std::cerr << "CLatConnect constructor:: Invalid HRU ID: " << _source_hru_id << std::endl;
    ExitGracefully("CLatConnect constructor:: Invalid HRU ID", BAD_DATA);
  }
  if (_recipient_hru_id < 0) {
    std::cerr << "CLatConnect constructor:: Invalid Connected HRU ID: " << _recipient_hru_id << std::endl;
    ExitGracefully("CLatConnect constructor:: Invalid Connected HRU ID", BAD_DATA);
  }
  if (_weight < 0.0) {
    std::cerr << "CLatConnect constructor:: Invalid weight: " << _weight << std::endl;
    ExitGracefully("CLatConnect constructor:: Invalid weight", BAD_DATA);
  }
}

///////////////////////////////////////////////////////////////////
/// \brief Destructor for CLatConnect
//
CLatConnect::~CLatConnect()
{
  // Destructor logic if needed
}

///////////////////////////////////////////////////////////////////
/// \brief Returns the source HRU ID
/// \return  HRU ID
//
long long int CLatConnect::GetSourceHRUID() const
{
  return _source_hru_id;
}

///////////////////////////////////////////////////////////////////
/// \brief Returns the recipient HRU ID
/// \return  recipient HRU ID
//
long long int CLatConnect::GetRecipientHRUID() const
{
  return _recipient_hru_id;
}

///////////////////////////////////////////////////////////////////
/// \brief Returns the weight
/// \return Double weight
//
double CLatConnect::GetWeight() const
{
  return _weight;
}


///////////////////////////////////////////////////////////////////
/// \brief Checks if weights for each HRU add up to 1.0 within tolerance
/// \param connections [in] Array of lateral connections
/// \param nConnections [in] Number of connections
/// \param pModel [in] Pointer to model
/// \param tolerance [in] Acceptable deviation from 1.0 (default 0.05)
/// \return true if all HRU connection weights sum to approximately 1.0
//
bool CLatConnect::CheckConnectionWeights(const CLatConnect* const connections[],
                                         const int nConnections,
                                         const CModel* pModel,
                                         const double tolerance)
{
  if (connections == NULL || nConnections <= 0 || pModel == NULL) {
    return false;
  }

  // Get the number of HRUs in the model
  int nHRUs = pModel->GetNumHRUs();

  // Create arrays to store the sum of weights for each HRU
  double* sumWeights = new double[nHRUs];
  bool* hasConnections = new bool[nHRUs];

  // Initialize arrays
  for (int k = 0; k < nHRUs; k++) {
    sumWeights[k] = 0.0;
    hasConnections[k] = false;
  }

  // Calculate the sum of weights for each source HRU
  for (int i = 0; i < nConnections; i++)
{
    int sourceHRU = pModel->GetHRUByID(connections[i]->GetSourceHRUID())->GetGlobalIndex();
    double weight = connections[i]->GetWeight();

    sumWeights    [sourceHRU] += weight;
    hasConnections[sourceHRU] = true;
  }

  // Check if weights sum to approximately 1.0 for each HRU that has connections
  bool allValid = true;
  for (int k = 0; k < nHRUs; k++) {
    if (hasConnections[k]) {
      bool enabled = pModel->GetHydroUnit(k)->IsEnabled();

      // Check if weights sum deviates significantly from 1.0
      if ((fabs(sumWeights[k] - 1.0) > 0.0001) && enabled) {
        allValid = false;
      }
    }
  }

  // Clean up
  delete[] sumWeights;
  delete[] hasConnections;

  return allValid;
}
