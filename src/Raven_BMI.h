/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2023 the Raven Development Team
  ----------------------------------------------------------------*/

// The Basic Model Interface (BMI)  class for Raven.

#ifndef BMI_RAVEN
#define BMI_RAVEN

#include "BMI.h"
#include "Model.h"

class CRavenBMI : public bmixx::Bmi
{
  private:
    CModel     *pModel;
    optStruct   Options;
    time_struct tt;

    // Internal functions for reading the YAML config file.
    void ReadConfigFile(std::string config_file);
    std::vector<char *> SplitLine(std::string line);
    void ProcessConfigFileArgument(std::string config_key, std::string config_value);

  public:
    CRavenBMI();
    ~CRavenBMI();

    // Model control functions.
    void Initialize(std::string config_file);
    void Update();
    void UpdateUntil(double time);
    void Finalize();

    // Model information functions.
    std::string GetComponentName();
    int GetInputItemCount();
    int GetOutputItemCount();
    std::vector<std::string> GetInputVarNames();
    std::vector<std::string> GetOutputVarNames();

    // Variable information functions
    int GetVarGrid(std::string name);
    std::string GetVarType(std::string name);
    std::string GetVarUnits(std::string name);
    int GetVarItemsize(std::string name);
    int GetVarNbytes(std::string name);
    std::string GetVarLocation(std::string name);

    double GetCurrentTime();
    double GetStartTime();
    double GetEndTime();
    std::string GetTimeUnits();
    double GetTimeStep();

    // Variable getters
    void GetValue(std::string name, void *dest);
    void *GetValuePtr(std::string name);
    void GetValueAtIndices(std::string name, void *dest, int *inds, int count);

    // Variable setters
    void SetValue(std::string name, void *src);
    void SetValueAtIndices(std::string name, int *inds, int count, void *src);

    // Grid information functions
    int GetGridRank(const int grid);
    int GetGridSize(const int grid);
    std::string GetGridType(const int grid);

    void GetGridShape(const int grid, int *shape);
    void GetGridSpacing(const int grid, double *spacing);
    void GetGridOrigin(const int grid, double *origin);

    void GetGridX(const int grid, double *x);
    void GetGridY(const int grid, double *y);
    void GetGridZ(const int grid, double *z);

    int  GetGridNodeCount(const int grid);
    int  GetGridEdgeCount(const int grid);
    int  GetGridFaceCount(const int grid);

    void GetGridEdgeNodes(const int grid, int *edge_nodes);
    void GetGridFaceEdges(const int grid, int *face_edges);
    void GetGridFaceNodes(const int grid, int *face_nodes);
    void GetGridNodesPerFace(const int grid, int *nodes_per_face);
};


////////////////////////////////////////////////////////////////////
/// NextGen expects two externalized functions, one to create a new
/// instance of the model and one to destroy/free an instance.
/// More: https://github.com/NOAA-OWP/ngen/blob/master/doc/BMI_MODELS.md#bmi-c-model-as-shared-library-1

extern "C"
{
  //////////////////////////////////////////////////////////////////
  /// \brief Create a new instance of the model as expected by NextGen.
  /// \return A pointer to the newly allocated instance.
  //
	CRavenBMI *bmi_model_create()
	{
		return new CRavenBMI();
	}

  //////////////////////////////////////////////////////////////////
  /// \brief Destroy/free an instance created with @see bmi_model_create
  /// \param ptr A pointer to the instance to be destroyed.
  //
	void bmi_model_destroy(CRavenBMI *ptr)
	{
		delete ptr;
	}
}


#endif
