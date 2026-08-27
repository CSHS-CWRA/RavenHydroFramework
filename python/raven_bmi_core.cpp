#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
//#include "RavenInclude.h"
//#include "BMI.h"
//#include "Raven_BMI.h"
#include <RavenInclude.h>
#include <BMI.h>
#include <Raven_BMI.h>

#ifdef _RVNETCDF_
const bool    __HAS_NETCDF__ = true;
#endif
#ifndef _RVNETCDF_
const bool    __HAS_NETCDF__ = false;
#endif

namespace py = pybind11;

template <typename T>
using py_array_cont = py::array_t<T, py::array::c_style | py::array::forcecast>;

//PYBIND11_MODULE(libraven, m) {
PYBIND11_MODULE(ravenbmi_pycore, m) {
    m.doc() =
      R"pbdoc(A thin Python wrapper around the BMI C++ interface to the hydrologic modelling framework Raven.)pbdoc";

    m.attr("__raven_version__") = __RAVEN_VERSION__;
    m.attr("__netcdf__") = __HAS_NETCDF__;

    py::class_<bmixx::Bmi, py::smart_holder> abc(m, "_BMI");

    pybind11::class_<CRavenBMI, bmixx::Bmi, py::smart_holder> cl(m, "CRavenBMI", "");
    cl.def(pybind11::init( [](CRavenBMI const &o){ return new CRavenBMI(o); } ) );
    cl.def(py::init<>());
    cl.def("Initialize", (void (CRavenBMI::*)(std::string)) &CRavenBMI::Initialize, "C++: CRavenBMI::Initialize(std::string) --> void", pybind11::arg("config_file"));
    cl.def("Update", (void (CRavenBMI::*)()) &CRavenBMI::Update, "C++: CRavenBMI::Update() --> void");
    cl.def("UpdateUntil", (void (CRavenBMI::*)(double)) &CRavenBMI::UpdateUntil, "C++: CRavenBMI::UpdateUntil(double) --> void", pybind11::arg("time"));
    cl.def("Finalize", (void (CRavenBMI::*)()) &CRavenBMI::Finalize, "C++: CRavenBMI::Finalize() --> void");
    cl.def("GetComponentName", (std::string (CRavenBMI::*)()) &CRavenBMI::GetComponentName, "C++: CRavenBMI::GetComponentName() --> std::string");
    cl.def("GetInputItemCount", (int (CRavenBMI::*)()) &CRavenBMI::GetInputItemCount, "C++: CRavenBMI::GetInputItemCount() --> int");
    cl.def("GetOutputItemCount", (int (CRavenBMI::*)()) &CRavenBMI::GetOutputItemCount, "C++: CRavenBMI::GetOutputItemCount() --> int");
    cl.def("GetInputVarNames", (class std::vector<std::string > (CRavenBMI::*)()) &CRavenBMI::GetInputVarNames, "C++: CRavenBMI::GetInputVarNames() --> class std::vector<std::string >");
    cl.def("GetOutputVarNames", (class std::vector<std::string > (CRavenBMI::*)()) &CRavenBMI::GetOutputVarNames, "C++: CRavenBMI::GetOutputVarNames() --> class std::vector<std::string >");
    cl.def("GetVarGrid", (int (CRavenBMI::*)(std::string)) &CRavenBMI::GetVarGrid, "C++: CRavenBMI::GetVarGrid(std::string) --> int", pybind11::arg("name"));
    cl.def("GetVarType", (std::string (CRavenBMI::*)(std::string)) &CRavenBMI::GetVarType, "C++: CRavenBMI::GetVarType(std::string) --> std::string", pybind11::arg("name"));
    cl.def("GetVarUnits", (std::string (CRavenBMI::*)(std::string)) &CRavenBMI::GetVarUnits, "C++: CRavenBMI::GetVarUnits(std::string) --> std::string", pybind11::arg("name"));
    cl.def("GetVarItemsize", (int (CRavenBMI::*)(std::string)) &CRavenBMI::GetVarItemsize, "C++: CRavenBMI::GetVarItemsize(std::string) --> int", pybind11::arg("name"));
    cl.def("GetVarNbytes", (int (CRavenBMI::*)(std::string)) &CRavenBMI::GetVarNbytes, "C++: CRavenBMI::GetVarNbytes(std::string) --> int", pybind11::arg("name"));
    cl.def("GetVarLocation", (std::string (CRavenBMI::*)(std::string)) &CRavenBMI::GetVarLocation, "C++: CRavenBMI::GetVarLocation(std::string) --> std::string", pybind11::arg("name"));
    cl.def("GetCurrentTime", (double (CRavenBMI::*)()) &CRavenBMI::GetCurrentTime, "C++: CRavenBMI::GetCurrentTime() --> double");
    cl.def("GetStartTime", (double (CRavenBMI::*)()) &CRavenBMI::GetStartTime, "C++: CRavenBMI::GetStartTime() --> double");
    cl.def("GetEndTime", (double (CRavenBMI::*)()) &CRavenBMI::GetEndTime, "C++: CRavenBMI::GetEndTime() --> double");
    cl.def("GetTimeUnits", (std::string (CRavenBMI::*)()) &CRavenBMI::GetTimeUnits, "C++: CRavenBMI::GetTimeUnits() --> std::string");
    cl.def("GetTimeStep", (double (CRavenBMI::*)()) &CRavenBMI::GetTimeStep, "C++: CRavenBMI::GetTimeStep() --> double");
    //cl.def("GetValue", (void (CRavenBMI::*)(std::string, void *)) &CRavenBMI::GetValue, "C++: CRavenBMI::GetValue(std::string, void *) --> void", pybind11::arg("name"), pybind11::arg("dest"));
    cl.def("GetValue",
      [](CRavenBMI &r, std::string name, py_array_cont<double> dest) -> py_array_cont<double> {
        int grid_size = r.GetGridSize(r.GetVarGrid(name));
        if (dest.size() != grid_size) {
          throw std::runtime_error("Size of output array (" + std::to_string(dest.size()) + ") does not match"
            +" size of variable '" + name +"' (" + std::to_string(grid_size) + ")."
          );
        }
        r.GetValue(name, dest.mutable_data());
        return dest;
      },
      "C++: CRavenBMI::GetValue(std::string, py::array_t<double>) --> py::array_t<double>",
      pybind11::arg("name"),
      pybind11::arg("dest")
    );
    //cl.def("GetValuePtr", (void * (CRavenBMI::*)(std::string)) &CRavenBMI::GetValuePtr, "C++: CRavenBMI::GetValuePtr(std::string) --> void *", pybind11::return_value_policy::automatic, pybind11::arg("name"));
    //
    cl.def("GetValuePtr",
      [](CRavenBMI &r, std::string name) -> py::array_t<double> {
        int num_elements = r.GetGridSize(r.GetVarGrid(name));
        double *data = static_cast<double *>(r.GetValuePtr(name));
        // CRavenBMI::GetValuePtr does not return a pointer to internal storage,
        // instead it allocates new memory on the heap and returns a pointer to that.
        // Python needs to take ownership and free the memory.
        return py::array_t<double>(
          { static_cast<size_t>(num_elements) },  // shape
          { sizeof(double) },                      // strides
          data,                                    // data pointer
          py::capsule(data, [](void *p) { delete[] static_cast<double *>(p); }));
      },
      "C++: CRavenBMI::GetValuePtr(std::string) --> py::array_t<double>",
      pybind11::arg("name")
    );
    //cl.def("GetValueAtIndices", (void (CRavenBMI::*)(std::string, void *, int *, int)) &CRavenBMI::GetValueAtIndices, "C++: CRavenBMI::GetValueAtIndices(std::string, void *, int *, int) --> void", pybind11::arg("name"), pybind11::arg("dest"), pybind11::arg("inds"), pybind11::arg("count"));
    cl.def("GetValueAtIndices",
      [](CRavenBMI &r, std::string name, py_array_cont<double> dest, py_array_cont<int> inds) -> py_array_cont<double> {
        int grid_size = r.GetGridSize(r.GetVarGrid(name));
        size_t num_elements = inds.size();
        if (num_elements != dest.size()) {
          throw std::runtime_error("Size of output array (" + std::to_string(dest.size()) + ") does not match"
            +" size of index array (" + std::to_string(num_elements) + ")."
          );
        }
        const int *const ind_data = inds.data();
        // bound check supplied index values
        for (size_t i = 0; i < num_elements; i++) {
          int p = ind_data[i];
          if (p < 0 || p >= grid_size) {
            throw std::out_of_range("Index " + std::to_string(i) + " has value " + std::to_string(p)
              + ", which is outside the valid range [0, " + std::to_string(grid_size - 1) + "] for variable '" + name + "'."
            );
          }
        }
        r.GetValueAtIndices(name, dest.mutable_data(), inds.mutable_data(), static_cast<int>(num_elements));
        return dest;
      },
      "C++: CRavenBMI::GetValueAtIndices(std::string, py::array_t<double>, py::array_t<int>) --> py::array_t<double>",
      pybind11::arg("name"),
      pybind11::arg("dest"),
      pybind11::arg("inds")
    );
    //cl.def("SetValue", (void (CRavenBMI::*)(std::string, void *)) &CRavenBMI::SetValue, "C++: CRavenBMI::SetValue(std::string, void *) --> void", pybind11::arg("name"), pybind11::arg("src"));
    cl.def("SetValue",
      [](CRavenBMI &r, std::string name, py_array_cont<double> src) {
        int grid_size = r.GetGridSize(r.GetVarGrid(name));
        if (src.size() != grid_size) {
          throw std::runtime_error("Source array has wrong size with size " + std::to_string(src.size()) + "."
            + " Variable '" + name + "' needs to set exactly " + std::to_string(grid_size) + " values."
          );
        }
        r.SetValue(name, src.mutable_data());
      },
      "C++: CRavenBMI::SetValue(std::string, py::array_t<double>) --> void",
      pybind11::arg("name"),
      pybind11::arg("src")
    );
    //cl.def("SetValueAtIndices", (void (CRavenBMI::*)(std::string, int *, int, void *)) &CRavenBMI::SetValueAtIndices, "C++: CRavenBMI::SetValueAtIndices(std::string, int *, int, void *) --> void", pybind11::arg("name"), pybind11::arg("inds"), pybind11::arg("count"), pybind11::arg("src"));
    cl.def("SetValueAtIndices",
      [](CRavenBMI &r, std::string name, py_array_cont<int> inds, py_array_cont<double> src) {
        int grid_size = r.GetGridSize(r.GetVarGrid(name));
        size_t num_elements = inds.size();
        if (num_elements != src.size()) {
          throw std::runtime_error("Size of source array (" + std::to_string(src.size()) + ") does not match"
            +" size of index array (" + std::to_string(num_elements) + ")."
          );
        }
        const int *const ind_data = inds.data();
        // bound check supplied index values
        for (size_t i = 0; i < num_elements; i++) {
          int p = ind_data[i];
          if (p < 0 || p >= grid_size) {
            throw std::out_of_range("Index " + std::to_string(i) + " has value " + std::to_string(p)
              + ", which is outside the valid range [0, " + std::to_string(grid_size - 1) + "] for variable '" + name + "'."
            );
          }
        }
        r.SetValueAtIndices(name, inds.mutable_data(), static_cast<int>(num_elements), src.mutable_data());
      },
      "C++: CRavenBMI::SetValueAtIndices(std::string, py::array_t<int>, py::array_t<double>) --> void",
      pybind11::arg("name"),
      pybind11::arg("inds"),
      pybind11::arg("src")
    );
    cl.def("GetGridRank", (int (CRavenBMI::*)(const int)) &CRavenBMI::GetGridRank, "C++: CRavenBMI::GetGridRank(const int) --> int", pybind11::arg("grid"));
    cl.def("GetGridSize", (int (CRavenBMI::*)(const int)) &CRavenBMI::GetGridSize, "C++: CRavenBMI::GetGridSize(const int) --> int", pybind11::arg("grid"));
    cl.def("GetGridType", (std::string (CRavenBMI::*)(const int)) &CRavenBMI::GetGridType, "C++: CRavenBMI::GetGridType(const int) --> std::string", pybind11::arg("grid"));
    cl.def("GetGridShape", (void (CRavenBMI::*)(const int, int *)) &CRavenBMI::GetGridShape, "C++: CRavenBMI::GetGridShape(const int, int *) --> void", pybind11::arg("grid"), pybind11::arg("shape"));
    cl.def("GetGridSpacing", (void (CRavenBMI::*)(const int, double *)) &CRavenBMI::GetGridSpacing, "C++: CRavenBMI::GetGridSpacing(const int, double *) --> void", pybind11::arg("grid"), pybind11::arg("spacing"));
    cl.def("GetGridOrigin", (void (CRavenBMI::*)(const int, double *)) &CRavenBMI::GetGridOrigin, "C++: CRavenBMI::GetGridOrigin(const int, double *) --> void", pybind11::arg("grid"), pybind11::arg("origin"));
    cl.def("GetGridX", (void (CRavenBMI::*)(const int, double *)) &CRavenBMI::GetGridX, "C++: CRavenBMI::GetGridX(const int, double *) --> void", pybind11::arg("grid"), pybind11::arg("x"));
    cl.def("GetGridY", (void (CRavenBMI::*)(const int, double *)) &CRavenBMI::GetGridY, "C++: CRavenBMI::GetGridY(const int, double *) --> void", pybind11::arg("grid"), pybind11::arg("y"));
    cl.def("GetGridZ", (void (CRavenBMI::*)(const int, double *)) &CRavenBMI::GetGridZ, "C++: CRavenBMI::GetGridZ(const int, double *) --> void", pybind11::arg("grid"), pybind11::arg("z"));
    cl.def("GetGridNodeCount", (int (CRavenBMI::*)(const int)) &CRavenBMI::GetGridNodeCount, "C++: CRavenBMI::GetGridNodeCount(const int) --> int", pybind11::arg("grid"));
    cl.def("GetGridEdgeCount", (int (CRavenBMI::*)(const int)) &CRavenBMI::GetGridEdgeCount, "C++: CRavenBMI::GetGridEdgeCount(const int) --> int", pybind11::arg("grid"));
    cl.def("GetGridFaceCount", (int (CRavenBMI::*)(const int)) &CRavenBMI::GetGridFaceCount, "C++: CRavenBMI::GetGridFaceCount(const int) --> int", pybind11::arg("grid"));
    cl.def("GetGridEdgeNodes", (void (CRavenBMI::*)(const int, int *)) &CRavenBMI::GetGridEdgeNodes, "C++: CRavenBMI::GetGridEdgeNodes(const int, int *) --> void", pybind11::arg("grid"), pybind11::arg("edge_nodes"));
    cl.def("GetGridFaceEdges", (void (CRavenBMI::*)(const int, int *)) &CRavenBMI::GetGridFaceEdges, "C++: CRavenBMI::GetGridFaceEdges(const int, int *) --> void", pybind11::arg("grid"), pybind11::arg("face_edges"));
    cl.def("GetGridFaceNodes", (void (CRavenBMI::*)(const int, int *)) &CRavenBMI::GetGridFaceNodes, "C++: CRavenBMI::GetGridFaceNodes(const int, int *) --> void", pybind11::arg("grid"), pybind11::arg("face_nodes"));
    cl.def("GetGridNodesPerFace", (void (CRavenBMI::*)(const int, int *)) &CRavenBMI::GetGridNodesPerFace, "C++: CRavenBMI::GetGridNodesPerFace(const int, int *) --> void", pybind11::arg("grid"), pybind11::arg("nodes_per_face"));
}
