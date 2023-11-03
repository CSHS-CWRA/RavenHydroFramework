#include <pybind11/pybind11.h>
#include "RavenInclude.h"

#ifdef _RVNETCDF_
const bool    __HAS_NETCDF__ = true;
#endif
#ifndef _RVNETCDF_
const bool    __HAS_NETCDF__ = false;
#endif

namespace py = pybind11;

PYBIND11_MODULE(libraven, m) {
    m.doc() =
      R"pbdoc(A Python wrapper to the hydrologic modelling framework Raven.)pbdoc";

    m.attr("__version__") = __RAVEN_VERSION__;
    m.attr("__netcdf__") = __HAS_NETCDF__;
}
