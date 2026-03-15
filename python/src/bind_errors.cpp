#include <nanobind/nanobind.h>
#include "goddard/error.hpp"

namespace nb = nanobind;

void bind_errors(nb::module_& m) {
    nb::exception<Goddard::ConvergenceError>(m, "ConvergenceError");

    nb::register_exception_translator([](const std::exception_ptr& p, void*) {
        try {
            std::rethrow_exception(p);
        } catch (const Goddard::NotImplementedError& e) {
            PyErr_SetString(PyExc_NotImplementedError, e.what());
        }
    });
}
