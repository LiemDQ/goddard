#include <nanobind/nanobind.h>
#include "goddard/error.hpp"

namespace nb = nanobind;

void bind_errors(nb::module_& m) {
    // ConvergenceError subclasses RuntimeError and carries the solver's iteration count,
    // tolerance and last residual as attributes (-1 when the solver does not report them).
    nb::exception<Goddard::ConvergenceError> convergence_error(m, "ConvergenceError",
                                                               PyExc_RuntimeError);

    // Translators registered later are tried first, so this one takes over from the
    // message-only translator of nb::exception.
    nb::register_exception_translator([](const std::exception_ptr& p, void* payload) {
        try {
            std::rethrow_exception(p);
        } catch (const Goddard::ConvergenceError& e) {
            nb::handle type(static_cast<PyObject*>(payload));
            nb::object instance = type(e.what());
            instance.attr("iterations") = e.iterations();
            instance.attr("tolerance") = e.tolerance();
            instance.attr("residual") = e.residual();
            PyErr_SetObject(type.ptr(), instance.ptr());
        } catch (const Goddard::NotImplementedError& e) {
            PyErr_SetString(PyExc_NotImplementedError, e.what());
        }
    }, convergence_error.ptr());
}
