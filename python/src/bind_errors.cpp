#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <format>

#include "cantera/base/ctexceptions.h"
#include "goddard/error.hpp"
#include "goddard/solve_guard.hpp"
#include "goddard/status.hpp"
#include "bind_status.hpp"
#include "goddard_docstrings.h"

namespace nb = nanobind;
using namespace nb::literals;

namespace {

// The exception types, for the translator below. Borrowed: the module holds the references.
PyObject* solve_error_type = nullptr;
PyObject* convergence_error_type = nullptr;

/** Raise `type(what)` with the status attached as attributes. */
void raise_solve_error(PyObject* type, const char* what, const Goddard::SolveStatus& status,
                       bool with_iteration_data) {
    nb::object instance = nb::handle(type)(what);
    instance.attr("status") = nb::cast(status);
    if (with_iteration_data) {
        instance.attr("iterations") = status.iterations;
        instance.attr("tolerance") = status.tolerance;
        instance.attr("residual") = status.residual;
    }
    PyErr_SetObject(type, instance.ptr());
}

} // namespace

void bind_errors(nb::module_& m) {
    nb::enum_<Goddard::FailureReason>(m, "FailureReason", DOC(Goddard, FailureReason))
        .value("NONE", Goddard::FailureReason::NONE, DOC(Goddard, FailureReason, NONE))
        .value("NOT_CONVERGED", Goddard::FailureReason::NOT_CONVERGED,
               DOC(Goddard, FailureReason, NOT_CONVERGED))
        .value("NO_SOLUTION", Goddard::FailureReason::NO_SOLUTION,
               DOC(Goddard, FailureReason, NO_SOLUTION))
        .value("OUT_OF_DATA_RANGE", Goddard::FailureReason::OUT_OF_DATA_RANGE,
               DOC(Goddard, FailureReason, OUT_OF_DATA_RANGE))
        .value("THERMO_FAILURE", Goddard::FailureReason::THERMO_FAILURE,
               DOC(Goddard, FailureReason, THERMO_FAILURE))
        .value("INCOMPLETE", Goddard::FailureReason::INCOMPLETE,
               DOC(Goddard, FailureReason, INCOMPLETE))
        .value("NOT_COMPUTED", Goddard::FailureReason::NOT_COMPUTED,
               DOC(Goddard, FailureReason, NOT_COMPUTED));

    // SolveStatus: read-only, produced by the solvers.
    nb::class_<Goddard::SolveStatus>(m, "SolveStatus", DOC(Goddard, SolveStatus))
        .def_ro("reason", &Goddard::SolveStatus::reason, DOC(Goddard, SolveStatus, reason))
        .def_ro("message", &Goddard::SolveStatus::message, DOC(Goddard, SolveStatus, message))
        .def_ro("stage", &Goddard::SolveStatus::stage, DOC(Goddard, SolveStatus, stage))
        .def_ro("iterations", &Goddard::SolveStatus::iterations,
                DOC(Goddard, SolveStatus, iterations))
        .def_ro("tolerance", &Goddard::SolveStatus::tolerance,
                DOC(Goddard, SolveStatus, tolerance))
        .def_ro("residual", &Goddard::SolveStatus::residual, DOC(Goddard, SolveStatus, residual))
        .def_ro("warnings", &Goddard::SolveStatus::warnings, DOC(Goddard, SolveStatus, warnings))
        .def_prop_ro("ok", &Goddard::SolveStatus::ok, DOC(Goddard, SolveStatus, ok))
        .def("describe", &Goddard::SolveStatus::describe, DOC(Goddard, SolveStatus, describe))
        .def("raise_if_failed", &Goddard::SolveStatus::raise_if_failed,
             DOC(Goddard, SolveStatus, raise_if_failed))
        .def("__repr__", [](const Goddard::SolveStatus& self) {
            if (self.ok()) {
                return self.warnings.empty()
                    ? std::string("<SolveStatus ok>")
                    : std::format("<SolveStatus ok, {} warning(s)>", self.warnings.size());
            }
            return std::format("<SolveStatus {}>", self.describe());
        });

    // SolveError subclasses RuntimeError; ConvergenceError subclasses SolveError. Both carry the
    // failure as `status`; ConvergenceError also has the iteration data as attributes.
    nb::exception<Goddard::SolveError> solve_error(m, "SolveError", PyExc_RuntimeError);
    nb::exception<Goddard::ConvergenceError> convergence_error(m, "ConvergenceError",
                                                               solve_error.ptr());
    solve_error_type = solve_error.ptr();
    convergence_error_type = convergence_error.ptr();

    // error.hpp is not scanned for docstrings, so the class docs are written here.
    solve_error.attr("__doc__") =
        "A numerical solve failed.\n\n"
        "Raised by functions that have no result to carry a status: the state-changing methods\n"
        "of ``Gas`` and scalar-valued solves. Solvers that return a result report a failure as\n"
        "the result's ``status`` instead; ``SolveStatus.raise_if_failed()`` raises it as this\n"
        "exception. Errors of the thermodynamic backend have reason ``THERMO_FAILURE``.\n\n"
        "Attributes:\n"
        "    status: The failure, as a ``SolveStatus``.";
    convergence_error.attr("__doc__") =
        "An iterative solve did not converge: a ``SolveError`` with reason ``NOT_CONVERGED``.\n\n"
        "Attributes:\n"
        "    status: The failure, as a ``SolveStatus``.\n"
        "    iterations: Iterations performed before giving up; -1 if not reported.\n"
        "    tolerance: Convergence tolerance the solver applied; NaN if not reported.\n"
        "    residual: Last residual or step size; NaN if not reported.";

    // Translators registered later are tried first, so this one takes over from the
    // message-only translators of nb::exception.
    nb::register_exception_translator([](const std::exception_ptr& p, void*) {
        try {
            std::rethrow_exception(p);
        } catch (const Goddard::ConvergenceError& e) {
            raise_solve_error(convergence_error_type, e.what(), e.status(), true);
        } catch (const Goddard::SolveError& e) {
            raise_solve_error(solve_error_type, e.what(), e.status(), false);
        } catch (const Goddard::NotImplementedError& e) {
            PyErr_SetString(PyExc_NotImplementedError, e.what());
        } catch (const Cantera::CanteraError&) {
            // Safety net: the thermodynamic backend is an implementation detail, so an error of
            // its own that reaches Python is reported as a Goddard one.
            const Goddard::SolveStatus status = Goddard::status_from_current_exception({});
            raise_solve_error(solve_error_type, status.describe().c_str(), status, false);
        }
    }, nullptr);
}
