#include <nanobind/nanobind.h>
#include <format>
#include <functional>
#include <string>
#include <vector>
#include <nanobind/stl/function.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/eigen/dense.h>
#include "goddard/thermoarray.hpp"
#include "goddard_docstrings.h"

namespace nb = nanobind;

namespace {

using Goddard::ThermoArray;

/** The array's shape as a Python tuple. */
nb::tuple shape_tuple(const ThermoArray& self) {
    nb::list extents;
    for (long extent : self.shape()) {
        extents.append(extent);
    }
    return nb::tuple(extents);
}

/** One value per entry in storage order, as a numpy array of the array's shape. */
nb::object shaped(const ThermoArray& self, const Eigen::ArrayXd& values) {
    return nb::cast(values).attr("reshape")(shape_tuple(self));
}

nb::object shaped(const ThermoArray& self, const Eigen::Array<bool, Eigen::Dynamic, 1>& values) {
    return nb::cast(values).attr("reshape")(shape_tuple(self));
}

/** One row per entry, as a numpy array of shape `shape + (columns,)`. */
nb::object shaped_rows(const ThermoArray& self, const Eigen::ArrayXXd& rows) {
    nb::list extents;
    for (long extent : self.shape()) {
        extents.append(extent);
    }
    extents.append(rows.cols());
    return nb::cast(rows).attr("reshape")(nb::tuple(extents));
}

/**
 * Flat location of a Python index: an int for a 1-D array, or a tuple with one int per
 * dimension. Negative indices count from the end, as in numpy.
 */
int location(const ThermoArray& self, const nb::handle& index) {
    std::vector<long> indices;
    if (nb::isinstance<nb::tuple>(index)) {
        for (nb::handle item : nb::borrow<nb::tuple>(index)) {
            indices.push_back(nb::cast<long>(item));
        }
    } else {
        indices.push_back(nb::cast<long>(index));
    }
    if (indices.size() != self.shape().size()) {
        throw nb::index_error(std::format("ThermoArray of {} dimensions indexed with {} indices; "
            "give one integer per dimension.", self.shape().size(), indices.size()).c_str());
    }
    for (size_t dim = 0; dim < indices.size(); dim++) {
        const long extent = self.shape()[dim];
        const long given = indices[dim];
        if (indices[dim] < 0) {
            indices[dim] += extent;
        }
        if (indices[dim] < 0 || indices[dim] >= extent) {
            throw nb::index_error(std::format("ThermoArray index {} out of range for dimension {} "
                "of extent {}.", given, dim, extent).c_str());
        }
    }
    return self.flat_index(indices);
}

} // namespace

void bind_thermoarray(nb::module_& m) {
    using Goddard::CompositionModel;

    nb::class_<ThermoArray>(m, "ThermoArray", DOC(Goddard, ThermoArray))
        .def_prop_ro("shape", &shape_tuple, "Extent of each dimension, as a tuple.")
        .def_prop_ro("size", &ThermoArray::size, DOC(Goddard, ThermoArray, size))
        .def_prop_ro("ndim", &ThermoArray::ndim, DOC(Goddard, ThermoArray, ndim))
        .def_prop_ro("species_names", &ThermoArray::species_names,
             DOC(Goddard, ThermoArray, species_names))
        .def_prop_ro("num_condensed", &ThermoArray::num_condensed,
             DOC(Goddard, ThermoArray, num_condensed))
        .def_prop_ro("condensed_species_names", &ThermoArray::condensed_species_names,
             DOC(Goddard, ThermoArray, condensed_species_names))
        .def("__getitem__", [](const ThermoArray& self, nb::handle index) {
                return self.gas_at(location(self, index));
            }, nb::arg("index"),
            "Independent Gas in the state of the entry at `index` (an int for a 1-D array, a "
            "tuple of ints otherwise). Raises SolveError for a failed entry.")
        .def("get_state", [](const ThermoArray& self, nb::handle index) {
                return self.get_state(location(self, index));
            }, nb::arg("index"), DOC(Goddard, ThermoArray, get_state))
        .def("get_condensed_moles", [](const ThermoArray& self, nb::handle index) {
                return self.get_condensed_moles(location(self, index));
            }, nb::arg("index"), DOC(Goddard, ThermoArray, get_condensed_moles))
        .def("snapshot", [](const ThermoArray& self, nb::handle index, CompositionModel model) {
                return self.snapshot(location(self, index), model);
            }, nb::arg("index"), nb::arg("model") = CompositionModel::FROZEN,
            DOC(Goddard, ThermoArray, snapshot))
        .def("evaluate", [](const ThermoArray& self, const std::function<double(const Goddard::Gas&)>& property) {
                return shaped(self, self.evaluate(property));
            }, nb::arg("property"),
            "Evaluate `property(gas)` for the Gas of every entry and return the values as an "
            "array of the array's shape; failed entries are NaN. The Gas passed to `property` "
            "is reused between entries and must not be kept.")
        .def_prop_ro("temperature", [](const ThermoArray& self) {
            return shaped(self, self.temperature()); }, DOC(Goddard, ThermoArray, temperature))
        .def_prop_ro("pressure", [](const ThermoArray& self) {
            return shaped(self, self.pressure()); }, DOC(Goddard, ThermoArray, pressure))
        .def_prop_ro("density", [](const ThermoArray& self) {
            return shaped(self, self.density()); }, DOC(Goddard, ThermoArray, density))
        .def_prop_ro("enthalpy_mass", [](const ThermoArray& self) {
            return shaped(self, self.enthalpy_mass()); }, DOC(Goddard, ThermoArray, enthalpy_mass))
        .def_prop_ro("entropy_mass", [](const ThermoArray& self) {
            return shaped(self, self.entropy_mass()); }, DOC(Goddard, ThermoArray, entropy_mass))
        .def_prop_ro("internal_energy_mass", [](const ThermoArray& self) {
            return shaped(self, self.internal_energy_mass()); },
            DOC(Goddard, ThermoArray, internal_energy_mass))
        .def_prop_ro("cp_mass", [](const ThermoArray& self) {
            return shaped(self, self.cp_mass()); }, DOC(Goddard, ThermoArray, cp_mass))
        .def_prop_ro("cv_mass", [](const ThermoArray& self) {
            return shaped(self, self.cv_mass()); }, DOC(Goddard, ThermoArray, cv_mass))
        .def_prop_ro("mean_molecular_weight", [](const ThermoArray& self) {
            return shaped(self, self.molecular_weight()); },
            DOC(Goddard, ThermoArray, molecular_weight))
        .def_prop_ro("mixture_molecular_weight", [](const ThermoArray& self) {
            return shaped(self, self.mixture_molecular_weight()); },
            DOC(Goddard, ThermoArray, mixture_molecular_weight))
        .def_prop_ro("gas_mass_fraction", [](const ThermoArray& self) {
            return shaped(self, self.gas_mass_fraction()); },
            DOC(Goddard, ThermoArray, gas_mass_fraction))
        .def_prop_ro("gamma_s", [](const ThermoArray& self) {
            return shaped(self, self.gamma_s(CompositionModel::FROZEN)); },
            "Isentropic exponent at frozen composition [-], as Gas.gamma_s.")
        .def_prop_ro("equilibrium_gamma_s", [](const ThermoArray& self) {
            return shaped(self, self.gamma_s(CompositionModel::EQUILIBRIUM)); },
            "Isentropic exponent at shifting equilibrium [-], as Gas.equilibrium_gamma_s.")
        .def_prop_ro("speed_of_sound", [](const ThermoArray& self) {
            return shaped(self, self.speed_of_sound(CompositionModel::FROZEN)); },
            "Speed of sound at frozen composition [m/s], as Gas.speed_of_sound.")
        .def_prop_ro("equilibrium_speed_of_sound", [](const ThermoArray& self) {
            return shaped(self, self.speed_of_sound(CompositionModel::EQUILIBRIUM)); },
            "Speed of sound at shifting equilibrium [m/s], as Gas.equilibrium_speed_of_sound.")
        .def_prop_ro("mole_fractions", [](const ThermoArray& self) {
            return shaped_rows(self, self.mole_fractions()); },
            "Mole fractions of the gas-phase species [-], of shape `shape + (len(species_names),)`.")
        .def_prop_ro("mass_fractions", [](const ThermoArray& self) {
            return shaped_rows(self, self.mass_fractions()); },
            "Mass fractions of the gas-phase species [-], of shape `shape + (len(species_names),)`.")
        .def_prop_ro("mixture_mass_fractions", [](const ThermoArray& self) {
            return shaped_rows(self, self.mixture_mass_fractions()); },
            "Mass fractions of the whole mixture [-], as Gas.mixture_mass_fractions, of shape "
            "`shape + (len(species_names) + num_condensed,)`.")
        .def_prop_ro("status", &ThermoArray::status, DOC(Goddard, ThermoArray, status))
        .def("status_at", [](const ThermoArray& self, nb::handle index) {
                return self.status_at(location(self, index));
            }, nb::arg("index"), DOC(Goddard, ThermoArray, status_at))
        .def_prop_ro("num_failed", &ThermoArray::num_failed, DOC(Goddard, ThermoArray, num_failed))
        .def_prop_ro("ok", [](const ThermoArray& self) { return shaped(self, self.ok()); },
             "Whether each entry is ok, as a boolean array of the array's shape.")
        .def("__repr__", [](const ThermoArray& self) {
            std::string shape;
            for (long extent : self.shape()) {
                shape += (shape.empty() ? "" : ", ") + std::to_string(extent);
            }
            return std::format("<ThermoArray shape=({}) {} of {} entries failed>", shape,
                               self.num_failed(), self.size());
        });
}
