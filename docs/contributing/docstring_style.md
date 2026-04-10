# Docstring Style Guide

Goddard uses a pipeline that automatically propagates C++ documentation comments
into Python docstrings. This guide explains the conventions so that every public
symbol appears correctly in the rendered API reference.

## How the pipeline works

```
C++ header (Doxygen comment)
  -> pybind11_mkdoc (libclang)
  -> goddard_docstrings.h (DOC() macro)
  -> nanobind binding (.def with DOC())
  -> compiled _core.so
  -> nanobind stubgen -> _core.pyi
  -> mkdocstrings renders API page
```

Because of this pipeline, **C++ headers are the single source of truth** for
documentation. You should never duplicate a docstring between a C++ header and a
binding file.

---

## C++ headers (Doxygen format)

Every public symbol (class, method, free function, struct) that should appear in
the Python API docs **must** have a Doxygen comment in the C++ header.

### Comment style

Use `/** ... */` block comments immediately above the declaration.  libclang
requires this exact placement to associate the comment with the symbol.

```cpp
/** Set state from temperature [K] and pressure [Pa]. */
void set_state_TP(double T, double P);
```

For multi-line comments:

```cpp
/**
 * Stagnation (total) enthalpy h0 = h + v^2/2.
 *
 * @param velocity flow velocity [m/s]
 * @return stagnation enthalpy [J/kg]
 */
double stagnation_enthalpy(double velocity) const;
```

### Structure

- **First sentence** is the brief summary. pybind11_mkdoc uses it automatically
  -- do not add `@brief`.
- Separate the brief from the body with a blank line.
- Use `@param name description` for parameters.
- Use `@return description` for return values.
- Use `@note`, `@warning`, `@throws` for callouts.

### Units

Goddard is unit-sensitive. **Always** document units inline:

```cpp
/** Set state from specific enthalpy [J/kg] and pressure [Pa]. */
void set_state_HP(double H, double P);
```

### What NOT to do

```cpp
// Bad: wrong comment style -- libclang ignores // comments
// Sets the temperature
void set_state_TP(double T, double P);

// Bad: comment not immediately above the declaration
/** Set the temperature. */

void set_state_TP(double T, double P);

// Bad: too vague, no units
/** Setter. */
void set_state_TP(double T, double P);
```

---

## Binding files (nanobind + DOC() macro)

When a binding file includes `goddard_docstrings.h`, it can reference
documentation with the `DOC()` macro:

```cpp
#include "goddard_docstrings.h"

// Class-level doc:
nb::class_<Goddard::Gas>(m, "Gas", DOC(Goddard, Gas))

// Method-level doc:
.def("set_state_TP", &Goddard::Gas::set_state_TP,
     "T"_a, "P"_a, DOC(Goddard, Gas, set_state_TP))
```

### Overloads

pybind11_mkdoc assigns a suffix to each overload in declaration order:

```cpp
// First overload  -> DOC(Goddard, Gas, isenthalpic_velocity)
double isenthalpic_velocity(double H_stagnation) const;
// Second overload -> DOC(Goddard, Gas, isenthalpic_velocity, 2)
double isenthalpic_velocity() const;
```

### Keep argument annotations

Even when using generated docstrings, nanobind still needs `"param"_a` annotations
for keyword-argument support. These are orthogonal to `DOC()`:

```cpp
.def("set_state_TP", &Goddard::Gas::set_state_TP,
     "T"_a, "P"_a, DOC(Goddard, Gas, set_state_TP))
//   ^^^^^^^^^^^^  still required
```

### Fallback behavior

When `goddard_BUILD_DOCSTRINGS=OFF` (the default outside the `docs` pixi
environment), `DOC(...)` expands to `""`. The binding compiles normally; Python
methods just have empty docstrings. This keeps the build working without
libclang.

---

## Pure-Python functions (Google style)

Functions in `python/goddard/convenience.py` and other pure-Python modules use
**Google-style** docstrings. mkdocstrings auto-detects this format:

```python
def gas_from_yaml(yaml_file, phase_name="", species=None,
                  chemistry=GasChemistry.FROZEN):
    """Create a Gas from a YAML thermodynamic data file.

    Args:
        yaml_file: Path to a Cantera YAML thermodynamic data file.
        phase_name: Name of the phase in the YAML file. Defaults to the
            first phase found.
        species: Optional list of species names to include. If None, all
            species in the phase are used.
        chemistry: GasChemistry mode for derived property calculations.

    Returns:
        Gas object initialized from the YAML data.

    Example:
        >>> gas = gas_from_yaml("h2o2.yaml", phase_name="gas")
        >>> gas.set_state_TP(3000.0, 101325.0)
    """
```

---

## What to document

- **Public classes** and their public methods
- **Free functions** exposed in the Python API
- **Struct fields** (use a Doxygen comment on the member)
- **Enum values** when the meaning isn't obvious from the name

## What to skip

- Private members (`m_*`)
- Internal helpers not exposed through bindings
- Trivial getters where the name is self-documenting and there are no units
  (e.g. `name()`, `num_species()`)
