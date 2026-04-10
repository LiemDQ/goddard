# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Goddard is a C++/Python rocket engine simulation toolkit that uses Cantera for chemical kinetics and thermodynamics. It aims to replicate and extend NASA CEA (Chemical Equilibrium with Applications) functionality with a modern API. The project is structured as a hybrid C++/Python codebase with CMake build system.

## Core Architecture

### Features
Goddard performs 5 main types of computations:
1. General thermodynamic and compressible flow properties: speed of sound, stagnation pressure, stagnation enthalpy, heat capacity ratio, etc, integrated with Cantera's thermodynamic solvers.
2. Combustion/equilibrium calculations. 
3. 1D nozzle flow: frozen, equilibrium and kinetic chemistries.
4. Method of characteristics for 2D or axisymmetric supersonic flow fields for design and analysis.
5. Shock properties (incident, reflected, and oblique)


### C++ Core Components
- **Gas** (`gas.hpp/cpp`): Core primitive for querying thermodynamic properties. Building block for the main solvers. 
- **Equilibrium** (`equilibrium.hpp/cpp`): Chemical equilibrium calculations and thermodynamic derivatives
- **Combustor** (`combustor.hpp/cpp`): Isobaric combustion reaction handling with support for infinite area, finite mass flux, and finite contraction ratio modes
- **Nozzle** (`nozzle.hpp/cpp`): Nozzle flow calculations with inheritance hierarchy:
  - `NozzleBase`: Abstract base class
  - `EquilibriumNozzle`: Chemical equilibrium nozzle flow
  - `FrozenNozzle`: Frozen composition nozzle flow
- **ThermoArray** (`thermoarray.hpp/cpp`): Batch thermodynamic property calculations
- **MoC** (`moc.hpp/cpp`, `characteristics.hpp/cpp`, `prandtlmeyer.hpp/cpp`): 2D supersonic nozzle flow via Method of Characteristics
  - `MocNozzle`: Solver class supporting design (minimum-length) and analysis modes
  - Planar and axisymmetric flow with perfect gas, frozen, or equilibrium chemistry
  - `CharacteristicNet`/`CharacteristicPoint`: Flow field mesh and point data
  - `NozzleProfile`: Wall contour representation with CSV I/O
  - `PrandtlMeyerTable`: Precomputed isentropic expansion data for non-ideal gas
  - `compute_thrust_coefficient()`: Exit plane integration for thrust performance
- **KineticNozzle** (`kinetic_nozzle.hpp/cpp`): 1D supersonic nozzle with finite-rate chemistry via Cantera's `IdealGasMoleReactor`. Standalone class — does not inherit `NozzleBase` because its spatially-resolved output is incompatible with the discrete area-ratio interface. Instead uses composition: owns a `NozzleBase`-derived object internally for throat conditions only. Output is `KineticNozzleResults` containing a `ThroatCondition` and a vector of `KineticNozzleStation` (x, velocity, Mach, area_ratio, thermo state, per-species Damköhler numbers). The `NozzleChemistryType` constructor parameter selects the throat model (EQUILIBRIUM or FROZEN); KINETIC is not valid as a throat model.
- **Shocks** (`shocks.hpp/cpp`): Normal shock relations (Rankine-Hugoniot) with perfect-gas and Cantera-state variants. `ShockResult` and `ObliqueShockResult` structs. Oblique shock API is declared but not yet fully implemented.

### Python Bindings (`python/`)
- Built with **nanobind** (`python/src/bind_*.cpp`), exposed as `goddard._core` extension module
- One binding file per C++ domain (enums, structs, problem, combustor, nozzle, thermoarray, mixture_ratio, equilibrium, errors, moc, kinetic_nozzle)
- Each file defines a `void bind_X(nb::module_& m)` function called from `bind_main.cpp`
- Pure-Python convenience layer in `python/goddard/` (`__init__.py` re-exports, `convenience.py` has factory functions)
- No Cantera internals such as `Cantera::Solution` should be exposed in the Python API. Only Goddard and STL types (which are converted to corresponding Python types).
- Dev workflow: `pixi run compile` builds `_core.cpython-*.so` and copies it to `python/goddard/` automatically; use `PYTHONPATH=python` to import
- Install workflow: `pip install -e . --no-build-isolation` via scikit-build-core

### Key Dependencies
- **Cantera**: Primary backend for chemical kinetics and thermodynamics
- **Eigen3**: Linear algebra operations
- **nanobind**: Python bindings (C++ → Python bridge)
- **scikit-build-core**: Python package build backend for CMake projects
- **Python >=3.11**: For Python bindings and examples
- **GTest**: Unit testing framework

## Build System


### Build Commands
This project uses the Pixi package manager to manage dependencies. The goddard conda environment needs to be active before building the project. It can be activated with `pixi shell`. There are also pixi tasks that implement basic functionality (consult `pixi.toml` for details.)

```bash
# activate environment
pixi shell -e default
# clean project
pixi run clean
# configure build
pixi run configure Debug
# compile project
pixi run compile
# test project
pixi run ctest
```
Alternatively, CLI tools can be used directly as long as the correct environment is active.
```bash
pixi shell -e default
# Configure build (from project root)
mkdir -p build && cd build
cmake ..

# Build main executable and library
cmake --build .

# Build specific targets
cmake --build . --target goddard_main     # Main executable
cmake --build . --target goddard_lib      # Shared library
```
To update the Python package, use
```bash
pixi run pip-install
```
### C++ tests
```bash
# Build and run tests
cmake --build . --target goddardTests
./test/goddardTests

# Run specific test suites
ctest -R "goddard"
```

### Python binding tests
```bash
# Install as editable package
pixi run -e test pip-install

# Run Python tests
pixi run -e test test-python
```

## Data Files Structure

The `data/` directory contains:
- `nasa9.dat`: NASA polynomial thermodynamic data
- `*.yaml`: Cantera-format species and reaction data files
- `cea_results/`: Reference CEA output files for validation

Data file paths are configured via CMake and accessible through `DATA_DIR` macro in `config.h`.

## Development Patterns

### Architecture & User API
Users start by constructing a `Gas` object by specifying thermodynamic data. `Gas` is the fundamental construct and should be the most feature-rich API. Next, `Gas` is passed to different solver objects like `MocNozzle`, `Nozzle`, `ShockSolver` along with configuration options, which are then solved.  


### Formatting
- Snake case for functions and variable names, Pascal case for types.
- Prefer descriptive variable names even if this makes them a bit longer.
- Private class members start with `m_`. Public class members are named normally. 

### State Management
- Cantera `Solution` objects are wrapped in `std::shared_ptr` for safe state management
- Initial states are preserved using `saveState()/restoreState()` pattern
- Thermodynamic calculations modify underlying Cantera state objects

### Error Handling
- Custom error handling through `error.hpp/cpp`
- Convergence checking with configurable tolerances via `reltol`/`abstol` parameters
- Boolean `converged` flags in result structures

### Batch Operations
- `ThermoArray` class handles vectorized thermodynamic calculations
- Eigen arrays used for efficient numerical operations on multiple states
- Temperature/pressure arrays processed in parallel where possible

### Miscellaneous

- As the scope of the library is relatively constrained, there is little need to write 'modular' and 'extensible' code except in cases where modularity is clearly needed (e.g. output report formatting).
- Long functions are OK if it makes sense for everything in them to be computed together, and intermediate results aren't needed.


## Python Integration

### Binding architecture
- Bindings live in `python/src/bind_*.cpp` — one file per C++ header domain
- Adding a new struct field to bindings = one `.def_rw(...)` line in the corresponding `bind_*.cpp`
- `goddard_warnings` is linked PRIVATE to `goddard_lib` so consumers (bindings, tests) don't inherit `-Werror` and strict warning flags
- `python/CMakeLists.txt` additionally suppresses warnings from Python C API headers via `-Wno-*` flags
- `ThermoArray` is bound read-only (getters only); Eigen arrays auto-convert to numpy via `nanobind/eigen/dense.h`
- `RocketProblemResults.cases` accessed via `get_case(name)` / `case_names()` methods (returns references) rather than exposing the raw `unordered_map`
- `moc_design()` / `moc_analysis()` convenience functions for common MoC workflows

### Examples
Python examples and notebooks in `examples/` demonstrate:
- Basic rocket engine calculations (`basic_rocket.ipynb`)
- Chemical kinetics analysis (`kinetics.ipynb`)
- API usage patterns (`example_api.ipynb`)

## Testing Strategy

Unit tests focus on:
- Individual component functionality (equilibrium, combustor, nozzle)
- Numerical accuracy against CEA reference data
- State management and error handling
- Located in `test/` with naming pattern `test_*.cpp`
- MoC solver validation against Anderson Ch. 11 reference data and analytical solutions
- Python smoke tests and integration tests in `python/tests/`

## Code style

In general, the paradigm is "data-oriented with algebraic data types", reminiscent of Rust.

### Code organization

- Classes are used primarily as a way to manage state and organize functions.
- Standalone functions should be "pure", i.e. don't mutate state.
- Structs are used as "dumb" data containers.

### Class design
Class members are public by default unless directly modifying them can violate an invariant of the class (e.g. size attribute in a container).

Most classes are standalone. In cases where it is needed  the inheritance hierarchies are generally "shallow" (at most 1-2 layers of inheritance) and flat.

### Function dispatch and polymorphism

Dispatch is primarily achieved with enum classes (AKA discriminated unions) and switch cases. 

Create enum classes liberally, and avoid:
* Raw C enums because of namespace pollution. 
* `std::variant`

Polymorphism (virtual functions) is used in the rare event where "2D" dispatch is needed when it keeps the code simpler than multiple levels of switch cases.

### Other language features
Keep things simple and sparingly use newer language features (post C++17). `auto` is fine in moderation but prefer to be explicit with type declarations if they aren't too verbose.

However, external facing APIs should use a minimum of post-C++11 features to keep language interop and bindings simple.


## Documentation

### Documentation pipeline

Goddard uses MkDocs + Material for MkDocs + mkdocstrings to generate its documentation site. The docs source lives in `docs/` and the site config is `mkdocs.yml` at the project root.

C++ Doxygen comments are the single source of truth for API docstrings. They flow into the Python API via this pipeline:

1. CMake runs `pybind11_mkdoc` (via `scripts/generate_docstrings.py`) to parse C++ headers with libclang and emit `python/src/goddard_docstrings.h` containing `DOC(Namespace, Class, method)` macros. When `goddard_BUILD_DOCSTRINGS=OFF` (default), a stub header is emitted instead so bindings still compile.
2. nanobind binding files `#include "goddard_docstrings.h"` and pass `DOC(...)` as the docstring argument to `.def()`.
3. The `docs-stubs` pixi task runs `nanobind.stubgen` to emit `python/goddard/_core.pyi` with embedded docstrings. This runs outside CMake (as a pixi task after compilation) because stubgen needs to import the compiled module, which requires the full Python environment.
4. mkdocstrings reads the stubs and renders the Python API reference.

Both `goddard_docstrings.h` and `_core.pyi` are generated artifacts (gitignored). The `docs` pixi environment enables `goddard_BUILD_DOCSTRINGS=ON` automatically.

### Quick commands

```bash
# Build the full docs site (configure + compile + copy notebooks + mkdocs build)
pixi run -e docs docs

# Live-preview server with hot reload
pixi run -e docs docs-serve

# Deploy a versioned build to gh-pages
pixi run -e docs docs-deploy 0.1
```

### Writing docstrings

When adding or modifying public C++ APIs or Python bindings, follow the conventions in `docs/contributing/docstring_style.md`. Key points:

- Use `/** ... */` Doxygen block comments on C++ headers. Always document units.
- In binding files, use `DOC(Goddard, Class, method)` instead of literal strings.
- Pure-Python functions use Google-style docstrings.
- Do not duplicate docstrings between C++ headers and binding files.

## Instructions

### File reads

DO NOT read any files in the data/ folder unless explicitly asked to.