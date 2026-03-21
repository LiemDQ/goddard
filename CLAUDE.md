# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Goddard is a C++/Python rocket engine simulation toolkit that uses Cantera for chemical kinetics and thermodynamics. It aims to replicate and extend NASA CEA (Chemical Equilibrium with Applications) functionality with a modern API. The project is structured as a hybrid C++/Python codebase with CMake build system.

## Core Architecture

### C++ Core Components
- **Equilibrium** (`equilibrium.hpp/cpp`): Chemical equilibrium calculations and thermodynamic derivatives
- **Combustor** (`combustor.hpp/cpp`): Isobaric combustion reaction handling with support for infinite area, finite mass flux, and finite contraction ratio modes
- **Nozzle** (`nozzle.hpp/cpp`): Nozzle flow calculations with inheritance hierarchy:
  - `NozzleBase`: Abstract base class
  - `EquilibriumNozzle`: Chemical equilibrium nozzle flow
  - `FrozenNozzle`: Frozen composition nozzle flow
- **ThermoArray** (`thermoarray.hpp/cpp`): Batch thermodynamic property calculations
- **MixtureRatio** (`mixture_ratio.hpp/cpp`): Fuel/oxidizer mixture ratio handling

### Python Bindings (`python/`)
- Built with **nanobind** (`python/src/bind_*.cpp`), exposed as `goddard._core` extension module
- One binding file per C++ domain (enums, structs, problem, combustor, nozzle, thermoarray, mixture_ratio, equilibrium, errors)
- Each file defines a `void bind_X(nb::module_& m)` function called from `bind_main.cpp`
- Pure-Python convenience layer in `python/goddard/` (`__init__.py` re-exports, `convenience.py` has factory functions)
- `Cantera::Solution` exposed as opaque `shared_ptr` handle (`SolutionHandle`) — no Cantera internals in the Python API
- Cantera Python interop via state reconstruction (`from_cantera()` in `convenience.py`), not pointer extraction (Cython bindings don't expose `shared_ptr`)
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

### Testing
```bash
# Build and run tests
cmake --build . --target goddardTests
./test/goddardTests

# Run specific test suites
ctest -R "goddard"
```

### Python Bindings
```bash
# Build (included in normal compile)
pixi run compile

# Install as editable package
pixi run pip-install

# Run Python tests
pixi run test-python
```

## Data Files Structure

The `data/` directory contains:
- `nasa9.dat`: NASA polynomial thermodynamic data
- `*.yaml`: Cantera-format species and reaction data files
- `cea_results/`: Reference CEA output files for validation

Data file paths are configured via CMake and accessible through `DATA_DIR` macro in `config.h`.

## Development Patterns

### Architecture
Goddard is rocket science, so its codebase shouldn't be. 
- The code should be as simple and readable as possible even at the cost of some repetition. 
- Avoid overabstraction and excessive use of inheritance.
- Prefer composition over inheritance. 
- As the scope of the library is relatively constrained, there is little need to write 'modular' and 'extensible' code except in cases where modularity is clearly needed (e.g. output report formatting).
- Long functions are OK if it makes sense for everything in them to be computed together, and intermediate results aren't needed.

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

## Python Integration

### Binding architecture
- Bindings live in `python/src/bind_*.cpp` — one file per C++ header domain
- Adding a new struct field to bindings = one `.def_rw(...)` line in the corresponding `bind_*.cpp`
- `goddard_warnings` is linked PRIVATE to `goddard_lib` so consumers (bindings, tests) don't inherit `-Werror` and strict warning flags
- `python/CMakeLists.txt` additionally suppresses warnings from Python C API headers via `-Wno-*` flags
- `ThermoArray` is bound read-only (getters only); Eigen arrays auto-convert to numpy via `nanobind/eigen/dense.h`
- `RocketProblemResults.cases` accessed via `get_case(name)` / `case_names()` methods (returns references) rather than exposing the raw `unordered_map`

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

## Code style

In general, the paradigm is "data-oriented with algebraic data types", reminiscent of Rust.

### Code organization

- Classes are used primarily as a way to manage state and organize functions.
- Standalone functions are "pure", i.e. don't mutate state.
- Structs are used as "dumb" data containers.

### Class design
Class members are public by default unless directly modifying them can violate an invariant of the class (e.g. size attribute in a container).

Minimize inheritance. Most classes are standalone. In cases where it is needed  the inheritance hierarchies are generally "shallow" (at most 1-2 layers of inheritance) and flat.

### Function dispatch and polymorphism

Dispatch is primarily achieved with enum classes (AKA discriminated unions) and switch cases. 

Create enum classes liberally, and avoid:
* Raw C enums because of namespace pollution. 
* `std::variant`: it's bloated, makes compilation and debugging harder and [the performance can be mediocre depending on the stdlib implementation](https://stackoverflow.com/questions/57726401/stdvariant-vs-inheritance-vs-other-ways-performance). `std::visit` often results in verbose and hard-to-understand code. Using it properly often requires defining custom functors and overload helper classes 

Polymorphism (virtual functions) is used in the rare event where "2D" dispatch is needed as it keeps the code simpler than multiple levels of switch cases, in my opinion.

### Generics
Templates are fine in moderation, and are the right tool in many contexts, particularly library code. Using heavily templated libraries like Boost or Eigen is acceptable when appropriate.

But in practice maintaining 2-3 overloaded function signatures is often a more pragmatic choice and it's rare that you'll need to support more types than that unless you're writing e.g. a numerics library.

Using C++20 Concepts in library code is a good idea if appropriate. 

The STL is mostly fine, but stick to a handful of things that are useful and practical: `vector`, `string`, `memory`, `utility` (mostly `std::move`), and `iostream`. `unordered_map` if a quick-and-dirty hashmap is needed and performance is unimportant. Other things are outside of this are acceptable on occasion if it makes a lot of sense and saves a lot of time.

### Memory management
Raw pointers should rarely be used. Modifying state is achieved via member functions. Avoiding copies is achieved via const references. 
Having to modify two different objects in a single function call is usually a code smell.

Using `shared_ptr` everywhere is an anti pattern. I find that I rarely need to use `shared_ptr` or really direct heap allocations in general. When I do, `unique_ptr` is often sufficient. In my opinion, many instances of heap allocation are due to a need for polymorphic dispatch and that simply doesn't happen often when enum-based dispatch is the default.

In practice, most dynamic memory allocation use cases are covered by containers like `vector`.

### Other language features
Keep things simple and sparingly use newer language features (post C++17). `auto` is fine in moderation but prefer to be explicit with type declarations if they aren't too verbose.

However, external facing APIs should use a minimum of post-C++11 features to keep language interop and bindings simple.


## Instructions

### File reads

DO NOT read any files in the data/ folder unless explicitly asked to.