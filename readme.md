# Goddard

Goddard is a Python simulation toolkit for designing rocket engines. It features a wrapper around [Cantera](https://cantera.org/), a chemical kinetics and thermodynamics simulation suite for analyzing rocket combustion. Goddard attempts to replicate and extend most of the functionality of NASA's [Chemical Equilibrium with Applications](https://www1.grc.nasa.gov/research-and-engineering/ceaweb/) software, but with a native Python interface and greater flexibility for ease of integration into other programmatic analyses, and without requiring detailed knowledge of the underlying simulation libraries. See also [RocketCEA](https://rocketcea.readthedocs.io/en/latest/) for a Python wrapper around the Fortran CEA code itself.

## Features

### Roadmap
- [ ] Support for analyzing hybrid engine performance
- [ ] Support for analyzing liquid engine performance
- [ ] Support for analyzing solid engine performance
- [ ] Export motors in a .eng format 
- [ ] Integration with RocketPy

## Cantera vs. CEA
Cantera is a very popular suite of tools used to simulate and analyze combustion, and comes with its own object-oriented Python interface.

Some of the challenges adapting Cantera to rocket engine analysis are discussed in this [set of notes](https://kyleniemeyer.github.io/rocket-propulsion/thermochemistry/cea_cantera.html) by Kyle Niemeyer.

## Contact
Quest