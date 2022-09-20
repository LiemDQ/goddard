# goddard

goddard is a wrapper around [Cantera](https://cantera.org/), a chemical kinetics and thermodynamics simulation suite, to adapt it for rocket engine combustion simulations. goddard attempts to replicate most of the functionality of NASA's [Chemical Equilibrium with Applications](https://www1.grc.nasa.gov/research-and-engineering/ceaweb/) software, but with a native Python interface and greater flexibility for ease of integration into other programmatic analyses. See also [RocketCEA](https://rocketcea.readthedocs.io/en/latest/) for a Python wrapper around the Fortran CEA code itself.

## Features


## Why Cantera?
Cantera is a very popular suite of tools used to simulate and analyze combustion, and comes with its own Python interface.  

Some of the challenges adapting Cantera to rocket engine analysis are discussed in this [set of notes](https://kyleniemeyer.github.io/rocket-propulsion/thermochemistry/cea_cantera.html) by Kyle Niemeyer. 