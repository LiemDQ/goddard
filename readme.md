# Goddard

Goddard is a Python simulation toolkit for simulating rocket engines. It uses [Cantera](https://cantera.org/), a chemical kinetics and thermodynamics simulation suite as a backend. Goddard attempts to replicate and extend most of the 'rocket' functionality of NASA's [Chemical Equilibrium with Applications](https://www1.grc.nasa.gov/research-and-engineering/ceaweb/) (CEA), but with a native Python interface and greater flexibility for ease of integration into other programmatic analyses, and without requiring detailed knowledge of the underlying simulation libraries. See also [RocketCEA](https://rocketcea.readthedocs.io/en/latest/) for a Python wrapper around the Fortran CEA code itself.

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
Questions? Reach out by opening a discussion topic, or via email at [liem.dam-quang@mail.mcgill.ca](mailto:liem.dam-quang@mail.mcgill.ca). Please keep in mind that I work on Goddard in my spare time, and I may take some time to respond.


## License
LGPL v3.0

Copyright © 2023 Liem Dam-Quang

Goddard is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 3.0 of the License, or (at your option) any later version.

Goddard is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.