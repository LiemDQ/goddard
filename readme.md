# Goddard

Goddard is a software library for simulating chemically reactive flows, particularly at the high temperatures and supersonic velocities found in rocket engines. At a basic level, Goddard provides similar functionality as the "rocket" module of NASA's [Chemical Equilibrium with Applications](https://www1.grc.nasa.gov/research-and-engineering/ceaweb/) (CEA), but with a modern API and more flexibility in problem specification. 


Goddard is comprised of a C++ core, which is built on top of the [Cantera](https://cantera.org/) C++ API, and a Python wrapper. Cantera provides very similar results as NASA CEA in most situations, but if you prefer to use CEA, consider [RocketCEA](https://rocketcea.readthedocs.io/en/latest/) or [cea](https://nasa.github.io/cea/) for wrappers around the CEA Fortran code.



## Contact
Questions? Reach out by opening a discussion topic, or via email at [dq@liem.ca](mailto:dq@liem.ca). Please keep in mind that I work on Goddard in my spare time, so updates may be sporadic.


## License
Goddard is subject to the terms of the Mozilla Public License, v. 2.0.

Copyright © 2026 Liem Dam-Quang