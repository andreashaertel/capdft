<!--
SPDX-FileCopyrightText: 2008-2021 Andreas Härtel <http://andreashaertel.anno1982.de/>

SPDX-License-Identifier: CC-BY-SA-4.0
-->

# capDFT

The capdft project provides a C++ library and several applications to perform 
calculations in the framework of classical density functional theory.

The project started in 2008 with a research project by Andreas Härtel in the context of a Diploma thesis. 

Further details will be found on capdft.org soon. 


## Install

First, check that the required libraries are installed:
1. GNU Scientific Library, e.g. `gsl-devel`
2. FFTW3, e.g. `fftw3-devel`
3. Doxygen
4. Parameter handler from https://github.com/andreashaertel/parameter_handler.git
(The latter should be located in the same folder as the capdft project.)

To download and install the project, you need to clone the project from github and to run `make init` and `make`: 
```bash
git clone https://github.com/andreashaertel/capdft.git
cd capdft/
make init
make
```

To make the documentation, you need to run make in the doc/ subfolder: 
```bash
cd doc/
make
```

See doc/html/index.html for further details. 


## Source code

The source code contains implementations for multiple excess free energy
functionals that can be used for density functional theory calculations
in different contexts, a Picard iteration class for the numerical solution
of the Euler-Lagrange equation based on those functionals, and many
related classes and functions which are employed on the way from the definition of the system to the 
analysis of the calculated density profiles.

Depending on the application of interest, choose specific classes for your calculations:
1. For different spatial symmetries of the system, different coordinate systems
might be chosen. This library contains implementations for effectively 1d
systems (corresponding files and classes labeled by 'planar'), cartesian 3d coordinates ('cartesian') and spherical coordinates
('spherical' or 'radial').
2. The functionals can treat hard sphere systems by the use of fundamental measure theory ('fmt')
or electrically charged particles ('es') by mean field approximations ('mf' or 'delta'). These can also
be combined for treating charged hard spheres.
3. You can modify the boundary positions or introduce additional surfaces / fixed objects
into the system by using the PoissonSolvers with appendix 'any' and BoundarySurface objects.

## Examples

Extensive examples are given in the examples directory. Code snippets demonstrating the
usage of single classes and functions and tests that can be used to check their
correct functioning can be found in the tests directory.


## Maintainers

- Andreas Härtel - <andreas.haertel@anno1982.de>
- Moritz Bültmann - <mbueltmann_uni@posteo.net>


## Contribute

Any pull requests or suggestions are welcome at 
<https://github.com/andreashaertel/capdft> or via e-mail to one of the maintainers. 


## License

This work is licensed under multiple licences. Because keeping this section 
up-to-date is challenging, here is a brief summary as of May 2021: 

- All original source code is licensed under either LGPL-3.0-or-later or GPL-3.0-or-later. 
- All documentation is licensed under CC-BY-SA-4.0. 
- Some configuration and data files are licensed under CC0-1.0. 
- Some code borrowed from 
  [anno1982/parameter...](https://github.com/... !!PARAMETER FILE READER) is licensed under
  LGPL-3.0-or-later. 

For more accurate information, check the individual files.



