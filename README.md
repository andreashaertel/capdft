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

To download and install the project, you need to clone from github and to run make install and make: 
```bash
git clone git@github.com:andreashaertel/capdft.git
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



## Examples

Examples are given in the examples directory. 




## Maintainers

- Andreas Härtel - <andreas.haertel@anno1982.de>



## Contribute

Any pull requests or suggestions are welcome at 
<https://github.com/andreashaertel/capdft> or via e-mail to one of the maintainers. 

Explain how to contribute ... . 



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



