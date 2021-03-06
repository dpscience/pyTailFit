Support this project and keep always updated about recent software releases, bug fixes and major improvements by [following on github](https://github.com/dpscience?tab=followers).

![badge-followers](https://img.shields.io/github/followers/dpscience?style=social)
![badge-stars](https://img.shields.io/github/stars/dpscience/pyTailFit?style=social)
![badge-forks](https://img.shields.io/github/forks/dpscience/pyTailFit?style=social)

# pyTailFit

![badge-language](https://img.shields.io/badge/language-Python-blue)
![badge-license](https://img.shields.io/badge/license-GPL-blue)

Copyright (c) 2021 Danny Petschke (danny.petschke@uni-wuerzburg.de). All rights reserved.<br><br>
<b>pyTailFit</b> - A simple Python program enabling tail-fitting for the analysis of lifetime spectra using least-square optimzation.

# Quickstart Guide (see example.py)

* import the 'pyTailFit' module

```python
import pyTailFit as ptf
```

* import your data (or generate it synthetically) and apply the fit using pyTailFit ...

```python
__,spectrum  = np.loadtxt('.../test-data.dat', delimiter='\t', skiprows=5, unpack=True, dtype='float')
        
time,data,fit=ptf.tail_fit(spectrum=spectrum[:],no_of_expected_decays=1,no_chn_right_of_peak=400,bin_width_in_ps=8.)
```

# How to cite this Software?

* <b>You should at least cite the applied version of this program in your study.</b><br>

You can cite all versions by using the <b>DOI 10.5281/zenodo.5520010</b>. This DOI represents all versions, and will always resolve to the latest one.<br>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.408955036.svg)](https://doi.org/10.5281/zenodo.5520010)

## ``v1.x``
<b>pyTailFit v1.0</b><br>[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.408955036.svg)](https://doi.org/10.5281/zenodo.5520011)<br>
 
 # License of pyTailFit (GNU General Public License) 
 Copyright (c) 2021 Danny Petschke (danny.petschke@uni-wuerzburg.de) All rights reserved.<br>

<p align="justify">This program is free software: you can redistribute it and/or modify<br>
it under the terms of the GNU General Public License as published by<br>
the Free Software Foundation, either version 3 of the License, or<br>
(at your option) any later version.<br><br>

This program is distributed in the hope that it will be useful,<br>
but WITHOUT ANY WARRANTY; without even the implied warranty of<br>
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.<br><br></p>

For more details see [GNU General Public License v3](https://www.gnu.org/licenses/gpl-3.0)
