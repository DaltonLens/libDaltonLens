# libDaltonLens

[![CMake Build and Test](https://github.com/DaltonLens/libDaltonLens/actions/workflows/cmake_build_and_test.yml/badge.svg)](https://github.com/DaltonLens/libDaltonLens/actions/workflows/cmake_build_and_test.yml)

Single-file, public domain implementation of color blindness / color vision
deficiency (CVD) simulation in C.

There are no dependencies, just include the C header and source file in your
project and call:

```
void dl_simulate_cvd (enum DLDeficiency deficiency, 
                      float severity,
                      unsigned char *srgba_image, 
                      size_t width, 
                      size_t height, 
                      size_t bytesPerRow);
```

Two algorithms are implemented:

- _Computerized simulation of color appearance for dichromats_ by Brettel, H.
  and Viénot, F. and Mollon, J. D. (1997). This is the most accurate approach,
  but it's slightly more expensive, so it's only chosen by default for
  tritanopia.

- _Digital video colourmaps for checking the legibility of displays by
  dichromats_ by Viénot, F. and Brettel, H. and Mollon, J. D. (1999). Their
  followup paper is a bit faster for protanopia and tritanopia and has a similar
  accuracy, so it's used by default for these deficiencies.

This implementation is meant to be an easy-to-read reference that can be easily
copy/pasted/adapted/optimized in other languages or contexts. 

It is unit-tested against
[DaltonLens-Python](https://github.com/DaltonLens/DaltonLens-Python), which is
itself unit tested against external references. The precomputed matrices /
values are also generated from [a notebook in
DaltonLens-Python](https://github.com/DaltonLens/DaltonLens-Python/blob/master/research/for-desktop/precomputed_matrices.ipynb),
so you can easily change the LMS model by running it with different parameters.

For a discussion about which CVD simulation algorithms are the most accurate see
our [Review of Open Source Color Blindness
Simulations](https://daltonlens.org/opensource-cvd-simulation/).

For more information about the math of the chosen algorithms see our article
[Understanding CVD
Simulation](https://daltonlens.org/understanding-cvd-simulation/).

This library is part of the [DaltonLens project](https://daltonlens.org).
## SVG Filters

A port as SVG filter is available in the svg subfolder. It includes an
implementation of Brettel et al., which is non-trivial. [Live demo](https://daltonlens.github.io/libDaltonLens/svg/cvd_svg_filters.html).
