# garma - R package for estimation of Gegenbauer Seasonal/Cyclical long memory processes.
[![CRAN Version](https://img.shields.io/cran/v/garma)](https://img.shields.io/cran/v/garma)
[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

## Overview & Introduction
This package fits a GARMA model (refer documentation) to a univariate time series.

GARMA models are extensions of ARIMA models which allow for both fractional differencing (like "fracdiff") but also allow that to happen at a non-zero frequency in the spectrum.

This package will estimate that frequency (which is known for technical reasons as the "Gegenbauer" frequency).

At time of writing several estimation methods are supports as well as a number of (non-linear) optimisation routines.

  ## Installation.
Ensure you have the "devtools" package installed:

```s
> install.packages('devtools')
```

After this you can install this package by typing:
```s
> remotes::install_github('rlph50/garma')
```

## Documentation
An Introduction to the "garma" packages is available [here](https://github.com/rlph50/garma/blob/master/inst/vignette_introduction.pdf), and the reference documentation is available [here](https://github.com/rlph50/garma/blob/master/inst/garma_0.9.0.pdf).
