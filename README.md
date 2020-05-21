# garma - R package for estimation of Gegenbauer Seasonal/Cyclical long memory processes.
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

## Overview & Introduction
This package fits a GARMA model (refer documentation) to a univariate time series.

GARMA models are extensions of ARIMA models which allow for both fractional differencing (like "fracdiff") but also allow thta to happen at a non-zero frequency in the spectrum.

This package will estimate that frequency.

At time of writing several estimation methods are supports as well as a number of (non-linear) optimisation routines.

However only k=1 models may be fit (ie only a single Gegenbauer factor) and also only non-seasonal integer differencing of d=0 or d=1 is supported (however you can always manually implement further differencing yourself prior to calling these routines).

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
Unfortunately the only documentation currently available is in the help files - after installing the package, type:
```s
> library('garma')
> help('garma')
```
