# garma - R package for support for of Gegenbauer Seasonal/Cyclical long memory processes.
## Overview & Introduction
This package fits a GARMA model (refer documentation) to a univariate time series.

GARMA models are extensions of ARIMA models which allow for both fractional differencing (like "fracdiff") but also allow thta to happen at a non-zero frequency in the spectrum.

This package will estimate that frequency.

At time of writing several estimation methods are supports as well as a number of (non-linear) optimisation routines.

However only k=1 models may be fit (ie only a single Gegenbauer factor) and also only non-seasonal integer differencing of d=0 or d=1 is supported (however you can always manually implement further differencing yourself prior to calling these routines).

## Installation.
Ensure you have the "devtools" package installed:

```
> install.packages('devtools')
```

After this you can install this package by typing:
```
> install_github('rlph50/garma')
```

## Documentation
Unfortunately the only documentation currently available is in the help files - after installing the package, type:
```
> library('garma')
> help('garma')
```
