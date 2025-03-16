# Check Output

    Code
      fit <- garma(dap, order = c(9, 0, 0), k = 1)
      print(fit)
    Output
      
      Call:
      garma(x = dap, order = c(9, 0, 0), k = 1)
      
      Mean term was fitted.
      No Drift (trend) term was fitted.
      
    Condition
      Warning in `internal_print_garma_model()`:
      model estimates are not Stationary! Forecasts may become unbounded.
    Output
      
      Coefficients:
            intercept       u1      fd1        ar1         ar2     ar3         ar4
             0.009440  0.49364  0.28067  2.277e-05  -2.157e-05  0.5000  -1.518e-05
      s.e.   0.008911  0.02842  0.08278  5.898e+01   5.898e+01  0.1028   2.949e+01
                   ar5         ar6        ar7         ar8      ar9
            -1.933e-05  -3.622e-05  2.893e-06  -4.449e-05  0.50002
      s.e.   2.949e+01   8.371e-02  2.950e+01   2.950e+01  0.07065
      
                              Factor1
      Gegenbauer frequency:    0.1678
      Gegenbauer Period:       5.9583
      Gegenbauer Exponent:     0.2807
      
      
      sigma^2 estimated as 0.0012: log likelihood = -140.051963, aic = 304.103927
      

