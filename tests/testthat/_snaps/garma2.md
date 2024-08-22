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
            intercept       u1      fd1        ar1        ar2     ar3         ar4
             0.009440  0.49364  0.28067  1.113e-19  1.368e-16  0.5000  -7.658e-17
      s.e.   0.008911  0.02839  0.08278  6.834e-02  5.555e-02  0.1023   8.147e-02
                  ar5         ar6         ar7        ar8      ar9
            4.235e-17  -3.587e-17  -1.089e-16  3.520e-17  0.50000
      s.e.  6.963e-02   8.312e-02   5.490e-02  5.389e-02  0.06725
      
                              Factor1
      Gegenbauer frequency:    0.1678
      Gegenbauer Period:       5.9583
      Gegenbauer Exponent:     0.2807
      
      
      sigma^2 estimated as 0.0012: log likelihood = -140.053316, aic = 304.106632
      

