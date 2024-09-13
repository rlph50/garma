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
            intercept       u1      fd1      ar1      ar2     ar3      ar4      ar5
             0.009440  0.49364  0.28067  0.00000  0.00000  0.5000  0.00000  0.00000
      s.e.   0.008911  0.02839  0.08278  0.06834  0.05555  0.1023  0.08147  0.06963
                ar6     ar7      ar8      ar9
            0.00000  0.0000  0.00000  0.50000
      s.e.  0.08312  0.0549  0.05389  0.06725
      
                              Factor1
      Gegenbauer frequency:    0.1678
      Gegenbauer Period:       5.9583
      Gegenbauer Exponent:     0.2807
      
      
      sigma^2 estimated as 0.0012: log likelihood = -140.053316, aic = 304.106632
      

