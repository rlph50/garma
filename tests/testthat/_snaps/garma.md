# garma xreg

    Code
      print(fit)
    Output
      
      Call:
      garma(x = x, xreg = m)
      
      Mean term was fitted.
      No Drift (trend) term was fitted.
      
      
      Coefficients:
            intercept       R1         R2      u1      fd1
              0.45728  0.14706  -0.009936  0.2487  0.05572
      s.e.    0.05234  0.07003   0.069077  0.1066  0.04333
      
                              Factor1
      Gegenbauer frequency:    0.2100
      Gegenbauer Period:       4.7619
      Gegenbauer Exponent:     0.0557
      
      
      sigma^2 estimated as 0.0120: log likelihood = -1112.626391, aic = 2231.252783
      
    Code
      predict(fit, n.ahead = 10, newdata = m.ahead)
    Output
      $pred
      Time Series:
      Start = 201 
      End = 210 
      Frequency = 1 
              1         2         3         4         5         6         7         8 
      0.5549148 0.4559084 0.5172535 0.4832164 0.5877493 0.4816626 0.5738843 0.5097813 
              9        10 
      0.4972928 0.5312529 
      

