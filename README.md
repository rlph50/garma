# garma - R code for support for of Gegenbauer Seasonal/Cyclical long memory processes
## Estimate the parameters of a GARMA model.

The garma function is the main function for the garma package. Depending on the parameters it will
calculate the parameter estimates for the GARMA process, and if available the standard errors (se's)
for those parameters.

The GARMA model is specified as
\begin{equation*}
\phi(B)\prod_{i=1}^{k}(1-2u_{i}B+B^{2})^{d_{i}}(X_{t}-\mu)= \theta(B) \epsilon _{t}}}{\prod(i=1 to k) (1-2u(i)B+B^2)^d(i) \phi(B) X(t) = \theta(B) \epsilon(t)}
\end{equation*}

where
* \eqn{\phi(B)}{\phi(B)} represents the short-memory Autoregressive component of order p,
* \eqn{\theta(B)}{\theta(B)} represents the short-memory Moving Average component of order q,
* \eqn{(1-2u_{i}B+B^{2})^{d_{i}}}{(1-2u(i)B+B^2)^d(i)} represents the long-memory Gegenbauer component (there may in general be k of these),
* \eqn{X_{t}}{X(t)} represents the observed process,
* \eqn{\epsilon_{t}}{\epsilon(t)} represents the random component of the model - these are assumed to be uncorrelated but identically distributed variates.
      Generally the routines in this package will work best if these have an approximate Gaussian distribution.
* \eqn{B}{B} represents the Backshift operator, defined by \eqn{B X_{t}=X_{t-1}}{B X(t) = X(t-1)}.

when k=0, then this is just a short memory model as fit by the stats "arima" function.
