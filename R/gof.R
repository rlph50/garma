## ' Goodness-of-Fit test for a garma_model.
## This function commented out for now. Does not appear to work well.
## '
## ' Provides a goodness-of-fit test for a GARMA Model, using the method of Delgado, Hidalgo and Velasco (2005).
## '
## ' This routine provides a test for White noise using eqn (6) of Delgado, Hidalgo and Velasco (2005). The
## ' statistic calculated is \eqn{\alpha_n^0=\frac{1}{\sqrt{\tilde{n}}}\left(\frac{G_n^0(\lambda_j)}{G_T^0(\pi)}-\frac{\lambda_j}{\pi}\right)}
## ' which is asymptotically distributed as a Brownian Bridge,
## ' where \eqn{G_n^0(\lambda_j)=\frac{2\pi}{\tilde{n}}\sum_{j=1}^{[\tilde{n}\lambda_j/\pi]}I_\epsilon(\lambda_j)}
## '
## ' So any interval \eqn{[0,\lambda]}
## ' will have an asymptotically \eqn{N\left(0,\frac{\lambda}{\pi}\left(1-\frac{\lambda}{\pi}\right)\right)} distribution.
## '
## ' @param object (garma_model) The garma_model to test.
## ' @param gof.lag (int) max lag to test.
## ' @return None.
## ' @export
# gof<-function(object,gof.lag=10) {
#   r <- as.numeric(residuals(object))
#   n <- length(r)
#   tilde_n <- as.integer(n/2)
#   if (gof.lag>tilde_n) stop('ERROR: number of lags not supported by length of data.\n')
#
#   I_eps <- spec.pgram(r,taper=0,fast=FALSE,demean=FALSE,detrend=FALSE,plot=FALSE)
#   G_lambda <- sum(I_eps$spec[1:gof.lag])
#   G_pi <- sum(I_eps$spec)
#   # in next line the (2*gof.lag/n) is a simplification:
#   # lambda_{gof.lag}/pi = (2*pi*gof.lag/n)/pi = (2*gof.lag/n)
#   s <- 1/sqrt(tilde_n)*(G_lambda/G_pi-(2*gof.lag/n))
#   var_s <- (2*gof.lag/n)*(1-2*gof.lag/n)
#   cat(sprintf('lambda/pi %f\nG_lambda %f\nG_pi %f\ns %f\nvar_s %f\n',2*gof.lag/n,G_lambda,G_pi,s,var_s))
#   s2 <- s/sqrt(var_s)
#   p <- 1-pnorm(s2)
#
#   cat(sprintf('\nTest for H0: Residuals up to lag %d are White Noise.\napprox z-statistic: %0.4f\np-value:     %0.4f\n',gof.lag,s2,p))
# }
