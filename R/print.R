#' summarise a garma_model object.
#'
#' The summary function provides a summary of a "garma_model" object.
#' @param object (garma_model) The garma_model from which to print the values.
#' @param ... Other arguments. Ignored.
#' @examples
#' data(AirPassengers)
#' ap  <- as.numeric(diff(AirPassengers,12))
#' mdl <- garma(ap,order=c(9,1,0),k=0,method='CSS',include.mean=FALSE)
#' summary(mdl)
#' @export
summary.garma_model<-function(object,...) {
  .print_garma_model(object,verbose=FALSE)
}

#' print a garma_model object.
#'
#' The print function prints a summary of a "garma_model" object.
#' @param x (garma_model) The garma_model from which to print the values.
#' @param ... Other arguments. Ignored.
#' @examples
#' data(AirPassengers)
#' ap  <- as.numeric(diff(AirPassengers,12))
#' mdl <- garma(ap,order=c(9,1,0),k=0,method='CSS',include.mean=FALSE)
#' print(mdl)
#' @export
print.garma_model<-function(x,...) {
  .print_garma_model(x,verbose=TRUE)
}

.print_garma_model<-function(mdl,verbose=TRUE) {
  cat("\nCall:", deparse(mdl$call, width.cutoff = 75L), "", sep = "\n")
  if (!mdl$include.mean) cat('No ')
  cat('Mean term was fitted.\n')
  if (!mdl$include.drift) cat('No ')
  cat('Drift (trend) term was fitted.\n\n')
  if (verbose) {
    with(mdl,
         cat(sprintf('Summary of a Gegenbauer Time Series model.\n\nFit using %s method.\nOrder=(%d,%d,%d) k=%d %s\n\nOptimisation.\nMethod:  %s\nMaxeval: %d\n',
                     method,order[1],order[2],order[3],k,ifelse(mdl$method=='QML',sprintf('QML Truncation at %d',mdl$m_trunc),''),
                     paste(mdl$opt_method,collapse=', '),mdl$maxeval))
    )
    if (mdl$opt_method[[1]]=='best') cat(sprintf('Best optimisation method selected: %s\n',mdl$opt_method_selected))
    cat(sprintf('Convergence Code: %d\nOptimal Value found: %0.8f\n\n',mdl$convergence,mdl$obj_value))
  }
  if (mdl$opt_method[[1]]=='solnp'&mdl$convergence!=0) cat('ERROR: Convergence not achieved. Please try another method.\n')
  if (mdl$convergence<0) cat(sprintf('Model did not converge.\n\n',mdl$conv_message))
  else {
    if (mdl$convergence>0)
      cat(sprintf('WARNING: Only partial convergence achieved!\n%s reports: %s (%d)\n\n',
                  ifelse(mdl$opt_method=='best',mdl$opt_method_selected,mdl$opt_method),mdl$conv_message,mdl$convergence))
    cat('Coefficients:\n')
    print.default(mdl$coef, print.gap=2, digits=4)
    cat('\n')

    if (mdl$k>0) print(mdl$ggbr_factors)

    if (mdl$sigma2>0) {
      cat(sprintf('\n\nsigma^2 estimated as %0.4f',mdl$sigma2))
      if (mdl$method %in% c('CSS','QML','Whittle')) cat (': ')
    }
    if (mdl$method=='CSS') cat(sprintf('part log likelihood = %f',mdl$loglik))
    if (mdl$method=='QML') cat(sprintf('log likelihood = %f',mdl$loglik))
    if (mdl$method=='Whittle') cat(sprintf('log likelihood = %f, aic = %f',mdl$loglik, mdl$aic))
    cat('\n')
  }
}


