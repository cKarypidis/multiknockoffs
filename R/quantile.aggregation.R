#' P-values quantile aggregation.
#'
#' This function aggregates p-values by the quantile aggregation proposed by
#' Meinshausen et al. (2009).
#'
#' @param pvals B x p matrix of p-values for all variables.
#' @param B number of p-values for one variable.
#' @param gamma value between (0,1) which defines the quantile value used for the aggregation.
#'               If \code{gamma = NULL}, the adaptive search by Meinshausen et al. (2009) is used.
#'               Default: \code{NULL}.
#'
#' @return Aggregated p-values.
#'
#' @details
#'
#' The function is a helper function that it used in the p-value knockoff filter.
#'
#' @references
#'   Meinshausen, Meier and Buehlmann (2009). \emph{p-Values for High-Dimensional Regression}.
#'   Journal of the American Statistical Association 104(488), 1671-1681.
#'
#' @export
quantile.aggregation <- function(pvals, B, gamma = NULL){

  #Input checks
  #Validate input checks
  if (!is.matrix(pvals)){
    stop("Input pvals must be a matrix")
  }

  if(!B == round(B)){
    stop("B must be an integer")
  }


  #If no gamma is provides, compute sequence
  if(is.null(gamma)){ gamma <- seq(ceiling(0.05 * B) / B, 1 - 1 / B, by = 1/B)}

  p <- ncol(pvals)
  pvals.aggr <- numeric(p)

  for(j in 1:p){#loop over all features
    Qj <- quantile(pvals[,j], gamma) / gamma
    penalty <- if(length(gamma) > 1) (1 - log(min(gamma))) else 1
    pvals.pre <- min(Qj) * penalty
    pvals.aggr[j] <- pmin(pvals.pre, 1)
  }
  return(pvals.aggr)
}
