#' @title ADAGES knockoff filter
#'
#' @description
#' This function runs the whole ADAGES procedure in the multiple knockoff setting, i.e.
#' it generates multiple knockoff matrices, estimates the score functions and the selection sets
#' of multiple knockoff runs, which are then aggregated by ADAGES to obtain the final selection set.
#'
#' @param X n x p matrix or data frame of original variables.
#' @param y response vector of length n.
#' @param knockoffs function for the knockoff construction. It must take the n x p matrix as input
#'                  and it must return a n x p knockoff matrix. Either choose a knockoff sampler of
#'                  the \code{knockoff} package or define it manually. Default: \code{create.second_order} (see below).
#' @param statistic function that computes the score vector \eqn{W} of length p. It must take the data matrix,
#'                  knockoff matrix and response vector as input and outputs a vector of computed
#'                  scores. Either choose one score statistic from the \code{knockoff} package or
#'                  define it manually. Default: \code{stat.glmnet_coefdiff} (see below).
#' @param q nominal level for the FDR control. Default: 0.2.
#' @param K number of knockoff runs. Default: 5.
#' @param offset either 0 (knockoff) or 1 (knockoff+). Default: 1.
#' @param type either \code{"ADAGES"} (default) or \code{"ADAGES.mod"} (see below).
#' @param sets logical argument if the K selection sets of each knockoff run
#'            should be returned. Default: \code{FALSE}.
#'
#' @return A list containing following components:
#'  \item{Shat}{aggregated selection set.}
#'  \item{c}{optimal threshold.}
#'  \item{K}{number of knockoff runs.}
#'  \item{sets}{if specified, individual selection sets of each knockoff run.}
#'
#' @details
#'
#' This function requires the installation of the \code{knockoff} package prior to its execution.
#'
#' The default knockoff sampler \code{create.second_order} is the second-order Gaussian knockoff construction from
#' the \code{knockoff} package.
#'
#' The default score function \code{stat.glmnet_coefdiff} is from the \code{knockoff} package.
#' It fits a  Lasso regression where the regularization parameter \eqn{\lambda} is tuned by cross-validation.
#' Then, the score is computed as the difference between
#'   \deqn{W_j = |Z_j| - |\tilde{Z}_j|}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the coefficient estimates for the
#' jth variable and its knockoff, respectively.
#'
#' \code{ADAGES} applies the minimization of the complexity ratio as a criterion to determine
#' the optimal threshold.
#'
#' \code{ADAGES.mod} minimizes the trade-off between the threshold and the model complexity \eqn{c |S|}
#' to determine the optimal threshold.
#'
#' @references
#'   Gui (2020). \emph{ADAGES: adaptive aggregation with stability for distributed feature selection}.
#'   Proceedings of the 2020 ACM-IMS on Foundations of Data Science Conference.
#'   \url{https://arxiv.org/pdf/2007.10776.pdf}
#'
#' @examples
#' n <- 400; p <- 200; s_0 <- 30
#' amplitude <- 1; mu <- rep(0,p); rho <- 0.25
#' Sigma <- toeplitz(rho^(0:(p-1)))
#'
#' X <- MASS::mvrnorm(n, mu, Sigma)
#' nonzero <- sample(p, s_0)
#' beta <- amplitude * (1:p %in% nonzero)
#' y <- X %*% beta + rnorm(n)
#'
#' # Basic usage with default arguments
#' res.ADAGES <- run.ADAGES(X, y, sets = TRUE)
#' res.ADAGES
#'
#' # Advanced usage with customized knockoff construction (equi-correlated)
#' equi.knock <- function(X) create.second_order(X, method = "equi")
#' res.ADAGES <- run.ADAGES(X, y, knockoffs = equi.knock, sets = TRUE)
#' res.ADAGES
#'
#' @export
run.ADAGES <- function(X, y,
                      knockoffs = create.second_order,
                      statistic = stat.glmnet_coefdiff,
                      q = 0.2, K = 5, offset = 1,
                      type = "ADAGES", sets = FALSE){

  library(knockoff)

  #Validate input checks
  if (is.data.frame(X)) {
    X.names = names(X)
    X = as.matrix(X, rownames.force = F)
  } else if (is.matrix(X)) {
    X.names = colnames(X)
  } else {
    stop('Input X must be a numeric matrix or data frame')
  }
  if (!is.numeric(X)) stop('Input X must be a numeric matrix or data frame')


  if (!is.factor(y) && !is.numeric(y)) {
    stop('Input y must be either of numeric or factor type')
  }
  if(is.numeric(y)) y = as.vector(y)


  n = nrow(X); p = ncol(X)
  stopifnot(length(y) == n)


  if(offset!=1 && offset!=0) {
    stop('Input offset must be either 0 or 1')
  }

  if(q < 0 | q > 1) {
    stop('q must be between 0 and 1')
  }

  #For New part Checks
  if (!is.function(knockoffs)) stop('Input knockoffs must be a function')
  if (!is.function(statistic)) stop('Input statistic must be a function')

  if(!K == round(K)){
    stop("K must be an integer")
  }

  #Knockoff construction
  Xk <- multi.knockoffs(X = X, K = K, knockoffs = knockoffs)

  #Knockoff filters
  multi.filter <- multi.knockfilter(X, Xk, y, q = q, offset = offset,
                                  statistic = statistic)


  Shat.list <- multi.filter$Shat.list

  if(type == "ADAGES"){
    res <- agg.ADAGES(Shat.list = Shat.list, p = p)
    }
  else if(type == "ADAGES.mod"){
    res <- agg.ADAGES.mod(Shat.list = Shat.list, p = p)
    }
  else{
    stop("Unknown type of ADAGES method")
    }

  if(sets == TRUE){
    res$sets <- Shat.list
  }

  return(res)
}
