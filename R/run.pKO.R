#' @title P-value knockoffs
#'
#' @description This function runs the whole p-value knockoff procedure, i.e. it generates multiple knockoff matrices,
#'              estimates the scores, computes intermediate p-values and aggregates them before applying
#'              Benjamini-Hochberg or Benjamini-Yekutieli in the last step to obtain the final selection set.
#'
#' @param X n x p matrix of original variables.
#' @param y response vector of length n.
#' @param knockoffs function for the knockoff construction. It must take the n x p matrix as input
#'                  and it must return a n x p knockoff matrix. Either choose a knockoff sampler of
#'                  the \code{knockoff} package or define it manually. Default: \code{create.second_order} (see below).
#' @param statistic function that computes the score vector \eqn{W} of length p. It must take the data matrix,
#'                  knockoff matrix and response vector as input and outputs a vector of computed
#'                  scores. Either choose one score statistic from the \code{knockoff} package or
#'                  define it manually. Default: \code{stat.glmnet_coefdiff} (see below).
#' @param q nominal level for the FDR control. Default: 0.2.
#' @param B number of knockoff runs. Default: 25.
#' @param gamma value between (0,1) which defines the quantile value used for the aggregation.
#'               If \code{gamma = NULL}, the adaptive search by Meinshausen et al. (2009) is used. Default: 0.3.
#' @param offset either 0 (knockoff) or 1 (knockoff+). Default: 1.
#' @param method the FDR controlling method in the last step. Either \code{"BH"} (default) or \code{"BY"}.
#' @param pvals logical argument if the aggregated p-values should be reported. Default: \code{FALSE}.
#'
#' @return A list containing following components:
#'  \item{Shat}{aggregated selection set.}
#'  \item{B}{number of knockoff matrices.}
#'  \item{pvals}{if specified, vector of aggregated p-values.}
#'
#' @details
#'
#' This function requires the installation of the \code{knockoff} package prior to its execution.
#'
#' The default knockoff sampler \code{create.second_order} is the second-order Gaussian knockoff construction from
#' the \code{knockoff} package.
#'
#' Although the default knockoff sampler is based on ASDP, we recommend using the equi-correlated
#' construction within \code{create.second_order} because it performs significantly better. See the example below
#' that shows how the user can change the knockoff sampler \code{create.second_order} to create equi-correlated knockoffs.
#'
#' The default score function \code{stat.glmnet_coefdiff} is from the \code{knockoff} package.
#' It fits a  Lasso regression where the regularization parameter \eqn{\lambda} is tuned by cross-validation.
#' Then, the score is computed as the difference between
#'   \deqn{W_j = |Z_j| - |\tilde{Z}_j|}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the coefficient estimates for the
#' jth variable and its knockoff, respectively.
#'
#'
#' @references
#'   Benjamini and Hochberg (1995). \emph{Controlling the False Discovery Rate: A Practical
#'   and Powerful Approach to Multiple Testing}. Journal of the Royal Statistical Society. Series B (Methodological) 57(1), 289-300.
#'
#'   Benjamini and Yekutieli (2001). \emph{The control of the false discovery rate in multiple testing under dependency}.
#'   The Annals of Statistics 29(4), 1165-1188.
#'
#'   Meinshausen, Meier and Buehlmann (2009). \emph{p-Values for High-Dimensional Regression}.
#'   Journal of the American Statistical Association 104(488), 1671-1681.
#'
#'   Nguyen, Chevalier, Thirion and Arlot (2020). \emph{Aggregation of Multiple Knockoffs}.
#'   Proceedings of the 37th International Conference on Machine Learning.
#'   \url{https://arxiv.org/abs/2002.09269}
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
#' res.pKO <- run.pKO(X, y, pvals = T)
#' res.pKO
#'
#' # Advanced usage with customized knockoff construction (equi-correlated)
#' equi.knock <- function(X) create.second_order(X, method = "equi")
#' res.pKO <- run.pKO(X, y, knockoffs = equi.knock, pvals = T)
#' res.pKO
#'
#' @export
run.pKO <-  function(X, y,
                     knockoffs = create.second_order,
                     statistic = stat.glmnet_coefdiff,
                     q = 0.2, B = 25, gamma = 0.3,
                     offset = 1, method = "BH", pvals = F, ...){

  library(knockoff)

  #Validate input checks
  if (!is.matrix(X)){
    stop("Input X must be a matrix")
  }


  if (!is.factor(y) && !is.numeric(y)) {
    stop('Input y must be either of numeric or factor type')
  }
  if( is.numeric(y) ) y = as.vector(y)

  n = nrow(X); p = ncol(X)
  stopifnot(length(y) == n)

  if(offset!=1 && offset!=0) {
    stop('Input offset must be either 0 or 1')
  }

  #For New part Checks
  if (!is.function(knockoffs)) stop('Input knockoffs must be a function')
  if (!is.function(statistic)) stop('Input statistic must be a function')

  if(!B == round(B)){
    stop("B must be an integer")
  }

  #Knockoff construction
  Xk <- multi.knockoffs(X = X, K = B, knockoffs = knockoffs)

  #Knockoff filters
  multi.filter <- multi.knockfilter(X, Xk, y, q = q, offset = offset,
                                  statistic = statistic)


  W.list <- multi.filter$W.list


  pKO.res <- agg.pKO(W.list = W.list, q=q, gamma = gamma, offset = offset, method = method,
                 pvals = pvals)


  return(pKO.res)

}
