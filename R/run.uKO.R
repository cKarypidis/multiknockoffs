#' Union knockoff filter
#'
#' This function runs the whole union knockoff procedure, i.e. it generates multiple
#' knockoff matrices, estimates the score functions and the selection sets of multiple knockoff
#' runs, which are then aggregated by their union to obtain the final selection set.
#'
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
#' @param qk  sequence of nominal levels. Either choose \code{"decseq"} (default) for \eqn{q_{k} = q/2^{k-1}}
#'            or \code{"ave"} for \eqn{q_{k} = q/K}.
#' @param q nominal level for the FDR control. Default: 0.2.
#' @param K number of knockoff runs. Default: 5.
#' @param q_seq  manual sequence of nominal level which has to match in length
#'               with the number of knockoff runs \code{K}. If this argument is specified, \code{qk} and \code{q} are
#'               ignored.
#' @param offset either 0 (knockoff) or 1 (knockoff+). Default: 1.
#' @param sets logical argument if the K selection sets of each knockoff run
#'            should be returned. Default: \code{FALSE}.
#'
#' @return A list containing following components:
#'  \item{Shat}{aggregated selection set.}
#'  \item{K}{number of knockoff runs.}
#'  \item{FDRbound}{theoretical FDR bound.}
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
#' The user has to specify either \code{qk} together with \code{q} to apply one of the pre-defined
#' nominal levels or has to define the argument \code{q_seq} for an own sequence of nominal levels.
#'
#' @references
#'   Xie and Lederer (2021). \emph{Aggregating Knockoffs for False Discovery Rate Control with an Application to Gut Microbiome Data.}
#'   Entropy 23(2), 230.
#'   \url{https://www.mdpi.com/1099-4300/23/2/230/xml}
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
#' res.uKO <- run.uKO(X, y, sets = T)
#' res.uKO
#'
#' # Advanced usage with customized knockoff construction (equi-correlated)
#' equi.knock <- function(X) create.second_order(X, method = "equi")
#' res.uKO <- run.uKO(X, y, knockoffs = equi.knock, sets = T)
#' res.uKO
#'
#' @export
run.uKO <- function(X, y,
                    knockoffs = create.second_order,
                    statistic = stat.glmnet_coefdiff,
                    qk = "decseq", q = 0.2, K = 5, q_seq = NULL,
                    offset = 1, sets = F){

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

  if(q < 0 | q > 1) {
    stop('q must be between 0 and 1')
  }

  if (!is.function(knockoffs)) stop('Input knockoffs must be a function')
  if (!is.function(statistic)) stop('Input statistic must be a function')

  if(!K == round(K)){
    stop("K must be an integer")
    }

  #If own sequence not specified: authors sequences
  if(is.null(q_seq)){

    if(qk == "decseq"){
      k_seq <- 1:K
      q_seq <- q/(2^(k_seq-1))
    }

    else if(qk == "ave"){
      q_seq <- q/rep(K,K)
    }

    else  stop("Unknown type specified for qk")
  }

  #Validate equal length of q_seq and K.
  stopifnot(length(q_seq) == K)

  ##end checks

  #Knockoff construction
  Xk <- multi.knockoffs(X = X, K = K, knockoffs = knockoffs)

  #Knockoff filters
  multi.filter <- multi.knockfilter(X, Xk, y, q = q_seq, offset = offset,
                                  statistic = statistic)


  q_bound <- sum(multi.filter$q)
  Shatk <- multi.filter$Shat.list

  if(K == 1)
  {S_hat.final = Shatk[[1]]
  }else{
    S_hat.final <- sort(Reduce(union, Shatk))  #Take Union
  }


  #Return Part
  if(sets == T){
    return(list(Shat=S_hat.final, K = K, FDRbound = q_bound, sets = Shatk))
  }else{
    return(list(Shat=S_hat.final, K = K, FDRbound = q_bound))
  }
}
