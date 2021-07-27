#' Construction multiple knockoff matrices
#'
#' This function generates a specified number of knockoff matrices given the data
#' and a knockoff sampling function.
#'
#' @param X n x p matrix of original variables.
#' @param K number of knockoff matrices.
#' @param knockoffs function for the knockoff construction. It must take the n x p matrix as input
#'                  and it must return a n x p knockoff matrix. Either choose a knockoff sampler of
#'                  the \code{knockoff} package or define it manually. Default: \code{create.second_order} (see below).
#'
#' @return A list that contains the K knockoff matrices.
#'
#' @details
#'
#' This function requires the installation of the \code{knockoff} package prior to its execution.
#'
#' The default knockoff sampler \code{create.second_order} is the second-order Gaussian knockoff construction from
#' the \code{knockoff} package.
#'
#' @references
#'   Candes, Fan, Janson, and Lv (2018). \emph{Panning for gold. model-X knockoffs for high
#'   dimensional controlled variable selection}. Journal of the Royal Statistical Society:
#'   Series B (Statistical Methodology) 80(3), 551-577.
#'
#' @examples
#' n <- 400; p <- 200
#' mu <- rep(0,p); rho <- 0.25
#' Sigma <- toeplitz(rho^(0:(p-1)))
#' X <- MASS::mvrnorm(n, mu, Sigma)
#'
#' # Basic usage with default arguments
#' Xk <- mult.knockoffs(X, K = 5)
#' Xk
#'
#' # Advanced usage with customized knockoff construction (equi-correlated)
#' equi.knock <- function(X) create.second_order(X, method = "equi")
#' Xk <- mult.knockoffs(X, K = 5, knockoffs = equi.knock)
#'
#'
#' @export
mult.knockoffs <- function(X, K, knockoffs = create.second_order){

  library(knockoff)

  #Input checks
  if (!is.matrix(X)){
    stop("Input X must be a matrix")
  }

  if (!is.function(knockoffs)) stop('Input knockoffs must be a function')

  if(!K == round(K)){
    stop("K must be an integer")
  }

  #End checks

  # Validate input dimensions
  knock.list <- list()

  for(k in 1:K){

  # Create knockoff variables
  knock_variables = knockoffs(X)

  #Check for correct function outout
  if (is(knock_variables,"matrix")){
    Xk = knock_variables
    rm(knock_variables)
  } else {
    stop('Knockoff variables of incorrect type')
  }
  knock.list[[k]] <- Xk
  }
  names(knock.list) <- paste("Xk", sep = "", 1:K)
  return(knock.list)
}


#' Multiple knockoff filter
#'
#' This function estimates multiple knockoff runs with different knockoff matrices.
#'
#' @param X n x p matrix of original variables.
#' @param Xk list with K elements containing the n x p knockoff matrices.
#' @param y response vector of length n.
#' @param q nominal level for the FDR control. Default: 0.2.
#' @param offset either 0 (knockoff) or 1 (knockoff+). Default: 1.
#' @param statistic function that computes the score vector \eqn{W} of length p. It must take the data matrix,
#'                  knockoff matrix and response vector as input and outputs a vector of computed
#'                  scores. Either choose one score statistic from the \code{knockoff} package or
#'                  define it manually. Default: \code{stat.glmnet_coefdiff} (see below).
#'
#' @return A list containing following components:
#'  \item{W.list}{the K score vectors of each knockoff run.}
#'  \item{Shat.list}{the K selection sets of each knockoff run.}
#'  \item{q}{the nominal level of each knockoff run.}
#'
#' @details
#'
#' This function requires the installation of the \code{knockoff} package prior to its execution.
#'
#' The default score function \code{stat.glmnet_coefdiff} is from the \code{knockoff} package.
#' It fits a  Lasso regression where the regularization parameter \eqn{\lambda} is tuned by cross-validation.
#' Then, the score is computed as the difference between
#'   \deqn{W_j = |Z_j| - |\tilde{Z}_j|}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the coefficient estimates for the
#' jth variable and its knockoff, respectively.
#'
#' The function should be used in combination with \code{\link{mult.knockoffs}} (see example).
#'
#' @references
#'   Candes, Fan, Janson, and Lv (2018). \emph{Panning for gold. model-X knockoffs for high
#'   dimensional controlled variable selection}. Journal of the Royal Statistical Society:
#'   Series B (Statistical Methodology) 80(3), 551-577.
#'
#' @examples
#' n <- 400; p <- 200; s_0 <- 30
#' amplitude <- 1; mu <- rep(0,p); rho <- 0.25
#' Sigma = toeplitz(rho^(0:(p-1)))
#'
#' X <- MASS::mvrnorm(n, mu, Sigma)
#' nonzero <- sample(p, s_0)
#' beta <- amplitude * (1:p %in% nonzero)
#' y <- X %*% beta + rnorm(n)
#'
#' # Construction of K knockoff matrices
#' Xk <- mult.knockoffs(X, K = 5)
#'
#' # Basic usage with default arguments
#' mult.res <- mult.knockfilter(X, Xk, y)
#'
#'
#' @export
mult.knockfilter <- function(X, Xk, y, q = 0.2, offset = 1, statistic = stat.glmnet_coefdiff){

  library(knockoff)

  # Input checks.
  if (is.data.frame(X)) {
    X.names = names(X)
    X = as.matrix(X, rownames.force = F)
  } else if (is.matrix(X)) {
    X.names = colnames(X)
  } else {
    stop('Input X must be a numeric matrix or data frame')
  }
  if (!is.numeric(X)) stop('Input X must be a numeric matrix or data frame')

  #Checks
  if (!is.list(Xk)) stop('Input Xk must be a list')
  if (!is.function(statistic)) stop('Input statistic must be a function')

  if(offset!=1 && offset!=0) {
    stop('Input offset must be either 0 or 1')
  }

  K <- length(Xk)
  if(length(q) == 1){q <- rep(q, K)}
  if (length(q) != K) stop('Vector q must be of length K')


  #Knockoff filter for 1 to K
  W.list <- list()
  Shat.list <- list()

  for(l in 1:K){
    W.list[[l]] <- statistic(X, Xk[[l]], y)
    thres <- knockoff.threshold(W.list[[l]], fdr = q[l], offset = offset)
    Shat.k <- which(W.list[[l]] >= thres)
    if (!is.null(X.names)){names(Shat.k) = X.names[Shat.k]}
    Shat.list[[l]] <-  Shat.k
  }
  names(W.list) <- paste("W", sep = "", 1:K)
  names(Shat.list) <- paste("S", sep = "", 1:K)


  res <- list(W.list = W.list, Shat.list = Shat.list, q = q)
  return(res)
}
