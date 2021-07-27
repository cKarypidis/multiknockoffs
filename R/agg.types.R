#' Union of selection sets
#'
#' This function aggregates a list of selection sets by their union.
#'
#' @param Shat.list list of K elements containing the selection sets.
#'
#' @return Aggregated union set of selected variables.
#'
#' @details
#'
#' The function can be used in combination with \code{\link{mult.knockoffs}} and
#' \code{\link{mult.knockfilter}}.
#'
#' The function does not have to be used in the context of multiple knockoffs.
#' It can also be used to aggregate the selection sets of multiple FDR controlling procedures
#' in general.
#'
#'
#' @examples
#' #Multiple knockoff example (union knockoffs)
#' n <- 400; p <- 200; s_0 <- 30
#' amplitude <- 1; mu <- rep(0,p); rho <- 0.25
#' Sigma <- toeplitz(rho^(0:(p-1)))
#' X <- MASS::mvrnorm(n, mu, Sigma)
#' nonzero <- sample(p, s_0)
#' beta <- amplitude * (1:p %in% nonzero)
#' y <- X %*% beta + rnorm(n)
#'
#' K <- 5
#' Xk <- mult.knockoffs(X, K = K)
#' q <- 0.2/2^((1:K)-1)
#' mult.res <- mult.knockfilter(X, Xk, y, q = q)
#' agg.union(mult.res$Shat.list)
#'
#'
#' #General example (selection sets with indices between 1 and 30)
#' Shat.list <- list(s1 = c(2,4,3,1,20,30), s2 = c(3,30,23,1,4,8),
#'                  s3 = c(3,4,5,13,15,12, 23:29, 30), s4 = c(1:10, 13:15, 17),
#'                  s5 = c(15:20, 23))
#' agg.union(Shat.list)
#'
#'
#' @export
agg.union <- function(Shat.list){

  if (!is.list(Shat.list)) stop('Input S.hat must be a list')

  K <- length(Shat.list)

  if(K == 1){
    S_hat.final = Shat.list[[1]]
  } else{
    S_hat.final <- sort(Reduce(union, Shat.list))  #Take Union
  }
  return(S_hat.final)
}

#' @title Aggregation by p-value knockoffs
#'
#' @description
#' This function only runs the aggregation. Given a list containing B elements
#' with the score vectors of length p, it computes aggregated p-values and applies
#' either BH or BY in the last step.
#'
#' @param W.list list with B elements containing the vectors of scores with length p each.
#' @param q nominal level for the FDR control. Default: 0.2.
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
#' This function should be used in combination with \code{\link{mult.knockoffs}} and
#' \code{\link{mult.knockfilter}} (see example).
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
#' # Construction of K knockoff matrices
#' equi.knock <- function(X) create.second_order(X, method = "equi")
#' Xk <- mult.knockoffs(X, K = 20, knockoffs = equi.knock)
#'
#' #Multiple knockoff filter
#' mult.res <- mult.knockfilter(X, Xk, y)
#'
#' pKO.res <- agg.pKO(mult.res$W.list)
#' pKO.res
#'
#' @export
agg.pKO <-  function(W.list, q = 0.2, gamma = 0.3, offset = 1, method = "BH", pvals = F){

  #Input checks
  if (!is.list(W.list)) stop('Input W.list must be a list')

  if(offset!=1 && offset!=0) {
    stop('Input offset must be either 0 or 1')
  }

  # Input dimensions
  p <- length(W.list[[1]])
  B <- length(W.list)

  #Store of all B x p p-values
  int.pval <- matrix(0, nrow = B, ncol = p)

  #Step1: Compute Intermediate p-values
  for(b in 1:B){#Loop over all b

    W <- W.list[[b]]

    #Intermediate p-values
    for(j in 1:p){
      if(W[j] <= 0){
        int.pval[b,j] <- 1
      }
      else{
        int.pval[b,j] <- (1/p)* (offset + sum(W <= -W[j]))
      }
    }
  }

  #Step2: Aggregated p-values
  agg.pval <- quantile.aggregation(int.pval, B, gamma)

  #Step3: Compute BH/BY
  if(method %in% c("BH", "BY")){
    p_corr <-  p.adjust(agg.pval,method = method)
  } else  stop("Unknown method specified")


  #Step 4: Selection Set
  Shat <- as.numeric(which(p_corr <= q))
  res <- list(Shat = Shat, B = B)

  if(pvals == T){res <- list(Shat=Shat, B = B, pvals = agg.pval)}

  return(res)
}


#' Aggregation by ADAGES
#'
#' This function aggregates a list of selection sets by ADAGES.
#'
#' @param Shat.list list of K elements containing the selection sets.
#' @param p number of original variables of the model.
#'
#' @return A list containing following components:
#'  \item{Shat}{aggregated selection set.}
#'  \item{c}{optimal threshold.}
#'  \item{K}{number of aggregated sets.}
#'
#' @details
#'
#' The function can be used in combination with \code{\link{mult.knockoffs}} and
#' \code{\link{mult.knockfilter}}.
#'
#' The function does not have to be used in the context of multiple knockoffs.
#' It can also be used to aggregate the selection sets of any FDR controlling procedure
#' in general.
#'
#' Applies the minimization of the complexity ratio as a criterion to determine the optimal threshold.
#'
#' @references
#'   Gui (2020). \emph{ADAGES: adaptive aggregation with stability for distributed feature selection}.
#'   Proceedings of the 2020 ACM-IMS on Foundations of Data Science Conference.
#'   \url{https://arxiv.org/pdf/2007.10776.pdf}
#'
#' @examples
#' #Multiple knockoff example
#' n <- 400; p <- 200; s_0 <- 30
#' amplitude <- 1; mu <- rep(0,p); rho <- 0.25
#' Sigma <- toeplitz(rho^(0:(p-1)))
#' X <- MASS::mvrnorm(n, mu, Sigma)
#' nonzero <- sample(p, s_0)
#' beta <- amplitude * (1:p %in% nonzero)
#' y <- X %*% beta + rnorm(n)
#'
#' Xk <- mult.knockoffs(X, K = 5)
#' mult.res <- mult.knockfilter(X, Xk, y)
#' agg.ADAGES(mult.res$Shat.list, p = p)
#'
#'
#' #General example (selection sets with indices between 1 and 30)
#' Shat.list <- list(s1 = c(2,4,3,1,20,30), s2 = c(3,30,23,1,4,8),
#'                  s3 = c(3,4,5,13,15,12, 23:29, 30), s4 = c(1:10, 13:15, 17),
#'                  s5 = c(15:20, 23))
#' agg.ADAGES(Shat.list, p = 30)
#'
#'
#' @export
agg.ADAGES <- function(Shat.list, p){

  #Input check
  if (!is.list(Shat.list)) stop('Input Shat.list must be a list')

  K <- length(Shat.list)

  #Step1: Compute mj
  m <- sapply(seq_len(p),
              function(x) {
                sum(x == unlist(Shat.list))
              })

  #Step2: Compute complexity ratio eta
  S_c_card = numeric(K)
  eta = numeric(K-1)

  for (c in 1:K){
    S_c <- which(m >= c)
    S_c_card[c] <- length(S_c)
    if(c > 1){
      eta[c-1] <- (S_c_card[c-1]+1)/(S_c_card[c]+1) #surrogate
    }
  }

  #Step3: Compute sbar
  S_card <- lengths(Shat.list)  #Cardinality of each Set/length of all sets
  s_bar <- mean(S_card)  #s bar

  #Step4: Compute sequence of c: |S_c| > s_bar
  c_seq <- which(S_c_card >= s_bar)

  #Step5: Find c0 and compute c*
  ll <- min(K-1, length(c_seq))

  if(K <= 2){c_op <- 2}else{
    c_op <-  max(c_seq[eta[1:ll] == min(eta[1:ll])])
    #this notation makes sense if we have two equal etas
  }

  #c0 <- max(which(S_c_card >= s_bar))
  #which.min(eta[1:c0])

  #Step6: Aggrgated selection set
  Shat <- which(m >= c_op)

  return(list(Shat = Shat, c = c_op, K = K))
}


#' Aggregation by modified ADAGES
#'
#' This function aggregates a list of selection sets by modified ADAGES.
#'
#' @param Shat.list list of K elements containing the selection sets.
#' @param p number of original variables of the model.
#'
#' @return A list containing following components:
#'  \item{Shat}{aggregated selection set.}
#'  \item{c}{optimal threshold.}
#'  \item{K}{number of aggregated sets.}
#'
#' @details
#'
#' The function can be used in combination with \code{\link{mult.knockoffs}} and
#' \code{\link{mult.knockfilter}}.
#'
#' The function does not have to be used in the context of multiple knockoffs.
#' It can also be used to aggregate the selection sets of any FDR controlling procedure
#' in general.
#'
#' Applies the minimization of the trade-off between the threshold and the model complexity \eqn{c |S|}.
#'
#' @references
#'   Gui (2020). \emph{ADAGES: adaptive aggregation with stability for distributed feature selection}.
#'   Proceedings of the 2020 ACM-IMS on Foundations of Data Science Conference.
#'   \url{https://arxiv.org/pdf/2007.10776.pdf}
#'
#' @examples
#' #Multiple knockoff example
#' n <- 400; p <- 200; s_0 <- 30
#' amplitude <- 1; mu <- rep(0,p); rho <- 0.25
#' Sigma <- toeplitz(rho^(0:(p-1)))
#' X <- MASS::mvrnorm(n, mu, Sigma)
#' nonzero <- sample(p, s_0)
#' beta <- amplitude * (1:p %in% nonzero)
#' y <- X %*% beta + rnorm(n)
#'
#' Xk <- mult.knockoffs(X, K = 5)
#' mult.res <- mult.knockfilter(X, Xk, y)
#' agg.ADAGES.mod(mult.res$Shat.list, p = p)
#'
#'
#' #General example (selection sets with indices between 1 and 30)
#' Shat.list <- list(s1 = c(2,4,3,1,20,30), s2 = c(3,30,23,1,4,8),
#'                  s3=c(3,4,5,13,15,12, 23:29, 30), s4 = c(1:10, 13:15, 17),
#'                  s5 = c(15:20, 23))
#' agg.ADAGES.mod(Shat.list, p = 30)
#'
#'
#' @export
agg.ADAGES.mod <- function(Shat.list, p){

  #Input check
  if (!is.list(Shat.list)) stop('Input Shat.list must be a list')

  K <- length(Shat.list)

  #Step1: Compute mj
  m <- sapply(seq_len(p),
              function(x) {
                sum(x == unlist(Shat.list))
              })

  #Compute Cardinalities
  S_c_card = numeric(K)
  for(c in 1:K){
    S_c <- which(m >= c)
    S_c_card[c] <- length(S_c)
  }

  #Compute c0
  S_card <- lengths(Shat.list)
  s_bar <- mean(S_card)
  #Compute sequence of c: |S_c| > s_bar
  c_seq <- which(S_c_card >= s_bar)
  c0 <- max(c_seq)

  #Compute criterion and selection set
  obj <- S_c_card[1:c0] * (1:c0)
  c_op <- (1:c0)[obj == min(obj)]
  Shat <- which(m >= c_op)
  return(list(Shat = Shat, c = c_op, K = K))
}
