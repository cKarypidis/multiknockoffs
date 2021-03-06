% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiknockfilter.R
\name{multi.knockfilter}
\alias{multi.knockfilter}
\title{Multiple knockoff filter}
\usage{
multi.knockfilter(
  X,
  Xk,
  y,
  q = 0.2,
  offset = 1,
  statistic = stat.glmnet_coefdiff
)
}
\arguments{
\item{X}{n x p matrix or data frame of original variables.}

\item{Xk}{list with K elements containing the n x p knockoff matrices.}

\item{y}{response vector of length n.}

\item{q}{either a scalar or vector of nominal levels. If a scalar is supplied,
then the same nominal level is used for each knockoff run. Default: 0.2.}

\item{offset}{either 0 (knockoff) or 1 (knockoff+). Default: 1.}

\item{statistic}{function that computes the score vector \eqn{W} of length p. It must take the data matrix,
knockoff matrix and response vector as input and outputs a vector of computed
scores. Either choose one score statistic from the \code{knockoff} package or
define it manually. Default: \code{stat.glmnet_coefdiff} (see below).}
}
\value{
A list containing following components:
\item{W.list}{the K score vectors of each knockoff run.}
\item{Shat.list}{the K selection sets of each knockoff run.}
\item{q}{the nominal level of each knockoff run.}
}
\description{
This function estimates multiple knockoff runs with different knockoff matrices.
}
\details{
This function requires the installation of the \code{knockoff} package prior to its execution.

The default score function \code{stat.glmnet_coefdiff} is from the \code{knockoff} package.
It fits a  Lasso regression where the regularization parameter \eqn{\lambda} is tuned by cross-validation.
Then, the score is computed as the difference between
\deqn{W_j = |Z_j| - |\tilde{Z}_j|}
where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the coefficient estimates for the
jth variable and its knockoff, respectively.

The function should be used in combination with \code{\link{multi.knockoffs}} (see example).
}
\examples{
n <- 400; p <- 200; s_0 <- 30
amplitude <- 1; mu <- rep(0,p); rho <- 0.25
Sigma = toeplitz(rho^(0:(p-1)))

X <- MASS::mvrnorm(n, mu, Sigma)
nonzero <- sample(p, s_0)
beta <- amplitude * (1:p \%in\% nonzero)
y <- X \%*\% beta + rnorm(n)

# Construction of K knockoff matrices
Xk <- multi.knockoffs(X, K = 5)

# Basic usage with default arguments
multi.res <- multi.knockfilter(X, Xk, y)


}
\references{
Candes, Fan, Janson, and Lv (2018). \emph{Panning for gold. model-X knockoffs for high
dimensional controlled variable selection}. Journal of the Royal Statistical Society:
Series B (Statistical Methodology) 80(3), 551-577.
}
