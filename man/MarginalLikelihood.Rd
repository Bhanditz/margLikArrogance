\name{MarginalLikelihood}
\alias{MarginalLikelihood}
\title{Estimate the marginal likelihood of a model}

\description{This is the main function of the \code{margLikArrogance}
  package.  It estimates the marginal likelihood of a Bayesian model
  from samples of the model's posterior parameter distribution.
}

\usage{
MarginalLikelihood(theta, ll, cov.prob=0.5,
                   bounds=matrix(c(-Inf, Inf), ncol=ncol(theta), nrow=2))
}

\arguments{
  \item{theta}{Samples from the posterior parameter distribution}
  \item{ll}{Log-likelihood of data and samples from theta}
  \item{cov.prob}{Coverage probability of histogram}
  \item{bounds}{Matrix specifying cube where log-likelihood positive}
}

\details{

  This function implements an ''arrogance sampling'' method.  For
  details, see the files \code{margLikArrogance.pdf} and
  \code{vignette.pdf} (or the corresponding Sweave files) in this
  module's installation directory.  They explain the method and provide
  several examples of the module's use.

  The main arguments are \code{theta} and \code{ll}.  Suppose you have
  just observed some data \eqn{x} and are evaluating model \eqn{T}.
  Then \code{theta} contains samples \eqn{\theta_i} from the posterior
  parameter distribution \eqn{p(\theta_i|x, T)}.  \code{ll} contains the
  value \eqn{log(p(x|\theta_i, T)p(\theta_i|T))} for each
  \eqn{\theta_i}.

  The size of the histogram is selected so that \code{cov.prob} of the
  points lie within it.  By default, this is set to 0.5, so that the
  histogram will be postive for half the samples \eqn{\theta}.

  The constructed histogram must lie entirely in a region where
  \eqn{p(\theta|x, T)} is positive.  If it reaches 0, the resulting
  estimate may be wrong or have infinite variance.  The \code{bounds}
  matrix can be used to check for this.  The \code{bounds} matrix should
  have two rows, be of the same dimensionality as \code{theta}, and
  specify the corners of a hypercube where \eqn{p(\theta|x, T)} is known
  to be positive.  If the histogram falls outside the cube, an error
  will be raised.  The user can also check the constructed histogram
  manually.
  
}

\value{This function returns a list with these elements:
  \describe{
    \item{mll}{The actual estimate of marginal likelihood}
    \item{width}{The selected width of the histogram}
    \item{conf.interval}{The 95\% confidence interval of the estimate}
    \item{hist}{The histogram built from the samples}
  }
}

\examples{
# Find the normalizing constant for a two-dimensional normal distribution
theta <- cbind(rnorm(1000), rnorm(1000))
ll <- 12345 + dnorm(theta[, 1], log=TRUE) + dnorm(theta[, 2], log=TRUE)
ml.out <- MarginalLikelihood(theta, ll)
ml.out$mll
ml.out$conf.interval

# See the vignette for more examples
}