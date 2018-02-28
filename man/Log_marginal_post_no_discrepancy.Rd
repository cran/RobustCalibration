\name{Log_marginal_post_no_discrepancy}
\alias{Log_marginal_post_no_discrepancy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Natural Logorithm of the posterior with no discrepancy function. }
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function compute the natural Logorithm of the posterior assuming no discrepancy function.
}
\usage{
Log_marginal_post_no_discrepancy(param,  output,  p_theta,   X, have_mean, 
inv_output_weights,  cm_obs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param}{
%%     ~~Describe \code{param} here~~
Current parameters in the MCMC. 
}
  \item{output}{
Experimental observations.
}

  \item{p_theta}{
Number of calibration parameters.
}

  \item{X}{
Number of mean discrepancy parameters.
}
  \item{have_mean}{
Whether the mean discrepancy is zero or not.
}

  \item{inv_output_weights}{
Inverse of the weights of the outputs}

  \item{cm_obs}{
Outputs from the mathematical model.
}



}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Natural logorithm of the posterior assuming no discrepancy function. }
\references{
A. O'Hagan and M. C. Kennedy (2001), \emph{Bayesian calibration of computer models}, \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology}, \bold{63}, 425-464.

Mengyang Gu. (2016). Robust Uncertainty Quantification and Scalable Computation for Computer Models with Massive Output. Ph.D. thesis. Duke University.

M. Gu and L. Wang (2017) \emph{An improved approach to Bayesian computer model calibration and prediction}. arXiv preprint arXiv:1707.08215.

}
\author{
%%  ~~who you are~~
\packageAuthor{RobustCalibration}

Maintainer: \packageMaintainer{RobustCalibration}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
% \code{\link{log_marginal_lik}}.
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")

\keyword{internal}