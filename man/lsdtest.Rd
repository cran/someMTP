\name{lsd.test}
\alias{lsd.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
Multivariate Left Spherically Distributed (LSD) linear scores test.
}
\description{
It performs the multivariate Left Spherically Distributed linear scores test of L\"auter et al. (The Annals of Statistics, 1998) (see also details below).
}

\usage{
lsd.test(resp, alternative = rep(1, dim(as.matrix(resp))[2]), null = NULL, D = function(resp) t(as.matrix(diag(t(resp) \%*\% resp))))
}

\arguments{
  \item{resp}{
  is a n x p matrix of spherical responses.
}
  \item{alternative}{
  is the n x k matrix of predictors to be tested
}
  \item{null}{
  is the n x h matrix of predictors not to be tested (i.e. covariates) 
  }
  \item{D}{
  is a function with argument \code{resp} prducing the q x p transformation matrix or is the trasformation matrix itself.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{F }{is the F statistic}
  \item{p }{is the associated p-value}
}
\references{

J. L\"auter, E. Glimm and S. Kropf (1998) Multivariate test based on Left-Spherically Distributed Linear Scores. 
The Annals of Statistics, Vol. 26, No. 5, 1972-1988

Finos (2010) A note on Left-Spherically Distributed Test with covariates. University of Padua, Dept Statistical Science. Working Paper #8
}
\author{
Livio Finos and Dario Basso
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
set.seed(1)
#Standard multivariate LSD test for one sample case
X=matrix(rnorm(50),5,10)+5
lsd.test(resp=X,alternative=rep(1,5))

#Standard multivariate LSD test for two sample case
X2=X+matrix(c(0,0,1,1,1),5,10)*10
lsd.test(resp=X2,null=rep(1,5),alternative=c(0,0,1,1,1))

#General multivariate LSD test for linear predictor with covariates
lsd.test(resp=X2,null=cbind(rep(1,5),c(0,0,1,1,1)),alternative=1:5)
}
\keyword{htest}
