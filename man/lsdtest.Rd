\name{lsd.test}
\alias{lsd.test}
\alias{lsd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
Multivariate Left Spherically Distributed (LSD) linear scores test.
}
\description{
It performs the multivariate Left Spherically Distributed linear scores test of L\"auter et al. (The Annals of Statistics, 1998) (see also details below).
}

\usage{
lsd.test(resp, alternative = 1, null = NULL, D = NULL, data=NULL)
}

\arguments{
 \item{resp}{The response vector of the regression model. May be
      supplied as a vector or as a \code{\link[stats:formula]{formula}}
      object. In the latter case, the right hand side of \code{Y}
      is passed on to \code{alternative} if that argument is missing, or
      otherwise to \code{null}.}

    \item{alternative}{The part of the design matrix corresponding to
      the alternative hypothesis. The covariates of the null model do
      not have to be supplied again here. May be given as a half
      \code{\link[stats:formula]{formula}} object (e.g. \code{~a+b}). In
      that case the intercept is always suppressed.}

    \item{null}{The part of the design matrix corresponding to the null hypothesis. May be given as a design matrix or as a half \code{\link[stats:formula]{formula}} object (e.g. \code{~a+b}). The default for \code{Z} is \code{~1}, i.e. only an intercept. This intercept may be suppressed, if desired, with \code{Z = ~0}.}
    
	\item{data}{Only used when \code{Y}, \code{X}, or \code{Z} is given in formula form. An optional data frame, list or environment containing the variables used in the formulae. If the variables in a formula are not found in \code{data}, the variables are taken from environment(formula), typically the environment from which \code{gt} is called.}
	
  \item{D}{
  is q x p matrix or it is a function with arguments \code{resp} and \code{null}
  returning the q x p transformation matrix.
  When \code{D = NULL}, then  \code{D =  diag(t(resp)\%*\%IP0\%*\%resp)} with \code{IP0 = diag(n) - null\%*\%solve(t(null)\%*\%null)\%*\%t(null)}
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{The function returns an object of class \code{lsd.object}.
  \item{F }{the test statistic}
  \item{df }{the degrees of freedom of F}
  \item{p }{the associated p-value}
  \item{D }{the matrix used in the test (it provide information on the influence of columns in \code{resp} to the test)}
  \item{\code{call}:}{The matched call to \code{\link{lsd.test}}.}
}

\references{

J. Laeuter, E. Glimm and S. Kropf (1998) Multivariate test based on Left-Spherically Distributed Linear Scores. 
The Annals of Statistics, Vol. 26, No. 5, 1972-1988

L. Finos (2011). A note on Left-Spherically Distributed Test with
covariates, Statistics and Probabilty Letters, Volume 81, Issue 6, June 2011, Pages 639-641 
}

\author{
Livio Finos
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
X=matrix(rnorm(50),5,10)+2
lsd.test(resp=X,alternative=~1)

#Standard multivariate LSD test for two sample case
X2=X+matrix(c(0,0,1,1,1),5,10)*10
lsd.test(resp=X2,null=~1,alternative=c(0,0,1,1,1))

#General multivariate LSD test for linear predictor with covariates
lsd.test(resp=X2,null=cbind(rep(1,5),c(0,0,1,1,1)),alternative=1:5)

}
\keyword{htest}
