\name{lsd.object class}
\docType{class}
\alias{lsd.object}
\alias{lsd.object-class}
\alias{show,lsd.object-method}
\alias{p.value,lsd.object-method}
\alias{summary,lsd.object-method}
\alias{summary}
\alias{p.value}


\title{Class "lsd.object" for storing the result of the function lsd}

\description{The class lsd.object is the output of a call to \code{\link{lsd.test}}}

\section{Slots}{
  
  \describe{
    \item{\code{F }:}{the test statistic}
    \item{\code{df }:}{the degrees of freedom of F}
    \item{\code{globalP}:}{the associated p-value}
    \item{\code{D}:}{the matrix used in the test (it provide information on the influence of columns in \code{resp} to the test)}     
    \item{\code{call}:}{The matched call to \code{\link{lsd}}.}
	\item{\code{MTP}:}{The procedure used ("fdrOrd", "kfweOrd" or others).}
  }
}

\section{Methods}{
  \describe{
    \item{p.value}{(gt.object): Extracts the p-values.}
    \item{show}{lsd.object: Prints the test results: p-value, test statistic, expected value of the test statistic under the null hypothesis, standard deviation of the test statistic under the null hypothesis, and number of covariates tested.}
    \item{summary}{lsd.object: Prints the test results: p-value, test statistic, expected value of the test statistic under the null hypothesis, standard deviation of the test statistic under the null hypothesis, and number of covariates tested. }
   }
}

\author{Livio Finos: \email{livio@stat.unipd.it}}

\seealso{\code{\link{lsd}}}

\keyword{methods}


\examples{
    # Simple examples with random data here
    set.seed(1)
	#Standard multivariate LSD test for one sample case
	X=matrix(rnorm(50),5,10)+5
	res <- lsd.test(resp=X,alternative=rep(1,5))
	print(res)
	p.value(res)
    summary(res,showD=TRUE)
} 