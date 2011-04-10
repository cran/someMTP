\name{fdrOrd}
\alias{fdrOrd}
\alias{kfweOrd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Controlling the False Discovery Rate in ordered Test}
\description{
Ordinal procedure controlling the FDR and the Generalized FWER
}
\usage{
fdrOrd(p, q = .01, ord = NULL, GD=FALSE)
kfweOrd(p, k = 1, alpha = 0.01, ord = NULL, alpha.prime = alpha, J = qnbinom(alpha, k, alpha.prime), GD = FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{ vector of p-values
%%     ~~Describe \code{p} here~~
}
  \item{ord}{ Values on the basis of which the procedure select the hypotheses (following decreasing order). The vector have the same length of \code{p}. If \code{NULL} the natural ordering is considered.
%%     ~~Describe \code{ord} here~~
}
  \item{q}{ average FDR level
%%     ~~Describe \code{alpha} here~~
}
  \item{alpha}{ global significance level
%%     ~~Describe \code{alpha} here~~
}
  \item{k}{ number of allowed errors in kFWE controls
%%     ~~Describe \code{k} here~~
}
  \item{J}{ number of allowed jumps befor stopping
%%     ~~Describe \code{J} here~~
}
  \item{alpha.prime}{ univariate alpha for single step Guo and Romano procedure
%%     ~~Describe \code{alpha} here~~
}
  \item{GD}{ Logic value. Should the correction for general dependence be applied?
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{The function returns an object of class \code{someMTP.object}.
    \item{\code{rej}:}{a logical vector indicating whenever the related hypotesis have been rejected. }
	\item{\code{p}:}{the vector of p-values used in the call}
    \item{\code{ord}:}{The vector used to sort the p-values (decrasing).}
    \item{\code{MTP}:}{"fdrOrd" or "kfweOrd"}
    \item{\code{GD}:}{A logical value incating if the correction for General Dependence have been used or not.}
    \item{\code{q}:}{The level of controlled FDR.}
	\item{\code{alpha}:}{The level of controlled k-FWER}
	\item{\code{alphaprime}:}{The significance level of individual tests}
	\item{\code{k}:}{Number of allowed Errors}
	\item{\code{J}:}{Number of allowed Jumps}
}

\references{
L. Finos, A. Farcomeni (2010). k-FWER Control without p-value Adjustment, with Application to Detection of Genetic Determinants of Multiple Sclerosis in Italian Twins. Biometrics.

}
\author{ L. Finos and A. Farcomeni
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{draw}}
}
\examples{


	set.seed(17)
	x=matrix(rnorm(60),3,20)
	x[,1:10]=x[,1:10]+2 ##variables 1:10 have tests under H1
	ts=apply(x,2,function(x) t.test(x)$statistic)
	ps=apply(x,2,function(x) t.test(x)$p.value) #compute p-values
	m2=apply(x^2,2,mean)           #compute ordering criterion

	pOrd <- fdrOrd(ps,q=.05,ord=m2)   #ordinal Procedure
	pOrd
	draw(pOrd)	
	sum(p.adjust(ps,method="BH")<=.05)  #rejections with BH
	
	kOrd <- kfweOrd(ps,k=5,ord=m2)                          #ordinal procedure
	kOrd
	kOrdGD <- kfweOrd(ps,k=5,ord=m2,GD=TRUE)                  #ordinal procedure (any dependence)
	kOrdGD


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
