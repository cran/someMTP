\name{fdrOrd}
\alias{fdrOrd}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Controlling the False Discovery Rate in ordered Test}
\description{
Ordinal procedure controlling the FDR
}
\usage{
fdrOrd(p, q = .01, ord = NULL, GD=FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{ vector of p-values
%%     ~~Describe \code{p} here~~
}
  \item{q}{ average FDR level
%%     ~~Describe \code{alpha} here~~
}
  \item{ord}{ Values on the basis of which the procedure select the hypotheses (following decreasing order). The vector have the same length of \code{p}. If \code{NULL} the natural ordering is considered.
%%     ~~Describe \code{ord} here~~
}
  \item{GD}{ Logic value. Should the correction for general dependence be applied?
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{The function returns an object of class \code{someMTP.object}.
    \item{\code{rej}:}{a logical vector indicating whenever the related hypotesis have been rejected. }
    \item{\code{ord}:}{The vector used to sort the p-values (decrasing).}
    \item{\code{MTP}:}{"fdrOrd"}
    \item{\code{GD}:}{A logical value incating if the correction for General Dependence have been used or not.}
    \item{\code{q}:}{The level of contrelled FDR.}
}

\references{
%%  ~~who you are~~
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
