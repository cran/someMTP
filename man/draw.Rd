\name{draw}
\alias{draw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots results of fdrOrd()
}
\description{
Plots results of fdrOrd()

}
\usage{
draw(object, what = c("all", "ordVsP", "stepVsR"), pdfName = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ a \code{someMTP.object} resulting from fdrOrd()
}
  \item{what}{
  what to plot; \code{"all"} is the default
}
  \item{pdfName}{ it is the pdf filename where the plot will be saved. If \code{pdfNane} is null (the default) the plot will show as window.
}
}
\value{
No value is returned
}
\author{
Livio Finos
}
\seealso{
See Also \code{\link{fdrOrd}}.
}
\examples{
set.seed(17)
	x=matrix(rnorm(60),3,20)
	x[,1:10]=x[,1:10]+2 ##variables 1:10 have tests under H1
	ts=apply(x,2,function(x) t.test(x)$statistic)
	ps=apply(x,2,function(x) t.test(x)$p.value)
	m2=apply(x^2,2,mean)
	pOrd <- fdrOrd(ps,q=.05,ord=m2)	
	draw(pOrd)
}
