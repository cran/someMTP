\name{someMTP.object class}
\docType{class}
\alias{someMTP.object}
\alias{length-method}
\alias{length}
\alias{summary-method}
\alias{someMTP.object-class}
\alias{show,someMTP.object-method}
\alias{summary,someMTP.object-method}
\alias{length,someMTP.object-method}
\alias{names,someMTP.object-method}
\alias{names<-,someMTP.object-method}
\alias{sort,someMTP.object-method}


\title{Class "someMTP.object" for storing the result of the function fdrOrd}

\description{The class someMTP.object is the output of a call to
\code{\link{fdrOrd}}. It also stores the information needed for related plots. }

\section{Slots}{
  \describe{
    \item{\code{rej}:}{a logical vector indicating whenever the related hypotesis have been rejected. }
	\item{\code{p}:}{The vector of (raw) p-values used in the procedure.}
    \item{\code{ord}:}{The vector used to sort the p-values (decreasing).}
	\item{\code{idOrd}:}{The vector of indices used in sorting.}
    \item{\code{MTP}:}{The type of procedure used.}
    \item{\code{GD}:}{A logical value incating if the correction for General Dependence have been used or not.}
    \item{\code{q}:}{The level of contrelled FDR when MTP=="fdrOrd".}
    \item{\code{k}:}{The number of false rejection when MTP=="kfweOrd"}
	\item{\code{J}:}{The number of allowed Jumps when MTP=="kfweOrd"}
    \item{\code{alpha}:}{The significance level when MTP=="kfweOrd"}
    \item{\code{alphaprime}:}{The significance level of individual tests.}
  }
}

\section{Methods}{
  \describe{
    \item{show}{someMTP.object: Prints the test results.}
    \item{summary}{someMTP.object: Prints the test results (as \code{show}).}
	\item{draw}{someMTP.object: Plots results; \code{what = c("all","ordVsP", "stepVsR")}  }
    \item{sort}{\code{signature(x = "someMTP.object")}: Sorts the p-values to decreasing order of \code{ord}.}
    \item{length}{\code{signature(x = "someMTP.object")}: The number of tests performed.}
    \item{names}{\code{signature(x = "someMTP.object")}:  Extracts the row names of the results matrix.}
    \item{names<-}{\code{signature(x = "someMTP.object")}: Changes the row names of the results matrix. Duplicate names are not allowed, but see \code{alias}.}
    }
}

\author{Livio Finos: \email{livio@stat.unipd.it}}

\seealso{\code{\link{someMTP.object}}}

\keyword{methods}


\examples{
    # Simple examples with random data
    set.seed(17)
	x=matrix(rnorm(60),3,20)
	x[,1:10]=x[,1:10]+2 ##variables 1:10 have tests under H1
	ts=apply(x,2,function(x) t.test(x)$statistic)
	ps=apply(x,2,function(x) t.test(x)$p.value)
	m2=apply(x^2,2,mean)
	pOrd <- fdrOrd(ps,q=.05,ord=m2)
	pOrd
    length(pOrd)
	names(pOrd) <- paste("V",1:20,sep="")
	names(pOrd)
	
} 