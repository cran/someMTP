\name{someMTP-package}
\alias{someMTP-package}
\alias{someMTP}
\docType{package}
\title{ Some Multiple Testing Procedures
}
\description{
It is a collection of functions for Multiplicty Correction and Multiple Testing.
}
\details{
\tabular{ll}{
Package: \tab someMTP\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-07-06\cr
License: \tab GPL (>= 2) \cr
LazyLoad: \tab yes\cr
}

}
\author{
livio finos

Maintainer: <livio@stat.unipd.it>
}
\references{
For weighted methods:

Benjamini, Hochberg (1997). Multiple hypotheses testing with weights. Scand. J. Statist. 24, 407-418.

Finos,  Salmaso (2007). FDR- and FWE-controlling methods using data-driven weights. Journal of Statistical Planning and Inference, 137,12, 3859-3870.



For LSD test:

J. L\"auter, E. Glimm and S. Kropf (1998). Multivariate test based on Left-Spherically Distributed Linear Scores. 
The Annals of Statistics, Vol. 26, No. 5, 1972-1988

Finos (2010). A note on Left-Spherically Distributed Test with covariates. University of Padua, Dept Statistical Science. Working Paper #8
}
\keyword{ package }
%\seealso{ \code{\link[someKfwer:kfwe-package.Rd]{someKfwer}} }

\examples{
set.seed(13)
y <- matrix(rnorm(5000),5,1000) #create toy data
y[,1:100] <- y[,1:100]+3 #create toy data

p <- apply(y,2,function(y) t.test(y)$p.value) #compute p-values
M2 <- apply(y^2,2,mean) #compute ordering criterion

fdr   <- p.adjust(p,method="BH") #(unweighted) procedure, fdr control
 sum(fdr<.05)
fdr.w <- p.adjust.w(p,method="BH",w=M2) #weighted procedure, weighted fdr control
 sum(fdr.w<.05)
 
fwer   <- p.adjust(p,method="holm") #(unweighted) procedure, fwer control
 sum(fwer<.05)
fwer.w <- p.adjust.w(p,method="BHfwe",w=M2) #weighted procedure, weighted fwer (=fwer) control
sum(fwer.w<.05)

plot(M2,-log10(p))
}
