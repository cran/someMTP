\name{p.adjust.w}
\alias{p.adjust.w}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adjust P-values for Multiple Comparisons
}
\description{
Given a set of p-values, returns p-values adjusted using one of several (weighted) methods.
It extends the method of \code{p.adjust\{stats\}}
}
\usage{
p.adjust.w(p, method = c("bonferroni","holm","BHfwe","BH","BY"), n = length(p),w=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
vector of p-values (possibly with NAs)
}
  \item{method}{
correction method
}
  \item{n}{
	number of comparisons, must be at least length(p); only set this (to non-default) when you know what you are doing!
}
  \item{w}{
  weigths to be used. \code{p.adjust.w(..., rep(1,length(p)))} produces the same results as in \code{p.adjust(...)} (i.e. the unweighted counterpart).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector of corrected p-values (same length as p) having two attributes: \code{attributes(...)$w} is the vecotr of used weights and \code{attributes(...)$method} is the method used.
}
\references{

Benjamini, Hochberg (1997). Multiple hypotheses testing with weights. Scand. J. Statist. 24, 407-418.

Finos,  Salmaso (2007). FDR- and FWE-controlling methods using data-driven weights. Journal of Statistical Planning and Inference, 137,12, 3859-3870.
}
\author{
Livio Finos
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\seealso{
 \code{\link{p.adjust}}
}
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
\keyword{htest}
