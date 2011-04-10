\name{p.adjust.w}
\alias{p.adjust.w}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
p.adjust.w(p, method = c("holm", "BHfwe", "BH", "BY"), n = length(p), w = rep(1, length(p)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{w}{
%%     ~~Describe \code{w} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (p, method = c("holm", "BHfwe", "BH", "BY"), n = length(p), 
    w = rep(1, length(p))) 
{
    w = w/sum(w) * length(w)
    method <- match.arg(method)
    if (method == "fdr") 
        method <- "BH"
    p0 <- p
    if (all(nna <- !is.na(p))) 
        nna <- TRUE
    p <- as.vector(p[nna])
    stopifnot(n >= length(p))
    if (n <= 1) 
        return(p0)
    p0[nna] <- switch(method, bonferroni = pmin(1, p * w), holm = {
        o <- order(p/w)
        ro <- order(o)
        i <- cumsum(w[o])
        pmin(1, cummax(n - i + w[o] * p[o]))[ro]
    }, BHfwe = {
        o <- order(p)
        ro <- order(o)
        i <- cumsum(w[o])
        pmin(1, cummax(n - i + w[o] * p[o]))[ro]
    }, BH = {
        o <- order(p, decreasing = TRUE)
        ro <- order(o)
        i <- cumsum(w[o[length(w):1]])[length(w):1]
        pmin(1, cummin(n/i * p[o]))[ro]
    }, BY = {
        o <- order(p, decreasing = TRUE)
        ro <- order(o)
        q <- sum(1/(cumsum(w[o[length(w):1]])))
        pmin(1, cummin(q * n/i * p[o]))[ro]
    }, none = p)
    attributes(p0)$w = w
    attributes(p0)$method = method
    p0
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
