\name{fdrOrd}
\alias{fdrOrd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fdrOrd(p, q = 0.01, ord = NULL, GD = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{q}{
%%     ~~Describe \code{q} here~~
}
  \item{ord}{
%%     ~~Describe \code{ord} here~~
}
  \item{GD}{
%%     ~~Describe \code{GD} here~~
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
function (p, q = 0.01, ord = NULL, GD = FALSE) 
{
    m <- length(p)
    alphaprime = q
    if (!is.null(ord)) 
        o <- order(ord, decreasing = T)
    else {
        o <- 1:length(p)
        ord = o[length(p):1]
    }
    if (GD) 
        alphaprime = alphaprime/(sum(1/((1 + (1:length(p)))/(2 - 
            alphaprime))))
    else alphaprime = q
    u <- cumsum(p[o] > alphaprime)
    h <- rep(0, m)
    names(h) <- names(p)
    w = min(c(which((u * alphaprime/(1 - alphaprime))/((1:m) - 
        u) > q) - 1, m))
    h[o[1:w]] <- 1
    h[p > alphaprime] <- 0
    out <- new("someMTP.object")
    out@rej = h == 1
    out@p = p
    out@ord = ord
    out@idOrd = o
    out@MTP = "fdrOrd"
    out@GD = GD
    out@q = q
    out@k = NULL
    out@alpha = NULL
    out@alphaprime = alphaprime
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
