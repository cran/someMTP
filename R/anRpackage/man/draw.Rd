\name{draw}
\alias{draw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
draw(object, what = c("all", "ordVsP", "stepVsR"), pdfNane = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{what}{
%%     ~~Describe \code{what} here~~
}
  \item{pdfNane}{
%%     ~~Describe \code{pdfNane} here~~
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
function (object, what = c("all", "ordVsP", "stepVsR"), pdfNane = NULL) 
{
    if (missing(what)) 
        what <- "all"
    else what <- match.arg(what, c("all", "ordVsP", "stepVsR"))
    if (!is.null(pdfNane)) 
        pdf(pdfNane, width = 10)
    if (what == "all") 
        par(mfrow = c(1, 2))
    if (what \%in\% c("all", "ordVsP")) {
        par(cex = 1.5)
        plot(object@ord, -log10(object@p), xlab = "Ordering values", 
            ylab = "-log10(p-values)", pch = 20, axes = TRUE, 
            main = "Ordering Criterion vs -log10(p-values)", 
            col = object@rej + 1)
        abline(-log10(object@q), 0, col = "gray", lwd = 2)
        vline = length(object@rej) - which.max(which(cumsum(object@rej[object@idOrd[length(object@idOrd):1]]) == 
            0)) + 1
        abline(v = ifelse(length(vline) > 0, object@ord[object@idOrd[vline]], 
            0), col = "gray", lwd = 2)
    }
    if (what \%in\% c("all", "stepVsR")) {
        par(cex = 1.5)
        plot(cumsum(object@p[object@idOrd] <= object@alphaprime), 
            xlab = "Steps", ylab = "# of rejections", lwd = 2, 
            axes = TRUE, main = "Step vs Number of Rejections", 
            type = "l")
        abline(0, 1/(2 - object@q), col = "red", lwd = 2)
        abline(0, 1, col = "gray", lwd = 2)
        legend(x = 0.05, y = sum(object@p <= object@alphaprime), 
            legend = c("Maximum", "Stop if below"), lty = 1, 
            col = c("gray", "red"), lwd = 2, bty = "n", cex = 1)
    }
    if (!is.null(pdfNane)) 
        dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
