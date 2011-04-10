\name{step.adj}
\alias{step.adj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
step.adj(object, MC = 1000, scope = NULL, scale = 0, direction = c("both", "backward", "forward"), trace = 0, keep = NULL, steps = 1000, k = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{MC}{
%%     ~~Describe \code{MC} here~~
}
  \item{scope}{
%%     ~~Describe \code{scope} here~~
}
  \item{scale}{
%%     ~~Describe \code{scale} here~~
}
  \item{direction}{
%%     ~~Describe \code{direction} here~~
}
  \item{trace}{
%%     ~~Describe \code{trace} here~~
}
  \item{keep}{
%%     ~~Describe \code{keep} here~~
}
  \item{steps}{
%%     ~~Describe \code{steps} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
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
function (object, MC = 1000, scope = NULL, scale = 0, direction = c("both", 
    "backward", "forward"), trace = 0, keep = NULL, steps = 1000, 
    k = 2) 
{
    family = object$family
    if (is.character(family)) 
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    object$family = family
    TEST = ifelse((object$family[1] == "gaussian") | (is(object, 
        "lm")), "F", "Ch")
    direction <- match.arg(direction)
    p = rep(NA, MC)
    formula0 = object$formula
    formula0[[3]] = 1
    mod0 = glm(formula = formula0, family = object$family, control = object$control, 
        data = object$data)
    mod.step = step(object, scope = scope, scale = scale, direction = direction, 
        trace = trace, keep = keep, steps = steps, k = k)
    anova.sw = anova(mod0, mod.step, test = TEST)
    p[MC] = anova.sw[[6]][2]
    redo <- function(mod, mod0, TEST, scope, scale, direction, 
        keep, steps, k) {
        mod$data[, 1] = sample(mod$data[, 1])
        mod = update(mod, data = mod$data)
        mod.step = step(mod, scope = scope, scale = scale, direction = direction, 
            trace = 0, keep = keep, steps = steps, k = k)
        anova.swr = anova(mod0, mod.step, test = TEST)[[5 + (TEST == 
            "F")]][2]
        p = min(anova.swr, 1, na.rm = T)
    }
    p[-MC] = replicate(MC - 1, redo(object, mod0, TEST, scope, 
        scale, direction, keep, steps, k))
    attributes(anova.sw)$heading[1] = paste("Analysis of Deviance Table on Selected Model \nwith Multiplicity Correction (weak FWE-control)\nbased on", 
        MC, "Random permutations\n")
    anova.sw[[" sw.adj P"]] = c(NA, mean(p <= p[MC]))
    return(anova.sw = anova.sw)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
