\name{step.adj}
\alias{step.adj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multipicity correction for Stepwise Selected models
}
\description{
Corrects the false discoveries due to model selection. It works with models of class \code{glm} and selected with function \code{step \{stats\)}.
}
\usage{
step.adj(object, MC = 1000, scope = NULL, scale = 0, direction = c("both", "backward", "forward"), trace = 0, keep = NULL, steps = 1000, k = 2) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{object of class \code{glm}. Note that
  \code{formula} have to write by variables name like \code{y~var1+var2+var3},
 \code{data} is a data.frame (see example below),
\code{offset} is not yet implemented, avoid its use,
\code{glm(formula, data, family=gaussian)} produce the same result of lm(formula, data), then linear model can be allways performed
}
  \item{MC}{
number of random permutations for the dependent variable
}
  \item{scope}{as in function \code{step}}
  \item{scale}{as in function \code{step}}
  \item{direction}{as in function \code{step}}
  \item{trace}{as in function \code{step}}
  \item{keep}{as in function \code{step}}
  \item{steps}{as in function \code{step}}
  \item{k}{as in function \code{step}, other arguments are not implemented yet.}
}


\details{
It performs anova function (stats library) on the model selected by function step vs the null model with the only intercept 
and it corrects for multiplicity.
For \code{lm} models and gaussian \code{glm} models it computes a F-test, form other models it uses Chisquare-test (see also \code{anova.glm} and \code{anova.lm} help).
}
\value{
An \code{anova} table with an extra column reporting the corrected p-value
}
\references{
L. Finos, C. Brombin, L. Salmaso (2010). Adjusting stepwise p-values in generalized linear models. Communications in Statistics - Theory and Methods. 
}
\author{
Livio Finos and Chiara Brombin
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{glm}}, \code{\link{anova}}
}
\examples{
set.seed(17)
y=rnorm(10)
x=matrix(rnorm(50),10,5)
#define a data.frame to be used in the glm function
DATA=data.frame(y,x)
#fit the model on a toy dataset
mod=glm(y~X1+X2+X3+X4+X5,data=DATA)

#select the model using function step
mod.step=step(mod, trace=0)
#test the selected model vs the null model
anova(glm(y~1, data=DATA),mod.step,test="F")

#step.adj do the same, but it also provides multiplicity control
step.adj(mod,MC=101, trace=0)
}
\keyword{htest}

