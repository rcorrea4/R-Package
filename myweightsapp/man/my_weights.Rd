\name{my_weights}
\alias{my_weights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Re-Weights Data before Applying Regression Discontinuity Designs
%%  ~~function to do ... ~~
}
\description{
my_weights re-weightes data that will be used in a regression discontinuity design that the sample have manipulated data and doesn't aprove the McCrary's test. The sample is re-weighted as proposed by Garcia (2020) in order to be able to make interpretations about the variable of interest after having applied the regression discontinuity design.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
my_weights(r, x, c, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{numeric vector r. Represents the assignment variable.
%%     ~~Describe \code{r} here~~
}
  \item{x}{numeric matrix x. Represents the observables variables.
%%     ~~Describe \code{x} here~~
}
  \item{c}{scalar c. Represents the cutoff point that indicates from which value the sample receives or does not receive treatment.
%%     ~~Describe \code{c} here~~
}
  \item{n}{optional scalar n. Represents the number of rows of r and x.
%%     ~~Describe \code{n} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{my_weights returns a list of weighters for each row of the sample. Multiplying these values with the original data, a continuous density sample is obtained. With this sample it is posible to make interpretations of the variable of interest after the discontinuous regression.
  \item{w}{list of the weighters values}
  \item{h_x}{bandwith x}
  \item{h_r}{bandwith r}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Garcia, A. (2020) "A Reweighting Approach for Regression     Discontinuity Designs with Discontinuous Density of the Running Variable", Universidad de Los Andes, Chile.
%% ~put references to the literature/web site here ~
}
\author{
Raimundo Correa D. and Vicente Lehmann C. under the supervision of Andres Garcia E.
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

x=cbind(rnorm(100),rnorm(100,2))
r=rnorm(100,0,6)
c=0
n=100

my_weights(r,x,c)
my_weights(r,x,c,n)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
