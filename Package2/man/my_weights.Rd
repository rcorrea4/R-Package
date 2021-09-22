\name{my_weights}
\alias{my_weights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Re-weights of data for regression discontinuity designs
%%  ~~function to do ... ~~
}
\description{
Re-weights the data before appliying regression discontinuity designs
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
my_weights(r,x,c,n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{vector r
%%     ~~Describe \code{x} here~~
}
  \item{x}{vector x}
  \item{c}{bandwith c}
  \item{n}{length of vectors r and x}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{w}{weighter}
  \item{h_x}{rellenar}
  \item{h_r}{rellenar}
  \item{N_ef_w}{rellenar}

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{García, A. (2020) "A Reweighting Approach for Regression Discontinuity Designs with Discontinuous Density of the Running Variable", Universidad de Los Andes, Chile.
%% ~put references to the literature/web site here ~
}
\author{
Raimundo Correa D. and Vicente Lehmann C. under the supervision of Andres García E.
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
my_weights(c(1,2),c(1,2),2,2)
function (x)
{
  }
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
