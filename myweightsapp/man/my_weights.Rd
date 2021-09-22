\name{my_weights}
\alias{my_weights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Re-weights of data for regression discontinuity designs
%%  ~~function to do ... ~~
}
\description{
Re-weights the data before appliying regression discontinuity designs
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
my_weights(r, x, c, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{vector r
%%     ~~Describe \code{r} here~~
}
  \item{x}{vector x
%%     ~~Describe \code{x} here~~
}
  \item{c}{bandwith c
%%     ~~Describe \code{c} here~~
}
  \item{n}{length n
%%     ~~Describe \code{n} here~~
}
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
\references{
García, A. (2020) "A Reweighting Approach for Regression     Discontinuity Designs with Discontinuous Density of the Running Variable", Universidad de Los Andes, Chile.
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
function (r, x, c, n)
{
    verificator_list = list(verify_is_vector(r), verify_is_vector(x),
        verify_vector_numeric(r), verify_vector_numeric(x), same_length(r,
            x), verify_c(c))
    error_list = list("R must be a vector with length greater than one",
        "X must be a vector with length greater than one", "R must be a numeric vector",
        "X must be a numeric vector", "R and X must be the same length",
        "C must be numeric of length one")
    problem = FALSE
    counter = 1
    for (i in verificator_list) {
        if (i == FALSE) {
            message(error_list[counter])
            problem = TRUE
        }
        counter = counter + 1
    }
    if (problem == TRUE) {
        message("There is an input problem")
        return("")
    }
    if (missing(n)) {
    }
    else {
        if (length(r) == n) {
            problem = FALSE
        }
        else {
            problem = TRUE
            message("N no es igual al largo de los vectores")
        }
    }
    if (problem == TRUE) {
        message("There is an input problem")
        return("")
    }
    bw_xr = npcdensbw(x ~ r, cykertype = "epanechnikov", cxkertype = "epanechnikov")
    h_x = bw_xr$ybw
    h_r = (bw_xr$xbw)
    n = length(r)
    w = array(NaN, dim = c(n, 1))
    N_ef_w = array(NaN, dim = c(1, 1))
    condition1 = (abs((r - c)/h_r) < sqrt(5))
    w[condition1 == FALSE] = 1
    N_ef_w = sum(condition1)
    condition2 = (abs((r - c)/h_r) <= 2 * sqrt(5))
    id = seq(1:n)
    nece = data.frame(id = id[condition2], r = r[condition2],
        x = x[condition2], estos = condition1[condition2])
    w_si = rep(NaN, sum(condition2))
    for (ii in (which.max(nece$estos)):(sum(nece$estos) + which.max(nece$estos) -
        1)) {
        condition3 = (abs((nece$r - nece$r[ii])/h_r) < sqrt(5))
        ker_x = (0.75 * (1 - 0.2 * ((nece$x[condition3] - nece$x[ii])/h_x)^2)/sqrt(5)) *
            (abs((nece$x[condition3] - nece$x[ii])/h_x) < sqrt(5)) *
            (1/h_x)
        R = cbind(rep(1, sum(condition3)), (nece$r[condition3] -
            nece$r[ii]))
        ker_r = (0.75 * (1 - 0.2 * ((nece$r[condition3] - nece$r[ii])/h_r)^2)/sqrt(5))
        GinvC = chol2inv(chol(crossprod(sqrt(ker_r) * R)))
        num = c(1, 0) \%*\% GinvC \%*\% crossprod(R * ker_r, ker_x)
        if (nece$r[ii] < c) {
            ker_rL = (0.75 * (1 - 0.2 * ((nece$r[condition3] -
                nece$r[ii])/h_r)^2)/sqrt(5)) * (nece$r[condition3] <
                c)
            GinvL = chol2inv(chol(crossprod(sqrt(ker_rL) * R)))
            den2 = c(1, 0) \%*\% GinvL \%*\% crossprod(R * ker_rL,
                ker_x)
            den = ifelse(abs(den2) < 1e-04, 1e-04, den2)
        }
        else {
            ker_rR = (0.75 * (1 - 0.2 * ((nece$r[condition3] -
                nece$r[ii])/h_r)^2)/sqrt(5)) * (nece$r[condition3] >=
                c)
            GinvR = chol2inv(chol(crossprod(sqrt(ker_rR) * R)))
            den2 = c(1, 0) \%*\% GinvR \%*\% crossprod(R * ker_rR,
                ker_x)
            den = ifelse(abs(den2) < 1e-04, 1e-04, den2)
        }
        w_si[ii] = num/den
    }
    send_w = data.frame(id = nece$id, w = w_si)
    w[(nece$id[which.max(nece$estos)]):(nece$id[sum(nece$estos) +
        which.max(nece$estos) - 1])] = send_w$w[!is.na(send_w$w)]
    out = list(w = w, h_x = h_x, h_r = h_r, N_ef_w = N_ef_w)
    return(out)
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
