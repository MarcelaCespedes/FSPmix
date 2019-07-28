\name{spec_sens}
\alias{spec_sens}
\title{
Calculate specificity and sensitivity on FSPmix simulation study predictions
}
\description{
Compute confusion table
}
\usage{
spec_sens(real.dat, pred.dat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{real.dat}{Simulated data }
  \item{pred.dat}{
%%     ~~Describe \code{pred.dat} here~~
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
function (real.dat, pred.dat)
{
    real.dat <- real.dat[order(real.dat$ppl), ]
    pred.dat <- pred.dat[order(pred.dat$ppl), ]
    dat.sol <- subset(real.dat, ppl \%in\% pred.dat$ppl)
    tot.obs <- dim(dat.sol)[1]
    tot.obs
    sol.A <- subset(dat.sol, group == "A")
    tot.A = dim(sol.A)[1]
    tot.A
    p.A <- subset(pred.dat, Pred == "Pred.A")
    pred.A <- dim(p.A)[1]
    pred.A
    tp = length(intersect(sol.A$ppl, p.A$ppl))
    tp
    pred.A.sol.B = abs(pred.A - tp)
    pred.A.sol.B
    sol.B <- subset(dat.sol, group == "B")
    tot.B = dim(sol.B)[1]
    tot.B
    p.B <- subset(pred.dat, Pred == "Pred.B")
    pred.B <- dim(p.B)[1]
    pred.B
    tn <- length(intersect(sol.B$ppl, p.B$ppl))
    tn
    pred.B.sol.A = abs(pred.B - tn)
    pred.B.sol.A
    conf.mat <- data.frame(pred.A = c(NA, NA), pred.B = c(NA,
        NA))
    conf.mat[1, 1] <- tp
    conf.mat[2, 2] <- tn
    conf.mat[1, 2] <- pred.B.sol.A
    conf.mat[2, 1] <- pred.A.sol.B
    row.names(conf.mat) <- c("sol.A", "sol.B")
    conf.mat.ratio <- conf.mat
    conf.mat.ratio[1, ] <- conf.mat.ratio[1, ]/tot.A
    conf.mat.ratio[2, ] <- conf.mat.ratio[2, ]/tot.B
    return(list(tot.obs = tot.obs, conf.mat = conf.mat, conf.mat.ratio = round(conf.mat.ratio,
        3)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
