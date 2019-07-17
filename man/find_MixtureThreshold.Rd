\name{find_MixtureThreshold}
\alias{find_MixtureThreshold}
Derive the intersection between two univariate Gaussian mixture models
\title{
Single bootstrap iteration of application of mixture model and threshold value estimation on selected feature. To be applied on real data, as known classification of groups is unknown.
}
\description{
Bootstrap two component univariate mixture models on feature of interest and find real intersection (valley between the two mixture model means).
}
\usage{
find_MixtureThreshold(dat, boot.size, method=c('diff', 'intersect'),
                          apply.all.dat = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{vector of selected features.}
  \item{method}{Either 'diff' or 'intersection'. Former argument estimates the difference in the posterior probability. Later argument calls on function findInt() to find the real value intersection between two Guassian mixture means.}
  \item{apply.all.dat}{boolean: FALSE implies take sample with replacement on dat and is the default for bootstrapping; TRUE denotes otherwise and no sampling of dat is performed}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{mix.threshold}{Indicator variable denoting if }
 return(list(mix.threshold=mix.threshold$rr, boot.samp = samp.dat, mix.means = mix.means1,
              sw = sw))
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Marcela Cespedes and Amy Chan
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(mixtools)
library(reshape2)
source("SimulationStudy2.r")

# Simulation study for 20 Genes
dat<- SimulationStudy2()
head(dat)

op<- find_MixtureThreshold(dat = dat[,1],
                        boot.size = 800, method='intersection')

str(op)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
