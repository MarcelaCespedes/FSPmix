\name{find_MixtureThreshold_Simulation}
\alias{find_MixtureThreshold_Simulation}
\title{
Single bootstrap sample of mixture model and threshold on selected feature for simulation study.
}
\description{
Function used to apply a single bootstrap sample of two component univariate mixture models on feature of interest and find real intersection (valley between the two mixture model means). This function is intended for simulation studies, where the underlying classification of each feature is known; to be called from inside a for-loop for no.bootstrap iterations.
}
\usage{
find_MixtureThreshold(dat, gRoup, boot.size, method=c('diff', 'intersect'),
                          apply.all.dat = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{vector of selected features.}
  \item{gRoup}{Known classification of features (groups A or B).}
  \item{method}{Either 'diff' or 'intersection'. Former argument estimates the difference in the posterior probability. Later argument calls on function findInt() to find the real value intersection between two Guassian mixture means.}
  \item{apply.all.dat}{boolean: FALSE implies take sample with replacement on dat and is the default for bootstrapping; TRUE denotes otherwise and no sampling of dat is performed}
}
\details{
Default arguments are method = 'intersection', apply.all.dat = FALSE in order to implement FSPmix algorithm.
}
\value{
  \item{mix.threshold}{Real value of mixture intersection}
  \item{boot.samp}{Bootstrap sample of data set/ vector}
  \item{mix.mean}{Vector of mixture model means}
  \item{sw}{boolean: Indicator variable denoting if mix.threshold was found (sw == 1); if threhsold was not found sw == 0}
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
dat<- SimulationStudy1()
head(dat)

op<- find_MixtureThreshold_Simulation(dat = dat[,1], gRoup =as.character(dat$group),
                        boot.size = 800, method='intersection')

str(op)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Sample with replacement}% use one of  RShowDoc("KEYWORDS")
