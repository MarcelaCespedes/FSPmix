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
  \item{dat}{Vector of selected features}
  \item{method}{Either 'diff' or 'intersection'. Former argument estimates the difference in the posterior probability. Later argument calls on function findInt() to find the real value intersection between two Guassian mixture means.}
  \item{apply.all.dat}{boolean: FALSE implies take sample with replacement on dat and is the default for bootstrapping; TRUE denotes otherwise and no sampling of dat is performed}
}
\value{
  \item{mix.threshold}{Indicator variable denoting if }
 return(list(mix.threshold=mix.threshold$rr, boot.samp = samp.dat, mix.means = mix.means1,
              sw = sw))
}
\author{
Marcela Cespedes and Amy Chan
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
