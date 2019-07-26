\name{SimulationStudy1}
\alias{SimulationStudy1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Synthetic normalised gene expression data for 20 genes
}
\description{
As described in xx paper
}
\usage{
SimulationStudy1(sEEd = 999, no.ppl= 1000)
}

\arguments{
  \item{sEED}{Positive integer; denoting the random generator state for number random generation (see ?set.seed).}
  \item{no.ppl}{Positive integer; Number of participants to simulate normalised gene expression data.}
}
\details{
First simulation for FSPmix as described in xx.
}
\value{
  \item{dat}{data.frame object of the 20 simulated gene expressions (columns) for no.ppl (rows). Additional columns include participant identifier and latent classification as group.}
  \item{plot.ColouredCoded.genes}{ggplot object of the density of the 20 features coloured coded into 2 classifications.}
  \item{plot.WhatFSPmixSees}{ggplot object of uncoloured densities of the distributions that FSPmix will categorise into either 'A' or 'B' groups.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Marcela Cespedes
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
op<- SimulationStudy1(sEEd = 487, no.ppl = 1000)

head(op$dat)

# View 20 simulated genes (features) coloured coded into groups 'A' and 'B'
op$plot.ColouredCoded.genes

# View 20 simulated genes (features) densities
op$plot.WhatFSPmixSees
}
\keyword{ Simulation }% use one of  RShowDoc("KEYWORDS")
