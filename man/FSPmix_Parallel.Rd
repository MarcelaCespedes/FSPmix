\name{FSPmix_Parallel}
\alias{FSPmix_Parallel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Implements FSPmix in parallelised manner (splits load across multiple CPUs)
}
\description{
Implements FSPmix using parallelisation. Prior to calling FSPmix_Parallel(), user must load required parallelisation packages and allocate the number of cores. FSPmix will be implemented with workload split between the number of allocated CPUs.
}
\usage{
FSPmix_Parallel(dat, boot.size = NULL, no.bootstrap=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{Data frame (required) of the number of features to analyse as columns, number of observations as rows}
  \item{boot.size}{Integer, size of bootstrap sample. If boot.size = NULL, then by default the bootstrap sample becomes floor(0.8*dim(dat)[1])}
  \item{no.bootstrap}{Integer, number of times to perform bootstrapping per feature. If no.bootstrap = NULL, then by default no.bootstrap = 500}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item()
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
# Generate data from simulation  study
no.ppl = 500
op<- SimulationStudy1(sEEd = 948575, no.ppl = no.ppl)

#x11()
#op$plot.ColourCoded.genes
#op$plot.WhatFSPmixSees

# set up simulation variables for FSPmix_Sim
feature.dat<- op$dat[, 1:20]
boot.size<- round(no.ppl*0.8) # set bootstrap size sample
no.bootstrap<- 300 # 500 as per paper - this takes a while


library(doParallel)
# Create and registed cluster. See ?makeCluster for further details
# Below we allocate 2 CPU's for

cl<-makeCluster(2, outfile="MulticoreLogging.txt")
cl
registerDoParallel(cl)

# Run FSPmix
sim.op<- FSPmix_Parallel(dat=feature.dat,
                         boot.size = boot.size,
                         no.bootstrap= no.bootstrap)

length(sim.op)
class(sim.op)
names(sim.op)

# Plot output
# multiplot(plotlist = sim.op$all.plots, cols = 5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
