\name{FSPmix_Parallel}
\alias{FSPmix_Parallel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Implements FSPmix in parallelised manner (splits load across multiple CPUs)
}
\description{
Implements FSPmix using parallelisation. Prior to calling FSPmix_Parallel(), user must load required parallelisation packages and allocate the number of cores. FSPmix will be implemented with workload split between the number of allocated CPUs (workers). Input requires the data and parameters required to enable FSPmix algorithm, this function returns a list of the combined output from all workers.
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
\value{
  \item{}{List with length of the number of features containing the following objects}
  \item{two.groups}{Boolean variable denoting if two groups were found, two.groups = TRUE, or otherwise; two.groups = FALSE}
  \item{summ.op}{Data frame summary of the feature, including indicator variable if two groups were found (two.groups), summary statistics of the seperation interval, as well as bootstrapped means of the mixture model components}
  \item{SampDat_Store}{Data frame with bootstrapped samples}
  \item{Plot}{ggplot object which shows the overlapped bootstrapped densities in black, bootstrapped mixture model component means (in gray) as well as the threshold interval (values within the blue lines).}
  \item{Classification.Pred}{Data frame containing the classification of each observation (or person) if two.groups = TRUE. Classification categories are Pred.A and Pred.B according to their feature values determined by FSPmix; and Pred.C denoting an observation remains unclassified as their feature value lies within the seperation interval}
}
\author{
Marcela Cespedes
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
no.bootstrap<- 300


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

length(sim.op) # length of number of Features to analyse
class(sim.op)

# Show which features FSPmix found two groups
plot.all<- list()
for(i in 1:20){
  tmp.op<- sim.op[[i]]
  plot.all[[i]]<- tmp.op$Plot
  print(tmp.op$two.groups)
}

multiplot(plotlist = plot.all, cols = 5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Parallelisation}% use one of  RShowDoc("KEYWORDS")
\keyword{Mixture Model}% __ONLY ONE__ keyword per line
\keyword{Feature selection}
\keyword{Prediction}
\keyword{Classification}
