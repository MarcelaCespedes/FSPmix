\name{FSPmix}
\alias{FSPmix}

\title{
Feature selection and prediction (FSPmix) algorithm
}
\description{
Implements the FSPmix algorithm on real data in a serial manner, that is, using a single CPU. For parallelised version use FSPmix_Parallel().
}
\usage{
FSPmix(dat, boot.size = NULL, no.bootstrap=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{Data frame containing features of interest (columns) of each observation (rows)}
  \item{boot.size}{Integer, size of bootstrap sample. If boot.size = NULL, then by default the bootstrap sample becomes floor(0.8*dim(dat)[1])}
  \item{no.bootstrap}{Integer, number of times to perform bootstrapping per feature. If no.bootstrap = NULL, then by default no.bootstrap = 500}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{}{List with length of the number of features containing the following objects}
  \item{two.groups}{Boolean variable denoting if two groups were found, two.groups = TRUE, or otherwise; two.groups = FALSE}
  \item{summ.op}{Data frame summary of the feature, including indicator variable if two groups were found (two.groups), summary statistics of the seperation interval, as well as bootstrapped means of the mixture model components}
  \item{SampDat_Store}{Data frame with bootstrapped samples}
  \item{Plot}{ggplot object which shows the overlapped bootstrapped densities in black, bootstrapped mixture model component means (in gray) as well as the threshold interval (values within the blue lines).}
  \item{Classification.Pred}{Data frame containing the classification of each observation (or person) if two.groups = TRUE. Classification categories are Pred.A and Pred.B according to their feature values determined by FSPmix; and Pred.C denoting an observation remains unclassified as their feature value lies within the seperation interval}
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
library(FSPmix)
library(ggplot2)
library(reshape2)

## ***************************************************
## Process real data
# Apply on real data: Downloaded from https://milxview.csiro.au/MilxXplore/Demo/xplorer_study/AIBL/Demo
# Monday 13 August 2019
pet.ROI<- read.csv(file = "PiB_uptake.csv", header=TRUE, sep = ",")
# Remove Vermis, cerebellum and other columns which are not needed
pet.ROI<- pet.ROI[, -c(5,6, 97:123)]
# combine left and right ROIs together to increase sample size
# Currently have 16 observations per ROI - by combining the left and right, we double to 32 obs per ROI
left.ROIs<- pet.ROI[, seq(from =5, to=93, by=2)]
colnames(left.ROIs)<- substr(colnames(left.ROIs), start=1, stop = nchar(colnames(left.ROIs))-5)
colnames(left.ROIs)[7]<- "Inferior_frontal_gyrus.triangular_part"
right.ROIs<- pet.ROI[, seq(from =6, to=94, by=2)]
colnames(right.ROIs)<- substr(colnames(right.ROIs), start=1, stop = nchar(colnames(right.ROIs))-6)

feature.dat<- rbind(left.ROIs, right.ROIs)
##
## Number of features
dim(feature.dat)[2]  # 45

## View the raw data
plot.fd<- melt(feature.dat)
ggplot(plot.fd, aes(x = value)) + geom_density(alpha = 0.8) +
  facet_wrap(.~variable)+ #, scales="free") +
  theme_bw() +
  ggtitle("Combined Left & Right ROIs densities")

## **********************************************
## Run FSPmix
no.ppl<- dim(feature.dat)[1]
boot.size<- round(no.ppl*0.8)
no.bootstrap<- 300

sim.op<- FSPmix(dat=feature.dat,
                boot.size = boot.size,
                no.bootstrap= no.bootstrap)

length(sim.op) # length of number of Features to analyse
class(sim.op)

# Show which features FSPmix found two groups
plot.all<- list()
for(i in 1:45){
  tmp.op<- sim.op[[i]]
  plot.all[[i]]<- tmp.op$Plot
  print(tmp.op$two.groups)
}
multiplot(plotlist = plot.all, cols = 7)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Mixture Model}% __ONLY ONE__ keyword per line
\keyword{Feature selection}
\keyword{Prediction}
\keyword{Classification}
