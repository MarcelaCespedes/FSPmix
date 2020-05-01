\name{PiB_uptake}
\alias{PiB_uptake}
\docType{data}
\title{PET PiB AIBL demonstration data set}
\description{
Demonstration data set of beta-amyloid detected in the human brain measured via Positron Emission Topography imaging with Pittsburgh compound B (PiB) tracer. Data set comprises of PET PiB measurements for five participants from the AIBL study who have between 2-4 observations each. PET PiB readings pertain to brain beta-amyloid measurements parcellated onto 116 brain regions of interest (ROIs).

Data set was downloaded on Monday 13 August 2019, from \url{https://milxview.csiro.au/MilxXplore/Demo/xplorer_study/AIBL/Demo}.
}
\usage{pet.ROI<- read.csv(file = "PiB_uptake.csv", header=TRUE, sep = ",")}
\format{
  Data frame with PET PiB brain measurements for five AIBL participants. Human brain was parcellated onto 116 ROIs as defined by the Automated Anatomical Atlas (AAL). Column names of data set are described below
  \describe{
    \item{\code{CaseID}}{Participant identification number}
    \item{\code{TP}}{Time point at imaging, baseline is denoted by T0, and subsequent follow-ups are T1-T3}
    \item{\code{Diagnosis}}{Clinical diagnosis of participant at each TP. Participants are either healthy (Normal control), mild cognitive impaired (MCI patient) or diagnosed with Alzheimer's disease (AD patient)}
    \item{\code{Age.at.MRI}}{In order to parcellate the brain into regions, both PET and Magnetic Resonance Imaging (MRI) is required. Age.at.MRI denotes the age of the participant at time of MRI scan}
    \item{\code{Filename}}{Name of imaging file}
    \item{\code{Neocortex}}{PET PiB uptake of the neocortex}
    \item{\code{Columns 7-123}}{Denote PET PiB update on 116 ROIs of the left and right brain hemisphere.}
    \item{\code{SUV.in.Reference.Region}}{Total Standard Uptake Value (SUV) in reference brain region}
  }
}
\details{
MilxXplore \href{https://milxview.csiro.au/MilxXplore/}{website} for further details, Biomedical Informatics Team who developed MilxXplore \href{https://aehrc.com/biomedical-informatics/}{website}, and AIBL longitudinal study of ageing website \href{https://aibl.csiro.au/}{website}. Brain regions for AAL see
Tzourio-Mazoyer, N., Landeau, B., Papathanassiou, D., Crivello, F., Etard, O., Delcroix, N., Mazoyer, B. and Joliot, M., 2002. Automated anatomical labeling of activations in SPM using a macroscopic anatomical parcellation of the MNI MRI single-subject brain. Neuroimage, 15(1), pp.273-289.
}
\source{
Bourgeat, P., Dore, V., Villemagne, V.L., Rowe, C.C., Salvado, O. and Fripp, J., 2013. MilxXplore: a web-based system to explore large imaging datasets. Journal of the American Medical Informatics Association, 20(6), pp.1046-1052
}
\examples{
library(FSPmix)
library(reshape2)
library(ggplot2)

pet.ROI<- read.csv(file = "PiB_uptake.csv", header=TRUE, sep = ",")
# Remove Vermis, cerebellum and other brain regions which are not needed
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
  ggtitle("Combined Left & Right ROIs densities") +
  xlab("Amyloid (detected by PiB uptake)")
}
\keyword{Beta-amyloid}% use one of  RShowDoc("KEYWORDS")
\keyword{PET}
\keyword{MRI}
\keyword{Neuroimaging}
\keyword{Alzheimer's disease}
\keyword{ROI}
