\name{PiB_uptake}
\alias{PiB_uptake}
\docType{data}
\title{
PET PiB tracer data
}
\description{
Positron Emission Topography (PET) with PiB (Pittsburgh Compound-B) tracer standard uptake value ratio (SUVR) brain region measures.
}
\usage{data("PiB_uptake")}
\format{
  A data frame with 17 observations on the following 123 variables; four colums of demographic information and generic imaging file description, as well as measures on 116 brain regions.
  \describe{
    \item{\code{CaseID}}{a numeric vector: participant ID}
    \item{\code{TP}}{a numeric vector: time point}
    \item{\code{Diagnosis}}{character vector: clinical diagnosis, healthy denoted as Healthy Control, mild cognitive impairment as MCI patient, and Alzheimer's disease denoted as AD patient}
    \item{\code{Age.at.MRI}}{a numeric vector: age of participant at MRI scan}
    \item{\code{Filename}}{a character vector: information of DICOM image file location}
  }
}
\details{
SUVR measures for 116 brain regions, parcellated according to the Automated Anatomical Atlas (AAL).
}

\examples{
data(PiB_uptake)
}
\keyword{datasets}
