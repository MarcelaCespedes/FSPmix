\name{spec_sens}
\alias{spec_sens}
\title{
Calculate specificity and sensitivity on FSPmix simulation study predictions
}
\description{
Compute 2x2 confusion matrix as a ratio as well as with numerical values.
}
\usage{
spec_sens(real.dat, pred.dat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{real.dat}{Vector of simulated data which was fed to FSPmix algorithm for a particular feature. Feature classified into groups A or B.}
  \item{pred.dat}{Vector of predictions for a particular feature, classified into either groups A or B.}
}
\value{
  \item{tot.obs}{Total number of observations classified}
  \item{conf.mat}{Sensitivity and specificity values on a 2 x 2 confusion matrix}
  \item{con.mat.ratio}{Sensitivity and specificity percentages, 2 x 2 matrix}
}
\author{
Marcela Cespedes
}
\examples{
set.seed(123)

real.data<- data.frame(ppl = 1:100,
                      group = sample(c("A", "B"), size = 100,prob = c(0.5, 0.5), replace=TRUE ))

pred.data<- data.frame(ppl = 1:100,
                       Pred = sample(c("Pred.A", "Pred.B"), size = 100,prob = c(0.5, 0.5), replace=TRUE ))

spec_sens(real.dat = real.data, pred.dat = pred.data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Specificity}% use one of  RShowDoc("KEYWORDS")
\keyword{Sensitivity}% __ONLY ONE__ keyword per line
