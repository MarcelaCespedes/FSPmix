\name{findInt}
\alias{findInt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Finds the intersection between the two fitted Gaussians.
}
\description{
Find the real value intersection (deepest part of the valley) between the means of two univariate Guassian distributions
}
\usage{
findInt(mu1, mu2, sd1, sd2, p1, p2, filter=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu1, mu2}{fitted means of two components}
  \item{sd1, sd2}{fitted standard deviations of two components}
  \item{p1, p2}{fitted proportions for the two components}
  \item{filter}{boolean: if TRUE, try to find the "middle" & real-valued intersection}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{rr}{Real value of the intersection between the means of the two Guassian distributions}
  \item{sw}{boolean: value 0 if real continuous intersection value was found; otherwise if no real (and imaginary) value was found sw = 0}
}
\references{
Similar to this (but with mixing proportions)
https://stats.stackexchange.com/questions/311592/how-to-find-the-point-where-two-normal-distributions-intersect
}
\author{
Amy Chan
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
y <- rnormmix(1e3, lambda=c(0.2, 0.8), mu=c(-0.5, 0.5), sigma=c(0.1, 0.3))
hist(y)
m <- normalmixEM(y, k=2)
thr <- findInt(m$mu[1], m$mu[2], m$sigma[1], m$sigma[2], m$lambda[1], m$lambda[2])
thr
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Intersection}% use one of  RShowDoc("KEYWORDS")
\keyword{Mixture model}% __ONLY ONE__ keyword per line
