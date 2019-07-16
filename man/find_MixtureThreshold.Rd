\name{find_MixtureThreshold}
\alias{find_MixtureThreshold}
Derive the intersection between two univariate Gaussian mixture models
\title{
Bootstrap mixture model and threshold on selected feature.
}
\description{
Bootstrap two component univariate mixture models on feature of interest and find real intersection (valley between the two mixture model means).
}
\usage{
find_MixtureThreshold(dat, gRoup, boot.size, method=c('diff', 'intersect'),
                          apply.all.dat = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{Data frame of }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
