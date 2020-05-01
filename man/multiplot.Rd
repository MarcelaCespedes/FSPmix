\name{multiplot}
\alias{multiplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multiplot various ggplot objects into a single frame
}
\description{
Multiplot function to plot multiplt ggplot objects onto one plot window. This code is freely available from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
}
\usage{
multiplot(x)
multiplot(plotlist = List.Of.Plots, cols = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ggplot object}
}
\references{
Cookbook for R
}
\author{
Cookbook for R
}

\examples{
# Example taken from website above.
library(ggplot2)
p1 <- ggplot(ChickWeight, aes(x=Time, y=weight, colour=Diet, group=Chick)) +
    geom_line() +
    ggtitle("Growth curve for individual chicks")

# Second plot
p2 <- ggplot(ChickWeight, aes(x=Time, y=weight, colour=Diet)) +
    geom_point(alpha=.3) +
    geom_smooth(alpha=.2, size=1) +
    ggtitle("Fitted growth curve per diet")

# Third plot
p3 <- ggplot(subset(ChickWeight, Time==21), aes(x=weight, colour=Diet)) +
    geom_density() +
    ggtitle("Final weight, by diet")

# Fourth plot
p4 <- ggplot(subset(ChickWeight, Time==21), aes(x=weight, fill=Diet)) +
    geom_histogram(colour="black", binwidth=50) +
    facet_grid(Diet ~ .) +
    ggtitle("Final weight, by diet") +
    theme(legend.position="none")        # No legend (redundant in this graph)

multiplot(p1, p2, p3, p4, cols = 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ggplot2}
\keyword{multiplot}
