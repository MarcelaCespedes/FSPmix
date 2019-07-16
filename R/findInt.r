# relies on this function
# @@ Alternative method - just find the intersection between the two fitted Gaussians.
#' Find intersection(s) between two Gaussians, possibly with mixing proportions.
#' @param m1, m2: fitted means of two components
#' @param sd1, sd2: fitted standard deviations of two components
#' @param p1, p2: fitted proportions for the two components
#' @param filter {boolean}: if TRUE, try to find the "middle" & real-valued intersection.
#' There can be more than one intersection between two Gaussians,
#' We try to find the Real one that is between the two means.
#' i.e. the "middle" one.
#' I have the working .. somewhere .., I dug this out of my code.
#' Similar to this:: (but with mixing proportions)
#' https://stats.stackexchange.com/questions/311592/how-to-find-the-point-where-two-normal-distributions-intersect
#' @examples
#' library(mixtools)
#' y <- rnormmix(1e3, lambda=c(0.2, 0.8), mu=c(-0.5, 0.5), sigma=c(0.1, 0.3))
#' hist(y)
#' m <- normalmixEM(y, k=2)
#' thr <- findInt(m$mu[1], m$mu[2], m$sigma[1], m$sigma[2], m$lambda[1], m$lambda[2])
#' DenPlot(m, BiomThresh=thr, BiomName='y', BinWidth=0.1, BiomRange=waiver())
findInt <- function (m1, m2, sd1, sd2, p1=1, p2=1, filter=T) {
  a <- 1/(2*sd1^2) - 1/(2*sd2^2)
  b <- m2/(sd2^2) - m1/(sd1^2)
  c <- m1^2 /(2*sd1^2) - m2^2 / (2*sd2^2) - log(sd2/sd1)
  c <- m1^2 /(2*sd1^2) - m2^2 / (2*sd2^2) - log((sd2*p1)/(sd1*p2))
  rr <- polyroot(c(c,b,a))

  sw = 0
  # find the one between the means
  if (filter) {
    rr <- Re(rr[sapply(Im(rr), function (x) isTRUE(all.equal(x, 0)))])
    rr <- rr[rr >= min(m1, m2) & rr <= max(m1, m2)]
    #rr <- rr[rr >= m1 & rr <= m2]
    if (length(rr) != 1){
      sw = 1
      #warning(sprintf("Failed to find intersection between the means (%i candidates)",
      #             length(rr))) # supress warnings, throw out switch instead
      #stop(sprintf("Failed to find intersection between the means (%i candidates)", length(rr)))
    }
  }
  return(list(rr = rr, sw = sw))
}
