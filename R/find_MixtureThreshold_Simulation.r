#' Derive the intersection between two univariate Gaussian mixture models
#' dat: vector of continuous observations to be classified
#' gRoup: vector of string characters - known membership of observations
#' bootzise: integer, size of sub-sample of data
#' method: two options as to how to find threshold

find_MixtureThreshold_Simulation<- function(dat, gRoup, boot.size, method=c('diff', 'intersect'),
                          apply.all.dat = FALSE){

  method <- match.arg(method)
  sw = 0
  #if(!is.null(boot.size)){  # set default bootstrap sample size
  #  boot.size = floor(0.8*dim(dat)[1])
  #}
  full.dat<- samp.dat<- data.frame(dat=dat, gRoup=gRoup)

  if(apply.all.dat == FALSE){
    # Not proper bootstrap if I do sampling WITHOUT replacement
    #samp.dat<- sample_n(full.dat, boot.size, replace=FALSE)
    samp.dat<- sample_n(full.dat, boot.size, replace=TRUE)
  }

  # supress the output - methinks with capture.output() : BUT this turns the
  # whole output to a character vector
  # @@ can invisible(capture.output(mix.mod <- ...)) & it will work
  invisible(capture.output(mix.mod<- normalmixEM(samp.dat[,1], k=2, epsilon = 0.05) ))
  #str(mix.mod)

  mix.threshold <- NA
  if (method == 'intersect') { # intersect two gaussians
    mix.threshold <- findInt(mix.mod$mu[1], mix.mod$mu[2],
                             mix.mod$sigma[1], mix.mod$sigma[2],
                             mix.mod$lambda[1], mix.mod$lambda[2])

    sw = mix.threshold$sw
  } else { # diff in posterior probability
    # note - if speed becomes an issue I can do some things here.
    post.df<- data.frame(cbind(x=mix.mod$x, mix.mod$posterior))
    post.df2 <- post.df[order(post.df$x),]
    post.df2$diff <- post.df2$comp.2- post.df2$comp.1

    rownames(post.df2)<- 1:dim(samp.dat)[1]
    # Threshold estimated via mixture model: T_e
    mix.threshold <- post.df2[which.min(abs(post.df2$diff)),1]
    mix.threshold

  }

  ### Attempt to correct label switching
  mix.means<- c(mix.mod$mu[1], mix.mod$mu[2])
  mix.means1 = mix.means[order(mix.means)]

  #x11()
  #DenPlot(MixData = mix.mod, BiomRange = c(1,2,3,4),
  #        BiomName = "Test 1st Var dat.1",
  #        BiomThresh = round(mix.threshold,3),BinWidth = 0.1)

  return(list(mix.threshold=mix.threshold$rr, boot.samp = samp.dat, mix.means = mix.means1,
              sw = sw))
}


