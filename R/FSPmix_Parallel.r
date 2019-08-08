#' FSPmix_Parallel.r implements the FSPmix algorithm in a parallelised manner
#' utilising multiple CPU's (as specified by the user)
#' Input
#' dat: data.frame of features for FSPmix to search. Each column represents a feature and each row denotes an observation
#' boot.size: positive integer, size of boot strap sample
#' no.bootstrap: positive integer, number of times to bootstrap

FSPmix_Parallel<- function(dat,
                           boot.size = NULL, no.bootstrap=NULL){

  if(is.null(boot.size)){
    boot.size<- round(dim(dat)[1]*0.8) # size of bootstrap sample is 80% of participant size
  }

  if(is.null(no.bootstrap)){
    no.bootstrap<- 100 #500  # as per manuscript (takes a while)
  }

  library(ggplot2)
  library(reshape2)
  library(dplyr)
  library(mixtools)
  library(grDevices)
  fade <- function(colors,alpha) {  # <-- to plot simulation densities
    rgbcols <- col2rgb(colors)
    rgb(rgbcols[1,],rgbcols[2,],rgbcols[3,],alpha/100*255,max=255)
  }

  ## ****************************************************************
  ## Begin FSPmix algorithm
  ##

  no.gene<- dim(dat)[2]
  final.op <- foreach(p=1:no.gene, .errorhandling='pass') %dopar% {

    library(dplyr)
    library(reshape2)
    library(ggplot2)
    library(mixtools)

    ## Include all functions within the for-loop: in order to make them accessible
    ## for all nodes
    find_MixtureThreshold_Simulation<- function(dat, boot.size,
                                                method=c('diff', 'intersect'),
                                                apply.all.dat = FALSE){

      method <- match.arg(method)
      sw = 0

      full.dat<- samp.dat<- data.frame(dat=dat)
      if(apply.all.dat == FALSE){
        samp.dat<- sample_n(full.dat, boot.size, replace=TRUE)
      }

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

      ### Correct label switching
      mix.means<- c(mix.mod$mu[1], mix.mod$mu[2])
      mix.means1 = mix.means[order(mix.means)]

      return(list(mix.threshold=mix.threshold$rr, boot.samp = samp.dat, mix.means = mix.means1,
                  sw = sw))
    }

    findInt <- function (m1, m2, sd1, sd2, p1=1, p2=1, filter=T) {
      a <- 1/(2*sd1^2) - 1/(2*sd2^2)
      b <- m2/(sd2^2) - m1/(sd1^2)
      c <- m1^2 /(2*sd1^2) - m2^2 / (2*sd2^2) - log(sd2/sd1)
      c <- m1^2 /(2*sd1^2) - m2^2 / (2*sd2^2) - log((sd2*p1)/(sd1*p2))
      rr <- polyroot(c(c,b,a))

      sw = 0
      # find the one between the means
      if (filter) {
        # take only the Real roots
        rr <- Re(rr[sapply(Im(rr), function (x) isTRUE(all.equal(x, 0)))])
        if (length(rr) == 0) {
          sw <- 1
        } else {
          inBetween <- rr >= min(m1, m2) & rr <= max(m1, m2)
          if (sum(inBetween) == 1) # we found exactly one
            rr <- rr[inBetween]
          else
            rr <- sample(rr, 1) # take one randomly
        }
      }
      return(list(rr = rr, sw = sw))
    }

    ######################################################
    # Start FSPmix implementation

    # prep variables to store output
    rownames(dat)<- ppl<- 1:dim(dat)[1]
    THRESHOLD.METHOD <- 'intersect'

    SampDat_Store<- data.frame(dat = NA, boot.str = NA)
    Te_mu_store<- data.frame(Thresh = rep(NA, no.bootstrap),
                             mu1 = rep(NA, no.bootstrap),
                             mu2 = rep(NA, no.bootstrap))

    two.groups<- st.dev.T_e<- NA

    ##
    ## Conduct bootstrap
    for(j in 1:no.bootstrap){
      op<- find_MixtureThreshold_Simulation(dat = dat[,p],
                                            boot.size = boot.size, method=THRESHOLD.METHOD)

      count = 1
      sw = op$sw  # if sw = 0, solution found
      sw
      while(sw ==  1 | count == 10){ # <-- means above op threw an error
        op<- find_MixtureThreshold_Simulation(dat = dat[,p],
                                              boot.size = boot.size, method=THRESHOLD.METHOD)
        count = count + 1
        sw = op$sw
      }

      Te_mu_store$Thresh[j]<- op$mix.threshold
      Te_mu_store$mu1[j]<- op$mix.means[1]
      Te_mu_store$mu2[j]<- op$mix.means[2]

      temp.d<- op$boot.samp
      temp.d$boot.str<- rep(paste("Boot.", j, sep = ""),boot.size )
      SampDat_Store = rbind(SampDat_Store, temp.d)
    }

    ##
    ## Criterion to determine if there are two groups in the data
    mean.T_e<- mean(Te_mu_store$Thresh)
    mean.T_e

    sd.T_e<- sd(Te_mu_store$Thresh)
    sd.T_e

    mean.mu<- apply(Te_mu_store[, 2:3], 2, mean)
    mean.mu
    #st.dev.T_e[p]<- sd.T_e
    #all.mu_Store[[p]]<- mu_Store
    #all.mu_summary[[p]]<- as.data.frame(apply(Te_mu_store[, 2:3], 2, summary))

    SampDat_Store<- SampDat_Store[-1,]
    rownames(SampDat_Store)<- 1:dim(SampDat_Store)[1]

    # two groups found?
    if(mean.mu[1] < (mean.T_e - sd.T_e) & (mean.T_e + sd.T_e) < mean.mu[2]){
      two.groups<- TRUE
      interval.T_e<- c(mean.T_e - sd.T_e, mean.T_e + sd.T_e)
    }else{
      two.groups<- FALSE
      interval.T_e<- c(NA,NA)
    }

    interval.T_e

    # determined by the range of T_e
    #if(mean(mu_Store[,1]) < rangeT_e[1] & rangeT_e[2] < mean(mu_Store[,2])){
    #  two.groups[p]<- TRUE
    #}else(two.groups[p]<- FALSE)


    #### plot genes  --------------------------------------------
    SampDat_Store$boot.str<- factor(SampDat_Store$boot.str)

    #x11()
    p.Feature<- ggplot(SampDat_Store, aes(x = dat, group = boot.str)) +
      #geom_density(alpha = 0.1) +
      geom_density(colour = fade("black",20))+
      ggtitle(paste("Feature ", p,sep = "") ) +
      geom_vline(xintercept = interval.T_e, colour = "blue", size=1) +
      geom_vline(xintercept = mean(mean.mu[1]), colour = "gray55", size=1) +
      geom_vline(xintercept = mean(mean.mu[2]), colour = "gray55", size=1) +
      #annotate('text', x = Inf, y = Inf, hjust = 1.2, vjust = 2,
      #         label = paste("Thresh: ",round(interval.T_e[1],2), ",",
      #                       round(interval.T_e[2],2), sep = "") ,
      #         size=3.5, colour = "blue") +
      annotate('text', x = Inf, y = Inf, hjust = 1, vjust = 1.5,
               label = paste("Thresh: (", round(interval.T_e[1],2), ",",
                             round(interval.T_e[2],2), ")", sep = "") ,
               size=3, colour = "blue")+
      theme_bw() +
      annotate('text', x = Inf, y = Inf, hjust = 1, vjust = 3,
               label = paste("mu1: ", round(mean.mu[1],2), sep = "") ,
               size=3, colour = "gray55") +
      annotate('text', x = Inf, y = Inf, hjust = 1, vjust = 4.5,
               label = paste("mu2: ", round(mean.mu[2],2), sep = "") ,
               size=3, colour = "gray55") +
      theme(text = element_text(size=12))+
      xlab("Feature value")

    #x11()
    #p.Feature

    ##
    ## Identify the group (A/B) for those individuals who the algorithm
    ## identified two groups over all genes

    if(two.groups){
      sub.d<- data.frame(Feature = dat[,p], ppl = ppl)

      sub.d<- mutate(sub.d,
                     Pred = factor(ifelse(Feature < interval.T_e[1], "Pred.A",
                                          ifelse(Feature > interval.T_e[2], "Pred.B", "Pred.C")),
                                   levels = c("Pred.A", "Pred.B", "Pred.C")),
                     id = rep(colnames(dat)[p], dim(sub.d)[1]))

      sub.d

      # How many times did the algorithm detect 2 groups?
      summ.op<-data.frame(Feature.no = p, two.groups=two.groups,
                          mean.Te = round(mean.T_e,2),
                          sd.Te = round(sd.T_e,2),
                          mean.mu1 = round(mean.mu[1],2),
                          mean.mu2 = round(mean.mu[2],2) )
    }

    #summ.op

    ###
    if(two.groups){
      op<- list(Classification.Pred=sub.d,
                two.groups=two.groups,
                summ.op = summ.op,
                SampDat_Store=SampDat_Store,
                Plot = p.Feature)
    }else{
      op<- list(two.groups=two.groups,
                summ.op = NULL,
                SampDat_Store=SampDat_Store,
                Plot = p.Feature)
    }


  }



  stopCluster(cl)
  return(final.op)
}


