#' FSPmix_Sim.r implements the FSPmix algorithm with simulated data (known solution)
#' Input
#' dat: data.frame of features for FSPmix to search. Each column represents a feature and each row denotes an observation
#' class: vector with classification A, B or C solution for each feature
#' boot.size: positive integer, size of boot strap sample
#' no.bootstrap: positive integer, number of times to bootstrap

FSPmix_Sim<- function(dat, class,
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

  #source("find_MixtureThreshold_Simulation.r")

  library(grDevices)
  fade <- function(colors,alpha) {  # <-- to plot simulation densities
    rgbcols <- col2rgb(colors)
    rgb(rgbcols[1,],rgbcols[2,],rgbcols[3,],alpha/100*255,max=255)
  }

  # prep variables to store output
  rownames(dat)<- ppl<- 1:dim(dat)[1]
  no.genes<- dim(dat)[2]
  THRESHOLD.METHOD <- 'intersect'

  all.Te_mu_Store<- list()
  all.SampDat_Store<- all.mu_summary<- all.pred.st<- all.plots<- list()
  two.groups<- st.dev.T_e<- rep(NA, no.genes)

  pb <- txtProgressBar(min = 0, max = no.genes, style = 3)
  for(p in 1:no.genes){

    SampDat_Store<- data.frame(dat = NA, gRoup = NA, boot.str = NA)
    Te_mu_store<- data.frame(Thresh = rep(NA, no.bootstrap),
                             mu1 = rep(NA, no.bootstrap),
                             mu2 = rep(NA, no.bootstrap))

    ##
    ## Conduct bootstrap
    for(j in 1:no.bootstrap){
      op<- find_MixtureThreshold_Simulation(dat = dat[,p], gRoup =as.character(class),
                          boot.size = boot.size, method=THRESHOLD.METHOD)

      count = 1
      sw = op$sw  # if sw = 0, solution found
      sw
      while(sw ==  1 | count == 10){ # <-- means above op threw an error
        # try to find threshold again
        #op<- find_thresh_v3(dat = dat.2[,p], gRoup =as.character(dat.2$group),
        #                    boot.size = boot.size, method=THRESHOLD.METHOD)
        op<- find_MixtureThreshold_Simulation(dat = dat[,p], gRoup =as.character(class),
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

    all.Te_mu_Store[[p]]<- Te_mu_store
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
    all.SampDat_Store[[p]]<- SampDat_Store

    # two groups found?
    if(mean.mu[1] < (mean.T_e - sd.T_e) & (mean.T_e + sd.T_e) < mean.mu[2]){
      two.groups[p]<- TRUE
      interval.T_e<- c(mean.T_e - sd.T_e, mean.T_e + sd.T_e)
    }else{
      two.groups[p]<- FALSE
      interval.T_e<- c(NA,NA)
    }

    interval.T_e

    # determined by the range of T_e
    #if(mean(mu_Store[,1]) < rangeT_e[1] & rangeT_e[2] < mean(mu_Store[,2])){
    #  two.groups[p]<- TRUE
    #}else(two.groups[p]<- FALSE)


    #### plot genes  --------------------------------------------

    SampDat_Store$gRoup<- factor(SampDat_Store$gRoup, levels = c("A", "B"))
    SampDat_Store$boot.str<- factor(SampDat_Store$boot.str)

    #x11()
    p.prot<- ggplot(SampDat_Store, aes(x = dat, fill = gRoup,group = boot.str)) +
      #geom_density(alpha = 0.1) +
      geom_density(colour = fade("black",20))+
      ggtitle(paste("Sim Gene ", p,sep = "") ) +
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
      xlab("Simulated gene expression")


      #annotate('text', x = Inf, y = Inf, hjust = 2, vjust = 2,
      #         label = paste("mu means: ",round(mean.mu[1],2), ",",
     #                       round(mean.mu[2],2), sep = "") ,
     #          size=3.5, colour = 'blue')

    #x11()
    #p.prot
    all.plots[[p]]<- p.prot

    ##
    ## Identify the group (A/B) for those individuals who the algorithm
    ## identified two groups over all genes

    if(two.groups[p]){
      sub.d<- data.frame(Gene = dat[,p], ppl = ppl)

      sub.d<- mutate(sub.d,
                     Pred = factor(ifelse(Gene < interval.T_e[1], "Pred.A",
                                          ifelse(Gene > interval.T_e[2], "Pred.B", "Pred.C")),
                                   levels = c("Pred.A", "Pred.B", "Pred.C")),
                     id = rep(colnames(dat)[p], dim(sub.d)[1]))

      all.pred.st[[p]] <- sub.d

    }

    setTxtProgressBar(pb, p)
  }
  close(pb)

  # How many times did the algorithm detect 2 groups?
  summ.op<-data.frame(gene.no = 1:no.genes, two.groups=two.groups,
                      mean.Te = round(unlist(lapply(all.Te_mu_Store, function(x) mean(x$Thresh))),2),
                      sd.Te = round(unlist(lapply(all.Te_mu_Store, function(x) sd(x$Thresh))),2),
                      mean.mu1 = round(unlist(lapply(all.Te_mu_Store, function(x) mean(x$mu1))),2),
                      mean.mu2 = round(unlist(lapply(all.Te_mu_Store, function(x) mean(x$mu2))),2) )

  #summ.op

  ###
  op<- list(all.pred.st=all.pred.st,
            two.groups=two.groups,
            summ.op = summ.op,
            all.Te_mu_Store=all.Te_mu_Store,
            all.SampDat_Store=all.SampDat_Store,
            all.mu_summary=all.mu_summary,
            all.plots = all.plots)

  return(op)
}


