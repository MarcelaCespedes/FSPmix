#' Simulation Study 1 - simulation of 20 genes
#' Generate gene expression data with two latent groups (A and B)
#' TODO Intention is to modify this for a more realistic real-life scenario
#' input
#' sEED: random seed generator
#' no.ppl: Even interger - number of observations (or participants) to simulate
SimulationStudy1<- function(sEEd = 123, no.ppl = 1000){

  library(ggplot2)
  library(reshape2)
  library(dplyr)
  library(mixtools)

  set.seed(sEEd)
  #set.seed(789)  # NOTE: different seeds lead to different results

  no.genes = 20
  dat.2<- data.frame(V1 = rep(0, no.ppl))

  ## **********************************************************
  ## simulate 10 RIGHT skewed data
  no.gene.rs = 10

  from = -1.5; to= -1  # <-- group A
  mean1 <- seq(from = from, to=to, by = ((to - from)/(no.gene.rs - 1)))
  var1<- runif(no.gene.rs, min = 0.5, max =0.7)


  from = 3; to= 2  # <--- group B
  mean2<- seq(from = from, to=to, by=((to - from)/(no.gene.rs - 1)))
  from=5; to=10
  var2<- seq(from = from, to= to, by = ((to - from)/(no.gene.rs - 1)))

  dat.2<- data.frame(V1 = rep(NA, no.ppl))
  for(i in 1:no.gene.rs){
    dat.2[1:(no.ppl/2), i]<- rnorm(no.ppl/2, mean = mean1[i], sd = sqrt(var1[i]) )
    dat.2[(no.ppl/2 +1):no.ppl, i]<- rnorm(no.ppl/2, mean = mean2[i], sd = sqrt(var2[i]) )
  }

  dat.2r = mutate(dat.2, ppl = 1:no.ppl,
                  group = factor(c(rep("A", no.ppl/2), rep("B", no.ppl/2)), levels = c("A", "B")))

  gene.list<- paste("gene.", 1:no.gene.rs, sep = "")
  colnames(dat.2r)[1:no.gene.rs]<- gene.list

  dat.right.skewed<- dat.2r
  dat.2ar<- melt(dat.2r, id.vars = c("ppl", "group"))



  ## ******************************************************************
  ## simulate 10 LEFT skewed data
  no.gene.ls = 10

  from = 1.5; to= 1  # <-- group
  mean1<- seq(from = from, to=to, by=((to - from)/(no.gene.rs - 1)))
  var1<- runif(no.gene.rs, min = 0.5, max =0.7)

  from = -3; to= -2  # <-- group
  mean2 <- seq(from = from, to=to, by = ((to - from)/(no.gene.rs - 1)))
  from=5; to=10
  var2<- seq(from = from, to= to, by = ((to - from)/(no.gene.rs - 1)))


  dat.2<- data.frame(V1 = rep(NA, no.ppl))
  for(i in 1:no.gene.ls){
    #dat.2[1:(no.ppl/2), i]<- rnorm(no.ppl/2, mean = mean1[i], sd = sqrt(var1[i]) )
    #dat.2[(no.ppl/2 +1):no.ppl, i]<- rnorm(no.ppl/2, mean = mean2[i], sd = sqrt(var2[i]) )
    dat.2[(no.ppl/2 +1):no.ppl, i]<- rnorm(no.ppl/2, mean = mean1[i], sd = sqrt(var1[i]) )
    dat.2[1:(no.ppl/2), i]<- rnorm(no.ppl/2, mean = mean2[i], sd = sqrt(var2[i]) )
  }

  dat.2l = mutate(dat.2, ppl = 1:no.ppl,
                  group = factor(c(rep("A", no.ppl/2), rep("B", no.ppl/2)), levels = c("A", "B")))

  gene.list<- paste("gene.",(no.gene.ls+1):(2*no.gene.ls), sep = "")
  colnames(dat.2l)[1:no.gene.ls]<- gene.list

  dat.left.skewed<- dat.2l
  dat.2al<- melt(dat.2l, id.vars = c("ppl", "group"))


  ###
  ###
  ### Combine both data sets

  # check
  #all.equal(dat.right.skewed$group, dat.left.skewed$group)
  #all.equal(dat.right.skewed$ppl, dat.left.skewed$ppl)

  dat<- cbind(dat.right.skewed[, 1:10], dat.left.skewed)
  dat.2all<- melt(dat, id.vars = c("ppl", "group"))

  colnames(dat)[1:20]<- paste("Feature.", 1:20, sep = "")
  dat.2all<- melt(dat, id.vars = c("ppl", "group"))
  # check allocation of groups A and B
  p1<- ggplot(dat.2all, aes(x = value, fill = group)) + geom_density(alpha = 0.2) +
    facet_wrap(~variable, scales="free", ncol = 5) +
    theme_bw() + ggtitle("Simulated synthetic feature data: colour coded into groups") + xlab("Simulated normalised feature") +
    theme(legend.position="bottom") +
    theme(strip.text.x = element_text(size = 10))

  #x11()
  #p1

  p2<- ggplot(dat.2all, aes(x = value)) + geom_density(alpha = 0.2) +
    facet_wrap(~variable, scales="free", ncol = 10) +
    theme_bw() + ggtitle("Simulated left & right skewed data: what FSPmix sees") + xlab("Simulated normalised Gene expression")

  #x11()
  #p2

  #################################
  op<-list(dat = dat,
           plot.ColourCoded.genes = p1,
           plot.WhatFSPmixSees = p2)

  #################################
  return(op)
}
