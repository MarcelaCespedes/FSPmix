#' Function to evaluate the specificity and sensitivity
#' of FSPmix based on simulated data
#' input:
#' real.dat: data.frame which consists of all individuals applied to find_MixtureThreshold_Simulation.r to have the following columns
#'      ppl: ID of participants
#'      Protein.1: data that was classified
#'      group: "A" and "B" classifications
#'
#' output
#' pred.dat: combined predicted data.frame - to have column names
#'      ppl: ID of participants
#'      Pred: two levels of classification - "Pred.A" and "Pred.B"
#'      Protein.1: data that was classified

spec_sens<- function(real.dat, pred.dat){

  real.dat<- real.dat[order(real.dat$ppl),]
  pred.dat<- pred.dat[order(pred.dat$ppl),]

  ## need to remove those individuals which had no predicted value
  dat.sol<- subset(real.dat, ppl %in% pred.dat$ppl)
  tot.obs<- dim(dat.sol)[1]
  tot.obs

  ## True positive: A
  sol.A<- subset(dat.sol, group == "A")
  tot.A = dim(sol.A)[1]
  tot.A

  p.A<- subset(pred.dat, Pred == "Pred.A")
  pred.A<- dim(p.A)[1]
  pred.A

  tp = length(intersect(sol.A$ppl, p.A$ppl))  # those individuals in A sol && prediction: correct classification
  tp  # true positive

  pred.A.sol.B = abs(pred.A - tp)  # any individual classed as A - but they are actually B: incorrect classification
  pred.A.sol.B  # Type I error

  ## for group B
  sol.B<- subset(dat.sol, group == "B")
  tot.B = dim(sol.B)[1]
  tot.B

  p.B<- subset(pred.dat, Pred == "Pred.B")
  pred.B<- dim(p.B)[1]
  pred.B

  tn<- length(intersect(sol.B$ppl, p.B$ppl))
  tn

  pred.B.sol.A = abs(pred.B - tn)
  pred.B.sol.A

  conf.mat<- data.frame(pred.A = c(NA, NA), pred.B = c(NA, NA) )
  conf.mat[1,1]<- tp
  conf.mat[2,2]<- tn
  conf.mat[1,2]<- pred.B.sol.A
  conf.mat[2,1]<- pred.A.sol.B

  row.names(conf.mat)<- c("sol.A", "sol.B")
  ########

  conf.mat.ratio<- conf.mat
  conf.mat.ratio[1,] <- conf.mat.ratio[1,]/tot.A
  conf.mat.ratio[2,] <- conf.mat.ratio[2,]/tot.B

  #########################################
  return(list(tot.obs = tot.obs, conf.mat = conf.mat, conf.mat.ratio = round(conf.mat.ratio,3)) )

}

