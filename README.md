# FSPmix
 Fast unsupervised feature selection using mixture models
 
## Install FSPmix

```{r}
library(devtools)
library(roxygen2)
install_github('MarcelaCespedes/FSPmix')
library(FSPmix)
```

## Simple example

```{r}
library(ggplot2)

source("SimulationStudy1.r")
source("FSPmix_Sim.r")
source("spec_sens.r")
source("multiplot.r")  # <-- once package is loaded am not sure if I should include all these source() calls
```

Simulation study 1: simulate 20 gene data

```{r}
no.ppl = 500
op<- SimulationStudy1(sEEd = 948575, no.ppl = no.ppl)
#x11()
#op$plot.ColourCoded.genes
#op$plot.WhatFSPmixSees

dat<- op$dat
head(dat)
dim(dat)
```

Set up simulation variables for FSPmix_Sim
```{r}
feature.dat<- op$dat[, 1:20]
class<- op$dat$group
boot.size<- round(no.ppl*0.8) # size of boot size sample
no.bootstrap<- 100 #500 # as per paper - this takes a while

sim.op<- FSPmix_Sim(boot.size = boot.size,
                    no.bootstrap = no.bootstrap,
                    dat = feature.dat, 
                    class = class)

length(sim.op)
sim.op$summ.op # <-- view summary of output
```

View simulation results
```{r}
x11()
multiplot(plotlist = sim.op$all.plots, cols = 5)
```

Assess prediction performance
```{r}
length(sim.op$all.pred.st)
# Summary of predictions for each simulated Gene
summary.pred <- sapply(sim.op$all.pred.st, function (x) table(x$Pred))
colnames(summary.pred)<- paste("Gene.",1:20, sep = "")
summary.pred

pred.A<- pred.B<- pred.C<- matrix(NA, nrow = 1, ncol = 4)
colnames(pred.A)<- colnames(pred.B)<- colnames(pred.C)<- c("Gene", "ppl", "Pred", "id")

Specificity.Sensitivity.Res<- list()
Fig1.dat<- data.frame(Feature = paste("Feature.", 1:20),
                      Specificity = rep(NA, 20),
                      Sensitivity = rep(NA, 20))

for(i in 1:20){
  temp.d<- sim.op$all.pred.st[[i]]
  
  tA<- subset(temp.d, Pred == "Pred.A")
  tB<- subset(temp.d, Pred == "Pred.B")
  tC<- subset(temp.d, Pred == "Pred.C")
  
  ## Get their spcificity/ sensitivity
  pred.prot.i = rbind(tA, tB)
  temp.dat<- subset(dat, ppl %in% pred.prot.i$ppl) # need to remove all those in limbo
  
  performance = spec_sens(real.dat = temp.dat, pred.dat = pred.prot.i)
  Specificity.Sensitivity.Res[[i]] <- performance$conf.mat.ratio
  Fig1.dat$Sensitivity[i]<- performance$conf.mat.ratio[1,1]
  Fig1.dat$Specificity[i]<- performance$conf.mat.ratio[2,2]
  ###
  pred.A<- rbind(pred.A, tA)
  pred.B<- rbind(pred.B, tB)
  pred.C<- rbind(pred.C, tC)
}

pred.A<- pred.A[-1,]
pred.B<- pred.B[-1,]
pred.C<- pred.C[-1,]
```

View specificity and sensitivity output for all simulated genes
```{r}
names(Specificity.Sensitivity.Res)<- paste("Gene.",1:20, sep = "")
Specificity.Sensitivity.Res
```

Replicate Figure 1 in manuscript
```{r}
Fig1.dat<- melt(Fig1.dat, id.vars = "Feature")
Fig1.dat$Feature<- factor(Fig1.dat$Feature, levels = unique(Fig1.dat$Feature) )
Fig1.dat$Performance<- factor(Fig1.dat$Performance, levels = c("Sensitivity", "Specificity"))
colnames(Fig1.dat)[2]<- "Performance"

x11()
ggplot(Fig1.dat, aes(x = Feature, y = value, colour = Performance)) +
  geom_point(size=4) +
  geom_hline(yintercept = 0.9) +
  theme_bw() +
  xlab("") + ylab("Performance value (Sensitivity/Specificity)")+
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 70,hjust = 1, size = 14) ) 
```
