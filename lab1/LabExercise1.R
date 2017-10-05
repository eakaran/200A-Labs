# Elizabeth Karan
# Lab Exercise 1: Introduction to estimating diversification rates

# 1. Lineage through time plots

library(phytools)

setwd("~/Dropbox/200A_Evolution/Labs/lab1")
# load tree
darter.tree <- read.tree("etheostoma_percina_chrono.tre")
# plot tree
plotTree(darter.tree,ftype="i",fsize=0.4,type="fan",lwd=1)
# generate LTT plot
obj<-ltt(darter.tree,log.lineages=FALSE)

# check tree
darter.tree
# a fully birfucating tree with N tips will have N-1 internal nodes
is.binary(darter.tree)
# comes back false. we must fix this
darter.tree<-multi2di(darter.tree)
darter.tree
# now it has 201 tips and 200 nodes
is.binary(darter.tree)
# we have resolved the inner nodes

# create a new object storing tree data
obj<-ltt(darter.tree,plot=FALSE)
obj
# make a new LTT plot
plot(obj,log.lineages=FALSE,main="LTT plot for darters")
# plot tree and LTT plot together
plot(obj,show.tree=TRUE,log.lineages=FALSE,main="LTT plot for darters")
# number of linages rises exponentially through time with the pure birth model
# we want to show a linear relationship
# plot our pure-birth trees on a semi-log scale
# vertical axis is on a log-scale and horizontal axis is on a linear scale.
plot(obj,log.lineages=FALSE,log="y",main="LTT plot for darters",
     ylim=c(2,Ntip(darter.tree)))
## we can overlay the pure-birth prediction:
h<-max(nodeHeights(darter.tree))
x<-seq(0,h,by=h/100)
b<-(log(Ntip(darter.tree))-log(2))/h
lines(x,2*exp(b*x),col="red",lty="dashed",lwd=2)

# compare observed LTT to simulated LTTs assuming a pure-birth process of the 
# same duration & resulting in the same total number of species
# use phytools function pbtree
trees<-pbtree(b=b,n=Ntip(darter.tree),t=h,nsim=100,method="direct",
              quiet=TRUE)
obj<-ltt(trees,plot=FALSE)
plot(obj,col="grey",main="LTT of darters compared to simulated LTTs")
lines(c(0,h),log(c(2,Ntip(darter.tree))),lty="dashed",lwd=2,col="red")
## now let's overlay our original tree
ltt(darter.tree,add=TRUE,lwd=2)

# ltt95 gives a (1-α)% CI for the LTT based on a set of trees
ltt95(trees,log=TRUE)
title(main="LTT of darters compared to simulated LTTs")
ltt(darter.tree,add=TRUE,log.lineages=FALSE,col="red",lwd=2)


# 2. Fitting pure birth and birth-death models to trees
library(phytools)
fitbd <- birthdeath(darter.tree)
fitbd
# use a phytools function bd() to get b and d rate separately
bd(fitbd)


# 3. The γ-statistic

# γ-statistic has astandard normal distribution for trees generated under a pure-birth speciation model
# test this with our 100 simulated pure-birth phylogenies
g<-sapply(trees,function(x) ltt(x,plot=FALSE)$gamma)
hist(g,main=expression(paste("Distribution of ",gamma," from simulation")))
mean(g)
var(g)

# test hypotheses about γ automatically with ltt
obj<-ltt(darter.tree,plot=FALSE)
print(obj)

# simulate a tree under a different model of lineage accumulation - the coalescent - 
# and see what the result is
# this should result in a significantly positive γ
coal.tree<-rcoal(n=100)
plotTree(coal.tree,ftype="off")
# and it does
obj<-ltt(coal.tree,log.lineages=FALSE,log="y")
darter.gamma <- obj$gamma # ltt returns a gamma value as one of its elements

# compare the mean gamma value of 200 pb trees to gamma from the coalescent tree
trees<-pbtree(n=100,nsim=200,scale=max(nodeHeights(coal.tree)))
ltt95(trees,log=TRUE)
title(main="Simulated coalescent trees compared to pure-birth LTTs")
ltt(coal.tree,add=TRUE,log.lineages=FALSE,col="red",lwd=2,lty="dashed")


# 3. Incomplete sampling

# the authors only sampled 201 of 216 species
# run MCCR test
library(geiger)
age <- 25.91862
richness <- 216
darterbirth =  (log(richness) - log(2))/age
darterbirth

richness <- 216
missing <- 15
#this simulates gamma values when trees are undersampled.
#we will grow trees with n=34 and prune them down to 13 taxa

num_simulations<-200 #number of simulations
g1_null<-numeric(num_simulations) #g1_null will hold the simulated gamma values
for(i in 1:num_simulations) {
  sim.bdtree(darterbirth, d=0, stop = "taxa", n=richness)->sim_tree 
  drop.random(sim_tree, missing)->prune # prune down to the # of taxa in the phylogeny
  gammaStat(prune)->g1_null[i]
}
# create a histogram of the null distribution
hist(g1_null)

#arrow indicates where the observed gamma falls in the null you just generated
arrows(darter.gamma, 40, darter.gamma, 0, col="red", lwd=2) 

# Which of the null values are smaller (more negative) than the data?
smallerNull<-g1_null<=darter.gamma
# How many TRUEs are there?
count<-sum(smallerNull)

# finally, what is the p-value?
mccr_pval<-(count+1)/(num_simulations+1)
mccr_pval
