# Elizabeth Karan
# Lab Exercise 1: Introduction to estimating diversification rates

## FIRST: A TUTORIAL
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
par(mfrow=c(1,1), mar=c(4, 4, 3, 2)) 
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
par(mfrow=c(1,1), mar=c(4, 4, 3, 2))
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

## END OF TUTORIAL


# EXERCISES

# Calculate the gamma statistic for this phylogeny of homalopsid snakes from 
# Alfaro et al., 2008: snake.tre.

library(phytools)
setwd("~/Dropbox/200A_Evolution/Labs/lab1")
snake_phy <- read.tree("homalops.phy")

# look at tree
plotTree(snake_phy,ftype="i",fsize=0.4,type="fan",lwd=1)
plotTree(snake_phy,ftype="i",fsize=0.4,type="cladogram",lwd=1)
# I like to look at it this way
plotTree(snake_phy,ftype="i",fsize=0.4,type="phylogram",lwd=1)

# get an LTT plot
snakeLTT <- ltt(snake_phy, plot = FALSE); snakeLTT
# gives a "gamma" statistic of -3.2411
# A negative gamma statistic suggests that the rate of speciation was more rapid 
# earlier in the tree

# Given that the crown age of the snake radiation is 22 MY and the total richness of 
# the clade is 34 species, determine whether the observed gamma could be due to the 
# amount of incomplete sampling in the empirical tree. On the basis of the MCCR test 
# what can you conclude about the tempo of homalopsid snake diversification?
library(geiger)
age <- 22
# use a purebirth estimate of the the tree based upon their total age and richness
richness <- 34
missing <- 13
snakebirth = (log(richness) - log(2))/age

snake_gamma <- snakeLTT$gamma


#this simulates gamma values when trees are undersampled.
#we will grow trees with n=34 and prune them down to 13 taxa

num_simulations<-200 #number of simulations
g1_null<-numeric(num_simulations) #g1_null will hold the simulated gamma values
for(i in 1:num_simulations) {
  sim.bdtree(snakebirth, d=0, stop = "taxa", n=richness)->sim_tree 
  drop.random(sim_tree, missing)->prune # prune down to the # of taxa in the phylogeny
  gammaStat(prune)->g1_null[i]
}
# create a histogram of the null distribution
hist(g1_null)

#arrow indicates where the observed gamma falls in the null you just generated
arrows(snake_gamma, 40, snake_gamma, 0, col="red", lwd=2) 

# Which of the null values are smaller (more negative) than the data?
smallerNull<-g1_null<=snake_gamma
# How many TRUEs are there?
count<-sum(smallerNull); count

# finally, what is the p-value?
mccr_pval<-(count+1)/(num_simulations+1)
mccr_pval

# What is the birth rate and death rate of the homalopsid tree?
fitbd <- birthdeath(snake_phy)
bd(fitbd)


# Pomacentridae phylogeny with 216 tips (Frédérich, B. et al. 2013)
poma_tree <- read.nexus("damsels260miltime.tree.txt")
# take a look at it
poma_tree
plotTree(poma_tree,ftype="i",fsize=0.4,type="phylogram",lwd=1)
# remove the outgroups from tree
tip <- c("OUG_Centropyge_bicolor", "OUG_Cymatogaster_aggregata", "OUG_Cypho_purpurascens", 
         "OUG_Dicentrarchus_labrax", "OUG_Embiotoca_jacksoni", "OUG_Ptychochromis_oligacanthus", 
         "OUG_Semicossyphus_pulcher", "OUG_Thorichthys_meeki")
poma_only_tree <- drop.tip(poma_tree, tip)
# take a look at it
poma_only_tree
# tree with 208 tips only with species from family Pomacentridae
plotTree(poma_only_tree,ftype="i",fsize=0.4,type="phylogram",lwd=1)
# fit a birth-death model
poma_bd <- birthdeath(poma_only_tree)
# report b and b
bd(poma_bd)

# perform MCCR test

# get an LTT plot
pomaLTT <- ltt(poma_only_tree, plot = FALSE); pomaLTT
# age of clade in MY
age <- 50
# use a purebirth estimate of the the tree based upon their total age and richness
richness <- 386
missing <- 178
poma_birth = (log(richness) - log(2))/age
poma_gamma <- pomaLTT$gamma



num_simulations<-200 #number of simulations
gamma_null<-numeric(num_simulations) #g1_null will hold the simulated gamma values
for(i in 1:num_simulations) {
  sim.bdtree(poma_birth, d=0, stop = "taxa", n=richness)->sim_tree 
  drop.random(sim_tree, missing)->prune # prune down to the # of taxa in the phylogeny
  gammaStat(prune)->gamma_null[i]
}
# create a histogram of the null distribution
hist(gamma_null)

#arrow indicates where the observed gamma falls in the null you just generated
arrows(poma_gamma, 40, poma_gamma, 0, col="red", lwd=2) 

# Which of the null values are smaller (more negative) than the data?
smallerNull<-gamma_null<=poma_gamma
# How many TRUEs are there?
count<-sum(smallerNull); count

# finally, what is the p-value?
mccr_pval<-(count+1)/(num_simulations+1)
mccr_pval
