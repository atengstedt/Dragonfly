########################################################
############# (A) Test migration events ################
########################################################

library(OptM)
setwd("/faststorage/project/Coregonus/Aja/Dragonfly/TreeMix_LD-thinned/")

folder <- file.path(path="./noRoot2_k100_bs500_m1-10_rep10/test_migrations/")                     #path to files of TreeMix replicates with different migration edges (m) to test
test.linear = optM(folder, method = "linear", tsv="linear.txt")   #test m: produces a list of which the $out dataframe suggests optimum m based on multiple linear models
plot_optM(test.linear, method = "linear")                         #shows changes in log likelihood for different m, and suggests optimum m as 'change points
savePlot("noRoot_k100_bs500_m1-10_rep10/plot2_optM.png","png")  #save plot

#Choose optimum number of m and continue with step 3 in the TreeMix pipeline

########################################################
################## (B) Plot Tree #######################
########################################################

source("/faststorage/project/Coregonus/Aja/Coregonus/WGS/treemix/TreeMix_functions.R") #path to required functions for this analysis

## 1. From the final runs, compare tree likelihoods, select tree with highest likelihood, remove duplicates and retain tree(s) with unique topology. 
#Adapted from R functions written by Zecca, Labra and Grassi (2019).

setwd("/faststorage/project/Coregonus/Aja/Dragonfly/TreeMix_LD-thinned/noRoot_k100_bs500_m1-10_rep10/final_runs/")                            #folder with all TreeMix outputs from the final runs
maxLL("mytree_1m_finalrun_", nt=10)                                         #first argument is stem of TreeMix output files, nt = number of runs
                                                                  #shows ML trees and highest likelihood, as well as tree(s) with unique topology. Outputs "TreeLLs.txt" into workign directory

#If n of unique trees = 1, continue with step #2, if n > 1, you might want to create a consensus tree. 
#Note that bootstrap and migration values will not be available for consensus tree, thus you could also choose one ML tree 

cfTrees("mytree_1m_finalrun_", nt=10, p=1, m='PH85')                        #m is the method to calculate pairwise distances between unique trees (default is Robinson-Foulds distance)
                                                #p (number from 0.5 to 1)-proportion for a clade to be represented in the consensus tree. Default (1) is strict tree, 0.5 for majority-rule
                                                                  #plots consensus tree and saves "Consensus.newick" into working directory

## 2. Now plot and save unique tree with highest likelihood:

source("/faststorage/project/Coregonus/Aja/Coregonus/WGS/treemix/BITE/R/treemix_bootstrap.R")
source("/faststorage/project/Coregonus/Aja/Coregonus/WGS/treemix/BITE/R/plotting_funcs_bite.R")
source("/faststorage/project/Coregonus/Aja/Coregonus/WGS/treemix/BITE/R/newick_split.R")

png(file="../TreeMix.png", units="in", width = 5, height = 3.2, res=300)

par(mar=c(2.5,0.5,0.2,0.2)+0.1, mgp=c(1.2,0.1,0), xpd=TRUE,tcl=NA,cex=0.7)                                
treemix.bootstrap("mytree_1m_finalrun_1", out.file = "tmp",                 #stem of TreeMix files (as above) + number of run with highest likelihood and unique topology
                  phylip.file = "../mytree_constree.newick",           #consensus tree in newick format (from the bootstrap procedure generated with PHYLIP)    
                  nboot = 500, fill = TRUE,                       #nboot is the number of bootstraps used
                  pop.color.file = "../../pop_color.txt",                     #specify colours with a tab delimited pop.color.file - first column is pop name, second the colour
                  boot.legend.location = "topright", inset.vector=c(0,0),
                  plotmig=FALSE, mbar=FALSE,
                  cex=1, lwd=1.5,boot.cex=2, boot.legend.cex = 1)
dev.off()


#tree with 1 migration edge

png(file="../TreeMix_m1.png", units="cm", width = 7, height = 9.8, res=300)

par(mar=c(2.5,0.5,0.2,0.2)+0.1, mgp=c(1.2,0.1,0), xpd=TRUE,tcl=NA,cex=0.7)                                
treemix.bootstrap("mytree_1m_finalrun_1", out.file = "tmp",                 #stem of TreeMix files (as above) + number of run with highest likelihood and unique topology
                  phylip.file = "../mytree_constree.newick",           #consensus tree in newick format (from the bootstrap procedure generated with PHYLIP)    
                  nboot = 500, fill = TRUE,                       #nboot is the number of bootstraps used
                  pop.color.file = "../../pop_color.txt",                     #specify colours with a tab delimited pop.color.file - first column is pop name, second the colour
                  boot.legend.location = "topright", 
                  plotmig=TRUE, mbar=TRUE,
                  cex=0.8, lwd=1.5,boot.cex=2, boot.legend.cex = 1)
dev.off()

#tree with 2 migration edges

png(file="../TreeMix_m2.png", units="cm", width = 7, height = 9.8, res=300)

par(mar=c(2.5,0.5,0.2,0.2)+0.1, mgp=c(1.2,0.1,0), xpd=TRUE,tcl=NA,cex=0.7)                                
treemix.bootstrap("mytree_2m_finalrun_1", out.file = "tmp",                 #stem of TreeMix files (as above) + number of run with highest likelihood and unique topology
                  phylip.file = "../mytree_constree.newick",           #consensus tree in newick format (from the bootstrap procedure generated with PHYLIP)    
                  nboot = 500, fill = TRUE,                       #nboot is the number of bootstraps used
                  pop.color.file = "../../pop_color.txt",                     #specify colours with a tab delimited pop.color.file - first column is pop name, second the colour
                  boot.legend.location = "bottomright", 
                  plotmig=TRUE, mbar=TRUE,
                  cex=1, lwd=1.5,boot.cex=2, boot.legend.cex = 1)
dev.off()

