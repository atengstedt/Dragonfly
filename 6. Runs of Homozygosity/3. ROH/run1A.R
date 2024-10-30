#read libraries
library("detectRUNS")
library("ggplot2")

#set working directory
setwd("/faststorage/project/Coregonus/Aja/Dragonfly/ROH/")

#read data
genotypeFilePath <- "/faststorage/project/Coregonus/Aja/Dragonfly/ROH/input_files/Dragonfly_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.MappabilityMask.RepeatMask.miss1.ped"
mapFilePath <- "/faststorage/project/Coregonus/Aja/Dragonfly/ROH/input_files/Dragonfly_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.MappabilityMask.RepeatMask.miss1.map"

### DETECT RUNS

slidingRuns <- slidingRUNS.run(
  genotypeFile = genotypeFilePath, 
  mapFile = mapFilePath, 
  minSNP = 100,
  threshold = 0.05,	
  windowSize = 100,
  maxOppWindow = 1, 	# maximum number of heterozygous SNPs allowed in each window
  maxMissWindow = 1,	# maximum number of missing SNPs allowed in each window
  maxGap = 280000, 	# the maximum distance between two consecutive SNPs to still be considered a potential ROH
  minLengthBps = 10000, 	# minimum length of a ROH in base pairs (bps)
  minDensity = 1/50000, 		# minimum acceptable number of SNPs per bps 
)

# save output
write.table(x=slidingRuns, file = "./run1A/slidingRuns.txt", quote=F, row.names=F)

# read runs file
#slidingRuns=readExternalRuns(inputFile="./run1A/slidingRuns.txt", program="detectRUNS")

#summary statistics
summaryList <- summaryRuns(
  runs = slidingRuns, mapFile = mapFilePath, genotypeFile = genotypeFilePath, 
  Class = 0.1, snpInRuns = TRUE)

#have a look at some of the 10 items in the returned summary list
write.csv(summaryList$summary_ROH_count_chr, file="./run1A/summary_ROH_count_chr.csv")
write.csv(summaryList$summary_ROH_percentage_chr, file="./run1A/summary_ROH_percentage_chr.csv")
write.csv(summaryList$summary_ROH_count, file="./run1A/summary_ROH_count.csv")
write.csv(summaryList$summary_ROH_percentage, file="./run1A/summary_ROH_percentage.csv")
write.csv(summaryList$summary_ROH_mean_chr, file="./run1A/summary_ROH_mean_chr.csv")
write.csv(summaryList$summary_ROH_mean_class, file="./run1A/summary_ROH_mean_class.csv")
write.csv(summaryList$result_Froh_genome_wide, file="./run1A/result_Froh_genome_wide.csv")
write.csv(summaryList$result_Froh_chromosome_wide, file="./run1A/result_Froh_chromosome_wide.csv")
write.csv(summaryList$result_Froh_class, file="./run1A/result_Froh_class.csv")

#produce basic plot of all runs detected in each individual against their position along the chromosome
plot_Runs(runs = slidingRuns, savePlots = TRUE, outputName = "./run1A/Runs")

### F_ROH: ROH-based inbreeding

#calculate individual inbreeding/consanguinity
#chromosome by chromosome or on all genome depending on genome_wide()
Froh_inbreeding(slidingRuns, mapFilePath, genome_wide=TRUE)

#plot inbreeding levels 
#remember to choose style (one type or ALL)
plot_InbreedingChr(slidingRuns, mapFilePath, groupSplit = TRUE,
                   style =  c("All"),
                   outputName = "./run1A/InbreedingChr", plotTitle = NULL, savePlots = TRUE)

