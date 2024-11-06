#read libraries
library("detectRUNS")
library("ggplot2")

#set working directory
setwd("/faststorage/project/Coregonus/Aja/Dragonfly/ROH/")

#read data
genotypeFilePath <- "/faststorage/project/Coregonus/Aja/Dragonfly/ROH/input_files/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.MappabilityMask.RepeatMask.miss1.ped"
mapFilePath <- "/faststorage/project/Coregonus/Aja/Dragonfly/ROH/input_files/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.MappabilityMask.RepeatMask.miss1.map"

### DETECT RUNS

slidingRuns <- slidingRUNS.run(
  genotypeFile = genotypeFilePath, 
  mapFile = mapFilePath, 
  minSNP = 116,
  threshold = 0.05,	
  windowSize = 116,
  maxOppWindow = 0, 	# maximum number of heterozygous SNPs allowed in each window
  maxMissWindow = 0,	# maximum number of missing SNPs allowed in each window
  maxGap = 280000, 	# the maximum distance between two consecutive SNPs to still be considered a potential ROH
  minLengthBps = 10000, 	# minimum length of a ROH in base pairs (bps)
  minDensity = 1/10000, 		# minimum acceptable number of SNPs per bps 
)

# save output
write.table(x=slidingRuns, file = "./run0/slidingRuns.txt", quote=F, row.names=F)

# read runs file
#slidingRuns=readExternalRuns(inputFile="./run0/slidingRuns.txt", program="detectRUNS")

#summary statistics
summaryList <- summaryRuns(
  runs = slidingRuns, mapFile = mapFilePath, genotypeFile = genotypeFilePath, 
  Class = 0.5, snpInRuns = TRUE)

#have a look at some of the 10 items in the returned summary list
write.csv(summaryList$summary_ROH_count_chr, file="./run0/summary_ROH_count_chr.csv")
write.csv(summaryList$summary_ROH_percentage_chr, file="./run0/summary_ROH_percentage_chr.csv")
write.csv(summaryList$summary_ROH_count, file="./run0/summary_ROH_count.csv")
write.csv(summaryList$summary_ROH_percentage, file="./run0/summary_ROH_percentage.csv")
write.csv(summaryList$summary_ROH_mean_chr, file="./run0/summary_ROH_mean_chr.csv")
write.csv(summaryList$summary_ROH_mean_class, file="./run0/summary_ROH_mean_class.csv")
write.csv(summaryList$result_Froh_genome_wide, file="./run0/result_Froh_genome_wide.csv")
write.csv(summaryList$result_Froh_chromosome_wide, file="./run0/result_Froh_chromosome_wide.csv")
write.csv(summaryList$result_Froh_class, file="./run0/result_Froh_class.csv")

#produce basic plot of all runs detected in each individual against their position along the chromosome
#plot_Runs(runs = slidingRuns, savePlots = TRUE, outputName = "./run0/Runs")

### F_ROH: ROH-based inbreeding

#calculate individual inbreeding/consanguinity
#chromosome by chromosome or on all genome depending on genome_wide()
Froh_inbreeding(slidingRuns, mapFilePath, genome_wide=TRUE)

#plot inbreeding levels 
#remember to choose style (one type or ALL)
plot_InbreedingChr(slidingRuns, mapFilePath, groupSplit = TRUE,
                   style =  c("All"),
                   outputName = "./run0/InbreedingChr", plotTitle = NULL, savePlots = TRUE)

#produce plot of all runs detected in each individual
library(ggplot2)
library(stringr)
require(gridExtra)
library(dplyr)
library(tidyr)

#read popmap
popmap<-read.table("/faststorage/project/Coregonus/Aja/Dragonfly/PCA_LD-thinned/popmap2_region.txt",sep="\t",stringsAsFactors=F)

# Create a list of vectors - one for each pop
vector_list <- lapply(split(popmap$V1, popmap$V2), function(x) as.character(x))

#create rows of data to make sure entire chromosomes are plotted regardless of presence of data points
poplist<-c("NWJ","EJ","WJ","SJ","NEZ","SZ","B")

for (pop in poplist){
individuals <- vector_list[[pop]]  #assign vector from list to variable "files"

for (ind in individuals){

row1<-c(pop,ind,"Contig9961","1","0","1","1")
row2<-c(pop,ind,"Contig94","1","0","1","1")
row3<-c(pop,ind,"Contig6623","1","0","1","1")
row4<-c(pop,ind,"Contig6610","1","0","1","1")
row5<-c(pop,ind,"Contig6679","1","0","1","1")
row6<-c(pop,ind,"Contig8990","1","0","1","1")
row7<-c(pop,ind,"Contig8857","1","0","1","1")
row8<-c(pop,ind,"Contig6617","1","0","1","1")
row9<-c(pop,ind,"Contig6697","1","0","1","1")
row10<-c(pop,ind,"Contig8852","1","0","1","1")
row11<-c(pop,ind,"Contig6842","1","0","1","1")
row12<-c(pop,ind,"Contig389","1","0","1","1")
row13<-c(pop,ind,"Contig6654","1","0","1","1")
row14<-c(pop,ind,"Contig9963","1","0","1","1")
row15<-c(pop,ind,"Contig6827","1","0","1","1")
row16<-c(pop,ind,"Contig9967","1","0","1","1")
row17<-c(pop,ind,"Contig583","1","0","1","1")
row18<-c(pop,ind,"Contig8991","1","0","1","1")
row19<-c(pop,ind,"Contig725","1","0","1","1")
row20<-c(pop,ind,"Contig49","1","0","1","1")
row21<-c(pop,ind,"Contig6776","1","0","1","1")
row22<-c(pop,ind,"Contig6930","1","0","1","1")
row23<-c(pop,ind,"Contig6795","1","0","1","1")
row24<-c(pop,ind,"Contig6708","1","0","1","1")
row25<-c(pop,ind,"Contig9098","1","0","1","1")
row26<-c(pop,ind,"Contig6692","1","0","1","1")
row27<-c(pop,ind,"Contig6845","1","0","1","1")
row28<-c(pop,ind,"Contig6876","1","0","1","1")
row29<-c(pop,ind,"Contig6941","1","0","1","1")
row30<-c(pop,ind,"Contig39","1","0","1","1")
row31<-c(pop,ind,"Contig6806","1","0","1","1")
row32<-c(pop,ind,"Contig6612","1","0","1","1")
row33<-c(pop,ind,"Contig6811","1","0","1","1")
row34<-c(pop,ind,"Contig6922","1","0","1","1")
row35<-c(pop,ind,"Contig1362","1","0","1","1")
row36<-c(pop,ind,"Contig6785","1","0","1","1")
row37<-c(pop,ind,"Contig6852","1","0","1","1")
row38<-c(pop,ind,"Contig6670","1","0","1","1")
row39<-c(pop,ind,"Contig6656","1","0","1","1")
row40<-c(pop,ind,"Contig7404","1","0","1","1")
row41<-c(pop,ind,"Contig6627","1","0","1","1")
row42<-c(pop,ind,"Contig7027","1","0","1","1")
row43<-c(pop,ind,"Contig7279","1","0","1","1")
row44<-c(pop,ind,"Contig7114","1","0","1","1")
row45<-c(pop,ind,"Contig8850","1","0","1","1")
row46<-c(pop,ind,"Contig7270","1","0","1","1")
row47<-c(pop,ind,"Contig7472","1","0","1","1")
row48<-c(pop,ind,"Contig6666","1","0","1","1")
row49<-c(pop,ind,"Contig6999","1","0","1","1")
row50<-c(pop,ind,"Contig7255","1","0","1","1")
row51<-c(pop,ind,"Contig7252","1","0","1","1")
row52<-c(pop,ind,"Contig6699","1","0","1","1")
row53<-c(pop,ind,"Contig6830","1","0","1","1")
row54<-c(pop,ind,"Contig7264","1","0","1","1")
row55<-c(pop,ind,"Contig1835","1","0","1","1")
row56<-c(pop,ind,"Contig230","1","0","1","1")
row57<-c(pop,ind,"Contig657","1","0","1","1")
row58<-c(pop,ind,"Contig7578","1","0","1","1")
row59<-c(pop,ind,"Contig6735","1","0","1","1")
row60<-c(pop,ind,"Contig6801","1","0","1","1")
row61<-c(pop,ind,"Contig7361","1","0","1","1")
row62<-c(pop,ind,"Contig9973","1","0","1","1")
row63<-c(pop,ind,"Contig9984","1","0","1","1")
row64<-c(pop,ind,"Contig6747","1","0","1","1")
row65<-c(pop,ind,"Contig2026","1","0","1","1")
row66<-c(pop,ind,"Contig6873","1","0","1","1")
row67<-c(pop,ind,"Contig7030","1","0","1","1")
row68<-c(pop,ind,"Contig7193","1","0","1","1")
row69<-c(pop,ind,"Contig7396","1","0","1","1")
row70<-c(pop,ind,"Contig9976","1","0","1","1")
row71<-c(pop,ind,"Contig7083","1","0","1","1")
row72<-c(pop,ind,"Contig6787","1","0","1","1")
row73<-c(pop,ind,"Contig7035","1","0","1","1")
row74<-c(pop,ind,"Contig8812","1","0","1","1")
row75<-c(pop,ind,"Contig1279","1","0","1","1")
row76<-c(pop,ind,"Contig7648","1","0","1","1")
row77<-c(pop,ind,"Contig7228","1","0","1","1")
row78<-c(pop,ind,"Contig8807","1","0","1","1")
row79<-c(pop,ind,"Contig2357","1","0","1","1")
row80<-c(pop,ind,"Contig9988","1","0","1","1")
row81<-c(pop,ind,"Contig2368","1","0","1","1")
row82<-c(pop,ind,"Contig7010","1","0","1","1")
row83<-c(pop,ind,"Contig7638","1","0","1","1")
row84<-c(pop,ind,"Contig1257","1","0","1","1")
row85<-c(pop,ind,"Contig8847","1","0","1","1")
row86<-c(pop,ind,"Contig6837","1","0","1","1")
row87<-c(pop,ind,"Contig6921","1","0","1","1")
row88<-c(pop,ind,"Contig6668","1","0","1","1")
row89<-c(pop,ind,"Contig7510","1","0","1","1")
row90<-c(pop,ind,"Contig6597","1","0","1","1")
row91<-c(pop,ind,"Contig2289","1","0","1","1")
row92<-c(pop,ind,"Contig7103","1","0","1","1")
row93<-c(pop,ind,"Contig6997","1","0","1","1")
row94<-c(pop,ind,"Contig7052","1","0","1","1")
row95<-c(pop,ind,"Contig8809","1","0","1","1")
row96<-c(pop,ind,"Contig9996","1","0","1","1")
row97<-c(pop,ind,"Contig6847","1","0","1","1")
row98<-c(pop,ind,"Contig7111","1","0","1","1")
row99<-c(pop,ind,"Contig7078","1","0","1","1")
row100<-c(pop,ind,"Contig7299","1","0","1","1")
row101<-c(pop,ind,"9961","1","20551762","20551763","1")
row102<-c(pop,ind,"94","1","14767039","14767040","1")
row103<-c(pop,ind,"6623","1","14361220","14361221","1")
row104<-c(pop,ind,"6610","1","13784832","13784833","1")
row105<-c(pop,ind,"6679","1","13435427","13435428","1")
row106<-c(pop,ind,"8990","1","12344626","12344627","1")
row107<-c(pop,ind,"8857","1","11989409","11989410","1")
row108<-c(pop,ind,"6617","1","11246234","11246235","1")
row109<-c(pop,ind,"6697","1","11209931","11209932","1")
row110<-c(pop,ind,"8852","1","10803280","10803281","1")
row111<-c(pop,ind,"6842","1","9848567","9848568","1")
row112<-c(pop,ind,"389","1","9686657","9686658","1")
row113<-c(pop,ind,"6654","1","8795943","8795944","1")
row114<-c(pop,ind,"9963","1","8657779","8657780","1")
row115<-c(pop,ind,"6827","1","8573604","8573605","1")
row116<-c(pop,ind,"9967","1","8329309","8329310","1")
row117<-c(pop,ind,"583","1","7810297","7810298","1")
row118<-c(pop,ind,"8991","1","7297351","7297352","1")
row119<-c(pop,ind,"725","1","6893062","6893063","1")
row120<-c(pop,ind,"49","1","6861040","6861041","1")
row121<-c(pop,ind,"6776","1","6785176","6785177","1")
row122<-c(pop,ind,"6930","1","6716959","6716960","1")
row123<-c(pop,ind,"6795","1","6551247","6551248","1")
row124<-c(pop,ind,"6708","1","6254755","6254756","1")
row125<-c(pop,ind,"9098","1","6165118","6165119","1")
row126<-c(pop,ind,"6692","1","6159818","6159819","1")
row127<-c(pop,ind,"6845","1","6109702","6109703","1")
row128<-c(pop,ind,"6876","1","6060377","6060378","1")
row129<-c(pop,ind,"6941","1","5889816","5889817","1")
row130<-c(pop,ind,"39","1","5759186","5759187","1")
row131<-c(pop,ind,"6806","1","5565633","5565634","1")
row132<-c(pop,ind,"6612","1","5341650","5341651","1")
row133<-c(pop,ind,"6811","1","5330790","5330791","1")
row134<-c(pop,ind,"6922","1","5051912","5051913","1")
row135<-c(pop,ind,"1362","1","5051128","5051129","1")
row136<-c(pop,ind,"6785","1","5018331","5018332","1")
row137<-c(pop,ind,"6852","1","4968962","4968963","1")
row138<-c(pop,ind,"6670","1","4927287","4927288","1")
row139<-c(pop,ind,"6656","1","4879140","4879141","1")
row140<-c(pop,ind,"7404","1","4846614","4846615","1")
row141<-c(pop,ind,"6627","1","4815769","4815770","1")
row142<-c(pop,ind,"7027","1","4789712","4789713","1")
row143<-c(pop,ind,"7279","1","4715109","4715110","1")
row144<-c(pop,ind,"7114","1","4706046","4706047","1")
row145<-c(pop,ind,"8850","1","4705212","4705213","1")
row146<-c(pop,ind,"7270","1","4694112","4694113","1")
row147<-c(pop,ind,"7472","1","4675307","4675308","1")
row148<-c(pop,ind,"6666","1","4574083","4574084","1")
row149<-c(pop,ind,"6999","1","4531221","4531222","1")
row150<-c(pop,ind,"7255","1","4530421","4530422","1")
row151<-c(pop,ind,"7252","1","4450808","4450809","1")
row152<-c(pop,ind,"6699","1","4440508","4440509","1")
row153<-c(pop,ind,"6830","1","4382150","4382151","1")
row154<-c(pop,ind,"7264","1","4336239","4336240","1")
row155<-c(pop,ind,"1835","1","4306919","4306920","1")
row156<-c(pop,ind,"230","1","4302506","4302507","1")
row157<-c(pop,ind,"657","1","4285162","4285163","1")
row158<-c(pop,ind,"7578","1","4214734","4214735","1")
row159<-c(pop,ind,"6735","1","4102920","4102921","1")
row160<-c(pop,ind,"6801","1","4095184","4095185","1")
row161<-c(pop,ind,"7361","1","4093523","4093524","1")
row162<-c(pop,ind,"9973","1","4067942","4067943","1")
row163<-c(pop,ind,"9984","1","4048754","4048755","1")
row164<-c(pop,ind,"6747","1","4017073","4017074","1")
row165<-c(pop,ind,"2026","1","4005502","4005503","1")
row166<-c(pop,ind,"6873","1","3992084","3992085","1")
row167<-c(pop,ind,"7030","1","3973253","3973254","1")
row168<-c(pop,ind,"7193","1","3938723","3938724","1")
row169<-c(pop,ind,"7396","1","3937937","3937938","1")
row170<-c(pop,ind,"9976","1","3930512","3930513","1")
row171<-c(pop,ind,"7083","1","3856834","3856835","1")
row172<-c(pop,ind,"6787","1","3844804","3844805","1")
row173<-c(pop,ind,"7035","1","3771056","3771057","1")
row174<-c(pop,ind,"8812","1","3738770","3738771","1")
row175<-c(pop,ind,"1279","1","3729699","3729700","1")
row176<-c(pop,ind,"7648","1","3710199","3710200","1")
row177<-c(pop,ind,"7228","1","3696238","3696239","1")
row178<-c(pop,ind,"8807","1","3672132","3672133","1")
row179<-c(pop,ind,"2357","1","3668220","3668221","1")
row180<-c(pop,ind,"9988","1","3667233","3667234","1")
row181<-c(pop,ind,"2368","1","3661826","3661827","1")
row182<-c(pop,ind,"7010","1","3648646","3648647","1")
row183<-c(pop,ind,"7638","1","3647478","3647479","1")
row184<-c(pop,ind,"1257","1","3637579","3637580","1")
row185<-c(pop,ind,"8847","1","3606059","3606060","1")
row186<-c(pop,ind,"6837","1","3585165","3585166","1")
row187<-c(pop,ind,"6921","1","3552548","3552549","1")
row188<-c(pop,ind,"6668","1","3544806","3544807","1")
row189<-c(pop,ind,"7510","1","3468267","3468268","1")
row190<-c(pop,ind,"6597","1","3447950","3447951","1")
row191<-c(pop,ind,"2289","1","3399061","3399062","1")
row192<-c(pop,ind,"7103","1","3362221","3362222","1")
row193<-c(pop,ind,"6997","1","3340223","3340224","1")
row194<-c(pop,ind,"7052","1","3339609","3339610","1")
row195<-c(pop,ind,"8809","1","3323069","3323070","1")
row196<-c(pop,ind,"9996","1","3320920","3320921","1")
row197<-c(pop,ind,"6847","1","3311793","3311794","1")
row198<-c(pop,ind,"7111","1","3301879","3301880","1")
row199<-c(pop,ind,"7078","1","3237682","3237683","1")
row200<-c(pop,ind,"7299","1","3214576","3214577","1")

runs <- rbind(slidingRuns, row1, row2, row3, row4, row5, row6, row7, row8, row9, row10, 
              row11, row12, row13, row14, row15, row16, row17, row18, row19, row20, 
              row21, row22, row23, row24, row25, row26, row27, row28, row29, row30, 
              row31, row32, row33, row34, row35, row36, row37, row38, row39, row40, 
              row41, row42, row43, row44, row45, row46, row47, row48, row49, row50, 
              row51, row52, row53, row54, row55, row56, row57, row58, row59, row60, 
              row61, row62, row63, row64, row65, row66, row67, row68, row69, row70, 
              row71, row72, row73, row74, row75, row76, row77, row78, row79, row80, 
              row81, row82, row83, row84, row85, row86, row87, row88, row89, row90, 
              row91, row92, row93, row94, row95, row96, row97, row98, row99, row100, 
              row101, row102, row103, row104, row105, row106, row107, row108, row109, row110, 
              row111, row112, row113, row114, row115, row116, row117, row118, row119, row120, 
              row121, row122, row123, row124, row125, row126, row127, row128, row129, row130, 
              row131, row132, row133, row134, row135, row136, row137, row138, row139, row140, 
              row141, row142, row143, row144, row145, row146, row147, row148, row149, row150, 
              row151, row152, row153, row154, row155, row156, row157, row158, row159, row160, 
              row161, row162, row163, row164, row165, row166, row167, row168, row169, row170, 
              row171, row172, row173, row174, row175, row176, row177, row178, row179, row180, 
              row181, row182, row183, row184, row185, row186, row187, row188, row189, row190, 
              row191, row192, row193, row194, row195, row196, row197, row198, row199, row200)
}
}

runs$chrom<-as.character(runs$chrom)
runs$from<-as.integer(runs$from)
runs$to<-as.integer(runs$to)

source("/faststorage/project/Coregonus/Aja/Dragonfly/ROH/plot_runs_custom.R")

#plot final figure
plot_Runs_custom(runs = runs, savePlots = TRUE, outputName = "./run0/custom_Runs")





