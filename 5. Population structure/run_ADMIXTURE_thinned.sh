##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

#-------------Generate ADMIXTURE input file
input="/faststorage/project/Coregonus/Aja/Dragonfly/VCF/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.miss1.LD-thinned.recode.vcf"
output="/faststorage/project/Coregonus/Aja/Dragonfly/ADMIXTURE_LD-thinned"

mkdir $output

sbatch -A Coregonus -t 12:00:00 --job-name plink --wrap\
 "plink --vcf ${input} --allow-extra-chr --make-bed --out ${output}/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.miss1.LD-thinned"
 
# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column with 0
file="/faststorage/project/Coregonus/Aja/Dragonfly/ADMIXTURE_LD-thinned/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.miss1.LD-thinned.bim"

awk '{$1="0";print $0}' ${file} > ${file}.tmp
#mv ${file}.tmp ${file}

#===========Run ADMIXTURE
folder="/faststorage/project/Coregonus/Aja/Dragonfly/ADMIXTURE_LD-thinned/"

cd ${folder}

for i in {2..10}
do
sbatch -A Coregonus -t 24:00:00 --mem 16G --job-name admixture --wrap\
 "admixture --cv=10 -B100 ${folder}/Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.miss1.LD-thinned.bed $i > log${i}.out"
done

#rm admixture-*.out

#----------Identify best number of K
grep -h CV log*.out>cross_validation2.txt

#plot cross validation (R)
library(stringr)
library(ggplot2)
library(dplyr)

cv <- read.table("cross_validation2.txt")

cv$K <-gsub("[\\(\\)]", "", regmatches(cv$V3, gregexpr("\\(.*?\\)", cv$V3)))
CV <- select(cv, V4,K)

colnames(CV) <- c("CV","K")

CV$K <- factor(CV$K, levels = c("K=2","K=3","K=4","K=5","K=6","K=7","K=8","K=9","K=10"))

png(file="Admixture_cross-validation2.png", width = 16.9, height =7, units="cm",res=200)

y_title="Cross-validation error"
graph_1<-ggplot(CV,aes(x=K,y=CV,group=1))
graph_1+geom_line()+
  labs(y=y_title)+labs(x=element_blank())+
  theme(axis.text.x=element_text(colour="black"))+
  theme(legend.title=element_blank())+
  theme(axis.text.y=element_text(colour="black",size=8))+
  theme(axis.text.x=element_text(colour="black",size=8))+
  theme(panel.border = element_rect(colour="black", fill=NA, size=1),
        axis.title=element_text(size=9,colour="black"))

dev.off()

#==========Plot ADMIXTURE (R)
library("tidyverse")
require(gridExtra)

samplelist <- read_table("/faststorage/project/Coregonus/Aja/Dragonfly/ADMIXTURE_LD-thinned/popmap2_region.txt",
                       col_names = c("sample","pop"))

all_data <- tibble(sample=character(),
                   pop=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())

for (k in 2:10){
  data <- read_delim(paste0("Dragonfly2_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.miss1.LD-thinned.",k,".Q"),
                  col_names = paste0("Q",seq(1:k)),
                  delim=" ")
  data$sample <- samplelist$sample
  data$pop <- samplelist$pop
  data$k <- k
  
  #This step converts from wide to long.
  data %>% gather(Q, value, -sample,-pop,-k) -> data
  all_data <- rbind(all_data,data)
}

#plot K2
plot1 <- all_data %>%
  filter(k == 2) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack",width=0.85) +
  facet_grid(~factor(pop, levels=c("NWJ","EJ","WJ","SJ","NEZ","SZ","B")),switch="x",scales="free", space="free_x") +
  theme_bw() +
  ylab("K=2") +
  theme(legend.position="none", 
    panel.border = element_rect(colour = "black", fill=NA, size=0.1), 
    panel.spacing.x = unit(0.2,"mm")) +
  theme(axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title = element_text(size=7),
    axis.ticks.x = element_blank(),
    plot.title = element_text(size=6, hjust=0.5, margin=margin(1,0,1,0)),
    strip.text.x = element_text(margin = margin(b = 0.1, t = 0.1)),
    axis.ticks.length=unit(0, "cm")) +
  theme(strip.background = element_blank(),
    strip.text.x = element_blank()) +
  theme(plot.margin = margin(0.04,0.02,0.04,0.04, "cm")) #+
#  scale_x_continuous(expand = c(0, 0)) + 
#  scale_y_continuous(expand = c(0, 0),breaks=c(0.0,0.5,1.0))

#plot K3
plot2 <- all_data %>%
  filter(k == 3) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack",width=0.85) +
  facet_grid(~factor(pop, levels=c("NWJ","EJ","WJ","SJ","NEZ","SZ","B")),switch="x",scales="free", space="free_x") +
  theme_bw() +
  ylab("K=3") +
  theme(legend.position="none", 
    panel.border = element_rect(colour = "black", fill=NA, size=0.1), 
    panel.spacing.x = unit(0.2,"mm")) +
  theme(axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title = element_text(size=7),
    axis.ticks.x = element_blank(),
    plot.title = element_text(size=6, hjust=0.5, margin=margin(1,0,1,0)),
    strip.text.x = element_text(margin = margin(b = 0.1, t = 0.1)),
    axis.ticks.length=unit(0, "cm")) +
  theme(strip.background = element_blank(),
    strip.text.x = element_blank()) +
  theme(plot.margin = margin(0.04,0.02,0.04,0.04, "cm")) #+
#  scale_x_continuous(expand = c(0, 0)) + 
#  scale_y_continuous(expand = c(0, 0),breaks=c(0.0,0.5,1.0))

#plot K4
plot3 <- all_data %>%
  filter(k == 4) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack",width=0.85) +
  facet_grid(~factor(pop, levels=c("NWJ","EJ","WJ","SJ","NEZ","SZ","B")),switch="x",scales="free", space="free_x") +
  theme_bw() +
  ylab("K=4") +
  theme(legend.position="none", 
    panel.border = element_rect(colour = "black", fill=NA, size=0.1), 
    panel.spacing.x = unit(0.2,"mm")) +
  theme(axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title = element_text(size=7),
    axis.ticks.x = element_blank(),
    plot.title = element_text(size=6, hjust=0.5, margin=margin(1,0,1,0)),
    strip.text.x = element_text(margin = margin(b = 0.1, t = 0.1)),
    axis.ticks.length=unit(0, "cm")) +
  theme(strip.background = element_blank(),
    strip.text.x = element_blank()) +
  theme(plot.margin = margin(0.04,0.02,0.04,0.04, "cm")) #+
#  scale_x_continuous(expand = c(0, 0)) + 
#  scale_y_continuous(expand = c(0, 0),breaks=c(0.0,0.5,1.0))

#plot K5
plot4 <- all_data %>%
  filter(k == 5) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack",width=0.85) +
  facet_grid(~factor(pop, levels=c("NWJ","EJ","WJ","SJ","NEZ","SZ","B")),switch="x",scales="free", space="free_x") +
  theme_bw() +
  ylab("K=5") +
  theme(legend.position="none", 
    panel.border = element_rect(colour = "black", fill=NA, size=0.1), 
    panel.spacing.x = unit(0.2,"mm")) +
  theme(axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title = element_text(size=7),
    axis.ticks.x = element_blank(),
    plot.title = element_text(size=6, hjust=0.5, margin=margin(1,0,1,0)),
    strip.text.x = element_text(margin = margin(b = 0.1, t = 0.1)),
    axis.ticks.length=unit(0, "cm")) +
  theme(strip.background = element_blank(),
    strip.text.x = element_blank()) +
  theme(plot.margin = margin(0.04,0.02,0.04,0.04, "cm")) #+
#  scale_x_continuous(expand = c(0, 0)) + 
#  scale_y_continuous(expand = c(0, 0),breaks=c(0.0,0.5,1.0))


#plot K6
plot5 <- all_data %>%
  filter(k == 6) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack",width=0.85) +
  facet_grid(~factor(pop, levels=c("NWJ","EJ","WJ","SJ","NEZ","SZ","B")),switch="x",scales="free", space="free_x") +
  theme_bw() +
  ylab("K=6") +
  theme(legend.position="none", 
    panel.border = element_rect(colour = "black", fill=NA, size=0.1), 
    panel.spacing.x = unit(0.2,"mm")) +
  theme(axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title = element_text(size=7),
    axis.ticks.x = element_blank(),
    plot.title = element_text(size=6, hjust=0.5, margin=margin(1,0,1,0)),
    strip.text.x = element_text(margin = margin(b = 0.1, t = 0.1)),
    axis.ticks.length=unit(0, "cm")) +
  theme(strip.background = element_blank(),
    strip.text.x = element_text(size=6)) +
  theme(plot.margin = margin(0.04,0.02,0.04,0.04, "cm")) #+
#  scale_x_continuous(expand = c(0, 0)) + 
# scale_y_continuous(expand = c(0, 0), breaks=c(0.0,0.5,1.0))

png(file="/faststorage/project/Coregonus/Aja/Dragonfly/ADMIXTURE_LD-thinned/ADMIXTURE_K2-K6.region.thinned2.png", units="cm", width = 15, height = 10, res=300)

lay <- rbind(c(1),
             c(2),
             c(3),
             c(4),
             c(5))
grid.arrange(plot1,plot2,plot3,plot4,plot5, heights=c(1,1,1,1,1.16), layout_matrix=lay)
dev.off()

