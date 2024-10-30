#R script originally from https://jujumaan.com/2017/07/15/linkage-disequilibrium-decay-plot/

# assign file names
input<-commandArgs(trailingOnly=T)
file<-input[1]
plot<-input[2]

#run R
pop <- read.table(file, header=TRUE)

# 2x number of individuals
n = 58

# Change fifth column name and add dist column
names(pop)[5] <- "rsq"
pop$dist <- ((pop$POS2 - pop$POS1)/1000000)

Cstart <- c(C=0.1)
# let's fit a non linear model using the arbitrary C value
modelC <- nls(rsq ~ ((10+C*dist)/((2+C*dist)*(11+C*dist)))*(1+((3+C*dist)*(12+12*C*dist+(C*dist)^2))/(n*(2+C*dist)*(11+C*dist))), data=pop, start=Cstart, control=nls.control(maxiter=100))
# extract rho, the recombination parameter, 4Nr
rho <- summary(modelC)$parameters[1]
# feed in the new value of rho to obtain LD values adjusted for their distances along the chromosome/genome
newrsq <- ((10+rho*pop$dist)/((2+rho*pop$dist)*(11+rho*pop$dist)))*(1+((3+rho*pop$dist)*(12+12*rho*pop$dist+(rho*pop$dist)^2))/(n*(2+rho*pop$dist)*(11+rho*pop$dist)))
newfile <- data.frame(pop$dist, newrsq)
#maxld <- max(pop$rsq) #using max LD value from initial input file
maxld <- max(newfile$newrsq) #using max LD value from adjusted data
halfdecay = maxld*0.5
halfdecaydist <- newfile$pop.dist[which.min(abs(newfile$newrsq-halfdecay))]
threshold<-0.1
thresholddist <- newfile$pop.dist[which.min(abs(newfile$newrsq-threshold))]
newfile <- newfile[order(newfile$pop.dist),]

# create the plot window
tiff(plot, units="mm", width=100, height=80, res=300)
par(mar=c(5.1,4.1,3.1,1.1), mgp=c(2.5,1,0))

# plot
plot(pop$dist, pop$rsq, main="pop", pch=".", cex=2, xlab="Distance (Mbp)", ylab=expression("LD (r"^2*")"))
lines(newfile$pop.dist, newfile$newrsq, col="blue", lwd=3)
#abline(v=halfdecaydist, col="red")
#mtext(round(halfdecaydist,2), side=3, line=0.05, at=halfdecaydist, cex=0.75, col="red")
abline(v=thresholddist, col="red")
mtext(round(thresholddist,2), side=3, line=0.05, at=thresholddist, cex=0.75, col="red")
#side=3 above graph

dev.off()