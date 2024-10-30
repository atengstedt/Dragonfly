# Rescale the ith iteration result of PSMC, and make ready for plotting.
# Modified from old script by Shenglin Liu, Mar 25, 2019.

# i.iteration: the ith iteration
# mu: mutation rate
# s: bin size
# g: years per generation
psmc.result<-function(file,i.iteration=25,mu=2.5e-8,s=100,g=3.5)
{
	X<-scan(file=file,what="",sep="\n",quiet=TRUE)
	START<-grep("^RD",X)
	END<-grep("^//",X)
	X<-X[START[i.iteration+1]:END[i.iteration+1]]
	
	TR<-grep("^TR",X,value=TRUE)
	RS<-grep("^RS",X,value=TRUE)
	
	theta0<-as.numeric(strsplit(TR,"\t")[[1]][2])
	N0<-theta0/4/mu/s
	
	a<-t(as.data.frame(strsplit(RS,"\t")))
	Time<-2*N0*as.numeric(a[,3])*g
	Ne<-N0*as.numeric(a[,4])
	
	n<-length(Ne)
	Time<-c(as.numeric(rbind(Time[-n],Time[-1])),Time[n])
	Ne<-c(as.numeric(rbind(Ne[-n],Ne[-n])),Ne[n])
	
	data.frame(Time,Ne)
}