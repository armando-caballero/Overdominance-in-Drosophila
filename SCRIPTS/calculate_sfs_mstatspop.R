#Estimate Variability
deltak <- function(i,j) {
  if(i==j) return(1)
  return(0)
}
weight.watt.folded <-function(nsam) {
  w <- array(0,dim=c(floor(nsam/2)))
  for(i in 1:length(w)) {
    w[i] <- nsam/(i*(nsam-i)*(1+deltak(i,nsam-i)))
  }
  w
}
weight.taj.folded <-function(nsam) {
  w <- array(0,dim=c(floor(nsam/2)))
  for(i in 1:length(w)) {
    w[i] <- nsam/(1+deltak(i,nsam-i))
  }
  w
}
weight.fuli.folded <-function(nsam) {
  w <- array(0,dim=c(floor(nsam/2)))
  w[1] <- nsam
  w
}
phi.i <- function(nsam) {
  w <- array(0,dim=c(floor(nsam/2)))
  for(i in 1:length(w)) { 
    w[i] <- (nsam/(i*(nsam-i))) / (1+deltak(i,nsam-i))
  }
  w
}
Calc.Theta.folded <- function(sfs,w,phi) {
  th <- 0
  for(i in 1:length(sfs)) {
    th <- th + w[i] * sfs[i] * 1/phi[i]
  }
  th <- th/(sum(w))
  return(th)
}

#READ DATA
args = commandArgs(trailingOnly=TRUE)
if(length(args)<1) {
  args[1]<-"Drosophila_Arm"
}
data_namefile <- args[1]

pathr <- c(".") #c("/Users/sramos/Desktop/Armando-Overdominance2025/SeqData-plink/VariabilityAnalysis")
SIT <- c("TOTAL","Silent","Syn","NSyn")

theta <- matrix(0,ncol=3,nrow=4)
colnames(theta) <- c("Watt","Pi","FuLi")
rownames(theta) <- SIT
theta
theta.no2L <- theta
theta.no2L

SFS.all <- NULL

pdf(file=sprintf("%s/SFS_mstatspop.pdf",pathr))

for(sites in SIT) {
  data <- read.table(file=sprintf("%s/%s.tfa_%s_Variability_o2.COLS.txt",pathr,data_namefile,sites),header=T)
  head(data)
  dim(data)
  #obtain nsam
  nsam <- data$nsam.0..
  nsam
  nsam <- nsam[1]
  nsam
  #obtain lengths of chromosomes
  len <- round(data$Eff_length2_pop.0..,0)
  len
  #obtain folded SFS
  SFS <- data[,c(22:(22+51-1))]
  SFS
  
  #plot
  
  plot(x=c(1:(nsam/2)),y=SFS[1,],type="b",xlab="SFS folded",ylab="freq",main=sprintf("Folded SFS %s per chr. arm",sites))
  lines(x=c(1:(nsam/2)),SFS[2,],type="b",col="blue")
  lines(x=c(1:(nsam/2)),SFS[3,],type="b",col="red")
  lines(x=c(1:(nsam/2)),SFS[4,],type="b",col="green")
  legend("topright",lwd=c(1,1,1,1),col=c("black","blue","red","green"),legend = c(data$scaffold_name.)[-5])
  
  SFS.tot <- apply(SFS[c(1:4),],2,sum)
  plot(x=c(1:(nsam/2)),SFS.tot,type="b",xlab="SFS folded",ylab="freq",main=sprintf("Folded SFS %s",sites))
  SFS.tot.no2L <- apply(SFS[c(2:4),],2,sum)
  plot(x=c(1:(nsam/2)),SFS.tot.no2L,type="b",xlab="SFS folded",ylab="freq",main=sprintf("Folded SFS %s NO-2L",sites))
  
  #SFS.tot <- SFS.tot.no2L
  L.tot <- sum(len[1:4])
  L.tot.no2L <- sum(len[2:4])
  
  SFS.all <- rbind(SFS.all,c(nsam,L,SFS.tot))
  SFS.all.no2L <- rbind(SFS.all,c(nsam,L,SFS.tot.no2L))
  
  #Calculate theta
  phi <- phi.i(nsam)
  w.wat <- weight.watt.folded(nsam)
  theta.w <- Calc.Theta.folded(SFS.tot,w.wat,phi)
  theta.w/L.tot
  theta.w.no2L <- Calc.Theta.folded(SFS.tot.no2L,w.wat,phi)
  theta.w.no2L/L.tot.no2L
  
  w.taj <- weight.taj.folded(nsam)
  theta.t <- Calc.Theta.folded(SFS.tot,w.taj,phi)
  theta.t/L.tot #PI
  theta.t.no2L <- Calc.Theta.folded(SFS.tot.no2L,w.taj,phi)
  theta.t.no2L/L.tot.no2L #PI
  
  w.fl <- weight.fuli.folded(nsam)
  theta.fl <- Calc.Theta.folded(SFS.tot,w.fl,phi)
  theta.fl/L.tot
  theta.fl.no2L <- Calc.Theta.folded(SFS.tot.no2L,w.fl,phi)
  theta.fl.no2L/L.tot.no2L
  
  theta[which(sites==SIT),] <- c(round(theta.w/L.tot,5),round(theta.t/L.tot,5),round(theta.fl/L.tot,5))
  theta.no2L[which(sites==SIT),] <- c(round(theta.w.no2L/L.tot.no2L,5),round(theta.t.no2L/L.tot.no2L,5),round(theta.fl.no2L/L.tot.no2L,5))
}
dev.off()

write.table(x="Theta/nt",file=sprintf("%s/Theta.txt",pathr), quote = F, sep="\t", append=F, col.names = F,row.names = F,eol = "\t")
write.table(x=theta,file=sprintf("%s/Theta.txt",pathr),quote = F, sep="\t", append=T)

write.table(x="Theta.no2L/nt",file=sprintf("%s/Theta.no2L.txt",pathr), quote = F, sep="\t", append=F, col.names = F,row.names = F,eol = "\t")
write.table(x=theta.no2L,file=sprintf("%s/Theta.no2L.txt",pathr),quote = F, sep="\t", append=T)

rownames(SFS.all) <- SIT
write.table(x=c("Type","nsam","length",sprintf("SFS.%s",c(1:(nsam/2)))),file=sprintf("%s/SFS.all.txt",pathr), quote = F, sep="\t", append=F, col.names = F,row.names = F,eol = "\t")
write.table(x="",file=sprintf("%s/SFS.all.txt",pathr), quote = F,append=T,sep="",row.names = F,col.names = F)
write.table(x=SFS.all,file=sprintf("%s/SFS.all.txt",pathr),quote = F, sep="\t", append=T,col.names = F)
