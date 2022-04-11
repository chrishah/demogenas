#setwd("/home/lv71312/hahnc/Monogeniac/samples/Gsalv_Tyri/kmc")
args = commandArgs(trailingOnly=TRUE)

file = args[1]
data = read.table(file, sep="\t", col.names=c("kmer_cov","count"))

sample=args[2]
#kmer_size=47
#cum_read_length=100000000000
#avg_read_length=145
kmer_size=as.integer(args[3])
cum_read_length=as.numeric(args[4])
avg_read_length=as.integer(args[5])

minsearch=5
maxsearch=200

frequency <- data$count/sum(as.numeric(data$count))

start <- which(data$kmer_cov == minsearch)-1
end <- which(data$kmer_cov == maxsearch)
cat("searching from",start,"to",end,fill=T)
for (i in start:end+1){
#  cat(i, data$count[i], max(data$count[i:end]), data$kmer_cov[i],frequency[i],fill=T)
  if ( max(data$count[i:end]) != data$count[i]){
#    cat("below max", fill=T)
#    if ( data$count[i-1] > data$count[i] & data$count[i] < data$count[i+1] & max(data$count[i:end])/data$count[i] >= 2 ){
    if ( data$count[i-1] > data$count[i] & data$count[i] < data$count[i+1] ){
#      cat("candidate",fill=T)
      #check if valley and next peak are at least 5 steps apart
      if ( (which(data$count == max(data$count[i:end])) - i) >= 5 ){
        break
      }
    }
#  } else {
#    cat("still at max",fill=T)
  }
}

#if ( which(data$count == max(data$count[minsearch:maxsearch])) > minsearch ){
if (i <= end){
  cat("found valley at kmercoverage:",i, fill=T)
  #find peak
#  valley <- i-1
#  index <- which(data$count == max(data$count[valley:maxsearch]))
  index <- which(data$count == max(data$count[i:end]))
  mode <- data$kmer_cov[index]
#  cat("index:",index,", mode:",mode,fill=T)
  
#  xmax <- round(mode*3,-1)
#  cat("X axis max:", xmax, fill=T)
#  cat("Y axis max:", round(frequency[index]*2,3), frequency[index], fill=T)
  peak_kmer_coverage = mode
  cat("found peak at kmer coverage:",peak_kmer_coverage, fill=T)
  read_coverage_peak=peak_kmer_coverage*avg_read_length/(avg_read_length-kmer_size+1)
  read_coverage_peak
  est_genome_size=cum_read_length/read_coverage_peak
  est_genome_size
  
  pdf(paste(sample, "-k", kmer_size, "-distribution.pdf", sep=""))
  plot(data$kmer_cov[0:end], frequency[0:end], ylim=c(0,round(frequency[index]*2,3)), xlim=c(0,end), pch=20, 
       cex=0.5, xlab=paste(kmer_size,"-mer coverage",sep=""), ylab="frequency", main=sample, axes=F)
  lines(data$kmer_cov[0:end], frequency[0:end], type="b")
  
  ymax <- 0
  ymax_decimals <- 0
  while (ymax == 0) {
    ymax_decimals = ymax_decimals+1
    ymax <- round(frequency[index]*2,ymax_decimals)
    #  print(ymax)
  }
  
  yticks <- 0
  yticks_decimals <- 0
  while (yticks == 0) {
    yticks_decimals = yticks_decimals+1
    yticks <- round(frequency[index]*2/4,yticks_decimals)
    #  print(yticks)
  }
  
  axis(1,at=seq(0,round(mode*3)+10,10),labels=T,las=2)
  axis(2,at=seq(0,round(frequency[index]*2,ymax_decimals),round((frequency[index]*2)/4,yticks_decimals)),labels=T)
  
  abline(v=peak_kmer_coverage,lty=2)
  legend(x="topright", legend=c(paste(
    " total data:", round(cum_read_length/1e+09,1), "Gbp\n",
    "peak", kmer_size, "-mer coverage =", mode, "x\n",
    "peak read coverage =", round(read_coverage_peak, 2), "x\n",
    "est. genome size =", round(est_genome_size/1e+06,1), "Mb"))
    ,bty="n")
  box(which = "plot")
  dev.off()
  
  pdf(paste(sample, "-k", kmer_size, "-distribution-full.pdf", sep=""))
  plot(data$kmer_cov[0:end], frequency[0:end], pch=20, 
       cex=0.5, xlab=paste(kmer_size,"-mer coverage",sep=""), ylab="frequency", main=sample, axes=T)
  lines(data$kmer_cov[0:end], frequency[0:end], type="b")
  abline(v=peak_kmer_coverage,lty=2)
  legend(x="topright", legend=c(paste(
    " peak", kmer_size, "-mer coverage =", mode, "x\n",
    "peak read coverage =", round(read_coverage_peak, 2), "x\n",
    "est. genome size =", round(est_genome_size/1e+06,1), "Mb"))
    ,bty="n")
  dev.off()
  
} else {
  cat("found no peak", fill=T)
  pdf(paste(sample, "-k", kmer_size, "-distribution-full.pdf", sep=""))
  plot(data$kmer_cov[0:end], frequency[0:end], pch=20,
       cex=0.5, xlab=paste(kmer_size,"-mer coverage",sep=""), ylab="frequency", main=sample, axes=T)
  lines(data$kmer_cov[0:end], frequency[0:end], type="b")
  legend(x="topright", legend=c(paste(" Found no peak between coverage", minsearch, "and", maxsearch, "x ")))
  dev.off()
}

