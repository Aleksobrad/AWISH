## TCR analysis for ITN030ST samples

## For clonal analysis, read tables into R, then create cd4 and cd8 objects with columns as below and use "run" function
# cd4= data[,c(1,2,3)] #pre-Tx unstim,CFSE-low,pre- or post-Tx sample
# cd8= data[,c(4,5,6)] #pre-Tx unstim,CFSE-low,pre- or post-Tx sample
# definitions for run function: cd4 as above, cd8 as above, fold and freqstim are fold expansion CFSE-low vs unstim and minimum frequency 
# in CFSE-low required, respectively, for identification of donor-reactive clones; cd4mincellcount and cd8mincellcount are minimum number of 
# productive templates of CD4 and CD8 unstimulated samples for subject, respectively; ambiguityRatio is used for minimizing sorting error 
# in identification of donor-reactive clones CD4 vs CD8

## For list of alloreactive clones, use listAlloreactive function where cd4 and cd8 objects have column 1 is unstim and column 2 is CFSE-low sample
# rownames of cd4 and cd8 must be defined as nucleotide sequence

## For JSD analysis, use jsdReport function where data is the columns of interest and topN are the number of clones to be compared between samples.

run <- function(cd4, cd8, fold=5, freqstim=0.0001, cd4mincellcount=Inf, cd8mincellcount=Inf, ambiguityRatio=5, filename="Default") {
  write.table(paste("\n Number of unique clones in postTx (from CD4 table):", length(which(cd4[,3]>0)),  sep=","),  quote=F, col.names = F, row.names=F, file=filename, append=T)
  write.table(paste("Number of unique clones in postTx (from CD8 table):", length(which(cd8[,3]>0)),  sep=","),  quote=F, col.names = F, row.names=F, file=filename, append=T)
 
  cd4counts<-cd4
  cd8counts<-cd8
  
  cd4 = normalize(cd4)
  cd8 = normalize(cd8)
  
  rows1 = cleanup(cd4[,1], cd8[,1], ratio = ambiguityRatio)
  rows2 = cleanup(cd4[,2], cd8[,2], ratio = ambiguityRatio)
  
  cd4 = exclude(cd4counts, rows1[,1] | rows2[,1])
  cd8 = exclude(cd8counts, rows1[,2] | rows2[,2])
  
  cd4freq=380000/cd4mincellcount*0.00001
  cd8freq=380000/cd8mincellcount*0.00001
  
  cat(paste("\n",colnames(cd4)[3], ": cd4\n"), file=filename, append=T)
  clonefrac(cd4,fold=fold,freqstim=freqstim, frequnstim=cd4freq, filename=filename)
  clonefreq(cd4,fold=fold,freqstim=freqstim, frequnstim=0, filename=filename)
  turnover(cd4,fold=fold,freqstim=freqstim, frequnstim=0, filename=filename)
  absoluteexpansion(cd4,fold=fold,freqstim=freqstim, frequnstim=cd4freq, filename=filename)
  
  cat(paste("\n",colnames(cd8)[3], ": cd8\n"), file=filename, append=T)
  
  clonefrac(cd8,fold=fold,freqstim=freqstim, frequnstim=cd8freq, filename=filename)
  clonefreq(cd8,fold=fold,freqstim=freqstim, frequnstim=0, filename=filename)
  turnover(cd8,fold=fold,freqstim=freqstim, frequnstim=0, filename=filename)
  absoluteexpansion(cd8,fold=fold,freqstim=freqstim, frequnstim=cd8freq, filename=filename)
}

reactiveClonescount <- function(datacounts, data, fold=5, freq=1e-5) {
  rclones = datacounts[data[,2] > freq & data[,2] > data[,1] * fold, ]
  return(rclones)
}

reactiveClones <- function(data, fold=5, freq=1e-5) {
  rclones = data[data[,2] > freq & data[,2] > data[,1] * fold, ]
  return(rclones)
}


clonefrac<-function(data,fold=5,freqstim=.0001, frequnstim=0, filename="Default"){
  data = normalize(data)

  allPre = data[data[,1] > frequnstim,]
  allBx = data[data[,3] > frequnstim, ]
  
  HvGinPre = reactiveClones(allPre, fold=fold, freq=freqstim)
  HvGinBx = reactiveClones(allBx, fold = fold, freq = freqstim)
  
  npre= nrow(allPre)
  nhvgpre= nrow(HvGinPre)
  nBx= nrow(allBx)
  nhvgpost= nrow(HvGinBx)
  
  # Expansion analysis (default): Is a detectable clone (above threshold freq), more likely to be alloreactive in postTx than preTx?
  #(i.e.) is the alloreactivity rate higher in postTx than preTx? 
  # allopre, Npre-allopre
  # allopost, Npost-allopost
  write.table("clone fraction", quote=F, col.names=F, row.names=F, file=filename, append=T)
  write.table(paste("Samples","detectable_clones","detectable_allo_clones","clone_frac", sep=","),quote=F,col.names=F,row.names=F, file=filename, append=T)
  write.table(paste("preTx", npre, nhvgpre,round(nhvgpre/npre,5),  sep=","),  quote=F, col.names = F, row.names=F, file=filename, append=T)
  write.table(paste("postTx",  nBx, nhvgpost,round(nhvgpost/nBx,5),  sep=","),  quote=F, col.names = F, row.names=F, file=filename, append=T)
  
  x = c(nhvgpost, nBx)
  y = c(nhvgpre, npre)
  if(sum(x)>0 && sum(y)>0){
    ft = fisher.test(cbind(x,y))
    write.table(paste("Fisher'sExactTest:\tp-value=\t",  signif(ft$estimate, 3),signif(ft$p.value,3),   "\t",  signif(ft$conf.int[1], 3), "-",signif(ft$conf.int[2], 3), sep="," ),  sep="," , quote=F, col.names = F, row.names=F, file=filename, append=T)
  }
}


clonefreq<-function(data,fold=5,freqstim=.0001, frequnstim=0, filename="Default"){
  datacount=data
  data = normalize(data)
  
  allPrecount = datacount[data[,1] > frequnstim,]
  allBxcount = datacount[data[,3] > frequnstim, ]
  allPre = data[data[,1] > frequnstim,]
  allBx = data[data[,3] > frequnstim, ]
  
  HvGinPrecount = reactiveClonescount(allPrecount, allPre, fold=fold, freq=freqstim)
  HvGinBxcount = reactiveClonescount(allBxcount, allBx, fold = fold, freq = freqstim)
  
  freqpre= sum(allPrecount[,1])
  hfreqpre= sum(HvGinPrecount[,1])
  freqBx= sum(allBxcount[,3])
  hfreqpost= sum(HvGinBxcount[,3])  
  
  write.table("cumulative frequency", quote=F, col.names=F, row.names=F, file=filename, append=T)
  write.table(paste("Samples","total_cells","allo_cells","allo_cum_freq", sep=","),quote=F,col.names=F,row.names=F, file=filename, append=T)
  write.table(paste("preTx", freqpre, hfreqpre, round(hfreqpre/freqpre, 5),  sep=","),  quote=F, col.names = F, row.names=F, file=filename, append=T)
  write.table(paste("postTx", freqBx, hfreqpost, round(hfreqpost/freqBx, 5),  sep=","),  quote=F, col.names = F, row.names=F, file=filename, append=T)
  
  x = c(hfreqpost, freqBx)
  y = c(hfreqpre, freqpre)
  if(sum(x)>0 && sum(y)>0){
    ft = fisher.test(cbind(x,y))
    write.table(paste("Fisher'sExactTest:\tp-value=", signif(ft$estimate, 3),signif(ft$p.value,3), " (95%CI:",  signif(ft$conf.int[1], 3), "-",signif(ft$conf.int[2], 3), ")", sep="," ),  sep="," , quote=F, col.names = F, row.names=F, file=filename, append=T)
  }
}


turnover<-function(data,fold=5,freqstim=.0001,frequnstim=0,filename="Default"){
  data=normalize(data)
  
  allPre = data[data[,1] > frequnstim, ]
  Detected = data[data[,1]>frequnstim & data[,3] > frequnstim, ]
  
  try = data[data[,3]>frequnstim,] #Looking in post-Tx
  
  allHvG = reactiveClones(data, fold=fold, freq=freqstim) #Total donor-reactive
  HvGDetected = reactiveClones(try, fold = fold, freq = freqstim) #D-reactive detected postTx
  
  nhvg= nrow(allHvG)
  npre=nrow(allPre)
  ndetected=nrow(Detected)
  nhvgdetected=nrow(HvGDetected)
  
  # repertoire turnover analysis: are preTx donor-reactive clones more likely to be persist than total preTx clones
  # N1,N0  (detected preTx clones, undetected preTx clones)
  # D1,D0 (donor-reactive clones, undetected donor-reactive clones)
    write.table("Turnover", quote=F, col.names=F, row.names=F, file=filename, append=T)
    write.table(paste("clone sets","in_postTx","not_in_postTx","persistent_ratio", sep=","),quote=F,col.names=F,row.names=F, file=filename, append=T)
    write.table(paste("preTx unstim", ndetected, npre-ndetected, ndetected/npre, sep=","),  quote=F, col.names = F, row.names=F, file=filename, append=T)
    write.table(paste("allo clones",  nhvgdetected, nhvg-nhvgdetected, nhvgdetected/nhvg, sep=","),  quote=F, col.names = F, row.names=F, file=filename, append=T)
    x = c(ndetected, npre-ndetected)
    y = c(nhvgdetected,nhvg-nhvgdetected) 
    if(sum(x)>0 && sum(y)>0){
      ft = fisher.test(cbind(y,x))
      write.table(paste("Fisher'sExactTest:", signif(ft$estimate, 3),signif(ft$p.value,3), " (95%CI:",  signif(ft$conf.int[1], 3), "-",signif(ft$conf.int[2], 3), ")", sep="," ),  sep="," , quote=F, col.names = F, row.names=F, file=filename, append=T)
    }
}


absoluteexpansion<-function(data,fold=5,freqstim=.0001, frequnstim=0, filename="Default"){
 data=normalize(data)
 
 allPre = data[data[,1] > frequnstim,]
 allBx = data[data[,3] > frequnstim, ]
 
 allHvG = reactiveClones(data, fold=fold, freq=freqstim)
 HvGinPre = reactiveClones(allPre, fold=fold, freq=freqstim)
 HvGinBx = reactiveClones(allBx, fold = fold, freq = freqstim)
 
 nhvgpre= nrow(HvGinPre)
 nhvgpost= nrow(HvGinBx)
 nhvg=nrow(allHvG)
  
  #absolute expansion analysis: are alloreative clones more likely to be detected in pre or post-transplant sample?
  # pre, N-pre
  # post, N-post   where N is the total number of clones termed alloreactive
    write.table("Absolute Expansion", quote=F, col.names=F, row.names=F, file=filename, append=T)
    write.table(paste("clone sets","allo_detected","allo_undetected","detection_rate", sep=","),quote=F,col.names=F,row.names=F, file=filename, append=T)
    write.table(paste("preTx", nhvgpre, nhvg-nhvgpre, nhvgpre/nhvg, sep=","),  quote=F, col.names = F, row.names=F, file=filename, append=T)
    write.table(paste("postTx",  nhvgpost, nhvg-nhvgpost, nhvgpost/nhvg, sep=","),  quote=F, col.names = F, row.names=F, file=filename, append=T)
    x = c(nhvgpre, nhvg-nhvgpre)
    y = c(nhvgpost, nhvg-nhvgpost) 
    if(sum(x)>0 && sum(y)>0){
      ft = fisher.test(cbind(y,x))
      write.table(paste("Fisher'sExactTest:\tp-value=\t", signif(ft$estimate, 3),signif(ft$p.value,3), "\t",  signif(ft$conf.int[1], 3), "-",signif(ft$conf.int[2], 3), sep="," ),  sep="," , quote=F, col.names = F, row.names=F, file=filename, append=T)
    }
}


normalize <- function(data) {
  nc = ncol(data)
  for (i in 1:nc) {
    data[,i] = data[,i ] / sum(data[,i])
  }
  return(data)
}


cloneCal <- function(array) {
  x = array[array >0 ] / sum(array)
  #  x = sort(array, decreasing=T)
  l = length(x)
  entropy = sum(x * -1 * log2(x))
  maxentropy = -log2(1/l)
  return(signif(1 - entropy / maxentropy, 2))
}

entropyCal<-function(array) {
  x = array[array >0 ] / sum(array)
  entropy = sum(x * -1 * log2(x))
  return(signif(entropy,2))
}

simpsonCal<-function(array) {
  x = array[array >0 ] / sum(array)
  simpson = sum(x * x)
  return(signif(simpson,2))
}

cleanup <- function(cd4, cd8, ratio = 2) {
  ## remove contaminated clones  
  ambi = (cd4 > 0 & cd8 > 0 & cd4 / cd8 > 1/ratio & cd4 / cd8 < ratio) 
  cd4exclude = (cd4 > 0 & cd8 > 0 & cd4 / cd8 <= 1/ratio ) 
  cd8exclude = (cd4 > 0 & cd8 > 0 & cd4 / cd8 >= ratio )  
  return(cbind(ambi | cd4exclude, ambi | cd8exclude))
}

exclude <- function(data, excluderows) {
  data = data[excluderows == F, ]
  return(data)
}

#before running, must do rownames(data)=data[,1]
listAlloreactive<-function(cd4,cd8, fold=5, freq1=.00001, ambiguityRatio=5) # rows must be indexed by clone ID
{
  # need only column 1 unstim and column 2 stim
  cd4 = normalize(cd4)
  cd8 = normalize(cd8)
  
  rows1 = cleanup(cd4[,1], cd8[,1], ratio = ambiguityRatio)
  rows2 = cleanup(cd4[,2], cd8[,2], ratio = ambiguityRatio)
  
  cd4 = exclude(cd4, rows1[,1] | rows2[,1])
  cd8 = exclude(cd8, rows1[,2] | rows2[,2])
  
  allHvGcd4 = reactiveClones(normalize(cd4), fold=fold, freq=freq1)
  allHvGcd8 = reactiveClones(normalize(cd8), fold=fold, freq=freq1)
  
  return(list(rownames(allHvGcd4), rownames(allHvGcd8)))
}

jsdReport<-function(data, topN=-1)
{
  out<-matrix(nrow=ncol(data),ncol=ncol(data), dimnames=list(colnames(data), colnames(data)))
  for(i in 1:ncol(data)){
    for(j in 1:ncol(data)){
      if(topN==-1){
        out[i,j]<-jensen_shannon(data[,i], data[,j])
      }
      else{
        a<-order(data[,i],decreasing=TRUE)[1:topN]
        b<-order(data[,j],decreasing=TRUE)[1:topN]
        z<-data[union(a,b),]
        out[i,j]<-jensen_shannon(z[,i],z[,j])
      }
    }
  }
  return(out)
}

shannon.entropy <- function(p)
{
  if (min(p) < 0 || sum(p) <= 0)
    return(NA)
  p.norm <- p[p>0]/sum(p)
  -sum(log2(p.norm)*p.norm)
}

cloneCal <- function(p) {
  x = p[p>0] / sum(p)
  l = length(x)
  entropy = shannon.entropy(p)
  maxentropy = -log2(1/l)
  return(signif(1 - entropy / maxentropy, 2)) 
}

jensen_shannon <- function(p, q){
  ## JSD = H(0.5 *(p +q)) - 0.5H(p) - 0.5H(q)
  # H(X) = \sum x_i * log2(x_i)
  #  p = p[p >0 & q >0]
  #  q = q[p>0 & q>0]
  p = p / sum(p)
  q = q / sum(q)
  Hj = shannon.entropy(0.5 *(p+q)) 
  Hp = shannon.entropy(p) 
  Hq = shannon.entropy(q)
  
  jsd = Hj - 0.5*(Hp+Hq)
  #	cat(Hj, Hp, Hq, jsd, "\n")
  return(jsd)
}

