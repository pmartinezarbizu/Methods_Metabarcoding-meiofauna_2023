#script to remove ASV which produce "no hit" in blast.
#
# copyright Pedro Martinez Arbizu & Sahar Khodami (2023)
# pmartinez@senckenberg.de
# 
pathblast <- "/Dataexchange/blastn/output/blast/"

blast <- read.csv(paste(pathblast,'all.otus.V1V2.dada2.final.blast.tab.txt',sep=''),header=FALSE,sep='')
otutab <- read.csv(paste(pathblast,'all.otutab.V1V2.dada2.tab.txt',sep=''),header=FALSE,sep='')
otuseq <- read.csv(paste(pathblast,'otuseqs.V1V2.dada2.txt',sep=''),header=FALSE,sep='')

otunew <- c()
otuseqnew <-c()

for (elem in 1:nrow(blast)){
otunew <- rbind(otunew,otutab[blast[elem,1],])
otuseqnew <- rbind(otuseqnew, as.character(otuseq[blast[elem,1],]))
}

otunew <- cbind(blast,otunew)

write.table(otunew,paste(pathblast,'ASV_Table.txt',sep=''),quote=FALSE,row.names=FALSE,col.names=FALSE)

write.table(otuseqnew,paste(pathblast,'otuseqnew.txt',sep=''),quote=FALSE,row.names=FALSE,col.names=FALSE)

