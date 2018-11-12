# Setup For Analysis
library(dplyr)
library(reshape2)
library(ggplot2)
library(affy)
wd<-"C:/Users/el3548/OneDrive - DuPont/Never_Late_Pipeline/Never_Late_Pipeline/Null_Ontology_Distribution"
setwd(wd)

# Iterate Through
files<-list.files(wd)[grep("*output.Rda", list.files(wd))]
for(i in files){
  load(i)
  if(exists('DavidSummary')){
    DavidSummary<-rbind(DavidSummary, simDavidSummary)
  }
  else{
    DavidSummary<-simDavidSummary
  }
  if(exists('DavidTerms')){
    DavidTerms<-rbind(DavidTerms, simDavidTerms[order(simDavidTerms$Benjamini),])
  }
  else{
    DavidTerms<-simDavidTerms[order(simDavidTerms$Benjamini),]
  }
  
  rm(list=c('simDavidSummary', 'simDavidTerms'))
}


# Validate Queries: Confirm that the top 5 (Benjamini) Terms are different for each dataset
Sets<-unique(DavidTerms$dataFile)
if(exists('checkDavidTerms')){rm(checkDavidTerms)}
for(i in Sets){
  if(exists('checkDavidTerms')){
    tmp<-data.frame(
      term1 = DavidTerms[DavidTerms$dataFile == i, ][1,'Term'],
      term2 = DavidTerms[DavidTerms$dataFile == i, ][2,'Term'],
      term3 = DavidTerms[DavidTerms$dataFile == i, ][3,'Term'],
      term4 = DavidTerms[DavidTerms$dataFile == i, ][4,'Term'],
      term5 = DavidTerms[DavidTerms$dataFile == i, ][5,'Term']
    )
    row.names(tmp)[1]<-i
    checkDavidTerms <-rbind(checkDavidTerms, tmp)
    rm(tmp)
  }
  else{
    checkDavidTerms<-data.frame(
      term1 = DavidTerms[DavidTerms$dataFile == i, ][1,'Term'],
      term2 = DavidTerms[DavidTerms$dataFile == i, ][2,'Term'],
      term3 = DavidTerms[DavidTerms$dataFile == i, ][3,'Term'],
      term4 = DavidTerms[DavidTerms$dataFile == i, ][4,'Term'],
      term5 = DavidTerms[DavidTerms$dataFile == i, ][5,'Term']
    )
    row.names(checkDavidTerms)[1]<-i
  }
  
}

checkDavidTerms$check <-0
for( i in Sets){
  x<-as.character(checkDavidTerms[i,1:5])
  ck<-0
  for(j in setdiff( Sets, i)){
    y<-as.character(checkDavidTerms[j,1:5])
    if(length(intersect(x,y)) > 3){
      ck<-1
    }
  }
  checkDavidTerms[i,'check']<-ck
}



# Plot Distribution of number of significant hits
opar = par()
png(filename = "Term_Counts",   # Opens a PNG device that writes to the named file
    width = 8.8,
    height = 5.5,
    units = "in",
    res = 96                          # In pixels per inch, determines maximum scaling range (and file size!)
)

par(mfrow = c(1,2))
hist(DavidSummary$Count_FAT_Sig,
     xlab = "Number of Terms detected with Benjamini < 0.05",
     ylab = "Frequency in 1000 Random Gene Lists",
     main = "Counts for David Category \nGOTERM_BP_FAT",
     ylim = c(0, 800)
)

hist(DavidSummary$Count_DIRECT_Sig,
     xlab = "Number of Terms detected with Benjamini < 0.05",
     ylab = "Frequency in 1000 Random Gene Lists",
     main = "Counts for David Category \nGOTERM_BP_DIRECT",
     ylim = c(0,1000),
     yaxp = c(0, 1000, 5))

# Tabulate Significant Term Counts

dev.off()
table(DavidSummary$Count_DIRECT_Sig)
FAT_Counts<-table(
  cut(DavidSummary$Count_FAT_Sig, 
      breaks = c(0,50, 100, 150, 200, 250), include.lowest = T)
)

DIRECT_Counts<-table(DavidSummary$Count_DIRECT_Sig)

# Plot Distributions of minimum observed Bennjamini
png(filename = "Min_Benjamini",   # Opens a PNG device that writes to the named file
    width = 8.8,
    height = 5.5,
    units = "in",
    res = 96                          # In pixels per inch, determines maximum scaling range (and file size!)
)
par(mfrow=c(1,2))
hist(0 - log10(DavidSummary$Min_FAT_Benjamini),
     xlab = "Negative Log 10 Min. Benjamini",
     ylab = "Frequency in 1000 Random Gene Lists",
     main = "Minimum Benjamini:\n GOTERM_BP_FAT",
     ylim = c(0, 200))
segments(0-log10(0.05), 0,0-log10(0.05), 200, col = "red")

hist(0 - log10(DavidSummary$Min_DIRECT_Benjamini),
     xlab = "Negative Log 10 Min. Benjamini",
     ylab = "Frequency in 1000 Random Gene Lists",
     main = "Minimum Benjamini:\n GOTERM_BP_DIRECT",
     ylim = c(0, 700))
segments(0-log10(0.05), 0,0-log10(0.05), 700, col = "red")
dev.off()




# Terms with a Benjamini < 0.05 and Detection Frequency In the Category GOTERM_BP_FAT
SigDavidTerms<-DavidTerms[DavidTerms$Benjamini < 0.05,]
length(SigDavidTerms[SigDavidTerms$Category == "GOTERM_BP_FAT", 'Term'])
TermFrequencies<-table(SigDavidTerms[SigDavidTerms$Category == "GOTERM_BP_FAT", 'Term'])

TF<-data.frame(
	Term = names(TermFrequencies),
	Frequency = as.numeric(TermFrequencies),
	stringsAsFactors = F
)


TabulatedTerms<-data.frame(
  Term = names(TermFrequencies),
  Count = as.numeric(TermFrequencies),
  Percentage = 100 * as.numeric(TermFrequencies) / length(SigDavidTerms[SigDavidTerms$Category == "GOTERM_BP_FAT", 'Term'])
)
write.table(
  TabulatedTerms[order(TabulatedTerms$Percentage, decreasing = T),],
  file = "Recurring_FAT_Terms.txt",
  quote = F, 
  sep = "\t",
  row.names = F
)


term<-"GO:0000003~reproduction"
term<-"GO:0012501~programmed cell death"
hist(0-log10(SigDavidTerms[SigDavidTerms$Term == term & SigDavidTerms$Category == "GOTERM_BP_FAT",'Benjamini']))

TF$MinEase<-sapply(TF$Term, function(x){
min(SigDavidTerms[SigDavidTerms$Term == x & SigDavidTerms$Category == "GOTERM_BP_FAT",'Benjamini'])})
head(TF)


min(SigDavidTerms[SigDavidTerms$Term == term & SigDavidTerms$Category == "GOTERM_BP_FAT",'Benjamini'])