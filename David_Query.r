# File: David_Query.r
# Purpose: Repeat limma analysis for all possible sample groupings
# Created: 2/18/18

# Setup
library('RDAVIDWebService')


# Please Edit the following information. 
workingDirectory<-"/home/adam/Desktop/Masters_Thesis_2017/Never_Late_Pipeline"  # Directory on your computer
dataFile <- "randomGenes_1-100.Rda"                                             # Name of the file you are processing
email <- "abf@udel.edu"                                                         # email address you registered with DAVID

# # # # # # # # # # # # # # # #Do not edit below this line # # # # # # # # # # # # # # # # #
setwd(workingDirectory)
davidUrl<-"https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/"
load(dataFile)
randomGenes<-get(dataFile)
rm(list=dataFile)

if(exists('simDavidSummary')){rm(simDavidSummary)}
if(exists('simDavidTerms')){rm(simDavidTerms)}

for (i in 1:ncol(randomGenes)){
  t<-Sys.time()
  print(i)
  david<-DAVIDWebService$new(email=email, url = davidUrl)
  gl<-randomGenes[,i]
  result<-addList(david, gl, idType="AFFYMETRIX_3PRIME_IVT_ID", listName="rand", listType="Gene")
  setAnnotationCategories(david, categories = c("GOTERM_BP_DIRECT", "GOTERM_BP_FAT"))
  chart<-getFunctionalAnnotationChart(david)
  chart$Perm<-i

  if(exists('simDavidTerms')){
    simDavidTerms<-rbind(simDavidTerms, chart)
  }
  
  else{
    simDavidTerms<-chart
  }

  if(exists('simDavidSummary')){
    tmp<-data.frame(
      Perm = i,
      Min_FAT_Benjamini = min(chart[chart$Category == "GOTERM_BP_FAT", "Benjamini"]),
      Min_DIRECT_Benjamini = min(chart[chart$Category == "GOTERM_BP_DIRECT", "Benjamini"]),
      Count_FAT_Sig = nrow(chart[chart$Category == "GOTERM_BP_FAT" & chart$Benjamini < 0.05, ]),
      Count_DIRECT_Sig = nrow(chart[chart$Category == "GOTERM_BP_DIRECT" & chart$Benjamini < 0.05, ])
    )
    simDavidSummary<-rbind(simDavidSummary, tmp)
    rm(tmp)
  }
  else{
    simDavidSummary<-data.frame(
      Perm = i,
      Min_FAT_Benjamini = min(chart[chart$Category == "GOTERM_BP_FAT", "Benjamini"]),
      Min_DIRECT_Benjamini = min(chart[chart$Category == "GOTERM_BP_DIRECT", "Benjamini"]),
      Count_FAT_Sig = nrow(chart[chart$Category == "GOTERM_BP_FAT" & chart$Benjamini < 0.05, ]),
      Count_DIRECT_Sig = nrow(chart[chart$Category == "GOTERM_BP_DIRECT" & chart$Benjamini < 0.05, ])
    )
  }
  while((t - Sys.time()) > -30){}
}

save(list = c('simDavidSummary', 'simDavidTerms'), file = paste(email, '_out.Rda'))



