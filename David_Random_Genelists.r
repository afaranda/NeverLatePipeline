# File: David_Random_Genelists.r
# Purpose: Repeat limma analysis for all possible sample groupings
# Created: 2/11/18

# Setup
library('affy')
library('limma')
library('RDAVIDWebService')
library(hgu133a.db)
library('dplyr')

wd<-"/home/adam/Desktop/Masters_Thesis_2017/Never_Late_Pipeline"
setwd(wd)
davidUrl<-"https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/"

# Generate Random DAVID distribution -- randomly genereated gene lists tested for enrichment
aProbes<-mappedkeys(hgu133aACCNUM)
aProbes<-aProbes[-grep("AFFX", aProbes)]

set.seed(15926)
randomGenes<-data.frame(row.names = 1:500)
for (i in 1:1000){
  randomGenes[paste('set', i, sep="_")]<-sample(aProbes, 500, replace = F)
}

# Split random lists into 10 subsets
for(i in 1:10){
  m <- (100*i)-99
  n <- 100*i
  
  fn<-paste("randomGenes_",m,"-",n,".Rda", sep = "")
  assign(fn, randomGenes[m:n])
  save(list = fn, file = fn)
  
}


# Test Rda Files
load("randomGenes_1-100.Rda")
load("randomGenes_101-200.Rda")
load("randomGenes_201-300.Rda")
load("randomGenes_301-400.Rda")
load("randomGenes_401-500.Rda")
load("randomGenes_501-600.Rda")
load("randomGenes_601-700.Rda")
load("randomGenes_701-800.Rda")
load("randomGenes_801-900.Rda")
load("randomGenes_901-1000.Rda")

for(i in 1:10){
  m <- (100*i)-99
  n <- 100*i
  check <- get(paste("randomGenes_",m,"-",n,".Rda", sep = ""))
  print(setdiff(randomGenes[,m], check[,1]))
  
}

# Please Edit the following information. 
workingDirectory<-"/home/adam/Desktop/Masters_Thesis_2017/Never_Late_Pipeline"
dataFile <- "randomGenes_1-100.Rda"
email <- "abf@udel.edu"

# Do not edit below this line
davidUrl<-"https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/"
load(dataFile)
randomGenes<-get(dataFile)
rm(list=dataFile)

if(exists('simDavidSummary')){rm(simDavidSummary)}
if(exists('simDavidTerms')){rm(simDavidTerms)}

for (i in 1:10){
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






