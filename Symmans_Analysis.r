# Process Demographics and assign samples to groups
# Setup
wd<-"/home/adam/Desktop/Masters_Thesis_2017/Never_Late_Pipeline/Symmans_Analysis"
celDir<-"/mnt/Storage/adam/kgsCelFiles/"
softDir<-"/mnt/Storage/adam/kgsSoftFiles"
studyName<-"Symmans"
setwd(wd)
library(limma)
library(affy)
library(GEOquery)
library(dplyr)
library(VennDiagram)
library(hgu133acdf)
library(hgu133a.db)
library("annotate")
library(arrayQualityMetrics)
source("/home/adam/Desktop/Masters_Thesis_2017/Never_Late_Pipeline/DataProcessingScripts/KeyValueExtract.r")
source("/home/adam/Desktop/Masters_Thesis_2017/Never_Late_Pipeline/DataProcessingScripts/Covariate_Descriptive_Stats.r")
options(download.file.method.GEOquery = "libcurl") # Required to properly fetch GPLXX platforms

# Read covariates from GSE
seriesMatrix<-getGEO("GSE17705", destdir = softDir)[[1]]     # GEO gene expression data set

# Process GEO Demographics Data from the soft file and assign samples to categories
phenotype<-keysToTable(pData(seriesMatrix)[,grep("characteristics_ch1",names(pData(seriesMatrix)))])

# Apply standardized column names, edit the search pattern in the 'grep' statement. 
names(phenotype)[grep("estrogen receptor (er) status", names(phenotype), fixed =T)] <-"ER_STATUS"
names(phenotype)[grep("event time (years)", names(phenotype), fixed =T)] <-"DMFS_TIME"
names(phenotype)[grep("distant relapse (1=dr, 0 censored)", names(phenotype), fixed =T)] <-"EVENT_DMFS"
names(phenotype)[grep("patient id", names(phenotype), fixed =T)] <-"Patient_ID"

# Convert time and event columns to numeric, and add grouping column
phenotype$DMFS_TIME<-as.numeric(phenotype$DMFS_TIME)
phenotype$EVENT_DMFS<-as.numeric(phenotype$EVENT_DMFS)
phenotype$`nodal status (0=negative, 1=positive, na=not applicable)`<-as.factor(phenotype$`nodal status (0=negative, 1=positive, na=not applicable)`)

phenotype$group<-factor(nrow(phenotype),levels=c("LATE", "NEVER"))


# Merge in filenames and re-assign row names.
phenotype$geo_accession<-row.names(phenotype)
phenotype<-merge(
  x=phenotype,
  y=pData(seriesMatrix)[,c("geo_accession", "supplementary_file")],
  by='geo_accession'
)


phenotype$supplementary_file<-sapply(
  strsplit(
    as.character(
      phenotype$supplementary_file), 
    "/", fixed=T), function(x){x[9]})
row.names(phenotype)<-phenotype$supplementary_file

# Get last event
lastEvent<-max(phenotype[phenotype$EVENT_DMFS == 1, "DMFS_TIME"], na.rm = T)

# Assign Groups and filter on unused samples
phenotype[
  !is.na(phenotype$DMFS_TIME) &
  !is.na(phenotype$EVENT_DMFS) &
  phenotype$DMFS_TIME > 5 &
  phenotype$DMFS_TIME <= lastEvent &
  phenotype$EVENT_DMFS == 1,
  "group"
]<-"LATE"

phenotype[
  !is.na(phenotype$EVENT_DMFS) &
    !is.na(phenotype$DMFS_TIME) &
    phenotype$DMFS_TIME>lastEvent &
    phenotype$EVENT_DMFS == 0,
  "group"
  ]<-"NEVER"


phenotype<-phenotype[phenotype$group %in% c("LATE", "NEVER"),]

# Drop samples that have been identified as outliers in previous runs.
# Samples were analyzed at two different labs. Four LATE samples analzyed
# at lab "JBI" were excluded
outlier<-c("GSM441627",
           "GSM441691",
           "GSM441692",
           "GSM441705"
           )
phenotype<-phenotype[!(phenotype$geo_accession %in% outlier), ]
phenotype<-phenotype[c("group", names(phenotype)[names(phenotype)!="group"])]
# Prepare Summary Statistics
descTable<-statsByType(phenotype)
writeDescTables(descTable)

statsTable<-calculatePairwise(phenotype)
writePairwiseStats(statsTable)

figuresByType(phenotype)

# Calculate Fishers for node manually, due to one NA
# tbl<-table(
#   phenotype[phenotype$`nodal status (0=negative, 1=positive, na=not applicable)` != "NA",
#           c('group', 'nodal status (0=negative, 1=positive, na=not applicable)')]
# )
# tbl<-tbl[,1:2]
# fisher.test(tbl)

# Build AffyBatch
sym.bat<-ReadAffy(
  filenames = row.names(phenotype),
  phenoData = phenotype,
  celfile.path = celDir
)
sym.eset<-rma(sym.bat)
output.gct(sym.eset, paste(studyName, "NormalizedExpression", sep = "_"))
output.cls(pData(sym.eset), "group", paste(studyName, "ClassLables", sep = "_"))

setwd(celDir)
for (i in phenotype$supplementary_file){
  cmd<-paste("cp", i, "ZhangCels")
  system(cmd)
}

# Increment outdir number if outliers found on first run eg. Quality2 etc...
# Un-Comment if 
# arrayQualityMetrics(expressionset = sym.eset,
#                     outdir = paste(studyName, "Quality",sep="_" ),
#                     intgroup = "group",
#                     force = T)


# Run Limma Analysis
source("/home/adam/Desktop/Masters_Thesis_2017/Never_Late_Pipeline/DataProcessingScripts/Limma_DAVID_Analysis.r")
