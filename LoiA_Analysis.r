# Process Demographics and assign samples to groups
# Setup
wd<-"/home/adam/Desktop/Masters_Thesis_2017/Never_Late_Pipeline/LoiA_Analysis"
celDir<-"/mnt/Storage/adam/kgsCelFiles/"
softDir<-"/mnt/Storage/adam/kgsSoftFiles"
studyName<-"LoiA"
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
library(pheatmap)
library(ArrayTools)
source("/home/adam/Desktop/Masters_Thesis_2017/Never_Late_Pipeline/DataProcessingScripts/KeyValueExtract.r")
source("/home/adam/Desktop/Masters_Thesis_2017/Never_Late_Pipeline/DataProcessingScripts/Covariate_Descriptive_Stats.r")
options(download.file.method.GEOquery = "libcurl") # Required to properly fetch GPLXX platforms

# Fetch  SOFT file from Gene expression omnibus
# seriesMatrix<-getGEO("GSE6532", destdir = softDir)     # GEO gene expression data set
# seriesMatrix<-seriesMatrix[[2]]

# Load SOFT file from local copy
setwd(softDir)
seriesMatrix<-getGEO(filename = "GSE6532-GPL96_series_matrix.txt.gz")
setwd(wd)

# Process GEO Demographics Data from the soft file and assign samples to categories
phenotype<-read.table("GSE6532_LUMINAL_demo.txt", header = T, sep="\t")

# Apply standardized column names, edit the search pattern in the 'grep' statement. 
names(phenotype)[grep("t.dmfs", names(phenotype), fixed =T)] <-"DMFS_TIME"
names(phenotype)[grep("e.dmfs", names(phenotype), fixed =T)] <-"EVENT_DMFS"
names(phenotype)[grep("t.rfs", names(phenotype), fixed =T)] <-"DFS_TIME"
names(phenotype)[grep("e.rfs", names(phenotype), fixed =T)] <-"EVENT_DFS"
names(phenotype)[grep("series", names(phenotype), fixed =T)] <-"LAB"
names(phenotype)[grep("er", names(phenotype), fixed =T)] <-"ER_STATUS"
names(phenotype)[grep("LAB", names(phenotype), fixed =T)] <-"series"
names(phenotype)[grep("samplename", names(phenotype), fixed =T)] <-"Patient_ID"
names(phenotype)[grep("geo_accn_hg.u133a", names(phenotype), fixed =T)] <-"geo_accession"



# Convert metadata columns to appropriate type (numeric or factor), and add grouping column
phenotype$DMFS_TIME<-as.numeric(phenotype$DMFS_TIME)
phenotype$EVENT_DMFS<-as.numeric(phenotype$EVENT_DMFS)
phenotype$DFS_TIME<-as.numeric(phenotype$DFS_TIME)
phenotype$EVENT_DFS<-as.numeric(phenotype$EVENT_DFS)
phenotype$geo_accession<-as.character(phenotype$geo_accession)
phenotype$Patient_ID<-as.character(phenotype$Patient_ID)
phenotype$node<-as.factor(phenotype$node)
phenotype$pgr<-as.factor(phenotype$pgr)
phenotype$age<-as.numeric(phenotype$age)
phenotype$size<-as.numeric(phenotype$size)
phenotype$grade<-as.factor(phenotype$grade)

phenotype$group<-factor(nrow(phenotype),levels=c("LATE", "NEVER"))

# Get last event -- 4906 days
lastEvent<-max(phenotype[phenotype$EVENT_DFS == 1, "DFS_TIME"], na.rm = T)

# Merge in filenames and re-assign row names -- rows from the Plus2 Chip fall away. 
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



# Assign Groups and filter on unused samples
phenotype[
  !is.na(phenotype$DFS_TIME) &
    !is.na(phenotype$DMFS_TIME) &
    !is.na(phenotype$EVENT_DFS) &
    !is.na(phenotype$EVENT_DMFS) &
    !is.na(phenotype$ER_STATUS) &
    phenotype$ER_STATUS == 1 &
    phenotype$EVENT_DFS == 1 &
    phenotype$DFS_TIME > 1825 &
    phenotype$EVENT_DMFS == 1 &
    phenotype$DMFS_TIME > 1825 &
    phenotype$DMFS_TIME <= lastEvent,
  "group"
  ] <-"LATE"

phenotype[
  !is.na(phenotype$DFS_TIME) &
    !is.na(phenotype$DMFS_TIME) &
    !is.na(phenotype$EVENT_DFS) &
    !is.na(phenotype$EVENT_DMFS) &
    !is.na(phenotype$ER_STATUS) &
    phenotype$EVENT_DFS == 0 &
    phenotype$DFS_TIME > lastEvent &
    phenotype$EVENT_DMFS == 0 &
    phenotype$ER_STATUS == 1 &
    phenotype$DMFS_TIME > lastEvent,
  "group"
  ]<-"NEVER"

phenotype<-phenotype[phenotype$group %in% c("LATE", "NEVER"),]


# Drop samples that have been identified as outliers in previous runs.
# Samples were analyzed at two different labs. Four LATE samples analzyed
# at lab "JBI" were excluded
outlier<-c()
phenotype<-phenotype[!(phenotype$geo_accession %in% outlier), ]
phenotype<-phenotype[c("group", names(phenotype)[names(phenotype)!="group"])]


# Prepare Summary Statistics for final data set.
descTable<-statsByType(phenotype)
writeDescTables(descTable)

statsTable<-calculatePairwise(phenotype)
writePairwiseStats(statsTable)

figuresByType(phenotype)

# Build AffyBatch
sym.bat<-ReadAffy(
  filenames = row.names(phenotype),
  phenoData = phenotype,
  celfile.path = celDir
)
sym.eset<-rma(sym.bat)

# Write out normalized expression matrix to a GCT file. 
output.gct(sym.eset, paste(studyName, "NormalizedExpression", sep = "_"))
output.cls(pData(sym.eset), "group", paste(studyName, "ClassLables", sep = "_"))

setwd(celDir)
system("mkdir LoiA")
for (i in phenotype$supplementary_file){
  cmd<-paste("cp", i, "ZhangCels")
  system(cmd)
}

# Increment outdir number if outliers found on first run eg. Quality2 etc...
# Un-Comment if the analysis is to be repeated with outliers removed
# arrayQualityMetrics(expressionset = sym.eset,
#                     outdir = paste(studyName, "Quality1",sep="_" ),
#                     intgroup = "group",
#                     force = T)


# Run Limma Analysis
source("/home/adam/Desktop/Masters_Thesis_2017/Never_Late_Pipeline/DataProcessingScripts/Limma_DAVID_Analysis.r")
