# Process Demographics and assign samples to groups
# Setup
wd<-"/home/adam/Desktop/Masters_Thesis_2017/Never_Late_Pipeline/Zhang_Analysis/Zhang_Full_QC"
celDir<-"/mnt/Storage/adam/kgsCelFiles/"
softDir<-"/mnt/Storage/adam/kgsSoftFiles"
studyName<-"Zhang"
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
library(ArrayTools)
source("/home/adam/Desktop/Masters_Thesis_2017/Never_Late_Pipeline/DataProcessingScripts/KeyValueExtract.r")
options(download.file.method.GEOquery = "libcurl") # Required to properly fetch GPLXX platforms

# Read covariates from GSE
#seriesMatrix<-getGEO("GSE12093", destdir = softDir)[[1]]     # GEO gene expression data set
setwd(softDir)   # GEO gene expression data set
seriesMatrix<-getGEO(filename = "GSE12093_series_matrix.txt.gz")
setwd(wd)
# Process GEO Demographics Data from the soft file and assign samples to categories
phenotype<-keysToTable(pData(seriesMatrix)[,grep("characteristics_ch1",names(pData(seriesMatrix)))])

# Apply standardized column names, edit the search pattern in the 'grep' statement. 
names(phenotype)[grep("DFS.time", names(phenotype), fixed =T)] <-"DMFS_TIME"
names(phenotype)[grep("DFS.status", names(phenotype), fixed =T)] <-"EVENT_DMFS"
phenotype$Patient_ID<-sub('breast tumor sample ', '',as.character(pData(seriesMatrix)[row.names(phenotype),"title"]), fixed=T)

# Convert time and event columns to numeric, and add grouping column
phenotype$DMFS_TIME<-as.numeric(phenotype$DMFS_TIME)
phenotype$EVENT_DMFS<-as.numeric(phenotype$EVENT_DMFS)
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
  phenotype$DMFS_TIME > 60 &
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



# Drop samples that was excluded in GDOC analysis.
outlier<-c("GSM305135", "GSM305264", "GSM305151", "GSM305172", "GSM305226", "GSM305203", "GSM305212")
phenotype<-phenotype[!(phenotype$geo_accession %in% outlier), ]
phenotype<-phenotype[c("group", names(phenotype)[names(phenotype)!="group"])] # Reorder for arrayQualityMetrics

# Build AffyBatch
sym.bat<-ReadAffy(
  filenames = row.names(phenotype),
  phenoData = phenotype,
  celfile.path = celDir
)
sym.eset<-rma(sym.bat)
output.gct(sym.eset, paste(studyName, "NormalizedExpression", sep = "_"))
output.cls(pData(sym.eset), "group", paste(studyName, "ClassLables", sep = "_"))

# Export Cel Files to a separate directory for Processing on ArrayAnalysis.org
setwd(celDir)
system("mkdir ZhangCels")
for (i in phenotype$supplementary_file){
  cmd<-paste("cp", i, "ZhangCels")
  system(cmd)
}

# Increment outdir number if outliers found on first run eg. Quality2 etc...
# Un-Comment if this is required
arrayQualityMetrics(expressionset = sym.eset,
                  outdir = paste(studyName, "Quality_5",sep="_" ),
                  intgroup = "group",
                  force = T)


# Run Limma Analysis
source("/home/adam/Desktop/Masters_Thesis_2017/Never_Late_Pipeline/DataProcessingScripts/Limma_DAVID_Analysis.r")

