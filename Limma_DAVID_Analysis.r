# Normalize CEL files obtained from GEO
# Setup
setwd(wd)
library(limma)
library(affy)
library(dplyr)
library(VennDiagram)
library(hgu133acdf)
library(hgu133a.db)
library(hgu133plus2.db)
library("annotate")
library(RColorBrewer)
library(RDAVIDWebService)
library(arrayQualityMetrics)
library(gplots)
library(pheatmap)
source("/home/adam/Desktop/Masters_Thesis_2017/Never_Late_Pipeline/DataProcessingScripts/GeneFilters.r")
options(download.file.method.GEOquery = "libcurl")




# Annotate Expression Set
ID<-featureNames(sym.eset)
Symbol<-getSYMBOL(ID,"hgu133a.db")
Name<-as.character(lookUp(ID, "hgu133a.db", "GENENAME"))
Entrez<-as.character(lookUp(ID,"hgu133a.db", "ENTREZID"))

tmp<-data.frame(ID=ID, Symbol=Symbol, Name=Name, Entrez=Entrez, stringsAsFactors = F)
fData(sym.eset)<-tmp


normTitle = paste(studyName, " RMA Normalized Probe Intensities", sep = " ")
palette(brewer.pal(3,'Set2'))         # RColorBrewer provides ~25 pre defined palettes
png(filename = gsub(" ", "_",paste(normTitle, ".png", sep="")),   # Opens a PNG device that writes to the named file
    width = 8.8,
    height = 5.5,
    units = "in",
    res = 96                          # In pixels per inch, determines maximum scaling range (and file size!)
    )
par(mfrow=c(1,1),                     # mfrow specifies the number / arrangement of plot panels in the device
    mar=c(9,5,3,2),                   # Margins (clearly  not set in inches. distance from device border to blot)
    cex =1                            # Global text scaling factor.
    )
boxplot(
  exprs(sym.eset),                    # Data to plot, arranged in columns by sample.
  # mat[20:50, 1:50],                  
  las = 2,                            # Force vertical orientation of axis index.
  col=pData(sym.eset)$group,          # group assignment determines box color: factor level = palette item
  main = normTitle,                   # Main title
  # xlab = "Sample",                    # X axis Label
  ylab = "Log 2 Probe Intensity",     # Y axis label
  cex.main = 1.5,                     # Magnify size of main title text
  cex.lab = 1.5,                      # Magnify size of axis label text
  xaxt="n"                            # Suppress automatic x axis index
)


fill=as.factor(levels(pData(sym.eset)$group))
legend('bottomright',                 # Inset only works with named positions
       legend=c("LATE", "NEVER"),     # Legend items 
       horiz = T,                     # Display Legend Items horizontally
       fill = fill,                   # Pass the levels of "group" to fill as a factor
       xpd =T,                        # When TRUE, allows plot items to display outside plot margins. 
       inset=c(0.38,-0.2)             # Off set legend position by (x,y) in proportion to device size?
)
              
dev.off()

# Filter for A Chip Probes.
nrow(exprs(sym.eset))
old.eset<-sym.eset                                  # For validation, 
aProbes<-mappedkeys(hgu133aACCNUM)
missing<-setdiff(aProbes, row.names(exprs(sym.eset)))
Probes<-setdiff(aProbes, missing)
exprs(sym.eset)<-exprs(sym.eset)[Probes,]
fData(sym.eset)<-fData(sym.eset)[Probes,]


# Verify that no probe mangling has occured. Check rows in the original matrix against rows in the new oen
if(exists('badProbe')){rm('badProbe')}
badProbe<-c()
for( i in 1:nrow(exprs(sym.eset))){
  probe<-row.names(exprs(sym.eset))[i]
  if( any(exprs(sym.eset)[probe,] - exprs(old.eset)[probe,] != 0)){
    msg = paste(probe, sum(exprs(sym.eset)[probe,] - exprs(old.eset)[probe,]))
    print(msg)
    badProbe<-c(badProbe, probe)
  }
}

if(length(badProbe)>0){
  print("Error found after filtering for A chip probes")
}
rm('old.eset')



# Calculate Differential Expression and filter probes for fold change, raw p-value
# Code adapted from GEO2R implementation
design <- model.matrix(~ group + 0, sym.eset)
colnames(design) <- levels(sym.eset$group)
fit <- lmFit(sym.eset, design)
cont.matrix <- makeContrasts(NEVER - LATE, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
top.sym<-topTable(fit2, coef = 1, lfc = log2(1.5), n = Inf )
top.sym<-top.sym[top.sym$P.Value < 0.05,]


# QC Sample assignments. 
qc<-as.data.frame(cbind(row.names(design), 
                        phenotype$supplementary_file, 
                        colnames(exprs(sym.eset)), 
                        as.character(phenotype$group), 
                        design[,1:2])
                  )
names(qc)<-c("Names_Design", "Name_Phenotype", "Name_Expression_Matrix", "Group_Assignment", "LATE Design", "NEVER Design")
write.table(
  qc,
  file = paste(studyName, 'GroupAssignment_QC.txt', sep = "_"),
  sep = "\t",
  quote = F, 
  row.names = F
)

# Verify Differential Expression Directionality and T-test results. 
myTt<-function(x, g1, g2){
  t.test(x[g1], x[g2], paired = F, var.equal = T)$statistic
}

myTp<-function(x, g1, g2){
  t.test(x[g1], x[g2], paired = F, var.equal = T)$p.value
}

myFc<-function(x, g1, g2){
  mean(x[g1], na.rm = T) - mean(x[g2], na.rm=T)
}

late<-as.character(pData(sym.eset)[pData(sym.eset)$group == "LATE", "supplementary_file"])
never<-as.character(pData(sym.eset)[pData(sym.eset)$group == "NEVER", "supplementary_file"])

tstat.ord<- fit2$coef / fit2$stdev.unscaled / fit2$sigma
top.sym$t.ord<-sapply(row.names(top.sym), function(x){tstat.ord[x,1]})
top.sym$myTt<-apply(exprs(sym.eset)[row.names(top.sym),], 1, myTt, never, late)
top.sym$myTp<-apply(exprs(sym.eset)[row.names(top.sym),], 1, myTp, never, late)
top.sym$myFc<-apply(exprs(sym.eset)[row.names(top.sym),], 1, myFc, never, late)
top.sym$LateAvg<-apply(exprs(sym.eset)[row.names(top.sym), late], 1, mean)
top.sym$NeverAvg<-apply(exprs(sym.eset)[row.names(top.sym), never], 1, mean)

# Get Affymetrix_QC Probes (if any)
qcProbes<-top.sym[grep("AFFX", top.sym$ID),]
qcFileName<-paste(studyName, 'QC_Probes.txt', sep='_')
write.table(qcProbes, file=qcFileName, sep = "\t", row.names = F)

# Merge in Genes identified in DAVID
david<-DAVIDWebService$new(email='abf@udel.edu', 
                           url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
result<-addList(david,top.sym$ID,idType="AFFYMETRIX_3PRIME_IVT_ID", listName=studyName, listType="Gene")
davidGenes<-getGeneListReport(david)
checkDavidMerge<-top.sym[,c('ID', 'logFC')]
top.sym<-merge(x=top.sym, y=davidGenes, by="ID", all.x=TRUE)
names(top.sym)[grep("Name.x", names(top.sym), fixed =T)] <-"Name_hgu133a"
names(top.sym)[grep("Name.y", names(top.sym), fixed =T)] <-"Name_DAVID"
checkDavidMerge$new.logFC<-sapply(checkDavidMerge$ID, function(x){top.sym[top.sym$ID == x, 'logFC']})
cor(checkDavidMerge$logFC, checkDavidMerge$new.logFC)


# Add DAVID Symbols
symExtract<-function(x){
  f<-gregexpr('(', x, fixed=T)
  l<-gregexpr(')', x, fixed=T)
  s<-substr(x, 
            f[[1]][length(f[[1]])]+1,
            l[[1]][length(l[[1]])]-1
            )
  s
}
top.sym$DAVID_Symbol<-sapply(top.sym$Name_DAVID, symExtract)


# Filter out bi-directional genes
top.filt<-simpleDupFilter(top.sym, idCol = "ID", symCol = "Name_DAVID", fcCol="logFC")
dupGenes<-getDups(top.sym)

# Get Number of genes meeting various criteria
degSummary<-data.frame(
  Total_Probes = nrow(top.sym),
  Total_Genes= nrow(top.filt),
  Up_Count = nrow(top.filt %>% filter(logFC > 0)),
  Down_Count = nrow(top.filt %>% filter(logFC < 0)),
  Two_Fold_Up = nrow(top.filt %>% filter(logFC > 1)),
  Two_Fold_Down = nrow(top.filt %>% filter(logFC < -1)),
  min_P_Value = min(top.filt$P.Value),
  min_FDR = min(top.filt$adj.P.Val),
  checkFC = cor(top.filt$logFC, top.filt$myFc),
  NumConflict = nrow(dupGenes[dupGenes$Count_UP >0 & dupGenes$Count_DN >0,])
)

write.table(
  degSummary, 
  file = paste(studyName, "DEG_Summary.txt", sep = "_"),
  sep = "\t",
  row.names = F,
  quote = F
)

write.table(
  dupGenes,
  file = paste(studyName, "Multiple_Probes.txt", sep = "_"),
  sep = "\t",
  row.names = F,
  quote = F
)



# Write top "Up" genes -- by Fold Change
nUp<-degSummary$Up_Count
nDn<-degSummary$Down_Count
write.table(
  top.filt %>%
    filter(logFC > 0) %>%
    arrange(desc(logFC)) %>%
    dplyr::select(ID, Symbol, Name = Name_DAVID, logFC, P.Value, adj.P.Val) %>%  #select masked by 'annotate'
    top_n(min(nUp, 20)), 
  file = paste(studyName, "_Top_UP.txt", sep = ""),
  sep = "\t",
  row.names = F,
  quote = F
)


# Write top "Down" genes
write.table(
  top.filt %>%
    filter(logFC < 0) %>%
    arrange(desc(logFC)) %>%
    dplyr::select(ID, Symbol, Name = Name_DAVID, logFC, P.Value, adj.P.Val) %>%
    top_n(min(nDn, 20)), 
  file = paste(studyName, "_Top_DN.txt", sep = ""),
  sep = "\t",
  row.names = F,
  quote = F
)

# Write David Genes, all genes
degFileName<-paste(studyName, "Never_vs_Late_DEG.txt", sep="_")
write.table(top.filt[order(top.filt$logFC),], file=degFileName, row.names = F, sep ="\t", quote = F)
write.table(top.sym[is.na(top.sym$DAVID_Symbol),],"NoGene.txt", row.names = F, sep ="\t", quote = F)
write.table(top.sym[order(top.sym$logFC),], file=paste(studyName, "NvL_All_DE_Probes.txt"), row.names = F, sep ="\t", quote = F)

# Draw Heatmap of 2 fold Genes

mat<-exprs(sym.eset)[top.filt[abs(top.filt$logFC) > log2(1.5),"ID"],]
colnames(mat)<-paste(phenotype$group, phenotype$Patient_ID)
row.names(mat)<-top.filt[top.filt$ID %in% row.names(mat), "DAVID_Symbol"]
hmTitle <-paste(studyName, "2 Fold Differentially Expressed Genes")
annotation_col<-phenotype['group']
row.names(annotation_col)<-colnames(mat)
png(filename = gsub(" ", "_",paste(hmTitle, ".png", sep="")),
  width =12,
  height = 0.15*nrow(mat),
  units = 'in',
  res = 96)
# x11(width = 9, height = 6)
par(mar=c(20,5,4,1))
# heatmap.2(mat,
#         ColSideColors = ifelse(phenotype$group =="LATE", "green", "red"),
#         main = hmTitle,
#         key = T, dendrogram = "column",
#         margins = c(6,8)
#         )
pheatmap(mat, annotation_col = annotation_col)
dev.off()

# Get DAVID Functional Annotation Charts
bpLevels<-c("FAT", "DIRECT", "ALL", as.character(1:5))
davidCat<-c(paste("GOTERM_BP", bpLevels, sep = "_"), "KEGG_PATHWAY", "REACTOME_PATHWAY")
setAnnotationCategories(david, categories =  davidCat)
getFunctionalAnnotationChartFile(david, fileName = paste(studyName, "GOBP_ALL.txt", sep = "_"))


# Save data objects
deg = topTable(fit2, coef = 1, n=Inf)
assign(paste(studyName, "deg", sep ="."), deg)
assign(paste(studyName, "eset", sep="."), sym.eset)
save(
  list = c(
    paste(studyName, "deg", sep ="."),
    paste(studyName, "eset", sep=".")
  ),
  file = paste(studyName, "RMA_Limma.Rda", sep = "_")
)


