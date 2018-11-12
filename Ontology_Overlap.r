# Extract Ontology terms detected in 3 out of 4 data sets. 
# Setup
mainDir<-"C:/Users/el3548/OneDrive - DuPont/Never_Late_Pipeline/Never_Late_Pipeline"
wd<-"C:/Users/el3548/OneDrive - DuPont/Never_Late_Pipeline/Never_Late_Pipeline/Ontology_Overlap"
setwd(wd)
library(dplyr)
library(VennDiagram)
library(RDAVIDWebService)
library(ggplot2)
library(venn)
library(reshape2)
library(GOplot)
source("C:/Users/el3548/OneDrive - DuPont/Never_Late_Pipeline/Never_Late_Pipeline/DataProcessingScripts/DAVID_FunctionalAnnotationClusters.r")

# Import Enrichment tables -- Traverse Analysis Directories and assemble a master table of enrichment results. 
getStudyName<-function(fullName){
  Study<-tail(strsplit(fullName, '/', fixed = T)[[1]], 1)
  Study<-sub("_Analysis", "", Study)
  Study
}

termFiles<-data.frame(
  analysisDirs=list.dirs(
    mainDir, recursive = F)[grep("_Analysis", list.dirs(mainDir, recursive = F))], stringsAsFactors = F)

termFiles$Study<-sapply(termFiles[,1], getStudyName)
termFiles$GOFile<-paste(termFiles$Study, "_GOBP_ALL.txt", sep = "")
termFiles$DEGFile<-paste(termFiles$Study, "_Never_vs_Late_DEG.txt", sep="")
termFiles$SummaryFile<-paste(termFiles$Study, "_DEG_Summary.txt", sep = "")
termFiles$FATClusterFile<-paste(termFiles$Study, "_B05_Clusters.txt", sep="")

# Iterate through each analysis directory, and concatenate summary tables
if(exists('summaryTable')){rm(summaryTable)}
for(i in termFiles$analysisDirs){
  if(exists("summaryTable")){
    tmp<-read.table(
      file = paste(i, termFiles[termFiles$analysisDirs == i, "SummaryFile"], sep = "/"),
      header = T, 
      sep = "\t",
      quote = "",
      stringsAsFactors = F
    )
    tmp$Study<-termFiles[termFiles$analysisDirs == i, "Study"]
    summaryTable<-rbind(summaryTable, tmp)
  }
  
  else{
    summaryTable<-read.table(
      file = paste(i, termFiles[termFiles$analysisDirs == i, "SummaryFile"], sep = "/"),
      header = T, 
      sep = "\t",
      quote = "",
      stringsAsFactors = F
    )
    summaryTable$Study<-termFiles[termFiles$analysisDirs == i, "Study"]
  }
}

write.table(
  summaryTable, 
  file = "DegSummary.txt",
  quote = F, 
  sep = "\t",
  row.names = F
)

# Iterate through each analysis directory, and concatenate term-tables
if(exists('termTable')){rm(termTable)}
for(i in termFiles$analysisDirs){
  if(exists("termTable")){
    tmp<-read.table(
    file = paste(i, termFiles[termFiles$analysisDirs == i, "GOFile"], sep = "/"),
    header = T, 
    sep = "\t",
    quote = "",
    stringsAsFactors = F
  )
  tmp$Study<-termFiles[termFiles$analysisDirs == i, "Study"]
  termTable<-rbind(termTable, tmp)
  }

  else{
   termTable<-read.table(
      file = paste(i, termFiles[termFiles$analysisDirs == i, "GOFile"], sep = "/"),
      header = T, 
      sep = "\t",
      quote = "",
      stringsAsFactors = F
    )
   termTable$Study<-termFiles[termFiles$analysisDirs == i, "Study"]
  }
}

# For each category select terms with a Benjamini < 0.05 in at least three studies,convert Study to factor
sigTerms<-termTable[termTable$Benjamini < 0.05,]
sigTerms$Study<-as.factor(sigTerms$Study)
sigTerms$Category<-as.factor(sigTerms$Category)

getHitCount<-function(x, termTable = sigTerms){
  nrow(termTable[termTable$Term == x,])
}
for (i in levels(sigTerms$Category)){
  sigTerms[sigTerms$Category == i, "FoundIn"]<-sapply(sigTerms[sigTerms$Category == i, "Term"], 
                                                      getHitCount, 
                                                      termTable = sigTerms[sigTerms$Category == i,]
                                                      )
}

# Write Terms that hit in all 4 studies, and in 3 of 5 studies
write.csv(
  dcast(sigTerms[sigTerms$FoundIn >3,c("Term","Category","Study", "Benjamini")], Term + Category ~ Study),
  file = "Ontology_Overlap_4_Studies.csv"
)

write.csv(
  dcast(sigTerms[sigTerms$FoundIn >2,c("Term", "Category", "Study", "Benjamini")], Term + Category ~ Study),
  file = "Ontology_Overlap_3_Studies.csv"
)

dcast(sigTerms[sigTerms$FoundIn >2,c("Term", "Category", "Study", "Benjamini")], Term + Category ~ Study)


# Prepare dataset for REVIGO
# Rank Terms by: Median Benjamini, divided by the number of studies that hit. 
scoreFunction<-function(x, termTable = sigTerms){
  y<-termTable[termTable$Term == x, c("Term","Benjamini")]
  if(nrow(y) > 1){
    termScore<-median(y$Benjamini)/nrow(y)
  }
  else{
    termScore<-y$Benjamini
  }
}
attach(sigTerms)
sigTerms[Category == "GOTERM_BP_FAT","Score"]<-sapply(sigTerms[Category == "GOTERM_BP_FAT", "Term"], scoreFunction)
revi.input<-distinct(sigTerms[Category == "GOTERM_BP_FAT" &FoundIn > 2, c("Term", "Score")], Term, Score)
revi.input$Term <- sapply(strsplit(revi.input$Term, '~'), function(x){x[1]})
write.csv(revi.input, "Revigo_Input_Data_2.csv")


# Select Interesting terms for Pathway analysis

interestingTerms<-c(
  "GO:0030334~regulation of cell migration",
  "GO:0048870~cell motility",
  "GO:0012501~programmed cell death",
  "GO:0043549~regulation of kinase activity",
  "GO:0060326~cell chemotaxis",
  "GO:0097190~apoptotic signaling pathway",
  "GO:0006915~apoptotic process",
  "GO:0043410~positive regulation of MAPK cascade",
  "GO:0008285~negative regulation of cell proliferation",
  "GO:0008284~positive regulation of cell proliferation",
  "GO:0098609~cell-cell adhesion",
  "GO:0045785~positive regulation of cell adhesion",
  "GO:0007155~cell adhesion"
)

revi.input$Term[grep("adhesion",revi.input$Term)]
revi.input$Term[grep("proliferation",revi.input$Term)]
revi.input$Term[grep("MAPK",revi.input$Term)]
sigTerms[sigTerms$Term %in% interestingTerms, -6 ]



# Assemble Master Table of differentially expressed genes (all studies) -- degTable
if(exists("degTable")){rm(degTable)}
for(i in termFiles$analysisDirs){
  if(exists("degTable")){
    tmp<-read.table(
      file = paste(i, termFiles[termFiles$analysisDirs == i, "DEGFile"], sep = "/"),
      header = T, 
      sep = "\t",
      stringsAsFactors = F,
      quote = ""
    )
    tmp$Study<-termFiles[termFiles$analysisDirs == i, "Study"]
    degTable<-rbind(degTable, tmp)
  }
  else{
    degTable<-read.table(
      file = paste(i, termFiles[termFiles$analysisDirs == i, "DEGFile"], sep = "/"),
      header = T, 
      sep = "\t",
      stringsAsFactors = F,
      quote = ""
    )
    degTable$Study<-termFiles[termFiles$analysisDirs == i, "Study"]
  }
}

if(exists("degList")){rm(degList)}
for(i in termFiles$analysisDirs){
  if(exists('degList')){
    tmp<-read.table(
      file = paste(i, termFiles[termFiles$analysisDirs == i, "DEGFile"], sep = "/"),
      header = T, 
      sep = "\t",
      stringsAsFactors = F,
      quote = ""
    )
    degList[[termFiles[termFiles$analysisDirs == i, "Study"]]]<-tmp
  }
  else{
    degList = list()
    tmp<-read.table(
      file = paste(i, termFiles[termFiles$analysisDirs == i, "DEGFile"], sep = "/"),
      header = T, 
      sep = "\t",
      stringsAsFactors = F,
      quote = ""
    )
    degList[[termFiles[termFiles$analysisDirs == i, "Study"]]]<-tmp
  }
}





# Look at genes accross studies
bigProbeTable<-dcast(degTable, ID + DAVID_Symbol ~ Study, value.var = "logFC")
bigProbeTable[is.na(bigProbeTable)]<-0
write.table(bigProbeTable, "bigProbeTable.txt", sep = "\t", quote = F, row.names = F)

# Organize probe sets for each of the interesting terms, in each study)
attach(sigTerms)
if(exists("termGenes")){rm(termGenes)}
for(j in levels(sigTerms$Study)){
  for (i in levels(sigTerms$Category)){
    x<-(sigTerms %>% 
          filter(Study == j, Category == i) %>%
          dplyr::select(Genes))$Genes
    names(x)<-(sigTerms %>% 
                 filter(Study == j, Category == i) %>%
                 dplyr::select(Term))$Term
    x<-strsplit(tolower(x), ", ")
    if(exists("y")){
      y<-append(y, list(x))
      names(y)[[length(y)]]<-i
    }
    else{
      y<-list(x)
      names(y)[[length(y)]]<-i
    }
  }
  if(exists("termGenes")){
    termGenes<-append(termGenes, list(y))
    names(termGenes)[[length(termGenes)]]<-j
  }
  else{
    termGenes<-list(y)
    names(termGenes)[[length(termGenes)]]<-j
  }
  rm("y")
}
detach(sigTerms)

# Validate Gene Lists for each Term
# for(i in levels(sigTerms$Study)){
#   for(j in levels(sigTerms$Category)){
#     for(t in unique(sigTerms[Study == i & Category == j, "Term"])){
#       x<-(sigTerms %>% 
#             filter(Study == i, Category == j, Term == t) %>%
#             dplyr::select(Genes))$Genes
#       x<-tolower(strsplit(x, ", ")[[1]])
#       y<-termGenes[[i]][[j]][[t]]
#       print(setdiff(x,y))
#       if(length(setdiff(x,y)) >1){stop("ERROR HERE")}
#     }
#   }
# }


# How many Benjamini < 0.05 terms in each category, By individual and by overlap. 
if(exists('x')){rm(x)}
if(exists('termCountTable')){rm(termCountTable)}
termCountTable<-data.frame(
  Category = levels(sigTerms$Category)
)

# Get Term Counts by Category for individual studies
for(i in levels(sigTerms$Study)){
  for(j in termCountTable$Category){
    if(exists('x')){
      x[length(x)+1]<-nrow(sigTerms[sigTerms$Study == i & sigTerms$Category == j, ])
    }
    else{
      x<-nrow(sigTerms[sigTerms$Study == i & sigTerms$Category == j, ])
    }
  }
  termCountTable<-cbind(termCountTable, x)
  names(termCountTable)[ncol(termCountTable)]<-i
  rm(x)
}

# Get Counts for term overlap, by category. 
termCountTable<-merge(
  termCountTable, as.data.frame(sigTerms %>% 
                                  dplyr::select(Term, Category, FoundIn) %>% 
                                  distinct(Term, Category, FoundIn) %>%
                                  filter(FoundIn > 1) %>% 
                                  arrange(Term) %>% 
                                  arrange(Category) %>%
                                  group_by(Category) %>%
                                  summarise( 
                                    At_Least_Two_Studies = n()
                                  )), 
  by = "Category", 
  all.x = T)

termCountTable<-merge(
  termCountTable, as.data.frame(sigTerms %>% 
                                  dplyr::select(Term, Category, FoundIn) %>% 
                                  distinct(Term, Category, FoundIn) %>%
                                  filter(FoundIn > 2) %>% 
                                  arrange(Term) %>% 
                                  arrange(Category) %>%
                                  group_by(Category) %>%
                                  summarise( 
                                    At_Least_Three_Studies = n()
                                  )), 
  by = "Category", 
  all.x = T)

termCountTable<-merge(
  termCountTable, as.data.frame(sigTerms %>% 
                                  dplyr::select(Term, Category, FoundIn) %>% 
                                  distinct(Term, Category, FoundIn) %>%
                                  filter(FoundIn > 3) %>% 
                                  arrange(Term) %>% 
                                  arrange(Category) %>%
                                  group_by(Category) %>%
                                  summarise( 
                                    All_Studies = n()
                                  )), 
  by = "Category", 
  all.x = T)
termCountTable[is.na(termCountTable)] <-0
write.csv(termCountTable, "TermCounts_by_David_Category.csv")
# Validate termCountTable Construction. 
# attach(sigTerms)
# length(unique(sigTerms[Category == "GOTERM_BP_1"  & FoundIn > 2,"Term" ]))
# length(unique(sigTerms[Category == "GOTERM_BP_2"  & FoundIn > 2,"Term" ]))
# length(unique(sigTerms[Category == "GOTERM_BP_3"  & FoundIn > 2,"Term" ]))
# length(unique(sigTerms[Category == "GOTERM_BP_4"  & FoundIn > 2,"Term" ]))
# length(unique(sigTerms[Category == "GOTERM_BP_5"  & FoundIn > 2,"Term" ]))
# length(unique(sigTerms[Category == "GOTERM_BP_ALL"  & FoundIn > 2,"Term" ]))
# length(unique(sigTerms[Category == "GOTERM_BP_DIRECT"  & FoundIn > 2,"Term" ]))
# length(unique(sigTerms[Category == "GOTERM_BP_FAT"  & FoundIn > 2,"Term" ]))
# 
# length(unique(sigTerms[Category == "GOTERM_BP_1"  & FoundIn > 3,"Term" ]))
# length(unique(sigTerms[Category == "GOTERM_BP_2"  & FoundIn > 3,"Term" ]))
# length(unique(sigTerms[Category == "GOTERM_BP_3"  & FoundIn > 3,"Term" ]))
# length(unique(sigTerms[Category == "GOTERM_BP_4"  & FoundIn > 3,"Term" ]))
# length(unique(sigTerms[Category == "GOTERM_BP_5"  & FoundIn > 3,"Term" ]))
# length(unique(sigTerms[Category == "GOTERM_BP_ALL"  & FoundIn > 3,"Term" ]))
# length(unique(sigTerms[Category == "GOTERM_BP_DIRECT"  & FoundIn > 3,"Term" ]))
# length(unique(sigTerms[Category == "GOTERM_BP_FAT"  & FoundIn > 3,"Term" ]))
# detach(sigTerms)


# Generate Venn Diagram of overlapping genes 
probeList<-sapply(degList, function(x){x[,'DAVID_Symbol']})
png(
  filename = 'Gene_Overlap.png'
)
venn::venn(probeList, zcolor='style', cexil=1.5, cexsn = 1.5)
dev.off()

if(exists('matchTable')){rm(matchTable)}
if(exists('studyList')){rm(studyList)}

matchTable<-data.frame(
  Study1 = character(0),
  Study2 = character(0),
  Pair = character(0),
  Total = numeric(0),
  Match = numeric(0),
  Mismatch = numeric(0)
)


studyList<-termFiles$Study
for(i in termFiles$Study){
  studyList<-studyList[studyList != i]
  for(j in studyList){
    # Get intersection and match table for ith and jth study
    ips<-intersect(degList[[i]][,'DAVID_Symbol'], degList[[j]][,'DAVID_Symbol'])
    x <-degList[[i]][degList[[i]][,'DAVID_Symbol'] %in% ips, c('DAVID_Symbol', 'logFC')]
    y <-degList[[j]][degList[[j]][,'DAVID_Symbol'] %in% ips, c('DAVID_Symbol', 'logFC')]
    z<-merge(x,y,by = 'DAVID_Symbol')
    z$match<-apply(z[,2:3], 1, function(x){abs(sign(x[1]) - sign(x[2]))/2})
    match<-table(z$match)['0']
    miss<-table(z$match)['1']
    tmp<-data.frame(
      Study1 = i,
      Study2 = j,
      Pair = paste(i, " - ", j, sep=""),
      Total = sum(match, miss),
      Match = match,
      Mismatch = miss
    )
    matchTable<-rbind(matchTable, tmp)
  }

}
write.table(matchTable, file = "Table_Match_Missmatch.txt", sep="\t", quote=F, row.names = F)

# Write out Individual Term x Benjamini Lists for Revigo
cat='GOTERM_BP_FAT'
for(i in sigTerms$Study){
  out<-data.frame(
    termId = sapply(
      strsplit(sigTerms[sigTerms$Study == i & sigTerms$Category == cat, 'Term'], '~'),
      function(x){x[1]}
    ),
    PValue = sigTerms[sigTerms$Study == i & sigTerms$Category == cat, 'PValue']
  )
  filename <-paste(i,cat,'RevigoInput.txt', sep='_')
  write.table(out, file = filename, sep ='\t', quote = F, row.names = F)
}

# Write out individual probelists for Manual DAVID Analysis. 
cut<-1.5
for(i in names(degList)){
  filename<-paste(i,"_fc-",cut, "_Affymetrix_Probes.txt", sep="")
  x<-degList[[i]][abs(degList[[i]]$logFC) > log2(cut),c('logFC','ID')]['ID']
  names(x)<-i
  write.table(
    x,
    file = filename, 
    sep = '\t',
    quote = F,
    row.names = F
  )
  print(nrow(x))
}

# Get p-values associated with max Benjamini for each data set. 
cat<-"GOTERM_BP_FAT"
if(exists('maxPtable')){rm(maxPtable)}
maxPtable<-data.frame(
  study<-character(0),
  maxPvalue<-numeric(0),
  maxBenjamini<-numeric(0)
)

for(i in levels(sigTerms$Study) ){
  maxPtable<-rbind(
    maxPtable,
    data.frame(
      study = i,
      maxPvalue = max(sigTerms[sigTerms$Study == i & sigTerms$Category == cat, 'PValue']),
      maxBenjamini = max(sigTerms[sigTerms$Study == i & sigTerms$Category == cat, 'Benjamini'])
    )
  )
}


# Assemble List of cluster files for plotting. 
if(exists('clusterList')){rm(clusterList)}
clusterList<-list()
for(i in termFiles$analysisDirs){
  fn<-paste(i, termFiles[termFiles$analysisDirs == i, 'FATClusterFile'], sep ="/")
  clusterList[[ termFiles[termFiles$analysisDirs == i, 'Study']]] <-DAVIDTermCluster(fn)
}

# Iterate through clusters and generate tables. 
for (study in names(clusterList)){
  clust<-clusterList[[study]]
  fn<-paste(study, "_FAT_Cluster_SummaryTable.txt", sep="")
  write.table(
    clusterSummary(clust, genelist=degList[[study]], study=study, n=3),
    file = fn, 
    quote = F, 
    sep = "\t",
    row.names = F
  )
}

# Iterate through clusters and generate plots
for (study in names(clusterList)){
  clust<-clusterList[[study]]
  plotClusterHeats(clust, genelist = degList[[study]], study = study)
}
# Fix Poorly Plotted Clusters
printGO






#- - - - - - - -  Working Below This Line- - - - - - - - - - - - - - - - - - - #
termsAllStudies<-dcast(sigTerms[sigTerms$FoundIn >2,c("Term","Category","Study", "Benjamini")], Term + Category ~ Study)
termsAllStudies %>% filter(Category == "KEGG_PATHWAY")

termsThreeStudies<-dcast(sigTerms[sigTerms$FoundIn > 2,c("Term","Category","Study", "Benjamini")], Term + Category ~ Study)
termsThreeStudies %>% filter(Category == "GOTERM_BP_FAT")


# Print a Gene-count table of the Reactome Pathways enriched in 2 or more studies
termsTwoStudies<-dcast(sigTerms[sigTerms$FoundIn > 1,c("Term","Category","Study", "Count")], Term + Category ~ Study)
termsTwoStudies[is.na(termsTwoStudies)]<-0
write.csv(
  termsTwoStudies %>% filter(Category == "REACTOME_PATHWAY"),
  "Reactome_in_two_studies.csv"
)




#Get Union of genes annotated with a term in each study.
if(exists('x')){rm(x)}
for(i in levels(sigTerms$Study)){
  if(exists("x")){
    x<-union(x,as.character(termGenes[[i]][['GOTERM_BP_FAT']][['GO:0098609~cell-cell adhesion']]))
  }
  else{
    x<-as.character(termGenes[[i]][['GOTERM_BP_FAT']][['GO:0098609~cell-cell adhesion']])
  }
}



# Get Intersection of genes annotated with a term in each study. 
if(exists('x')){rm(x)}
for(i in c("LoiA", "LoiP")){
  if(exists("x")){
    x<-intersect(x,as.character(termGenes[[i]][['GOTERM_BP_FAT']][['GO:0098609~cell-cell adhesion']]))
  }
  else{
    x<-as.character(termGenes[[i]][['GOTERM_BP_FAT']][['GO:0098609~cell-cell adhesion']])
  }
}

#kegg<-read.csv("Kegg_hsa04514_Cell_Adhesion_Nodes.csv")
kegg<-read.csv("hsa04510_NodeTable.csv")
for (i in levels(sigTerms$Study)){
  cols<-c("DAVID_Symbol", "logFC")
  kegg<-merge(kegg, degTable[degTable$Study == i,cols], by.x = "KEGG_NODE_LABEL_LIST_FIRST", by.y = "DAVID_Symbol", all.x = T)
  names(kegg)[grep('logFC', names(kegg))]<-i
}
kegg[is.na(kegg)]<-0
write.csv(kegg[,c("shared.name", unique(degTable$Study))], "hsa04510_logFoldChanges.csv")

















termsTwoStudies<-dcast(sigTerms[sigTerms$FoundIn > 2,c("Term","Category","Study", "Benjamini")], Term + Category ~ Study)
termsTwoStudies[is.na(termsTwoStudies)]<-0
write.csv(
  termsTwoStudies %>% filter(Category == "GOTERM_BP_FAT"),
  "BP_FAT_in_three_studies.csv"
)




termsTwoStudies %>% filter(Category == "GOTERM_BP_FAT")


min(termTable[
  termTable$Study == "LoiA" &
    termTable$Category == "GOTERM_BP_DIRECT",
  "Benjamini"
  ])

min(termTable[
  termTable$Study == "LoiP" &
    termTable$Category == "GOTERM_BP_DIRECT",
  "Benjamini"
  ])

min(termTable[
  termTable$Study == "Zhang" &
    termTable$Category == "GOTERM_BP_DIRECT",
  "Benjamini"
  ])

min(termTable[
  termTable$Study == "Symmans" &
    termTable$Category == "GOTERM_BP_DIRECT",
  "Benjamini"
  ])




termTable %>% 
  filter(Term == "GO:0030198~extracellular matrix organization") %>%
  dplyr::select(Study, Category, Term, Benjamini)


write.csv(
	dcast(sigTerms[sigTerms$FoundIn >2,c("Term", "Category", "Study", "Benjamini")], Term + Category ~ Study) %>%
		filter(Category == "GOTERM_BP_FAT") %>%
		mutate(GO_ID = sapply(strsplit(Term, '~'), function(x){x[1]})) %>%
		mutate(TermName = sapply(strsplit(Term, '~'), function(x){x[2]})) %>%
		inner_join(revi.input, by=c("GO_ID" = "Term")) %>%
		left_join(TF, by="Term") %>%
		arrange(Score) %>%
		dplyr::select(TermName, Frequency, MinEase, LoiA, LoiP, Symmans, Zhang),
	"Overlap_Terms.csv"
)
	
getwd()





