#    File: KeyValueExtract.r
#    Purpose: Extract phenotype data stored as key-value pairs in GEO series matrix files
#    Created: 9/23/17
#    Author:  Adam Faranda

# Get a list of all keys in the data frame
scanKVvars<-function(df, delim=": "){
  kv<-character()
  for (i in 1:ncol(df)){
    kv<-c(kv, sapply(strsplit(as.character(df[,i]), delim), function(x){if(!is.na(x[2])){return(x[1])}else {return()}}))
    kv<-unique(kv)
  }
  kv<-kv[kv !="NULL"]
  kv
}

# Make a data frame of character vectors from a list of variable names
makeCharDF<-function(x){
  df<-data.frame()
  for(i in x){
    newCol<-data.frame(character(), stringsAsFactors = F)
    names(newCol)[1]<-i
    df<-cbind(df,newCol)
  }
  df
}

# Extract delimited values and store them in a new data frame
keysToTable<-function(df, delim=": "){
  keyList<-scanKVvars(df, delim)
  extracted<-makeCharDF(keyList)
  for (i in 1:nrow(df)){
    for(j in 1:ncol(df)){
      for(k in keyList){
        if(grepl(paste(k,delim,sep=""), df[i,j], fixed=T)){
          extracted[i,k]<-sapply(strsplit(as.character(df[i,j]), delim), function(x){x[2]})
        }
      }
    }
    rownames(extracted)[i]<-rownames(df)[i]
  }
  extracted
}
