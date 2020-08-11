library(plyr)
library(dplyr)
library(ggplot2)
library(plotly)
library(boot)
library(pscl)
library(RColorBrewer)
library(gplots)

##############################
########### STEP 1 ###########
##############################

### import Alex extended database "KEGG pathways" with all metabolic pathways ###
### select the first column and make a list ###
allpathways <- KEGG_pathways_alex_tab
allpathwayslist <- list(allpathways[,1])

## make a list without the rownames (first column), than delete the missings##
for(i in 1:nrow(allpathways)){
  v <- allpathways[i,2:ncol(allpathways)]
  allpathwayslist[i] <- list(v[!is.na(v)])                
}

## place back the rownames ##
names(allpathwayslist) <- rownames(allpathways)

##############################
########### STEP 2 ###########
##############################

### import all metabolites that were identified in 200HIV (1986 metabolites in total) ###
hivmet <- `200hiv_metab`
hivmet <- as.matrix(hivmet)

### make an empty matrix ###
### the number of rows is the same as the number of rows with a KEGG-id (no empty second column) ###
### the number of columns is the same as the number of columns in the original table ###
hivmetkegempty <- matrix(NA, length(which(hivmet[,2]!="")), ncol(hivmet))

### remove metabolites without KEGG-id ###
### 827 metabolites with at least 1 KEGGid remain ###
j <- 1
for(i in 1:nrow(hivmet)){
  if(hivmet[i,2]!=""){  ### if the second column of a row is not empty ###
    hivmetkegempty[j,] <- hivmet[i,]  ### than place this row of position j in the empty matrix ###
    j <- j + 1    ### j = (1, 1+1, 2+1, enz) ###
  }}
hivmetkeg <- hivmetkegempty

row.names(hivmetkeg) <- hivmetkeg[,1] ### make rownames out of column one ###
hivmetkeg <- hivmetkeg[,-1]   ### than delete column 1 ###

#### Check if there are KEGG-ids without a role in human metabolic pathways ####
#### Create thee function "idinpathway" for this purpose ####
#### m is the matrix with KEGG-ids for the corresponding m/z ratios; ####
#### l is the list of metabolic pathways with corresponding KEGGs ####
idinpathway <- function(m, l){       
  mat <- matrix(NA, nrow(m), ncol(m))     ### create an empty matrix with an identical number of rows and columns as hivmetkeg (= m) ###
  u <- unique(unlist(l))                 ### create a vector with only unique elements of allpathwayslist (= l) ###
  for(i in 1:nrow(m)){                   ### for each row of hivmetkeg###
    for(j in 1:ncol(m)){                  ### for each column of hivmetkeg ###
      if(m[i,j]!="" && m[i,j] %in% u)   ### if the cell is not empty and is included in the vector allpathways ###
        mat[i,j] <- m[i,j]              ### than add to the empty matrix  ###
    }}
  return(mat)
}

hivpathkeg <- idinpathway(hivmetkeg, allpathwayslist) ### execute function for m= hivmetkeg and l=allpathwayslist###
row.names(hivpathkeg) <- row.names(hivmetkeg) ### add rownames ###

### delete rows with only NAs, pathways without a metabolite with a KEGG-id in 200HIV will be deleted    ##
onlynalines <- hivpathkeg %>% is.na() %>% apply(MARGIN = 1, FUN = all) ### make a vector for all rows with only NAs ###
hpknonaline <- hivpathkeg[!onlynalines,]   ### filter out all rows with only NAs ###

### delete columns with only NAs, metabolites without a KEGG-id in 1 of the pathways will be deleted ###
onlynacol <- hpknonaline %>% is.na() %>% apply(MARGIN = 2, FUN = all) ### make a vector for all columns with only NAs ###
hivpathkegcompl <- hpknonaline[, !onlynacol]   ### filter out all columns with only NAs ###

######## create hivpathkegcompl: a vector with only 1 kegg-id per metabolite ##############
vecthiv <- c(hivpathkegcompl[,1])
names(vecthiv) <- rownames(hivpathkegcompl)
for(i in 1:length(vecthiv)){
  if(is.na(vecthiv[i])==T){
    f <- which(!is.na(hivpathkegcompl[i,]))
    vecthiv[i] <- hivpathkegcompl[i,f[1]]
  }
} 
vecthiv <- as.data.frame(vecthiv, row.names=names(vecthiv))

vecthiv[,1] <- as.character(vecthiv[ ,1])

##############################
########### STEP 3 ###########
##############################

# 3.1 ### Import data with metabolite measurements, keep only metabolites that are included in human pathways  ####
metabmeasure <- read.table("~/Google Drive/PhD/CCR5 paper/CCR5 18 maart 2020/metaboloom_200hiv_R_data_en_script/200HIV_DESeq2_norm_13.09.2018.txt", sep="\t", header = T, row.names = 1, quote = "")
hpmetabmeasure <- matrix(NA, nrow(hivpathkegcompl), ncol(metabmeasure)) ### make an empty matrix ###
hpmetabmeasure <- metabmeasure[rownames(hivpathkegcompl),] ### import the rows from hivdesqnorm with rownames that correspond with the metabolite measurements ####
hpmetabmeasure <- as.matrix(hpmetabmeasure)
###filter out all columns with only NAs (429 to 394 rows) ###
thpmetabmeasure <- t(hpmetabmeasure)
thpmetabmeasure <- thpmetabmeasure[ ,colSums(is.na(thpmetabmeasure)) != nrow(thpmetabmeasure)]

## I had a problem with importing the metabolites measurements data (problem: all numerical data were characters, with dots being thousand seperators) so i had to execute the following commands to fix this:
thpmetabmeasure2 <- as.data.frame(thpmetabmeasure)
thpmetabmeasure3 <- data.frame(lapply(thpmetabmeasure2, function(x) as.numeric(gsub(",", ".", gsub("\\.", "", x)))))
rownames(thpmetabmeasure3) <- rownames(thpmetabmeasure2)
thpmetabmeasure4 <- as.matrix(thpmetabmeasure3)

### import clinical data ###
clinicaldata <- HIV.dataset3 ## I created a new dataset (HIV.dataset3) with rownames (record IDs) that are similar to the rownames used in other datasets

### make sure that metabolimics data and clinical data have the same number of rows ###
clincommon <- clinicaldata[rownames(thpmetabmeasure),]

### make sure that all patients are both in the metabolomics and the clinical data, i had to execute the following commands ### 
clincommon <- clincommon[rownames(thpmetabmeasure4),]
clincommon <- clincommon[ rowSums(is.na(clincommon)) != ncol(clincommon) , ]
thpmetabmeasure4 <- thpmetabmeasure4[rownames(clincommon),]
thpmetabmeasure4 <- thpmetabmeasure4[ rowSums(is.na(thpmetabmeasure4)) != ncol(clincommon) , ]

##############################################
##########STEP 4 #############################
#### regression analysis for CCR5 data #######
##############################################

outcome <- "Pgated_CCR5P_CD45"
pccr5 <- c()
bccr5 <- c()
confounders <- c("clincommon$Age_at_visit", "clincommon$Sex", "clincommon$BMI")

for(j in 1:ncol(thpmetabmeasure4)){
  yvar <- 'clincommon[,outcome]'
  xvar <- "thpmetabmeasure4[,j]"
  xvars <- confounders
  my.formula <- paste( yvar, '~', xvar, '+', paste(xvars, collapse=' + ' ))
  my.formula <- as.formula(my.formula)
  model <- lm(my.formula)
  bccr5[j] <- summary(model)$coefficients[2,1]                    ###b-waarde: “estimate”####
  pccr5[j] <- summary(model)$coefficients[2,4]                      ###p-waarde: “Pr(>|t|)”####
}
names(pccr5) <- colnames(thpmetabmeasure4)
names(bccr5) <- colnames(thpmetabmeasure4)

sortedpccr5 <- sort(pccr5)
pbelow5ccr5 <- sortedpccr5[sortedpccr5 < 0.05]
bpbelow5ccr5 <- bccr5[names(pbelow5ccr5)]
pfccr5 <- p.adjust(pccr5, method = "fdr") ###FDR correction ###
sortedpfccr5 <- sort(pfccr5)
pfbelow5ccr5 <- sortedpfccr5[sortedpfccr5 < 0.05]

keggspbelow5ccr5 <-c()
for(l in 1:nrow(vecthiv)){
  if(rownames(vecthiv)[l] %in% names(pbelow5ccr5)){
    v <- c(vecthiv[l,1])
    keggspbelow5ccr5 <- c(keggspbelow5ccr5,v)
  }
}
keggspbelow5ccr5 <-  keggspbelow5ccr5[!is.na(keggspbelow5ccr5)]
as.data.frame(keggspbelow5ccr5)

w <- c()
keggsbpbelow5ccr5 <- c()
for(h in 1:length(bpbelow5ccr5)){     ##### b-values of metabolites with p <0.05######
  for(s in 1:nrow(vecthiv)){
    if(names(bpbelow5ccr5[h])==rownames(vecthiv)[s]){
      v <- c(vecthiv[s,])
      v <- v[!is.na(v)]
      w <- rep(bpbelow5ccr5[h], length(v))
      names(w)<- v
      keggsbpbelow5ccr5 <- c(bpbelow5ccr5, w)  ########## add b-values to KEGG-list ####
    }
  }}

a <- c()
keggsppbelow5ccr5 <- c()
for(h in 1:length(pbelow5ccr5)){ ######## p-values of metabolites with p <0.05####
  for(s in 1:nrow(vecthiv)){
    if(names(pbelow5ccr5[h])==rownames(vecthiv)[s]){
      v <- c(vecthiv[s,])
      v <- v[!is.na(v)]
      a <- rep(pbelow5ccr5[h], length(v))
      names(a)<- v
      keggsppbelow5ccr5 <- c(pbelow5ccr5, a)   ####### add p-values to the KEGG list #####
    }
  }}

resultskeggpbccr5<- cbind( vecthiv[ intersect(rownames(vecthiv), names(pbelow5ccr5)), ])
rownames(resultskeggpbccr5) <- intersect(rownames(vecthiv), names(pbelow5ccr5))
resultskeggpbccr5 <- cbind(resultskeggpbccr5, pbelow5ccr5[rownames(resultskeggpbccr5)], bpbelow5ccr5[rownames(resultskeggpbccr5)])
colnames(resultskeggpbccr5) <- c("KEGG", "p-value", "b-value")
resultskeggpbccr5 <- resultskeggpbccr5[order(resultskeggpbccr5[,ncol(resultskeggpbccr5)-1]),]

write.csv(keggspbelow5ccr5)
write.csv(pfbelow5ccr5)
write.csv(resultskeggpbccr5)

