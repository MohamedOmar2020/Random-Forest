#############################################################################
#############################################################################
## Mohamed Omar
## 10/4/2019
## Goal: Getting and organizing the data for the next step (Creating the K-TSP classifier)
############################################################################
rm(list = ls())

setwd("/Users/mohamedomar/Documents/Research/Projects/New_KTSP")


library(GEOquery)
library(Biobase)
library(sampling)
library(limma)


### Getting the data 
dataset1 <- getGEO("GSE57813", GSEMatrix = TRUE, AnnotGPL = TRUE)
dataset1 <- dataset1$GSE57813_series_matrix.txt.gz

#dataset2 <- getGEO("GSE39282", GSEMatrix = TRUE, AnnotGPL = TRUE)
#dataset2 <- dataset2$`GSE39282-GPL10150_series_matrix.txt.gz`
#dataset3 <- dataset2$`GSE39282-GPL15793_series_matrix.txt.gz`

dataset4 <- getGEO("GSE37817", GSEMatrix = TRUE, AnnotGPL = TRUE)
dataset4 <- dataset4$`GSE37817-GPL6102_series_matrix.txt.gz`
#dataset5 <- dataset4$`GSE37817-GPL8490_series_matrix.txt.gz`

#dataset6 <- getGEO("GSE19915", GSEMatrix = TRUE, AnnotGPL = TRUE)
#dataset6 <- dataset6$`GSE19915-GPL3883_series_matrix.txt.gz`
#dataset7 <- dataset6$`GSE19915-GPL4723_series_matrix.txt.gz`
#dataset8 <- dataset6$`GSE19915-GPL5186_series_matrix.txt.gz`

##################
## Getting old data
### Old objects
load("./Data/Old/allExprsData.rda")
load("./Data/Old/allPhenoData.rda")

### New from GEO
load("./Data/allGEOexprs.rda")
load("./Data/allGEOpheno.rda")


## Expression
expr1 <- exprs(dataset1)
#expr2 <- exprs(dataset2)
#expr3 <- exprs(dataset3)
expr4 <- exprs(dataset4)
#expr5 <- exprs(dataset5)
#expr6 <- exprs(dataset6)
#exprs7 <- exprs(dataset7)
#expr8 <- exprs(dataset8)
expr9 <- read.csv(file = "/Users/mohamedomar/Desktop/UROMOL_gene_fpkm_gtf.txt", sep = "\t")


## Feature data
featData1 <- fData(dataset1)
#featData2 <- fData(dataset2)
#featData3 <- fData(dataset3)
featData4 <- fData(dataset4)
#featData5 <- fData(dataset5)
#featData6 <- fData(dataset6)
#featData7 <- fData(dataset7)
#featData8 <- fData(dataset8)



## Annotation
rownames(expr1) <- featData1$Symbol
summary(is.na(rownames(expr1)))
rownames(expr1) <- gsub("-","", rownames(expr1))
rownames(expr1) <- gsub("_","",rownames(expr1))
sel <- which(apply(expr1, 1, function(x) all(is.finite(x)) ))
expr1 <- expr1[sel, ]
expr1 <- expr1[!is.na(rownames(expr1)),]
dim(expr1)
#rownames(expr2) <- featData2$GENE_SYMBOL
# Remove empty rows (no gene symbol)
#expr2 <- expr2[!(rownames(expr2) == ""), ]

#rownames(expr3) <- featData3$Gene
#expr3 <- expr3[!(rownames(expr3) == "-"), ]

rownames(expr4) <- featData4$`Gene symbol`
expr4 <- expr4[!(rownames(expr4) == ""), ]
rownames(expr4) <- gsub("-","", rownames(expr4))
rownames(expr4) <- gsub("_","",rownames(expr4))


#rownames(expr5) <- featData5$Symbol
#expr5 <- expr5[!(rownames(expr5) == ""), ]
# Remove NAs
#sel <- which(apply(expr5, 1, function(x) all(is.finite(x)) ))
#expr5 <- expr5[sel, ]

#rownames(expr6) <- featData6$`Gene symbol`
#expr6 <- expr6[!(rownames(expr6) == ""), ]
# Remove NAs
#sel <- which(apply(expr6, 1, function(x) all(is.finite(x)) ))
#expr6 <- expr6[sel, ]
#rownames(expr6) <- gsub("-","", rownames(expr6))
#rownames(expr6) <- gsub("_","",rownames(expr6))
#rownames(expr6) <- gsub(".+\\///","",rownames(expr6))

#rownames(expr8) <- featData8$V3.0.3_gene_symbol
#expr8 <- expr8[!(rownames(expr8) == ""), ]
# Remove NAs
#sel <- which(apply(expr8, 1, function(x) all(is.finite(x)) ))
#expr8 <- expr8[sel, ]
#rownames(expr8) <- gsub("-","", rownames(expr8))
#rownames(expr8) <- gsub("_","",rownames(expr8))

expr9 <- expr9[!duplicated(expr9$gene.name), ]
rownames(expr9) <- expr9$gene.name
expr9 <- expr9[,-c(1:4)]
table(is.na(expr9))
expr9 <- expr9[!(rownames(expr9) == ""), ]
expr9 <- expr9[!duplicated(gsub("-","", rownames(expr9))),]
rownames(expr9) <- gsub("_","",rownames(expr9))
# use voom on expr9
expr9 <- voom(expr9)
expr9 <- expr9$E

### Phenotype
pheno1 <- pData(dataset1)
#pheno2 <- pData(dataset2)
#pheno3 <- pData(dataset3)
pheno4 <- pData(dataset4)
#pheno5 <- pData(dataset5)
#pheno6 <- pData(dataset6)
#pheno8 <- pData(dataset8)
pheno9 <- read.csv(file = "/Users/mohamedomar/Desktop/E-MTAB-4321.sdrf.txt", sep = "\t")

### Modify the pheno

pheno1$PROGRESSION <- as.factor(pheno1$`liver sample group:ch1`)
levels(pheno1$PROGRESSION) <- c("NoProgression", "Progression")

# pheno2/pheno3 ## Not good (GSE39282)

pheno4$PROGRESSION <- as.factor(pheno4$`progression:ch1`)
levels(pheno4$PROGRESSION) <- c("NoProgression", "Progression")
# Remove normal samples
pheno4 <- pheno4[-c(1:6), ]
# Modify expr4
expr4 <- expr4[, colnames(expr4) %in% rownames(pheno4)]
all(rownames(pheno4) == colnames(expr4))

# pheno5 # Not good

##Remove Ta from pheno6
#pheno6 <- pheno6[which(!grepl("Ta", pheno6$`tumor stage:ch1`)), ]
#pheno6$PROGRESSION <- as.factor(pheno6$`tat1 progress^to minv/cystectomy/dead:ch1`)
#levels(pheno6$PROGRESSION) <- c("NoProgression", "Progression", "Progression", "Progression", "Progression")
#pheno6 <- pheno6[!is.na(pheno6$PROGRESSION), ]
# Modify expr6
#expr6 <- expr6[, colnames(expr6) %in% rownames(pheno6)]
#all(rownames(pheno6) == colnames(expr6))


## Remove Ta from pheno8
#pheno8 <- pheno8[which(!grepl("Ta", pheno8$`tumor stage:ch1`)), ]
#pheno8$PROGRESSION <- as.factor(pheno8$`tat1 progress^to minv/cystectomy/dead:ch1`)
#levels(pheno8$PROGRESSION) <- c("NoProgression", "Progression", "Progression", "Progression")
#pheno8 <- pheno8[!is.na(pheno8$PROGRESSION), ]
# Modify expr8
#expr8 <- expr8[, colnames(expr8) %in% rownames(pheno8)]
#all(rownames(pheno8) == colnames(expr8))


## Remove Ta from pheno9
pheno9 <- pheno9[which(!grepl("Ta", pheno9$Factor.Value.disease.staging.)), ]
pheno9$PROGRESSION <- as.factor(pheno9$Factor.Value.progression.to.T2..)
levels(pheno9$PROGRESSION) <- c("NoProgression", "Progression")
rownames(pheno9) <- pheno9$Source.Name
# Modify expr9
expr9 <- expr9[,colnames(expr9) %in% rownames(pheno9)]
all(rownames(pheno9) == colnames(expr9))


############

allpheno1 <- list(pheno1, pheno4, pheno9)
names(allpheno1) <- c("GSE57813", "GSE37817-GPL6102", "E-MTAB-4321")

allExpr1 <- list(expr1, expr4, expr9)
names(allExpr1) <- c("GSE57813", "GSE37817-GPL6102", "E-MTAB-4321")


### Filter phenotype information for the required samples
groupProgression1 <- mapply(x=allpheno1, FUN=function(x) {
  x <- x[,"PROGRESSION"]
  out <- factor(x, levels=c("Progression", "NoProgression"))
  out
})
###################################################################################
## Process old data

### Combine expression and phenotypes (OLD DATA)
allExprs2 <- c(allExprsData, allGEOexprs)
allPheno2 <- c(allPhenoData, allGEOpheno)


### Prepare objects
allPheno2 <- lapply(allPheno2, function(x) {
    colnames(x) <- toupper(colnames(x))
    colnames(x) <- gsub("\\.$", "", gsub("\\.\\.", ".", colnames(x)))
    colnames(x) <- gsub("\\.", "_", colnames(x))
    x
})

### Retain only datasets with survival in caucasians
keep <- c("GSE13507", "pmid15930337", "GSE32894")

### Drop unwanted datasets and format survival time
allPheno2 <- allPheno2[keep]
allExprs2 <- allExprs2[keep]

### Process information for consistency
allPheno2 <- lapply(allPheno2, function(x) {
    x[ , grep("_MO", colnames(x)) ] <- sapply(x[ , grep("_MO", colnames(x)) ], as.numeric)
    x[ , grep("DATEFACT", colnames(x)) ] <- sapply(x[ , grep("DATEFACT", colnames(x)) ], as.numeric)
    x[ , sapply(x, class) == "factor" ] <- sapply(x[ , sapply(x, class) == "factor"], as.character)
    x
})

### Make IDs as rownames
rownames(allPheno2$GSE13507) <- allPheno2$GSE13507$GSM
rownames(allPheno2$pmid15930337) <- allPheno2$pmid15930337$IDS
rownames(allPheno2$GSE32894) <- allPheno2$GSE3289$GEO_ACCESSION

###########################################################################
### Keep only relevant samples
keepSamples <- sapply(allPheno2, function(x) {
  if (length(grep("STAGE$", colnames(x))) > 0) {
    out <- x[ , grep("STAGE$", colnames(x))] %in% c("T1", "STAGE:T1N0M0") &
      !is.na(x[ , grep("STAGE$", colnames(x)) ]) &
      !is.na(x[ , grep("PROGR", colnames(x)) ])
  } else {
    out <- rep(TRUE, nrow(x))
  }
  out
})


### Filter phenotype information for the required samples
groupProgression2 <- mapply(x=allPheno2, y=keepSamples, FUN=function(x, y) {
  nms <- rownames(x)[y]
  x <- x[y, grep("PROGRESSION", colnames(x))]
  names(x) <- nms
  x[grep("YES", x)] <- "YES"
  x[x %in% c("1", "NO", "NonProg")] <- "NoProgression"
  x[x %in% c("2", "YES", "Prog")] <- "Progression"
  out <- factor(x, levels=c("Progression", "NoProgression"))
  out
})

### Modify allExprs2 to include only the relevant samples
allExprs2$GSE13507 <- allExprs2$GSE13507[,colnames(allExprs2$GSE13507)%in% names(groupProgression2$GSE13507)]
allExprs2$pmid15930337 <- allExprs2$pmid15930337[, colnames(allExprs2$pmid15930337)%in% names(groupProgression2$pmid15930337)]
allExprs2$GSE32894 <- allExprs2$GSE32894[, colnames(allExprs2$GSE32894) %in% names(groupProgression2$GSE32894)]

###################################################################
### Combine old and new
allPheno <- c(allPheno2, allpheno1)
allExpr <- c(allExprs2, allExpr1)
groupProgression_ALL <- c(groupProgression2, groupProgression1)



### Keep good expression data and phenotypes
allExpr <- allExpr[ names(groupProgression_ALL) ]

### Find commom subset of genes
commonGenes <- Reduce("intersect", lapply(allExpr, rownames))


### Filter expression for the required samples
exprsProgression <- mapply(x=allExpr, FUN=function(x, gns) {
  x <- x[ gns ,]
}, MoreArgs=list(gns=commonGenes))

### Check
all(names(exprsProgression) == names(groupProgression_ALL))

### Check order
all(rownames(allPheno$GSE57813) == colnames(allExpr$GSE57813))
all(rownames(allPheno$`GSE37817-GPL6102`) == colnames(allExpr$`GSE37817-GPL6102`))
all(rownames(allPheno$`E-MTAB-4321`) == colnames(allExpr$`E-MTAB-4321`))


###################################################################
#############
## Combine for KTSP analysis
train <- c("GSE13507","pmid15930337","GSE32894")
test <- c("GSE57813", "GSE37817-GPL6102", "E-MTAB-4321")

## Training
trainMat <- do.call("cbind", exprsProgression[train])
trainGroup <- factor(do.call("c", groupProgression_ALL[train]))
levels(trainGroup) <- c("Progression", "NoProgression")

## Testing
testMat <- do.call("cbind", exprsProgression[test])
testGroup <- factor(do.call("c", groupProgression_ALL[test]))
levels(testGroup) <- c("Progression", "NoProgression")

##################################################################
#####################
## All combined
allMat <- do.call("cbind", exprsProgression)
allGroup <- unlist(groupProgression_ALL)
allStudies <- gsub("\\..+", "", names(allGroup))

all(colnames(allMat) == names(allStudies))


### Prepare to extract additional info
myPheno1 <- mapply(function(x, y) x[y,], allPheno[1:3], keepSamples)
myPheno <- c(myPheno1, allpheno1)

#############################################################
### Grade

## allpheno$GSE57813$Grade <- allpheno$GSE57813$.. ##No grade
myPheno$`GSE37817-GPL6102`$GRADE <- as.factor(myPheno$`GSE37817-GPL6102`$`grade:ch1`)
myPheno$`E-MTAB-4321`$GRADE <- as.factor(myPheno$`E-MTAB-4321`$Factor.Value.tumor.grading.)



### Covariates of relevance select complete cases: GRADE
allGrade <- lapply(myPheno, function(x) {
  i <- grep("GRADE", colnames(x))
  if (length(i) == 0) out <- factor(rep("", nrow(x)))
  else x <- factor(x[, i  ])
})


levels(allGrade$GSE13507) <- c("low", "high")
levels(allGrade$pmid15930337) <- c("low", "high")
levels(allGrade$GSE32894) <- c("low", "high", "high")
levels(allGrade$`GSE37817-GPL6102`) <- c("high", "low")
levels(allGrade$`E-MTAB-4321`) <- c("high", "low")


allGrade <- factor(unlist(allGrade))
levels(allGrade)[levels(allGrade) == "GX"] <- NA
allGrade <- ordered(allGrade, levels=c("low", "high"))

################################################################################
########################
### Recurrence

myPheno$`GSE37817-GPL6102`$RECURRENCE <- as.factor(myPheno$`GSE37817-GPL6102`$`recurrence:ch1`)
levels(myPheno$`GSE37817-GPL6102`$RECURRENCE) <- c("NoRecurrence", "Recurrence")

### Covariates of relevance select complete cases: RECURRENCE
allRecurrence <- lapply(myPheno, function(x) {
  i <- grep("RECURRENCE", colnames(x))
  if (length(i) == 0) out <- factor(rep("", nrow(x)))
  else x <- factor(x[, i  ])
})
levels(allRecurrence$GSE13507) <- c("Recurrence", "NoRecurrence") # Value 1=Rec, value 2=NoRec

allRecurrence <- factor(unlist(allRecurrence))
levels(allRecurrence)[levels(allRecurrence) == ""] <- NA

###############################################################################
############################
### Intravesical therapy

myPheno$`GSE37817-GPL6102`$INTRAVESICAL <- as.factor(myPheno$`GSE37817-GPL6102`$`intravesical therapy:ch1`)
myPheno$`E-MTAB-4321`$INTRAVESICAL <- as.factor(myPheno$`E-MTAB-4321`$Factor.Value.BCG.treatment.)


### Covariates of relevance select complete cases: INTRAVESCICAL
allRX <- lapply(myPheno, function(x) {
  i <- grep("INTRAVESICAL", colnames(x))
  if (length(i) == 0) out <- factor(rep("", nrow(x)))
  else x <- factor(x[, i  ])
})

levels(allRX$GSE13507) <- c("YES", "NO") # Value 1=YES, value 2=NO
levels(allRX$`GSE37817-GPL6102`) <- c("NO", "YES")
levels(allRX$`E-MTAB-4321`) <- c("NO", "YES")

allRX <- factor(unlist(allRX))
levels(allRX)[levels(allRX) == ""] <- NA

#################################################################################
###############################
### Age
myPheno$`GSE37817-GPL6102`$AGE <- as.character(myPheno$`GSE37817-GPL6102`$`AGE:ch1`)
myPheno$`E-MTAB-4321`$AGE <- as.character(myPheno$`E-MTAB-4321`$Factor.Value.age.)

### Covariates of relevance select complete cases: AGE
allAGE <- lapply(myPheno, function(x) {
  i <- grep("^AGE$", colnames(x))
  if (length(i) == 0) out <- rep(NA, nrow(x))
  else x <- as.numeric(x[, i  ])
})
allAGE <- unlist(allAGE)

##################################################################################
################################
### Sex
myPheno$`GSE37817-GPL6102`$GENDER <- as.factor(myPheno$`GSE37817-GPL6102`$`SEX:ch1`)
myPheno$`E-MTAB-4321`$GENDER <- as.factor(myPheno$`E-MTAB-4321`$Factor.Value.sex.)
myPheno$GSE13507$GENDER <- myPheno$GSE13507$SEX

### Covariates of relevance select complete cases: SEX
allGENDER <- lapply(myPheno, function(x) {
  i <- grep("GENDER", colnames(x))
  if (length(i) == 0) out <- rep(NA, nrow(x))
  else x <- factor(x[, i  ])
})
levels(allGENDER$`E-MTAB-4321`) <- c("F", "M")
allGENDER <- factor(unlist(allGENDER))
levels(allGENDER)[levels(allGENDER) == ""] <- NA


#########################################################################

### Assemble in one data.frame and turn numeric
covs <- data.frame(STUDIES=allStudies,
                   GRADE=allGrade,RECURRENCE=allRecurrence, RX=allRX,
                   GENDER=allGENDER, AGE=allAGE)

### Prepare vocs for sampling
covs <- sapply(covs , function(x) as.numeric(factor(paste(x))) )


###########################################################################
###SAMPLING

### Balanced stratification
set.seed(333)
trainingOrTesting <- balancedstratification(
  covs[ , , drop=FALSE], strata=1*(allGroup == "Progression"),
  pik=inclusionprobabilities(1:nrow(covs), nrow(covs)/2),
  comment=TRUE, method=1)

### Show
apply(covs[, -ncol(covs),drop=FALSE], 2, table, allGroup, trainingOrTesting)

### Subset Training
mixTrainMat <- allMat[ , trainingOrTesting == 1]
mixTrainGroup <- allGroup[ trainingOrTesting == 1]
mixTrainStudy <- allStudies[ trainingOrTesting == 1]

### Subset Testing
mixTestMat <- allMat[ , trainingOrTesting == 0]
mixTestGroup <- allGroup[ trainingOrTesting == 0]
mixTestStudy <- allStudies[ trainingOrTesting == 0]


###########################################################################
### Save
save(exprsProgression, trainMat, trainGroup, testMat, testGroup,
     mixTrainMat, mixTrainGroup, mixTrainStudy,
     mixTestMat, mixTestGroup, mixTestStudy,
     file="./objs/progressionDataGood.rda")

#########################################################################
#########################################################################
########################################################################









