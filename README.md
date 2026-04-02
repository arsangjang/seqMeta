# seqMeta/MetaSKAT
## Performing a survival meta-analysis of rare variants from sequencing data using seqMeta/MetaSKAT

### Instalation SKAT/ MetaSKAT/seqMeta


-   Use [R](https://cran.r-project.org/) Version \< 4.0.2 if you got a instalation "had non-zero exit status" error & and fix C++/g++/ rtools.

### Required packages:

```{r}
setwd("D:/folder")
install_github("brentp/skatMeta")

library(CompQuadForm)
library(twostageGWASsurvival)
library(seqMeta)
library(MASS)
library(writexl)
library(Matrix)
```

### Required data files:

1- SNP information (SNPInfo): contains SNP & gene names (establishing a connection between the genetic loci and their corresponding genes)

[![](SNP data format.PNG)](https://drive.google.com/file/d/1nOs2OyH6q_GUuTImJf9dgHdifO2chzSQ/view?usp=drive_link)

2- Demographic & clinical data file (data): Includes age, sex, Bp, survival time, outcome variable , ...for each individual.

3- Z matrix (Genotype Matrix): This matrix represents the genotype data for the genetic variants. The rows of this matrix correspond to individual samples, while the columns indicate the SNP names (genetic loci). The entries of the matrix represent the alleles present at each SNP for each individual. It is generally a binary matrix, where 0 denotes one allele and 1 denotes the other allele.

[![](Z matrix.PNG)](https://drive.google.com/file/d/10t80BN7hGG4__dXEF_IGA_uQyzgflGpU/view?usp=drive_link)

To create Z we need allele data set contains individuals ID and position:

               allel_binary_data.txt

| Pos | Ref | Alter | sample 1 | sample 2 | sample 3 |
|-----|-----|-------|----------|----------|----------|
| 45  | A   | G     | 1        | 1        | 0        |
| 158 | TG  | G     | 0        | 1        | 0        |
| 86  | T   | C     | 0        | 0        | 0        |

```{r}
genomatrix = read.table("allel_binary_data.txt",header=TRUE)
snpInfo = read.csv('SNPInfo.csv',header=TRUE)
names(genomatrix)

I = match(snpInfo$POS, genomatrix$POS)
genomatrix = genomatrix[I,]
genomatrix = t(genomatrix[,-c(1:3)])
colnames(genomatrix) = snpInfo$SNP
``` 

### Runing survival model

```{r}
cox_model<-prepCox(genomatrix, Surv(time,binary_outcome) ~ Age +  Sex +  intervention, 
             SNPInfo = SNPInfo, data = my_data, verbose = FALSE)

output_cox <- skatMeta(cox_model, SNPInfo = SNPInfo)
head(output_cox)
```


## SKAT Optimal Test

SKAT O Test introduced to test weighted averages of SKAT & burden tests

```{r}
cox_model1<- prepScores(genomatrix_a, outcome ~ sex+age+ intervention, 
             SNPInfo = SNPInfo, data = my_data, verbose = FALSE)

cox_model2<-prepScores(genomatrix_b, outcome ~ sex+age+ intervention, 
             SNPInfo = SNPInfo, data = my_data_b, verbose = FALSE)


Out_combined <- skatOMeta(cox_model1,cox_model2, SNPInfo=SNPInfo, method="int")

head(Out_combined)

```
##Survival data

```{r}
#survival model 1
Cox_model1<-prepCox(genomatrixS, Surv(time,rel) ~ sex+age+ intervention,
                    SNPInfo = SNPInfo, data = my_data, verbose = FALSE)

out.Cox_model1 <- skatOMeta(Cox_model1, SNPInfo = SNPInfo)
head(out.Cox_model1)

#survival model 2
Cox_model2<-prepCox(genomatrix_b, Surv(time,rel) ~   sex+age+ intervention, 
             SNPInfo = SNPInfo, data = my_data_b, verbose = FALSE)

out_Cox_model2 <- skatMeta(Cox_model2, SNPInfo = SNPInfo)
head(out.Cox_model2)


overall_Cox <- skatOMeta(Cox_model1, Cox_model2, SNPInfo = SNPInfo)
head(overall_Cox)

```

## Results from the SKAT test

| gene         | p                   | pmin        | rho      | cmaf                              | nmiss             | nsnps              | errflag             |
|--------|--------|--------|-------|-------------|-----------|---------|-----------|
| Name of Gene | P-value: is \<0.05? | Min P-value | 𝜌 (0, 1) | Cumulative minor allele frequency | n of missing SNPs | n SNPs in the gene | Inaccurate p-values |






rm(list = ls())

## Load libraries
library(tidyverse)
library(seqMeta)
library(CompQuadForm)
library(twostageGWASsurvival)
library(MASS)
library(writexl)
library(readxl)
library(Matrix)
library(rmarkdown)
library(knitr)
library(data.table)
library(gtsummary)
library(dplyr)

# Set directory
setwd("C:/Users/Sarsangjang/Desktop/GVHD Project/Results_SKAT-o_March2026_JCO_reply")


#### 1-a) Prepare data for main analysis (*only one time*)


# Read validation dataset- NCI
genomatrixV = read.table("C:/Users/Sarsangjang/Desktop/GVHD Project/Results_SKAT-o_March2026_JCO_reply/MDS295_2942_binary.txt",header=TRUE)
snpInfoV = read.csv('C:/Users/Sarsangjang/Desktop/GVHD Project/Results_SKAT-o_March2026_JCO_reply/SNV_info_excluded_non-coding.csv',header=TRUE)
pheno_fileV = read_excel("C:/Users/Sarsangjang/Desktop/GVHD Project/Results_SKAT-o_March2026_JCO_reply/MDS295_clinical_R3.xlsx")

pheno_fileV$event_aGvhD24[pheno_fileV$event_aGvhD24==98] <- NA
pheno_fileV$event_aGvhD34[pheno_fileV$event_aGvhD34==98] <- NA
pheno_fileV$cgvhd[pheno_fileV$cgvhd== 99] <- NA


pheno_fileV$event_aGvhD24<-as.numeric(pheno_fileV$event_aGvhD24)
pheno_fileV$event_aGvhD34<-as.numeric(pheno_fileV$event_aGvhD34)
pheno_fileV$intxagvhd24<-as.numeric(pheno_fileV$intxagvhd24)
pheno_fileV$intxagvhd34<-as.numeric(pheno_fileV$intxagvhd34)
pheno_fileV$intxcgvhd<-as.numeric(pheno_fileV$intxcgvhd)
pheno_fileV$cgvhd<-as.numeric(pheno_fileV$cgvhd)
pheno_fileV$prophx<-as.numeric (pheno_fileV$prophx)

# Remove rows where prophx and cgvhd is NA
pheno_fileV <- pheno_fileV[!is.na(pheno_fileV$prophx), ]
pheno_fileV <- pheno_fileV[!is.na(pheno_fileV$cgvhd), ]

nrow(pheno_fileV) # 285

##### 1-c) Create genomatrix for cox model to merge snpInfo and genomatrix info- recipients


IV = match(snpInfoV$POS, genomatrixV$POS)

genomatrixV = genomatrixV[IV,]

genomatrixV = t(genomatrixV[,-c(1:3)])

colnames(genomatrixV) = snpInfoV$SNP


samplesV = row.names(genomatrixV)
IV = match(as.character(pheno_fileV$subject_id),samplesV)
genomatrixV = genomatrixV[IV,]

rownames(genomatrixV) <- pheno_fileV$subject_id



aggregate_byV = 'Gene18b'





# SeqMeta Glm model

Fit_model_cgvhdn <- prepScores(
  Z = genomatrixV,
  formula = cgvhd ~ age +donorage+donorgp+atgcampathgp_r +prophx +pc1+ pc2,
  family = binomial(),
  SNPInfo = snpInfoV,
  snpNames = "SNP",
  aggregateBy = aggregate_byV,
  data = pheno_fileV,
  verbose = TRUE
)

Fit_skato_cgvhn  <- skatOMeta(
  Fit_model_cgvhdn,
  SNPInfo = snpInfoV,
  aggregateBy = aggregate_byV,
  snpNames = "SNP",
  verbose = TRUE
)



Ref.

[seqMeta](http://cran.nexr.com/web/packages/seqMeta/seqMeta.pdf)

[seqMeta: an R Package for meta-analyzing region-based tests of rare DNA variants](https://rdrr.io/cran/seqMeta/f/inst/doc/seqMeta.pdf)

[Meta-MultiSKAT: Multiple phenotype meta-analysis for region-based association test](https://www.biorxiv.org/content/10.1101/593814v1.full)

