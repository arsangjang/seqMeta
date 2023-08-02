

# Performing a survival meta-analysis of rare variants from sequencing data using seqMeta/MetaSKAT

## Instalation SKAT/ MetaSKAT/seqMeta


setwd("D:/folder")
install_github("brentp/skatMeta")

library(CompQuadForm)
library(twostageGWASsurvival)
library(seqMeta)
library(MASS)
library(writexl)
library(Matrix)



#### Z genomatrix
genomatrix = read.table("allel_binary_data.txt",header=TRUE)
snpInfo = read.csv('SNPInfo.csv',header=TRUE)
names(genomatrix)

I = match(snpInfo$POS, genomatrix$POS)
genomatrix = genomatrix[I,]
View(genomatrix)
genomatrix = t(genomatrix[,-c(1:3)])
colnames(genomatrix) = snpInfo$SNP


### Running survival model


cox_model<-prepCox(genomatrix, Surv(time,binary_outcome) ~ Age +  Sex +  intervention, 
             SNPInfo = SNPInfo, data = my_data, verbose = FALSE)

outacox_model <- skatMeta(cox_model, SNPInfo = SNPInfo)
head(outacox_model)
print(outacox_model)


######################################################################
## SKAT Optimal Test
#SKAT O Test introduced to test weighted averages of SKAT & burden tests
######################################################################

model1<- prepScores(genomatrix_a, outcome ~ sex+age+ intervention, 
             SNPInfo = SNPInfo, data = my_data, verbose = FALSE)

model2<-prepScores(genomatrix_b, outcome ~ sex+age+ intervention, 
             SNPInfo = SNPInfo, data = my_data_b, verbose = FALSE)


Out_combined <- skatOMeta(model1,model2, SNPInfo=SNPInfo, method="int")
head(Out_combined)


##########################################################
##Survival data
##########################################################

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
print(overall_Cox)
