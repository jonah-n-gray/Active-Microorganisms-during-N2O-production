### this file contains the code used to run the mEQO analysis.  




library(mEQO)
library(slam)
library(gurobi)
library(tidyverse)
library(ggplot2)
library(dplyr)
palette_14<- c( "#1F78B4","#E31A1C",  "#FB9A99", "#33A02C", "#FF7F00",  "#FDBF6F", "#CAB2D6",
                "gold",  "#EB05DB","#FFFF99","#EDCEEB","#1A635A","#9AD6CE" , "#969696")

#reading in data
setwd("C:/Users/murra/OneDrive/Desktop/Research/Experiments/Act_N2O/Microbial data/Sequence_data_raw/250324_M07914_0225_000000000-M258W")
melt.active=read.csv("melt.active.csv")
meta.active=read.csv("meta.active.csv")
absolute.abundance=read.csv("absolute.abundance.csv")
genus.abundance=read.csv("genus.abundance.csv")
asv.abundance=read.csv("ASV.abundance.csv")

#reorganizing dataframe to be complient with mEQO inputs
meqo_microbiome=asv.abundance[asv.abundance$Abundance !=0,]
meqo_microbiome <- reshape(
  meqo_microbiome[, c("Sample", "OTU", "abundance_per_gdw")],
  timevar = "OTU",
  idvar = "Sample",
  direction = "wide"
)

#reordering samples.  For some reason, they keep being read in out of order.
samplelist=meqo_microbiome$Sample
samplelist=c(samplelist[21:25],samplelist[1:20])
meqo_microbiome=meqo_microbiome[order(meqo_microbiome$Sample), ]
meqo_microbiome=meqo_microbiome[ ,-1 ]
#set NA to 0
meqo_microbiome[is.na(meqo_microbiome)] = 0
meqo_microbiome=meqo_microbiome[-11, ]

#reading in N2O data
N2O_trait=meta.active$ug_N2O_N_per_gdw_hr
N2O_trait=N2O_trait[-11]

#filtering out rare asvs (only show up in 10/24 samples)
meqo_microbiome_counts=colSums(meqo_microbiome>0)
meqo_microbiome_filtered=meqo_microbiome[,meqo_microbiome_counts>9]

saveRDS(meqo_microbiome_filtered, "meqo_microbiome_filtered.rds")

#########
#Finding the best model
BLSaic=sapply(1:20,function(N){
  assemblage=(EQO_bls(meqo_microbiome_filtered,N2O_trait, Nmax=N, K=10))
  return(list(AIC=AIC(lm(N2O_trait~assemblage$abundance)), members=length(assemblage$members), y=assemblage$y))
})
aic.table=data.frame(BLSaic)

#running the model
#8 ASV model is best.  comparing to 1 asv model

BLS_8=EQO_bls(meqo_microbiome_filtered,N2O_trait, Nmax=8, K=10)
BLS_1=EQO_bls(meqo_microbiome_filtered,N2O_trait, Nmax=1, K=10)


#write.csv(BLS_8, "BLS_8_output.csv")