# code used to extract 16S rRNA sequences from 
library(phyloseq)
library(tidyverse)

setwd("C:/Users/murra/OneDrive/Desktop/Research/Experiments/Act_N2O/Microbial data/Sequence_data_raw/250324_M07914_0225_000000000-M258W")
seqtab.nochim=readRDS("C:/Users/murra/OneDrive/Desktop/Research/Experiments/Act_N2O/Microbial data/Sequence_data_raw/250324_M07914_0225_000000000-M258W/filtered/sequencetable_nochimeras.rds")
ps=readRDS("C:/Users/murra/OneDrive/Desktop/Research/Experiments/Act_N2O/Microbial data/Sequence_data_raw/250324_M07914_0225_000000000-M258W/ps.rds")
ps.pruned=readRDS("C:/Users/murra/OneDrive/Desktop/Research/Experiments/Act_N2O/Microbial data/Sequence_data_raw/250324_M07914_0225_000000000-M258W/ps_pruned.rds")
##ordered(sample_sums(ps.pruned))
BLS_8=read.csv("BLS_8_output.csv")
melt.active=read.csv("melt.active.csv")

ps.1 <- ps

dna <- Biostrings::DNAStringSet(taxa_names(ps.1))
names(dna) <- taxa_names(ps.1)
ps.1 <- merge_phyloseq(ps.1, dna)
taxa_names(ps.1) <- paste0("ASV", seq(ntaxa(ps.1)))
ps.1

asv_sequences <- phyloseq::refseq(ps.1)

BLS_ASV1="ASV109"
BLS_ASV4=c("ASV1317", "ASV1742", "ASV2532", "ASV2544")
BLS_ASV10=c("ASV1278", "ASV1305", "ASV1550", "ASV1688", "ASV1799", "ASV1883", "ASV2209", "ASV2544", "ASV2556", "ASV2767")
BLS_ASV8=BLS_8$members

asv_id1=which(asv_sequences@ranges@NAMES %in% BLS_ASV1)
asv_id4=which(asv_sequences@ranges@NAMES %in% BLS_ASV4)
asv_id10=which(asv_sequences@ranges@NAMES %in% BLS_ASV10)
asv_id8=which(asv_sequences@ranges@NAMES %in% BLS_ASV8)

BLS_seq1=dna@ranges@NAMES[asv_id1]
BLS_seq4=dna@ranges@NAMES[asv_id4]
BLS_seq10=dna@ranges@NAMES[asv_id10]
BLS_seq8=dna@ranges@NAMES[asv_id8]


meqo.taxonomy1=melt.active[melt.active$OTU%in%BLS_ASV1,]%>% arrange(Sample)
meqo.taxonomy4=melt.active[melt.active$OTU%in%BLS_ASV4,]%>% arrange(Sample)
meqo.taxonomy10=melt.active[melt.active$OTU%in%BLS_ASV10,]
meqo.taxonomy8=melt.active[melt.active$OTU%in%BLS_ASV8,]%>% arrange(Sample)


write.csv(meqo.taxonomy1, "meqo_taxonomy_1.csv")
write.csv(meqo.taxonomy8, "meqo_taxonomy_8.csv")


firmicutes=melt.active[melt.active$Phylum=="Firmicutes",]
firmicutes.asv=unique(firmicutes$OTU)
firmicutes.id=which(asv_sequences@ranges@NAMES %in% firmicutes.asv)

firmicutes.seq=data.frame(dna@ranges@NAMES[firmicutes.id])
firmicutes.asv=na.omit(firmicutes.asv)
firmicutes.taxonomy=melt.active[melt.active$OTU%in%firmicutes.asv,]%>% arrange(Sample)
firmicutes.taxonomy=firmicutes.taxonomy[!duplicated(firmicutes.taxonomy$OTU),]%>% arrange(Abundance)
write.csv(firmicutes.taxonomy, "firmicutes_taxonomy.csv")

topasv=c("ASV1","ASV2")
firmicutes.id=which(asv_sequences@ranges@NAMES %in% topasv)
firmicutes.seq=dna@ranges@NAMES[firmicutes.id]
