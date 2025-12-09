###Jonah Gray
##This file contains code used to make the figures in the manuscript. It also contains some additional data proceessing steps not included in the 16S rRNA markdown file.
##These steps include calculating relative and absolute abundances, and organizing data.
##Not all plots included in this file are in the manuscript.  
##exported figures were made publication ready in Adobe Illustrator.


library(ggplot2)
library(dada2)
library(ShortRead)
library(Biostrings)
library(phyloseq)
library(corncob)
library(DESeq2)
library(microbiome)
library(DECIPHER)
library(phangorn)
library(tibble)
library(lme4)
library(lmerTest)
library(ggplot2)
library(car)
library(vegan)
library(RVAideMemoire)
library(emmeans)
library(tidyverse)

palette_50 <- c(
  "#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231",
  "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe",
  "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000",
  "#aaffc3", "#808000", "#ffd8b1", "#000080", "#808080",
  "#FFFFFF", "#000000", "#e41a1c", "#377eb8", "#4daf4a",
  "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf",
  "#999999", "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
  "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3", "#1b9e77",
  "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02",
  "#a6761d", "#666666", "#7fc97f", "#beaed4", "#fdc086"
)
palette_14<- c( "purple","#E31A1C",  "#FB9A99", "#33A02C", "#FF7F00", "#9AD6CE" , "#CAB2D6",
                "gold",  "#8B4513","blue","#f032e6","#1A635A", "#FDBF6F", "#969696")

palette_6=c("black", "#228B22", "#696969", "#8B4513", "red", "#377eb8")

setwd("C:/Users/murra/OneDrive/Desktop/Research/Experiments/Act_N2O/Microbial data/Sequence_data_raw/250324_M07914_0225_000000000-M258W")
seqtab.nochim=readRDS("C:/Users/murra/OneDrive/Desktop/Research/Experiments/Act_N2O/Microbial data/Sequence_data_raw/250324_M07914_0225_000000000-M258W/filtered/sequencetable_nochimeras.rds")
ps=readRDS("C:/Users/murra/OneDrive/Desktop/Research/Experiments/Act_N2O/Microbial data/Sequence_data_raw/250324_M07914_0225_000000000-M258W/ps.rds")
ps.pruned=readRDS("C:/Users/murra/OneDrive/Desktop/Research/Experiments/Act_N2O/Microbial data/Sequence_data_raw/250324_M07914_0225_000000000-M258W/ps_pruned.rds")
ordered(sample_sums(ps.pruned))

ps.1 <- ps
dna <- Biostrings::DNAStringSet(taxa_names(ps.1))
names(dna) <- taxa_names(ps.1)
ps.1 <- merge_phyloseq(ps.1, dna)
taxa_names(ps.1) <- paste0("ASV", seq(ntaxa(ps.1)))
ps.1


alpha.div<-estimate_richness(ps.pruned, measures=c("Shannon", "Observed"))
set.seed(500)
ps.pruned=rarefy_even_depth(ps.pruned)

ps.perc <- transform_sample_counts(ps.pruned, function(x) x / sum(x)) 

ps.active= subset_samples(ps.perc, sample_type=="a")
ps.inactive =subset_samples(ps.perc, sample_type=="i")
ps.all =subset_samples(ps.perc, sample_type=="a+i")
ps.bulk =subset_samples(ps.perc, sample_type=="bulk")

meta <- read.csv("metadata.csv", header = TRUE, row.names = 1)
alpha.div2<-estimate_richness(ps.pruned, measures=c("Shannon", "Observed"))
even <- evenness(ps.pruned, 'pielou')
meta=meta[-c(43,77),]
meta$Shannon <- paste(alpha.div$Shannon)
meta$Observed <- paste(alpha.div$Observed)
meta$Evenness <- even$pielou
meta$Observed <- as.numeric(meta$Observed)
meta$Shannon <- as.numeric(meta$Shannon)
meta$originalID<- as.numeric(meta$originalID)
meta$time_point<- as.numeric(meta$time_point)
meta$cell_count<- as.numeric(meta$cell_count)
meta$sample_type <- as.factor(meta$sample_type)

meta %>% rownames_to_column(var = "Sample") %>% filter(sample_type=="bulk") -> meta_bulk
meta %>% rownames_to_column(var = "Sample") %>% filter(sample_type=="a") ->meta_active
meta %>% rownames_to_column(var = "Sample") %>% filter(sample_type=="i") ->meta_inactive
meta %>% rownames_to_column(var = "Sample") %>% filter(sample_type=="a+i") ->meta_all
ps.phyla.perc.active = tax_glom(ps.active, "Phylum")
ps.phyla.perc.inactive = tax_glom(ps.inactive, "Phylum")
melt.active=psmelt(ps.active)
melt.inactive=psmelt(ps.inactive)

melt.phylum.active <- psmelt(ps.phyla.perc.active)
melt.phylum.inactive <- psmelt(ps.phyla.perc.inactive)
### setting relative abundance to 1
phyl.means.active <- aggregate(Abundance~Sample+Phylum, melt.phylum.active, FUN=mean)
phyl.means.inactive <- aggregate(Abundance~Sample+Phylum, melt.phylum.inactive, FUN=mean)
#finds relative abundance for each sample so we can calculate absolute abundances
ASV.means.active <- aggregate(Abundance~Sample+OTU, melt.active, FUN=mean)
ASV.means.active$Sample=as.numeric(ASV.means.active$Sample)
meta_active$Sample=as.numeric(meta_active$Sample)
ASV.abs.abun.active= ASV.means.active %>% left_join(meta_active, by = "Sample") %>% mutate(absolute_abundance=cell_count*Abundance) %>% mutate(abundance_per_gdw=Abundance*count_per_gdw)

ASV.abs.abun.active %>%
  filter(time_point == 5,
         OTU %in% c("ASV1", "ASV2")) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE),sd= sd(Abundance))

ASV.abs.abun.active %>%
  filter(time_point == 5) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE), sd= sd(Abundance))

ASV.abs.abun.active %>%
  filter(Sample == "101") %>%
  summarise(total_abundance = sum(Abundance))

phyl.means.active$Sample <- as.character(phyl.means.active$Sample)
meta_active$Sample <- as.character(meta_active$Sample)
phyl.means.inactive$Sample <- as.character(phyl.means.inactive$Sample)
meta_inactive$Sample <- as.character(meta_inactive$Sample)

phyl.abs.abun.active= phyl.means.active %>% left_join(meta_active, by = "Sample") %>% mutate(absolute_abundance=cell_count*Abundance) %>% mutate(abundance_per_gdw=Abundance*count_per_gdw)
phyl.abs.abun.inactive= phyl.means.inactive %>% left_join(meta_inactive, by = "Sample") %>% mutate(absolute_abundance=cell_count*Abundance) %>% mutate(abundance_per_gdw=Abundance*count_per_gdw)
phyl.abs.abun.active$Phylum[phyl.abs.abun.active$Abundance<.01] <- "xother"
phylum=unique(phyl.abs.abun.active$Phylum)
phyl.abs.abun.inactive$Phylum[!(phyl.abs.abun.inactive$Phylum %in% phylum)] <- "xother"

phyl.abs.abun.active %>%
  filter(time_point == 4,
         Phylum %in% "Firmicutes") %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE),sd= sd(Abundance))

phyl.abs.abun.active %>%
  filter(time_point == 5,
         Phylum %in% "Firmicutes") %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE),sd= sd(Abundance))

alpha.div2<-estimate_richness(ps.pruned, measures=c("Chao1"))
meta$Chao1 <- paste(alpha.div2$Chao1)
meta %>% rownames_to_column(var = "Sample") %>% filter(sample_type=="a") ->meta_active

####data analysis

###observed
one_way=aov(Observed~as.factor(sample_type), meta)
tukey=TukeyHSD(one_way)
summary(one_way)
tukey=as.data.frame(tukey$`as.factor(sample_type)`)
write.csv(tukey, "tukey_observed_sample_output.csv")
anova=as.data.frame(summary(one_way)[[1]])
write.csv(anova, "anova_observed_sample_output.csv")

one_way=aov(Observed~as.factor(time_point), meta_active)
tukey=TukeyHSD(one_way)
summary(one_way)
tukey=as.data.frame(tukey$`as.factor(time_point)`)
write.csv(tukey, "tukey_observed_active_output.csv")
anova=as.data.frame(summary(one_way)[[1]])
write.csv(anova, "anova_observed_active_output.csv")

one_way=aov(Observed~as.factor(time_point), meta_inactive)
tukey=TukeyHSD(one_way)
summary(one_way)
tukey=as.data.frame(tukey$`as.factor(time_point)`)
write.csv(tukey, "tukey_observed_inactive_output.csv")
anova=as.data.frame(summary(one_way)[[1]])
write.csv(anova, "anova_observed_inactive_output.csv")

one_way=aov(Observed~as.factor(time_point), meta_all)
tukey=TukeyHSD(one_way)
summary(one_way)
tukey=as.data.frame(tukey$`as.factor(time_point)`)
write.csv(tukey, "tukey_observed_actinact_output.csv")
anova=as.data.frame(summary(one_way)[[1]])
write.csv(anova, "anova_observed_actinact_output.csv")

### evenness
one_way=aov(Evenness~as.factor(sample_type), meta)
tukey=TukeyHSD(one_way)
summary(one_way)
tukey=as.data.frame(tukey$`as.factor(sample_type)`)
write.csv(tukey, "tukey_Evenness_sample_output.csv")
anova=as.data.frame(summary(one_way)[[1]])
write.csv(anova, "anova_evenness_sample_output.csv")

one_way=aov(Evenness~as.factor(time_point), meta_active)
tukey=TukeyHSD(one_way)
summary(one_way)
tukey=as.data.frame(tukey$`as.factor(time_point)`)
write.csv(tukey, "tukey_Evenness_active_output.csv")
anova=as.data.frame(summary(one_way)[[1]])
write.csv(anova, "anova_evenness_active_output.csv")

one_way=aov(Evenness~as.factor(time_point), meta_inactive)
tukey=TukeyHSD(one_way)
summary(one_way)
tukey=as.data.frame(tukey$`as.factor(time_point)`)
write.csv(tukey, "tukey_Evenness_inactive_output.csv")
anova=as.data.frame(summary(one_way)[[1]])
write.csv(anova, "anova_evenness_inactive_output.csv")

one_way=aov(Evenness~as.factor(time_point), meta_all)
tukey=TukeyHSD(one_way)
summary(one_way)
tukey=as.data.frame(tukey$`as.factor(time_point)`)
write.csv(tukey, "tukey_Evenness_actinact_output.csv")
anova=as.data.frame(summary(one_way)[[1]])
write.csv(anova, "anova_evenness_actinact_output.csv")


###permanova on bray-curtis
bc.asv.active <- phyloseq::distance(otu_table(ps.active), "bray")
bc.asv.inactive <- phyloseq::distance(otu_table(ps.inactive), "bray")
bc.asv.all= phyloseq::distance(otu_table(ps.perc), "bray")

#permanova for time and sample type
permanova_bray_all=adonis2(bc.asv.all~sample_type * time_point, data= meta, permutations = 999, method = "bray")
summary(permanova_bray_all)
permanova_bray_all
permanova=permanova_bray_all
write.csv(permanova, "permanova_permanova_bray_all_output.csv")


## 21% of the variation is explained by whether the sample was active or inactive

perm_active=adonis2(bc.asv.active~time_point, data= meta_active, permutations = 999, method = "bray")
permanova=perm_active
write.csv(permanova, "permanova_perm_active_output.csv")

perm_inactive=adonis2(bc.asv.inactive~ time_point, data= meta_inactive, permutations = 999, method = "bray")
permanova=perm_inactive
write.csv(permanova, "permanova_perm_inactive_output.csv")

perm_active #25% of variation is explained by time
perm_inactive # only 5% of the variation is explained by time

###Microbiome data plots

plot(Observed ~ time_point, data = meta_active)
plot(Chao1 ~ time_point, data = meta_active)
plot(Shannon ~ time_point, data = meta_active)
plot(Evenness ~ time_point, data = meta_active)

ggplot(meta_active, aes(x = factor(time_point), y = as.numeric(Observed),group = factor(time_point))) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray") +   # remove outliers so they don't duplicate
  geom_jitter(width = 0.2, size = 3, alpha = 0.7, colour = "darkblue") +  # scatter points
  theme_classic()+xlab("Time group") + ylab("Observed ASV")
ggplot(meta_active, aes(x = factor(time_point), y = as.numeric(Evenness),group = factor(time_point))) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray") +   # remove outliers so they don't duplicate
  geom_jitter(width = 0.2, size = 3, alpha = 0.7, colour = "darkblue") +  # scatter points
  theme_classic()+xlab("Time group") + ylab("Evenness")


observed= ggplot(meta, aes(x = factor(time_point), y = as.numeric(Observed),group = factor(time_point))) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray") +   # remove outliers so they don't duplicate
  geom_jitter(width = 0.2, size = 3, alpha = 0.7, colour = "darkblue") +  # scatter points
  facet_grid(~ sample_type, scales = "free_x", space = "free_x")+
  theme_classic()+xlab("Time group") + ylab("Observed ASV")
evenness= ggplot(meta, aes(x = factor(time_point), y = as.numeric(Evenness),group = factor(time_point))) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray") +   # remove outliers so they don't duplicate
  geom_jitter(width = 0.2, size = 3, alpha = 0.7, colour = "darkblue") +  # scatter points
  facet_grid(~ sample_type, scales = "free_x", space = "free_x")+
  theme_classic()+xlab("Time group") + ylab("Evenness")
ggsave("observed.pdf", observed, width = 5000, height = 2000, units = "px")
ggsave("evenness.pdf", evenness, width = 5000, height = 2000, units = "px")

ggplot(meta_active, aes(x = factor(time_point), y = as.numeric(Chao1),group = factor(time_point))) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray") +   # remove outliers so they don't duplicate
  geom_jitter(width = 0.2, size = 3, alpha = 0.7, colour = "darkblue") +  # scatter points
  scale_y_continuous(breaks = seq(600, 1400, by = 100),limits = c(600, 1400)) +
  theme_classic()

ggplot(meta_active, aes(x = factor(time_point), y = as.numeric(Shannon),group = factor(time_point))) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray") +   # remove outliers so they don't duplicate
  geom_jitter(width = 0.2, size = 3, alpha = 0.7, colour = "darkblue") +  # scatter points
  #cale_y_continuous(breaks = seq(600, 1400, by = 100),limits = c(600, 1400)) +
  theme_classic()

sample_data(ps.perc)$time_point <- as.factor(sample_data(ps.perc)$time_point)
PCOA_BRAY=plot_ordination(ps.perc, ordinate(ps.perc, "PCoA", "bray"),
                          color = "sample_type", shape = "time_point") +
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(size = 10)),
         shape = guide_legend(override.aes = list(size = 10))) +
  geom_point(size = 7)+
  scale_shape_manual(values = c(1,2,3,4,7,8))+theme(axis.text = element_text(size = 20))+
  theme(axis.title.x = element_text(size = 30), 
        axis.title.y = element_text(size = 30))
PCOA_BRAY_sample=PCOA_BRAY+ stat_ellipse(geom = "polygon", type="norm", linetype= 2, alpha=0.0, aes(group= sample_type, fill=sample_type))
ggsave("PCOA_BRAY_sample.pdf", PCOA_BRAY_sample, width = 4000, height = 3000, units = "px")

PCOA_BRAY_sample_time=PCOA_BRAY+ stat_ellipse(geom = "polygon", type="norm", linetype= 2, alpha=0.0, aes(fill=sample_type))
ggsave("PCOA_BRAY_sample_time.pdf", PCOA_BRAY_sample_time, width = 4000, height = 3000, units = "px")

###N2O and activity plots and models
N2O_v_activity_lm= lm(meta_active$ug_N2O_N_per_gdw_hr~meta_active$percent_gated)
summary(N2O_v_activity_lm)

N2O_v_fluor_lm= lm(meta_active$ug_N2O_N_per_gdw_hr~log(meta_active$af_median))
summary(N2O_v_fluor_lm)

activity_v_time_lm= lm(meta_active$percent_gated~meta_active$t)
summary(activity_v_time_lm)

fluor_v_time_lm= lm(log(meta_active$af_median)~meta_active$t)
summary(fluor_v_time_lm)


N2O_v_activity= ggplot(meta_active, aes(x = meta_active$percent_gated, y = meta_active$ug_N2O_N_per_gdw_hr)) +
  geom_smooth(method = "lm", colour = "blue", fill = "gray", size=3) +
  geom_point(colour = "darkblue", size=4) +  # Add points
  # Set x and y limits
  scale_y_continuous(breaks = seq(0, 5.5, by = 1),limits = c(0, 5.5)) +
  theme_classic()+xlab("Percent activity") + ylab("ug N2O-N gdw-1 hr-1") +  # Set axis labels
  ggtitle("N2O-N production rate vs Percent activity (active cells/total cells)") +  # Set plot title
  theme(title = element_text(size = 20),
        axis.text = element_text(size = 20))+
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20))  # Use a minimal theme, you can change this to other themes if you prefer

N2O_v_fluor= ggplot(meta_active, aes(x = log(meta_active$af_median), y = meta_active$ug_N2O_N_per_gdw_hr)) +
  geom_smooth(method = "lm", colour = "blue", fill = "gray", size=1) +
  geom_point(colour = "darkblue", size= 3) +  # Add points
  # Set x and y limits
  scale_y_continuous(breaks = seq(0, 5, by = 1),limits = c(0, 5)) +
  theme_classic()+xlab("Log(median fluorescence)") + ylab("ug N2O-N gdw-1 hr-1") +  # Set axis labels
  ggtitle("Level of activity vs N2O-N production rate") +  # Set plot title
  theme(title = element_text(size = 20),
        axis.text = element_text(size = 20))+
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20))# Use a minimal theme, you can change this to other themes if you prefer

ggsave("N2O_v_fluor.pdf", N2O_v_fluor, width = 3000, height = 3000, units = "px")

activity_v_fluor= ggplot(meta_active, aes(x = log(meta_active$af_median), y = meta_active$percent_gated)) +
  geom_smooth(method = "lm", colour = "blue", fill = "gray", size=3) +
  geom_point(colour = "darkblue", size= 4) +  # Add points
  # Set x and y limits
  theme_classic()+xlab("Log(median fluorescence)") + ylab("Percent activity") +  # Set axis labels
  ggtitle("Level of activity vs Percent activity") +  # Set plot title
  theme(title = element_text(size = 20),
        axis.text = element_text(size = 20))+
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20))# Use a minimal theme, you can change this to other themes if you prefer

activity_v_time= ggplot(meta_active, aes(x = meta_active$t, y = meta_active$percent_gated, group=t)) +
  geom_boxplot(colour = "darkblue",fill = "gray")+
  geom_jitter(width = 0.2, color="blue", size=3)+
  scale_x_continuous(breaks = seq(0, 15, by =3),limits = c(0, 18)) +
  theme_classic()+xlab("Time (hr)") + ylab("Pecent of organisms active") +  # Set axis labels
   # Set plot title
  theme(title = element_text(size = 20),
        axis.text = element_text(size = 20))+
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20))# Use a minimal theme, you can change this to other themes if you prefer

ggsave("activity_v_time.pdf", activity_v_time, width = 3000, height = 2000, units = "px")
activity_v_time


  
  ggplot(meta_active, aes(x = meta_active$Sample, y = meta_active$percent_gated, group=t)) +
    geom_boxplot(colour = "darkblue",fill = "gray")+
    geom_jitter(width = 0.2, color="blue", size=3)
  
  
### mirror plot moment
act_inact=rbind(phyl.abs.abun.active,phyl.abs.abun.inactive)
act_inact$Sample=as.numeric(act_inact$Sample)
act_inact$Sample=c((act_inact$Sample[1:1100]+100),act_inact$Sample[1101:2112])
act_inact$abundance_per_gdw=c(act_inact$abundance_per_gdw[1:1100],(act_inact$abundance_per_gdw[1101:2112]*-1))


ggplot(act_inact, aes(x = Sample, y = abundance_per_gdw, fill = Phylum)) + theme_classic() +
  geom_bar(stat = "identity",position = "stack" ) + scale_fill_manual(values=palette_14) + 
  facet_grid(~ time_point, scales = "free_x", space = "free_x") +  theme(axis.text.x = element_text(angle = 45, hjust = 1))

mirror_lower= ggplot(act_inact, aes(x = Sample, y = abundance_per_gdw, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = palette_14) +
  ylim(-1.7e+8,0)+
  facet_grid(~ time_point, scales = "free_x", space = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("mirror_lower.pdf", mirror_lower, width = 3000, height = 2000, units = "px")

mirror_upper=ggplot(act_inact, aes(x = Sample, y = abundance_per_gdw, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = palette_14) +
  ylim(0,1e+6)+
  facet_grid(~ time_point, scales = "free_x", space = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("mirror_upper.pdf", mirror_upper, width = 3000, height = 2000, units = "px")


####BLS and meqo plots
setwd("C:/Users/murra/OneDrive/Desktop/Research/Experiments/Act_N2O/Microbial data/Sequence_data_raw/250324_M07914_0225_000000000-M258W")
melt.active=read.csv("melt.active.csv")
meta.active=read.csv("meta.active.csv")
meqo_microbiome_filtered=readRDS("meqo_microbiome_filtered.rds")
BLS_8=read.csv("BLS_8_output.csv")
BLS_8$members=paste0("abundance_per_gdw.", BLS_8$members)


N2O_trait=meta.active$ug_N2O_N_per_gdw_hr
N2O_trait=N2O_trait[-11]

samplelist=seq(101,125,by=1)
samplelist=samplelist[-11]
meqoplot_members8=data.frame(samplelist,meqo_microbiome_filtered[,BLS_8$members],N2O_trait,time=meta.active$time_point[-11])
meqoplot_members8=meqoplot_members8[,c(0:9,26,27)]
meqoplot_members8=meqoplot_members8%>%
  pivot_longer(
    cols= -c(samplelist,N2O_trait,time),
    names_to = "ASV",
    values_to = "Abundance per gdw"
  ) 
##
coeff <- .005
blsplot_8=ggplot(meqoplot_members8, aes(x= samplelist))+ 
  geom_bar(aes(y=`Abundance per gdw` , fill=ASV), stat= "identity")+
  scale_fill_manual(values=palette_14)+
  geom_point(aes(y=N2O_trait/coeff),colour = "black", size =3)+ 
  scale_y_continuous(name= "viable cell abundance g soil-1", limits =c(0,800),sec.axis = sec_axis(~.*coeff, name= "ug N2O-N g soil-1 hr-1"))+ 
  facet_grid(~ time, scales = "free_x", space = "free_x") +  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Replicates")+
  theme_classic()+
  ggtitle("") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20))
blsplot_8
ggsave("blsplot_8.pdf", blsplot_8, width = 3000, height = 2000, units = "px")


###Nitrate and ammonia data
setwd("C:/Users/murra/OneDrive/Desktop/Research/Experiments/Act_N2O/Soil Data (water, N)")
t=c(3,5.67,8.42,12,15)
time=rep(t, each=5)
dat_no3_nh4=read.csv("Final NO3 NH4 data.csv")
dat_no3_nh4=dat_no3_nh4[,c(1,3,4)]
dat=as.data.frame(cbind(time,dat_no3_nh4))
dat=dat[-9,] # removing an outlier
data_all=data.frame(Time = dat$time, mg_NH4_N = dat$Final.conc..NH4.N..mg.kg.1.soil.,mg_NO3_N = dat$Final.conc..NO3.N..mg.kg.1.soil.)

no3_nh4=ggplot(data_all, aes(x = Time)) +
  geom_smooth(aes(y = mg_NO3_N), method = "lm", colour = "blue", fill = "gray", size=3) +
  geom_point(aes(y = mg_NO3_N), colour = "darkblue", size = 4) +  # Add points
  geom_smooth(aes(y = mg_NH4_N), method = "lm", colour = "maroon", fill = "gray", size = 3) +
  geom_point(aes(y = mg_NH4_N), colour = "darkred", size =4) +  # Add points
  scale_x_continuous(breaks = seq(0, 15, by = 3),limits = c(3, 15)) +
  xlab("Time (hr)") + ylab("mg N/kg soil") +  # Set axis labels
  ggtitle("Soil NH4-N and NO3-N") +  # Set plot title
  theme_classic()+
  theme(
    axis.text = element_text(size = 20),
    axis.title.x = element_text(size = 25), 
    axis.title.y = element_text(size = 25))  # Use a minimal theme, you can change this to other themes if you prefer


###Nitrous oxide data
setwd("C:/Users/murra/OneDrive/Desktop/Research/Experiments/Act_N2O/Licor data")
N2O_flux_midpoint=read.csv("N2O_flux_midpoint.csv")
dat=readRDS("n2o_dat.rds")
data2=data.frame(Time=dat$time, N2O_N=dat$ug_N2O_N_per_g_dry_hr)

N2O_flux=ggplot(N2O_flux_midpoint, aes(x = midpoint.hr, y = ug_N2O_N_per_g_dryhr)) +
  geom_smooth(method = "lm", colour = "blue", fill = "gray") +
  geom_point(colour = "darkblue", size =3) +  # Add points
  xlim(0, 16) + ylim(0, 4.0) +  # Set x and y limits
  xlab("Time (hr)") + ylab("ug N2O-N gdw-1 hr-1") +  # Set axis labels
  ggtitle("N2O-N production rate") +  # Set plot title
  theme_classic()+ theme(title = element_text(size = 25),
                         axis.text = element_text(size = 20))+
  theme(axis.title.x = element_text(size = 25), 
        axis.title.y = element_text(size = 25))



coeff=100
N2O_flux_nh4_no3<- ggplot() +
  # Second dataset (24 rows): NO3
  geom_point(data = data_all,
             aes(x = Time, y = mg_NO3_N),
             colour = "darkgreen", size = 3) +
  geom_smooth(data = data_all,
              aes(x =Time, y = mg_NO3_N),
              method = "lm", colour = "green", fill = "gray", lwd = 1) +
  # Second dataset (24 rows): NH4
  geom_point(data = data_all,
             aes(x = Time, y = mg_NH4_N),
             colour = "darkred", size = 3) +
  geom_smooth(data = data_all,
              aes(x = Time, y = mg_NH4_N),
              method = "lm", colour = "maroon", fill = "gray", lwd = 1) +
  # First dataset (49 rows)
  geom_point(data = N2O_flux_midpoint,
             aes(x = midpoint.hr, y = ug_N2O_N_per_g_dryhr*coeff/2.5),
             colour = "darkblue", size = 3) +
  geom_smooth(data = N2O_flux_midpoint,
              aes(x = midpoint.hr, y = ug_N2O_N_per_g_dryhr*coeff/2.5),
              method = "lm", colour = "blue", fill = "gray",lwd = 1) +
  scale_y_continuous(name = "mg N kg soil-1", limits =c(0,200),sec.axis = sec_axis(~ . / coeff*2.5, name = "ug N2O-N g soil-1 hr-1")) +
  scale_x_continuous(breaks = seq(0, 15, by = 3),limits = c(0, 15)) +
  xlab("Time (hr)") +
  ggtitle("") +
  theme_classic() +
  theme(title = element_text(size = 20),
        axis.text = element_text(size = 20))+
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20))
N2O_flux_nh4_no3

ggsave("N2O_flux_nh4_no3.pdf", N2O_flux_nh4_no3, width = 3000, height = 3000, units = "px")



