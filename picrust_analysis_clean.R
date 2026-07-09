#this file is for picrust data analysis using the ggpicrust2 package
#J. Gray, 4/29/26
#remotes::install_version("ggplot2", version = "3.4.0")
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(GGally)
library(KEGGREST)
library(ALDEx2)
library(data.table)
library(phyloseq)

#read in phyloseq objects for finding mEQO ASVs.  Duplicate code from 16S data analysis.
ps=readRDS("C:/Users/murra/OneDrive/Desktop/Research/Experiments/Act_N2O/Microbial data/Sequence_data_raw/250324_M07914_0225_000000000-M258W/ps.rds")
ps.pruned=readRDS("C:/Users/murra/OneDrive/Desktop/Research/Experiments/Act_N2O/Microbial data/Sequence_data_raw/250324_M07914_0225_000000000-M258W/ps_pruned.rds")
ordered(sample_sums(ps.pruned))

ps.1 <- ps
dna <- Biostrings::DNAStringSet(taxa_names(ps.1))
names(dna) <- taxa_names(ps.1)
ps.1 <- merge_phyloseq(ps.1, dna)
taxa_names(ps.1) <- paste0("ASV_", seq(ntaxa(ps.1)))
ps.1

set.seed(500)
ps.pruned=rarefy_even_depth(ps.pruned)

ps.perc <- transform_sample_counts(ps.pruned, function(x) x / sum(x)) 

ps.active= subset_samples(ps.perc, sample_type=="a")
ps.inactive =subset_samples(ps.perc, sample_type=="i")
ps.all =subset_samples(ps.perc, sample_type=="a+i")
ps.bulk =subset_samples(ps.perc, sample_type=="bulk")
ps.phyla.perc.active = tax_glom(ps.active, "Phylum")
ps.phyla.perc.inactive = tax_glom(ps.inactive, "Phylum")
melt.active=psmelt(ps.active)
melt.inactive=psmelt(ps.inactive)

melt.phylum.active <- psmelt(ps.phyla.perc.active)
melt.phylum.inactive <- psmelt(ps.phyla.perc.inactive)

phyl.means.active <- aggregate(Abundance~Sample+Phylum, melt.phylum.active, FUN=mean)
phyl.means.inactive <- aggregate(Abundance~Sample+Phylum, melt.phylum.inactive, FUN=mean)
#finds relative abundance for each sample so we can calculate absolute abundances
phyl.means.active$Sample <- as.character(phyl.means.active$Sample)
#meta_active$Sample <- as.character(meta_active$Sample)
phyl.means.inactive$Sample <- as.character(phyl.means.inactive$Sample)
#meta_inactive$Sample <- as.character(meta_inactive$Sample)

phyl.abs.abun.active= phyl.means.active #%>% left_join(meta_active, by = "Sample") %>% mutate(absolute_abundance=cell_count*Abundance) %>% mutate(abundance_per_gdw=Abundance*count_per_gdw)
phyl.abs.abun.inactive= phyl.means.inactive #%>% left_join(meta_inactive, by = "Sample") %>% mutate(absolute_abundance=cell_count*Abundance) %>% mutate(abundance_per_gdw=Abundance*count_per_gdw)
phyl.abs.abun.active$Phylum[phyl.abs.abun.active$Abundance<.01] <- "xother"
phylum=unique(phyl.abs.abun.active$Phylum)
phyl.abs.abun.inactive$Phylum[!(phyl.abs.abun.inactive$Phylum %in% phylum)] <- "xother"


#Start picrust analysis
setwd("C:/Users/murra/OneDrive/Desktop/Research/Experiments/Act_N2O/Microbial data/Sequence_data_raw/250324_M07914_0225_000000000-M258W")
meta <- read.csv("metadata.csv", header = TRUE, row.names = 1)
meta$time_point<- as.numeric(meta$time_point)
meta$cell_count<- as.numeric(meta$cell_count)
meta$sample_type <- as.factor(meta$sample_type)
setwd("C:/Users/murra/OneDrive/Desktop/Research/Experiments/Act_N2O/Picrust")


#combined_KO_predicted=fread("combined_KO_predicted.tsv")
#matching_cols <- colnames(combined_KO_predicted) %in% paste0("ko:", denit_KO$KO, sep="")
#denit_ASV=combined_KO_predicted[ , ..matching_cols]
#denit_ASV=data_frame(combined_KO_predicted[,1], denit_ASV)
#denit_ASV=denit_ASV[rowSums(denit_ASV[, 2:12]) != 0, ]
#saveRDS(denit_ASV, "denit_ASV.rds")

denit_ASV=readRDS("denit_ASV.rds")

#asv_interest=read.csv("ASV_interest.csv", header = FALSE)
#ASV_KO=combined_KO_predicted= combined_KO_predicted[ combined_KO_predicted$sequence %in% asv_interest$V1,]
#saveRDS(ASV_KO, "ASV_KO.rds")

ASV_KO=readRDS("ASV_KO.rds")

abundance_file=read_tsv("pred_metagenome_unstrat.tsv")
abundance_file$`function` <- sub("^ko:", "", abundance_file$`function`)
colnames(abundance_file) [1] = "#NAME"
abundance_file <- abundance_file[, !(colnames(abundance_file) %in% c("217", "4"))]

picrust_metadata=read_delim(
  "picrust_metadata.txt",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)
picrust_metadata <- picrust_metadata[-c(43,77),]

####


### just active and inactive models
#relative abundance of denit genes across samples
denit_KO=read.csv("denit_KO.csv")

picrust_metadata_trimmed <- picrust_metadata[
  !(picrust_metadata$sample_type %in% c("a+i", "bulk")),]

abundance_file_trimmed <- abundance_file[, colnames(abundance_file) %in% picrust_metadata_trimmed$Sample]

picrust_active<- picrust_metadata[
  !(picrust_metadata$sample_type %in% c("a+i", "bulk","i")),]
abundance_active<- abundance_file[, colnames(abundance_file) %in% picrust_active$Sample]

####

#just active
col_totals <- colSums(abundance_active, na.rm = TRUE)

rel_abun_act <- sweep(abundance_active, 2,
                ifelse(col_totals == 0, NA, col_totals),
                FUN = "/")
rel_abun_act= cbind(abundance_file[,1],rel_abun_act)
rel_abun_act_denit = rel_abun_act[rel_abun_act$`#NAME` %in% denit_KO$KO, ]
rel_abun_act_denit <- rel_abun_act_denit[match(denit_KO$KO, rel_abun_act_denit$`#NAME`), ]
rel_abun_act_denit= data.frame(gene= denit_KO[,1], rel_abun_act_denit, check.names = FALSE)
rel_abun_act_denit=rel_abun_act_denit[-c(7,9),]
plot_df <- rel_abun_act_denit %>%
  pivot_longer(
    cols = 3:27,
    names_to = "sample",
    values_to = "rel_abundance"
  )

ggplot(plot_df, aes(x = gene, y = rel_abundance, fill = gene)) +
  geom_col() +
  facet_wrap(~ sample) +   # fixed y-axis across facets
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#active averages
act_denit_avg=rel_abun_act_denit[,-c(1,2)]

samples_by_time <- split(
  picrust_metadata_trimmed$Sample,
  picrust_metadata_trimmed$`time point`
)

for(tp in names(samples_by_time)) {
  
  sample_cols <- intersect(
    samples_by_time[[tp]],
    colnames(act_denit_avg)
  )
  
  act_denit_avg[[paste0("mean_", tp)]] <-
    rowMeans(act_denit_avg[, sample_cols, drop = FALSE], na.rm = TRUE)
  
  act_denit_avg[[paste0("sd_", tp)]] <-
    apply(act_denit_avg[, sample_cols, drop = FALSE],
          1,
          sd,
          na.rm = TRUE)
}

act_denit_avg=data.frame(gene=rel_abun_act_denit[,c(1)], act_denit_avg[,c(26:35)])

plot_df <- act_denit_avg |>
  pivot_longer(
    cols = -gene,
    names_to = c(".value", "timepoint"),
    names_sep = "_"
  )

# bar plot
ggplot(plot_df,
       aes(x = gene,
           y = mean,
           fill = gene)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean - sd,
                    ymax = mean + sd),
                width = 0.2) +
  facet_wrap(~timepoint) +
  labs(x = "Gene",
       y = "Mean Relative Abundance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))


#active and inactive
col_totals <- colSums(abundance_file_trimmed, na.rm = TRUE)

rel_abun_act_in <- sweep(abundance_file_trimmed, 2,
                      ifelse(col_totals == 0, NA, col_totals),
                      FUN = "/")
rel_abun_act_in= cbind(abundance_file[,1],rel_abun_act_in)
rel_abun_act_in_denit = rel_abun_act_in[rel_abun_act_in$`#NAME` %in% denit_KO$KO, ]
rel_abun_act_in_denit <- rel_abun_act_in_denit[match(denit_KO$KO, rel_abun_act_in_denit$`#NAME`), ]
rel_abun_act_in_denit= data.frame(gene= denit_KO[,1], rel_abun_act_in_denit, check.names = FALSE)
rel_abun_act_in_denit=rel_abun_act_in_denit[-c(7,9),]

plot_df2 <- rel_abun_act_in_denit %>%
  pivot_longer(
    cols = 3:50,
    names_to = "sample",
    values_to = "rel_abundance"
  )

ggplot(plot_df2, aes(x = gene, y = rel_abundance, fill = gene)) +
  geom_col() +
  facet_wrap(~ sample) +   # fixed y-axis across facets
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

##inactive averages
picrust_inactive<- picrust_metadata[
  !(picrust_metadata$sample_type %in% c("a+i", "bulk","a")),]
abundance_inactive<- abundance_file[, colnames(abundance_file) %in% picrust_inactive$Sample]
col_totals <- colSums(abundance_inactive, na.rm = TRUE)

rel_abun_in <- sweep(abundance_inactive, 2,
                         ifelse(col_totals == 0, NA, col_totals),
                         FUN = "/")
rel_abun_in= cbind(abundance_file[,1],rel_abun_in)
rel_abun_in_denit = rel_abun_in[rel_abun_in$`#NAME` %in% denit_KO$KO, ]
rel_abun_in_denit <- rel_abun_in_denit[match(denit_KO$KO, rel_abun_in_denit$`#NAME`), ]
rel_abun_in_denit= data.frame(gene= denit_KO[,1], rel_abun_in_denit, check.names = FALSE)
rel_abun_in_denit=rel_abun_in_denit[-c(7,9),]

inact_denit_avg=rel_abun_in_denit[,-c(1,2)]

samples_by_time <- split(
  picrust_inactive$Sample,
  picrust_inactive$`time point`
)




for(tp in names(samples_by_time)) {
  
  sample_cols <- intersect(
    samples_by_time[[tp]],
    colnames(inact_denit_avg)
  )
  
  inact_denit_avg[[paste0("mean_", tp)]] <-
    rowMeans(inact_denit_avg[, sample_cols, drop = FALSE], na.rm = TRUE )
  
  inact_denit_avg[[paste0("sd_", tp)]] <-
    apply(inact_denit_avg[, sample_cols, drop = FALSE],
          1,
          sd,
          na.rm = TRUE)
}

inact_denit_avg=data.frame(gene=rel_abun_act_in_denit[,c(1)], inact_denit_avg[,c(24:33)])

#plotting
# active dataframe
act_long <- act_denit_avg |>
  pivot_longer(
    cols = -gene,
    names_to = c(".value", "timepoint"),
    names_sep = "_"
  ) |>
  mutate(group = "Active")

# inactive dataframe
inact_long <- inact_denit_avg |>
  pivot_longer(
    cols = -gene,
    names_to = c(".value", "timepoint"),
    names_sep = "_"
  ) |>
  mutate(group = "Inactive")

# combine
plot_df <- bind_rows(act_long, inact_long)

# plot
plot_df <- plot_df |>
  mutate(group = factor(group, levels = c("Inactive", "Active")))
x=c(1,2,3,4,5)
y=c(3,5.67,8.42,12,15)

t=data.frame(timepoint=x, t=y)
plot_df$t <- t$t[
  match(plot_df$timepoint, t$timepoint)
]



plot_df$gene <- factor(
  plot_df$gene,
  levels = c(
    "nifH",
    "amoB",
    "amoC",
    "hao",
    "nrfA",
    "napA",
    "narG",
    "nirK",
    "nirS",
    "norB",
    "nosZ"))

ggplot(plot_df,
       aes(x = t,
           y = mean,
           color = group,
           group = interaction(group, gene))) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.2) +
  geom_errorbar(
    aes(ymin = mean - sd,
        ymax = mean + sd),
    width = 0.15,
    alpha = 0.5 ) +
  facet_wrap(~gene) +
  scale_color_manual(
    values = c(
      Inactive = "steelblue",
      Active = "blue"
    ),
    name = "State") +
  labs(
    x = "Time Point",
    y = "Predicted Mean Relative Abundance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

rel_abun_long <- rel_abun_act_in_denit %>%
  pivot_longer(
    cols = -c(gene,`#NAME`),
    names_to = "Sample",
    values_to = "rel_abundance"
  )
rel_abun_long$time_point=rep(picrust_metadata_trimmed$`time point`, times=11)
rel_abun_long$sample_type=rep(picrust_metadata_trimmed$sample_type, times=11)

picrust_metadata_trimmed$Sample <- as.character(picrust_metadata_trimmed$Sample)
rel_abun_long$Sample <- as.character(rel_abun_long$Sample)

x=c(1,2,3,4,5)
y=c(3,5.67,8.42,12,15)

t=data.frame(time_point=x, t=y)
rel_abun_long$t <- t$t[
  match(rel_abun_long$time_point, t$time_point)
]


#### Models
rel_abun_in_denit
rel_abun_act_denit

stat_df=rel_abun_long
stat_df$sample_type <- factor(
  stat_df$sample_type,
  levels = c("i", "a")
)
pseudo <- min(stat_df$rel_abundance[stat_df$rel_abundance > 0]) / 2

napA=lm(log(rel_abundance+pseudo) ~ sample_type * t, data = subset(stat_df, gene == "napA"))
narG=lm(log(rel_abundance+pseudo) ~ sample_type * t, data = subset(stat_df, gene == "narG"))
nifH=lm(log(rel_abundance+pseudo) ~ sample_type * t, data = subset(stat_df, gene == "nifH"))
nirS=lm(log(rel_abundance+pseudo) ~ sample_type * t, data = subset(stat_df, gene == "nirS"))
nirK=lm(log(rel_abundance+pseudo) ~ sample_type * t, data = subset(stat_df, gene == "nirK"))
norB=lm(log(rel_abundance+pseudo) ~ sample_type * t, data = subset(stat_df, gene == "norB"))
nosZ=lm(log(rel_abundance+pseudo) ~ sample_type * t, data = subset(stat_df, gene == "nosZ"))
amoB=lm(log(rel_abundance+pseudo) ~ sample_type * t, data = subset(stat_df, gene == "amoB"))
amoC=lm(log(rel_abundance+pseudo) ~ sample_type * t, data = subset(stat_df, gene == "amoC"))
hao=lm(log(rel_abundance+pseudo) ~ sample_type * t, data = subset(stat_df, gene == "hao"))
nrfA=lm(log(rel_abundance+pseudo) ~ sample_type * t, data = subset(stat_df, gene == "nrfA"))

summary(napA)
summary(narG)
summary(nifH)
summary(nirS)
summary(nirK)
summary(norB)
summary(nosZ)
summary(amoB)
summary(amoC)
summary(hao)
summary(nrfA)


###ASV of interest
#
ASV_KO2KEGG= ASV_KO |>
  column_to_rownames(var = colnames(ASV_KO)[1]) |>
  t() |>
  as.data.frame() |>
  rownames_to_column(var = colnames(ASV_KO)[1])
ASV_KO2KEGG=ko2kegg_abundance(ASV_KO2KEGG)
#
ASV_KO=readRDS("ASV_KO.rds")

colnames(ASV_KO) = sub("^ko:", "", colnames(ASV_KO))

ASV_KO= ASV_KO |>
  column_to_rownames(var = colnames(ASV_KO)[1]) |>
  t() |>
  as.data.frame() |>
  rownames_to_column(var = colnames(ASV_KO)[1])


ASV_KO2=ASV_KO[,-1]
col_totals <- colSums(ASV_KO2, na.rm = TRUE)

rel_abun_ASV <- sweep(ASV_KO2, 2,
                      ifelse(col_totals == 0, NA, col_totals),
                      FUN = "/")
rel_abun_ASV= cbind(ASV_KO[,1],rel_abun_ASV)

rel_abun_ASV_denit = rel_abun_ASV[rel_abun_ASV$`ASV_KO[, 1]` %in% denit_KO$KO, ]

rel_abun_ASV_denit <- rel_abun_ASV_denit[match(denit_KO$KO, rel_abun_ASV_denit$`ASV_KO[, 1]`), ]
rel_abun_ASV_denit= data.frame(gene= denit_KO[,1], rel_abun_ASV_denit, check.names = FALSE)
rel_abun_ASV_denit=rel_abun_ASV_denit[-c(7,9),]

plot_df <- rel_abun_ASV_denit %>%
  pivot_longer(
    cols = 3:11,
    names_to = "sample",
    values_to = "rel_abundance"
  )

ggplot(plot_df, aes(x = gene, y = rel_abundance, fill = gene)) +
  geom_col() +
  facet_wrap(~ sample) +   # fixed y-axis across facets
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#across all the predicted KOs per ASV, this is the relative abundance of this KO 







