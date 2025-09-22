#Script written by Meriam van Os
#Used for downstream analysis of metagenomic shotgun data from kuri (dog) palaeofaeces
#Uploaded 18/09/2025

library(scales)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(patchwork) 

setwd("~/Documents/Rstudio/r-tidyverse")

###########
### SFS ###
###########

SFS <- read.csv("SFS_per_site.csv") %>%
  select("Frequency", "All", "Pre_contact") %>% 
  pivot_longer(cols = c("All", "Pre_contact"),
               names_to = "Dataset",
               values_to = "Number.of.SNPs")

colnames(SFS)[1] <- "Minor.Allele.Count" # Rename first column

### Calculate percentages per dataset
SFS <- SFS %>%
  group_by(Dataset) %>%
  mutate(Total_SNPs = sum(Number.of.SNPs, na.rm = TRUE),
    Percentage = 100 * Number.of.SNPs / Total_SNPs) %>%
  ungroup()

### Proportions 
SFS %>% group_by(Dataset) %>% 
  summarise( total_all = sum(Number.of.SNPs), 
             total_sub = sum(Number.of.SNPs[Minor.Allele.Count > 0]))
# A tibble: 2 × 3
#Dataset     total_all total_sub
#<chr>           <int>     <int>
#  1 All         211471754   4681508
#  2 Pre_contact 138738633   3040505

proportions <- SFS %>%
  group_by(Dataset) %>%
  summarise(
    prop = sum(Number.of.SNPs[Minor.Allele.Count > 1], na.rm = TRUE) /
      sum(Number.of.SNPs[Minor.Allele.Count > 0], na.rm = TRUE))
print(proportions)
# A tibble: 2 × 2
#Dataset       prop
#<chr>        <dbl>
#  1 All         0.0421
#  2 Pre_contact 0.0419

####################
### SFS plot #######
####################

# Subset to remove non-variant sites
SFS.sub <- filter(SFS, Minor.Allele.Count > 0)

### Expected values, using n = 16 for "All" and n = 13 for Pre_contact
expected_list <- SFS.sub %>%
  distinct(Dataset) %>%
  rowwise() %>%
  mutate(n = ifelse(Dataset == "All", 16, 13)) %>%
  do({ data.frame(
         Dataset = .$Dataset,
         Minor.Allele.Count = 1:(.$n - 1),
         Expected = 1 / (1:(.$n - 1))) }) %>%
  ungroup()

# Scale expected values to match each dataset’s observed total
scale_factors <- SFS.sub %>%
  group_by(Dataset) %>%
  summarise(scale_factor = sum(Percentage) / sum(expected_list$Expected[expected_list$Dataset == unique(Dataset)]))

expected_scaled <- expected_list %>%
  left_join(scale_factors, by = "Dataset") %>%
  mutate(Expected = Expected * scale_factor) %>%
  filter(Minor.Allele.Count <= 8)   # Keep same filtering as before

# Restrict observed to same Minor Allele Counts
SFS.sub <- filter(SFS.sub, Minor.Allele.Count <= 8)

ggplot(SFS.sub, aes(x = Minor.Allele.Count, y = Percentage, fill = Dataset)) +
  geom_col() +
  geom_point(data = expected_scaled,
             aes(x = Minor.Allele.Count, y = Expected),
             inherit.aes = FALSE,
             color = "black", size = 3) +
  scale_x_continuous(breaks = expected_scaled$Minor.Allele.Count) +
  scale_y_continuous(labels = scales::comma) +
  facet_wrap(~Dataset, ncol = 2) +
  labs(title = "Frequency of minor variant with expected values",
       x = "Minor Allele Count", y = "Frequency (%)") +
  scale_fill_manual(
    values = c("Pre_contact" = "purple4", "All" = "skyblue"),
    labels = c("Pre_contact" = "pre-contact samples (geno 0.4)",
               "All" = "all samples (geno 0.5)")) +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        legend.position = "bottom",  
        legend.direction = "horizontal")

ggsave("SFS_expected_percentage.png", width = 8, height = 4)

ggplot(SFS.sub, aes(x = Minor.Allele.Count, y = Number.of.SNPs, fill = Dataset)) +
  geom_col() +
  scale_x_continuous(breaks = expected_scaled$Minor.Allele.Count) +
  scale_y_continuous(labels = scales::comma) +
  facet_wrap(~Dataset, ncol = 2) +
  labs(title = "Frequency of minor variant with expected values",
       x = "Minor Allele Count", y = "Counts") +
  scale_fill_manual(
    values = c("Pre_contact" = "purple4", "All" = "skyblue"),
    labels = c("Pre_contact" = "pre-contact samples (geno 0.4)",
               "All" = "all samples (geno 0.5)")) +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        legend.position = "bottom",  
        legend.direction = "horizontal")

ggsave("SFS_counts.png", width = 8, height = 4)

###############################
### Nucleotide diversity Pi ###
###############################

# Load covered data for both
covered_04 <- read_tsv("positions_counts_per_10kb_window.tsv")
sum(covered_04$N_SITES) # [1] 138,738,633
covered_05 <- read_tsv("positions_counts_per_10kb_window_all_samples.tsv")
sum(covered_05$N_SITES) # [1] 211,471,754

# Load pi data
pi_04 <- read.table("angsd_bam_trimmed_SE_no_post_04_haploid.windowed.pi", header = TRUE)
pi_05 <- read.table("angsd_bam_trimmed_SE_all_05_haploid.windowed.pi", header = TRUE)

# Process in a function
process_pi <- function(pi, covered, label) {
  pi.sub <- subset(pi, N_VARIANTS > 0) %>%
    mutate(BIN_START = BIN_START - 1)
  
  pi.cov <- left_join(pi.sub,
    covered %>% select(CHROM, START, N_SITES),
    by = c("CHROM", "BIN_START" = "START"))
  
  pi.covered <- subset(pi.cov, N_SITES > 100)
  pi.covered$Dataset <- label
  return(pi.covered) }

pi_04.proc <- process_pi(pi_04, covered_04, "no_post_04")
mean(pi_04.proc$PI) # [1] 0.0001332216
median(pi_04.proc$PI) # [1] 0.000128307
pi_05.proc <- process_pi(pi_05, covered_05, "all_05")
mean(pi_05.proc$PI) # [1] 0.0001331843
median(pi_05.proc$PI) # [1] 0.000131035

pi.all.proc <- bind_rows(pi_04.proc, pi_05.proc)

# Plot side by side with facets
ggplot(pi.all.proc, aes(x = PI)) +
  geom_histogram(bins = 20, fill = "skyblue", color = "black") +
  theme_minimal() +
  scale_y_continuous(labels = comma) +
  labs(title = "Histogram of nucleotide diversity (π)",
       x = "Nucleotide diversity (π)", y = "Frequency") +
  facet_wrap(~Dataset, ncol = 2) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 10), 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 14))

################
#### TajimaD ###
################

tajD_04 <- read.table("angsd_bam_trimmed_SE_no_post_04_haploid.Tajima.D", header = TRUE)
tajD_05 <- read.table("angsd_bam_trimmed_SE_all_05_haploid.Tajima.D", header = TRUE)

process_tajD <- function(tajD, covered, label) {
  tajD.sub <- subset(tajD, N_SNPS > 0)
  tajD.cov <- left_join(tajD.sub,covered %>% 
                          select(CHROM, START, N_SITES),
                        by = c("CHROM", "BIN_START" = "START"))
  tajD.covered <- subset(tajD.cov, N_SITES > 50)
  tajD.covered$Dataset <- label
  return(tajD.covered) }

tajD_04.proc <- process_tajD(tajD_04, covered_04, "no_post_04")
mean(tajD_04.proc$TajimaD) # [1] -1.080473
median(tajD_04.proc$TajimaD) # [1] -1.11908
tajD_05.proc <- process_tajD(tajD_05, covered_05, "all_05")
mean(tajD_05.proc$TajimaD) # [1] -0.9422591
median(tajD_05.proc$TajimaD) # [1] -0.9700975

tajD.all.proc <- bind_rows(tajD_04.proc, tajD_05.proc)

ggplot(tajD.all.proc, aes(x = TajimaD)) +
  geom_histogram(bins = 20, fill = "skyblue", color = "black") +
  theme_minimal() +
  scale_y_continuous(labels = comma) +
  labs(title = "Histogram of Tajima's D",
       x = "Tajima's D", y = "Frequency") +
  facet_wrap(~Dataset, ncol = 2) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 10), 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 14))

#####################
### Combined plot ###
#####################

df_all <- bind_rows(
  pi.all.proc %>% 
    transmute(Value = PI, Statistic = "Pi", Dataset),
  tajD.all.proc %>% 
    transmute(Value = TajimaD, Statistic = "TajimaD", Dataset))

ggplot(df_all, aes(x = Dataset, y = Value, fill = Dataset)) +
  geom_violin(outlier.shape = NA, alpha = 0.7) +
  facet_wrap(~Statistic, scales = "free_y") +
  scale_fill_manual(
    values = c("no_post_04" = "purple4", "all_05" = "skyblue"),
    labels = c("no_post_04" = "pre-contact samples (geno 0.4)",
               "all_05" = "all samples (geno 0.5)")) +
  theme_minimal() +
  labs(x = "Dataset", y = "Value") +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 16, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = "bottom",  
        legend.direction = "horizontal")

#ggsave("Violinplots_TajD_Pi_two_datasets.png", width = 8, height = 4)
