#Script written by Meriam van Os
#Used for downstream analysis of metagenomic shotgun data from kuri (dog) palaeofaeces
#Uploaded 18/09/2025

library(scales)
library(ggplot2)
library(tidyverse)
library(dplyr)

setwd("~/Documents/Rstudio/r-tidyverse")

###########
### SFS ###
###########
 
SFS <- read.csv("SFS_per_site.csv") %>%
        select("Frequency", "All") #or select site "Whenua_Hou", "Long_Bay", "Kahukura"
colnames(SFS) <- c("Minor.Allele.Count", "Number.of.SNPs")

SFS <- SFS %>%
  mutate(Percentage = 100 * Number.of.SNPs / sum(Number.of.SNPs))

SFS.sub <- subset(SFS, Minor.Allele.Count > 0)

cols <- colnames(SFS)[-1] # Get column names (excluding the Frequency column)

# Initialize result vector
proportions <- sapply(cols, function(col) {
  sum_gt1 <- sum(SFS[SFS$Frequency > 1, col], na.rm = TRUE)
  sum_gt0 <- sum(SFS[SFS$Frequency > 0, col], na.rm = TRUE)
  sum_gt1 / sum_gt0 })

print(proportions)

sum(SFS$Number.of.SNPs) #138,738,633
sum(SFS.sub$Number.of.SNPs) #3,040,505

ggplot(SFS.sub, aes(x = Minor.Allele.Count, y = Percentage)) +
  geom_col(fill = "purple4") +
#  scale_x_continuous(breaks = SFS.sub$Minor.Allele.Count) +
  scale_x_continuous(breaks = sort(unique(c(0, SFS.sub$Minor.Allele.Count)))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Frequency of minor variant",
       x = "Minor Allele Count", y = "Percentage") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 12), 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14))

#ggsave("SFS_1-6.png", width = 8, height = 4)

### Map expected values ###

n <- 13 # n = 13 for the pre-contact kuri only

expected <- data.frame(
  Minor.Allele.Count = 1:(n - 1),
  Expected = 1 / (1:(n - 1)))

# Scale to match observed total SNPs
scale_factor <- sum(SFS.sub$Percentage) / sum(expected$Expected)
expected$Expected <- expected$Expected * scale_factor

all_x <- sort(unique(as.numeric(c(SFS.sub$Minor.Allele.Count, expected$Percentage))))

SFS.sub <- SFS.sub %>% filter(Minor.Allele.Count <= 8)
expected <- expected %>% filter(Minor.Allele.Count <= 8)

ggplot(SFS.sub, aes(x = Minor.Allele.Count, y = Percentage)) +
  geom_col(fill = "skyblue") +
  geom_point(data = expected, aes(x = Minor.Allele.Count, y = Expected),
             color = "black", size = 3) +
  scale_x_continuous(breaks = expected$Minor.Allele.Count) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Frequency of minor variant",
       x = "Minor Allele Count", y = "Frequency (%)") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 12), 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14))

#ggsave("SFS_1-8_percentage_expected.png", width = 6, height = 4)

############################
### Nucleotide diversity ###
############################

covered <- read_tsv("positions_counts_per_10kb_window.tsv")
sum(covered$N_SITES) #138,738,633
head(covered)
# A tibble: 6 × 4
# CHROM START END N_SITES

### Pi 
pi.all <- read.table("angsd_bam_trimmed_SE_no_post_04_haploid.windowed.pi", header = TRUE)

sum(pi.all$N_VARIANTS) #3,040,505

pi.all.sub <- subset(pi.all, N_VARIANTS > 0)
#pi.all.sub <- subset(pi.all, PI > 0.0015)

pi.sub.edit <- pi.all.sub %>%
  mutate(BIN_START = BIN_START - 1)

pi.cov <- left_join(pi.sub.edit, 
                    covered %>% select(CHROM, START, N_SITES),
                    by = c("CHROM", "BIN_START" = "START"))

pi.covered <- subset(pi.cov, N_SITES > 250)
# pi.cov.pos <- subset(pi.covered, PI > 0)

mean(pi.covered$PI)
median(pi.covered$PI)
mean(pi.covered$N_VARIANTS)
median(pi.covered$N_VARIANTS)

#hist(pi.all.sub$PI,br=20)

boxplot(pi.covered$PI,ylab="diversity")

ggplot(pi.covered, aes(x = PI)) +
  geom_histogram(bins = 20, fill = "skyblue", color = "black") +
  theme_minimal() +
  scale_y_continuous(labels = comma) +
  labs(title = "Histogram of nucleotide diversity (π) values across 10kp windows",
       x = "Nucleotide diversity (π)",
       y = "Frequency") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 10), 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 14))

#ggsave("Histogram_pi.png", width = 4, height = 4)

################
#### TajimaD ###
################

tajD.all <- read.table("angsd_bam_trimmed_SE_no_post_04_haploid.Tajima.D", header = TRUE)

sum(tajD.all$N_SNPS) #3,040,505

tajD.all.sub <- subset(tajD.all, N_SNPS > 0)

tajD.cov <- left_join(tajD.all.sub, 
                      covered %>% select(CHROM, START, N_SITES),
                      by = c("CHROM", "BIN_START" = "START"))

tajD.covered <- subset(tajD.cov, N_SITES > 250)
#tajD.cov.pos <- subset(tajD.covered, TajimaD > 0)

mean(tajD.covered$TajimaD)
median(tajD.covered$TajimaD)
mean(tajD.covered$N_SNPS)
median(tajD.covered$N_SNPS)

#hist(tajD.all.sub$TajimaD,br=20)

ggplot(tajD.covered, aes(x = TajimaD)) +
  geom_histogram(bins = 20, fill = "skyblue", color = "black") +
  theme_minimal() +
  scale_y_continuous(labels = comma) +
  labs(title = "Histogram of Tajima's D values across 10kp windows - Autosomes",
       x = "Tajima's D",
       y = "Frequency") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 10), 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 14))

#ggsave("Histogram_tajimaD.png", width = 4, height = 4)

