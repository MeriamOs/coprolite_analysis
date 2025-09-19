#Script written by Meriam van Os
#Used for downstream analysis of metagenomic shotgun data from kuri (dog) palaeofaeces
#Uploaded 18/09/2025


library(tibble)
library(dplyr)
library(readr)
library(forcats)
library(ggplot2)
library(plotly)
library(stringr)

###################
### GLOBAL DOGS ###
###################

setwd("~/Documents/Rstudio/r-tidyverse")

# Read the pca file
data_pca <- read_csv("merged_kuri_souilmi_dingo_rstudio.csv")
data_pca <- read_csv("merged_kuri_souilmi_dingo_nowolves_rstudio.csv") 

# Read the eigenvalues file
eigs <- scan("merged_kuri_souilmi_dingo.eval")
eigs <- scan("merged_kuri_souilmi_dingo_nowolves.eval")

total_var <- sum(eigs) # Calculate total variance

percent_var <- eigs / total_var * 100 # Calculate percent variance explained

# Put into a tidy dataframe
explained <- tibble(
  PC = paste0("PC", seq_along(eigs)),
  Eigenvalue = eigs,
  Percent = percent_var,
  Cumulative = cumsum(percent_var))

# Print first 10 PCs
print(head(explained, 10))

cols(Individual = col_character(), PC1 = col_double(), 
     PC2 = col_double(), PC3 = col_double(), PC4 = col_double(), 
     Population = col_character(), population2 = col_character())

data_meta <- read_csv("dogs_pca_population_list.csv")
cols(Population = col_character(), colorNr = col_character(), symbolNr2 = col_double(), alpha = col_double())

data_combined <- left_join(data_pca, data_meta) %>% print()
data_combined <- data_combined %>% mutate(Population = as_factor(Population)) %>% print()

list_colour <- data_combined %>% 
  dplyr::select(Population, colorNr) %>% 
  deframe()

list_shape <- data_combined %>% 
  dplyr::select(Population, symbolNr2) %>% 
  deframe()

list_alpha <- data_combined %>% 
  dplyr::select(Population, alpha) %>% 
  deframe()

ggplot(data = data_combined %>% arrange(Population == "BorneoVillage"), 
       aes(x = PC1, y = PC3, colour = Population, shape = Population,
                             alpha = Population)) +
          scale_alpha_manual(values = list_alpha) +
          scale_colour_manual(values = list_colour,
            # limits = c(
            #            "Wolf", 
            #            "Village", "Asia Village", "Breed", 
            #            "Dingo_ancient",
            #            "Dingo", "NewGuineaSingingDog", "Kuri"),
            # labels = c(
            #            "Wolf", 
            #            "Village dogs", "Asian Village", "Modern breeds", 
            #            "Dingo ancient",
            #            "Dingo","NG Singing dog", "Kuri")) 
            )+
        scale_shape_manual(values = list_shape,
            # limits = c(
            #            "Wolf", 
            #              "Village", "Asia Village", "Breed", 
            #            "Dingo_ancient",
            #              "Dingo", "NewGuineaSingingDog", "Kuri"),
            # labels = c(
            #            "Wolf", 
            #              "Village dogs", "Asian Village", "Modern breeds", 
            #            "Dingo ancient",
            #              "Dingo","NG Singing dog", "Kuri")
            ) +
#  guides(colour = guide_legend(ncol = 2)) +
  xlab(paste0("PC1 (", round(explained$Percent[explained$PC == "PC1"], 1), "%)")) +
#  ylab(paste0("PC2 (", round(explained$Percent[explained$PC == "PC2"], 1), "%)")) +
  ylab(paste0("PC3 (", round(explained$Percent[explained$PC == "PC3"], 1), "%)")) +
  geom_point(size = 3) +
  theme(legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        legend.position = "bottom",  # Move the legend to the bottom
        legend.direction = "horizontal",
#        plot.title = element_text(size = 16PCA_nuclear_DNA_subset_dogs_shrink_lsq, face = "bold", , hjust = 0.5),
        axis.text.x = element_text(size = 16), 
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 16))

#ggsave("PCA_nuclear_DNA_global_dogs_no_wolves.png", width = 10, height = 6)
ggsave("PCA_nuclear_DNA_global_dogs_ancient.png", width = 10, height = 6)

##########################
### SUBSET GLOBAL DOGS ###
##########################

data_pca <- read_csv("merged_kuri_souilmi_dingo_subset_shrink2_rstudio.csv") %>% 
#  filter(Population != "Dingo_ancient") %>%
  mutate(Population = if_else(str_detect(Population, "^China") | 
                              Population == "ChongqingDog", "ChinaVillage", Population))

cols(Individual = col_character(), PC1 = col_double(), 
     PC2 = col_double(), PC3 = col_double(), 
     PC4 = col_double(), Population = col_character())

eigs <- scan("merged_kuri_souilmi_dingo_subset_shrink.eval")

total_var <- sum(eigs) # Calculate total variance

percent_var <- eigs / total_var * 100 # Calculate percent variance explained

# Put into a tidy dataframe
explained <- tibble(
  PC = paste0("PC", seq_along(eigs)),
  Eigenvalue = eigs,
  Percent = percent_var,
  Cumulative = cumsum(percent_var))

# Print first 10 PCs
print(head(explained, 10))

data_meta <- read_csv("dogs_pca_population_list_subset.csv")
cols(Population = col_character(), colorNr = col_character(), symbolNr = col_double())

data_combined <- left_join(data_pca, data_meta) %>% print()
data_combined <- data_combined %>% mutate(Population = as_factor(Population)) %>% print()

#data_combined <- data_combined %>% filter(Population != "Dingo_ancient")

list_colour <- data_combined %>% 
  dplyr::select(Population, colorNr) %>% 
  deframe()

list_shape <- data_combined %>% 
  dplyr::select(Population, symbolNr) %>% 
  deframe()

ggplot(data_pca, aes(x = PC1, y = PC3, colour = Population)) +
  geom_point()

ggplot(data = data_combined %>% 
         arrange(Population == "VietnamVillage"), #to superpose Vietnam samples
              aes(x = PC1, y = PC4, colour = Population, shape = Population)) + 
       scale_colour_manual(values = list_colour, 
                 limits = c("Kuri", "Dingo",
                            "Dingo_ancient",
                            "NGSingingDog", "BorneoVillage", "VietnamVillage",
                            "ChinaVillage", "TaiwainVillage", "ChowChow", 
                            "XiasiDog", "Jindo"),
                 labels = c("Kuri", "Dingo",
                           "Dingo ancient", 
                           "New Guinea singing dog", "Village dog Borneo", 
                           "Village dog Vietnam", "Village dog China", "Village Taiwan", 
                           "Chow chow", "Xiasi", "Jindo")) + 
      scale_shape_manual(values = list_shape,
                 limits = c("Kuri", "Dingo",
                            "Dingo_ancient",
                            "NGSingingDog", "BorneoVillage", "VietnamVillage",
                            "ChinaVillage", "TaiwainVillage", "ChowChow", 
                            "XiasiDog", "Jindo"),
                 labels = c("Kuri", "Dingo",
                           "Dingo ancient", 
                            "New Guinea singing dog", "Village dog Borneo", 
                            "Village dog Vietnam", "Village dog China", "Village Taiwan", 
                            "Chow chow", "Xiasi", "Jindo")) + 
  xlab(paste0("PC1 (", round(explained$Percent[explained$PC == "PC1"], 1), "%)")) +
  ylab(paste0("PC4 (", round(explained$Percent[explained$PC == "PC4"], 1), "%)")) +
  geom_point(size = 5) +
  theme(legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(size = 16), 
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 16),
        legend.position = "bottom",  # Move the legend to the bottom
        legend.direction = "horizontal")

ggsave("PCA_nuclear_DNA_subset_dogs_shrink_lsq_PC1_PC4.png", width = 10, height = 6)

######################
### KURI ONLY PLOT ###
######################

### smartPCA ###
data_pca <- read_csv("angsd_bam_trimmed_SE_geno08_maf0.25_no_damage_rstudio.csv")
data_pca <- read_csv("angsd_bam_trimmed_SE_geno08_maf0.33_no_damage_rstudio.csv")
data_pca <- read_csv("angsd_bam_trimmed_SE_geno08_maf0.4_no_damage_rstudio.csv")
data_pca <- read_csv("angsd_bam_trimmed_SE_geno08_no_damage_shared_sites_rstudio.csv")
data_pca <- read_csv("pileupcaller.single.maf0.01_rstudio.csv")

cols(Individual = col_character(), 
     PC1 = col_double(), PC2 = col_double(), 
     PC3 = col_double(), PC4 = col_double(), Population = col_character())

data_meta <- read_csv("kuri_pca_population_list.csv")
cols(Population = col_character(), colorNr = col_character(), symbolNr = col_double())

data_combined <- left_join(data_pca, data_meta) %>% print()
data_combined <- data_combined %>% mutate(Population = as_factor(Population)) %>% print()

samples <- read.csv("angsd_numbers.csv", header = TRUE)
data_combined$Sex <- samples$Sex

list_colour <- data_combined %>% 
  dplyr::select(Population, colorNr) %>% 
  deframe()

list_shape <- data_combined %>% 
  dplyr::select(Population, symbolNr) %>% 
  deframe()

#pca_kuri <- 
ggplot(data = data_combined, # %>% arrange(desc(Population)), 
       aes(x = PC1, y = PC2, colour = Population, shape = Sex,
           text = paste("Sample:", Individual))) + 
        scale_colour_brewer(palette = "Set1",
#                      limits = c("LB", "KH", 
#                                 "WH_pre", 
#                                 "WH_post"),
                      limits = c("Long_Bay", "Kahukura", 
                                 "Whenua_Hou_pre", 
                                 "Whenua_Hou_post"),
                      labels = c("Long Bay pre-contact (n=4)", 
                                 "Kahukura pre-contact (n=4)", 
                                 "Whenua Hou pre-contact (n=5)", 
                                 "Whenua Hou post-contact (n=3)")) +
        geom_point(size = 5) +
        labs(title = "PCA of kuri nuclear DNA - smartpca pileupcaller") +
        theme(plot.title = element_text(size = 16, face = "bold", , hjust = 0.5),
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14), 
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14),
              legend.title = element_text(size = 16, face = "bold"),
              legend.text = element_text(size = 16),
              legend.position = "bottom",  
              legend.direction = "horizontal") +
        guides(colour = guide_legend(ncol = 2))

#interactive_plot <- ggplotly(pca_kuri, tooltip = "text")
#interactive_plot

ggsave("PCA_nuclear_DNA_kuri_smartPCA_pileupcaller.png", width = 10, height = 8)
ggsave("PCA_nuclear_DNA_kuri_smartPCA_angsd_no_damage_maf0.25.png", width = 10, height = 8)
ggsave("PCA_nuclear_DNA_kuri_smartPCA_angsd_no_damage_shared.png", width = 10, height = 8)


### PCAngsd ###
file <- "angsd_bam_trimmed_SE_geno08_maf0.4_no_damage"
file <- "angsd_bam_trimmed_SE_geno08_maf0.33_no_damage"
file <- "angsd_bam_trimmed_SE_geno08_maf0.25_no_damage"
file <- "angsd_bam_trimmed_SE_geno08_no_damage_shared_sites"
file <- "pileupcaller.single.pcangsd"
file <- "angsd_bam_trimmed_GL"

C <- as.matrix(read.table(paste0(file, ".cov"))) # Reads estimated covariance matrix
D <- as.matrix(read.table(paste0(file, ".selection"))) # Reads PC based selection statistics

samples <- read.csv("angsd_numbers.csv", header = TRUE)
samples_sorted <- samples[order(samples$angsd), ]

# Plot PCA plot
e <- eigen(C)
pca_data <- as.data.frame(e$vectors[,1:2])  # Extract the first two principal components

#pileupcaller order
pca_data$Sample <- samples$Sample
pca_data$Site <- samples$Site
pca_data$Sex <- samples$Sex

#angsd order
pca_data$Sample <- samples_sorted$Sample
pca_data$Site <- samples_sorted$Site
pca_data$Sex <- samples_sorted$Sex

#pca_kuri <- 
ggplot(pca_data, aes(x = V1, y = V2, colour = Site,
                     shape = Sex, text = paste("Sample:", Sample))) +
  geom_point(size = 5, alpha = 0.7) +
#  labs(title = "PCA of kuri nuclear DNA - PCAngsd maf 0.25", 
  labs(title = "PCA of kuri nuclear DNA - PCAngsd Genotype likelihoods", 
#  labs(title = "PCA of kuri nuclear DNA - PCAngsd pileupcaller", 
       x = "PC1", y = "PC2") +
  scale_colour_brewer(palette = "Set1",
                      limits = c("Long Bay", "Kahukura", 
                                 "Whenua Hou - PC", 
                                 "Whenua Hou - H"),
                      labels = c("Long Bay pre-contact (n=4)", 
                                 "Kahukura pre-contact (n=4)", 
                                 "Whenua Hou pre-contact (n=5)", 
                                 "Whenua Hou post-contact (n=3)")) +
  geom_point(size = 5) +
  theme(plot.title = element_text(size = 16, face = "bold", , hjust = 0.5),
                      axis.title.x = element_text(size = 16),
                      axis.text.x = element_text(size = 14), 
                      axis.title.y = element_text(size = 16),
                      axis.text.y = element_text(size = 14),
                      legend.title = element_text(size = 16, face = "bold"),
                      legend.text = element_text(size = 16),
                      legend.position = "bottom",  
                      legend.direction = "horizontal") +
  guides(colour = guide_legend(ncol = 2)) +
  guides(shape = guide_legend(ncol = 1))

#interactive_plot <- ggplotly(pca_kuri, tooltip = "text")
#interactive_plot

ggsave("PCA_nuclear_DNA_kuri_pileupcaller.png", width = 10, height = 8)
#ggsave("PCA_nuclear_DNA_kuri_PCAngsd_no_damage_maf0.25.png", width = 10, height = 8)
#ggsave("PCA_nuclear_DNA_kuri_PCAngsd_no_damage_shared.png", width = 10, height = 8)
ggsave("PCA_nuclear_DNA_kuri_PCAngsd_GL.png", width = 10, height = 8)


### Barplot of covered SNPs per sample for SNP panel dataset
## All samples have ~4.5% of non reference SNPs

SNP_stats <- read_csv("SNP_stats.csv")

ggplot(SNP_stats, aes(x = Sample, y = CoveredSites)) +
  geom_col() +
  geom_col(aes(x = Sample, y = NonReferenceSNPs, 
               fill = "Non reference SNPs")) +
  theme(plot.title = element_text(size = 16, 
                                  face = "bold", , hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 90), 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.position = "bottom",  
        legend.direction = "horizontal")

ggsave("Non_reference_SNPs_stats_barplots.png", width = 10, height = 8)
