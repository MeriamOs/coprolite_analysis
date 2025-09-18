#Script written by Meriam van Os
#Used for downstream analysis of metagenomic shotgun data from kuri (dog) palaeofaeces
#Uploaded 18/09/2025

library(tibble)
library(dplyr)
library(readr)
library(ape)
library(RColorBrewer)
library(ggtree)
library(treeio)
library(phangorn)
library(ggtreeExtra)
library(ggplot2)

setwd("~/Documents/Rstudio/r-tidyverse")

############################
### global dog phylogeny ###
############################

nj_tree <- read.newick("merged_kuri_souilmi_ibs_BOOSTER_tree.nwk")

# use only subset os samples for better visualisation
samples_to_keep <- read.csv("dogs_to keep_tree.csv", header = FALSE)[[1]]

samples_to_remove <- setdiff(nj_tree$tip.label, samples_to_keep)  

pruned_tree <- drop.tip(nj_tree, samples_to_remove)

metadata <- read_csv("Metadata_NJ_tree_global_dogs.csv")

tip_df <- data.frame(Sample = pruned_tree$tip.label, stringsAsFactors = FALSE)

tree_data <- dplyr::full_join(tip_df, metadata, by = "Sample")

groups <- unique(tree_data$Group)

# Generate a palette for all groups
palette <- setNames(
  brewer.pal(max(length(groups), 3), "Set3")[1:length(groups)],
  groups)

# Ensure colours are consistent
palette["Jackal"] <- "black"
palette["Wolf"] <- "darkgray"
palette["Basal breed"] <- "#FFFFB3"
palette["Modern breed"] <- "#8DD3C7"
palette["Kuri"] <- "#80B1D3"
palette["Dingo"] <- "#FDB462"
palette["New Guinea Singing Dog"] <- "#FB8072" 
palette["Indigenous dog China"] <- "#FCCDE5"
palette["Indigenous dog Vietnam"] <- "#B3DE69"
palette["Village dog Borneo"] <- "#BC80BD"
palette["Village dog China"] <- "#BEBADA"
palette["Village dog other"] <- "#CCEBC5"

pruned_tree$edge.length <- NULL

ggtree(pruned_tree, layout = "circular") %<+% tree_data +
  geom_tiplab(colour = "black", size = 4, offset = 0.5) +  
  geom_tippoint(aes(colour = Group), size = 3) +
  scale_colour_manual(values = palette) +
#  geom_nodelab(aes(label = sprintf("%.2f", as.numeric(label))), size = 3) +
  theme_tree2() +
#  ggtitle("NJ tree of global dogs") +
  theme(
        legend.position = "none",
        # legend.position = "right",
        # legend.key.size = unit(1, "lines"),
        # legend.text = element_text(size = 20),
        # legend.title = element_text(size = 24, face = "bold"),
        # legend.direction = "vertical",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
  guides(color = guide_legend(ncol = 1, 
                              override.aes = list(size = 5))) +
  geom_nodepoint(aes(subset = (as.numeric(label) > 0.9)), colour = "black")

ggsave("NJ_tree_circular_global_dogs.png", width = 12, height = 8)
ggsave("NJ_tree_circular_global_dogs_text_no_legend.png", width = 12, height = 12)

#################
### kuri only ###
#################

# Construct Neighbor-Joining tree
kuri_tree <- read.newick("pileupcaller.single.maf0.01_BOOSTER_tree.nwk")
kuri_tree <- read.newick("angsd_bam_trimmed_SE_geno08_maf0.4_no_damage_BOOSTER_tree.nwk")
kuri_tree <- read.newick("angsd_bam_trimmed_SE_geno08_no_damage_shared_sites_BOOSTER_tree.nwk")

site_info <- read.csv("sample_site_info.csv", header=TRUE)
site_info <- site_info[1:16,]

ggtree(kuri_tree, layout = "rectangular") %<+% site_info +
  geom_tiplab(colour = "black", size = 6, offset = 0.005) +  
  geom_tippoint(aes(colour = Site), size = 6) +
  scale_colour_brewer(palette = "Set1",
                    limits = c("Long_Bay", "Kahukura", 
                               "Whenua_Hou_pre", 
                               "Whenua_Hou_post"),
                    labels = c("Long Bay pre-contact (n=4)", 
                               "Kahukura pre-contact (n=4)", 
                               "Whenua Hou pre-contact (n=5)", 
                               "Whenua Hou post-contact (n=3)")) +
  theme_tree2() +
  ggtitle("NJ tree of kuri - angsd data geno 0.4") +
#  ggtitle("NJ tree of kuri - angsd data shared sites") +
#  ggtitle("NJ tree of kuri - pileupcaller") +
  xlim(0, max(kuri_tree$edge.length) + 0.05) +
  geom_nodelab(aes(label = sprintf("%.2f", as.numeric(label))), size = 3,
               nudge_x = 0.005, colour = "red4") +
  theme(
        legend.key.size = unit(1, "lines"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24, face = "bold"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "bottom",  # Move the legend to the bottom
        legend.direction = "horizontal") +
  guides(color = guide_legend(ncol = 1, 
                              override.aes = list(size = 5))) +
  geom_nodepoint(aes(subset = (as.numeric(label) > 0.9)), colour = "black", shape = 18, size = 3) +
  guides(colour = guide_legend(ncol = 2))

#ggsave("Midpoint_rooted_NJ_tree_kuri_pileupcaller.png", width = 10, height = 6)
#ggsave("Midpoint_rooted_NJ_tree_kuri_angsd_geno0.4.png", width = 10, height = 6)
#ggsave("Midpoint_rooted_NJ_tree_kuri_angsd_shared_sites.png", width = 10, height = 6)
ggsave("Midpoint_rooted_NJ_tree_angsd_kuri_angsd_no_damage.png", width = 10, height = 8)

### outgroup kuri ###

kuri_tree <- read.newick("kuri_outgroup_BOOSTER_tree.nwk")

site_info <- read.csv("sample_site_info.csv", header=TRUE)
site_info <- site_info[1:16,]

outgroup_tip <- which(kuri_tree$tip.label == "SiberianHusky01")

# Find the edge leading to that tip
outgroup_edge <- which(kuri_tree$edge[, 2] == outgroup_tip)

# Shrink that branch length (e.g., to 20% of its original length)
kuri_tree$edge.length[outgroup_edge] <- kuri_tree$edge.length[outgroup_edge] * 0.2

ggtree(kuri_tree, layout = "rectangular") %<+% site_info +
  geom_tiplab(colour = "black", size = 6, offset = 0.005) +  
  geom_tippoint(aes(colour = Site), size = 6) +
  scale_colour_brewer(palette = "Set1",
                      limits = c("Long_Bay", "Kahukura", 
                                 "Whenua_Hou_pre", 
                                 "Whenua_Hou_post"),
                      labels = c("Long Bay pre-contact (n=4)", 
                                 "Kahukura pre-contact (n=4)", 
                                 "Whenua Hou pre-contact (n=5)", 
                                 "Whenua Hou post-contact (n=3)")) +
  theme_tree2() +
#  ggtitle("NJ tree of kuri - angsd data geno 0.4") +
  xlim(0, max(kuri_tree$edge.length) + 0.125) +
  geom_nodelab(aes(label = sprintf("%.2f", as.numeric(label))), size = 3,
               nudge_x = 0.008, colour = "red4") +
  theme(
    legend.key.size = unit(1, "lines"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 24, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom",  # Move the legend to the bottom
    legend.direction = "horizontal") +
  guides(color = guide_legend(ncol = 1, 
                              override.aes = list(size = 5))) +
  geom_nodepoint(aes(subset = (as.numeric(label) > 0.8)), colour = "black", shape = 18, size = 3) +
  guides(colour = guide_legend(ncol = 2))

ggsave("Outgroup_NJ_tree_kuri_pileupcaller.png", width = 10, height = 6)
