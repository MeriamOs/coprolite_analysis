#Script written by Meriam van Os
#Used for downstream analysis of metagenomic shotgun data from kuri (dog) palaeofaeces
#Uploaded 18/09/2025

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(stringr)
library(scales)

### Dog reads

read_length_dog <- readr::read_csv(file = 'final_merged_dog_dna.csv')
dim(read_length_dog)
# [1] 82 17

# Reshape data
read_long <- read_length_dog %>%
  pivot_longer(cols = starts_with("MS"),    
               names_to = "Sample",         
               values_to = "number")  

head(read_long)

ggplot(read_long, aes(x = read_length, 
                      y = number, color = Sample)) +
  geom_line() +
  labs(x = "Read length (bp)", y = "Number of reads", 
       title = "Read length distribution of host mapped reads across samples") +
  scale_x_continuous(breaks = seq(30, 110, by = 10), limits = c(30, 110)) +
  scale_y_continuous(labels = scales::comma) +
  theme(
    plot.title = element_text(size = 18, face = "bold", , hjust = 0.5),
    legend.title = element_text(size = 16), 
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14))

#ggsave("dog_reads_fragment_length_plot.png", width = 12, height = 8)

### Microbial reads

read_length_microbial <- readr::read_csv(file = 'final_microbial.csv')

read_long_m <- read_length_microbial %>%
  pivot_longer(cols = starts_with("MS"),  # Gather all columns that start with "MS"
               names_to = "Sample",           # Create a new column "MS" to hold names
               values_to = "number")        # Create a new column "coverage" to hold values

# View reshaped data
head(read_long_m)

ggplot(read_long_m, aes(x = read_length, 
                      y = number, color = Sample)) +
  geom_line() +
  labs(x = "Read length (bp)", y = "Number of reads", 
       title = "Read length distribution of unmapped reads across samples") +
  scale_x_continuous(breaks = seq(30, 110, by = 10), limits = c(30, 110)) +
  scale_y_continuous(labels = scales::comma) +
  theme(
    plot.title = element_text(size = 18, face = "bold", , hjust = 0.5),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12),
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12))

#ggsave("microbial_reads_fragment_length_plot.png", width = 12, height = 6)


# Combine results into a new dataframe
df_frac <- data.frame(read_length = read_length_dog$read_length)

# Add summed values for each sample and calculate fractions
for (sample in colnames(read_length_dog)[-1]) {  
  
  # Add the sum for the current sample from df1 and df2
  df_frac[[paste0("sum_", sample)]] <- read_length_dog[[sample]] + read_length_microbial[[sample]]
  
  # Calculate the fraction of df1 for the current sample
  df_frac[[paste0("fraction_df1_", sample)]] <- read_length_dog[[sample]] / df_frac[[paste0("sum_", sample)]]
  
  # Calculate the fraction of df2 for the current sample
  df_frac[[paste0("fraction_df2_", sample)]] <- read_length_microbial[[sample]] / df_frac[[paste0("sum_", sample)]]
}

df_final <- dplyr::select(df_frac, read_length, contains("fraction_"))

#############################
### fraction of dog reads ###
#############################

frac_long_dog <- df_final %>%
  pivot_longer(cols = starts_with("fraction_df1"),  
               names_to = "sample",          
               values_to = "number")     

# View reshaped data
head(frac_long_dog)

ggplot(frac_long_dog, aes(x = read_length, 
                      y = number, color = sample)) +
  geom_line() +
  labs(x = "Read length (bp)", y = "Fraction of total number of reads", 
       title = "Fraction of host mapped fragment length distribution") +
  scale_x_continuous(breaks = seq(30, 110, by = 10), limits = c(30, 110)) +
  scale_color_discrete(labels = c("MS11669", "MS11670", "MS11673", "MS11674",
             "MS11675", "MS11676", "MS11677", "MS11678",
             "MS11679", "MS11683", "MS11684", "MS11686", 
             "MS11770", "MS11771", "MS11774", "MS11775")) +
  theme(
    plot.title = element_text(size = 20, face = "bold", , hjust = 0.5),
    legend.title = element_text(size = 16), 
    legend.text = element_text(size = 16),
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16))

#ggsave("dog_reads_fraction_plot.png", width = 12, height = 7)

###################################
### fraction of microbial reads ###
###################################

frac_long_mic <- df_final %>%
  pivot_longer(cols = starts_with("fraction_df2"),  # Gather all columns that start with "MS"
               names_to = "sample",           # Create a new column "MS" to hold names
               values_to = "number")        # Create a new column "coverage" to hold values

# View reshaped data
head(frac_long_mic)

ggplot(frac_long_mic, aes(x = read_length, 
                          y = number, color = sample)) +
  geom_line() +
  labs(x = "Read length", y = "Fraction", 
       title = "Fraction of unmapped read length distribution") +
  scale_x_continuous(breaks = seq(30, 110, by = 10), limits = c(30, 110)) +
  scale_color_discrete(labels = c("MS11669", "MS11670", "MS11673", "MS11674",
                                  "MS11675", "MS11676", "MS11677", "MS11678",
                                  "MS11679", "MS11683", "MS11684", "MS11686", 
                                  "MS11770", "MS11771", "MS11774", "MS11775")) +
  theme()

#ggsave("microbial_read_fraction_plot.png", width = 12, height = 6)

###################################
### total number of read length ###
###################################

df_sum <- dplyr::select(df_frac, read_length, contains("sum_"))

sum_long <- df_sum %>%
  pivot_longer(cols = starts_with("sum"),  # Gather all columns that start with "MS"
               names_to = "sample",           # Create a new column "MS" to hold names
               values_to = "number")        # Create a new column "coverage" to hold values

# View reshaped data
head(sum_long)

ggplot(sum_long, aes(x = read_length, 
                          y = number, color = sample)) +
  geom_line() +
  labs(x = "Read length", y = "Number of reads total", 
       title = "Total number of read length distribution") +
  scale_x_continuous(breaks = seq(30, 110, by = 10), limits = c(30, 110)) +
  scale_color_discrete(labels = c("MS11669", "MS11670", "MS11673", "MS11674",
                                  "MS11675", "MS11676", "MS11677", "MS11678",
                                  "MS11679", "MS11683", "MS11684", "MS11686", 
                                  "MS11770", "MS11771", "MS11774", "MS11775")) +
  scale_y_continuous(labels = scales::comma) +
  theme()

#ggsave("total_reads_length_distribution_plot.png", width = 12, height = 6)
