# Software Tools - Anopheles Classification Activity ----

##***************************
## Software Tools - Anopheles Classification Activity with Random Forest
## Script: anopheles_scrambled.R
##
## Student: Stephanie Saab
## Author: Karl Cottenie
## Updated: October 23rd, 2025 by Stephanie Saab
## Description: This script fetches Anopheles COI gene sequences from GenBank,
## and builds a Random Forest classifier to predict species based on their 
## dinucleotide frequencies. Data are split into training (80%) and validation
## (20%) datasets to assess the model's generalization.
##***************************

# PART 1: LOAD PACKAGES & SETUP ----
library(tidyverse)
library(viridis)
library(Biostrings)
library(randomForest)
library(rentrez)
library(conflicted)

# Adjust for conflicts
conflict_prefer("filter", "dplyr")

# PART 2: DOWNLOAD & IMPORT DATA ----

# Get the GenBank search hits based on the search terms of interest (database = nucleotide, organism = anopheles, gene = COI)
df_anopheles_search <- entrez_search(db = "nuccore", term = "anopheles [ORGN] AND COI [gene]", use_history = T)
df_anopheles_search # Check that it resulted in 10730 hits, 20 IDs and a web_history object

# Read in the sequence data as a DNA StringSet from data GenBank database
st_anopheles <- readDNAStringSet("../data/sequence.fasta")
st_anopheles # Check that it is a DNAStringSet object of length 10730

#Convert to dataframe and check it has both titles of entries and sequences
df_anopheles <- data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
head(df_anopheles)

# PART 3: CALCULATE DINUCLEOTIDE FREQUENCIES ----

# Add column with nucleotides as a DNAStringSet class to use Biostrings package functions
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2)
head(df_anopheles) # Check that the nucleotides2 column is class 'DNAStringSet'

# Calculate the dinucleotide frequencies and append them to our dataframe,
df_anopheles <- cbind(
  df_anopheles,
  as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE))
)
head(df_anopheles)

# Remove the redundant nucleotides2 column and convert df to tibble
df_anopheles <- df_anopheles %>%
  select(-nucleotides2) %>%
  as_tibble()

# PART 4: CLEAN DATA & CREATE TRAIN-VALIDATION SPLITS ----

# Extract species name (2nd and 3rd words) and unique ID (1st word) from the sequence title, to group and filter sequences later.
df_anopheles <- df_anopheles %>%
  mutate(species_name = word(title, 2L, 3L)) %>%
  mutate(unique_id = word(title, 1L))
df_anopheles # Check that updated tibble has 10730 rows x 20 columns

# Get the unique entries
df_anopheles$unique_id %>%
  unique() %>%
  length() # Check that it is 10730 rows

# Create a vector list of frequently occurring species names ( > 200 times)
ls_vectors <- df_anopheles %>%
  group_by(species_name) %>%
  count() %>%
  arrange(desc(n)) %>% # View()
  filter(n > 200) %>%
  pull(species_name)
ls_vectors #Check that it has 14 unique species names

# Create a df of entries with only frequently occurring species (n>200)
df_anopheles_vector <-
  df_anopheles %>%
  filter(species_name %in% ls_vectors)

# Check that the species counts all > 200, expect a tibble of 14 rows x 2 columns
df_anopheles_vector %>%
  group_by(species_name) %>%
  count()

# Remove unnecessary objects
rm(df_anopheles, st_anopheles)

# Get the smallest species count to prevent oversampling later
smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample # Check it returns 220

# Create a balanced validation dataset by sampling about 20% of the smallest species group from each species.
# Use set.seed() for reproducibility. This avoids oversampling or biasing species representation by taking the same number of rows per species.
set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation

# Check that the validation dataset is balanced, each species should have the same number of samples (expected 44)
df_anopheles_validation %>%
  group_by(species_name) %>%
  count()

# Create a balanced training dataset that takes about 80% of the smallest species group from each species.
# Exlude rows that have a unique_id in the validation dataset to prevent overlap.
set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!unique_id %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

# Check that the sampling is balanced, count should be the same for each species (expect 176)
df_anopheles_training %>%
  group_by(species_name) %>%
  count()

# PART 5: TRAIN CLASSIFICATION MODEL: SPECIES ----

# Find which columns that have dinucleotide frequencies
names(df_anopheles_training)

# Building a species classifier.
# Use the dinucleotide frequencies (col 3-18) as predictors, use species name as response variable. Create 100 trees as we have a fairly large dataset.
anopheles_class <- randomForest(
  x = df_anopheles_training[, 3:18],
  y = base::as.factor(df_anopheles_training$species_name),
  ntree = 100, importance = TRUE
)

anopheles_class # Display the Random Forest model summary (i.e. OOB error, confusion matrix, number of trees)
plot(anopheles_class) # Plot OOB error rate vs. number of trees to assess model stability
varImpPlot(anopheles_class) # Plot variable importance to identify which dinucleotide frequencies contribute most to classification

# PART 6: VALIDATE CLASSIFICATION MODEL: SPECIES ----
# Use the trained model to predict the species names for samples in our unseen data (validation dataset) based on dinucleotide frequencies.
predict_validation <- predict(
  anopheles_class,
  df_anopheles_validation[, c(19, 3:18)]
)
predict_validation # View the predicted species for each sample

# Create a confusion matrix of the observed species names vs. the predicted species names to evaluate how the model did on unseen data
table(
  observed = df_anopheles_validation$species_name,
  predicted = predict_validation
)
