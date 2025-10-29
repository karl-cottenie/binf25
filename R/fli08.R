## ***************************
# Title: Anopheles COI Classification using Dinucleotide Frequencies
# Date: <2025-10-23>
#
# Purpose:
#   This script retrieves Anopheles COI gene sequences, processes
#   them into a tidy format, computes dinucleotide frequency features,
#   and trains a Random Forest classifier to distinguish species.
#
## ***************************

## ---- 1. Setup ----

# Packages used
library(tidyverse)
library(stringr)
library(Biostrings)
library(randomForest)
library(rentrez)

set.seed(100)  # global seed

# Search GenBank for Anopheles COI sequences
df_anopheles_search = entrez_search(db = "nuccore", 
                                    term = "anopheles [ORGN] AND COI [gene]",
                                    use_history = T)
df_anopheles_search


## ---- 2. Read FASTA ----
# Example sequence (not exist)
st_anopheles = readDNAStringSet("../data/sequence.fasta")
st_anopheles 


## ---- 3. Create a dataframe from FASTA data ----
# Convert the DNAStringSet object into a dataframe for easier handling.
df_anopheles = data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
head(df_anopheles)


## ---- 4. Add DNAStringSet column ----
# Store each sequence again as a DNAStringSet object so we can calculate
# dinucleotide frequencies using Biostrings functions.
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2)
head(df_anopheles)


## ---- 5. Calculate dinucleotide frequencies ----
# Compute the proportion of each of the 16 possible dinucleotides
# (AA, AC, AG, ..., TT) in every sequence. These serve as numerical
# features for machine learning.
df_anopheles = cbind(df_anopheles, 
                     as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, 
                                                         as.prob = TRUE)))
head(df_anopheles)


## ---- 6. Clean up data format ----
# Remove the large DNAStringSet column and convert to tibble for easier handling.
df_anopheles = df_anopheles %>% 
  select(-nucleotides2) %>% 
  as_tibble()

## ---- 7. Extract metadata from FASTA headers ----
# Parse the FASTA title string to extract:
# - species_name: combines 2nd and 3rd words (e.g., "Anopheles gambiae")
# - unique_id: uses the first word (record identifier)
df_anopheles = df_anopheles %>% 
  mutate(species_name = word(title, 2L, 3L)) %>% 
  mutate(unique_id = word(title, 1L)) 
df_anopheles


## ---- 8. Count unique sequence IDs ----
# Quick check to ensure all unique identifiers were parsed correctly.
df_anopheles$unique_id %>% 
  unique() %>% 
  length()


# ---- 9. Identify well-represented species ----
# Count how many sequences each species has, sort descending,
# and keep only species with >200 sequences (for balanced model training).
ls_vectors = df_anopheles %>% 
  group_by(species_name) %>% 
  count() %>% 
  arrange(desc(n)) %>% #View()
  filter(n > 200) %>% 
  pull(species_name)
ls_vectors


## ---- 10. Filter dataset to keep only selected species ----
# Keep only the species with sufficient sequence counts.
df_anopheles_vector = 
  df_anopheles %>% 
  filter(species_name %in% ls_vectors)

## ---- 11. Check how many sequences per species remain ----
# Sanity check class balance in the validation set:
# Count sequences per species to confirm the stratified split worked as intended.
df_anopheles_vector %>% 
  group_by(species_name) %>% 
  count()


## ---- 12. Free up memory ----
# Remove large intermediate objects that are no longer needed.
rm(df_anopheles, st_anopheles)


## ---- 13. Determine smallest sample size among species ----
# This ensures balanced sampling across all species.
smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample


## ---- 14. Create validation dataset ----
# Randomly sample 20% of the smallest class size for each species.
set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation


## ---- 15. Check class distribution in validation set ----
# Confirm that each species contributes the same number of samples.
df_anopheles_validation %>% 
  group_by(species_name) %>% 
  count()


## ---- 16. Create training dataset ----
# From the remaining data, randomly sample 80% of the smallest class size per species.
set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

df_anopheles_training %>% 
  group_by(species_name) %>% 
  count()

## ---- 17. Inspect feature columns ----
names(df_anopheles_training)


## ---- 18. Train Random Forest model ----
# Use the 16 dinucleotide frequencies as predictors (columns 3:18).
# The response variable is species_name. We use 100 trees for stability.
anopheles_class <- randomForest(x = df_anopheles_training[, 3:18],
                                y = as.factor(df_anopheles_training$species_name), 
                                ntree = 100, importance = TRUE)

## ---- 19. View and visualize the trained model ----
# anopheles_class: shows overall model accuracy (OOB error) and per-class errors.
# plot(): shows error rate vs. number of trees.
# varImpPlot(): ranks which dinucleotides are most informative.
anopheles_class
plot(anopheles_class)
varImpPlot(anopheles_class)


## ---- 20. Predict on validation dataset ----
# Use the trained model to classify the held-out validation sequences.
predict_validation <- predict(anopheles_class, 
                              df_anopheles_validation[, c(19, 3:18)])
predict_validation


## ---- 21. Confusion matrix ----
# Compare predicted vs. true species names to evaluate model performance.
table(observed = df_anopheles_validation$species_name, 
      predicted = predict_validation)

