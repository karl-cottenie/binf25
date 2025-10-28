# MBINF 6210
# Matei Dan-Dobre
# October 23rd, 2025

#### LOAD PACKAGES####

library(rentrez)
library(Biostrings)
library(tidyverse)
library(randomForest)
# conflict_prefer("filter", "dplyer)

#### OBTAIN DATA####
# Find on NCBI
df_anopheles_search <- entrez_search(db = "nuccore", term = "anopheles [ORGN] AND COI [gene]", use_history = T)
df_anopheles_search

# Read data into R
st_anopheles <- readDNAStringSet("../data/sequence.fasta")
st_anopheles

# Create a dataframe
df_anopheles <- data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
head(df_anopheles)
#### TODO: Explore dataframe more thoroughly, make histogram as well

#### GET DINUCLEOTIDE FREQUENCIES####
# Make column nucleotides2
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2)
head(df_anopheles)

# Calculate dinucleotide frequencies and add them to the dataframe as columns
df_anopheles <- cbind(
  df_anopheles,
  as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE))
)
head(df_anopheles)

# Add species name with second and third word from title, and unique_id with first word from title
df_anopheles <- df_anopheles %>%
  mutate(species_name = word(title, 2L, 3L)) %>%
  mutate(unique_id = word(title, 1L))
df_anopheles

#### PREPARE DATA FOR DATASET CREATION####
# Remove column nucleotides2 and make the dataframe into a tibble
df_anopheles <- df_anopheles %>%
  select(-nucleotides2) %>%
  as_tibble()

# Determine number of unique ids
df_anopheles$unique_id %>%
  unique() %>%
  length()

# Create an object containing species names, organized by the number of samples in descending order, revoming all those less than 200 samples
ls_vectors <- df_anopheles %>%
  group_by(species_name) %>%
  count() %>%
  arrange(desc(n)) %>% # View()
  filter(n > 200) %>%
  pull(species_name)
ls_vectors

# Filter by species names from ls_vectors ONLY
df_anopheles_vector <-
  df_anopheles %>%
  filter(species_name %in% ls_vectors)

# Remove uneccessary objects in the environment
rm(df_anopheles, st_anopheles)

# Count number of each species
df_anopheles_vector %>%
  group_by(species_name) %>%
  count()

# Find minimum sample size for later creation of validation data
smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample

#### CREATE VALIDATION DATASET####
# Validation dataset created with each species only having 44 records sampled
set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation

# Check to make sure each species has 44 records sampled
df_anopheles_validation %>%
  group_by(species_name) %>%
  count()
#### TODO: Make test set named df_anopheles_test, assign it a portion of the species names

#### CREATE TRAINING DATASET####
# Training dataset created excluding all samples that appear in the validation dataset, with each species only having 176 records sampled
set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

# CHeck to make sure they all have 176 records sampled
df_anopheles_training %>%
  group_by(species_name) %>%
  count()

# Check names in training dataset
names(df_anopheles_training)

#### CREATE RANDOM FOREST MODEL####
# Use columns including the dinucleotide frequencies, with species name as a factor
anopheles_class <- randomForest(
  x = df_anopheles_training[, 3:18],
  y = as.factor(df_anopheles_training$species_name),
  ntree = 100, importance = TRUE
)

# Check output of random forest by visualizing confusion matrix, plotting the error as more trees are explored, and plotting the mean decrease accuracy to determine the importance of each variable, and the mean decrease Gini to determine how useful each variable is to creating clean splits.
anopheles_class
plot(anopheles_class)
varImpPlot(anopheles_class)
#### TODO: Use previously created test set to test output of the random forest

#### VALIDATE PREDICTION####
# Run prediction
predict_validation <- predict(
  anopheles_class,
  df_anopheles_validation[, c(19, 3:18)]
)
predict_validation

# Check output in a table
table(
  observed = df_anopheles_validation$species_name,
  predicted = predict_validation
)
#### TODO: Visualize the accuracy of the model, represent graphically and write output to a jpeg in the fig folder of the file tree for this R script
