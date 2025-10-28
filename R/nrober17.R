#### Software Tools - Random Forest Model for Anopheles Species Classification using Dinucleotide frequencies from COI sequences ----
# Nadira Robertson
# 23/10/25

#### Load packages ----

library(rentrez)
library(seqinr)
library(Biostrings)
library(tidyverse)
library(randomForest)

#### Part 1: AQUIRE FILES ----

# Searches NCBI for anopheles COI gene sequences

df_anopheles_search <- entrez_search(db = "nuccore", term = "anopheles [ORGN] AND COI [gene]", use_history = T)
df_anopheles_search

# Loads FASTA File for data

st_anopheles <- readDNAStringSet("../data/sequence.fasta")
st_anopheles

#### Part 2: PREPARE AND FORMAT FILES ----

# Create data frame for anopheles with titles and sequences

df_anopheles <- data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
head(df_anopheles)


# Create nucleotide2 column in DNAStringSet format using sequence column

df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2)
head(df_anopheles)

# Calculate the dinucleotide frequency and combine with anopheles file

df_anopheles <- cbind(
  df_anopheles,
  as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE))
)
head(df_anopheles)

# Set nucleotides column as tibble vector and remove nucleotides2

df_anopheles <- df_anopheles %>%
  select(-nucleotides2) %>%
  as_tibble()

# Create species name column using FASTA file headers and ID column

df_anopheles <- df_anopheles %>%
  mutate(species_name = word(title, 2L, 3L)) %>%
  mutate(unique_id = word(title, 1L))
df_anopheles

# Counts the number of unique IDs

df_anopheles$unique_id %>%
  unique() %>%
  length()

#### Part 3: CREATE TRAINING AND VALIDATION DATA SETS ----

# Collect species with sequence counts over 200 and create vector

ls_vectors <- df_anopheles %>%
  group_by(species_name) %>%
  dplyr::count() %>%
  arrange(desc(n)) %>% # View()
  filter(n > 200) %>%
  pull(species_name)
ls_vectors

# Create data set with only vector species

df_anopheles_vector <-
  df_anopheles %>%
  filter(species_name %in% ls_vectors)

# Lists the number of sequences each species in the vector has

df_anopheles_vector %>%
  group_by(species_name) %>%
  dplyr::count()

# Remove uneeded files

rm(df_anopheles, st_anopheles)

# Creates vector of smallest group value in anopheles vector

smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample

# Set seed and create validation data set with 20% of sequences

set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation

# Double check the count for validation file for each species

df_anopheles_validation %>%
  group_by(species_name) %>%
  dplyr::count()

# Set seed and create training data using vector with 80% of sequences, exclude those used

set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

# Double check the count for training data for each species

df_anopheles_training %>%
  group_by(species_name) %>%
  dplyr::count()

# Check names of columns in training file

names(df_anopheles_training)

#### Part 4: RUNNING RANDOM FOREST MODEL ----

# Use random forest and training set to train model on dinucleotide frequencies

anopheles_class <- randomForest(
  x = df_anopheles_training[, 3:18],
  y = as.factor(df_anopheles_training$species_name),
  ntree = 100, importance = TRUE
)

#### Part 5: PLOT RESULTS AND MODEL PREFORMANCE ----

# View results summary, plot random forest error rate results, and plot dinucleotide features

anopheles_class
plot(anopheles_class)
varImpPlot(anopheles_class)

# Predict species in validation data set

predict_validation <- predict(
  anopheles_class,
  df_anopheles_validation[, c(19, 3:18)]
)
predict_validation

# Compare predicted and observed species labels

table(
  observed = df_anopheles_validation$species_name,
  predicted = predict_validation
)
