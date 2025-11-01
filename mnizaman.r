##### BINF 6210 Software Tools  ----
##
##
##
## Mahnoor Nizamani
##
## 2025-10-23
##
##### ---

## Load Packages -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
conflicted::conflicts_prefer(dplyr::count())
library(viridis)
library(Biostrings)
library(rentrez)
library(seqinr)
library(randomForest)


# Get data -------


# Search the NCBI Nucleotide database for *Anopheles* COI gene sequences.
df_anopheles_search <- entrez_search(db = "nuccore", term = "anopheles [ORGN] AND COI [gene]", use_history = T)
df_anopheles_search

# Read in the downloaded fasta file containing Anopheles sequences.
st_anopheles <- readDNAStringSet("../data/sequence.fasta")
st_anopheles

# Prepare data for analysis -------

# Create a dataframe containing the sequence titles and raw DNA sequence data.
df_anopheles <- data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
head(df_anopheles)

# Convert nucleotides to a DNAStringSet for further analysis with Biostrings package.
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2)
head(df_anopheles)

# This calculates dinucleotide frequencies for each sequence.
df_anopheles <- cbind(df_anopheles,
  as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE)))
head(df_anopheles)

# This removes the Biostring column nucleotides2 and create a tibble.
df_anopheles <- df_anopheles %>%
  select(-nucleotides2) %>%
  as_tibble()

# Extract species names and unique Ids from the fasta titles, saving it to a new column.
df_anopheles <- df_anopheles %>%
  mutate(species_name = word(title, 2L, 3L)) %>%
  mutate(unique_id = word(title, 1L))
df_anopheles

# check what length of unique sequences there are in the data frame (checks how many there are)
df_anopheles$unique_id %>%
  unique() %>%
  length()

# Filter data set -------

# The next block of code filters the data set and only keeps the species that are considered well represented in their sequences, in this data set.

# this pipe only keeps species that have over 200 sequences present in the data set and extracts the filtered species name column as a vector.
ls_vectors <- df_anopheles %>%
  group_by(species_name) %>%
  count() %>%
  arrange(desc(n)) %>%
  filter(n > 200) %>%
  pull(species_name)
ls_vectors

# This takes the full data set and keeps only those rows whose species name matches with one in ls_vectors.
df_anopheles_vector <-
  df_anopheles %>%
  filter(species_name %in% ls_vectors)

# verify number of sequences per species
df_anopheles_vector %>%
  group_by(species_name) %>%
  count()

# remove unneeded
rm(df_anopheles, st_anopheles)

# Split data into training and validation sets -------

# determines smallest sample size among filtered species.
smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample

# this creates a validation set with 20% of the smallest sample size per species, by grouping them and taking the bottom 20%.
# set seed ensures reproducibility
set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation


# Create training set with the rest of the data
set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

# ensure the sample sizes in the set
df_anopheles_training %>%
  group_by(species_name) %>%
  count()

names(df_anopheles_training)

# Train Random Forest -------

# Train random forest model using dinucleotide frequencies from columns 3 through 18 as predictor variables and species names as respose labels/identifications.
anopheles_class <- randomForest(
  x = df_anopheles_training[, 3:18],
  y = as.factor(df_anopheles_training$species_name),
  ntree = 100, importance = TRUE)

# model summary
anopheles_class

# plots model error rate over trees
plot(anopheles_class)

# visualizes which dinucleotides are best at differentiating between species.
varImpPlot(anopheles_class)

# Evaluate model performance by first using the trained model to prefict species identities in the validation set
# Evaluate Performance -------

predict_validation <- predict(anopheles_class,
  df_anopheles_validation[, c(19, 3:18)])
predict_validation

# Then create a table comparing observed vs. predicted labels/identifications of species names.
table(observed = df_anopheles_validation$species_name,
  predicted = predict_validation)
