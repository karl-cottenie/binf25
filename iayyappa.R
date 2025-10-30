## Anopheles COI randomForest

## _ Packages used -------
library(rentrez)
library(tidyverse)
library(randomForest)
library(Biostrings)
library(conflicted)
conflicted::conflicts_prefer(base::as.factor)

#### 1. Finding COI records in NCBI ----
# This shows the search pattern to find relevant COI records.
df_anopheles_search <- entrez_search(db = "nuccore", term = "anopheles [ORGN] AND COI [gene]", use_history = T)
df_anopheles_search

#### 2. Read data and build dataframe ----
# Load data from the folder
st_anopheles <- readDNAStringSet("../data/sequence.fasta")
st_anopheles

# Create a simple data frame with sequence headers and sequence strings.
df_anopheles <- data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
head(df_anopheles)

#### 3. DNAStringSet Column ----
# Keep a DNAStringSet for computing
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2)
head(df_anopheles)

# Compute dinucleotide frequencies and bind to df
df_anopheles <- cbind(
  df_anopheles,
  as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE))
)
head(df_anopheles)

# Drop nucleotides column and convert to tibble
df_anopheles <- df_anopheles %>%
  select(-nucleotides2) %>%
  as_tibble()

#### 4. Cleaning and parsing ----
# Extract only required headers (species name and unique id)
df_anopheles <- df_anopheles %>%
  mutate(species_name = word(title, 2L, 3L)) %>%
  mutate(unique_id = word(title, 1L))
df_anopheles

# Check unique records
df_anopheles$unique_id %>%
  unique() %>%
  length()

# Identify well-sampled species
ls_vectors <- df_anopheles %>%
  group_by(species_name) %>%
  count() %>%
  arrange(desc(n)) %>% # View()
  filter(n > 200) %>%
  pull(species_name)
ls_vectors

# Store only the records for those well-sampled species
df_anopheles_vector <-
  df_anopheles %>%
  filter(species_name %in% ls_vectors)

# Class size overview after filtering
df_anopheles_vector %>%
  group_by(species_name) %>%
  count()

# Cleanup
rm(df_anopheles, st_anopheles)

#### 5. Smallest class size ----
smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample

#### 6. Validation data set — 20% of class ----
set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation

# Check counts in validation
df_anopheles_validation %>%
  group_by(species_name) %>%
  count()

#### 7. Make training set — 80% of class ----
set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

# Verifying training counts
df_anopheles_training %>%
  group_by(species_name) %>%
  count()

# Verifying column names
names(df_anopheles_training)

#### 8. Random Forest classifier ----
anopheles_class <- randomForest(
  x = df_anopheles_training[, 3:18],
  y = as.factor(df_anopheles_training$species_name),
  ntree = 100, importance = TRUE
)

#### 9. Analysis and plots ----
anopheles_class
plot(anopheles_class)
varImpPlot(anopheles_class)


#### 10. Prediction on validation set ----
predict_validation <- predict(
  anopheles_class,
  df_anopheles_validation[, c(19, 3:18)]
)
predict_validation

# Final evaluation (observed vs predicted count)
table(
  observed = df_anopheles_validation$species_name,
  predicted = predict_validation
)