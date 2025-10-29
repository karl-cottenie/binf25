#####
## Script: anopheles_RF_dinuc_features.R
## Purpose: Classify Anopheles species from COI sequences using
##          dinucleotide composition and a Random Forest model.
##
## Author: Kexin Gong
#####

#### 1 - QUERY NCBI: Entrez search handle (history) ----
# Pull a search history for Anopheles COI from the 'nuccore' db.
# Note: Requires an NCBI key (rentrez::set_entrez_key()) if rate-limited.
df_anopheles_search = entrez_search(db = "nuccore", 
                                    term = "anopheles [ORGN] AND COI [gene]",
                                    use_history = T)
df_anopheles_search

#### 2 - LOAD SEQUENCES: FASTA -> DNAStringSet ----
# Read a local FASTA of sequences selected for this analysis.
# Adjust the relative path to your environment.
st_anopheles = readDNAStringSet("../data/sequence.fasta")
st_anopheles 

#### 3 - BUILD BASE DATA FRAME ----
# Convert DNAStringSet to a data.frame with 'title' and string 'sequence'.
df_anopheles = data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
head(df_anopheles)

#### 3.1 - CAST TO DNAStringSet & CALCULATE DINUCLEOTIDES ----
# Keep a DNAStringSet column for feature extraction.
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2)
head(df_anopheles)

# Compute dinucleotide frequencies (as probabilities); yields 16 features.
df_anopheles = cbind(df_anopheles,
                     as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE)))
head(df_anopheles)

#### 3.2 - CLEAN TEMP DNA COLUMN; COERCE TO TIBBLE ----
# Drop helper DNAStringSet column and standardize to tibble for dplyr ops.
df_anopheles = df_anopheles %>% 
  select(-nucleotides2) %>% 
  as_tibble()

#### 4 - ADD METADATA FROM TITLES ----
# Extract species binomial (2nd–3rd tokens) and a unique id (1st token).
df_anopheles = df_anopheles %>% 
  mutate(species_name = word(title, 2L, 3L)) %>% 
  mutate(unique_id = word(title, 1L)) 
df_anopheles

# Quick sanity check: count of unique accessions / records.
df_anopheles$unique_id %>% 
  unique() %>% 
  length()

#### 5 - SELECT VECTORS (ABUNDANT SPECIES) ----
# Keep species with >200 sequences to balance classes for classification.
ls_vectors = df_anopheles %>% 
  group_by(species_name) %>% 
  count() %>% 
  arrange(desc(n)) %>% #View()
  filter(n > 200) %>% 
  pull(species_name)
ls_vectors

# Subset to those species only.
df_anopheles_vector = 
  df_anopheles %>% 
  filter(species_name %in% ls_vectors)

# Class counts per species after filtering.
df_anopheles_vector %>% 
  group_by(species_name) %>% 
  count()

# (Optional) free memory from very large objects.
rm(df_anopheles, st_anopheles)

#### 6 - STRATIFIED TRAIN/VALIDATION SPLIT ----
# Use the smallest class size to define per-class sample sizes.
smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample

# Validation: 20% of the smallest class per species (stratified).
set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation

# Check validation class balance.
df_anopheles_validation %>% 
  group_by(species_name) %>% 
  count()

# Training: remaining records; then sample 80% of smallest class per species.
# Note: Filtering by unique_id ensures no leakage if titles contain ids.
set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

# Confirm training class balance.
df_anopheles_training %>% 
  group_by(species_name) %>% 
  count()

# Peek at columns: first 2 are metadata; 3:18 are 16 dinucleotide features.
names(df_anopheles_training)

#### 7 - MODEL: Random Forest training & importance ----
# Fit RF classifier with 100 trees on dinucleotide features.
# x: columns 3–18 (16 features), y: species factor.
anopheles_class <- randomForest(x = df_anopheles_training[, 3:18],
                                y = as.factor(df_anopheles_training$species_name), 
                                ntree = 100, importance = TRUE)

# Model summary and learning curve + variable importance plots.
anopheles_class
plot(anopheles_class)
varImpPlot(anopheles_class)

#### 8 - EVALUATE: Prediction & confusion matrix ----
# Predict on validation set. Ensure feature columns align.
predict_validation <- predict(anopheles_class, 
                              df_anopheles_validation[, c(19, 3:18)])
predict_validation

# Confusion matrix of observed vs. predicted species.
table(observed = df_anopheles_validation$species_name, 
      predicted = predict_validation)
