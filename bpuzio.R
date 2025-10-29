# This R script performs an analysis of mitochondrial COI gene sequences from Anopheles mosquito species. It includes:
#Retrieving COI sequences from the NCBI Nucleotide database.
#Reading and formatting FASTA data into a structured dataframe.
#Calculating dinucleotide frequencies
#Extracting species identifiers and metadata from sequence titles.
#Selecting species with sufficient representation
#Splitting the dataset into balanced training and validation sets (80/20).
#Training a Random Forest classifier to predict species identity.
#Evaluating model performance
# OUTPUT:
# Random forest model summary and variable importance plots.
# Confusion matrix comparing predicted vs. observed species.
 
library(rentrez)
library(Biostrings)
library(dplyr)
library(seqinr)
library(randomForest)

#### 1. DOWNLOAD SEQUENCES FROM NCBI ####
# Search the NCBI Nucleotide database for sequences
df_anopheles_search = entrez_search(
  db = "nuccore",
  term = "anopheles [ORGN] AND COI [gene]",
  use_history = TRUE
)
df_anopheles_search  # inspect search summary 

#### 2. READ FASTA SEQUENCES FROM FILE ####
.
st_anopheles = readDNAStringSet("../data/sequence.fasta")
st_anopheles  # shows how many sequences

#### 3. BUILD A DATA FRAME ####

# Convert the DNAStringSet into a regular data frame with:
# title: sequence header from FASTA
# sequence: the actual nucleotide sequence as plain text
df_anopheles = data.frame(
  title = names(st_anopheles),
  sequence = paste(st_anopheles)
)
head(df_anopheles)

# Create a DNAStringSet column from the character sequence. We need this object type to calculate dinucleotide frequencies.
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)

class(df_anopheles$nucleotides2)  # confirm it's DNAStringSet
head(df_anopheles)

#### 4. CALCULATE DINUCLEOTIDE FREQUENCIES ####

# dinucleotideFrequency(as.prob = TRUE) returns the frequencies of all 16 possible dinucleotides  for each sequence.We cbind() those frequency columns onto df_anopheles so that each row now has sequence-level features we can later use as predictors.
df_anopheles = cbind(
  df_anopheles,
  as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE))
)

head(df_anopheles)  # confirm new feature columns are present

#### 5. CLEAN UP STRUCTURE AND ADD METADATA ####

# We no longer need the DNAStringSet column 'nucleotides2' itself
# the numeric features are what we care about for modeling.
# Convert to tibble for nicer dplyr behavior.
df_anopheles = df_anopheles %>% 
  select(-nucleotides2) %>% 
  as_tibble()

# Parse identifiers out of the FASTA title line.
# species_name: 2nd and 3rd words
# unique_id:    1st word
df_anopheles = df_anopheles %>% 
  mutate(species_name = word(title, 2L, 3L)) %>% 
  mutate(unique_id   = word(title, 1L)) 

df_anopheles  # view full tibble with features

# Count how many unique sequence IDs exist.
df_anopheles$unique_id %>% 
  unique() %>% 
  length()

#### 6. FIND WELL-SAMPLED SPECIES ####

# Count how many sequences we have per species_name, arrange from most to least, keep only species with >200 sequences and pull just those species names.
ls_vectors = df_anopheles %>% 
  group_by(species_name) %>% 
  count() %>% 
  arrange(desc(n)) %>%
  filter(n > 200) %>% 
  pull(species_name)

ls_vectors  # list of high-representation species


#### 7. SUBSET TO ONLY THOSE HIGHLY REPRESENTED SPECIES ####

# Filter the main dataset to include only species found in ls_vectors
df_anopheles_vector = df_anopheles %>% 
  filter(species_name %in% ls_vectors)

# Sanity check: how many sequences per species after filtering
df_anopheles_vector %>% 
  group_by(species_name) %>% 
  count()

# Drop intermediate objects we don't need anymore
rm(df_anopheles, st_anopheles)

# Find the smallest class size among the kept species.
smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample

#### 8. CREATE VALIDATION SET (20% PER SPECIES) ####

set.seed(51)
# For each species_name group:
# - sample 20% of the "smaller_sample" size to form the validation set.
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))

df_anopheles_validation

# Count the number of validation sequences per species to verify that each species is equally represented in the validation subset.
df_anopheles_validation %>% 
  group_by(species_name) %>% 
  count()

#### 9. CREATE TRAINING SET (REMAINING 80%) ####

# Create the training dataset by:
#   1. Excluding any sequences already used in the validation set 
#   2. Grouping by species to maintain class balance.
#   3. Sampling approximately 80% of the smallest class size for each species.
# Finally, count how many training sequences belong to each species to confirm that the training data is evenly balanced.
set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

df_anopheles_training %>% 
  group_by(species_name) %>% 
  count()

names(df_anopheles_training)  # inspect column order

#### 10. TRAIN RANDOM FOREST CLASSIFIER ####

# Train a random forest to classify species_name based on the dinucleotide frequency features.

anopheles_class <- randomForest(
  x = df_anopheles_training[, 3:18],
  y = as.factor(df_anopheles_training$species_name), 
  ntree = 100,
  importance = TRUE
)

anopheles_class        # model summary 
plot(anopheles_class)  # error rate vs number of trees
varImpPlot(anopheles_class)  # which dinucleotides matter most

#### 11. VALIDATE MODEL ON DATA ####

# Use the trained model to predict species_name for the validation set.

predict_validation <- predict(
  anopheles_class, 
  df_anopheles_validation[, c(19, 3:18)]
)
predict_validation  # predicted class labels

# observed (true) vs predicted species_name.
table(
  observed  = df_anopheles_validation$species_name, 
  predicted = predict_validation
)