##########################################################
# Title: Anopheles COI Gene Classification using Random Forest
# Author: Vian Lelo
# Description:
#   This script retrieves Anopheles COI sequences from NCBI,
#   processes nucleotide data, extracts dinucleotide frequencies,
#   prepares training/validation datasets,
#   and builds a Random Forest classifier to predict species.
##########################################################

## 1. Load Required Libraries -----

# install.packages(c("tidyverse", "Biostrings", "randomForest", "rentrez", "stringr")) #TODO confirm versions before running

library(tidyverse) # Data wrangling
library(Biostrings) # DNA sequence operations
library(randomForest) # Random Forest model building
library(rentrez) # Retrieve sequence data from NCBI
library(stringr) # Text manipulation


## 2. Retrieve Sequence Data -----
## Goal: Query NCBI for all Anopheles COI gene sequences.
## Justification: COI gene is a common DNA barcode for species classification (#==> ensures comparable data source)

df_anopheles_search <- entrez_search(
  db = "nuccore",
  term = "anopheles [ORGN] AND COI[gene]",
  use_history = T
)

df_anopheles_search # (#==> returns search metadata)


## 3. Load FASTA Sequence Data -----
## Justification: Using readDNAStringSet() ensures consistent parsing of FASTA-formatted nucleotide data.

st_anopheles <- readDNAStringSet("../data/sequence.fasta")
st_anopheles # (#==> object of class DNAStringSet; contains sequence headers and strings)


## 4. Convert Sequences into a Data Frame -----
## Goal: Convert DNA sequences into a tidy format for downstream processing
#- "title": sequence header (metadata line)
#- "sequence": the actual nucleotide sequence
df_anopheles <- data.frame(
  title = names(st_anopheles),
  sequence = paste(st_anopheles)
)

head(df_anopheles) # (#==> preview data: shows sequence IDs and raw nucleotide strings)


## 5. Convert Sequences to DNAStringSet for Analysis -----
## Justification: The DNAStringSet class allows efficient computation of base composition and frequencies.

df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2) # (#==> should return "DNAStringSet")
head(df_anopheles)


## 6. Calculate Dinucleotide Frequencies -----
## Justification: Dinucleotide frequencies (AA, AT, CG, etc.) serve as molecular descriptors for classification.

df_anopheles <- cbind(
  df_anopheles,
  as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE))   # "as.prob = TRUE" gives proportions instead of raw counts.
)

head(df_anopheles) # (#==> new columns correspond to normalized dinucleotide probabilities)


## 7. Clean and Prepare Metadata -----
## Justification: Metadata extraction creates clean identifiers for grouping and modeling.

# Remove nucleotides2 (unnecessary large object) and convert to tibble
df_anopheles <- df_anopheles %>%
  select(-nucleotides2) %>%
  as_tibble()

# Extract useful info from sequence titles:
#   species_name → words 2 and 3 of title (usually "Anopheles gambiae")
#   unique_id → word 1 of title (unique accession or identifier)
df_anopheles <- df_anopheles %>%
  mutate(species_name = word(title, 2L, 3L)) %>%
  mutate(unique_id = word(title, 1L))

df_anopheles # (#==> shows species names and unique identifiers ready for grouping)


## 8. Explore the Dataset -----
## Justification: Determine dataset balance and sequence representation per species.
## Ref: Data inspection step follows tidyverse best practices for exploratory data analysis (EDA).

df_anopheles$unique_id %>%
  unique() %>%
  length() # (#==> number of unique sequences retrieved)

ls_vectors <- df_anopheles %>%
  group_by(species_name) %>%
  count() %>%
  arrange(desc(n)) %>%
  filter(n > 200) %>%
  pull(species_name)

ls_vectors # (#==> list of species with >200 sequences retained for balanced analysis)


## 9. Filter for Common Species -----
## Justification: Maintain species with enough data points for reliable machine learning results.
## Reference: Similar thresholding used in population genomic datasets to avoid sampling bias.

df_anopheles_vector <- df_anopheles %>%
  filter(species_name %in% ls_vectors)

df_anopheles_vector %>%
  group_by(species_name) %>%
  count() # (#==> balanced dataset ready for sampling)


## 10. Clean Up Environment -----
## Justification: Remove large unused objects to optimize memory.
rm(df_anopheles, st_anopheles) # (#==> cleanup complete)


## 11. Prepare Training and Validation Sets -----
## Justification: Split data into train/test for unbiased model evaluation (80/20 split).
## Ref: Kuhn & Johnson, “Applied Predictive Modeling” (Springer, 2013).
## (*_*) Corresponds to “Model training and validation” subsection in manuscript.

smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample # (#==> smallest group size among selected species)

set.seed(51) # (#==> ensures reproducibility)

df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation # (#==> validation subset sampled)

df_anopheles_validation %>%
  group_by(species_name) %>%
  count() # (#==> verify per-species balance in validation set)

set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

df_anopheles_training %>%
  group_by(species_name) %>%
  count() # (#==> confirm balanced training dataset)

names(df_anopheles_training) # (#==> check input feature names)


## 12. Train Random Forest Classifier -----
## Justification: Random Forest is robust to correlated predictors (dinucleotide frequencies).
## Ref: Breiman (2001) “Random Forests,” Machine Learning 45(1):5–32.
## (*_*) Will produce Figure 3: Feature importance visualization.

anopheles_class <- randomForest(
  x = df_anopheles_training[, 3:18],
  y = as.factor(df_anopheles_training$species_name),
  ntree = 100,
  importance = TRUE
)

anopheles_class # (#==> model summary: OOB error rate, confusion matrix)
plot(anopheles_class) # (*_*) Figure 2: Random Forest training curve
varImpPlot(anopheles_class) # (*_*) Figure 3: Variable importance plot


## 13. Validate Model on Test Set -----
## Justification: Assess generalization accuracy using unseen validation data.
## Ref: Validation approach aligns with standard supervised learning frameworks.

predict_validation <- predict(
  anopheles_class,
  df_anopheles_validation[, c(19, 3:18)]
)

predict_validation # (#==> predicted species labels)

table(
  observed = df_anopheles_validation$species_name,
  predicted = predict_validation
) # (#==> confusion matrix; used for Table 2 in manuscript)

## (#==> Model achieved classification accuracy to be reported in “Results” section)
## #TODO compute metrics like precision, recall, F1-score in next iteration

##########################################################
# End of Script
##########################################################
