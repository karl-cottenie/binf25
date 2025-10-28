##### BINF6120 2025-10-23 =======
## anopheles script puzzle

# === Packages used =======
library(rentrez)
library(Biostrings)
library(tidyverse)
library(randomForest)
# Startup ends here

# Searching Nucleotide database for Anopheles records with COI gene
df_anopheles_search = entrez_search(db = "nuccore", term = "anopheles [ORGN] AND COI [gene]", use_history = T)
df_anopheles_search

# Reading in sequence data
st_anopheles = readDNAStringSet("../data/sequence.fasta")
st_anopheles

# creating dataframe with title and sequence columns
df_anopheles = data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
head(df_anopheles)

# Creating nucleotides column from sequence data
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2)
head(df_anopheles)

# Calculating dinucleotide frequency and joining to dataframe
df_anopheles = cbind(df_anopheles,
                     as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE)))
head(df_anopheles)

# Removing nucleotides column and setting every remaining column as tibble
df_anopheles = df_anopheles %>% 
  select(-nucleotides2) %>% 
  as_tibble()

# Creating species name and unique id columns 
df_anopheles = df_anopheles %>% 
  mutate(species_name = word(title, 2L, 3L)) %>% 
  mutate(unique_id = word(title, 1L)) 
df_anopheles

# Validating number of records (n = 10730)
df_anopheles$unique_id %>% 
  unique() %>% 
  length()

# Creating list of species with more than 200 records
ls_vectors = df_anopheles %>% 
  group_by(species_name) %>% 
  count() %>% 
  arrange(desc(n)) %>% #View()
  filter(n > 200) %>% 
  pull(species_name)
ls_vectors

# Create vector only including species in ls_vectors
df_anopheles_vector = 
  df_anopheles %>% 
  filter(species_name %in% ls_vectors)

# Validating counts of species (all should be more than 200)
df_anopheles_vector %>% 
  group_by(species_name) %>% 
  count()

rm(df_anopheles, st_anopheles)

# Taking the smallest number of records to ensure no bias
smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample

# Creating validation set using 20% of minimum sample size
set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation

# Checking sample size for each species in validation set (n = 44)
df_anopheles_validation %>% 
  group_by(species_name) %>% 
  count()

# Creating training set using 80% of minimum sample size and excluding samples in validation set
set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

# Validating sample size for each species in training set (n = 176)
df_anopheles_training %>% 
  group_by(species_name) %>% 
  count()

# Checking which columns to use in classifier
names(df_anopheles_training)

# Creating classifier with dinucleotide frequencies as predictor and species name as response
anopheles_class <- randomForest(x = df_anopheles_training[, 3:18],
                                y = as.factor(df_anopheles_training$species_name), 
                                ntree = 100, importance = TRUE)

# Checking results of classifier after training (error, accuracy, Gini)
anopheles_class
plot(anopheles_class)
varImpPlot(anopheles_class)

# validating classifier with validation set
predict_validation <- predict(anopheles_class, 
                              df_anopheles_validation[, c(19, 3:18)])
predict_validation

table(observed = df_anopheles_validation$species_name, 
      predicted = predict_validation)