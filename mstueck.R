#Call necessary libraries
library(Biostrings)
library(tidyverse)
library(rentrez)
library(dplyr)
library(randomForest)

### Gather and Clean Data ------------------------------------------------------

#Use rentrez nuccore to search for COI gene sequences in Anopheles species
#Searches GenBank's nucleotide database (nuccore)
#Stores search history for later data retrieval
df_anopheles_search = entrez_search(db = "nuccore", term = "anopheles [ORGN] AND COI [gene]", use_history = T)
df_anopheles_search

#Read in sequence data from local FASTA file (downloaded from NCBI)
st_anopheles = readDNAStringSet("../data/sequence.fasta")
st_anopheles 

#Convert FASTA object to dataframe
df_anopheles = data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
head(df_anopheles)

#Create DNAStringSet column for sequence manipulation
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2)
head(df_anopheles)

#Get dinucleotide frequencies from nucleotide sequences
df_anopheles = cbind(df_anopheles,
                     as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE)))

#Clean up unnecessary columns and convert to tibble for readability
df_anopheles = df_anopheles %>% 
  select(-nucleotides2) %>% 
  as_tibble()

head(df_anopheles)

#Separate busy columns - unpack title column into separate variables 
#Extract the species name and unique ID into new columns
df_anopheles = df_anopheles %>% 
  mutate(species_name = word(title, 2L, 3L)) %>% 
  mutate(unique_id = word(title, 1L)) 
df_anopheles

#Quality check: count number of unique sequence IDs
df_anopheles$unique_id %>% 
  unique() %>% 
  length()


### Prepare data for ML --------------------------------------------------------

#Create vector of species with at least 200 sequences
ls_vectors = df_anopheles %>% 
  group_by(species_name) %>% 
  count() %>% 
  arrange(desc(n)) %>% #View()
  filter(n > 200) %>% 
  pull(species_name)
ls_vectors

#Subset dataset to include only selected species
df_anopheles_vector = 
  df_anopheles %>% 
  filter(species_name %in% ls_vectors)

#Check how many samples per species
df_anopheles_vector %>% 
  group_by(species_name) %>% 
  count()

#Identify smallest sample size across selected species (to balance training and validation sets)
smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample

### Split Data into training and validation sets -------------------------------

#Create validation set (20% of each species’ sequences)
set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation

#Check species distribution in validation set
df_anopheles_validation %>% 
  group_by(species_name) %>% 
  count()


#Create training set (remaining 80% of samples)
#Ensures no overlap with validation IDs
set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

#Check class balance in training set
df_anopheles_training %>% 
  group_by(species_name) %>% 
  count()


names(df_anopheles_training)


### Train Random Forest --------------------------------------------------------

#Train Random Forest classifier
#x: dinucleotide frequency features (columns 3:18)
#y: target species labels
#ntree: number of trees to grow
#importance: enables variable importance calculation

anopheles_class <- randomForest(x = df_anopheles_training[, 3:18],
                                y = as.factor(df_anopheles_training$species_name), 
                                ntree = 100, importance = TRUE)


#Visualize model results: Confusion matrix, error rate per tree, and feature importance plots
anopheles_class
plot(anopheles_class)
varImpPlot(anopheles_class)

### Evaluate model success -----------------------------------------------------

#Predict species on validation set
predict_validation <- predict(anopheles_class, 
                              df_anopheles_validation[, c(19, 3:18)])
predict_validation

#Compare observed vs predicted species (confusion matrix)
table(observed = df_anopheles_validation$species_name, 
      predicted = predict_validation)
