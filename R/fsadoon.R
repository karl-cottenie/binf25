##***************************
## Software Tools for Bioinformatics - Lecture 14 - In-Class Activity
##
## Farah Sadoon
##
## 2025-10-23
##
##***************************

### Load Libraries
library(tidyverse)
library(rentrez)
library(Biostrings)
library(randomForest)

## Load Data and Create Data Frame ----
# load fasta file
st_anopheles = readDNAStringSet("../data/sequence.fasta")
st_anopheles

# Check to make sure that the search here gives you the same number of hits as what you get 
df_anopheles_search = entrez_search(db = "nuccore", term = "anopheles [ORGN] AND COI [gene]", use_history = T)
df_anopheles_search

# create df_anopheles by selecting names and sequences from fasta file
df_anopheles = data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
head(df_anopheles)

## Add Columns with Dinucleotide Frequencies, Species Names and IDs ----
# create DNA string set from sequences in df_anopheles
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2)
head(df_anopheles)

# determine dinucleotide frequency for the sequences
df_anopheles = cbind(df_anopheles,
                     as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE)))
head(df_anopheles)

# extract names and unique ids from title column and create new columns
df_anopheles = df_anopheles %>% 
  mutate(species_name = word(title, 2L, 3L)) %>% # from the title field select the second and third words (identified by space delimiters) to extract species names
  mutate(unique_id = word(title, 1L)) # use the first word in the title field as the unique_ids for each record
df_anopheles

# check to make sure the number of unique ids matches the number of records in the data frame
df_anopheles$unique_id %>% 
  unique() %>% 
  length()

# make df_anopheles a tibble, and remove the nucleotides2 (as it is a repeat of the sequences, )
df_anopheles = df_anopheles %>% 
  select(-nucleotides2) %>% 
  as_tibble()

## Prepare Data Frame for Random Forest Model ----
# create a list of vectors to extract from df_anopheles
ls_vectors = df_anopheles %>% 
  group_by(species_name) %>% # grouping by species
  count() %>% # count number of records for each species
  arrange(desc(n)) %>% #View() # put the list in order by number of records descending
  filter(n > 200) %>% # filter out the values where there are less than 200 records
  pull(species_name) # create vectors where each species is the key (each vector will contain information for each species)
ls_vectors

# create a data frame from df_anopheles, filter it out to make sure the the species names exist in the list of vectors
df_anopheles_vector = 
  df_anopheles %>% 
  filter(species_name %in% ls_vectors)

# view the number of records for each species in df_anopheles_vector (each record contains sequences for a specific species), how many times is each species represented
df_anopheles_vector %>% 
  group_by(species_name) %>% 
  count()

# remove objects that we don't need anymore
rm(df_anopheles, st_anopheles)

# create a smaller sample of data for creating and testing model - we're limiting the data size to the lowest common denominator of number of species
smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample

## Create Validation and Training Data Sets for Random Forest Model ----
# Create validation set of data to remove from original data and test on it later
set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation

# Look at the validation data and the counts
df_anopheles_validation %>% 
  group_by(species_name) %>% 
  count()

# Create training data
set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

# Look at the training data to make sure the number for each species is the same
df_anopheles_training %>% 
  group_by(species_name) %>% 
  count()

# Look at the names of species in the training dataset 
names(df_anopheles_training)

## Build and Apply Random Forest Model ----
# Create random forest model for the training data using dinucleotide frequencies as the predictors, create 100 decision trees
anopheles_class <- randomForest(x = df_anopheles_training[, 3:18],
                                y = as.factor(df_anopheles_training$species_name), 
                                ntree = 100, importance = TRUE)

# look at the results of the random forest model 
anopheles_class # prints the output of the random forest model
plot(anopheles_class) # plots the error rate of the model as you add more trees
varImpPlot(anopheles_class) # variable importance plot - plots which variables are the most important predictors in the model

# apply the model to the validation data, using dinucelotide frequencies
predict_validation <- predict(anopheles_class, 
                              df_anopheles_validation[, c(19, 3:18)])
predict_validation

## Compare Observed vs. Predicted Results from Original Data and Random Forest Model ----
# Look at and compare the differences between what the model predicted and what was actually observed in the validation data
table(observed = df_anopheles_validation$species_name, 
      predicted = predict_validation)