# ---- Dependencies --------------------
rm(list=ls())
library(rentrez)
library (Biostrings) 
library(stringr)
library(tidyverse)
library(randomForest)

# ---- Data Import and Formatting --------------------

#search NCBI (nuccore) for COI gene sequences from Anopheles species.
df_anopheles_search = entrez_search(db = "nuccore", term = "anopheles [ORGN] AND COI [gene]", use_history = T)
df_anopheles_search

#reads DNA sequences from a local FASTA file.
st_anopheles = readDNAStringSet("../data/sequence.fasta")
st_anopheles 

#converts to data frame: the title is the seq headers from FASTA file, sequence is sequence (as a string) 
df_anopheles = data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
head(df_anopheles)

# TODO: Explore this dataframe more thoroughly — check sequence lengths, base composition, and make a histogram of sequence lengths or GC content to understand data quality.

# ---- Measure Sequence Patterns --------------------

#re-formatted string sequences into DNA biostrings and calculates dinucleotide frequencies (e.g., AA, AT, etc.), stored as proportions (is this adding a column?) 
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2)
head(df_anopheles)

df_anopheles = cbind(df_anopheles,
                     as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE)))
head(df_anopheles)

#---- Add Species Information  --------------------

#gather info from seq titles. 
df_anopheles = df_anopheles %>% 
  mutate(species_name = word(title, 2L, 3L)) %>% ## take 2nd–3rd words in header (Genus species) 
  mutate(unique_id = word(title, 1L)) #take 1st token in header as a simple record id
df_anopheles #show the updated table

#Counts the number of unique sequences.
df_anopheles$unique_id %>% 
  unique() %>% 
  length()

#Removes the Biostrings column.
df_anopheles = df_anopheles %>% 
  dplyr::select(-nucleotides2) %>% 
  as_tibble() #for dplyr compatibility

#Keeps only species with at least 200 sequences, pulls them into a list (ls_vectors) for filtering.
ls_vectors = df_anopheles %>% 
  group_by(species_name) %>% 
  count() %>%  # per-species counts
  arrange(desc(n)) %>% #View() # arranges by largest first
  filter(n > 200) %>% # keep species having > 200 sequences
  pull(species_name) # return as a character vector
ls_vectors # print the retained species names

#Filters the full dataset to only sequences from those top vector species.
df_anopheles_vector = 
  df_anopheles %>% 
  filter(species_name %in% ls_vectors)

df_anopheles_vector %>% 
  group_by(species_name) %>% 
  count() # check counts after filtering

rm(df_anopheles, st_anopheles) #removed bcs they wont be used again

#---- Split the Data for Training and Testing --------------------
##Gets the size of the smallest group (species) to balance training/testing datasets.
smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample # size of the rarest class among retained species

#Randomly samples 20% of the smallest group size from each species. These samples are held out entirely for model validation.
set.seed(51) # reproducible sampling
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample)) # take the same number from each class
df_anopheles_validation

#check that the validation dataset has the same number of seq per species, which confirms that your earlier sample step worked properly and that the validation set is evenly balanced (ie has the same count for each species (should be 44)).
df_anopheles_validation %>% 
  group_by(species_name) %>% 
  count() # verify equal validation rows per species

# TODO: Create a separate *test set* from df_anopheles_vector, assigning some species names to it for final evaluation after training and validation

set.seed(40) # set a random seed
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>% # removes any sequences that were used for validation
  group_by(species_name) %>% # groups remaining data by species
  sample_n(ceiling(0.8 * smaller_sample)) # randomly picks 80% of the smallest group size for training, keeping the same number per species ()

df_anopheles_training %>% 
  group_by(species_name) %>% 
  count() #should be the same count for every species all the way through (should be 176) 
names(df_anopheles_training)

# ---- Train and Test the Model --------------------

#trains the model to tell different Anopheles species apart based on their DNA pattern frequencies
anopheles_class <- randomForest(x = df_anopheles_training[, 3:18], # trains a random forest model and saves it as 'anopheles_class'. uses columns 3 to 18 (the dinucleotide values) as the input features. 
                                y = as.factor(df_anopheles_training$species_name),  #uses the species name as the label to predict
                                ntree = 100, importance = TRUE) #builds 100 decision trees in the forest. importance = true tells the model to measure which features are most useful

# check how well the random forest model works (how accurate the model is and which DNA base pair patterns are most useful for telling the species apart)
anopheles_class #numeric summary of model performance
plot(anopheles_class) #visual check of how error changes with more trees
#MeanDecreaseAccuracy: how much accuracy drops if that variable is randomized (higher = more important)
#MeanDecreaseGini: how useful the variable is for making clean splits between species (higher = better)
varImpPlot(anopheles_class) #visual check of which dinucleotide patterns were most informative 

# TODO: Once test set is created, use it here to evaluate how well the trained model performs on completely unseen data 

# use the trained random forest model to make predictions
predict_validation <- predict(anopheles_class, #the model we trained earlier
                              df_anopheles_validation[, c(19, 3:18)]) #the validation data (should include the same feature columns used for training)
predict_validation #shows the predicted species for each validation sequence

#test how well the model works by predicting species on the validation data and comparing those predictions to the actual species names
table(observed = df_anopheles_validation$species_name,  # the real species names (actual labels)
      predicted = predict_validation) # the species names the model guessed (predictions)

# TODO: Visualize the accuracy results — plot a bar chart or heatmap of the confusion matrix






