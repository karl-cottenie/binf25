#### - Load libraries----

library(rentrez)
library(Biostrings)
library(stringr)
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(randomForest)

#### - Load Data and Create Data Frame----

#searcheing the nucleotide NCBI database for the COI gene found in Anopheles
df_anopheles_search = entrez_search(db = "nuccore", term = "anopheles [ORGN] AND COI [gene]", use_history = T)
df_anopheles_search

#loading sequences in a fasta file
st_anopheles = readDNAStringSet("../data/sequence.fasta")
st_anopheles 

#creating a data frame with "title" and "sequence" columns
df_anopheles = data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
head(df_anopheles)

#### - Add Additional Columns to the Data Frame----

#creating a column titled "nucleotides2" containing sequences in a Biostrings format
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2)
head(df_anopheles)

#calculating the proportions of each dinucleotide and adding dinucleotide frequencies as new columns into the data frame
df_anopheles = cbind(df_anopheles,
                     as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE)))
head(df_anopheles)

#deleting the column "nucleotides2" and converting to a tibble for tidyverse
df_anopheles = df_anopheles %>% 
  select(-nucleotides2) %>% 
  as_tibble()

#creating a "species_name" column by taking the 2nd and 3rd word from the title column and creating a "unique_id" column by taking the first word from the title column
df_anopheles = df_anopheles %>% 
  mutate(species_name = word(title, 2L, 3L)) %>% 
  mutate(unique_id = word(title, 1L)) 
df_anopheles

#counting the number of unique id's
df_anopheles$unique_id %>% 
  unique() %>% 
  length()

#### - Manipulate Data Frame for Random Forest Model----

#creating a list that identifies species with more than 200 COI sequences in the database, and sorting from most to least abundant
ls_vectors = df_anopheles %>% 
  group_by(species_name) %>% 
  dplyr::count() %>% 
  arrange(desc(n)) %>% #View()
  filter(n > 200) %>% 
  pull(species_name)
ls_vectors

#creating new dataframe that filters for abundant species only (those that have at least 200 sequences)
df_anopheles_vector = 
  df_anopheles %>% 
  filter(species_name %in% ls_vectors)

#checking counts of how many sequences each species has in the new filtered data frame
df_anopheles_vector %>% 
  group_by(species_name) %>% 
  dplyr::count()

#removes the variables "df_anopheles" and "st_anopheles"
rm(df_anopheles, st_anopheles)

#### - Create Validation and Training Groups for Random Forest Model----

#identifying the smallest group size of species (species with the fewest sequences) and storing that number (220)
smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample

#creating a validation group by randomly selecting 20% of the smallest group size (20% of 220)
set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation

#checking validation set to make sure each species contributes the same number of sequences (44), which is 20% of the smaller_sample
df_anopheles_validation %>% 
  group_by(species_name) %>% 
  dplyr::count()

#creating a training group that uses species names NOT found in the validation group and randomly selects 80% of the smallest group size (80% of 220)
set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

#checking training set to make sure each species contributes the same number of sequences (176), which is 80% of the smaller_sample
df_anopheles_training %>% 
  group_by(species_name) %>% 
  dplyr::count()

#viewing names of all the columns in the training set
names(df_anopheles_training)

#### - Build and Run Random Forest Model----

#training the random forest classifier by inputting lines 3-18 (dinucleotide frequencies) to learn to predict which anopheles species a sequence belongs to based on the dinucleotide composition. The model uses 100 decision trees
anopheles_class <- randomForest(x = df_anopheles_training[, 3:18],
                                y = base::as.factor(df_anopheles_training$species_name), 
                                ntree = 100, importance = TRUE)

#### - View Results of Random Forest Model----

#viewing model summary (with a confusion matrix and error summary) and overall model performance on the training set, plotting to visualizes how the error changes with the number of trees and plotting to see which dinucleotides are the most accurate for species classification
anopheles_class
plot(anopheles_class)
varImpPlot(anopheles_class)

#creating a predictor variable based on the validation group to validate model - outputs the predicted species for each sequence
predict_validation <- predict(anopheles_class, 
                              df_anopheles_validation[, c(19, 3:18)])
predict_validation

#### - Compare Observed Vs. Predicted Results from Random Forest Model----

#generating a confusion matrix to view the predicted species vs the true (observed) species to evaluate model performance and accuracy
table(observed = df_anopheles_validation$species_name, 
      predicted = predict_validation)










