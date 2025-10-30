##load packages
library(rentrez)
library(Biostrings)
library(tidyverse)
library(stringr)
library(randomForest)

##Load Data and Create Data Frame -----
##search for records in GenBank using search terms anopheles for organism and COI for marker
df_anopheles_search = entrez_search(db = "nuccore", 
                                    term = "anopheles [ORGN] AND COI [gene]",
                                    use_history = T)
df_anopheles_search

##read in downloaded fasta file 
st_anopheles = readDNAStringSet("../data/sequence.fasta")
st_anopheles 

##create dataframe using GenBank string label and the sequence
df_anopheles = data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
head(df_anopheles)

## Add Dinucleotide Frequencies, Species Name and ID -----
##create a variable (nucleotides2) that is a DNA string set of the sequence column
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2) ##"Biostrings"
head(df_anopheles)

##calculate nucleotide frequencies and add these columns to the data frame
df_anopheles = cbind(df_anopheles,
                     as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE)))
head(df_anopheles)

##transform the dataframe into a tibble and remove the nucleotides2 column
df_anopheles = df_anopheles %>% 
  select(-nucleotides2) %>% 
  as_tibble()

##create 2 new columns (species name & unique_id) using the title column where binomial species name is generated using the 2nd term (genus) and 3rd term (species) and the unique_id is generated using the 1st term (identifier) 
df_anopheles = df_anopheles %>% 
  mutate(species_name = word(title, 2L, 3L)) %>% 
  mutate(unique_id = word(title, 1L)) 
df_anopheles

##retrieve the total count (length) of unique ids
df_anopheles$unique_id %>% 
  unique() %>% 
  length()

##Prepare data for Random Forest Classifier Model
##create a list (ls_vectors) of species name that are arranged in descending order by count where only counts greater than 200 are included
ls_vectors = df_anopheles %>% 
  group_by(species_name) %>% 
  count() %>% 
  arrange(desc(n)) %>% #View()
  filter(n > 200) %>% 
  pull(species_name)
ls_vectors

##create a variable (df_anopheles_vector) that filters the species names in df_anopheles through the ls_vectors list
df_anopheles_vector = 
  df_anopheles %>% 
  filter(species_name %in% ls_vectors)

##retrieve the count of each species name
df_anopheles_vector %>% 
  group_by(species_name) %>% 
  count()

##remove variables that will not be used anymore
rm(df_anopheles, st_anopheles)

##create a variable (smaller_subsample) to retrieve the lowest species name count
smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample

## Create Training and Validation Data Sets for Classifier -----
##set the seed (for randomness), then create the validation test set using sample size as 20% of the lowest species count (use floor to retrieve nearest integer below value)
set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation

##retrieve the count of each species name in the validation test set; should all be the same as each other and equivalent to 20% of smaller_sample
df_anopheles_validation %>% 
  group_by(species_name) %>% 
  count()

##set seed (for randomness), the create the training test set using 80% of the lowest species count (use ceiling to retrieve nearest integer above value)
set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

##retrieve the count of each species name in the training set; should all be the same as each other and equivalent to 80% of the smaller_sample
df_anopheles_training %>% 
  group_by(species_name) %>% 
  count()

##see column names
names(df_anopheles_training)

## Train and Validate Random Forest Model Classifer, and Visualize -----
#train classifier using random forest method with training data set using nucleotide frequencies as predictor variables and species name as response variable generating 100 trees
anopheles_class <- randomForest(x = df_anopheles_training[, 3:18],
                                y = as.factor(df_anopheles_training$species_name), 
                                ntree = 100, importance = TRUE)

##plot the tree results again classification error rate, plot the model accuracy as predictor variable permutes for each dinucleotide, and plot the predictor variable contribution to node purity for each dinucleotide
anopheles_class
plot(anopheles_class)
varImpPlot(anopheles_class)

##validate the classifier using the validation test data using classification error and the nucleotide frequences
predict_validation <- predict(anopheles_class, 
                              df_anopheles_validation[, c(19, 3:18)])
predict_validation

## Generate Confusion Matrix to Compare Observed vs Predicted Results -----
##display classification prediction results in a confusion matrix
table(observed = df_anopheles_validation$species_name, 
      predicted = predict_validation)