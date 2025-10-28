## _ Packages used -------
library(tidyverse)
library(viridis)
# + scale_color/fill_viridis(discrete = T/F)
theme_set(theme_light())
library(randomForest)
library(Biostrings)
library(rentrez)
library(seqinr)
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflicts_prefer(dplyr::count)
conflicted::conflicts_prefer(base::as.factor)


# Get data from the NCBI nuccore database online using the same search as the function below. Verify data has same number of entires as the count of the code below, download as "Complete Record" and "FASTA" file, place in data folder in same parent directory as the working directory's parent.
df_anopheles_search = entrez_search(db = "nuccore", term = "anopheles [ORGN] AND COI [gene]", use_history = T)
df_anopheles_search


#Loading the Anopheles fasta data from above into environment
st_anopheles = readDNAStringSet("../data/sequence.fasta")
st_anopheles


#Storing the Anopheles data as a dataframe.
df_anopheles = data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
head(df_anopheles)


#Extracting the species name and unique sample ID from the dataframe, then adding it as unique column for ease of reference.
df_anopheles = df_anopheles %>% 
  mutate(species_name = word(title, 2L, 3L)) %>% 
  mutate(unique_id = word(title, 1L)) 
df_anopheles


#Validation, confirming that number of unique ID's equals dataset entries.
df_anopheles$unique_id %>% 
  unique() %>% 
  length()


#Create a dedicated column where the original sequences are converted from characters to Biostrings data for compatiblity and analysis.
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2)
head(df_anopheles)

#Create columns of dinucleotide frequency for each sequence.
df_anopheles = cbind(df_anopheles,
                     as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE)))
head(df_anopheles)

#Converts the rest of the Anopheles dataframe (except nucleotides2) to a tibble for compatibility.
df_anopheles = df_anopheles %>% 
  select(-nucleotides2) %>% 
  as_tibble()


#Create a set of Anopheles species with more than 200 entries in the original database, allowing for a robust testing poplation for statistical analysis.
ls_vectors = df_anopheles %>% 
  group_by(species_name) %>% 
  count() %>% 
  arrange(desc(n)) %>% #View()
  filter(n > 200) %>% 
  pull(species_name)
ls_vectors

#Create a dedicated subset of the original dataframe matching species within the above list.
df_anopheles_vector = 
  df_anopheles %>% 
  filter(species_name %in% ls_vectors)

#Sort the above and check if it worked.
df_anopheles_vector %>% 
  group_by(species_name) %>% 
  count()

#Remove unused variables once data processing using them is done for memory efficiency and general tidiness.
rm(df_anopheles, st_anopheles)

#Find the smallest species population number in the above list, use that number for test pools to avoid oversampling distortion of the model later.
smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample

#Set seed for reproduciblity, then pull a random selection of entries for each species amounting to 20% of the smallest species entry population size, preserving them for model validation.
set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation

#Checking that the above produced a consistent number of entries for each species.
df_anopheles_validation %>% 
  group_by(species_name) %>% 
  count()

#Set seed for reproducibility, then pull a random selection of entries amounting to 80% of the smallest species entry population size, to be used for model training.
set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

#Checking that the above produced a consistent number of entries for each species.
df_anopheles_training %>% 
  group_by(species_name) %>% 
  count()

#Confirm that our training data kept the dinucleotide columns for model training, and which columns correlate to them.
names(df_anopheles_training)


#Train a random forest model from our test data, predicting species name from dinucleotide frequency (dropping the columns containing redundant information).
anopheles_class <- randomForest(x = df_anopheles_training[, 3:18],
                                y = as.factor(df_anopheles_training$species_name), 
                                ntree = 100, importance = TRUE)

#Check if the model worked on the training data, plot model architecture and variable use.
anopheles_class
plot(anopheles_class)
varImpPlot(anopheles_class)

#Apply model to the validation model data, attempting to predict species from dinucleotide data.
predict_validation <- predict(anopheles_class, 
                              df_anopheles_validation[, c(19, 3:18)])
predict_validation


#Check if the model accurately predicted the species.
table(observed = df_anopheles_validation$species_name, 
      predicted = predict_validation)


#Code improvement thoughts 1: Dinucleotide comparison is created for the full dataset, but is only actually used for the high-population species subset. It could save marginal processing time and memory storage to extract that information only from the population subset to-be used, which may become non-incidental for larger datasets.

#Code improvement thoughts 2: No outlier processing is performed on this data for sequences which are abnormally short/long. This is a possible source of error if the model is applied to other datasets, which could be resolved by population quartile filtering.
