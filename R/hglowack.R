####
#Hannah Glowacki
#Quiz 2
#October 26, 2025
####

#1. Load libraries----
library("rentrez")
library("tidyverse")
conflicted::conflict_prefer("filter", "dplyr")
library("randomForest")
library("Biostrings")

#2. Search for and download data  ----

#Search for Anopheles and COI gene on NCBI website 
df_anopheles_search = entrez_search(db = "nuccore", 
                                    term = "anopheles [ORGN] AND COI [gene]",
                                    use_history = T)
df_anopheles_search

#Upload to R studio from downloaded Biostring DNA string fasta data file on desktop
st_anopheles = readDNAStringSet("../data/sequence.fasta")
st_anopheles 

#3. Convert and calculate data sequence features----
#Create dataframe from DNAstring set fasta file, containing only the title and nucleotide sequences
df_anopheles = data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))

#Display the first part of the dataframe, "df_anopheles".

head(df_anopheles)

#Converting the sequence column to a DNAStringSet (class) called 'nucleotides2' so that we can use functions from the Biostrings package.
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)

#Confirm class of nucleotide data (should be Biostring)
class(df_anopheles$nucleotides2)

#Display the first part of the dataframe, "df_anopheles".
head(df_anopheles)

#Calculate dinucleotide frequency (k-mers of length of 2) using the data in "nucleotides2" column) - useful as future predictor when creating a classifier.
df_anopheles = cbind(df_anopheles,
                     as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE)))

#Display the first part of the dataframe, "df_anopheles".
head(df_anopheles)

#Convert the regular dataframe, 'df_anopheles' to a tibble dataframe. Remove the column called 'nucleotides2'. 
df_anopheles = df_anopheles %>% 
  select(-nucleotides2) %>% 
  as_tibble()

#4. Data filtering by creating more dataframes and lists ####

#Based on the title column, create the additional following 2 columns: i) create a column called "species_name" that contains the 2nd and 3rd word from the title column, ii) create a column called "unique_id" that contains the 1st word from the title column
df_anopheles = df_anopheles %>% 
  mutate(species_name = word(title, 2L, 3L)) %>% 
  mutate(unique_id = word(title, 1L)) 

#View the tibble dataframe.
df_anopheles

#From df_anopheles, specifically the column_id: find the unique id numbers and count how many unique nucleotide sequences
df_anopheles$unique_id %>% 
  unique() %>% 
  length()

#Create creates a list of species names that appear more than 200 times in the dataset 
ls_vectors = df_anopheles %>% 
  group_by(species_name) %>% 
  count() %>% 
  arrange(desc(n)) %>% #View - commented out but would open the data in a viewer if uncommented
  filter(n > 200) %>% 
  pull(species_name)

#View the list of Anopheles species that is most represented in the dataset in descending order
ls_vectors

#Create another dataframe called 'df_anopheles_vector' that only contains the names of species found within the list, 'ls_vector' (species that appear more than 200 times in the df_anopheles dataframe).
df_anopheles_vector = 
  df_anopheles %>% 
  filter(species_name %in% ls_vectors)

#Take the just created dataframe, 'df_anopheles_vector' and group it based on species name and then count the frequency each name appears 
df_anopheles_vector %>% 
  group_by(species_name) %>% 
  count()

#Remove dataframe, 'df_anaopheles', and DNAStringSet, 'st_anopheles'. 
rm(df_anopheles, st_anopheles)

#5. Training Classification Model ####
#From the dataframe, 'df_anopheles_vector', find the smallest number of times a species name is recorded (we know that 200 times is the cutoff but we want to make sure that there is sample balance - important for validation and training dataset functionality).

smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample


#The smallest number of time a species is mentioned in the 'df_anopheles_vector' is 220 times.  

#Seed ensures reproducible random sampling.
#Validation dataset is created: Groups by species, then samples 20% of the smallest group size from each species

set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation

#View the validation dataset - confirm that 20% of the smallest group size from each species is taken - should be 44 samples 
df_anopheles_validation %>% 
  group_by(species_name) %>% 
  count()

##Seed ensures reproducible random sampling.
#Training dataset is created. It is important that the training dataset does not overlap with the validation dataset - this is acheived by asking for species names that are not (!) in the validation set using %in%. Second, among records remaining that are not in the validation set, 80% of the smallest group size from each species are chosen.

set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

#View the training dataset - confirm that 80% of the smallest group size from each species is taken - should be 176 samples (220-44 = 176)

df_anopheles_training %>% 
  group_by(species_name) %>% 
  count()

#List the column names in the dataframe, 'df_anopheles_training'
names(df_anopheles_training)

#6. Create Classifier, Train Classifer, & Running Analysis----
#Build random Forest Model.
#Create classifier to separate species using columns #3-18 (dinucleotide frequency) as predictors. Builds 100 decision tress in the forest.
anopheles_class <- randomForest(x = df_anopheles_training[, 3:18],
                                y = as.factor(df_anopheles_training$species_name), 
                                ntree = 100, importance = TRUE)

#View model summary
anopheles_class

#Plot model error - helps to determine if 100 trees is sufficient  
plot(anopheles_class)

#Variable Importance plot - Visualizes which variables (columns 3-18) are most important for species classification
varImpPlot(anopheles_class)

#Makes predictions on validation data.
predict_validation <- predict(anopheles_class, 
                              df_anopheles_validation[, c(19, 3:18)])
predict_validation

#Create confusion matrix. Compares actual species names vs. predicted species names in a cross-tabulation matrix to evaluate model accuracy.
table(observed = df_anopheles_validation$species_name, 
      predicted = predict_validation)

