# Script Name: anopheles_coi_classification.R
# Lishita Rowjee
# Date:  29th October 2025

# Section 1:Install required packages ----
#(BiocManager,tidyverse, Biostrings, rentrez, randomForest)


# Section 2: Load libraries----

library(randomForest)
library(rentrez)
library(Biostrings)
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())


# Section 3: Search for COI sequences in NCBI----

df_anopheles_search = entrez_search(db = "nuccore", term = "anopheles [ORGN] AND COI [gene]", use_history = T)
df_anopheles_search

#The number of hits from NCBI is the same (10730) as that on the NCBI website after searching for "anopheles [ORGN] AND COI [gene]" in the nucleotide database.


#Section 4: Download and load sequence data----
st_anopheles = readDNAStringSet("../data/sequence.fasta")
st_anopheles 

#Section 5: Convert DNAStringSet into a standard dataframe ----
#This is for easier manipulation
df_anopheles = data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
head(df_anopheles) #Look at first few rows

#Convert text sequences back into DNAStringSet  for analysis
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2) #confirms class is DNAStringSet
head(df_anopheles)


#Section 6: Calculate dinucleotide frequencies----

#Dinucleotide frequencies: the frequency of every possible dinucleotide is computed for each sequence in the dataframe.

df_anopheles = cbind(df_anopheles,
                     as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE)))
head(df_anopheles)

#Section 7: Clean and organize dataframe----

#clean-up dataframe by removing the nucleotides2 column and converting the dataframe to a tibble for easier viewing
df_anopheles = df_anopheles %>% 
  select(-nucleotides2) %>% 
  as_tibble()

#Section 8: Extract species and Unique IDs----

#Extracting the species name and unique identifier from the titles of the different sequences. The following code says to extract the 2nd and 3rd word as the species_name and 1st word as the unique_id.

df_anopheles = df_anopheles %>% 
  mutate(species_name = word(title, 2L, 3L)) %>% 
  mutate(unique_id = word(title, 1L)) 
df_anopheles #view the new data frame

#Counting the number of unique sequences in the data frame
df_anopheles$unique_id %>% 
  unique() %>% 
  length()
# 10730 obtained confirming the number of unique sequences in data set

# Section 9: Filter species > 200 species----

#Only selecting species having greater than 200 sequences from the data set

ls_vectors = df_anopheles %>% 
  group_by(species_name) %>% 
  count()  %>% 
  arrange(desc(n)) %>% #View() #count how many sequences belongs to each species and arrange in descending order
  filter(n > 200) %>% 
  pull(species_name)
ls_vectors #provide list of selected species name

#Filter the data set to only include the species having > 200 sequences
df_anopheles_vector = 
  df_anopheles %>% 
  filter(species_name %in% ls_vectors)

#checking the number of sequences per species after filtering to ensure only species with >200 sequences have been kept
df_anopheles_vector %>% 
  group_by(species_name) %>% 
  count()

#remove the large objects as they are no longer necessary since df_anopheles_vector has been created already
rm(df_anopheles, st_anopheles)

#Determining the smallest number of sequences available across the selected species to ensure balanced sampling

smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample #smallest number of sequences is 220

#Section 10: Create validation and trining dataset----

#Validation dataset is created by set.seed(51) to ensure reproducibility first and then randomly sampling 20% of each species for validation

set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation

#After creating validation data set, check that each species has the same number of validation samples
df_anopheles_validation %>% 
  group_by(species_name) %>% 
  count()

# The training data set is created with the remaining 80% of species and set.seed(40) is to ensure reproducibility
set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

#confirm that there is balanced sampling per species
df_anopheles_training %>% 
  group_by(species_name) %>% 
  count()

#checking the structure of the data structure
names(df_anopheles_training) #visualized the difference column names


#Section 11: Training a random Forest classifier----

#Dinucleotide frequencies (columns 3–18) as predictors (X)and species_name as the target variable (Y). Also R was told to make 100 decision trees

anopheles_class <- randomForest(x = df_anopheles_training[, 3:18],
                                y = as.factor(df_anopheles_training$species_name), 
                                ntree = 100, importance = TRUE)

#A conflict error message was presented because "as.factor" is found in 3 packages. It got solved using the code below.
conflicted::conflicts_prefer(IRanges::as.factor)

# Section 12: Evaluate and visualize model performance----

anopheles_class #Provides summary of model performance and OOB error rate
plot(anopheles_class) 
# Error rate is plotted as the number of trees increases
varImpPlot(anopheles_class) #a plot is generated to show which dinucleotide is the most important to accurately determine species classification

#Section 13: Validate model on test data ----

#Validating the model on the test set (20%) by using the trained model to predict species on the validation data
predict_validation <- predict(anopheles_class, 
                              df_anopheles_validation[, c(19, 3:18)])
predict_validation

#Comparing observed species vs predicted species
table(observed = df_anopheles_validation$species_name, 
      predicted = predict_validation)