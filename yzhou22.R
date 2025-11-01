#
# Title: Read and Analyze GenBank Sequences 
# Author: XXX
# Date: 2025/10/25
#
# Description:
# This script connects to NCBI’s GenBank to retrieve nucleotide sequences of Anopheles # #species. 
#
# Input: GenBank species keywords.
#
# Output: FASTA sequence files and summary tables.
#
# Libraries:
# - rentrez, tidyverse,dplyr
############################################################

##outlines:

#1.load libraries
#2. Retrieve GenBank Dataand load/read local fasta file 
#3. convert Sequences to Dataframe
#4. Calculate Dinucleotide Frequencies
#5. filter out dataframe: Extract Species Name and Unique ID from title column
#6. filter Species for Sample Size
#7. Split Data into Training and Validation Sets. validation set: 20% of smaller sample; Training set: remaining 80%
#8. Random Forest Classifier; plot(anopheles_class), varImpPlot(anopheles_class)
#9. Evaluate Model on Validation Set; generate Confusion matrix
####


## _ Packages used -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(viridis)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())
library(Biostrings)
library(rentrez)
library(tidyverse)
library(Biostrings)
library(stringr)

#Search GenBank for Anopheles COI sequences
df_anopheles_search = entrez_search(db = "nuccore", 
                                    term = "anopheles [ORGN] AND COI [gene]",
                                    use_history = T)
df_anopheles_search


#Read FASTA file of sequences from local folder
st_anopheles = readDNAStringSet("../data/sequence.fasta")
st_anopheles 

#Convert sequences to dataframe
df_anopheles = data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
## preview first few rows
head(df_anopheles)

#Convert strings to DNAStringSet 
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2)
## preview first few rows
head(df_anopheles)

#Calculate dinucleotide frequencies
df_anopheles = cbind(df_anopheles,
                     as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE)))
## preview first few rows
head(df_anopheles)

#Clean the dataframe and convert to tibble(), and remove unnecessary columns 
df_anopheles = df_anopheles %>% 
  select(-nucleotides2) %>% 
  as_tibble()

#Extract species_name and unique_ID from title
##df_anopheles data frame will have two new columns: unique_id and species_name
#unique column extracts the first word from the title column
#species_name extracts 2nd the 3rd from the title column
df_anopheles = df_anopheles %>% 
  mutate(species_name = word(title, 2L, 3L)) %>% 
  mutate(unique_id = word(title, 1L)) 
df_anopheles


#Count number of unique sequences
df_anopheles$unique_id %>% 
  unique() %>% 
#total number of the sequences
  length()

#sort and Identify species with >200 sequences for analysis;return species names as a vector
ls_vectors = df_anopheles %>% 
  group_by(species_name) %>% 
  count() %>% 
  arrange(desc(n)) %>% #View()
  filter(n > 200) %>% 
  pull(species_name)
ls_vectors

#Filter dataframe to only include selected species; and to verify the count() at the end
df_anopheles_vector = 
  df_anopheles %>% 
  filter(species_name %in% ls_vectors)

df_anopheles_vector %>% 
  group_by(species_name) %>% 
  count()

rm(df_anopheles, st_anopheles)

#Determine the smallest sample size among selected species
smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample

#Split data into validation dataset (20% of the smaller sample); and check counts 
set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation

#returns a table showing how many entries belong to each species_name
df_anopheles_validation %>% 
  group_by(species_name) %>% 
  count()

#load the rest of the 80% into training set and check id and counts 
set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
## exclude validation IDs
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))


df_anopheles_training %>% 
  group_by(species_name) %>% 
  count()

#recall and verify 
names(df_anopheles_training)


#Train Random Forest model to classify species, based on dinucleotide frequencies, response variable and tree numbers; variable improtance active
anopheles_class <- randomForest(x = df_anopheles_training[, 3:18],
                                y = as.factor(df_anopheles_training$species_name), 
                                ntree = 100, importance = TRUE)


#showes the summary of anopheles_class model 
#plots the ramdom forest model 

#showes the variable importance plot of which section contributes most to the model
#model performance 
anopheles_class

## OOB error vs number of trees
plot(anopheles_class)
varImpPlot(anopheles_class)


#Predict species on validation set; with 3-18 columns for prediction
predict_validation <- predict(anopheles_class, 
                              df_anopheles_validation[, c(19, 3:18)])
predict_validation


# Confusion matrix table: observed vs predicted
table(observed = df_anopheles_validation$species_name, 
      predicted = predict_validation)

