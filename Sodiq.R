library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(viridis)
install.packages("rentrez")
library(rentrez)
library(Biostrings)
library(dplyr)
library(randomForest)
library(tidyr)
library(stringr)

#
#search for the anopheles in the ncbi database#
df_anopheles_search = entrez_search(db = "nuccore", term = "anopheles [ORGN] AND COI [gene]", use_history = T)
df_anopheles_search

#read the data from the directory#
st_anopheles <- readDNAStringSet("../data/sequence.fasta")
st_anopheles 

#rename the st_anopheles to become a data frame of df_anopheles)
df_anopheles <- data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
head(df_anopheles)
#Convert sequences to DNAStringSet to enable us know the amount#
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2)
head(df_anopheles)
#compute dinucleotides as probabilities from the nucleotides2 in the data frame#
df_anopheles <- cbind(df_anopheles,
                     as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE)))
head(df_anopheles)

#Extract two new columns and give them these names; species_name and unique_id#
df_anopheles <- df_anopheles %>% 
  mutate(species_name = word(title, 2L, 3L)) %>% 
  mutate(unique_id = word(title, 1L)) 
df_anopheles

#Check the newly created column unique_id for its length# 
df_anopheles$unique_id %>% 
  unique() %>% 
  length()
#Remove the nucleotides2 column
df_anopheles <- df_anopheles %>% 
  select(-nucleotides2) %>% 
  as_tibble()

#listing vectors only in the species column with 200 as the size#
ls_vectors <- df_anopheles %>% 
  group_by(species_name) %>% 
  count() %>% 
  arrange(desc(n)) %>% #View()
  filter(n > 200) %>% 
  pull(species_name)
ls_vectors

#create a new data framwe that keeps only the species_name#
df_anopheles_vector <-  df_anopheles %>% 
  filter(species_name %in% ls_vectors)
#making a head count for the each specie in the data frame#
df_anopheles_vector %>% 
  group_by(species_name) %>% 
  count()

#select the specie with the minimum value#
smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample

#Create a new data frame from data frame anopheles vector with 20% of the smallest size#
set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation

#count the number of species present in the dataframe#
df_anopheles_validation %>% 
  group_by(species_name) %>% 
  count()

#Create a new data frame by filtering out the sequences already in the validation set with 80% of the smallest size#
set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

#count the number of species present in the dataframe#
df_anopheles_training %>% 
  group_by(species_name) %>% 
  count()
names(df_anopheles_training)

# Train Random Forest model to classify Anopheles species b7 using dinucleotide frequencies in column 3 to 18
anopheles_class <- randomForest(x = df_anopheles_training[, 3:18],
                                y = base::as.factor(df_anopheles_training$species_name), 
                                ntree = 100, importance = TRUE)
#view the new dataframe that has the classification of anopheles species#
anopheles_class

#plot a graph showing anopheles class#
plot(anopheles_class)
#this plots the dinucleotide frequency that has the most contributions#
varImpPlot(anopheles_class) 

#predict the species on the validation dataframe#
predict_validation <- predict(anopheles_class, 
                              df_anopheles_validation[, c(19, 3:18)])
predict_validation

#make a table showing the observed and predicted#
table(observed = df_anopheles_validation$species_name, 
      predicted = predict_validation)

#remove data frame#
rm(df_anopheles, st_anopheles)







