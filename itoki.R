## _ Packages used -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter, dplyr::combine, base::as.factor)
library(viridis)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())
library(Biostrings)
library(rentrez)
library(seqinr)
library(dplyr)
library(randomForest)

#Code comments below each block of code

##1. Search and read in data-----
df_anopheles_search = entrez_search(db = "nuccore", term = "anopheles [ORGN] AND COI [gene]", use_history = T)
df_anopheles_search
#Search for data on entrez for taxonomic class Anopheles and COI marker gene gives us 10730 hits shows us what to expect

st_anopheles = readDNAStringSet("../data/sequence.fasta")
st_anopheles 
#Read in anopheles data downloaded from NCBI with the previously used search parameters 

##2. Data Transformation for analysis----

df_anopheles = data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
head(df_anopheles)
#Creating dataframe from the DNAStringset 


df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2)
head(df_anopheles)
#Creating nucleotides2 column as a DNAstringset to function in biostrings 


df_anopheles = cbind(df_anopheles,
                     as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE)))
head(df_anopheles)
#Using the biostrings package to generate dinucleotide frequency from nucleotides2 column

df_anopheles = df_anopheles %>% 
  select(-nucleotides2) %>% 
  as_tibble()
#remove the nucleotides2 column and convert our dataframe to tibble to use on tidyverse

df_anopheles = df_anopheles %>% 
  mutate(species_name = word(title, 2L, 3L)) %>% 
  mutate(unique_id = word(title, 1L)) 
df_anopheles
#Creating columns unique_id and species_name from the Title column

df_anopheles$unique_id %>% 
  unique() %>% 
  length()
#Crosschecking our dataset to see our sample unique IDs match what we have from our search


ls_vectors = df_anopheles %>% 
  group_by(species_name) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  filter(n > 200) %>% 
  pull(species_name)
ls_vectors
#creating a list of vectors by filtering for  species with more than 200 samples

df_anopheles_vector = 
  df_anopheles %>% 
  filter(species_name %in% ls_vectors)
#Keeping only samples for our vectors in our new dataframe

df_anopheles_vector %>% 
  group_by(species_name) %>% 
  count()
#Inspecting to see that our data transformation worked correctly and it did since none of the values are < or = 200



rm(df_anopheles, st_anopheles)
#Remove dataframes we no longer need


##3. Random forest Training and validation dataset creation-------

smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample
#Smallest number of samples for a species is 220

set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation
#Creating Validation data set with 20%(44) of our smallest  specie group for each specie    sample group this will be separated from the training data for our prediction model.

df_anopheles_validation %>% 
  group_by(species_name) %>% 
  count()
#Confirming that we have 44 samples for each group


set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!unique_id %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))
#Next we create our training data set by first removing all samples present in our validation data set, group by species name and take 80%(176) of our smallest specie sample group for each specie 

#Changed species name to unique_id in filter I think it works better

df_anopheles_training %>% 
  group_by(species_name) %>% 
  count()
# Confirming that we have 176 samples for each group.

##4. Random forest Prediction model-----

names(df_anopheles_training) 
#ensuring all columns are present in our training data set and checking column position


anopheles_class <- randomForest(x = df_anopheles_training[, 3:18],
                                y = as.factor(df_anopheles_training$species_name), 
                                ntree = 100, importance = TRUE)
#Generating a random forest with our training data that has  100 decision trees using dinucleotide frequency(columns 3-18) as our predictor and Species name as our response variable 

anopheles_class
plot(anopheles_class)
varImpPlot(anopheles_class)
#Data exploration of our random forest


predict_validation <- predict(anopheles_class, 
                              df_anopheles_validation[, c(19, 3:18)])
predict_validation

#Testing out our random forest on the validation dataset



table(observed = df_anopheles_validation$species_name, 
      predicted = predict_validation)

#Observing the accuracy of our analysis, 98% accuracy.

















