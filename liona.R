####1. Required packages----
library("rentrez")
library("randomForest")
library("tidyverse")
conflicted::conflicts_prefer(dplyr::count)
conflicted::conflicts_prefer(dplyr::filter)
library("Biostrings")

####2. Loading our data, checking, and filtering our data----
df_anopheles_search <- entrez_search(db = "nuccore", term = "anopheles [ORGN] AND COI [gene]", use_history = T)

df_anopheles_search

#Reading in the file, has 10730 elements which matches the count from df_anopheles search
st_anopheles <- readDNAStringSet("../data/sequence.fasta")
st_anopheles 

#Generating a dataframe with different columns for names and sequence
df_anopheles <- data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))

#Checking out the first few rows of df_anopheles
head(df_anopheles)

####3. Sequence features, adding species name and IDs, and generating dinucleotide frequencies ----
#Storing the sequence column (currently character class) to DNA string set sequence
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2)
head(df_anopheles)

#Generating dinucleotide kmers and making separate columns for each unique kmer
df_anopheles <- cbind(df_anopheles,
                      as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE)))

#Checking the headers to ensure that the kmers columns were made
head(df_anopheles)

#Removing the nucleotides sequence column and making tibble class
df_anopheles <- df_anopheles %>% 
  select(-nucleotides2) %>% 
  as_tibble()

#Generating two columns from the title column, species name (1st word) and unique ID (2nd and third word)
df_anopheles <- df_anopheles %>% 
  mutate(species_name = word(title, 2L, 3L)) %>% #
  mutate(unique_id = word(title, 1L)) 
df_anopheles

#Checking to see if the numbers of unique ids were created properly 
df_anopheles$unique_id %>% 
  unique() %>% 
  length()

####4. Generating training and validation data ----

#Grouping by species name, counting how many species show up and filtering and pulling
#names that have ones greater than 200 counts.
ls_vectors <- df_anopheles %>% 
  group_by(species_name) %>% 
  count() %>% 
  arrange(desc(n)) %>% #View()
  filter(n > 200) %>% 
  pull(species_name)
ls_vectors

#Generating a dataframe filtering species name that are found in both species name 
#of the anopheles vector and ls vectors
df_anopheles_vector <- df_anopheles %>% 
  filter(species_name %in% ls_vectors)

#Counting the number of occurences for each species name
df_anopheles_vector %>% 
  group_by(species_name) %>% 
  count()

#Finds the lowest count from the species name, ensures validation and training 
#data is equally represented.
smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample

#Set seed to ensure reproducibility
set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample)) #floor rounds down integers
df_anopheles_validation

# Counting the number of records for each species name present
df_anopheles_validation %>% 
  group_by(species_name) %>% 
  count()

#set seed to make our code reproducible
set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

# Counting the number of records for each species present
df_anopheles_training %>% 
  group_by(species_name) %>% 
  count()

#Looks at header names of df_anopheles_training
names(df_anopheles_training)

####5. Training classification model using random forests and plotting data ----
conflicted::conflicts_prefer(base::as.factor)
anopheles_class <- randomForest(x = df_anopheles_training[, 3:18],
                                y = as.factor(df_anopheles_training$species_name), 
                                ntree = 100, importance = TRUE)

#Checking random forest model and checking its performance by plotting
#error rate plot and variable importance plot
anopheles_class
plot(anopheles_class) 
varImpPlot(anopheles_class)

#Error rates plateaus as more trees are added, we had 100 trees so error rate is about less than 0.05

#Checks to see if our prediction classifier works
predict_validation <- predict(anopheles_class, 
                              df_anopheles_validation[, c(19, 3:18)])
predict_validation

#Compares the predictions from classifier model to the actual labels for each species
table(observed = df_anopheles_validation$species_name, 
      predicted = predict_validation)
#Anopheles triannulatus for example, was classified 43 times as Anopheles triannulatus,
# missaclassfied 1 time as Anopheles maculipennis, and 0 times for other species.
#Indicating the classifier model is working well.

#no longer needed so we remove
rm(df_anopheles, st_anopheles)
