##***************************
## Sabrina Saiphoo
## 2025-10-23
##***************************

## --- Packages  -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
conflicted::conflicts_prefer(dplyr::count)
conflicted::conflicts_prefer(base::as.factor)
library(viridis)
theme_set(theme_light())
library(Biostrings)
library(rentrez)
library(randomForest)
library(seqinr)

# ------ Determine Search Term -------------

# 1 - Shows how many records are found with the specific search term on NCBI Nucleotide
df_anopheles_search = entrez_search(db = "nuccore", term = "anopheles [ORGN] AND COI [gene]", use_history = T)
df_anopheles_search

# ------ Read/Format FASTA File -------------

#2 - Loading the FASTA files using readDNAStringSet, stores them as a DNAStringSet class
st_anopheles = readDNAStringSet("../data/sequence.fasta")
st_anopheles 

#3 - Turning the FASTA files into a dataframe
df_anopheles = data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
head(df_anopheles)

# ------ Calculate Dinucleotide Frequency -------------

#4 - Adding the nucleotides2 column that includes the sequence data back as a DNAStringSet class, once this is done the we now have the sequence data in the right format to run the dinucleotideFrequency function
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2)
head(df_anopheles)

#5 - Here we are calculating the dinucleotide frequencies of each sequence and binding them to our df_anopheles dataframe
df_anopheles = cbind(df_anopheles,
                     as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE)))
head(df_anopheles)

#6 - Here we are using select to remove the DNAStringSet column nucelotides2 as this it is unneeded now that we have the dinucleotide frequencies
df_anopheles = df_anopheles %>% 
  select(-nucleotides2) %>% 
  as_tibble()

# ------ Formatting Species & Unique IDs -------------

#7 - Here we are adding in two new columns to simplify the species name and also the individual id number, these columns will be added to the very right of our anopheles dataframe. More specifically, the word(title, 2L, 3L) is taking the value in the title column and pulling the second and third word, which makes up the species name. This new species name will be held in the mutated species_name column. The command, word(title, 1L) takes the first word from the value in the title column. The first word happens to be the unique id for the sequence and this will be the new value of the unique_id column.
df_anopheles = df_anopheles %>% 
  mutate(species_name = word(title, 2L, 3L)) %>% 
  mutate(unique_id = word(title, 1L)) 
df_anopheles

#8 - This shows how many unique ids are present in the dataset
df_anopheles$unique_id %>% 
  unique() %>% 
  length()

# ------ Feature Vector Preparation -------------

#9 - Here we are grouping by species and counting how many of each species there are. We are only keeping the species that have over 200 samples and then pull only the species_names. This will give us a set of valid classes for our RandomForest model.
ls_vectors = df_anopheles %>% 
  group_by(species_name) %>% 
  count() %>% 
  arrange(desc(n)) %>% #View()
  filter(n > 200) %>% 
  pull(species_name)
ls_vectors

#10 - We are now creating the data frame that will house both the features and the labels. We are starting by filtering the df_anopheles to only include species that are a part of the ls_vectors
df_anopheles_vector = 
  df_anopheles %>% 
  filter(species_name %in% ls_vectors)

#11 - Here we are displaying the count per species, where we can check to ensure all the counts are over 200
df_anopheles_vector %>% 
  group_by(species_name) %>% 
  count()

#12 - Here we are removing these so not get confused
rm(df_anopheles, st_anopheles)

# ------ Validation Data Preparation -------------

#13 - We are finding the smallest count of samples from a specific species, this is so we can stratify our dataset for the model. This will help us account for any imbalances.
smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample

#14 - Here we are setting a seed so that the results can be reproducible. We are creating our validation set which will be used at the end to validate our model. We are filtering out samples from each of the species to set aside. The number of samples that are set aside are 20% of the lowest count of samples for a species. We take this value to account for imbalances.
set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation

#15 - This specific line is a checkpoint for the previous one. In the previous step, we created the df_anopheles_validation dataset by taking a portion of each species’ samples that had over 200 entries and setting them aside. The portion that was taken was determined based on the species with the smallest sample count, where 20% of that smallest value was taken from each species to ensure balance across the validation set. This line of code checks that all count values are the same number and confirms that they are equivalent to 0.20 * smaller_sample.
df_anopheles_validation %>% 
  group_by(species_name) %>% 
  count()

# ------ Training Data Preparation -------------

#16 - In this step, we are creating our training dataset that does not include any samples found in the validation set. This helps remove potential bias, because if the model has already seen the validation samples during training, it could influence the results. Keeping the validation and training sets separate ensures no bias. The training dataset is created by taking 80% of the smallest species sample size from each group to maintain balance.(Not sure if the filter is filtering correctly as species_name and unique_id don't have the same values, should both be unique_id to properly remvove the validation data samples)
set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

#17 - Ensure that the number is the same for all
df_anopheles_training %>% 
  group_by(species_name) %>% 
  count()

#18 - Showing that the column headings of the training set include the features (the dinucleotide frequencies) and the labels (the species names)
names(df_anopheles_training)

# ------ RandomForest Classifier -------------

#19 - Training the classifier on the x = features, and providing the y = labels. Here we are making 100 decision trees. Columns 3:18 are considered the features (dinucleotide frequencies). The species names are the labels.
anopheles_class <- randomForest(x = df_anopheles_training[, 3:18],
                                y = as.factor(df_anopheles_training$species_name), 
                                ntree = 100, importance = TRUE)

# ------ Classifier Analysis -------------

#20 - This part of the code displays the Anopheles RandomForest classifier statistics. It shows different analysis plots and outputs that help us understand how well the model performed. The classifier was created using 100 decision trees and had a 1.7% OOB error rate, showing that the model performed well. It also shows a confusion matrix that demonstrates how many times one species was classified as another. When plotting anopheles_class, we can see the error rate compared to the number of trees per species, where many of the species begin to plateau between 20–40 trees. It is clear that as the number of trees increases, the error rate starts to decrease. Lastly, the varImpPlot is a variable importance plot that shows which dinucleotide frequencies were the most useful in helping the model make accurate classifications.
anopheles_class
plot(anopheles_class)
varImpPlot(anopheles_class)


# ------ Validating the Classifier -------------

#21 - Here we are testing out the model on the validation dataset. This is to see how well it predicts unseen data. Calling on the feature columns and the labels column.
predict_validation <- predict(anopheles_class, 
                              df_anopheles_validation[, c(19, 3:18)])
predict_validation

# ------ Classifier Results -------------

#22 - This is essentially a confusion matrix showing how well each of the observed species were classified. It shows if the observed species was the same as the predicted. If the observed was the same as the predicted then the matrix coordinate will increase by 1. Example: Anopheles triannulatus was mis-classified as Anopheles albimanus once. 
table(observed = df_anopheles_validation$species_name, 
      predicted = predict_validation)




