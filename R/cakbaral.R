## CLASSIFICATION OF ANOPHELES COI GENE DATABASE SEQUENCES BY SPECIES----

#The purpose of this script was to build a classifier model that could accurately predict and classify COI gene sequences of Anopheles, retrieved from the NCBI Nucleotide database, by their species name.

#_ R packages used and conflict preferences----
library(tidyverse)
library(viridis)
library(dplyr)
library(Biostrings)
library(rentrez)
library(conflicted)
library(randomForest)
conflict_prefer("filter", "dplyr")
conflicts_prefer(base::as.factor)

##1: NCBI DATABASE SEARCH FOR ANOPHELES AND READING IN COI GENE SEQUENCE DATA----

#Step 1: The rentrez package used the entrez_search() command to check that the total search results from NCBI matched all sequences recorded in the database, based on the search terms given for COI genes in Anopheles, along with the database used (Nucleotide).
df_anopheles_search = entrez_search(db = "nuccore", term = "anopheles [ORGN] AND COI [gene]", use_history = T)
df_anopheles_search


#Step 2: The FASTA file containing all 10730 total hits were read in as a DNAStringSet using the Biostrings package, and the following summarizes the information about the nature of the data.
st_anopheles = readDNAStringSet("../data/sequence.fasta")
st_anopheles


##2: CONSTRUCTING THE ANOPHELES COI SEQUENCE DATA FRAME AND CLEANING THE DATA----

#Step 3: All sequence hits were summarized in a data frame to visualize the information about the Anopheles species samples and their corresponding COI gene sequences from the NCBI search.
df_anopheles = data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
head(df_anopheles)

#Step 4: To calculate frequencies using the Biostrings pacakge, all sequences needed to be converted into DNAStringSet data.
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2)
head(df_anopheles)

#Step 5: Using the sequence data from the DNAStringSet, all dinucleotide frequencies for each sample were calculated and appended to the data frame.
df_anopheles = cbind(df_anopheles,
                     as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE)))
head(df_anopheles)

#Step 6: Now that the dinucleotide frequencies were calculated, the DNAStringSet sequences were removed and reassigned as a tibble.
df_anopheles = df_anopheles %>% 
  select(-nucleotides2) %>% 
  as_tibble()

#Step 7: In this step, the species names and unique species IDs were extracted from the NCBI entry titles.  Using the mutate() command, two new columns were appended to the data frame in order to visualize the list of all species and IDs.
df_anopheles = df_anopheles %>% 
  mutate(species_name = word(title, 2L, 3L)) %>% 
  mutate(unique_id = word(title, 1L)) 
df_anopheles


#Step 8: By including this check, this ensures that the data frame still includes all total results from the NCBI database.
df_anopheles$unique_id %>% 
  unique() %>% 
  length()

#Step 9: The names of all Anopheles species were listed out to count how many different ones were included, ensuring that only species with a sufficient number of entry samples (>200) were considered for downstream analysis.
ls_vectors = df_anopheles %>% 
  group_by(species_name) %>% 
  count() %>% 
  arrange(desc(n)) %>% #View()
  filter(n > 200) %>% 
  pull(species_name)
ls_vectors

#Step 10: To clean the data and only include these species that are part of the aforementioned list, Anopheles species that contained an insufficient number of entry samples (<=200) were removed from the data frame.
df_anopheles_vector = 
  df_anopheles %>% 
  filter(species_name %in% ls_vectors)

#Step 11: With this command, the sample counts for each type of species were listed. This also helps to ensure that sample counts were all >200, which, in this case, is TRUE.
df_anopheles_vector %>% 
  group_by(species_name) %>% 
  count()

#Step 12: With the data cleaned for classifier training, since the original data frame and DNAString from the FASTA file were no longer needed for downstream analysis, they were removed from the environment.
rm(df_anopheles, st_anopheles)

#Step 13: The smaller_sample variable listed the minimum entry sample size amongst these counts. This also doubles as a secondary check to ensure that even the Anopheles species with the lowest sample size was >200.
smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample


##3: ANOPHELES CLASSIFIER MODEL TRAINING----

#Step 14: The random seed was set, and the validation data set was created.  The number of samples taken for each class of Anopheles was 20% of the minimum sample size as defined by "smaller_sample", which was (0.20 * 220) = 44.
set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation

#Step 15: It was essential that all Anopheles species had equal sample sizes in the validation set to avoid class imbalance, which could result in bias by majority. By running this command, this list of counts checked that all sample sizes across the classes were indeed equal to 44.
df_anopheles_validation %>% 
  group_by(species_name) %>% 
  count()

#Step 16: Next, the seed was set again, and training dataset was retrieved by sampling 80% of the minimum sample size for each class of Anopheles (0.80 * 220 = 176). This data frame should not overlap with the validation, which was accomplished by reverse matching entries from its data frame (i.e. the training dataset contains samples NOT from the validation dataset).
set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

#Step 17: Once more, this check ensured that all Anopheles species contained equal sample sizes for the training dataset, which were indeed all equal to 176.
df_anopheles_training %>% 
  group_by(species_name) %>% 
  count()

#Step 18: This check listed all of the categorical names of the training dataset to ensure all relevant information was included, such as title, sequence, species name and ID, and dinucleotide frequencies.
names(df_anopheles_training)


#Step 19: The classifer model was constructed to sort the COI sequences by the species of Anopheles, using the "species_name" column of the training data frame. The number of decision trees in the forest was set to 100.
anopheles_class <- randomForest(x = df_anopheles_training[, 3:18],
                                y = as.factor(df_anopheles_training$species_name), 
                                ntree = 100, importance = TRUE)

#Step 20: Now that the classifier model has been defined, it was trained using the sample training data set. The first resulting plot compared how the error rate was affected based on the number of decision trees, and the variable importance plot illustrated how different variables were prioritized in classifying the COI sequences accurately.
anopheles_class
plot(anopheles_class)
varImpPlot(anopheles_class)


##4: PERFORMANCE EVALUATION OF ANOPHELES COI CLASSIFIER MODEL ON UNKNOWN DATA SETS----

##Step 21: Using the validation dataset as an example, this variable helped to confirm that the model can work as a proper method for predicting the classifications of any unknown dataset of Anopheles COI gene sequences.
predict_validation <- predict(anopheles_class, 
                              df_anopheles_validation[, c(19, 3:18)])
predict_validation


#Step 22: Lastly for this analysis, a confusion matrix was constructed to compare the observed and predicted values, serving as a final evaluation of the model's performance in case further adjustments or amendments are required.
table(observed = df_anopheles_validation$species_name, 
      predicted = predict_validation)