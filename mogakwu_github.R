##### BINF 6890: Software Tools: Anopheles corrected script file with comments
##*
## Maryanne Ogakwu 1395098
##
## 2025-10-24
##
##***************************

## _ Packages used -------
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(viridis)
# + scale_color/fill_viridis(discrete = T/F)
theme_set(theme_light())
library("vegan")
library(Biostrings)
install.packages("BiocManager")
BiocManager::install("GenomeInfoDb")
`force = TRUE`
library(randomforest)
install.packages("rentrez")
library(rentrez)
install.packages("seqinr")
library(seqinr)

# Setting working directory
getwd()
setwd("C:/Users/marya/OneDrive/Dokumenty/Masters/binf_6210/Oct 23 Malaria Vector Classifier/R script/")

# Startup ends here

# PART 1: OBTAINING THE DATABASE INFO ------

# This line is to search through the NCBI GenBank database to obtain the Anopheles genus and COI gene data, and prints the search results in a summary. This is an easier alternative method to obtain data, rather that searching through the NCBI database manually.
df_anopheles_search <- entrez_search(
  db = "nuccore",
  term = "anopheles [ORGN] AND COI [gene]",
  use_history = T
)
df_anopheles_search

# PART 2: READING THE OBTAINED DATA INTO A DATASET ---------------

# Now that the data has beem uptained, we use the readDNAStringSet to read the the FASTA file into R as a dataset.
st_anopheles <- readDNAStringSet("../data/sequence.fasta")
# Print a summary of the dataset to show the number of sequences that were loaded
st_anopheles

# PART 3: CREATING A DATAFRAME -------

# To be able to view the data properly, manipulate and analyze its columns, it needs to be loaded into a dataframe.So two coulns will be produced: the sequence name from the FASTA file header and the DNA sequences.
df_anopheles <- data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
# Displays the first few rows of the dataframe, and allows us to have an idea of how the dataframe is arranged without opening it.
head(df_anopheles)

## PART 4: INSPECTING THE DNA SEQUENCES -------
# For the DNA sequences to be read and analyzed, it first has to be converted to a format that Biostrings ca analyze, so we are craeating a new column where only the sequences are stored as a DNAStringSet.The new column with the DNA sequences is called nucleotides2.
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
# Checking the class of the new column in the dataframe
class(df_anopheles$nucleotides2)
# Checking the first few rows of the dataframe.
head(df_anopheles)

# PART 5: OBTAINING DINUCLEOTIDE SEQUENCES ---------

# We are calculating the frequency of different dinucleotide pairs in each of the sequences, as.prob = TRUE is used to display the frequencies as numeric values instead of raw counts. cbind is used to combine the frequency results, this leads to each sequence having a numeric feature that describes it's nucleotide composition.
df_anopheles <- cbind(
  df_anopheles,
  as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE))
)
# Checks the first few rows of the dataframe
head(df_anopheles)

# Since we have obtained our nucleotide frequencies for each of the sequences, we no longer need the nucleotides2 column, so we are deleting it from the dataframe. as_tibble will convert the data frame into a tidyverse tibble.
df_anopheles <- df_anopheles %>%
  select(-nucleotides2) %>%
  as_tibble()

## PART 6:  EXTRACTING SPECIES AND ID (ASSENCION) INFO ------

# To get the species names and IDs from the column, were extracting the second word (2L) and third word (3L) from the main column and adds thes to new columns. The assecion ID  will be extracted from the first word (1l). This done using the mutate function.
df_anopheles <- df_anopheles %>%
  mutate(species_name = word(title, 2L, 3L)) %>%
  mutate(unique_id = word(title, 1L))
# Prints then rsults of the muated dataset
df_anopheles

# As we've extracted the unique IDs, we need gto find out how many of them are unique, we'll use the code belowto count how many of the IDs (asscession number) are unique while avoiding duplicates.
df_anopheles$unique_id %>%
  unique() %>%
  length()

# Counting how many sequences appear for each species, and grouping the results by species.We are also filtering to keeps only the sequences that appear more than 200 times. These results are saved to new list called ls_vectors.

## PART 7: BUILDING DATASETS WITH ENOUGH SAMPLES --------
ls_vectors <- df_anopheles %>%
  group_by(species_name) %>%
  count() %>%
  arrange(desc(n)) %>% # View()
  filter(n > 200) %>%
  pull(species_name)
# Prints the objects in the ls_vectors list to allow us to view them.
ls_vectors

# This filters our df_anopheles dataset to only keep the frequently occuring species and assigns it to another dataframe called df_anopheles_vactor.
df_anopheles_vector <-
  df_anopheles %>%
  filter(species_name %in% ls_vectors)
# Counts how many of the sequences remain for each species.
df_anopheles_vector %>%
  group_by(species_name) %>%
  count()

# To ensure our environment isn't cluttered, we remove excess objects that we have derived our data from. So these are no longer needed.
rm(df_anopheles, st_anopheles)

## PART 8: PREPARING TRAINING DATASET AND VALIDATION DATASET -------

# This finds the species with the least occuring number of sequences in the df_anopheles_vector dataset, and assigns the results to the object called smaller_sample.
smaller_sample <- min(table(df_anopheles_vector$species_name))
# Prints the smaller sample data
smaller_sample

# We are setting a random seed to ensure the anaysis is reproducible and we get the same reults each time this line of code is run.
set.seed(51)
# Randomly selects 20% of sequences for every species in the df_anopheles_vector dataset, to be used for a validation test.Selecting only 20% makes sure that every species is equally represented in the dataset asto avoid skewed and biased results.
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
# Prints the data in the validation set.
df_anopheles_validation

# This ensures that each species is well represented in the training data, and that each species have has the same number of sequences. This is because most of software tools tend to have a bias towards larger datasets, ultimately and indiscriminately favoring them in analysis results due to the sample size.
df_anopheles_validation %>%
  group_by(species_name) %>%
  count()

# Sets another random seed for tking samples out of the training dta, once again making sure that the results are reproducible
set.seed(40)
# Samples 80% of the sequences in the df_anopheles_vector dataset to create a training data. Filtering the species name in the dataframe by the unique ID and grouping the reults by the species name.
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

## PART 9: TRAINING RANDOM FOREST MODEL ----
# This checks the data generated in the df_anopheles_training set is balanced and has the correct number of samples, by grouping the data by species name and counting the results.
df_anopheles_training %>%
  group_by(species_name) %>%
  count()
# Checks the column names
names(df_anopheles_training)

## PART 10: RANDOM FOREST MODEL EVALUATION ---------

# Here we are building a random forest model using data we derived and saved in our training dataset.
anopheles_class <- randomForest(
  x = df_anopheles_training[, 3:18],
  y = as.factor(df_anopheles_training$species_name),
  ntree = 100, importance = TRUE
)

## PART 11: EVALUATION OF MODEL PERFORMANCE ------

# Prints the model summary to the console
anopheles_class
# creates a plot of the random forest error rate as the trees increases.
plot(anopheles_class)
# Creates a plot that shows which dinucleotide frequencies are the most important when it comes to classifying species.
varImpPlot(anopheles_class)

# Makes predictions on the validation set using the trained model we created earlier. Mkaes use of the dinucleotide frequency columns that exist in the datasets

## PART 12: VALIDATION OF MODEL PREDICTIONS ------
predict_validation <- predict(
  anopheles_class,
  df_anopheles_validation[, c(19, 3:18)]
)
# Prints the predicted species name sfor each sequence in the console
predict_validation

# Makes a table that compare the true "observed" species and the mode predicted species to measure  accuracy of the results.
table(
  observed = df_anopheles_validation$species_name,
  predicted = predict_validation
)

usethis::create_from_github(
  "https://github.com/karl-cottenie/binf25.git",
  # replace the above github url with the green Code > Local > Clone using the web URL
  # you need to update this for each repository you want to contribute to
  destdir = "./binf25/", 
  # or whatever directory location you want
  # this destdir can be the same across different repositories if you want all your github projects to be in the same location
  # each repository will be in a different subdirectory in this location
  fork = TRUE
)
usethis::git_remotes()
usethis::pr_init(branch = "Maryanne_binf6210_github")
usethis::git_sitrep()
#Trying to push again
usethis::pr_push()
