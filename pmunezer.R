##### Software Tools (BINF*6210) – Fall 2025 – Quiz #2 ----
##***************************
##  QUIZ 2 - RANDOMFOREST - ANOPHELES - COMMENT THE CODES
##
## Pierre Celestin Munezero (1378625)
##
## Instructor: Karl Cottenie
##
## 2025-10-26
##
##***************************
## Main tasks:
##
## PART 1: OVERVIEW
## PART 2: FETCHING RECORDS OF DNA SEQUENCES OF ANOPHERES
## PART 3: DATA PARSING AND CLEANING
## PART 4: TRAINING THE CLASSIFICATION MODEL
## PART 5: VERIFICATION OF MODEL GENERALIZATION

# _ Comment codes ------

# Coding explanation (#)
# Solution/result/interpretation (#==>)

## _ LOADING PACKAGES -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::count)
conflicted::conflicts_prefer(base::as.factor)
library(viridis)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())
library(Biostrings)
library(rentrez)
library(seqinr)
library(randomForest)

# _ PART 1: OVERVIEW -------

# Anopheles genus includes more than 460 species across 7 subgenera, with Anopheles, Cellia, Kerteszia, and Nyssorhynchus being the main malaria vectors (https://doi.org/10.5772/54695). These mosquitoes, particularly Anopheles, show remarkable diversity that shapes malaria transmission dynamics. Some molecular markers, such as COI and ITS2, are critical for accurate identification of cryptic species (https://doi.org/10.1186/s13071-014-0592-5). Machine learning is an important tool that can be used to classify Anopheles species. In the following workflow, COI is used to make a classification model for Anopheles species using DNA sequences (nucleotides) from NCBI.

# _ PART 2: FETCHING RECORDS OF DNA SEQUENCES OF ANOPHERES -------

# Let's search for DNA sequences from nucleotide databse (nuccore) for "anopheles" as the organism and "COI" as important gene for our prediction. The argument "use_history" is used to extract a large amount of hits from NCBI. We use entrez_search(), a function from rentrez R package which will allow to make a query to NCBI Entrez database (nuccore). 
df_anopheles_search = entrez_search(db = "nuccore", 
                                    term = "anopheles [ORGN] AND COI [gene]",
                                    use_history = T)
df_anopheles_search
#==> The output is 10730 hits. These are available results in NCBI corresponding to our query.

#Let's use readDNAStringSet() function from Biostrings package to read DNA sequences from a fasta file (format) and store them in DNAStrinngSet object, which is more efficient for large dataset.
st_anopheles = readDNAStringSet("../data/sequence.fasta")
st_anopheles 

# _ PART 3: DATA PARSING AND CLEANING ---------

# Let's convert DNAStringSet format into dataframe format. The data.frame() function is used to convert the nucleotide sequences into character strings that can be easly stored in dataframe format as plain text. DNAStringSet is an S4 object that can't help see DNA sequences as plain text.
df_anopheles = data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
head(df_anopheles) # Display the first DNA sequences in the console

# Let's convert back our dataframe into DNAStringSet. We do in place conversion because we already checked out the function and R needs the column "nucleotides2" to be a DNAStringSet, which allows Biostring function to recognize and analyse. This is done just for immediate check of the class of the column, nucleotides2. In place the DNAStringSet try to create an object containing all sequences at once, and not one per row. This will not change the class of df_anopheles.
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2)
head(df_anopheles)
#==> R stores the whole DNAStringSet object in a single complex entity (S4 object) and the remaining rows are filled with NA.

# Let's calculate the frequency of dinucleotides (2-mers). cbind() will allow to append new columns to the existing data frame, df_anopheles. First R calculate the dinucleotide frequencies for each DNA sequence, then convert the result into a data frame. Finally, it combines those newly created columns with the existing data frame, df_anopheles. The argument "as.prob"  returns frequencies as relative proportion. 
df_anopheles = cbind(df_anopheles,
                     as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE)))
head(df_anopheles)
#==> From 3 columns to 19 columns. Additional 16 variables come from 4 nucleotides and 2-mers, which make 4^2 = 16 (dinucleotides' proportions).

# There are two similar observations in columns "sequences" and "nucleotides2". The former has been manupilated to culculate dincleotides, we choose to cancel it and remain with "sequences" column. Then convert our object to tibble for tidyverse friendly handling.
df_anopheles = df_anopheles %>% 
  select(-nucleotides2) %>% 
  as_tibble()
#==> Now, we remain with 18 variables.

# The column "title" is too complex, containing different information, including unique id (1L) and species names (2L and 3L). Let's transform our tibble by creating new columns named "species_name" and "unique_id" using mutate() function. word() function will help access "title" column and extract specific words (e.g., 2L to mean 2nd word ...). "species_name" will contain the 2nd and 3rd words, while "unique_id" will contain the 1st word.
df_anopheles = df_anopheles %>% 
  mutate(species_name = word(title, 2L, 3L)) %>% 
  mutate(unique_id = word(title, 1L)) 
df_anopheles
#==> Two new columns ("species_name" and "unique_id") are appended to the tibble.

# Let's check the length or number of rows with unique id (using tidyverse). This will help know if the observations have distinct unique_id.
df_anopheles$unique_id %>% 
  unique() %>% 
  length()
#==> This returns 10730 unique id, the same number of hits we started with.

# Create a vector containing species_name (group_by()) and sort those species names from most to least frequent. Then maintain those having more than 200 nucleotides.
ls_vectors = df_anopheles %>% 
  group_by(species_name) %>% 
  count() %>% 
  arrange(desc(n)) %>% #View()
  filter(n > 200) %>% 
  pull(species_name)
ls_vectors
#==> The output is 14 species

# Create a new tibble containing only rows (in df_anopheles) whose  species_name appears in the vector created above, ls_vectors
df_anopheles_vector = 
  df_anopheles %>% 
  filter(species_name %in% ls_vectors)
#==> Only 4778 observations remain. These observations have distinct id and are grouped by species names. Each observation contains a sequence of more than 200 nucleotides.

# Still using tidyverse, let's find the records (rows) for each unique species.
df_anopheles_vector %>% 
  group_by(species_name) %>% 
  count()
#==> Records range from 609 to 220

# Let's remove all objects that are not needed in subsequent lines
rm(df_anopheles, st_anopheles)

# _ PART 4: TRAINING THE CLASSIFICATION MODEL ------

# To create the model, we need to know the lowest records across species. These records will serve as sample size for every species to avoid biasing the model toward the most represented species. Let's find the least records among all species using a frequency table. 
smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample
#==> A species with least records has 220 rows. The smallest group contains a number of records that can be considered as sample size for each group.

# Let's create a validation dataset, which is separate from training data set in order to learn real pattern in the data. This dataset will help evaluate externally how well the model generalizes on unseen data. We first set seed to keep our results reproducible. Then, we group species names and randomly sample 20% of the sample size for each species name.
set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation
#==> R returns 616 observations (check in the environment). This was calculated from: 220 rows (number of records in each species) * 14 (number of species) * 20% or simply 220 rows (the smallest record) * 20% = 44 rows from each species name (see the next code).

# A validation dataset has been created. The function group_by() will group species by species_name, then count() function  will count the number of records (rows) per each species_name in the validation dataset
df_anopheles_validation %>% 
  group_by(species_name) %>% 
  count()
#==> 44 rows are randomly sampled from each species to make a validation dataset of 616 rows (observations) in total.

# Let's create a training dataset with the remaing rows to the least frequent species (80%). Here we need to exclude rows of species_name already composing the validation dataset. Thus, "%in%" operator will allow to find any species_name appearing in the validation dataset, then "!" logical operator will negates the result from "%in%" to keep records which are not figuring in the validation dataset. ceiling() function will round up the resulting numerical value. Training dataset will teach the model how to recognize the patterns or relationships between features (dinucleotides) and response variable (species_nname) and make predictions on new data.
set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

# let's check records (rows) for each species in the training dataset
df_anopheles_training %>% 
  group_by(species_name) %>% 
  count()
#==> Exacly 176 rows per each species to make 2464 observations in total.

# Let's check the indexes of the 20 variables. These indexes will be used in randomForest() to locate our features.
names(df_anopheles_training)

# Now, let's build a classifier to species using dinucleotideFrequencies from column 3 to 18 as features and species_name as response or target variable. The latter are first converted into categorical factors to allow the randomForest build decision trees on separate classes. The randomForest() function allows to train a model by identifying features that best predict the species_name.  Argument "ntree" help define the number of trees to be generated. Here, the argument "importance" measures how each feature contributes to the classification.
anopheles_class <- randomForest(x = df_anopheles_training[, 3:18],
                                y = as.factor(df_anopheles_training$species_name), 
                                ntree = 100, importance = TRUE)

# Let's check the results of our classification model
anopheles_class
plot(anopheles_class)
varImpPlot(anopheles_class)
#==> The results show split: 4. This is the number of features used in the construction of each tree. 4 is a default calculated as the square root of the number of features (16). The randomForest validates itself internally using OOB estimate error rate. During training, a bootstrap sample (observations selected randomly with replacement from the training dataset) is created, and the samples left out are referred to as out-of-box observations. In internal validation, trees that have never seen that OOB sample predict the label (species_name) of those observations. In other words, the OOB estimate error rate of 1.7% indicates that the model classifies correctly 98.3% of samples it has not seen during training. Plot() shows that the OOB errors reduce drastically as the number of trees increases and start flattening at around 20 trees. This means that the model has converged and no additional trees can significantly improve accuracy. The varImpPlot of our classifier shows the importance of features in the classification model. In general, TT and CA are the most important predictors of our model, as removing any of them from the model can cause the greatest loss in accuracy of more than 23 and the largest drop in GINI impurity of more than 200 (in making pure splits in the decision tree) (refer to https://doi.org/10.1890/07-0539.1 for the interpretation of MeanDecreaseAccuracy).

# _ PART 5: VERIFICATION OF MODEL GENERALIZATION -----
# By using our validation dataset, we can check the generalization of our model outside training dataset.
predict_validation <- predict(anopheles_class, 
                              df_anopheles_validation[, c(19, 3:18)])
predict_validation

# Let's build our own confusion matrix. This matrix will show how samples in validation dataset were predicted using the model already created.
table(observed = df_anopheles_validation$species_name, 
      predicted = predict_validation)
#==> Results show that the least predicted were correctly classified in Anopheles maculipennis at 93.2%, while the model managed to correctly classify samples in 12 species names at 100%. This results confirm the generalization of our model.

### Thank you for this helpful exercise, which helped me learn alots about machine learning.