##Required libraries -----
library(rentrez)
library(Biostrings)
library(dplyr)
library(stringr)
library(randomForest)
library(tibble)

##Query NCBI for Anopheles COI sequences -----
df_anopheles_search = entrez_search(db = "nuccore",
                                    term = "anopheles [ORGN] AND COI [gene]",
                                    use_history = T)
df_anopheles_search

##Import DNA sequences from FASTA dataset -----
st_anopheles = readDNAStringSet("../data/sequence.fasta")
st_anopheles 

##Create data frame from sequence titles and bases -----
df_anopheles = data.frame(title = names(st_anopheles),
                          sequence = paste(st_anopheles))
head(df_anopheles)

##Convert raw sequences to DNAStringSet objects -----
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2)
head(df_anopheles)

##Find dinucleotide composition frequencies -----
df_anopheles = cbind(df_anopheles,
                     as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2,
                                                         as.prob = TRUE)))
head(df_anopheles)

##Remove DNA objects and convert to tibble -----
df_anopheles = df_anopheles %>%
  select(-nucleotides2) %>%
  as_tibble()

#Extracts species names and unique identifiers from FASTA titles into new columns.
df_anopheles = df_anopheles %>%
  mutate(species_name = word(title, 2L, 3L)) %>%
  mutate(unique_id = word(title, 1L))
df_anopheles

##Count unique sequence IDs -----
df_anopheles$unique_id %>%
  unique() %>%
  length()

##Summarize sequence counts per species (>200 only) -----
ls_vectors = df_anopheles %>%
  group_by(species_name) %>%
  count() %>%
  arrange(desc(n)) %>%
  filter(n > 200) %>%
  pull(species_name)
ls_vectors

##Filter dataset to include selected vector species -----
df_anopheles_vector =
  df_anopheles %>%
  filter(species_name %in% ls_vectors)

##Check counts per species in filtered dataset -----
df_anopheles_vector %>%
  group_by(species_name) %>%
  count()

##Remove un-needed objects -----
rm(df_anopheles, st_anopheles)

##Identify smallest species group size -----
smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample

##Split data into validation and training sets -----
set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation

#Count validation sequences per species 
df_anopheles_validation %>%
  group_by(species_name) %>%
  count()

set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

##Count training sequences per species -----
df_anopheles_training %>%
  group_by(species_name) %>%
  count()

##Inspect training dataframe structure -----
names(df_anopheles_training)

##Train Random Forest classifier model -----
anopheles_class <- randomForest(x = df_anopheles_training[, 3:18],
                                y = as.factor(df_anopheles_training$species_name),
                                ntree = 100,
                                importance = TRUE)

#Examine trained Random Forest model and visualize variable importance 
anopheles_class
plot(anopheles_class)
varImpPlot(anopheles_class)

##Validate model predictions -----
predict_validation <- predict(anopheles_class,
                              df_anopheles_validation[, c(19, 3:18)])
predict_validation

##Compare predicted vs observed species classifications -----
table(observed = df_anopheles_validation$species_name,
      predicted = predict_validation)
