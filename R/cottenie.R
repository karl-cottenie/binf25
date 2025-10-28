library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflicts_prefer(dplyr::desc)
conflicted::conflicts_prefer(base::as.factor)
library(rentrez)
library(Biostrings)
library(randomForest)
source("Entrez_Functions.R")

df_anopheles_search = entrez_search(db = "nuccore", 
                                    term = "anopheles [ORGN] AND COI [gene]",
                                    use_history = T)
df_anopheles_search

FetchFastaFiles(searchTerm = "anopheles [ORGN] AND COI [gene]", 
                                     fastaFileName = "../data/dl_anopheles")
df_anopheles = MergeFastaFiles(filePattern = "../data/dl_anopheles*")
df_anopheles
str(df_anopheles)
files_to_move = list.files(".", pattern = "\\.fasta$", full.names = TRUE)
file.rename(from = files_to_move,
            to = file.path("../data/", 
                           basename(files_to_move)))
names(df_anopheles) = c("title", "sequence") #jump to line w/nucleotides2

st_anopheles = readDNAStringSet("../data/sequence.fasta")
st_anopheles 

df_anopheles = data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
head(df_anopheles)

df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2)
head(df_anopheles)

df_anopheles = cbind(df_anopheles,
                     as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE)))
head(df_anopheles)

df_anopheles = df_anopheles %>% 
  select(-nucleotides2) %>% 
  as_tibble()

df_anopheles = df_anopheles %>% 
  mutate(species_name = word(title, 2L, 3L)) %>% 
  mutate(unique_id = word(title, 1L)) 
df_anopheles

df_anopheles$unique_id %>% 
  unique() %>% 
  length()

ls_vectors = df_anopheles %>% 
  group_by(species_name) %>% 
  count() %>% 
  arrange(desc(n)) %>% #View()
  filter(n > 200) %>% 
  pull(species_name)
ls_vectors

df_anopheles_vector = 
  df_anopheles %>% 
  filter(species_name %in% ls_vectors)
  
df_anopheles_vector %>% 
  group_by(species_name) %>% 
  count()

rm(df_anopheles, st_anopheles)

smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample

set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation

df_anopheles_validation %>% 
  group_by(species_name) %>% 
  count()

set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

df_anopheles_training %>% 
  group_by(species_name) %>% 
  count()

names(df_anopheles_training)

anopheles_class <- randomForest(x = df_anopheles_training[, 3:18],
                                y = as.factor(df_anopheles_training$species_name), 
                                ntree = 100, importance = TRUE)

anopheles_class
plot(anopheles_class)
varImpPlot(anopheles_class)

predict_validation <- predict(anopheles_class, 
                              df_anopheles_validation[, c(19, 3:18)])
predict_validation

table(observed = df_anopheles_validation$species_name, 
      predicted = predict_validation)

