# Title: Anopheles COI classification pipeline — Quiz 2
# Name: Yazan Alqawlaq
# Student number: 1021214

# ---- Packages -----
# What: load the tools we use across the 3 scripts.
# Why: each package handles a piece of the workflow (wrangle data, read FASTA, model).
library(tidyverse)    
library(stringr)      
library(Biostrings)  
library(randomForest) 
library(rentrez)
library(dplyr)

# ---- NCBI search (context/provenance) ----
# What: show how we'd find Anopheles COI records on NCBI (hits, web_history).
# Why: documents data provenance
df_anopheles_search = entrez_search(db = "nuccore", 
                                    term = "anopheles [ORGN] AND COI [gene]",
                                    use_history = T)
df_anopheles_search    

# ---- Read input sequences (FASTA → DNAStringSet) ----
# What: load the COI sequences from disk into R.
# Why: this is the raw data we’ll convert into numeric features for the model.
st_anopheles = readDNAStringSet("../data/sequence.fasta")
st_anopheles              

# ---- Build a basic data.frame (title + raw sequence) ----
# What: put the FASTA headers and sequences into a simple table.
# Why: gives us a tidy place to attach features and labels.
df_anopheles = data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
head(df_anopheles)

# ---- Prep DNA object for calculation ----
# What: add a DNAStringSet column from the raw sequence text.
# Why: dinucleotideFrequency() needs DNAStringSet input to compute features.
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2)
head(df_anopheles)   

# ---- Compute dinucleotide features ----
# What: turn each sequence into probabilities for all 16 dinucleotides (AA, AC, …, TT).
# Why: numeric features are what the classifier learns from.
df_anopheles = cbind(df_anopheles,
                     as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE)))
head(df_anopheles)

# ---- convert to tibble ----
# What: remove the heavy DNA object and keep a tidy tibble.
# Why: tibbles require simple vectors; avoids tibble/S4 issues later.
df_anopheles = df_anopheles %>% 
  select(-nucleotides2) %>% 
  as_tibble()

# ---- Parse labels from FASTA titles ----
# What: extract species_name (Genus species) and a unique_id from the FASTA header.
# Why: species_name is the target class; unique_id helps with bookkeeping.
df_anopheles = df_anopheles %>% 
  mutate(species_name = word(title, 2L, 3L)) %>% 
  mutate(unique_id = word(title, 1L)) 
df_anopheles         

# ---- ID check ----
# What: confirm how many unique records we have.
# Why: catches obvious duplication problems before modeling.
df_anopheles$unique_id %>% 
  unique() %>% 
  length()           

# ---- Keep species with enough data (n > 200) ----
# What: list species with > 200 sequences.
# Why: ensures each class has enough examples to train/evaluate fairly.
ls_vectors = df_anopheles %>% 
  group_by(species_name) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  filter(n > 200) %>% 
  pull(species_name)
ls_vectors            

# ---- Subset to well-represented species ----
# What: keep only rows from those species.
# Why: reduces class imbalance and focuses the model.
df_anopheles_vector = 
  df_anopheles %>% 
  filter(species_name %in% ls_vectors)

# ---- Check class sizes after filtering ----
# What: count rows per species in the working set.
# Why: verify we actually have multiple classes and enough samples.
df_anopheles_vector %>% 
  group_by(species_name) %>% 
  count()            

# ---- Clean up big objects we no longer need ----
# What: drop intermediates to free memory.
# Why: keeps the session lighter; good habit in longer pipelines.
rm(df_anopheles, st_anopheles)

# ---- Pick a balanced sample size ----
# What: use the smallest class size as our per-class cap.
# Why: makes splits balanced so the model isn’t biased by class size.
smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample        

# ---- Validation split (20% per class) ----
# What: hold out ~20% from each species for evaluation.
# Why: unbiased check of generalization after training.
set.seed(51)
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation  

# ---- Confirm validation balance ----
# What: count per class in the validation set.
# Why: ensure the holdout is balanced (fair comparison).
df_anopheles_validation %>% 
  group_by(species_name) %>% 
  count()             

# ---- Training split (80% per class) ----
# What: sample the training set from the remaining data (balanced).
# Why: provide the model with enough data from each species.
set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

# ---- Confirm training balance ----
# What: count per class in the training set.
# Why: check we’re actually training on balanced classes.
df_anopheles_training %>% 
  group_by(species_name) %>% 
  count()             

# ---- Peek at column names ----
# What: verify the feature columns match the indices used below.
# Why: prevents indexing mistakes when fitting/predicting.
names(df_anopheles_training)

# ---- Train Random Forest (100 trees) ----
# What: fit a classifier on dinucleotide features (cols 3:18).
# Why: RF is a strong baseline for tabular features; handles nonlinearity well.
anopheles_class <- randomForest(x = df_anopheles_training[, 3:18],
                                y = as.factor(df_anopheles_training$species_name), 
                                ntree = 100, importance = TRUE)

# ---- Model quick look & diagnostics ----
# What: print OOB error and show diagnostics/importance.
# Why: validation of model quality and which features mattered.
anopheles_class      
plot(anopheles_class)
varImpPlot(anopheles_class)

# ---- Predict on validation ----
# What: generate class predictions for the holdout set.
# Why: measure performance on unseen data.
predict_validation <- predict(anopheles_class, 
                              df_anopheles_validation[, c(19, 3:18)])
predict_validation    

# ---- Confusion table (observed vs predicted) ----
# What: compare true vs predicted labels.
# Why: see accuracy and any mix-ups between species.
table(observed = df_anopheles_validation$species_name, 
      predicted = predict_validation)   
