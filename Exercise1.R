# Load required packages
library(Biostrings)  # For sequence manipulation
library(ORFik)       # For finding ORFs
library(ape)         # For reading GenBank files
library(pwalign)     # For sequence alignments

# Hint: set working directory with setwd()
setwd("C:/Users/paolo/Desktop/Università/Magistrale/2° Anno/Secondo Semestre/High-throughput Bioinformatics/Lab/Exercise1/Exercise1")
getwd()

# PART 1: Analyzing sequence composition and open reading frames ----------

# 1. GC skew --------------------------------------------------------------

# Hint: DNAString(), alphabetFrequency(), for loop

# Retrieve SARS-CoV-2 genome sequence from GenBank
sars_cov2_genome <- read.GenBank("NC_045512", as.character = T)

# Convert the sequence uppercase character vector
sars_cov2_sequence <- toupper(paste(sars_cov2_genome[[1]], collapse = ""))

# Display the sequence
sars_cov2_sequence
sars_cov2_sequence <- strsplit(sars_cov2_sequence, "")[[1]]
sars_cov2_sequence
n <- length(sars_cov2_sequence) 
n

# Create windows and vectors to store results
window_size <- 500
step_size <- 50
GC_skew_values <- c()
base_position <- c()

# For loop to compute GC skew values
for (start in seq(1, n - window_size + 1, by = step_size)) {
  end <- start + window_size - 1
  window <- sars_cov2_sequence[start:end]
  
  g <- sum(window == "G")
  c <- sum(window == "C")
  
  GC_skew <- (g - c) / (g + c)
  
  base_position <- c(base_position, start)
  GC_skew_values <- c(GC_skew_values, GC_skew)
}

# Visualize first 10 values
head(base_position, 10)
head(GC_skew_values, 10)

# Plot
png("Figures/GC_skew_plot.png", width=800, height=600)
plot(base_position, GC_skew_values, type = "l",
     xlab = "Genomic Position (bp)", ylab = "GC Skew",
     main = "GC Skew of SARS-CoV-2 Genome",
     col = "blue", lwd = 2)
dev.off()


# 2. Identify ORFs --------------------------------------------------------

# Hint: findORFs(), which.max(), as.data.frame() for easier inspection of results

# Convert the sequence into a DNAString object
sars_cov2_dna <- DNAStringSet(paste(sars_cov2_sequence, collapse = ""))
sars_cov2_dna

# Extraction of all the codons from the sequence
orfs <- findORFs(sars_cov2_dna, startCodon = "ATG", stopCodon = "TAA|TAG|TGA", longestORF = FALSE)
orfs

# Convert ORFs to data frame
orfs_df <- as.data.frame(orfs)
head(orfs_df)

# Identify longest ORF
longest_orf <- orfs_df[which.max(orfs_df$width), ]
longest_orf

sink("Figures/longest_ORF.txt")
cat("Longest ORF in SARS-CoV-2 genome:\n")
cat("Start position: ", longest_orf$start, "\n")
cat("End position:   ", longest_orf$end, "\n")
cat("Length:         ", longest_orf$width, "bp\n")
sink()


# 3. Permutation test -----------------------------------------------------

# Define a function for the permutation test
perform_permutation_test <- function(sequence, num_iterations = 1000) {
  # Initialize a vector to store randomized ORF lengths
  randomized_ORF_lengths <- numeric()
  
  # Iterate over num_iterations, for loop
  for (i in 1:num_iterations) {
    set.seed(i)
    # In each iteration: 
    # Shuffle the sequence, useful R functions: strsplit(), sample(), paste()
    shuffled_seq <- sample(sequence)
    
    # Find ORFs in the randomized sequence and record their lengths
    shuffled_seq_dna <- DNAStringSet(paste(shuffled_seq, collapse = ""))
    shuffled_orfs <- findORFs(shuffled_seq_dna, startCodon = "ATG", stopCodon = "TAA|TAG|TGA", longestORF = FALSE)
    randomized_ORF_lengths <- c(randomized_ORF_lengths, width(unlist(shuffled_orfs))) #unlist() takes a list in R and turns it into a “flat” vector
    
  }
  # Return randomized ORF lengths
  return(randomized_ORF_lengths)
}

randomized_ORF_lengths <- perform_permutation_test(sars_cov2_sequence, 1000)
randomized_ORF_lengths

# Loop over the ORFs identified in the previous part and compute p-values based
# on the ORF lengths

iterations <- length(orfs_df$width)
iterations
p_values <- numeric(iterations)
for (i in 1:iterations) {
  p_values[i] <- sum(randomized_ORF_lengths >= orfs_df$width[i]) / length(randomized_ORF_lengths)
}
p_values

# Count how many ORFs have p-value < 0.01
num_significant_orfs <- sum(p_values < 0.01)
num_significant_orfs

sink("Figures/Number_significant_ORF.txt")
cat("Number of significant orfs:",  num_significant_orfs)
sink()

png("Figures/ORF_length_distribution.png", width = 800, height = 600)
hist(
  randomized_ORF_lengths,
  breaks = 50,
  prob = TRUE,
  col = "lightblue",
  main = "Distribution of ORF lengths: Randomized vs Observed",
  xlab = "ORF length",
  ylab = "Frequency"
)

# Overlay observed ORF lengths (use abline())
abline(v = orfs_df$width, col = "red", lwd = 0.5)

legend(
  "topright",
  legend = c("Randomized ORFs", "Observed ORFs"),
  pch = c(15, NA), 
  lty = c(NA, 1),
  col = c("lightblue", "red"),
  lwd = c(NA, 0.5),
  pt.cex = 1.5,     
  cex = 0.8,         
  bty = "o"         
)
dev.off()

# 4. Multiple testing correction ------------------------------------------

# Hint: Use p.adjust() for multiple testing correction (e.g., Bonferroni, BH).

# Apply Bonferroni correction
p_values_bonferroni <- p.adjust(p_values, method = "bonferroni")
p_values_bonferroni

# Apply Benjamini-Hochberg (BH) correction
p_values_bh <- p.adjust(p_values, method = "BH")
p_values_bh

# Significant ORF after correction
num_significant_bonferroni <- sum(p_values_bonferroni < 0.01)
num_significant_bh <- sum(p_values_bh < 0.01)
num_significant_bonferroni
num_significant_bh


sink("Figures/multiple_testing_results.txt")
cat("Multiple testing correction results:\n")
cat("Bonferroni correction:\n")
cat("Number of significant ORFs (p < 0.01): ", num_significant_bonferroni, "\n")
cat("Benjamini-Hochberg correction:\n")
cat("Number of significant ORFs (p < 0.01): ", num_significant_bh, "\n")
sink()

# PART 2: Identifying the spike protein encoding ORF ----------------------

# 5. Train Markov models on known spike sequences -------------------------

# Load spike protein sequences
spike_sequences <- readDNAStringSet("data/coronaviruses_spike.fasta")
spike_sequences

# Convert to character vector
spike_sequences <- as.character(spike_sequences)
spike_sequences

# Complete the function train a first-order markov model 
train_markov_model <- function(sequences) {
  
  # Initialize transition matrix
  transition_matrix <- matrix(0, nrow = 4, ncol = 4, 
                              dimnames = list(c("A", "T", "G", "C"), 
                                              c("A", "T", "G", "C")))
  
  
  # Count transitions across all sequences
  for (sequence in sequences) {
    for (i in 1:(nchar(sequence) - 1)) {
      
      # Hint: Use substring() to extract current and next nucleotides and update
      # the transition matrix
      current_nuc <- substring(sequence, i, i)
      next_nuc <- substring(sequence, i + 1, i + 1)
      #adding +1 to the corresponding cell in the matrix
      transition_matrix[current_nuc, next_nuc] <- transition_matrix[current_nuc, next_nuc] + 1
      
    }
  }
  
  
  # Remember to normalize the transition matrix by rows
  transition_matrix <- transition_matrix / rowSums(transition_matrix)
  
  return(transition_matrix) 
}

markov_model <- train_markov_model(spike_sequences)
markov_model


# 6. Compute the likelihoods for each candidate ORF -----------------------

# Function to calculate likelihoods
calculate_likelihood <- function(sequence, model) {
  
  # Initialize log-likelihood; assuming the initial state is given, i.e. the
  # initial probability is 1
  log_likelihood <- log(1)
  
  # Traverse the sequence and calculate the likelihood using the given model
  for (i in 1:(nchar(sequence) - 1)) {
    
    # Hint: Use substring() to extract current and following nucleotides and
    # update the likelihood based on the transition matrix
    current_nuc <- substring(sequence, i, i)
    next_nuc <- substring(sequence, i + 1, i + 1)
    prob <- model[current_nuc, next_nuc]
    log_likelihood <- log_likelihood + log(prob)
    
  }
  
  # Hint: remember to normalize the likelihood with the sequence length
  normalized_log_likelihood <- log_likelihood / (nchar(sequence) - 1)
  
  return(normalized_log_likelihood)
}


# Loop through each significant ORF and compute their log-likelihood
# Get all the significant ORFs after Bonferroni correction
significant_orfs <- orfs_df[p_values_bonferroni < 0.01, ]
significant_orfs

# Extract the sequences of ORF from dataframe
log_likelihoods <- numeric(nrow(significant_orfs))

for (i in 1:nrow(significant_orfs)) {
  orf <- significant_orfs[i,] # extract each row to get start and end
  sequence <- paste(sars_cov2_sequence[orf$start:orf$end], collapse = "")
  log_likelihoods[i] <- calculate_likelihood(sequence, markov_model)
} 

log_likelihoods
index <- which.max(log_likelihoods)
index
best_orf <- significant_orfs[index, ]
best_orf
sequence_best_orf <- paste(sars_cov2_sequence[best_orf$start:best_orf$end], collapse = "")
sequence_best_orf

sink("Figures/most_likely_ORF.txt")
cat("Most likely sequence to be the one coding for the spike protein in SARS-CoV-2: \n")
cat("Start position:", best_orf$start, "\n")
cat("End position:", best_orf$end, "\n")
cat("Length:", best_orf$width, "\n")
cat("Normalized Log-Likelihood:", log_likelihoods[index], "\n")
cat("Nucleotide sequence:", sequence_best_orf)
sink()

# PART 3: Sequence alignment ----------------------------------------------

# Load nucleotide sequences coding for spike proteins from other coronaviruses
spike_sequences <- readDNAStringSet("data/coronaviruses_spike.fasta")
spike_sequences

# Load spike protein coding sequence from SARS-CoV-2
sars_cov2_spike_sequence <- readDNAStringSet("data/SARS-CoV-2_spike.fasta")
sars_cov2_spike_sequence

# 7. Translating the nucleotide sequences into amino acid sequences -------
# Hint: translate()
spike_sequences_translated <- translate(spike_sequences)
spike_sequences_translated

for (i in 1:length(spike_sequences_translated)) {
  seq <- as.character(spike_sequences_translated[i])
  
  if (substr(seq, 1, 1) == "M") {
  } else {
    stop(paste("Sequence", names(spike_sequences_translated)[i], "does not start with Methionine"))
  }
}

# 8. Global alignment -----------------------------------------------------
# Hint: pairwiseAlignment()
sars_cov2_spike_sequence_translated <- translate(sars_cov2_spike_sequence)
sars_cov2_spike_sequence_translated

gap_opening_global <- 10 # Cost to start a new gap. A high penalty discourages opening too many gaps.
gap_extension_global <- 2 # Additional cost to extend an already opened gap A lower penalty allows long gaps without making the score too negative.

# Substitution matrix
data(BLOSUM62)

global_alignment_scores <- numeric(length(sars_cov2_spike_sequence_translated))

for (i in 1:length(spike_sequences_translated)) {
  target = spike_sequences_translated[i]
  
  alignment <- pairwiseAlignment(
    pattern = sars_cov2_spike_sequence_translated[[1]], # Reference sequence
    subject = target, # Other spike protein sequence from other coronavirus we want to compare to the refence sequence
    substitutionMatrix = BLOSUM62,
    gapOpening = gap_opening_global,
    gapExtension = gap_extension_global,
    type = "global" 
  )
  global_alignment_scores[i] <- score(alignment)
}

global_alignment_scores  
global_index_scores <- which.max(global_alignment_scores)
global_index_scores

sink("Figures/most_similar_spike_global.txt")
cat("Human coronavirus most similar globally to SARS-CoV-2 spike protein:\n")
cat("Name:", names(spike_sequences_translated[global_index_scores]), "\n")
cat("Global alignment score:", global_alignment_scores[global_index_scores], "\n")
cat("Amino acid sequence:", as.character(spike_sequences_translated[[global_index_scores]]),"\n")
sink()

# 9. Local alignment ------------------------------------------------------
# Hint: pairwiseAlignment()
gap_opening_local <- 5 
gap_extension_local <- 1 

local_alignment_scores <- numeric(length(sars_cov2_spike_sequence_translated))

for (i in 1:length(spike_sequences_translated)) {
  target = spike_sequences_translated[i]
  
  alignment <- pairwiseAlignment(
    pattern = sars_cov2_spike_sequence_translated[[1]], 
    subject = target, 
    substitutionMatrix = BLOSUM62,
    gapOpening = gap_opening_local,
    gapExtension = gap_extension_local,
    type = "local" 
  )
  local_alignment_scores[i] <- score(alignment)
}

local_alignment_scores  
local_index_scores <- which.max(local_alignment_scores)
local_index_scores

sink("Figures/most_similar_spike_local.txt")
cat("Human coronavirus most similar to SARS-CoV-2 spike protein:\n")
cat("Name:", names(spike_sequences_translated[local_index_scores]), "\n")
cat("Local alignment score:", local_alignment_scores[local_index_scores], "\n")
cat("Amino acid sequence:", as.character(spike_sequences_translated[[local_index_scores]]),"\n")
sink()

# 10. Blast results ------------------------------------------------------
# Export the translated SARS-CoV-2 spike protein sequence in FASTA format
writeXStringSet(
  sars_cov2_spike_sequence_translated,
  "data/SARS-CoV-2_spike_protein.fasta")

sink("Figures/BLAST_bat_coronavirus_results.txt")
cat("BLASTp analysis of SARS-CoV-2 spike protein against bat coronaviruses\n")

cat("Top BLAST hit:\n")
cat("Description: Spike glycoprotein\n")
cat("Accession code: UAY13217.1\n")

cat("Alignment score:\n")
cat("E-value: 0.0\n")
cat("Percentage identity: 98.43%\n")
cat("Host species: Rhinolophus malayanus\n\n")

sink()




