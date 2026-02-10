# Load required packages
library(ape)        # For reading GenBank files
library(HMM)        # For building simple HMMs
library(Biostrings) # For sequence manipulation
library(msa)        # For performing multiple sequence alignment
library(aphid)      # For building profile HMMs
library(igraph)     # For graph analysis
library(PRROC)      # For computing precision-recall curves and associated metrics

# if in doubt of the usage of a function, check documentation
# Hint: set working directory with setwd()
setwd("C:/Users/paolo/Desktop/Università/Magistrale/2° Anno/Secondo Semestre/High-throughput Bioinformatics/Lab/Exercise2/Exercise2")
getwd()

# Retrieve Olfactory receptor sequence from GenBank

# Read sequences, Hint: readAAStringSet
seq1 <- readAAStringSet("data/NP_001001957.2.faa")
seq2 <- readAAStringSet("data/NP_149420.4.faa")
seq3 <- readAAStringSet("data/WLF82657.1.faa")
seq1
seq2
seq3


# PART 1: Viterbi decoding---------------

# Define transition matrix
T_matrix <- matrix(c(0.8, 0.2, 0.05, 0.95), byrow=TRUE, 
            nrow=2, ncol=2, 
            dimnames=list(c("IN", "OUT"), c("IN", "OUT")))
T_matrix

# Define emission frequency matrix
in_freqs <- c(15, 11, 10, 9, 12, 8, 4, 12, 8, 29, 36, 8, 13, 24, 15, 34, 20, 1, 12, 31)
out_freqs <- c(10, 16, 12, 5, 17, 18, 4, 12, 11, 20, 45, 12, 17, 20, 19, 26, 20, 7, 6, 15)
aa_symbols <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

states <- c("IN", "OUT")
F_matrix <- matrix(c(in_freqs, out_freqs), nrow = 2, ncol = 20, byrow = TRUE, dimnames = list(states, aa_symbols))
F_matrix

# Convert F to an emission probability matrix, E
#Hint: colSums or rowSums()
row_totals <- rowSums(F_matrix)
E_matrix <- F_matrix / row_totals
E_matrix

# 1. Viterbi decoding algorithm from scratch
viterbi_decoding <- function(s, T_matrix, E_matrix){
	# Transform T and E to their log2 equivalents
	logT <- log2(T_matrix)
	logE <- log2(E_matrix)
	s <- (strsplit(as.character(s), "")[[1]])
	n <- length(s)
	num_states <- nrow(T_matrix)

	# Initialize viterbi and pointer matrices
	# of dimension: Num(states) x (n+1)
	V <- matrix(NA, nrow = num_states, ncol = n+1) # Log probability matrix
	# of dimension: Num(states) x (n)
	P <- matrix(NA, nrow = num_states, ncol = n) # Backpointer matrix
	
	# initialize the first column of V
	V[,1] <- log2(1/num_states)

	for (i in 2:(n+1)){
		for (l in 1:num_states){
			# populate V and P matrices using the recursive relation
		  log_probs <- V[, i-1] + logT[, l] + logE[l, match(s[i-1], colnames(logE))]
		  V[l,i] <- max(log_probs)
		  P[l,i-1] <- which.max(log_probs)
		}
	}

	# Use P to trace the most probable path
	path <- numeric(n)
	path[n] <- which.max(V[,n+1])
	for(i in (n-1):1){
	  path[i] <- P[path[i+1], i+1]
	}
	# Map numeric path to state labels
	states <- c("IN", "OUT")
	path <- states[path]
	                      
	#initialize a list or dictionary to store the results
	result <- list()
	result$V <- V # store viterbi matrix
	result$P <- P # store pointer matrix
	result$path <- path # store probable path
	result$log_prob <- max(V[,n+1]) # log probability of the probable path,

	return(result)
}

viterbi_result <- viterbi_decoding(seq1, T_matrix, E_matrix) #check result of the decoding
viterbi_result

# Piecewise plot Hint: ifelse() to map the states to 0/1
path_numeric <- ifelse(viterbi_result$path == "IN", 0, 1)

# plot(..., type="o", ...)
png("Figures/Most_probable_path.png", width = 800, height = 600)
plot(1:length(path_numeric), path_numeric, 
     type = "o",
     xlab = "Sequence Position", 
     ylab = "Hidden State",
     main = "Most Probable Hidden State Path",
     col = "darkblue",
     yaxt = "n")

# Changing the values on y-axis to state labels
axis(2, at = c(0, 1), labels = c("IN", "OUT"))

legend("bottomleft", 
       legend = c("IN (Hydrophobic)", "OUT (Hydrophilic)"),
       col = c("darkblue", "darkblue"))
dev.off()


# 2. Sequence scoring via viterbi decoding
score1 <- viterbi_decoding(seq1, T_matrix, E_matrix)$log_prob
score2 <- viterbi_decoding(seq2, T_matrix, E_matrix)$log_prob
score3 <- viterbi_decoding(seq3, T_matrix, E_matrix)$log_prob

scores <- c(score1, score2, score3)
names(scores) <- c("NP 001001957.2", "NP 149420.4", "WLF82657.1")
print(names(scores)[which.max(scores)])
print(max(scores))

sink("Figures/Most_likely_sequence_HMM.txt")
cat("Most probable sequence according to the 2-state HMM: \n")
cat("Most likely sequence: ", names(scores)[which.max(scores)], "\n")
cat("With a log-probability score of: ", max(scores))
sink()


# 3. HMM training

# initialize 2-state HMM
hmm <- initHMM(States = states, Symbols = aa_symbols, startProbs = c(0.5, 0.5),
               transProbs = T_matrix, emissionProbs = E_matrix)

# Conversion of the AAStringSet to a character vector
seq1_vector <- strsplit(as.character(seq1), "")[[1]]

trainedHMM <- baumWelch(hmm, seq1_vector, maxIterations=100, delta = 1e-9)$hmm
#or 
trainedHMM <- viterbiTraining(hmm, seq1_vector, maxIterations=100, delta = 1e-9)$hmm

# infer most likely hidden state sequence 
seq1_path <- viterbi(trainedHMM, seq1_vector)  
seq1_path

T_final <- trainedHMM$transProbs
sink("Figures/Trained_HMM_Transition_Matrix.txt")
cat("Transition Matrix initial:", T_matrix ,"\n")
cat("Transition Matrix final: ", T_final, "\n")
sink()

# plot hidden state sequence as piecewise 
path_numeric_training <- ifelse(seq1_path == "IN", 0, 1)
png("Figures/Most_probable_path_training.png", width = 800, height = 600)
plot(1:length(path_numeric_training), path_numeric_training, 
     type = "o",
     xlab = "Sequence Position", 
     ylab = "Hidden State",
     main = "Most Probable Hidden State Path",
     col = "darkblue",
     yaxt = "n")

# Changing the values on y-axis to state labels
axis(2, at = c(0, 1), labels = c("IN", "OUT"))

legend("bottomleft", 
       legend = c("IN (Hydrophobic)", "OUT (Hydrophilic)"),
       col = c("darkblue", "darkblue"))
dev.off()



# PART 2: Probabilistic modeling of protein families with profile HMMs---------------

# Read sequences, Hint: readAAStringSet
protein_family <- readAAStringSet("data/rats_olfactory_receptors.fasta")
protein_family

# Compute MSA
msa_result <- msa(protein_family, method = "ClustalW")
msa_result

# plot columns 50-165, Hint: msaPrettyPrint()
msaPrettyPrint(msa_result, 
               output="pdf", 
               file="MSA_segment.pdf", 
               showNames="none",    
               showLogo="none",
               y=c(50, 165))  


## Unfortunately, I was unable to experimentally validate this section of the analysis 
## due to hardware compatibility constraints. The aphid R package, which is essential 
## for deriving Profile HMMs, relies on native dependencies optimized for x86 architectures. 
## As my current workstation utilizes an ARM64 architecture (Windows on ARM), the installation 
## process consistently fails due to the absence of compatible pre-compiled binaries and the 
## inability to compile the source code for this specific instruction set

# Derive HMM based on given protein family
# convert msa to AAbin (binary data format)
msa_bin <- as(msa_result, "AAbin") 

profile_HMM <- derivePHMM(msa_result, residues="", pseudocount="Laplace", 
	alignment=TRUE)

# Score sequences
scoreHMM_1 <- Viterbi(profile_HMM, seq1)$score
scoreHMM_2 <- Viterbi(profile_HMM, seq2)$score
scoreHMM_3 <- Viterbi(profile_HMM, seq3)$score
scoresHMM <- c(scoreHMM_1, scoreHMM_2, scoreHMM_3)
names(scoresHMM) <- c("NP 001001957.2", "NP 149420.4", "WLF82657.1")
print(names(scoresHMM)[which.max(scoresHMM)])
print(max(scoresHMM))

sink("Figures/Most_likely_member_HMM.txt")
cat("Most probable member of the GPCR family: \n")
cat("Most likely member: ", names(scoresHMM)[which.max(scoresHMM)], "\n")
cat("With a log-probability score of: ", max(scoresHMM))
sink()



# PART 4: Post-hoc consistency correction of a model’s predictions---------------

# Define graph edges
edges <- c(
  "a","b",
  "a","c",
  "b","d",
  "b","e",
  "c","f",
  "c","g",
  "c","h",
  "f","i",
  "g","i"
)

# Create DAG graph, Hint: igraph::make_graph
dag <- graph(edges = edges, directed = TRUE)


# Enumerate all consistent subgraphs in dag

enumerate_subgraphs <- function(g){
	# Find all descendants of a given node (root node preferably)
	# Hint: igraph's subcomponent()
  in_degrees <- degree(g, mode="in")
  root_node <- V(g)[in_degrees == 0]$name
  
  # Get all descendants from the root (includes root itself)
  descendants <- V(g)[subcomponent(g, root_node, mode="out")]$name
  
	# initialize list to store consistent subgraphs found
  consistent_subgraphs <- list()
  
	# for each descendant
	for (d in 1:length(descendants)){
		# enumerate all combinations of nodes (i.e. random subgraphs)
		for (nodes in combn(descendants, d, simplify=FALSE)){
			# generate the subgraph induced by the set of nodes
			# Hint: igraph's subgraph
		  			subg <- induced_subgraph(g, vids = nodes)
			# node names can be accessed by: V(g)$name
		  			current_node_names <- V(subg)$name
		  			is_consistent <- TRUE
		  			for (n in current_node_names) {
		  			  parents <- neighbors(g, n, mode = "in")$name
		  			  if (length(parents) > 0 && !all(parents %in% current_node_names)) {
		  			    is_consistent <- FALSE
		  			    break
		  			  }
		  			}
		  			if (is_consistent) {
		  			  consistent_subgraphs[[length(consistent_subgraphs) + 1]] <- subg
			# Reflect on what a consistent graph is
			# Check if induced subgraph meets the consistency criteria

			# store the node set forming the consistent subgraph in the list 
		  			}
		}
	}
	
	return(consistent_subgraphs)
}

list_all_consistent_subgraphs <- enumerate_subgraphs(dag)

# number of consistent subgraphs
# number of subgraphs having 7 nodes
total_subgraphs <- length(list_all_consistent_subgraphs)
node_counts <- sapply(list_all_consistent_subgraphs, vcount)
subgraphs_7_nodes <- sum(node_counts == 7)
print(total_subgraphs)
subgraphs_7_nodes

sink("Figures/Consistent_subgraphs.txt")
cat("Consistent Subgraphs\n")
cat("Number of subgraphs having 7 nodes:", subgraphs_7_nodes, "\n\n")
cat("Number of consistent subgraphs:", total_subgraphs, "\n")
sink()

# Post-hoc consistency correction

consistency_correction <- function(g, node_scores){
	# Check lecture notes,
	# Useful function: topo_sort(), rev(), neighbors(), max()
  node_ordering <- topo_sort(g, mode = "out")
  reversed_ordering <- rev(node_ordering)
  corrected_scores <- node_scores
  
  for (v_idx in reversed_ordering) {
    v_name <- V(g)[v_idx]$name
    v_score <- corrected_scores[v_name]
    
    parents <- neighbors(g, v_name, mode = "in")$name
    
    for (p in parents) {
      prev_score <- corrected_scores[p]
      corrected_scores[p] <- max(prev_score, v_score)
    }
  }

	return(corrected_scores)
}

#Thresholding, Hint: sapply(),
initial_scores <- c(
  "a" = 0.55,
  "b" = 0.32,
  "c" = 0.95,
  "d" = 0.70,
  "e" = 0.49,
  "f" = 0.23,
  "g" = 0.63,
  "h" = 0.85,
  "i" = 0.40
)

corrected_scores <- consistency_correction(dag, initial_scores)

binary_predictions <- ifelse(corrected_scores >= 0.5, 1, 0)
nodes_in_corrected_subgraph <- names(binary_predictions[binary_predictions == 1])
corrected_scores
nodes_in_corrected_subgraph

sink("Figures/Corrected_Subgraphs_Scores.txt")
cat("Final Scores after Correction:\n")
cat(corrected_scores)
cat("\nNodes in Corrected Subgraph (Threshold >= 0.5):\n")
cat(nodes_in_corrected_subgraph, sep=", ")
sink()

#plot pr-curves
# Define groundtruth and predicted vectors
# corresponds to an alphabetical ordering of node labels/names
ytrue <- c(1, 0, 1, 0, 0, 1, 1, 0, 1); 
# assuring that the ordering of predicted scores matches the ordering of ytrue
ordered_names <- sort(names(initial_scores))
yhat_before_correction <- initial_scores[ordered_names]
yhat_after_correction  <- corrected_scores[ordered_names]


pr_curve1 <- pr.curve(scores.class0 = yhat_before_correction, 
                      weights.class0 = ytrue , 
                      curve = TRUE)


pr_curve2 <- pr.curve(scores.class0 = yhat_after_correction, 
                      weights.class0 = ytrue, 
                      curve = TRUE)

#retrieve AUC values:
auc_before_correction <- pr_curve1$auc.integral 
auc_after_correction  <- pr_curve2$auc.integral 
auc_before_correction
auc_after_correction

sink("Figures/AUC_Values.txt")
cat("AUC M1 (Original):", auc_before_correction, "\n")
cat("AUC M2 (Correct):", auc_after_correction, "\n")
sink()

png("Figures/Roc_curve.png", width = 800, height = 600)
# plot pr-curve 1
plot(pr_curve1, color = "red", auc.main = FALSE, main = "Precision-Recall Curves")

# add curve 2 onto the same plot with different line color
plot(pr_curve2, color = "blue", add = TRUE, auc.main = FALSE)

# use legend() to label them appropriately
legend("bottomleft", 
       legend = c(paste("M1 (Inital Scores) AUC:", round(auc_before_correction, 3)), 
                  paste("M2 (Corrected) AUC:", round(auc_after_correction, 3))), 
       col = c("red", "blue"), 
       lty = 1,
       cex = 0.8)
dev.off()


