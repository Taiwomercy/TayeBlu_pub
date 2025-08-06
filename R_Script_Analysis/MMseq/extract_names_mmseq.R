setwd("") #set working directory

# Read the MMseqs2 clustering results
mmseqs <- read.table("mmseq_unchar_seq_0.3_cov_0.8.tsv", 
                     sep = "\t", 
                     header = FALSE, 
                     stringsAsFactors = FALSE)

# Assign column names (typical MMseqs2 output format)
colnames(mmseqs) <- c("Cluster_Rep", "Cluster_Member")

# Extract phage names from protein identifiers
# For IDs like "AZP_TayeBlu_0006", we want "AZP_TayeBlu"
# For IDs like "NC_071011_CDS_0038", we want "NC_071011_CDS"
mmseqs$PhageName <- sub("_[^_]+$", "", mmseqs$Cluster_Member)

# Get all unique phages and clusters
all_phages <- unique(mmseqs$PhageName)
all_clusters <- unique(mmseqs$Cluster_Rep)

# Create a list of phages for each cluster
cluster_phages <- list()
for (cluster in all_clusters) {
  cluster_members <- mmseqs$Cluster_Member[mmseqs$Cluster_Rep == cluster]
  cluster_phages[[cluster]] <- unique(sub("_[^_]+$", "", cluster_members))
}

# Create a reorganized dataframe with TayeBlu as the representative when present
reorganized_clusters <- data.frame(
  Cluster_ID = integer(),
  Original_Rep = character(),
  New_Rep = character(),
  Phage_Count = integer(),
  Members = character(),
  stringsAsFactors = FALSE
)

# Process each cluster
cluster_id <- 1
for (cluster in all_clusters) {
  # Get members directly
  members <- mmseqs$Cluster_Member[mmseqs$Cluster_Rep == cluster]
  phage_count <- length(cluster_phages[[cluster]])
  
  # Check if TayeBlu is present in this cluster
  # FIXED: Changed from "^TayeBlu_" to "TayeBlu" to match "AZP_TayeBlu_"
  tayeblu_members <- members[grepl("TayeBlu", members)]
  
  if (length(tayeblu_members) > 0) {
    # Use the first TayeBlu member as the representative
    new_rep <- tayeblu_members[1]
  } else {
    # Keep the original representative
    new_rep <- cluster
  }
  
  # Count unique proteins per phage in this cluster
  proteins_per_phage <- table(sub("_[^_]+$", "", members))
  proteins_per_phage_str <- paste(names(proteins_per_phage), proteins_per_phage, sep=":", collapse="|")
  
  # Add to the reorganized dataframe
  reorganized_clusters <- rbind(reorganized_clusters, data.frame(
    Cluster_ID = cluster_id,
    Original_Rep = cluster,
    New_Rep = new_rep,
    Phage_Count = phage_count,
    Members = paste(members, collapse = "|"),
    Proteins_Per_Phage = proteins_per_phage_str,
    stringsAsFactors = FALSE
  ))
  
  cluster_id <- cluster_id + 1
}

# Sort by phage count (descending) and then by whether TayeBlu is the representative
reorganized_clusters <- reorganized_clusters[order(-reorganized_clusters$Phage_Count, 
                                                   grepl("TayeBlu", reorganized_clusters$New_Rep), 
                                                   decreasing = TRUE), ]

# Write the reorganized clusters to a TSV file
write.table(reorganized_clusters, "reorganized_phage_clusters.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Create a simplified version with only TayeBlu representatives
tayeblu_clusters <- reorganized_clusters[grepl("TayeBlu", reorganized_clusters$New_Rep), ]
tayeblu_simplified <- data.frame(
  New_Rep = tayeblu_clusters$New_Rep,
  Phage_Count = tayeblu_clusters$Phage_Count,
  Members = tayeblu_clusters$Members,
  stringsAsFactors = FALSE
)

# Write the TayeBlu clusters to a TSV file
write.table(tayeblu_simplified, "tayeblu_clusters.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Print summary statistics
cat("\nNumber of clusters:", nrow(reorganized_clusters), "\n")
cat("Clusters with TayeBlu representatives:", nrow(tayeblu_clusters), "\n")
cat("Clusters with all phages (if 9 phages):", sum(reorganized_clusters$Phage_Count == 9), "\n")

# Calculate total unique proteins per phage across all clusters
all_proteins_by_phage <- list()
for (phage in all_phages) {
  # Get all proteins for this phage
  phage_proteins <- mmseqs$Cluster_Member[mmseqs$PhageName == phage]
  all_proteins_by_phage[[phage]] <- length(unique(phage_proteins))
}

cat("\nTotal unique proteins per phage:\n")
for (phage in names(all_proteins_by_phage)) {
  cat(phage, ":", all_proteins_by_phage[[phage]], "\n")
}

