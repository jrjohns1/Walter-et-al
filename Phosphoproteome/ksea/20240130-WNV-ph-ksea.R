library(ksea)
library(reshape2)
library(ggplot2)
library(stringr)
library(readxl)



# Set working directory
setwd("~/OneDrive - The Mount Sinai Hospital/Data/Collaborations/RamageLab/20240130-WNV-reanalysis/20240130-WNV-ph/ksea/")


# Set prefix for naming output files
prefix <- "20240130-WNV-ph"



# Set database files
dbfiles <- c("~/OneDrive - The Mount Sinai Hospital/DB/20191017.SwissProt.Hsapiens.WNVNY99.fasta")



# Get the number of databases
numdb <- length(dbfiles)

# Loop through each database file
for (i in 1:numdb) {
  # If it's the first database then read lines into dblines
  if (i == 1) {
    dblines <- readLines(con = dbfiles[i])
  } else {
    # If it's not the first database then read lines into tmp
    # Then concatenate tmp onto the end of dblines
    tmp <- readLines(con = dbfiles[i])
    dblines <- c(dblines, tmp)
  }
}

# Grab header lines starting with > character
headerlines <- dblines[grep(">", dblines)]

# Extract protein accession and protein name from header lines
protein.accession <- str_match(headerlines, ">\\S\\S\\|(\\S+)\\|(\\S+.*)")[,2]
protein.description <- str_match(headerlines, ">\\S\\S\\|(\\S+)\\|(\\S+.*)")[,3]
protein.name <- str_match(headerlines, ">\\S\\S\\|(\\S+)\\|(\\S+)_.*")[,3]

# Create a named vector to lookup protein descriptions by accession
getProteinDescription <- protein.description
names(getProteinDescription) <- protein.accession
getProteinName <- protein.name
names(getProteinName) <- protein.accession



# Load ProtMapper kinase substrate table
protmapper <- as.data.frame(read_excel(path = "~/OneDrive - The Mount Sinai Hospital/DB/20230322-ProtMapper-kinase-substrate-table.xlsx"))


# Construct a phacc column
protmapper$PH_ACC <- paste(protmapper$TARGET_UP_ID, "_ph", protmapper$TARGET_POS, sep = "")

# Filter for belief = 1
protmapper <- protmapper[which(protmapper$BELIEF == 1),]

# Filter for kinases
protmapper <- protmapper[which(protmapper$CTRL_IS_KINASE == TRUE),]

# Make a table of kinases and substrates
kinsub <- protmapper[, c("CTRL_ID", "PH_ACC")]

# Change column names
colnames(kinsub) <- c("KINASE", "SUB_PH_ACC")

# Create a list to store kinase-substrate annotations
kinases_gsea <- list()

# Get list of kinases
kinases <- unique(kinsub$KINASE)

# Make a lookup table from kinase name to gene name
getGeneName <- protmapper$CTRL_GENE_NAME
names(getGeneName) <- protmapper$CTRL_ID


# Loop through each kinase and load into the list
for (i in 1:length(kinases)) {
  
  # Get the substrates of the current kinase
  #substrates <- unique(kinsub_allph[which(kinsub_allph$KINASE == kinases[i]), "PhAcc"])
  
  substrates <- unique(kinsub[which(kinsub$KINASE == kinases[i]), "SUB_PH_ACC"])
  
  # Store list of substrates in the kinases_gsea list
  kinases_gsea[[i]] <- substrates
  
}

# Name each element in kinases_gsea list with KINASE name
names(kinases_gsea) <- kinases



# Read in combined results
results <- read.delim(file = "../msstats/20240130-WNV-ph-results-ann.txt", quote = "", stringsAsFactors = FALSE)

# Get a list of comparisons
comparisons <- unique(results$Label)

# Create a data frame to store GSEA results
all_ksea_results <- data.frame(
  comparison = character(),
  kinase = character(),
  ES = double(),
  pval = double(),
  nsub = double()
)

# Loop through comparisons
for (i in 1:length(comparisons)) {
  
  # Get current comparison
  current_comparison <- comparisons[i]
  
  # Print current comparison
  print(current_comparison)
  
  # Get l2fc data for current comparison
  current_l2fc <- results[which((results$Label == current_comparison) & (is.finite(results$log2FC))), "log2FC"]
  
  # Get names
  current_l2fc_names <- results[which((results$Label == current_comparison) & (is.finite(results$log2FC))), "Protein"]
  
  # Rank order l2fc
  current_l2fc_ranked <- current_l2fc[order(current_l2fc, decreasing = TRUE)]
  
  # Order corresponding protein names by same ranking
  names(current_l2fc_ranked) <- current_l2fc_names[order(current_l2fc, decreasing = TRUE)]
  
  
  #ksea_result <- ksea_batchKinases(ranking = names(current_l2fc_ranked), logvalues = current_l2fc_ranked, regulons = kinases_gsea, trial = 1000)
  
  # Loop through each kinase
  for (j in 1:length(kinases)) {
    
    # Get the current kinase
    current_kinase <- kinases[j]
    
    # Print the current kinase
    print(current_kinase)
    
    # Check that we have measurements for at least 3 substrates
    
    nsub <- length(intersect(names(current_l2fc_ranked), kinases_gsea[[current_kinase]]))
    
    # Run GSEA for this comparison if nsub > 3
    if (nsub > 3) {
      
      # Run ksea
      ksea_result <- ksea(names(current_l2fc_ranked), current_l2fc_ranked, kinases_gsea[[current_kinase]], trial=1000, significance = TRUE)
      
      # Print ksea result
      print (paste(current_comparison, j, "of", length(kinases), current_kinase, ksea_result$ES, ksea_result$p.value))
      
      # Add ksea result to table of all results
      all_ksea_results <- rbind(
        data.frame(
          comparison = current_comparison,
          kinase = current_kinase,
          ES = ksea_result$ES,
          pval = ksea_result$p.value,
          nsub = nsub
        ), all_ksea_results
      )
      
    }
  }
}

# Add gene name to ksea results
all_ksea_results$kinase_gene <- getGeneName[all_ksea_results$kinase]

# Adjust p-values
all_ksea_results$adj.pval <- p.adjust(p = all_ksea_results$pval, method = "fdr")

# Write to a file
write.table(x = all_ksea_results, file = paste(prefix, "-ksea.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)




# Make a volcano plot of ksea results

# Filter ksea results for no NAs in kinase gene column
all_ksea_results <- all_ksea_results[which(!(is.na(all_ksea_results$kinase_gene))),]


# Add LogP column
all_ksea_results$LogP <- -1 * log(all_ksea_results$pval) / log(10)

# Get hits with pval < 0.05 to label on volcano plot
sig_ksea <- all_ksea_results[which(all_ksea_results$pval < 0.05),]

# Plot
ggplot(data = all_ksea_results, aes(x = all_ksea_results$ES, y = all_ksea_results$LogP)) +
  geom_point(aes(size = all_ksea_results$nsub)) +
  theme_bw() +
  geom_point(data = sig_ksea, aes(x = sig_ksea$ES, y = sig_ksea$LogP, size = sig_ksea$nsub), color = "red") +
  geom_text(data = sig_ksea, aes(x = sig_ksea$ES, y = sig_ksea$LogP, label = sig_ksea$kinase_gene), vjust = 1.5, hjust = -0.05) +
  scale_size(name = "# substrates") +
  xlab(label = "KSEA Enrichment Score") +
  ylab(label = "-Log10(p-value)") +
  theme(aspect.ratio = 1)
ggsave(filename = paste(prefix, "-ksea-volcano.pdf", sep = ""), width = 6, height = 6)





