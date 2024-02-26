library(reshape2)
library(stringr)
library(seqinr)
library(ggseqlogo)


# Set working directory
setwd("~/OneDrive - The Mount Sinai Hospital/Data/Collaborations/RamageLab/20240130-WNV-reanalysis/20240130-WNV-ph/motif-x/")


# Set prefix for naming output files
prefix <- "202401130-WNV-ph"


# Read in results
results <- read.delim(file = "../msstats/20240130-WNV-ph-results-ann.txt", quote = "", stringsAsFactors = FALSE)





# Set db files
dbfiles <- "~/OneDrive - The Mount Sinai Hospital/DB/20191017.SwissProt.Hsapiens.WNVNY99.fasta"

numdb <- length(dbfiles)

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

# Set a flag to zero for when to read in sequences
grabsequence <- 0

# Declare lists to store sequences and accessions
sequencelist <- list()
accessionlist <- list()

# Loop through databases line by line
for (i in 1:length(dblines)) {
  # If the line contains a ">" then this is a description line
  if(grepl(pattern = ">", x = dblines[i])) {
    # If this is not the first line then append the sequence and accession to their respective lists
    if (i > 1) {
      sequencelist <- append(sequencelist, sequence)
      accessionlist <- append(accessionlist, accession)
    }
    # Extract accession based on the pattern: ">sp|P12345|ASDF_HUMAN blah blah", where P12345 is accession
    accession <- str_match(string = dblines[i], pattern = ">\\S\\S\\|([^\\|]*)\\|.*")[2]
    
    # Set grab sequence flag to 1
    grabsequence <- 1
  } else {
    if (grabsequence == 1) {
      # If this is not a header line and grabsequence flag is set then read line into sequence
      sequence <- dblines[i]
      # And set the grabsequence flag back to zero
      grabsequence <- 0
    } else {
      # Otherwise append the sequence to the existing sequqence
      sequence <- paste(sequence, dblines[i], sep = '', collapse = '')
    }
  }
}

# Grab the last ones in the file
sequencelist <- append(sequencelist, sequence)
accessionlist <- append(accessionlist, accession)

# Store sequences in a named vector for quick retrieval by accession
getSequence <- sequencelist
names(getSequence) <- accessionlist




# Get phaccs
phaccs <- as.data.frame(unique(results$Protein))
colnames(phaccs) <- "Accession"

# Remove any lines with _phNA
phaccs <- phaccs[which(!(str_detect(string = phaccs$Accession, pattern = "_phNA"))),]

phaccs <- data.frame(Accession = phaccs)

###################

# Get flanking sequences



# Loop through list of phaccs
for (i in 1:nrow(phaccs)) {
  
  # Check if this is an unambiguous ID
  if (str_detect(string = phaccs$Accession[i], ";")) {
    
    # Split by semi colons
    tmp <- str_split(string = phaccs$Accession[i], pattern = ";")[[1]][1]
    
    # Grab the first one to use to find flankseq
    phacc_split <- str_split(string = tmp, pattern = "_ph")
  } else {
    # Otherwise just split ph accession on _ph characters
    phacc_split <- str_split(string = phaccs$Accession[i], pattern = "_ph")
    
  }
  
  
  # Get protein accession
  protein.accession <- phacc_split[[1]][1]
  
  # Loop through sites
  for (j in 2:length(phacc_split[[1]])) {
    
    # Get sequence for this protein
    sequence <- getSequence[protein.accession]
    
    # Get the length of the sequence
    seqlength <- str_length(sequence)
    
    # Get phosphosite position
    position <- as.numeric(phacc_split[[1]][j])
    
    # Calculate start and end of positions in string to grab
    start <- position - 7
    end <- position + 7
    
    # Check that there are at least 7 AAs before and after position to grab
    if (start < 1) {
      start <- 1
    }
    if (end > seqlength) {
      end = seqlength
    }
    
    # First grab the sequence up to and including the phosphosite
    flankseq <- str_sub(string = sequence, start = start, end = position)
    
    # Add a * to denote the phospho position
    flankseq <- paste(flankseq, "*", sep = "", collapse = "")
    
    # Extract sequence
    flankseq <- paste(flankseq, str_sub(string = sequence, start = position + 1, end = end), sep = "", collapse = "")
    
    
    
    # If this is a multiply phosphorylated site then concatenate by | characters
    if (j > 2) {
      flankseq_out <- paste(flankseq_out, flankseq, sep = "|")
    } else {
      flankseq_out <- flankseq
    }
    
    
    
  }
  # Add flankseq to phacc dataframe
  phaccs[i, "Flankseq"] <- flankseq_out
  
  
  
  
}


# Create a named vector of flankseqs
getFlankseq <- phaccs$Flankseq
names(getFlankseq) <- phaccs$Accession



# Add Flankseq to results
results$Flankseq <- getFlankseq[results$Protein]

# Write to a file
write.table(x = results, file = paste(prefix, "-results-flankseq.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)







# Filter out flankseq with length other than 16
results <- results[which(str_length(results$Flankseq) == 16),]

# Remove *s in flankseq
results$Flankseq <- str_replace(string = results$Flankseq, pattern = "\\*", replacement = "")


# Set values for diff abundance
l2fc_min <- -1
l2fc_max <- 1
pval_max <- 0.05
adjpval_max <- 0.05



# Get a list of comparisons
comparisons <- unique(results$Label)


# Loop through comparisons and write out files for motif-x analysis
for (i in 1:length(comparisons)) {
  
  # Get the current comparison
  current_comparison <- comparisons[i]
  
  # Get list of upreg flankseq
  upreg_flankseq <- results[which(
    (results$Label == current_comparison) &
      (results$log2FC > l2fc_max) &
      (results$pvalue < pval_max) &
      (results$adj.pvalue < adjpval_max)), "Flankseq"]
  
  # Get names
  upreg_names <- results[which(
    (results$Label == current_comparison) &
      (results$log2FC > l2fc_max) &
      (results$pvalue < pval_max) &
      (results$adj.pvalue < adjpval_max)), "Protein"]
  
  # Write out as fasta file
  write.fasta(sequences = as.list(upreg_flankseq), names = upreg_names, file.out = paste(current_comparison, "-upreg.fasta", sep = ""))
  
  # Get list of downreg flankseq
  downreg_flankseq <- results[which(
    (results$Label == current_comparison) &
      (results$log2FC < l2fc_min) &
      (results$pvalue < pval_max) &
      (results$adj.pvalue < adjpval_max)), "Flankseq"]
  
  # Get names
  downreg_names <- results[which(
    (results$Label == current_comparison) &
      (results$log2FC < l2fc_min) &
      (results$pvalue < pval_max) &
      (results$adj.pvalue < adjpval_max)), "Protein"]
  
  # Write out as fasta file
  write.fasta(sequences = as.list(downreg_flankseq), names = downreg_names, file.out = paste(current_comparison, "-downreg.fasta", sep = ""))
  
  
  
  
  
  # Get list of background flankseq
  all_flankseq <- results[which(results$Label == current_comparison), "Flankseq"]
  
  # GEt all names
  all_names <- results[which(results$Label == current_comparison), "Protein"]
  
  # Write out as fasta file
  write.fasta(sequences = as.list(all_flankseq), names = all_names, file.out = paste(current_comparison, "-all.fasta", sep = ""))
  
}




