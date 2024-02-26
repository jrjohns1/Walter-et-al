# Load required libraries
library(MSstats)
library(ggplot2)
library(reshape2)
library(stringr)
library(gplots)
library(progress)


# Set working directory
setwd("~/OneDrive - The Mount Sinai Hospital/Data/Collaborations/RamageLab/20240130-WNV-reanalysis/20240130-WNV-ph/msstats/")


# Set prefix for naming output files
prefix <- "20240130-WNV-ph"


# Set database files
dbfiles <- c("~/OneDrive - The Mount Sinai Hospital/DB/20191010.SwissProt.Hsapiens.fasta", "~/OneDrive - The Mount Sinai Hospital/DB/20191010.SwissProt.WNVNY99.fasta")




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


# Read in evidence file
raw <- read.delim(file = "101719-wnv-ph-evidence.txt", quote = "", stringsAsFactors = FALSE)






# Filter for only peptides containing phospho mods
raw_ph <- raw[grep(pattern = "\\(Phospho \\(STY\\)\\)", x = raw$Modified.sequence),]

# Collapse to unique accessions
raw_ph <- unique(raw_ph[, c("Proteins", "Modified.sequence", "Sequence")])

# Extract the text within () brackets in each modified peptide sequence (these are mods)
raw_ph$mods <- str_match_all(string = raw_ph$Modified.sequence, pattern = "\\(([^(]+\\([^)]+\\))\\)")

# Extract the positions of [] brackets in each modified peptide sequence
raw_ph$mod_positions <- str_locate_all(string = raw_ph$Modified.sequence, pattern = "\\(([^(]+\\([^)]+\\))\\)")

# Count the number of [] brackets in each modified peptide sequence
raw_ph$mod_count <- str_count(string = raw_ph$Modified.sequence, pattern = "\\(([^(]+\\([^)]+\\))\\)")

# Split out accessions containing semicolon characters
splitaccessions <- str_split(string = raw_ph$Proteins, pattern = ";")

# Create a list to store phosphosite-specific accessions by row
ph_accessions <- list()

# Calculate number of rows in data
numrows <- nrow(raw_ph)



pb <- progress_bar$new(total = numrows)


# Now loop through each row pf data one by one
for (i in 1:numrows) {
  
  pb$tick()
  
  # Offset is used to determine the position in the unmodified peptide sequence that is modified
  offset <- 0
  
  # Set a counter for phosphosites within each modified peptide sequence
  phcount <- 0
  
  # Create an empty list to store the positions of phopsho modifications
  ph_positions <- list()
  
  # Loop through each modification in each modified peptide sequence
  for (j in 1:raw_ph[i,"mod_count"]) {
    if (raw_ph$mods[[i]][j,2] == "Phospho (STY)") {
      # If it's a phospho modification then calculate the position in the peptide and append to the list of phosphorylation positions
      ph_positions <- append(ph_positions, unname(raw_ph$mod_positions[[i]][j,1]) - 1 - offset)
    }
    
    # Increase the offset by the length of the modification and then go on to the next one
    offset <- offset + str_length(raw_ph$mods[[i]][j,1])
  }
  
  # Create an empty list to store the peptide positions within protein sequence
  proteinpositions <- list()
  
  # Loop through each accession split by semicolons
  for (j in 1:length(splitaccessions[[i]])) {
    # Find peptide sequence position within protein sequence
    position <- str_locate(pattern = as.character(raw_ph[i, "Sequence"]), string = getSequence[splitaccessions[[i]][j]])[1]
    
    # Appen protein position to the end of the list of positions
    proteinpositions <- append(x = proteinpositions, values = position)
    
    if (j == 1) {
      # If this is the first accession then store it as the output accession
      outaccession <- splitaccessions[[i]][j]
    } else {
      # Otherwise, append to other output accession separated by semicolon
      outaccession <- paste(outaccession, ";", splitaccessions[[i]][j], sep = '', collapse = '')
    }
    
    # Loop through all the phosphosites on this peptide
    for (k in 1:length(ph_positions)) {
      # Append to uniprot accession as _ph23, where 23 is the position of the phosphosite within the protein sequence
      outaccession <- paste(outaccession, "_ph", ph_positions[[k]]+proteinpositions[[j]]-2, collapse = '', sep = '')
    }
    
    # Finally, store the assembled accession in the list of accessions
    ph_accessions[[i]] <- outaccession
  }
  
  # Clear the ph_positions list for the next row
  rm(ph_positions)
  
}



# Reassign protein accession to phosphospecific accessions in Spectronaut report
raw_ph$Proteins <- unlist(ph_accessions)

# Make a lookup vector to get phacc from modseq
getPhAccession <- raw_ph$Proteins
names(getPhAccession) <- raw_ph$Modified.sequence

# Replace ProteinAccessions and ProteinGroups in raw with ph accessions
raw$Proteins <- getPhAccession[raw$Modified.sequence]


# Remove rows with NA for protein accessions (not phosphorylated)
raw <- raw[which(!(is.na(raw$Proteins))),]

# Write to a file (without extract 3 columns added during processing)
write.table(x = raw, file = paste(prefix, "-report-modacc.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)




# Aggregate intensities by protein/peptide/charge/run
evidence <- aggregate(Intensity ~ Proteins + Modified.sequence + Charge + Raw.file, raw, sum)




# Read in keys file
keys <- read.delim(file = "20231215-WNV-ph-keys.txt", quote = "", stringsAsFactors = FALSE)



# Create lookup vector to add conditions and bioreplicate to evidence
getCondition <- keys$Condition
names(getCondition) <- keys$FileName
getBioReplicate <- keys$BioReplicate
names(getBioReplicate) <- keys$FileName

# Add condition and bioreplicate to evidence
evidence$Condition <- getCondition[evidence$Raw.file]
evidence$BioReplicate <- getBioReplicate[evidence$Raw.file]




# Rename columns to match MSstats quant format
colnames(evidence) <- c("ProteinName", "PeptideSequence", "PrecursorCharge", "Run", "Intensity", "Condition", "BioReplicate")

# Add FragmentIon and ProductCharge columns and fill with NAs
evidence$FragmentIon <- NA
evidence$ProductCharge <- NA

# Add IsotopeLabelType and fill with Ls
evidence$IsotopeLabelType <- "L"

# Rename to quant
quant <- evidence


# Remove REV__ and CON__ entries
quant <- quant[which(!(str_detect(string = quant$ProteinName, pattern = "CON__"))),]
quant <- quant[which(!(str_detect(string = quant$ProteinName, pattern = "REV__"))),]



# Add a count column
quant$Count <- 1

# Add a feature column
quant$Feature <- paste(quant$PeptideSequence, quant$PrecursorCharge, sep = "_")

# Count runs per features
runsPerFeature <- aggregate(Count ~ Feature, quant, sum)

# Get features in <3 runs to exclude from data
exclude <- runsPerFeature[which(runsPerFeature$Count < 3), "Feature"]

# Filter out low count features
quant <- quant[which(!(quant$Feature %in% exclude)),]


# Run data process function on quant
processed <- dataProcess(raw = quant, normalization = "equalizeMedians", summaryMethod = "TMP", cutoffCensored = "minFeature", censoredInt = "NA", MBimpute = TRUE, maxQuantileforCensored = 0.999, featureSubset = "topN", n_top_feature = 25)



# Get subquant data
subquant <- quantification(processed)

# WRite subquant to a file
write.table(x = subquant, file = paste(prefix, "-subquant.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)




# Extract numerical data from subquant matrix for sample correlation heatmap
subquantmatrix <- as.matrix(subquant[,2:ncol(subquant)])

# Calculate pairwise correlations for all samples
cormat <- round(cor(subquantmatrix, use="complete.obs"),2)

# Row names are the same as column names
row.names(cormat) <- colnames(cormat)

# Set color scale from white to blue from 0 to 1 in 0.1 increments
breaks <- seq(0, 1, 0.1)
colorpalette <- colorRampPalette(c("white", "blue"))

# Plot heatmap of correlations and save to pdf
pdf(paste(prefix, "-sample-correlation-heatmap.pdf", collapse = '', sep = ''), width = 8, height = 8)
heatmap.2(x=cormat, col = colorpalette(100), trace="none", density.info = "none", margins = c(10,10))
dev.off()




# Convert subquant to T/F based on NA values
subquant_tf <- is.na(subquant)[,2:ncol(subquant)]

# Sum each column
nacounts <- colSums(subquant_tf)

# Subtract na count from total count to get protien counts
proteincounts <- nrow(subquant_tf) - nacounts

# Convert to a data frame
proteincounts <- as.data.frame(proteincounts)

# Rename column of data frame
colnames(proteincounts) <- "Count"


# Add column of data frame to label Sample names
proteincounts$Sample <- row.names(proteincounts)

# Plot protein counts by samples
ggplot(data = proteincounts, aes(x = proteincounts$Sample, y = proteincounts$Count)) + 
  geom_bar(stat = "identity", fill = "gray75", color = "black") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_y_continuous(limits = c(0,10000), expand = c(0,0)) + 
  xlab(label = "Sample") + 
  ylab(label = "Num. Proteins Quantified")
ggsave(filename = paste(prefix, "-proteincounts-bysample.pdf", sep = ""), width = 6, height = 4)




#Melt subquant data to plot intensity box plots by sample
subquant_melted <- melt(data = subquant, id.vars = "Protein", measure.vars = colnames(subquant[,2:ncol(subquant)]), variable.name = "Sample", value.name = "Log2Intensity")

# Calculate medians by sample
subquantmedians <- aggregate(Log2Intensity ~ Sample, subquant_melted, median)

# Rename column header
colnames(subquantmedians) <- c("Sample", "Log2Intensity")


# Plot

ggplot(data = subquant_melted, aes(x = subquant_melted$Sample, y = subquant_melted$Log2Intensity)) + geom_boxplot() + theme_bw() + theme(axis.text.x = element_text(angle = 90))  + xlab(label = NULL) + ylab(label = "Log2Intensity") + geom_text(data = subquantmedians, aes(x = subquantmedians$Sample, y = subquantmedians$Log2Intensity, label = round(subquantmedians$Log2Intensity, 2)), vjust = -0.5, size = 2)
ggsave(filename = paste(prefix, "-sample-intensities.pdf", sep = ""), width = 6, height = 4)





# Extract numerical data from subquant matrix for sample correlation heatmap
subquantmatrix <- as.matrix(subquant[,2:ncol(subquant)])

# Calculate principle components based on subquant matrix
pca <- prcomp(na.omit(subquantmatrix), center = TRUE, scale = TRUE)

# Get the summary
pca_summary <- summary(pca)$importance

# Get the proportion of variance of PCs
pc1_propvar <- pca_summary[2,1]
pc2_propvar <- pca_summary[2,2]

# Extract PC1 and PC2 for plotting
pcaplot <- as.data.frame(cbind(pca$rotation[,1], pca$rotation[,2]))


# Parse drug from column names
condition <- str_match(string = colnames(subquantmatrix), pattern = "(\\S+)_.*")[,2]


# Plot pca and color by cell type
ggplot(data=pcaplot, aes(x=pcaplot[,1], y=pcaplot[,2])) + geom_point(aes(color = condition), size = 3) + xlab(label = paste("Principle Component 1 (", pc1_propvar*100, "%)", sep = "")) + ylab(label = paste("Principle Component 2 (", pc2_propvar*100, "%)", sep = "")) + theme_bw() + theme(aspect.ratio = 1) + scale_color_discrete(name = "Condition") + theme(aspect.ratio = 1)
ggsave(filename = paste(prefix, "-pca.pdf", sep = ""), width = 4, height = 4)




# Read in comparisons file
comparisons <- read.delim(file = "20231215-WNV-ab-comparisons.txt", quote = "")


# Grab condition columns from comparisons
comparisonmatrix <- comparisons[,2:ncol(comparisons)]

# Reorder conditions in alphabetical order (MSstast will fail if you don't do this)
comparisonmatrix <- as.matrix(comparisonmatrix[,order(names(comparisonmatrix))])

# Add names of comparisons as row names
row.names(comparisonmatrix) <- comparisons[,1]

# Run comparison function in MSstats
compared <- groupComparison(contrast.matrix = comparisonmatrix, data = processed)

# Store MSstats results (log2fold-changes and adjusted p-values) in results data frame
results <- compared$ComparisonResult



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

# Split accession column by semi-colons
split_acc <- str_split(string = results$Protein, pattern = ";")

# Go through row by row and extract protein accessions from site accessions
# Look up protein names and descriptions by lookup vectors
# Add Protein.Accession, Protein.Name, and Protein.Description columns to results table
# This section is a little slow but I can't figure out a better way to do this
for (i in 1:length(split_acc)) {
  protein_accession <- str_match(split_acc[[i]], "([A-Z0-9]+).*")[,2]
  proteinAccessions <- paste(protein_accession, collapse = ";")
  results[i,"Protein.Accession"] <- proteinAccessions
  lookup_function <- function(x) unname(getProteinName[x])
  proteinNames <- paste(lapply(protein_accession, lookup_function), collapse = ";")
  results[i,"Protein.Name"] = proteinNames
  lookup_function <- function(x) unname(getProteinDescription[x])
  proteinDescriptions <- paste(lapply(protein_accession, lookup_function), collapse = ";")
  results[i,"Protein.Description"] <- proteinDescriptions
}

# Write results file to text
write.table(x = results, file = paste(prefix, "-results-ann.txt", sep = '',  collapse = ''), sep = "\t", quote = FALSE, row.names = FALSE)

# Melt and cast results to wide format and write to file
results_melted <- melt(data = results, id.vars = c("Protein", "Protein.Description", "Protein.Name", "Label"), measure.vars = c("log2FC", "pvalue", "adj.pvalue"))
results_wide <- dcast(data = results_melted, formula = Protein + Protein.Description + Protein.Name ~ Label + variable, value.var = "value")
write.table(x = results_wide, file = paste(prefix, "-results-wide.txt", collapse = '', sep = ''), quote = FALSE, sep = "\t", row.names = FALSE)




# Impute infinity values 

# Melt down subquant
subquant_melted <- melt(data = subquant, id.vars = "Protein", value.name = "Log2Intensity", variable.name = "Sample")

# Parse out condition from Sample name
subquant_melted$Condition <- str_match(string = subquant_melted$Sample, pattern = "([^_]+)_.*")[,2]

# Add a count column and fill with zeros
subquant_melted$Count <- 0

# Change count to 1 if not an NA value
subquant_melted[which(!(is.na(subquant_melted$Log2Intensity))), "Count"] <- 1

# Aggregate count by protein and condition
protein_nsample <- aggregate(Count ~ Protein + Condition, subquant_melted, sum)

# Cast to wide format on count
protein_nsample_wide <- dcast(data = protein_nsample, formula = Protein ~ Condition, value.var = "Count")

# Get proteins in +Inf set that are detected in all WNV and no mock samples
pos_inf_proteins <- protein_nsample_wide[which(
  (protein_nsample_wide$Mock == 0) & (protein_nsample_wide$WNV == 4)), "Protein"]

# Get proteins in -Inf set that are detected in all mock and no WNV samples
neg_inf_proteins <- protein_nsample_wide[which(
  (protein_nsample_wide$Mock == 4) & (protein_nsample_wide$WNV == 0)), "Protein"]


# Replaces -Inf and +Inf values with NAs in results
results[which(results$log2FC == -Inf), "log2FC"] <- NA
results[which(results$log2FC == Inf), "log2FC"] <- NA


# Replace all/nothing +/-Inf proteins with log2FC of 10 or -10
results[which(results$Protein %in% pos_inf_proteins), "log2FC"] <- Inf
results[which(results$Protein %in% neg_inf_proteins), "log2FC"] <- -Inf

# Write to a file
write.table(x = results, file = paste(prefix, "-results-infImpute.txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")






# Plot volcano plots

# Add logp column to results
results$LogP <- -1 * log(results$pvalue) / log(10)

# Get sig up and down
results_up <- results[which(
  (results$log2FC > 1) & (results$adj.pvalue < 0.05)),]
results_down <- results[which(
  (results$log2FC < -1) & (results$adj.pvalue < 0.05)),]

# Count # up/down phosphsites
n_up <- nrow(results_up)
n_down <- nrow(results_down)

# Plot
ggplot(data = results, aes(x = results$log2FC, y = results$LogP)) +
  geom_point() +
  geom_point(data = results_up, aes(x = results_up$log2FC, y = results_up$LogP), color = "red") +
  geom_point(data = results_down, aes(x = results_down$log2FC, y = results_down$LogP), color = "blue") +
  theme_bw() +
  scale_x_continuous(limits = c(-5.5,5.5)) +
  scale_y_continuous(limits = c(0, 11), expand = c(0,0)) +
  xlab(label = "Log2FC (WNV/Mock)") +
  ylab(label = "-Log10(p-value)") 
ggsave(filename = paste(prefix, "-volcano.pdf", sep = ""), width = 5, height = 5)



# PLot ab vs ph fold changes

# Read in ph and ab results
results_ab <- read.delim(file = "../../20240130-WNV-ab/msstats/20240130-WNV-ab-results-ann.txt", quote = "", stringsAsFactors = FALSE)



# Make a lookup vector for ab log2fc
getAbL2FC <- results_ab$log2FC
names(getAbL2FC) <- results_ab$Protein


# Add ab l2fc to ph results
results$log2FC_ab <- getAbL2FC[results$Protein.Accession]



# Plot ph vs ab l2fc
ggplot(data = results, aes(x = results$log2FC_ab, y = results$log2FC)) + 
  geom_point(size = 0.5) +
  scale_x_continuous(limits = c(-2, 4)) +
  theme_bw() +
  xlab(label = "Log2FC (WNV/Mock) Protein") +
  ylab(label = "Log2FC (WNV/Mock) Phosphosite")
ggsave(filename = paste(prefix, "-vs-ab-scatter.pdf", sep = ""), width = 5, height = 5)
