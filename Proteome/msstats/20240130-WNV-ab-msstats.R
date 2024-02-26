library(MSstats)
library(reshape2)
library(stringr)
library(ggplot2)
library(gplots)



# Set working directory 
setwd("~/OneDrive - The Mount Sinai Hospital/Data/Collaborations/RamageLab/20240130-WNV-reanalysis/20240130-WNV-ab/msstats/")


# Set prefix fo rnaming output files
prefix <- "20240130-WNV-ab"


# Read in evidence file
evidence <- read.delim(file = "101019-wnv-abun-evidence.txt", quote = "", stringsAsFactors = FALSE)


# Aggregate intensities by protein/peptide/charge/run
evidence <- aggregate(Intensity ~ Proteins + Modified.sequence + Charge + Raw.file, evidence, sum)


# Read in keys file
keys <- read.delim(file = "122516-hr-406-421-abun-keys.txt", quote = "", stringsAsFactors = FALSE)

# Replace infected with WNV and uninfected with mock
keys$Condition <- str_replace(string = keys$Condition, pattern = "Infected", replacement = "WNV")
keys$Condition <- str_replace(string = keys$Condition, pattern = "Uninfected", replacement = "Mock")



# Create lookup vector to add conditions and bioreplicate to evidence
getCondition <- keys$Condition
names(getCondition) <- keys$RawFile
getBioReplicate <- keys$BioReplicate
names(getBioReplicate) <- keys$RawFile

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


# Get unique features and proteins
proteinFeatures <- unique(quant[, c("ProteinName", "Feature", "Count")])

# Count features per protein
featuresPerProtein <- aggregate(Count ~ ProteinName, proteinFeatures, sum)

# Get proteins with <2 feature to exclue
exclude <- featuresPerProtein[which(featuresPerProtein$Count < 2), "ProteinName"]

# Filter these out 
quant <- quant[which(!(quant$ProteinName %in% exclude)),]



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
  scale_y_continuous(limits = c(0,6000), expand = c(0,0)) + 
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






# Get a list of isgs
isgs <- read.delim(file = "~/OneDrive - The Mount Sinai Hospital/DB/ISGs.txt", quote = "", stringsAsFactors = FALSE)

# Filter results for isgs
results_isgs <- results[which(results$Protein %in% isgs$Entry),]

# Add a logp column
results_isgs$logp <- -1 * log(results_isgs$pvalue) / log(10)

# Create a column to lbael upreg ISgs
results_isgs$Text <- ""
results_isgs[which(results_isgs$log2FC > 0.75), "Text"] <- results_isgs[which(results_isgs$log2FC > 0.75), "Protein.Name"]

# Copy upreg ISGs to a separate df to plot in red
results_isgs_upreg <- results_isgs[which(results_isgs$log2FC > 0.75), ]


# Plot as a volcano plot
ggplot(data = results_isgs, aes(x = results_isgs$log2FC, y = results_isgs$logp)) +
  geom_point() +
  geom_point(data = results_isgs_upreg, aes(x = results_isgs_upreg$log2FC, y = results_isgs_upreg$logp), color = "red") +
  geom_text(aes(label = results_isgs$Text), vjust = 1.5) +
  theme_bw() + 
  xlab(label = "Log2FC (WNV/Mock)") +
  ylab(label = "-Log10(p-value)")
ggsave(filename = paste(prefix, "-isg-volcano.pdf", sep = ""), width = 4, height = 4)





# Get a list of upreg ISGs to plot as faceted intensity boxplots
isgs_upreg <- results_isgs_upreg$Protein

# Melt down subquant
subquant_melted <- melt(data = subquant, id.vars = "Protein", variable.name = "Sample", value.name = "Log2Intensity")

# Parse out condition from sample name
subquant_melted$Condition <- str_match(string = subquant_melted$Sample, pattern = "(\\S+)_.*")[,2]


# Filter subquant for the list
subquant_list <- subquant_melted[which(subquant_melted$Protein %in% isgs_upreg),]

# Get protein names
subquant_list$Protein.Name <- getProteinName[as.character(subquant_list$Protein)]

# Plot
ggplot(data = subquant_list, aes(x = subquant_list$Condition, y = subquant_list$Log2Intensity)) +
  geom_boxplot() +
  facet_wrap(Protein.Name ~ .) +
  theme_bw() +
  xlab(label = "Condition") + 
  ylab(label = "Log2Intensity") 
ggsave(filename = paste(prefix, "-upregISGs-boxplots.pdf", sep = ""), width = 4, height = 4)



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
  scale_x_continuous(limits = c(-4, 4)) +
  scale_y_continuous(limits = c(0, 10), expand = c(0,0)) +
  xlab(label = "Log2FC (WNV/Mock)") +
  ylab(label = "-Log10(p-value)") 
ggsave(filename = paste(prefix, "-volcano.pdf", sep = ""), width = 5, height = 5)
