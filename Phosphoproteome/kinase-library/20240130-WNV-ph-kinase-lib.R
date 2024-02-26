library(ggplot2)
library(reshape2)




# Set working directory
setwd("~/OneDrive - The Mount Sinai Hospital/Data/Collaborations/RamageLab/20240130-WNV-reanalysis/20240130-WNV-ph/kinase-library//")


# Set prefix for naming output files
prefix <- "20240130-WNV-ph"


# Read in kinase library reuslts
results <- read.delim(file = "enrichment-analysis-result-table.txt", quote = "", stringsAsFactors = FALSE)

# Get sig results with pval < 0.001
sig_results <- results[which(results$dominant_adjusted_p_value_log10_abs > 5),]

# GEt results for PAK kinase
PAKs <- results[which(results$kinase %in% c("PAK2", "PAK3")),]


# Plot as a volcano
ggplot(data = results, aes(x = results$dominant_enrichment_value_log2, y = results$dominant_p_value_log10_abs)) +
  geom_point() +
  geom_point(data = sig_results, aes(x = sig_results$dominant_enrichment_value_log2, y = sig_results$dominant_p_value_log10_abs), color = "red") +
  geom_text(data = sig_results, aes(x = sig_results$dominant_enrichment_value_log2, y = sig_results$dominant_p_value_log10_abs, label = sig_results$kinase), color = "red") +
  theme_bw() +
  xlab(label = "Log2(Enrichment Value)") +
  ylab(label = "-Log10(p-value)") +
  scale_x_continuous(limits = c(-2.5,2.5)) +
  theme(aspect.ratio = 1) +
  geom_hline(yintercept = 1.31, linetype = "dashed")
ggsave(filename = paste(prefix, "-kinase-library-volcano.pdf", sep = ""), width = 6, height = 6)
