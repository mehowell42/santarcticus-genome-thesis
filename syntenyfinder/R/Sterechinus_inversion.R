# R script for metazoan synteny analysis -- HAPLOTYPE (Sa1 vs Sa2) VERSION,
# ribbons and dot plots coloured by inferred inversion
# Adapted from Isabel Jiah-Yih Liao, Tom Lewin, Yi-Jyun Luo's synteny_v5.5.R
# Symbiosis Genomics & Evolution Lab 2024
#
# CHANGES FROM YOUR PASTED VERSION (colouring only):
#   1. Removed the chromosome-letter palette and the "untraceable" filter.
#      Colour now comes from the inversion classification, so those rows no
#      longer need to be dropped and every single-copy gene is kept.
#   2. Added one console line per pair reporting how many genes were called
#      colinear / inverted / unassessed, so you can sanity-check the colouring.
# Everything else is unchanged.
#
# HOW THE COLOURING WORKS
#   Gene strand is not kept in SyntenyFinder's output, so inversions are
#   inferred from gene order, as in an Oxford dot plot. For each chromosome
#   pair the dominant orientation is taken from the overall Spearman rho.
#   A sliding window (+/- 5 genes) then asks whether local gene order runs
#   the same way or the opposite way.
#     grey  (BBBBBB) = colinear with the dominant orientation
#     purple   (B19CD7) = local order runs opposite (candidate inversion)
#     pale  (EEEEEE) = too few genes in that chromosome pair to assess
#   This is a heuristic, not formal inversion calling. Inversions shorter than
#   about 6 genes will not show, and translocations can also appear red.

# Setting the working directory
setwd("C:/Users/maddi/OneDrive - James Cook University/Thesis/R/syntenyfinder/sterechinus_inversion/")

# Load required packages
library(RIdeogram)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(ape)

# List of species
species_list <- c("Sne", "San")

# Assign source and target species (overwritten by the loop below; only
# matters if you run this manually outside the loop, same as the original)
source <- "Sne"
target <- "San"

# Initialize a matrix to store mixing rate results
mixing_rate_matrix <- matrix(NA, nrow = length(species_list), ncol = length(species_list),
                             dimnames = list(species_list, species_list))

# Initialize a dataframe to store mixing rate results
chr_mixing_rate_data <- data.frame(Species1 = character(), Species2 = character(), OrthoChromosome = character(),
                                   MixingRate = character(), stringsAsFactors = FALSE)

classify_orientation <- function(df, window = 5, min_genes = 8) {
  df$orient <- "unassessed"
  pairs <- unique(df[, c("Species_1", "Species_2")])
  for (p in seq_len(nrow(pairs))) {
    idx <- which(df$Species_1 == pairs$Species_1[p] & df$Species_2 == pairs$Species_2[p])
    if (length(idx) < min_genes) next
    idx <- idx[order(df$Start_1[idx])]
    x <- df$Start_1[idx]; y <- df$Start_2[idx]
    rho <- suppressWarnings(cor(x, y, method = "spearman"))
    if (is.na(rho)) next
    dominant <- ifelse(rho >= 0, 1, -1)
    n <- length(idx)
    for (i in seq_len(n)) {
      lo <- max(1, i - window); hi <- min(n, i + window)
      if (hi - lo + 1 < 4) next
      r <- suppressWarnings(cor(x[lo:hi], y[lo:hi], method = "spearman"))
      if (is.na(r)) next
      df$orient[idx[i]] <- ifelse(ifelse(r >= 0, 1, -1) == dominant, "colinear", "inverted")
    }
  }
  df
}

# Loop over each pair of species
for (source in species_list) {
  for (target in species_list) {
    if (source != target) {  # Skip if source and target are the same
      
      # Setup filter for jumping genes (default = 5)
      jumping_genes = 5
      
      # If you would like to reorder the chromosomes, put the chromosome order here.
      # Otherwise, set reorder to FALSE
      # Each vector must be a permutation of 1:N (check with
      # sort(X_order) == 1:length(X_order)).
      reorder <- TRUE
      Sne_order <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
      San_order <- c(1,2,3,10,11,12,13,7,6,9,16,14,5,18,20,15,4,17,19,8)
      
      source_order <- get(paste0(source, "_order"))
      target_order <- get(paste0(target, "_order"))
      
      # List any chromosomes in either the source or the target genome which should be flipped
      Sne_flip <- c()
      San_flip <- c()
      source_flip <- get(paste0(source, "_flip"))
      target_flip <- get(paste0(target, "_flip"))
      
      ########################################
      ### 1. Set karyotypes
      ########################################
      
      source_karyotype <- read.table(paste0("input/", source, "_karyotype.txt"),
                                     sep = "\t", header = TRUE, stringsAsFactors = F)
      target_karyotype <- read.table(paste0("input/", target, "_karyotype.txt"),
                                     sep = "\t", header = TRUE, stringsAsFactors = F)
      
      source_karyotype <- source_karyotype[order(source_karyotype$End,
                                                 decreasing = TRUE), ]
      target_karyotype <- target_karyotype[order(target_karyotype$End,
                                                 decreasing = TRUE), ]
      
      source_mapping <- setNames(1:nrow(source_karyotype), source_karyotype$Chr)
      target_mapping <- setNames(1:nrow(target_karyotype), target_karyotype$Chr)
      
      source_karyotype$orig_ids <- source_karyotype$Chr
      target_karyotype$orig_ids <- target_karyotype$Chr
      
      source_karyotype$Chr <- source_mapping[source_karyotype$Chr]
      target_karyotype$Chr <- target_mapping[target_karyotype$Chr]
      
      if (reorder == TRUE) {
        source_karyotype <- source_karyotype[order(match(source_karyotype$Chr, source_order)),]
        target_karyotype <- target_karyotype[order(match(target_karyotype$Chr, target_order)),]
        source_remap <- setNames(1:nrow(source_karyotype), source_order)
        target_remap <- setNames(1:nrow(target_karyotype), target_order)
      }
      
      karyotype_all <- rbind(source_karyotype, target_karyotype,
                             make.row.names=FALSE)
      karyotype_all$orig_ids <- NULL
      karyotype_all$fill <- "cccccc"
      
      karyotype_all <- karyotype_all[, c(1,2,3,7,4,5,6)]
      
      karyotype_all <- as.data.frame(karyotype_all, stringsAsFactors = F )
      head(karyotype_all)
      
      ########################################
      ### 2. Set orthologous gene coordinates
      ########################################
      
      source_coords <- read.table(paste0("input/", source, "_coordinates.tsv"),
                                  sep = "\t", header = FALSE, stringsAsFactors = F)
      target_coords <- read.table(paste0("input/", target, "_coordinates.tsv"),
                                  sep = "\t", header = FALSE, stringsAsFactors = F)
      
      nrow(source_coords)
      nrow(target_coords)
      
      source_coords <- source_coords[source_coords$V3 %in% source_karyotype$orig_ids,]
      target_coords <- target_coords[target_coords$V3 %in% target_karyotype$orig_ids,]
      
      nrow(source_coords)
      nrow(target_coords)
      
      head(source_coords)
      head(target_coords)
      
      source_coords$V3 <- source_mapping[source_coords$V3]
      target_coords$V3 <- target_mapping[target_coords$V3]
      
      for (chromosome in source_flip){
        chromosome_length = source_karyotype[source_karyotype$Chr == chromosome, 'End']
        new_coords <- (chromosome_length - source_coords[source_coords$V3 == chromosome, c("V4", "V5")])
        source_coords[source_coords$V3 == chromosome, c("V4", "V5")] <- new_coords
      }
      for (chromosome in target_flip){
        chromosome_length = target_karyotype[target_karyotype$Chr == chromosome, 'End']
        new_coords <- (chromosome_length - target_coords[target_coords$V3 == chromosome, c("V4", "V5")])
        target_coords[target_coords$V3 == chromosome, c("V4", "V5")] <- new_coords
      }
      
      if (reorder == TRUE) {
        source_coords$V3 <- source_remap[as.character(source_coords$V3)]
        target_coords$V3 <- target_remap[as.character(target_coords$V3)]
      }
      
      # Create connections
      coords_all <- (merge(source_coords, target_coords, by = "V1")
                     [, c(-1, -2, -6, -7)])
      head(coords_all)
      
      # Relabel columns (the last column holds the chromosome-trace label; it
      # is not used for colouring in this version)
      colnames(coords_all) <- c("Species_1", "Start_1", "End_1", "Species_2", "Start_2", "End_2", "ALG")
      
      # Filter out the jumping genes first, so sparse chromosome pairs don't get classified
      coords_all <- coords_all %>%
        group_by(Species_1, Species_2) %>%
        filter(n() > jumping_genes) %>%
        ungroup()
      coords_all <- as.data.frame(coords_all)
      
      # Classify each gene as colinear or inverted relative to its chromosome pair's dominant orientation
      coords_all <- classify_orientation(coords_all)
      
      # Report the classification so the colouring can be sanity-checked
      orient_counts <- table(coords_all$orient)
      cat(source, "vs", target, "--",
          paste(names(orient_counts), orient_counts, sep = ": ", collapse = ", "), "\n")
      
      orient_hexmap <- c("unassessed" = "EEEEEE", "colinear" = "BBBBBB", "inverted" = "B19CD7")
      coords_all$fill <- orient_hexmap[coords_all$orient]
      
      # Draw order: unassessed (bottom), colinear, then inverted on top
      sort_key <- match(coords_all$orient, c("unassessed", "colinear", "inverted"))
      coords_all <- coords_all[order(sort_key), ]
      
      # Select the 7 columns by NAME. RIdeogram reads them by position, and we now
      # have extra columns (ALG, orient) that must not be passed through.
      synteny_coord <- as.data.frame(coords_all[, c("Species_1","Start_1","End_1",
                                                    "Species_2","Start_2","End_2","fill")])
      head(synteny_coord)
      head(karyotype_all)
      
      ########################################
      ### 3. Create synteny plot
      ########################################
      
      ideogram(karyotype = karyotype_all, synteny = synteny_coord,
               output = paste0(source, "_", target, "_synteny.svg"))
      
      convertSVG(paste0(source, "_", target, "_synteny.svg"),
                 file = paste0(source, "_", target, "_synteny.pdf"),
                 device = "pdf")
      
      ########################################
      ### 4. Generate Oxford dot plot
      ########################################
      
      source_karyotype$Cumulative_End <- 0
      target_karyotype$Cumulative_End <- 0
      
      for (i in 2:nrow(source_karyotype)) {
        source_karyotype$Cumulative_End[i] <- sum(source_karyotype$End[1:(i-1)])
      }
      for (i in 2:nrow(target_karyotype)) {
        target_karyotype$Cumulative_End[i] <- sum(target_karyotype$End[1:(i-1)])
      }
      
      source_cumsum <- source_karyotype[source_karyotype$species == source, "Cumulative_End"]/1e6
      target_cumsum <- target_karyotype[target_karyotype$species == target, "Cumulative_End"]/1e6
      
      synteny_coord$adjusted_start_1 <- synteny_coord$Start_1 + source_karyotype$Cumulative_End[synteny_coord$Species_1]
      synteny_coord$adjusted_start_2 <- synteny_coord$Start_2 + target_karyotype$Cumulative_End[synteny_coord$Species_2]
      
      # Midpoint of each chromosome (in Mb), in the displayed (reordered) order
      source_mid <- (source_karyotype$Cumulative_End + source_karyotype$End / 2) / 1e6
      target_mid <- (target_karyotype$Cumulative_End + target_karyotype$End / 2) / 1e6
      
      synteny_coord$fill <- paste0("#", synteny_coord$fill)
      dp <- ggplot(synteny_coord, aes(x=adjusted_start_1/1e6, y=adjusted_start_2/1e6, color=fill)) +
        scale_color_identity() +
        geom_point(size=0.5) +
        scale_x_continuous(breaks = source_mid, labels = source_karyotype$Chr, expand = c(0, 0)) +
        scale_y_continuous(breaks = target_mid, labels = target_karyotype$Chr, expand = c(0, 0)) +
        theme_classic() +
        theme(aspect.ratio = 1,
              axis.text = element_text(size = 7))
      
      for (x in source_cumsum) {
        dp <- dp + geom_vline(xintercept=x, color="gray", linewidth=0.3)
      }
      for (y in target_cumsum) {
        dp <- dp + geom_hline(yintercept=y, color="gray", linewidth=0.3)
      }
      dp <- dp + labs(x = paste(source, "chromosome"), y = paste(target, "chromosome")) +
        theme(legend.position = "none")
      
      synteny_coord$fill <- gsub("#", "", synteny_coord$fill)
      
      ggsave(paste0(source, "_", target, "_oxford_plot.pdf"), dp,
             device = "pdf", width = 7, height = 7)
      
      ########################################
      ### 5. Calculate genome-wide synteny
      ###    mixing rate (MACROSYNTENY)
      ########################################
      
      synteny_coord$rank_1 <- rank(synteny_coord$adjusted_start_1)
      synteny_coord$rank_2 <- rank(synteny_coord$adjusted_start_2)
      
      spearman_coefficient <- cor(synteny_coord$adjusted_start_1, synteny_coord$adjusted_start_2, method = "spearman")
      
      mixing_rate <- 1 - abs(spearman_coefficient)
      
      mixing_rate_matrix[source, target] <- mixing_rate
      
      ########################################
      ### 6. Calculate synteny mixing rate for
      ###    orthologous chromosome pairs (MICROSYNTENY)
      ########################################
      
      chromosomes <- sort(unique(synteny_coord$Species_1))
      
      for (chr in chromosomes) {
        chr_data <- subset(synteny_coord, Species_1 == chr)
        
        chr_data$rank_1 <- rank(chr_data$adjusted_start_1)
        chr_data$rank_2 <- rank(chr_data$adjusted_start_2)
        
        spearman_coefficient <- cor(chr_data$adjusted_start_1, chr_data$adjusted_start_2, method = "spearman")
        
        chr_mixing_rate <- 1 - abs(spearman_coefficient)
        
        new_row <- data.frame(
          Species1 = source,
          Species2 = target,
          OrthoChromosome = chr,
          MixingRate = chr_mixing_rate,
          stringsAsFactors = FALSE
        )
        chr_mixing_rate_data <- rbind(chr_mixing_rate_data, new_row)
      }
    }
  }
}

# End of loop.

diag(mixing_rate_matrix) <- 0

########################################################################
### Export the mixing rate results
########################################################################

write.csv(mixing_rate_matrix, "mixing_rate_genome_wide.csv")
write.csv(chr_mixing_rate_data, "mixing_rate_per_chromosome.csv", row.names = FALSE)

########################################################################
### (a) Create a heatmap visualizing the genome-wide synteny mixing rate
########################################################################

head(mixing_rate_matrix)

mixing_rate_data <- melt(mixing_rate_matrix, varnames = c("Source", "Target"), value.name = "Mixing_rate")

heatmap1 <- ggplot(mixing_rate_data, aes(x = Source, y = Target, fill = Mixing_rate)) +
  geom_tile() +
  scale_fill_gradient2(low = "royalblue", mid = "white",high = "darkorange",
                       midpoint = 0.3, limit = c(0, 0.6)) +
  theme_minimal() +
  coord_fixed(ratio = 1) +
  xlab("Source Species") +
  ylab("Target Species") +
  ggtitle("Heatmap of syteny mixing rate")

print(heatmap1)

my_palette <- colorRampPalette(c("royalblue", "white", "darkorange"))(n = 64)
heatmap2 <- pheatmap(mixing_rate_matrix,
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",
                     cluster_rows=TRUE, col=my_palette, scale="none",
                     key=TRUE, symkey=FALSE, density.info="none", trace="none",
                     cexCol=1, cexRow=1, cellwidth=30, cellheight=30)

print(heatmap2)

ggsave("mixing_rate_heatmap_1.pdf", heatmap1, device = "pdf")
ggsave("mixing_rate_heatmap_2.pdf", heatmap2, device = "pdf")
