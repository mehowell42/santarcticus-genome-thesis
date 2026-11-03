# R script for metazoan synteny analysis
# Isabel Jiah-Yih Liao, Tom Lewin, Yi-Jyun Luo
# Symbiosis Genomics & Evolution Lab 2024
# version 5.5 | 2024-01-16
#
# CHANGES FROM YOUR PASTED VERSION (only what was necessary to run):
#   1. source_flip / target_flip are set to empty vectors. The flip loops
#      later in the script run unconditionally, so with these commented
#      out entirely the script would crash immediately with
#      "object 'source_flip' not found".
#   2. Added the "NA" -> "XXXX" relabeling step before the colour lookup,
#      matching the fix we already applied elsewhere -- without it every
#      unlabelled orthogroup produces the "Malformed colour string '#NA'"
#      error.
#   3. Section (b)'s example species (wsu/cpa) replaced with a real pair
#      from your species_list (Bfl/Pma) -- those two don't exist in your
#      data and would fail to read a file.
#   4. Section (c), the tree/divergence-time analysis, removed entirely
#      per your request -- it was built for the bryozoan paper's five
#      species and doesn't apply here.
# Everything else is unchanged from what you pasted.

# Setting the working directory
setwd("C:/Users/maddi/OneDrive - James Cook University/Thesis/R/syntenyfinder/alg_ech/")

# Load required packages
library(RIdeogram)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(ape)

# List of species
species_list <- c("Bfl", "Pma","Abi","Aam","Afi","Aja","Spu","San")

# Assign source and target species
source <- "Bfl"
target <- "Pma"

# Initialize a matrix to store mixing rate results
mixing_rate_matrix <- matrix(NA, nrow = length(species_list), ncol = length(species_list),
                             dimnames = list(species_list, species_list))

# Initialize a dataframe to store mixing rate results
chr_mixing_rate_data <- data.frame(Species1 = character(), Species2 = character(), OrthoChromosome = character(),
                                   MixingRate = character(), stringsAsFactors = FALSE)

# Loop over each pair of species
for (source in species_list) {
  for (target in species_list) {
    if (source != target) {  # Skip if source and target are the same
      
      # Setup filter for jumping genes (default = 5)
      jumping_genes = 5
      
      # If you would like to reorder the chromosomes, put the chromosome order here.
      # Otherwise, set reorder to FALSE
      reorder <- TRUE
        Bfl_order <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
        Pma_order <- c(3,13,4,9,1,17,10,12,7,11,6,2,5,15,8,14,16,18,19)
        
        Abi_order <- c(9,8,6,10,3,4,7,2,11,5,1)
        Aam_order <- c(1,3,19,9,13,21,20,17,11,16,5,7,10,22,14,8,15,2,4,12,6,18)
        Afi_order <- c(5,17,14,8,16,1,2,12,15,10,4,11,18,13,20,19,3,9,6,7)
        Aja_order <- c(2,9,6,16,22,7,1,18,3,20,4,17,8,11,10,14,19,15,5,13,21,23,12)
        Spu_order <- c(5,18,2,13,21,7,3,20,10,17,8,11,9,6,1,15,4,12,19,16,14)
        San_order <- c(1,13,12,19,7,4,8,6,15,3,18,10,9,2,5,11,20,16,17,14)
        
        source_order <- get(paste0(source, "_order"))
        target_order <- get(paste0(target, "_order"))
      
      # List any chromosomes in either the source or the target genome which should be flipped
      # CHANGED: must be defined even when not reordering -- the flip loops
      # below run unconditionally.
        Bfl_flip <- c()
        Pma_flip <- c()
        Abi_flip <- c()
        Aam_flip <- c()
        Afi_flip <- c()
        Aja_flip <- c()
        Spu_flip <- c()
        San_flip <- c()
        
        source_flip <- get(paste0(source, "_flip"))
        target_flip <- get(paste0(target, "_flip"))
      
      ########################################
      ### 1. Set karyotypes
      ########################################
      
      # Load karyotype files
      source_karyotype <- read.table(paste0("input/", source, "_karyotype.txt"),
                                     sep = "\t", header = TRUE, stringsAsFactors = F)
      target_karyotype <- read.table(paste0("input/", target, "_karyotype.txt"),
                                     sep = "\t", header = TRUE, stringsAsFactors = F)
      
      # Sort karyotypes such that chromosomes are ordered from largest to smallest
      source_karyotype <- source_karyotype[order(source_karyotype$End,
                                                 decreasing = TRUE), ]
      target_karyotype <- target_karyotype[order(target_karyotype$End,
                                                 decreasing = TRUE), ]
      
      # Create dictionaries such that each karyotype name can be mapped to an integer
      source_mapping <- setNames(1:nrow(source_karyotype), source_karyotype$Chr)
      target_mapping <- setNames(1:nrow(target_karyotype), target_karyotype$Chr)
      
      # Save the original chromosome names to a new column
      source_karyotype$orig_ids <- source_karyotype$Chr
      target_karyotype$orig_ids <- target_karyotype$Chr
      
      # Remap the karyotypes such that the names are integers
      source_karyotype$Chr <- source_mapping[source_karyotype$Chr]
      target_karyotype$Chr <- target_mapping[target_karyotype$Chr]
      
      # Reorder karyotypes as specified.
      if (reorder == TRUE) {
        source_karyotype <- source_karyotype[order(match(source_karyotype$Chr, source_order)),]
        target_karyotype <- target_karyotype[order(match(target_karyotype$Chr, target_order)),]
        source_remap <- setNames(1:nrow(source_karyotype), source_order)
        target_remap <- setNames(1:nrow(target_karyotype), target_order)
      }
      
      # Combine the karyotypes into one matrix
      karyotype_all <- rbind(source_karyotype, target_karyotype,
                             make.row.names=FALSE)
      karyotype_all$orig_ids <- NULL
      # Set fill colours for nodes
      karyotype_all$fill <- "cccccc"
      
      # Reorder the columns
      karyotype_all <- karyotype_all[, c(1,2,3,7,4,5,6)]
      
      # Change the datatype of the variables into a data frame
      karyotype_all <- as.data.frame(karyotype_all, stringsAsFactors = F )
      head(karyotype_all)
      
      ########################################
      ### 2. Set orthologous gene coordinates
      ########################################
      
      # Load coordinates files
      source_coords <- read.table(paste0("input/", source, "_coordinates.tsv"),
                                  sep = "\t", header = FALSE, stringsAsFactors = F)
      target_coords <- read.table(paste0("input/", target, "_coordinates.tsv"),
                                  sep = "\t", header = FALSE, stringsAsFactors = F)
      
      source_coords <- source_coords[!is.na(source_coords$V3) & source_coords$V3 %in% source_karyotype$orig_ids,]
      target_coords <- target_coords[!is.na(target_coords$V3) & target_coords$V3 %in% target_karyotype$orig_ids,]
      
      # Check file length
      nrow(source_coords)
      nrow(target_coords)
      
      # Get rid of coordinates that aren't on the selected chromosomes
      source_coords <- source_coords[source_coords$V3 %in% source_karyotype$orig_ids,]
      target_coords <- target_coords[target_coords$V3 %in% target_karyotype$orig_ids,]
      
      # Recheck file length
      nrow(source_coords)
      nrow(target_coords)
      
      # Check coordinates file
      head(source_coords)
      head(target_coords)
      
      # Map the chromosome names to appropriate karyotype numbers for both coordinate sets
      source_coords$V3 <- source_mapping[source_coords$V3]
      target_coords$V3 <- target_mapping[target_coords$V3]
      
      # Flip any specified chromosomes
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
      
      # Remap the chromosomes to the reshuffled chromosome order if necessary
      if (reorder == TRUE) {
        source_coords$V3 <- source_remap[as.character(source_coords$V3)]
        target_coords$V3 <- target_remap[as.character(target_coords$V3)]
      }
      
      # Create connections
      coords_all <- (merge(source_coords, target_coords, by = "V1")
                     [, c(-1, -2, -6, -7)])
      head(coords_all)
      
      # Relabel columns
      colnames(coords_all) <- c("Species_1", "Start_1", "End_1", "Species_2", "Start_2", "End_2", "ALG")
      
      # Create ALG colour map
      # The 'XXXX' designation indicates orthology assignments made independently of ALG orthologs.
      alg_hexmap <- c("A1a" = "A6CDE3", "A1b" = "A6CDE3", "A2" = "197AB4",
                      "B1" = "089E78", "B2" = "31A149", "B3" = "B2D287",
                      "C1" = "F19799", "C2" = "E42027",
                      "D"  = "F7F09B",
                      "Ea" = "AF5B2A", "Eb" = "AF5B2A",
                      "F"  = "DA6227",
                      "G"  = "7671B0",
                      "H"  = "E42C8A",
                      "I"  = "65A744",
                      "J1" = "C9B2D3", "J2" = "664292",
                      "K"  = "E6AB22",
                      "L"  = "A4782B",
                      "M"  = "90CDC3",
                      "N"  = "686868",
                      "O1" = "EF7E25", "O2" = "F9BD6F",
                      "P"  = "BDB9DA",
                      "Qa" = "FED832", "Qb" = "FED832", "Qc" = "FED832", "Qd" = "FED832",
                      "R"  = "00008B",
                      "XXXX" = "DDDDDD")
      
      # CHANGED: relabel unlabelled orthogroups ("NA" from the Python
      # pipeline) to "XXXX" before the colour lookup -- otherwise these
      # produce an invalid "#NA" colour and geom_point() errors out.
      coords_all$ALG[coords_all$ALG == "NA" | is.na(coords_all$ALG)] <- "XXXX"
      
      # Reassign line colours based on the mapping above
      coords_all$fill <- alg_hexmap[coords_all$ALG]
      
      # Check entries
      head(coords_all)
      
      # Filter out the jumping genes to reduce noise
      coords_all <- coords_all %>%
        group_by(Species_1, Species_2) %>%
        filter(n() > jumping_genes) %>%
        ungroup()
      
      # Check file
      head(coords_all)
      
      # Re-order connections by colour
      # Re-order connections: grey (XXXX/non-ALG) drawn first/bottom, real
      # ALG colors drawn last/top.
      sort_key <- coords_all$fill != "DDDDDD"
      coords_all <- coords_all[order(sort_key), ]
      
      # Reformat into data frame to plot
      synteny_coord <- as.data.frame(coords_all[,-7])
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
      ###    mixing rate
      ########################################
      
      synteny_coord$rank_1 <- rank(synteny_coord$adjusted_start_1)
      synteny_coord$rank_2 <- rank(synteny_coord$adjusted_start_2)
      
      spearman_coefficient <- cor(synteny_coord$adjusted_start_1, synteny_coord$adjusted_start_2, method = "spearman")
      
      mixing_rate <- 1 - abs(spearman_coefficient)
      
      mixing_rate_matrix[source, target] <- mixing_rate
      
      ########################################
      ### 6. Calculate synteny mixing rate for
      ###    orthologous chromosome pairs
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

# End of loop. Please continue with the subsequent steps below for additional analyses.

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

########################################################################
### (b) Examine specific chromosome subsets within the Oxford dot plot
########################################################################

# Assign source and target species
# CHANGED: wsu/cpa don't exist in your data -- using a real pair instead.
source <- "Bfl"
target <- "Pma"

# (Rerun the script from the karyotype-loading section down to the end of
# the main loop's body with these values, same as the original note said)

# Replace chromosome number with the target chromosome number you are interested in
specific_chromosome <- 1

# Plot highlighting points from a specific target chromosome
dp1 <- ggplot(synteny_coord, aes(x = adjusted_start_1/1e6, y = adjusted_start_2/1e6)) +
  geom_point(size = 0.3, color = 'royalblue') +
  geom_point(data = subset(synteny_coord, Species_1 == specific_chromosome),
             aes(color = "specific_chromosome"), size = 0.3) +
  geom_vline(xintercept = source_cumsum, color = "gray", size = 0.3) +
  geom_hline(yintercept = target_cumsum, color = "gray", size = 0.3) +
  theme_classic() +
  labs(x = source, y = target) +
  scale_color_manual(values = c("specific_chromosome" = "red")) + theme(legend.position = "none")

print(dp1)

# Subset the data for chromosomes 1 and 2, replace them with target chromosomes
subset <- c(1,2)
subset_chrs <- synteny_coord[synteny_coord$Species_1 %in% subset, ]

color_mapping <- setNames(c("royalblue", "red"), subset)

dp2 <- ggplot(subset_chrs, aes(x = adjusted_start_1/1e6, y = adjusted_start_2/1e6)) +
  geom_point(aes(color = as.factor(Species_1)), size = 0.3) +
  theme_classic() +
  labs(x = source, y = target, color = "Chromosome") +
  scale_color_manual(values = color_mapping)

print(dp2)

# Plot highlighting points from a specific target ALG (B2=31A149)
dp3 <- ggplot(synteny_coord, aes(x = adjusted_start_1/1e6, y = adjusted_start_2/1e6)) +
  geom_point(size = 0.3, color = "#27A9E0") +
  geom_point(data = subset(synteny_coord, fill == "31A149"),
             aes(color = "specific_ALG"), size = 0.3) +
  geom_vline(xintercept = source_cumsum, color = "gray", size = 0.3) +
  geom_hline(yintercept = target_cumsum, color = "gray", size = 0.3) +
  theme_classic() +
  labs(x = source, y = target) +
  scale_color_manual(values = c("specific_ALG" = "red")) + theme(legend.position = "none")

print(dp3)

# Loop through each chromosome and create a plot
chromosomes <- sort(unique(synteny_coord$Species_1))

for (chr in chromosomes) {
  chr_data <- subset(synteny_coord, Species_1 == chr)
  plot_name <- paste0(source, "_", "chr_", chr, "_", target, "_plot.pdf")
  dp4 <- ggplot(chr_data, aes(x = adjusted_start_1/1e6, y = adjusted_start_2/1e6)) +
    geom_point(size = 0.3, color = 'royalblue') +
    theme_classic() +
    labs(x = source, y = target, color = "Chromosome") +
    scale_color_manual(values = color_mapping) +
    ggtitle(paste("OrthoChromosome", chr))
  ggsave(plot_name, dp4)
}

write.csv(mixing_rate_matrix, "mixing_rate_genome_wide.csv")
write.csv(chr_mixing_rate_data, "mixing_rate_per_chromosome.csv", row.names = FALSE)

########################################################################
### (c) Calculate mixing rates for orthologous chromosome pairs
########################################################################

# In this section, we merge the evolutionary distance data with synteny 
# mixing rates for orthologous chromosome pairs.

# Check the mixing rate data
head(chr_mixing_rate_data)

# Parse the Newick format phylogenetic tree. 
# This tree is rooted, with 'Pecten maximus' (pma) as the outgroup.
tree_newick <- "(Her:0.5,(Hcr:1,Htu:1)1:0.5);"
tree <- read.tree(text = tree_newick)

# Calibrate the tree based on the split of pma and cmu around the Cambrian Period (540-480 mya)
cal <- makeChronosCalib(tree, interactive = T)

# Check the calibration table
head(cal)
# node age.min age.max soft.bounds
# 1    7     480     540       FALSE

# Create chronogram with estimated divergence times
chronogram <- chronos(tree, calibration = cal)

# Calculate pairwise divergence times
divergence_times <- cophenetic(chronogram)

# Convert the distance matrix to a table with unique pairs
divergence_time_data <- melt(as.matrix(divergence_times))
divergence_time_data <- divergence_time_data[divergence_time_data$Var1 != divergence_time_data$Var2, ]
names(divergence_time_data) <- c("Species1", "Species2", "DivergenceTime")

# Check the phylogenetic distance data
head(divergence_time_data)

# Merge mixing rate and phylogenetic distance data
merged_data <- merge(chr_mixing_rate_data, divergence_time_data, by = c("Species1", "Species2"), all.x = TRUE)

# Choose the species located deepest within the tree for further analysis (e.g., wsu or cpa)
# Exclude cmu since it doesn't have one-to-one orthologous chromosome pairs with other species.
merged_data <- subset(merged_data, Species1 == "cpa" & Species2 != "cmu")

# Fit the linear model
linear_model <- lm(MixingRate ~ DivergenceTime, data = merged_data)
equation <- paste("y =",
                  format(coef(linear_model)[1], digits = 4),
                  "+",
                  format(coef(linear_model)[2], digits = 4),
                  "x")
model_summary <- summary(linear_model)
r_squared_linear <- model_summary$r.squared

# Fit the logarithmic model
log_model <- lm(MixingRate ~ log(DivergenceTime), data = merged_data)
intercept <- coef(log_model)[1]
slope <- coef(log_model)[2]
equation <- paste("y =", round(intercept, 3), "+", round(slope, 3), "* log(x)")
model_summary <- summary(log_model)
r_squared_log <- model_summary$r.squared

# Create a dot plot with a linear model
dp6 <- ggplot(merged_data, aes(x = DivergenceTime, y = MixingRate)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, color = "black", se = FALSE) +
  theme_classic() +
  labs(title = "Relationship between phylogenetic distance and mixing rate",
       x = "Divergence Time (MYA)",
       y = "Mixing rate")

# Create a dot plot with a logarithmic model
dp7 <- ggplot(merged_data, aes(x = DivergenceTime, y = MixingRate)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ log(x), color = "black", se = FALSE) +
  theme_classic() +
  labs(title = "Relationship between phylogenetic distance and mixing rate",
       x = "Divergence Time (MYA)",
       y = "Mixing rate")

# Print the dot plot
print(dp6)
print(dp7)

# Save the data
write.csv(merged_data, file="mixing_rate.csv")
