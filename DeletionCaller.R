## file = path to rds object output from mgatk
#cell_heteroplasmy_file = path to heteroplasmy csv, cell_heteroplasmic_df.tsv.gz
# mtDNA_file = path to mtDNA reference , chrM_refAllele.txt

DeletionDetection <- function(sample, o_call = 'tenx_call', min_readcount = 1, 
                              file, cell_heteroplasmy_file, mtDNA_file) {
  library(data.table)
  library(Matrix)
  
  # Load the object and variant table
  object <- readRDS(file)
  variants <- fread(cell_heteroplasmy_file)
  variants <- as.data.frame(variants)
  rownames(variants) <- variants$V1
  variants <- variants[, -1, drop = FALSE]
  
  # Handle coverage (sparse matrix)
  coverage_matrix <- object@assays@data$coverage  # dgCMatrix
  coverage <- Matrix::rowSums(coverage_matrix)
  
  coverage_df <- coverage_matrix
  colnames(coverage_df) <- object@colData$sample
  
  # Load mtDNA reference
  mtdna <- read.csv(mtDNA_file, header = FALSE, sep = '\t')
  
  # Detect potential deletions based on coverage changes
  diff <- coverage[2:length(coverage)] - coverage[1:(length(coverage) - 1)]
  zscores <- (diff - mean(diff)) / sd(diff)
  potential <- which(abs(zscores) > 3.5)
  
  opp_sign <- closeness <- c()
  if (length(potential) > 1) {
    for (i in 1:(length(potential) - 1)) {
      closeness[i] <- potential[i + 1] - potential[i]
      opp_sign[i] <- (zscores[potential[i]] < 0) & (zscores[potential[i + 1]] > 0)
    }
    evaluation <- which(closeness == 1 & opp_sign)
  } else {
    evaluation <- integer(0)
  }
  
  if (length(evaluation) > 0) {
    cat('Potential deletions found\n')
    potential <- potential[evaluation]
    
    for (i in seq_along(potential)) {
      position <- potential[i] + 1
      ref_allele <- toupper(mtdna[position, 2])
      assays <- paste0(ref_allele, c('_counts_fw', '_counts_rev'))
      
      # Estimate cell coverage
      positions_idx <- c(max(position - 5, 1):(position - 1),
                         (position + 1):min(position + 5, nrow(coverage_df)))
      est_cov <- Matrix::colMeans(coverage_df[positions_idx, rownames(variants), drop = FALSE])
      
      # Estimate reference allele coverage
      ref_allele_coverage <- rep(0, length(est_cov))
      for (k in seq_along(assays)) {
        if (assays[k] %in% names(object@assays@data)) {
          x <- object@assays@data[[assays[k]]]
          colnames(x) <- object@colData$sample
          ref_allele_coverage <- ref_allele_coverage + x[position, rownames(variants)]
        }
      }
      
      del_coverage <- (est_cov - ref_allele_coverage) / est_cov
      del_coverage[(!est_cov > min_readcount) | est_cov == 0] <- NA
      del_coverage[del_coverage < 0] <- 1
      name <- paste0(position, ref_allele, '>NULL')
      variants[name] <- del_coverage
    }
    cat('Done with adding deletions\n')
  }
  
  cat('Adding full data\n')
  FullData <- data.frame(row.names = object@colData$sample)
  
  for (i in 1:ncol(variants)) {
    name <- colnames(variants)[i]
    split_name <- strsplit(name, '')[[1]]
    loc <- which(split_name == '>')
    ref <- split_name[loc - 1]
    alt <- split_name[loc + 1]
    position <- as.numeric(paste0(split_name[1:(loc - 2)], collapse = ''))
    
    if (alt == 'N') {
      # Deletion
      assays <- paste0(ref, c('_counts_fw', '_counts_rev'))
      positions_idx <- c(max(position - 5, 1):(position - 1),
                         (position + 1):min(position + 5, nrow(coverage_df)))
      est_cov <- Matrix::colMeans(coverage_df[positions_idx, rownames(FullData), drop = FALSE])
      
      ref_allele_coverage <- rep(0, length(est_cov))
      for (k in seq_along(assays)) {
        if (assays[k] %in% names(object@assays@data)) {
          x <- object@assays@data[[assays[k]]]
          colnames(x) <- object@colData$sample
          ref_allele_coverage <- ref_allele_coverage + x[position, rownames(FullData)]
        }
      }
      
      del_coverage <- (est_cov - ref_allele_coverage) / est_cov
      del_coverage[(!est_cov > min_readcount) | est_cov == 0] <- NA
      del_coverage[del_coverage < 0] <- 1
      FullData[name] <- del_coverage
    } else {
      # SNP
      assays <- paste0(alt, c('_counts_fw', '_counts_rev'))
      est_cov <- coverage_df[position, rownames(FullData)]
      
      ref_allele_coverage <- rep(0, length(est_cov))
      for (k in seq_along(assays)) {
        if (assays[k] %in% names(object@assays@data)) {
          x <- object@assays@data[[assays[k]]]
          colnames(x) <- object@colData$sample
          ref_allele_coverage <- ref_allele_coverage + x[position, rownames(FullData)]
        }
      }
      
      ref_allele_coverage <- ref_allele_coverage / est_cov
      ref_allele_coverage[(!est_cov > min_readcount) | est_cov == 0] <- NA
      FullData[name] <- ref_allele_coverage
    }
  }
  
  return(list(FullData = FullData, variants = variants))
}



# Run the original function
results <- DeletionDetection(
  sample = "samplename",
  o_call = "tenx_call",
  min_readcount = 1,
  file = "/Volumes/4T1_data/output/mgatk_out/p0/final/mgatk.rds",
  cell_heteroplasmy_file = "/Volumes/4T1_data/output/mgatk_out/p0/final/mgatk.cell_heteroplasmic_df.tsv.gz",
  mtDNA_file = "/Volumes/4T1_data/output/mgatk_out/p0/final/chrM_refAllele.txt"
)

write.csv(
  results[["variants"]],
  "/Volumes/4T1_data/Rdraft/FMI/FMI_Deletion_Heteroplasmy_010726_variants.csv",
  row.names = TRUE
)

