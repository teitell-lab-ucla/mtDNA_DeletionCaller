# mtDNA DeletionCaller

A tool for detection of small deletions in mitochondrial DNA from single-cell sequencing data.

## Overview

mtDNA DeletionCaller is designed to identify potential deletions in mitochondrial DNA from single-cell sequencing data processed with the Mitochondrial Genome Analysis Toolkit (mgatk). This tool analyzes coverage patterns to identify significant drops that may indicate deletion events, and quantifies heteroplasmy levels of these events across cells.

## Features

- Detection of potential mtDNA deletions based on significant coverage drops
- Calculation of deletion heteroplasmy levels per cell
- Integration with mgatk output files
- Customizable minimum read count thresholds

## Requirements

- R (version 3.5+)
- Required R packages:
  
  - data.table (for efficient data loading with fread)
  - Matrix (for sparse coverage/allele matrix operations)
  - S4 object system support (included in base R) — the mgatk RDS is a SummarizedExperiment-style object accessed directly via its slots, so no additional Bioconductor packages need to be attached

## Input Files

mtDNA IndelCaller requires the following input files:

1. **mgatk output RDS file**: Contains coverage and allele count information
2. **Cell heteroplasmy file**: TSV/CSV containing heteroplasmy values per cell
3. **mtDNA reference file**: Text file with mitochondrial reference alleles

## Usage

```R
source("DeletionCaller.R")

results <- DeletionDetection(
  sample = "your_sample_name",
  o_call = "tenx_call",
  min_readcount = 1,
  file = "path/to/mgatk_output.rds",
  cell_heteroplasmy_file = "path/to/cell_heteroplasmic_df.tsv.gz",
  mtDNA_file = "path/to/chrM_refAllele.txt"
)
```

### Run Time 

mtDNA DeletionCaller will take approximately 10-15 mins to run. 

### Output

The function returns a list containing two data frames:

1. `FullData`: Heteroplasmy levels for all variants (SNVs and indels) across all cells
2. `variants`: Heteroplasmy levels for the newly detected deletion variants

## How It Works

1. The tool reads in coverage data from the mgatk output
2. It calculates differences in coverage between adjacent positions
3. It identifies patterns where one position has a significant drop in coverage (negative z-score > 3.5) followed by an adjacent position with a significant increase in coverage (positive z-score > 3.5)
4. These patterns indicate potential deletion boundaries where:
   - The first position marks where coverage suddenly drops (start of deletion)
   - The adjacent position marks where coverage suddenly increases again (end of deletion)
5. For identified deletion sites, it:
   - Calculates expected coverage based on surrounding positions (-5 to -1 and +1 to +5 around the site)
   - Estimates reference allele coverage at the position
   - Calculates deletion heteroplasmy as (expected - observed)/expected

The algorithm specifically looks for adjacent positions with significant opposite-signed z-scores. Based on the implementation, the tool is optimized for detecting small deletions, as it focuses on closely adjacent positions rather than scanning for larger ranges. The current implementation is most effective at detecting deletions between adjacent bases, but could potentially detect slightly larger deletions if they create the characteristic coverage pattern at their boundaries.

## Notes and Limitations

- The tool is primarily designed to detect small deletions between adjacent positions in mtDNA
- The current implementation detects deletions by identifying adjacent positions where coverage suddenly drops and then increases
- The algorithm uses a specific pattern of coverage changes to identify deletion events
- Larger deletions spanning multiple bases may be detected if they create the characteristic coverage pattern at their boundaries
- Very large deletions might be missed or improperly characterized if they don't fit the detection pattern
- Performance depends on sequencing depth and quality
- Results should be validated with alternative methods for critical analyses
- The tool assumes that the input data has been properly processed using mgatk
- The tool does not currently detect insertions, which would require different detection algorithms

## Citation

If you use this tool in your research, please cite:
[Paper Info]


## Contact

For questions or support, please contact [bpeng01@g.ucla.edu/eainthoneyaungwin@mednet.ucla.edu]
