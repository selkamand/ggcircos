#' Create a base circos plot
#'
#' @param reference a 2-column data.frame describing the name and length of each reference sequence.
#' Can create from .fai index file using [fai_to_reference()]
#' @param col_name name of column in reference data.frame that describes sequence name
#' @param col_length name of column in reference data.frame that describes sequence length
#'
#'
#' @returns A ggplot canvas with karyotype information
#' @export
#'
#' @examples
ggcircos <- function(reference, col_name = "chromosome", col_length = "length"){

  # Perform assertions and create a 3 column data.frame with columns: name, start, end
  reference <- standardise_reference_regions(reference, col_name, col_length)


  # Output basic ggcircos object
  list(
    reference = reference
  )
}


#' Parse fai index to reference data.frame
#'
#' Parse fai index to reference data.frame expected by ggcircos
#'
#' @param fai Path to a .fai index file created by running `samtools faidx <reference.fasta>`.
#'
#' @returns A 2-column data.frame with fields 'chromosome' and 'length' that can be parsed by ggcircos
#' @export
#'
#' @examples
#' fai <- system.file("hg38.fai", package = "ggcircos")
#' fai_to_reference(fai)
fai_to_reference <- function(fai){
  assertions::assert_file_exists(fai)

  # Read Data
  df_fai <- read.csv(fai, header = FALSE, sep = "\t")

  # Assert its a real fai file.
  assertions::assert(ncol(df_fai) == 5, msg = paste0("fai indexes are expected to have 5 columns. Are you sure [", fai, "] is an index file created by `samtools faidx` from a reference fasta file?"))

  # Set colnames
  colnames(df_fai) <- c("chromosome", "length", "offset", "linebases", "linewidth")

  # Just return chromosome and length
  reference <- dfx::bselect(df_fai, columns = c("chromosome", "length"))

  # reference <- standardise_reference_regions(df_fai, col_name = "chromosome", col_length = "length")

  # Return Standardised Data
  return(reference)
}

# Standardisation -----------------------------------------------------------------

#' Standardise reference regions
#'
#' @param reference a 2-column data.frame describing the name and length of each reference sequence.
#' Can create from .fai index file using [fai_to_reference()]
#' @param col_name name of column in reference data.frame that describes sequence name
#' @param col_length name of column in reference data.frame that describes sequence length
#'
#' @returns 4-column data.frame with columns 'name, start, end, refindex'.
#' Intervals are 1-based & both-end inclusive (\[start, end\]).
#' This means a \[1, 1000\] reference entry is a 1000 bp-long region.
#' refindex is a unique number from 1 to nrow that is unique for each reference chromosome name and is used
#' by downstream functions to convert chrom-pos positions into a single coord position on x-axis. Based on order in input reference data.frame.
#'
#' @export
#'
#' @examples
#' df = data.frame(chromosome = c("chr1", "chr2"), length = c(1000, 2000))
#' reference <- standardise_reference_regions(df)
#' print(reference)
#'
standardise_reference_regions <- function(reference, col_name = "chromosome", col_length = "length"){

  # Assertions
  assertions::assert_dataframe(reference)
  assertions::assert_string(col_name)
  assertions::assert_string(col_name)
  assertions::assert_names_include(reference, names = c(col_name, col_length))
  assertions::assert_no_missing(reference[[col_name]])
  assertions::assert_no_missing(reference[[col_length]])
  assertions::assert_no_duplicates(reference[[col_name]])
  assertions::assert_length_greater_than(reference[[1]], length = 0, msg = "reference data.frame must contain at least one row")
  assertions::assert_numeric(reference[[col_length]])
  assertions::assert_all_greater_than(reference[[col_length]], minimum = 1)
  assertions::assert_character(reference[[col_name]])

  # Standardise names so reference data.frame has cols 'name, start, end'
  reference <- dfx::brename(
    reference,
    namemap = c(
      "end" = col_length,
      "name" = col_name
    ))

  # Add start coord
  reference$start <- 1

  # Only return columns of interest
  reference <- dfx::bselect(reference, c("name", "start", "end"))
  reference$index <- seq_along(reference$name)

  # Note at this point reference is 1-based, both-end inclusive format (a.k.a [start, end])
  # This means a [1, 1000] reference entry is a 1000 bp-long region
  return(reference)
}


#' Create Coordinate Transformer
#'
#' Once we have a reference karyotype set, we to transform chrom + pos genomic positions to a single x-coordinate position
#' that tells us where on the x axis to plot.
#' For example a variant at chr1:100 we might want to plot at position 100 on the x axis, whereas chr2:50 we might want to plot at
#' x = chr1_length + 50.
#' If we want our x coord axis to have gaps between each chromosome, we might need to take into account additional buffer space.

#' @param reference standardised reference data.frame created from [standardise_reference_regions()]
#' @param inter_chrom_gap proportion of total genome size that will be taken up by gaps between chromosomes
#'
#' @returns A function that takes two arguments 'chromosome' and 'position' and returns a single number representing position on the new x axis
#' @export
#'
#' @examples
transform_position_to_reference <- function(reference, inter_chrom_gap = 0){


  reference$length <- reference$end

  total_length <- sum(reference$length)
  buffer_size = round(total_length * inter_chrom_gap, digits = 0)

  reference$realstart <- reference$start
}
