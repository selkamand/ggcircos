
# Exported Plot builders --------------------------------------------------

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
#' ## Filepaths
#' fai <- system.file("hg38.fai", package = "ggcircos")
#' mutations_tsv <- system.file("hg38.mutations.fai", package = "ggcircos")
#'
#' # Read in data
#' hg38 <- fai_to_reference(fai)
#' mutations <- read.csv(mutations_tsv, sep = "\t")
#'
#' # Create Circos Plot
#' multiomic_circos <- ggcircos(hg38) |>
#'   add_mutations(mutations)
ggcircos <- function(reference, col_name = "chromosome", col_length = "length", inter_chrom_spacing = 0.1) {
  # Perform assertions and create a 4 column data.frame with columns: name, start, end, refindex
  reference <- standardise_reference_regions(reference, col_name, col_length)

  # Output basic ggcircos object
  Circos(
    reference = reference,
    inter_chrom_spacing = inter_chrom_spacing
  )
}

#' Add mutations to circos plot
#'
#' @param circos a circos object created using `ggcircos()`
#' @param data a data.frame where each row represents a mutations.
#' @param col_chrom Name of column describing the chromosome of the mutation
#' @param col_pos Name of column describing the 1-based position of the mutation.
#' @param col_type Name of column describing the type of mutation. Used to colour mutations. If NULL mutations will not be coloured.
#' @param col_vaf Name of column
#'
#' @returns a Circos class object with mutations added
#' @export
#'
#' @examples
#' ## Filepaths
#' fai <- system.file("hg38.fai", package = "ggcircos")
#' mutations_tsv <- system.file("hg38.mutations.fai", package = "ggcircos")
#'
#' # Read in data
#' hg38 <- fai_to_reference(fai)
#' mutations <- read.csv(mutations_tsv, sep = "\t")
#'
#' # Create Circos Plot
#' multiomic_circos <- ggcircos(hg38) |>
#'   add_mutations(mutations)
add_mutations <- function(circos, mutations, col_chrom = "chromosome", col_pos = "position", col_vaf = "vaf", col_type = NULL, col_tooltip = NULL){

  # Assertions
  assertions::assert_class(circos, class = "ggcircos::Circos")
  assertions::assert_dataframe(mutations)
  assertions::assert_names_include(mutations, c(col_chrom, col_pos, col_vaf, col_type, col_tooltip))

  if(!is.null(col_type) || !is.null(col_tooltip)) {
    stop("col_type and col_tooltip have not been implemented yet")
  }

  # Standardise mutation data.frame into format for loading into the 'Track' class
  # Note we leave optional columns col_type and col_tooltip named whatever they were before for now
  rename_map <- list(
    name = col_chrom,
    position = col_pos,
    vaf = col_vaf
  )

  mutations = dfx::brename(mutations, rename_map)

  # Create Track object (at present col_type and col_tooltip are ignored)
  mutation_track <- Track(data = mutations, type = "SmallMutations")

  # Add mutation track to circos
  circos@tracks <- append(circos@tracks, list(mutation_track))

  return(circos)
}

generate <- function(circos){
  assertions::assert_class(circos, "ggcircos::Circos")

  # Grab reference Data
  refdata <- circos@reference_normalised_data

  # Count tracks (rings we need to plot)
  ntracks <- circos@n_tracks


  #smallmutdata <- circos@snv_normalised_data

  ggcirc <- ggplot2::ggplot(refdata) +
    ggiraph::geom_rect_interactive(
      # TODO change ymin and ymax based on other data in plot
      ggplot2::aes(
        xmin = start, xmax = end, ymin = 90, ymax = 100, tooltip = name,
        data_id = name
        ),
      color="black", fill = "white"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = start + (end-start)/2, label = name, y = 115)
      ) +
    ggplot2::scale_y_continuous(limits = c(1, NA)) +
    ggplot2::scale_x_continuous() +
    ggplot2::ylab(NULL) +
    ggplot2::xlab(NULL) +
    ggplot2::coord_radial(expand = FALSE) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank()
    )

  # Add SNV ring
  if(ntracks > 1){
    muts <- circos@tracks[[2]]@data
    muts$name <- paste0("chr", muts$name)
    #browser()
    muts$position_normalised <- purrr::map2_dbl(.x = muts$name, .y = muts$position, .f = circos@project_point)
    #browser()
    ggcirc <- ggcirc +
      ggiraph::geom_point_interactive(
        data = muts,
        ggplot2::aes(
          x = .data[["position_normalised"]],
          tooltip = paste0(.data[["name"]], ":", .data[["position"]]),
          data_id = paste0(.data[["name"]], ":", .data[["position"]])
          ),
        shape=1,
        size=0.5,
        y = 60,
      )
  }

  return(ggcirc)

}



make_interactive <- function(ggcircos, ...){
  ggiraph::girafe(ggcircos, ...)
}
# chr_to_start_and_end_modifier <- function(stdref, inter_chrom_spacing){
#
#
#   reflength <- sum(stdref$end)
#   nseqs <- nrow(stdref)
#   chroms <- stdref$name
#   chromlengths <- stdref$end
#
#   if(nseqs == 1) inter_chrom_spacing <- 0
#
#   total_spaces_in_bases <- inter_chrom_spacing * reflength
#   space_per_gap <- total_spaces_in_bases / (nseqs-1)
#
#
#   cumlength <- cumsum(stdref$end)
#
#   # Create gapless xmin & xmax
#   xmin_nogap <- stdref$start + dplyr::lag(cumlength, default = 0)
#   xmax_nogap <- xmin_nogap + chromlengths
#
#   # Slide min and max both along by 'gap' * refindex
#   gapvec <- space_per_gap * stdref$refindex-1
#   xmin_with_gap <- xmin_nogap + gapvec
#   xmax_with_gap <- xmax_nogap + gapvec
#
#   # Normalise so values are between 1 and 100
#   min_x <- min(xmin_with_gap)
#   max_x <- max(xmax_with_gap)
#   start_norm <- scales::rescale(xmin_with_gap, from = c(min_x, max_x), to = c(1, 100))
#   end_norm <- scales::rescale(xmax_with_gap, from = c(min_x, max_x), to = c(1, 100))
#
#   chrommap = lapply(seq_along(start_norm), FUN = function(i){
#     list(
#       origina
#       new = c("start_normalised" = start_norm, "start_normalised" = end_norm)
#     )
#   })
#   names(chrommap) <- chroms
#   return(chrommap)
#   #
#   # # Refdata
#   # newdata <- data.frame(
#   #   name = chroms,
#   #   start = stdref$start,
#   #   end = stdref$end,
#   #   start_normalised = start_norm,
#   #   end_normalised = end_norm
#   #   #start_normalised = start_norm,
#   #   #end_normalised = end_norm
#   # )
#   #
#   # list(
#   #   data = newdata,
#   #   rescale_point <- function(name, position){
#   #     idx = match(data$name, name)
#   #     range = data$
#   #   }
#   # )
#   # }
# }

get_normfactor <- function(x, range = c(1, 100)){
  to[1] + (x - old_min) * (to[2] - to[1]) / (old_max - old_min)
}

rescale_minmax <- function(x, to = c(1, 100), na.rm = TRUE) {
  old_min <- min(x, na.rm = na.rm)
  old_max <- max(x, na.rm = na.rm)

  if (old_min == old_max) {
    return(rep(mean(to), length(x)))
  }

  to[1] + (x - old_min) * (to[2] - to[1]) / (old_max - old_min)
}

#' Parse fai index to reference data.frame
#'
#' Parse fai index to reference data.frame expected by ggcircos
#'
#' @param fai Path to a .fai index file created by running `samtools faidx <reference.fasta>`.
#' @param exclude Names of chromosomes to exclude.
#' Many references will include decoy contigs that you probably don't want plotted in your circos.
#' It is also common to exclude mitochondrial contigs ('chrM' / 'MT') in eukaryotic genomes since they are too small
#' compared to nuclear chromosomes to see well in circos.
#' @returns A 2-column data.frame with fields 'chromosome' and 'length' that can be parsed by ggcircos
#' @export
#'
#' @examples
#' fai <- system.file("hg38.fai", package = "ggcircos")
#' fai_to_reference(fai)
fai_to_reference <- function(fai, exclude = NULL) {
  assertions::assert_file_exists(fai)

  # Read Data
  df_fai <- read.csv(fai, header = FALSE, sep = "\t")

  # Assert its a real fai file.
  assertions::assert(ncol(df_fai) == 5, msg = paste0("fai indexes are expected to have 5 columns. Are you sure [", fai, "] is an index file created by `samtools faidx` from a reference fasta file?"))

  # Set colnames
  colnames(df_fai) <- c("chromosome", "length", "offset", "linebases", "linewidth")

  # Just return chromosome and length
  reference <- dfx::bselect(df_fai, columns = c("chromosome", "length"))

  # Turn length into 'numeric'
  # (since if parsed as 'integer' type you can end up with integer overflow when we cumsum)
  reference$length <- as.numeric(reference$length)

  # Exclude reference sequences with names specified in `exclude`
  if(!is.null(exclude)){
    rows_to_include <- !reference$chromosome %in% exclude
    reference <- reference[rows_to_include,,drop=FALSE]
  }

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
#' df <- data.frame(chromosome = c("chr1", "chr2"), length = c(1000, 2000))
#' reference <- standardise_reference_regions(df)
#' print(reference)
#'
standardise_reference_regions <- function(reference, col_name = "chromosome", col_length = "length") {
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
    )
  )

  # Add start coord
  reference$start <- 1

  # Only return columns of interest
  reference <- dfx::bselect(reference, c("name", "start", "end"))
  reference$refindex <- seq_along(reference$name)

  # Ensure data.types are sensible
  reference$start <- as.numeric(reference$start)
  reference$end <- as.numeric(reference$end)

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
transform_position_to_reference <- function(reference, inter_chrom_gap = 0) {
  reference$length <- reference$end

  total_length <- sum(reference$length)
  buffer_size <- round(total_length * inter_chrom_gap, digits = 0)

  reference$realstart <- reference$start
}
