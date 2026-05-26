#' Simulate mutations based on a reference
#'
#' Simulates single base mutations based on a reference genome.
#'
#' Simulation is based on uniform distribution of mutations across the full length of the geome.
#' VAFs are normally distributed around 0.5
#' Type column is based on random sampling aginst
#'
#' @param reference a 2-column data.frame describing the name and length of each reference sequence.
#' Can create from .fai index file using [fai_to_reference()]
#' @param col_name name of column in reference data.frame that describes sequence name
#' @param col_length name of column in reference data.frame that describes sequence length
#' @param nmuts number of mutations to simulate
#' @param types a vector of possible mutation 'types'
#' @param type_probabilty a vector describing the probabiliy of assigning a mutation to each type.
#' Should be the same length as types vector.
#'
#'
#' @returns A mutations data.frame compatible with the [add_mutations()] function. Contains the columns:
#' 1. chromosome
#' 2. position
#' 3. vaf
#' 4. type
#' @export
#'
#' @examples
#' ## Get Reference Based on fai index
#' fai <- system.file("hg38.fai", package = "ggcircos")
#' hg38 <- fai_to_reference(fai)
#'
#' # Simulate 100 mutations
#' mutations <- simulate_mutations(hg38, nmuts = 100)
#'
#'
#' Plot simulated mutations
#' multiomic_circos <- ggcircos(hg38) |>
#'   add_mutations(mutations)
simulate_mutations <- function(reference, col_name = "chromosome", col_length = "length", nmuts = 100, types = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"), type_probability = rep(1/length(types), times = length(types))){

  # Assertions
  assertions::assert_no_duplicates(types)
  assertions::assert_number(nmuts)
  assertions::assert_equal(length(types), length(type_probability))

  # Simulate mutations from reference (uniformly distributed)
  l <- as.numeric(reference[[col_length]])
  r <- reference[[col_name]]

  # Total length including all chroms
  total_l <- sum(l)
  prop_l <- l/total_l
  cum_l <- cumsum(l) # Cumulative total length
  names(cum_l) <- r
  n_to_subtract <- dplyr::lag(cum_l, n=1, default = 0)

  # Mut positions (where position is on global coord)
  sim_global_positions <- round(runif(n = nmuts, min = 1, max = total_l), digits = 0)

  # Map back to chr + local positions
  sim_chr <- cut(sim_global_positions, breaks = c(1, cum_l), labels = r, right = TRUE, include.lowest = TRUE)
  sim_chr <- as.character(sim_chr)
  sim_position <- sim_global_positions - n_to_subtract[match(sim_chr, r)]
  unname(sim_position)

  # Get Vaf
  sim_vafs <- rnorm(mean = 0.5, n = nmuts, sd = 0.1)

  # Include random ref and alt
  type <- sample(types, replace = TRUE, size = nmuts, prob = type_probability)

  # Create Output Dataframe
  muts_unsorted <- data.frame(
    "chromosome" =  sim_chr,
    "position" = sim_position, #1-based
    "vaf" = sim_vafs,
    "type" = type
  )

  # Sort based on reference
  chr_order <- match(sim_chr, r)
  muts_unsorted[order(chr_order, muts_unsorted$position),,drop=FALSE]
}
