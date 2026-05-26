#' Project points onto circos mega-interval
#'
#' Note target range is defined by \code{target_range} argument of \code{create_normalisation_map}
#'
#' @param chr chromosome name
#' @param pos position in original coord system (1 based)
#' @param map produced by [create_normalisation_map()]. A list mapping chromosomes to original & new ranges.
#'
#' @returns a single number representing position in unified (mega-interval) coord system.
project_points <- function(chr, pos, map){

  # Assertions
  if(length(chr) != length(pos)) stop("`chr` and `pos` must be the same length")
  if(!is.numeric(pos)) stop("`pos` must be numeric")

  chr_in_ref <- chr %in% names(map)
  newpos <- Map(chr, pos, !chr_in_ref, f = function(c, p, missing){
     if(missing) return(NA)
     scales::rescale(x = p, from = map[[c]]$original, map[[c]]$new)
  })
  unname(unlist(newpos))
}

#' Project intervals onto circos mega-interval
#'
#' Note target range is defined by \code{target_range} argument of \code{create_normalisation_map}
#'
#' @param chr chromosome name
#' @param start start position in original coord system (1 based)
#' @param end end position in original coord system (1 based)
#' @param map produced by [create_normalisation_map()]. A list mapping chromosomes to original & new ranges.
#' @param return_type Output format. Either `"list"` to return a list of 2-element numeric vectors, or `"data.frame"` to return a data frame with columns `chr`, `start`, and `end`.
#'
#' @returns If `return_type = "list"`, a list of 2-element numeric vectors representing projected start and end positions. If `return_type = "data.frame"`, a data frame with columns `chr`, `start`, and `end`.
project_intervals <- function(chr, start, end, map, return_type = c("list", "data.frame")){
  if(length(chr) != length(start)) stop("`chr` and `start` must be the same length")
  if(length(chr) != length(end)) stop("`chr` and `end` must be the same length")
  return_type <- match.arg(return_type, several.ok = FALSE)
  new_start <- project_points(chr = chr, pos = start, map = map)
  new_end <- project_points(chr = chr, pos = end, map = map)

  if(return_type == "list"){
    Map(new_start, new_end, f = function(s, e) { c(s, e)})
  }
  else{
    data.frame("chr" = chr, "start" = new_start, "end" = new_end)
  }
}

#' Create a data structure mapping chr + original positions to a new mega-interval
#'
#' Take the original start & end for each chromosome, line them up into a single mega interval (with \code{inter_chrom_spacing}).
#' Returns a data structure mapping each chromosome to a list of two ranges,
#' one representing the original range, the other representing the new range.
#' This makes it easy to rescale a chr + postion using [project_points()]
#'
#' @param stdref data.frame with columns 'name', 'start', 'end' where intervals are 1-based both-end inclusive.
#' @param inter_chrom_spacing how much space is between each chromosome
#' @param target_range c(1, 100)
#'
#' @returns a e.g.
#' list(
#'  chr1 = list(
#'     original = c(start=start, end=end),
#'     new = c(start=start, end=end),
#'  ),
#'  chr2 = list(
#'     original = c(start=start, end=end),
#'     new = c(start=start, end=end),
#'  )
#' )
#'
create_normalisation_map <- function(stdref, inter_chrom_spacing, target_range = c(1, 100)){

  # Assertions
  required_names <- c("name", "start", "end")
  observed_names <- colnames(stdref)
  missing_names <- setdiff(required_names, observed_names)
  if(length(missing_names) > 0) stop("`stdref` missing expected names [", toString(missing_names), "]")

  reflength <- sum(stdref$end)
  nseqs <- nrow(stdref)
  chroms <- stdref$name
  chromlengths <- stdref$end

  if(nseqs == 1) inter_chrom_spacing <- 0

  total_spaces_in_bases <- inter_chrom_spacing * reflength
  space_per_gap <- total_spaces_in_bases / (nseqs-1)
  if(nseqs == 1) space_per_gap <- 0

  cumlength <- cumsum(stdref$end)

  # Create gapless xmin & xmax
  xmin_nogap <- stdref$start + dfx::lag(cumlength, default = 0)
  xmax_nogap <- xmin_nogap + chromlengths

  # Slide min and max both along by 'gap' * refindex
  gapvec <- space_per_gap * (stdref$refindex-1)
  xmin_with_gap <- xmin_nogap + gapvec
  xmax_with_gap <- xmax_nogap + gapvec

  # Normalise so values are between 1 and 100
  min_x <- min(xmin_with_gap)
  max_x <- max(xmax_with_gap)

  start_norm <- scales::rescale(xmin_with_gap, from = c(min_x, max_x), to = target_range)
  end_norm <- scales::rescale(xmax_with_gap, from = c(min_x, max_x), to = target_range)

  chrommap = lapply(seq_along(start_norm), FUN = function(i){
    list(
      original = c("start" = stdref$start[i], "end" = stdref$end[i]),
      new = c("start" = start_norm[i], "end" = end_norm[i])
    )
  })
  names(chrommap) <- chroms
  return(chrommap)
}
