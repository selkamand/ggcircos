# Create normalisation map ------------------------------------------------

test_that("create_normalisation_map maps a single chromosome to the full target range", {
  stdref <- data.frame(
    name = "chr1",
    start = 1,
    end = 100,
    refindex = 1
  )

  result <- create_normalisation_map(stdref, inter_chrom_spacing = 0.1)

  expect_named(result, "chr1")
  expect_equal(result$chr1$original, c(start = 1, end = 100))
  expect_equal(result$chr1$new, c(start = 1, end = 100))
})


test_that("create_normalisation_map preserves chromosome names and original ranges", {
  stdref <- data.frame(
    name = c("chr1", "chr2", "chr3"),
    start = c(1, 1, 1),
    end = c(100, 200, 300),
    refindex = c(1, 2, 3)
  )

  result <- create_normalisation_map(stdref, inter_chrom_spacing = 0)

  expect_named(result, c("chr1", "chr2", "chr3"))

  expect_equal(result$chr1$original, c(start = 1, end = 100))
  expect_equal(result$chr2$original, c(start = 1, end = 200))
  expect_equal(result$chr3$original, c(start = 1, end = 300))
})

test_that("create_normalisation_map maps chromosomes monotonically into target range", {
  stdref <- data.frame(
    name = c("chr1", "chr2", "chr3"),
    start = c(1, 1, 1),
    end = c(100, 200, 300),
    refindex = c(1, 2, 3)
  )

  result <- create_normalisation_map(stdref, inter_chrom_spacing = 0)

  starts <- vapply(result, function(x) x$new[["start"]], numeric(1))
  ends <- vapply(result, function(x) x$new[["end"]], numeric(1))

  expect_equal(starts[["chr1"]], 1)
  expect_equal(ends[["chr3"]], 100)

  expect_true(all(starts < ends))
  expect_true(all(diff(starts) > 0))
  expect_true(all(diff(ends) > 0))
})

test_that("create_normalisation_map creates gaps when inter_chrom_spacing is positive", {
  stdref <- data.frame(
    name = c("chr1", "chr2", "chr3"),
    start = c(1, 1, 1),
    end = c(100, 100, 100),
    refindex = c(1, 2, 3)
  )

  no_gap <- create_normalisation_map(stdref, inter_chrom_spacing = 0)
  with_gap <- create_normalisation_map(stdref, inter_chrom_spacing = 0.1)

  no_gap_chr1_end <- no_gap$chr1$new[["end"]]
  no_gap_chr2_start <- no_gap$chr2$new[["start"]]

  with_gap_chr1_end <- with_gap$chr1$new[["end"]]
  with_gap_chr2_start <- with_gap$chr2$new[["start"]]

  expect_equal(no_gap_chr1_end, no_gap_chr2_start)
  expect_gt(with_gap_chr2_start, with_gap_chr1_end)
})

test_that("create_normalisation_map applies no gap before the first chromosome", {
  stdref <- data.frame(
    name = c("chr1", "chr2"),
    start = c(1, 1),
    end = c(100, 100),
    refindex = c(1, 2)
  )

  result <- create_normalisation_map(stdref, inter_chrom_spacing = 0.1)

  expect_equal(result$chr1$new[["start"]], 1)
})


# Project Points ----------------------------------------------------------

test_that("project_points projects a single point", {
  map <- list(
    chr1 = list(
      original = c(start = 1, end = 100),
      new = c(start = 1, end = 50)
    )
  )

  result <- project_points("chr1", 50.5, map)

  expect_type(result, "double")
  expect_length(result, 1)
  # browser()
  expect_equal(result, 25.5)
})


test_that("project_points projects multiple points from the same chromosome", {
  map <- list(
    chr1 = list(
      original = c(start = 1, end = 100),
      new = c(start = 1, end = 50)
    )
  )

  result <- project_points(
    chr = c("chr1", "chr1", "chr1"),
    pos = c(1, 50.5, 100),
    map = map
  )

  expect_type(result, "double")
  expect_length(result, 3)
  expect_equal(result, c(1, 25.5, 50))
})

test_that("project_points projects points across multiple chromosomes", {
  map <- list(
    chr1 = list(
      original = c(start = 1, end = 100),
      new = c(start = 1, end = 40)
    ),
    chr2 = list(
      original = c(start = 1, end = 200),
      new = c(start = 50, end = 100)
    )
  )

  result <- project_points(
    chr = c("chr1", "chr2"),
    pos = c(100, 200),
    map = map
  )

  expect_type(result, "double")
  expect_length(result, 2)
  expect_equal(result, c(40, 100))
})

test_that("project_points preserves input order", {
  map <- list(
    chr1 = list(
      original = c(start = 1, end = 100),
      new = c(start = 1, end = 40)
    ),
    chr2 = list(
      original = c(start = 1, end = 200),
      new = c(start = 50, end = 100)
    )
  )

  result <- project_points(
    chr = c("chr2", "chr1", "chr2"),
    pos = c(1, 100, 200),
    map = map
  )

  expect_equal(result, c(50, 40, 100))
})

test_that("project_points returns a numeric vector, not a list", {
  map <- list(
    chr1 = list(
      original = c(start = 1, end = 10),
      new = c(start = 100, end = 200)
    )
  )

  result <- project_points(
    chr = c("chr1", "chr1"),
    pos = c(1, 10),
    map = map
  )

  expect_true(is.numeric(result))
  expect_false(is.list(result))
  # browser()
  expect_equal(result, c(100, 200))
})

test_that("project_points errors when chr and pos have different lengths", {
  map <- list(
    chr1 = list(
      original = c(start = 1, end = 100),
      new = c(start = 1, end = 50)
    )
  )

  expect_error(
    project_points(
      chr = c("chr1", "chr1"),
      pos = 1,
      map = map
    ),
    "`chr` and `pos` must be the same length",
    fixed = TRUE
  )
})

test_that("project_points errors when pos is not numeric", {
  map <- list(
    chr1 = list(
      original = c(start = 1, end = 100),
      new = c(start = 1, end = 50)
    )
  )

  expect_error(
    project_points("chr1", "50", map),
    "`pos` must be numeric",
    fixed = TRUE
  )
})

test_that("project_points returns NA for chromosomes missing from the map", {
  map <- list(
    chr1 = list(
      original = c(start = 1, end = 100),
      new = c(start = 1, end = 50)
    )
  )

  expect_true(is.na(project_points("chr2", 10, map)))
})


# Project Intervals -------------------------------------------------------
test_that("project_intervals projects a single interval", {
  map <- list(
    chr1 = list(
      original = c(start = 1, end = 100),
      new = c(start = 1, end = 50)
    )
  )

  result <- project_intervals(
    chr = "chr1",
    start = 1,
    end = 100,
    map = map
  )

  expect_type(result, "list")
  expect_length(result, 1)

  expect_true(is.numeric(result[[1]]))
  expect_length(result[[1]], 2)
  expect_equal(result[[1]], c(1, 50))
})

test_that("project_intervals projects multiple intervals from the same chromosome", {
  map <- list(
    chr1 = list(
      original = c(start = 1, end = 100),
      new = c(start = 1, end = 50)
    )
  )

  result <- project_intervals(
    chr = c("chr1", "chr1", "chr1"),
    start = c(1, 25.75, 50.5),
    end = c(25.75, 50.5, 100),
    map = map
  )

  expect_type(result, "list")
  expect_length(result, 3)

  expect_equal(result[[1]], c(1, 13.25))
  expect_equal(result[[2]], c(13.25, 25.5))
  expect_equal(result[[3]], c(25.5, 50))
})

test_that("project_intervals projects intervals across multiple chromosomes", {
  map <- list(
    chr1 = list(
      original = c(start = 1, end = 100),
      new = c(start = 1, end = 40)
    ),
    chr2 = list(
      original = c(start = 1, end = 200),
      new = c(start = 50, end = 100)
    )
  )

  result <- project_intervals(
    chr = c("chr1", "chr2"),
    start = c(1, 1),
    end = c(100, 200),
    map = map
  )

  expect_type(result, "list")
  expect_length(result, 2)

  expect_equal(result[[1]], c(1, 40))
  expect_equal(result[[2]], c(50, 100))
})

test_that("project_intervals preserves input order", {
  map <- list(
    chr1 = list(
      original = c(start = 1, end = 100),
      new = c(start = 1, end = 40)
    ),
    chr2 = list(
      original = c(start = 1, end = 200),
      new = c(start = 50, end = 100)
    )
  )

  result <- project_intervals(
    chr = c("chr2", "chr1", "chr2"),
    start = c(1, 1, 100.5),
    end = c(200, 100, 200),
    map = map
  )

  expect_equal(result[[1]], c(50, 100))
  expect_equal(result[[2]], c(1, 40))
  expect_equal(result[[3]], c(75, 100))
})

test_that("project_intervals returns one two-element numeric vector per interval", {
  map <- list(
    chr1 = list(
      original = c(start = 1, end = 10),
      new = c(start = 100, end = 200)
    )
  )

  result <- project_intervals(
    chr = c("chr1", "chr1"),
    start = c(1, 5.5),
    end = c(5.5, 10),
    map = map
  )

  expect_type(result, "list")
  expect_length(result, 2)

  expect_true(all(vapply(result, is.numeric, logical(1))))
  expect_true(all(vapply(result, length, integer(1)) == 2L))

  expect_equal(result[[1]], c(100, 150))
  expect_equal(result[[2]], c(150, 200))
})

test_that("project_intervals errors when chr and start have different lengths", {
  map <- list(
    chr1 = list(
      original = c(start = 1, end = 100),
      new = c(start = 1, end = 50)
    )
  )

  expect_error(
    project_intervals(
      chr = c("chr1", "chr1"),
      start = 1,
      end = c(10, 20),
      map = map
    ),
    "`chr` and `start` must be the same length",
    fixed = TRUE
  )
})

test_that("project_intervals errors when chr and end have different lengths", {
  map <- list(
    chr1 = list(
      original = c(start = 1, end = 100),
      new = c(start = 1, end = 50)
    )
  )

  expect_error(
    project_intervals(
      chr = c("chr1", "chr1"),
      start = c(1, 10),
      end = 20,
      map = map
    ),
    "`chr` and `end` must be the same length",
    fixed = TRUE
  )
})

test_that("project_intervals errors when start is not numeric", {
  map <- list(
    chr1 = list(
      original = c(start = 1, end = 100),
      new = c(start = 1, end = 50)
    )
  )

  expect_error(
    project_intervals(
      chr = "chr1",
      start = "1",
      end = 100,
      map = map
    ),
    "`pos` must be numeric",
    fixed = TRUE
  )
})

test_that("project_intervals errors when end is not numeric", {
  map <- list(
    chr1 = list(
      original = c(start = 1, end = 100),
      new = c(start = 1, end = 50)
    )
  )

  expect_error(
    project_intervals(
      chr = "chr1",
      start = 1,
      end = "100",
      map = map
    ),
    "`pos` must be numeric",
    fixed = TRUE
  )
})


test_that("project_intervals can return a data.frame", {
  map <- list(
    chr1 = list(
      original = c(start = 1, end = 100),
      new = c(start = 1, end = 50)
    )
  )

  result <- project_intervals(
    chr = c("chr1", "chr1"),
    start = c(1, 50.5),
    end = c(50.5, 100),
    map = map,
    return_type = "data.frame"
  )

  expect_s3_class(result, "data.frame")
  expect_named(result, c("chr", "start", "end"))
  expect_equal(nrow(result), 2L)
  expect_equal(ncol(result), 3L)
})

test_that("project_intervals data.frame output contains projected starts and ends", {
  map <- list(
    chr1 = list(
      original = c(start = 1, end = 100),
      new = c(start = 1, end = 50)
    )
  )

  result <- project_intervals(
    chr = c("chr1", "chr1"),
    start = c(1, 50.5),
    end = c(50.5, 100),
    map = map,
    return_type = "data.frame"
  )

  expect_equal(result$chr, c("chr1", "chr1"))
  expect_equal(result$start, c(1, 25.5))
  expect_equal(result$end, c(25.5, 50))
})
test_that("project_intervals data.frame output preserves input order across chromosomes", {
  map <- list(
    chr1 = list(
      original = c(start = 1, end = 100),
      new = c(start = 1, end = 40)
    ),
    chr2 = list(
      original = c(start = 1, end = 200),
      new = c(start = 50, end = 100)
    )
  )

  result <- project_intervals(
    chr = c("chr2", "chr1", "chr2"),
    start = c(1, 1, 100.5),
    end = c(200, 100, 200),
    map = map,
    return_type = "data.frame"
  )

  expected <- data.frame(
    chr = c("chr2", "chr1", "chr2"),
    start = c(50, 1, 75),
    end = c(100, 40, 100)
  )

  expect_equal(result, expected)
})

test_that("project_intervals defaults to list output", {
  map <- list(
    chr1 = list(
      original = c(start = 1, end = 100),
      new = c(start = 1, end = 50)
    )
  )

  result <- project_intervals(
    chr = "chr1",
    start = 1,
    end = 100,
    map = map
  )

  expect_type(result, "list")
  expect_equal(result[[1]], c(1, 50))
})

test_that("project_intervals errors for invalid return_type", {
  map <- list(
    chr1 = list(
      original = c(start = 1, end = 100),
      new = c(start = 1, end = 50)
    )
  )

  expect_error(
    project_intervals(
      chr = "chr1",
      start = 1,
      end = 100,
      map = map,
      return_type = "matrix"
    )
  )
})
