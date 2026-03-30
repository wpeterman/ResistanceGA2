test_that("To.From.ID returns data frame with correct columns and rows", {
  id <- To.From.ID(5)
  expect_s3_class(id, "data.frame")
  expect_named(id, c("pop1", "pop2"))
  expect_equal(nrow(id), choose(5, 2))  # 10
})

test_that("To.From.ID pop1 and pop2 are factors", {
  id <- To.From.ID(4)
  expect_s3_class(id$pop1, "factor")
  expect_s3_class(id$pop2, "factor")
})

test_that("To.From.ID all pairs are distinct (pop1 != pop2)", {
  id <- To.From.ID(6)
  expect_true(all(as.character(id$pop1) != as.character(id$pop2)))
})

test_that("To.From.ID correct row count for several n values", {
  for (n in 3:8) {
    id <- To.From.ID(n)
    expect_equal(nrow(id), choose(n, 2),
                 info = paste("n =", n))
  }
})

test_that("To.From.ID with pop_n returns individual-level data frame", {
  id <- To.From.ID(sampled_pops = 3, pop_n = c(1L, 2L, 1L))
  # 4 individuals total; 5 unique inter-population pairs (within-pop self-pairs excluded)
  expect_equal(nrow(id), 5L)
  expect_true("pop1.ind" %in% names(id))
  expect_true("pop1.pop" %in% names(id))
})

test_that("To.From.ID stops when spLoc supplied without nb", {
  spLoc <- terra::vect(cbind(c(1, 2, 9), c(1, 2, 9)), type = "points")
  expect_error(To.From.ID(3, spLoc = spLoc), "nb")
})

test_that("To.From.ID accepts SpatVector for spLoc", {
  spLoc <- terra::vect(cbind(c(1, 2, 9, 10), c(1, 2, 9, 10)), type = "points")
  id <- To.From.ID(4, spLoc = spLoc, nb = 3)
  expect_named(id, c("pop1", "pop2", "corr_", "cor.grp"))
  expect_equal(nrow(id), 6L)
})

test_that("To.From.ID errors on coordinate matrix for spLoc", {
  spLoc <- cbind(c(1, 2, 9, 10), c(1, 2, 9, 10))
  expect_error(To.From.ID(4, spLoc = spLoc, nb = 3), "SpatVector")
})
