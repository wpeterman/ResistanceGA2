test_that("lower() extracts correct number of values from 4x4 matrix", {
  m <- matrix(1:16, nrow = 4)
  expect_length(lower(m), 6L)  # choose(4, 2) = 6
})

test_that("lower() values match base lower.tri()", {
  m <- matrix(1:9, nrow = 3)
  expect_equal(lower(m), m[lower.tri(m)])
})

test_that("lower() is consistent across matrix sizes", {
  for (n in 3:7) {
    m <- matrix(seq_len(n^2), nrow = n)
    expect_length(lower(m), choose(n, 2))
  }
})

test_that("lower() warns on non-square matrix", {
  m <- matrix(1:6, nrow = 2, ncol = 3)
  expect_warning(lower(m))
})

test_that("lower() warns on vector input", {
  expect_warning(lower(1:5))
})
