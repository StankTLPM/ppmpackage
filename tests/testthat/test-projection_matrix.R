test_that("numeric pre-breeding projection_matrix returns expected matrix", {
  t <- list(c(2, 1, 0.5), c(3, 2, 0.6))  # t_1_2 = 0.5, t_2_3 = 0.6
  r <- list(c(1, 3, 1.2))               # r_3_1 = 1.2

  res <- projection_matrix(3, t, r, mode = "num", timing = "pre")

  expect_type(res, "list")
  expect_named(res, c("matrix", "growth_rate", "stable_distribution"))
  expect_equal(dim(res$matrix), c(3, 3))
  expect_true(all(res$matrix >= 0))
  expect_equal(sum(res$stable_distribution), 1, tolerance = 1e-6)
})

test_that("numeric post-breeding projection_matrix works", {
  t <- list(c(2, 1, 0.5), c(3, 2, 0.6))
  r <- list(c(1, 3, 1.2))

  res <- projection_matrix(3, t, r, mode = "num", timing = "post")

  expect_type(res, "list")
  expect_equal(dim(res$matrix), c(3, 3))
  expect_true(res$growth_rate > 0)
})

test_that("literal pre-breeding returns symbolic matrix", {
  t <- list(c(2, 1, "t12"), c(3, 2, "t23"))
  r <- list(c(1, 3, "r31"))

  res <- projection_matrix(3, t, r, mode = "lit", timing = "pre")

  expect_true(is.matrix(res))
  expect_equal(dim(res), c(3, 3))
  expect_true(all(sapply(res, is.character)))
  expect_true(any(grepl("t12", res)))
  expect_true(any(grepl("r31", res)))
})

test_that("literal post-breeding returns symbolic matrix", {
  t <- list(c(2, 1, "t12"), c(3, 2, "t23"))
  r <- list(c(1, 3, "r31"))

  res <- projection_matrix(3, t, r, mode = "lit", timing = "post")

  expect_equal(dim(res), c(3, 3))
  expect_true(any(grepl("t23", res)))
  expect_true(any(grepl("r31", res)))
})
