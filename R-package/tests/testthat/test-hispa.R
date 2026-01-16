context("HiSpaR basic functionality tests")

test_that("convolute_contacts works correctly", {
  # Create a simple test matrix
  n <- 10
  mat <- matrix(rpois(n*n, 5), n, n)
  mat <- (mat + t(mat)) / 2  # Make symmetric
  
  # Convolute
  result <- convolute_contacts(mat, half_k = 2)
  
  # Tests
  expect_equal(dim(result), c(n, n))
  expect_true(isSymmetric(result))
  expect_true(all(result >= 0))
})

test_that("hispa_analyze handles invalid inputs", {
  # Test with non-matrix input
  expect_error(hispa_analyze(c(1, 2, 3), "output"))
  
  # Test with non-symmetric matrix (should symmetrize with warning)
  mat <- matrix(1:9, 3, 3)
  expect_warning(
    hispa_analyze(mat, tempdir(), mcmc_iterations = 10)
  )
})

test_that("hispa_analyze returns correct structure", {
  skip_if_not_installed("HiSpaR")
  
  # Create small test matrix
  n <- 20
  mat <- matrix(rpois(n*n, 10), n, n)
  mat <- (mat + t(mat)) / 2
  
  result <- hispa_analyze(
    contact_matrix = mat,
    output_dir = tempdir(),
    mcmc_iterations = 100,
    mcmc_burn_in = 10,
    verbose = FALSE
  )
  
  # Check structure
  expect_type(result, "list")
  expect_true("position_matrix" %in% names(result))
  expect_true("beta0" %in% names(result))
  expect_true("beta1" %in% names(result))
  expect_true("log_likelihood" %in% names(result))
  
  # Check dimensions
  expect_equal(nrow(result$position_matrix), n)
  expect_equal(ncol(result$position_matrix), 3)
  
  # Check class
  expect_s3_class(result, "hispa_result")
})
