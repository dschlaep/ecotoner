context("Mathematical utilities")


x1 <- c(rep(NA, 3), 1:4)
x2 <- c(rep(1/3, 3), 1:4, rep(NA, 3))


test_that("Mathematical utilities", {
	expect_equal(majority(x1), NA_real_)
	expect_equal(majority(x1, return_single = FALSE, na.rm = TRUE), 1:4)
	expect_equal(majority(x1, use_mode = FALSE, return_single = FALSE, na.rm = TRUE), 1:4)
	expect_equal(majority(x1, use_mode = FALSE, return_single = FALSE, na.rm = FALSE), NA_real_)
	expect_equal(majority(x2, return_single = FALSE), c(NA, 1/3))
	expect_true(majority(x1, na.rm = TRUE) %in% na.exclude(x1))
})

test_that("Mathematical utilities with a random component", {
	expect_true(majority(x1, na.rm = TRUE) %in% na.exclude(x1))
	expect_true(majority(x1, use_mode = FALSE, na.rm = TRUE) %in% na.exclude(x1))
})
