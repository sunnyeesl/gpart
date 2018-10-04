context("test-bigld")

test_that("intervalCliqueList", {
 test <- intervalCliqueList(list(c(1,3), c(2,3,4)))
 exptest <- rbind(c(1, 3), c(1, 4), c(2, 4))
 expect_equivalent(test[,1:2], exptest)
})

test_that("findMaximumIndept", {
  test <- findMaximumIndept(rbind(c(1,3), c(4, 6)), c(1, 1))
  test2 <- findMaximumIndept(rbind(c(1,5), c(4, 6)), c(2, 1))
  expect_equivalent(test[[1]], c(1,2))
  expect_equivalent(test2[[1]], c(1))
  test3 <- findMaximumIndept(rbind(c(1,5), c(4, 8), c(6, 10)), c(1, 3, 1))
  expect_equivalent(test3[[1]], c(2))
})
