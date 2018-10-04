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


test_that("BigLD", {
  testmat <- matrix(1, 10, 10)
  testmat[1:3,]<-0
  testmat[7:10, ]<-2
  diag(testmat)<-0
  testSNP <- data.frame(chrN = 1,
                        rsID = paste("SNP", seq_len(10), sep = ""),
                        bp = seq(1, 100, 10), stringsAsFactors = FALSE)
  test1 <- BigLD(geno = testmat, SNPinfo = testSNP)
  test1exp <- data.frame(chr = 1, start.index = 1, end.index = 7,
                         start.rsID = "SNP1", end.rsID = "SNP7",
                         start.bp = 1, end.bp = 61)
  expect_equal(test1, test1exp)
})
