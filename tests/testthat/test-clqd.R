context("test-clqd")

test_that("newSplitCliques", {
  test <- newSplitCliques(cliques.bp = list(c(1, 500, 510)), gapdist = 20)
  test2 <- newSplitCliques(cliques.bp = list(c(1, 500, 510)), gapdist = 100)
  expect_equal(test, test2)
})

test_that("genoDp", {
  testmat <- matrix(1, 10, 10)
  testmat[1:3,]<-0
  testmat[7:10, ]<-2
  diag(testmat)<-0
  test1 <- genoDp2(testmat[,1,drop=FALSE], testmat[,2,drop=FALSE], strLD = FALSE)[1,1]
  expect_equivalent(test1, 1)
  test1 <- genoDp2(testmat[,1,drop=FALSE], testmat[,6,drop=FALSE], strLD = FALSE)[1,1]
  expect_equivalent(test1, 1)
})


test_that("CLQD", {
  testmat <- matrix(1, 10, 10)
  testmat[1:3,]<-0
  testmat[7:10, ]<-2
  diag(testmat)<-0
  testSNP <- data.frame(chrN = 1,
                        rsID = paste("SNP", seq_len(10), sep = ""),
                        bp = seq(1, 100, 10), stringsAsFactors = FALSE)
  test1 <- CLQD(geno = testmat, SNPinfo = testSNP)
  test1exp <- c(rep(1, 7), rep(NA, 3))
  expect_equal(test1, test1exp)

  test2 <- CLQD(geno = testmat, SNPinfo = testSNP, clstgap = 3)
  expect_true(all(is.na(test2)))
})
