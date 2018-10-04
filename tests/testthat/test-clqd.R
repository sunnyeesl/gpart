context("test-clqd")

test_that("newSplitCliques", {
  test <- newSplitCliques(cliques.bp = list(c(1, 500, 510)), gapdist = 20)
  test2 <- newSplitCliques(cliques.bp = list(c(1, 500, 510)), gapdist = 100)
  expect_equal(test, test2)
})
