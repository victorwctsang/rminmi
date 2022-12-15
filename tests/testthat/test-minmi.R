test_that("Delta model point estimate works", {
  theta.true <- 10000
  K <- 18000
  n <- 30
  B <- 100
  q <- 0.5

  set.seed(0)
  ages <- runif(n, theta.true, K)
  u <- runif(B, 0 , 1)

  EXPECTED_DELTA_MODEL_POINT_ESTIMATE <- K - q^(-1/n)*(K-min(ages))

  expect_equal(estimate_quantile.minmi(K=K, W=ages, u=u, eps.mean=0, eps.sigma=0, q=q), EXPECTED_DELTA_MODEL_POINT_ESTIMATE)
})
