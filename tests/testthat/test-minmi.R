test_that("Delta model point estimate works", {
  q <- 0.5
  theta.true <- 10000
  K <- 18000
  n <- 30
  B <- 100

  set.seed(0)
  ages <- runif(n, theta.true, K)
  u <- runif(B, 0 , 1)

  EXPECTED_DELTA_MODEL_POINT_ESTIMATE <- K - q^(-1/n)*(K-min(ages))

  expect_equal(estimate_quantile.minmi(K=K, W=ages, u=u, eps.sigma=0, q=q), EXPECTED_DELTA_MODEL_POINT_ESTIMATE)
})

test_that("Epsilon model (with same standard deviations) point estimate works", {
  q <- 0.5
  K <- 18000
  n <- 30
  sd <- rep(200, n)
  theta.true <- 10000
  K <- 18000
  n <- 30
  sd <- rep(200, n)
  theta.true <- 10000

  set.seed(0); W <- runif(n, theta.true, K) + sd
  set.seed(1); u <- matrix(runif(n*100, 0, 1), ncol=100)

  EXPECTED_EPSILON_SAME_SD_POINT_ESTIMATE <- 10151.6242

  expect_equal(estimate_quantile.minmi(K=K, W=W, u=u, eps.sigma=rep(200, 30), q=q), EXPECTED_EPSILON_SAME_SD_POINT_ESTIMATE)
})
