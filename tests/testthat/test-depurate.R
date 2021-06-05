## tests for depurate
set.seed(3)
crds <- data.frame(
  x = runif(20, min = 232400, max = 232500),
  y = runif(20, min = 6366600, max = 6366700)
)
crds <- crds[order(crds$x, crds$y),]

r <- sample(c(-1, 0, rnorm(17, 2.5), rnorm(1, 4)))
mapa <- data.frame(crds, r)

mapa <- sf::st_as_sf(mapa, coords = c('x', 'y'), crs = 32721)


test_that("edges removes fine", {
  dep <- depurate(mapa,
           'r',
           'edges',
           buffer = -10)

  expect_equal(sum(is.na(dep$condition)), 4)
  expect_equal(sum(dep$condition == 'border', na.rm = T), 16)
  expect_equal(sum(is.na(dep$condition)) +
                 sum(dep$condition == 'border', na.rm = T), nrow(mapa))
})

# expect_equal(sum(is.na(dep$condition)), 7)
# expect_equal(sum(dep$condition == 'border', na.rm = T), 16)
# expect_equal(sum(is.na(dep$condition)) +
#                sum(dep$condition == 'border', na.rm = T), nrow(mapa))


