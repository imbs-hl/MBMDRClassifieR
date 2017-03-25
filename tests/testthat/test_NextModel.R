
context("Function NextModel")

test_that("NextModel is integer vector", {

  k <- 10L
  order <- 2L
  model <- c(1L, 2L)

  expect_is(MBMDRClassifieR:::NextModel(model = model,
                                        k = k,
                                        order = order,
                                        j = order),
            class = "integer")
})

test_that("NextModel increases last position", {

  k <- 10L
  order <- 2L
  model <- c(1L, 2L)

  expect_equal(MBMDRClassifieR:::NextModel(model = model,
                                           k = k,
                                           order = order,
                                           j = order),
               c(1L, 3L))
})

test_that("NextModel switches position to increase", {

  k <- 10L
  order <- 2L
  model <- c(1L, 10L)

  expect_equal(MBMDRClassifieR:::NextModel(model = model,
                                           k = k,
                                           order = order,
                                           j = order),
               c(2L, 3L))

})

test_that("NextModel works for arbitrary orders", {

  k <- sample(10L:100L, 1)

  order <- 1L
  model <- 1L:order

  expect_is(MBMDRClassifieR:::NextModel(model = model,
                                        k = k,
                                        order = order,
                                        j = order),
            class = "integer")

  order <- 2L
  model <- 1L:order

  expect_is(MBMDRClassifieR:::NextModel(model = model,
                                        k = k,
                                        order = order,
                                        j = order),
            class = "integer")

  order <- 3L
  model <- 1L:order

  expect_is(MBMDRClassifieR:::NextModel(model = model,
                                        k = k,
                                        order = order,
                                        j = order),
            class = "integer")

  order <- 4L
  model <- 1L:order

  expect_is(MBMDRClassifieR:::NextModel(model = model,
                                        k = k,
                                        order = order,
                                        j = order),
            class = "integer")

  order <- 5L
  model <- 1L:order

  expect_is(MBMDRClassifieR:::NextModel(model = model,
                                        k = k,
                                        order = order,
                                        j = order),
            class = "integer")

})

test_that("NextModel stops", {

  k <- 10L
  order <- 2L
  model <- c(9L, 10L)

  expect_false(MBMDRClassifieR:::NextModel(model = model,
                                           k = k,
                                           order = order,
                                           j = order))

})
