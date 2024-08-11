data("exampleData1")
Data <- exampleData1[-c(83:138), ]
markers <- Data[, -1]
status <- factor(Data$group, levels = c("not_needed", "needed"))

load("result_data/mayo.rda")
Data2 <- mayo[-c(42:119), ]
markers2 <- Data2[, 3:4]
status2 <- factor(Data2[, 2], levels = c(1, 0))

Data3 <-
  read.csv(
    "result_data/wdbc.data.txt",
    header = FALSE
  )
Data3 <- Data3[-c(121:262), ]
markers3 <- Data3[, 4:5]
status3 <- factor(Data3[, 2], levels = c("B", "M"))

###############################################################################

status4 <- factor(Data3[, 2], levels = c("B", "M", "C"))
status4[[9]] <- "C"

test_that("mlComb functions ...", {
  expect_error(
    mlComb(
      markers = Data3[, 4:5],
      status = status4,
      event = "M",
      method = "knn",
      direction = "<",
      cutoff.method = "Youden"
    ),
    "the number of status levels should be 2"
  )

  expect_error(
    mlComb(
      markers = Data3[, 4:6],
      status = status3,
      event = "M",
      method = "svmLinear",
      direction = "<",
      cutoff.method = "roc01"
    ),
    "the number of markers should be 2"
  )
})

test_that("mlComb functions ...", {
  expect_error(
    mlComb(
      markers = markers,
      status = status,
      event = "needed",
      direction = "<",
      cutoff.method = "Youden"
    ),
    "The response given method is not available for mlComb function. See availableMethods function for the list of methods available."
  )

  expect_error(
    mlComb(
      markers = markers3,
      status = status3,
      event = "M",
      method = "asaddsa",
      direction = "auto",
      cutoff.method = "Youden"
    ),
    "The response given method is not available for mlComb function. See availableMethods function for the list of methods available."
  )

  expect_error(
    mlComb(
      markers = markers2,
      status = status2,
      event = "1",
      method = "rf",
      direction = "asdada",
      cutoff.method = "Youden"
    ),
    "direction should be one of 'auto', '<', '>'"
  )

  expect_error(
    mlComb(
      markers = markers2,
      status = status2,
      event = "1",
      method = "knn",
      direction = "auto",
      cutoff.method = "sadda"
    ),
    "The entered cutoff.method is invalid"
  )
})

###############################################################################

markers3[44, 1:2] <- "assay"

test_that("mlComb functions ...", {
  expect_error(
    mlComb(
      markers = markers3,
      status = status3,
      event = "M",
      method = "rf",
      direction = "<",
      cutoff.method = "Youden"
    ),
    "at least one variable is not numeric"
  )
  expect_error(
    mlComb(
      markers = markers,
      status = status,
      event = "C",
      method = "glm",
      direction = "<",
      cutoff.method = "Youden"
    ),
    "status does not include event"
  )
})
