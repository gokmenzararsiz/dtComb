library(APtools)
library(usethis)

data("exampleData1")
Data <- exampleData1[-c(83:138),]
markers <- Data[, -1]
status <- factor(Data$group, levels = c("not_needed", "needed"))

data(mayo)
Data2 <- mayo[-c(42:119),]
markers2 <- Data2[, 3:4]
status2 <- factor(Data2[, 2], levels = c(1, 0))

Data3 <-
  read.csv(
    "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data",
    header = FALSE
  )
Data3 <- Data3[-c(121:262),]
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
      cutoff.method = "youden"
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
      cutoff.method = "youden"
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
      cutoff.method = "youden"
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
      cutoff.method = "youden"
    ),
    "direction should be one of “auto”, “<”, “>”"
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
    "cutoff.method should be one of “youden”, “roc01”"
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
      cutoff.method = "youden"
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
      cutoff.method = "youden"
    ),
    "status does not include event"
  )
})

