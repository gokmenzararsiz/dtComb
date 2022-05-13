library(APtools)
library(usethis)



data("exampleData1")
Data <- exampleData1[-c(83:138), ]
markers <- Data[,-1]
status <- factor(Data$group, levels = c("not_needed", "needed"))


data(mayo)
Data2 <- mayo[-c(42:119), ]
markers2 <- Data2[, 3:4]
status2 <- factor(Data2[, 2], levels = c(1, 0))


Data3 <-
  read.csv(
    "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data",
    header = FALSE
  )
Data3 <- Data3[-c(121:262), ]
markers3 <- Data3[, 4:5]
status3 <- factor(Data3[, 2], levels = c("B", "M"))

load("result_data/test_mathComb.rda")

#comb.score, AUC, SEN, SPE ve Cutoff kontrolü

for (distance in c("lorentzian",
                   "avg",
                   "taneja",
                   "kumar-johnson")) {
  set.seed(14042022)
  res <- mathComb(
    markers = markers,
    status = status,
    event = "needed",
    method = "distance",
    distance = distance,
    direction = "<",
    cutoff.method = "youden"
  )
  
  test_that("mathComb functions ...", {
    expect_length(res, 11)
    expect_equal(as.numeric(res$AUC_table$AUC[[3]]),  r$AUC[r$Distance == distance][1], tolerance =
                   0.01)
    expect_equal(as.numeric(res$DiagStatCombined$detail$sp[[1]]),
                 r$SPE[r$Distance == distance][1],
                 tolerance = 0.01)
    expect_equal(as.numeric(res$DiagStatCombined$detail$se[[1]]),
                 r$SENS[r$Distance == distance][1],
                 tolerance = 0.01)
    expect_equal(as.numeric(res$ThresholdCombined), r$Cutoff[r$Distance == distance][1], tolerance =
                   0.01)
  })
}

#mayo Datası ile
for (method in c("add",
                 "multiply",
                 "divide",
                 "subtract",
                 "baseinexp",
                 "expinbase"
                 )) {
  set.seed(14042022)
  res <- mathComb(
    markers = markers2,
    status = status2,
    event = "1",
    method = method,
    direction = "<",
    cutoff.method = "youden"
  )
  
  test_that("mathComb functions ...", {
    expect_length(res, 11)
    expect_equal(as.numeric(res$AUC_table$AUC[[3]]),  r$AUC[r$Method == method][1], tolerance =
                   0.01)
    expect_equal(as.numeric(res$DiagStatCombined$detail$sp[[1]]),
                 r$SPE[r$Method == method][1],
                 tolerance = 0.01)
    expect_equal(as.numeric(res$DiagStatCombined$detail$se[[1]]),
                 r$SENS[r$Method == method][1],
                 tolerance = 0.01)
    expect_equal(as.numeric(res$ThresholdCombined), r$Cutoff[r$Method == method][1], tolerance =
                   0.01)
  })
}

#WDBC Datası ile
for (distance in c("kulczynski_d",
                 "euclidean",
                 "manhattan",
                 "chebyshev")) {
  set.seed(14042022)
  res <- mathComb(
    markers = markers3,
    status = status3,
    event = "M",
    method = "distance",
    distance = distance,
    direction = "<",
    cutoff.method = "youden"
  )
  
  test_that("mathComb functions ...", {
    expect_length(res, 11)
    expect_equal(as.numeric(res$AUC_table$AUC[[3]]),  r$AUC[r$Distance == distance][1], tolerance =
                   0.01)
    expect_equal(as.numeric(res$DiagStatCombined$detail$sp[[1]]),
                 r$SPE[r$Distance == distance][1],
                 tolerance = 0.01)
    expect_equal(as.numeric(res$DiagStatCombined$detail$se[[1]]),
                 r$SENS[r$Distance == distance][1],
                 tolerance = 0.01)
    expect_equal(as.numeric(res$ThresholdCombined), r$Cutoff[r$Distance == distance][1], tolerance =
                   0.01)
  })
}


#Hata kontrolü
status4 <- factor(Data3[, 2], levels = c("B", "M", "C"))
status4[[9]] <- "C"

test_that("mathComb functions ...", {
  expect_error(
    mathComb(
      markers = Data3[, 4:5],
      status = status4,
      event = "M",
      method = "multiply",
      direction = direction,
      standardize = "zScore",
      cutoff.method = cutoff.method
    ),
    "the number of status levels should be 2"
  )
  
  expect_error(
    mathComb(
      markers = Data3[, 4:6],
      status = status3,
      event = "M",
      method = "multiply",
      direction = direction,
      standardize = "zScore",
      cutoff.method = cutoff.method
    ),
    "the number of markers should be 2"
  )
})



test_that("mathComb functions ...", {
  expect_error(
    mathComb(
      markers = markers3,
      status = status3,
      event = "M",
      direction = "<",
      standardize = "none",
      cutoff.method = "youden"
    ),
    "method should be one of “add”, “multiply”, “divide”, “subtract” ,“distance”, “baseinexp”, “expinbase”"
  )
  
  expect_error(
    mathComb(
      markers = markers3,
      status = status3,
      event = "M",
      method = "adsadad",
      direction = "auto",
      standardize = "none",
      cutoff.method = "youden"
    ),
    "method should be one of “add”, “multiply”, “divide”, “subtract” ,“distance”, “baseinexp”, “expinbase”"
  )
  
  expect_error(
    mathComb(
      markers = markers3,
      status = status3,
      event = "M",
      method = "add",
      direction = "auto",
      standardize = "asdada",
      cutoff.method = "youden"
    ),
    "standardize should be one of “range”, “zScore”, “tScore”, “mean”, “deviance”"
  )
  
  expect_error(
    mathComb(
      markers = markers2,
      status = status2,
      event = "1",
      method = "baseinexp",
      direction = "asdada",
      standardize = "none",
      cutoff.method = "youden"
    ),
    "direction should be one of “auto”, “<”, “>”"
  )
  
  expect_error(
    mathComb(
      markers = markers2,
      status = status2,
      event = "1",
      method = "expinbase",
      direction = "auto",
      standardize = "none",
      cutoff.method = "sadda"
    ),
    "cutoff.method should be one of “youden”, “roc01”"
  )

  expect_warning(
    mathComb(
      markers = markers,
      status = status,
      event = "needed",
      method = "distance",
      distance = "kumar-johnson",
      direction = "<",
      standardize = "zScore",
      cutoff.method = "youden"
    ),
    "Infinity or NaNs values generated in markers, standardization changed to 'none'."
  )
  
  expect_warning(
    mathComb(
      markers = markers,
      status = status,
      event = "needed",
      method = "expinbase",
      transform = "exp",
      direction = "<",
      cutoff.method = "youden"
    ),
    "Infinity or NaNs values generated in markers, transformation changed to 'none'."
  )
  expect_error(
    mathComb(
      markers = markers3,
      status = status3,
      event = "M",
      method = "distance",
      direction = "auto",
      standardize = "none",
      cutoff.method = "youden"
    ),
    "distance should be one of “euclidean”, “manhattan”, “chebyshev”, “kulczynski_d”, “lorentzian”, “avg”, “taneja”, “kumar-johnson”"
  )
  
  expect_error(
    mathComb(
      markers = markers3,
      status = status3,
      event = "M",
      method = "distance",
      distance = "manhattan",
      direction = "auto",
      standardize = "none",
      transform = "adsd",
      cutoff.method = "youden"
    ),
    "transforms should be one of “none”, “log”, “exp”, “sin”, “cos”"
  )
  expect_error(
    mathComb(
      markers = markers3,
      status = status3,
      event = "M",
      method = "distance",
      distance = "dscsdcs",
      direction = "auto",
      standardize = "none",
      cutoff.method = "youden"
    ),
    "distance should be one of “euclidean”, “manhattan”, “chebyshev”, “kulczynski_d”, “lorentzian”, “avg”, “taneja”, “kumar-johnson”"
  )
  
})

# Markers için numeric Kontrolü ve event statüsü içeriyor mu?
markers3[44, 1:2] <- "assay"

test_that("mathComb functions ...", {
  expect_error(
    mathComb(
      markers = markers3,
      status = status3,
      event = "M",
      method = "add",
      direction = "<",
      standardize = "zScore",
      cutoff.method = "youden"
    ),
    "at least one variable is not numeric"
  )
  expect_error(
    mathComb(
      markers = markers,
      status = status,
      event = "C",
      method = "add",
      direction = "<",
      standardize = "zScore",
      cutoff.method = "youden"
    ),
    "status does not include event"
  )
})

# NA Kontrolü
markers3 <- Data3[, 4:5]
status3[[12]] <- NA



test_that("mathComb functions ...", {
  expect_warning(
    mathComb(
      markers = markers3,
      status = status3,
      event = "M",
      method = "divide",
      direction = "<",
      standardize = "zScore",
      cutoff.method = "youden"
    ),
    "Rows with NA removed from the dataset since status include NA"
  )
})


markers3[44, 1:2] <- NA
status3 <- factor(Data3[, 2], levels = c("B", "M"))

test_that("mathComb functions ...", {
  expect_warning(
    mathComb(
      markers = markers3,
      status = status3,
      event = "M",
      method = "add",
      direction = "<",
      standardize = "zScore",
      cutoff.method = "youden"
    ),
    "Rows with NA removed from the dataset since markers include NA"
  )
})

