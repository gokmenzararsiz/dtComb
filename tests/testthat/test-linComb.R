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


#comb.score, AUC, SEN, SPE ve Cutoff kontrolü

load("result_data/test_linComb.rda")

for (method in c("TS",
                 "minimax")) {
  set.seed(14042022)
  rf <- linComb(
    markers = markers,
    status = status,
    event = "needed",
    method = method,
    resample = "none",
    direction = "<",
    cutoff.method = "youden"
  )
  
  test_that("linComb functions ...", {
    expect_length(rf, 11)
    expect_equal(as.numeric(rf$AUC_table$AUC[[3]]),  r$AUC[r$Method == method][1], tolerance =
                   0.01)
    expect_equal(as.numeric(rf$DiagStatCombined$detail$sp[[1]]),
                 r$SPE[r$Method == method][1],
                 tolerance = 0.01)
    expect_equal(as.numeric(rf$DiagStatCombined$detail$se[[1]]),
                 r$SENS[r$Method == method][1],
                 tolerance = 0.01)
    expect_equal(as.numeric(rf$ThresholdCombined), r$Cutoff[r$Method == method][1], tolerance =
                   0.01)
  })
}

#mayo Datası ile
for (method in c("logistic",
                 "SL",
                 "scoring")) {
  set.seed(14042022)
  rf <- linComb(
    markers = markers2,
    status = status2,
    event = "1",
    method = method,
    resample = "none",
    direction = "<",
    cutoff.method = "youden"
  )
  
  test_that("linComb functions ...", {
    expect_length(rf, 11)
    expect_equal(as.numeric(rf$AUC_table$AUC[[3]]),  r$AUC[r$Method == method][1], tolerance =
                   0.01)
    expect_equal(as.numeric(rf$DiagStatCombined$detail$sp[[1]]),
                 r$SPE[r$Method == method][1],
                 tolerance = 0.01)
    expect_equal(as.numeric(rf$DiagStatCombined$detail$se[[1]]),
                 r$SENS[r$Method == method][1],
                 tolerance = 0.01)
    expect_equal(as.numeric(rf$ThresholdCombined), r$Cutoff[r$Method == method][1], tolerance =
                   0.01)
  })
}

#WDBC Datası ile 
for (method in c("PCL",
                 "PT",
                 "minmax"
                 )) {
  set.seed(14042022)
  rf <- linComb(
    markers = markers3,
    status = status3,
    event = "M",
    method = method,
    resample = "none",
    standardize = "range",
    direction = "<",
    cutoff.method = "youden"
  )
  
  test_that("linComb functions ...", {
    expect_length(rf, 11)
    expect_equal(as.numeric(rf$AUC_table$AUC[[3]]),  r$AUC[r$Method == method][1], tolerance =
                   0.01)
    expect_equal(as.numeric(rf$DiagStatCombined$detail$sp[[1]]),
                 r$SPE[r$Method == method][1],
                 tolerance = 0.01)
    expect_equal(as.numeric(rf$DiagStatCombined$detail$se[[1]]),
                 r$SENS[r$Method == method][1],
                 tolerance = 0.01)
    expect_equal(as.numeric(rf$ThresholdCombined), r$Cutoff[r$Method == method][1], tolerance =
                   0.01)
  })
}




#Hata kontrolü
status4 <- factor(Data3[, 2], levels = c("B", "M", "C"))
status4[[9]] <- "C"

test_that("linComb functions ...", {
  expect_error(
    linComb(
      markers = Data3[, 4:5],
      status = status4,
      event = "M",
      method = "scoring",
      direction = direction,
      standardize = "zScore",
      cutoff.method = cutoff.method
    ),
    "the number of status levels should be 2"
  )
  
  expect_error(
    linComb(
      markers = Data3[, 4:6],
      status = status3,
      event = "M",
      method = "PT",
      direction = direction,
      standardize = "zScore",
      cutoff.method = cutoff.method
    ),
    "the number of markers should be 2"
  )
})

test_that("linComb functions ...", {
  expect_error(
    linComb(
      markers = markers3,
      status = status3,
      event = "M",
      direction = "<",
      standardize = "none",
      cutoff.method = "youden"
    ),
    "method should be one of “scoring”, “SL”, “logistic”, “minmax”, “PT”, “PCL”, “minimax”, “TS”"
  )
  
  expect_error(
    linComb(
      markers = markers3,
      status = status3,
      event = "M",
      method = "asaddsa",
      direction = "auto",
      standardize = "none",
      cutoff.method = "youden"
    ),
    "method should be one of “scoring”, “SL”, “logistic”, “minmax”, “PT”, “PCL”, “minimax”, “TS”"
  )
  
  expect_error(
    linComb(
      markers = markers3,
      status = status3,
      event = "M",
      method = "minmax",
      direction = "auto",
      resample = "cv",
      standardize = "asdada",
      cutoff.method = "youden"
    ),
    "standardize should be one of “range”, “zScore”, “tScore”, “mean”, “deviance”"
  )
  
  expect_error(
    linComb(
      markers = markers2,
      status = status2,
      event = "1",
      method = "minimax",
      direction = "asdada",
      standardize = "none",
      cutoff.method = "youden"
    ),
    "direction should be one of “auto”, “<”, “>”"
  )
  
  expect_error(
    linComb(
      markers = markers2,
      status = status2,
      event = "1",
      method = "SL",
      direction = "auto",
      standardize = "tScore",
      cutoff.method = "sadda"
    ),
    "cutoff.method should be one of “youden”, “roc01”"
  )
  expect_error(
    linComb(
      markers = markers2,
      status = status2,
      event = "1",
      method = "scoring",
      resample = "sada",
      standardize = "range",
      direction = "<",
      cutoff.method = "youden"
    ),
    "resample should be one of “none“, “cv“, “repeatedcv“, “boot“"
  )
  
})

# Markers için numeric Kontrolü ve event statüsü içeriyor mu?
markers3[44, 1:2] <- "assay"

test_that("linComb functions ...", {
  expect_error(
    linComb(
      markers = markers3,
      status = status3,
      event = "M",
      method = "PCL",
      direction = "<",
      standardize = "zScore",
      cutoff.method = "youden"
    ),
    "at least one variable is not numeric"
  )
  expect_error(
    linComb(
      markers = markers,
      status = status,
      event = "C",
      method = "PCL",
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



test_that("linComb functions ...", {
  expect_warning(
    linComb(
      markers = markers3,
      status = status3,
      event = "M",
      method = "TS",
      direction = "<",
      standardize = "zScore",
      cutoff.method = "youden"
    ),
    "Rows with NA removed from the dataset since status include NA"
  )
})


markers3[44, 1:2] <- NA
status3 <- factor(Data3[, 2], levels = c("B", "M"))

test_that("linComb functions ...", {
  expect_warning(
    linComb(
      markers = markers3,
      status = status3,
      event = "M",
      method = "TS",
      direction = "<",
      standardize = "zScore",
      cutoff.method = "youden"
    ),
    "Rows with NA removed from the dataset since markers include NA"
  )
})