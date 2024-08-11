data("exampleData1")
Data <- exampleData1[-c(83:138), ]
markers <- Data[, -1]
status <- factor(Data$group, levels = c("not_needed", "needed"))

test <- exampleData1[c(83:138), ]

load("result_data/mayo.rda")
Data2 <- mayo[-c(42:119), ]
markers2 <- Data2[, 3:4]
status2 <- factor(Data2[, 2], levels = c(1, 0))

Data3 <-
  utils::read.csv(
    "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data",
    header = FALSE
  )
Data3 <- Data3[-c(121:262), ]
markers3 <- Data3[, 4:5]
status3 <- factor(Data3[, 2], levels = c("B", "M"))

###############################################################################

load("result_data/test_linComb.rda")

for (method in c(
  "TS",
  "minimax"
)) {
  set.seed(14042022)
  res <- linComb(
    markers = markers,
    status = status,
    event = "needed",
    method = method,
    resample = "none",
    direction = "<",
    cutoff.method = "Youden"
  )

  test_that("linComb functions ...", {
    expect_length(res, 15)
    expect_equal(as.numeric(res$CombScore), r$Comb.score[r$Method == method],
      tolerance =
        0.1
    )
    expect_equal(as.numeric(res$AUC_table$AUC[[3]]), r$AUC[r$Method == method][1],
      tolerance =
        0.01
    )
    expect_equal(as.numeric(res$DiagStatCombined$detail[4, 2]),
      r$SPE[r$Method == method][1],
      tolerance = 0.01
    )
    expect_equal(as.numeric(res$DiagStatCombined$detail[3, 2]),
      r$SENS[r$Method == method][1],
      tolerance = 0.01
    )
    expect_equal(as.numeric(res$ThresholdCombined), r$Cutoff[r$Method == method][1],
      tolerance =
        0.01
    )
  })
}

###############################################################################

for (method in c(
  "logistic",
  "SL",
  "scoring"
)) {
  set.seed(14042022)
  res <- linComb(
    markers = markers2,
    status = status2,
    event = "1",
    method = method,
    resample = "none",
    direction = "<",
    cutoff.method = "Youden"
  )

  test_that("linComb functions ...", {
    expect_length(res, 15)
    expect_equal(as.numeric(res$CombScore), r$Comb.score[r$Method == method],
      tolerance =
        0.1
    )
    expect_equal(as.numeric(res$AUC_table$AUC[[3]]), r$AUC[r$Method == method][1],
      tolerance =
        0.01
    )
    expect_equal(as.numeric(res$DiagStatCombined$detail[4, 2]),
      r$SPE[r$Method == method][1],
      tolerance = 0.01
    )
    expect_equal(as.numeric(res$DiagStatCombined$detail[3, 2]),
      r$SENS[r$Method == method][1],
      tolerance = 0.01
    )
    expect_equal(as.numeric(res$ThresholdCombined), r$Cutoff[r$Method == method][1],
      tolerance =
        0.01
    )
  })
}

###############################################################################

for (method in c(
  "PCL",
  "PT",
  "minmax"
)) {
  set.seed(14042022)
  res <- linComb(
    markers = markers3,
    status = status3,
    event = "M",
    method = method,
    resample = "none",
    standardize = "range",
    direction = "<",
    cutoff.method = "Youden"
  )

  test_that("linComb functions ...", {
    expect_length(res, 15)
    expect_equal(as.numeric(res$CombScore), r$Comb.score[r$Method == method],
      tolerance =
        0.1
    )
    expect_equal(as.numeric(res$AUC_table$AUC[[3]]), r$AUC[r$Method == method][1],
      tolerance =
        0.01
    )
    expect_equal(as.numeric(res$DiagStatCombined$detail[4, 2]),
      r$SPE[r$Method == method][1],
      tolerance = 0.01
    )
    expect_equal(as.numeric(res$DiagStatCombined$detail[3, 2]),
      r$SENS[r$Method == method][1],
      tolerance = 0.01
    )
    expect_equal(as.numeric(res$ThresholdCombined), r$Cutoff[r$Method == method][1],
      tolerance =
        0.01
    )
  })
}

###############################################################################

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
      cutoff.method = "Youden"
    ),
    "method should be one of 'scoring', 'SL', 'logistic', 'minmax', 'PT', 'PCL', 'minimax', 'TS'"
  )

  expect_error(
    linComb(
      markers = markers3,
      status = status3,
      event = "M",
      method = "asaddsa",
      direction = "auto",
      standardize = "none",
      cutoff.method = "Youden"
    ),
    "method should be one of 'scoring', 'SL', 'logistic', 'minmax', 'PT', 'PCL', 'minimax', 'TS'"
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
      cutoff.method = "Youden"
    ),
    "standardize should be one of 'range', 'zScore', 'tScore', 'mean', 'deviance'"
  )

  expect_error(
    linComb(
      markers = markers2,
      status = status2,
      event = "1",
      method = "minimax",
      direction = "asdada",
      standardize = "none",
      cutoff.method = "Youden"
    ),
    "direction should be one of 'auto', '<', '>'"
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
    "The entered cutoff.method is invalid"
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
      cutoff.method = "Youden"
    ),
    "resample should be one of 'none', 'cv', 'repeatedcv', 'boot'"
  )
})

###############################################################################

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
      cutoff.method = "Youden"
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
      cutoff.method = "Youden"
    ),
    "status does not include event"
  )
})

###############################################################################

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
      cutoff.method = "Youden"
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
      cutoff.method = "Youden"
    ),
    "Rows with NA removed from the dataset since markers include NA"
  )
})

test_that("linComb functions ...", {
  expect_warning(
    linComb(
      markers = markers,
      status = status,
      event = "needed",
      method = "PCL",
      direction = "<",
      cutoff.method = "Youden"
    ),
    "The used combination method requires range standardization. All biomarker values are standardized to a range between 0 and 1."
  )
})

test_that("linComb functions ...", {
  expect_warning(
    linComb(
      markers = markers,
      status = status,
      event = "needed",
      method = "PT",
      direction = "<",
      cutoff.method = "Youden"
    ),
    "The used combination method requires zScore standardization."
  )
})
