data("laparoscopy")
Data <- laparoscopy[-c(83:138), ]
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

load("result_data/test_mathComb.rda")

###############################################################################

for (distance in c(
  "lorentzian",
  "avg",
  "taneja",
  "kumar-johnson"
)) {
  set.seed(14042022)
  res <- mathComb(
    markers = markers,
    status = status,
    event = "needed",
    method = "distance",
    distance = distance,
    direction = "<",
    cutoff.method = "Youden"
  )

  test_that("mathComb functions ...", {
    expect_length(res, 15)
    expect_equal(as.numeric(res$CombScore), r$Comb.score[r$Distance == distance],
      tolerance =
        0.1
    )
    expect_equal(as.numeric(res$AUC_table$AUC[[3]]), r$AUC[r$Distance == distance][1],
      tolerance =
        0.01
    )
    expect_equal(as.numeric(res$DiagStatCombined$detail[4, 2]),
      r$SPE[r$Distance == distance][1],
      tolerance = 0.01
    )
    expect_equal(as.numeric(res$DiagStatCombined$detail[3, 2]),
      r$SENS[r$Distance == distance][1],
      tolerance = 0.01
    )
    expect_equal(as.numeric(res$ThresholdCombined), r$Cutoff[r$Distance == distance][1],
      tolerance =
        0.01
    )
  })
}

###############################################################################

for (method in c(
  "add",
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
    cutoff.method = "Youden"
  )

  test_that("mathComb functions ...", {
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

for (distance in c(
  "kulczynski_d",
  "euclidean",
  "manhattan",
  "chebyshev"
)) {
  set.seed(14042022)
  res <- mathComb(
    markers = markers3,
    status = status3,
    event = "M",
    method = "distance",
    distance = distance,
    direction = "<",
    cutoff.method = "Youden"
  )

  test_that("mathComb functions ...", {
    expect_length(res, 15)
    expect_equal(as.numeric(res$CombScore), r$Comb.score[r$Distance == distance],
      tolerance =
        0.1
    )
    expect_equal(as.numeric(res$AUC_table$AUC[[3]]), r$AUC[r$Distance == distance][1],
      tolerance =
        0.01
    )
    expect_equal(as.numeric(res$DiagStatCombined$detail[4, 2]),
      r$SPE[r$Distance == distance][1],
      tolerance = 0.01
    )
    expect_equal(as.numeric(res$DiagStatCombined$detail[3, 2]),
      r$SENS[r$Distance == distance][1],
      tolerance = 0.01
    )
    expect_equal(as.numeric(res$ThresholdCombined), r$Cutoff[r$Distance == distance][1],
      tolerance =
        0.01
    )
  })
}

###############################################################################

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
      cutoff.method = "Youden"
    ),
    "method should be one of 'add', 'multiply', 'divide', 'subtract' ,'distance', 'baseinexp', 'expinbase'"
  )

  expect_error(
    mathComb(
      markers = markers3,
      status = status3,
      event = "M",
      method = "adsadad",
      direction = "auto",
      standardize = "none",
      cutoff.method = "Youden"
    ),
    "method should be one of 'add', 'multiply', 'divide', 'subtract' ,'distance', 'baseinexp', 'expinbase'"
  )

  expect_error(
    mathComb(
      markers = markers3,
      status = status3,
      event = "M",
      method = "add",
      direction = "auto",
      standardize = "asdada",
      cutoff.method = "Youden"
    ),
    "standardize should be one of 'range', 'zScore', 'tScore', 'mean', 'deviance'"
  )

  expect_error(
    mathComb(
      markers = markers2,
      status = status2,
      event = "1",
      method = "baseinexp",
      direction = "asdada",
      standardize = "none",
      cutoff.method = "Youden"
    ),
    "direction should be one of 'auto', '<', '>'"
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
    "The entered cutoff.method is invalid"
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
      cutoff.method = "Youden"
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
      cutoff.method = "Youden"
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
      cutoff.method = "Youden"
    ),
    "distance should be one of 'euclidean', 'manhattan', 'chebyshev', 'kulczynski_d', 'lorentzian', 'avg', 'taneja', 'kumar-johnson'"
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
      cutoff.method = "Youden"
    ),
    "transforms should be one of 'none', 'log', 'exp', 'sin', 'cos'"
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
      cutoff.method = "Youden"
    ),
    "distance should be one of 'euclidean', 'manhattan', 'chebyshev', 'kulczynski_d', 'lorentzian', 'avg', 'taneja', 'kumar-johnson'"
  )
})

###############################################################################

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
      cutoff.method = "Youden"
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
      cutoff.method = "Youden"
    ),
    "status does not include event"
  )
})

###############################################################################

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
      cutoff.method = "Youden"
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
      cutoff.method = "Youden"
    ),
    "Rows with NA removed from the dataset since markers include NA"
  )
})
