# library(usethis)

data("laparoscopy")
Data <- laparoscopy[-c(83:138), ]
markers <- Data[, -1]
status <- Data$group

newmarkers <- laparoscopy[c(83:138), -1]

load("result_data/mayo.rda")
Data2 <- mayo[-c(42:119), ]
markers2 <- Data2[, 3:4]
status2 <- factor(Data2[, 2])

newmarkers2 <- mayo[c(42:119), 3:4]

Data3 <-
  read.csv(
    "result_data/wdbc.data.txt",
    header = FALSE
  )
Data3 <- Data3[-c(121:262), ]
markers3 <- Data3[, 4:5]
status3 <- factor(Data3[, 2], levels = c("B", "M"))

newmarkers3 <- Data3[c(121:262), 4:5]

load("result_data/test_predComb.rda")

###############################################################################

for (method in c(
  "scoring",
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
  pred <- predict(res, newmarkers)
  test_that("linComb functions ...", {
    expect_length(pred, 2)
    expect_equal(as.numeric(pred$comb.score), r$Comb.score[r$Method == method], tolerance = 0.01)
    expect_equal(pred$labels,
      r$Labels[r$Method == method],
      tolerance = 0.01
    )
  })
}

###############################################################################

for (method in c(
  "logistic",
  "SL",
  "TS"
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

  pred <- predict(res, newmarkers2)
  test_that("linComb functions ...", {
    expect_length(pred, 2)
    expect_equal(as.numeric(pred$comb.score), r$Comb.score[r$Method == method], tolerance = 0.01)
    expect_equal(pred$labels,
      r$Labels[r$Method == method],
      tolerance = 0.01
    )
  })
}

###############################################################################

for (method in c(
  "PCL",
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

  pred <- predict(res, newmarkers3)
  test_that("linComb functions ...", {
    expect_length(pred, 2)
    expect_equal(as.numeric(pred$comb.score), r$Comb.score[r$Method == method], tolerance = 0.01)
    expect_equal(pred$labels,
      r$Labels[r$Method == method],
      tolerance = 0.01
    )
  })
}

set.seed(14042022)
res <- linComb(
  markers = markers3,
  status = status3,
  event = "M",
  method = "PT",
  resample = "none",
  standardize = "range",
  direction = "<",
  cutoff.method = "Youden"
)

pred <- predict(res, newmarkers3)
test_that("linComb functions ...", {
  expect_length(pred, 2)
  expect_equal(as.numeric(pred$comb.score), r$Comb.score[r$Method == "PT"], tolerance = 0.01)
  expect_equal(pred$labels,
    r$Labels[r$Method == "PT"],
    tolerance = 0.01
  )
})
