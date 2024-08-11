data("exampleData1")
Data <- exampleData1[-c(83:138), ]
markers <- Data[, -1]
status <- factor(Data$group, levels = c("not_needed", "needed"))

load("result_data/test_std.train.rda")

for (standardize in c(
  "range",
  "zScore",
  "tScore",
  "mean",
  "deviance"
)) {
  res <- std.train(markers, standardize = standardize)

  test_that("std.train functions ...", {
    expect_length(res, 2)
    expect_equal(as.numeric(res$data$ddimer), r$ddimer[r$standardize == standardize],
      tolerance =
        0.01
    )
    expect_equal(as.numeric(res$data$log_leukocyte), r$log_leukocyte[r$standardize == standardize],
      tolerance =
        0.01
    )
  })
}

###############################################################################

test <- exampleData1[c(83:138), -1]

load("result_data/test_std.test.rda")

for (standardize in c(
  "range",
  "zScore",
  "tScore",
  "mean",
  "deviance"
)) {
  res <- linComb(
    markers = markers,
    status = status,
    event = "needed",
    method = "SL",
    resample = "none",
    standardize = standardize,
    direction = "<",
    cutoff.method = "Youden"
  )

  r.std.test <- std.test(test, res)

  test_that("std.test functions ...", {
    expect_length(r.std.test, 2)
    expect_equal(as.numeric(r.std.test$ddimer), r.test$ddimer[r.test$standardize == standardize],
      tolerance =
        0.01
    )
    expect_equal(as.numeric(r.std.test$log_leukocyte), r.test$log_leukocyte[r.test$standardize == standardize],
      tolerance =
        0.01
    )
  })
}
