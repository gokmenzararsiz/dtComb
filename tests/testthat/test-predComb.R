library(APtools)
data("exampleData1")
Data <- exampleData1[-c(83:138),]
markers <- Data[, -1]
status <- factor(Data$group, levels = c("not_needed", "needed"))

newmarkers <- exampleData1[c(83:138), -1]



data(mayo)
Data2 <- mayo[-c(42:119), ]
markers2 <- Data2[, 3:4]
status2 <- factor(Data2[, 2], levels = c(1, 0))

newmarkers2 <- mayo[c(42:119), 3:4]


Data3 <-
  read.csv(
    "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data",
    header = FALSE
  )
Data3 <- Data3[-c(121:262),]
markers3 <- Data3[, 4:5]
status3 <- factor(Data3[, 2], levels = c("B", "M"))

newmarkers3 <- Data3[c(121:262), 4:5]


r <- read.table("C:/Users/ilayd/Desktop/projects/test/p.result.txt", header = TRUE)

##################################

for (method in c("scoring",
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
 
  pred <- predComb(rf, newmarkers)
    test_that("linComb functions ...", {
      expect_length(pred, 2)
      expect_equal(as.numeric(pred$comb.score),  r$Comb.score[r$Method == method], tolerance = 0.01)
      expect_equal(pred$labels,
                   r$Labels[r$Method == method],
                   tolerance = 0.01)
    })
}

#mayo Datası ile
for (method in c("logistic",
                 "SL",
                 "TS")) {
  set.seed(14042022)
  rf <- nonlinComb(
    markers = markers2,
    status = status2,
    event = "1",
    method = "nsgam",
    resample = "none",
    direction = "<",
    cutoff.method = "youden"
  )
  
  pred <- predComb(rf, newmarkers2)
  test_that("linComb functions ...", {
    expect_length(pred, 2)
    expect_equal(as.numeric(pred$comb.score),  r$Comb.score[r$Method == method], tolerance = 0.01)
    expect_equal(pred$labels,
                 r$Labels[r$Method == method],
                 tolerance = 0.01)
  })
}

#WDBC Datası ile 
for (method in c("PCL",
                 "PT",
                 "minmax")) {
  set.seed(14042022)
  rf <- linComb(
    markers = markers3,
    status = status3,
    event = "M",
    method = method,
    resample = "none",
    direction = "<",
    cutoff.method = "youden"
  )
  
  pred <- predComb(rf, newmarkers3)
  test_that("linComb functions ...", {
    expect_length(pred, 2)
    expect_equal(as.numeric(pred$comb.score),  r$Comb.score[r$Method == method], tolerance = 0.01)
    expect_equal(pred$labels,
                 r$Labels[r$Method == method],
                 tolerance = 0.01)
  })
}

