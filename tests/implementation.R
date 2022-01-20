#linComb

data(exampleData1)

markers <- exampleData1[, -1]
status <- factor(exampleData1$group, levels = c("not_needed", "needed"))
event <- "needed"

score3 <- linComb(markers = markers, status = status, event = event,
              method = "logistic", direction = "<", cutoff.method = "youden")

score2 <- linComb(markers = markers, status = status, event = event,
                 method = "minmax", direction = "<", 
                  cutoff.method = "youden")

score1 <- linComb(markers = markers, status = status, event = event,
                  method = "scoring", ndigits = 0, standardize = "zScore", direction = "<", 
                  cutoff.method = "youden")


#mlComb

data(exampleData1)
markers <- exampleData1[, -1]
status <- factor(exampleData1$group, levels = c("not_needed", "needed"))
method <- "glm"
train_method <- "cv"
preProcess = c("center","scale")
train_repeats <- 1
event <- "needed"
train_number <- 10
direction <- "<"
cutoff.method <- "youden"

method <- "glm"

model <- mlComb(markers = markers, status = status, event = event,
                method = method,
                preProcess = preProcess,
                direction = direction, cutoff.method = cutoff.method)

method <- "glm"

model <- mlComb(markers = markers, status = status, event = event,
                method = method,
                direction = direction, cutoff.method = cutoff.method)

method <- "glm"
train_method <- "cv"
train_number <- 10

model <- mlComb(markers = markers, status = status, event = event,
                          method = method,
                          train_method = train_method, train_number = train_number,
                          preProcess = preProcess,
                          direction = direction, cutoff.method = cutoff.method)

method <- "glm"
train_method <- "repeatedcv"
train_repeats <- 5
train_number <- 10

model <- mlComb(markers = markers, status = status, event = event,
                method = method,
                train_method = train_method, train_number = train_number,
                train_repeats = train_repeats,
                preProcess = preProcess,
                direction = direction, cutoff.method = cutoff.method)


method <- "ordinalRF"

model <- mlComb(markers = markers, status = status, event = event,
                method = method,
                preProcess = preProcess,
                direction = direction, cutoff.method = cutoff.method)

method <- "gbm"
train_method <- "cv"
train_number <- 10

model <- mlComb(markers = markers, status = status, event = event,
                method = method,
                train_method = train_method, train_number = train_number,
                preProcess = preProcess,
                direction = direction, cutoff.method = cutoff.method)

method <- "gbm"
train_method <- "repeatedcv"
train_repeats <- 5
train_number <- 10

model <- mlComb(markers = markers, status = status, event = event,
                method = method,
                train_method = train_method, train_number = train_number,
                train_repeats = train_repeats,
                preProcess = preProcess,
                direction = direction, cutoff.method = cutoff.method)

method <- "bagFDA"
train_method <- "cv"
train_number <- 10

model <- mlComb(markers = markers, status = status, event = event,
                method = method,
                train_method = train_method, train_number = train_number,
                preProcess = preProcess,
                direction = direction, cutoff.method = cutoff.method)

method <- "bagFDA"
train_method <- "cv"
train_number <- 10

model <- mlComb(markers = markers, status = status, event = event,
                method = method,
                train_method = train_method, train_number = train_number,
                preProcess = preProcess, B=25,
                direction = direction, cutoff.method = cutoff.method)


method <- "bagFDA"
train_method <- "repeatedcv"
train_repeats <- 5
train_number <- 10

model <- mlComb(markers = markers, status = status, event = event,
                method = method,
                train_method = train_method, train_number = train_number,
                train_repeats = train_repeats,
                preProcess = preProcess, B=15,
                direction = direction, cutoff.method = cutoff.method)

method <- "bag" #will give errors
train_method <- "cv"
train_number <- 10

model <- mlComb(markers = markers, status = status, event = event,
                method = method,
                train_method = train_method, train_number = train_number,
                preProcess = preProcess,
                direction = direction, cutoff.method = cutoff.method)

method <- "bag" #parameter added to eliminate error
train_method <- "cv"
train_number <- 10

model <- mlComb(markers = markers, status = status, event = event,
                method = method,
                train_method = train_method, train_number = train_number,
                preProcess = preProcess,
                direction = direction, cutoff.method = cutoff.method, 
                bagControl = bagControl(fit = ldaBag$fit, predict = ldaBag$pred,
                                        aggregate = ldaBag$aggregate))



#mathComb

 data(exampleData1)
 markers <- exampleData1[, -1]
 status <- factor(exampleData1$group, levels = c("not_needed", "needed"))
 event <- "needed"
 direction <- "<"
 cutoff.method <- "youden"

 score1 <- mathComb(markers = markers, status = status, event = event,
 method = "distance", distance ="euclidean", direction = direction,
 cutoff.method = cutoff.method)

 score2 <- mathComb(markers = markers, status = status, event = event,
 method = "add", log.transform = TRUE, direction = direction, cutoff.method = cutoff.method)

 score2 <- mathComb(markers = markers, status = status, event = event,
                    method = "subtract", power.transform = TRUE, direction = direction,
                    cutoff.method = cutoff.method)
