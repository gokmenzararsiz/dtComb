data("exampleData1")
data <- exampleData1

markers <- data[, -1]
data$group <- factor(data$group, levels = c("not_needed", "needed"))
event <- "needed"

df <- 4 


gamfit1 = gam::gam(group ~ gam::s(ddimer, df = df) + gam::s(log_leukocyte, df = df), 
              data = data, family = binomial)

pred1 = predict(gamfit1,newdata=markers,type="response")

plot.score <- pROC::roc(data$group ~ pred1, plot = TRUE, print.auc = TRUE,
                        col = "purple", lwd = 4, legacy.axes = TRUE,
                        main = "ROC Curve for Combined Diagnostic Test 1")



gamfit2 = gam::gam(group ~ gam::s(ddimer, df = df) + log_leukocyte, 
                   data = data, family = binomial)

pred2 = predict(gamfit2,newdata=markers,type="response")

plot.score <- pROC::roc(data$group ~ pred2, plot = TRUE, print.auc = TRUE,
                        col = "red", lwd = 4, legacy.axes = TRUE,
                        main = "ROC Curve for Combined Diagnostic Test 2")


gamfit3 = gam::gam(group ~ ddimer +  gam::s(log_leukocyte, df = df), 
                   data = data, family = binomial)

pred3 = predict(gamfit3,newdata=markers,type="response")

plot.score <- pROC::roc(data$group ~ pred3, plot = TRUE, print.auc = TRUE,
                        col = "blue", lwd = 4, legacy.axes = TRUE,
                        main = "ROC Curve for Combined Diagnostic Test 3")

############################################################################################

gamfit4 = gam::gam(group ~ splines::ns(ddimer, df = df) + splines::ns(log_leukocyte, df = df), 
                   data = data, family = binomial)

pred4 = predict(gamfit4,newdata=markers,type="response")

plot.score <- pROC::roc(data$group ~ pred4, plot = TRUE, print.auc = TRUE,
                        col = "green", lwd = 4, legacy.axes = TRUE,
                        main = "ROC Curve for Combined Diagnostic Test 4")



gamfit5 = gam::gam(group ~ splines::ns(ddimer, df = df) + log_leukocyte, 
                   data = data, family = binomial)

pred5 = predict(gamfit5,newdata=markers,type="response")

plot.score <- pROC::roc(data$group ~ pred5, plot = TRUE, print.auc = TRUE,
                        col = "orange", lwd = 4, legacy.axes = TRUE,
                        main = "ROC Curve for Combined Diagnostic Test 5")


gamfit6 = gam::gam(group ~ ddimer +  splines::ns(log_leukocyte, df = df), 
                   data = data, family = binomial)

pred6 = predict(gamfit6,newdata=markers,type="response")

plot.score <- pROC::roc(data$group ~ pred6, plot = TRUE, print.auc = TRUE,
                        col = "pink", lwd = 4, legacy.axes = TRUE,
                        main = "ROC Curve for Combined Diagnostic Test 6")

###################################################################################

