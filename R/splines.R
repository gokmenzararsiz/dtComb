data("exampleData1")
data <- exampleData1

markers <- data[, -1]
data$group <- factor(data$group, levels = c("not_needed", "needed"))
event <- "needed"

df <- 4 
degree <- 1

bspline = glm(group ~ splines::bs(ddimer,degree = degree, df = df) + 
                splines::bs(log_leukocyte,degree = degree, df = df), 
              data = data, family = binomial)

bpred = predict(bspline,newdata=markers,type="response")

plot.score <- pROC::roc(data$group ~ bpred, plot = TRUE, print.auc = TRUE,
                        col = "purple", lwd = 4, legacy.axes = TRUE,
                        main = "ROC Curve for Combined Diagnostic Test")

#############################################################################

nspline = glm(group ~ splines::ns(ddimer,df=df) + splines::ns(log_leukocyte,df=df),
              data=data,family=binomial)

npred = predict(nspline,newdata=markers,type="response")

plot.score <- pROC::roc(data$group ~ npred, plot = TRUE, print.auc = TRUE,
                        col = "red", lwd = 4, legacy.axes = TRUE,
                        main = "ROC Curve for Combined Diagnostic Test")
