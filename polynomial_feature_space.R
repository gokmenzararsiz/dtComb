
data("exampleData1")
data <- exampleData1

markers <- data[, -1]
data$group <- factor(data$group, levels = c("not_needed", "needed"))
event <- "needed"

#data <- cbind(status,markers)
degree <- 2 

# plot(data$ddimer, data$log_leukocyte, col=data$group, xlab = "ddimer", 
#      ylab = "log_leukocyte", las = 1, xlim = c(0, 20))


interact <- markers[,1] * markers[,2]

polyFitin <- glm(group ~ poly(markers[, 1], degree) + poly(markers[,2], degree) 
                 + interact, data = data, family = binomial(link = "logit"))

polyFit <- glm(group ~ poly(ddimer, degree) + poly(log_leukocyte, degree), 
               data = data, family = binomial(link = "logit"))
summary(polyFit)
summary(polyFitin)

pred<-predict(polyFit,newdata=data,type="response")
predin<-predict(polyFitin,newdata=data,type="response")


plot.score <- pROC::roc(data$group ~ pred, plot = TRUE, print.auc = TRUE,
                        col = "purple", lwd = 4, legacy.axes = TRUE,
                        main = "ROC Curve for Combined Diagnostic Test")

plot.score <- pROC::roc(data$group ~ predin, plot = TRUE, print.auc = TRUE,
                        col = "red", lwd = 4, legacy.axes = TRUE,
                        main = "ROC Curve for Combined Diagnostic Test")




