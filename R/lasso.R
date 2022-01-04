data("exampleData1")
data <- exampleData1

markers <- data[, -1]
data$group <- factor(data$group, levels = c("not_needed", "needed"))
event <- "needed"

#data <- cbind(status,markers)
degree <- 5 

interact <- markers[,1] * markers[,2]

space <- cbind.data.frame(poly(markers[, 1], degree), poly(markers[,2], degree))


set.seed(123) 
cv.lasso <- glmnet::cv.glmnet(x = as.matrix(space), y = data$group, alpha = 1, family = "binomial")
# Fit the final model on the training data
model <- glmnet::glmnet(x = space, y = data$group, alpha = 1, family = "binomial",
                lambda = cv.lasso$lambda.min)
# Display regression coefficients
coef(model)


spacein <- cbind.data.frame(poly(markers[, 1], degree), poly(markers[,2], degree), interact)


set.seed(123) 
cv.lassoin <- glmnet::cv.glmnet(x = as.matrix(spacein), y = data$group, alpha = 1, family = "binomial")
# Fit the final model on the training data
modelin <- glmnet::glmnet(x = spacein, y = data$group, alpha = 1, family = "binomial",
                lambda = cv.lassoin$lambda.min)
# Display regression coefficients
coef(modelin)

pred <- predict(model, newx = as.matrix(space), type="response")
predin <- predict(modelin, newx = as.matrix(spacein), type="response")


plot.score <- pROC::roc(data$group ~ pred, plot = TRUE, print.auc = TRUE,
                        col = "purple", lwd = 4, legacy.axes = TRUE,
                        main = "ROC Curve for Combined Diagnostic Test")

plot.score <- pROC::roc(data$group ~ predin, plot = TRUE, print.auc = TRUE,
                        col = "red", lwd = 4, legacy.axes = TRUE,
                        main = "ROC Curve for Combined Diagnostic Test")

###############################################################################

space <- cbind.data.frame(poly(markers[, 1], degree), poly(markers[,2], degree))


set.seed(123) 
cv.lasso <- glmnet::cv.glmnet(x = as.matrix(space), y = data$group, alpha = 1, family = "binomial")
# Fit the final model on the training data
model <- glmnet::glmnet(x = space, y = data$group, alpha = 1, family = "binomial",
                lambda = cv.lasso$lambda.1se)
# Display regression coefficients
coef(model)


spacein <- cbind.data.frame(poly(markers[, 1], degree), poly(markers[,2], degree), interact)


set.seed(123) 
cv.lassoin <- glmnet::cv.glmnet(x = as.matrix(spacein), y = data$group, alpha = 1, family = "binomial")
# Fit the final model on the training data
modelin <- glmnet::glmnet(x = spacein, y = data$group, alpha = 1, family = "binomial",
                  lambda = cv.lassoin$lambda.1se)
# Display regression coefficients
coef(modelin)

pred <- predict(model, newx = as.matrix(space), type="response")
predin <- predict(modelin, newx = as.matrix(spacein), type="response")


plot.score <- pROC::roc(data$group ~ pred, plot = TRUE, print.auc = TRUE,
                        col = "green", lwd = 4, legacy.axes = TRUE,
                        main = "ROC Curve for Combined Diagnostic Test")

plot.score <- pROC::roc(data$group ~ predin, plot = TRUE, print.auc = TRUE,
                        col = "blue", lwd = 4, legacy.axes = TRUE,
                        main = "ROC Curve for Combined Diagnostic Test")

