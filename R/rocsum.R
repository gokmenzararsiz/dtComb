rocsum <- function(markers = NULL, comb.score = NULL, status = NULL, event = NULL, 
                    direction = c("auto", "<", ">"), conf.level = 0.95, 
                    cutoff.method = c("youden", "roc01")){
  
  roc.m1 <- suppressMessages(pROC::roc(status ~ markers[, 1], plot = TRUE, print.auc = TRUE, 
                      direction = direction, col = "#619CFF", lwd = 4, legacy.axes = TRUE, 
                      percent = FALSE, main = "ROC Curves for Combined Diagnostic Test"))
  
  roc.m2 <- suppressMessages(pROC::roc(status ~ markers[, 2], plot = TRUE, print.auc = TRUE, 
                      direction = direction, print.auc.y = 0.40, col = "#00BA38", lwd = 4, 
                      legacy.axes = TRUE, percent = FALSE, add = TRUE,
                      main = "ROC Curves for Combined Diagnostic Test"))
  
  roc.c <- suppressMessages(pROC::roc(status ~ as.numeric(comb.score), plot = TRUE, print.auc = TRUE, 
                     direction = direction, print.auc.y = 0.30, col = "#F8766D", lwd = 4, 
                     legacy.axes = TRUE, percent = FALSE, add = TRUE,
                     main = "ROC Curves for Combined Diagnostic Test"))
  
  legend("bottomright", legend = c(colnames(markers)[1], colnames(markers)[2], 
                                   "Combined Score"), col = c("#619CFF","#00BA38","#F8766D"), lwd = 4)
  
  
  # ROC coordinates
  coord.m1 <- pROC::coords(roc.m1)
  coord.m2 <- pROC::coords(roc.m2)
  coord.c <- pROC::coords(roc.c)
  coord.names <- c(rep(colnames(markers)[1], nrow(coord.m1)),
                   rep(colnames(markers)[2], nrow(coord.m2)),
                   rep("Combined", nrow(coord.c)))
  ROC_coordinates <- data.frame(coord.names, rbind(coord.m1, coord.m2, coord.c))
  colnames(ROC_coordinates) <- c("Marker", "Threshold", "Specificity", "Sensitivity")
  
  
  # Area under the ROC Curves
  auc.m1 <- pROC::ci.auc(roc.m1, method = "delong")
  auc.m2 <- pROC::ci.auc(roc.m2, method = "delong")
  auc.c <- pROC::ci.auc(roc.c, method = "delong")
  AUC_table <- rbind(auc.m1, auc.m2, auc.c)
  
  var.m1 <- pROC::var(roc.m1, method = "delong")
  var.m2 <- pROC::var(roc.m2, method = "delong")
  var.c <- pROC::var(roc.c, method = "delong")
  
  std.err <- c(sqrt(var.m1), sqrt(var.m2), sqrt(var.c))
  z.stat <- (AUC_table[, 2] - .5) / std.err
  p.val <- 2 * pt(-abs(z.stat), df = Inf)
  
  AUC_table <- cbind(AUC_table[, 2], std.err, AUC_table[, 1], AUC_table[, 3], 
                     z.stat, p.val)
  rownames(AUC_table) <- c(colnames(markers)[1], colnames(markers)[2], "Combined")
  colnames(AUC_table) <- c("AUC", "SE.AUC", "LowerLimit", "UpperLimit", "z", "p-value")
  AUC_table <- data.frame(AUC_table)
  
  
  # Comparison of the Area under the ROC curves
  roccm1 <- pROC::roc.test(roc.c, roc.m1, method = "delong")
  roccm2 <- pROC::roc.test(roc.c, roc.m2, method = "delong")
  rocm1m2 <- pROC::roc.test(roc.m1, roc.m2, method = "delong")
  
  MultComp_table <- matrix(0, 3, 6)
  MultComp_table[1, ] <- cbind(roccm1$estimate[1], 
                               roccm1$estimate[2], abs(roccm1$estimate[1] - roccm1$estimate[2]), 
                               abs(roccm1$estimate[1] - roccm1$estimate[2]) / roccm1$statistic, 
                               roccm1$statistic, roccm1$p.value)
  MultComp_table[2, ] <- cbind(roccm2$estimate[1], 
                               roccm2$estimate[2], abs(roccm2$estimate[1] - roccm2$estimate[2]), 
                               abs(roccm2$estimate[1] - roccm2$estimate[2]) / roccm2$statistic, 
                               roccm2$statistic, roccm2$p.value)
  MultComp_table[3, ] <- cbind(rocm1m2$estimate[1], 
                               rocm1m2$estimate[2], abs(rocm1m2$estimate[1] - rocm1m2$estimate[2]), 
                               abs(rocm1m2$estimate[1] - rocm1m2$estimate[2]) / rocm1m2$statistic, 
                               rocm1m2$statistic, rocm1m2$p.value)
  
  comp.names <- cbind(c("Combined", "Combined", colnames(markers)[1]), 
                      c(colnames(markers)[1], colnames(markers)[2], colnames(markers)[2]))
  
  
  MultComp_table <- data.frame(comp.names, MultComp_table)
  colnames(MultComp_table) <- c("Marker1 (A)", "Marker2 (B)", "AUC (A)", "AUC (B)",
                                "|A-B|", "SE(|A-B|)", "z", "p-value")
  
  # Diagnostic Statistics for Each Marker and Combined Function
  best.m1 <- pROC::coords(roc.m1, "best", ret = c("threshold","tp","tn","fp","fn"), 
                          best.method = cutoff.method)
  best.m2 <- pROC::coords(roc.m2, "best", ret = c("threshold","tp","tn","fp","fn"), 
                          best.method = cutoff.method)
  best.c <- pROC::coords(roc.c, "best", ret = c("threshold","tp","tn","fp","fn"), 
                         best.method = cutoff.method)
  best.m1.tbl <- as.table(matrix(c(best.m1$tp, best.m1$fp, best.m1$fn, best.m1$tn), 
                                 nrow = 2, byrow = TRUE))
  best.m2.tbl <- as.table(matrix(c(best.m2$tp, best.m2$fp, best.m2$fn, best.m2$tn), 
                                 nrow = 2, byrow = TRUE))
  best.c.tbl <- as.table(matrix(c(best.c$tp, best.c$fp, best.c$fn, best.c$tn), 
                                nrow = 2, byrow = TRUE))
  DiagStatMarker1 <- epiR::epi.tests(best.m1.tbl, conf.level = conf.level)
  DiagStatMarker2 <- epiR::epi.tests(best.m2.tbl, conf.level = conf.level)
  DiagStatCombined <- epiR::epi.tests(best.c.tbl, conf.level = conf.level)
  
  allres <- list(ROC_coordinates = ROC_coordinates, 
                 AUC_table = AUC_table,
                 MultComp_table = MultComp_table,
                 DiagStatMarker1 = DiagStatMarker1,
                 DiagStatMarker2 = DiagStatMarker2,
                 DiagStatCombined = DiagStatCombined,
                 ThresholdCombined = best.c$threshold)
  
  return(allres)
}
