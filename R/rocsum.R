#' @title Generate ROC curves and related  statistics for the given markers and
#' Combination score.
#'
#' @description The \code{rocsum} function returns the ROC curves with
#' coordinates, Area Under the Curves of markers and combination score, Area Under
#' the Curve comparison of markers and combination score, Confusion matrices for both
#' markers and combination score with the cutoff values derived from the ROC Curves.
#'
#' @param markers a \code{numeric} data frame that includes two diagnostic tests
#' results
#'
#' @param comb.score a matrix of \code{numeric} combination scores calculated
#' according to the given method
#'
#' @param status a \code{factor} vector that includes the actual disease
#'  status of the patients
#'
#' @param event a \code{character} string that indicates the event in the status
#' to be considered as positive event
#'
#' @param direction a \code{character} string determines in which direction the
#' comparison will be made.  “>”: if the predictor values for the control group
#' are higher than the values of the case group (controls > cases).
#' “<”: if the predictor values for the control group are lower or equal than
#' the values of the case group (controls < cases).
#'
#' @param conf.level a \code{numeric} values determines the confidens interval
#' for the ROC curve(0.95, default).
#'
#' @param cutoff.method  a \code{character} string determines the cutoff method
#' for the ROC curve.
#'
#' @param show.plot a \code{logical}. If TRUE, a ROC curve is plotted.
#' Default is FALSE.
#'
#' @return A list of \code{numeric} ROC Curves, AUC statistics and Confusion
#' matrices.
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec, Necla Kochan,
#' Gozde Erturk Zararsiz, Selcuk Korkmaz, Gokmen Zararsiz
#'

rocsum <- function(markers = NULL, comb.score = NULL, status = NULL, event = NULL,
                   direction = c("auto", "<", ">"), conf.level = 0.95,
                   cutoff.method = c(
                     "CB", "MCT", "MinValueSp", "MinValueSe", "ValueSp",
                     "ValueSe", "MinValueSpSe", "MaxSp", "MaxSe",
                     "MaxSpSe", "MaxProdSpSe", "ROC01", "SpEqualSe",
                     "Youden", "MaxEfficiency", "Minimax", "MaxDOR",
                     "MaxKappa", "MinValueNPV", "MinValuePPV", "ValueNPV",
                     "ValuePPV", "MinValueNPVPPV", "PROC01", "NPVEqualPPV",
                     "MaxNPVPPV", "MaxSumNPVPPV", "MaxProdNPVPPV",
                     "ValueDLR.Negative", "ValueDLR.Positive", "MinPvalue",
                     "ObservedPrev", "MeanPrev", "PrevalenceMatching"
                   ), show.plot = show.plot) {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  graphics::par(pty = "s")

  roc.m1 <- suppressMessages(pROC::roc(status ~ markers[, 1],
    plot = show.plot, print.auc = TRUE,
    direction = direction, col = "#619CFF", lwd = 6, legacy.axes = TRUE, grid = TRUE,
    percent = FALSE, main = "ROC Curves for Combination Diagnostic Test"
  ))

  roc.m2 <- suppressMessages(pROC::roc(status ~ markers[, 2],
    plot = show.plot, print.auc = TRUE,
    direction = direction, print.auc.y = 0.40, col = "#00BA38", lwd = 6,
    legacy.axes = TRUE, percent = FALSE, add = TRUE,
    main = "ROC Curves for Combination Diagnostic Test"
  ))

  roc.c <- suppressMessages(pROC::roc(status ~ as.numeric(comb.score),
    plot = show.plot, print.auc = TRUE,
    direction = direction, print.auc.y = 0.30, col = "#F8766D", lwd = 6,
    legacy.axes = TRUE, percent = FALSE, add = TRUE,
    main = "ROC Curves for Combination Diagnostic Test"
  ))

  if (show.plot == TRUE) {
    legend("bottomright", legend = c(
      colnames(markers)[1], colnames(markers)[2],
      "Combination Score"
    ), col = c("#619CFF", "#00BA38", "#F8766D"), lwd = 4)
  }
  coord.m1 <- pROC::coords(roc.m1)
  coord.m2 <- pROC::coords(roc.m2)
  coord.c <- pROC::coords(roc.c)
  coord.names <- c(
    rep(colnames(markers)[1], nrow(coord.m1)),
    rep(colnames(markers)[2], nrow(coord.m2)),
    rep("Combination", nrow(coord.c))
  )
  ROC_coordinates <- data.frame(coord.names, rbind(coord.m1, coord.m2, coord.c))
  colnames(ROC_coordinates) <- c("Marker", "Threshold", "Specificity", "Sensitivity")
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

  AUC_table <- cbind(
    AUC_table[, 2], std.err, AUC_table[, 1], AUC_table[, 3],
    z.stat, p.val
  )
  rownames(AUC_table) <- c(colnames(markers)[1], colnames(markers)[2], "Combination")
  colnames(AUC_table) <- c("AUC", "SE.AUC", "LowerLimit", "UpperLimit", "z", "p-value")
  AUC_table <- data.frame(AUC_table)
  roccm1 <- pROC::roc.test(roc.c, roc.m1, method = "delong")
  roccm2 <- pROC::roc.test(roc.c, roc.m2, method = "delong")
  rocm1m2 <- pROC::roc.test(roc.m1, roc.m2, method = "delong")

  MultComp_table <- matrix(0, 3, 6)
  MultComp_table[1, ] <- cbind(
    roccm1$estimate[1],
    roccm1$estimate[2], abs(roccm1$estimate[1] - roccm1$estimate[2]),
    abs(roccm1$estimate[1] - roccm1$estimate[2]) / roccm1$statistic,
    roccm1$statistic, roccm1$p.value
  )
  MultComp_table[2, ] <- cbind(
    roccm2$estimate[1],
    roccm2$estimate[2], abs(roccm2$estimate[1] - roccm2$estimate[2]),
    abs(roccm2$estimate[1] - roccm2$estimate[2]) / roccm2$statistic,
    roccm2$statistic, roccm2$p.value
  )
  MultComp_table[3, ] <- cbind(
    rocm1m2$estimate[1],
    rocm1m2$estimate[2], abs(rocm1m2$estimate[1] - rocm1m2$estimate[2]),
    abs(rocm1m2$estimate[1] - rocm1m2$estimate[2]) / rocm1m2$statistic,
    rocm1m2$statistic, rocm1m2$p.value
  )

  comp.names <- cbind(
    c("Combination", "Combination", colnames(markers)[1]),
    c(colnames(markers)[1], colnames(markers)[2], colnames(markers)[2])
  )


  MultComp_table <- data.frame(comp.names, MultComp_table)
  colnames(MultComp_table) <- c(
    "Marker1 (A)", "Marker2 (B)", "AUC (A)", "AUC (B)",
    "|A-B|", "SE(|A-B|)", "z", "p-value"
  )
  data <- data.frame(markers, status, comb.score)

  cutoff.m1 <- OptimalCutpoints::optimal.cutpoints(colnames(data[1]), colnames(data[3]),
    tag.healthy = min(status),
    methods = cutoff.method, data = data
  )
  cutoff.m2 <- OptimalCutpoints::optimal.cutpoints(colnames(data[2]), colnames(data[3]),
    tag.healthy = min(status),
    methods = cutoff.method, data = data
  )
  cutoff.mc <- OptimalCutpoints::optimal.cutpoints(colnames(data[4]), colnames(data[3]),
    tag.healthy = min(status),
    methods = cutoff.method, data = data
  )

  threshold.m1 <- cutoff.m1[[cutoff.method]]$Global$optimal.cutoff
  TP.m1 <- cutoff.m1[[cutoff.method]]$Global$measures.acc$n$d - threshold.m1$FN
  TN.m1 <- cutoff.m1[[cutoff.method]]$Global$measures.acc$n$h - threshold.m1$FP

  threshold.m2 <- cutoff.m2[[cutoff.method]]$Global$optimal.cutoff
  TP.m2 <- cutoff.m2[[cutoff.method]]$Global$measures.acc$n$d - threshold.m2$FN
  TN.m2 <- cutoff.m2[[cutoff.method]]$Global$measures.acc$n$h - threshold.m2$FP

  threshold.mc <- cutoff.mc[[cutoff.method]]$Global$optimal.cutoff
  TP.mc <- cutoff.mc[[cutoff.method]]$Global$measures.acc$n$d - threshold.mc$FN
  TN.mc <- cutoff.mc[[cutoff.method]]$Global$measures.acc$n$h - threshold.mc$FP

  best.m1.tbl <- as.table(matrix(c(TP.m1[1], threshold.m1$FP[1], threshold.m1$FN[1], TN.m1[1]),
    nrow = 2, byrow = TRUE
  ))
  best.m2.tbl <- as.table(matrix(c(TP.m2[1], threshold.m2$FP[1], threshold.m2$FN[1], TN.m2[1]),
    nrow = 2, byrow = TRUE
  ))
  best.c.tbl <- as.table(matrix(c(TP.mc[1], threshold.mc$FP[1], threshold.mc$FN[1], TN.mc[1]),
    nrow = 2, byrow = TRUE
  ))

  DiagStatMarker1 <- epiR::epi.tests(best.m1.tbl, conf.level = conf.level)
  DiagStatMarker2 <- epiR::epi.tests(best.m2.tbl, conf.level = conf.level)
  DiagStatCombined <- epiR::epi.tests(best.c.tbl, conf.level = conf.level)

  allres <- list(
    ROC_coordinates = ROC_coordinates,
    AUC_table = AUC_table,
    MultComp_table = MultComp_table,
    DiagStatMarker1 = DiagStatMarker1,
    DiagStatMarker2 = DiagStatMarker2,
    DiagStatCombined = DiagStatCombined,
    ThresholdMarker1 = threshold.m1$cutoff,
    ThresholdMarker2 = threshold.m2$cutoff,
    ThresholdCombined = threshold.mc$cutoff,
    Criterion.m1 = cutoff.m1[[cutoff.method]]$Global$optimal.criterion,
    Criterion.m2 = cutoff.m2[[cutoff.method]]$Global$optimal.criterion,
    Criterion.c = cutoff.mc[[cutoff.method]]$Global$optimal.criterion,
    CombScore = comb.score,
    Cuttoff.method = cutoff.method
  )

  class(allres) <- "dtComb"
  graphics::par(mfrow = c(2, 2))
  return(allres)
}
