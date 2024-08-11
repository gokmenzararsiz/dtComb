#' @title Plot the combination scores using the training model
#'
#' @description The \code{plotComb} a function that generates plots from the
#' training model. The function takes argument model. The outputs of the
#' function are three different plots generated from the combination scores.
#'
#' @param model a \code{list} object where the parameters from the training
#' model are saved.
#'
#' @param status a \code{factor} vector that includes the actual disease
#' status of the patients
#'
#' @return A \code{data.frame} plots
#'
#' @author Serra Ilayda Yerlitas, Serra Bersan Gengec, Necla Kochan,
#' Gozde Erturk Zararsiz, Selcuk Korkmaz, Gokmen Zararsiz
#'
#' @examples
#'
#' # call data
#' data(exampleData1)
#'
#' # define the function parameters
#' markers <- exampleData1[, -1]
#' status <- factor(exampleData1$group, levels = c("not_needed", "needed"))
#' event <- "needed"
#'
#' score1 <- linComb(
#'   markers = markers, status = status, event = event,
#'   method = "scoring", resample = "none",
#'   standardize = "none", direction = "<", cutoff.method = "Youden"
#' )
#'
#' plotComb(score1, status)
#'
#' score2 <- nonlinComb(
#'   markers = markers, status = status, event = event,
#'   method = "nsgam", resample = "cv", include.interact = FALSE, direction = "<",
#'   standardize = "zScore", cutoff.method = "Youden"
#' )
#'
#' plot.score2 <- plotComb(score2, status)
#'
#' score3 <- mathComb(
#'   markers = markers, status = status, event = event,
#'   method = "distance", distance = "euclidean", direction = "auto",
#'   standardize = "tScore", cutoff.method = "Youden"
#' )
#'
#' plot.score3 <- plotComb(score3, status)
#'
#' @export

plotComb <- function(model, status) {
  # density graph

  combinedData <- cbind(as.data.frame(model$CombScore), as.data.frame(status))

  names(combinedData) <- c("CombinationScore", "Labels")

  plotDensity <- ggplot2::ggplot(combinedData, ggplot2::aes(
    x = CombinationScore,
    colour = Labels
  )) +
    ggplot2::xlab("Combination Score") +
    ggplot2::ylab("Density") +
    ggplot2::geom_density(size = 2) +
    ggplot2::ggtitle("Distribution Plot") +
    ggplot2::geom_vline(xintercept = model$ThresholdCombined, linetype = "dotted") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 22, face = "bold")) +
    ggplot2::theme(text = ggplot2::element_text(size = 20)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 20)) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 20)) +
    ggplot2::theme(legend.position = "bottom")


  # distribution scatter

  combinedData <- cbind(as.data.frame(model$CombScore), as.data.frame(status))

  names(combinedData) <- c("CombinationScore", "Labels")

  plotScatter <- ggplot2::ggplot(combinedData, ggplot2::aes(
    x = Labels, y = CombinationScore,
    color = Labels
  )) +
    ggplot2::ylab("Combination Score") +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_jitter(width = 0.40) +
    ggplot2::geom_point() +
    ggplot2::ggtitle("Scatter Plot") +
    ggplot2::geom_hline(yintercept = model$ThresholdCombined, linetype = "dotted") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 22, face = "bold")) +
    ggplot2::theme(text = ggplot2::element_text(size = 20)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 20)) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 20)) +
    ggplot2::theme(legend.position = "bottom")


  # sens & spec curve

  results <- model$ROC_coordinates

  coord <- results[results[, "Marker"] == "Combination", ]

  colors <- c("Sensitivity" = "#f8766d", "Specificity" = "#00bfc4")

  plotSensSpec <- ggplot2::ggplot(coord, ggplot2::aes(x = Threshold)) +
    ggplot2::geom_line(ggplot2::aes(y = Sensitivity, color = "Sensitivity"), show.legend = TRUE, size = 2) +
    ggplot2::geom_line(ggplot2::aes(y = Specificity, color = "Specificity"), show.legend = TRUE, size = 2) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 22, face = "bold")) +
    ggplot2::theme(text = ggplot2::element_text(size = 20)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 20)) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 20)) +
    ggplot2::ggtitle("Sensitivity&Specificity Plot") +
    ggplot2::labs(y = "Value", x = "Combination Score", color = "Labels") +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::geom_vline(xintercept = model$ThresholdCombined, linetype = "dotted") +
    ggplot2::theme(legend.position = "bottom")


  plotDensity <- ggpubr::ggarrange(plotDensity,
    ncol = 1, nrow = 1
  )
  plotScatter <- ggpubr::ggarrange(plotScatter,
    ncol = 1, nrow = 1
  )
  plotSensSpec <- ggpubr::ggarrange(plotSensSpec + ggpubr::rremove("x.text"),
    ncol = 1, nrow = 1
  )
  all <- ggpubr::ggarrange(plotDensity, plotScatter, plotSensSpec + ggpubr::rremove("x.text"),
    ncol = 3, nrow = 1
  )
  allPlots <- list(plotDensity = plotDensity, plotScatter = plotScatter, plotSensSpec = plotSensSpec, all = all)
  return(allPlots)
}
