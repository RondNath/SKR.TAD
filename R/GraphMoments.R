## Step 3. Graphical representation functions ----
### a. Distribution moments representation (mean, variance, skewness and kurtosis) ----

#' @title Graph: Distribution moments
#' @description Generate and save plot with moments (mean, variance, skewness & kurtosis) and SES moments values (observations compared to randomizations)
#' @param MomentsDataFrame Moments data frame (mean, variance, skewness, kurtosis)
#' @param SESMomentsDataFrame SES of the Moments data frame and significance compared to null model
#' @param statisticsFactorName column of data use for colors discrimination
#' @param statisticsFactorNameBreaks vector of factor levels of the statisticsFactorName, same dimension than statisticsFactorNameCol
#' @param statisticsFactorNameCol vector of colors, same dimension than statisticsFactorNameBreaks
#' @param path_GraphMoments The path to save the graph of the moments (mean, variance, skewness and kurtosis)
#' @returns Plot of the moments in two panels up) the moments, bottom) the SES moments values
#' @export
#' @examples
#' head(MomentsDataFrame)
#' head(SESMomentsDataFrame)
#' GraphMoments(
#' MomentsDataFrame = MomentsDataFrame,
#' SESMomentsDataFrame = SESMomentsDataFrame,
#' statisticsFactorName = c("Treatment"),
#' statisticsFactorNameBreaks = c("Mown_Unfertilized", "Mown_NPK"),
#' statisticsFactorNameCol = c("#1A85FF", "#D41159"),
#' path_GraphMoments = NULL
#' )

GraphMoments <- function(
    MomentsDataFrame,
    SESMomentsDataFrame,
    statisticsFactorName = NULL,
    statisticsFactorNameBreaks = NULL,
    statisticsFactorNameCol = grDevices::palette(),
    path_GraphMoments = NULL
) {
  GraphMoments <- ggpubr::ggarrange(
    ggplot2::ggplot()+
      ggplot2::geom_boxplot(data = MomentsDataFrame |>
                              dplyr::filter(Number > 0),
                            ggplot2::aes(x = "Mean", y = mean),
                            col = "black", fill = "lightgrey", alpha = 0.4)+
      ggplot2::geom_point(data = MomentsDataFrame |>
                            dplyr::filter(Number == 0),
                          ggplot2::aes(x = "Mean", y = mean, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                          shape = 21, size = 4, alpha = 0.4, position = "jitter")+
      ggplot2::scale_fill_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
      ggplot2::scale_color_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
      ggplot2::theme_bw()+
      ggplot2::labs(y = "Moments")+
      ggplot2::theme(plot.title = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(size = 10),
                     axis.title.y = ggplot2::element_text(size = 15),
                     axis.text.x = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank(),
                     legend.title = ggplot2::element_text(size = 15, face = "bold"),
                     legend.text = ggplot2::element_text(size = 15),
                     legend.key.size = ggplot2::unit(1, 'cm')),
    ggplot2::ggplot()+
      ggplot2::geom_boxplot(data = MomentsDataFrame |>
                              dplyr::filter(Number > 0),
                            ggplot2::aes(x = "Variance", y = variance),
                            col = "black", fill = "lightgrey", alpha = 0.4)+
      ggplot2::geom_point(data = MomentsDataFrame |>
                            dplyr::filter(Number == 0),
                          ggplot2::aes(x = "Variance", y = variance, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                          shape = 21, size = 4, alpha = 0.4, position = "jitter")+
      ggplot2::scale_fill_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
      ggplot2::scale_color_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
      ggplot2::theme_bw()+
      ggplot2::labs(y = "Moments")+
      ggplot2::theme(plot.title = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(size = 10),
                     axis.title.y = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank(),
                     legend.title = ggplot2::element_text(size = 15, face = "bold"),
                     legend.text = ggplot2::element_text(size = 15),
                     legend.key.size = ggplot2::unit(1, 'cm')),
    ggplot2::ggplot()+
      ggplot2::geom_boxplot(data = MomentsDataFrame |>
                              dplyr::filter(Number > 0),
                            ggplot2::aes(x = "Skewness", y = skewness),
                            col = "black", fill = "lightgrey", alpha = 0.4)+
      ggplot2::geom_point(data = MomentsDataFrame |>
                            dplyr::filter(Number == 0),
                          ggplot2::aes(x = "Skewness", y = skewness, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                          shape = 21, size = 4, alpha = 0.4, position = "jitter")+
      ggplot2::scale_fill_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
      ggplot2::scale_color_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
      ggplot2::theme_bw()+
      ggplot2::labs(y = "Moments")+
      ggplot2::theme(plot.title = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(size = 10),
                     axis.title.y = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank(),
                     legend.title = ggplot2::element_text(size = 15, face = "bold"),
                     legend.text = ggplot2::element_text(size = 15),
                     legend.key.size = ggplot2::unit(1, 'cm')),
    ggplot2::ggplot()+
      ggplot2::geom_boxplot(data = MomentsDataFrame |>
                              dplyr::filter(Number > 0),
                            ggplot2::aes(x = "Kurtosis", y = kurtosis),
                            col = "black", fill = "lightgrey", alpha = 0.4)+
      ggplot2::geom_point(data = MomentsDataFrame |>
                            dplyr::filter(Number == 0),
                          ggplot2::aes(x = "Kurtosis", y = kurtosis, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                          shape = 21, size = 4, alpha = 0.4, position = "jitter")+
      ggplot2::scale_fill_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
      ggplot2::scale_color_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
      ggplot2::theme_bw()+
      ggplot2::labs(y = "Moments")+
      ggplot2::theme(plot.title = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(size = 10),
                     axis.title.y = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank(),
                     legend.title = ggplot2::element_text(size = 15, face = "bold"),
                     legend.text = ggplot2::element_text(size = 15),
                     legend.key.size = ggplot2::unit(1, 'cm')),
    ggplot2::ggplot()+
      ggplot2::geom_abline(intercept = 0, slope = 0, color = "grey", linewidth = 1, linetype = "dashed")+
      ggplot2::geom_point(data = SESMomentsDataFrame |>
                            dplyr::filter(significancemean == "TRUE"),
                          ggplot2::aes(x = "Mean", y = SESmean, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                          shape = 21, size = 4, alpha = 0.8, position = "jitter")+
      ggplot2::geom_point(data = SESMomentsDataFrame |>
                            dplyr::filter(significancemean == "FALSE"),
                          ggplot2::aes(x = "Mean", y = SESmean, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                          shape = 21, size = 4, alpha = 0.2, position = "jitter")+
      ggplot2::scale_fill_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
      ggplot2::scale_color_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
      ggplot2::theme_bw()+
      ggplot2::labs(y = paste0("SES ", "Moments"))+
      ggplot2::theme(plot.title = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(size = 10),
                     axis.title.y = ggplot2::element_text(size = 15),
                     axis.text.x = ggplot2::element_text(size = 10),
                     axis.title.x = ggplot2::element_blank(),
                     legend.title = ggplot2::element_text(size = 15, face = "bold"),
                     legend.text = ggplot2::element_text(size = 15),
                     legend.key.size = ggplot2::unit(1, 'cm')),
    ggplot2::ggplot()+
      ggplot2::geom_abline(intercept = 0, slope = 0, color = "grey", linewidth = 1, linetype = "dashed")+
      ggplot2::geom_point(data = SESMomentsDataFrame |>
                            dplyr::filter(significancevariance == "TRUE"),
                          ggplot2::aes(x = "Variance", y = SESvariance, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                          shape = 21, size = 4, alpha = 0.8, position = "jitter")+
      ggplot2::geom_point(data = SESMomentsDataFrame |>
                            dplyr::filter(significancevariance == "FALSE"),
                          ggplot2::aes(x = "Variance", y = SESvariance, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                          shape = 21, size = 4, alpha = 0.2, position = "jitter")+
      ggplot2::scale_fill_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
      ggplot2::scale_color_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
      ggplot2::theme_bw()+
      ggplot2::labs()+
      ggplot2::theme(plot.title = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(size = 10),
                     axis.title.y = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(size = 10),
                     axis.title.x = ggplot2::element_blank(),
                     legend.title = ggplot2::element_text(size = 15, face = "bold"),
                     legend.text = ggplot2::element_text(size = 15),
                     legend.key.size = ggplot2::unit(1, 'cm')),
    ggplot2::ggplot()+
      ggplot2::geom_abline(intercept = 0, slope = 0, color = "grey", linewidth = 1, linetype = "dashed")+
      ggplot2::geom_point(data = SESMomentsDataFrame |>
                            dplyr::filter(significanceskewness == "TRUE"),
                          ggplot2::aes(x = "Skewness", y = SESskewness, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                          shape = 21, size = 4, alpha = 0.8, position = "jitter")+
      ggplot2::geom_point(data = SESMomentsDataFrame |>
                            dplyr::filter(significanceskewness == "FALSE"),
                          ggplot2::aes(x = "Skewness", y = SESskewness, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                          shape = 21, size = 4, alpha = 0.2, position = "jitter")+
      ggplot2::scale_fill_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
      ggplot2::scale_color_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
      ggplot2::theme_bw()+
      ggplot2::labs()+
      ggplot2::theme(plot.title = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(size = 10),
                     axis.title.y = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(size = 10),
                     axis.title.x = ggplot2::element_blank(),
                     legend.title = ggplot2::element_text(size = 15, face = "bold"),
                     legend.text = ggplot2::element_text(size = 15),
                     legend.key.size = ggplot2::unit(1, 'cm')),
    ggplot2::ggplot()+
      ggplot2::geom_abline(intercept = 0, slope = 0, color = "grey", linewidth = 1, linetype = "dashed")+
      ggplot2::geom_point(data = SESMomentsDataFrame |>
                            dplyr::filter(significancekurtosis == "TRUE"),
                          ggplot2::aes(x = "Kurtosis", y = SESkurtosis, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                          shape = 21, size = 4, alpha = 0.8, position = "jitter")+
      ggplot2::geom_point(data = SESMomentsDataFrame |>
                            dplyr::filter(significancekurtosis == "FALSE"),
                          ggplot2::aes(x = "Kurtosis", y = SESkurtosis, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                          shape = 21, size = 4, alpha = 0.2, position = "jitter")+
      ggplot2::scale_fill_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
      ggplot2::scale_color_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
      ggplot2::theme_bw()+
      ggplot2::labs()+
      ggplot2::theme(plot.title = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(size = 10),
                     axis.title.y = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(size = 10),
                     axis.title.x = ggplot2::element_blank(),
                     legend.title = ggplot2::element_text(size = 15, face = "bold"),
                     legend.text = ggplot2::element_text(size = 15),
                     legend.key.size = ggplot2::unit(1, 'cm')),
    ncol = 4,
    nrow = 2,
    common.legend = T,
    legend = "bottom"
  )

  if (!is.null(path_GraphMoments)) {
    message("Save Moments plot", path_GraphMoments)
    ggplot2::ggsave(
      path_GraphMoments,
      GraphMoments,
      dpi = 600,
      width = 15,
      height = 8
    )
  }

  GraphMoments
}
