## Step 3. Graphical representation functions ----
### b. SKR representation ----

#' @title Graph: SKR
#' @description Generate and save plot with the Skewness-Kurtosis Relationship (SKR: Kurtosis ~ slope x Skewness² + intercept)
#' @param moments moments data frame (mean, variance, skewness, kurtosis)
#' @param statisticsFactorName column of data use for colors discrimination
#' @param statisticsFactorNameBreaks vector of factor levels of the statisticsFactorName, same dimension than statisticsFactorNameCol
#' @param statisticsFactorNameCol vector of colors, same dimension than statisticsFactorNameBreaks
#' @param slope_ref_TADs slope of the theoretical distribution law (default: slope = 1 intercept = 1.86 skew-uniform)
#' @param intercept_ref_TADs intercept of the theoretical distribution law (default: slope = 1 intercept = 1.86 skew-uniform)
#' @param path_GraphSKR The path to save the graph of the SKR
#' @returns Plot of the SKR for observed and randomized communities in two panels left) kurtosis ~ skewness, right) kurtosis ~ skewness²
#' @export
#' @examples
#' GraphSKR(moments = SKR.TAD::MomentsDataFrame,
#' statisticsFactorName = c("Treatment"),
#' statisticsFactorNameBreaks = c("Mown_Unfertilized", "Mown_NPK"),
#' statisticsFactorNameCol = c("#1A85FF", "#D41159"),
#' slope_ref_TADs = 1,
#' intercept_ref_TADs = 1.86,
#' path_GraphSKR = "./SKR.png")

GraphSKR <- function(
    moments,
    statisticsFactorName = NULL,
    statisticsFactorNameBreaks = NULL,
    statisticsFactorNameCol = grDevices::palette(),
    slope_ref_TADs = 1,
    intercept_ref_TADs = 1.86,
    path_GraphSKR = NULL
) {
  ggplot2::ggsave(path_GraphSKR,
                  ggpubr::ggarrange(
                    ggplot2::ggplot()+
                      ggplot2::stat_function(fun = \(x) intercept_ref_TADs + slope_ref_TADs*x^2,
                                             linewidth = 1,
                                             linetype = 2,
                                             col = "black") +
                      ggdensity::geom_hdr(data = moments|>
                                            dplyr::filter(Number > 0),
                                          ggplot2::aes(x = skewness,
                                                       y = kurtosis,
                                                       col = "grey",
                                                       fill = "lightgrey"),
                                          alpha = 0.5,
                                          probs = .9)+
                      ggplot2::geom_point(data = moments|>
                                            dplyr::filter(Number == 0),
                                          ggplot2::aes(x = skewness, y = kurtosis,
                                                       col = !!rlang::sym(statisticsFactorName),
                                                       fill = !!rlang::sym(statisticsFactorName)),
                                          alpha = 0.6,
                                          size = 5,
                                          shape = 21)+
                      ggplot2::scale_fill_manual(limits = statisticsFactorNameBreaks, values = statisticsFactorNameCol)+
                      ggplot2::scale_color_manual(limits = statisticsFactorNameBreaks, values = statisticsFactorNameCol)+
                      ggplot2::theme_bw() +
                      ggplot2::labs(x = "Skewness", y = "Kurtosis")+
                      ggplot2::theme(legend.position = "bottom",
                                     axis.text.y = ggplot2::element_text(size = 20),
                                     axis.title.y = ggplot2::element_text(size = 30),
                                     axis.text.x = ggplot2::element_text(size = 20),
                                     axis.title.x = ggplot2::element_text(size = 30)),
                    ggplot2::ggplot() +
                      ggplot2::geom_point(data = moments |>
                                            dplyr::filter(Number > 0),
                                          ggplot2::aes(x = skewness**2, y = kurtosis),
                                          shape = 21, size = 2, alpha = 0.4, col = "#D3D3D3", fill = "#D3D3D3")+
                      ggplot2::geom_smooth(data = moments |>
                                             dplyr::filter(Number > 0),
                                           ggplot2::aes(x = skewness**2, y = kurtosis, group = Number),
                                           col = "#D3D3D3", fill = "#D3D3D3", se = F, method = "lm", formula = y ~ x, linetype = 1, linewidth = 0.5, alpha = 0.1)+
                      ggplot2::geom_abline(intercept = intercept_ref_TADs, slope = slope_ref_TADs, linetype = "dashed", linewidth = 2) +
                      ggplot2::geom_point(data = moments |>
                                            dplyr::filter(Number == 0),
                                          ggplot2::aes(x = skewness**2, y = kurtosis, fill = !!rlang::sym(statisticsFactorName)),
                                          shape = 21, size = 6, alpha = 0.4)+
                      ggplot2::geom_smooth(data = moments |>
                                             dplyr::filter(Number == 0),
                                           ggplot2::aes(x = skewness**2, y = kurtosis, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                                           se = F, method = "lm", formula = y ~ x, linetype = 1, linewidth = 2, alpha = 0.1)+
                      ggpubr::stat_regline_equation(data = moments |>
                                                      dplyr::filter(Number == 0),
                                                    ggplot2::aes(skewness**2, y = kurtosis, label =  paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~"), col = !!rlang::sym(statisticsFactorName)),
                                                    alpha = 1, size = 8)+
                      ggplot2::scale_fill_manual(limits = statisticsFactorNameBreaks, values = statisticsFactorNameCol)+
                      ggplot2::scale_color_manual(limits = statisticsFactorNameBreaks, values = statisticsFactorNameCol)+
                      ggplot2::theme_bw()+
                      ggplot2::labs(x = "Skewness²", y = "Kurtosis")+
                      ggplot2::theme(legend.position = "none",
                                     plot.title = ggplot2::element_blank(),
                                     axis.text.y = ggplot2::element_text(size = 20),
                                     axis.title.y = ggplot2::element_text(size = 30),
                                     axis.text.x = ggplot2::element_text(size = 20),
                                     axis.title.x = ggplot2::element_text(size = 30)),

                    ncol = 2,
                    nrow = 1,
                    common.legend = T,
                    legend = "bottom"),
                  dpi = 600,
                  height = 10,
                  width = 15
  )
}
