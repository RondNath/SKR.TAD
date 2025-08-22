## Step 3. Graphical representation functions ----
### b. SKR representation ----

#' @title Graph: SKR
#' @description Generate and save plot with the Skewness-Kurtosis Relationship (SKR: Kurtosis ~ slope x Skewness² + intercept)
#' @param MOM moments data frame (mean, variance, skewness, kurtosis)
#' @param statisticsFactorName column of data use for colors discrimination
#' @param statisticsFactorNameBreaks vector of factor levels of the statisticsFactorName, same dimension than statisticsFactorNameCol
#' @param statisticsFactorNameCol vector of colors, same dimension than statisticsFactorNameBreaks
#' @param slope_speTADs slope of the theoretical distribution law (default: slope = 1 intercept = 1.86 skew-uniform)
#' @param intercept_speTADs intercept of the theoretical distribution law (default: slope = 1 intercept = 1.86 skew-uniform)
#' @param saveGraphSKR The path to save the graph
#' @returns Plot of the SKR and parameters for observed and randomized communities
#' @export
#' @examples
#'
#' Example of how to use the function for grassland plant communities
#' under contrasting management practices.
#'
#' SKR.TAD::GraphSKR(
#' MOM = readRDS("./MomentsDataFrame.RDS"),
#' statisticsFactorName = c("Treatment"),
#' statisticsFactorNameBreaks = c("Mown_Unfertilized", "Mown_NPK"),
#' statisticsFactorNameCol = c("#1A85FF", "#D41159"),
#' slope_speTADs = 1,
#' intercept_speTADs = 1.86,
#' saveGraphSKR = "./SKR.png"
#' )

GraphSKR <- function(
    MOM,
    statisticsFactorName,
    statisticsFactorNameBreaks = NULL,
    statisticsFactorNameCol = palette(),
    slope_speTADs = 1,
    intercept_speTADs = 1.86,
    saveGraphSKR
) {
  ggplot2::ggsave(saveGraphSKR,
                  ggplot2::ggplot() +
                    ggplot2::geom_point(data = MOM |>
                                          dplyr::filter(Number > 0),
                                        ggplot2::aes(x = skewness**2, y = kurtosis),
                                        shape = 21, size = 2, alpha = 0.4, col = "#D3D3D3", fill = "#D3D3D3")+
                    ggplot2::geom_smooth(data = MOM |>
                                           dplyr::filter(Number > 0),
                                         ggplot2::aes(x = skewness**2, y = kurtosis, group = Number),
                                         col = "#D3D3D3", fill = "#D3D3D3", se = F, method = "lm", formula = y ~ x, linetype = 1, linewidth = 0.5, alpha = 0.1)+
                    ggplot2::geom_abline(intercept = intercept_speTADs, slope = slope_speTADs, linetype = "dashed", linewidth = 2) +
                    ggplot2::geom_point(data = MOM |>
                                          dplyr::filter(Number == 0),
                                        ggplot2::aes(x = skewness**2, y = kurtosis, fill = !!rlang::sym(statisticsFactorName)),
                                        shape = 21, size = 6, alpha = 0.4)+
                    ggplot2::geom_smooth(data = MOM |>
                                           dplyr::filter(Number == 0),
                                         ggplot2::aes(x = skewness**2, y = kurtosis, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                                         se = F, method = "lm", formula = y ~ x, linetype = 1, linewidth = 2, alpha = 0.1)+
                    ggpubr::stat_regline_equation(data = MOM |>
                                                    dplyr::filter(Number == 0),
                                          ggplot2::aes(skewness**2, y = kurtosis, label =  paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~"), col = !!rlang::sym(statisticsFactorName)),
                                          alpha = 1, size = 8)+
                    ggplot2::scale_fill_manual(limits = statisticsFactorNameBreaks, values = statisticsFactorNameCol)+
                    ggplot2::scale_color_manual(limits = statisticsFactorNameBreaks, values = statisticsFactorNameCol)+
                    ggplot2::xlim(0, 10)+
                    ggplot2::ylim(0, 20)+
                    ggplot2::theme_bw()+
                    ggplot2::labs(x = "Skewness²", y = "Kurtosis")+
                    ggplot2::theme(legend.position = "bottom",
                                   plot.title = ggplot2::element_blank(),
                                   axis.text.y = ggplot2::element_text(size = 35),
                                   axis.title.y = ggplot2::element_text(size = 40),
                                   axis.title.x = ggplot2::element_text(size = 40),
                                   axis.text.x = ggplot2::element_text(size = 35)),
                  dpi = 600,
                  height = 10,
                  width = 10
  )
}
