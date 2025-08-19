## Step 3. Graphical representation functions ----
### a. Distribution moments representation (mean, variance, skewness and kurtosis) ----

#' @title Graph: Distribution moments
#' @param MOM Moments data frame (mean, variance, skewness, kurtosis)
#' @param SESMOM SES of the Moments data frame and significance compared to null model
#' @param statisticsFactorName column of data use for colors discrimination
#' @param statisticsFactorNameBreaks vector of factor levels of the statisticsFactorName, same dimension than statisticsFactorNameCol
#' @param statisticsFactorNameCol vector of colors, same dimension than statisticsFactorNameBreaks
#' @param saveGraphMoments The path to save the graph
#' @export

GraphMoments <- function(MOM,
                         SESMOM,
                         statisticsFactorName,
                         statisticsFactorNameBreaks = NULL,
                         statisticsFactorNameCol = palette(),
                         saveGraphMoments) {
  ggplot2::ggsave(saveGraphMoments,
                  ggpubr::ggarrange(
                    ggplot2::ggplot()+
                      ggplot2::geom_boxplot(data = MOM %>%
                                              dplyr::filter(Number > 0),
                                            ggplot2::aes(x = "Mean", y = mean),
                                            col = "black", fill = "lightgrey", alpha = 0.4)+
                      ggplot2::geom_point(data = MOM %>%
                                            dplyr::filter(Number == 0),
                                          ggplot2::aes(x = "Mean", y = mean, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                                          shape = 21, size = 4, alpha = 0.4, position = "jitter")+
                      ggplot2::scale_fill_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::scale_color_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::theme_bw()+
                      ggplot2::labs(y = "Moments")+
                      ggplot2::theme(plot.title = ggplot2::element_blank(),
                                     axis.text.y = ggplot2::element_text(size = 20),
                                     axis.title.y = ggplot2::element_text(size = 30),
                                     axis.text.x = ggplot2::element_blank(),
                                     axis.title.x = ggplot2::element_blank(),
                                     legend.title = ggplot2::element_text(size = 30, face = "bold"),
                                     legend.text = ggplot2::element_text(size = 30),
                                     legend.key.size = ggplot2::unit(1.5, 'cm')),
                    ggplot2::ggplot()+
                      ggplot2::geom_boxplot(data = MOM %>%
                                              dplyr::filter(Number > 0),
                                            ggplot2::aes(x = "Variance", y = variance),
                                            col = "black", fill = "lightgrey", alpha = 0.4)+
                      ggplot2::geom_point(data = MOM %>%
                                            dplyr::filter(Number == 0),
                                          ggplot2::aes(x = "Variance", y = variance, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                                          shape = 21, size = 4, alpha = 0.4, position = "jitter")+
                      ggplot2::scale_fill_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::scale_color_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::theme_bw()+
                      ggplot2::labs(y = "Moments")+
                      ggplot2::theme(plot.title = ggplot2::element_blank(),
                                     axis.text.y = ggplot2::element_text(size = 20),
                                     axis.title.y = ggplot2::element_blank(),
                                     axis.text.x = ggplot2::element_blank(),
                                     axis.title.x = ggplot2::element_blank(),
                                     legend.title = ggplot2::element_text(size = 30, face = "bold"),
                                     legend.text = ggplot2::element_text(size = 30),
                                     legend.key.size = ggplot2::unit(1.5, 'cm')),
                    ggplot2::ggplot()+
                      ggplot2::geom_boxplot(data = MOM %>%
                                              dplyr::filter(Number > 0),
                                            ggplot2::aes(x = "Skewness", y = skewness),
                                            col = "black", fill = "lightgrey", alpha = 0.4)+
                      ggplot2::geom_point(data = MOM %>%
                                            dplyr::filter(Number == 0),
                                          ggplot2::aes(x = "Skewness", y = skewness, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                                          shape = 21, size = 4, alpha = 0.4, position = "jitter")+
                      ggplot2::scale_fill_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::scale_color_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::theme_bw()+
                      ggplot2::labs(y = "Moments")+
                      ggplot2::theme(plot.title = ggplot2::element_blank(),
                                     axis.text.y = ggplot2::element_text(size = 20),
                                     axis.title.y = ggplot2::element_blank(),
                                     axis.text.x = ggplot2::element_blank(),
                                     axis.title.x = ggplot2::element_blank(),
                                     legend.title = ggplot2::element_text(size = 30, face = "bold"),
                                     legend.text = ggplot2::element_text(size = 30),
                                     legend.key.size = ggplot2::unit(1.5, 'cm')),
                    ggplot2::ggplot()+
                      ggplot2::geom_boxplot(data = MOM %>%
                                              dplyr::filter(Number > 0),
                                            ggplot2::aes(x = "Kurtosis", y = kurtosis),
                                            col = "black", fill = "lightgrey", alpha = 0.4)+
                      ggplot2::geom_point(data = MOM %>%
                                            dplyr::filter(Number == 0),
                                          ggplot2::aes(x = "Kurtosis", y = kurtosis, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                                          shape = 21, size = 4, alpha = 0.4, position = "jitter")+
                      ggplot2::scale_fill_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::scale_color_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::theme_bw()+
                      ggplot2::labs(y = "Moments")+
                      ggplot2::theme(plot.title = ggplot2::element_blank(),
                                     axis.text.y = ggplot2::element_text(size = 20),
                                     axis.title.y = ggplot2::element_blank(),
                                     axis.text.x = ggplot2::element_blank(),
                                     axis.title.x = ggplot2::element_blank(),
                                     legend.title = ggplot2::element_text(size = 30, face = "bold"),
                                     legend.text = ggplot2::element_text(size = 30),
                                     legend.key.size = ggplot2::unit(1.5, 'cm')),
                    ggplot2::ggplot()+
                      ggplot2::geom_abline(intercept = 0, slope = 0, color = "grey", linewidth = 1, linetype = "dashed")+
                      ggplot2::geom_point(data = SESMOM %>%
                                            dplyr::filter(significanceMean == "TRUE"),
                                          ggplot2::aes(x = "Mean", y = standardizedObservedMean, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                                          shape = 21, size = 4, alpha = 0.8, position = "jitter")+
                      ggplot2::geom_point(data = SESMOM %>%
                                            dplyr::filter(significanceMean == "FALSE"),
                                          ggplot2::aes(x = "Mean", y = standardizedObservedMean, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                                          shape = 21, size = 4, alpha = 0.2, position = "jitter")+
                      ggplot2::scale_fill_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::scale_color_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::theme_bw()+
                      ggplot2::labs(y = paste0("SES ", "Moments"))+
                      ggplot2::theme(plot.title = ggplot2::element_blank(),
                                     axis.text.y = ggplot2::element_text(size = 20),
                                     axis.title.y = ggplot2::element_text(size = 30),
                                     axis.text.x = ggplot2::element_text(size = 20),
                                     axis.title.x = ggplot2::element_blank(),
                                     legend.title = ggplot2::element_text(size = 30, face = "bold"),
                                     legend.text = ggplot2::element_text(size = 30),
                                     legend.key.size = ggplot2::unit(1.5, 'cm')),
                    ggplot2::ggplot()+
                      ggplot2::geom_abline(intercept = 0, slope = 0, color = "grey", linewidth = 1, linetype = "dashed")+
                      ggplot2::geom_point(data = SESMOM %>%
                                            dplyr::filter(significanceVariance == "TRUE"),
                                          ggplot2::aes(x = "Variance", y = standardizedObservedVariance, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                                          shape = 21, size = 4, alpha = 0.8, position = "jitter")+
                      ggplot2::geom_point(data = SESMOM %>%
                                            dplyr::filter(significanceVariance == "FALSE"),
                                          ggplot2::aes(x = "Variance", y = standardizedObservedVariance, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                                          shape = 21, size = 4, alpha = 0.2, position = "jitter")+
                      ggplot2::scale_fill_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::scale_color_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::theme_bw()+
                      ggplot2::labs()+
                      ggplot2::theme(plot.title = ggplot2::element_blank(),
                                     axis.text.y = ggplot2::element_text(size = 20),
                                     axis.title.y = ggplot2::element_blank(),
                                     axis.text.x = ggplot2::element_text(size = 20),
                                     axis.title.x = ggplot2::element_blank(),
                                     legend.title = ggplot2::element_text(size = 30, face = "bold"),
                                     legend.text = ggplot2::element_text(size = 30),
                                     legend.key.size = ggplot2::unit(1.5, 'cm')),
                    ggplot2::ggplot()+
                      ggplot2::geom_abline(intercept = 0, slope = 0, color = "grey", linewidth = 1, linetype = "dashed")+
                      ggplot2::geom_point(data = SESMOM %>%
                                            dplyr::filter(significanceSkewness == "TRUE"),
                                          ggplot2::aes(x = "Skewness", y = standardizedObservedSkewness, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                                          shape = 21, size = 4, alpha = 0.8, position = "jitter")+
                      ggplot2::geom_point(data = SESMOM %>%
                                            dplyr::filter(significanceSkewness == "FALSE"),
                                          ggplot2::aes(x = "Skewness", y = standardizedObservedSkewness, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                                          shape = 21, size = 4, alpha = 0.2, position = "jitter")+
                      ggplot2::scale_fill_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::scale_color_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::theme_bw()+
                      ggplot2::labs()+
                      ggplot2::theme(plot.title = ggplot2::element_blank(),
                                     axis.text.y = ggplot2::element_text(size = 20),
                                     axis.title.y = ggplot2::element_blank(),
                                     axis.text.x = ggplot2::element_text(size = 20),
                                     axis.title.x = ggplot2::element_blank(),
                                     legend.title = ggplot2::element_text(size = 30, face = "bold"),
                                     legend.text = ggplot2::element_text(size = 30),
                                     legend.key.size = ggplot2::unit(1.5, 'cm')),
                    ggplot2::ggplot()+
                      ggplot2::geom_abline(intercept = 0, slope = 0, color = "grey", linewidth = 1, linetype = "dashed")+
                      ggplot2::geom_point(data = SESMOM %>%
                                            dplyr::filter(significanceKurtosis == "TRUE"),
                                          ggplot2::aes(x = "Kurtosis", y = standardizedObservedKurtosis, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                                          shape = 21, size = 4, alpha = 0.8, position = "jitter")+
                      ggplot2::geom_point(data = SESMOM %>%
                                            dplyr::filter(significanceKurtosis == "FALSE"),
                                          ggplot2::aes(x = "Kurtosis", y = standardizedObservedKurtosis, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)),
                                          shape = 21, size = 4, alpha = 0.2, position = "jitter")+
                      ggplot2::scale_fill_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::scale_color_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::theme_bw()+
                      ggplot2::labs()+
                      ggplot2::theme(plot.title = ggplot2::element_blank(),
                                     axis.text.y = ggplot2::element_text(size = 20),
                                     axis.title.y = ggplot2::element_blank(),
                                     axis.text.x = ggplot2::element_text(size = 20),
                                     axis.title.x = ggplot2::element_blank(),
                                     legend.title = ggplot2::element_text(size = 30, face = "bold"),
                                     legend.text = ggplot2::element_text(size = 30),
                                     legend.key.size = ggplot2::unit(1.5, 'cm')),
                    ncol = 4,
                    nrow = 2,
                    common.legend = T,
                    legend = "bottom"
                  ),
                  dpi = 600,
                  width = 15,
                  height = 8
  )
}
