## Step 3. Graphical representation functions ----
### c. Parameters of the SKR ----

#' @title Graph: parameters of the SKR
#' @description Generate and save plot with the SES SKR parameters values (observations compared to randomizations)
#' @param SESSKRparam SES of SKR parameters data frame (SES and significance)
#' @param statisticsFactorName column of data use for colors discrimination
#' @param statisticsFactorNameBreaks vector of factor levels of the statisticsFactorName, same dimension than statisticsFactorNameCol
#' @param statisticsFactorNameCol vector of colors, same dimension than statisticsFactorNameBreaks
#' @param path_GraphparamSKR The path to save the graph
#' @param slope_reference_TADs slope of a specific SKR used as a baseline (default: slope_reference_TADs = 1; skew-uniform slope)
#' @param intercept_reference_TADs intercept of a specific SKR used as a baseline (default: intercept_reference_TADs = 1.86; skew-uniform intercept)
#' @returns Plot of the SES SKR parameters values and significance
#' @export
#' @examples
#'
#' Example of how to use the function for grassland plant communities
#' under contrasting management practices.
#'
#' SKR.TAD::GraphparamSKR(
#' SESSKRparam = readRDS("./SES_SKRDataFrame.RDS"),
#' statisticsFactorName = c("Treatment"),
#' statisticsFactorNameBreaks = c("Mown_Unfertilized", "Mown_NPK"),
#' statisticsFactorNameCol = c("#1A85FF", "#D41159"),
#' slope_reference_TADs = 1,
#' intercept_reference_TADs = 1.86,
#' path_GraphparamSKR = "./paramSKR.png"
#' )

GraphparamSKR <- function(SESSKRparam,
                          statisticsFactorName,
                          statisticsFactorNameBreaks = NULL,
                          statisticsFactorNameCol = palette(),
                          slope_reference_TADs = 1,
                          intercept_reference_TADs = 1.86,
                          path_GraphparamSKR) {

  title_dist_reference_TADs <- paste0("Distance from reference TADs:
K = ", slope_reference_TADs, " x S² + ", intercept_reference_TADs)

  ggplot2::ggsave(
    path_GraphparamSKR,
    ggplot2::ggplot()+
      ggplot2::geom_abline(intercept = 0, slope = 0, color = "grey", linewidth = 1, linetype = "dashed")+
      ggplot2::geom_point(data = SESSKRparam |>
                   dplyr::filter(significanceSlope == TRUE),
                   ggplot2::aes(x =  "Slope", y = SESSlope, fill = !!rlang::sym(statisticsFactorName)),
                 alpha = 0.8, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SESSKRparam |>
                   dplyr::filter(significanceSlope == FALSE),
                   ggplot2::aes(x =  "Slope", y = SESSlope, fill = !!rlang::sym(statisticsFactorName)),
                 alpha = 0.2, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SESSKRparam |>
                   dplyr::filter(significanceIntercept == TRUE),
                   ggplot2::aes(x =  "Intercept", y = SESIntercept, fill = !!rlang::sym(statisticsFactorName)),
                 alpha = 0.8, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SESSKRparam |>
                   dplyr::filter(significanceIntercept == FALSE),
                   ggplot2::aes(x =  "Intercept", y = SESIntercept, fill = !!rlang::sym(statisticsFactorName)),
                 alpha = 0.2, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SESSKRparam |>
                   dplyr::filter(significanceRsquare == TRUE),
                   ggplot2::aes(x =  "R²", y = SESRsquare, fill = !!rlang::sym(statisticsFactorName)),
                 alpha = 0.8, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SESSKRparam |>
                   dplyr::filter(significanceRsquare == FALSE),
                   ggplot2::aes(x =  "R²", y = SESRsquare, fill = !!rlang::sym(statisticsFactorName)),
                 alpha = 0.2, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SESSKRparam |>
                   dplyr::filter(significancedistance_reference_TADs == TRUE),
                   ggplot2::aes(x = title_dist_reference_TADs, y = SESdistance_reference_TADs, fill = !!rlang::sym(statisticsFactorName)),
                 alpha = 0.8, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SESSKRparam |>
                   dplyr::filter(significancedistance_reference_TADs == FALSE),
                   ggplot2::aes(x = title_dist_reference_TADs, y = SESdistance_reference_TADs, fill = !!rlang::sym(statisticsFactorName)),
                 alpha = 0.2, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SESSKRparam |>
                   dplyr::filter(significanceCV_distance_reference_TADs == TRUE),
                   ggplot2::aes(x = paste0("CV ", title_dist_reference_TADs), y = SESCV_distance_reference_TADs, fill = !!rlang::sym(statisticsFactorName)),
                 alpha = 0.8, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SESSKRparam |>
                   dplyr::filter(significanceCV_distance_reference_TADs == FALSE),
                   ggplot2::aes(x = paste0("CV ", title_dist_reference_TADs), y = SESCV_distance_reference_TADs, fill = !!rlang::sym(statisticsFactorName)),
                 alpha = 0.2, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SESSKRparam |>
                   dplyr::filter(significancedistance_predicted_TADs == TRUE),
                   ggplot2::aes(x =  "Distance from predicted TADs", y = SESdistance_predicted_TADs, fill = !!rlang::sym(statisticsFactorName)),
                 alpha = 0.8, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SESSKRparam |>
                   dplyr::filter(significancedistance_predicted_TADs == FALSE),
                   ggplot2::aes(x =  "Distance from predicted TADs", y = SESdistance_predicted_TADs, fill = !!rlang::sym(statisticsFactorName)),
                 alpha = 0.2, size = 6, color = "black", shape = 21)+
      ggplot2::scale_x_discrete(limits = c("Slope", "Intercept", "R²", "Distance from predicted TADs", title_dist_reference_TADs, paste0("CV ", title_dist_reference_TADs)))+
      ggplot2::scale_fill_manual(limits = statisticsFactorNameBreaks, values = statisticsFactorNameCol)+
      ggplot2::scale_color_manual(limits = statisticsFactorNameBreaks, values = statisticsFactorNameCol)+
      ggplot2::theme_bw()+
      ggplot2::labs(title = paste0("Parameters of the SKR"), y = "SES")+
      ggplot2::theme(legend.position = "bottom",
            plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
            axis.text.y = ggplot2::element_text(size = 10),
            axis.title.y = ggplot2::element_text(size = 12, face = "bold"),
            axis.title.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(size = 10, face = "bold")),
    dpi = 600,
    width = 15,
    height = 5
  )
}
