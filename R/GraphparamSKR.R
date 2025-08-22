## Step 3. Graphical representation functions ----
### c. Parameters of the SKR ----

#' @title Graph: parameters of the SKR
#' @description Generate and save plot with the SES SKR parameters values (observations compared to randomizations)
#' @param data SES of SKR parameters data frame (SES and Significance)
#' @param statisticsFactorName column of data use for colors discrimination
#' @param statisticsFactorNameBreaks vector of factor levels of the statisticsFactorName, same dimension than statisticsFactorNameCol
#' @param statisticsFactorNameCol vector of colors, same dimension than statisticsFactorNameBreaks
#' @param saveGraphparamSKR The path to save the graph
#' @param slope_speTADs slope of a specific SKR used as a baseline (default: slope_speTADs = 1; skew-uniform slope)
#' @param intercept_speTADs intercept of a specific SKR used as a baseline (default: intercept_speTADs = 1.86; skew-uniform intercept)
#' @returns Plot of the SES SKR parameters values and significance
#' @export
#' @examples
#'
#' Example of how to use the function for grassland plant communities
#' under contrasting management practices.
#'
#' SKR.TAD::GraphparamSKR(
#' SKRparam = readRDS("./SES_SKRDataFrame.RDS"),
#' statisticsFactorName = c("Treatment"),
#' statisticsFactorNameBreaks = c("Mown_Unfertilized", "Mown_NPK"),
#' statisticsFactorNameCol = c("#1A85FF", "#D41159"),
#' slope_speTADs = 1,
#' intercept_speTADs = 1.86,
#' saveGraphparamSKR = "./paramSKR.png"
#' )

GraphparamSKR <- function(SKRparam,
                          statisticsFactorName,
                          statisticsFactorNameBreaks = NULL,
                          statisticsFactorNameCol = palette(),
                          slope_speTADs = 1,
                          intercept_speTADs = 1.86,
                          saveGraphparamSKR) {

  title_dist_speTADs <- paste0("Distance from specific TADs:
K = ", slope_speTADs, " x S² + ", intercept_speTADs)

  ggplot2::ggsave(
    saveGraphparamSKR,
    ggplot2::ggplot()+
      ggplot2::geom_abline(intercept = 0, slope = 0, color = "grey", linewidth = 1, linetype = "dashed")+
      ggplot2::geom_point(data = SKRparam |>
                   dplyr::filter(Slope_Signi == TRUE),
                   ggplot2::aes(x =  "Slope", y = Slope_SES, fill = !!rlang::sym(statisticsFactorName)),
                 alpha = 0.8, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SKRparam |>
                   dplyr::filter(Slope_Signi == FALSE),
                   ggplot2::aes(x =  "Slope", y = Slope_SES, fill = !!rlang::sym(statisticsFactorName)),
                 alpha = 0.2, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SKRparam |>
                   dplyr::filter(Intercept_Signi == TRUE),
                   ggplot2::aes(x =  "Intercept", y = Intercept_SES, fill = !!rlang::sym(statisticsFactorName)),
                 alpha = 0.8, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SKRparam |>
                   dplyr::filter(Intercept_Signi == FALSE),
                   ggplot2::aes(x =  "Intercept", y = Intercept_SES, fill = !!rlang::sym(statisticsFactorName)),
                 alpha = 0.2, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SKRparam |>
                   dplyr::filter(Rsquare_Signi == TRUE),
                   ggplot2::aes(x =  "R²", y = Rsquare_SES, fill = !!rlang::sym(statisticsFactorName)),
                 alpha = 0.8, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SKRparam |>
                   dplyr::filter(Rsquare_Signi == FALSE),
                   ggplot2::aes(x =  "R²", y = Rsquare_SES, fill = !!rlang::sym(statisticsFactorName)),
                 alpha = 0.2, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SKRparam |>
                   dplyr::filter(distance_specific_TADs_Signi == TRUE),
                   ggplot2::aes(x = title_dist_speTADs, y = distance_specific_TADs_SES, fill = !!rlang::sym(statisticsFactorName)),
                 alpha = 0.8, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SKRparam |>
                   dplyr::filter(distance_specific_TADs_Signi == FALSE),
                   ggplot2::aes(x = title_dist_speTADs, y = distance_specific_TADs_SES, fill = !!rlang::sym(statisticsFactorName)),
                 alpha = 0.2, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SKRparam |>
                   dplyr::filter(CV_distance_specific_TADs_Signi == TRUE),
                   ggplot2::aes(x = paste0("CV ", title_dist_speTADs), y = CV_distance_specific_TADs_SES, fill = !!rlang::sym(statisticsFactorName)),
                 alpha = 0.8, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SKRparam |>
                   dplyr::filter(CV_distance_specific_TADs_Signi == FALSE),
                   ggplot2::aes(x = paste0("CV ", title_dist_speTADs), y = CV_distance_specific_TADs_SES, fill = !!rlang::sym(statisticsFactorName)),
                 alpha = 0.2, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SKRparam |>
                   dplyr::filter(distance_predicted_TADs_Signi == TRUE),
                   ggplot2::aes(x =  "Distance from predicted TADs", y = distance_predicted_TADs_SES, fill = !!rlang::sym(statisticsFactorName)),
                 alpha = 0.8, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SKRparam |>
                   dplyr::filter(distance_predicted_TADs_Signi == FALSE),
                   ggplot2::aes(x =  "Distance from predicted TADs", y = distance_predicted_TADs_SES, fill = !!rlang::sym(statisticsFactorName)),
                 alpha = 0.2, size = 6, color = "black", shape = 21)+
      ggplot2::scale_x_discrete(limits = c("Slope", "Intercept", "R²", "Distance from predicted TADs", title_dist_speTADs, paste0("CV ", title_dist_speTADs)))+
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
