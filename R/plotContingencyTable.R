#' @title plots the contingency table of two clustering results
#'
#' @author ranjanb
#'
#' @param cluster_labels_1 First vector of cluster labels for each cell
#' @param cluster_labels_2 Second vector of cluster labels for each cell
#' @param automateConsensus Boolean indicating whether automated consensus should be returned
#' @param minClustSize Minimum number of cells in a cluster. Default = 10.
#' @param filename name of contingency table file
#'
#' @return consensus cluster labels vector
#'
#' @export
#'
plotContingencyTable <- function(cluster_labels_1 = NULL, cluster_labels_2 = NULL, automateConsensus = T, minClustSize = 10, filename = "Contingency_Table.pdf") {

    if(is.null(cluster_labels_1) | is.null(cluster_labels_2)) {
        stop("Incomplete parameters provided.")
    }

    ctg_table <- table(cluster_labels_1, cluster_labels_2)
    ctg_df <- as.data.frame(ctg_table)
    ctg_df <- reshape2::dcast(data = ctg_df, formula = cluster_labels_1~cluster_labels_2, fun.aggregate = sum, value.var = "Freq")
    rownames(ctg_df) <- ctg_df$cluster_labels_1
    ctg_df$cluster_labels_1 <- NULL
    ctg_matrix <- as.matrix(ctg_df)


    pdf(filename, width=15, height=15)

    colorScheme <-
        circlize::colorRamp2(
            seq(-max(abs(
                ctg_matrix
            )), max(abs(
                ctg_matrix
            )), length.out = 5),
            c(
                "cyan",
                "#7FFF7F",
                "yellow",
                "#FF7F00",
                "red"
            )
        )

    ht <- ComplexHeatmap::Heatmap(matrix = ctg_matrix,
                                  col = colorScheme,
                                  column_title = "Cluster labels 2",
                                  row_title = "Cluster Labels 1",
                                  column_title_side = "top",
                                  column_title_gp = grid::gpar(fontsize = 30, fontface = "bold"),
                                  column_names_side = "top",
                                  column_names_gp = grid::gpar(fontsize = 20, fontface = "bold"),
                                  show_column_dend = FALSE,
                                  row_title_side = "left",
                                  row_title_gp = grid::gpar(fontsize = 30, fontface = "bold"),
                                  row_names_side = "left",
                                  row_names_gp = grid::gpar(fontsize = 20, fontface = "bold"),
                                  show_row_dend = FALSE,
                                  name = "Contingency Table",
                                  show_heatmap_legend = FALSE,
                                  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                                      grid::grid.text(ctg_matrix[i, j], x, y, gp = grid::gpar(fontsize = 20, fontface = "bold", col = "black"))
                                  })
    ComplexHeatmap::draw(ht)
    dev.off()

    if(automateConsensus) {
        if(length(unique(cluster_labels_1)) > length(unique(cluster_labels_2))){
            consensusClusterLabels = cluster_labels_1
            remainderLabels = cluster_labels_2
        } else if(length(unique(cluster_labels_1)) < length(unique(cluster_labels_2))) {
            consensusClusterLabels = cluster_labels_2
            remainderLabels = cluster_labels_1
        } else {
            if(median(table(cluster_labels_1)) > median(table(cluster_labels_2))) {
                consensusClusterLabels = cluster_labels_1
                remainderLabels = cluster_labels_2
            } else {
                consensusClusterLabels = cluster_labels_2
                remainderLabels = cluster_labels_1
            }
        }

        r = nrow(ctg_matrix)
        c = ncol(ctg_matrix)
        if (c > r) {
            ctg_matrix <- t(ctg_matrix)
            r = nrow(ctg_matrix)
            c = ncol(ctg_matrix)
        }
        else if(c == r) {
            if(length(intersect(unique(consensusClusterLabels), rownames(ctg_matrix)))!= r) {
                ctg_matrix <- t(ctg_matrix)
                r = nrow(ctg_matrix)
                c = ncol(ctg_matrix)
            }
        }


        for(i in 1:r) {

            row <- ctg_matrix[i, ]
            percent_row <- 100*(row/sum(row))

            for(j in 1:c) {
                if((percent_row[j] >= 10) && (row[j] > minClustSize)) {
                    consensusClusterLabels[which((consensusClusterLabels == rownames(ctg_matrix)[i]) & (remainderLabels == names(row)[j]))] <- paste(rownames(ctg_matrix)[i], names(row)[j], sep = "_")
                }
            }

        }
        return(consensusClusterLabels)
    }
}
