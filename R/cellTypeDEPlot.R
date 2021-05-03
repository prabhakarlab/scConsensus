#'@title Plots DE Heatmap showing how consensus cluster labels are refined by re-clustering over pair-wise DE genes
#'
#'@author ranjanb
#'
#'@param dataMatrix the log-transformed and normalized scRNAseq genes x cells matrix
#'@param nodg vector of NODG to visualize cell sequencing differences
#'@param cellTree tree obtained from hierarchical clustering
#'@param clusterLabels consensus cluster labels
#'@param dynamicColorsList labels of newly defined clusters on DE genes
#'@param geneLabels labels for the DE genes, here NULL.
#'@param colScheme c("green", "blue", "lightblue")
#'@param filename name of DE heatmap file
#'
#'@export
#'
#'
cellTypeDEPlot <- function(dataMatrix,
                           nodg = NULL,
                           cellTree,
                           clusterLabels,
                           dynamicColorsList,
                           geneLabels = NULL,
                           colScheme = "green",
                           filename = "DE_Heatmap") {

    # Initialize heatmap input
    heatmapIn <- dataMatrix


    ### Compute number of detected genes (NODG) for each cell
    if (is.null(nodg)) {
        nodg <- c()
        for (i in 1:length(colnames(heatmapIn))) {
            nodg[i] <- length(which(heatmapIn[, i] > 0))
        }
    }

    ### Create black and white annotation bars for each cluster

    # extract the unique cluster labels
    uniqueClusters <- unique(clusterLabels)

    # initialize empty dataframes to store cluster colors and labels
    allClusterColorDf <-
        data.frame(id = ncol(heatmapIn))
    allClusterDf <- data.frame(id = ncol(heatmapIn))

    # initialize empty vectors to store color and binary values for each cluster
    clusterColorVec <- c()
    clusterBinVec <- c()

    # for each cluster
    for (i in 1:length(uniqueClusters)) {
        # set cluster color as black where cluster label matches cluster i
        clusterColorVec[clusterLabels == uniqueClusters[i]] <-
            "black"

        # set cluster binary value 1 as black where cluster label matches cluster i
        clusterBinVec[clusterLabels == uniqueClusters[i]] <-
            "1"

        # set cluster color as white where cluster label matches cluster i
        clusterColorVec[clusterLabels != uniqueClusters[i]] <-
            "white"

        # set cluster binary value 1 as black where cluster label matches cluster i
        clusterBinVec[clusterLabels != uniqueClusters[i]] <-
            "0"

        # initialize data frame storing cluster colors
        clusterColorDf <- as.data.frame(clusterColorVec)

        # name data frame with cluster name
        names(clusterColorDf) <- uniqueClusters[i]

        # initialize data frame storing cluster binary values
        clusterDf <- as.data.frame(clusterBinVec)

        # name data frame with cluster name
        names(clusterDf) <- uniqueClusters[i]

        # append cluster color data frame to data frame of complete cluster set
        allClusterColorDf <-
            cbind(allClusterColorDf, clusterColorDf)

        # append cluster binary value data frame to data frame of complete cluster set
        allClusterDf <- cbind(allClusterDf, clusterDf)

        # reset cluster specific vectors and data frames
        clusterColorVec <- c()
        clusterColorDf <- NULL
        clusterBinVec <- c()
        clusterDf <- NULL

    }

    # reset id values to remove id column
    allClusterColorDf$id <- NULL
    allClusterDf$id <- NULL

    # save color list for annotation
    allClusterColors <- as.list(allClusterColorDf)

    ### Create annotation color list for cluster labels
    # initialise annotation color list with complete cluster color list
    colorList <- allClusterColors

    # append dynamic colors to annotation color list

    colorList <- append(colorList, dynamicColorsList)

    print(str(dynamicColorsList))

    # name each color element for HeatmapAnnotation

    for (i in 1:ncol(heatmapIn)) {
        for (j in 1:length(uniqueClusters)) {
            # if annotation color is black for cluster j, name the element 1
            if (colorList[[j]][i] == "black") {
                names(colorList[[j]])[i] <- "1"
            }
            # else name the element 0
            else {
                names(colorList[[j]])[i] <- "0"
            }
        }

        # name dynamic colors with color names
        k = 0
        for (ds in names(dynamicColorsList)) {
            names(colorList[[((k + 1) + length(uniqueClusters))]])[i] <-
                dynamicColorsList[[ds]][i]
            k = k + 1
        }

    }

    ### Create data frame of annotations
    # initialise annotation data frame with complete cluster set data frame
    annotationDf <- allClusterDf

    # append dynamic colors to annotation data frame

    dynamicColorsDf <- data.frame(dynamicColorsList)
    names(dynamicColorsDf) <- names(dynamicColorsList)
    annotationDf <-
        cbind(annotationDf, dynamicColorsDf)

    print(str(annotationDf))

    ### Plot Heatmap
    # initialize column color bar annotations along with NODG bar plot
    columnColorBar <-
        ComplexHeatmap::HeatmapAnnotation(df = annotationDf,
                          col = colorList,
                          NODG = anno_barplot(
                              nodg,
                              gp = grid::gpar(fill = "#777777", col = "#777777"),
                              axis = TRUE,
                              axis_param = list(side = "right")
                          ),
                          show_annotation_name = TRUE,
                          annotation_name_side = "left",
                          gap = grid::unit(5, "mm"),
                          which = "column"
        )

    # Z-transform rows for easy reading
    # heatmapIn <- t(apply(heatmapIn, 1, function(x){
    #     (x - mean(x))/sd(x)
    # }))

    # Select color scheme for heatmap ("blue",  "green" OR "lightblue")
    if (colScheme == "blue") {
        # "blue": color scheme varies from minimum value of data input to max value of data input
        # since most values are near the minimum value of input, majority of the heatmap will appear blue
        colorScheme <-
            circlize::colorRamp2(
                seq(min(heatmapIn), max(heatmapIn), length.out = 9),
                c(
                    "#00007F",
                    "blue",
                    "#007FFF",
                    "cyan",
                    "#7FFF7F",
                    "yellow",
                    "#FF7F00",
                    "red",
                    "#7F0000"
                )
            )
    } else if (colScheme == "green") {
        # "green": color scheme varies from negative of the maximum absolute value of data input to maximum absolute value of data input
        # since most values are near the middle of this range, majority of the heatmap will appear green
        colorScheme <-
            circlize::colorRamp2(
                seq(-max(abs(heatmapIn)), max(abs(heatmapIn)), length.out = 9),
                c(
                    "#00007F",
                    "blue",
                    "#007FFF",
                    "cyan",
                    "#7FFF7F",
                    "yellow",
                    "#FF7F00",
                    "red",
                    "#7F0000"
                )
            )
    } else if (colScheme == "violet") {
        # "violet": color scheme starts from light blue as the lowest value and goes through white to red and deep red
        # since most values are near the lower-mid range, the background will be violet
        colorScheme <-
            circlize::colorRamp2(
                seq(min(abs(heatmapIn)), max(abs(heatmapIn)), length.out = 5),
                c("#7777FF",
                  "white",
                  "red",
                  "#7F0000",
                  "#2F0000")
            )
    }

    # initialize heatmap object
    ht <- ComplexHeatmap::Heatmap(
        matrix = heatmapIn,
        col = colorScheme,

        cluster_columns = as.dendrogram(cellTree),
        cluster_rows = TRUE,

        column_dend_side = "top",
        row_dend_side = "left",

        show_row_dend = TRUE,
        row_dend_width = grid::unit(30, "mm"),

        show_column_dend = TRUE,
        column_dend_height = grid::unit(100, "mm"),
        column_dend_reorder = FALSE,

        show_column_names = TRUE,

        show_row_names = TRUE,
        row_names_gp = grid::gpar(fontsize = 5),

        top_annotation = columnColorBar,

        heatmap_legend_param = list(title = "legend", color_bar = "continuous"),
        use_raster = TRUE,
        raster_device = "png",
        raster_quality = 1
    )

    # create pdf object to hold heatmap
    pdf(paste0(filename, ".pdf"),
        width = 50,
        height = 50)

    if (!is.null(geneLabels)) {
        # Create row annotation for geneLabels
        geneColors <- labels2colors(geneLabels)
        names(geneColors) <- geneLabels

        row_list <- list(Gene_Groups = geneColors)

        row_hm_df <- data.frame(Gene_Groups = geneLabels)

        row_color_bar_anno <-
            ComplexHeatmap::HeatmapAnnotation(
                row_hm_df,
                col = row_list,
                show_annotation_name = TRUE,
                annotation_name_side = "top",
                which = "row"
            )

        # draw heatmap in pdf
        ComplexHeatmap::draw(ht + row_color_bar_anno,
             heatmap_legend_side = "left",
             annotation_legend_side = "left")

    } else {
        # draw heatmap in pdf
        ComplexHeatmap::draw(ht,
             heatmap_legend_side = "left",
             annotation_legend_side = "left")
    }


    # shut down device
    dev.off()
}
