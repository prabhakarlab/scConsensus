#' @title Recluster consensus clusters using DE gene analysis
#'
#' @author ranjanb, schmidtf
#'
#' @param dataMatrix the log-transformed and normalized scRNAseq genes x cells matrix
#' @param consensusClusterLabels consensus cell type labels for each cell
#' @param method Method used for DE gene statistical test
#' @param meanScalingFactor scale of the mean gene expression across the gene expression matrix to set a minimum threshold of average cluster expression for a gene to be considered a DE gene.
#' @param qValThrs maximum q-value threshold for statistical test
#' @param fcThrs minimum fold-change threshold for DE gene criterion (natural log following Seurat convention)
#' @param deepSplitValues vector of WGCNA tree cutting deepsplit parameters
#' @param minClusterSize specifies the type of minimum cluster size
#' @param filename name of DE gene object file
#' @param plotName name of DE Heatmap file
#'
#' @return list containing vector of DE genes, clustering tree and dynamic color list.
#'
#' @export
#'
reclusterDEConsensus <- function(dataMatrix,
                                 consensusClusterLabels,
                                 method = "Wilcoxon",
                                 meanScalingFactor = 5,
                                 qValThrs,
                                 fcThrs,
                                 deepSplitValues = 1:4,
                                 minClusterSize = 10,
                                 filename = "de_gene_object.rds",
                                 plotName = "DE_Heatmap") {

    ### Initialise data input
    dataIn <- as.matrix(dataMatrix)


    ### Initialize mean expression threshold
    meanExprsThrs = meanScalingFactor * mean(expm1(dataIn))

    ### Use only clusters with number of cells > minimum cluster size for DE gene calling
    which(table(consensusClusterLabels) > minClusterSize)
    colorCounts <-
        table(consensusClusterLabels)

    ### Extract unique cluster labels
    uniqueClusters <-
        names(colorCounts[colorCounts > minClusterSize])

    # ignore the grey cluster because it represents unclustered cells
    uniqueClusters <-
        uniqueClusters[!grepl("grey", uniqueClusters)]

    ### If the cluster label vector is unnamed, name it in the order of the data matrix columns
    if (is.null(names(consensusClusterLabels))) {
        names(consensusClusterLabels) <- colnames(dataIn)
    }

    ### Initialize nested lists to store q values, log-normalized fold change values and de genes
    qValueList <-
        rep(list(list()), length(uniqueClusters))
    logFCList <-
        rep(list(list()), length(uniqueClusters))
    deGeneList <- rep(list(list()), length(uniqueClusters))

    ### Initialise number of comparisons as n(n-1)/2
    numComparisons <-
        (length(uniqueClusters) * (length(uniqueClusters) - 1)) / 2

    ### Conduct pairwise cluster comparison to obtain q-values, log-normalized fold changes and
    ### DE genes for each comparison
    for (i in 1:(length(uniqueClusters) - 1)) {
        for (j in (i + 1):length(uniqueClusters)) {
            ### Get the cell names cell data for cluster i
            cellNamesi <-
                names(consensusClusterLabels)[which(consensusClusterLabels == uniqueClusters[i])]
            cellDatai <- dataIn[, cellNamesi]

            ### Get the cell names cell data for cluster j
            cellNamesj <-
                names(consensusClusterLabels)[which(consensusClusterLabels == uniqueClusters[j])]
            cellDataj <- dataIn[, cellNamesj]

            deCellData <- cbind(cellDatai, cellDataj)

            ### Declare variables to store p-values, q-values and log2 fold change values for each pairwise gene comparison
            pval <- NA
            qval <- NA
            logfc <- NA

            if (method == "Wilcoxon") {
                ### For each gene, conduct Wilcoxon Rank Sum test to obtain q-value and log2-normalized fold change
                for (k in 1:nrow(cellDatai)) {

                    # extract data of gene k for cluster i
                    geneDatai <- cellDatai[k, ]

                    # extract data of gene k for cluster j
                    geneDataj <- cellDataj[k, ]

                    # run wilcoxon test
                    wcTestOut <-
                        wilcox.test(geneDatai, geneDataj)

                    # initialize q value and fold change for each gene
                    pval[k] <- wcTestOut$p.value
                    logfc[k] <-
                        mean(geneDatai) - mean(geneDataj)
                }

                # check if gene satisfies mean expression threshold
                meanExprsLogicalVector <-
                    apply(deCellData, 1, function(row) {
                        mean(row[cellNamesi]) > log(meanExprsThrs) |
                            mean(row[cellNamesj]) > log(meanExprsThrs)
                    })

                # adjust p-value using Benjamini-Hochberg adjustment
                qval <-
                    p.adjust(
                        p = pval,
                        method = "BH",
                        n = nrow(cellDatai)
                    )

            } else if (method == "edgeR") {

                # check if gene satisfies mean expression threshold
                meanExprsLogicalVector <-
                    apply(deCellData, 1, function(row) {
                        mean(row[cellNamesi]) > log(meanExprsThrs) |
                            mean(row[cellNamesj]) > log(meanExprsThrs)
                    })

                # Create DGE object
                dgeObj <- edgeR::DGEList(counts=deCellData, group = c(rep(1, ncol(cellDatai)), rep(-1, ncol(cellDataj))))

                # Estimate common and tag-wise dispersion for the pair of clusters
                dgeObj <- edgeR::estimateCommonDisp(dgeObj)
                dgeObj <- edgeR::estimateTagwiseDisp(dgeObj)

                # Calculate norm factors for data
                dgeObj <- edgeR::calcNormFactors(dgeObj, method = "none")

                # Perform exact test
                etObj <- edgeR::exactTest(object = dgeObj)

                # Save p-value results
                pval <- etObj$table$PValue

                # adjust p-value using Benjamini-Hochberg adjustment
                qval <-
                    p.adjust(
                        p = pval,
                        method = "BH",
                        n = nrow(cellDatai)
                    )

                log2fc <- etObj$table$logFC

            } else {
                ### exit from function if no method is chosen
                print("Incorrect method chosen.")
                return(NULL)
            }

            ### determine if a gene is a DE gene based on thresholds
            deGeneLogicalVector <-
                qval < qValThrs &
                abs(logfc) > log(fcThrs)

            # filter de genes by mean expression
            deGeneLogicalVector <- deGeneLogicalVector & meanExprsLogicalVector

            print(paste0(
                uniqueClusters[i],
                ", ",
                uniqueClusters[j],
                " DE genes: ",
                sum(deGeneLogicalVector)
            ))

            ### store q-values, log-normalized fold change and DE genes
            qValueList[[i]][[j]] <- qval
            logFCList[[i]][[j]] <- logfc
            deGeneList[[i]][[j]] <-
                rownames(cellDatai)[deGeneLogicalVector]

            if (sum(deGeneLogicalVector) <= 1)
                next

            de_genes <- deGeneList[[i]][[j]]


        }
    }

    ### Name the q-value, log2 fold change and DE gene lists by cluster names
    names(qValueList) <- uniqueClusters
    names(logFCList) <- uniqueClusters
    names(deGeneList) <- uniqueClusters

    saveRDS(object = qValueList, file = "qValueList.rds")
    saveRDS(object = logFCList, file = "logFCList.rds")
    saveRDS(object = deGeneList, file = "deGeneList.rds")

    ### Obtain union of all de genes
    # initialize empty DE gene union vector
    deGeneUnion <- c()

    # For each pair of clusters
    for (i in 1:(length(uniqueClusters) - 1)) {
        for (j in (i + 1):length(uniqueClusters)) {

            # Rank the DE genes for each comparison by fold change and take the top 30

            de_logfc <- logFCList[[i]][[j]][which(rownames(dataIn) %in% deGeneList[[i]][[j]])]
            names(de_logfc) <- rownames(dataIn)[which(rownames(dataIn) %in% deGeneList[[i]][[j]])]
            sorted_de_logfc <- sort(x = abs(de_logfc), decreasing = T)

            if(length(sorted_de_logfc) > 30) {
                de_genes <- names(sorted_de_logfc)[1:30]
            } else {
                de_genes <- names(sorted_de_logfc)
            }

            deGeneUnion <-
                union(deGeneUnion, de_genes)
        }
    }

    print(str(deGeneUnion))

    saveRDS(deGeneUnion, file = "deGeneUnion.rds")

    ### Compute PCA + euclidean distance of cells based on the union of DE genes
    pca.data <- irlba::prcomp_irlba(x = t(dataIn[deGeneUnion, ]), n = min(length(deGeneUnion), 15), center = TRUE, scale. = FALSE)$x

    d <- dist(pca.data, method = "euclidean")

    ### Compute 2-level correlation distance of cells based on the union of DE genes
    # d = as.dist(1 - cor(dataIn[deGeneUnion, ], method = "pearson"))

    ### Build dendrogram using distance matrix d
    if(require(fastcluster)){
        cellTree = fastcluster::hclust(d, method = "ward.D2")
    }else{
        cellTree = stats::hclust(d, method = "ward.D2")
    }

    # Initialize empty dynamic colors list
    dynamicColorsList <- list()

    ### For all values of deepsplit specified
    for (dsv in deepSplitValues) {
        ### Use WGCNA dynamic tree cutting to obtain clusters
        dynamicGroups = dynamicTreeCut::cutreeDynamic(
            dendro = cellTree,
            distM = as.matrix(d),
            deepSplit = dsv,
            pamStage = FALSE,
            minClusterSize = minClusterSize
        )
        dynamicColors = WGCNA::labels2colors(dynamicGroups)

        dynamicColorsList[[paste("deepsplit:", dsv)]] <-
            dynamicColors

    }

    ### Name dynamic colors list
    names(dynamicColorsList) <- paste("deepsplit:", deepSplitValues)

    ### Compute number of detected genes for each cell in the input data
    nodg <- c()
    for (i in 1:length(colnames(dataIn))) {
        nodg[i] <- length(which(dataIn[, i] > 0))
    }

    ### Create and initialise object to return to function caller
    returnObj = list(
        "deGeneUnion" = deGeneUnion,
        "cellTree" = cellTree,
        "dynamicColors" = dynamicColorsList
    )

    ### Save DE object
    saveRDS(object = returnObj, file = filename)

    ### Plot DE Gene Plot
    cellTypeDEPlot(
        dataMatrix = dataIn[deGeneUnion,],
        nodg = nodg,
        cellTree = cellTree,
        clusterLabels = consensusClusterLabels,
        dynamicColorsList = dynamicColorsList,
        colScheme = "violet",
        filename = plotName
    )

    return(returnObj)
}
