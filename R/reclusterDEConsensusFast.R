#' @title Recluster consensus clusters using improved DE gene analysis
#'
#' @author ranjanb, schmidtf
#'
#' @param dataMatrix the log-transformed and normalized scRNAseq genes x cells matrix
#' @param consensusClusterLabels consensus cell type labels for each cell
#' @param method Method used for DE gene statistical test
#' @param qValThrs maximum q-value threshold for statistical test
#' @param fcThrs minimum fold-change threshold for DE gene criterion (natural log as by Seurat convention)
#' @param deepSplitValues vector of WGCNA tree cutting deepsplit parameters
#' @param minClusterSize specifies the type of minimum cluster size
#' @param minPerCent Minium percentage of cells expressing a gene within a cluster
#' @param filename name of DE gene object file
#' @param plotName name of DE Heatmap file
#' @param NumbertopDEGenes number of top DE genes to consider for final clustering
#' @param nCores nuber of cores to use for DE gene calculation
#'
#' @return list containing vector of DE genes, clustering tree and dynamic color list.
#'
#' @export
#'
reclusterDEConsensusFast <- function(dataMatrix,
                                     consensusClusterLabels,
                                     method = "wilcox",
                                     qValThrs=0.1,
                                     logFCThrs=0.5,
                                     deepSplitValues = 1:4,
                                     minClusterSize = 10,
                                     minPerCent = 20,
                                     filename = "de_gene_object.rds",
                                     plotName = "DE_Heatmap",
                                     NumbertopDEGenes = 30,
                                     nCores = 1) {
    #libraries required for parallel execution
    require(foreach)
    require(doParallel)
    require(dplyr)

    ### Use only clusters with number of cells > minimum cluster size for DE gene calling
    which(table(consensusClusterLabels) > minClusterSize)
    colorCounts <- table(consensusClusterLabels)

    ### Extract unique cluster labels
    uniqueClusters <- names(colorCounts[colorCounts > minClusterSize])

    # ignore the grey cluster because it represents unclustered cells
    uniqueClusters <- uniqueClusters[!grepl("grey", uniqueClusters)]
    numberUniqueClusters <- length(uniqueClusters)

    ### If the cluster label vector is unnamed, name it in the order of the data matrix columns
    if (is.null(names(consensusClusterLabels))) {
        names(consensusClusterLabels) <- colnames(dataMatrix)
    }


    ### Initialise number of comparisons as n(n-1)/2
    numComparisons <- (numberUniqueClusters* (numberUniqueClusters - 1)) / 2

    ### Conduct pairwise cluster comparison to obtain q-values, log-normalized fold changes and
    ### DE genes for each comparison
    cl <- makeCluster(nCores) #not to overload your computer
    registerDoParallel(cl)

    #Parallel outer loop
    deGenes <- foreach (i = 1:(numberUniqueClusters - 1),.combine=rbind) %dopar% {
        library(tidyverse)
        library(ROCR)
        library(dplyr)
        library(tidyr)
        library(ggpubr)
        library(cluster)
        #Function definition instead of export
        ComputePairWiseDE <- function(object, cells.1 = NULL, cells.2 = NULL, features = NULL,
                                      logfc.threshold = log(1.5), test.use = "wilcox", min.pct = 0.25,
                                      min.diff.pct = -Inf, verbose = TRUE, only.pos = FALSE, max.cells.per.ident = Inf,
                                      random.seed = 1, min.cells.group = 3, pseudocount.use = 1, MeanExprsThrs = 0, p.adjust.methods = "BH"){
            ## for Wilcox test
            WilcoxDETest <- function(data.use, cells.1, cells.2, verbose = TRUE){
                group.info <- data.frame(row.names = c(cells.1, cells.2))
                group.info[cells.1, "group"] <- "Group1"
                group.info[cells.2, "group"] <- "Group2"
                group.info[, "group"] <- factor(x = group.info[, "group"])
                data.use <- data.use[, rownames(x = group.info), drop = FALSE]
                p_val <- sapply(
                    X = 1:nrow(x = data.use),
                    FUN = function(x) {
                        return(wilcox.test(data.use[x, ] ~ group.info[, "group"])$p.value)
                    }
                )
                return(data.frame(p_val, row.names = rownames(x = data.use)))
            }

            bimodLikData <- function(x, xmin = 0) {
                x1 <- x[x <= xmin]
                x2 <- x[x > xmin]
                xal <- MinMax(
                    data = length(x = x2) / length(x = x),
                    min = 1e-5,
                    max = (1 - 1e-5)
                )
                likA <- length(x = x1) * log(x = 1 - xal)
                if (length(x = x2) < 2) {
                    mysd <- 1
                } else {
                    mysd <- sd(x = x2)
                }
                likB <- length(x = x2) *
                    log(x = xal) +
                    sum(dnorm(x = x2, mean = mean(x = x2), sd = mysd, log = TRUE))
                return(likA + likB)
            }
            DifferentialLRT <- function(x, y, xmin = 0) {
                lrtX <- bimodLikData(x = x)
                lrtY <- bimodLikData(x = y)
                lrtZ <- bimodLikData(x = c(x, y))
                lrt_diff <- 2 * (lrtX + lrtY - lrtZ)
                return(pchisq(q = lrt_diff, df = 3, lower.tail = F))
            }
            DiffExpTest <- function(data.use, cells.1,  cells.2, verbose = TRUE) {
                p_val <- unlist(
                    x = sapply(
                        X = 1:nrow(x = data.use),
                        FUN = function(x) {
                            return(DifferentialLRT(
                                x = as.numeric(x = data.use[x, cells.1]),
                                y = as.numeric(x = data.use[x, cells.2])
                            ))
                        }
                    )
                )
                to.return <- data.frame(p_val, row.names = rownames(x = data.use))
                return(to.return)
            }
            ## for ROC test
            DifferentialAUC <- function(x, y) {
                prediction.use <- prediction(
                    predictions = c(x, y),
                    labels = c(rep(x = 1, length(x = x)), rep(x = 0, length(x = y))),
                    label.ordering = 0:1
                )
                perf.use <- performance(prediction.obj = prediction.use, measure = "auc")
                auc.use <- round(x = perf.use@y.values[[1]], digits = 3)
                return(auc.use)
            }
            AUCMarkerTest <- function(data1, data2, mygenes, print.bar = TRUE) {
                myAUC <- unlist(x = lapply(
                    X = mygenes,
                    FUN = function(x) {
                        return(DifferentialAUC(
                            x = as.numeric(x = data1[x, ]),
                            y = as.numeric(x = data2[x, ])
                        ))
                    }
                ))
                myAUC[is.na(x = myAUC)] <- 0
                iterate.fxn <- ifelse(test = print.bar, yes = pblapply, no = lapply)
                avg_diff <- unlist(x = iterate.fxn(
                    X = mygenes,
                    FUN = function(x) {
                        return(
                            ExpMean(
                                x = as.numeric(x = data1[x, ])
                            ) - ExpMean(
                                x = as.numeric(x = data2[x, ])
                            )
                        )
                    }
                ))
                toRet <- data.frame(cbind(myAUC, avg_diff), row.names = mygenes)
                toRet <- toRet[rev(x = order(toRet$myAUC)), ]
                return(toRet)
            }
            MarkerTest <- function(data.use, cells.1, cells.2, verbose = TRUE){
                to.return <- AUCMarkerTest(
                    data1 = data.use[, cells.1, drop = FALSE],
                    data2 = data.use[, cells.2, drop = FALSE],
                    mygenes = rownames(x = data.use),
                    print.bar = verbose
                )
                to.return$power <- abs(x = to.return$myAUC - 0.5) * 2
                return(to.return)
            }
            # Differential expression testing using Student's t-test
            # Identify differentially expressed genes between two groups of cells using the Student's t-test
            DiffTTest <- function(data.use, cells.1, cells.2, verbose = TRUE) {
                p_val <- unlist(
                    x = sapply(
                        X = 1:nrow(data.use),
                        FUN = function(x) {
                            t.test(x = data.use[x, cells.1], y = data.use[x, cells.2])$p.value
                        }
                    )
                )
                to.return <- data.frame(p_val,row.names = rownames(x = data.use))
                return(to.return)
            }


            features <- features %||% rownames(x = object)
            # error checking
            if (length(x = cells.1) == 0) {
                stop("Cell group 1 is empty - identity of group 1 need to be defined ")
            }
            else if (length(x = cells.2) == 0) {
                stop("Cell group 2 is empty - identity of group 2 need to be defined ")
                return(NULL)
            }
            else if (length(x = cells.1) < min.cells.group) {
                stop("Cell group 1 has fewer than ", min.cells.group, " cells")
            }
            else if (length(x = cells.2) < min.cells.group) {
                stop("Cell group 2 has fewer than ", min.cells.group, " cells")
            }
            else if (any(!cells.1 %in% colnames(x = object))) {
                bad.cells <- colnames(x = object)[which(x = !as.character(x = cells.1) %in% colnames(x = object))]
                stop(
                    "The following cell names provided to cells.1 are not present: ",
                    paste(bad.cells, collapse = ", ")
                )
            } else if (any(!cells.2 %in% colnames(x = object))) {
                bad.cells <- colnames(x = object)[which(x = !as.character(x = cells.2) %in% colnames(x = object))]
                stop(
                    "The following cell names provided to cells.2 are not present: ",
                    paste(bad.cells, collapse = ", ")
                )
            }

            # feature selection (based on percentages)
            thresh.min <- 0
            pct.1 <- round(
                x = 100*Matrix::rowSums(x = object[features, cells.1, drop = FALSE] > thresh.min) /
                    length(x = cells.1),
                digits = 16
            )
            pct.2 <- round(
                x = 100*Matrix::rowSums(x = object[features, cells.2, drop = FALSE] > thresh.min) /
                    length(x = cells.2),
                digits = 16
            )
            object.alpha <- cbind(pct.1, pct.2)
            colnames(x = object.alpha) <- c("pct.1", "pct.2")
            alpha.min <- apply(X = object.alpha, MARGIN = 1, FUN = max)
            names(x = alpha.min) <- rownames(x = object.alpha)
            features <- names(x = which(x = alpha.min > min.pct))
            if (length(x = features) == 0) {
                warning("No features pass min.pct threshold")
                return(NULL)
            }
            alpha.diff <- alpha.min - apply(X = object.alpha, MARGIN = 1, FUN = min)
            features <- names(
                x = which(x = alpha.min > min.pct & alpha.diff > min.diff.pct)
            )
            if (length(x = features) == 0) {
                warning("No features pass min.diff.pct threshold")
                return(NULL)
            }

            # feature selection (based on average difference)
            mean.fxn <- function(x) {
                return(log(x = mean(x = expm1(x = x)) + pseudocount.use))
            }
            object.1 <- apply(
                X = object[features, cells.1, drop = FALSE],
                MARGIN = 1,
                FUN = mean.fxn
            )
            object.2 <- apply(
                X = object[features, cells.2, drop = FALSE],
                MARGIN = 1,
                FUN = mean.fxn
            )
            total.diff <- (object.1 - object.2)

            # feature selection (based on mean exprssion threshold)
            features = names(x = which(x = expm1(object.1) > MeanExprsThrs  | expm1(object.2) > MeanExprsThrs ))
            if (length(x = features) == 0) {
                warning("No features pass log mean exprssion threshold")
                return(NULL)
            }

            # feature selection (based on logfc threshold)
            features.diff <- if (only.pos) {
                names(x = which(x = total.diff > logfc.threshold))
            } else {
                names(x = which(x = abs(x = total.diff) > logfc.threshold))
            }
            features <- intersect(x = features, y = features.diff)
            if (length(x = features) == 0) {
                warning("No features pass logfc.threshold threshold")
                return(NULL)
            }

            # sampling cell for DE computation
            if (max.cells.per.ident < Inf) {
                set.seed(seed = random.seed)
                # Should be cells.1 and cells.2?
                if (length(x = cells.1) > max.cells.per.ident) {
                    cells.1 <- sample(x = cells.1, size = max.cells.per.ident)
                }
                if (length(x = cells.2) > max.cells.per.ident) {
                    cells.2 <- sample(x = cells.2, size = max.cells.per.ident)
                }
            }

            # perform DE
            de.results <- switch(
                EXPR = test.use,
                'wilcox' = WilcoxDETest(
                    data.use = object[features, c(cells.1, cells.2), drop = FALSE],
                    cells.1 = cells.1,
                    cells.2 = cells.2,
                    verbose = verbose
                ),
                'bimod' = DiffExpTest(
                    data.use = object[features, c(cells.1, cells.2), drop = FALSE],
                    cells.1 = cells.1,
                    cells.2 = cells.2,
                    verbose = verbose
                ),
                'roc' = MarkerTest(
                    data.use = object[features, c(cells.1, cells.2), drop = FALSE],
                    cells.1 = cells.1,
                    cells.2 = cells.2,
                    verbose = verbose
                ),
                't' = DiffTTest(
                    data.use = object[features, c(cells.1, cells.2), drop = FALSE],
                    cells.1 = cells.1,
                    cells.2 = cells.2,
                    verbose = verbose
                ),
                stop("Unknown test: ", test.use)
            )

            diff.col <- "avg_logFC"
            de.results[, diff.col] <- total.diff[rownames(x = de.results)]
            de.results <- cbind(de.results, object.alpha[rownames(x = de.results), , drop = FALSE])

            if (only.pos) {
                de.results <- de.results[de.results[, diff.col] > 0, , drop = FALSE]
            }

            if (test.use == "roc") {
                de.results <- de.results[order(-de.results$power, -de.results[, diff.col]), ]
            } else {
                de.results <- de.results[order(de.results$p_val, -de.results[, diff.col]), ]
                de.results$p_val_adj = p.adjust(
                    p = de.results$p_val,
                    method = p.adjust.methods
                )
            }
            de.results$Gene<-row.names(de.results)
            row.names(de.results)<-NULL
            return(de.results)
        }

        marker.genes=c()
        #Sequential inner loop
        for (j in (i + 1):numberUniqueClusters) {
            ### Get the cell names cell data for cluster i
            cellNamesi <- names(consensusClusterLabels)[which(consensusClusterLabels == uniqueClusters[i])]
            cellDatai <- dataMatrix[, cellNamesi]
            ### Get the cell names cell data for cluster j
            cellNamesj <- names(consensusClusterLabels)[which(consensusClusterLabels == uniqueClusters[j])]
            cellDataj <- dataMatrix[, cellNamesj]

            tmp = cbind(Cluster1=uniqueClusters[i],Cluster2=uniqueClusters[j],
                        ComputePairWiseDE(object = as.matrix(dataMatrix),
                                          cells.1 = cellNamesi,
                                          cells.2 = cellNamesj,
                                          logfc.threshold = logFCThrs,
                                          test.use = method,
                                          min.pct = minPerCent))

            #Filter for q-value, fc filtering internally
            if (dim(tmp)[1]>1){
                marker.genes<-rbind(marker.genes,tmp[tmp$p_val_adj < qValThrs,])

            }
        }
        marker.genes
    }

    stopCluster(cl)

    deGenes$avg_logFC<-abs(deGenes$avg_logFC)
    # saveRDS(deGenes,"PairwiseDifferentiallyExpressedGenes.RDS")

    ### Obtain union of all de genes
    # initialize empty DE gene union vector
    tmp<-deGenes %>% group_by(Cluster1,Cluster2) %>% top_n(n=NumbertopDEGenes,wt=avg_logFC)
    deGeneUnion<-unique(tmp$Gene)
    print(str(deGeneUnion))

    # saveRDS(deGeneUnion, file = "deGeneUnion.rds")

    ### Compute PCA + euclidean distance of cells based on the union of DE genes
    pca.data <- irlba::prcomp_irlba(x = t(dataMatrix[deGeneUnion, ]), n = min(length(deGeneUnion), 15), center = TRUE, scale. = FALSE)$x

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
    
    #DataFrame holding clusterinformation depending on deepsplit values
    deepSplitInfo<-c()

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
        
        deepSplitInfo<-rbind(deepSplitInfo,cbind("DeepSplit"=dsv,"NumbersOfClusters"=length(unique(dynamicGroups)),"SI"=mean(summary(cluster::silhouette(dynamicGroups,dmatrix=as.matrix(d)))$clus.avg.widths)))
    }

    ### Name dynamic colors list
    names(dynamicColorsList) <- paste("deepsplit:", deepSplitValues)

    ### Compute number of detected genes for each cell in the input data
    nodg <- c()
    for (i in 1:length(colnames(dataMatrix))) {
        nodg[i] <- length(which(dataMatrix[, i] > 0))
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
        dataMatrix = dataMatrix[deGeneUnion,],
        nodg = nodg,
        cellTree = cellTree,
        clusterLabels = consensusClusterLabels,
        dynamicColorsList = dynamicColorsList,
        colScheme = "violet",
        filename = plotName
    )
    


    return(returnObj)
}

