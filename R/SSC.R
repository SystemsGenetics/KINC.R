#' Uses the edge sample string to return the list of sample indexes in the edge.
#'
#' @param i
#'   The index of the edge in the network
#' @param net
#'   A network dataframe containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @return
#'   An array of sample indexes that are part of the cluster. Because these
#'   indexes are derived from the edge sample string they indexes should
#'   correspond to the order of samples in the expression matrix.
#' @export
getEdgeSamples <- function(i, net) {
  # Convert the sample string to a numerical vector.
  edge = net[i,]
  sample_str = edge$Samples
  edge_samples = as.numeric(strsplit(sample_str, "")[[1]])
  sample_indexes = which(edge_samples == 1)
  return (sample_indexes)
}
#' Uses the edge sample string to return the list of sample indexes not in the edge.
#'
#' @param i
#'   The index of the edge in the network
#' @param net
#'   A network dataframe containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @return
#'   An array of sample indexes that are NOT part of the cluster. Because these
#'   indexes are derived from the edge sample string they indexes should
#'   correspond to the order of samples in the expression matrix.
#'
#' @export
getNonEdgeSamples <- function(i, net) {
  # Convert the sample string to a numerical vector.
  edge = net[i,]
  sample_str = edge$Samples
  edge_samples = as.numeric(strsplit(sample_str, "")[[1]])
  sample_indexes = which(edge_samples != 1)
  return (sample_indexes)
}
#' Uses the edge sample string to return the list of sample that are missing.
#'
#' Note that a edge will never have missing samples. These missing samples
#' are from the pairwise comparision between the two genes of the edge.
#'
#' @param i
#'   The index of the edge in the network
#' @param net
#'   A network dataframe containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @return
#'   An array of sample indexes that are NOT part of the cluster. Because these
#'   indexes are derived from the edge sample string they indexes should
#'   correspond to the order of samples in the expression matrix.
#'
#' @export
#'
getMissingEdgeSamples <- function(i, net) {
  # Convert the sample string to a numerical vector.
  edge = net[i,]
  sample_str = edge$Samples
  edge_samples = as.numeric(strsplit(sample_str, "")[[1]])
  sample_indexes = which(edge_samples == 9)
  return (sample_indexes)
}

#' Converts and edge sample string into an array.
#'
#' @param i
#'   The index of the edge in the network
#' @param net
#'   A network dataframe containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @return
#'   An array whose values are  1 if the sample is in the cluster,
#'   0 if the sample is in another cluster, 9 if the sample is missing and
#'   7 if the sample was removed as a pairwise outlier or 8 if the sample
#'   was removed as a cluster outlier.
#' @export
getSampleStringArray <-function(i, net){
  edge = net[i,]
  sample_str = edge$Samples
  edge_samples = as.numeric(strsplit(sample_str, "")[[1]])
  return(edge_samples)
}
#' Creates an Edge Expression Matrix (EEM) from a network data frame.
#'
#' The EEM is an n x m*2 array were n is the number of samples in the annotation matrix and
#' m is the number of edges in the network. Each edge is represented in the matrix by two columns
#' with each containing gene expression of the source and target genes respectively.
#'
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param osa
#'   The ordered sample annotation data frame.
#' @param ematrix
#'   The expression matrix data frame.
#'
#' @return
#'    An nxm*2 dataframe contianing the EEM
#' @export
getEdgeExpressionMatrix <- function(net, osa, ematrix) {
  samples = osa$Sample
  num_edges = dim(net)[1]
  sm = data.frame(matrix(NA, nrow=length(samples), ncol=num_edges*2))
  row.names(sm) = samples
  cnames = vector(mode="character", length=num_edges*2)
  for (i in 1:num_edges) {
    print(paste(i,'of',num_edges))

    # Get the sample string in array form.
    edge_samples = getSampleStringArray(i, net)

    # Get the gene expression for the source and target genes of the edge.
    source = net[i, 'Source']
    target = net[i, 'Target']
    sourceExp = ematrix[source,]
    targetExp = ematrix[target,]

    # Create names for the EEM columns.
    cnames[i*2-1] = paste('E', i, "_", source, sep="")
    cnames[i*2] = paste('E', i, "_", target, sep="")

    # Remove expression from non-edge samples
    sourceExp[which(edge_samples != 1)] = 0
    targetExp[which(edge_samples != 1)] = 0

    # Add the new columns.
    sm[,i*2-1] = t(sourceExp)
    sm[,i*2] = t(targetExp)
  }
  colnames(sm) = cnames
  return(sm)
}
#' Performs hierarchical clustering of edges in a network based on their sample compositions.
#'
#' This function uses the dist function to calucate a distance and the
#' hclust function to generate the dendrogram.
#'
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param distMethod
#'   The method to be provided to the dist function.
#' @param hclustMethod
#'   The method to be provided to the hclust method.
#'
#' @return
#'    An object of class hclust which describes the tree produced by the
#'    clustering process.
#'
clusterEdges = function(net, distMethod = "manhattan", hclustMethod = "ward.D2") {

  # Convert the samples strings into a matrix.
  num_samples = nchar(net$Samples[1])
  sample_strs = net$Samples
  samples = sapply(sample_strs, FUN=function(x) {
    s = as.numeric(strsplit(x, "")[[1]])
    s[which(s != 1)] = 0
    return(s)
  })
  colnames(samples) = c()
  samples=t(samples)

  # Perform clustering on the sample tree
  sample_dist  = dist(samples, method = distMethod)
  tree = hclust(sample_dist, method = hclustMethod)

  return(tree)
}

#' Generates a matrix containing edge sample strings digits.
#'
#' The returned matrix contains as many rows as there are edges
#' in the network and as many columns as their are samples. Each
#' column corresponds to a sample in the same order as the
#' sample strings from the network.
#'
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @return
#'   A matrix
#'
#' @export
getSampleMatrix = function(net) {
  # Convert the samples strings into a matrix.
  num_samples = nchar(net$Samples[1])
  sample_strs = net$Samples
  samples = sapply(sample_strs, FUN=function(x) {
    s = as.numeric(strsplit(x, "")[[1]])
    s[which(s != 1)] = 0
    return(s)
  })
  colnames(samples) = c()
  samples = t(samples)

  return(samples)
}

#' Draws a heatmap with dendrogram with clusters identified.
#'
#' @param sampleMatrix
#'   The sample matrix as created by the getSampleMatrix() function.
#' @param tree
#'   An instance of an hclust object as created by the clusterEdges() function.
#' @param osa
#'   The sample annotation matrix as created by the loadSampleAnnotations()
#'   function.
#' @param fieldOrder
#'   A vector containing a list of sample attribute names for reordering of
#'   of the samples. Each element of this vector should be the name of
#'   a column in the osa matrix.  The sorting occurs first by the first
#'   element, then by the second, etc.
#' @param outfile_prefix
#'   An output file prefix to add to the beginnging of the output heatmap image file.
#'   If this argument is not provided then the figure will be plotted to the
#'   default devide rather than to a file.
drawEdgeTreeHeatMap = function(sampleMatrix, tree, osa, fieldOrder, outfile_prefix = NA) {

  # Reorder samples in the sampleMatrix according to the fieldOrder argument.
  sample_order = eval(parse(text=paste('order(osa$', paste(fieldOrder, collapse=' ,osa$'), ")", sep="")))
  sampleMatrix2 = sampleMatrix[, c(sample_order)]

  categories = do.call(paste, list(c(osa[, fieldOrder[1]]), sep="-"))
  num_categories = length(unique(categories))
  osa$hmap_categories = categories
  tColors = data.frame(
    Field = unique(categories),
    Color = rgb(runif(num_categories), runif(num_categories), runif(num_categories))
  )
  colColors = as.character(merge(osa[sample_order,], tColors, by.x="hmap_categories", by.y="Field", sort=FALSE)$Color)

  if (!is.na(outfile_prefix)) {
    outfile = paste(outfile_prefix, paste(fieldOrder, collapse="-"), num_clusters, "png", sep=".")
    png(filename = outfile, width=3000, height=21000, res=300)
  }
  heatmap.2(sampleMatrix2,
            Rowv=as.dendrogram(tree),
            Colv=FALSE,
            dendrogram = 'row',
            col = c("blue", "green"),
            breaks = c(-1, 0, 1),
            trace = 'none',
            key = FALSE,
            ColSideColors = colColors,
            labRow = NULL,
            labCol = NULL
  )
  if (!is.na(outfile_prefix)) {
    dev.off()
  }
}

#' Draws a heatmap with dendrogram of a module in the netework.
#'
#' @param edge_indexes
#'   The index of the edges in the network that comprise the module.
#' @param osa
#'   The sample annotation matrix. One column must contain the header 'Sample'
#'   and the remaining colums correspond to an annotation type.  The rows
#'   of the anntation columns should contain the annotations.
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param fieldOrder
#'   A vector containing a list of sample attribute names for reordering of
#'   of the samples. Each element of this vector should be the name of
#'   a column in the osa matrix.  The sorting occurs first by the first
#'   element, then by the second, etc.
#' @export
drawEdgeListHeatMap = function(edge_indexes, osa, net, fieldOrder) {
  ce.tree = clusterEdges(net[edge_indexes,])
  sm = getSampleMatrix(net[edge_indexes,])
  drawEdgeTreeHeatMap(sm, ce.tree, osa, fieldOrder)
}


#' Removes edges after KINC has been run that are supsected to be biased.
#'
#' This function tests for two modes of bias: collinearity and missingness.  In
#' either case, an edge can incorrectly appear to be signficantly associated
#' with an experimental condition.
#'
#' To correct for collinearity, a Welch one-way ANOVA test is performed for each
#' gene. This tests if the distribution of the edge samples is different from the
#' distribution of non-edge samples in both genes.  Collinearity can be detected
#' when at least one of the gene shows no signficant difference between the two.
#'
#' Missigness can bias an edge towards a condition if one gene that is highly
#' condition-specific has missing values for samples that did not experience the
#' condition but the other gene is not.  In order for KINC to perform correlation
#' analysis it must remove any samples that are missing in either gene.  To
#' correct for this, a paired t-test is performed.
#'
#' @param net
#'   A network dataframe containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param ematrix
#'   The expression matrix data frame.
#' @param th
#'   The signficance threshold for bias testing.  Edges with P-values below this
#'   value, for any of the tests performed, are kept.
#' @param progressBar
#'   Set to FALSE to repress progress bar
#'
#' @export
filterBiasedEdges <- function(net, ematrix, th = 1e-3, progressBar = TRUE) {

  # Intialize a vector that we'll use to indicate which rows to keep.
  keep = rep(FALSE, dim(net)[1])
  # net$SourceDiff = NA
  # net$TargetDiff = NA
  # net$MissingDiff = NA

  # Intialize the progress bar
  if (progressBar){
    pb = txtProgressBar(min = 0, max = nrow(net), style = 3)
  }

  for (i in 1:dim(net)[1]) {
    gene1 = net$Source[i]
    gene2 = net$Target[i]
    cluster_samples = getEdgeSamples(i, net)
    non_cluster_samples = getNonEdgeSamples(i,net)
    missing_samples = getMissingEdgeSamples(i,net)

    if (progressBar){
      setTxtProgressBar(pb, i)
    }

    # First, performs a Welch's ANOVA test for both genes of every network edge.
    # The Welch's Anova test is used for heteroscedastic data (non-equal variance)
    # and does not assume that both categories have the same distribution.  This
    # test is used for each edge of the network to test if the mean of the
    # in-cluster and out-cluster samples of each gene are the same.
    # TODO: we should do a power-analysis to see how many samples we must have
    # to do thie Welch ANOVA test. For now we'll leave it at 15.
    if (length(non_cluster_samples) - length(missing_samples) >= 10 &
        length(cluster_samples) >= 10) {

      # Get the expression vectors of the samples not in the cluster
      # and remove outliers. We don't need to do this for the in-cluster
      # samples because clusters should be normal and have outliers removed
      # by KINC.
      ncs1 = as.vector(t(ematrix[gene1, non_cluster_samples]))
      ncs2 = as.vector(t(ematrix[gene2, non_cluster_samples]))
      b1 = boxplot(ncs1, plot = FALSE)
      b2 = boxplot(ncs2, plot = FALSE)
      # remove the outliers from the original expression data for the gene
      g1 = as.vector(t(ematrix[gene1,]))
      g2 = as.vector(t(ematrix[gene2,]))
      g1[which(g1 %in% b1$out)] = NA
      g2[which(g2 %in% b2$out)] = NA

      # Create an array of which samples are in/out of the cluster
      position = rep('Out',dim(ematrix)[2])
      position[cluster_samples] = 'In'

      # Perform the Welch one-way ANOVA test on in/out of gene1
      source = data.frame(Position = position, Expression = g1)
      w1 = oneway.test(Expression ~ Position, data=source, var.equal=FALSE)
      # net$SourceDiff[i] = w1$p.value

      # Perform the Welch one-way ANOVA test on in/out of gene2
      target = data.frame(Position = position, Expression = g2)
      w2 = oneway.test(Expression ~ Position, data=target, var.equal=FALSE)
      # net$TargetDiff[i] = w2$p.value

      # If the welch one-way ANOVA threshold is not significant then that
      # means that the distributions cannot be determined to be different
      # and we should not keep this edge.
      if (w1$p.value < th & w2$p.value < th) {
        keep[i] = TRUE
      }
    }

    # If this looks like a good edge, let's do one more check to make
    # sure there isn't bias in the missigness.
    if (keep[i] == TRUE) {

      # Convert the missigness of both genes to a feature where 1 == missing
      # and 0 == not missing
      g1m = as.vector(t(ematrix[gene1,]))
      g1m[which(!is.na(g1m))] = 0
      g1m[which(is.na(g1m))] = 1
      g2m = as.vector(t(ematrix[gene2,]))
      g2m[which(!is.na(g2m))] = 0
      g2m[which(is.na(g2m))] = 1

      # Perform a paired t-test to see how similar the missingness is. If the
      # p-value is signficant then that means the pattern of missing
      # values is different and there may be bias.
      # We should not keep this edge.
      t = t.test(g1m, g2m, paired = TRUE, alternative = "two.sided")
      if (!is.na(t)) {
        # net$MissingDiff[i] = t$p.value
        if (t$p.value < 0.05) {
           keep[i] = FALSE
        }
      }
    }
  }

  if (progressBar){
    close(pb)
  }

  qs = qvalue(net[[colname]], fdr.level = fdr.level, pi0 = 1)
  newcol = sub('_pVal$', '_qVal', colname)
  net[[newcol]] = qs$qvalues
  sig[[newcol]] = qs$significant
  colorder[length(colorder)+1] = newcol

  return(net[which(keep == TRUE), ])
}

#' Calculates QValues from Pvalues in the network file.
#'
#' @param net
#'   The KINC v3 network data frame.
#' @param fdr.level
#'   A level at which to control the FDR. Must be in (0,1].
#'
#' @export
#'
filterNetFDR <- function(net, fdr.level = 0.001) {
  colnames = colnames(net)
  colorder = c()
  sig = data.frame(index = 1:dim(net)[1])
  for (colname in colnames) {
    colorder[length(colorder)+1] = colname
    # Is this a p-value column?
    if (length(grep('_pVal$', c(colname), perl=TRUE)) == 1) {
      qs = qvalue(net[[colname]], fdr.level = fdr.level, pi0 = 1)
      newcol = sub('_pVal$', '_qVal', colname)
      net[[newcol]] = qs$qvalues
      sig[[newcol]] = qs$significant
      colorder[length(colorder)+1] = newcol
    }
  }
  return (net[which(apply(sig, 1, any)),colorder])
}

#' Filters a KINC v3 network data frame by p-value, q-value or r-squared value.
#'
#' Any edge where all applicable scores fail to pass the threshold will
#' be excluded.  If any score meets the requirement within the edge then
#' the edge is kept.
#'
#' @param net
#'   The KINC v3 network data frame.
#' @param p_th
#'   The p-value threshold. Any p-value column with a value below will
#'   result in the edge being retained.
#' @param q_th
#'   The q-value threshold. Any p-value column with a value below will
#'   result in the edge being retained. Q-values are added to the
#'   network data frame using the calculateNetQValues function.
#' @param r_th
#'   The p-value threshold. Any p-value column with a value below will
#'   result in the edge being retained.
#'
#'@export
filterInsignficantEdges <- function(net, p_th = NA, q_th = 1e-3, r_th = NA) {
  colnames = colnames(net)
  keep = rep(FALSE, dim(net)[1])
  for (i in 1:dim(net)[1]) {
    for (colname in colnames) {
      # Is this a p-value column?
      if (!is.na(p_th) &
          length(grep('_pVal$', c(colname), perl=TRUE)) == 1 &
          net[i,colname] <= p_th) {
        keep[i] = TRUE
      }
      if (!is.na(q_th) &
          length(grep('_qVal$', c(colname), perl=TRUE)) == 1 &
          net[i,colname] <= q_th) {
        keep[i] = TRUE
      }
      if (!is.na(r_th) &
          length(grep('_RSqr$', c(colname), perl=TRUE)) == 1 &
          net[i,colname] >= r_th) {
        keep[i] = TRUE
      }
    }
  }
  return (net[which(keep==TRUE),])
}

