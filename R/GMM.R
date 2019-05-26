#' Performs pair-wise Gaussian Mixture Models (GMM) on a pair of genes
#'
#' @param source
#'   The name of the source gene.
#' @param target
#'   The name of the target gene.
#' @param ematrix
#'   The expression matrix data frame.
#' @param method
#'   The correlation method to perform (see the cor function for valid options)
#' @param th
#'   The similarity value threshold. Any correlation value greater than or equal to the
#'   absolute value of this number is kept as an edge. Derfault = 0.5.
#' @param minc
#'   The minimum number of samples that must be present to identify a cluster. Default = 30.
#' @param plot
#'   Boolean indicating if the Rmixmod clusters scatterplot should be drawn. Default = FALSE
#' @return
#'   Returns a data frame of edges in the same format as the network returned by
#'   the loadNetwork function.
#' @export
#'
getPairGMMEdges = function(source, target, ematrix, minc = 30,
  method="spearman", th=0.5, plot=FALSE) {

  # Get the source and target genes from the expression matrix
  x =  t(ematrix[source,])
  y =  t(ematrix[target,])

  # We'll use this data frame to calculate the sample strings.
  S = data.frame(x=x, y=y, stype=0, index=1:length(x), cluster=NA)
  colnames(S) = c('x', 'y', 'stype', 'index', 'cluster')

  # Remove missing values
  x2 = x[!(is.na(x) | is.na(y))]
  y2 = y[!(is.na(x) | is.na(y))]

  # Missing values get a number 9.
  S[which(is.na(x) | is.na(y)), 'stype'] = 9

  # Remove outliers
  bx = boxplot(x2, plot = FALSE)
  by = boxplot(y2, plot = FALSE)
  x_outliers = which(x2 %in% bx$out)
  y_outliers = which(y2 %in% by$out)
  outliers = c(x_outliers, y_outliers)
  if (length(outliers) > 0) {
    x2 = as.matrix(x[-outliers])
    y2 = as.matrix(y[-outliers])
  }

  # Pair-wise outliers get a number 7.
  S[which(S$x %in% bx$out), 'stype'] = 7
  S[which(S$y %in% by$out), 'stype'] = 7

  # Perform the GMM clustering.
  X = S[which(S$stype == 0),]
  xem = mixmodCluster(X[c("x", "y")], nbCluster=1:5, criterion="ICL")

  # Get the best set of clusters.
  b=xem['bestResult']
  if (plot) {
    plotCluster(b, X, xlab=source, ylab=target, cex.lab=1.5, cex.axis=1.5)
  }

  # Keep track of the cluster that each sample belongs to.
  X$cluster = b@partition
  S[X$index, 'cluster'] = X$cluster

  # Build the edges data frame, we'll return this when the function completes.
  edges = data.frame(Source = character(), Target = character(), sc = numeric(),
                     Interaction = character(), Cluster = numeric(),
                     Num_Clusters = numeric(), Cluster_Samples = numeric(), Missing_Samples = numeric(),
                     Cluster_Outliers = numeric(), Pair_Outliers = numeric(), Too_Low = numeric(),
                     Samples = character(), stringsAsFactors = FALSE)

  # Itearate through the clusters, remove outliers, perform correlation
  # analyis and if the correlation value is > then the threshodl, add the
  # cluster as an edge.
  for (ci in 1:b@nbCluster) {

    # Get the cluster points.
    cx = X$x[which(X$cluster == ci)]
    cy = X$y[which(X$cluster == ci)]

    # First, remove outliers from the cluster
    cbx = boxplot(cx, plot = FALSE)
    cby = boxplot(cy, plot = FALSE)
    cx_outliers = which(cx %in% cbx$out)
    cy_outliers = which(cy %in% cby$out)
    coutliers = c(cx_outliers, cy_outliers)
    if (length(coutliers) > 0) {
      cx = as.matrix(cx[-coutliers])
      cy = as.matrix(cy[-coutliers])
    }

    # Skip clusters that are too small.
    if (length(cx) < minc) {
      next
    }

    # Perform correlation.
    r = cor(cx, cy, method=method)

    # If the correlation value is > the threshold add the edge.
    if (abs(r) >= th) {

      # Create the sample string
      sample_str = S$stype
      sample_str[which(S$cluster == ci)] = 1
      sample_str[coutliers] = 8
      sample_str = paste(sample_str, sep='', collapse='')

      # Create a new edge dataframe and merge it into our edges data frame.
      edge = data.frame(Source = source, Target = target, sc = r, Interaction = 'co', Cluster = ci,
                       Num_Clusters = b@nbCluster, Cluster_Samples = length(cx),
                       Missing_Samples = length(which(S$stype == 9)),
                       Cluster_Outliers = length(which(S$stype == 8)),
                       Pair_Outliers = length(which(S$stype == 7)),
                       Too_Low = 0, Samples = sample_str, stringsAsFactors=FALSE)
      edges = rbind(edges, edge)
    }
  }
  return (edges)
}


countEmatrixModes = function(ematrix, minc = 30, plot=FALSE, progressBar = TRUE) {

  num_genes = dim(ematrix)[1]

  # Intalize the progress bar.
  if (progressBar){
    pb <- txtProgressBar(min = 0, max = num_genes, style = 3)
  }

  # Initalize the counts array that we will return.
  counts <- vector(mode="numeric", length=num_genes)

  # Perform GMM on each row of the GMM and we'll calculate the
  # modes of expression for each gene/transcript.
  for (i in 1:num_genes)  {

    # Set the progress bar
    if (progressBar){
      setTxtProgressBar(pb, i)
    }

    # Get the gene expression and initalize its counts.
    gene = t(ematrix[i,])
    gene_name = colnames(gene)[1]

    # Remove missing values
    gene = data.frame(gene_name = gene[complete.cases(gene)])

    # Only perform GMMs for rows with enough samples.
    if (dim(gene)[1] < minc) {
      next
    }

    # Perform GMMs on this edge
    xem = mixmodCluster(gene, nbCluster=1:5, criterion="ICL")

    # Get the best set of clusters.
    b=xem['bestResult']
    if (plot) {
      hist(xem)
    }

    counts[i] = b@nbCluster
  }

  # Close down the progress bar.
  if (progressBar){
    close(pb)
  }

  return(counts)
}
