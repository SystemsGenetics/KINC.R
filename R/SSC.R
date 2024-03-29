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
  edge_samples = as.numeric(strsplit(net[i, 'Samples'], "")[[1]])
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
  edge_samples = as.numeric(strsplit(net[i, 'Samples'], "")[[1]])
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
  edge_samples = as.numeric(strsplit(net[i, 'Samples'], "")[[1]])
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


#' Calcualtes the distance matrix between edges.
#'
#' This function uses the dist function to calucate a distance
#'
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param distMethod
#'   The method to be provided to the dist function.
#' @param strict
#'   If TRUE then all values in the netwrok sample string other than 1 are
#'   converted to 0.  This means that numeric values indicating
#'   missingness or removal as outliers are considered as 0.
#'
#' @return
#'    returns an object of class "dist" from the dist function.
#'
#' @export
getSampleDistance = function(net, distMethod = "manhattan", strict=TRUE) {

  # Convert the samples strings into a matrix.
  samples=getSampleMatrix(net, strict = strict)

  # Calcualte the distance matrix.
  sample_dist  = dist(samples, method = distMethod)

  return(sample_dist)
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
#' @param strict
#'   If TRUE then all values in the netwrok sample string other than 1 are
#'   converted to 0.  This means that numeric values indicating
#'   missingness or removal as outliers are considered as 0.
#' @return
#'    An object of class hclust which describes the tree produced by the
#'    clustering process.
#'
#' @export
clusterEdges = function(net, distMethod = "manhattan",
                        hclustMethod = "ward.D2", strict=TRUE) {

  # Calcualte the distance matrix.
  sample_dist  = getSampleDistance(net, distMethod, strict=strict)

  # Perform clustering on the sample dendrogram tree
  dendro = hclust(sample_dist, method = hclustMethod)

  return(dendro)
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
#'
#' @param strict
#'   If TRUE then all values in the netwrok sample string other than 1 are
#'   converted to 0.  This means that numeric values indicating
#'   missingness or removal as outliers are considered as 0.
#' @return
#'   An n x m matrix where n is the number of edges and m is
#'   the number of samples. The values are the numbers from
#'   the sample strings in the network.
#'
#' @export
getSampleMatrix = function(net, strict=TRUE) {
  # Convert the samples strings into a matrix.
  num_samples = nchar(net$Samples[1])
  sample_strs = net$Samples
  samples = sapply(sample_strs, FUN=function(x) {
    s = as.numeric(strsplit(x, "")[[1]])
    if (strict) {
      s[which(s != 1)] = 0
    }
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
#'   A vector containing a list of sample names for reordering of
#'   of the samples. Each element of this vector should be the name of
#'   a column in the osa matrix.  The sorting occurs first by the first
#'   element, then by the second, etc.
#' @param image_name
#'   The filename for saving the image. If no name is provided then
#'   the image is not saved.
#'
#' @export
drawEdgeTreeHeatMap = function(sampleMatrix, tree, osa, fieldOrder, image_name = NA) {

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

  if (!is.na(image_name)) {
    png(filename = image_name, width=21000, height=3000, res=300)
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
  if (!is.na(image_name)) {
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
#' This function tests for two modes of bias in the GMM cluster: lack of differential
#' cluster expression (DCE) and inconsistent missingness between the two genes of
#' the edge.  In either case, an edge can incorrectly appear to be signficantly
#' associated with an experimental condition.
#'
#' To identify lack of DCE, a Welch one-way ANOVA test is performed for each gene.
#' This tests if the distribution of the edge samples is different from the
#' distribution of non-edge samples in both genes.  Both genes must be
#' differentially expressed for the edge samples from the non-edge samples
#' to not be biased by one or the other.
#'
#' Missigness can bias an edge towards a condition if one gene that is highly
#' condition-specific has missing values for samples that did not experience the
#' condition but the other gene is not.  In order for KINC to perform correlation
#' analysis it must remove any samples that are missing in either gene.  To
#' correct for this, a paired t-test is performed.
#'
#' This function returns a new network dataframe with only edges that passed
#' the tests. It also adds three new columns to the input network that contains
#' the p-values of the three tests: SourceWAnova, TargetWAnova, MissingTtest.
#' These can be filtered later if needed.
#'
#' @param net
#'   A network dataframe containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param ematrix
#'   The expression matrix data frame. The expression matrix must have the samples
#'   ordered in the same order that was provided to KINC and which matches the
#'   'Samples" string column in the network.
#' @param threads
#'   The number of computational threads to use for parallel processing. By
#'   default, all but 2 cores available on the local machine will be used.
#' @param wa_th
#'   The signficance threshold for performing the Welch Anova test. This test
#'   checks for differential expression of the GMM cluster vs non-cluster for each gene,
#'   thus, two tests are performed: one for each gene in the edge.  Edges with both tests
#'   having P-values below this value, are kept.
#' @param wa_samples
#'   The list of numeric column indexes in the expression matrix that can be considered
#'   "control" samples. When multiple conditions are present in an experiment we will
#'   want to perform the Welch Anova test on the cluster vs the control (or base) rather
#'   than on all of the non-cluster samples. This is because a gene may be differentically
#'   expressed in one gene in a condition other than the one being tested and this can
#'   result in a cluster showing as signficantly different because the variance of the
#'   outgroup is wider in the gene with differential expression in another condition. By
#'   only performing Welch's Anova on the "control" or base group it prevents this bias.
#'   There must be at least 10 samples in this list for comparision.
#' @param mtt_th
#'   The signficance threshold for performing the paired T-test for missingness. This test
#'   checks for signicant difference in missingness between the two genes of an edge.
#'   This is important because one gene with a high level of missginess will bias
#'   the relationship if that missigness is condition-specific. Only edges whose
#'   genes have the same pattern of missingness should be considered.  Those with
#'   p-values greater than the threshold are considered non-different.
#' @export
filterBiasedEdges <- function(net, ematrix, threads=0, wa_th = 1e-3, mtt_th = 0.1,
                              wa_samples = c()) {

  min_samples = max(c(30, min(net$Cluster_Size)))
  # Use all but two CPU cores. If there aren't more
  # than 2 then just default to using 1.
  if (threads == 0) {
    ncores = detectCores()
    threads = ncores - 2
    if (threads < 1) {
      threads = 1
    }
  }
  cl = makeCluster(threads)

  clusterExport(cl, c("ematrix", "wa_samples", "performBiasTests", "getEdgeSamples",
                      "getMissingEdgeSamples", "getNonEdgeSamples"))
  results = pbapply(net, 1, performBiasTests, min_samples=min_samples, cl=cl)
  stopCluster(cl)

  results = as.data.frame(t(results), stringsAsFactors=FALSE)

  # Fix the network dataframe to have the correct types.
  results$Similarity_Score = as.numeric(results$Similarity_Score)
  results$Cluster_Index = as.numeric(results$Cluster_Index)
  results$Cluster_Size = as.numeric(results$Cluster_Size)
  results$WAnova_Max = as.numeric(results$WAnova_Max)
  results$WAnova_Min = as.numeric(results$WAnova_Min)
  results$MissingTtest = as.numeric(results$MissingTtest)

  # Handle differences between tidy and not tidy networks.
  if ('Test_Name' %in% colnames(net)) {
    results$p_value = as.numeric(results$p_value)
    results$r_squared = as.numeric(results$r_squared)
  } else {
    for (i in 8:dim(net)[2]) {
      results[,i] = as.numeric(results[,i])
    }
  }

  # Keep any edges that pass the filters or have an NA for each test.  The
  # NA for a test indicates the test could not be run and therefore the
  # edge could not be excluded.
  keep = which((results$WAnova_Max <= wa_th | is.na(results$WAnova_Max)) &
               (results$MissingTtest >= mtt_th | is.na(results$MissingTtest)))

  return (results[keep,])
}

#' A helper function for the filterBiasedEdges.
#'
#' This function should not be called directly.
#'
#' @param row
#'   The row from the network dataframe.
#'
#' @export
performBiasTests <- function(row, min_samples=30) {
  gene1 = row['Source']
  gene2 = row['Target']
  row['WAnova_Max'] = NA
  row['WAnova_Min'] = NA
  row['MissingTtest'] = NA
  cluster_samples = getEdgeSamples(1, t(as.data.frame(row)))
  non_cluster_samples = getNonEdgeSamples(1, t(as.data.frame(row)))
  missing_samples = getMissingEdgeSamples(1, t(as.data.frame(row)))

  # Remove the missing samples from the non_cluster_samples
  non_cluster_samples = non_cluster_samples[which(!non_cluster_samples %in% missing_samples)]

  # If the user provided the list of control or base samples then
  # we want to restrict the non cluster samples to only that list.
  if (length(wa_samples) > 0) {
    non_cluster_samples = non_cluster_samples[which(non_cluster_samples %in% wa_samples)]
  }

  # First, perform a Welch's ANOVA test for both genes of every network edge.
  # The Welch's Anova test is used for heteroscedastic data (non-equal variance)
  # and does not assume that both categories have the same distribution.  This
  # test is used for each edge of the network to test if the mean of the
  # in-cluster and out-cluster samples of each gene are the same.

  # Only perform the test if we have at least 10 samples in and
  # out of the cluster.
  if (length(non_cluster_samples) >= 10 &
      length(cluster_samples) >= 10) {

    # Get the expression vectors of the samples not in the cluster
    # and remove outliers. We don't need to do this for the in-cluster
    # samples because clusters should be normal and have outliers removed
    # by KINC.
    ncs1 = as.vector(t(ematrix[gene1,]))
    ncs2 = as.vector(t(ematrix[gene2,]))
    ncs1[!seq(1:dim(ematrix)[2]) %in% non_cluster_samples] = NA
    ncs2[!seq(1:dim(ematrix)[2]) %in% non_cluster_samples] = NA
    b1 = boxplot(ncs1, plot = FALSE)
    b2 = boxplot(ncs2, plot = FALSE)
    ncs1[which(ncs1 %in% b1$out)] = NA
    ncs2[which(ncs2 %in% b2$out)] = NA

    # Get the expression of the genes for only the cluster samples
    # and the non_cluster samples.
    g1 = as.vector(t(ematrix[gene1,]))
    g2 = as.vector(t(ematrix[gene2,]))
    g1[!seq(1:dim(ematrix)[2]) %in% cluster_samples] = NA
    g2[!seq(1:dim(ematrix)[2]) %in% cluster_samples] = NA
    g1[non_cluster_samples] = ncs1[non_cluster_samples]
    g2[non_cluster_samples] = ncs2[non_cluster_samples]

    # We'll do bootstrapping so that we can subsample and
    # ensure comparable p-values.
    source_anova = 0
    target_anova = 0
    num_repeats = 15
    for (i in 1:num_repeats) {

      # Create an array of which samples are in/out of the cluster/
      # We will limit the sample to only 10 items so that p-values
      # between clusters are comparable.
      inout = rep(NA,dim(ematrix)[2])
      inout[sample(cluster_samples, min_samples, replace=TRUE)] = 'In'
      inout[sample(non_cluster_samples, min_samples, replace=TRUE)] = 'Out'

      # Perform the Welch one-way ANOVA test on in/out of gene1
      source = na.omit(data.frame(InOut = inout, Expression = g1))
      w1 = oneway.test(Expression ~ InOut, data=source, var.equal=FALSE)
      source_anova = source_anova + w1$p.value

      # Perform the Welch one-way ANOVA test on in/out of gene2
      target = na.omit(data.frame(InOut = inout, Expression = g2))
      w2 = oneway.test(Expression ~ InOut, data=target, var.equal=FALSE)
      target_anova = target_anova + w2$p.value
    }
    source_anova = source_anova / num_repeats
    target_anova = target_anova / num_repeats
    row['WAnova_Max'] = max(source_anova, target_anova)
    row['WAnova_Min'] = min(source_anova, target_anova)
  }

  # Let's do one more check to make sure there isn't bias in the missigness.
  # Convert the missigness of both genes to a feature where 1 == missing
  # and 0 == not missing
  g1m = as.vector(t(ematrix[gene1,]))
  g1m[which(!is.na(g1m))] = 0
  g1m[which(is.na(g1m))] = 1
  g2m = as.vector(t(ematrix[gene2,]))
  g2m[which(!is.na(g2m))] = 0
  g2m[which(is.na(g2m))] = 1

  # If the two vectors have identical missing then don't do the test.
  if (sum(abs(g1m - g2m)) != 0) {
    # Perform a paired t-test to see how similar the missingness is. If the
    # p-value is signficant then that means the pattern of missing
    # values is different and there may be bias.
    # We should not keep this edge.
    t = t.test(g1m, g2m, paired = TRUE, alternative = "two.sided")
    row['MissingTtest'] = t$p.value
  }

  # Perform multiple testing correction.
  #row['WAnova_Max'] = p.adjust(row['WAnova_Max'], 'BH')
  #row['WAnova_Min'] = p.adjust(row['WAnova_Min'], 'BH')
  #row['MissingTtest'] = p.adjust(row['MissingTtest'], 'BH')

  return(row)
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

#' Generates a matrix containing edge sample strings digits.
#'
#' The returned matrix contains as many rows as there are edges
#' in the network and as many columns as their are samples. Each
#' column corresponds to a sample in the same order as the
#' sample strings from the network.
#'
#' @param net
#'   A network data frame in Tidy format, containing the KINC-produced network.
#'   The loadNetwork function imports a dataframe in the correct format for this function.
#'
#' @return
#'   A list containing the following keys and values.
#'   - modules: a matrix with the list of modules in order of the edges
#'     present in the network provided by the net argument.
#'   - sampleMatrix: the sample matrix used to caculate the similarity matrix
#'     using the network sample strings.
#'   - dendro: An object of class hclust which describes the tree produced by the
#'     clustering process.
#'   - scores: a data frame containing the scores representing the strength
#'     of the association with the test names
#'
#' @export
getSSLinkModules = function(net, k=25, min_cluster_size = 10, show_plots = TRUE) {

  # Set p-values at 0 to a very low value so our
  # scoring works.
  net[which(net$p_value == 0), 'p_value'] = 1e-100

  # Create the dendrograms of sample strings in both strict and non-strict modes.
  dendro = clusterEdges(net, strict=TRUE)
  sampleMatrix = getSampleMatrix(net, strict=TRUE)

  # Get the list of modules from the tree.
  modules = cutree(dendro, k=k)

  # Create a dataframe for storing the scores.
  scores = data.frame(Module = as.character(), Size = as.numeric(),
                      Test_Name = as.character(), Score = as.numeric())

  # Create an empty array for storing the cluster membership.
  membership = rep('', dim(net)[1])

  # Iterate through the modules of the dendrogram.
  for (m in unique(modules)) {

    # Get the edges for this module and convert to iGraph object.
    m_indexes = (which(modules == m))

    mNet = netT[m_indexes,]
    g = KINCtoIgraph(mNet)
    g = set_edge_attr(g, 'Indexes', value=m_indexes)

    # skip modules smaller than the min size.
    min_verticies = 3
    if (gorder(g) <= min_verticies) {
      next
    }

    # Break the module into it's connected sub parts with
    # sub modules no smaller than 3 edges.
    subg = decompose(g,  min.vertices = min_verticies);

    # Iterate through the sub graphs and calculate a score for how well
    # overall represent the tests.
    if (length(subg) == 0) {
      subg[[1]] = g
    }

    for (gi in 1:length(subg)) {
      indexes = edge_attr(subg[[gi]], 'Indexes')
      size = length(indexes)
      subnet = netT[indexes, ]

      # Get some stats about the modules.
      netByMod = subnet %>% group_by(Test_Name)
      modScores = netByMod %>% summarise(p_value = mean(log10(p_value)) / log10(min(p_value)))
      modScores = as.data.frame(modScores)
      colnames(modScores) = c('Test_Name', 'Score')
      mod_name = paste0('SSM', sprintf('%04d', m), 'C', sprintf('%03d', gi))
      modScores['Module'] = mod_name
      modScores['Size'] = size
      if (size >= min_cluster_size) {
        membership[indexes] = mod_name
      }
      scores = rbind(scores, modScores)
    }
  }

  if (show_plots) {
    par(mfrow=c(2,2))
    hist(table(modules), xlab='Module Size', main=k)
    plot(density(table(modules)))
    boxplot(scores$Score, ylim=c(0,1), xlab="Score")
    plot(scores$Score, scores$Size, xlim=c(0,1))
  }

  return(list(
    'modules' = membership,
    'sampleMatrix' = sampleMatrix,
    'dendro' = dendro,
    'scores' = scores
  ))
}

#' A helper function for the filterPhasedEdges
#'
#' This function should not be called directly.
#'
#' @param row
#'   The row from the network dataframe.
#' @param test_column
#'   The name of the column in the osa matrix that should be tested for
#'   differential expression within the edge.
#' @param sample_col
#'   The name of the column in the osa matrix that contains the ssample name. The
#'   default is 'Sample'
#' @param min_cluster_size
#'   The minimum size of the cluster for each label in the test_column.  The
#'   default is 15.
#' @export
performHotellingTest = function(row, test_column, sample_col, min_cluster_size) {
  source = row["Source"]
  target = row["Target"]

  esi = getEdgeSamples(1, t(as.data.frame(row)))
  eamx = amx[esi,  c(sample_col, test_column)]
  eemx = t(emx[c(source,target),][esi])

  edata = merge(eamx,eemx,by="row.names")[,c(3,4,5)]
  colnames(edata) = c('Label', 'Gene1', 'Gene2')

  num_labels = table(edata$Label)
  test_labels_i = which(num_labels > min_cluster_size)
  if (length(test_labels_i) == 2) {
    test_labels = row.names(num_labels)[test_labels_i]
    edata = edata[which(edata$Label %in% test_labels),]
    #fit = hotelling.test(.~Label ,edata, perm=TRUE, B=30, progBar=FALSE)
    fit = hotelling.test(.~Label ,edata)
    return(fit$pval)
  }
  return (NA)
}

#' This function identifies "phased" edges of a network.
#'
#' A phased edge is one in which another categorical column other than the one
#' identified for the edge has differential expression within the GMM cluster
#' underlying the edge but did not have sufficient 'uniqueness" within the
#' cluster to identify a separate label-specific edge.
#'
#' This function currently only works when there are only two labels from the
#' test_column of min_cluster_size or greater in the edge. If fewer or more
#' are present then the edge is not tested.
#'
#' @param net
#'   A network dataframe containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param ematrix
#'   The expression matrix data frame as created by the loadGEM() function.
#' @param osa
#'   The sample annotation matrix as created by the loadSampleAnnotations()
#'   function.
#' @param test_column
#'   The name of the column in the osa matrix that should be tested for
#'   differential expression within the edge.
#' @param sample_col
#'   The name of the column in the osa matrix that contains the ssample name. The
#'   default is 'Sample'
#' @param min_cluster_size
#'   The minimum size of the cluster for each label in the test_column.  The
#'   default is 15.
#' @param threads
#'   The number of computational threads to use for parallel processing. By
#'   default, all but 2 cores available on the local machine will be used.
#'
#' @return
#'   The original network dataframes with an extra column added
#'   containing the p-value of the test.
#'
#' @export
performEdgeDiff = function(net, emx, amx, test_column, sample_col="Sample",
                             min_cluster_size = 15, threads = 0 ) {

  # Use all but two CPU cores. If there aren't more
  # than 2 then just default to using 1.
  if (threads == 0) {
    ncores = detectCores()
    threads = ncores - 2
    if (threads < 1) {
      threads = 1
    }
  }

  cl = makeCluster(threads)
  clusterExport(cl, c("net", "emx", "amx", "getEdgeSamples", "performHotellingTest", "hotelling.test"))
  results = pbapply(net, 1, performHotellingTest,  test_column, sample_col, min_cluster_size, cl=cl)
  stopCluster(cl)

  net['hotelling_p_value'] = as.numeric(results)

  return (net)

}
