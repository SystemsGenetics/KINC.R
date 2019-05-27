#' Converts an edge sample string into an array of indexes of samples belonging to the cluseter.
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

#' Performs a fisher's exact test on a set of samples.
#'
#' An annotation field and a specific category are provided. This
#' test reports if the category is enriched within the cluster.
#' It does not indiciate if the category is singificantly more
#' prominent in the cluster. For that, try the sampleClusterBTest
#' function.
#' @param category
#'   The annotation category to be used for testing. It must be a valid
#'   category in the field if the osa identified by the field argument.
#' @param field
#'   The field in the osa variable on which enrichment will be performed.
#' @param i
#'   The index of the edge in the network
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param osa
#'   The sample annotation matrix. One column must contain the header 'Sample'
#'   and the remaining colums correspond to an annotation type.  The rows
#'   of the anntation columns should contain the annotations.
#' @param ematrix
#'   The expression matrix data frame.
#' @param samples
#'   Limit the analysis to only the samples indexes provided.
performClusterFishersTest = function(category, field, i, net, osa, ematrix, verbose=FALSE) {

  # Convert the sample string to a numerical vector.
  cluster_samples = getEdgeSamples(i, net)

  # Get osa elements for samples that belong to the cluster and those that don't.
  osa_cat_indexes = which(osa$Sample %in% names(ematrix)[cluster_samples])
  osa_out_indexes = which(!osa$Sample %in% names(ematrix)[cluster_samples])

  #  Contingency matrix for each category in a cluster:
  #
  #                   Is Cat  Not Cat    Totals
  #                  ------------------
  #  In Cluster      |  n11   |   n12   | n1p
  #  Not in Cluster  |  n21   |   n22   | n2p
  #                  ------------------
  #  Totals             np1       np2     npp
  #
  n11 = length(which(osa[[field]][osa_cat_indexes] == category))
  n12 = length(which(osa[[field]][osa_cat_indexes] != category))
  n21 = length(which(osa[[field]][osa_out_indexes] == category))
  n22 = length(which(osa[[field]][osa_out_indexes] != category))

  contmatrix = matrix(
    as.numeric(c(n11, n21, n12, n22)),
    nr=2,
    dimnames = list(
      In_Cluster = c("Yes", "No"),
      Is_Category = c("Yes", "No")
    )
  )

  if (verbose) {
    print(contmatrix)
  }
  res = fisher.test(contmatrix, alternative = "greater")
  p.value = res$p.value
  return(p.value)
}

#' Performs a fisher's exact test on a set of samples.
#'
#' An annotation field and a specific category are provided. This
#' test reports if the category is enriched within the cluster.
#' It does not indiciate if the category is singificantly more
#' prominent in the cluster. For that, try the sampleClusterBTest
#' function.
#' @param category
#'   The annotation category to be used for testing. It must be a valid
#'   category in the field if the osa identified by the field argument.
#' @param field
#'   The field in the osa variable on which enrichment will be performed.
#' @param i
#'   The index of the edge in the network
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param osa
#'   The sample annotation matrix. One column must contain the header 'Sample'
#'   and the remaining colums correspond to an annotation type.  The rows
#'   of the anntation columns should contain the annotations.
#' @param ematrix
#'   The expression matrix data frame.
#' @param t1_psucc
#'   The probabilty of success for the first binomial test.  This value
#'   indicates the percentage of samples in the cluster that must
#'   belong to the category.
#' @param t2_psucc
#'   The probability of success for the second binomial test. This value
#'   indicates the percentage of category samples that must belong
#'   to the cluster.
performClusterBinomialTest = function(category, field, i, net, osa, ematrix,
                                      t1_psucc = 0.25, t2_psucc = 0.75, verbose=FALSE) {

  # Convert the sample string to a numerical vector.
  cluster_samples = getEdgeSamples(i, net)

  # Get the cluster size
  cluster_size = length(cluster_samples)

  # How many samples belong to this category
  num_category = length(which(osa[[field]] == category))
  num_not_category = length(which(osa[[field]] != category))

  # Get the indexes in the OSA of the samples that are in the cluster and not in
  # the cluster.
  cluster_osa_indexes = which(osa$Sample %in% names(ematrix)[cluster_samples])
  non_cluster_osa_indexes = which(!osa$Sample %in% names(ematrix)[cluster_samples])

  # Get the number of samples of the category that are in and not in the cluster.
  num_category_in_cluster = length(which(osa[[field]][cluster_osa_indexes] == category))
  num_category_not_in_cluster = num_category - num_category_in_cluster

  # Get the number of samples not of the category that are in and not in the cluster.
  num_other_in_cluster = length(which(osa[[field]][cluster_osa_indexes] != category))
  num_other_not_in_cluster = num_not_category - num_other_in_cluster

  # Test #1
  # successes = number of non-category samples in the cluster
  # failures = number of non-category samples not in the cluster
  # Ho: successes >= 0.15
  # Ha: successes < 0.15
  p1.pval = 0
  if (num_not_category > 0) {
    p1 = binom.test(num_other_in_cluster, num_not_category, p=t1_psucc, alternative='less')
    p1.pval = p1$p.value
  }

  # Test #1b
  # Perform the same test as #1 but for each other category separately.

  # Test #2
  # successes = number of category samples in the cluster
  # failures = number of category samples out of the cluster
  # Ho: successes = 0.85
  # Ha: successes > 0.85
  p2.pval = 0
  if (num_category > 0) {
    p2 = binom.test(num_category_in_cluster, num_category, p=t2_psucc, alternative='greater')
    p2.pval = p2$p.value
  }

  # Return the maximum p-value
  # TODO: perhaps we should employ a method of combining p-values rather
  # than just using the minimum.
  # https://academic.oup.com/bioinformatics/article/32/17/i430/2450768
  p.value = max(p1.pval, p2.pval)

  if (verbose) {
    print(p1)
    print(p2)
  }

  return(p.value)
}

#' Performs significant testing of a single edge in the network for a set of annotation categories.
#'
#' @param i
#'   The index of the edge in the network
#' @param osa
#'   The sample annotation matrix. One column must contain the header 'Sample'
#'   and the remaining colums correspond to an annotation type.  The rows
#'   of the anntation columns should contain the annotations.
#' @param net
#'   A network dataframe containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param ematrix
#'   The expression matrix data frame.
#' @param field
#'   The field in the osa variable on which enrichment will be performed.
#' @param test
#'   The statistical test to perform: 'fishers' or 'binomial'
#' @param category
#'   By default this function will cylce through all categories of the
#'   'field' argument. But to limit the analysis to just one category set it here.
#' @param t1_psucc
#'   For the bionomial test only.
#'   The probabilty of success for the first binomial test.  This value
#'   indicates the percentage of samples in the cluster that must
#'   belong to the category.
#' @param t2_psucc
#'   For the binomial test only.
#'   The probability of success for the second binomial test. This value
#'   indicates the percentage of category samples that must belong
#'   to the cluster.
#' @export
#'
#' @examples
#'
analyzeEdgeCat = function(i, osa, net, ematrix, field, test = 'binomial',
                          category = NA, t1_psucc = 0.25, t2_psucc = 0.75) {

  sample_types = as.character(osa[[field]])
  num_samples = length(sample_types)

  categories = unique(osa[[field]])
  if (!is.na(category)) {
    categories = c(category)
  }

  if (test == 'fishers') {
    pvals = sapply(categories, performClusterFishersTest, field, i, net, osa, ematrix, FALSE)
    names(pvals) = categories
    return(pvals);
  }
  if (test == 'binomial') {
    pvals = sapply(categories, performClusterBinomialTest, field, i, net, osa, ematrix, t1_psucc, t2_psucc, FALSE)
    names(pvals) = categories
    return(pvals);
  }

}

#' Performs significant testing of each edge in the network for a set of annotation categories.
#'
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param osa
#'   The ordered sample annotation data frame.
#' @param ematrix
#'   The expression matrix data frame.
#' @param field
#'   The field in the osa variable on which enrichment will be performed.
#' @param test
#'   The statistical test to perform.  This is 'fishers' for an enrichment
#'   test.
#' @param correction.
#'   The method to apply for multiple testing correction. Valid values are identical
#'   to those available to the p.adjust function.  The default is to
#'   apply 'hochberg' correction.
#' @param progressBar
#'   Set to FALSE to repress progress bar
#' @param category
#'   By default this function will cylce through all categories of the
#'   'field' argument. But to limit the analysis to just one category set it here.
#' @param t1_psucc
#'   For the bionomial test only.
#'   The probabilty of success for the first binomial test.  This value
#'   indicates the percentage of samples in the cluster that must
#'   belong to the category.
#' @param t2_psucc
#'   For the binomial test only.
#'   The probability of success for the second binomial test. This value
#'   indicates the percentage of category samples that must belong
#'   to the cluster.
#' @export
#'
analyzeNetCat = function(net, osa, ematrix, field, test = 'binomial',
                         correction = 'hochberg', progressBar = TRUE,
                         category = NA, t1_psucc = 0.25, t2_psucc=0.75) {

  # Get the list of categories to analyze.
  sample_types = as.character(osa[[field]])
  categories = unique(sample_types)
  if (!is.na(category)) {
    categories = c(category)
  }

  # Add in new columns for each category.
  net2 = net
  for (label in categories) {
    subname = paste(field, label, sep='_')
    net2[subname] = NA
  }

  # Intialize the progress bar
  if (progressBar){
    pb <- txtProgressBar(min = 0, max = nrow(net), style = 3)
  }

  # Iterate through the rows of the network
  for (i in 1:nrow(net)) {
    if (progressBar){
      setTxtProgressBar(pb, i)
    }

    # Calculate the probability that this edge is a result of any specific
    # known experimental condition (i.e. category) for the given field.
    p.vals = analyzeEdgeCat(i, osa, net, ematrix, field, test = test,
                            category, t1_psucc, t2_psucc)
    for (label in names(p.vals)) {
      subname = paste(field, label, sep='_')
      net2[i, subname] = p.vals[label]
    }
  }
  if (progressBar){
    close(pb)
  }

  # Perform multiple testing correction on the p-values. This may not
  # be necessary because it overly reduces p-values from the binomial
  # test, but for backwards compatiblity, let's leave it for fishers:
  if (test == 'fishers') {
    for (category in categories) {
      subname = paste(field, category, sep='_')
      net2[subname] = p.adjust(as.numeric(unlist(net2[subname])), method=correction)
    }
  }

  return(net2);
}
#' Performs linear regression of a quantitative traits against a single edge in the network.
#'
#' @param i
#'   The index of the edge in the network
#' @param osa
#'   The sample annotation matrix. One column must contain the header 'Sample'
#'   and the remaining colums correspond to an annotation type.  The rows
#'   of the anntation columns should contain the annotations.
#' @param net
#'   A network dataframe containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param field
#'   The field in the osa variable on which enrichment will be performed.
#' @param samples
#'   Limit the analysis to only the samples indexes provided.
#'
#' @export
#'
#' @examples
#'
analyzeEdgeQuant = function(i, osa, net, field, samples = c()) {
  source = net[i, 'Source']
  target = net[i, 'Target']

  # Convert the sample string to a numerical vector.
  edge_samples = getEdgeSamples(i, net)
  if (length(samples) > 0) {
    edge_samples = samples
  }

  # Build the expression vectors using only the samples provided.
  x = t(ematrix[source, edge_samples])
  y = t(ematrix[target, edge_samples])
  z = as.numeric(as.factor(osa[edge_samples, field]))

  # If our z-dimension only has one value then there's no need to
  # perform a linear regression, so just return.
  if (dim(table(z)) == 1) {
    return(list(
      p = NA,
      roccm = NA,
      model = NA
    ))
  }

  # Use linear regression to obtain a p-value for the association.
  model = lm(y + x ~ z, data=data.frame(x=x, y=y, z=z))
  s = summary(model)
  pval = NA
  roccm = NA
  if (dim(s$coefficients)[1] > 1) {
    pval = s$coefficients[2,4]
    # Get the rate of change of the mean of y + x
    roccm = s$coefficients[2,1]
  }

  return(list(
    p = pval,
    roccm = roccm,
    model = model
  ))
}

#' Performs significance association for every edge in the network.
#'
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param osa
#'   The ordered sample annotation data frame.
#' @param field
#'   The field in the osa variable on which enrichment will be performed.
#' @param correction.
#'   The method to apply for multiple testing correction. Valid values are identical
#'   to those available to the p.adjust function.  The default is to
#'   apply 'hochberg' correction.
#' @param samples
#'   Limit the annalysis to only the samples indexes provided.
#' @param progressBar
#'   Set to FALSE to repress progress bar
#'
#' @export
#'
analyzeNetQuant = function(net, osa, field, correction = 'hochberg',
                           samples = c(), progressBar = TRUE) {

  # Add in new columns for each category.
  net2 = net
  net2[field] = NA
  net2[paste(field, '_roccm', sep='')] = NA

  if(progressBar){pb <- txtProgressBar(min = 0, max = nrow(net), style = 3)}
  for (i in 1:nrow(net)) {
    if (progressBar){setTxtProgressBar(pb, i)}
    result = analyzeEdgeQuant(i, osa, net, field, samples)
    net2[i, field] = result[['p']]
    net2[i, paste(field, '_roccm', sep='')] = result[['roccm']]
  }
  if (progressBar){close(pb)}

  # Perform multiple testing correction on the p-values
  net2[field] = p.adjust(as.numeric(unlist(net2[field])), method=correction)

  return(net2);
}

analyzeNetDiffQuant = function(net, osa, field, model_samples = NA, test_samples,
                               progressBar = FALSE) {

  # Add in new columns for each category.
  net2 = net
  column_name = paste(field, 'diff', paste(test_samples, collapse='_'), sep="-")
  net2[column_name] = NA

  if(progressBar){pb <- txtProgressBar(min = 0, max = nrow(net), style = 3)}
  for (i in 1:nrow(net)) {
    if (progressBar){setTxtProgressBar(pb, i)}
    result = analyzeEdgeDiffQuant(i, osa, net, field, model_samples, test_samples)
    net2[i, column_name] = median(result)
  }
  if (progressBar){close(pb)}

  return(net2);
}

analyzeEdgeDiffQuant = function(i, osa, net, field, model_samples = c(), test_samples = c()) {
  source = net[i, 'Source']
  target = net[i, 'Target']

  # Convert the sample string to a numerical vector.
  edge = net[i,]
  edge_samples = getEdgeSamples(edge$Samples)
  if (length(samples) > 0) {
    edge_samples = samples
  }

  # Build the expression vectors using only the edges of the edge.
  x = t(ematrix[source, ])
  y = t(ematrix[target, ])
  z = as.numeric(as.factor(osa[, field]))

  # Use linear regression to obtain a model for the samples the model samples group.
  model = lm(y + x ~ z, data=data.frame(x=x, y=y, z=z), subset = edge_samples)

  # Now get the difference for the samples in the test_sample group
  obs = t(ematrix[source, test_samples]) + t(ematrix[target, test_samples])
  colnames(obs) = 'y + x'
  exp = predict(model, data.frame(z = z[test_samples]), interval = 'predict')
  res = (obs - exp[, 'lwr']) / (exp[, 'upr'] - exp[, 'lwr']) * 2 - 1
  return (res)
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
