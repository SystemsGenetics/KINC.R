
#' Perfoms heirarchical clustering of edges in a network based on their sample compositions.
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

#' Performs a fisher's exact test on a sample string.
#'
#' The sample annotation matrix consists of multiple columns of annotation
#' types (or fields).  Each annotation type may have multiple distinct
#' values. For example the type may be 'sex' and the categories may be
#' 'male' and 'female'.  This function performs a Fisher's test on
#' a single category amongst all the other categories within the annotation
#' type.
#'
#' @param category
#'   The annotation category.
#' @param sample_types
#'   The ordered array of annotations for all samples (i.e. a full column
#'   of the sample annotation matrix).
#' @param rep_samples
#'   An array of numeric values indicating the strength of the sample in
#'   the edge(s).
#' @param correction.
#'   The method to apply for multiple testing correction. Valid values are identical
#'   to those available to the p.adjust function.  The default is to
#'   apply no correction.
#' @param alternative
#'   indicates the alternative hypothesis and must be one of "two.sided",
#'   "greater" or "less". You can specify just the initial letter.
#'
fishers_test = function(category, sample_types, sample_ref, correction = NA, alternative = 'greater') {
  #  Contingency matrix for each category  in a module:
  #
  #                   Is Type  Not Type Totals
  #                  ------------------
  #  Module          |  n11   |   n12   | n1p
  #  Not in Mod      |  n21   |   n22   | n2p
  #                  ------------------
  #  Totals             np1       np2     npp
  #
  n11 = round(sum(sample_ref[which(sample_types == category)]))
  n12 = round(sum(sample_ref[which(sample_types != category)]))
  np1 = length(which(sample_types == category))
  np2 = length(which(sample_types != category))
  n21 = np1 - n11
  n22 = np2 - n12

  contmatrix = matrix(
    as.numeric(c(n11, n21, n12, n22)),
    nr=2,
    dimnames = list(
      In_Module = c("Yes", "No"),
      Is_Type = c("Yes", "No")
    )
  )
  #print(category)
  #print(contmatrix)
  res = fisher.test(contmatrix, alternative = alternative)
  p.value = res$p.value
  if (!is.na(correction)) {
    p.value = p.adjust(p.value, method=correction, n=length(unique(sample_types)))
  }
  return(p.value)
} # end fisher’s test function

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

#' Performs enrichment analysis of traits against a a single edge in the network.
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
#' @param correction.
#'   The method to apply for multiple testing correction. Valid values are identical
#'   to those available to the p.adjust function.  The default is to
#'   apply 'bonferroni' correction.
#' @export
#'
#' @examples
#'
analyzeEdgeCat = function(i, osa, net, field, correction = 'bonferroni') {

  sample_types = as.character(osa[[field]])
  num_samples = length(sample_types)

  sample_types.prob = table(sample_types) / num_samples
  nona.idx = which(!is.na(sample_types))
  categories = sort(unique(sample_types[nona.idx]))

  # Convert the sample string into an integer array.
  edge = net[i,]
  sample_str = edge$Samples
  samples = as.numeric(strsplit(sample_str, "")[[1]])
  samples[which(samples != 1)] = 0

  pvals = sapply(categories, fishers_test, sample_types, samples, correction)
  return(pvals);
}

#' Performs linear regression of a quantitative traits against a a single edge in the network.
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
#' @param correction.
#'   The method to apply for multiple testing correction. Valid values are identical
#'   to those available to the p.adjust function.  The default is to
#'   apply 'bonferroni' correction.
#' @param samples
#'   Limit the annalysis to only the samples indexes provided.
#'
#' @export
#'
#' @examples
#'
analyzeEdgeQuant = function(i, osa, net, field, samples = c()) {
  source = net[i, 'Source']
  target = net[i, 'Target']

  # Convert the sample string to a numerical vector.
  edge = net[i,]
  sample_str = edge$Samples
  edge_samples = as.numeric(strsplit(sample_str, "")[[1]])
  num_samples = length(edge_samples)

  # If the caller provided a set of samples to focus on then
  # we want to ignore any that are in the edge but not in the list.
  if (length(samples) > 0) {
    edge_samples[which(!(1:num_samples %in% samples))] = 0
  }

  # Set any value not equal to 1 to be zero.
  edge_samples[which(edge_samples != 1)] = 0

  # Convert the edge_sample list back to a listo f indexes.
  edge_samples = which(edge_samples == 1)

  # Build the expression vectors using only the edges of the edge.
  x = t(ematrix[source, edge_samples])
  y = t(ematrix[target, edge_samples])
  z = as.factor(osa[edge_samples, field])

  model = lm(y + x ~ as.numeric(z), data=data.frame(x=x, y=y, z=z))
  s = summary(model)
  pval = s$coefficients[2,4]
  return(pval)
}

#' Performs significance association of a subnetwork with a field.
#'
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param edge_indexes
#'   The index of the edges in the network that comprise the module.
#' @param osa
#'   The ordered sample annotation data frame.
#' @param field
#'   The field in the osa variable on which enrichment will be performed.
#' @param correction.
#'   The method to apply for multiple testing correction. Valid values are identical
#'   to those available to the p.adjust function.  The default is to
#'   apply 'bonferroni' correction.
#' @param alternative
#'   indicates the alternative hypothesis and must be one of "two.sided",
#'   "greater" or "less". You can specify just the initial letter.
#' @export
#'
analyzeEdgeListCat = function(edge_indexes, osa, net, field, correction = 'bonferroni', alternative = 'greater') {

  sample_types = as.character(osa[[field]])
  sample_types[which(sample_types == "null")] = NA
  sample_types[which(sample_types == "notreported")] = NA
  num_samples = length(sample_types)

  sample_types.prob = table(sample_types) / num_samples
  nona.idx = which(!is.na(sample_types))
  categories = sort(unique(sample_types[nona.idx]))

  # Calculate the relative frequency of each sample in the module.
  mod_edges = net[edge_indexes,]
  num_edges = nrow(mod_edges)
  sample_strs = mod_edges$Samples
  mod_ref = vector(mode="numeric", length=num_samples)
  for (sample_str in sample_strs) {
    samples = as.numeric(strsplit(sample_str, "")[[1]])
    samples[which(samples != 1)] = 0
    mod_ref = mod_ref + samples
  }
  mod_ref = mod_ref / num_edges

  results = sapply(categories, fishers_test, sample_types, mod_ref, correction, alternative)
  names(results) = categories
  return(results)
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
#'   apply 'bonferroni' correction.
#' @export
#'
analyzeNetCat = function(net, osa, field, correction = 'bonferroni') {

  sample_types = as.character(osa[[field]])

  categories = unique(sample_types)

  # Add in new columns for each category.
  net2 = net
  for (category in categories) {
    subname = paste(field, category, sep='_')
    net2[subname] = NA
  }

  pb <- txtProgressBar(min = 0, max = nrow(net), style = 3)
  for (i in 1:nrow(net)) {
    setTxtProgressBar(pb, i)
    p.vals = analyzeEdgeCat(i, osa, net, field, correction)
    for (category in names(p.vals)) {
      subname = paste(field, category, sep='_')
      net2[i, subname] = p.vals[category]
    }
  }
  close(pb)
  return(net2);
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
#' @param samples
#'   Limit the annalysis to only the samples indexes provided.
#'
#' @export
#'
analyzeNetQuant = function(net, osa, field, samples = NA) {

  # Add in new columns for each category.
  net2 = net
  net2[field] = NA

  pb <- txtProgressBar(min = 0, max = nrow(net), style = 3)
  for (i in 1:nrow(net)) {
    setTxtProgressBar(pb, i)
    p.val = analyzeEdgeQuant(i, osa, net, field, samples)
    net2[i, field] = p.val
  }
  close(pb)
  return(net2);
}
