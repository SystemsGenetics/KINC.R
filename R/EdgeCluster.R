#' Imports sample annotations.
#'
#' The format of the sample annotation file is tab delimited where the
#' first column is the sample name and every other column contains
#' the annotation.  The file must contain a header with the annotation
#' types.
#'
#' @param annotation_file
#'   The path to the file containing the annotations.
#' @param smaple_order
#'   A file containing the order of samples in the sample strings of
#'   the network.  The file should contain the list of samples each on
#'   a separate line.
#'
#' @return
#'   A data frame containing the annotations in the order of the samples.
#'
#' @export
loadSampleAnnotations = function (annotation_file, sample_order) {
  sample_annots = read.table(annotation_file, sep="\t", header=TRUE, row.names=NULL)
  sample_order = read.table(sample_order, colClasses=c('character'),
    col.names=c('Sample'))
  osa = merge(sample_order, sample_annots, by = "Sample", sort=FALSE)

  return(osa)
}

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
#' @export
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
#' @export
drawNetHeatMap = function(sampleMatrix, tree, osa, fieldOrder, outfile_prefix = NA) {

  # Reorder samples in the sampleMatrix according to the fieldOrder argument.
  sample_order = eval(parse(text=paste('order(osa$', paste(fieldOrder, collapse=' ,osa$'), ")", sep="")))
  sampleMatrix2 = sampleMatrix[, c(sample_order)]

  # Get the members of each cluster.
  members = cutree(tree, k = num_clusters)

  #
  categories = do.call(paste, list(c(osa[, fieldOrder]), sep="-"))
  num_categories = length(unique(categories))
  osa$hmap_categories = categories
  tColors = data.frame(
    Field = unique(categories),
    Color = rgb(runif(num_categories), runif(num_categories), runif(num_categories))
  )
  mColors = data.frame(
    Cluster = unique(members),
    color = rgb(runif(num_clusters), runif(num_clusters), runif(num_clusters))
  )
  colColors = as.character(merge(osa[sample_order,], tColors, by.x="hmap_categories", by.y="Field", sort=FALSE)$Color)
  rowColors = as.character(factor(members, labels=mColors$color))

  if (!is.na(outfile_prefix)) {
    outfile = paste(outfile_prefix, paste(fieldOrder, collapse="-"), num_clusters, "png", sep=".")
    png(filename = outfile, width=3000, height=21000, res=300)
  }
  heatmap.2(sampleMatrix2,
    Rowv=as.dendrogram(tree),
    Colv=FALSE,
    dendrogram = 'row',
    col = c("green", "red"),
    breaks = c(-1, 0, 1),
    trace = 'none',
    key = TRUE,
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
drawModuleHeatMap = function(edge_indexes, osa, net, fieldOrder) {
  ce.tree = clusterEdges(net[edge_indexes,])
  sm = getSampleMatrix(net[edge_indexes,])
  drawNetHeatMap(sm, ce.tree, osa, fieldOrder)
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
analyzeEdge = function(i, osa, net, field, correction = 'bonferroni') {

  sample_types = as.character(osa[[field]])
  sample_types[which(sample_types == "null")] = NA
  sample_types[which(sample_types == "notreported")] = NA
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
analyzeModule = function(edge_indexes, osa, net, field, correction = 'bonferroni', alternative = 'greater') {

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
analyzeNet = function(net, osa, field, correction = 'bonferroni') {

  sample_types = as.character(osa[[field]])
  sample_types[which(sample_types == "null")] = NA
  sample_types[which(sample_types == "notreported")] = NA

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
    p.vals = analyzeEdge(i, osa, net, field, correction)
    for (category in names(p.vals)) {
      subname = paste(field, category, sep='_')
      net2[i, subname] = p.vals[category]
    }
  }
  close(pb)
  return(net2);
}

#' Performs enrichment analysis of traits against a network dendrogram
#'
#' Iterates through the tree created by the clusterEdges() function and
#' performs a Fisher's test on each of the
#'
#' @param tree
#'   An instance of an hclust object.
#' @param osa
#'   The sample annotation matrix. One column must contain the header 'Sample'
#'   and the remaining colums correspond to an annotation type.  The rows
#'   of the anntation columns should contain the annotations.
#' @param net
#'   A network dataframe containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param field
#'   The field in the osa variable on which enrichment will be performed.
#' @param min_cluster_size
#'   The minimum number of network edges that must be present in a cluster
#'   In order for that cluster to be analyzed.
#'
#' @export
#'
#' @examples
#'
analyzeEdgeTree = function(tree, osa, net, field, min_cluster_size = 3) {

  height_array = unique(tree$height)
  height_array = sort(height_array, decreasing = TRUE)

  # Initialize seen_checksums empty vector
  seen_checksums <- vector(mode="character", length=0)
  prev_indexes = list()

  # Create the list that will hold the results.
  final_results = list();

  pb <- txtProgressBar(min = 0, max = length(height_array), style = 3)

  # Iterate through the dendrogram at increasing heights and perform
  # enrichment analysis
  depth = 1
  for (height in height_array) {
    setTxtProgressBar(pb, depth)
    members = cutree(tree, h = height)
    num_clusters = length(unique(members))
    #drawDendro(osa, members, samples, depth, 'Treatment')

    prev_indexes[[depth]] = list()
    for(i in unique(members)) {
      # print(paste("Depth: ", depth, ", Num Clusters: ", num_clusters, ", Cluster: ", i, sep=""))

      cluster_indexes = as.integer(which(members == i))
      clust_chksm = digest(cluster_indexes)
      prev_indexes[[depth]][[clust_chksm]] = cluster_indexes

      # if the checksum of the cluster index array is not in seen_checksums
      if (!is.element(clust_chksm, seen_checksums)) {

        # Add it to seen checksums
        seen_checksums = c(seen_checksums, clust_chksm)

        # Find the parent
        parent = NA
        if (depth > 1) {
          for (pchksm in names(prev_indexes[[depth-1]])) {
            if (all(cluster_indexes %in% prev_indexes[[depth-1]][[pchksm]])) {
              parent = pchksm
            }
          }
        }

        # Calculate the average degree
        cluster_nodes = unique(c(
          net[which(members == i), c('Source')],
          net[which(members == i), c('Target')]
        ))
        avg_degree = (length(which(members == i)) * 2) / length(cluster_nodes)

        # Run enrichment if the cluster size meets the minimum size limit.
        if (length(cluster_indexes) >= min_cluster_size) {
          enrichment = analyzeModule(cluster_indexes, osa, net, field)
          n = length(final_results)
          final_results[[n+1]] = list(
            "id" = clust_chksm,
            "parent" = parent,
            "depth" = depth,
            "height" = height,
            "num_clusters" = num_clusters,
            "cluster" = i,
            "cluster_size" = length(cluster_indexes),
            "avg_degree" = avg_degree,
            "field" = field,
            "categories" = unique(osa[,c(field)]),
            "p_vals" = enrichment,
            "nodes" = cluster_nodes,
            "edge_indexes" = cluster_indexes
          );
        } # end if (length(cluster_indexes
      } # end if (!is.element(clust_chksm,
    } # end for(i in unique(members)) {
    depth = depth + 1
  } # end for (height in height_array) {
  close(pb)
  return(final_results)
}

#' Writes to a file the results from the analyzeEdgeTree function.
#'
#' @param results
#'   A list containing the results from the analyzeEdgeTree function.
#' @param outfile
#'   The name of the file where output results are stored.
#' @param field
#'   The field in the osa variable on which enrichment will be performed.
#' @export
#'
writeEdgeTreeAnalysis = function(results, outfile = 'output.scc.txt') {
  output = file(outfile, "w")
  write(paste("Cluster_ID", "Parent_ID", "Depth", "Height", "Clusters", "Cluster_Num", "Size",
              "Avg Degree", "Type", "Categories", "Enrichment", "Nodes", "Edge Index",
              sep="\t"), file=output, append=FALSE)
  pb <- txtProgressBar(min = 0, max = length(results), style = 3)

  for(i in 1:length(results)) {
    setTxtProgressBar(pb, i)
    row = results[[i]]
    write(
      paste(
        row$id,
        row$parent,
        row$depth,
        row$height,
        row$num_clusters,
        row$cluster,
        row$cluster_size,
        row$avg_degree,
        row$field,
        paste(row$categories, collapse=","),
        paste(row$"p-pvals", collapse=","),
        paste(row$nodes, collapse=","),
        paste(row$edge_indexes, collapse=","),
        sep="\t"),
      file = output,
      append = TRUE
    )
  }
  close(output)
  close(pb)
}

#' Finds the best subgraph for the given field and category.
#'
#' Searches through the results from the analyzeEdgeTree for the cluster whose
#' p-value is best given a specific field and category.
#'
#' @param results
#'   A list containing the results from the analyzeEdgeTree function.
#' @param field
#'   One of the header names from the sample annnotation matrix.
#' @param category
#'   For categorical data in the sample annotation matrix, this is one of the
#'   categories found in the field column.
#' @param alpha
#'   The alpha value used for significance.
#'
#' @return
#'   A dataframe containing the set of subgraphs (heirarchical clusters)
#'   that have a p-value less than alpha.
#'
#' @export
findSubgraphs = function(results, field, category, alpha = 0.001) {
  # The best matching row index and p-value in the results will be stored here.
  best = list()

  for(i in 1:length(results)) {
    row = results[[i]]
    if (row$field == field) {
      cindex = which(row$categories == category)
      if (cindex == 0) {
        next;
      }
      sig.pvals = which(row$p_vals < 0.001)
      if (length(sig.pvals) == 1 && cindex %in% sig.pvals) {
        p = row$p_vals[cindex]
        best[[length(best) + 1]] = list(
          index = i,
          p = p,
          size = row$cluster_size,
          avg_degree = row$avg_degree
        )
      }
    }
  }
  matches = data.frame(matrix(unlist(best), nrow=length(best), ncol=4, byrow=T))
  names(matches) = c('index', 'p' ,'size', 'avg_degree')
  return(matches)
}
