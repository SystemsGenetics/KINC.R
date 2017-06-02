library(gplots)
library(plotrix)
library(digest)
 
 
fishers_test = function(category, sample_types, rep_samples) {
  #  Contingency matrix for each category  in a module:
  #       
  #                   Is Type  Not Type Totals
  #                  ------------------
  #  Module              |  n11   |   n12   | n1p
  #  Not in Mod        |  n21   |   n22   | n2p
  #                  ------------------
  #  Totals               np1       np2       npp
  #
  n11 = length(which(sample_types[rep_samples] == category))
  n1p = length(which(!is.na(sample_types[rep_samples])))
  np1 = length(which(sample_types == category))
  np2 = length(which(sample_types != category & !is.na(sample_types)))
  n21 = np1 - n11
  n12 = n1p - n11
  n22 = np2 - n12
 
  contmatrix = matrix(
    as.numeric(c(n11, n12, n21, n22)),
    nr=2,
    dimnames = list(
      Is_Type = c("Yes", "No"),
      In_Module = c("Yes", "No")
    )
  )
  #print(contmatrix)
  res = fisher.test(contmatrix, alternative="greater")
  return(res$p.value)
} # end fisher’s test function
 
# Params
# @net:  the network
# @edge_indexes:  the index of the edges in the network that comprise the module.
# @osa: the ordered sample annotation data frame
# @field: the field in the osa variable on which enrichment will be performed.
# @min_presence: the percentage of edges in the module for which a sample
#   must be present in order to be counted.
test_module = function(net, edge_indexes, osa, field, min_presence = 0.95) {
 
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
 
  # Get the list of samples that are present > min_presence% in the module.
  rep_samples = which(mod_ref >= min_presence)
  num_rep_samples = length(rep_samples)
  if (num_rep_samples == 0) {
    results = rep(NA, length(categories))
    names(results) = categories
    return(results)
  }
 
  # get the frequencies of each category.
  rep_sample_types = as.character(na.omit(sample_types[rep_samples]))
  mod.freq = table(rep_sample_types)
  #print(mod.freq)  
 
  results = sapply(categories, fishers_test, sample_types, rep_samples)
  names(results) = categories
  return(list(results, table(rep_sample_types)))
}
 
drawNetDendro = function(osa, tree, height, samples, depth, field, tOrder) {
 
  members = cutree(tree, h = height)
 
  num_categories = length(unique(osa[, c(field)]))
  num_clusters = length(unique(members))
  #tOrder = order(osa[, c(field)])
  samples2 = samples[,tOrder]
  tColors = data.frame(
    Field = unique(osa[, c(field)]), 
    Color = rgb(runif(num_categories), runif(num_categories), runif(num_categories))
  )
  mColors = data.frame(
    Cluster = unique(members),
    color = rgb(runif(num_clusters), runif(num_clusters), runif(num_clusters))
  )
  colColors = as.character(merge(osa, tColors, by.x=field, by.y="Field",    
    sort=FALSE)$Color)
  rowColors = as.character(factor(members, labels=mColors$color))
 
  png(filename=paste("heatmap","-", field, "-", depth, ".png", sep=""), 
    width=3000, height=21000, res=300)
  heatmap.2(samples2, Rowv=as.dendrogram(tree), Colv=FALSE, dendrogram = 'row', 
    col = c("red", "green"), 
    breaks = c(-1, 0, 1), 
    trace = 'none', 
    key = FALSE, 
    RowSideColors = rowColors,
    ColSideColors = colColors
  )
  dev.off()
}

#' Performs enrichment analysis of traits against a network dendrogram
#'
#' TODO: add documentation here.
#'
#' 
#' @param tree
#'   An instance of an hclust object.
#' @param osa
#'   The sample annotation matrix. One column must contain the header 'Sample'
#'   and the remaining colums correspond to an annotation type.  The rows
#'   of the anntation columns should contain the annotations.
#' @param net
#'   A network-style representation of the similarity matrix. This is a 
#'   data.frame that must have a column named 'Source', another named 
#'   'Target' and a third named 'Similarity'.  Any other columns are allowed
#'   but ignored.
#' @export
#' @examples
#'   
analyzeNetDendro = function(tree, osa, net, alpha = 0.001, min_presence = 0.95,
  min_cluster_size = 3) {
  height_array = unique(tree$height)
  height_array = sort(height_array, decreasing = TRUE)
  seen = list()
 
  # Initialize seen_checksums empty vector
  seen_checksums <- vector(mode="character", length=0)
  prev_indexes = list()
 
  # Iterate through the dendrogram at increasing heights and perform
  # enrichment analysis
  output = file("output.txt", "w")
  write(paste("Cluster_ID", "Parent_ID", "Depth", "Height", "Clusters", "Cluster_Num", "Size", 
    "Avg Degree", "Type", "Counts", "Categories", "Enrichment", "Nodes", "Edge Index", 
    sep="\t"), file=output, append=TRUE)
  depth = 1
  for (height in height_array) {
    members = cutree(tree, h = height)
    num_clusters = length(unique(members))
    #drawDendro(osa, members, samples, depth, 'Treatment')
  
    prev_indexes[[depth]] = list() 
    for(i in unique(members)) {
       print(paste("Depth: ", depth, ", Num Clusters: ", num_clusters, ", Cluster: ", i, sep=""))
 
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
 
         # Run enrichment
         if (length(cluster_indexes) >= min_cluster_size) {
           for (field in fields) {
             results = test_module(net, cluster_indexes, osa, field, min_presence)
             enrichment = results[[1]]
             counts = results[[2]]
             #if (length(which(enrichment < alpha)) > 0) {
              write(
                 paste(
                   clust_chksm,
                   parent,
                   depth,
                   height,
                   num_clusters,
                   i,
                   length(cluster_indexes),
               avg_degree,
                   field,
                   paste(counts, collapse=","),
                   paste(unique(osa[,c(field)]), collapse=","),
                   paste(enrichment, collapse=","), 
                   paste(cluster_nodes, collapse=","), 
                   paste(cluster_indexes, collapse=","),
                   sep="\t"),
                  file = output, 
                 append = TRUE
               )
             #} # end if(length(which(enrichment….
              } # end for (field in fields)
        } # end if (length(cluster_indexes
        else {
          print("Skipping, cluster too small")
        }
      } # end if (!is.element(clust_chksm, 
      else {
        print("Skipping, already seen")
      }
    } # end for(i in unique(members)) {
    depth = depth + 1
  } # end for (height in height_array) {
  close(output)
}

