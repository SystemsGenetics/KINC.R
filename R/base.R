#' Imports a network file produced by KINC into a data frame.
#'
#' @param network_file
#'   The full path to the network file on the file system.
#'
#' @return
#'    A dataframe containing the network file contents
#'
#' @export
loadKINCNetwork = function(network_file, KINC_version = '1.0') {
  # Load in the KINC network file
  colNames = c(
    "Source", "Target", "sc", "Interaction", "Cluster",
    "Num_Clusters", "Cluster_Samples", "Missing_Samples", "Cluster_Outliers",
    "Pair_Outliers", "Too_Low", "Samples"
  )
  colClasses = c(
    "character", "character", "numeric", "character", "numeric",
    "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
    "character"
  )
  net = read.table(network_file, header = TRUE,
    sep="\t", colClasses=colClasses, col.names=colNames)
  #net = read.big.matrix(network_file, sep="\t", header = TRUE,
  #                      col.names=colNames, type=colClasses)

}

#' Exports a network data frame as a KINC compatible file.
#'
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param network_file
#'   The path to which the network file will be saved.
#'
#' @export
saveKINCNetwork = function(net, network_file) {
  write.table(net, file= network_file, sep = "\t", na="NA", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

#' Imports sample annotations.
#'
#' The format of the sample annotation file is tab delimited where the
#' first column is the sample name and every other column contains
#' the annotation.  The file must contain a header with the annotation
#' types and the first column of the header should be the word 'Sample'
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
loadSampleAnnotations = function (annotation_file, sample_order_file) {
  sample_annots = read.table(annotation_file, sep="\t", header=TRUE, row.names=NULL)
  sample_order = read.table(sample_order_file, colClasses=c('character'),
                            col.names=c('Sample'))
  osa = merge(sample_order, sample_annots, by = "Sample", sort=FALSE)
  row.names(osa) = osa$Sample

  return(osa)
}
#' Draws the graph of the entire network
#'
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param osa
#'   The sample annotation matrix. One column must contain the header 'Sample'
#'   and the remaining colums correspond to an annotation type.  The rows
#'   of the anntation columns should contain the annotations.
#' @export
graphNet = function(net, osa = data.frame(), sim_col = 'sc') {
  edge_indexes = c(1:length(net$Source));
  g = graphEdgeList(edge_indexes, net, sim_col, osa)
  return(g)
}

#' Draws the graph of a set of edges in the network
#'
#' @param edge_indexes
#'   The index of the edges in the network that comprise the module.
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param osa
#'   The sample annotation matrix. One column must contain the header 'Sample'
#'   and the remaining colums correspond to an annotation type.  The rows
#'   of the anntation columns should contain the annotations.
#'
#' @export
graphEdgeList = function(edge_indexes, net, sim_col = 'sc', osa = data.frame()) {

  g = graph.edgelist(as.matrix(net[edge_indexes, c('Source', 'Target')]), directed = F)
  E(g)$weight = abs(net[edge_indexes, sim_col])

  # Don't show node labels if we have too many nodes.
  vlc = 0.5
  vs = 10
  if (length(edge_indexes) > 100) {
    vlc = 0.25
    vs = 5
  }
  if (length(edge_indexes) > 1000) {
    vlc = 0.001
    vs = 2
  }

  # Plot the graph.
  l = layout_with_kk(g)
  plot(g,
       vertex.label.color="black", vertex.label.cex = vlc, vertex.size = vs,
       edge.color='black', vertex.color='cyan',
       layout = l)

  return(list(
    graph = g,
    layout = l
  ))
}

#' Colors the edges of a given edge list on the currently plotted graph.
#'
#'
#' @param edge_indexes
#'   The index of the edges in the network that comprise the module.
#' @param graph
#'   The igraph object as returned from graphNet or graphEdgeList.
#' @param layout
#'   The layout object as returned from graphNet or graphEdgeList
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param color
#'   The color that should be used to highlith the edges
#'
#' @export
highlightEdgeList = function(edge_indexes, graph, layout, net, color) {
   edges = net[edge_indexes, c('Source', 'Target')]
   verts = names(V(graph))
   layout = norm_coords(layout)
   for (i in 1:nrow(edges)) {
     source = edges[i, c('Source')]
     target = edges[i, c('Target')]
     sindex = which(verts == source)
     tindex = which(verts == target)
     x0 = layout[sindex, 1]
     y0 = layout[sindex, 2]
     x1 = layout[tindex, 1]
     y1 = layout[tindex, 2]
     segments(x0, y0, x1, y1, col = color)
   }
}

#' Plots a 3D scatterplot for a given pair of nodes.
#'
#' @param i
#'   The index of the first node in the expression matrix.
#' @param j
#'   The index of the second node in the expression matrix.
#' @param osa
#'   The sample annotation matrix. One column must contain the header 'Sample'
#'   and the remaining colums correspond to an annotation type.  The rows
#'   of the anntation columns should contain the annotations.
#' @param ematrix
#'   The expression matrix.
#' @param field
#'   The field in the sample annotation matrix by which to split the points
#'   into separate layers.
#' @param colors
#'   An array of colors per sample.
#' @export
plot3DPair = function(i, j, osa, ematrix, field, colors = NA, samples = NA, highlight = c()) {
  rgl.open()
  rgl.bg(color = "white")

  if (!is.na(samples)) {
    x = as.factor(osa[samples, field])
    y = t(ematrix[i, samples])
    z = t(ematrix[j, samples])
  }
  else {
    x = as.factor(osa[, field])
    y = t(ematrix[i, ])
    z = t(ematrix[j, ])
  }

  size = rep(0.25, length(x))
  if (length(highlight) > 0) {
    size[highlight] = 0.5
  }

  plot3d(x, y, z, type = 's', size = size, col = colors)
}

#' Plots a 3D scatterplot for a given list of edges.
#'
#' @param edge_indexes
#'   The index of the edge in the network data frame.
#' @param osa
#'   The sample annotation matrix. One column must contain the header 'Sample'
#'   and the remaining colums correspond to an annotation type.  The rows
#'   of the anntation columns should contain the annotations.
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param ematrix
#'   The expression matrix.
#' @param field
#'   The field in the sample annotation matrix by which to split the points
#'   into separate layers.
#' @param colors
#'   An array of colors per sample.
#' @param samples
#'   Limit the plot to only the samples indexes provided.
#' @param highlight
#'   If set to TRUE then the samples in the edge cluster are drawn larger than the
#'   other samples in the plot.
#' @export
plot3DEdgeList = function(edge_indexes, osa, net, ematrix, field, label_field = field,
                          colors = NA, samples = NA, highlight = TRUE) {
  rgl.open()
  rgl.bg(color = "white")

  # Set the row names to be the sample names.
  row.names(osa) = osa$Sample

  # Loop until the user decides to leave.
  j = 1
  done = FALSE;
  while (!done) {
    i = edge_indexes[j]

    source = net[i, 'Source']
    target = net[i, 'Target']
    sample_indexes = c()
    if ('Samples' %in% names(net)) {
      sample_indexes = getEdgeSamples(i, net)
    }

    if (!is.na(samples)) {
      x = t(ematrix[source, samples])
      y = t(ematrix[target, samples])
      z = as.factor(osa[[field]][samples])
      colors = colors[samples]
    }
    else {
      x = t(ematrix[source, ])
      y = t(ematrix[target, ])
      z = as.factor(osa[, field])
    }

    size = rep(0.25, length(x))
    if (highlight & length(sample_indexes) > 0) {
      size[sample_indexes] = 0.5
    }

    main = paste(source, 'vs', target, 'Edge:', i)
    plot3d(x, y, z, type = 's', main = main, size = size, col = colors, xlab=source, ylab=target, zlab=field, axes=FALSE)
    z_labels = unique(sort(as.factor(osa[[field]])))
    box3d()
    axis3d('x')
    axis3d('y')
    axis3d('z', nticks = length(unique(osa[, field])) - 1, labels = z_labels)
    if (!is.na(samples)) {
      text3d(x, y, z, osa[names(ematrix[samples]), ][[label_field]], adj=c(1, 1))
    }
    else {
      text3d(x, y, z, osa[names(ematrix), ][[label_field]], adj=c(1, 1))
    }

    # add legend
    #legend3d("topright", legend = paste('Type', c('A', 'B', 'C')), pch = 16, col = rainbow(3), cex=1, inset=c(0.02))

    # Use the key press to navigate through the images.
    val = readline(prompt="Type control keys the [enter]. Keys: n > forward, b > backwards, q > quit.")
    if (val == 'n') {
      if (j < length(edge_indexes)) {
        j = j + 1;
      }
    }
    if (val == 'b') {
      if (j > 1) {
        j = j - 1
      }
    }
    if(val == 'q') {
      done = TRUE
    }
  }
}

#' Plots a scatterplot of gene expression across a set of samples.
#'
#' Samples will be ordered by the values of field.
#'
#' @param gene
#'   The name of the gene to plot
#' @param ematrix
#'   The expression matrix.
#' @param osa
#'   The sample annotation matrix. One column must contain the header 'Sample'
#'   and the remaining colums correspond to an annotation type.  The rows
#'   of the anntation columns should contain the annotations.
#' @param field
#'   The name of the field by which to order the points along the
#'   x-axis.
#' @param xlab
#'   The label for the x-axis of the plot.
#' @export
plotGene = function(gene, ematrix, osa, field, xlab = field) {
  condition = osa[c(field)]
  expdata = merge(t(ematrix[gene,]), condition, by="row.names")
  colnames(expdata) = c('Sample', 'y', 'x')
  expdata = expdata[complete.cases(expdata),]

  expplot = ggplot(expdata, aes(x, y)) +
    geom_point(size=1) +
    xlab(xlab) + ylab('Expression Level')
  print(expplot)
  return(expplot)
}

#' Plots a 2D scatterplot for a given list of edges.
#'
#' @param edge_indexes
#'   The index of the edge in the network data frame.
#' @param osa
#'   The sample annotation matrix. One column must contain the header 'Sample'
#'   and the remaining colums correspond to an annotation type.  The rows
#'   of the anntation columns should contain the annotations.
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param ematrix
#'   The expression matrix.
#' @param field
#'   The field in the sample annotation matrix by which to split the points
#'   into separate layers.
#' @param samples
#'   An array of sample names for which the scatterplot should be included in the plot.
#'   If no value is provided then all samples with values are included.
#' @param highlight
#'   If set to TRUE, makes the samples belonging to the cluster that
#'   underlies the edge larger in the plot. Default = TRUE.
#' @param title
#'   A title to give the plot. Default = NULL
#' @export
plot2DEdgeList = function(edge_indexes, osa, net, ematrix,
                          field = NA, samples = NA, highlight = TRUE,
                          title = NULL) {
  j = 1
  done = FALSE;
  while (!done) {
    i = edge_indexes[j]
    source = net[i, 'Source']
    target = net[i, 'Target']
    sample_indexes = c()
    if ('Samples' %in% names(net)) {
      sample_indexes = getEdgeSamples(i, net)
    }

    if (!is.na(samples)) {
      x = t(ematrix[source, samples])
      y = t(ematrix[target, samples])
      colors = colors[samples]
      condition = osa[[field]][samples]
    }
    else {
      x = t(ematrix[source, ])
      y = t(ematrix[target, ])
      condition = osa[[field]]
    }

    size = rep(0.25, length(x))
    if (highlight & length(sample_indexes) > 0) {
      size[sample_indexes] = 1.5
    }

    coexpdata = data.frame(source = x, target = y, category = condition, size = size)
    colnames(coexpdata) = c('x', 'y', 'category', 'size')
    coexpdata = coexpdata[complete.cases(coexpdata),]
    coexpplot = ggplot(coexpdata, aes(x, y, color=category)) +
      geom_point(size=coexpdata$size) +
      xlab(source) + ylab(target) + labs(colour=field)
    if (!is.numeric(condition[0])) {
      coexpplot = coexpplot + scale_color_brewer(palette="Set1")
    }
    if (!is.null(title)) {
      coexpplot = coexpplot + ggtitle(title)
    }
    print(coexpplot)

    if (length(edge_indexes) == 1) {
      return(coexpplot)
    }
    # Use the key press to navigate through the images.
    val = readline(prompt=paste("Edge: ", i, '. ', j, " of ", length(edge_indexes), ". Type control keys the [enter]. Keys: n > forward, b > backwards, q > quit.", sep=""))
    if (val == 'n') {
      if (j < length(edge_indexes)) {
        j = j + 1;
      }
    }
    if (val == 'b') {
      if (j > 1) {
        j = j - 1
      }
    }
    if(val == 'q') {
      done = TRUE
    }
  }
}

#' Finds Linked communities in the network.
#'
#' This function generates three output files that it writes to the current
#' working directory.
#'
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param file_prefix
#'   A prefix to add to the beginning of each file name.
#' @param module_prefix
#'   A prefix to add to the beginning of the module names. By deafult this is simply
#'   the letter 'M'.
#' @param sim_col
#'   The column name in the network data frame for the similarity score. Default is 'sc',
#' @param hcmethod
#'   A character string naming the hierarchical clustering method to use. Can be one of
#'   "ward.D", "single", "complete", "average", "mcquitty", "median", or "centroid".
#'   Defaults to "single".
#' @param meta
#'   Indicates if modules should be collapsed into meta-modules. If set to TRUE then
#'   the linked communities returned are meta modules. Defaults to FALSE.
#'
#' @return
#'   The linked communities object.
#' @export
#'
findLinkedCommunities = function(net, file_prefix, module_prefix = 'M', sim_col = 'sc',
                                 hcmethod = 'complete', meta = TRUE, ignore_inverse = TRUE) {
  new_net = net
  new_net$Module = NA

  # If the user requested to ignore inverse correlation edges we'll take those out.
  if (ignore_inverse) {
    g = graph.edgelist(as.matrix(net[which(net[[sim_col]] > 0), c('Source', 'Target')]), directed = F)
    E(g)$weight = abs(net[which(net[[sim_col]] > 0), sim_col])
  }
  else {
    g = graph.edgelist(as.matrix(net[, c('Source', 'Target')]), directed = F)
    E(g)$weight = abs(net[, sim_col])
  }

  # Remove duplicated edges or it drammatically slows down cluster merging.
  #g = unique(g)

  # Linked community method doesn't do well with disconnected graphs. So, we want
  # to decompose the graph into its disjointed subgraphs.  If the callee doesn't want
  # modules to span inverse edges then remove those.
  subg = decompose(g,  min.vertices = 30);

  lcs = list()

  # Iterate through the subgraphs and find the communities.
  for (gi in 1:length(subg)) {

    print(paste('Working on subgraph', gi, 'of', length(subg), '...', sep=" "))

    # We need a large enough subnetwork to find clusters.
    subnet = as_edgelist(subg[[gi]])

    # Create the Link Community clusters.
    lc = getLinkCommunities(subnet, hcmethod = hcmethod, removetrivial = FALSE, plot = FALSE)

    # If the meta analysis was selected. Rather than use the cluster merging
    # function that comes with LCM we perform our own. This is because the
    # heirarchical clustering method of LCM improperly merges some clusters that
    # are not connected.
    if (meta) {
      r = mergeCommunities(lc)
    }
    else {
      r = list(
        cedges = lc$clusters,
        cnodes = lc$nodeclusters
      )
    }

    lcs[[gi]] = r

    # Iterate through all of the edges in the graph and add the module to the
    # network
    for (i in 1:length(r$cedges)) {
      for (j in 1:length(r$cedges[[i]])) {
        edge = subnet[r$cedges[[i]][j],]
        node1 = edge[1]
        node2 = edge[2]
        cluster = i

        ei = which(new_net$Source == node1 & new_net$Target == node2)
        if (length(ei) > 0) {
          module = sprintf("SG%02dM%04d", gi, cluster)
          new_net[ei,]$Module = module
        }
        ei = which(new_net$Target == node1 & new_net$Source == node2)
        if (length(ei) > 0) {
          module = sprintf("SG%02dM%04d", gi, cluster)
          new_net[ei,]$Module = module
        }
      }
    }
  }


  # Convert all of the arrays above into a data frame for printing the modules
  # list.
  write.table(new_net, file=paste(file_prefix, "coexpnet.edges.lcm.txt", sep=".") ,sep="\t", row.names=FALSE, append=FALSE, quote=FALSE)

  # Write out the node list
  node_list = data.frame(Node=c(as.character(new_net$Source), as.character(new_net$Target)), Cluster=c(new_net$Module, new_net$Module))
  write.table(unique(node_list), file=paste(file_prefix, "coexpnet.edges.lcmByNodes.txt", sep=".") ,sep="\t", row.names=FALSE, append=FALSE, quote=FALSE, col.names=FALSE)

  # convert the edge and cluster into a report of modules by edge which can
  # be used for Cytoscape.
  new_net_edges = cbind(paste(new_net$Source, "(co)", new_net$Target, sep=" "), new_net$Module)
  colnames(new_net_edges) = c('Edge', 'Module')
  new_net_edges = new_net_edges[which(!is.na(new_net_edges[,2])),]
  write.table(new_net_edges, file=paste(file_prefix, "coexpnet.edges.lcm.cys.txt", sep=".") ,sep="\t", row.names=FALSE, append=FALSE, quote=FALSE)

  return(lcs)
}

#' Merges linked community clusters that have the most similar sets of nodes.
#'
#' This merging function recursively iterates through all of the clusters
#' and performs a pair-wse Jaccard comparision between all clusters.
#' Those with the highest Jaccard score that are above the given
#' threshold (i.e. deafult of 0.5) are candidates for merging.
#'
#' This is a helper function for the findLinkedCommunities() function and is
#' not meant to be called on its own.
#'
#' @param lc
#'   A linkcomm object.
#'
#' @return
#'   A list were each element of the list is the
#'   set of nodes and edges of the merged clusters.

mergeCommunities = function(lc){

  cedges = lc$clusters
  cnodes = lc$nodeclusters
  cnodes$cluster = as.integer(cnodes$cluster)
  r = mergeClusters(cedges, cnodes)

  return(r)
}

#' The recursive merging function called by mergeCommunities().
#'
mergeClusters = function(cedges, cnodes, th = 0.5) {
  nclusters = length(cedges)

  # Create a dataframe for storing the best pairwise Jaccard similarity scores.
  best = data.frame(i = 0, j = 0, ji = 0)

  for (i in 1:nclusters) {
    if (best$ji[1] == 1) next
    for (j in 1:nclusters) {
      if (best$ji[1] == 1) next
      if (j >= i) next
      A = as.character(cnodes$node[which(cnodes$cluster == i)])
      B = as.character(cnodes$node[which(cnodes$cluster == j)])
      overlap_i = length(intersect(A, B))/length(A)
      overlap_j = length(intersect(A, B))/length(B)

      if (overlap_i > best$ji[1]) {
        best$i[1] = i
        best$j[1] = j
        best$ji[1] = overlap_i
      }
      if (overlap_j > best$ji[1]) {
        best$i[1] = i
        best$j[1] = j
        best$ji[1] = overlap_j
      }
    }
  }

  # Now merge the best two clusters as long as the jaccard score is above
  # our given threshold.
  if (best$ji[1] > th) {
    i = best$i[1]
    j = best$j[1]
    cat(paste("Merging", i, 'and', j, 'from', nclusters, 'similarity:', sprintf("%.02f", best$ji[1]), "\r", sep=" "))
    cedges[[i]] = union(cedges[[i]], cedges[[j]])
    cedges[[j]] = NULL
    cnodes$cluster[which(cnodes$cluster == j)] = i
    cnodes$cluster[which(cnodes$cluster > j)] = cnodes$cluster[which(cnodes$cluster > j)] - 1
    cnodes = cnodes[which(!duplicated(cnodes)),]
    r = mergeClusters(cedges, cnodes)
    return(r)
  }

  return(list(
    cedges = cedges,
    cnodes = cnodes
  ))
}


