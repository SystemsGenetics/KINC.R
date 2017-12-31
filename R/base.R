#' Imports a network file produced by KINC into a data frame.
#'
#' @param network_file
#'   The full path to the network file on the file system.
#'
#' @return
#'    A dataframe containing the network file contents
#'
#' @export
loadNetwork = function(network_file, KINC_version = '1.0') {
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
graphNet = function(net, osa = data.frame()) {
  edge_indexes = c(1:length(net$Source));
  g = graphEdgeList(edge_indexes, net, osa)
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
graphEdgeList = function(edge_indexes, net, osa = data.frame(), field = NA) {

  edges = as.vector(unlist(t(net[edge_indexes, c('Source', 'Target')])))
  g = graph(edges, directed = F)

  # Don't show node labels if we have too many nodes.
  vlc = 1
  vs = 10
  if (length(edge_indexes) > 100) {
    vlc = 0.25
    vs = 5
  }
  if (length(edge_indexes) > 1000) {
    vlc = 0.001
    vs = 2
  }

  if (!is.na(field)) {
    sample_types = as.character(osa[[field]])
    sample_types[which(sample_types == "null")] = NA
    sample_types[which(sample_types == "notreported")] = NA

    categories = sort(unique(sample_types[nona.idx]))
    max = lapply(categories, FUN=function(x) {
      subname = paste(field, x, sep='_')
      return(min(net[, c(subname)], na.rm = TRUE))
    })
  }

  # Plot the graph.
  l = layout_with_kk(g)
  plot(g,
       # Verticies
       vertex.label.color="black", vertex.label.cex = vlc, vertex.size = vs,
       # Edges
       edge.color=rep("#AAAAAA", length(edge_indexes)), vertex.color='cyan',
       # Layout
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
plot3DPair = function(i, j, osa, ematrix, field, colors = NA, samples = NA) {
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

  plot3d(x, y, z, type = 's', size = 0.5, col = colors)
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
#' @export
plot3DEdgeList = function(edge_indexes, osa, net, ematrix, field, colors = NA, samples = NA) {
  rgl.open()
  rgl.bg(color = "white")

  for (i in edge_indexes) {
    source = net[i, 'Source']
    target = net[i, 'Target']

    if (!is.na(samples)) {
      x = t(ematrix[source, samples])
      y = t(ematrix[target, samples])
      z = as.factor(osa[samples, field])
      colors = colors[samples]
    }
    else {
      x = t(ematrix[source, ])
      y = t(ematrix[target, ])
      z = as.factor(osa[, field])
    }

    plot3d(x, y, z, type = 's', size = 0.5, col = colors, xlab=source, ylab=target, zlab=field)
    text3d(x, y, z, osa[[field]], adj=c(1, 1))
    val = readline(prompt="Press enter to continue to next plot. Press 'q' and enter to quit.")

    if(val == 'q') {
      return()
    }
  }
}

#' Plots a 2D scatterplot for a given list of edges.
#'
#' @param edge_indexes
#'   The index of the edge in the network data frame.
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param ematrix
#'   The expression matrix.
#' @param cor_col
#'   Then name of the column in the net argument that provides the similarity score.
#' @param colors
#'   An array of colors per sample.
#' @param samples
#'   An array of sample names for which the scatterplot should be included in the plot.
#'   If no value is provided then all samples with values are included.
#' @export
plotEdgeList = function(edge_indexes, net, ematrix, cor_col = 'sc', colors = NA, samples = NA) {
  for (i in edge_indexes) {
    source = net[i, 'Source']
    target = net[i, 'Target']

    if (!is.na(samples)) {
      x = t(ematrix[source, samples])
      y = t(ematrix[target, samples])
    }
    else {
      x = t(ematrix[source, ])
      y = t(ematrix[target, ])
    }

    # Set a default point color of black.
    colors = rep('#000000', length(x))

    # If there is a sample string then set the colors to match that.
    if ('Samples' %in% names(net)) {
      sample_str = net[i, 'Samples']
      samples_vec = as.numeric(strsplit(sample_str, "")[[1]])
      samples_vec[which(samples_vec != 1)] = 0
      colors[which(samples_vec == 1)] = '#FF0000'
    }

    main = paste('Cor: ', net[i, cor_col]);
    if ('Missing_Samples' %in% names(net)) {
      main = paste(main, '; ', 'Missing: ', net[i, 'Missing_Samples'], sep="")
    }
    if ('Cluster_Samples' %in% names(net)) {
      main = paste(main, '; ', 'Size: ', net[i, 'Cluster_Samples'], sep="")
    }
    if ('Num_Clusters' %in% names(net)) {
      main = paste(main, '; ', 'Clusters: ', net[i, 'Num_Clusters'], sep="")
    }

    plot(x, y, col = colors, main=main)
    val = readline(prompt="Press enter to continue to next plot. Press 'q' and enter to quit.")

    if(val == 'q') {
      return()
    }
  }
}
