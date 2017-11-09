#' Imports a network file produced by KINC into a data frame.
#'
#' @param network_file
#'   The full path to the network file on the file system.
#'
#' @return
#'    A dataframe containing the network file contents
#'
#' @export
loadNetwork = function(network_file) {
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
graphNet = function(net, osa) {
  edge_indexes = c(1:length(net$Source));
  g = graphEdgeList(edge_indexes, osa, net)
  return(g)
}

#' Draws the graph of a set of edges in the network
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
#'
#' @export
graphEdgeList = function(edge_indexes, osa, net) {
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

  # Plot the graph.
  l = layout_with_kk(g)
  plot(g,
       # Verticies
       vertex.label.color="black", vertex.label.cex = vlc, vertex.size = vs,
       # Edges
       edge.color="#AAAAAA", vertex.color='cyan',
       # Layout
       layout = l)
  return(list(
    graph = g,
    layout = l
  ))
}

#' Colors the edges of a given edge list on the currentlly plotted graph.
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
