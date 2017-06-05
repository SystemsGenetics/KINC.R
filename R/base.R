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
