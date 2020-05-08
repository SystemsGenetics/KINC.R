#' Prior to importing a KINC network, this function checks its Size.
#'
#' For large networks it may be useful to preserve system memory
#' by reading in the network in subsets.  Prior to calling the
#' loadKINCNetwork function this function can be used to
#' determine the size of the network in terms of number of
#' edges.  Once the size is know the nrows and skip arguments
#' of the loadKINCNetwork function can be used to load
#' subsets of the network at a time.
#'
#' @param network_file
#'   The path to the network file on the file system.
#'
#' @return
#'   An integer indicating the number of edges in the
#'   network.
#'
#' @export
getKINCNetworkSize = function (network_file) {
  num_lines = countLines(network_file)
  # Subtract 1 from the total to account for the header.
  return (num_lines[1] - 1)
}

#' Imports a network file produced by KINC into a data frame.
#'
#' @param network_file
#'   The path to the network file on the file system.
#' @param nrows
#'   For large network files it may be useful to preserve system
#'   memory be reading in subsets of large network files.  This
#'   argument specifies the number of rows of the network that
#'   should be read.  Use this in combination with the skip
#'   argument to cycle through the entire network file
#' @param skip
#'   For large network files it may be useful to preserve system
#'   memory be reading in subsets of large network files.  This
#'   argument specifies the number of rows to skip before reading.
#'   Use this in combination with the nrows argument to cycle
#'   through the entire network file
#'
#' @return
#'    A dataframe containing the network file contents
#'
#' @export
loadKINCNetwork = function(network_file, nrows=-1, skip=0) {

  # First get the headers so we can tell if this is full or tidy network file.
  type = 'text'
  headers = read.table(file=network_file, header = TRUE, sep="\t", nrows=1, skip=0, stringsAsFactors = FALSE)
  if ('Test_Name' %in% colnames(headers)) {
    type = 'tidy'
  }

  # If we are skipping rows then we need to tell the loader not to try to
  # retrieve a header.
  header = TRUE
  if (skip > 0) {
    header = FALSE
  }

  if (type == 'text') {
    # Get the number of columns in the file and set the column classes.
    # The first 7 columns are set by KINC.
    ncols = max(count.fields(network_file, sep = "\t"))
    colClasses = c(c(
      "character", "character", "numeric", "character", "numeric",
      "numeric", "character"
    ), rep("numeric", ncols-7))

    net = read.table(file=network_file, header = header, sep="\t", colClasses=colClasses, nrows=nrows, skip=skip)
  }
  if (type == 'tidy') {
    # Get the number of columns in the file and set the column classes.
    # The tidy format has 10 columns.
    colClasses = c(
      "character", "character", "numeric", "character", "numeric",
      "numeric", "character", "character", "numeric", "numeric")
    net = read.table(file=network_file, header = header, sep="\t", colClasses=colClasses, na.strings="nan", nrows=nrows, skip=skip)
  }

  colnames(net) = colnames(headers)
  return (net)
}

#' Imports a network file produced by KINC into a data frame.
#'
#' @param GEM_file
#'   The path to the GEM file on the file system. This file should
#'   contain gene expression (or abundance data) where the rows
#'   are genes (or compounds) and the columns are samples.  The first
#'   column should countain the gene name.  The first row contain a
#'   a header with only the sample names. Thus the header row has
#'   one less values then every other row.
#'
#' @return
#'    A dataframe containing the GEM.
#'
#' @export
loadGEM = function(GEM_file) {
  emx = read.table(GEM_file, header=TRUE, sep="\t")
  return(emx)
}
#' Exports a network data frame as a KINC compatible file.
#'
#' @param net
#'   A network data frame int tidy format containing the KINC-produced network.
#'   The loadKINCNetwork function can be used to import the network.
#' @param network_file
#'   The path to which the network file will be saved.
#' @param append
#'   For large network files it may be useful to preserve system
#'   memory be reading in subsets of large network files, processing
#'   and then saving the processed subset before moving on to the
#'   next. Set this to TRUE to append to an existing network file.
#'   For the header to appear correctly, be sure to set append
#'   to FALSE when saving the first subset.
#'
#' @export
saveKINCNetwork = function(net, network_file, append=FALSE) {
  col.names = TRUE
  if (append == TRUE) {
    col.names = FALSE
  }

  num_cols = dim(net)[2]
  if (num_cols > 12) {
    for (i in 13:num_cols) {
      if (is.numeric(net[,i])) {
        net[,i] = format(net[,i], digits=4, nsmall=1)
      }
    }
  }
  write.table(format(net, digits=4), file = network_file,
              sep = "\t", na="", quote=FALSE,
              row.names=FALSE, col.names=col.names, append=append)
}
#' Exports a condition-specific network as a Tidy KINC compatible file.
#'
#' This function will save as many condition-specific network files as directed.
#'
#' @param net
#'   A network data frame int tidy format containing the KINC-produced network.
#'   The loadKINCNetwork function can be used to import the network.
#' @param network_file
#'   The path to which the network file will be saved. The token %condition%
#'   should be present in the filename. The condition name will be substituted
#'   in place of this token resulting in unique files for each condition. The
#'   token %filter% can also be provided. If present the value provided to
#'   the filter argument will be substituted.
#' @param conditions
#'   An array of condition names. These must be present in the Test_Name
#'   column of the network data frame. If no conditions are provided then
#'   all of the conditions in the Test_Name column will be used.
#' @param filter
#'   The type of condition-specific network to create: 'label', 'class' or 'full'.
#'   If 'label' is provided only edges that have no other signficant tests other
#'   than the one specifid are included.  If 'class' is provided only edges in one
#'   type of class will be included. For example, if the Test_Name value is
#'   'Subspecies__Japonica' the class is 'Subspecies'.  If no value is provided
#'   then the filter defaults to 'full'.
#' @export
  saveConditionKINCNetwork = function(net, network_file, conditions=c(), filter='label') {
  if (!'Test_Name' %in% colnames(net)) {
    message('ERROR: Please provide a network dataframe in Tidy format to save condition-specific networks.')
    return(NA)
  }

  if (length(conditions) == 0) {
    conditions = unique(net$Test_Name)
  }

  if (filter == 'class') {
    message(paste0('Filtering for edges with unique conditional classes...'))
    net$Class = str_replace(net$Test_Name, '__.*$','')
    edge_grp = net %>% group_by(Source,Target)
    net = edge_grp %>% filter(length(unique(Class)) == 1) %>% as.data.frame()
    for (class in unique(net$Class)) {
      message(paste0('Writing the condition-specific network for the class: ', class, '...'))
      csGCN = net[which(net$Class == class),]
      if (dim(csGCN)[1] == 0) {
        message(paste0('The class ', class, ' is not present in the network dataframe. Skipping.'))
        next
      }
      csGCN_file = str_replace(network_file, '%condition%', class)
      csGCN_file = str_replace(csGCN_file, '%filter%', 'unique_class')
      saveKINCNetwork(csGCN, csGCN_file)
    }

  } else if (filter == 'label') {
    message(paste0('Filtering for edges with unique conditional labels...'))
    edge_grp = net %>% group_by(Source,Target)
    net = edge_grp %>% filter(n() == 1) %>% as.data.frame()
    for (condition in conditions) {
      message(paste0('Writing the condition-specific network for the condition: ', condition, '...'))
      csGCN = net[which(net$Test_Name == condition),]
      if (dim(csGCN)[1] == 0) {
        message(paste0('The condition ', condition, ' is not present in the network dataframe. Skipping.'))
        next
      }
      csGCN_file = str_replace(network_file, '%condition%', condition)
      csGCN_file = str_replace(csGCN_file, '%filter%', 'unique_label')
      saveKINCNetwork(csGCN, csGCN_file)
    }

  } else {
    for (condition in conditions) {
      message(paste0('Writing the condition-specific network for the condition: ', condition, '...'))
      csGCN = net[which(net$Test_Name == condition),]
      if (dim(csGCN)[1] == 0) {
        message(paste0('The condition ', condition, ' is not present in the network dataframe. Skipping.'))
        next
      }
      csGCN_file = str_replace(network_file, '%condition%', condition)
      csGCN_file = str_replace(csGCN_file, '%filter%', 'all')
      saveKINCNetwork(csGCN, csGCN_file)
    }
  }
}

#' Converts the KINC network dataframe to an iGraph object.
#'
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @return
#'   An iGraph object.
#'
#' @export
KINCtoIgraph = function(net) {
  g = graph.edgelist(as.matrix(net[, c('Source', 'Target')]), directed = FALSE)
  return(g)
}

#' Creates all PNG figures describing p-value and r-squared distributions of the network.
#'
#' The network must be in tidy format and must include condition-specfic p-values.
#'
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param out_prefix
#'   The prefix of the output images files.
#'
#' @return dataframe
#'   A network dataframe with non-significant edges removed.
#'
#' @export
saveKINCplots = function(net, out_prefix = "KINC_network") {

  if (!'Test_Name' %in% colnames(net)) {
    message('ERROR: Please provide a network dataframe in Tidy format to generate plots.')
    return(NA)
  }

  has_r2 = TRUE
  if(sum(is.na(net$r_squared)) == dim(net)[1]) {
    has_r2 = FALSE
  }

  plotPvalPerCondition(net, out_prefix)
  plotPvalPerScore(net, out_prefix)
  plotPvalPerScoreAndCondition(net, out_prefix)
  plotPvalHistPerCondition(net, out_prefix)

  if (has_r2) {
    plotRsqrPerScoreAndCondition(net, out_prefix)
    plotRsqrPerScore(net, out_prefix)
    plotRsqrHistPerCondition(net, out_prefix)
    plotRsqrPerCondition(net, out_prefix)
  }

  plotEdgeCountPerCondition(net, out_prefix)

}

#'Plots the distribution of p-values per condition.
#'
#'@param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#'@param out_prefix
#'   The prefix of the output images files.
#'
#'@export
plotPvalPerCondition = function(net, out_prefix=NA) {
  if (!'Test_Name' %in% colnames(net)) {
    message('ERROR: Please provide a network dataframe in Tidy format to generate plots.')
    return(NA)
  }

  # Explore the distribution of p-values & rSqr for each trait
  if (!is.na(out_prefix)) {
    print(paste0("Saving: ", out_prefix, '.PvalPerCondition.png'))
    png(paste0(out_prefix, '.PvalPerCondition.png'), width=6 ,height=6, units="in", res=300)
  }
  plot = ggplot(net, aes(x=Test_Name, y=p_value, fill=Test_Name)) +
    geom_boxplot() +
    scale_y_log10(breaks=c(1e-3, 1e-5, 1e-10, 1e-15, 1e-20)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "none") +
    xlab("") + ylab(bquote(-log^10~(p)))
  print(plot)
  if (!is.na(out_prefix)) {
    dev.off()
  }
}
#'Plots the distribution of p-values per condition.
#'
#'@param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#'@param out_prefix
#'   The prefix of the output images files.
#'
#'@export
plotRsqrPerCondition = function(net, out_prefix=NA) {
  if (!'Test_Name' %in% colnames(net)) {
    message('ERROR: Please provide a network dataframe in Tidy format to generate plots.')
    return(NA)
  }

  net$r_squared = as.numeric(net$r_squared)

  if (!is.na(out_prefix)) {
    print(paste0("Saving: ", out_prefix, '.RsqrPerCondition.png'))
    png(paste0(out_prefix, '.RsqrPerCondition.png'), width=6 ,height=6, units="in", res=300)
  }
  plot = ggplot(net, aes(x=Test_Name, y=r_squared, fill=Test_Name)) +
    geom_boxplot() +
    scale_y_continuous(breaks=seq(-100,100, by=10)/100) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "none") +
    xlab("") + ylab(bquote(R^2))
  print(plot)
  if (!is.na(out_prefix)) {
    dev.off()
  }
}
#'Plots the distribution of p-values per similarity score
#'
#'@param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#'@param out_prefix
#'   The prefix of the output images files.
#'
#'@export
plotPvalPerScore = function(net, out_prefix=NA)  {
  if (!'Test_Name' %in% colnames(net)) {
    message('ERROR: Please provide a network dataframe in Tidy format to generate plots.')
    return(NA)
  }

  # Round off the similarity sccores
  mround = function(x, base) {
    base * floor(x/base)
  }

  net$Similarity_Score = as.factor(mround(net$Similarity_Score, 0.05))

  # Explore the distribution of p-values per correlation range.
  if (!is.na(out_prefix)) {
    print(paste0("Saving: ", out_prefix, '.PvalPerScore.png'))
    png(paste0(out_prefix, '.PvalPerScore.png'), width=6 ,height=6, units="in", res=300)
  }
  plot = ggplot(net, aes(x=Similarity_Score, y=p_value, fill=Similarity_Score)) +
    geom_boxplot() +
    scale_y_log10(breaks=c(1e-3, 1e-5, 1e-10, 1e-15, 1e-20)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "none") +
    xlab("Similarity Score") + ylab(bquote(-log^10~(p)))
  print(plot)
  if (!is.na(out_prefix)) {
    dev.off()
  }
}

#'Plots the distribution of R-squared values per similarity score
#'
#'@param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#'@param out_prefix
#'   The prefix of the output images files.
#'
#'@export
plotRsqrPerScore = function(net, out_prefix=NA)  {
  if (!'Test_Name' %in% colnames(net)) {
    message('ERROR: Please provide a network dataframe in Tidy format to generate plots.')
    return(NA)
  }

  # Round off the similarity sccores
  mround = function(x, base) {
    base * floor(x/base)
  }
  net$Similarity_Score = as.factor(mround(net$Similarity_Score, 0.05))
  net$r_squared = as.numeric(net$r_squared)

  if (!is.na(out_prefix)) {
    print(paste0("Saving: ", out_prefix, '.RsqrPerScore.png'))
    png(paste0(out_prefix, '.RsqrPerScore.png'), width=6 ,height=6, units="in", res=300)
  }
  plot = ggplot(net, aes(x=Similarity_Score, y=r_squared, fill=Similarity_Score)) +
    geom_boxplot() +
    scale_y_continuous(breaks=seq(-100,100, by=10)/100) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "none") +
    xlab("Similarity Score") + ylab(bquote(R^2))
  print(plot)
  if (!is.na(out_prefix)) {
    dev.off()
  }

}
#'Plots boxplots of the p-values per similarity score and condition.
#'
#'@param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#'@param out_prefix
#'   The prefix of the output images files.
#'
#'@export
plotPvalPerScoreAndCondition = function(net, out_prefix=NA) {
  if (!'Test_Name' %in% colnames(net)) {
    message('ERROR: Please provide a network dataframe in Tidy format to generate plots.')
    return(NA)
  }

  # Round off the similarity sccores
  mround = function(x, base) {
    base * floor(x/base)
  }
  net$Similarity_Score = as.factor(mround(net$Similarity_Score, 0.05))

  # Explore the distribution of p-values per correlation range and trait.
  if (!is.na(out_prefix)) {
    print(paste0("Saving: ", out_prefix, '.PvalPerScoreAndCondition.png'))
    png(paste0(out_prefix, '.PvalPerScoreAndCondition.png'), width=6 ,height=6, units="in", res=300)
  }
  plot = ggplot(net, aes(x=Similarity_Score, y=p_value, fill=Test_Name)) +
    geom_boxplot() +
    facet_wrap(~Test_Name) +
    scale_y_log10(breaks=c(1e-3, 1e-5, 1e-10, 1e-15, 1e-20)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(legend.position = "none") +
    xlab("Similarity Score") + ylab(bquote(-log^10~(p)))
  print(plot)
  if (!is.na(out_prefix)) {
    dev.off()
  }
}
#'Plots boxplots of the R-squared values per similarity score and condition.
#'
#'@param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#'@param out_prefix
#'   The prefix of the output images files.
#'
#'@export
plotRsqrPerScoreAndCondition = function(net, out_prefix=NA) {
  if (!'Test_Name' %in% colnames(net)) {
    message('ERROR: Please provide a network dataframe in Tidy format to generate plots.')
    return(NA)
  }

  # Round off the similarity sccores
  mround = function(x, base) {
    base * floor(x/base)
  }
  net$Similarity_Score = as.factor(mround(net$Similarity_Score, 0.05))
  net$r_squared = as.numeric(net$r_squared)

  if (!is.na(out_prefix)) {
    print(paste0("Saving: ", out_prefix, '.RsqrPerScoreAndCondition.png'))
    png(paste0(out_prefix, '.RsqrPerScoreAndCondition.png'), width=6 ,height=6, units="in", res=300)
  }
  plot = ggplot(net, aes(x=Similarity_Score, y=r_squared, fill=Test_Name)) +
    geom_boxplot() +
    facet_wrap(~Test_Name) +
    scale_y_continuous(breaks=seq(-100,100, by=10)/100) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(legend.position = "none") +
    xlab("Similarity Score") + ylab(bquote(R^2))
  print(plot)
  if (!is.na(out_prefix)) {
    dev.off()
  }
}
#'Plots a histogram of the p-values per condition
#'
#'@param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#'@param out_prefix
#'   The prefix of the output images files.
#'
#'@export
plotPvalHistPerCondition = function(net, out_prefix=NA) {
  if (!'Test_Name' %in% colnames(net)) {
    message('ERROR: Please provide a network dataframe in Tidy format to generate plots.')
    return(NA)
  }

  # How many rows does each trait contribute
  if (!is.na(out_prefix)) {
    print(paste0("Saving: ", out_prefix, '.PvalHistPerCondition.png'))
    png(paste0(out_prefix, '.PvalHistPerCondition.png'), width=6 ,height=6, units="in", res=300)
  }
  plot = ggplot(net, aes(x=p_value, color=Test_Name)) +
    geom_histogram(binwidth=0.5) +
    facet_wrap(~Test_Name) +
    scale_x_log10(breaks=c(1e-3, 1e-5, 1e-10, 1e-15, 1e-20)) +
    scale_y_log10() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    theme(legend.position = "none") +
    ylab("Count") + xlab(bquote(-log^10~(p)))
  print(plot)
  if (!is.na(out_prefix)) {
    dev.off()
  }
}
#'Plots a histogram of the R-squared values per condition
#'
#'@param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#'@param out_prefix
#'   The prefix of the output images files.
#'
#'@export
plotRsqrHistPerCondition = function(net, out_prefix=NA) {
  if (!'Test_Name' %in% colnames(net)) {
    message('ERROR: Please provide a network dataframe in Tidy format to generate plots.')
    return(NA)
  }

  # Round off the similarity sccores
  mround = function(x, base) {
    base * floor(x/base)
  }
  net$r_squared = as.numeric(net$r_squared)

  if (!is.na(out_prefix)) {
    print(paste0("Saving: ", out_prefix, '.RsqrHistPerCondition.png'))
    png(paste0(out_prefix, '.RsqrHistPerCondition.png'), width=6 ,height=6, units="in", res=300)
  }
  plot = ggplot(net, aes(x=r_squared, color=Test_Name)) +
    geom_histogram(binwidth=0.01) +
    facet_wrap(~Test_Name) +
    scale_y_continuous(breaks=seq(-100,100, by=10)/100) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    theme(legend.position = "none") +
    ylab("Count") + xlab(bquote(R^2))
  print(plot)
  if (!is.na(out_prefix)) {
    dev.off()
  }
}
#'Plots a bar graph of the number of edges per condition.
#'
#'@param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#'@param out_prefix
#'   The prefix of the output images files.
#'
#'@export
plotEdgeCountPerCondition = function(net, out_prefix=NA) {
  if (!'Test_Name' %in% colnames(net)) {
    message('ERROR: Please provide a network dataframe in Tidy format to generate plots.')
    return(NA)
  }

  # Plot the edge count per category
  if (!is.na(out_prefix)) {
    print(paste0("Saving: ", out_prefix, '.EdgeCountPerCondition.png'))
    png(paste0(out_prefix, '.EdgeCountPerCondition.png'), width=6 ,height=6, units="in", res=300)
  }
  plot = ggplot(net, aes(x=Test_Name, fill=Test_Name)) +
    geom_bar() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    theme(legend.position = "none") +
    ylab("Number of Edges") + xlab("")
  print(plot)
  if (!is.na(out_prefix)) {
    dev.off()
  }
}

#' Ranks the edges of a network and returns those ranks.
#'
#' Final ranks are determined by ordering p-values, R-squared and
#' similarity scores for each edge and assigning a rank to each edge
#' corresponding to the order of each value.  For example, edges with the
#' lowest p-value will receive a rank of 1.  The final rank is a weighted
#' sum of each of the three categories.
#'
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param by_condition
#'   If TRUE then ranks are calculated independently for each condition
#'   (i.e. value in the Test_Name column). If FALSE then ranks are
#'   calculated for the entire network.
#' @param sscore_weight
#'   The weight to apply to the Similarity Score in rankings.
#' @param pval_weight
#'   The weight to apply to the p-value in rankings.
#' @param rsqr_weight
#'   The weight to apply to the r-squared value in rankings.
#' @return
#'   An array of ranks. The order of the items in the ranks array corresponds
#'   to the order of edges in the provided network dataframe.
#'
#' @export
getEdgeRanks = function(net, by_condition = TRUE, sscore_weight = 1, pval_weight = 1, rsqr_weight = 1) {

  if (!'Test_Name' %in% colnames(net)) {
    message('ERROR: Please provide a network dataframe in Tidy format to generate plots.')
    return(NA)
  }

  if (by_condition == FALSE) {
    net = getRanks(net)
    valuations = (net$score_rank * sscore_weight) +
                 (net$pval_rank * pval_weight) +
                 (net$rsqr_rank * rsqr_weight)
    net$valuation = valuations

    # Now order the valuations
    unique_val = unique(valuations)
    ordered_val = unique_val[order(unique_val)]
    final_ranks = data.frame(valuation = ordered_val, rank = seq(1:length(ordered_val)))

    net = left_join(net, final_ranks, by=c('valuation' = 'valuation'))

  } else {
    for (condition in unique(net$Test_Name)) {
      csGCN = net[which(net$Test_Name == condition),]
      if ('rank' %in% colnames(csGCN)) {
        csGCN = csGCN[, -which(colnames(csGCN) %in% c("rank"))]
      }
      csGCN = getRanks(csGCN)

      # calculate the valuation for each edge.
      valuations = (csGCN$score_rank * sscore_weight) +
                  (csGCN$pval_rank * pval_weight) +
                  (csGCN$rsqr_rank * rsqr_weight)
      csGCN$valuation = valuations

      # Now order the valuations
      unique_val = unique(valuations)
      ordered_val = unique_val[order(unique_val)]
      final_ranks = data.frame(valuation = ordered_val, rank = seq(1:length(ordered_val)))

      csGCN = left_join(csGCN, final_ranks, by=c('valuation' = 'valuation'))
      net[which(net$Test_Name == condition), 'rank'] = csGCN$rank
    }
  }
  return(net$rank)
}

#' A helper function for the getEdgeRanks
#'
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @return
#'   A new network containing the p-value, r-squared and similarity score
#'   rankings.
getRanks = function(net) {
  net$ABS_Similarity_Score = abs(as.numeric(net$Similarity_Score))

  unique_scores = unique(net$ABS_Similarity_Score)
  ordered_scores = unique_scores[order(unique_scores, decreasing=TRUE)]
  score_ranks = data.frame(score = ordered_scores, score_rank = seq(1:length(ordered_scores)))

  unique_pvals = unique(net$p_value)
  ordered_pvals = unique_pvals[order(unique_pvals)]
  pval_ranks = data.frame(p_value = ordered_pvals, pval_rank = seq(1:length(ordered_pvals)))

  unique_rsqr = unique(net$r_squared)
  ordered_rsqr = unique_rsqr[order(unique_rsqr, decreasing=TRUE)]
  rsqr_ranks = data.frame(r_squared = ordered_rsqr, rsqr_rank = seq(1:length(ordered_rsqr)))

  net = left_join(net ,score_ranks, by=c('ABS_Similarity_Score' = 'score'))
  net = left_join(net, pval_ranks, by=c('p_value' = 'p_value'))
  net = left_join(net, rsqr_ranks, by=c('r_squared' = 'r_squared'))

  return (net)
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
#' @param sample_header
#'   The name of the column containing the sample names.
#'
#' @return
#'   A data frame containing the annotations in the order of the samples.
#'
#' @export
loadSampleAnnotations = function (annotation_file, sample_header="Sample") {
  # Read in the annotation file
  osa = read.table(annotation_file, sep="\t", header=TRUE, row.names=NULL, quote="", fill=TRUE)
  #sample_order = read.table(sample_order_file, colClasses=c('character'),
  #                          col.names=c('Sample'))
  #osa = merge(sample_order, sample_annots, by = "Sample", sort=FALSE)
  row.names(osa) = osa[,sample_header]

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
graphEdgeList = function(edge_indexes, net, osa = data.frame()) {

  g = graph.edgelist(as.matrix(net[edge_indexes, c('Source', 'Target')]), directed = F)
  E(g)$weight = abs(net[edge_indexes, 'Similarity_Score'])

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
  l = layout_nicely(g)
  plot(g, vertex.label.color="black", vertex.label.cex = vlc, vertex.size = vs,
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
#' @param gene1
#'   The name of the first gene in the pair.
#' @param gene2
#'   The name of the second gene in the pair.
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
plot3DPair = function(gene1, gene2, osa, ematrix, field,
                      colors = NA, samples = NA, highlight = c()) {

  rgl.open()
  rgl.bg(color = "white")


  if (!is.na(samples)) {
    x = as.factor(osa[samples, field])
    y = t(ematrix[gene1, ])
    z = t(ematrix[gene2, ])
  }
  else {
    x = as.factor(osa[, field])
    y = t(ematrix[gene1, ])
    z = t(ematrix[gene2, ])
  }

  size = rep(0.25, length(x))
  if (length(highlight) > 0) {
    size[highlight] = 0.5
  }

  if (is.na(colors) == TRUE) {
    plot3d(x, y, z, type = 's', size = size)
  } else {
    plot3d(x, y, z, type = 's', size = size, col = colors)
  }
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
#' @param samples
#'   Limit the plot to only the samples indexes provided.
#' @param highlight
#'   If set to TRUE then the samples in the edge cluster are drawn larger than the
#'   other samples in the plot.
#' @export
plot3DEdgeList = function(edge_indexes, osa, net, ematrix, field, label_field = NA,
                          samples = NA, highlight = TRUE) {


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
      z = osa[[field]][samples]
    }
    else {
      x = t(ematrix[source, ])
      y = t(ematrix[target, ])
      z = osa[, field]
    }

    size = rep(0.3, length(x))
    if (highlight & length(sample_indexes) > 0) {
      size[sample_indexes] = 0.75
    }

    # Build a color vector for the data.
    zcolors = rep('blue', length(z))
    if (class(z) != 'factor' & class(z) != 'numeric') {
      z = as.factor(z)
    }
    if (class(z) == 'factor') {
      if (nlevels(z) > 2 & nlevels(z) <= 8) {
        colors = data.frame(colors = brewer.pal(nlevels(z), "Dark2"), zval = unique(z), stringsAsFactors=FALSE)
        for (i in 1:length(z)) {
          zcolors[i] = colors$colors[which(colors$zval == z[i])]
        }
      }
      if (nlevels(z) > 8) {
        colfunc <- colorRampPalette(c("darkblue", "cyan"))
        colors = data.frame(colors = colfunc(length(sort(unique(z)))),
                            zval = sort(unique(z)),
                            stringsAsFactors=FALSE)
        for (i in 1:length(z)) {
          if (is.na(z[i])) {
            zcolors[i] = '#000000'
          } else {
            zcolors[i] = colors$colors[which(colors$zval == z[i])]
          }
        }
      }
      if (nlevels(z) == 2) {
        colors = data.frame(colors = brewer.pal(3, "Dark2")[1:2],
                            zval = unique(z),
                            stringsAsFactors=FALSE)
        for (i in 1:length(z)) {
          if (is.na(z[i])) {
            zcolors[i] = '#000000'
          } else {
            zcolors[i] = colors$colors[which(colors$zval == z[i])]
          }
        }
      }
    }
    if (class(z) == 'numeric') {
      # Use a blue gradient
      colfunc <- colorRampPalette(c("darkblue", "cyan"))
      colors = data.frame(colors = colfunc(length(sort(unique(z)))),
                          zval = sort(unique(z)),
                          stringsAsFactors=FALSE)
      for (i in 1:length(z)) {
        if (is.na(z[i])) {
          zcolors[i] = '#000000'
        }
        else {
          zcolors[i] = colors$colors[which(colors$zval == z[i])]
        }
      }
    }

    main = paste(source, 'vs', target, 'Edge:', i)
    plot3d(x, y, z, type = 's', main = main, size = size, col = zcolors, xlab=source, ylab=target, zlab=field, axes=FALSE)
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
#'   The label for the x-axis
#' @param colfield
#'   The field in the OSA that should be used to color the points. Each
#'   category in the field recieves a unique color.
#' @param show.legend
#'   Set to TRUE to have a legend appear on the figure. Defaults to TRUE.
#' @param fig_title
#'   The text to use for the title of the figure.
#' @param highlight
#'   A vector of sample indexes to highlight in the plot.
#' @export
plotGene = function(gene, ematrix, osa, field, xlab = field, colfield = field,
                    show.legend=TRUE, fig_title = NA, highlight=c()) {

  condition = data.frame(c = osa[[field]], c2 = osa[[colfield]])
  row.names(condition) = row.names(osa)
  expdata = merge(t(ematrix[gene,]), condition, by="row.names")
  colnames(expdata) = c('Sample', 'y', 'x', 'z')
  row.names(expdata) = make.names(row.names(expdata))

  expdata$size = 1
  if (length(highlight) > 0) {
    expdata$size = 0.25
    sample_names = colnames(ematrix)[highlight]
    expdata$size[which(expdata$Sample %in% sample_names)] = 1
  }

  expplot = ggplot(expdata, aes(x, y, color=z)) +
    geom_point(size=expdata$size, show.legend = show.legend) +
    xlab(xlab) + ylab(gene) +
    theme(legend.title = element_blank())
  if (!is.numeric(expdata$z)) {
    expplot = expplot + scale_color_brewer(palette="Dark2")
  }
  if (!is.null(fig_title)) {
    expplot = expplot + ggtitle(fig_title)
  }
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
#' @param fig_title
#'   A title to give the plot. Default = NULL
#' @export
plot2DEdgeList = function(edge_indexes, osa, net, ematrix,
                          field = NA, samples = NA, highlight = TRUE,
                          fig_title = NULL, show.legend=TRUE, legend.title = field) {
  j = 1
  done = FALSE;
  while (!done) {
    i = edge_indexes[j]
    source = net$Source[i]
    target = net$Target[i]
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
    if (!highlight) {
      size = rep(1, length(x))
    }
    if (highlight & length(sample_indexes) > 0) {
      size[sample_indexes] = 1.5
    }

    coexpdata = data.frame(source = x, target = y, category = condition, size = size)
    colnames(coexpdata) = c('x', 'y', 'category', 'size')
    coexpdata = coexpdata[complete.cases(coexpdata),]
    coexpplot = ggplot(coexpdata, aes(x, y, color=category)) +
      geom_point(size=coexpdata$size, show.legend = show.legend) +
      xlab(source) + ylab(target) + labs(colour = legend.title)
    if (!is.null(fig_title)) {
      coexpplot = coexpplot + ggtitle(fig_title)
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

#' Plots a 2D scatterplot for a given pair of genes
#'
#' @param gene1
#'   The name of the first gene in the pair.
#' @param gene2
#'   The name of the second gene in the pair.
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
#' @param fig_title
#'   A title to give the plot. Default = NULL
#' @export
plot2DPair = function(gene1, gene2, osa, net, ematrix,
                      field = NA, fig_title = NA, show.legend=TRUE) {

    # Get the gene expression and order it by the osa sample names.
    x = t(ematrix[gene1, ])
    y = t(ematrix[gene2, ])

    # Get the conditions of the field specified by the user.
    condition = osa[colnames(ematrix),field]
    size = 1

    # Build a datafame suitable for ggplot
    coexpdata = data.frame(source = x, target = y, category = condition, size = size)
    colnames(coexpdata) = c('x', 'y', 'category', 'size')
    coexpdata = coexpdata[complete.cases(coexpdata),]

    # If we have an edge in the network for this pair then we'll change
    # the size of the points to match
    edge = which((net$Source == gene1 & net$Target == gene2) | (net$Source == gene2 & net$Target == gene1))
    if (edge & 'Samples' %in% colnames(net)) {
      print(paste("Network Edge #", edge, sep=""))
      coexpdata$size = 0.25
      samples = getEdgeSamples(edge, net)
      coexpdata$size[which(row.names(coexpdata) %in% colnames(ematrix)[samples])] = 1
    }

    row.names(coexpdata) = make.names(row.names(coexpdata))

    coexpplot = ggplot(coexpdata, aes(x, y, color=category)) +
      geom_point(size=coexpdata$size, show.legend = show.legend) +
      xlab(gene1) + ylab(gene2) + labs(colour=field)
    if (!is.numeric(coexpdata$category)) {
      coexpplot = coexpplot + scale_color_brewer(palette="Dark2")
    }
    if (!is.null(fig_title)) {
      coexpplot = coexpplot + ggtitle(fig_title)
    }
    print(coexpplot)
    r = cor(coexpdata$x, coexpdata$y, method="spearman", use="complete.obs")
    print(r)
    return(coexpplot)
}

#' Plots four images for each edge to explore their relationship.
#'
#' @param gene1
#'   The name of the first gene in the pair.
#' @param gene2
#'   The name of the second gene in the pair.
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
#'   into separate layers in the pairwise scatterplot
#' @param field2
#'   The field in the sample annotation matrix by which to order points on
#'   the y-axis of the bottom two gene expression plots. The x-axis is
#'   ordered by expression level.
#' @export
plot2DEdgeReport <- function(edge_indexes, osa, net, ematrix,
                             field = NA, field2 = field, show.legend=TRUE) {
  j = 1
  done = FALSE;
  while (!done) {
    i = edge_indexes[j]
    source = net[i, 'Source']
    target = net[i, 'Target']

    plot2DPairReport(source, target, osa, net, ematrix, field, field2, show.legend)

    if (length(edge_indexes) == 1) {
      return
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

#' Plots four images for a given gene pair to explore their relationship.
#'
#' @param gene1
#'   The name of the first gene in the pair.
#' @param gene2
#'   The name of the second gene in the pair.
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
#'   into separate layers in the pairwise scatterplot
#' @param field2
#'   The field in the sample annotation matrix by which to order points on
#'   the y-axis of the bottom two gene expression plots. The x-axis is
#'   ordered by expression level.
#' @export
plot2DPairReport <-function(gene1, gene2, osa, net, ematrix,
                            field = NA, field2 = field, show.legend=TRUE) {

  # If we have an edge in the network for this pair then we'll change
  # the size of the points to match
  groups = osa[c(field2, field2)]
  edge = which((net$Source == gene1 & net$Target == gene2) | (net$Source == gene2 & net$Target == gene1))
  samples=c()
  if (edge & 'Samples' %in% colnames(net)) {
    samples = getEdgeSamples(edge, net)
    groups = rep('Out', length(colnames(ematrix)))
    groups[samples] = "In"
    groups = data.frame(groups)
    row.names(groups) = colnames(ematrix)
  }

  FigYa = plot2DPair(gene1, gene2, osa, net, ematrix, field, fig_title='(a)', show.legend=TRUE)
  FigYc = plotGene(gene1, ematrix, osa, field2, colfield=field, show.legend=TRUE, fig_title='(c)', highlight=samples)
  FigYd = plotGene(gene2, ematrix, osa, field2, colfield=field, show.legend=TRUE, fig_title='(d)', highlight=samples)



  expdata = merge(t(ematrix[c(gene1,gene2),]), groups, by="row.names")
  cols = c('Sample', gene1, gene2, 'Edge')
  colnames(expdata) = cols
  expdata = expdata[, cols]
  expdata = melt(expdata)
  colnames(expdata) = c('Sample', 'Edge', 'Gene', 'Expression')
  FigYb <- ggplot(arrange(expdata, Edge), aes(x=Edge, y=Expression, fill=Edge)) +
    geom_violin(trim=FALSE, show.legend=FALSE) + ggtitle('(b)') +
    geom_boxplot(width=0.1, show.legend=FALSE) +
    scale_fill_brewer(palette="Paired") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(. ~ Gene)

  layout <- matrix(seq(1,100), ncol = 10, nrow = 10)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
  print(FigYa, vp = viewport(layout.pos.row = c(1:6), layout.pos.col = c(1:5)))
  print(FigYb, vp = viewport(layout.pos.row = c(1:6), layout.pos.col = c(6:10)))
  print(FigYc, vp = viewport(layout.pos.row = c(7:10), layout.pos.col = c(1:5)))
  print(FigYd, vp = viewport(layout.pos.row = c(7:10), layout.pos.col = c(6:10)))
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
findLinkedCommunities = function(net, file_prefix="net", module_prefix = 'M',
                                 hcmethod = 'complete', meta = TRUE, ignore_inverse = TRUE) {
  new_net = net
  new_net$Module = NA
  sim_col = 'Similarity_Score'

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
  write.table(new_net, file=paste(file_prefix, "gcn.lcm.txt", sep=".") ,sep="\t", row.names=FALSE, append=FALSE, quote=FALSE)

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


