#' T-test Edge
#' 
#' Performs a weltch's t-test on the two genes of an edge. Returns 1 if there is not enough power to do test,
#' otherwise will return t-test stat for both genes
#' 
#' @param i
#'   The index of the edge in the network
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#'   This function is intended to be run on a network that has already been
#'   reduced by analyzeNetCat.
#' @param ematrix
#'   The expression matrix data frame.
#'   
#' @export
TtestEdgeCat = function(i, net, ematrix){

  # Gets the name of the Target and the Source gene
  Source = net[i,]$Source
  Target = net[i,]$Target
  pvals = c(1,1) # Assume that this fails, modify from there
  names(pvals) = c("Ttest_Source","Ttest_Target")
  k = 1

  # Convert the sample string to a numerical vector.
  cluster_indexes = getEdgeSamples(i, net)
  non_cluster_indexes = getNonEdgeSamples(i, net)

  min_samples = 2 # This parameter needs to be better thought out I think
  # It is possible that we could make it so the sample sizes must be kind of equal 
  # for it to allow the t-test to happen?
  
  # First check to make sure that there are at least 2
  # samples in each cluster, then Perform t-test for each gene
  if(length(cluster_indexes) >= min_samples & length(non_cluster_indexes) >= min_samples)
  {
    for(j in c(Source, Target))
    {
      # Test to see if gene breaks the null hypothesis
      Ttest_result = t.test(ematrix[j, cluster_indexes], ematrix[j, non_cluster_indexes],
                            mu=0,alt="two.sided",conf=0.95,var.eq=F,paired=F)
      pvals[k] = Ttest_result$p.value
      k = k + 1
    }
  }

  # return the p.values
  return(pvals);
}


  

#' T-test Net
#' 
#' Performs a weltch's t-test on a network
#' 
#' 
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#'   This function is intended to be run on a network that has already been
#'   reduced by analyzeNetCat.
#' @param ematrix
#'   The expression matrix data frame.
#' @param progressBar
#'   Set to FALSE to repress progress bar
#' 
#' @export
TtestNetCat = function(net, ematrix, progressBar = TRUE){
  # Add in new columns for the Ttest
  net2 = net
  net2["Ttest_Source"] = NA
  net2["Ttest_Target"] = NA
  
  # Intialize the progress bar
  if (progressBar){
    pb <- txtProgressBar(min = 0, max = nrow(net), style = 3)
  }
  
  # Iterate through the rows of the network
  for (i in 1:nrow(net)) {
    if (progressBar){
      setTxtProgressBar(pb, i)
    }
    # print(i)
    # Calculate the probability that each edge meets the assumption of the Ttest
    # for each gene.
    # The mean of the non-group should be independent of the group.
    p.value = TtestEdgeCat(i, net, ematrix)
    
    #print(p.value)
    net2[i, "Ttest_Source"] = p.value["Ttest_Source"]
    net2[i, "Ttest_Target"] = p.value["Ttest_Target"]
  }
  if (progressBar){
    close(pb)
  }
  
  return(net2);
}
