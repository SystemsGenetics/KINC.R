# Perform t-test

#' Performs a weltch's t-test
#' @param i
#'   The index of the edge in the network
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param ematrix
#'   The expression matrix data frame.
#' @param num_not_category
#'   The number of samples that are not in the cluster
#' @param num_category
#'   The number of samples that are in the cluster
#' @param cluster_osa_indexes
#'   indexes of the samples that are in the cluster 
#' @param non_cluster_osa_indexes
#'   indexes of the samples that are not in the cluster
performClusterTTest = function(i, net, ematrix, num_not_category, num_category, 
                               cluster_osa_indexes, non_cluster_osa_indexes) 
{
  # Gets the total samples that are in the cluster and non cluster
  total_samples = num_category + num_not_category
  
  # Gets the name of the Target and the Source gene
  Target = net[i,]$Target
  Source = net[i,]$Source
  
  #'
  #' Power Analysis
  #' Determine if we have the power to distinguish if there is a difference with our number of samples
  #' If we do not have the power for one or both of the genes, we can not determine if it is signicantly different
  #' return 1 if this is true
  #' 
  Source_cluster_mean = mean(as.numeric(ematrix[Source,cluster_osa_indexes]), na.rm = TRUE)
  Source_cluster_sd = sd(as.numeric(ematrix[Source,cluster_osa_indexes]), na.rm = TRUE)
  Source_non_cluster_mean = mean(as.numeric(ematrix[Source,non_cluster_osa_indexes]), na.rm = TRUE)
  Source_non_cluster_sd = sd(as.numeric(ematrix[Source,non_cluster_osa_indexes]), na.rm = TRUE)
  
  # This is to determine the power of the test with the number of samples:
  Source_Power = AB_t2n(N = (total_samples), percent_B = num_not_category/total_samples, mean_diff = abs(Source_cluster_mean - Source_non_cluster_mean)
                        , sd_A = Source_cluster_sd, sd_B = Source_non_cluster_sd, sig_level = 0.001, alternative = 'two_sided')
  
  
  Target_cluster_mean = mean(as.numeric(ematrix[Target,cluster_osa_indexes]), na.rm = TRUE)
  Target_cluster_sd = sd(as.numeric(ematrix[Target,cluster_osa_indexes]), na.rm = TRUE)
  Target_non_cluster_mean = mean(as.numeric(ematrix[Target,non_cluster_osa_indexes]), na.rm = TRUE)
  Target_non_cluster_sd = sd(as.numeric(ematrix[Target,non_cluster_osa_indexes]), na.rm = TRUE)
  
  # This is to determine the power of the test with the number of samples:
  Target_Power = AB_t2n(N = (total_samples), percent_B = num_not_category/total_samples, mean_diff = abs(Target_cluster_mean - Target_non_cluster_mean)
                        , sd_A = Target_cluster_sd, sd_B = Target_non_cluster_sd, sig_level = 0.001, alternative = 'two_sided')
  
  # Determine if we even have the power to attempt a t-test. If we do not have the power, then return 1
  # Also must check to make sure that we actually have a power value for each
  if (min(Target_Power$power, Source_Power$power) < 0.8 || is.na(Target_Power$power) || is.na(Source_Power$power))
  {
    return(1)
  }
  else
  {
    ###################################################
    # Check to see if min non-cluster size is satisfied
    # Get Names of Source and Target Gene for this edge
    Source = net[i,]$Source
    Target = net[i,]$Target
    
    # Test to see if Source or Target gene break the null hypothesis
    s = t.test(ematrix[Source,cluster_osa_indexes], ematrix[Source,non_cluster_osa_indexes],mu=0,alt="two.sided",conf=0.95,var.eq=F,paired=F)
    s.pval = s$p.value
    t = t.test(ematrix[Target,cluster_osa_indexes], ematrix[Target,non_cluster_osa_indexes],mu=0,alt="two.sided",conf=0.95,var.eq=F,paired=F)
    t.pval = t$p.value
    
    # Detemine higher p-value.
    p.value = max(s.pval, t.pval)
    
    # return the p.value
    return(p.value)
    
  }
}