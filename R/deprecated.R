#' Performs a fisher's exact test on a set of samples.
#'
#' An annotation field and a specific category are provided. This
#' test reports if the category is enriched within the cluster.
#' It does not indiciate if the category is singificantly more
#' prominent in the cluster. For that, try the sampleClusterBTest
#' function.
#' @param category
#'   The annotation category to be used for testing. It must be a valid
#'   category in the field if the osa identified by the field argument.
#' @param field
#'   The field in the osa variable on which enrichment will be performed.
#' @param i
#'   The index of the edge in the network
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param osa
#'   The sample annotation matrix. One column must contain the header 'Sample'
#'   and the remaining colums correspond to an annotation type.  The rows
#'   of the anntation columns should contain the annotations.
#' @param ematrix
#'   The expression matrix data frame.
#' @param samples
#'   Limit the analysis to only the samples indexes provided.
performClusterFishersTest = function(category, field, i, net, osa, ematrix, verbose=FALSE) {

  # Convert the sample string to a numerical vector.
  cluster_samples = getEdgeSamples(i, net)

  # Get osa elements for samples that belong to the cluster and those that don't.
  osa_cat_indexes = which(osa$Sample %in% names(ematrix)[cluster_samples])
  osa_out_indexes = which(!osa$Sample %in% names(ematrix)[cluster_samples])

  #  Contingency matrix for each category in a cluster:
  #
  #                   Is Cat  Not Cat    Totals
  #                  ------------------
  #  In Cluster      |  n11   |   n12   | n1p
  #  Not in Cluster  |  n21   |   n22   | n2p
  #                  ------------------
  #  Totals             np1       np2     npp
  #
  n11 = length(which(osa[[field]][osa_cat_indexes] == category))
  n12 = length(which(osa[[field]][osa_cat_indexes] != category))
  n21 = length(which(osa[[field]][osa_out_indexes] == category))
  n22 = length(which(osa[[field]][osa_out_indexes] != category))

  contmatrix = matrix(
    as.numeric(c(n11, n21, n12, n22)),
    nr=2,
    dimnames = list(
      In_Cluster = c("Yes", "No"),
      Is_Category = c("Yes", "No")
    )
  )

  if (verbose) {
    print(contmatrix)
  }
  res = fisher.test(contmatrix, alternative = "greater")
  p.value = res$p.value
  return(p.value)
}



#' Performs significant testing of a single edge in the network for a set of annotation categories.
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
#' @param ematrix
#'   The expression matrix data frame.
#' @param field
#'   The field in the osa variable on which enrichment will be performed.
#' @param test
#'   The statistical test to perform: 'fishers' or 'binomial'
#' @param category
#'   By default this function will cylce through all categories of the
#'   'field' argument. But to limit the analysis to just one category set it here.
#' @param t1_psucc
#'   For the bionomial test only.
#'   The probabilty of success for the first binomial test.  This value
#'   indicates the percentage of samples in the cluster that must
#'   belong to the category.
#' @param t2_psucc
#'   For the binomial test only.
#'   The probability of success for the second binomial test. This value
#'   indicates the percentage of category samples that must belong
#'   to the cluster.
#' @examples
#'
analyzeEdgeCat = function(i, osa, net, ematrix, field, test = 'binomial',
                          category = NA, t1_psucc = 0.25, t2_psucc = 0.75) {

  sample_types = as.character(osa[[field]])
  num_samples = length(sample_types)

  categories = unique(osa[[field]])
  if (!is.na(category)) {
    categories = c(category)
  }

  pvals = sapply(categories, performClusterBinomialTest, field, i, net, osa, ematrix, t1_psucc, t2_psucc, FALSE)
  names(pvals) = categories
  return(pvals);

}

#' Performs significant testing of each edge in the network for a set of annotation categories.
#'
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param osa
#'   The ordered sample annotation data frame.
#' @param ematrix
#'   The expression matrix data frame.
#' @param field
#'   The field in the osa variable on which enrichment will be performed.
#' @param test
#'   The statistical test to perform.  This is 'fishers' for an enrichment
#'   test.
#' @param correction.
#'   The method to apply for multiple testing correction. Valid values are identical
#'   to those available to the p.adjust function.  The default is to
#'   apply 'hochberg' correction.
#' @param progressBar
#'   Set to FALSE to repress progress bar
#' @param category
#'   By default this function will cylce through all categories of the
#'   'field' argument. But to limit the analysis to just one category set it here.
#' @param t1_psucc
#'   For the bionomial test only.
#'   The probabilty of success for the first binomial test.  This value
#'   indicates the percentage of samples in the cluster that must
#'   belong to the category.
#' @param t2_psucc
#'   For the binomial test only.
#'   The probability of success for the second binomial test. This value
#'   indicates the percentage of category samples that must belong
#'   to the cluster.
#' @param nthreads
#'   The number of threads (CPU cores/threads) used to analyze the network.
#'
analyzeNetCat = function(net, osa, ematrix, field, test = 'binomial',
                         correction = 'hochberg', progressBar = TRUE,
                         category = NA, t1_psucc = 0.25, t2_psucc=0.75) {

  # Get the list of categories to analyze.
  sample_types = as.character(osa[[field]])
  categories = unique(sample_types)
  if (!is.na(category)) {
    categories = c(category)
  }

  # Intialize the progress bar
  if (progressBar){
    pb = txtProgressBar(min = 0, max = nrow(net), style = 3)
  }

  # Iterate through the rows of the network
  for (i in 1:nrow(net)) {

    if (progressBar){
      setTxtProgressBar(pb, i)
    }

    # Add in new columns for each category.
    for (label in categories) {
      subname = paste(field, label, sep='_')
      net[i, subname] = NA
    }

    # Calculate the probability that this edge is a result of any specific
    # known experimental condition (i.e. category) for the given field.
    p.vals = analyzeEdgeCat(i, osa, net, ematrix, field, test = test,
                            category, t1_psucc, t2_psucc)

    # Add the p-values into the edge.
    for (label in names(p.vals)) {
      subname = paste(field, label, sep='_')
      net[i, subname] = p.vals[label]
    }
  }
  if (progressBar){
    close(pb)
  }

  # Perform multiple testing correction on the p-values. This may not
  # be necessary because it overly reduces p-values from the binomial
  # test, but for backwards compatiblity, let's leave it for fishers:
  if (test == 'fishers') {
    for (category in categories) {
      subname = paste(field, category, sep='_')
      net[subname] = p.adjust(as.numeric(unlist(net[subname])), method=correction)
    }
  }

  return(net);
}
#' Performs linear regression of a quantitative traits against a single edge in the network.
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
#' @param samples
#'   Limit the analysis to only the samples indexes provided.
#'
#'
#' @examples
#'
analyzeEdgeQuant = function(i, osa, net, field, samples = c()) {
  source = net[i, 'Source']
  target = net[i, 'Target']

  # Convert the sample string to a numerical vector.
  edge_samples = getEdgeSamples(i, net)
  if (length(samples) > 0) {
    edge_samples = samples
  }

  # Build the expression vectors using only the samples provided.
  x = as.numeric(t(ematrix[source, edge_samples]))
  y = as.numeric(t(ematrix[target, edge_samples]))
  z = as.numeric(as.factor(osa[edge_samples, field]))

  # If our z-dimension only has one value then there's no need to
  # perform a linear regression, so just return.
  if (dim(table(z)) == 1) {
    return(list(
      p = NA,
      roccm = NA,
      model = NA
    ))
  }

  # Use linear regression to obtain a p-value for the association.
  model = lm(z ~ x + y, data=data.frame(x=x, y=y, z=z))
  s = summary(model)
  pval = NA
  roccm = NA
  if (dim(s$coefficients)[1] > 1) {
    # Choose the smallest p-value of the x and y slopes
    pval = min(s$coefficients[2,4], s$coefficients[3,4])

    # Get the rate of change of the mean of y + x
    roccm = s$coefficients[2,1]
  }

  return(list(
    p = pval,
    roccm = roccm,
    model = model
  ))
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
#'   apply 'hochberg' correction.
#' @param samples
#'   Limit the annalysis to only the samples indexes provided.
#' @param progressBar
#'   Set to FALSE to repress progress bar
#'
analyzeNetQuant = function(net, osa, field, correction = 'hochberg',
                           samples = c(), progressBar = TRUE) {

  # Intialize the progress bar
  if (progressBar) {
    pb <- txtProgressBar(min = 0, max = nrow(net), style = 3)
  }
  for (i in 1:nrow(net)) {

    if (progressBar) {
      setTxtProgressBar(pb, i)
    }

    # Add in new columns for the category.
    net[i, field] = NA
    net[i, paste(field, '_roccm', sep='')] = NA

    result = analyzeEdgeQuant(i, osa, net, field, samples)
    net[i, field] = result[['p']]
    net[i, paste(field, '_roccm', sep='')] = result[['roccm']]
  }
  if (progressBar){
    close(pb)
  }


  # Perform multiple testing correction on the p-values
  #net[field] = p.adjust(as.numeric(unlist(net[field])), method=correction)

  return(net);
}


#' Performs a fisher's exact test on a set of samples.
#'
#' An annotation field and a specific category are provided. This
#' test reports if the category is enriched within the cluster.
#' It does not indiciate if the category is singificantly more
#' prominent in the cluster. For that, try the sampleClusterBTest
#' function.
#' @param category
#'   The annotation category to be used for testing. It must be a valid
#'   category in the field if the osa identified by the field argument.
#' @param field
#'   The field in the osa variable on which enrichment will be performed.
#' @param i
#'   The index of the edge in the network
#' @param net
#'   A network data frame containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param osa
#'   The sample annotation matrix. One column must contain the header 'Sample'
#'   and the remaining colums correspond to an annotation type.  The rows
#'   of the anntation columns should contain the annotations.
#' @param ematrix
#'   The expression matrix data frame.
#' @param t1_psucc
#'   The probabilty of success for the first binomial test.  This value
#'   indicates the percentage of samples in the cluster that must
#'   belong to the category.
#' @param t2_psucc
#'   The probability of success for the second binomial test. This value
#'   indicates the percentage of category samples that must belong
#'   to the cluster.
performClusterBinomialTest = function(category, field, i, net, osa, ematrix,
                                      t1_psucc = 0.25, t2_psucc = 0.75, verbose=FALSE) {

  # Convert the sample string to a numerical vector.
  cluster_samples = getEdgeSamples(i, net)

  # Get the cluster size
  cluster_size = length(cluster_samples)

  # How many samples belong to this category
  num_category = length(which(osa[[field]] == category))
  num_not_category = length(which(osa[[field]] != category))

  # Get the indexes in the OSA of the samples that are in the cluster and not in
  # the cluster.
  cluster_osa_indexes = which(osa$Sample %in% names(ematrix)[cluster_samples])
  non_cluster_osa_indexes = which(!osa$Sample %in% names(ematrix)[cluster_samples])

  # Get the number of samples of the category that are in and not in the cluster.
  num_category_in_cluster = length(which(osa[[field]][cluster_osa_indexes] == category))
  num_category_not_in_cluster = num_category - num_category_in_cluster

  # Get the number of samples not of the category that are in and not in the cluster.
  num_other_in_cluster = length(which(osa[[field]][cluster_osa_indexes] != category))
  num_other_not_in_cluster = num_not_category - num_other_in_cluster

  # Test #1
  # successes = number of non-category samples in the cluster
  # failures = number of non-category samples not in the cluster
  # Ho: successes >= 0.15
  # Ha: successes < 0.15
  p1.pval = 0
  if (num_not_category > 0) {
    p1 = binom.test(num_other_in_cluster, num_not_category, p=t1_psucc, alternative='less')
    p1.pval = p1$p.value
  }

  # Test #2
  # successes = number of category samples in the cluster
  # failures = number of category samples out of the cluster
  # Ho: successes = 0.85
  # Ha: successes > 0.85
  p2.pval = 0
  if (num_category > 0) {
    p2 = binom.test(num_category_in_cluster, num_category, p=t2_psucc, alternative='greater')
    p2.pval = p2$p.value
  }

  # Set the maximum p-value
  # TODO: perhaps we should employ a method of combining p-values rather
  # than just using the minimum.
  # https://academic.oup.com/bioinformatics/article/32/17/i430/2450768
  p.value = max(p1.pval, p2.pval)

  if (verbose) {
    print(p1)
    print(p2)
  }

  # If all tests pass then return the highest binomial
  # p-value.
  return(p.value)
}


analyzeNetDiffQuant = function(net, osa, field, model_samples = NA, test_samples,
                               progressBar = FALSE) {

  # Add in new columns for each category.
  net2 = net
  column_name = paste(field, 'diff', paste(test_samples, collapse='_'), sep="-")
  net2[column_name] = NA

  if(progressBar){pb <- txtProgressBar(min = 0, max = nrow(net), style = 3)}
  for (i in 1:nrow(net)) {
    if (progressBar){setTxtProgressBar(pb, i)}
    result = analyzeEdgeDiffQuant(i, osa, net, field, model_samples, test_samples)
    net2[i, column_name] = median(result)
  }
  if (progressBar){close(pb)}

  return(net2);
}

analyzeEdgeDiffQuant = function(i, osa, net, field, model_samples = c(), test_samples = c()) {
  source = net[i, 'Source']
  target = net[i, 'Target']

  # Convert the sample string to a numerical vector.
  edge = net[i,]
  edge_samples = getEdgeSamples(edge$Samples)
  if (length(samples) > 0) {
    edge_samples = samples
  }

  # Build the expression vectors using only the edges of the edge.
  x = t(ematrix[source, ])
  y = t(ematrix[target, ])
  z = as.numeric(as.factor(osa[, field]))

  # Use linear regression to obtain a model for the samples the model samples group.
  model = lm(y + x ~ z, data=data.frame(x=x, y=y, z=z), subset = edge_samples)

  # Now get the difference for the samples in the test_sample group
  obs = t(ematrix[source, test_samples]) + t(ematrix[target, test_samples])
  colnames(obs) = 'y + x'
  exp = predict(model, data.frame(z = z[test_samples]), interval = 'predict')
  res = (obs - exp[, 'lwr']) / (exp[, 'upr'] - exp[, 'lwr']) * 2 - 1
  return (res)
}
