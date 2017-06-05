#' Performs the RMT analysis on a gene similarity matrix.
#'
#' During Random Matrix Theory (RMT) analysis, a subset of the similarity 
#' matrix is tested to see how closely the nearest neighbor spacing 
#' distribution (NNSD) of its eigenvalues. A NNSD that appears
#' Guassian is indicative of a random matrix.  Naturally occuring networks
#' do not exhibit random behavior and their NNSD appears Poisson.  This 
#' function iterates through a successively decreasing similarity score
#' values, creates a subset of the similarity matrix containing only genes
#' with similarity values at or above the test score and performs a Chi-square
#' test on its NNSD.  Prior to testing, eigenvalues are spectrally unfolded
#' and then fit with a spline.  The NNSD is essentially a historgram with
#' 60 bins.  Therefore, the Chi-square test has 60 degrees of freedom and
#' the test statistic, 99.607, corresponding to a p-value of 0.001. Thus
#' A Chi-square test less than or equal to 99.607 indicates that the 
#' similarity submatrix appears Poisson.
#'
#' 
#' @param net
#'   A network dataframe containing the KINC-produced network.  The loadNetwork
#'   function imports a dataframe in the correct format for this function.
#' @param start
#'   The similarity score at which to start RMT analysis. The RMT will 
#'   iterate through successively lower similarity scores. By default
#'   it performs the ChiSquare test of a reduced matrix that only contains
#'   genes with correlation at or above this start value. Default is 0.99.
#' @param end
#'   The similarity score at which the RMT analysis will not proceed further.
#'   Default is 0.75
#' @param step
#'   The size that the similarity score is reduced at each iteration tested.
#'   Default is -0.001.
#' @param min.eigens
#'   To perform a Chi-squared test of the NNSD, the eigenvalues of the
#'   reduces similarity matrix must first be calculated.  This argument
#'   sets the limit for the number of eigenvalues that must be present in
#'   order to perform the test. This roughly corresponds to the number of
#'   genes in the reduced similarity matrix.  The actual number of 
#'   eigenvalues may be smaller than the number of genes if duplicates
#'   are present.  Duplicates are removed. Default is 100.
#' @param save.plot
#'   Set to TRUE to save a plot at each iteration. 
#' @export
#' @examples
#'   
RMT = function(net, start = 0.99, end = 0.75, step = -0.001, min.eigens = 100, save.plot = FALSE) {
  for (k in seq(start, end, step)) {
    sig = net[net$Similarity > k | net$Similarity < -k,]
    names=unique(c(as.vector(unique(sig$Source)), as.vector(unique(sig$Target))))
    n = length(names)
    num_rows = nrow(sig)

    if (n < 100) {
      print(paste("Skipping ", k, ". Not enough nodes: ", n, sep=""))
      next
    }

    # create the similarity matrix
    m = matrix(nrow=n, ncol=n)
    colnames(m) = names
    rownames(m) = names
    diag(m) = 1
    for (i in 1:num_rows) {
      cat(paste(i, "of", num_rows, "\r", sep=" "))
      g1 = as.character(sig[i,]$Source)
      g2 = as.character(sig[i,]$Target)
      m[g1,g2] = sig[i,]$Similarity
      m[g2,g1] = sig[i,]$Similarity
    }
    m[is.na(m)] = 0

    ##
    # Unfolding
    #
    # https://books.google.com/books?id=29UIY9EZTzYC&pg=PA715&lpg=PA715#v=onepage&q&f=false
    # Step 1: Calculate the eigenvalues of the matrix of order n
    # get the eigenvalues, remove duplicates, and order.
    e = eigen(m)$values
    # Set values less than 0.00001 to zero. (From KINC code)
    e[which(abs(e) < 0.000001)] = 0
    # Remove duplicates (those that are within 0.000001 of each other
    eo = unique(e[order(e)])
    v = eo[c(1, which(diff(eo) > 0.000001) + 1)]

    if (length(v) < min.eigens) {
      next;
    }

    # Calculate the average Chi-square
    chi=c()
    for (i in seq(40, 10, ,31)) {
      chi[length(chi)+1] = getNNSDPaceChiSquare(k, v, i, save.plot)
    }
    print(paste("th: ", sprintf("%.5f", k), ", Avg Chi: ", sprintf("%.2f", mean(chi)), sep=""));
  }
}

#
#
#
getNNSDPaceChiSquare = function(k, v, pace = 10, save.plot = FALSE) {
  # Calculate the GOE curve.
  goe=c()
  i=0
  for (s in seq(0,3, 0.01)) {
    i = i + 1
    goe[i] = (32/(pi^2))*s^2*exp(-(4/pi)*s^2)
  }

  # Calculate the Poisson curve.
  poisy = c()
  n = 60
  bin = 3 / n
  for (j in 1:n) {
    poisy[j] = (exp(-1 * j * bin) - exp(-1 * (j + 1) * bin))
  }
  poisy = poisy / max(poisy)    

    
  # Step 3: Fit a spline through the eigenvalues
  # Select elements evenly spaced accorindg to the pace
  oX = v[seq(1,length(v),, round(length(v)/pace))]
  oY = seq(0, 1,, round(length(v)/pace))

  has_error = 0
  tryCatch({
    sp = splinefun(oX, oY, method="hyman")
    ei = sp(v[2:(length(v)-1)])
  },
  error = function(err) {
    print(paste("ERROR:", err, sep=" "))
    has_error = 1 
  })
  if (has_error) {
    next
  } 

  ##
  # Calculate the NNSD
  #
  # Step 1: calculate the spacing levels.
  nsl = diff(ei)
  nsl = nsl * length(v)
    
  # Step 2: calculate the NNSD of the spacing levels
  # we will focus only on the region between 0 and 3.
  nsl = nsl[which(nsl >= 0 & nsl <= 3)]
  nnsd = hist(nsl, seq(0, 3,, 61), plot=FALSE)

  ##
  # Calculate the chi-square value
  chi = 0;
  n = 60
  bin = 3 / n
  exp = c();
  for (j in 1:n) {
    obs = nnsd$counts[j]
    expect = (exp(-1 * j * bin) - exp(-1 * (j + 1) * bin)) * length(v)
    exp[j] = expect
    if (expect != 0) {
      chi = chi + ((obs - expect)^2 / expect)
    }
  }   
  
  stats = paste("th:", round(k, digits=5), ", chi:", 
      round(chi, digits=2), ", pace:", pace, ", n:", length(v), sep=" ")

  # Plot the Density Distribution
  if (pace == 10 & save.plot == TRUE) {
    png(filename = paste("RMT-NNSD-pace10-", sprintf("%.5f", k), ".png", sep=""), width=960)
  }

  layout(matrix(c(1,2), 1, 2, byrow = TRUE))
  plot(seq(min(v),max(v),,5000), sp(seq(min(v), max(v),, 5000)), type="l", xlab="Ordered Eigenvalues", ylab="")
  lines(v, seq(0,1,,length(v)), col="green")
  points(oX, oY, col="blue", pch=19, lw=1)
  # points(oX2, oY2)
  # lines(q, seq(0,1,,71),  col="blue")
  # lines(q, c(0,yy,1), col="orange")
  # points(q, seq(0,1,,71), col="red")

  nnsd$counts = nnsd$counts / max(nnsd$counts)
  plot(nnsd, main=stats, xlim=c(0,3), xlab="Nearest Neighbor Spacing", ylab="Frequency")
  lines(seq(0,3,,n), poisy, col="blue", lw=2)
  lines(seq(0,3, 0.01), goe, col="red", lw=2)

  if (pace == 10 & save.plot == TRUE) {
    dev.off()
  }
  
  #print(stats)
  return(chi)
}
