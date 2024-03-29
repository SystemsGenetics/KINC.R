% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base.R
\name{mergeCommunities}
\alias{mergeCommunities}
\title{Merges linked community clusters that have the most similar sets of nodes.}
\usage{
mergeCommunities(lc, th = 0.5)
}
\arguments{
\item{lc}{A linkcomm object.}

\item{th}{Threshold to be used by `mergeClusters`. Lower threshold results in merging
being more common.}
}
\value{
A list were each element of the list is the
  set of nodes and edges of the merged clusters.
}
\description{
This merging function recursively iterates through all of the clusters
and performs a pair-wse Jaccard comparision between all clusters.
Those with the highest Jaccard score that are above the given
threshold (i.e. deafult of 0.5) are candidates for merging.
}
\details{
This is a helper function for the findLinkedCommunities() function and is
not meant to be called on its own.
}
