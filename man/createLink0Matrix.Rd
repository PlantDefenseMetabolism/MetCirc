% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linkMatrix.R
\name{createLink0Matrix}
\alias{createLink0Matrix}
\title{Create a link matrix}
\usage{
createLink0Matrix(similarityMatrix)
}
\arguments{
\item{similarityMatrix}{\code{matrix}, a similarity matrix that contains the 
NDP similarity measure between all precursors in the data set}
}
\value{
createLink0Matrix returns a \code{matrix} that gives per each row 
information on linked features
}
\description{
Create a link matrix which links every feature in similarity
matrix with another.
}
\details{
createLink0Matrix creates a \code{matrix} from a similarity 
matrix which includes all connections between features in the 
similarity matrix, but 
exclude links which have a similarity of exactly 0.
}
\examples{
data("binnedMSP", package = "MetCirc")
## truncate binnedMSP
binnedMSP <- binnedMSP[1:28,]
namesPrec <- rownames(binnedMSP)
similarityMat <- createSimilarityMatrix(binnedMSP)
link0Mat <- createLink0Matrix(similarityMatrix = similarityMat)
}
\author{
Thomas Naake, \email{thomasnaake@googlemail.com}
}
