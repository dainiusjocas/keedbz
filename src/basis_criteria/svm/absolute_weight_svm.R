###############################################################################
# Library to rank features by svm weights
###############################################################################

library(e1071)

# ENTRY POINT
# This method computes the weights of features of dataset. First, we build
#   a model of svm classifier, then we get weight values of every feature 
#   from that model. 
# NOTE: Only one iteration of building the svm model is done.
# input: dataset: rows - features, columns - tuples
# input: pos - indexes of normal patients in the dataset
# input: neg - indexes of patients with tumor in the dataset
# output: svm_weights - a vector of scores of features
GetAbsoluteWeightSVMFeatureWeights <- function(dataset, pos, neg)
{
  y <- c()
  y[pos] <- -1
  y[neg] <- 1
  # dataset should be transposed in order to compute weights
  data <- t(dataset) 
  svmModel = svm(data, y, cost = 10, cachesize=500,
                 type="C-classification", kernel="linear" )
  # For two-class situation
  svm_weights <- abs(t(svmModel$coefs)%*%svmModel$SV)
  # For multiclass situation
  # svm_weights <- abs(svm.weights(svmModel))
  return(as.vector(svm_weights))
}

# This method returns feature ranking according to svm weights
GetAbsoluteWeightSVMFeatureRanking <- function(dataset, pos, neg)
{
  ranking <- sort(GetAbsoluteWeightSVMFeatureWeights(dataset, pos, neg),
                  decreasing=T, index.return=T)$ix
  return(ranking)
}
