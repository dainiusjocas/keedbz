# This library contains functions to get feature ranking according to 
# SVM-RFE (Recursive Feature Ranking) algorithm

source('absolute_weight_svm.R', chdir=T)

GetSVMRFEFeatureRanking <- function(dataset, normal, tumor, E=0.2) {
  # Gets SVM-RFE feature ranking
  #
  # Vars:
  #   dataset: matrix where rows - features, columns - tuples
  #   normal: indexes of patients without tumor
  #   tumor: indexes of patients with tumor
  #
  # Returns:
  #   Vector with feature rankings (first position in vector contains an index 
  #   of best ranked feature).
  if (E < 0 || E > 1) {
    return(NaN)
  }
  ranking <- c()
  indexes = c(1:length(dataset[ , 1]))
  i <- length(dataset[ , 1])
  while (i > 2) {
    ranks <- GetAbsoluteWeightSVMFeatureRanking(dataset[indexes, ],
                                                normal,
                                                tumor)
    ranking <- c(ranking, indexes[tail(ranks, n=round(i * E))])
    indexes <- indexes[-tail(ranks, n=round(i * E))]
    i <- i - round(i * E)
  } 
  ranking <- c(ranking, indexes)
  return(rev(ranking))
}