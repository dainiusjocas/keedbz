###############################################################################
# Library to compute Score-Based Multicriterion Fusion
###############################################################################

source('../basis_criteria/fisher_score.R', chdir=T)
source('../basis_criteria/asymmetric_dependency_coefficient.R', chdir=T)
source('../basis_criteria/relief_measurement.R', chdir=T)
source('../basis_criteria/svm/absolute_weight_svm.R', chdir=T)

GetScoreBasedMulticriterionFusionFeatureRanking <- function(dataset, pos, neg) {
  # Gets score-based multicriterion fusion
  #
  # Args:
  #   dataset: matrix where rows - features, columns - tuples
  #   normal: indexes of patients without tumor
  #   tumor: indexes of patients with tumor
  #
  # Results:
  #   feature.ranking.according.to.fusioned.scores: vector where first element
  #     is an index of best feature and so on
  #
  list.of.scores <- GetListOfScores(dataset, pos , neg)
  fusioned.scores <- GetFusionedScores(list.of.scores)
  feature.ranking.according.to.fusioned.scores <- 
    sort(fusioned.scores, decreasing=T, index.return=T)$ix
  return(feature.ranking.according.to.fusioned.scores)
}

GetFusionedScores <- function(list.of.scores) {
  # Gets fusioned weights of a feature's scores.
  #
  # Args:
  #   list.of.scores: matrix where in every column are weights for all features
  #     gained by one feature weighting method.
  #
  # Results:
  #   fusion.of.scores: vector where every element is a fusioned weight of a
  #   feature
  #
  normalized.list.of.scores <- NormalizeListOfScores(list.of.scores)
  fusioned.scores <- ComputeFusionOfScores(normalized.list.of.scores)
  return(fusioned.scores)
}

NormalizeListOfScores <- function(list.of.scores) {
  # Normalizes the list of scores by formula:
  # U[i]' = (U[i] - U[i[min]]) / (U[i[max]] - U[i[min]]), where
  # U[i] is scores of one ranking method
  #
  # Args: 
  #  list.of.scores: rows - scores of feature, columns - scores by ranking
  #    methods
  # 
  # Results:
  #   normalized.list.of.scores - all values falls into [0..1]
  #
  normalized.list.of.scores <- list.of.scores
  for (i in 1:length(list.of.scores[1, ])) {
    new_score <- (list.of.scores[, i]) - min(list.of.scores[ , i])
    normalized.list.of.scores[ , i] <- 
      new_score / 
      (max(list.of.scores[ , i]) - 
      min(list.of.scores[ , i]))
  }
  return(normalized.list.of.scores)
}

ComputeFusionOfScores <- function(normalized.list.of.scores) {
  # This method computes fusion of scores by formula:
  # U = (1/m) * sum(U[i]', i=1..m), where m - number of methods used,
  #   U[i]' - normalized values of scores computed by one method
  #
  # Args:
  #   normalized.list.of.scores
  #   
  # Results: 
  #   fusion.of.scores - vector of values of fusioned feature scores
  #
  fusion.of.scores <- c()
  number.of.methods <- length(normalized.list.of.scores[1, ])
  for (i in 1:length(normalized.list.of.scores[, 1])) {
    fusion.of.scores[i] <- 
      sum(normalized.list.of.scores[i, ]) / 
      number.of.methods
  }
  return(fusion.of.scores)
}

GetListOfScores <- function(dataset, pos, neg) {
  # Constructs list of lists of feature scores according to multiple feature 
  #   scoring methods. One list (column) one scoring method.
  #
  # Args:
  #   dataset: matrix where rows - features, columns - tuples
  #   normal: indexes of patients without tumor
  #   tumor: indexes of patients with tumor
  #
  # Results:
  #   list.of.scores: matrix where columns weights, rows - features
  #
  fisher.scores <- GetFisherScores(dataset, pos, neg)
  relief.weights <- GetReliefWeights(dataset, pos, neg)
  adc.coefficients <- GetADCCoefficients(dataset, pos, neg)
  svm.absolute.weights <- GetAbsoluteWeightSVMFeatureWeights(dataset, pos, neg)
  list.of.scores <- NULL
  list.of.scores <- cbind(list.of.scores, fisher.scores)
  list.of.scores <- cbind(list.of.scores, relief.weights)
  list.of.scores <- cbind(list.of.scores, adc.coefficients)
  list.of.scores <- cbind(list.of.scores, svm.absolute.weights)
  return(list.of.scores)
}
