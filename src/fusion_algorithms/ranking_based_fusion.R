###############################################################################
# Library to compute Ranking-Based Multicriterion Fusion
###############################################################################

source('../basis_criteria/fisher_score.R', chdir=T)
source('../basis_criteria/asymmetric_dependency_coefficient.R', chdir=T)
source('../basis_criteria/relief_measurement.R', chdir=T)
source('../basis_criteria/svm/absolute_weight_svm.R', chdir=T)

GetRankingBasedMulticriterionFusionFeatureRanking <- function(dataset,
                                                              normal, tumor) {
  # Gets ranking-based multicriterion fusion feature ranking
  # 
  # Args:
  #   dataset: rows - features, columns - tuples
  #   normal: indexes of patiend without tumor
  #   tumor: indexes of patients with tumor
  # 
  # Results: 
  #   ranking: vector where first element is an index of best ranked feature
  #
  feature.ranking.scores <- GetFeatureRankingScores(dataset, normal, tumor)
  fusioned.ranking.scores <- 
    ComputeFusionedRankingScores(feature.ranking.scores)
  return(fusioned.ranking.scores)
}

ComputeFusionedRankingScores <- function(feature.ranking.scores) {
  # Computes the ranking-based fusion of features by formula:
  #   fusion[i] = sum(v[i], i=1..m), where m - number of features, v - vector 
  #   of ranking scores
  #
  # Args:
  #   feature.ranking.scores: matrix where columns - ranking scores, 
  #   rows - features
  #
  # Returns:
  #   fusioned.ranking.scores: vector where every element is fusioned ranking 
  #     score of a feature
  #
  fusioned.ranking.scores <- c()
  for (i in 1: length(feature.ranking.scores[ , 1])) {
    fusioned.ranking.scores[i] <- sum(feature.ranking.scores[i, ]) 
  }
  return(fusioned.ranking.scores)
}

GetFeatureRankingScores <- function(dataset, normal, tumor) {
  # Gets feature ranking from various feature ranking methods
  #
  # Args:
  #   dataset: rows - features, columns - tuples
  #   normal: patiens without tumor
  #   tumor: indexes of patients with tumor
  #
  # Results:
  #   rankings: matrix where columns are feature rankings, rows - feature rank
  #     scores. Ex. Feature with index 1, get 218 ranking points by fisher. 
  #     Column <<index>> is not returned!
  #   +-----------+----------------+----------------+
  #   | <<index>> | fisher.ranking | relief.ranking | ...
  #   +-----------+----------------+----------------+
  #   |     1     |     218        |     1200       | ...
  #   |    ...    |     ...        |     ....       | ...
  #
  feature.ranking.scores <- NULL
  # Get ranking scores out of fisher ranking
  fisher.ranking <- GetFisherFeatureRanking(dataset, normal, tumor)
  fisher.ranking.scores <- sort(rev(fisher.ranking), index.return=T)$ix
  feature.ranking.scores <- cbind(feature.ranking.scores, fisher.ranking.scores)
  # ADC ranking scores
  adc.ranking <- GetADCFeatureRanking(dataset, normal, tumor)
  adc.ranking.scores <- sort(rev(adc.ranking), index.return=T)$ix
  feature.ranking.scores <- cbind(feature.ranking.scores, adc.ranking.scores)
  # Relief ranking scores
  relief.ranking <- GetReliefFeatureRanking(dataset, normal, tumor)
  relief.ranking.scores <- sort(rev(relief.ranking), index.return=T)$ix
  feature.ranking.scores <- cbind(feature.ranking.scores, relief.ranking.scores)
  # AW-SVM ranking scores
  aw.svm.ranking <- GetAbsoluteWeightSVMFeatureRanking(dataset, normal, tumor)
  aw.svm.ranking.scores <- sort(rev(aw.svm.ranking), index.return=T)$ix
  feature.ranking.scores <- cbind(feature.ranking.scores, aw.svm.ranking.scores)
  return(feature.ranking.scores)
}
