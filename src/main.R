# This script is like main in the java - just an entry point
# Clean working space

rm(list = ls())

load('res/RData/colon.RData')

source('./src/basis_criteria/relief_measurement.R', chdir=T)
source('./src/basis_criteria/fisher_score.R', chdir=T)
source('./src/basis_criteria/asymmetric_dependency_coefficient.R', chdir=T)
source('./src/basis_criteria/svm/absolute_weight_svm.R', chdir=T)
source('./src/basis_criteria/svm/svm_rfe_ranking.R', chdir=T)


res1 <- GetADCFeatureRanking(nncolon, pos, neg)
res2 <- GetFisherFeatureRanking(nncolon, pos, neg)
res3 <- GetReliefFeatureRanking(nncolon, pos, neg)
res4 <- GetAbsoluteWeightSVMFeatureRanking(nncolon, pos, neg)
res5 <- GetSVMRFEFeatureRanking(nncolon, pos, neg)

print('done')