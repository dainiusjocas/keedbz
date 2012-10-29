###############################################################################
# Library for computing Asymetric Dependency Coefficients
###############################################################################

library(entropy)
library(multicore)
library(doMC)

# ENTRT POINT
# This method computes Asymetric Dependency Coefficients for the dimensions
#   of a given dataset. 
# DATASET SHOULD BE ALREADY NORMALIZED!!!
# input: dataset: rows - features, columns - tuples
# input: normal - indexes of normal patients in the dataset
# input: tumor - indexes of patients with tumor in the dataset
# output: list of Asymetric Dependency Coefficiets dimentions.
GetADCCoefficients <- function(dataset, normal, tumor) {
  total_number_of_items <- length(dataset[ , 1]) 
  class_probability <- c(length(normal) / total_number_of_items,
                         length(tumor)  / total_number_of_items)
  dataset_as_list <- c()
  for(i in 1:length(dataset[,1])) 
  {
    dataset_as_list[i] <- list(dataset[i, ])
  }
  # MI(Y, X)
  values_of_mutual_information <- mclapply(X=dataset_as_list,
                                           get_mutual_information,
                                           indexes_normal=normal,
                                           indexes_tumor=tumor,
                                           probability_normal=class_probability[1],
                                           probability_tumor=class_probability[2])
  result_as_vector <- c()
  for(i in 1:length(dataset[,1]))
  {
    result_as_vector[i] <- values_of_mutual_information[[i]]
  }
  # H(Y)
  entropy_of_classes <- entropy(c(length(normal), length(tumor))) 
  # ADC(Y, X) = MI(Y, X) / H(Y)
  adc <- result_as_vector / entropy_of_classes;
  return(adc)
}

# This method computes mutual information of the feature.
# input: values_of_feature - values of a specific feature
# input: indexes_normal - indexes in the list values_of_feature
#   of normal patients
# input: indexes_tumor - indexes in the list values_of_feature
#   of patients with tumor
# input: probability_normal - probability that patient is normal
# input: probability_tumor - probability thet patient has tumor
get_mutual_information <- function(values_of_feature,
                                   indexes_normal,
                                   indexes_tumor,
                                   probability_normal,
                                   probability_tumor)
{
  probability_of_normal <- density(values_of_feature[indexes_normal], kernel="epanechnikov",
                                   from=-12, to=12, n=512, adjust=0.3)$y
  probability_of_normal <- probability_of_normal / sum(probability_of_normal)
  probability_of_tumor <- density(values_of_feature[indexes_tumor], kernel="epanechnikov",
                                  from=-12, to=12, n=512, adjust=0.3)$y
  probability_of_tumor <- probability_of_tumor / sum(probability_of_tumor)
  joint_probability <- matrix(cbind(probability_of_normal * probability_normal,
                                    probability_of_tumor * probability_tumor), ncol=2)
  mutual_information <- mi.plugin(joint_probability)
  if (is.na(mutual_information)) { print(mutual_information) }
  return(mutual_information)
}

# This method returns ranking of features according to ADC features ranking
#   method. 
GetADCFeatureRanking <- function(dataset, normal, tumor) {
  ranking <- sort(GetADCCoefficients(dataset, normal, tumor), decreasing=T, index.return=T)$ix
  return(ranking)
}
