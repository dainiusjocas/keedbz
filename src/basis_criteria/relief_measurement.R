###############################################################################
# Implementation of RELIEF feature selection algorithm
###############################################################################

# input: dataset - rows -> features, columns -> tuples
# input: pos - indexes of healthy patients
# input: neg - indexes of not healthy patients
# input: m - how many times random selection of tuple from dataset will be done
# output: weights - quality scores of all features 
GetReliefWeights <- function(dataset, pos, neg, m=length(dataset[1, ]))
{
  distance_matrix <- make_distance_matrix(dataset)
  min_max_diff <- get_min_max_diff(dataset)
  weights <- c(seq(from=0, to=0, length.out=length(dataset[ , 1])))
  for(i in 1:m)
  {
    instance_index <- sample(x=c(1:length(dataset[1, ])), size=1) 
    nearest_hit <- 
      get_index_of_nearest_hit(distance_matrix, instance_index, pos, neg)
    nearest_miss <- 
      get_index_of_nearest_miss(distance_matrix, instance_index, pos, neg)
    for(j in 1:length(dataset[ , 1]))
    {
      weights[j] <- weights[j] - 
        get_diff(dataset, instance_index, j, nearest_hit, nearest_miss) /
        (m * min_max_diff[j]) 
    }
  }
  return(weights)
}

get_diff <- function(dataset, instance_index, feature_index,
                     nearest_hit, nearest_miss)
{
  difference <- 
    abs(dataset[feature_index, instance_index] -
    dataset[feature_index, nearest_hit]) -
    abs(dataset[feature_index, instance_index] -
    dataset[feature_index, nearest_miss])
}

# This method computes differences between max and min values of every feature
# input: dataset - rows -> features, columns -> tuples 
get_min_max_diff <- function(dataset)
{
  min_max_diff <- c()
  for(i in 1:length(dataset[ , 1]))
  {
    min_max_diff[i] <- max(dataset[i, ]) - min(dataset[i, ])
  }
  return(min_max_diff)
}

# This method gets an index of a nearest neighbour of current instance from 
#   the same class
# input: distance_matrix
# input: instance_index
# input: pos
# input: neg
# output: index_of_nearest_hit
get_index_of_nearest_hit <- function(distance_matrix, instance_index, pos, neg)
{
  if (any(instance_index == pos))
  {
    nearest_hit <- 
      which(distance_matrix[instance_index, pos] 
            == min(distance_matrix[instance_index, pos]))
    nearest_hit <- pos[nearest_hit]
  }
  else
  {
    nearest_hit <- 
      which(distance_matrix[instance_index, neg]
            == min(distance_matrix[instance_index, neg]))
    nearest_hit <- neg[nearest_hit]
  }
  return(head(nearest_hit))
}

# This method gets an index of a nearest neighbour of current instance from
#   other class.
# input: distance_matrix
# input: instance_index
# input: pos
# input: neg
# output: index_of_nearest_miss
get_index_of_nearest_miss <- function(distance_matrix, instance_index, pos, neg)
{
  if (any(instance_index == pos))
  {
    nearest_miss <- 
      which(distance_matrix[instance_index, neg] 
            == min(distance_matrix[instance_index, neg]))
    nearest_miss <- neg[nearest_miss]
  }
  else
  {
    nearest_miss <- 
      which(distance_matrix[instance_index, pos]
            == min(distance_matrix[instance_index, pos]))
    nearest_miss <- pos[nearest_miss]
  }
  return(head(nearest_miss))
}

# This method makes distance matrix which if full matrix and this matrix
#   is specific because in the diagonal there are Inf's (normally diag is 0)
# input: dataset - rows -> features, columns -> tuples
# output: distance_matrix - in the diagonal there is Inf
make_distance_matrix <- function(dataset)
{
  distance_matrix <- as.matrix(dist(t(dataset), 
                                    method = "euclidean",
                                    upper = TRUE))
  diag(distance_matrix) <-Inf
  distance_matrix <- matrix(distance_matrix, ncol=length(dataset[1, ]))
  return(distance_matrix)
}

GetReliefFeatureRanking <- function(dataset, normal, tumor, m=length(dataset[1, ])) 
{
  # Gets ranking of features according to relief weights based feature ranking
  #   method.
  #
  # Args:
  #   dataset: matrix where rows - features, columns - tuples
  #   normal: indexes of patients without tumor
  #   tumor: indexes of patients with tumor
  #   m: number of iterations. By default it is equal to number of tuples in the 
  #      dataset.
  #
  # Returns:
  #   Vector with feature rankings (first position in vector contains an index 
  #   of best ranked feature).
  ranking <- sort(GetReliefWeights(dataset, normal, tumor, m), decreasing=T,
                  index.return=T)$ix  
  return(ranking)
}
