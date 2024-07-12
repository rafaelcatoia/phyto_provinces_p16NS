
# distance_matrix <- matrix(c(
#   0, 7.76, 5.18, 5.27, 1.11, 1.39, 3.02, 4.77,
#   7.76, 0, 2.73, 8.78, 8.91, 3.54, 1.83, 1.90,
#   5.18, 2.73, 0, 2.45, 7.89, 0.87, 6.31, 2.06,
#   5.27, 8.78, 2.45, 0, 7.00, 4.78, 6.91, 2.48,
#   1.11, 8.91, 7.89, 7.00, 0, 1.78, 8.12, 4.14,
#   1.39, 3.54, 0.87, 4.78, 1.78, 0, 6.03, 3.96,
#   3.02, 1.83, 6.31, 6.91, 8.12, 6.03, 0, 8.75,
#   4.77, 1.90, 2.06, 2.48, 4.14, 3.96, 8.75, 0
# ), nrow = 8, ncol = 8, byrow = TRUE)
# 
# group_vector=c(1,1,2,2,2,3,3,4)

# Function to compute average mean distance to the medoid within groups and average distance between medoids of each group
compute_medoid_distances <- function(distance_matrix, group_vector) {
  unique_groups <- unique(group_vector)
  num_groups <- length(unique_groups)
  
  # Initialize storage for medoids and within-group distances
  medoids <- vector("list", num_groups)
  within_group_distances <- numeric(num_groups)
  
  for (i in 1:num_groups) {
    group <- unique_groups[i]
    group_indices <- which(group_vector == group)
    
    if(length(group_indices)>1){
    
    # Extract the submatrix for the group
    group_distance_matrix <- distance_matrix[group_indices, group_indices]
    
    # Compute the medoid
    pam_result <- cluster::pam(as.dist(group_distance_matrix), k = 1)
    medoid_index <- group_indices[pam_result$id.med]
    medoids[[i]] <- medoid_index
    
    # Compute average mean distance to the medoid within the group
    within_group_distances[i] <- mean(group_distance_matrix[upper.tri(group_distance_matrix)])
    }else{
      medoids[[i]] <- group_indices
      within_group_distances[i] <- 0
    }
  }
  
  # Compute the average distance between medoids of each group
  medoid_distance_matrix <- distance_matrix[unlist(medoids), unlist(medoids)]
  between_medoids_distances <- medoid_distance_matrix[upper.tri(medoid_distance_matrix)]
  average_between_medoids_distance <- mean(between_medoids_distances)
  
  result <- data.frame(
    average_within_group_distances = mean(within_group_distances),
    average_between_medoids_distance = average_between_medoids_distance,
    ratio_within_between = mean(within_group_distances)/average_between_medoids_distance
  )
  
  return(result)
}
