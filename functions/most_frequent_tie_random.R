most_frequent_tie_random <- function(matrix) {
  get_most_frequent <- function(row) {
    tab <- table(row)
    max_count <- max(tab)
    most_frequent_values <- names(tab[tab == max_count])
    selected_value <- sample(most_frequent_values, 1)
    return(c(as.numeric(selected_value), max_count))
  }
  
  result <- t(apply(matrix, 1, get_most_frequent))
  df <- data.frame(value = result[, 1], frequency = result[, 2])
  return(df)
}

## Example usage
#matrix <- matrix(c(1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5), nrow = 4, byrow = TRUE)
#print(matrix)
#print(most_frequent(matrix))
