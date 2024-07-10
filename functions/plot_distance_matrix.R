plot_distance_matrix <- function(distance_matrix) {
  # Convert the distance matrix to a data frame suitable for ggplot
  distance_df <- as.data.frame(as.table(as.matrix(distance_matrix)))
  
  # Plot using ggplot
  p <- ggplot(data = distance_df, aes(Var1, Var2, fill = Freq)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "white", 
      mid = "gray", 
      high = "purple4",
      na.value = 'red',
      midpoint = median(distance_df$Freq,na.rm=T))
    labs(x = 'Row', y = 'Column', fill = 'Distance') +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Print the plot
  print(p)
}
