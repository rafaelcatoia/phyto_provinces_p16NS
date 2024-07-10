### -------------------
### Abiotic Variability
### -------------------
library(dplyr) ; library(ggplot2) ; library(gridExtra)# ; library(grid)
#f_GA = df_geo_abiotics
#bioVariable = 'Temperature'
#
abioVariation <- function(df_GA,abioVariable){
  
  ######### -----------------------------------------------------------------
  df_GA$AbioFact = df_GA %>% select(any_of(abioVariable))
  df_GA = df_GA %>% select(Latitude,Depth,AbioFact)
  df_GA = df_GA %>% mutate(DistToEquador = distance_to_equator(Latitude)/100) %>% 
    mutate(DistToEquadorSplit2=cut(DistToEquador,breaks = seq(0,80,5),include.lowest = T)) %>% 
    mutate(DepthSplit2=cut(Depth,breaks = c(0,25,50,100,150,250,400,600),include.lowest = T)) %>% 
    mutate(DistToEquadorSplit = categorize_by_quantiles(DistToEquador)) %>%
    mutate(DepthSplit = categorize_by_quantiles(Depth))
    #mutate(abioQuantile = categorize_by_quantiles(AbioFact))
    
  abioDistEqPLot = df_GA %>% ggplot(aes(x=DistToEquador,y=AbioFact,color=Depth))+
    geom_point(size=3)+
    theme_minimal(base_size = 16)+
    theme(legend.position = 'bottom')
  
  abioDepthPLot = df_GA %>% ggplot(aes(x=Depth,y=AbioFact,color=DistToEquador))+
    geom_point(size=3)+
    theme_minimal(base_size = 16)+
    theme(legend.position = 'bottom')
  
  
  plt_points <- ggpubr::ggarrange(abioDistEqPLot,abioDepthPLot,ncol=2)
  plt_points <- ggpubr::annotate_figure(plt_points,top=textGrob(abioVariable))
  ######### -----------------------------------------------------------------
  df_qts = df_GA %>% 
    group_by(DistToEquadorSplit2) %>% 
    mutate(Eq_Qt0.1 = quantile(AbioFact,0.1),
           Eq_Qt0.9 = quantile(AbioFact,0.9)) %>% 
    ungroup() %>% 
    group_by(DepthSplit2) %>% 
    mutate(Dp_Qt0.1 = quantile(AbioFact,0.1),
           Dp_Qt0.9 = quantile(AbioFact,0.9)) %>% 
    ungroup()
    
  
  plt_quant_DistEq = df_qts  %>% select(DistToEquadorSplit2,Eq_Qt0.1,Eq_Qt0.9) %>%
    distinct() %>%
    ggplot()+
    geom_point(aes(x=DistToEquadorSplit2,y=Eq_Qt0.1,color='Quantile_10'),size=3)+
    geom_point(aes(x=DistToEquadorSplit2,y=Eq_Qt0.9,color='Quantile_90'),size=3)+
    scale_fill_brewer(palette="Spectral")+
    theme_minimal(base_size = 16)+
    theme(legend.position = 'bottom')+
    ylab('Quantiles')+
    scale_color_manual(values = c("Quantile_10" = 'deepskyblue3',"Quantile_90" = 'firebrick'))
  
  plt_quant_Depth = df_qts %>% select(DepthSplit2,Dp_Qt0.1,Dp_Qt0.9) %>%
    distinct() %>% ggplot()+
    geom_point(aes(x=DepthSplit2,y=Dp_Qt0.1,color='Quantile_10'),size=3)+
    geom_point(aes(x=DepthSplit2,y=Dp_Qt0.9,color='Quantile_90'),size=3)+
    scale_fill_brewer(palette="Spectral")+
    theme_minimal(base_size = 16)+
    theme(legend.position = 'bottom')+
    ylab('Quantiles')+ 
    scale_color_manual(values = c("Quantile_10" = 'deepskyblue3',"Quantile_90" = 'firebrick'))
  
  plt_quantiles <- ggpubr::ggarrange(plt_quant_DistEq,plt_quant_Depth,ncol=2)
  plt_quantiles <- ggpubr::annotate_figure(plt_quantiles,top=textGrob(abioVariable) )
  
  return(list(plt_points = plt_points,plt_quantiles=plt_quantiles))
}


distance_to_equator <- function(latitude_degrees){
  # Radius of the Earth in kilometers
  earth_radius <- 6371
  # Convert latitude from degrees to radians
  latitude_radians <- latitude_degrees * pi / 180
  # Latitude of the Equator in radians
  equator_latitude <- 0
  # Haversine formula
  distance <- earth_radius * acos(cos(latitude_radians) * cos(equator_latitude) +
                                    sin(latitude_radians) * sin(equator_latitude))
  return(distance)
}

categorize_by_quantiles <- function(numeric_vector, quantile_vector=seq(0.1,0.9,0.1)) {
  length_qtVector = length(quantile_vector)
  if(length(numeric_vector) == 0 || length(quantile_vector) == 0) {
    stop("Input vectors cannot be empty.")
  }
  if(any(quantile_vector < 0) || any(quantile_vector > 1)) {
    stop("Quantiles must be between 0 and 1.")
  }
  
  # Compute quantiles
  quantiles <- quantile(numeric_vector, quantile_vector)
  quantiles_names = paste0('qt',1:(length_qtVector+1)-1)
  # Categorize values based on quantiles
  categorized <- cut(numeric_vector,
                     breaks = c(-Inf, quantiles, Inf),
                     labels = quantiles_names,
                     include.lowest = TRUE)
  
  return(categorized)
}



#categorize_by_quantiles(numeric_vector = 1:100,quantile_vector = seq(0.1,0.9,0.1))
