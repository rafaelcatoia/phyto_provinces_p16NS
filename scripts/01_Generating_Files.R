########################################################################################################################################
## 01 - Generating Files:
## - data frame containing lat depth abiotics and long hurst provinces for each sample
## - list of distance matrices (Biotic, Abiotic, and Geographical)
## - list of distance matrices normalized (Biotic, Abiotic, and Geographical)
## - list containing perturbed versions of the Biotic component - to find k (number of clusters) and alpha (geographical structure)
########################################################################################################################################

### packages  ----------------------------------------------------------------------
library(dplyr) ; library(ggplot2)

## functions  ---------------------------------------------------------------------- 
### to filter ASVs
`%not_in%` <- purrr::negate(`%in%`)

### normalizing matrices
normalizeMatrix <- function(XX){
  normMat = norm(XX,type='2')
  return(XX/normMat)
}


### Getting data & path
root <- rprojroot::has_file(".git/index")
datadir = root$find_file("data")
savingdir = root$find_file("saved_files")
path_df = root$find_file('data/grump_phyto.csv')

## filtering only the cruises that we want and depth
grump_longer <- data.table::fread(path_df) %>%
  filter(Cruise %in% c('P16N','P16S')) %>% 
  filter(Depth<600)

# creating a key for the ASVs -- only valid for this data set. 
grump_longer = grump_longer  %>% 
  left_join(grump_longer  %>% 
              select(ASV_hash) %>% distinct() %>% 
              mutate(ID_ASV_Num = 1:n()) %>% 
              mutate(ID_ASV = ifelse(ID_ASV_Num<10,paste('000',ID_ASV_Num,sep=''),
                                     ifelse(ID_ASV_Num<100,paste('00',ID_ASV_Num,sep=''),
                                            ifelse(ID_ASV_Num<1000,paste('0',ID_ASV_Num,sep=''),
                                                   paste(ID_ASV_Num))))) %>% 
              mutate(ID_ASV = paste0('ASV_',ID_ASV)) %>% 
              select(ASV_hash,ID_ASV)) %>% 
  filter(!is.na(ID_ASV)) %>% 
  mutate(Raw.Sequence.Counts = Corrected_sequence_counts)

### -----------------------------------------------------------
## Creating abiotics dataframe  -------------------------------
### -----------------------------------------------------------
vet_abiotic = c("Temperature","Salinity","Oxygen",
  "Silicate","NO2","NO3","PO4"
)

## --- we must correct the missing on the abiotics !
df_geo_abiotics <- grump_longer %>%
  select(SampleID,one_of(vet_abiotic),Latitude,Longitude,Depth,Longhurst_Short) %>%
  distinct() %>% 
  arrange(SampleID)

## Here we have the sample with missing values on Oxygen and NO3
df_geo_abiotics %>% 
  filter(NO3<0 |Oxygen<0)


## First inputation = 212.9
df_geo_abiotics %>% filter(Latitude ==18) %>% arrange(Depth)
df_geo_abiotics$Oxygen[df_geo_abiotics$SampleID=='P16N-S40-N10'] <- 212.9

## Second inputation = 0.13 + 2.49 / 2 = 1.375
df_geo_abiotics %>% filter(Latitude ==-18.0) %>% arrange(Depth)
df_geo_abiotics$NO3[df_geo_abiotics$SampleID=='P16S-S05-N15'] <- 1.375

## Third inputation = 346.9
df_geo_abiotics %>% filter(Latitude ==-63.5) %>% arrange(Depth)
df_geo_abiotics$Oxygen[df_geo_abiotics$SampleID=='P16S-S96-N28'] <- 346.9

## looking at the samples
# df_geo_abiotics %>% ggplot(aes(x=Latitude,y=Depth))+ geom_point(size=2)+scale_y_reverse()

## This is how it looks after the imputation 
df_geo_abiotics %>% 
  tidyr::pivot_longer(cols = any_of(vet_abiotic),names_to = 'AbioticFactor') %>% 
  group_split(AbioticFactor) %>% 
  purrr::map(
    ~ggplot(., aes(x=Latitude,y=Depth, color = value)) + 
      geom_point(size = 3) +
      theme_minimal()+
      scale_y_reverse() +
      scale_colour_gradient2(
        low = "black", 
        mid = "gray", 
        high = "deepskyblue3",
        na.value = 'red',
        midpoint = median(.$value,na.rm=T))+
      ggtitle(unique(.$AbioticFactor))
    #facet_wrap(~ AbioticFactor, labeller = function(x) label_value(x, multi_line = FALSE))
  ) %>% 
  cowplot::plot_grid(plotlist = ., align = 'hv', ncol = 2)

## Saving the abiotics df
saveRDS(df_geo_abiotics,file = paste0(savingdir,'/','df_geo_abiotics'))

## Saving the filtered Grump stacked
saveRDS(grump_longer,file = paste0(savingdir,'/','dfGrump_longer_phyto'))

### --------------------------------------------------------------
### Creating the list containing distance matrices 
### Normalized an unnormalized
### --------------------------------------------------------------

###### First the normalized ones  ################################
geo_Dist_mirrored = df_geo_abiotics %>% 
  mutate(Latitude = abs(Latitude)) %>% 
  transmute(lat_scaled = (Latitude - mean(Latitude)) /sd(Latitude),
            depht_scaled = (Depth -mean(Depth)) /sd(Depth)) %>%
  as.matrix() %>% dist() %>% as.matrix() #%>% normalizeMatrix()

abioticDist = df_geo_abiotics %>% 
  transmute(Temperature = scale(Temperature),
            Salinity = scale(Salinity),
            Oxygen = scale(Oxygen),
            Silicate = scale(Silicate),
            NO2 = scale(NO2),
            NO2 = scale(NO3),
            PO4 = scale(PO4)) %>% 
  as.matrix() %>% dist() %>% as.matrix() #%>%  normalizeMatrix()


biotic_dist = grump_longer %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  select(SampleID,ID_ASV,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,ID_ASV) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = ID_ASV ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() #%>% normalizeMatrix()

list_abio_bio_geo_dist <- list(
  geoDist = geo_Dist_mirrored,
  abioticDist = abioticDist,
  bioticDist = biotic_dist
  )

list_abio_bio_geo_dist_normallized <- list(
  geoDist = geo_Dist_mirrored %>% normalizeMatrix(),
  abioticDist = abioticDist %>% normalizeMatrix(),
  bioticDist = biotic_dist %>% normalizeMatrix()
)

saveRDS(list_abio_bio_geo_dist,file = paste0(savingdir,'/','list_abio_bio_geo_dist'))
saveRDS(list_abio_bio_geo_dist_normallized,file = paste0(savingdir,'/','list_abio_bio_geo_dist_normallized'))


### ----------------------------------------------------------------------
### Creating the B replicates of aitchison distance -- To tune K and alpha
### ----------------------------------------------------------------------

## number of replicates 
B=500

## list with perturbed versions of the original data
list_AitDist = list()

## Creating linst of all ASVs
idASVs = grump_longer %>% select(ID_ASV) %>% distinct() %>% pull

## percentage of columns
pct_colSubSample = 0.5 

## minimal value to add in order to get the CLR transform / Aitchison distance
min_raw_count = grump_longer %>% select(Raw.Sequence.Counts) %>% min()
min_raw_count = min_raw_count/100

## initiating the counter we will use a while loop, because when removing some ASVs
## it is possible that some samples are going to be excluded, and we don't
## want that to happen.

accepted_sample <- 1 
seedI = 5757

## untill we have less then or equal to B replicates:
while(accepted_sample<=B){
  ## sampling the ASVs to be used
  asvSubset=sample(idASVs,size = pct_colSubSample*length(idASVs))
  ## creating the aitchison distance matrix 
  mat_repB <- grump_longer %>%
    ## 
    filter(ID_ASV %in% asvSubset) %>% 
    mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
    select(SampleID,ID_ASV,Raw.Sequence.Counts) %>% distinct() %>% 
    group_by(SampleID,ID_ASV) %>% 
    summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
    tidyr::pivot_wider(id_cols = SampleID,names_from = ID_ASV ,
                       values_from = Sum_RawCounts,
                       values_fill = min_raw_count) %>%
    data.frame() %>% 
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
    relocate(SampleID,sort(names(.))) %>% 
    arrange(SampleID) %>% 
    select(-SampleID) %>% 
    vegan::vegdist(method = 'aitchison') %>% as.matrix() %>% 
    normalizeMatrix()
  
  mat_evalB <- grump_longer %>%
    filter(ID_ASV %not_in% asvSubset) %>% 
    mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
    select(SampleID,ID_ASV,Raw.Sequence.Counts) %>% distinct() %>% 
    group_by(SampleID,ID_ASV) %>% 
    summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
    tidyr::pivot_wider(id_cols = SampleID,names_from = ID_ASV ,
                       values_from = Sum_RawCounts,
                       values_fill = min_raw_count) %>%
    data.frame() %>% 
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
    relocate(SampleID,sort(names(.))) %>% 
    arrange(SampleID) %>% 
    select(-SampleID) %>% 
    vegan::vegdist(method = 'aitchison') %>% as.matrix()
  
  if( sum(c(dim(mat_repB),dim(mat_evalB))!=nrow(df_geo_abiotics))  == 0){
    list_AitDist = rlist::list.append(list_AitDist,list(INS=mat_repB,OOS=mat_evalB))
    accepted_sample = accepted_sample+1
  }
  cat(paste0('Iteration number ------ ::',accepted_sample,'\n'))
}

saveRDS(list_AitDist,file = paste0(savingdir,'/','list_AitDist_IS_OOS'))

### -------------------------------------------------------------------------
### Creating the B replicates of Aitchison distance -- To find Provinces
### -------------------------------------------------------------------------
## number of replicates 
B=100

## list with perturbed versions of the original data
list_AitDist = list()

## Creating linst of all ASVs
idASVs = grump_longer %>% select(ID_ASV) %>% distinct() %>% pull

## percentage of columns
pct_colSubSample = 0.9

## minimal value to add in order to get the CLR transform / Aitchison distance
min_raw_count = grump_longer %>% select(Raw.Sequence.Counts) %>% min()
min_raw_count = min_raw_count/100

## initiating the counter we will use a while loop, because when removing some ASVs
## it is possible that some samples are going to be excluded, and we don't
## want that to happen.

accepted_sample <- 1 
seedI = 5757

while(accepted_sample<=B){
  ## sampling the ASVs to be used
  asvSubset=sample(idASVs,size = pct_colSubSample*length(idASVs))
  ## creating the aitchison distance matrix 
  mat_repB <- grump_longer %>%
    ## 
    filter(ID_ASV %in% asvSubset) %>% 
    mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
    select(SampleID,ID_ASV,Raw.Sequence.Counts) %>% distinct() %>% 
    group_by(SampleID,ID_ASV) %>% 
    summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
    tidyr::pivot_wider(id_cols = SampleID,names_from = ID_ASV ,
                       values_from = Sum_RawCounts,
                       values_fill = min_raw_count) %>%
    data.frame() %>% 
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
    relocate(SampleID,sort(names(.))) %>% 
    arrange(SampleID) %>% 
    select(-SampleID) %>% 
    vegan::vegdist(method = 'aitchison') %>% as.matrix() %>% 
    normalizeMatrix()
  
  if( sum(dim(mat_repB) != nrow(df_geo_abiotics))  == 0){
    list_AitDist = rlist::list.append(list_AitDist,mat_repB)
    accepted_sample = accepted_sample+1
  }
  cat(paste0('Iteration number ------ ::',accepted_sample,'\n'))
}

saveRDS(list_AitDist,file = paste0(savingdir,'/','list_AitDist_IS_provinces'))

### -------------------------------------------------------------------------
### Now, moove to script 02! 
### -------------------------------------------------------------------------