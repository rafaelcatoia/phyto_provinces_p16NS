########################################################
## Script to get a list containing k "Bootstrap" samples
########################################################
#
getwd()
grump_stacked <- data.table::fread('/Users/rafaelcatoia/Desktop/repos/clusterPhyto_16NS/phyto_provinces_p16NS/data/ex_grump_asv_long20240501.csv')

## Yubins code
df <- read.csv("/file/path/to/that/dropbox/link/file")
## Subset out the picocyanobacteria and diazotrophs from the main dataframe
df[, c(2:11,54)][is.na(df[, c(2:11,54)])] <- "Not Applicable" 
cyano <- rbind(df[df$Eco_relevant_plank_groups=="Prochlorococcus",],
               df[df$Eco_relevant_plank_groups=="Synechococcus",],
               df[grepl("Richelia", df$Genus),],
               df[grepl("Trichodesmium", df$Genus),],
               df[grepl("UCYN-A", df$Genus),],
               df[grepl("Crocosphaera", df$Genus),],
               df[grepl("UCYN-C", df$Genus),])
## Subset out the chloroplast 16S ASVs from the main dataframe
chl.16s <- df[df$Sequence_Type=="Chloroplast_16S",]
## Remove the following ASVs from that subset
chl.16s <- chl.16s[!chl.16s$Supergroup=="Rhizaria",]
chl.16s <- chl.16s[!chl.16s$Supergroup=="Excavata",]
chl.16s <- chl.16s[!chl.16s$Supergroup=="Alveolata",]
chl.16s <- chl.16s[!chl.16s$Division=="Rhodophyta",]
## Final "phytoplankton" community we will use for the time being
phyto.all <- rbind(cyano, chl.16s)
## From here, subset out the P16 S/N cruises to focus on that <- this is more or less the 
## transect that I have been using for a lot of other analysis



grump_stacked %>% colnames()

grump_stacked %>%  filter(Level_2 %in% c(“Cyanobacteria”, “Chloroplast_16S”, “Phytoplankton_18S”))

library(dplyr)

root <- rprojroot::has_file(".git/index")
datadir = root$find_file("data")
funsdir = root$find_file("functions")
savingdir = root$find_file("saved_files")

files_vec <- list.files(funsdir)

for( i in 1:length(files_vec)){
  source(root$find_file(paste0(funsdir,'/',files_vec[i])))
}

dat_tax = data.table::fread('https://raw.githubusercontent.com/rafaelcatoia/zoop_16N/main/treated_taxonomy_dat.csv') %>%
  as_tibble()

### Loading the new data 

## Use this if you are using the new GRUMP Data Set
datapath = root$find_file(paste0(datadir,'/','grump_asv_long.csv'))
dframe = data.table::fread(input = datapath) %>%
  filter(Cruise %in% c('P16N','P16S')) %>% 
  #filter(Raw.Sequence.Counts>0) %>% 
  #filter(Domain!='Unassigned') %>% 
  mutate(Raw.Sequence.Counts = Corrected_sequence_counts) %>% 
  #### New filter by Depth!
  filter(Depth<600)

#### -- now lets subset by only using the ASV that we had before
dframe <- dframe %>% filter(ASV_hash %in% dat_tax$ASV_ID) %>% 
  left_join(dframe %>% filter(ASV_hash %in% dat_tax$ASV_ID) %>% 
              select(ASV_hash) %>% distinct() %>% 
              mutate(ID_ASV_Num = 1:n()) %>% 
              mutate(ID_ASV = ifelse(ID_ASV_Num<10,paste('000',ID_ASV_Num,sep=''),
                                     ifelse(ID_ASV_Num<100,paste('00',ID_ASV_Num,sep=''),
                                            ifelse(ID_ASV_Num<1000,paste('0',ID_ASV_Num,sep=''),
                                                   paste(ID_ASV_Num))))) %>% 
              mutate(ID_ASV = paste0('ASV_',ID_ASV)) %>% 
              select(ASV_hash,ID_ASV)) %>% 
  filter(!is.na(ID_ASV)) 

### -----------------------------------------------------------
## Creating abiotics dataframe  -------------------------------
### -----------------------------------------------------------
vet_abiotic = c(
  "Temperature",
  "Salinity",
  "Oxygen",
  "Silicate",
  "NO2",
  "NO3",#this causes duplicates
  #"NH3",#this is empty
  "PO4"
)

## Filtering the abiotics per sample
df_geo_abiotics <- dframe %>%
  select(SampleID,one_of(vet_abiotic),Latitude,Longitude,Depth,Longhurst_Short) %>%
  mutate(Oxygen=ifelse(Oxygen<0,NA,Oxygen),
         NO3=ifelse(NO3<0,NA,NO3)) %>% 
  distinct() %>% arrange(SampleID)

## looking at the samples
library(ggplot2)
df_geo_abiotics %>% ggplot(aes(x=Latitude,y=Depth))+
  geom_point(size=5)+
  scale_y_reverse()

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
        midpoint = median(.$value,na.rm=T)
      )+
      ggtitle(unique(.$AbioticFactor))
      #facet_wrap(~ AbioticFactor, labeller = function(x) label_value(x, multi_line = FALSE))
  ) %>% 
  cowplot::plot_grid(plotlist = ., align = 'hv', ncol = 2)



## Solving the missing problem -------------------------------------------------

df_geo_abiotics <- dframe %>%
  select(SampleID,one_of(vet_abiotic),Latitude,Longitude,Depth,Longhurst_Short) %>%
  distinct() %>% 
  #mutate(Oxygen=ifelse(Oxygen<0,NA,Oxygen),
  #       NO3=ifelse(NO3<0,NA,NO3)) %>% 
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

## This is how it looks after the inputation 
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


### Inputing those values:: 

## Saving the abiotics df
saveRDS(df_geo_abiotics,file = paste0(savingdir,'/','df_geo_abiotics'))
dim(df_geo_abiotics)

## Saving the filtered Grump stacked
saveRDS(dframe,file = paste0(savingdir,'/','dfGrump_longer_filtered'))

### --------------------------------------------------------------
### Creating the B replicates of aitchison distances 
### --------------------------------------------------------------
B=250
list_AitDist = list()
idASVs = dframe %>% select(ID_ASV) %>% distinct() %>% pull
pct_colSubSample = 0.75 

min_raw_count = dframe %>% select(Raw.Sequence.Counts) %>% min()
min_raw_count = min_raw_count/1000

seedI = 5757
set.seed(seedI)
for(ii in 1:B){
  
  set.seed(seedI+ii)
  asvSubset=sample(idASVs,size = pct_colSubSample*length(idASVs))
  
  list_AitDist[[ii]] <- dframe %>%
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
}

### Sometimes, when we colsample the ASVs, we basically loose some samples. So we are only keeping the subset that 
length(list_AitDist)
list_AitDist = Filter(function(x) dim(x)[1] == 174, list_AitDist)
length(list_AitDist)
####### here we have a list of "bootstraped" normalized aitDist
####### now lets save this and move to the next step
saveRDS(list_AitDist,file = paste0(savingdir,'/','list_AitDist'))


#####################################################################################################

###### creating the distance matrices that we will use in the convex mixture in the future
df_geo_abiotics = readRDS(file = paste0(savingdir,'/','df_geo_abiotics'))
df_geo_abiotics = df_geo_abiotics %>% arrange(SampleID)

## First the normalized ones 
geo_Dist = df_geo_abiotics %>% 
  transmute(lat_scaled = (Latitude -mean(Latitude)) /sd(Latitude),
            depht_scaled = (Depth -mean(Depth)) /sd(Depth)) %>%
  as.matrix() %>% dist() %>% as.matrix() %>% normalizeMatrix()

geo_Dist_mirrored = df_geo_abiotics %>% 
  mutate(Latitude = abs(Latitude)) %>% 
  transmute(lat_scaled = (Latitude - mean(Latitude)) /sd(Latitude),
            depht_scaled = (Depth -mean(Depth)) /sd(Depth)) %>%
  as.matrix() %>% dist() %>% as.matrix() %>% normalizeMatrix()

abioticDist = df_geo_abiotics %>% 
  transmute(Temperature = scale(Temperature),
            Salinity = scale(Salinity),
            Oxygen = scale(Oxygen),
            Silicate = scale(Silicate),
            NO2 = scale(NO2),
            NO2 = scale(NO3),
            PO4 = scale(PO4)) %>% 
  as.matrix() %>% dist() %>% as.matrix() %>%  normalizeMatrix()


list_normalized_geo_abiotics_dists <- list(
  geo_Dist = geo_Dist,
  geo_Dist_mirrored = geo_Dist_mirrored,
  abioticDist = abioticDist
)


## Now the unnormalized
geo_Dist = df_geo_abiotics %>% 
  transmute(lat_scaled = (Latitude -mean(Latitude)) /sd(Latitude),
            depht_scaled = (Depth -mean(Depth)) /sd(Depth)) %>%
  as.matrix() %>% dist() %>% as.matrix()

geo_Dist_mirrored = df_geo_abiotics %>% 
  mutate(Latitude = abs(Latitude)) %>% 
  transmute(lat_scaled = (Latitude -mean(Latitude)) /sd(Latitude),
            depht_scaled = (Depth -mean(Depth)) /sd(Depth)) %>%
  as.matrix() %>% dist() %>% as.matrix()

abioticDist = df_geo_abiotics %>% 
  transmute(Temperature = scale(Temperature),
            Salinity = scale(Salinity),
            Oxygen = scale(Oxygen),
            Silicate = scale(Silicate),
            NO2 = scale(NO2),
            NO2 = scale(NO3),
            PO4 = scale(PO4)) %>% 
  as.matrix() %>% dist() %>% as.matrix()

min_raw_count = dframe %>% select(Raw.Sequence.Counts) %>% min()
min_raw_count = min_raw_count/1000

list_geo_abiotics_dists <- list(
  geo_Dist = geo_Dist,
  geo_Dist_mirrored = geo_Dist_mirrored,
  abioticDist = abioticDist,
  aitDist = dframe %>%
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
)


## So here we have two different list of distance matrices. 
saveRDS(list_geo_abiotics_dists,paste0(savingdir,'/','list_geo_abiotics_dists'))
saveRDS(list_normalized_geo_abiotics_dists,paste0(savingdir,'/','list_normalized_geo_abiotics_dists'))



