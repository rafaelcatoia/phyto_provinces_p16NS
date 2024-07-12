########################################################
## Script to filter only phyto ASVs
########################################################

### Getting data & path
root <- rprojroot::has_file(".git/index")
datadir = root$find_file("data")
savingdir = root$find_file("saved_files")
path_df = root$find_file('data/grump_asv_long.csv')
df <- data.table::fread(path_df)

## Running the code Yubin sent us to subset only the Phytoplaktons --------------------------------

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

## Saving file:: 
data.table::fwrite(x = phyto.all,file = paste0(datadir,'/','grump_phyto.csv'))

## Now move to Script 01 to create the other files.

