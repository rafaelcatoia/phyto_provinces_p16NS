## grum2comp
## inputs:
#### dflonger = stacked dataset
#### SampleIdVar = name of the variable containing the sample ID
#### bioGranVar = name of the variable that contains the biological unity to be used
################# i.e.   Family,Genus,Species,ASV...
#### CountingVar = name of the variable containing the counts
#### if fillPercentage = False : Each row will be a composition, 
######################## but if a certain bioGranVar is not present for a specific sample it will be filled with 0
#### if fillPercentage = True : After creating the composition, it adds the smallest percentage of the dataset * 1e-4
######################## to all places and than re-scale to get compositions in the row. 

grump2comp <- function(
    dfLonger,
    sampleIdVar='SampleID',
    bioGranVar ='ASV_name',
    countingVar='Raw.Sequence.Counts',
    proportion=T,
    rounding=T,
    fillPercentage=F,
    constant=1e-2){
  
  dfOutput = dfLonger %>% 
    select(dplyr::all_of(c(sampleIdVar,bioGranVar,countingVar))) %>% 
    group_by_at(c(sampleIdVar,bioGranVar)) %>% 
    summarise_at(vars(any_of(countingVar)),.funs = sum) %>% 
    tidyr::pivot_wider(id_cols = dplyr::any_of(sampleIdVar),
                       values_from = any_of(countingVar),
                       names_from= dplyr::any_of(bioGranVar),values_fill = 0) %>%
    ungroup()
  
  if(proportion){
    dfLonger = dfLonger %>% 
      mutate(across(where(is.numeric))/rowSums(across(where(is.numeric))))
    
    if(fillPercentage){
      vetAux = dfOutput[,-1] %>% unlist()
      vetAux=vetAux[vetAux!=0]
      addValue=min(vetAux)*constant
      
      dfOutput = dfOutput %>% 
        mutate(across(where(is.numeric),function(x) {x+addValue}))%>% 
        mutate(across(where(is.numeric))/rowSums(across(where(is.numeric))))
    }
  }else{
    if(rounding){
      dfOutput = dfOutput %>% 
        mutate_if(is.numeric,round)
    }
  }
  
  return(dfOutput)
}

# Example
#dfTest <- grump2comp(dfZoop,SampleIdVar = 'SampleID',bioGranVar = 'ASV_name',CountingVar = 'Raw.Sequence.Counts',fillPercentage =T)
#adist <- robCompositions::aDist(x = dfTest[,-1])
