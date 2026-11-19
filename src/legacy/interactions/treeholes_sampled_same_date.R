library(tidyverse)

library(robCompositions)
library(philentropy)





ST_metadata = read.table("Metadata_singleTree_wDirich-K10.tsv", header = T, sep = "\t")

ST_metadata_filter <- ST_metadata[ST_metadata$Location %in% c(1:10),]
ST_metadata_filter <- ST_metadata_filter[!is.na(ST_metadata_filter$Dirichlet),]


dat = c()
dir = c()
loc = c()

for (i in 1:10)
{
  
  if (length(unique(ST_metadata_filter[ST_metadata_filter$Dirichlet == i,"Location"]))>1)
  {
    ST_metadata_Location <- ST_metadata_filter %>% filter(Dirichlet == i)
    dates <- names(table(ST_metadata_Location$R.date)[table(ST_metadata_Location$R.date) > 1])
    
    ST_metadata_Date <- ST_metadata_Location %>%  filter(R.date %in% dates)
    
    
    for (j in dates)
    {
      locs = unique((ST_metadata_Location %>% filter(R.date == j))$Location)
      if(length(locs > 1))
      {
        dat = c(dat, j)
        dir = c(dir,i)
        loc = c(loc, paste(locs, collapse = "_"))
        
      }
      
      
    }
    
    
  }
    
}

all_d <- data.frame(Dir = dir, Loc = loc, Date = dat)

ST_metadata_filter %>% filter((Dirichlet %in% all_d$Dir) & (R.date %in% all_d$Date))
