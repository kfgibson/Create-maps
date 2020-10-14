#Last updated October 13, 2020
#Kate Gibson
#Purpose: Create maps showing the spatial-temporal distribution of different ASVs and taxonomic groups.
#Input: Phyloseq object with OTU table. 

########################
#load packages
library(phyloseq)
library(compositions)
library(sf)
library(sp)
library(raster)
library(fasterize)
library(tmap)
library(grid)
######################### 

#function to convert count data to compositional data
convert_compositional <- function(phyloseq_object){

  count_tab <- as.data.frame((t(as.matrix(otu_table(phyloseq_object)))))
  count_tab_composition <- count_tab

  for (col in colnames(count_tab)){
    for (row in 1:nrow(count_tab)){
      if (count_tab[row,col] == 0){
        count_tab_composition[row,col] <- (count_tab[row,col]+0.1)/sum(count_tab[,col])
      }
      else {
        count_tab_composition[row,col] <- count_tab[row,col]/sum(count_tab[,col])
      }
    }
  }
  return(count_tab_composition)
}

#clr-transform compositional data & update phyloseq object 
compositional_to_clr <- function(compositional_data, phyloseq_object){
  
  count_tab_clr <- as.data.frame(clr(compositional_data))
  otu_table(phyloseq_object) <- otu_table(count_tab_clr, taxa_are_rows = T)
  
  return(phyloseq_object)
  
}

#parse ASVs: clr-transform all ASVs & marge sample data and otu table 
parse_ASVs <- function(phyloseq_object){
  
  #convert to compositional data 
  count_tab_composition <- convert_compositional(phyloseq_object)
  #clr transformation 
  phyloseq_object <- compositional_to_clr(count_tab_composition, phyloseq_object)
  #merge sample data and otu table & return
  merged<-merge(as.data.frame(sample_data(phyloseq_object)), as.data.frame(t(otu_table(phyloseq_object))), by=0)
  row.names(merged) <- merged$Row.names
  merged<-within(merged, rm("Row.names"))
  return(merged)
  
}

#parse all taxonomic ranks: clr-transform & marge sample data and otu table
parse_ranks <- function(phyloseq_object, merged){
  
  for (rank in rank_names(phyloseq_object)[1:6]){
    
    #aglomerate to the rank level
    ASV_physeq_glom <- tax_glom(phyloseq_object, taxrank=rank)
    #convert to compositional data
    count_tab_composition <- convert_compositional(ASV_physeq_glom)
    #clr transform
    ASV_physeq_glom <- compositional_to_clr(count_tab_composition, ASV_physeq_glom)
    
    #replace taxa names with rank name 
    taxa_names(ASV_physeq_glom) <- paste0(rank, "_", as.character(data.frame(as(tax_table(ASV_physeq_glom), "matrix"))[,rank]))
    
    #merge with previous data frame 
    new_otu<-as.data.frame(t(otu_table(ASV_physeq_glom)))
    merged<-merge(merged, new_otu, by=0)
    row.names(merged) <- merged$Row.names
    merged<-within(merged, rm("Row.names"))
    
  }
  
  return(merged)
  
}

#function to create individual map 
map_data<-function(map_df, min, max, ASV, legend_title, title, plot_legend = F){
  
  map_df<-subset(map_df, !is.na(Latitude))
  plot_locations<-st_as_sf(map_df, coords=c('Latitude', 'Longitude'), crs=4326)
  
  attr(plot_locations$geometry, "bbox") = c(xmin=108, ymin=20, xmax=111, ymax=22)
  
  #make empty grid
  g = st_make_grid(marine, n=c(50,50))
  sf::st_crs(g) <- sf::st_crs(plot_locations)
  #idw interpolation
  int<-gstat::idw(as.formula(paste(ASV,  " ~ 1")), plot_locations, newdata=g)
  
  #convert to raster data
  r<-raster::raster(int, ncols=400,nrows=400)
  crs(r) <- CRS(projargs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  r<-fasterize(int, r, field='var1.pred')
  r.m<-raster::mask(r,marine)
  
  if (plot_legend == T){
    plot<-tm_shape(r.m) + tm_raster(n=10,palette="-RdBu", breaks=seq(min,max,by=0.5), legend.show = T, midpoint=NA, title=legend_title) +
      tm_shape(plot_locations) + tm_dots(size = 0.05) +
      tm_layout(legend.only = TRUE, scale=2, legend.text.size = 0.6, legend.title.size = 1) 
  }
  
  else {
    plot<-tm_shape(r.m) + tm_raster(n=10,palette="-RdBu", breaks=seq(min,max,by=0.5), legend.show = F, midpoint=NA) +
      tm_shape(plot_locations) + tm_dots(size = 0.05) +
      tm_layout(title, title.size=1, outer.margins=c(0.01, 0.01, 0.01, 0.01), title.position = c(0.3,0.95)) +
      tm_legend(legend.outside=TRUE)
  }
  
  return(plot)
}

#function to make grid of maps for an ASV or taxonomic group 
make_grid<-function(group_name, group_label, clr_df){
  
  minimum<-floor(min(clr_df[,group_name]))
  maximum<-ceiling(max(clr_df[,group_name]))  
  
  legend<-map_data(exped1_small, minimum, maximum, group_name, group_label, expression("January (0.22-10" ~ mu*"m)"), T)
  plot1<-map_data(exped1_small, minimum, maximum, group_name, group_label, expression("January (0.22-10" ~ mu*"m)"))
  plot2<-map_data(exped1_large, minimum, maximum, group_name, group_label, expression("January (10-200" ~ mu*"m)"))
  plot3<-map_data(exped2_small, minimum, maximum, group_name, group_label, expression("Feb-March (0.22-10" ~ mu*"m)"))
  plot4<-map_data(exped2_large, minimum, maximum, group_name, group_label, expression("Feb-March (10-200" ~ mu*"m)"))
  plot5<-map_data(exped3_small, minimum, maximum, group_name, group_label, expression("April (0.22-10" ~ mu*"m)"))
  plot6<-map_data(exped3_large, minimum, maximum, group_name, group_label, expression("April (10-200" ~ mu*"m)"))
  
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow=2,ncol=4, widths=c(2.66,2.66,2.66,2), heights=c(5,5))))
  print(plot1, vp=viewport(layout.pos.col = 1, layout.pos.row=1))
  print(plot3, vp=viewport(layout.pos.col = 2, layout.pos.row=1))
  print(plot5, vp=viewport(layout.pos.col = 3, layout.pos.row=1))
  print(legend, vp=viewport(layout.pos.col = 4, layout.pos.row=1:2))
  print(plot2, vp=viewport(layout.pos.col = 1, layout.pos.row=2))
  print(plot4, vp=viewport(layout.pos.col = 2, layout.pos.row=2))
  print(plot6, vp=viewport(layout.pos.col = 3, layout.pos.row=2))
  
}

#################################################
#################################################
#MAIN

#set working directory 
setwd("C:/Users/kfgib/Documents/sfuvault_sync/Thesis/Data/18s project/Beibu_gulf/16S/maps_github/")
#read in phyloseq object (saved as rds file)
ASV_physeq_filt <- readRDS("ASV_physeq_filt.rds")
#process at the ASV-level
clr_df_ASVs <- parse_ASVs(ASV_physeq_filt)
#process the taxonomic ranks 
clr_df_all_levels <- parse_ranks(ASV_physeq_filt, clr_df_ASVs)

#subset for each expedition and size 
exped1_small <- subset(clr_df_all_levels, Size == 0.2 & Expedition == 1)
exped1_large <- subset(clr_df_all_levels, Size == 10 & Expedition == 1)
exped2_small <- subset(clr_df_all_levels, Size == 0.2 & Expedition == 2)
exped2_large <- subset(clr_df_all_levels, Size == 10 & Expedition == 2)
exped3_small <- subset(clr_df_all_levels, Size == 0.2 & Expedition == 3)
exped3_large <- subset(clr_df_all_levels, Size == 10 & Expedition == 3)

#read in marine shape file 
marine <- st_read("ne_10m_geography_marine_polys.shp")
attr(marine$geometry, "bbox") = c(xmin=108, ymin=20, xmax=111, ymax=22)

#sample map
make_grid("Order_Flavobacteriales", "Flavobacteriales", clr_df_all_levels)

