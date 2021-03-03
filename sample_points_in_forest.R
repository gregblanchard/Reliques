
library(spatial)
library(sp)
library(sf)
library(biodivMapR)
library(reticulate)
library(rgee)
library(googledrive)
library(spatial)
library(rgdal)
library(raster)
library(rgdal)
library(UScensus2010)
library(rgeos)
library(maptools)
library(RSAGA)
library(SDMTools)

#####################################################################################################
#####################################################################################################
####  Analysis of the signature of forest fragmentation on spectral beta diversity (Sentinel 2)  ####
#####################################################################################################
#####################################################################################################

# #### get previous points sample #### 
# sample_pca_S2_NC_year <- readRDS(paste0(local_path, ".rds"))
# nrow(sample_pca_S2_NC_year)
# 
# #### plot smpled points #### 
# plot(sample_pca_S2_NC_year[1:100,])
# 
# # get NC map
# NC <- st_read(dsn = '/home/thesardfou/Documents/GIS/divers/contour_grande_terre/contour_grande_terre_clean.shp')
# NC <- st_transform(NC, crs = "+proj=longlat +datum=WGS84 +no_defs")
# 
# plot(NC$geometry, axes = T)
# plot(sample_pca_S2_NC_year$geometry, add = T)

#####################################################################################################
#####################################################################################################
####  Compute spatial metrics  ####
#####################################################################################################
#####################################################################################################


#####################################################################################################
####  load data  ####
#####################################################################################################

# forest polygones
# forest_map <- readOGR(
#   "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/forest_NC.shp")
# forest_map@data$id_tmp <- 1:dim(forest_map@data)[1]
# forest_map_utm <- spTransform(forest_map, CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# all_patch_areas <- gArea(forest_map_utm,  byid=TRUE)
# forest_map@data$area <- forest_map_utm@data$area <- all_patch_areas
# saveRDS(forest_map, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/forest_NC.rds")

forest_map <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/forest_NC.rds")
forest_map_utm <- spTransform(forest_map, CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# plot(forest_map)

forest_map_utm$id_tmp <- 1:nrow(forest_map_utm)

#####################################################################################################
####  Sample random points in forest  ####
#####################################################################################################
#### 20 000 pts ####
# points_spl <- spsample(forest_map_utm, 20000, type = "stratified")
# saveRDS(points_spl, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_spectrale_NC/points_sample/points_spl.rds")

#### 30 000 pts ####
 points_spl <- spsample(forest_map_utm, 30000, type = "stratified")
# saveRDS(points_spl, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_spectrale_NC/points_sample/points_spl.rds")

# 
# points_spl_utm <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_spectrale_NC/points_sample/points_spl.rds")

points_spl_utm <- SpatialPointsDataFrame(points_spl_utm, data.frame(ID=1:length(points_spl_utm)))

#####################################################################################################
####  Get environmental variables and spatial metrics from random points in forest  ####
#####################################################################################################

# points_spl <- sample_pca_S2_NC_year
# points_spl <- as(points_spl, "Spatial")
# points_spl_utm <- spTransform(points_spl, CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#### extract spatial variables ####

# ID patch
patchs_over_points_utm <- over(points_spl_utm,forest_map_utm) 
points_spl_utm$id_patch <- patchs_over_points_utm$id_tmp
# saveRDS(patchs_over_points_utm, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/patchs_over_points_utm.rds")

# patch area
points_spl_utm$area <- patchs_over_points_utm$area

# keep only points in forest polygons
points_spl_utm <- points_spl_utm[!is.na(points_spl_utm$id_patch),]

#### distance to edge #### 

cl <- makeCluster(3)
registerDoSNOW(cl)
clusterEvalQ(cl, library(sp))
clusterEvalQ(cl, library(rgeos))
pb <- txtProgressBar(max = , style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
system.time(
    dist_edge <- foreach(i = 1:n_pts, .combine = c, 
                         .options.snow = opts) %dopar%
        {
      pts_tmp <- points_spl_utm[i,]
      patch_tmp <- forest_map_utm[forest_map_utm$id_tmp <- pts_tmp$id_patch,]
      forest_edge <- as(patch_tmp, "SpatialLines") 
      dist_tmp <- gDistance(pts_tmp, forest_edge, byid=TRUE)
    }
)
close(pb)
stopCluster(cl) 
 
points_spl_utm$dist_edge <- dist_edge


###########################################################################
#####   landscape metrics #####  
###########################################################################

forest_map_raster <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/forest_raster/forest_raster_GT.tif")

##### only points in forest raster ####
points_spl <- spTransform(points_spl_utm, CRS("+proj=longlat +datum=WGS84 +no_defs"))
pts_in_forest_raster <- extract(forest_map_raster, points_spl)
pts_in_forest_raster[is.na(pts_in_forest_raster)] <- 0
points_spl <- points_spl[pts_in_forest_raster == 1 ,]
points_spl_utm <- spTransform(points_spl, CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs"))
#####   patch metrics #####  
mat = ConnCompLabel(forest_map_raster)
patch_stats = PatchStat(mat)
id_patch_in_grid = extract(mat,points_spl)
patch_stats$patchID
# !!!! in patch_stats, patchs id begin at 0 not 1 !!!!
pts_patch_stats = patch_stats[id_patch_in_grid+1,]

# saveRDS(pts_patch_stats,
#   "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/pts_patch_stats_metrics.rds")

#####  landscape metrics for diferent landscape buffers ##### 

# buffer values 
buffers = c(250,500,1000,2000)
# buffers = c(250,500)
n_pts <- nrow(points_spl_utm)
n_pts <- 10
library(doSNOW)


local_stats = c()
for (j in 1:length(buffers)){
  dist_tmp = buffers[j]
  tab_tmp = c()
  cl <- makeCluster(3)
  registerDoSNOW(cl)
  clusterEvalQ(cl, library(sp))
  clusterEvalQ(cl, library(rgeos))
  clusterEvalQ(cl, library(raster))
  clusterEvalQ(cl, library(SDMTools))
  pb <- txtProgressBar(max = , style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  system.time(
    tab_tmp <- foreach(i = 1:n_pts, .combine = rbind, 
                                                          .options.snow = opts) %dopar%
      {
        points_spl_buffer_tmp <- gBuffer(points_spl_utm[i,], width=dist_tmp, byid=TRUE)
        points_spl_buffer_tmp <- spTransform(points_spl_buffer_tmp, CRS("+proj=longlat +datum=WGS84 +no_defs"))
        local_habitat =  crop(forest_map_raster, points_spl_buffer_tmp)
        local_habitat =  mask(local_habitat, points_spl_buffer_tmp)
        local_stats_tmp = ClassStat(local_habitat)
        local_stats_tmp = local_stats_tmp[2,]
        }
  )
  close(pb)
  stopCluster(cl) 
  colnames(tab_tmp) <- paste0( colnames(tab_tmp), "_", buffers)
  if(is.null(local_stats)) local_stats <- matrix(nrow = nrow(tab_tmp), ncol = 0)
  local_stats <- cbind(local_stats, tab_tmp)
}


# saveRDS(local_stats,
# "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/landscape_metrics.rds")


#### extract environmental variables ####

# get values from gdal script (/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/gis_script/var_topo)

dem <- raster("/home/thesardfou/Documents/GIS/MNT/mnt10_GT.tif")
slope <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/slope_GT.tif")
aspect <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/aspect_GT.tif")
curv <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/curvature_GT.tif")
twi <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/twi_GT.tif")

pts_alt <- extract(dem, points_spl_utm)
pts_slope <- extract(slope, points_spl_utm)
pts_aspect <- extract(aspect, points_spl_utm)
pts_cur <- extract(curv, points_spl_utm)
pts_twi <- extract(twi, points_spl_utm)

points_spl_utm$elevation <- pts_alt
points_spl_utm$slope <- pts_slope
points_spl_utm$aspect <- pts_aspect
points_spl_utm$curvature <- pts_cur

# distance to river network ==> récuérer le reseau hydrographique sur le nas! 
# channel <- readOGR(".shp")
# all_dist_to_channel <- gDistance(points_spl, channel, byid=TRUE)
# dist_to_channel <- apply(points_spl, 2,min)


##### UM and non-UM areas ####

UM_GT <- readOGR("/home/thesardfou/Documents/GIS/divers/NC_GT_UM_areas.shp")
UM_GT_UTM <-  spTransform(UM_GT, CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs "))

pts_in_UM <- over( points_spl_utm, UM_GT_UTM)

substrat <- ifelse(is.na(pts_in_UM$gid), "non_UM", "UM")

points_spl_utm$substrat <- substrat

saveRDS(points_spl_utm,
  "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_spectrale_NC/data_frag/points_spl_utm_landscape_metrics_topo2.rds")

