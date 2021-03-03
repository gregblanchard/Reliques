
#####################################################################################################
#####################################################################################################
####  Spectral beta diversity from Sentinel 2  ####
#####################################################################################################
#####################################################################################################

#### packages ####

library(raster)
library(rgdal)
library(UScensus2010)
library(rgeos)
library(sp)
library(maptools)
library(sf)
library(RSAGA)
library(SDMTools)
library(biodivMapR)
library(mapview)



#####################################################################################################
####  load image from GEE  ####
#####################################################################################################

#### select clustered image of dry season ####

im_path <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_spectrale_NC/2_seasons"
im_name <- "kmeans_S2_NC_SS_50grp"

# kmeans groups
kmeans_NC_forest_SS <- raster( paste0(im_path, "/", im_name, ".tif") )

#### mask raster with forest #### 
# get forest raster
forest_map_raster <- raster(
  "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/forest_raster/forest_raster_GT.tif")

kmeans_NC_forest_SS <- mask(kmeans_NC_forest_SS, forest_map_raster)

writeRaster(kmeans_NC_forest_SS,  paste0(im_path, "/", im_name, "_forest.tif"),
            datatype="INT1U", options="COMPRESS=LZW")

#### select clustered image of humid season ####
im_name <- "kmeans_S2_NC_SH_50grp"
# kmeans groups
kmeans_NC_forest_SH <- raster( paste0(im_path, "/", im_name, ".tif") )

#### mask raster with forest #### 
# get forest raster
forest_map_raster <- raster(
  "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/forest_raster/forest_raster_GT.tif")

kmeans_NC_forest_SH <- mask(kmeans_NC_forest_SH, forest_map_raster)

writeRaster(kmeans_NC_forest_SH,  paste0(im_path, "/", im_name, "_forest.tif"),
            datatype="INT1U", options="COMPRESS=LZW")

#####################################################################################################
####  load image from GEE with forest mask  ####
#####################################################################################################
im_path <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_spectrale_NC/2_seasons"
im_name <- "kmeans_S2_NC_SS_50grp_forest"
# kmeans groups
kmeans_NC_forest_SS <- raster( paste0(im_path, "/", im_name, ".tif") )

#######################################
#### alpha diversity  ####
#######################################

# raster with dimension = X 10

res_S2 <- res(kmeans_NC_forest_SS)
r <- raster(ext = extent(kmeans_NC_forest_SS), res=c(res_S2*10,res_S2*10))
p <- rasterToPolygons(r) 


#######################################
#### get sampling points ####
#######################################
points_spl_utm <- readRDS(
  "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_spectrale_NC/data_frag/points_spl_utm_landscape_metrics_topo.rds")



#######################################
#### alpha an beta diversity maps ####
#######################################

#### File options ###
Input_Image_File=(paste0(im_path, "/biodivMapR_Convert_BIL/", im_name))
Input_Image_File
#Input_Mask_File   ="D:/Mes Donnees/Google_Cloud/mask_NC_sud"
Output_Dir        = paste0(im_path, "/RESULTS/", im_name)

#### Computing options ###
nbCPU         = 2
MaxRAM        = 0.5
nbclusters    = 50

window_size <- 100

Index_Alpha   = c('Shannon')



################################################################################
##              MAP FUNCTIONAL DIVERSITY METRICS FRic, FEve, FDiv             ##
##          (Villeger et al, 2008 https://doi.org/10.1890/07-1206.1)          ##
################################################################################
## read selected features from dimensionality reduction
Selected_Features <- read.table(Sel_PC)[[1]]
## path for selected components
map_functional_div(Original_Image_File = Input_Image_File, Functional_File = PCA_Files,
                   Selected_Features = Selected_Features, Output_Dir = Output_Dir,
                   window_size = window_size, nbCPU = nbCPU, MaxRAM = MaxRAM,TypePCA = TypePCA)
