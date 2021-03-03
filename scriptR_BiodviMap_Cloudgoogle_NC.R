#######Packages
library("raster")
library("devtools")
#devtools::install_github('jbferet/biodivMapR')
library("biodivMapR")
library("mapview")



## rasterOptions(tmpdir="D:/A_traiter/R_temp")
## system.file(package = "biodivMapR")
###############
############### 1/ Lire les donn?es 
###############

# note that there is a function in BiodivMapR to merge multiple images : build_image_from_list 


# S2 images
im_path <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_spectrale_NC/2_seasons/S2_images"
im_name <- "S2_NC_SH_int_crop"
im <- stack(paste0(im_path, "/", im_name, ".tif"))

# forest raster 

forest_map_raster <- raster(
  "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/forest_raster/forest_raster_GT.tif")

#conversion en BIL

#TRansfo en format ENVI pour passer convertir en bil apr?s
im_envi <- writeRaster(im, paste0(im_path, "/", im_name, ".tif"), format="ENVI", datatype = "INT2U", options="INTERLEAVE=BSQ", overwrite=TRUE) 

im_envi=stack(paste0(im_path, "/", im_name, ".envi"))
plot(im_envi)

# # convert the image using Convert.Raster2BIL if not in the proper format
Input_Image_File  <- raster2BIL(Raster_Path = paste0(im_path, "/", im_name, ".envi"),
                               Sensor = "unknown",
                               Convert_Integer = FALSE,
                               Multiplying_Factor= 10000,
                               #Multiplying_Factor_Last =1,
                               Output_Dir = FALSE)

im_bil=stack(paste0(im_path, "/biodivMapR_Convert_BIL/", im_name))

#### same for forest map ####

# GRD TO TIFF !!!! à faire !!!! 

im_envi <- writeRaster(im, paste0(im_path, "/", im_name, ".tif"), format="ENVI",  options="INTERLEAVE=BSQ", overwrite=TRUE) 

im_envi=stack(paste0(im_path, "/", im_name, ".envi"))
plot(im_envi)

# # convert the image using Convert.Raster2BIL if not in the proper format
Input_Image_File  <- raster2BIL(Raster_Path = paste0(im_path, "/", im_name, ".envi"),
                                Sensor = 'SENTINEL_2A',
                                Convert_Integer = FALSE,
                                Multiplying_Factor= 10000,
                                #Multiplying_Factor_Last =1,
                                Output_Dir = FALSE)

im_bil=stack(paste0(im_path, "/biodivMapR_Convert_BIL/", im_name))


Input_Image_File=(paste0(im_path, "/biodivMapR_Convert_BIL/", im_name))
Input_Image_File
#Input_Mask_File   ="D:/Mes Donnees/Google_Cloud/mask_NC_sud"
Output_Dir        = paste0(im_path, "/RESULTS/", im_name)
Input_Mask_File   = 
#
#red=raster('J:/Sentinel2/Gilles/biodivMapR_Convert_BIL/Saotome_S2', band=3)
#nir=raster('J:/Sentinel2/Gilles/biodivMapR_Convert_BIL/Saotome_S2', band=7)
#nir=calc(nir,function(x)ifelse(x>0,x,NA))
#range(values(nir), na.rm = T)
#ndvi=(nir-red)/(nir+red)
#toto=raster('J:/Sentinel2/Gilles/biodivMapR_Convert_BIL/Saotome_S2', band=7)
#mapview(toto)


window_size = 10
TypePCA = 'SPCA'
FilterPCA = TRUE

########Computing options
nbCPU         = 2
MaxRAM        = 0.5
nbclusters    = 80

# 1- MASK : Filter data in order to discard non vegetated / shaded / cloudy pixels
NDVI_Thresh = 0.53
Blue_Thresh = 500
NIR_Thresh  = 1500
print("PERFORstM RADIOMETRIC FILTERING")
Input_Mask_File = perform_radiometric_filtering(Input_Image_File,Input_Mask_File,Output_Dir,
                                                NDVI_Thresh = NDVI_Thresh, Blue_Thresh = Blue_Thresh,
                                                NIR_Thresh = NIR_Thresh)


############## PCA
print("PERFORM PCA ON RASTER")
PCA_Output        = perform_PCA(Input_Image_File, Input_Mask_File, Output_Dir,
                                TypePCA = TypePCA, FilterPCA = FilterPCA, nbCPU = nbCPU,MaxRAM = MaxRAM, Continuum_Removal = TRUE)


# Path for the PCA raster
PCA_Files         = PCA_Output$PCA_Files
# number of pixels used for each partition used for k-means clustering
Pix_Per_Partition = PCA_Output$Pix_Per_Partition
# number of partitions used for k-means clustering
nb_partitions     = PCA_Output$nb_partitions

# Path for the updated mask
Input_Mask_File   = PCA_Output$MaskPath
# parameters of the PCA model
PCA_model         = PCA_Output$PCA_model
# definition of spectral bands to be excluded from the analysis
SpectralFilter    = PCA_Output$SpectralFilter


print("Select PCA components for diversity estimations")
select_PCA_components(Input_Image_File, Output_Dir, PCA_Files, File_Open = TRUE)

print("MAP SPECTRAL SPECIES")
map_spectral_species(Input_Image_File, Output_Dir, PCA_Files, PCA_model, SpectralFilter, Input_Mask_File,
                     Pix_Per_Partition, nb_partitions, nbCPU=nbCPU, MaxRAM=MaxRAM, 
                     nbclusters = nbclusters, TypePCA = TypePCA, Continuum_Removal = TRUE)
###################################
####"#alpha an beta diversity maps

print("MAP ALPHA DIVERSITY")
# Index.Alpha   = c('Shannon','Simpson')
Index_Alpha   = c('Shannon')
map_alpha_div(Input_Image_File, Output_Dir, window_size, nbCPU=nbCPU,
              MaxRAM=MaxRAM, Index_Alpha = Index_Alpha, nbclusters = nbclusters)

print("MAP BETA DIVERSITY")
map_beta_div(Input_Image_File, Output_Dir, window_size, nb_partitions=nb_partitions,
             nbCPU=nbCPU, MaxRAM=MaxRAM, nbclusters = nbclusters)

