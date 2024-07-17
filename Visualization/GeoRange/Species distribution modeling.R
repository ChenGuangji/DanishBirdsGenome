## Visualizing the spatial ranges of Europe birds by ENMeval2.0 (Kass et al, 2021)
## Author: Xu Jiaqi; Email: Xujiaqi09@zju.edu.cn

########## Importing the Necessary R packages
library(ENMeval)
library(raster)
library(dplyr)
library(RStoolbox)
library(rasterVis)
library(viridis)
library(rgbif)
library(spocc)
library(taxize)

########## Get Occurances ##########
args <- commandArgs(T)
species <- args[1]
re_back <- name_backbone(gsub("_"," ",species))

## Get the accuracy scientific name on gbif
taxon_id <- ifelse(re_back$status == "SYNONYM", re_back$acceptedUsageKey, re_back$usageKey)
rename_species <- gbif_name_usage(taxon_id)$species

bv <- occ(rename_species, 'gbif', limit=10000, has_coords=TRUE)
occs <- as.data.frame(bv$gbif$data[[1]][,2:3])
occs <- occs[!duplicated(occs),]

result <- paste0("01.data/01.occs/", species, "_occs.txt")
write.table(occs, file = result, quote = F, sep = "\t", col.names = T, row.names = F)

########## Data Process ##########
## Set a random seed 
set.seed(48)

## Import the occurance data
occs <- read.table(paste0("01.data/01.occs/", species, "_occs.txt"), header = T)
# Limit the number of samples to avoid excessive computation time.
if(nrow(occs) > 3000){
  occs <- occs[sample(1:nrow(occs), 3000),]
}

## Filter the collinearity bioclimatical variables
need_bio <- c(4, 8, 10, 11, 15, 16, 18)

## Import the present climate data for modeling
envs.files <- paste0("01.data/02.climate/WC_2.5m/wc2.1_2.5m_bio_", need_bio, ".tif")
pre_envs <- stack(lapply(envs.files, raster))
names(pre_envs) <- paste0("bio_", need_bio)

## For masking somewhere incorrect
mask_envs <- raster("01.data/02.climate/BA_2.5m/bio_1.tif")

origin(pre_envs) <- c(0,0)
origin(mask_envs) <- c(0,0)

## Crop and mask the raster
envs_extent <- extent(pre_envs)
mask.envs_extent <- extent(mask_envs)
crop_buf <- extent(max(envs_extent[1], mask.envs_extent[1]), min(envs_extent[2], mask.envs_extent[2]),
                   max(envs_extent[3], mask.envs_extent[3]), min(envs_extent[4], mask.envs_extent[4]))

pre_envs <- crop(pre_envs, crop_buf)
mask_envs <- crop(mask_envs, crop_buf)
pre_envs <- mask(pre_envs, mask_envs)


########## Model Construction ##########
## Selection the range of raster according to the species occurances
bufsize <- 10
sp1 <- occs %>% na.omit() %>% SpatialPoints()
bb_buf2 <- extent(bbox(sp1)[1] - bufsize, bbox(sp1)[3] + bufsize, 
                  bbox(sp1)[2] - bufsize, bbox(sp1)[4] + bufsize)

envs <- crop(pre_envs, bb_buf2)
mask_envs <- crop(mask_envs, bb_buf2)

## Delete the duplicated occurances
occs.cells <- raster::extract(envs[[1]], occs, cellnumbers = TRUE)
occs.cellDups <- duplicated(occs.cells[,1])
occs <- occs[!occs.cellDups,]

## Select the background
occs.sf <- sf::st_as_sf(occs, coords = c("longitude", "latitude"), crs = raster::crs(envs))
eckertIV <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
occs.sf <- sf::st_transform(occs.sf, crs = eckertIV)

occs.buf <- sf::st_buffer(occs.sf, dist = 500000) %>% 
  sf::st_union() %>% 
  sf::st_sf() %>%
  sf::st_transform(crs = raster::crs(envs))

envs.bg <- crop(envs, occs.buf)
envs.bg <- mask(envs.bg, occs.buf)

bg <- dismo::randomPoints(envs.bg[[1]], n = 10000) %>% as.data.frame() ## select the background points randomly
colnames(bg) <- colnames(occs)

## Dimension reduction by rasterPCA()
envs_pca <- rasterPCA(envs, spca=TRUE, maskCheck=FALSE, nSamples=10000)
envs_pca$map <- subset(envs_pca$map, 1:6)

e.mx <- ENMevaluate(occs = occs, envs = envs_pca$map, bg = bg, 
                    algorithm = 'maxnet', partitions = 'randomkfold', 
                    tune.args = list(fc=c("L", "LQ", "LQH", "LQHP", "LQHPT"),
                                     rm=c(1:5)))
## Get the best model
res <- eval.results(e.mx)
opt.aicc <- res %>% filter(delta.AICc == 0)
opt.seq <- res %>% filter(or.10p.avg == min(or.10p.avg)) %>% filter(auc.val.avg == max(auc.val.avg))

mod_best <- e.mx@models[[opt.seq$tune.args]]

########## Niche Prediction ##########
## Crop the raster to Europe range
bufsize2 <- 5
xmin <- -10.34
xmax <- 66.1
ymin <- 34.55
ymax <- 81.15

envs_extent <- extent(pre_envs)
bb_buf2 <- extent(max(envs_extent[1], xmin-bufsize2), min(envs_extent[2], xmax+bufsize2),
                  max(envs_extent[3], ymin-bufsize2), min(envs_extent[4], ymax+bufsize2))

pred_envs <- crop(pre_envs, bb_buf2)

## Dimension reduction
pred_envs_pca <- rasterPCA(pred_envs, spca=TRUE, maskCheck=FALSE, nSamples=10000)
pred_envs_pca$map <- subset(pred_envs_pca$map, 1:6)

## Prediction
map_pred <- enm.maxnet@predict(mod_best, pred_envs_pca$map, list(pred.type = "logistic", doClamp = FALSE))
plot(map_pred)