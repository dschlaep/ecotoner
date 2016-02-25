library(raster)
#------Get data
#--Settings
transect_type <- 4
longlat <- FALSE

#--Directories
dir.gis <- "/Volumes/BookDuo_12TB/BigData/GIS/Data"
dir.dem <- file.path(dir.gis, "Topography", "NED_USA")
dir.ned1s <- file.path(dir.dem, "NED_1arcsec")
dir.veg <- file.path(dir.gis, "SpeciesVegetation", "GAP_2011_v2")
dir.gap <- file.path(dir.veg, "Nat_GAP_LandCover")

#--Vegetation data
gap <- raster::raster(temp <- file.path(dir.gap, "natgaplandcov_v2_1.img"))
stopifnot(all(res(gap) == get("cell_width_m", envir = etr_vars)))
stopifnot(all(isLonLat(gap) == longlat))
proj_NAD83 <- paste(projection(gap), "+datum=NAD83")

#gap.rat <- GDALinfo(temp, silent=TRUE, returnRAT=TRUE, returnCategoryNames=TRUE, returnStats=TRUE, returnColorTable=TRUE)	#RAT is empty --> read vat.dbf and link back to raster
gap.rat <- foreign::read.dbf(file=paste0(temp, ".vat.dbf"))

#Human impact on land cover:
#CL - NVC_CLASS
#07	- "Agricultural Vegetation": all
#08 - "Developed & Other Human Use" (Mining and developed areas): all
#10 - "Recently Disturbed or Modified": all (i.e., generic disturbance, logged or harvested, disturbed/successional)
#	except 'recently burned' (LEVEL3: 8301, 8302, 8303, 8304)
humanRatValue <- na.omit(gap.rat$Value[gap.rat$CL %in% c("07", "08", "10") & !(gap.rat$LEVEL3 %in% c(8301, 8302, 8303, 8304))])

#Big sagebrush ecosystems:
#Value - gap.rat$LEVEL3	- gap.rat$ECOLSYS_LU:
#489	- 5706	- Inter-Mountain Basins Big Sagebrush Shrubland
#490	- 5307	- Inter-Mountain Basins Big Sagebrush Steppe
#491	- 5308	- Inter-Mountain Basins Montane Sagebrush Steppe
bseRatValue <- c(489, 490, 491)
bse_EndToLeft <- FALSE

#Temperate forests:
tfRatValue <- na.omit(gap.rat$Value[gap.rat$NVC_FORM %in% c("Warm Temperate Forest", "Cool Temperate Forest")]) # i.e., excluding boreal, tropical and flooded forests
tf_EndToLeft <- TRUE

#Abutting cells
bseABUTtf <- raster(file.path(dir.veg, "Types", "gapv2_bseABUTtf.tif"))
raster::crs(bseABUTtf) <- proj_NAD83


#--Topography
elev <- raster(file.path(dir.ned1s, "ned_1s_westernUS_AEANAD83.tif"))
stopifnot(all(res(elev) == get("cell_width_m", envir = etr_vars)))
raster::crs(elev) <- proj_NAD83

if(transect_type == 4){
	#These rasters were calculated (/Users/drschlaep/Documents/drschlaepfer/3_BigData/200907_UofWyoming_PostDoc/GISonE/Data/Topography/NED_USA/NED_1arcsec/3_Calc_Terrain.R) assuming:
	#	- min_slope_with_aspect == 2 * pi/180
	# 	- bandTransect_width_cellN == 200
	asp201Mean <- raster(file.path(dir.ned1s, "terrain", "aspect_201Mean_ned_1s_westernUS_AEANAD83.tif"))
	raster::crs(asp201Mean) <- proj_NAD83
	asp201SD <- raster(file.path(dir.ned1s, "terrain", "aspect_201SD_ned_1s_westernUS_AEANAD83.tif"))
	raster::crs(asp201SD) <- proj_NAD83
}


#-------
dir.out <- "~/Dropbox/Work_Stuff/2_Research/200907_UofWyoming_PostDoc/Product17_EcotoneGradients/3_EcotoneGradient/4_Intermountain_v6/ecotoner/inst/extdata/raster"

ext_crop <- raster::extent(x = -1.52 * 1e6, xmax = -1.50 * 1e6, ymin = 2.03 * 1e6, ymax = 2.05 * 1e6)

veg_eg <- raster::crop(gap, ext_crop, filename = file.path(dir.out, "veg_eg.grd"), prj = TRUE, datatype = 'INT2U', overwrite = TRUE)
abutt_eg <- raster::crop(bseABUTtf, ext_crop, filename = file.path(dir.out, "abutt_eg.grd"), prj = TRUE, datatype = 'INT1S', overwrite = TRUE)
elev_eg <- raster::crop(elev, ext_crop)
elev_eg <- raster::calc(elev_eg, fun = function(x) round(x), filename = file.path(dir.out, "elev_eg.grd"), prj = TRUE, datatype = 'INT2S', overwrite = TRUE)
asp201Mean_eg <- raster::crop(asp201Mean, ext_crop, filename = file.path(dir.out, "asp201Mean_eg.grd"), prj = TRUE, datatype = 'FLT4S', overwrite = TRUE)
asp201SD_eg <- raster::crop(asp201SD, ext_crop, filename = file.path(dir.out, "asp201SD_eg.grd"), prj = TRUE, datatype = 'FLT4S', overwrite = TRUE)


save(
	bse_EndToLeft,
	bseRatValue,
	gap.rat,
	humanRatValue,
	longlat,
	tf_EndToLeft,
	tfRatValue,
	file = file.path(dir.out, "meta_eg.RData"))

#-------
grid_mask1 <- abutt_eg
grid_maskNA <- elev_eg
