path_demo <- system.file("extdata", "raster", package = "ecotoner")


load(file.path(path_demo, "meta_eg.RData"))

abutt_eg <- raster::raster(file.path(path_demo, "abutt_eg.grd"))
asp201Mean_eg <- raster::raster(file.path(path_demo, "asp201Mean_eg.grd"))
asp201SD_eg <- raster::raster(file.path(path_demo, "asp201SD_eg.grd"))
veg_eg <- raster::raster(file.path(path_demo, "veg_eg.grd"))
elev_eg <- raster::raster(file.path(path_demo, "elev_eg.grd"))


bseABUTtf <- abutt_eg
asp201Mean <- asp201Mean_eg
asp201SD <- asp201SD_eg
gap <- veg_eg
elev <- elev_eg
