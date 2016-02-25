context("Locate search points")

path_demo <- system.file("extdata", "raster", package = "ecotoner")
abutt_eg <- raster::raster(file.path(path_demo, "abutt_eg.grd"))
elev_eg <- raster::raster(file.path(path_demo, "elev_eg.grd"))

initfile1 <- "test_initpoints_reference_01.rds"
initfile2 <- "test_initpoints_reference_02.rds"
initfile3 <- "test_initpoints_reference_03.RData"

#--- Generate reference object
temp <- get_transect_search_points(N = 10, grid_mask1 = abutt_eg, grid_maskNA = elev_eg, initfile = initfile1, seed = 123, verbose = FALSE)
if (!file.exists(initfile3)) save(temp, file = initfile3)
# dput(temp)
temp_exp <- new("SpatialPoints"
    , coords = structure(c(-1502970, -1503180, -1501230, -1501350, -1507590, 
-1500360, -1507800, -1504530, -1501410, -1500630, 2039479.43228714, 
2041009.43228714, 2030719.43228714, 2030419.43228714, 2033719.43228714, 
2049289.43228714, 2032399.43228714, 2036509.43228714, 2030359.43228714, 
2045209.43228714), .Dim = c(10L, 2L), .Dimnames = list(c("19", 
"17", "16", "20", "14", "18", "13", "12", "11", "15"), c("x", 
"y")))
    , bbox = structure(c(-1507800, 2030359.43228714, -1500360, 2049289.43228714
), .Dim = c(2L, 2L), .Dimnames = list(c("x", "y"), c("min", "max"
)))
    , proj4string = new("CRS"
    , projargs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
)
)


test_that("Locate search points", {
	expect_is(get_transect_search_points(N = 10, grid_mask1 = abutt_eg, grid_maskNA = elev_eg, seed = 123, verbose = FALSE), "SpatialPoints")
	expect_equal(get_transect_search_points(N = 10, grid_mask1 = abutt_eg, grid_maskNA = elev_eg, seed = 123, verbose = FALSE), temp_exp)
	expect_equal(get_transect_search_points(N = 10, grid_mask1 = abutt_eg, grid_maskNA = elev_eg, initfile = initfile1, seed = 123, verbose = FALSE), readRDS(initfile1))
	expect_equal_to_reference(get_transect_search_points(N = 10, grid_mask1 = abutt_eg, grid_maskNA = elev_eg, seed = 123, verbose = FALSE), initfile2)
	expect_error(get_transect_search_points(N = 0, grid_mask1 = abutt_eg, grid_maskNA = elev_eg, seed = 123, verbose = FALSE), "N must be > 0")
	expect_error(get_transect_search_points(N = 10, grid_mask1 = NULL, grid_maskNA = elev_eg, seed = 123, verbose = FALSE), "not of class 'Raster'")
	expect_error(get_transect_search_points(N = 10, grid_mask1 = abutt_eg, grid_maskNA = elev_eg, initfile = initfile3, seed = 123, verbose = FALSE), "extension of argument 'initfile' must be '.rds'")
})
