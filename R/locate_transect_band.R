
#' Draws a rectangular band (of #cells=width) around transect defined by two points (start, end)
set_stband_polygon <- function(start, end, width_n, res_m, Longitude_WGS84 = 0) {
	#---Test input
	temp_prj <- raster::crs(start)
	if (!inherits(temp_prj, "CRS") || !inherits(raster::crs(end), "CRS")) {
		stop("ecotoner::set_stband_polygon(): argument(s) 'start' and/or 'end' lack a 'CRS'")
	}

	#---Begin function calculations
	#project grid onto a Cartesian surface first: UTM
	crs_UTM <- sp::CRS(paste0("+proj=utm +zone=", (floor((Longitude_WGS84 + 180)/6) %% 60) + 1, " +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")) #doesn't work for couple of exceptional areas (Svalbard and parts of Norway)
	
	#project transect-defining points
	pt.startUTM <- sp::spTransform(start, CRS = crs_UTM)
	pt.endUTM <- sp::spTransform(end, CRS = crs_UTM)
	
	#coordinates of originating transect around which band is laid
	startxy <- sp::coordinates(pt.startUTM)
	endxy <- sp::coordinates(pt.endUTM)
	
	#calculate vectors to new points
	r <- res_m * width_n/2
	vecX <- endxy[1] - startxy[1]
	vecY <- endxy[2] - startxy[2]
	declin <- atan2(vecY, vecX)
	if (isTRUE(all.equal(vecX, 0)) || isTRUE(all.equal(vecY, 0))) {
		declin <- pi/2 - declin
	}
	deltaX <- r * sin(declin)
	deltaY <- r * cos(declin)
	
	#calculate delineating points
	startxy1 <- sp::coordinates(data.frame(startxy[1] + deltaX, startxy[2] - deltaY))
	startxy2 <- sp::coordinates(data.frame(startxy[1] - deltaX, startxy[2] + deltaY))
	endxy1 <- sp::coordinates(data.frame(endxy[1] - deltaX, endxy[2] + deltaY))
	endxy2 <- sp::coordinates(data.frame(endxy[1] + deltaX, endxy[2] - deltaY))
	
	# check that opposite sides are of equal length, i.e., that band transect is a parallelogram
	if (any(!isTRUE(all.equal(sp::spDistsN1(startxy1, endxy1), sp::spDistsN1(startxy2, endxy2))),
			!isTRUE(all.equal(sp::spDistsN1(startxy1, endxy2), rect.length <- sp::spDistsN1(startxy, endxy))),
			!isTRUE(all.equal(sp::spDistsN1(startxy2, endxy1), rect.length)),
			!isTRUE(all.equal(sp::spDistsN1(startxy1, startxy2), rect.width <- 2*r)),
			!isTRUE(all.equal(sp::spDistsN1(endxy1, endxy2), rect.width)))) {
		stop("ecotoner::set_stband_polygon(): calculated band transect is not a parallelogram")
	}
	# check that all angles between sides are normal, i.e., that band transect is a rectangle
	tol <- res_m / 100	# this needs to be only a rectangle relative to the size of the cells and not absolute
	if (any(!isTRUE(all.equal(as.numeric(tcrossprod(startxy2 - startxy1, endxy2 - startxy1)), 0, tolerance = tol)),
			!isTRUE(all.equal(as.numeric(tcrossprod(startxy1 - startxy2, endxy1 - startxy2)), 0, tolerance = tol)),
			!isTRUE(all.equal(as.numeric(tcrossprod(startxy2 - endxy1, endxy2 - endxy1)), 0, tolerance = tol)),
			!isTRUE(all.equal(as.numeric(tcrossprod(endxy1 - endxy2, startxy1 - endxy2)), 0, tolerance = tol)))) {
		stop("ecotoner::set_stband_polygon(): calculated band transect ist not a rectangle")
	}
	
	#make polygon
	tpoints <- sp::rbind.SpatialPoints(sp::SpatialPoints(startxy1), sp::SpatialPoints(startxy2), sp::SpatialPoints(endxy1), sp::SpatialPoints(endxy2))
	raster::crs(tpoints) <- crs_UTM
	pts_corner <- sp::spTransform(tpoints, CRS = temp_prj)
	tpoly_UTM <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(sp::rbind.SpatialPoints(tpoints, tpoints[1,]), hole=FALSE)), ID=1)), proj4string = crs_UTM)
	tpoly_orig <- sp::spTransform(tpoly_UTM, CRS = temp_prj)

	list(tpoly_UTM = tpoly_UTM, tpoly_orig = tpoly_orig, crs_UTM = crs_UTM, center = startxy1, pts_corner = pts_corner, declin = declin)
}

#' @export
extract_tband_grids <- function(tpoly_UTM, width_n, crs_UTM, declin, center, res_rotate_grid = 8, grid1, grid2 = NULL, fun1, fun2 = NULL) {
	#the calls to rasterize() are the most time consuming calls of this function
	
	#project grid onto a Cartesian surface first: UTM
	to_UTM <- raster::projectExtent(grid1, crs = crs_UTM)
	raster::res(to_UTM) <- raster::res(grid1)
	
	#extract grid values
	grid1UTMr2 <- raster::disaggregate(to_UTM, fact = res_rotate_grid) #need more cells otherwise after rotation not all 30-m will be covered and would need to be interpolated (not good for categorical data)
	band0 <- raster::rasterToPoints(x = raster::rasterize(x = tpoly_UTM, y = grid1UTMr2, silent = TRUE), spatial = TRUE)
	
	#rotate band so that long side of rectangle lies flat with function elide(), elide requires SpatialPoints
	if (requireNamespace("maptools", quietly = TRUE)) {
		band1 <- maptools::elide(obj = band0, rotate = declin * 180/pi, center = center)
	} else {
		temp <- rotate_coords(sp::coordinates(band0), angle_rad = -declin, center = center)
		band1 <- sp::SpatialPointsDataFrame(coords = data.frame(x = temp[, 1], y = temp[, 2]), data = data.frame(layer = rep(1, nrow(temp))))
	}
	to_UTM2 <- raster::raster(raster::extent(band1), crs = crs_UTM)
	raster::res(to_UTM2) <- raster::res(grid1)
	raster::origin(to_UTM2) <- c(0, 0)
	
	#extract grid values for locations within band
	band0_orig <- sp::spTransform(band0, CRS = raster::crs(grid1))	#rows of band1 and band0_orig correspond to each other
	if (!is.null(grid2) && !is.null(fun2)) {
		griddata <- sapply(list(grid1, grid2), FUN = function(grid) raster::extract(grid, band0_orig))
		band.grid1 <- raster::rasterize(x = band1, y = to_UTM2, field = griddata[, 1], fun = fun1)
		if (nrow(band.grid1) < width_n) {
			stop("ecotoner::extract_tband_grids(): failed because extracted band transect is not wide enough: ", nrow(band.grid1), " instead of ", width_n, " cells")
		}
		band.grid2 <- raster::rasterize(x = band1, y = to_UTM2, field = griddata[, 2], fun = fun2)
	} else {
		band.grid1 <- raster::rasterize(x = band1, y = to_UTM2, field = raster::extract(grid1, band0_orig), fun = fun1)
		if (nrow(band.grid1) < width_n) {
			stop("ecotoner::extract_tband_grids(): failed because extracted band transect is not wide enough: ", nrow(band.grid1), " instead of ", width_n, " cells")
		}
		band.grid2 <- NULL
	}
	
	list(grid1 = band.grid1, grid2 = band.grid2)
}

#' @export
locate_candidate_transect_band <- function(stline, width_N, res_rotate_grid, elevCropped, elev, gap, bseRatValue, tfRatValue, seed = NULL) {	
	stband_polygon <- set_stband_polygon(start = stline$endPoints[1, ],
						end = stline$endPoints[2, ],
						width_n = width_N,
						res_m = raster::xres(elev),
						Longitude_WGS84 = sp::coordinates(stline$endPoints_WGS84[1, ])[1])

	#Guarantee that band transect is contained by grid_gradient
	stband_ext <- raster::extent(stband_polygon$tpoly_orig)
	temp <- raster::extent(elevCropped)
	
	if (identical(raster::union(temp, stband_ext), temp)) {

		#Extract a maximal band transect along elevation transect
		stband_grids <- extract_tband_grids(tpoly_UTM = stband_polygon$tpoly_UTM,
											width_n = width_N,
											crs_UTM = stband_polygon$crs_UTM,
											declin = stband_polygon$declin,
											center = stband_polygon$center,
											res_rotate_grid = res_rotate_grid,
											grid1 = raster::crop(elev, stband_ext), grid2 = raster::crop(gap, stband_ext),
											fun1 = mean, fun2 = function(x, ...) majority(x, seed = seed))
	
		stband_grids$DistAlongXaxis_m <- raster::xres(stband_grids$grid1) * ((1:dim(stband_grids$grid1)[2]) - 1/2)
		#Note: identical(	[clump(extract_tband_grids(gap) %in% valList)],
		#					[extract_tband_grids(clump(gap %in% valList))])
		#may or may not be true depending on valList and spatial configuration,
		#e.g., true for i=16 & valList=c(489, 490, 491) and false for i=16 & valList=489
		stband_grids$bse <- raster::calc(stband_grids$grid2, fun = function(x) ifelse(x %in% bseRatValue, x, NA))
		stband_grids$bse_dens <- count_y_at_each_x(stband_grids$bse)
		stband_grids$tf <- raster::calc(stband_grids$grid2, fun = function(x) ifelse(x %in% tfRatValue, x, NA))
		stband_grids$tf_dens <- count_y_at_each_x(stband_grids$tf)
	} else {
		# goto next candidate transect
		stband_grids <- NULL
	}

	list(polygon = stband_polygon, grids = stband_grids, ext = stband_ext)
}



#' Crop band transect to limits1_ZoneEcolBoundary
#'
#' @export
trim1_band_to_ecotone <- function(stband_grids, start_icol, end_icol) {
	ext_to_crop <- raster::extent(stband_grids$grid1, r1 = 1, r2 = nrow(stband_grids$grid1), c1 = start_icol, c2 = end_icol)
	ids <- start_icol:end_icol
	
	res <- list()
	res$bse <- raster::crop(stband_grids$bse, ext_to_crop)
	res$tf <- raster::crop(stband_grids$tf, ext_to_crop)
	res$bse_dens <- stband_grids$bse_dens[ids]
	res$tf_dens <- stband_grids$tf_dens[ids]
	
	res
}

			
#' @export	
trim2_band_to_ecotone <- function(tempBand1, tempBand2, tempBand3, tempPatches, tempHuman, start1_icol, start2_icol, end2_icol){
	ext1 <- raster::extent(tempBand2$grid1, r1 = 1, r2 = nrow(tempBand2$grid1), c1 = start2_icol, c2 = end2_icol) # extent to crop to
	
	res <- list()
	#grids
	res$elev <- raster::crop(tempBand2$grid1, ext1)
	ext2 <- raster::extent(res$elev) # extent to shift to zero origin
	res$elev <- set_zero_origin(res$elev, shift_x = -ext2@xmin, shift_y = -ext2@ymin)
	res$bse <- set_zero_origin(raster::crop(tempBand2$bse, ext1), shift_x = -ext2@xmin, shift_y = -ext2@ymin)
	res$tf <- set_zero_origin(raster::crop(tempBand2$tf, ext1), shift_x = -ext2@xmin, shift_y = -ext2@ymin)
	res$human <- set_zero_origin(raster::crop(tempHuman, ext1), shift_x = -ext2@xmin, shift_y = -ext2@ymin)
	res$Veg1patch8 <- set_zero_origin(raster::crop(tempPatches$grid1, ext1), shift_x = -ext2@xmin, shift_y = -ext2@ymin)
	res$Veg2patch8 <- set_zero_origin(raster::crop(tempPatches$grid2, ext1), shift_x = -ext2@xmin, shift_y = -ext2@ymin)
	
	#location
	res$declin <- tempBand1$declin
	crs_orig <- raster::crs(tempBand1$pts_corner)
	crs_UTM <- raster::crs(tempBand2$grid1)
	corners_elidedUTM <- sp::SpatialPoints(coords = data.frame(x = c(ext1@xmin, ext1@xmin, ext1@xmax, ext1@xmax), y = c(ext1@ymin, ext1@ymax, ext1@ymax, ext1@ymin)), proj4string = crs_UTM)
	if (requireNamespace("maptools", quietly = TRUE)) {
		corners_UTM <- maptools::elide(obj = corners_elidedUTM, rotate = -tempBand1$declin * 180/pi, center = tempBand1$center)
	} else {
		temp <- rotate_coords(sp::coordinates(corners_elidedUTM), angle_rad = tempBand1$declin, center = tempBand1$center)
		corners_UTM <- sp::SpatialPoints(coords = data.frame(x = temp[, 1], y = temp[, 2]))
	}
	raster::crs(corners_UTM) <- crs_UTM
	res$band.pts <- sp::spTransform(corners_UTM, crs_orig)
	#Veg densities
	#res$Veg1dens <- tempBand3$Veg1dens[ids2 <- (start2_icol:end2_icol - start1_icol + 1)]
	#res$Veg2dens <- tempBand3$Veg2dens[ids2]
	
	res
}
