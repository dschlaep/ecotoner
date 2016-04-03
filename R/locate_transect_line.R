
#Determine location of highest value within neighborhood
locate_max_value <- function(x, point, neighbors) {
	icelle <- raster::cellFromXY(x, point)
	if (raster::validCell(x, icelle)) {
		irowcole <- raster::rowColFromCell(x, icelle)
		n_neigh <- round((neighbors - 1) / 2)
		dims <- dim(x)
		trow <- min(max(1, irowcole[1] - n_neigh), dims[1] - 1)
		tcol <- min(max(1, irowcole[2] - n_neigh), dims[2] - 1)
		foc <- raster::getValuesBlock(x, row=trow, nrows=min(max(1, neighbors), dims[1] - trow), col=tcol, ncols=min(max(1, neighbors), dims[2] - tcol), format = 'matrix')
		maxe <- which(foc == max(foc, na.rm = TRUE), arr.ind = TRUE) #row/col location in local block
		imax <- sp::SpatialPoints(data.frame(x = xFromCol(x, tcol + maxe[2] - 1), y=yFromRow(x, trow + maxe[1] - 1)), proj4string = raster::crs(x))
	} else {
		imax <- NULL
	}
	
	imax
}


#Determine location with steepest slope within neighborhood
locate_steepest_slope <- function(x, point, neighbors_max, neighbors_min) {
	icelle <- raster::cellFromXY(x, point)
	if (raster::validCell(x, icelle)) {
		#Crop to neighbors
		irowcole <- raster::rowColFromCell(x, icelle)
		n_neigh <- round((neighbors_max - 1) / 2)
		dims <- dim(x)
		r1 <- min(max(1, irowcole[1] - n_neigh), dims[1] - 1)
		r2 <- r1 + min(max(1, neighbors_max), dims[1] - r1)
		c1 <- min(max(1, irowcole[2] - n_neigh), dims[2] - 1)
		c2 <- c1 + min(max(1, neighbors_max), dims[2] - c1)
		cexp <- 1.01 * raster::extent(x, r1, r2, c1, c2)
		xcrop <- raster::crop(x, cexp)
		
		#Calculate slope from point
		dist_fromPoint <- raster::distanceFromPoints(xcrop, point)
		point_elev <- raster::extract(xcrop, point)
		neigh_maxDist <- neighbors_max * raster::xres(x)
		neigh_minDist <- neighbors_min * raster::xres(x)
		slope_fromPoint <- raster::overlay(xcrop, dist_fromPoint, fun=function(e, d) ifelse(d >= neigh_minDist & d <= neigh_maxDist, abs(e - point_elev) / d, 0))

		#Locate point with largest slope(s)
		maxSlope <- raster::maxValue(slope_fromPoint)
		if(!is.numeric(maxSlope)){
			slope_fromPoint <- raster::setMinMax(slope_fromPoint)
			maxSlope <- raster::maxValue(slope_fromPoint)
		}
		maxs <- raster::calc(slope_fromPoint, fun=function(x) ifelse(!is.na(x) & x >= maxSlope - 10*sqrt(.Machine$double.neg.eps), 1, NA)) #changed x from 1 to 10 in 'maxSlope - x*sqrt(.Machine$double.neg.eps)'; for some points x=1 didn't work (why?)
		imax <- sp::SpatialPoints(data.frame(x=(temp <- raster::rasterToPoints(x=maxs, spatial=FALSE)[, 1:2, drop=FALSE])[, 1], y=temp[, 2]), proj4string=raster::crs(x)) #this function doesn't work for spatial=TRUE if cellStats(maxs, "sum") == 1
		if (length(imax) > 1) { #In case, there are several, choose the one closest to point
			mdist <- raster::extract(dist_fromPoint, imax)
			imax <- imax[which.min(mdist)]
		}
	} else {
		imax <- NULL
	}
	
	imax
}



#Determines transect between point and highLoc
elongate_linear_transect <- function(x, point1, point2, efac = 1, extend.dir12 = FALSE, extend.dir21 = TRUE) {
	efac <- max(1, round(efac))
	point1xy <- sp::coordinates(point1)
	point2xy <- sp::coordinates(point2)
	if (extend.dir21) {
		point3 <- sp::coordinates(data.frame(point1xy[1] + efac * (point1xy[1] - point2xy[1]), point1xy[2] + efac * (point1xy[2] - point2xy[2])))	#Get point at distance efac * (point2 - point1) in opposite direction from point1
	} else {
		point3 <- point1
	}
	if (extend.dir12) {
		point4 <- sp::coordinates(data.frame(point2xy[1] + efac * (point2xy[1] - point1xy[1]), point2xy[2] + efac * (point2xy[2] - point1xy[2])))	#Get point at distance efac * (point1 - point2) in opposite direction from point2
	} else {
		point4 <- point1
	}
	temp_crs <- raster::crs(point2)
	tpoints <- sp::rbind.SpatialPoints(sp::SpatialPoints(point3, proj4string=temp_crs),
										sp::SpatialPoints(point1, proj4string=temp_crs),
										point2,
										sp::SpatialPoints(point4, proj4string=temp_crs))
	tline <- sp::SpatialLines(list(sp::Lines(sp::Line(tpoints), ID=1)), proj4string=temp_crs)
	
	temp <- raster::rasterize(x = tline, y = x)
	path <- raster::rasterToPoints(x = temp, spatial = TRUE)
	
	sort_path_max_range(start = point3, path = path, grid_gradient = elev)
}

# Extract from the longest segment that lowest to hightest value piece of the transect, identified by a homogeneous aspect
subset_candidate_THA <- function(grid_gradient, init_spLine, init_rle){
	temp <- raster::rasterize(x = init_spLine, y = grid_gradient)
	path <- raster::rasterToPoints(x = temp, spatial = TRUE)
	
	ids <- cumsum(c(0, init_rle$lengths))
	id <- init_rle$isegments[which.max(init_rle$lengths[init_rle$isegments])]
	
	path <- path[(ids[id]+1):ids[id+1], ] #Cut to run of focal mean aspect of init point	
	
	sort_path_max_range(start = sp::coordinates(path[1, ]), path = path, grid_gradient = grid_gradient)
}


#' Locate candidate transect lines of (sufficient) homogeneous aspect THAs
#' 
#' @param n An integer. The size of initial random sample of starting points to locate candidate transect lines
#'						within areas of \code{abut} and where \code{asp201SD} <= \code{max_aspect_sd}.
#' @param grid_gradient A \code{\linkS4class{RasterLayer}} object representing the environmental gradient, e.g., elevation.
#' @param asp201Mean A \code{\linkS4class{RasterLayer}} object.
#' @param asp201SD A \code{\linkS4class{RasterLayer}} object.
#' @param abuts A \code{\linkS4class{RasterLayer}} object.
#' @param max_aspect_sd A numeric. Limits the search area for line-starting points to where \code{asp201SD} <= \code{max_aspect_sd}.
#' @param max_aspect_diff A numeric. Cuts the candidate lines to areas where \code{asp201Mean} differs from the mean aspect at the starting point differs by less than \code{max_aspect_diff}.
#' @param width_N An integer. The width of the band transect in number of grid cells which will be put around the transect line.
#' @param seed An integer value, \code{NA}, or \code{NULL} to set the random number generator. \code{NULL} follows the behavior of \code{\link{set.seed}}.
locate_candidate_THAs <- function(n, grid_gradient, max_neighborhood, asp201Mean, asp201SD, abuts, max_aspect_sd, max_aspect_diff, width_N, seed = NULL) {
							
	candidates <- list()

	#1. Calculate aspect: focal mean and SD of band transect width (calls to focal() are VERY slow for large w as here)
	asp201Mean_atlowSD <- raster::mask(asp201Mean, raster::calc(asp201SD, fun = function(x) ifelse(x <= max_aspect_sd, x, NA)))

	#2. Find potential init points among abutting locations in areas of low focal aspect RMSE, take a random sample of 20
	lowSDabuts <- raster::mask(asp201Mean_atlowSD, abuts)
	if ((raster::ncell(lowSDabuts) - raster::cellStats(lowSDabuts, 'countNA')) >= n) {
		lowSDabuts_aspm <- raster::rasterToPoints(lowSDabuts, spatial = TRUE)
		
		if (!is.na(seed)) set.seed(seed)
		lowSDabuts_aspm_spPoints <- lowSDabuts_aspm[sample(x = nrow(lowSDabuts_aspm), size = n, replace = FALSE), ]

		#3. Create lines for each potential init point in the direction of its focal mean aspect
		extend_m <- sqrt(2) / 2 * raster::xres(grid_gradient) * (max_neighborhood - 1)
		pcoords <- sp::coordinates(lowSDabuts_aspm_spPoints)
		extXY_m <- cbind(tempX <- ifelse(abs(tempA <- tan(pi / 2 - lowSDabuts_aspm_spPoints@data[, 1])) > 1, extend_m / tempA, extend_m), tempX * tempA)
		points_pos <- sp::SpatialPoints(sp::coordinates(temp <- data.frame(pcoords[, 1] + extXY_m[, 1], pcoords[, 2] + extXY_m[, 2])), proj4string = raster::crs(lowSDabuts))
		points_neg <- sp::SpatialPoints(sp::coordinates(data.frame(pcoords[, 1] - extXY_m[, 1], pcoords[, 2] - extXY_m[, 2])), proj4string = raster::crs(lowSDabuts))
		seq_n <- seq_len(n)
		extLinesAspm_spLine <- sp::SpatialLines(lapply(seq_n, FUN = function(ip) sp::Lines(sp::Line(sp::rbind.SpatialPoints(points_pos[ip, ], points_neg[ip, ])), ID = ip)), proj4string = raster::crs(lowSDabuts))

		#4. Extract focal mean aspect within low RMSE areas along those lines (slow call)
		extLinesAspm_aspm <- raster::extract(asp201Mean_atlowSD, extLinesAspm_spLine, method = "simple") # NA, if asp201SD > max_aspect_sd

		#5. Retain those lines that have a sufficient length of homogeneous aspect (i.e., !is.na) around the init points
##TODO(drs): should this be difference from mean aspect of init point (as implemented currently) or of mean of all points?
		runs <- lapply(seq_n, function(il) rle(ifelse(abs(extLinesAspm_aspm[[il]] - lowSDabuts_aspm_spPoints@data[il, 1]) < max_aspect_diff, 1, 0)))
		runs <- lapply(runs, function(r) {r[["isegments"]] <- which(!is.na(r$values) & r$values == 1); r})
		good_lines <- which(sapply(runs, function(r) length(r$isegments) > 0 && max(r$lengths[r$isegments]) > width_N))

		if (length(good_lines) > 0) for(i in seq_along(good_lines)) {
			candidates[[i]] <- subset_candidate_THA(grid_gradient = grid_gradient,
													init_spLine = extLinesAspm_spLine[good_lines[i], ],
													init_rle = runs[[good_lines[i]]])
		}
	}
		
	candidates
}


#' @export
locate_candidate_transect_lines <- function(transect_type, start, neighbors, max_neighborhood, grid_gradient, grid_flow = NULL, asp201MeanCropped = NULL, asp201SDCropped = NULL, bseABUTtfCropped = NULL, candidate_THA_N = 1, width_N, max_aspect_sd, max_aspect_diff, seed = NULL) {
	#---Test input
	if (!inherits(grid_gradient, "RasterLayer")) {
		stop("ecotoner::locate_candidate_transect_lines(): argument 'grid_gradient' is not of class 'Raster'")
	}
	if (transect_type == 3 && !inherits(grid_flow, "RasterLayer")) {
		stop("ecotoner::locate_candidate_transect_lines(): argument 'grid_flow' is not of class 'Raster'")
	}
	if (transect_type == 4 && !inherits(asp201MeanCropped, "RasterLayer") && !inherits(asp201SDCropped, "RasterLayer") && !inherits(bseABUTtfCropped, "RasterLayer")) {
		stop("ecotoner::locate_candidate_transect_lines(): argument(s) 'asp201MeanCropped', 'asp201SDCropped', and/or 'bseABUTtfCropped' are not of class 'Raster'")
	}

	#---Begin function calculations
	candidates <- list()
	
	#2a. Elevation transect
	if (transect_type %in% 1:3) {
		#Get a second point within neighborhood
		if (transect_type == 1 || transect_type == 3) { 
			# Get the highest point within neighborhood
			highLoc <- locate_max_value(x = grid_gradient,
										point = start,
										neighbors = neighbors)
		} else if(transect_type == 2) {
			highLoc <- locate_steepest_slope(x = grid_gradient,
											point = start,
											neighbors_max = neighbors,
											neighbors_min = width_N)
		}
		
		if (!is.null(highLoc)) {
			if (transect_type == 1) {
				# Approach: transect = direct route between highLoc and ipoint and extend straight by factor x, cut out section between highest and lowest elevation
				candidates[[1]] <- elongate_linear_transect(x = grid_gradient,
														point1 = start,
														point2 = highLoc,
														efac = 2,
														extend.dir12 = FALSE,
														extend.dir21 = TRUE)
			} else if (transect_type == 2) {
				# Approach: transect = direct route between highLoc and ipoint and extend straight in both directions by factor x, cut out section between highest and lowest elevation
				candidates[[1]] <- elongate_linear_transect(x = grid_gradient,
														point1 = start,
														point2 = highLoc,
														efac = 2,
														extend.dir12 = TRUE,
														extend.dir21 = TRUE)
			} else if (transect_type == 3) {
				# Approach: transect = flowpath, but calc_flow_path gets stuck easily on flat areas or in small bowls
				temp <- calc_flow_path(flowdir = grid_flow, point = highLoc)
				candidates[[1]] <- spdf_with_extract(points = temp, grid = grid_gradient)
			}
		}
	} else if (transect_type == 4) {
		# Approach: transect: sufficient length along homogeneous aspect
		candidates <- locate_candidate_THAs(n = candidate_THA_N,
									grid_gradient = grid_gradient,
									max_neighborhood = max_neighborhood,
									asp201Mean = asp201MeanCropped,
									asp201SD = asp201SDCropped,
									abuts = bseABUTtfCropped,
									max_aspect_sd = max_aspect_sd,
									max_aspect_diff = max_aspect_diff,
									width_N = width_N,
									seed = seed)
	}
	
	list(lines = candidates, line_N = length(candidates))
}
	

#' @export
orient_transect_line <- function(pts_tcand, longlat) {
	#---Test input
	if (!inherits(raster::crs(pts_tcand), "CRS")) {
		stop("ecotoner::orient_transect_line(): argument 'pts_tcand' lacks a 'CRS'")
	}

	#---Begin function calculations
	stline <- list()
	
	pts_N <- nrow(pts_tcand)
	if (pts_N > 1) { # At least two points required to define start and end of a line
		temp <- pts_tcand[c(1, pts_N), ]
		
		if (!isTRUE(all.equal(min(temp$elev), max(temp$elev)))) { # We need a gradient > 0
			#Test elevation transect and orient it from low (row 1) -> high (row 2)
			if (which.max(temp$elev) == 1 && which.min(temp$elev) == 2) {
				# need to flip transect so that start=first-index is at lowest elevation
				pts_tcand <- pts_tcand[pts_N:1, ]
			}

			#Get transect start and end point (end point with highest elev)
			stline$endPoints <- pts_tcand[c(1, pts_N), ]
			stline$endPoints_WGS84 <- sp::spTransform(stline$endPoints, sp::CRS("+proj=longlat +datum=WGS84"))

			#Get transect distances and sort transect line
			temp <- sort_path(start = stline$endPoints[1, ], path = pts_tcand, longlat)
			stline$pts <- temp$path
			stline$dist_m <- temp$dist_m
			stline$pts_N <- pts_N
		}
	}
	
	stline
}
	

#' @export
trim_line_to_ecotone <- function(stline_pts, distX, start_m, end_m, crs, longlat){
	startID <- findInterval(start_m, c(0, distX), all.inside = TRUE)
	endID <- findInterval(end_m, c(0, distX), all.inside = TRUE)
	ids <- startID:endID
	
	res <- list()
	res$pts <- stline_pts[ids,]
	res$pts_N <- nrow(res$pts)
	res$endPoints <- res$pts[c(1, res$pts_N), ]
	raster::crs(res$endPoints) <- crs
	res$endPoints_WGS84 <- sp::spTransform(res$endPoints, sp::CRS("+proj=longlat +datum=WGS84"))
	res$dist_m <- sp::spDistsN1(res$pts, res$endPoints[1, ], longlat = longlat) * if(longlat) 1000 else 1
	res$StartInElevationTransect_m <- start_m
		
	res
}
				
