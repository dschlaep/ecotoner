#' @export
crop_to_neighborhood <- function(pt_start, neighbor_n, grid_gradient, grid_veg = NULL, grid_veg_abut12 = NULL, grid_aspect_mean = NULL, grid_aspect_sd = NULL, transect_type) {
	icelle <- raster::cellFromXY(grid_gradient, pt_start)
	if (raster::validCell(grid_gradient, icelle)) {
		irowcole <- raster::rowColFromCell(grid_gradient, icelle)
		n_extend <- round((neighbor_n - 1) / 2)
		dims <- dim(grid_gradient)
		r1 <- min(max(1, irowcole[1] - n_extend), dims[1] - 1)
		r2 <- r1 + min(max(1, 2 * n_extend), dims[1] - r1)
		c1 <- min(max(1, irowcole[2] - n_extend), dims[2] - 1)
		c2 <- c1 + min(max(1, 2 * n_extend), dims[2] - c1)
		ext_cropped <- 1.05 * raster::extent(grid_gradient, r1, r2, c1, c2)
	
		grid_elev_cropped <- raster::crop(grid_gradient, ext_cropped)
		grid_veg_cropped <- raster::crop(grid_veg, ext_cropped)
		if (transect_type == 4) {
			grid_veg_abut12_cropped <- raster::crop(grid_veg_abut12, ext_cropped)
			grid_aspect_mean_cropped <- raster::crop(grid_aspect_mean, ext_cropped)
			grid_aspect_sd_cropped <- raster::crop(grid_aspect_sd, ext_cropped)
		} else {
			grid_veg_abut12_cropped <- NULL
			grid_aspect_mean_cropped <- NULL
			grid_aspect_sd_cropped <- NULL
		}
	} else {
		# goto next start point
		grid_elev_cropped <- NULL
		grid_veg_cropped <- NULL
		grid_veg_abut12_cropped <- NULL
		grid_aspect_mean_cropped <- NULL
		grid_aspect_sd_cropped <- NULL
	}
	
	list(grid_elev_cropped = grid_elev_cropped,
		grid_veg_cropped = grid_veg_cropped,
		grid_veg_abut12_cropped = grid_veg_abut12_cropped,
		grid_aspect_mean_cropped = grid_aspect_mean_cropped,
		grid_aspect_sd_cropped = grid_aspect_sd_cropped)
}
					


#' Determine the range of spatial autocorrelation
calc_variogram_range <- function(grid, maxX) {
	res <- NA
	
	if (requireNamespace("gstat", quietly = TRUE)) {
		vetemp <- as(grid, "SpatialGridDataFrame")
		if (!isTRUE(all.equal(vetemp@grid@cellsize[1], vetemp@grid@cellsize[2], tolerance=.Machine$double.eps))) {
			warning("'ecotoner': forcing cells of 'SpatialGridDataFrame' to be square")
			vetemp@grid@cellsize[1] <- vetemp@grid@cellsize[2] <- raster::xres(grid)
		}
	
		sample.var <- gstat::variogram(layer ~ 1, data = vetemp) #~ 1 means no regression, i.e. normal semivariogram
		fvar <- gstat::fit.variogram(sample.var, model = gstat::vgm(1, "Sph", 3000, 1000), debug.level = 0)
		if (fvar$range[2] > maxX) {
			fvar <- gstat::fit.variogram(sample.var, model = gstat::vgm(1, "Gau", 3000, 1000), debug.level = 0)
		}
		
		res <- fvar$range[2]
	} else {
		warning("Package 'gstat' not installed: the range of the variogram will not be determined by function 'calc_variogram_range'")
	}
	
	res
}


classify_grid <- function(grid1, grid2, grid3, k, seed = NULL) {
	#Adopted 2013/08/23: from 'Terrain Classification Experiment 2: GRASS, R, and the raster package': http://casoilresource.lawr.ucdavis.edu/drupal/node/993
	res <- NA
	
	if (requireNamespace("cluster", quietly = TRUE) && requireNamespace("randomForest", quietly = TRUE)) {	
		if (!is.na(seed)) set.seed(seed)
		s <- raster::stack(grid1, grid2, grid3)
		names(s) <- c('grid1', 'grid2', 'grid3')
		sr <- na.omit(as.data.frame(raster::sampleRegular(s, min(1e4, 0.1*ncell(elev))))) #randomForest doesn't like NAs
		sclara <- cluster::clara(sr, stand = TRUE, k = k, samples = 50, rngR = TRUE, pamLike = TRUE, medoids.x = FALSE, keep.data = FALSE)
		sr$cluster <- factor(sclara$clustering)
		rf <- randomForest::randomForest(cluster ~ grid1 + grid2 + grid3, data = sr, importance = TRUE, ntree = 201)
	
		res <- predict(s, rf, type = 'response')
	} else {
		warning("Package(s) 'cluster' and/or 'randomForest' not installed: the classification will not be determined by function 'classify_grid'")
	}
	
	res
}

#' export
count_y_at_each_x <- function(grid) {
	#maybe better use instead: sample density: y1means <- sapply(1:(N-1), FUN=function(n) sapply(1:(N-n), FUN=function(i) mean(y1[i:(i+n)]))
	
	apply(mat <- raster::as.matrix(grid), 2, FUN = function(x) sum(!is.na(x))) / dim(mat)[1]
}

mean_y_at_each_x <- function(grid) {
	apply(raster::as.matrix(grid), 2, mean, na.rm = TRUE)
}

sd_y_at_each_x <- function(grid) {
	apply(raster::as.matrix(grid), 2, sd, na.rm = TRUE)
}

calc_abutting <- function(Veg1, Veg2) {
	edge1 <- raster::boundaries(Veg1, type = "outer", classes = FALSE, directions = 4)
	edge2 <- raster::boundaries(Veg2, type = "outer", classes = FALSE, directions = 4)
	
	raster::overlay(edge1, edge2, Veg1, Veg2, fun = function(e1, v1, e2, v2) ifelse((e1 == 1 & v2 == 1) | (v1 == 1 & e2 == 1), 1, NA))	
}

set_zero_origin <- function(grid, shift_x = 0, shift_y = 0) {
	raster::shift(grid, x = shift_x, y = shift_y)
}

spdf_with_extract <- function(points, grid) {
	temp_prj <- if (inherits(raster::crs(points), "CRS")) raster::crs(points) else raster::crs(grid)
	sp::SpatialPointsDataFrame(coords = sp::coordinates(points),
								data = data.frame(elev = raster::extract(grid, points)),
								proj4string = temp_prj)
}


sort_path <- function(start, path, longlat) {
	#sort points along path by increasing distance along transect
	path_dist_m <- sp::spDistsN1(path, start, longlat = longlat) * if (longlat) 1000 else 1
	isort <- order(path_dist_m)
	
	list(path = path[isort, ], dist_m = path_dist_m[isort])
}


sort_path_max_range <- function(start, path, grid_gradient) {
	path <- spdf_with_extract(points = path, grid = grid_gradient)
	path <- path[!is.na(path$elev),]
	
	path <- sort_path(start, path, longlat = raster::isLonLat(grid_gradient))$path
	
	#cut out section between highest and lowest elevation
	path[sort(which.min(path$elev):which.max(path$elev)), ]

}

get_xyz <- function(spdf, field) {
	sp:::as.data.frame.SpatialPointsDataFrame(spdf)[c("x", "y", field)]
}

#' rasterToPoints
#'
#' Identical to \code{\link{raster::rasterToPoints}} if \code{not_convertible} is \code{NULL}.
#'
#' rasterToPoints does not convert cells with values that are contained as single values in not_convertible or that are within and including the range of min-max of the length two vectors of not_convertible.
#'
#' @inheritParams raster::rasterToPoints
#' @param not_convertible A list of single elements, e.g., NA, 10, and of vectors of length two, e.g., c(0, 10).
#'
#' @examples
#' if (requireNamespace("raster")) {
#'   r <- raster::raster(nrow = 5, ncol = 5)
#'   r[] <- 1:(raster::ncell(r))
#'   r[1:5] <- NA
#'   
#'   identical(rasterToPoints(r, not_convertible = NULL), raster::rasterToPoints(r))
#'   rasterToPoints(r, not_convertible = NA)
#'   rasterToPoints(r, not_convertible = list(NA, 25, c(7.1, 20.5)))
#' }
#' @export
rasterToPoints <- function(x, fun = NULL, spatial = FALSE, not_convertible = list(NA), ...) {
	if (!is.null(not_convertible)) {
		keep_NA <- FALSE
		temp_grid <- if (anyNA(not_convertible) && raster::cellStats(x, "countNA") > 0) {
							# we want cells with NA values to be converted to points
							# find a placeholder value
							x_vals <- sort(raster::unique(x))
							val_for_NA <- min(x_vals, na.rm = TRUE) - 1
							keep_NA <- TRUE
							raster::calc(x, function(x) ifelse(is.na(x), val_for_NA, x))
						} else x
		
		res <- raster::rasterToPoints(temp_grid, fun = fun, spatial = spatial, ...)
		
		for (i in seq_along(not_convertible)) {
			if (nrow(res) == 0) break
			if (!is.na(not_convertible[i])) {
				iremove <- FALSE
				if (length(not_convertible[i][[1]]) == 1) {
					iremove <- res[, 'layer'] == not_convertible[i][[1]]
				}
				if (length(not_convertible[i][[1]]) == 2) {
					iremove <- res[, 'layer'] >= not_convertible[i][[1]][1] & res[, 'layer'] <= not_convertible[i][[1]][2]
				}
				if (sum(iremove) > 0) res <- res[!iremove, , drop = FALSE]
			}
		}

		if (keep_NA) res[res[, 'layer'] == val_for_NA, 'layer'] <- NA
		
	} else {
		res <- raster::rasterToPoints(x, fun = fun, spatial = spatial, ...)
	}
	
	res
}

transect_to_long <- function(x, y, na.rm = TRUE) {
	if (inherits(x, "RasterLayer")) x <- rasterToPoints(x, not_convertible = NA)
	if (inherits(y, "RasterLayer")) y <- rasterToPoints(y, not_convertible = NA)
	
	if (inherits(x, "matrix") && inherits(y, "matrix") && nrow(x) == nrow(y) && ncol(x) >= 3 && ncol(y) >= 3) {
		if (na.rm) isnotna <- !(is.na(x[, 3]) | is.na(y[, 3])) else rep(TRUE, nrow(x))
		vy <- y[isnotna, 3]
		vx <- x[isnotna, 3]
		rows <- x[isnotna, 2]
	} else {
		vx <- vy <- rows <- NA
	}
	
	list(x = vx, y = vy, reps = rows)
}

