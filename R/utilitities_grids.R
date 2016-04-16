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
calc_variogram_range <- function(data, maxX, res_m, na.rm = FALSE) {
	res <- NA
	
	if (requireNamespace("gstat", quietly = TRUE)) {
		if (inherits(data, "SpatialPixelsDataFrame")) {
			if (na.rm) data <- data[complete.cases(data@data), ]
		} else {
			if (na.rm) {
				data <- as(data, "SpatialPointsDataFrame") # NAs are excluded from coercion to "SpatialPointsDataFrame" (but not if coerced to "SpatialGridDataFrame")
				sp::gridded(data) <- TRUE #upgrade to SpatialPixelDataFrame
			} else {
				data <- as(data, "SpatialGridDataFrame") 
			}
		}

		if (!isTRUE(all.equal(data@grid@cellsize[1], data@grid@cellsize[2], tolerance = .Machine$double.eps, check.attributes = FALSE))) {
			warning("ecotoner:calc_variogram_range(): forcing cells to be square")
			data@grid@cellsize[1] <- data@grid@cellsize[2] <- res_m
		}
				
		# gstat::variogram will fail with error "dimensions do not match" if data contains NAs
		sample.var <- gstat::variogram(layer ~ 1, data = data) #~ 1 means no regression, i.e. normal semivariogram
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

grid_to_NA1 <- function(grid, vals, ...) {
	raster::calc(grid, fun = function(x) ifelse(is.na(x) | !(x %in% vals), NA, 1), ...)
}

gridNA1_to_01 <- function(grid, ...) {
	raster::calc(grid, fun = function(x) ifelse(!is.na(x) & sapply(x, function(x2) isTRUE(all.equal(x2, 1))), 1, 0), ...)
}


fun_abut <- function(e1, v1, e2, v2) ifelse((e1 == 1 & v2 == 1) | (v1 == 1 & e2 == 1), 1, NA)

calc_abutting <- function(Veg1, Veg2, filename, asNA = TRUE, ...) {
	dots_edges <- dots <- list(...)
	
	if (is.null(dots_edges$type)) dots_edges$type <- "outer"
	if (is.null(dots_edges$classes)) dots_edges$classes <- FALSE
	if (is.null(dots_edges$directions)) dots_edges$directions <- 4
	dots_edges$asNA <- asNA
	dots_edges$x <- Veg1
	edge1 <- do.call(raster::boundaries, args = dots_edges)
	dots_edges$x <- Veg2
	edge2 <- do.call(raster::boundaries, args = dots_edges)
	
	
	if (!missing(filename)) dots$filename <- filename
	args_overlay <- c(list(x = edge1, y = edge2, Veg1, Veg2, fun = fun_abut), dots)
	grid_res <- do.call(raster::overlay, args = args_overlay)
	
	if (!asNA) {
		args_to01 <- c(list(grid = grid_res), dots)
		if (!missing(filename)) {
			if (is.null(args_to01$overwrite)) {
				args_to01 <- c(args_to01, list(overwrite = TRUE))
			} else {
				args_to01$overwrite <- TRUE
			}
		}
		grid_res <- do.call(gridNA1_to_01, args = args_to01)
	}
	
	grid_res
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
	i_min <- which(path$elev == min(path$elev))[1] # the first point with minimal elevation
	i_max <- tail(which(path$elev == max(path$elev)), n = 1) # the last point with maximal elevation
	path[sort(i_min:i_max), ]

}

get_xyz <- function(spdf, field) {
	sp:::as.data.frame.SpatialPointsDataFrame(spdf)[c("x", "y", field)]
}

#' rasterToPoints
#'
#' Identical to \code{\link{raster::rasterToPoints}} if \code{not_convertible} is \code{NA}.
#'
#' rasterToPoints does not convert cells with values that are contained as single values in not_convertible or that are within and including the range of min-max of the length two vectors of not_convertible.
#'
#' @inheritParams raster::rasterToPoints
#' @param not_convertible A list of single elements, e.g., NA, 10, and of vectors of length two, e.g., c(0, 10).
#'
#' @seealso \code{\link{raster::rasterToPoints}} which this function simplifies
#' @examples
#' if (requireNamespace("raster")) {
#'   r <- raster::raster(nrow = 5, ncol = 5)
#'   r[] <- 1:(raster::ncell(r))
#'   r[1:5] <- NA
#'   
#'   identical(rasterToPoints(r, not_convertible = NA), raster::rasterToPoints(r))
#'   rasterToPoints(r, not_convertible = NULL)
#'   rasterToPoints(r, not_convertible = 25)
#'   rasterToPoints(r, not_convertible = list(NA, 25, c(7.1, 20.5)))
#' }
#' @export
rasterToPoints <- function(x, fun = NULL, spatial = FALSE, not_convertible = list(NA), ...) {
	if (identical(not_convertible, NA)) {
		res <- raster::rasterToPoints(x, fun = fun, spatial = spatial, ...)
	} else {
		keep_NA <- FALSE
		temp_grid <- if (!anyNA(not_convertible) && raster::cellStats(x, "countNA") > 0) {
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
				if (any(iremove)) res <- res[!iremove, , drop = FALSE]
			}
		}

		if (keep_NA) res[res[, 'layer'] == val_for_NA, 'layer'] <- NA
		
	}
	
	res
}

transect_to_long <- function(x, y, na.rm = TRUE) {
	if (inherits(x, "RasterLayer")) x <- rasterToPoints(x, not_convertible = NULL)
	if (inherits(y, "RasterLayer")) y <- rasterToPoints(y, not_convertible = NULL)
	
	if (inherits(x, "matrix") && inherits(y, "matrix") && nrow(x) == nrow(y) && ncol(x) >= 3 && ncol(y) >= 3) {
		isnotna <- if (na.rm) !(is.na(x[, 3]) | is.na(y[, 3])) else rep(TRUE, nrow(x))
		res <- matrix(NA, nrow = sum(isnotna), ncol = 4, dimnames = list(NULL, c("x", "y", "cols", "rows")))
		res[, "x"] <- x[isnotna, 3] # values of first raster
		res[, "y"] <- y[isnotna, 3] # values of second raster
		res[, "cols"] <- x[isnotna, 1] # coordinates along length of transect
		res[, "rows"] <- x[isnotna, 2] # coordinates along width of transect
	} else {
		res <- matrix(NA, nrow = 0, ncol = 4, dimnames = list(NULL, c("x", "y", "cols", "rows")))
	}
	
	res
}


#' An implementation of a Simple Sequential Inhibition (SSI) process
#'
#' rSSI of cell-center points and avoiding spatstat. This function generates a random pattern of points selected from \code{win} all of which have at least \code{r} units distance to each other.
#'
#' This function implements the same process as \code{\link{spatstat::rSSI}}. \code{\link{spatstat::rSSI}} relies on \code{win} being of class \code{owin}, produces an object of class \code{ppp} - both very memory-intensive classes -, and copies such objects excessively. For instance, a matrix with 75 million rows (1.1 GB memory) was translated into an object of class \code{owin} of 16.8 GB memory usage. Running \code{\link{spatstat::rSSI}} on this 'owin' object failed after the process required more than than 84 GB memory - more than was available.
#' 
#' @inheritParams spatstat::rSSI
#' @seealso \code{\link{spatstat::rSSI}} which this function simplifies
#' @export
rSSI2 <- function(r, n, win, giveup = 1000, x.init = NULL, ...) {
    X <- if (inherits(x.init, "matrix") && ncol(x.init) == 2) {
				x.init[, c("x", "y")]
			} else {
				matrix(NA, nrow = 0, ncol = 2, dimnames = list(NULL, c("x", "y")))
			}
	
    r2 <- r^2
    N <- nrow(win)
    
    ntries <- 0    
    while (ntries < giveup && nrow(X) < n) {
        qq <- win[sample.int(N, size = 1L), c("x", "y")]

        if (all((qq["x"] - X[, "x"])^2 + (qq["y"] - X[, "y"])^2 > r2)) {
        	X <- rbind(X, qq)
        	ntries <- 0
        } else ntries <- ntries + 1
    }
    
	X
}


get_transect_grids_as_df <- function(i, et_desc = NULL, ecotoner_settings, migtypes, veg_types) {
	iflag <- flag_itransect(ecotoner_settings, i)
	
	do_measure <- TRUE
	if (file.exists(fname_etlocated(ecotoner_settings, iflag))) {
		load(fname_etlocated(ecotoner_settings, iflag)) #i, b, and etransect loaded
	} else {
		do_measure <- FALSE # no suitable transect located for search point i
	}
	
	if (do_measure && exists("etransect") && !is.null(etransect)) {
		res <- list()
		if (!is.null(et_desc)) idesc_remove <- which(colnames(et_desc) %in% c("Transect_ID", "Neighbor_Cells"))
		
		for (b in seq_len(neighbors_N(ecotoner_settings))) {
			for (im in seq_along(migtypes)) {
				for (iveg in seq_along(veg_types)) {
					# Interpret NAs as absences of iveg
					y <- raster::calc(etransect[["etbands"]][[b]]$Veg[[migtypes[im]]][[veg_types[iveg]]]$grid, function(x) ifelse(is.na(x), 0, x))
					temp <- transect_to_long(x = etransect[["etbands"]][[b]]$Env$elev$grid, y = y)
				
					if (im == 1 && iveg == 1) {
						mat <- matrix(NA, nrow = nrow(temp), ncol = 3 + length(migtypes) * length(veg_types))
						mat[, 1:4] <- temp[, c("x", "cols", "rows", "y")]
					} else {
						mat[, 3 + (im - 1) * length(veg_types) + iveg] <- temp[, "y"]
					}
				}
			}
			
			res[[b]] <- if (is.null(et_desc)) {
							data.frame(Transect_ID = i,
										Neighbor_Cells = neighborhoods(ecotoner_settings)[b],
										mat,
										row.names = NULL)
						} else {
							data.frame(Transect_ID = i,
										Neighbor_Cells = neighborhoods(ecotoner_settings)[b],
										mat,
										et_desc[et_desc[, "Transect_ID"] == i & et_desc[, "Neighbor_Cells"] == neighborhoods(ecotoner_settings)[b], -idesc_remove],
										row.names = NULL)
						}
		}
		
		do.call(rbind, res)
	} else NULL
}
	

# modified from raster:::.circular.weight
#' @param rs A vector of length two. The resolution of a raster object.
#' @param nbs A numerical value. The distance of the circular weights matrix in raster units.
focalWeight_inverse <- function (rs, nbs) {
	nx <- 1 + 2 * floor(nbs / rs[1])
    ny <- 1 + 2 * floor(nbs / rs[2])
    m <- matrix(ncol = nx, nrow = ny)
	m[ceiling(ny/2), ceiling(nx/2)] <- 1
    if (nx == 1 && ny == 1) return(m)

	d <- raster::raster(m, xmn = 0, xmx = nx * rs[1], ymn = 0, ymx = ny * rs[2], crs = "+proj=utm +zone=1 +datum=WGS84")
	d <- raster::distance(d)
	d <- raster::calc(d, fun = function(x) ifelse(is.na(x) | x > nbs | abs(x) < sqrt(.Machine$double.eps), 0, 1 / x))
	m <- raster::as.matrix(d)
    m / sum(m)
}
