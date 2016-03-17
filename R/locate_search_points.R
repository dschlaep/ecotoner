# Sample random points where grid==1 to initiate gradient search
sample_init_cells <- function(N, grid, gridNA = NULL, pts_init = NULL, seed = NULL, verbose = FALSE, iter_max = 7, iter = 1) {
	#---Test input
	if (inherits(pts_init, "SpatialPoints")) {
		pts_init <- sp::coordinates(pts_init)
	}

	#---Begin function calculations
	names(grid) <- "stratum"
	tempfile_init <- file.path(tempdir(), "ecotoner_temp_init_points.RData")
	
	if (!is.na(seed)) set.seed(seed)

	pts_new <- try(raster::sampleStratified(grid, size = N, xy = TRUE))
		#Error in ys[[i]] <- y : attempt to select less than one element
		#	- raster(file.path(dir.bse, "bseEdge_CroppedWyoming.tif"))
		#	- raster(file.path(dir.bse, "gapv2_bse_iedge.tif"))
	if (inherits(class(pts_new), "try-error")) {
		if (!is.na(seed)) set.seed(seed)
		pts_new <- raster::sampleRandom(grid, 3 * N, na.rm = TRUE, xy = TRUE)
	}
	
	if (!is.null(gridNA)) {
		# Do not use cells where corresponding gridNA value is NA
		pts_new[, "stratum"] <- ifelse(is.na(raster::extract(gridNA, sp::SpatialPoints(pts_new[, c("x", "y")]))), -1, pts_new[, "stratum"])
	}
	pts_new <- (temp <- as.data.frame(pts_new))[temp[, "stratum"] == 1, c("x", "y")]
	pts_init <- rbind(pts_init, pts_new)
	
	if (length(temp <- which(duplicated(pts_init))) > 0) {
		pts_init <- pts_init[-temp, ] # remove duplicated grid cells
	}
	
	if (nrow(pts_init) < N & iter <= iter_max) {
		save(N, pts_init, iter, file = tempfile_init)
		if(verbose) cat("'ecotoner' starting points: ", iter, " iteration of produced ", nrow(pts_init), " points at ", format(Sys.time(), format = ""), "\n", sep = "")
		pts_init <- Recall(N, grid, inhibit_dist, gridNA, pts_init, seed, verbose, iter_max, iter + 1)
	} else {
		if (!is.na(seed)) set.seed(seed)
		pts_init <- pts_init[sample(x = nrow(pts_init), size = min(nrow(pts_init), N)), ]
		if(verbose) cat("'ecotoner' starting points: concluded with ", nrow(pts_init), " points at ", format(Sys.time(), format = ""), "\n", sep = "")
	}
	
	if (file.exists(tempfile_init)) unlink(tempfile_init)
	
	sp::SpatialPoints(coords = pts_init, proj4string = raster::crs(grid))
}


# Sample rSSI points either from mywin or where grid==1 to initiate gradient search
sample_inhibited_init_cells <- function(N, grid, inhibit_dist, mywin = NULL, gridNA = NULL, pts_init = NULL, initwindowfile = NA_character_, seed = NULL, verbose = FALSE) {
	#---Test input
	if (is.null(mywin)) {
		if (!inherits(grid, "RasterLayer")) {
			stop("ecotoner::sample_inhibited_init_cells(): either 'mywin' must be a matrix with a 'x' and a 'y' column or 'grid' must be of class 'Raster' to generate a matrix of suitable coordinates from which to randomly sample")
		} else {
			if(verbose) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "'sample_inhibited_init_cells' creating a matrix of coordinates", "\n")
			
			# Do not use cells where corresponding gridNA value is NA
			if (!is.null(gridNA) && inherits(gridNA, "RasterLayer")) {
				do_mask <- if (raster::canProcessInMemory(gridNA)) {
								if (raster::cellStats(gridNA, "countNA") > 0) TRUE else FALSE
							} else TRUE
				if (do_mask) {
					if (!raster::compareRaster(grid, gridNA, stopiffalse = FALSE)) {
						ext_intersect <- raster::intersect(raster::extent(grid), raster::extent(gridNA))
						grid <- raster::crop(grid, ext_intersect)
						gridNA <- raster::crop(gridNA, ext_intersect)
					}
					grid <- raster::mask(grid, gridNA)
				}
			}
			
			# Do not use cells that are 0
			if (raster::freq(grid, value = 0) > 0) grid <- raster::calc(grid, function(x) ifelse(abs(x) < sqrt(.Machine$double.eps), NA, x))
						
			# create matrix
			mywin <- raster::rasterToPoints(grid)[, c("x", "y")]
			
			if (!is.na(initwindowfile)) saveRDS(mywin, file = initwindowfile)
		}
	} else if (inherits(mywin, "matrix")) {
		mywin <- mywin[, c("x", "y")]
	} else if (inherits(mywin, "data.frame")) {
		mywin <- as.matrix(mywin[, c("x", "y")])
	} else {
		stop("ecotoner::sample_inhibited_init_cells(): either 'mywin' must be a matrix with a 'x' and a 'y' column or 'grid' must be of class 'Raster' to generate a matrix of suitable coordinates from which to randomly sample")
	}
	
	if (inherits(pts_init, "SpatialPoints"))  pts_init <- sp::coordinates(pts_init)

	if (!is.na(seed)) set.seed(seed)

	#---Begin function calculations
	x.init <- if (nrow(pts_init) > 0) pts_init[, c("x", "y")] else NULL
	
	if(verbose) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "'sample_inhibited_init_cells' started to generate a random point pattern as a realisation of a simple sequential inhibition process", "\n")
	# spatstat::rSSI() used to much memory due to excessive copying and the fact that ppp and owin classes are not well-suited for sparse objects; rSSI2 is a downscaled version that reduces copying and operates on coordinates instead of point-patterns
	# Note: the owin object was 16.8 GB in memory
	pts_new <- rSSI2(r = inhibit_dist, n = N, win = mywin, giveup = max(1000L, as.integer(0.1 * N)), x.init = x.init)
	pts_init <- rbind(pts_init, pts_new)
	
	sp::SpatialPoints(coords = pts_init, proj4string = raster::crs(grid))
}

#' Sample randomly \code{N} points to initiate the search for suitable gradients
#'
#' @param N An integer value. The number of search points to generate.
#' @param grid_mask1 A \code{\linkS4class{RasterLayer}} object from which cells with the value of 1 will be sampled.
## \link[package:MyClass-class]{MyClass}
#' @param mywindow A \code{\linkS4class{owin}} of the package 'spatstat' or NULL. If inhibit_dist is not 0 or NULL, then an inihibited random point process generates points within this window or, if NULL, then one is created based on grid_mask1 and grid_maskNA.
#' @param inhibit_dist A numeric value or \code{NULL}. The inhibition distance in meters between randomly sampled points.
#' @param grid_maskNA A \code{\linkS4class{RasterLayer}} object whose cells with the value of NA will prevent sampling of corresponding cells in \code{grid_mask1}.
#' @param initfile A filename (and its path) or \code{NA}. If not \code{NA}, then the generated search points are stored to disk in \link[=readRDS]{.rds} format.
#' @param initwindowfile A filename (and its path) or \code{NA}. If not \code{NA} and if a spatstat::owin object is created for an inhibited random sample, then the generated 'owin' object is stored to disk in \link[=readRDS]{.rds} format.
#' @param seed An integer value or \code{NULL} to set the random number generator. \code{NULL} follows the behavior of \code{\link{set.seed}}.
#' @param verbose A logical value. If \code{TRUE}, then progress statements are printed out.
#'
#' @return A \code{\link[sp]{SpatialPoints}} object with \code{N} unique points. If a suitable \code{initfile} is provided, then the returned object is also stored to disk. 
#'
#' @examples
#' path_demo <- system.file("extdata", package = "ecotoner")
#' abutt_eg <- raster::raster(file.path(path_demo, "abutt_eg.grd"))
#' elev_eg <- raster::raster(file.path(path_demo, "elev_eg.grd"))
#'
#' (temp <- get_transect_search_points(N = 10, grid_mask1 = abutt_eg, grid_maskNA = elev_eg))
#'
#' raster::plot(elev_eg)
#' raster::image(raster::calc(abutt_eg, fun = function(x) ifelse(!is.na(x) & x == 1, 1, NA)), col = "black", add = TRUE)
#' points(temp, col = "red", lwd = 2)
#'
#' @export
get_transect_search_points <- function(N, grid_mask1, mywindow = NULL, inhibit_dist = NULL, grid_maskNA = NULL, initfile = NA_character_, initwindowfile = NA_character_, seed = NULL, verbose = TRUE) {
	#---Test input
	if (!inherits(grid_mask1, "RasterLayer") || (!is.null(grid_maskNA) && !inherits(grid_maskNA, "RasterLayer"))) {
		stop("ecotoner::get_transect_search_points(): argument(s) 'grid_mask1' and/or 'grid_maskNA' are not of class 'Raster'")
	}
	if (!is.null(grid_maskNA) && !raster::compareRaster(grid_mask1, grid_maskNA, extent = FALSE, rowcol = FALSE, crs = TRUE, res = TRUE, orig = TRUE, rotation = TRUE)) {
		stop("ecotoner::get_transect_search_points(): arguments 'grid_mask1' and 'grid_maskNA' must have identical CRS, resolution, origin, and rotation")
	}
	if (!is.na(initfile) && !identical(tolower(tail(strsplit(basename(initfile), split = ".", fixed = TRUE)[[1]], n = 1)), "rds")) {
		stop("ecotoner::get_transect_search_points(): extension of argument 'initfile' must be '.rds'")
	}
	N <- as.integer(N)
	if (N < 1L) {
		stop("ecotoner::get_transect_search_points(): N must be > 0")
	}
	
	#---Begin function calculations
	if (!is.na(initfile) && file.exists(initfile)) {
		initpoints <- readRDS(file = initfile)
		ip_N <- length(initpoints)
	} else {
		initpoints <- array(NA, dim = c(0, 2), dimnames = list(NULL, c("x", "y")))
		ip_N <- 0L
	}
	
	if (!is.na(seed)) set_default_seed(seed)

	if (ip_N < N) {
		dont_inhibit <- is.null(inhibit_dist) || inhibit_dist <= 0
		
		if (!dont_inhibit) {
			temp <- try(sample_inhibited_init_cells(N = N, grid = grid_mask1, inhibit_dist = inhibit_dist, mywin = mywindow, gridNA = grid_maskNA, pts_init = initpoints, initwindowfile = initwindowfile, seed = NA, verbose = verbose), silent = TRUE)
			if (inherits(temp, "try-error")) {
				dont_inhibit <- TRUE
			} else {	
				initpoints <- temp
			}
		}
		
		if (dont_inhibit) {
			initpoints <- sample_init_cells(N = N, grid = grid_mask1, gridNA = grid_maskNA, pts_init = initpoints, seed = NA, verbose = verbose)
		}
		
		if (!is.na(initfile)) saveRDS(initpoints, file = initfile)
	} else if (ip_N > N) {
		initpoints <- initpoints[sample(x = length(initpoints), size = N), ]
	}
	
	initpoints
}
