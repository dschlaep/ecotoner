# Sample random point where grid==1 to initiate gradient search
sample_init_cells <- function(N, grid, gridNA = NULL, pts_init = NULL, seed = NULL, verbose = FALSE, iter = 1) {
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
	
	if (nrow(pts_init) < N & iter <= 7) {
		save(N, pts_init, iter, file = tempfile_init)
		if(verbose) cat("'ecotoner' starting points: ", iter, " iteration of produced ", nrow(pts_init), " points at ", format(Sys.time(), format = ""), "\n", sep = "")
		pts_init <- Recall(N, grid, gridNA, pts_init, seed, verbose, iter + 1)
	} else {
		if (!is.na(seed)) set.seed(seed)
		pts_init <- pts_init[sample(x = nrow(pts_init), size = min(nrow(pts_init), N)), ]
		if(verbose) cat("'ecotoner' starting points: concluded with ", nrow(pts_init), " points at ", format(Sys.time(), format = ""), "\n", sep = "")
	}
	
	if (file.exists(tempfile_init)) unlink(tempfile_init)
	
	sp::SpatialPoints(coords = pts_init, proj4string = raster::crs(grid))
}


#' Sample randomly \code{N} points to initiate the search for suitable gradients
#'
#' @param N An integer value. The number of search points to generate.
#' @param grid_mask1 A \code{\linkS4class{RasterLayer}} object from which cells with the value of 1 will be sampled.
## \link[package:MyClass-class]{MyClass}
#' @param grid_maskNA A \code{\linkS4class{RasterLayer}} object whose cells with the value of NA will prevent sampling of corresponding cells in \code{grid_mask1}.
#' @param initfile A filename (and its path) or \code{NA}. If not \code{NA}, then the generated search points are stored to disk in \link[=readRDS]{.rds} format.
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
get_transect_search_points <- function(N, grid_mask1, grid_maskNA = NULL, initfile = NA_character_, seed = NULL, verbose = TRUE) {
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
	
	if (ip_N < N) {
		initpoints <- sample_init_cells(N = N, grid = grid_mask1, gridNA = grid_maskNA, pts_init = initpoints, seed = seed, verbose = verbose)
		if (!is.na(initfile)) saveRDS(initpoints, file = initfile)
	} else if (ip_N > N) {
		initpoints <- initpoints[sample(x = length(initpoints), size = N), ]
	}
	
	initpoints
}
