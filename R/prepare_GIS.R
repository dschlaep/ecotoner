#' Mosaic raster tiles
#'
#' This function can be used to mosaic the NED tiles downloaded from nationalmap.gov
#'
#' If 'fname_grid_ned' already exists and agree with origin, resolution, and coordinate reference system with the tiles, then tiles are added instead of creating a new raster object.
#' 
#' @param dir_tiles A character string. The path to the directory where the folders of each tile are stored on disk.
#' @param chunksize An integer value. The number of tiles that are mosaic-ed per function call.
#' @param fname_grid_ned A character string. The path with file name (with or without file extension) of the final NED raster.
#' @param format A character string or missing. The file type. For possible values see \code{\link{raster::writeFormats}}.
#' @param ... Additional arguments for writing raster files to disk as in \code{\link{raster::writeRaster}}.
#' @seealso This function calls \code{\link{raster::mosaic}}
#' @return A RasterLayer object stored on disk as \code{fname_grid_ned}
#' @export
mosaic_tiles <- function(dir_tiles, chunksize = 10L, fname_grid_ned, format, ...) {
	#---Test input
	dir_create(dirname(fname_grid_ned))
	chunksize <- as.integer(chunksize)

	filetype <- raster:::.filetype(format = format, filename = fname_grid_ned)
	if (!(filetype %in% raster::writeFormats()[, "name"]))
		stop("ecotoner::mosaic_NED_tiles(): file type '", filetype, "' is not supported")
	fname_grid_ned <- raster:::.getExtension(fname_grid_ned, filetype)

	
	#---Begin function calculations
	# Read the individual tiles to a list of RasterLayer objects
	grid_tiles <- lapply(list.dirs(path = dir_tiles, full.names = TRUE, recursive = FALSE), function(x) {
						ftemp1 <- file.path(x, paste0("grd", basename(x), "_1"))
						if (file.exists(ftemp1)) {
							raster::raster(ftemp1)
						} else {
							# Sometimes the grids don't have the "_1" appended
							ftemp2 <- file.path(x, paste0("grd", basename(x)))
							if (file.exists(ftemp2)) raster::raster(ftemp2) else NULL
						}
					})
	
	# Subset to only those tiles that have the same origin, resolution, and coordinate reference system as the first tile
	if (length(grid_tiles) > 0) {
		good_tiles <- sapply(grid_tiles, function(x) raster::compareRaster(grid_tiles[[1]], x, extent = FALSE, rowcol = FALSE, crs = TRUE, res = TRUE, orig = TRUE, rotation = TRUE, values = FALSE, stopiffalse = FALSE, showwarning = FALSE))
		if (any(!good_tiles)) {
			cat("ecotoner::mosaic_NED_tiles: ", any(!good_tiles), " tiles deviate in origin, resolution, and/or coordinate reference system from the first tile '", names(grid_tiles[[1]]), "' and will be removed prior to mosaic-ing", "\n", sep = "")
			grid_tiles <- grid_tiles[good_tiles]
		}
	}
	
	# Use the 'final' target file if it exists and fits the mosaic criteria
	if (file.exists(fname_grid_ned)) {
		temp <- raster::raster(fname_grid_ned)
		if (length(grid_tiles) > 0) {
			if (raster::compareRaster(temp, grid_tiles[[1]], extent = FALSE, rowcol = FALSE, crs = TRUE, res = TRUE, orig = TRUE, rotation = TRUE, values = FALSE, stopiffalse = FALSE, showwarning = FALSE)) {
				grid_ned <- temp
			} else {
				grid_ned <- grid_tiles[[1]]
				grid_tiles <- grid_tiles[-1]
			}
		} else {
			grid_ned <- temp
		}
	} else {
		grid_ned <- grid_tiles[[1]]
		grid_tiles <- grid_tiles[-1]
	}
	
	n <- length(grid_tiles)
	
	# Mosaic the tiles
	ftemp2 <- rep(NA_character_, 2)
	if (n > 0) {
		cat("ecotoner::mosaic_NED_tiles: ", format.POSIXct(Sys.time(), ""), "; n = ", n, " tiles will be mosaic-ed", "\n", sep = "")
		
		dir_temp <- dirname(fname_grid_ned)
		todos <- parallel:::splitList(seq_len(n), n / chunksize)
		
		i <- 1
		repeat {
			ftemp2 <- c(ftemp2[2], file.path(dir_temp, paste0(format.POSIXct(Sys.time(), "%Y%m%d_%H%M"), "_temp_raster.tif")))
			
			if (length(todos[[i]]) > 0) {
				args_to_mosaic <- if (length(todos[[i]]) > 1) {
										c(list(x = grid_ned, y = grid_tiles[[todos[[i]][1]]]), c(grid_tiles[todos[[i]][-1]], list(...)), list(fun = mean, filename = ftemp2[2]))
									} else {
										list(x = grid_ned, y = grid_tiles[[todos[[i]][1]]], fun = mean, filename = ftemp2[2], ...)
									}

				grid_ned <- do.call(getFromNamespace("mosaic", "raster"), args = args_to_mosaic)
			
				cat("ecotoner::mosaic_NED_tiles: ", format.POSIXct(Sys.time(), ""), "; ", length(unlist(todos[1:i])), " out of ", n, " tiles mosaic-ed", "\n", sep = "")
			}
			
			file.remove(ftemp2[1])
			if (i < length(todos)) i <- i + 1 else break
		}
	}
	
	done <- try(raster::writeRaster(grid_ned, filename = fname_grid_ned, ...), silent = TRUE)
	if (!inherits(done, "try-error")) file.remove(ftemp2[2])
	
	grid_ned
}


#' Project a raster
#' 
#' This function may be very SLOW and memory intensive!
#' 
#' @param grid_from A RasterLayer object. The raster to be projected
#' @param fname_grid_to A character string. The path with file name (with or without file extension) of the resulting projected raster.
#' @param res_to The raster cell resolution of the target.
#' @param crs_to The CRS of the target.
#' @param parallel_N An integer value. If larger than 1, then the re-projection will use a snow cluster, see \code{\link{raster::beginCluster}} for details
#' @param chunksize An integer value or missing. Maximum number of cells to read/write in a single chunk while processing (chunk by chunk) disk based Raster* objects (see \code{\link{raster::rasterOptions}} for details)
#' @param maxmemory An integer value or missing. Maximum number of cells to read into memory. I.e., if a Raster* object has more than this number of cells, \code{\link{raster::canProcessInMemory}} will return \code{FALSE} (see \code{\link{raster::rasterOptions}} for details)
#' @param ... Additional arguments for writing raster files to disk as in \code{\link{raster::writeRaster}}.
#' @seealso This function calls \code{\link{raster::projectRaster}}
#' @return A RasterLayer object stored on disk as \code{fname_grid_to}
#' @export
project_raster <- function(grid_from, fname_grid_to, res_to, crs_to, parallel_N = 1, ...) {
	parallel_N <- as.integer(parallel_N)
	cat("ecotoner::project_raster: ", format.POSIXct(Sys.time(), ""), ": raster will be projected using n = ", parallel_N, " cores; NOTE: this call may be very SLOW and memory intensive!", "\n", sep = "")

	# Prepare resulting grid
	grid_to <- raster::projectExtent(grid_from, crs = crs_to)
	raster::res(grid_to) <- res_to
	
	# Set up work zone
	dots <- list(...)
	op_old <- raster::rasterOptions(todisk = FALSE)	
	if (!is.null(dots[["chunksize"]])) {
		raster::rasterOptions(chunksize = dots[["chunksize"]])
		dots[["chunksize"]] <- NULL
		on.exit(raster::rasterOptions(op_old[["chunksize"]]), add = TRUE)
	}
	if (!is.null(dots[["maxmemory"]])) {
		raster::rasterOptions(maxmemory = dots[["maxmemory"]])
		dots[["maxmemory"]] <- NULL
		on.exit(raster::rasterOptions(op_old[["maxmemory"]]), add = TRUE)
	}
	if (parallel_N > 1) {
		raster::beginCluster(n = parallel_N)
		on.exit(raster::endCluster(), add = TRUE)
	}

	# Project with progress report
	dots$from <- grid_from
	dots$to <- grid_to
	dots$method <- "bilinear"
	dots$filename <- fname_grid_to
	dots$progress <- "text"
	grid_res <- do.call(getFromNamespace("projectRaster", "raster"), args = dots)
	
	grid_res
}



#' @export
determine_abutters <- function(grid, grid_veg1, grid_veg2, filename, asNA = TRUE, ...) {
	calc_abutting(Veg1 = grid_veg1, Veg2 = grid_veg2, filename, asNA = asNA, ...)
}

#' @export			
extract_vegetation <- function(grid, ids, filename, parallel_N, ...) {
	parallel_N <- as.integer(parallel_N)
	cat("ecotoner::extract_vegetation: ", format.POSIXct(Sys.time(), ""), ": vegetation will be extracted using n = ", parallel_N, " cores; NOTE: this call may be very SLOW and memory intensive!", "\n", sep = "")

	if (parallel_N > 1) {
		raster::beginCluster(n = parallel_N)
		on.exit(raster::endCluster(), add = TRUE)

		raster::clusterR(x = grid, fun = grid_to_NA1, args = list(vals = ids),
						filename = filename, ...)
	} else {
		grid_to_NA1(grid, ids, filename, ...)
	}
}


#' @export
terrain_slope <- function(grid_elev, filename, ...) {
	raster::terrain(x = grid_elev, opt = 'slope', unit = 'radians', neighbors = 8, filename = filename, ...) 
}

#' @export
terrain_aspect <- function(grid_elev, grid_slope, min_slope, parallel_N, filename, ...) {
	cat("ecotoner::terrain_aspect: ", format.POSIXct(Sys.time(), ""), ": aspect will be calculated using n = ", parallel_N, " cores; NOTE: this call may be very SLOW and memory intensive!", "\n", sep = "")

	grid_temp <- raster::terrain(x = grid_elev, opt = 'aspect', unit = 'radians', neighbors = 8, ...) 

	fun_asp_smin <- if (parallel_N > 1) {
						compiler::cmpfun(function(x) ifelse(x[2] >= min_slope, x[1], NA))
					} else {
						compiler::cmpfun(function(x, y) ifelse(y >= min_slope, x, NA))
					}
	
	if (parallel_N > 1) {
		raster::beginCluster(n = parallel_N)
		on.exit(raster::endCluster(), add = TRUE)

		grid_aspect <- raster::clusterR(x = raster::brick(grid_temp, grid_slope), fun = fun_asp_smin, filename = filename, ...)
	} else {
		grid_aspect <- raster::overlay(grid_temp, grid_slope, fun = fun_asp_smin, filename = filename, ...)
	}
	
	grid_aspect
}


#	asp201Mean <- focal(aspectCropped, w = ifelse(is_odd(width_N), width_N, width_N+1), fun = function(x, ...) circ_mean(x, int = 2*pi, na.rm = TRUE), pad = TRUE, padValue = NA, filename = file.path(dir_fig, "asp201MeanCroppedF.tif"))
#	asp201SD <- focal(aspectCropped, w = ifelse(is_odd(width_N), width_N, width_N+1), fun = function(x, ...) circ_sd(x, int = 2*pi, na.rm = TRUE), pad = TRUE, padValue = NA)

get_window_size <- function(N) {
	N <- as.integer(N)
	if (!is_odd(N)) N <- N + 1
	if (N > 0) N else 0
}

#' @export
homogenous_aspect <- function(grid_aspect, fun = c("mean", "sd"), window_N, filename, ...) {
	cat("ecotoner::homogenous_aspect: ", format.POSIXct(Sys.time(), ""), ": smoothed ", fun, " of aspect will be calculated; NOTE: this call may be very SLOW and memory intensive!", "\n", sep = "")

	fun <- match.arg(fun)
	window_N <- get_window_size(window_N)
	w <- raster::focalWeight(grid_aspect,
					d = get_window_size(window_N) * raster::xres(grid_aspect) / 2,
					type = "rectangle")
	
	fun_circ <- switch(fun,
						mean = function(x, ...) circ_mean(x, int = 2 * pi, na.rm = TRUE),
						sd = function(x, ...) circ_sd(x, int = 2 * pi, na.rm = TRUE))

	raster::focal(grid_aspect, w = w, fun = fun_circ, pad = TRUE, padValue = NA,
					filename = filename, progress = "text", ...)
}										


#' @export
characterize_veg_data <- function(ecotoner_settings, ecotoner_grids, initpoints, inhibit_dist) {
	ftemp <- file.path(dir_out(ecotoner_settings), "Metadata_vegetation.rds")
	
	if (file.exists(ftemp)) {
		res <- readRDS(ftemp)
	} else {

		# grid prevalence
		lgrids <- list(veg = grid_veg(ecotoner_grids), veg1 = grid_veg1(ecotoner_grids), veg2 = grid_veg2(ecotoner_grids), abut = grid_abut(ecotoner_grids))
		dat_grids <- sapply(lgrids, FUN = function(x)
								c(	ncell_total = temp <- raster::ncell(x), 
									ncell_narm = temp -  raster::cellStats(x, "countNA"),
									sum_values = raster::cellStats(x, "sum"),
									fsize_GB = 1e-9 * file.size(x@file@name)))
								
		prevalence <- dat_grids["ncell_narm", ] / dat_grids["ncell_narm", "veg"]
	
	
		# packing intensity
		row.names(initpoints) <- 1:length(initpoints)
		# method with gBuffer: is much slower, but seems to be more accurate
		initpoints_buffering <- rgeos::gBuffer(initpoints, byid = TRUE, id = 1:length(initpoints), width = inhibit_dist / 2)
		cells_buffering <- raster::extract(grid_abut(ecotoner_grids), y = initpoints_buffering, fun = sum, na.rm = TRUE, df = TRUE)

		cells_total <- dat_grids["ncell_narm", "abut"]
		tau_packing <- sum(cells_buffering[, "gapv2_bseABUTtf"]) / cells_total

		# method with extract(buffer): seems to be less accurate and slightly more inclusive
		if (FALSE) {
			cells_buffering2 <- raster::extract(grid_abut(ecotoner_grids), y = initpoints,
												buffer = inhibit_dist / 2, fun = sum, na.rm = TRUE, df = TRUE)
			plot(cells_buffering2[, "gapv2_bseABUTtf"], cells_buffering[, "gapv2_bseABUTtf"])
		
			tau_packing2 <- sum(cells_buffering2[, "gapv2_bseABUTtf"]) / cells_total
		}

		res <- list(dat_grids = dat_grids,
					prevalence = prevalence,
					tau_packing = tau_packing,
					initpoints_buffering = initpoints_buffering,
					cells_buffering = cells_buffering)
					
		saveRDS(res, file = ftemp)
	}

	res
}
