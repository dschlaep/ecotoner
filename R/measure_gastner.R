

#------Gastner, M., Oborny, B., Zimmermann, D.K. & Pruessner, G. (2009) Transition from connected to fragmented vegetation across an environmental gradient: scaling laws in ecotone geometry. The American Naturalist, 174, E23-E39.

#---Gastner et al. 2009: percolation neighborhoods based on Fig. 1A
get.stepMaskAndDistance <- function(steplength, res_m) {
	sn <- matrix(NA, ncol = 1+2*steplength, nrow = 1+2*steplength)
	sn[1+steplength, 1+steplength] <- 0 #center
	rsn <- raster(sn, xmn = 0, xmx = ncol(sn) * res_m, ymn = 0, ymx = nrow(sn) * res_m, crs = "+proj = utm +zone = 13 +ellps = GRS80 +towgs84 = 0,0,0,0,0,0,0 +units = m +no_defs")	
	dsn <- distance(rsn, doEdge = FALSE)
	rsn[] <- dsn[] <= steplength * res_m
	sn <- raster::as.matrix(rsn, maxpixels = raster::ncell(rsn))
	sn[sn == 0] <- NA
	sn[1+steplength, 1+steplength] <- 0
	
	list(sn = sn, dsn = raster::as.matrix(dsn, maxpixels = raster::ncell(dsn)))
}


#---Gastner et al. 2009: Appendix 'Determination of the Hull Edge by a Biased Walk'
calc_Gastner2010_hulledge <- function(i, steplength, veg, end_toLeft) {
	helper_CRS <- raster::crs(veg)
	res_m <- raster::xres(veg)
	
	#--- 0. Calculate window properties
	stepwindow_N <- 1 + 2 * steplength
	stepwindow <- array(NA, dim = c(stepwindow_N, stepwindow_N, 5), dimnames = list(NULL, NULL, c("mask", "azimuth", "distance", "deltaX", "deltaY")))
	temp <- get.stepMaskAndDistance(steplength, res_m = res_m)
	stepwindow[,, "mask"] <- temp$sn
	stepwindow[,, "distance"] <- temp$dsn #unit = meters
	stepwindow[,, "deltaX"] <- matrix(temp <- rep(seq_len(stepwindow_N) - (1 + steplength), times = stepwindow_N), ncol = stepwindow_N, nrow = stepwindow_N, byrow = TRUE)
	stepwindow[,, "deltaY"] <- matrix(temp, ncol = stepwindow_N, nrow = stepwindow_N, byrow = FALSE)
	#azimuth = clockwise angle (in radians) between possible steps and general walking direction (if end_toLeft then 'up' = North else 'down' = South)
 	stepwindow[,, "azimuth"] <- rotate_azimuth_matrix(as.matrix(temp <- raster::atan2(raster(ifelse(end_toLeft, -1, 1) * stepwindow[,, "deltaX"], xmn = 0, xmx = stepwindow_N * res_m, ymn = 0, ymx = stepwindow_N * res_m, crs = helper_CRS), raster(ifelse(end_toLeft, 1, -1) * stepwindow[,, "deltaY"], xmn = 0, xmx = stepwindow_N * res_m, ymn = 0, ymx = stepwindow_N * res_m, crs = helper_CRS)), maxpixels = raster::ncell(temp)), 0)
		
	#--- 1. Identify largest connected patch
	#veg <- calc(veg, fun = function(x) ifelse(is.na(x), NA, 1))
	vegclumps_rookIDs <- clump(veg, directions = 4)
	#Rook's case
	vegclumps_largestID <- which.max((listFreqs <- freq(vegclumps_rookIDs))[!is.na(listFreqs[, 1]), 2])
	#vegclumps_largest <- calc(vegclumps_rookIDs, fun = function(x) ifelse(x %in% vegclumps_largestID & !is.na(x), 1, NA))
	vegclumps_largest <- calc(vegclumps_rookIDs, fun = function(x) ifelse(match(x, vegclumps_largestID, nomatch = 0) > 0 & !is.na(x), 1, NA))										
	#steplength larger than Rook's case
	if (steplength > 1) {
		vegclumps_largest_uniqueID <- 10^floor(log(sum(!is.na(listFreqs[, 1])), base = 10)+1)
		#vegclumps_outlying <- calc(vegclumps_rookIDs, fun = function(x) ifelse(x %in% vegclumps_largestID | is.na(x), NA, x)) #raster with clumps except the largest clump
		vegclumps_outlying <- calc(vegclumps_rookIDs, fun = function(x) ifelse(match(x, vegclumps_largestID, nomatch = 0) > 0 | is.na(x), NA, x)) #raster with clumps except the largest clump
		ir <- 1
		repeat { #include outlying patches into largest clump if they are within reach by steplength
			vegclumps_largest_enlarged <- focal(vegclumps_largest, w = stepwindow[,, "mask"], fun = function(x, ...) sum(!is.na(x)), pad = TRUE, padValue = NA) #raster with largest clump made larger by steplength: this call is by ca. 5x the most time consuming within this loop
			vegclumps_addedIDs <- raster::overlay(vegclumps_largest_enlarged, vegclumps_outlying, fun = function(x, y) ifelse(!is.na(x) & x > 0, vegclumps_largest_uniqueID, 0) + ifelse(!is.na(y), y, 0)) #raster with outlying clumps marked that overlap with enlarged largest clump, i.e., they are reachable from largest clump with a steplength
			addedIDs_freq <- unique(vegclumps_addedIDs)
			outlyingIDs_ToAdd <- addedIDs_freq[addedIDs_freq > vegclumps_largest_uniqueID] - vegclumps_largest_uniqueID #extract IDs of clumps that are reachable from largest clump by steplength
			if (length(outlyingIDs_ToAdd) > 0 & ir < 100) {
				#vegclumps_largest <- raster::overlay(vegclumps_largest, vegclumps_outlying, fun = function(x, y) ifelse(!is.na(x) | y %in% outlyingIDs_ToAdd, 1, NA)) #add outlying reachable clumps to largest clump
				vegclumps_largest <- raster::overlay(vegclumps_largest, vegclumps_outlying, fun = function(x, y) ifelse(!is.na(x) | match(y, outlyingIDs_ToAdd, nomatch = 0) > 0, 1, NA)) #add outlying reachable clumps to largest clump
				#vegclumps_outlying <- calc(vegclumps_outlying, fun = function(x) ifelse(is.na(x) | x %in% outlyingIDs_ToAdd, NA, x)) #remove outlying added clumps from remaining outlying clumps
				vegclumps_outlying <- calc(vegclumps_outlying, fun = function(x) ifelse(is.na(x) | match(x, outlyingIDs_ToAdd, nomatch = 0) > 0, NA, x)) #remove outlying added clumps from remaining outlying clumps
				ir <- ir + 1
			} else {
				break
			}
		}
		rm(addedIDs_freq, vegclumps_outlying, outlyingIDs_ToAdd, vegclumps_addedIDs, vegclumps_largest_enlarged)
	}
	rm(vegclumps_largestID)
	
	#---2. Identify cells to start and end biased walk
	runners <- calc_Eppinga2013_frontrunners(vegclumps_largest, end_toLeft)
	if (end_toLeft) { #left biased walk and from grid bottom to top 
		start_row <- max(which(!is.na(runners)))
		end_row <- min(which(!is.na(runners)))
	} else { #right biased walk and from grid top to bottom
		start_row <- min(which(!is.na(runners)))
		end_row <- max(which(!is.na(runners)))
	}
	start_col <- runners[start_row]
	end_col <- runners[end_row]
	
	is_spanning <- (yspan <- abs(end_row - start_row) + 1) == nrow(veg) #is largest connected patch spanning across y-axis of grid?
	
	#---3. Biased walk
	#Init
	grid_hullEdge <- raster::raster(veg)
	mat_hullEdge <- matrix(NA, nrow = 0, ncol = 2, dimnames = list(NULL, c("x", "y")))
	spLine_hullEdge <- sp::SpatialLines(list(sp::Lines(sp::Line(rbind(temp <- c(raster::xFromCol(veg, start_col), raster::yFromRow(veg, start_row)), temp)), ID = 1)), proj4string = sp::CRS(helper_CRS))
	iLine <- NULL
	iazimuth <- 0 #clockwise angle between current step and general walking direction
	irow <- start_row #matrix row (y-axis) of current position
	icol <- start_col #matrix column (x-axis) of current position
	nextStep_pos <- rep(NA, 2)
	#Walk
	repeat {
		#-Take step
		grid_hullEdge[irow, icol] <- 1
		mat_hullEdge <- rbind(mat_hullEdge, c(icol, irow))
		spLine_hullEdge <- if (!is.null(iLine)) rgeos::gLineMerge(rbind(spLine_hullEdge, iLine)) else spLine_hullEdge
		if (irow == end_row && icol == end_col) break

		#-Decide on next step
		#account for edges of the grid
		clipw <- c(ifelse(irow > steplength, 1, steplength + 2 - irow),
					ifelse((rtemp <- nrow(veg) - irow) >= steplength, stepwindow_N, steplength + 1 + rtemp),
					ifelse(icol > steplength, 1, steplength + 2 - icol),
					ifelse((ctemp <- ncol(veg) - icol) >= steplength, stepwindow_N, steplength + 1 + ctemp)) #usable extent of window (row min + max, col min + max)
		clipb <- c(ifelse(clipw[1] == 1, irow - steplength, 1),
					1 + ifelse(clipw[1] == 1, steplength, irow - 1) + ifelse(clipw[2] == stepwindow_N, steplength, rtemp),
					ifelse(clipw[3] == 1, icol - steplength, 1),
					1 + ifelse(clipw[3] == 1, steplength, icol - 1) + ifelse(clipw[4] == stepwindow_N, steplength, ctemp)) #usable extent of block (row min, nrows, col min, ncols)
		#apply stepwindow's mask
		neigh_sw <- stepwindow[clipw[1]:clipw[2], clipw[3]:clipw[4], ]
		neigh_vals <- getValuesBlock(vegclumps_largest, row = clipb[1], nrows = clipb[2], col = clipb[3], ncols = clipb[4], format = 'matrix')
		neigh_available <- neigh_vals * neigh_sw[,, "mask"]
		neigh_available[neigh_available == 0] <- NA
		#get the leftmost available cell(s) as seen from current direction
		neigh_azimuths <- neigh_available * rotate_azimuth_matrix(neigh_sw[,, "azimuth"], iazimuth)
		#loop through the cells with increasing leftmost cells and pick the first cell that doesn't cross previously walked path
		nextStep_pos <- rep(NA, 2)
		repeat {
			leftmost_value <- min(neigh_azimuths, na.rm = TRUE)
			if(is.infinite(leftmost_value)) break #no possible next steps to take
			leftmost_pos <- which(neigh_azimuths - sqrt(.Machine$double.neg.eps) <= leftmost_value & neigh_azimuths + sqrt(.Machine$double.eps) >= leftmost_value, arr.ind = TRUE, useNames = FALSE)
			#choose the nearest leftmost cell
			if((temp <- nrow(leftmost_pos)) > 1){
				leftmost_pos <- leftmost_pos[which.min(sapply(seq_len(temp), FUN = function(i) neigh_sw[leftmost_pos[i, 1], leftmost_pos[i, 2], "distance"])), ]
			} else {
				leftmost_pos <- drop(leftmost_pos)
			}
			#check that this step would not cross previously walked path c(yFromRow(veg, start_row), xFromCol(veg, start_col)
			iLine <- sp::SpatialLines(list(Lines(Line(rbind(c(xFromCol(veg, icol), yFromRow(veg, irow)), c(xFromCol(veg, icol + neigh_sw[leftmost_pos[1], leftmost_pos[2], "deltaX"]), yFromRow(veg, irow + neigh_sw[leftmost_pos[1], leftmost_pos[2], "deltaY"])))), ID = 0)), proj4string = CRS(helper_CRS))
			iIntersections <- rgeos::gIntersection(spLine_hullEdge, iLine)
			#remove points of cells along walk from this intersection (i.e., walk is allowed to revisit cells, but not to cross)
			temp1 <- apply(if(identical(class(temp <- coordinates(iIntersections)), "list")) as.matrix(temp[[1]][[1]]) else as.matrix(temp), 1, FUN = function(x) paste(round(x), collapse = "_"))
			temp2 <- apply(rbind(as.matrix(coordinates(iLine)[[1]][[1]]), as.matrix(coordinates(spLine_hullEdge)[[1]][[1]])), 1, FUN = function(x) paste(round(x), collapse = "_"))
			if (any(!(temp1 %in% temp2))) {
				neigh_azimuths[neigh_azimuths == leftmost_value] <- NA #remove this leftmost value from the possible choices
			} else { #next step found
				nextStep_pos <- leftmost_pos
				break
			}
		}
		if (all(is.na(nextStep_pos))) break #no next step: i.e. walk ended (probably prematurely)
		
		#-Prepare this step
		irow <- irow + neigh_sw[nextStep_pos[1], nextStep_pos[2], "deltaY"]
		icol <- icol + neigh_sw[nextStep_pos[1], nextStep_pos[2], "deltaX"]
		iazimuth <- rotate_azimuth_matrix(neigh_sw[nextStep_pos[1], nextStep_pos[2], "azimuth"], pi)
	}
		
	list(grid_LargestPatch = vegclumps_largest, is_spanning = is_spanning,
		grid_hullEdge = grid_hullEdge,
		mat_hullEdge = mat_hullEdge,
		spLine_hullEdge = spLine_hullEdge)
}

#---Gastner et al. 2009: Appendix 'Width and Length of the Hull'
calc_Gastner2010_hulledge_Statistics <- function(mat_hullEdge, spLine_hullEdge, res_m) {
	pos <- mean(mat_hullEdge[, 1])	#eq. A2
	width <- sd(mat_hullEdge[, 1])	#eq. A3
	len <- rgeos::gLength(spLine_hullEdge)	#eq. A4
		
	list(position_cells = pos, position_m = res_m * pos, width_m = res_m * width, length_m = len)
}

calc_Gastner2010_hulledge_Gradient <- function(elev, slope, grid_hullEdge) {
	telev <- raster::overlay(elev, grid_hullEdge, fun = function(x, y) ifelse(!is.na(y), x, NA))
	elev_Mean <- raster::cellStats(telev, 'mean')
	elev_SD <- raster::cellStats(telev, 'sd')
	tslope <- raster::overlay(slope, grid_hullEdge, fun = function(x, y) ifelse(!is.na(y), x, NA))
	slope_Mean <- raster::cellStats(tslope, 'mean')
	slope_SD <- raster::cellStats(tslope, 'sd')
		
	
	list(Elev_Mean_m = elev_Mean, Elev_SD_m = elev_SD, Slope_Mean_rad = slope_Mean, Slope_SD_rad = slope_SD)
}

calc_Gastner2010_hulledge_DistanceToStruggleZone <- function(vegOther, mat_hullEdge, end_toLeft, width_N) {
	#distances between hullEdge and front runners of other vegetation; negative if front runners of other vegetation stop 'before' hull edge of vegetation (i.e, no interaction)
	runners <- calc_Eppinga2013_frontrunners(vegOther, !end_toLeft)
	
	ifelse(end_toLeft, -1, 1) * raster::xres(vegOther) * sapply(seq_len(width_N), FUN = function(irow) {temp <- mat_hullEdge[mat_hullEdge[, 2] == irow, 1]; return(ifelse(length(temp) > 0, ifelse(end_toLeft, min(temp, na.rm = TRUE), max(temp, na.rm = TRUE)), NA) - runners[irow])})
}

tabulate_Gastner2010_hulledge <- function(etable, b, data, steplength, width_N, flag_migtype) {

	place_veg <- function(etable, b, data, steplength, vegno, flag_migtype, width_N) {
		colnamesAdd <- paste0(flag_migtype, "_Gastner2009_HullEdge_Step", steplength, "_Veg", vegno, "_", 
						c("Location_m", "Width_m", "Length_m", "Elev_Mean_m", "Elev_SD_m", 
						"Slope_Mean_rad", "Slope_SD_rad", "LargestPatch_isSpanning_TF", "LargestPatch_Area_fractionOfTransect",
						"VegDensity_atMeanHullPosition",
						"MeanDistanceHullEdgeToOtherVegFrontRunners_m", "RowsOtherVegFrontRunnersCloserToLargestPatchOfHullEdge_Fraction", "HullEdgeAffectedByOtherVeg_TF"))
		res <- matrix(NA, nrow = max(1, nrow(etable)), ncol = length(colnamesAdd), dimnames = list(NULL, colnamesAdd))

		res[b, 1] <- data$stats$position_m
		res[b, 2] <- data$stats$width_m
		res[b, 3] <- data$stats$length_m
		res[b, 4] <- data$grad$Elev_Mean_m
		res[b, 5] <- data$grad$Elev_SD_m
		res[b, 6] <- data$grad$Slope_Mean_rad
		res[b, 7] <- data$grad$Slope_SD_rad
		res[b, 8] <- data$is_spanning
		res[b, 9] <- raster::cellStats(data$grid_LargestPatch, 'sum') / raster::ncell(data$grid_LargestPatch)
		res[b, 10] <- data$VegDensity_atMeanHullPosition
		res[b, 11] <- mean(data$HullEdgeDistToOtherVeg_m, na.rm = TRUE)
		res[b, 12] <- sum(data$HullEdgeDistToOtherVeg_m > 0, na.rm = TRUE) / width_N
		res[b, 13] <- res[b, 12] > 0
		
		cbind(etable, res)
	}
	
	etable <- place_veg(etable, b, data = data$Veg1, steplength, vegno = 1, flag_migtype = flag_migtype, width_N = width_N)
	etable <- place_veg(etable, b, data = data$Veg2, steplength, vegno = 2, flag_migtype = flag_migtype, width_N = width_N)
	
	etable
}



Gastner2010AmNat <- function(i, b, migtype, ecotoner_settings, etband, etmeasure, flag_bfig, copy_FromMig1_TF, do_figures, ...) {
	#3c. Gastner et al. 2010 Am.Nat.: Location of boundary as hull edge
	
	dots <- list(...)

	etmeasure$etable[b, "Transect_ID"] <- i
	etmeasure$etable[b, "Neighbor_Cells"] <- neighborhoods(ecotoner_settings)[b]
	
	step_lengths <- stepsHullEdge(ecotoner_settings)
	etmeasure$gETmeas[[b]][[migtype]] <- vector(mode = "list", length = length(step_lengths))
	step_names <- paste0("Step", step_lengths)
	names(etmeasure$gETmeas[[b]][[migtype]]) <- step_names
	
	for (i_step in seq_along(step_names)) {
		istep <- step_names[i_step]
		if (!copy_FromMig1_TF) {
			for (iveg in c("Veg1", "Veg2")) {
				iother <- setdiff(c("Veg1", "Veg2"), iveg)
				type_veg <- if (iveg == "Veg1") type_veg1(ecotoner_settings) else if (iveg == "Veg2") type_veg2(ecotoner_settings) else NULL
			
				etmeasure$gETmeas[[b]][[migtype]][[istep]][[iveg]] <- calc_Gastner2010_hulledge(i,
																steplength = step_lengths[i_step],
																veg = etband$Veg[[migtype]][[iveg]]$grid,
																end_toLeft = end_to_left(type_veg))
				etmeasure$gETmeas[[b]][[migtype]][[istep]][[iveg]]$Hull_density <- count_y_at_each_x(etmeasure$gETmeas[[b]][[migtype]][[istep]][[iveg]]$grid_hullEdge)
				etmeasure$gETmeas[[b]][[migtype]][[istep]][[iveg]]$stats <- calc_Gastner2010_hulledge_Statistics(etmeasure$gETmeas[[b]][[migtype]][[istep]][[iveg]]$mat_hullEdge,
																etmeasure$gETmeas[[b]][[migtype]][[istep]][[iveg]]$spLine_hullEdge,
																res_m = raster::xres(etband$Env$elev$grid))
				etmeasure$gETmeas[[b]][[migtype]][[istep]][[iveg]]$VegDensity_atMeanHullPosition <- etband$Veg[[migtype]][[iveg]]$density[round(etmeasure$gETmeas[[b]][[migtype]][[istep]][[iveg]]$stats$position_cells)]
				etmeasure$gETmeas[[b]][[migtype]][[istep]][[iveg]]$grad <- calc_Gastner2010_hulledge_Gradient(etband$Env$elev$grid,
																etband$Env$slope$grid,
																etmeasure$gETmeas[[b]][[migtype]][[istep]][[iveg]]$grid_hullEdge)
				etmeasure$gETmeas[[b]][[migtype]][[istep]][[iveg]]$HullEdgeDistToOtherVeg_m <- calc_Gastner2010_hulledge_DistanceToStruggleZone(vegOther = etband$Veg[[migtype]][[iother]]$grid,
																mat_hullEdge = etmeasure$gETmeas[[b]][[migtype]][[istep]][[iveg]]$mat_hullEdge,
																end_toLeft = end_to_left(type_veg),
																width_N = bandTransect_width_cellN(ecotoner_settings))
			}


			if (do_figures) plot_Gastner2010_hulledge(filename = file.path(dir_fig, paste0(flag_bfig, "Gastner2009_FittedHullEdge_step", step_lengths[i_step], "_", migtype, ".pdf")),
													eB_Env = etband$Env,
													eB_Veg = etband$Veg[[migtype]],
													datFit = etmeasure$gETmeas[[b]][[migtype]][[istep]])
		} 
		
		etmeasure$etable <- tabulate_Gastner2010_hulledge(etable = etmeasure$etable, b = b,
															data = etmeasure$gETmeas[[b]][[if (copy_FromMig1_TF) 1 else migtype]][[istep]],
															steplength = step_lengths[i_step],
															width_N = bandTransect_width_cellN(ecotoner_settings),
															flag_migtype = migtype)
	}
	
	etmeasure
}

