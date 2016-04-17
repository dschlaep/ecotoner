#------Gastner, M., Oborny, B., Zimmermann, D.K. & Pruessner, G. (2009) Transition from connected to fragmented vegetation across an environmental gradient: scaling laws in ecotone geometry. The American Naturalist, 174, E23-E39.

version_Gastner2010AmNat <- function() numeric_version("0.1.0")


#---Gastner et al. 2009: percolation neighborhoods based on Fig. 1A
get.stepMaskAndDistance <- function(steplength, res_m, pk = 100) {
	sn <- array(NA, dim = rep(1 + 2 * steplength, 2))
	sn[1 + steplength, 1 + steplength] <- 0 #center from which to compute distances
	
	dsn <- raster::raster(sn, xmn = 0, xmx = ncol(sn) * res_m, ymn = 0, ymx = nrow(sn) * res_m, crs = "+proj=utm +zone=13 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")	
	dsn <- raster::as.matrix(raster::distance(dsn, doEdge = FALSE))
	dsn <- round(dsn * pk / res_m) * res_m / pk # round to res_m / pk
	
	sn <- dsn <= {steplength * res_m + sqrt(.Machine$double.eps)}
	sn[!sn] <- NA
	mode(sn) <- "integer"
	sn[1 + steplength, 1 + steplength] <- 0 #center
	
	list(sn = sn, dsn = dsn)
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
 	rtemp1 <- rtemp2 <- raster::raster(stepwindow[,, "azimuth"],
 							xmn = 0, xmx = stepwindow_N * res_m,
 							ymn = 0, ymx = stepwindow_N * res_m,
 							crs = helper_CRS)
	rtemp1 <- raster::setValues(rtemp1, {if (end_toLeft) -1 else 1} * stepwindow[,, "deltaX"])
	rtemp2 <- raster::setValues(rtemp2, {if (end_toLeft) 1 else -1} * stepwindow[,, "deltaY"])
 	rtemp <- raster::atan2(rtemp1, rtemp2)
 	mat <- raster::as.matrix(rtemp, maxpixels = raster::ncell(rtemp))
 	stepwindow[,, "azimuth"] <- rotate_azimuth_matrix(mat, 0)
		
	#--- 1. Identify largest connected patch
	#veg <- calc(veg, fun = function(x) ifelse(is.na(x), NA, 1))
	vegclumps_rookIDs <- raster::clump(veg, directions = 4)
	#Rook's case
	vegclumps_largestID <- which.max((listFreqs <- raster::freq(vegclumps_rookIDs))[!is.na(listFreqs[, 1]), 2])
	vegclumps_largest <- raster::calc(vegclumps_rookIDs, fun = function(x) ifelse(match(x, vegclumps_largestID, nomatch = 0) > 0 & !is.na(x), 1, NA))										
	#steplength larger than Rook's case
	if (steplength > 1) {
		vegclumps_largest_uniqueID <- 10^floor(log(sum(!is.na(listFreqs[, 1])), base = 10)+1)
		vegclumps_outlying <- raster::calc(vegclumps_rookIDs, fun = function(x) ifelse(match(x, vegclumps_largestID, nomatch = 0) > 0 | is.na(x), NA, x)) #raster with clumps except the largest clump
		ir <- 1
		repeat { #include outlying patches into largest clump if they are within reach by steplength
			vegclumps_largest_enlarged <- raster::focal(vegclumps_largest, w = stepwindow[,, "mask"], fun = function(x, ...) sum(!is.na(x)), pad = TRUE, padValue = NA) #raster with largest clump made larger by steplength: this call is by ca. 5x the most time consuming within this loop
			vegclumps_addedIDs <- raster::overlay(vegclumps_largest_enlarged, vegclumps_outlying, fun = function(x, y) ifelse(!is.na(x) & x > 0, vegclumps_largest_uniqueID, 0) + ifelse(!is.na(y), y, 0)) #raster with outlying clumps marked that overlap with enlarged largest clump, i.e., they are reachable from largest clump with a steplength
			addedIDs_freq <- raster::unique(vegclumps_addedIDs)
			outlyingIDs_ToAdd <- addedIDs_freq[addedIDs_freq > vegclumps_largest_uniqueID] - vegclumps_largest_uniqueID #extract IDs of clumps that are reachable from largest clump by steplength
			if (length(outlyingIDs_ToAdd) > 0 & ir < 100) {
				vegclumps_largest <- raster::overlay(vegclumps_largest, vegclumps_outlying, fun = function(x, y) ifelse(!is.na(x) | match(y, outlyingIDs_ToAdd, nomatch = 0) > 0, 1, NA)) #add outlying reachable clumps to largest clump
				vegclumps_outlying <- raster::calc(vegclumps_outlying, fun = function(x) ifelse(is.na(x) | match(x, outlyingIDs_ToAdd, nomatch = 0) > 0, NA, x)) #remove outlying added clumps from remaining outlying clumps
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
	mat_hullEdge <- matrix(NA, nrow = 0, ncol = 3, dimnames = list(NULL, c("x", "y", "azimuth")))
	spLine_hullEdge <- sp::SpatialLines(list(sp::Lines(sp::Line(rbind(temp <- c(raster::xFromCol(veg, start_col), raster::yFromRow(veg, start_row)), temp)), ID = 1)), proj4string = helper_CRS)
	iLine <- NULL
	iazimuth <- 0 #clockwise angle between current step and general walking direction
	irow <- start_row #matrix row (y-axis) of current position
	icol <- start_col #matrix column (x-axis) of current position

	#Walk
	if (requireNamespace("rgeos", quietly = TRUE)) {
		repeat {
			#-Take step
			grid_hullEdge[irow, icol] <- 1
			mat_hullEdge <- rbind(mat_hullEdge, c(icol, irow, iazimuth))
			spLine_hullEdge <- if (!is.null(iLine)) rgeos::gLineMerge(sp::rbind.SpatialLines(spLine_hullEdge, iLine)) else spLine_hullEdge
			if (irow == end_row && icol == end_col) break

			#-Decide on next step
			#account for edges of the grid
			clipw <- c(if (irow > steplength) 1 else steplength + 2 - irow,
						if ({rtemp <- nrow(veg) - irow} >= steplength) stepwindow_N else {steplength + 1 + rtemp},
						if (icol > steplength) 1 else steplength + 2 - icol,
						if ({ctemp <- ncol(veg) - icol} >= steplength) stepwindow_N else {steplength + 1 + ctemp}) #usable extent of window (row min + max, col min + max)
			clipb <- c(if (clipw[1] == 1) irow - steplength else 1,
						1 + {if (clipw[1] == 1) steplength else {irow - 1}} +
							{if (clipw[2] == stepwindow_N) steplength else rtemp},
						if (clipw[3] == 1) {icol - steplength} else 1,
						1 + {if (clipw[3] == 1) steplength else {icol - 1}} +
							{if (clipw[4] == stepwindow_N) steplength else ctemp}) #usable extent of block (row min, nrows, col min, ncols)

			#apply stepwindow's mask
			neigh_sw <- stepwindow[clipw[1]:clipw[2], clipw[3]:clipw[4], ]
			neigh_vals <- raster::getValuesBlock(vegclumps_largest, row = clipb[1], nrows = clipb[2], col = clipb[3], ncols = clipb[4], format = 'matrix')
			neigh_available <- neigh_vals * neigh_sw[,, "mask"]
			neigh_available[neigh_available == 0] <- NA
			#get the leftmost available cell(s) as seen from current direction
			neigh_azimuths <- neigh_available * rotate_azimuth_matrix(neigh_sw[,, "azimuth"], iazimuth)
			
			#loop through the cells with increasing leftmost cells and pick the first cell that doesn't cross previously walked path
			take_nextStep <- FALSE
			if (any(!is.na(neigh_azimuths))) repeat {
				leftmost_value <- min(neigh_azimuths, na.rm = TRUE)
				if (is.infinite(leftmost_value)) break #no possible next steps to take
				leftmost_pos <- which(neigh_azimuths - sqrt(.Machine$double.neg.eps) <= leftmost_value & neigh_azimuths + sqrt(.Machine$double.eps) >= leftmost_value, arr.ind = TRUE, useNames = FALSE)
				#choose the nearest leftmost cell
				leftmost_pos <- if ((temp <- nrow(leftmost_pos)) > 1) {
									leftmost_pos[which.min(sapply(seq_len(temp), FUN = function(i) neigh_sw[leftmost_pos[i, 1], leftmost_pos[i, 2], "distance"])), ]
								} else {
									drop(leftmost_pos)
								}
				
				icol_next <- icol + neigh_sw[leftmost_pos[1], leftmost_pos[2], "deltaX"]
				irow_next <- irow + neigh_sw[leftmost_pos[1], leftmost_pos[2], "deltaY"]
				iazimuth_next <- rotate_azimuth_matrix(neigh_sw[leftmost_pos[1], leftmost_pos[2], "azimuth"], pi)
				
				#check that this step would not cross previously walked path c(yFromRow(veg, start_row), xFromCol(veg, start_col)
				has_next1 <- FALSE
				iLine <- sp::SpatialLines(list(sp::Lines(sp::Line(rbind(
									c(raster::xFromCol(veg, icol), raster::yFromRow(veg, irow)),
									c(raster::xFromCol(veg, icol_next), raster::yFromRow(veg, irow_next)))),
								ID = 0)), proj4string = helper_CRS)
				iIntersections <- rgeos::gIntersection(spLine_hullEdge, iLine)
				if (is.null(iIntersections)) {
					has_next1 <- TRUE
				} else {
					#remove points of cells along walk from this intersection (i.e., walk is allowed to revisit cells, but not to cross)
					temp1 <- apply(if (inherits(temp <- sp::coordinates(iIntersections), "list")) as.matrix(temp[[1]][[1]]) else as.matrix(temp), 1, FUN = function(x) paste(round(x), collapse = "_"))
					temp2 <- apply(rbind(as.matrix(sp::coordinates(iLine)[[1]][[1]]), as.matrix(sp::coordinates(spLine_hullEdge)[[1]][[1]])), 1, FUN = function(x) paste(round(x), collapse = "_"))
				
					if (any(!(temp1 %in% temp2))) {
						neigh_azimuths[neigh_azimuths == leftmost_value] <- NA #remove this leftmost value from the possible choices
					} else has_next1 <- TRUE
				}
				
				# check that this step would not retrace previously walked path in the same direction
				this_sequence <- rbind(mat_hullEdge[nrow(mat_hullEdge), ],
										c(icol_next, irow_next, iazimuth_next))
				sequences_unique <- sapply(1:nrow(mat_hullEdge), FUN = function(i) {
										if (i > 1) {
											sum(abs(mat_hullEdge[(i - 1):i, ] - this_sequence)) > sqrt(.Machine$double.eps)
										} else TRUE
									})
				has_next2 <- all(sequences_unique)
				
				# next step found
				if (has_next1 && has_next2) {
					take_nextStep <- TRUE
					break
				}
			}

			if (!take_nextStep) break #no next step: i.e. walk ended (probably prematurely)
		
			#-Prepare this step
			irow <- irow_next
			icol <- icol_next
			iazimuth <- iazimuth_next
		}
	} else {
		warning("Package 'rgeos' not installed: 'calc_Gastner2010_hulledge' will not performe the biased walk")
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
	if (requireNamespace("rgeos", quietly = TRUE)) {
		len <- rgeos::gLength(spLine_hullEdge)	#eq. A4
	} else {
		warning("Package 'rgeos' not installed: 'calc_Gastner2010_hulledge_Statistics' will not determine the length of the hull edge")
	}
		
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
	
	{if (end_toLeft) -1 else 1} * raster::xres(vegOther) *
		sapply(seq_len(width_N), FUN = function(irow) {
				temp <- mat_hullEdge[mat_hullEdge[, 2] == irow, 1]
				{if (length(temp) > 0) {
					if (end_toLeft) min(temp, na.rm = TRUE) else max(temp, na.rm = TRUE)
				} else NA} - runners[irow]
	})
}

tabulate_Gastner2010_hulledge <- function(etable, index, data, steplength, width_N) {

	place_veg <- function(etable, index, data, steplength, vegno, width_N) {
		colnamesAdd <- paste0("Gastner2009_HullEdge_Step",
						steplength, "_Veg", vegno, "_", 
						c("Location_m", "Width_m", "Length_m", "Elev_Mean_m", "Elev_SD_m", 
							"Slope_Mean_rad", "Slope_SD_rad", "LargestPatch_isSpanning_TF", "LargestPatch_Area_fractionOfTransect",
							"VegDensity_atMeanHullPosition",
							"MeanDistanceHullEdgeToOtherVegFrontRunners_m", "RowsOtherVegFrontRunnersCloserToLargestPatchOfHullEdge_Fraction", "HullEdgeAffectedByOtherVeg_TF"))
		res <- as.data.frame(t(rep(NA_real_, length(colnamesAdd))))
		colnames(res) <- colnamesAdd

		res[1] <- data$stats$position_m
		res[2] <- data$stats$width_m
		res[3] <- data$stats$length_m
		res[4] <- data$grad$Elev_Mean_m
		res[5] <- data$grad$Elev_SD_m
		res[6] <- data$grad$Slope_Mean_rad
		res[7] <- data$grad$Slope_SD_rad
		res[8] <- data$is_spanning
		res[9] <- raster::cellStats(data$grid_LargestPatch, 'sum') / raster::ncell(data$grid_LargestPatch)
		res[10] <- data$VegDensity_atMeanHullPosition
		res[11] <- mean(data$HullEdgeDistToOtherVeg_m, na.rm = TRUE)
		res[12] <- sum(data$HullEdgeDistToOtherVeg_m > 0, na.rm = TRUE) / width_N
		res[13] <- res[12] > 0
		
		tabulate_merge_into_etable(etable, index, res)
	}
	
	vegs <- c("Veg1", "Veg2")

	for (iveg in seq_along(vegs))
		etable <- place_veg(etable, index, data = data[[vegs[iveg]]], steplength, vegno = iveg, width_N)
	
	etable
}

#' @export
plot_Gastner2010_hulledge <- function(filename, eB_Env, eB_Veg, datFit) {
	pdf(width=7, height=10, file=filename)
	par_old <- par(mfrow=c(nr <- 2, 1), mar=c(0, 0.1, 1, 1), mgp=c(1.5, 0.5, 0), cex=cex <- 1.5)
	on.exit({par(par_old); dev.off()}, add = TRUE)


	#Panel a: map
	ext1 <- raster::extent(eB_Env$elev$grid)
	xlim <- c(-1000, ext1@xmax)
	ylim <- c(-1000, ext1@ymax)
	
	raster::image(eB_Env$elev$grid, col=gray(0:255/255), xlim=xlim, ylim=ylim, main="", xlab="", ylab="", asp=1, axes=FALSE)
	raster::image(eB_Veg$Veg1$grid, col=adjustcolor("red", alpha.f = 0.3), add=TRUE)
	raster::image(datFit$Veg1$grid_LargestPatch, col=adjustcolor("red", alpha.f = 0.3), add=TRUE)
	raster::image(eB_Veg$Veg2$grid, col=adjustcolor("darkgreen", alpha.f = 0.3), add=TRUE)
	raster::image(datFit$Veg2$grid_LargestPatch, col=adjustcolor("darkgreen", alpha.f = 0.3), add=TRUE)
	atx <- c((atx <- axTicks(1))[atx >= 0 & atx < ext1@xmax], ext1@xmax)
	axis(1, pos=0, at=atx)
	axis(2, pos=0, at=c(0, 2000, 4000, 6000))
	text(x=ext1@xmax/2, y=-strheight("0", units="user", cex=cex)*(0.5+2), labels="Transect length (m)", xpd=NA)
	text(x=-strwidth("0", units="user", cex=cex)*(0.5+2.5), y=ext1@ymax/2, labels="Transect width (m)", srt=90)

	lines(x=sp::coordinates(datFit$Veg1$spLine_hullEdge)[[1]][[1]], col="orange")
	lines(x=sp::coordinates(datFit$Veg2$spLine_hullEdge)[[1]][[1]], col="green")
	mtext(text="(a)", line=-1, cex=cex, adj=0.01)

	#Panel b: densities
	par(mar=c(2.5, 2.5, 1, 1.5))
	#Veg1
	plot(eB_Env$DistAlongXaxis_m, eB_Veg$Veg1$density, lty=1, col="red", type="l", xlim=c(0, ext1@xmax), ylim=c(0, 1), xlab="Transect length (m)", ylab="Density", xaxs="i", axes=FALSE)
	atx <- c((atx <- axTicks(1))[atx >= 0 & atx < ext1@xmax], ext1@xmax)
	axis(1, pos=0, at=atx)
	axis(2, pos=0)
	include <- datFit$Veg1$density > 0
	lines(eB_Env$DistAlongXaxis_m[include], datFit$Veg1$density[include], col="red", lwd=2)
	arrows(x0=pos <- datFit$Veg1$stats$position_m, y0=ytemp <- 1, x1=pos, y1=0, lwd=2, col="orange", length=0.03*par()$pin[1]/nr)
	arrows(x0=pos - datFit$Veg1$stats$width_m, y0=ytemp, x1=pos + datFit$Veg1$stats$width_m, y1=ytemp, lwd=2, col="orange", length=0.03*par()$pin[1]/nr, code=3)
	#Veg2
	lines(eB_Env$DistAlongXaxis_m, eB_Veg$Veg2$density, lty=1, col="darkgreen")
	include <- datFit$Veg2$density > 0
	lines(eB_Env$DistAlongXaxis_m[include], datFit$Veg2$density[include], col="darkgreen", lwd=2)
	arrows(x0=pos <- datFit$Veg2$stats$position_m, y0=ytemp-0.01, x1=pos, y1=0, lwd=2, col="green", length=0.03*par()$pin[1]/nr)
	arrows(x0=pos - datFit$Veg2$stats$width_m, y0=ytemp-0.01, x1=pos + datFit$Veg2$stats$width_m, y1=ytemp-0.01, lwd=2, col="green", length=0.03*par()$pin[1]/nr, code=3)

	mtext(text="(b)", line=-0.5, cex=cex, adj=0.01)

	invisible()
}




Gastner2010AmNat <- function(i, b, migtype, ecotoner_settings, etband, etmeasure, copy_FromMig1_TF, do_figures, ...) {
	#3c. Gastner et al. 2010 Am.Nat.: Location of boundary as hull edge
	
	dots <- list(...)
	seed <- if ("seed" %in% names(dots)) dots["seed"] else NULL
	if ("flag_bfig" %in% names(dots)) flag_bfig <- dots["flag_bfig"] else do_figures <- FALSE
	if ("dir_fig" %in% names(dots)) dir_fig <- dots["dir_fig"] else do_figures <- FALSE

	b_migtype <- (b - 1) * length(get("migtypes", envir = etr_vars)) + which(migtype == get("migtypes", envir = etr_vars))
	etmeasure$etable[b_migtype, "Transect_ID"] <- i
	etmeasure$etable[b_migtype, "Neighbor_Cells"] <- neighborhoods(ecotoner_settings)[b]
	etmeasure$etable[b_migtype, "Migration_Type"] <- migtype
	
	step_lengths <- stepsHullEdge(ecotoner_settings)
	etmeasure$gETmeas[[b]][[migtype]] <- vector(mode = "list", length = length(step_lengths))
	step_names <- paste0("Step", step_lengths)
	names(etmeasure$gETmeas[[b]][[migtype]]) <- step_names
	
	etmeasure$gETmeas[[b]][[migtype]]$meta <- list()
	etmeasure$gETmeas[[b]][[migtype]]$meta$method <- "Gastner2010AmNat"
	etmeasure$gETmeas[[b]][[migtype]]$meta$version <- version_Gastner2010AmNat()
	etmeasure$gETmeas[[b]][[migtype]]$meta$copy_FromMig1_TF <- copy_FromMig1_TF
	
	for (i_step in seq_along(step_names)) {
		istep <- step_names[i_step]
		if (!copy_FromMig1_TF) {
			for (iveg in c("Veg1", "Veg2")) {
				iother <- setdiff(c("Veg1", "Veg2"), iveg)
				type_veg <- if (iveg == "Veg1") type_veg1(ecotoner_settings) else if (iveg == "Veg2") type_veg2(ecotoner_settings) else NULL
				
				lvar <- calc_Gastner2010_hulledge(i,
						steplength = step_lengths[i_step],
						veg = etband$Veg[[migtype]][[iveg]]$grid,
						end_toLeft = end_to_left(type_veg))
				
				lvar$Hull_density <- count_y_at_each_x(lvar$grid_hullEdge)
				
				lvar$stats <- calc_Gastner2010_hulledge_Statistics(
						lvar$mat_hullEdge,
						lvar$spLine_hullEdge,
						res_m = raster::xres(etband$Env$elev$grid))
				
				lvar$VegDensity_atMeanHullPosition <- 
					etband$Veg[[migtype]][[iveg]]$density[round(lvar$stats$position_cells)]
				
				lvar$grad <- calc_Gastner2010_hulledge_Gradient(
						etband$Env$elev$grid,
						etband$Env$slope$grid,
						lvar$grid_hullEdge)
				
				lvar$HullEdgeDistToOtherVeg_m <- 
					calc_Gastner2010_hulledge_DistanceToStruggleZone(
						vegOther = etband$Veg[[migtype]][[iother]]$grid,
						mat_hullEdge = lvar$mat_hullEdge,
						end_toLeft = end_to_left(type_veg),
						width_N = bandTransect_width_cellN(ecotoner_settings))
						
				etmeasure$gETmeas[[b]][[migtype]][[istep]][[iveg]] <- lvar
			}


			if (do_figures) plot_Gastner2010_hulledge(
				filename = file.path(dir_fig, paste0(flag_bfig, "Gastner2009_FittedHullEdge_step", step_lengths[i_step], "_", migtype, ".pdf")),
				eB_Env = etband$Env,
				eB_Veg = etband$Veg[[migtype]],
				datFit = etmeasure$gETmeas[[b]][[migtype]][[istep]])
		} 
		
		etmeasure$etable <- tabulate_Gastner2010_hulledge(
			etable = etmeasure$etable, index = b_migtype,
			data = etmeasure$gETmeas[[b]][[if (copy_FromMig1_TF) "AllMigration" else migtype]][[istep]],
			steplength = step_lengths[i_step],
			width_N = bandTransect_width_cellN(ecotoner_settings))
	}
	
	etmeasure
}

