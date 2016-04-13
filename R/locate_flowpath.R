

calc.flowdir.uphill <- function(elev, seed = NULL) {
	#calc.flowdir.uphill returns the reverse 'flow direction' (of water), i.e.
	#	- the direction of the smallest rise in elevation
	#	- If two cells or more cells have the same smallest rise in elevation, a random cell is picked.
	#	- If all surrounding cells are lower, then top is reached and set to NA
	#
	#The are encoded as powers of 2 (0 to 7). The cell to the right of the focal cell 'x' is 1, the one below that is 2, and so on:
	#32	64	 128
	#16	 x	 1
	#8	 4	 2
	vals <- 2^matrix(c(5:7, 4, NA, 0, 3:1), nrow = 3, byrow = TRUE)
	
	uphill <- function(n3, ...){
		n3 <- matrix(n3, nrow = 3, ncol = 3, byrow = TRUE)
		res <- NA
		difs <- n3 - n3[2, 2]
		if (any(difs > 0, na.rm = TRUE)) { #there is a rise
			difs[difs <= 0] <- NA
			flowdir <- which(difs == min(difs, na.rm = TRUE), arr.ind = TRUE) #get the smallest rise(s)
			if (nrow(flowdir) > 1) {
				if (!is.na(seed)) set.seed(seed)
				flowdir <- flowdir[sample(nrow(flowdir), 1), ] #If two or more cells have the same smallest rise in elevation, a random cell is picked.
			}
			res <- vals[flowdir[1], flowdir[2]]
		}# else all directions are downslope or flat
		
		res
	}
	
	raster::focal(x = elev, w = matrix(1, nrow = 3, ncol = 3), fun = uphill, pad = TRUE, padValue = NA)
}

#Determines flow path from point
calc_flow_path <- function(flowdir, pathlimits = NULL, point) {
	#point = location to start flowpath calculation
	#pathlimits = raster with NAs and numbers >= 0; if flowpath reaches a cell with >= 0 then flowpath stops (limited_TF becomes TRUE and cell value is returned in limit_id)
	limited_TF <- FALSE
	limit_id <- NULL
	limitval <- NA
	res_m <- raster::xres(flowdir)
	
	#flowdir returns the downhill 'flow direction' (of water), i.e.
	#	- the direction of the greatest drop in elevation
	# 	- (or the smallest rise if all neighbors are higher).
	#	- If two cells have the same drop in elevation, a random cell is picked.
	#The are encoded as powers of 2 (0 to 7). The cell to the right of the focal cell 'x' is 1, the one below that is 2, and so on:
	#32	64	 128
	#16	 x	 1
	#8	 4	 2
	vals <- 2^(0:7)

	#Init
	fi <- 1
	projCRS <- raster::crs(flowdir)
	if ("SpatialPointsDataFrame" %in% (temp <- class(point))) {
		path <- sp::SpatialPoints(coords = point, proj4string = projCRS)
	} else if(identical("numeric", temp)){
		path <- sp::SpatialPoints(coords = data.frame(t(point)), proj4string = projCRS)
	} else if(identical("data.frame", temp)){
		path <- sp::SpatialPoints(coords = point[, 1:2], proj4string = projCRS)
	} else path <- point
	val_last <- val <- raster::extract(flowdir, path)
	if (!is.null(pathlimits)) limitval <- raster::extract(pathlimits, path)

	if (!is.na(limitval)) {
		limited_TF <- TRUE
		limit_id <- limitval
		path <- NULL
	} else if (!is.na(val)) repeat {	#Loop until out of raster area or a pathlimit is reached
		if (!(val %in% vals)) {#ArcGIS flow path adds up cell codes if sink and multiple options:
			#Get the different directions encoded
			tval <- val
			oval <- NULL
			iv <- length(vals)
			while (tval > 0 & iv > 0) {
				temp <- tval - vals[im]
				if (temp >= 0) {
					oval <- c(oval, vals[im])
					tval <- temp
					if (temp == 0) break
				}
				iv <- iv - 1
			}
			#Choose one direction: We attempt to follow the previous direction if possible, otherwise we choose a direction that is not opposite of the previous direction, otherwise we stop
			if (val_last %in% oval) {
				val <- val_last
			} else {
				vopp <- vals[(which(val_last == vals) + 0.5 * length(vals) - 1) %% length(vals) + 1]
				if (any(oval == vopp)) {
					oval <- oval[-which(vopp == oval)]
					if (length(oval) == 0) break
				}
				val <- oval[round(runif(1, min = 1, max = length(oval)))]
			}
		}
		
		#Calculate next cell based on value of flow direction; stop at raster boundary
		rc <- raster::rowColFromCell(flowdir, raster::cellFromXY(flowdir, path[fi, ]))
		newr <- rc[1] + switch(EXPR = paste0("d", val), d1 = 0, d2 = 1, d4 = 1, d8 = 1, d16 = 0, d32 = -1, d64 = -1, d128 = -1)
		newc <- rc[2] + switch(EXPR = paste0("d", val), d1 = 1, d2 = 1, d4 = 0, d8 = -1, d16 = -1, d32 = -1, d64 = 0, d128 = 1)
		if (!(raster::validCol(flowdir, newc) & raster::validRow(flowdir, newr))) break
		
		#Create next point
		newpoint <- raster::xyFromCell(flowdir, cell = raster::cellFromRowCol(flowdir, rownr = newr, colnr = newc), spatial = TRUE)
		#Get value of flow direction and pathlimits of next point
		val_last <- val
		val <- raster::extract(flowdir, newpoint)
		if (!is.null(pathlimits)) limitval <- raster::extract(pathlimits, newpoint)
		#Add next point to path if point is not already existing and valid and a limit is not reached
		if (nrow(sp::zerodist2(newpoint, path, zero = res_m/2)) > 0) break
		if (!is.na(limitval)) {
			limited_TF <- TRUE
			limit_id <- limitval
			break
		}
		path <- sp::rbind.SpatialPoints(path, sp::spTransform(newpoint, projCRS))
		fi <- fi + 1
		if (is.na(val)) break
	}
	
	list(path = path, limited_TF = limited_TF, limit_id = limit_id)
}


#' @export
calc.MigrationRoutes_EstimateFlowpaths <- function(elev, flowdir = NULL, patches, end_toLeft, seed = NULL) {
	#Categorize patches that are connected by a flowpath with the low x-border (good) or that connect with one of the transect sides (y-borders) (not good)
	#x-border = border where x = 0 or x = xmax (left or right side of transect)
	#y-border = border where y = 0 or y = ymax (bottom or top side)
	res_m <- raster::xres(elev)
	
	if(is.null(flowdir)) flowdir <- raster::terrain(elev, opt = "flowdir")
	#X and Y limits
	xlims <- c((temp <- raster::extent(elev))@xmin + 2 * res_m, temp@xmax - 2 * res_m)
	ylims <- c(temp@ymin + 2 * res_m, temp@ymax - 2 * res_m)
	#Inits
	patch_ids <- raster::unique(patches)
	minSample <- 3
	maxSample <- 6
	si <- 1 #Counter for flowpaths
	rpaths <- raster::raster(flowdir) #raster with established flowpaths (si): new flowpaths can connect and don't need to retrace already established downstream segments
	rpaths[] <- NA
	maxpathlength <- 10 ^ floor(1 + log(max(2 * ncol(flowdir), 2 * floor(sqrt(raster::ncell(flowdir)))), base = 10))
	#Result container
	tally <- matrix(0, nrow = length(patch_ids), ncol = 8, dimnames = list(NULL, c("PatchID", "GoodMigration_N", "BadMigration_N", "SuccessfulFlowpaths_N", "TestedFlowpaths_N", "PatchCoversGoodBorder_TF", "PatchCoversBadBorder_012", "PatchCells_N")))
		# tally[, "PatchCoversBadBorder_012"]: 0, does not reach y border; >= 1, does reach y border; 2, is completely within y border (i.e., 1 cell wide) and hence may not have a flowpath
	tally[, "PatchID"] <- patch_ids
	paths <- list()
	
	if (end_toLeft) {
		#We cannot search uphill for patch-connection; thus, we calculate flowpaths starting at each y-border cell and test whether any patch connects

		#Create points along lower and upper y-boundary
		xseq <- seq(from = xlims[1] - 1/2 * res_m, to = xlims[2] + 1/2 * res_m, by = res_m)
		yborderPoints <- sp::SpatialPoints(coord = rbind(cbind(xseq, ylims[1] - 1/2 * res_m), cbind(xseq, ylims[2] + 1/2 * res_m)), proj4str = raster::crs(elev))
		#Get flowpaths
		for (ip in seq_along(yborderPoints)) { #Loop through every y-border cell
			segment <- calc_flow_path(flowdir = flowdir, pathlimits = rpaths, point = yborderPoints[ip, ])
			if (!is.null(segment$path)) rpaths <- raster::overlay(rpaths, raster::rasterize(x = segment$path, y = rpaths, field = si * maxpathlength + seq_along(segment$path)), fun = function(x, y) ifelse(!is.na(y), y, x))
			if (segment$limited_TF) {
				connection_id <- segment$limit_id %/% maxpathlength
				connectingcell_id <- segment$limit_id - (connection_id * maxpathlength)
				addpath <- (temp <- paths[[connection_id]])[connectingcell_id:length(temp), ]
				if (!is.null(segment$path)) {
					path <- sp::rbind.SpatialPoints(segment$path, addpath)
				} else {
					path <- addpath
				}
			} else {
				path <- segment$path
			}
			paths[[si]] <- path
			si <- si + 1
		}
		paths_success_TF <- matrix(NA, nrow = length(paths), ncol = 3, dimnames = list(NULL, c("success", "out_good", "out_bad")))
		paths_success_TF[, "success"] <- TRUE
		
		#Calculate data for tally
		tally[, "TestedFlowpaths_N"] <- 1
		for (id in seq_along(patch_ids)) { #Loop through every patch
			i_vegPatch <- raster::calc(patches, fun = function(x) ifelse(!is.na(x) & x == patch_ids[id], 1, NA))
			tally[id, "PatchCells_N"] <- raster::cellStats(i_vegPatch, 'sum')
			tally[id, c("PatchCoversGoodBorder_TF", "GoodMigration_N", "SuccessfulFlowpaths_N")] <- as.numeric((extemp <- raster::extent(raster::trim(i_vegPatch)))@xmax >= xlims[2])
			tally[id, c("PatchCoversBadBorder_012", "BadMigration_N", "SuccessfulFlowpaths_N", "TestedFlowpaths_N")] <- as.numeric(any(extemp@ymin <= ylims[1]) || any(extemp@ymax >= ylims[2]))
			if (!tally[id, "PatchCoversGoodBorder_TF"]){ #Only trace flowpaths for patches that are not covering the upper x border
				i_yborderpath <- raster::mask(x = i_vegPatch, mask = rpaths)
				if (raster::cellStats(i_yborderpath, 'sum') > 0) { #this patch is connected with a y-border via a flowpath
					tally[id, "BadMigration_N"] <- tally[id, "SuccessfulFlowpaths_N"] <- 1
					ids_yborderpath <- raster::mask(x = rpaths, mask = i_vegPatch)
					ids_yborderpath <- unique(raster::unique(ids_yborderpath) %/% maxpathlength)
					paths_success_TF[ids_yborderpath, "out_bad"] <- TRUE
				} else if(raster::cellStats(raster::mask(x = i_vegPatch, mask = flowdir), 'sum') == 0 && (extemp@xmin > xlims[1] && extemp@xmax < xlims[2])){ #patch has cells only in y-border
					tally[id, "PatchCoversBadBorder_012"] <- 2
				} else {#not connected with any y-border
					tally[id, "GoodMigration_N"] <- tally[id, "SuccessfulFlowpaths_N"] <- 1
				}
			}
		}
	} else {
		paths_success_TF <- matrix(NA, nrow = maxSample * length(patch_ids), ncol = 3, dimnames = list(NULL, c("success", "out_good", "out_bad")))

		#We start in a random set of cells for each path and calculate flowpaths and test whether flowpaths leaves transect through a x- or y-border
		for (id in seq_along(patch_ids)) { #Loop through every patch
			i_vegPatch <- raster::calc(patches, fun = function(x) ifelse(!is.na(x) & x == patch_ids[id], 1, NA))
			tally[id, "PatchCells_N"] <- raster::cellStats(i_vegPatch, 'sum')
			tally[id, c("PatchCoversGoodBorder_TF", "GoodMigration_N", "SuccessfulFlowpaths_N", "TestedFlowpaths_N")] <- as.numeric((temp <- raster::extent(raster::trim(i_vegPatch)))@xmin <= xlims[1])
			tally[id, c("PatchCoversBadBorder_012", "BadMigration_N", "SuccessfulFlowpaths_N", "TestedFlowpaths_N")] <- as.numeric(any(temp@ymin <= ylims[1]) || any(temp@ymax >= ylims[2]))
			
			if (!tally[id, "PatchCoversGoodBorder_TF"]) { #Only trace flowpaths for patches that are not covering the lower x border
				i_vegPatch <- raster::mask(x = i_vegPatch, mask = flowdir)
				if ((ntemp <- raster::cellStats(i_vegPatch, 'sum')) > 0) { #test in case the entire patch lies on the y-border where flowdir cells are NA
					if (!is.na(seed)) set.seed(seed)
					i_vegPatch_sample <- raster::sampleRandom(i_vegPatch, size = min(ntemp, min(maxSample, max(minSample, round(0.05 * tally[id, "PatchCells_N"])))), sp = TRUE)
					if ((tally[id, "TestedFlowpaths_N"] <- length(i_vegPatch_sample)) > 0) {

						for (ip in seq_along(i_vegPatch_sample)) { #Loop through every sampled cell within patch
							segment <- calc_flow_path(flowdir = flowdir, pathlimits = rpaths, point = i_vegPatch_sample[ip, ])
							
							if (!is.null(segment$path)) rpaths <- raster::overlay(rpaths, raster::rasterize(x = segment$path, y = rpaths, field = si * maxpathlength + seq_along(segment$path)), fun = function(x, y) ifelse(!is.na(y), y, x))
							if (segment$limited_TF) {
								connection_id <- segment$limit_id %/% maxpathlength
								connectingcell_id <- segment$limit_id - (connection_id * maxpathlength)
								addpath <- (temp <- paths[[connection_id]])[connectingcell_id:length(temp), ]
								if (!is.null(segment$path)) {
									path <- sp::rbind.SpatialPoints(segment$path, addpath)
								} else {
									path <- addpath
								}
							} else {
								path <- segment$path
							}
							
							paths[[si]] <- path
							paths_success_TF[si, "out_good"] <- as.numeric(any((temp1 <- sp::coordinates(path[length(path)])[1]) <= xlims[1])) #only lower x-boundary is good
							paths_success_TF[si, "out_bad"] <- as.numeric(any(temp1 >= xlims[2]) || any((temp2 <- sp::coordinates(path[length(path)])[2]) <= ylims[1]) || any(temp2 >= ylims[2]))
							paths_success_TF[si, "success"] <- paths_success_TF[si, "out_good"] || paths_success_TF[si, "out_bad"]
							tally[id, "GoodMigration_N"] <- tally[id, "GoodMigration_N"] + paths_success_TF[si, "out_good"]
							tally[id, "BadMigration_N"] <- tally[id, "BadMigration_N"] + paths_success_TF[si, "out_bad"]
							tally[id, "SuccessfulFlowpaths_N"] <- tally[id, "SuccessfulFlowpaths_N"] + paths_success_TF[si, "success"]
							si <- si + 1
						}
					}
				} else {
					tally[id, "PatchCoversBadBorder_012"] <- 2
				}
			}
		}
		
		paths_success_TF <- paths_success_TF[seq_along(paths), ]
	}
	
	list(tally = tally, paths = paths, paths_success_TF = paths_success_TF)
}

#' @export
calc.Identify_GoodvsBadMigration <- function(patches4, tally, paths, paths_success_TF, end_toLeft) {
	patch_ids <- raster::unique(patches4)
	patchIDs_crossedByPaths <- lapply(paths, FUN = function(p) unique(na.omit(raster::extract(patches4, p)))) #entries in order of downstream patches4 crossed
	#patches4 that connect directly with the x-border
	patchIDs_BadMigration <- tally[tally[, "PatchCoversGoodBorder_TF"] == 1 | (tally[, "SuccessfulFlowpaths_N"] > 0 & tally[, "GoodMigration_N"] > 0 & tally[, "BadMigration_N"] == 0), "PatchID"]
#TODO(drs): shouldn't this be 'patchIDs_GoodMigration'??
		
	if (end_toLeft) {
		#Check network of flowpaths (originating from every y-border cell) and decide for each patch whether to be excluded or not
		#a. Do not exclude patches4 that are
		#	- on the upper x-border (PatchCoversGoodBorder_TF == 1)
		#	- or are not crossed by such a y-border cell flowpath
		#b. Exclude patches4 that are
		#	- not on the upper x-border and
		#	- directly downhill on a flowpath coming directly from a y-border cell,
		#	- or directly downhill of such a patch or chain of patches4, i.e., without crossing a x-patch

		#Helper functions: focus on target = first patch in flowpath
		get.targetPatchIDInPaths <- function(paths) unique(sapply(paths, FUN = function(p) ifelse(length(p) > 0, p[1], NA)))
		get.PathsWithTargetPatchID <- function(paths, targetID) sapply(paths, FUN = function(p) ifelse(length(p) > 0, targetID == p[1], FALSE))
		get.2ndPatchInPaths <- function(pathsDone, pathsTemp) unlist(lapply(pathsDone, FUN = function(ids) sapply(which(ids), FUN = function(id) ifelse(length(temp1 <- pathsTemp[id][[1]]) > 1, temp1[2], NA))))
		remove.targetPatchIDsFromPaths <- function(paths, removeIDs) lapply(paths, FUN = function(p) if(length(p) > 0 && p[1] %in% removeIDs) p[-1] else p)

	} else {
		#Check network of flowpaths (originating from every patch) and decide for each patch whether to be excluded or not
		#a. Do not exclude patches4 that are
		#	- on the lower x-border (PatchCoversGoodBorder_TF == 1)
		#	- or have flowpaths that do not leave transect directly on a y-border or in a y-only chain
		#		- (i.e.,flowpaths that end prematurely
		#		- or flowpaths that leave through a x-border
		#		- or flowpaths that connect to a x-patch)
		#b. Exclude patches4 that are
		#	- not on the lower x-border and
		#	- directly uphill of y-border even if patch contains flowpaths that leave transect through both x- and y-border,
		#	- or directly uphill of such a patch or chain of patches4, i.e., without crossing a x-patch

		#Helper functions: focus on target = last patch in flowpath
		get.targetPatchIDInPaths <- function(paths) unique(sapply(paths, FUN = function(p) ifelse(length(p) > 0, tail(p, n = 1), NA)))
		get.PathsWithTargetPatchID <- function(paths, targetID) sapply(paths, FUN = function(p) ifelse(length(p) > 0, targetID == tail(p, n = 1), FALSE))
		get.2ndPatchInPaths <- function(pathsDone, pathsTemp) unlist(lapply(pathsDone, FUN = function(ids) sapply(which(ids), FUN = function(id) ifelse((temp2 <- length(temp1 <- pathsTemp[id][[1]])) > 1, temp1[temp2-1], NA))))
		remove.targetPatchIDsFromPaths <- function(paths, removeIDs) lapply(paths, FUN = function(p) if(length(p) > 0 && tail(p, n = 1) %in% removeIDs) head(p, n = -1) else p)
	}

## TODO(drs): replace 'get.2ndPatchInPaths' with 'get.PatchPrevioustoIDInPaths'

	#patches4 that connect directly with a y-border
	lpatch0 <- tally[tally[, "PatchCoversBadBorder_012"] == 2, "PatchID"]
	lpatch1 <- tally[tally[, "PatchCoversGoodBorder_TF"] == 0 & (tally[, "PatchCoversBadBorder_012"] == 1 | (tally[, "SuccessfulFlowpaths_N"] > 0 & tally[, "PatchCoversGoodBorder_TF"] == 0 & tally[, "BadMigration_N"] > 0)), "PatchID"]
	lpatch2 <- get.targetPatchIDInPaths(patchIDs_crossedByPaths)
	lpatch3 <- lpatch1[lpatch1 %in% lpatch2]
	if (length(lpatch3) > 0) {
		patchIDs_GoodMigration <- get.targetPatchIDInPaths(patchIDs_crossedByPaths[apply(sapply(lpatch3, FUN = function(id) get.PathsWithTargetPatchID(patchIDs_crossedByPaths, id)), 1, any) & paths_success_TF[, "success"]])
	} else {
		patchIDs_GoodMigration <- NULL
	}
	patchIDs_GoodMigration <- unique(c(lpatch0, patchIDs_GoodMigration))
	#patches4 that connect in direct chain with a patch that connects directly with a y-border
	if (length(patchIDs_GoodMigration) > 0) {
		pathsTemp <- patchIDs_crossedByPaths
		repeat {
			lpatch1 <- lpatch2 <- get.targetPatchIDInPaths(pathsTemp)
			ilps <- lpatch1[lpatch1 %in% patchIDs_GoodMigration]
			if (length(ilps) > 0) {
				pathsYDone <- lapply(ilps, FUN = function(id) get.PathsWithTargetPatchID(pathsTemp, id))
				candidates_2ndPatchInPath_addPatchIDs <- get.2ndPatchInPaths(pathsYDone, pathsTemp)
				if (!is.null(candidates_2ndPatchInPath_addPatchIDs)) candidates_2ndPatchInPath_addPatchIDs <- na.exclude(candidates_2ndPatchInPath_addPatchIDs)
				if (length(candidates_2ndPatchInPath_addPatchIDs) > 0) {
					addPatchIDs <- candidates_2ndPatchInPath_addPatchIDs[!(candidates_2ndPatchInPath_addPatchIDs %in% patchIDs_BadMigration)]
					pathsTemp <- remove.targetPatchIDsFromPaths(pathsTemp, patchIDs_GoodMigration)
					patchIDs_GoodMigration <- unique(c(patchIDs_GoodMigration, addPatchIDs))
					lpatch1 <- lpatch1[!(lpatch1 %in% ilps)]
				}
			}
			if (identical(lpatch1, lpatch2)) break
		}
	} else {
		patchIDs_GoodMigration <- NULL
	}
	
	patchIDs_BadMigration <- patch_ids[!(patch_ids %in% patchIDs_GoodMigration)]
	
	list(patchIDs_BadMigration = patchIDs_BadMigration, patchIDs_GoodMigration = patchIDs_GoodMigration)
}

#' @export
calc.RemoveBadMigration_fromVeg <- function(patches4, patches8, patch4IDs_remove) {
	patched4 <- raster::calc(patches4, fun = function(x) ifelse(!is.na(x) & !(x %in% patch4IDs_remove), x, NA))
	patched8 <- raster::mask(patches8, mask = patched4)
	veg <- raster::calc(patched4, fun = function(x) ifelse(!is.na(x), 1, NA))
	
	list(grid = veg, patches8 = patched8, patches4 = patched4, patchID_removed = patch4IDs_remove)
}


