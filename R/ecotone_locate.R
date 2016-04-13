estimate_transect_homogeneity <- function(ecotoner_settings, ecotoner_grids, etband, temp_band4, aspC, aspL) {
	#2e2i. Quality/Homogeneity of environmental gradient along transect that is assumed in below methods
	#--Human footprint
	etband$Env$human$grid <- temp_band4$human
	etband$Env$human$density <- count_y_at_each_x(etband$Env$human$grid)
	if (sum(etband$Env$human$density) > 0) {
		humanClump8 <- raster::clump(etband$Env$human$grid, directions=8)
		i.largest_hC8 <- (temp <- (temp <- raster::freq(humanClump8))[complete.cases(temp), , drop=FALSE])[which.max(temp[, 2]), 1]
		largest_hC8 <- raster::trim(raster::calc(humanClump8, fun=function(x) ifelse(x == i.largest_hC8, 1, NA)))
		etband$Env$human$LargestClump_XExtent <- raster::xmax(largest_hC8) - raster::xmin(largest_hC8)
		etband$Env$human$LargestClump_YExtent <- raster::ymax(largest_hC8) - raster::ymin(largest_hC8)
		# rm(humanClump8, i.largest_hC8, largest_hC8)
	} else {
		etband$Env$human$LargestClump_XExtent <- etband$Env$human$LargestClump_YExtent <- 0
	}
	#--Quality of aspect: local (RMSE), patches (clumps), and global (mean) aspect (at slopes more than a minimum) should be about equal to transect gradient = get("transect_azimuth", envir = etr_vars) (from left to right)
	#'local' error of aspect
	resids <- circ_sum(aspC-transect_azimuth(ecotoner_settings), pi, int=2*pi) - pi
	etband$Env$aspect$Aspect_Cellwise_RMSE <- sqrt(circ_mean(resids^2, int=2*pi, na.rm = TRUE))
	etband$Env$aspect$Aspect_Cellwise_95PercQuantileAbsError <- quantile_95th_of_abs(resids, na.rm = TRUE)
	#'y-means' error of aspect
	resids <- circ_sum(etband$Env$aspect$YMeans_ForEachX-transect_azimuth(ecotoner_settings), pi, int=2*pi) - pi
	etband$Env$aspect$Aspect_YMeans_RMSE <- sqrt(circ_mean(resids^2, int=2*pi, na.rm = TRUE))
	etband$Env$aspect$Aspect_YMeans_95PercQuantileAbsError <- quantile_95th_of_abs(resids, na.rm = TRUE)
	#'patches' of aspects grouped into the four cardinal directions
	aspD <- raster::calc(aspL, fun=function(x) 1 + ((x + pi/4) %% (2*pi)) %/% (pi/2)) #North (325-45 deg, 1; 45-135 deg, 2; 135-225 deg, 3; 225-315 deg, 4); transect aspect is 270 deg
	aspC1 <- raster::clump(raster::calc(aspD, fun=function(x) ifelse(x==1, 1, 0)), directions=4)
	aspC2 <- raster::clump(raster::calc(aspD, fun=function(x) ifelse(x==2, 1, 0)), directions=4)
	aspC3 <- raster::clump(raster::calc(aspD, fun=function(x) ifelse(x==3, 1, 0)), directions=4)
	ncells <- raster::ncell(aspD)
	etband$Env$aspect$Aspect_Patch_FractionLargestDeviatingPatch <- max((temp1 <- raster::freq(aspC1))[!is.na(temp1[, 1]), 2], (temp2 <- raster::freq(aspC3))[!is.na(temp2[, 1]), 2], (temp3 <- raster::freq(aspC3))[!is.na(temp3[, 1]), 2]) / ncells #fraction of band transect that is largest patch of aspects with deviating cardinal directions
	etband$Env$aspect$Aspect_Patch_FractionAllDeviatingPatches <- (sum(temp1[!is.na(temp1[, 1]), 2]) + sum(temp2[!is.na(temp2[, 1]), 2]) + sum(temp3[!is.na(temp3[, 1]), 2])) / ncells #fraction of band transect that is made up by aspects of deviating cardinal directions
	#'overall' error of aspect
	etband$Env$aspect$Aspect_Overall_ME <- circ_mean(aspC, int=2*pi, na.rm = TRUE) - transect_azimuth(ecotoner_settings)
	# rm(aspL, aspC, aspD, aspC1, aspC2, aspC3)

	#--Quality of elevation:
	telev <- raster::as.matrix(etband$Env$elev$grid, maxpixels = raster::ncell(etband$Env$elev$grid))
	#'y/x ratios' of elevation
	yratio <- apply(telev, 2, FUN = function(x) max(x) - min(x))[-c(1, ncols <- ncol(telev))] / (tan(etband$Env$slope$YMeans_ForEachX[-c(1, ncols)])*bandTransect_width_cellN(ecotoner_settings)*res_m(specs_grid(ecotoner_grids))) #no slope in first and last column, because slope unidentifiable
	etband$Env$elev$Elevation_RatioYrangeToExpectXrange_RMSE <- root_mean_square(yratio, na.rm = TRUE)
	etband$Env$elev$Elevation_RatioYrangeToExpectXrange_95PercQuantileRatioAbsError <- quantile_95th_of_abs(yratio, na.rm = TRUE)
	etband$Env$elev$Elevation_RatioYrangeToExpectXrange_MaxAbsError <- max_of_abs(yratio, na.rm = TRUE)

	#'local' elevational error from simple linear model transect elevation gradient
	resids1 <- rep(etband$Env$elev$Ymeans_lm, each=nrow(telev)) - telev
	etband$Env$elev$Elevation_Cellwise_lmResiduals_RMSE <- root_mean_square(resids1, na.rm = TRUE)
	etband$Env$elev$Elevation_Cellwise_lmResiduals_95PercQuantileAbsError <- quantile_95th_of_abs(resids1, na.rm = TRUE)
	etband$Env$elev$Elevation_Cellwise_lmResiduals_MaxAbsError <- max_of_abs(resids1, na.rm = TRUE)
	#'local' elevational error from y-means transect elevation gradient
	resids2 <- rep(etband$Env$elev$YMeans_ForEachX, each=nrow(telev)) - telev
	etband$Env$elev$Elevation_Cellwise_YMeans_RMSE <- root_mean_square(resids2, na.rm = TRUE)
	etband$Env$elev$Elevation_Cellwise_YMeans_95PercQuantileAbsError <- quantile_95th_of_abs(resids2, na.rm = TRUE)
	etband$Env$elev$Elevation_Cellwise_YMeans_MaxAbsError <- max_of_abs(resids2, na.rm = TRUE)
	
	#global spatial autocorrelation in elevation residuals from y-means transect elevation gradient
	mresids <- sp::SpatialPointsDataFrame(coords = sp::coordinates(etband$Env$elev$grid),
											data = data.frame(layer = as.vector(t(resids2))), proj4string = raster::crs(etband$Env$elev$grid))
	sp::gridded(mresids) <- TRUE #upgrade to SpatialPixelDataFrame
	if (anyNA(resids2)) mresids <- mresids[!is.na(mresids$layer), ] # Remove NAs for 'moran.mc' and 'calc_variogram_range'
	
	nbw <- spdep::listw2U(spdep::nb2listw(spdep::knn2nb(spdep::knearneigh(sp::coordinates(mresids), k=8, longlat=FALSE)), style="W")) #neighbors: Queen's case; weights: row standardized ("Row standardization is recommended whenever the distribution of your features is potentially biased due to sampling design or an imposed aggregation scheme"); listw2U(): makes weight matrix symmetric
	mimc <- spdep::moran.mc(x = mresids$layer, listw = nbw, alternative = "greater", nsim = 199) #Null hypothesis: Moran's I = 0 = = no spatial autocorrelation (alternative default >0); spdep::moran.mc must be instructed how to deal with NAs
	etband$Env$elev$Elevation_Cellwise_YMeans_MoransI <- mimc$statistic
	etband$Env$elev$Elevation_Cellwise_YMeans_MoransIp <- mimc$p.value			
	
	#Gastner et al. 2009: we hypothesize that heterogeneity in the environmental background would not influence the model’s predictions unless the sizes of vegetation patches are smaller than the correlation length of the heterogeneity			
	etband$Env$elev$Elevation_YMeansResiduals_VariogramRange <- calc_variogram_range(mresids, max(etband$Env$DistAlongXaxis_m), res_m = res_m(specs_grid(ecotoner_grids))) #Spatial autocorrelation distance: range
	
	etband
}


#' Locate ecotone, prepare ecotone line and band transect, extract relevant information
#'
#' @export
establish_ecotone_transect <- function(x) x
#TODO(drs): devtools::document() doesn't like to export the real function below; this is a hack

establish_ecotone_transect <- function(i, b, etband, etable, ipoint, ecotoner_settings, ecotoner_grids, elevCropped, gapCropped, bseABUTtfCropped, asp201MeanCropped, asp201SDCropped, tempData = NULL, dir_fig, figBasename, verbose, do_figures, seed) {
	if	(verbose) {
		cat("'ecotoner' establishing: tr = ", i, "; neigh = ", b, "; start at ", format(Sys.time(), format = ""), "\n", sep = "")
		idh <- 0 #counter for debug 'here' statements
	}
	
## TODO(drs): assumptions used about 'gap.rat'; replace properly with object of class 'TypeInfo'
gap.rat <- df_veg(ecotoner_grids)
	
	#Result containers for GIS and similar data
	# "search-transect-line = temp_stline" vs. "ecotone-transect-line = temp_etline"
	# "search-transect-band = temp_stband" vs. "ecotone-transect-band = temp_etband"
	temp_stline <- temp_etline <- list()	#elevation (start), linear, and band transects
	status <- "searching"

	template_eBTable <- data.frame(matrix(NA, nrow=1, ncol=0))
	template_eBTable$Transect_ID <- i
	template_eBTable$Neighbor_Cells <- neighborhoods(ecotoner_settings)[b]

	
	#---1. Locate neighborhood and elevation gradient: define transect across elevation

	#---2. Locate candidate search-transect lines					
	if (!is.null(tempData) && file.exists(tempData)) {
		temp_stline_cands <- readRDS(tempData)
	} else {
		temp_stline_cands <- locate_candidate_transect_lines(transect_type = transect_type(ecotoner_settings),
								start = ipoint,
								neighbors = neighborhoods(ecotoner_settings)[b],
								max_neighborhood = max(neighborhoods(ecotoner_settings)),
								grid_gradient = elevCropped,
								grid_flow = grid_flow(ecotoner_grids),
								asp201MeanCropped = asp201MeanCropped,
								asp201SDCropped = asp201SDCropped,
								bseABUTtfCropped = bseABUTtfCropped,
								candidate_THA_N = candidate_THAs_N(ecotoner_settings),
								width_N = bandTransect_width_cellN(ecotoner_settings),
								max_aspect_sd = Aspect_SDof201Mean_Maximum(ecotoner_settings),
								max_aspect_diff = Aspect_DeviationAlongCandidateLines_Maximum(ecotoner_settings),
								seed = seed)
		if (!is.null(tempData)) saveRDS(temp_stline_cands, file = tempData)
	}
	
	if (temp_stline_cands$line_N > 0) {	
		if(verbose) cat("'ecotoner' establishing: tr = ", i, "; neigh = ", b, "; prog: ", idh <- idh + 1, "; # candidates = ", temp_stline_cands$line_N, "\n", sep = "")

		transect_success <- FALSE
		for (ic in seq_len(temp_stline_cands$line_N)) {
			
			# "search-transect-line = temp_stline" vs. "ecotone-transect-line = temp_etline"
			temp_stline <- orient_transect_line(pts_tcand = temp_stline_cands$lines[[ic]], longlat = longlat(specs_grid(ecotoner_grids)))
		
			if (length(temp_stline) == 0) {
				if(verbose) cat("'ecotoner' establishing: tr = ", i, "; neigh = ", b, "; prog: ", idh, "; cand = ", ic,
								": no transect located because transect line has 0 length after orienting along gradient and ordering by distance\n", sep = "")
				next #goto next candidate transect
			}
		
			#Store elevation transect data for analysis
			temp.etable <- get.ElevationTransect_tempTableData(temp_Table = template_eBTable, data = temp_stline)
		

			#2b. Extract band transect along elevation transect with a width of 200 cells
			#Calculate outline of maximal band transect along elevation transect
			# and how to transform from original coordinates to band transect coordinates
		
			temp_stband <- locate_candidate_transect_band(stline = temp_stline,
													width_N = bandTransect_width_cellN(ecotoner_settings),
													res_rotate_grid = res_rotate_grid(ecotoner_settings),
													elevCropped = elevCropped,
													elev = grid_env(ecotoner_grids), gap = grid_veg(ecotoner_grids),
													bseRatValue = type_ids(type_veg1(ecotoner_settings)),
													tfRatValue = type_ids(type_veg2(ecotoner_settings)),
													seed = seed)
		
			if (is.null(temp_stband$grids)) {
				if(verbose) cat("'ecotoner' establishing: tr = ", i, "; neigh = ", b, "; prog: ", idh, "; cand = ", ic,
								": no transect located because band transect falls outside neighborhood window\n", sep = "")			
				next # goto next candidate transect
			}				

			#2c. Apply definition of zone of ecological boundary to BSE and TF
			limits1_ZoneEcolBoundary <- determine_ecotone_linear_extent(dens_BSE=NULL, dens_Veg1=temp_stband$grids$bse_dens, end1_toLeft=end_to_left(type_veg1(ecotoner_settings)),
																		dens_TF=NULL, dens_Veg2=temp_stband$grids$tf_dens, end2_toLeft=end_to_left(type_veg2(ecotoner_settings)),
																		veg_density_low = vegDefiningDensityTransect_low(ecotoner_settings),
																		veg_density_high = vegDefiningDensityTransect_high(ecotoner_settings),
																		veg_density_extended_min = vegDefiningDensityTransectExtended_min(ecotoner_settings))
		
		
			if (anyNA(limits1_ZoneEcolBoundary)) {
				if(verbose) cat("'ecotoner' establishing: tr = ", i, "; neigh = ", b, "; prog: ", idh, "; cand = ", ic,
								": no transect located because there is no ecological boundary zone (step 1)\n", sep = "")			
				next #goto next candidate transect; there is no zone of the ecological boundary according to definition
			} else {
				if(verbose) cat("'ecotoner' establishing: tr = ", i, "; neigh = ", b, "; prog: ", idh, "; cand = ", ic,
								"; limits1 = ", paste(limits1_ZoneEcolBoundary, collapse=", "), "\n", sep = "")			
			}

			#Crop band transect to limits1_ZoneEcolBoundary
			temp_etband <- trim1_band_to_ecotone(stband_grids = temp_stband$grids,
													start_icol = limits1_ZoneEcolBoundary[1],
													end_icol = limits1_ZoneEcolBoundary[2])
			temp_etband$DistAlongXaxis_m <- res_m(specs_grid(ecotoner_grids)) * ((1:dim(temp_etband$bse)[2]) - 1/2)

			#Determine majority BSE and majority TF types within limits1_ZoneEcolBoundary
# TODO(drs): move to function; maybe create sub-lists for each veg type?
			#majority BSE
			bseFreq <- raster::freq(temp_etband$bse)
			temp_etband$ratValue_BSE <- bseFreq[bseFreq[, 1] %in% type_ids(type_veg1(ecotoner_settings)), 1]
			if (length(temp_etband$ratValue_BSE) == 0) {
				if(verbose) cat("'ecotoner' establishing: tr = ", i, "; neigh = ", b, "; prog: ", idh, "; cand = ", ic,
								": no transect located because there is no Veg1 within ecological boundary zone\n", sep = "")			
				break #no BSE in transect
			}
			temp_etband$level3_densities_BSE <- (dtemp <- bseFreq[itemp <- match(temp_etband$ratValue_BSE, bseFreq[, 1], nomatch = 0), 2]) / raster::ncell(temp_etband$bse)
			temp_etband$majorityVal_BSE <- bseFreq[itemp, 1][which.max(dtemp)]
			temp_etband$majorityType_BSE <- paste(sapply(gap.rat[gap.rat$Value == temp_etband$majorityVal_BSE, c("FRM", "DIV", "MACRO_CD", "LEVEL3", "ECOLSYS_LU")], FUN = as.character), collapse = "_")
			temp_etband$Veg1 <- raster::calc(temp_etband$bse, fun = function(x) ifelse(x == temp_etband$majorityVal_BSE, 1, NA))
			temp_etband$Veg1dens <- count_y_at_each_x(temp_etband$Veg1)
		
			#majority TF
			tfFreq <- raster::freq(temp_etband$tf)
			temp_etband$ratValue_TF <- tfFreq[tfFreq[, 1] %in% type_ids(type_veg2(ecotoner_settings)), 1]
			if (length(temp_etband$ratValue_TF) == 0) {
				if(verbose) cat("'ecotoner' establishing: tr = ", i, "; neigh = ", b, "; prog: ", idh, "; cand = ", ic,
								": no transect located because there is no Veg2 within ecological boundary zone\n", sep = "")			
				break #no TF in transect
			}
			temp_etband$level3_densities_TF <- tfFreq[match(temp_etband$ratValue_TF, tfFreq[, 1], nomatch = 0), 2] / raster::ncell(temp_etband$tf)
			temp_etband$macrocd_TF <- unique(temp <- gap.rat$MACRO_CD[match(temp_etband$ratValue_TF, gap.rat$Value, nomatch = 0)])
			temp_etband$macrocd_densities_TF <- aggregate(temp_etband$level3_densities_TF, by = list(temp), FUN = sum)[, 2]
			temp_etband$majorityMACROCD_TF <- temp_etband$macrocd_TF[which.max(temp_etband$macrocd_densities_TF)]
			temp_etband$majorityVals_TF <- temp_etband$ratValue_TF[temp == temp_etband$majorityMACROCD_TF]
			temp_etband$majorityType_TF <- paste(sapply(gap.rat[which(gap.rat$MACRO_CD == temp_etband$majorityMACROCD_TF)[1], c("FRM", "DIV", "MACRO_CD", "NVC_MACRO")], FUN = as.character), collapse = "_")
			temp_etband$Veg2 <- raster::calc(temp_etband$tf, fun = function(x) ifelse(x %in% temp_etband$majorityVals_TF, 1, NA))
			temp_etband$Veg2dens <- count_y_at_each_x(temp_etband$Veg2)

			#2c. Apply definition of zone of ecological boundary to majority of BSE and majority of TF within limits1_ZoneEcolBoundary and make sure there is not too much of remaining BSE or TF types around
			limits2_ZoneEcolBoundary <- determine_ecotone_linear_extent(dens_BSE=temp_etband$bse_dens, dens_Veg1=temp_etband$Veg1dens, end1_toLeft=end_to_left(type_veg1(ecotoner_settings)),
																		dens_TF=temp_etband$tf_dens, dens_Veg2=temp_etband$Veg2dens, end2_toLeft=end_to_left(type_veg2(ecotoner_settings)),
																		veg_density_low = vegDefiningDensityTransect_low(ecotoner_settings),
																		veg_density_high = vegDefiningDensityTransect_high(ecotoner_settings),
																		veg_density_extended_min = vegDefiningDensityTransectExtended_min(ecotoner_settings))

			if (anyNA(limits2_ZoneEcolBoundary)) {
				if(verbose) cat("'ecotoner' establishing: tr = ", i, "; neigh = ", b, "; prog: ", idh, "; cand = ", ic,
								": no transect located because there is no ecological boundary zone (step 2)\n", sep = "")			
				next #goto next candidate transect; there is no zone of the ecological boundary according to definition
			} else {
				if(verbose) cat("'ecotoner' establishing: tr = ", i, "; neigh = ", b, "; prog: ", idh, "; cand = ", ic,
								"; limits2 = ", paste(limits2_ZoneEcolBoundary, collapse=", "), "\n", sep = "")			
				transect_success <- TRUE
				break #leave loop through candidate transects
			}
		} #end of loop through candidate transects


		if (transect_success) {
			#end of testing
			limits2_ZoneEcolBoundary <- limits1_ZoneEcolBoundary[1] + limits2_ZoneEcolBoundary - 1

			#Plot map of density profiles
			if(do_figures) plot_transect_density_profiles(filename=file.path(dir_fig, paste0(figBasename, "DensityProfiles.pdf")),
													distX1=distX1 <- temp_stband$grids$DistAlongXaxis_m,
													distX2=(limits1_ZoneEcolBoundary[1] - 1) * res_m(specs_grid(ecotoner_grids)) + temp_etband$DistAlongXaxis_m,
													elev=mean_y_at_each_x(temp_stband$grids$grid1),
													dens_Veg1=temp_etband$Veg1dens, dens_BSE=temp_stband$grids$bse_dens,
													dens_Veg2=temp_etband$Veg2dens, dens_TF=temp_stband$grids$tf_dens,
													start1_m=distX1[limits1_ZoneEcolBoundary[1]], end1_m=distX1[limits1_ZoneEcolBoundary[2]],
													start2_m=distX1[limits2_ZoneEcolBoundary[1]], end2_m=distX1[limits2_ZoneEcolBoundary[2]])
	
			if(verbose) cat("'ecotoner' establishing: tr = ", i, "; neigh = ", b, "; prog: ", idh <- idh + 1, "\n", sep = "")

			#---Above was the last test about whether to accept a transect
			#Now copy temporary data table into final variable
			etable <- get.CopyTempTo_TableData(temp_Table=temp.etable, etable=etable, b=b)


			#2d. Crop the linear transect to the ecological boundary zone
			temp_etline <- trim_line_to_ecotone(stline_pts = temp_stline$pts,
												distX = temp_stline$dist_m,
												start_m = temp_stband$grids$DistAlongXaxis_m[limits2_ZoneEcolBoundary[1]],
												end_m = temp_stband$grids$DistAlongXaxis_m[limits2_ZoneEcolBoundary[2]],
												crs = crs(specs_grid(ecotoner_grids)),
												longlat = longlat(specs_grid(ecotoner_grids)))

			#Store linear transect data for analysis
			etable <- get.LinearTransect_TableData(etable=etable, b=b, data=temp_etline)
			if(verbose) cat("'ecotoner' establishing: tr = ", i, "; neigh = ", b, "; prog: ", idh <- idh + 1, "\n", sep = "")


			#2e. Crop the band transect to the ecological boundary zone
			#'Global' patch size distribution: calculate patches on gapCropped extent: hoping that this is as close as possible to global patches for the band transect extent
			Veg1Cropped <- raster::calc(gapCropped, fun = function(x) ifelse(x %in% temp_etband$majorityVal_BSE, 1, NA))
			Veg2Cropped <- raster::calc(gapCropped, fun = function(x) ifelse(x %in% temp_etband$majorityVals_TF, 1, NA))
			temp.patches1 <- list()
			temp.patches1$Veg1_patches8 <- raster::clump(Veg1Cropped, directions = 8)
			temp.patches1$freq_patches8_Veg1 <- data.frame(raster::freq(temp.patches1$Veg1_patches8, useNA = "no"))
			temp.patches1$Veg2_patches8 <- raster::clump(Veg2Cropped, directions = 8)
			temp.patches1$freq_patches8_Veg2 <- data.frame(raster::freq(temp.patches1$Veg2_patches8, useNA = "no"))
			temp.patches2 <- extract_tband_grids(tpoly_UTM = temp_stband$polygon$tpoly_UTM,
												width_n = bandTransect_width_cellN(ecotoner_settings),
												crs_UTM = temp_stband$polygon$crs_UTM,
												declin = temp_stband$polygon$declin,
												center = temp_stband$polygon$center,
												res_rotate_grid = res_rotate_grid(ecotoner_settings),
												grid1 = temp.patches1$Veg1_patches8,
												grid2 = temp.patches1$Veg2_patches8,
												fun1 = function(x, ...) majority(x, seed = seed),
												fun2 = function(x, ...) majority(x, seed = seed))

			#Human impacts
			HumanCropped <- raster::calc(gapCropped, fun = function(x) ifelse(x %in% type_ids(type_excl(ecotoner_settings)), 1, NA))
			temp.human <- extract_tband_grids(tpoly_UTM = temp_stband$polygon$tpoly_UTM,
												width_n = bandTransect_width_cellN(ecotoner_settings),
												crs_UTM = temp_stband$polygon$crs_UTM,
												declin = temp_stband$polygon$declin,
												center = temp_stband$polygon$center,
												res_rotate_grid = res_rotate_grid(ecotoner_settings),
												grid1 = HumanCropped,
												fun1 = function(x, ...) majority(x, seed = seed))

			#2e1. Crop band transect to zone of ecol boundary and set start of 'my' coordinate system
			temp_band4 <- trim2_band_to_ecotone(tempBand1 = temp_stband$polygon,
													tempBand2 = temp_stband$grids,
													tempBand3 = temp_etband,
													tempPatches = temp.patches2,
													tempHuman = temp.human$grid1,
													start1_icol = limits1_ZoneEcolBoundary[1],
													start2_icol = limits2_ZoneEcolBoundary[1],
													end2_icol = limits2_ZoneEcolBoundary[2])


			#2e2. Environmental gradient along transect
			#Copy temporary environmental data into final variable
			etband$Env$DistAlongXaxis_m <- res_m(specs_grid(ecotoner_grids)) * ((1:dim(temp_band4$elev)[2]) - 1/2)
			etband$Env$elev$grid <- temp_band4$elev
			etband$Env$band.pts <- temp_band4$band.pts
			etband$Env$declin <- temp_band4$declin

			#Extract environmental conditions along transect band
			etband$Env$slope$grid <- raster::terrain(etband$Env$elev$grid, opt = "slope", neighbors = 8)  #Holt 1981 algorithm
			etband$Env$slope$YMeans_ForEachX <- mean_y_at_each_x(grid = etband$Env$slope$grid)
			etband$Env$elev$YMeans_ForEachX <- mean_y_at_each_x(grid = etband$Env$elev$grid)
			etband$Env$elev$YSDs_ForEachX <- sd_y_at_each_x(grid = etband$Env$elev$grid)
			etband$Env$elev$Ymeans_lm <- predict(lm(etband$Env$elev$YMeans_ForEachX ~ etband$Env$DistAlongXaxis_m))
			etband$Env$aspect$grid <- raster::terrain(etband$Env$elev$grid, opt = "aspect", neighbors = 8) #it is faster to re-calculate slope and aspect than to load it from the master file and re-sample it to our band transect
			aspL <- raster::overlay(etband$Env$aspect$grid, etband$Env$slope$grid, fun = function(a, sl) ifelse(abs(sl) > min_slope_with_aspect(ecotoner_settings), a, ifelse(sl > 0, transect_azimuth(ecotoner_settings), 2*pi - transect_azimuth(ecotoner_settings))))
			aspC <- raster::as.matrix(aspL, maxpixels = raster::ncell(aspL))
			etband$Env$aspect$YMeans_ForEachX <- apply(aspC, 2, FUN = circ_mean, int = 2*pi, na.rm = TRUE)
			etband$Env$aspect$YSDs_ForEachX <- apply(aspC, 2, FUN = circ_sd, int = 2*pi, na.rm = TRUE)
	
			if(verbose) cat("'ecotoner' establishing: tr = ", i, "; neigh = ", b, "; prog: ", idh <- idh + 1, "\n", sep = "")
			
			# Quality/Homogeneity of environmental gradient along transect
			etband <- estimate_transect_homogeneity(ecotoner_settings, ecotoner_grids, etband, temp_band4, aspC, aspL)

			#Plot map of linear and band transect
			if(do_figures) map_transect(filename = file.path(dir_fig, paste0(figBasename, "TransectMap_v1.pdf")),
									stline_pts = temp_stline$pts, etline_pts = temp_etline$pts, gB_Env = etband$Env, 
									elevCropped = elevCropped,
									bse = bseCropped <- raster::calc(gapCropped, fun = function(x) ifelse(x %in% type_ids(type_veg1(ecotoner_settings)), x, NA)),
									tf = tfCropped <- raster::calc(gapCropped, fun = function(x) ifelse(x %in% type_ids(type_veg2(ecotoner_settings)), x, NA)),
									Veg1 = Veg1Cropped, Veg2 = Veg2Cropped,
									human = HumanCropped)
			if(do_figures) map_transect(filename = file.path(dir_fig, paste0(figBasename, "TransectMap_v2.pdf")),
									stline_pts = temp_stline$pts, etline_pts = temp_etline$pts, gB_Env = etband$Env, 
									elevCropped = raster::crop(grid_env(ecotoner_grids), temp_stband$ext),
									bse = raster::crop(bseCropped, temp_stband$ext), tf = raster::crop(tfCropped, temp_stband$ext),
									Veg1 = raster::crop(Veg1Cropped, temp_stband$ext), Veg2 = raster::crop(Veg2Cropped, temp_stband$ext),
									human = raster::crop(HumanCropped, temp_stband$ext))
			#rm(Veg1Cropped, Veg2Cropped, HumanCropped, bseCropped, tfCropped)
	
			#Store environmental band transect data for analysis
			etable <- get.BandTransect_TableDataEnv(etable = etable, b = b, data = etband$Env)

			if(verbose) cat("'ecotoner' establishing: tr = ", i, "; neigh = ", b, "; prog: ", idh <- idh + 1, "\n", sep = "")

			#2e3. Vegetation data along transect
			#Copy temporary vegetation data into final variable:
			#	- Veg$AllMigration: all vegetation in transect
			#	- Veg$OnlyGoodMigration: all vegetation minus patches whose flowpath connects to beyond y-border indicating non-x-axis migration origin

			#2e3i. Veg$AllMigration
			#Need to recalculate veg1, veg2, veg1density, and veg2density based on patches because identical([clumps(extract_tband_grids(gap) %in% valList)], [extract_tband_grids(clumps(gap %in% valList))]) is likely FALSE
			etband$Veg$AllMigration$Veg1$patches8 <- temp_band4$Veg1patch8
			etband$Veg$AllMigration$Veg1$patches8_freq <- temp.patches1$freq_patches8_Veg1
			etband$Veg$AllMigration$Veg1$grid <- raster::calc(etband$Veg$AllMigration$Veg1$patches8, fun = function(x) ifelse(!is.na(x), 1, NA)) #cell values: NA, no vegetation; 1, vegetation (edge or interior)
			etband$Veg$AllMigration$Veg1$density <- count_y_at_each_x(etband$Veg$AllMigration$Veg1$grid)
			etband$Veg$AllMigration$Veg1$patches4 <- raster::clump(etband$Veg$AllMigration$Veg1$grid, directions = 4)	#re-calc clumps (now Rook's case) because rotation of band transect may distort detected patches8
			etband$Veg$AllMigration$Veg1$ratValue_majority <- temp_etband$majorityVal_BSE
			etband$Veg$AllMigration$Veg1$type_majority <- temp_etband$majorityType_BSE
			etband$Veg$AllMigration$Veg2$patches8 <- temp_band4$Veg2patch8
			etband$Veg$AllMigration$Veg2$patches8_freq <- temp.patches1$freq_patches8_Veg2
			etband$Veg$AllMigration$Veg2$grid <- raster::calc(etband$Veg$AllMigration$Veg2$patches8, fun = function(x) ifelse(!is.na(x), 1, NA))
			etband$Veg$AllMigration$Veg2$density <- count_y_at_each_x(etband$Veg$AllMigration$Veg2$grid)
			etband$Veg$AllMigration$Veg2$patches4 <- raster::clump(etband$Veg$AllMigration$Veg2$grid, directions = 4)	#re-calc clumps (now Rook's case) because rotation of band transect may distort detected patches8
			etband$Veg$AllMigration$Veg2$ratValues_majority <- temp_etband$majorityVals_TF
			etband$Veg$AllMigration$Veg2$type_majority <- temp_etband$majorityType_TF

			#Recalculate majority BSE and majority TF type
			#majority BSE
			bseFreq <- raster::freq(temp_band4$bse)
			etband$Veg$AllMigration$BSE$ratValues <- bseFreq[bseFreq[, 1] %in% type_ids(type_veg1(ecotoner_settings)), 1]
			etband$Veg$AllMigration$BSE$level3 <- gap.rat$LEVEL3[match(etband$Veg$AllMigration$BSE$ratValues, gap.rat$Value, nomatch=0)]
			etband$Veg$AllMigration$BSE$level3_density <- bseFreq[match(etband$Veg$AllMigration$BSE$ratValues, bseFreq[, 1], nomatch=0), 2] / raster::ncell(temp_band4$bse)
			#majority TF
			tfFreq <- raster::freq(temp_band4$tf)
			etband$Veg$AllMigration$TF$ratValues <- tfFreq[tfFreq[, 1] %in% type_ids(type_veg2(ecotoner_settings)), 1]
			etband$Veg$AllMigration$TF$level3 <- gap.rat$LEVEL3[match(etband$Veg$AllMigration$TF$ratValues, gap.rat$Value, nomatch=0)]
			etband$Veg$AllMigration$TF$level3_density <- tfFreq[match(etband$Veg$AllMigration$TF$ratValues, tfFreq[, 1], nomatch=0), 2] / raster::ncell(temp_band4$tf)
			#rm(temp_stband, temp_etband, temp_band4, temp.patches1, temp.patches2, temp.human)

			#Quality: vegetation patch sizes: Gastner et al. 2009: we hypothesize that heterogeneity in the environmental background would not influence the model’s predictions unless the sizes of vegetation patches are smaller than the correlation length of the heterogeneity			
			#Sizes of vegetation patches
			etband$Veg$AllMigration$Veg1$MedianPatch_radius_m <- sqrt(median((temp <- raster::freq(etband$Veg$AllMigration$Veg1$patches8))[!is.na(temp[, 1]), 2]) * res_m(specs_grid(ecotoner_grids))^2 / pi)
			etband$Veg$AllMigration$Veg2$MedianPatch_radius_m <- sqrt(median((temp <- raster::freq(etband$Veg$AllMigration$Veg2$patches8))[!is.na(temp[, 1]), 2]) * res_m(specs_grid(ecotoner_grids))^2 / pi)
		} else { # else no transect successfully located
			status <- "notransect"
		}
	} else { # else: no candidates
		status <- "notransect"
	}

	# Return
	list(etband = etband, etable = etable, status = status)
}


#' Identify suitably migrated vegetation patches of the ecotone band transect, extract relevant information
#'
#' @export
identify_migration_patches <- function(i, b, ecotoner_settings, etband, etable, dir_fig, figBasename, verbose, do_figures, seed) {
	if(verbose) {
		cat("'ecotoner' migrating: tr = ", i, "; neigh = ", b, "; start at ", format(Sys.time(), format = ""), "\n", sep = "")
		idh <- 0 #counter for debug 'here' statements
	}
	
	#2e3ii. Veg$OnlyGoodMigration: Quality of x- vs y-axis migration routes: idea: identify patches that 'drain' out either at y=1 or at y=200; assumption: most likely migration route = flowpath
	#Use Rook's case re-calculated patches (patches4) and not rotates ones (patches8)
	flowdir <- raster::terrain(etband$Env$elev$grid, opt="flowdir")

	for (iveg in c("Veg1", "Veg2")) {
		if(verbose) cat("'ecotoner' migrating: tr = ", i, "; neigh = ", b, "; prog: ", idh <- idh + 1, "\n", sep = "")
		
		type_veg <- if (iveg == "Veg1") type_veg1(ecotoner_settings) else type_veg2(ecotoner_settings)
		migration <- calc.MigrationRoutes_EstimateFlowpaths(elev = etband$Env$elev$grid, flowdir = flowdir,
														patches = etband$Veg$AllMigration[[iveg]]$patches4,
														end_toLeft = end_to_left(type_veg), seed = seed)

		#Identify x vs y migration
		migration$direction <- calc.Identify_GoodvsBadMigration(patches4 = etband$Veg$AllMigration[[iveg]]$patches4,
												tally = migration$tally, paths = migration$paths,
												paths_success_TF = migration$paths_success_TF, end_toLeft = end_to_left(type_veg))

		#Get grid, patches4, patches8, and patchID_removed (patchID_removed == NULL if no patch removed)
		etband$Veg$OnlyGoodMigration[[iveg]] <- calc.RemoveBadMigration_fromVeg(patches4 = etband$Veg$AllMigration[[iveg]]$patches4,
																				patches8 = etband$Veg$AllMigration[[iveg]]$patches8,
																				patch4IDs_remove = migration$direction$patchIDs_GoodMigration)
		etband$Veg$OnlyGoodMigration[[iveg]]$Migration_Tally <- migration
		#rm(flowdir, migration)

		#Calculate missing variables for Veg$OnlyGoodMigration
		etband$Veg$OnlyGoodMigration[[iveg]]$density <- count_y_at_each_x(etband$Veg$OnlyGoodMigration[[iveg]]$grid)
	}

	etband$Veg$OnlyGoodMigration$diffs_amongMigTypes_TF <- !is.null(etband$Veg$OnlyGoodMigration$Veg1$patchID_removed) |
																	!is.null(etband$Veg$OnlyGoodMigration$Veg2$patchID_removed)

	#Store migration band transect data for analysis
	etable <- get.BandTransect_TableDataMigration(etable=etable, b=b, data=etband$Veg$OnlyGoodMigration)
	
	if(verbose) cat("'ecotoner' migrating: tr = ", i, "; neigh = ", b, "; prog: ", idh <- idh + 1, "\n", sep = "")

	if(do_figures) map_flowpaths(filename=file.path(dir_fig, paste0(figBasename, "FlowpathsApproximatingMigration.pdf")),
							elev=etband$Env$elev$grid, 
							bse=etband$Veg$AllMigration$Veg1$grid, bse_remain=etband$Veg$OnlyGoodMigration$Veg1$grid, paths_bse=etband$Veg$OnlyGoodMigration$Veg1$Migration_Tally$paths,
							tf=etband$Veg$AllMigration$Veg2$grid, tf_remain=etband$Veg$OnlyGoodMigration$Veg2$grid, paths_tf=etband$Veg$OnlyGoodMigration$Veg2$Migration_Tally$paths)
	
	list(etband = etband, etable = etable)
}



# Workhorse function
#' @export
detect_ecotone_transects_from_searchpoint <- function(i, initpoints, ecotoner_settings, ecotoner_grids, seed_streams = NULL, do_interim = TRUE, verbose = TRUE, do_figures = TRUE) {
	t1t <- Sys.time()
	if (verbose) {
		idh <- 0 #counter for debug statements
		cat("'ecotoner' detecting: tr = ", i, "; start at ", format(t1t, format = ""), "\n", sep = "")
	}
	
	# Clean up old temporary raster files; assuming that no call to this function took much longer than previous calls
	raster::removeTmpFiles(h = ceiling(get_max_timing(file_timing_locate(esets))))

	iflag <- flag_itransect(ecotoner_settings, i)

	# Result container
	etransect <- list(etbands = vector(mode = "list", length = neighbors_N(ecotoner_settings)), # list for 'etband' of each neighborhood
						status = rep("searching", neighbors_N(ecotoner_settings)), # "located", "notransect", "searching", "error"
						seeds = vector(mode = "list", length = neighbors_N(ecotoner_settings)),
						etable = list()) #container combining output table results from each neighborhood


	# Check if transect successfully located, none found, or error occurred, then do not continue with locating a transect from search point
	if (file.exists(fname_etlocated(ecotoner_settings, iflag))) {
		flag_status <- "located"
		load(fname_etlocated(ecotoner_settings, iflag)) # i, b, etransect
	
	} else if (file.exists(fname_etnone(ecotoner_settings, iflag))) {
		flag_status <- "notransect"
	
	} else if (file.exists(fname_etfailed(ecotoner_settings, iflag))) {
		flag_status <- "error"
	
	} else {
		flag_status <- "searching"
		
		ipoint <- initpoints[i, ]
		# Directory for temporary output
		dir_fig <- file.path(dir_out_fig(ecotoner_settings), iflag)
		dir_create(dir_fig)

		if(verbose) cat("'ecotoner' detecting: tr = ", i, "; prog: ", idh <- idh + 1, "\n", sep = "")
		
		#Crop rasters around ipoint to define linear transect: account for largest neighborhood
		temp <- crop_to_neighborhood(pt_start = ipoint,
									neighbor_n = max(neighborhoods(ecotoner_settings)),
									grid_gradient = grid_env(ecotoner_grids),
									grid_veg = grid_veg(ecotoner_grids),
									grid_veg_abut12 = grid_abut(ecotoner_grids),
									grid_aspect_mean = grid_aspect_mean(ecotoner_grids),
									grid_aspect_sd = grid_aspect_sd(ecotoner_grids),
									transect_type = transect_type(ecotoner_settings))
		
		elevCropped <- temp$grid_elev_cropped
		
		if (is.null(elevCropped)) {
			etransect[["status"]] <- rep("notransect", neighbors_N(ecotoner_settings))
			flag_status <- "notransect"
			b <- 0

		} else {
			gapCropped <- temp$grid_veg_cropped
			bseABUTtfCropped <- temp$grid_veg_abut12_cropped
			asp201MeanCropped <- temp$grid_aspect_mean_cropped
			asp201SDCropped <- temp$grid_aspect_sd_cropped
	
			for (b in seq_len(neighbors_N(ecotoner_settings))) { #if no transect can be established then proceed to next initpoint
				t1b <- Sys.time()
				do_timing <- FALSE
				
				flag_bfig <- flag_basename(ecotoner_settings, iflag, b)
				etransect[["status"]][b] <- "searching"
			
				# Check whether partial (debug) output was already generated, if so read and jump to continue
				tempData0 <- file.path(dir_fig, paste0("ecotoner_tmp0_candidates_tr", i, "_neigh", b, ".rds"))
				tempData1 <- file.path(dir_fig, paste0("ecotoner_tmp1_establish_tr", i, "_neigh", b, ".RData"))
				tempData2 <- file.path(dir_fig, paste0("ecotoner_tmp2_identify_tr", i, "_neigh", b, ".RData"))

				if (file.exists(tempData2)) {
					 # output from call to 'establish_ecotone_transect' available from previous run
					 # output from call to 'identify_migration_patches' available from previous run
					if(verbose) cat("'ecotoner' detecting: tr = ", i, "; prog: ", idh, ": 'do.tempData2' available and loaded\n", sep = "")
					
					load(file = tempData2)
					etransect[["etbands"]][[b]] <- etband
					etransect[["etable"]] <- etable
				} else {
					if (file.exists(tempData1)) {
						 # output from call to 'establish_ecotone_transect' available from previous run
						if(verbose) cat("'ecotoner' detecting: tr = ", i, "; prog: ", idh, ": 'do.tempData1' available and loaded\n", sep = "")
					
						load(file = tempData1)
						etransect[["etbands"]][[b]] <- etband
						etransect[["etable"]] <- etable
					} else {
						do_timing <- !file.exists(tempData0)

						# Set random seed and generator
						etransect[["seeds"]][[b]] <- if (is.null(seed_streams)) NULL else if (inherits(seed_streams, "list")) seed_streams[[(i - 1) * neighbors_N(ecotoner_settings) + b]] else NA
						set_RNG_stream(etransect[["seeds"]][[b]])
					
						temp <- try(establish_ecotone_transect(i, b, 
													etband = etransect[["etbands"]][[b]],
													etable = etransect[["etable"]],
													ipoint, ecotoner_settings, ecotoner_grids,
													elevCropped, gapCropped, bseABUTtfCropped, asp201MeanCropped, asp201SDCropped,
													tempData = if (do_interim) tempData0 else NULL,
													dir_fig, figBasename = flag_bfig,
													verbose, do_figures,
													seed = NA), silent = TRUE)
					
						if (!inherits(temp, "try-error")) {
							etransect[["status"]][b] <- temp[["status"]]
						
							if (temp[["status"]] == "searching") {
								etransect[["etbands"]][[b]] <- temp[["etband"]]
								etransect[["etable"]] <- temp[["etable"]]
							}
							rm(temp)
						} else {
							etransect[["status"]][b] <- "error"
							warning("'detect_ecotone_transects_from_searchpoint': ", temp, immediate. = TRUE)
						}
					
						if (do_interim) {
							#saveAll <- ls(all=TRUE, name=environment())
							#if(!identical(environment(), sys.frame())) saveAll <- c(saveAll, ls(all=TRUE, name=sys.frame()))
							etband <- etransect[["etbands"]][[b]]
							etable <- etransect[["etable"]]
							save(etband, etable, file = tempData1)
						}
					} # end do.tempData1
				
					if (etransect[["status"]][b] == "searching") {
						# Re-set random seed and generator (in case code reentry is here after interruption and tempData1 object already stored to disk)
						set_RNG_stream(etransect[["seeds"]][[b]])
					
						temp <- try(identify_migration_patches(i, b, ecotoner_settings,
																etband = etransect[["etbands"]][[b]],
																etable = etransect[["etable"]],
																dir_fig, flag_bfig, verbose, do_figures,
																seed = NA), silent = TRUE)
				
						if (!inherits(temp, "try-error")) {
							etransect[["etbands"]][[b]] <- temp[["etband"]]
							etransect[["etable"]] <- temp[["etable"]]
						} else {
							etransect[["status"]][b] <- "error"
							warning("'detect_ecotone_transects_from_searchpoint': ", temp, immediate. = TRUE)
						}

						if (do_interim) {
							etband <- etransect[["etbands"]][[b]]
							etable <- etransect[["etable"]]
							save(etband, etable, file = tempData2)
						}
					}
				} # end do.tempData2

				#----3. Loop through both versions of vegetation for applying the methods : Veg$AllMigration and Veg$OnlyGoodMigration
				if (etransect[["status"]][b] == "searching") {
					for (migtype in get("migtypes", envir = etr_vars)) {
						copy_FromMig1_TF <- if (migtype == "AllMigration") FALSE else !etransect$etbands[[b]]$Veg[["OnlyGoodMigration"]]$diffs_amongMigTypes_TF
						if (verbose) cat("'ecotoner' detecting: tr = ", i, "; neigh = ", b, ": prog: ", idh, "; mig-type: ", migtype, "\n", sep = "")
					
						if (!copy_FromMig1_TF) {
							#Extract abutting cells between Veg1 and Veg2
							etransect$etbands[[b]]$Veg[[migtype]]$Veg1ABUTVeg2$grid <- calc_abutting(Veg1=etransect$etbands[[b]]$Veg[[migtype]]$Veg1$grid, Veg2=etransect$etbands[[b]]$Veg[[migtype]]$Veg2$grid)
							etransect$etbands[[b]]$Veg[[migtype]]$Veg1ABUTVeg2$density <- count_y_at_each_x(etransect$etbands[[b]]$Veg[[migtype]]$Veg1ABUTVeg2$grid)

							#Store vegetation band transect data for analysis
							etransect$etable <- get.BandTransect_TableDataVeg(etable=etransect$etable, b=b, data=etransect$etbands[[b]]$Veg[[migtype]], flag_migtype = migtype)
						} else {
							etransect$etable <- get.BandTransect_TableDataVeg(etable=etransect$etable, b=b, data=etransect$etbands[[b]]$Veg[["AllMigration"]], flag_migtype = migtype)
						}
					}#end loop through migtypes

					#Temporarily save data to disk file
					temp <- write_ecotoner_row(data_row = etransect$etable[b, ],
												filename = file_etsummary_temp(ecotoner_settings),
												tag_fun = 'detect_ecotone_transects_from_searchpoint',
												tag_id = paste0(i, " 'etransect$etable[", b, ", ]'"))
					if (identical(temp, "error")) etransect[["status"]][b] <- "error"
					
					if(verbose) cat("'ecotoner' detecting: tr = ", i, "; neigh = ", b, ": prog: ", idh <- idh + 1, "\n", sep = "")
				}
				
				# If code advanced until here, then transect located and characterized successfully
				if (etransect[["status"]][b] == "searching") etransect[["status"]][b] <- "located"
				
				# Save data to RData file
				if (etransect[["status"]][b] == "error") {
					break
				} else if (b < neighbors_N(ecotoner_settings)) {
					save(i, b, etransect, file = fname_etsearching(ecotoner_settings, iflag))
				}

				# Timing information: add only if everything was freshly calculated
				if (do_timing) {
					t2b <- Sys.time()
					# Timing for search point and this neighborhood
					add_new_timing(i, b, time_h = difftime(t2b, t1b, units = "hours"), filename = file_timing_locate(esets))
					# Timing for search point and all neighborhoods
					add_new_timing(i, -1, time_h = difftime(t2b, t1t, units = "hours"), filename = file_timing_locate(esets))
				}
			} # end for-loop through neighborhoods


			# Save data to RData file and remove no longer used files
			remove_fsearching <- FALSE
			
			if (all(etransect[["status"]] == "located")) {
				save(i, b, etransect, file = fname_etlocated(ecotoner_settings, iflag))
				remove_fsearching <- TRUE				
				flag_status <- "located"
				unlink(x = c(file.path(dir_fig, "ecotoner_tmp0_candidates_tr*.rds"),
							 file.path(dir_fig, "ecotoner_tmp1_establish_tr*.RData"),
							 file.path(dir_fig, "ecotoner_tmp2_identify_tr*.RData")))
			
			} else {
				unlink(dir_fig, recursive = TRUE)
				
				if (any(etransect[["status"]] == "error")) {
					save(i, b, etransect, file = fname_etfailed(ecotoner_settings, iflag))
					remove_fsearching <- TRUE
					flag_status <- "error"
				} else if (any(etransect[["status"]] == "notransect")) {
					save(i, b, etransect, file = fname_etnone(ecotoner_settings, iflag))
					remove_fsearching <- TRUE
					flag_status <- "notransect"
				}
			}
			
			if (remove_fsearching){
				unlink(fname_etsearching(ecotoner_settings, iflag)) #Not deleting a non-existent file is not a failure
			}
		}
	}
	
	if (verbose) cat("'ecotoner' detecting: tr = ", i, "; completed with status '", flag_status, "' at ", format(t2t <- Sys.time(), format = ""), " after ", round(difftime(t2t, t1t, units = "hours"), 2), " hours\n", sep = "")

	etransect$etable #if unsuccessful then NULL else data.frame
}	
