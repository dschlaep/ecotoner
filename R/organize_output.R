#------------------------------------------------------------#
#Functions writing data in temporary table

get.ElevationTransect_tempTableData <- function(temp_Table, data) {
	temp_Table$Transect_Azimuth_deg <- if (requireNamespace("geosphere", quietly = TRUE)) {
											geosphere::bearing(data$endPoints_WGS84[1,], data$endPoints_WGS84[2,])
										} else {
											warning("Package 'geosphere' not installed: the azimuth of the transect will not be calculated by the function 'get.ElevationTransect_tempTableData'")
											NA
										}
	temp_Table$Transect_DirectionAEA_deg <- 180 / pi * atan2((temp <- sp::coordinates(data$endPoints))[1,2]-temp[2,2], temp[1,1]-temp[2,1])
	temp_Table$TransectElevation_Length_m <- 1000 * sp::spDists(data$endPoints_WGS84[1,], data$endPoints_WGS84[2,], longlat = TRUE)


	temp_Table[c("TransectElevation_StartPoint_X",
				"TransectElevation_StartPoint_Y",
				"TransectElevation_StartPoint_Elev_m")] <- get_xyz(data$endPoints[1, ], "elev")
	temp_Table[c("TransectElevation_EndPoint_X",
				"TransectElevation_EndPoint_Y",
				"TransectElevation_EndPoint_Elev_m")] <- get_xyz(data$endPoints[2, ], "elev")

	temp_Table[c("TransectElevation_StartPoint_WGS84_Long",
				"TransectElevation_StartPoint_WGS84_Lat",
				"TransectElevation_StartPoint_WGS84_Elev_m")] <- get_xyz(data$endPoints_WGS84[1, ], "elev")
	temp_Table[c("TransectElevation_EndPoint_WGS84_Long",
				"TransectElevation_EndPoint_WGS84_Lat",
				"TransectElevation_EndPoint_WGS84_Elev_m")] <- get_xyz(data$endPoints_WGS84[2, ], "elev")
	
	return(temp_Table)
}

tabulate_merge_into_etable <- function(etable, index, data) {
	cname <- if (is.vector(data)) names(data) else colnames(data)
	cname_exist <- match(cname, colnames(etable), nomatch = 0)
	
	icol <- cname_exist > 0
	if (any(icol)) etable[index, cname_exist[icol]] <- data[icol]

	icol <- cname_exist == 0
	if (any(icol)) {
		res <- as.data.frame(matrix(NA, nrow = max(index, nrow(etable)), ncol = sum(icol), dimnames = list(NULL, cname[icol])))
		res[index, ] <- data[icol]
		etable <- cbind(etable, res)
	}
	
	etable
}


get.CopyTempTo_TableData <- function(temp_Table, etable, b) {
	cname <- colnames(temp_Table)
	cname_exist <- match(cname, colnames(etable), nomatch = 0)
	if (any(icol <- (cname_exist > 0))) {
		etable[b, cname_exist[icol]] <- unlist(temp_Table[, icol])
	}
	if (any(icol <- (cname_exist == 0))) {
		res <- matrix(NA, nrow = max(1, nrow(etable)), ncol = sum(icol), dimnames = list(NULL, cname[icol]))
		res[b, ] <- unlist(temp_Table[, icol])
		etable <- if (length(etable) > 0) data.frame(cbind(etable, res)) else data.frame(res)
	}
	
	etable
}		

put.ListData1Level_TableData <- function(etable, b, data, colFlag, listname, flag_migtype){
	cname <- paste(flag_migtype, colFlag, listname, names(data), sep="_")
	cname_exist <- match(cname, colnames(etable), nomatch=0)
	if(any(icol <- (cname_exist > 0))){
		etable[b, cname_exist[icol]] <- unlist(data[icol])
	}
	if(any(icol <- (cname_exist == 0))){
		res <- matrix(NA, nrow=max(1, nrow(etable)), ncol=sum(icol), dimnames=list(NULL, cname[icol]))
		res[b, ] <- unlist(data[icol])
		etable <- cbind(etable, res)
	}
	return(etable)
}



get.LinearTransect_TableData <- function(etable, b, data) {
	etable$Transect_StartInElevationTransect_m[b] <- data$StartInElevationTransect_m

	etable$TransectLinear_Length_m[b] <- sp::spDists(data$endPoints_WGS84[1,], data$endPoints_WGS84[2,], longlat=TRUE) * 1000

	ids <- c("TransectLinear_XXXPoint_WGS84_Long", "TransectLinear_XXXPoint_WGS84_Lat", "TransectLinear_XXXPoint_WGS84_Elev_m")
	ids_temp <- gsub("XXX", "Start", ids)
	etable[, ids_temp] <- NA
	etable[b, ids_temp] <- get_xyz(data$endPoints_WGS84[1, ], "elev")
	ids_temp <- gsub("XXX", "End", ids)
	etable[, ids_temp] <- NA
	etable[b, ids_temp] <- get_xyz(data$endPoints_WGS84[2, ], "elev")

	etable$TransectLinear_ElevMax_m[b] <- max(data$pts$elev)
	etable$TransectLinear_ElevMin_m[b] <- min(data$pts$elev)
	etable$TransectLinear_ElevDiff_m[b] <- abs(etable$TransectLinear_EndPoint_WGS84_Elev_m[b] - etable$TransectLinear_StartPoint_WGS84_Elev_m[b])
	etable$TransectLinear_ElevDiffSum_m[b] <- sum(abs(data$pts$elev[-1] - data$pts$elev[-data$pts_N]))
	etable$TransectLinear_MeanSlope_mpm[b] <- etable$TransectLinear_ElevDiff_m[b] / max(data$dist_m)
	etable$TransectLinear_MeanSlope_deg[b] <- atan(etable$TransectLinear_MeanSlope_mpm[b]) * 180/pi
	etable$TransectLinear_SlopeRoughness_mpm[b] <- (etable$TransectLinear_ElevDiffSum_m[b] - etable$TransectLinear_ElevDiff_m[b]) / max(data$dist_m)

	etable
}

get.BandTransect_TableDataEnv <- function(etable, b, data){
	ncells <- raster::ncell(data$elev$grid)
	res_m <- raster::xres(data$elev$grid)

	etable$TransectBand_Length_m[b] <- ncol(data$elev$grid) * res_m
	etable$TransectBand_Width_m[b] <- nrow(data$elev$grid) * res_m
	etable$TransectBand_Area_m2[b] <- ncells * res_m^2

	etable$TransectBand_Quality_HumanFootprint_Area_Fraction[b] <- raster::cellStats(data$human$grid, 'sum') / ncells
	etable$TransectBand_Quality_HumanFootprint_Mode_XLocation_m[b] <- ifelse(any(data$human$density > 0), which.max(data$human$density) * res_m, 0)
	etable$TransectBand_Quality_HumanFootprint_LargestClump_XExtent_Fraction[b] <- data$human$LargestClump_XExtent / etable$TransectBand_Length_m[b]
	etable$TransectBand_Quality_HumanFootprint_LargestClump_YExtent_Fraction[b] <- data$human$LargestClump_YExtent / etable$TransectBand_Width_m[b]
	
	etable$TransectBand_Quality_Aspect_Cellwise_RMSE_rad[b] <- data$aspect$Aspect_Cellwise_RMSE
	etable$TransectBand_Quality_Aspect_Cellwise_95PercQuantileAbsError_rad[b] <- data$aspect$Aspect_Cellwise_95PercQuantileAbsError
	etable$TransectBand_Quality_Aspect_YMeans_RMSE_rad[b] <- data$aspect$Aspect_YMeans_RMSE
	etable$TransectBand_Quality_Aspect_YMeans_95PercQuantileAbsError_rad[b] <- data$aspect$Aspect_YMeans_95PercQuantileAbsError
	etable$TransectBand_Quality_Aspect_Patch_LargestDeviatingPatch_Fraction[b] <- data$aspect$Aspect_Patch_FractionLargestDeviatingPatch
	etable$TransectBand_Quality_Aspect_Patch_AllDeviatingPatches_Fraction[b] <- data$aspect$Aspect_Patch_FractionAllDeviatingPatches
	etable$TransectBand_Quality_Aspect_Overall_MeanError_rad[b] <- data$aspect$Aspect_Overall_ME
	
	etable$TransectBand_Quality_Elevation_Cellwise_lmResiduals_RMSE_m[b] <- data$elev$Elevation_Cellwise_lmResiduals_RMSE
	etable$TransectBand_Quality_Elevation_Cellwise_lmResiduals_95PercQuantileAbsError_m[b] <- data$elev$Elevation_Cellwise_lmResiduals_95PercQuantileAbsError
	etable$TransectBand_Quality_Elevation_Cellwise_lmResiduals_MaxAbsError_m[b] <- data$elev$Elevation_Cellwise_lmResiduals_MaxAbsError
	etable$TransectBand_Quality_Elevation_Cellwise_YMeans_RMSE_m[b] <- data$elev$Elevation_Cellwise_YMeans_RMSE
	etable$TransectBand_Quality_Elevation_Cellwise_YMeans_95PercQuantileAbsError_m[b] <- data$elev$Elevation_Cellwise_YMeans_95PercQuantileAbsError
	etable$TransectBand_Quality_Elevation_Cellwise_YMeans_MaxAbsError_m[b] <- data$elev$Elevation_Cellwise_YMeans_MaxAbsError
	etable$TransectBand_Quality_Elevation_Cellwise_YMeans_MoransI[b] <- data$elev$Elevation_Cellwise_YMeans_MoransI
	etable$TransectBand_Quality_Elevation_Cellwise_YMeans_MoransI_pvalue[b] <- data$elev$Elevation_Cellwise_YMeans_MoransIp
	etable$TransectBand_Quality_Elevation_RatioYrangeToExpectXrange_RMSE[b] <- data$elev$Elevation_RatioYrangeToExpectXrange_RMSE
	etable$TransectBand_Quality_Elevation_RatioYrangeToExpectXrange_95PercQuantileRatioAbsError[b] <- data$elev$Elevation_RatioYrangeToExpectXrange_95PercQuantileRatioAbsError
	etable$TransectBand_Quality_Elevation_RatioYrangeToExpectXrange_MaxAbsError[b] <- data$elev$Elevation_RatioYrangeToExpectXrange_MaxAbsError
	etable$TransectBand_Quality_Elevation_YMeansResiduals_VariogramRange_m[b] <- data$elev$Elevation_YMeansResiduals_VariogramRange

	etable$TransectBand_meanElevMax_m[b] <- max(data$elev$YMeans_ForEachX)
	etable$TransectBand_meanElevMin_m[b] <- min(data$elev$YMeans_ForEachX)
	etable$TransectBand_meanElevDiff_m[b] <- abs(data$elev$YMeans_ForEachX[length(data$elev$YMeans_ForEachX)] - data$elev$YMeans_ForEachX[1])
	etable$TransectBand_meanElevDiffSum_m[b] <- sum(abs(data$elev$YMeans_ForEachX[-1] - data$elev$YMeans_ForEachX[-length(data$elev$YMeans_ForEachX)]))
	etable$TransectBand_MeanSlope_mpm[b] <- etable$TransectBand_meanElevDiff_m[b] / etable$TransectBand_Length_m[b]
	etable$TransectBand_MeanSlope_deg[b] <- atan(etable$TransectBand_MeanSlope_mpm[b]) * 180/pi
	etable$TransectBand_SlopeRoughness_mpm[b] <- (etable$TransectBand_meanElevDiffSum_m[b] - etable$TransectBand_meanElevDiff_m[b]) / etable$TransectBand_Length_m[b]

	return(etable)
}

get.BandTransect_TableDataMigration <- function(etable, b, data){
	etable$TransectBand_Quality_Diffs_amongMigrationTypes_TF[b] <- data$diffs_amongMigTypes_TF

	etable$TransectBand_Quality_Veg1_Migration_CellsTested_Number[b] <- sum(data$Veg1$Migration_Tally$tally[, "TestedFlowpaths_N"])

	tally <- data$Veg1$Migration_Tally$tally
	tallyX <- tally[tally[, "PatchID"] %in% data$Veg1$Migration_Tally$direction$patchIDs_BadMigration, , drop=FALSE]
	tallyY <- tally[tally[, "PatchID"] %in% data$Veg1$Migration_Tally$direction$patchIDs_GoodMigration, , drop=FALSE]
	etable$TransectBand_Quality_Veg1_PatchesGoodMigration_Number[b] <- nrow(tallyX)
	etable$TransectBand_Quality_Veg1_PatchesGoodMigration_NumberFraction[b] <- ifelse(nrow(tally) > 0, nrow(tallyX) / nrow(tally), NA)
	etable$TransectBand_Quality_Veg1_PatchesGoodMigration_AreaFraction[b] <- ifelse(nrow(tally) > 0, sum(tallyX[, "PatchCells_N"]) / sum(tally[, "PatchCells_N"]), NA)
	etable$TransectBand_Quality_Veg1_PatchesBadMigration_Number[b] <- nrow(tallyY)
	etable$TransectBand_Quality_Veg1_PatchesBadMigration_NumberFraction[b] <- ifelse(nrow(tally) > 0, nrow(tallyY) / nrow(tally), NA)
	etable$TransectBand_Quality_Veg1_PatchesBadMigration_AreaFraction[b] <- ifelse(nrow(tally) > 0, sum(tallyY[, "PatchCells_N"]) / sum(tally[, "PatchCells_N"]), NA)

	tally <- data$Veg2$Migration_Tally$tally
	tallyX <- tally[tally[, "PatchID"] %in% data$Veg2$Migration_Tally$direction$patchIDs_BadMigration, , drop=FALSE]
	tallyY <- tally[tally[, "PatchID"] %in% data$Veg2$Migration_Tally$direction$patchIDs_GoodMigration, , drop=FALSE]
	etable$TransectBand_Quality_Veg2_PatchesGoodMigration_Number[b] <- nrow(tallyX)
	etable$TransectBand_Quality_Veg2_PatchesGoodMigration_NumberFraction[b] <- ifelse(nrow(tally) > 0, nrow(tallyX) / nrow(tally), NA)
	etable$TransectBand_Quality_Veg2_PatchesGoodMigration_AreaFraction[b] <- ifelse(nrow(tally) > 0, sum(tallyX[, "PatchCells_N"]) / sum(tally[, "PatchCells_N"]), NA)
	etable$TransectBand_Quality_Veg2_PatchesBadMigration_Number[b] <- nrow(tallyY)
	etable$TransectBand_Quality_Veg2_PatchesBadMigration_NumberFraction[b] <- ifelse(nrow(tally) > 0, nrow(tallyY) / nrow(tally), NA)
	etable$TransectBand_Quality_Veg2_PatchesBadMigration_AreaFraction[b] <- ifelse(nrow(tally) > 0, sum(tallyY[, "PatchCells_N"]) / sum(tally[, "PatchCells_N"]), NA)

	return(etable)
}

get.BandTransect_TableDataVeg <- function(etable, b, data, flag_migtype){
	ncells <- raster::ncell(data$Veg1$grid)

	put.VegData <- function(etable, b, data, flag_migtype, vegno){
		colnamesAdd <- paste0("TransectBand_", flag_migtype, "_Veg", vegno, "_", "rDensity")
		if(flag_migtype==get("migtypes", envir = etr_vars)[1]){
			colnamesAdd <- c(colnamesAdd, paste0("TransectBand_", flag_migtype, "_Veg", vegno, "_",
							c("Quality_MedianPatch_radius_m", "MajorityType")))
		}
		res <- data.frame(matrix(NA, nrow=max(1, nrow(etable)), ncol=length(colnamesAdd), dimnames=list(NULL, colnamesAdd)))

		res[b, 1] <- raster::cellStats(data$grid, 'sum', na.rm = TRUE) / ncells
		if(flag_migtype==get("migtypes", envir = etr_vars)[1]){
			res[b, 2] <- data$MedianPatch_radius_m
			res[b, 3] <- data$type_majority
		}
		
		etable <- cbind(etable, res)
		return(etable)
	}
	etable <- put.VegData(etable=etable, b=b, data=data$Veg1, flag_migtype=flag_migtype, vegno=1)
	etable <- put.VegData(etable=etable, b=b, data=data$Veg2, flag_migtype=flag_migtype, vegno=2)

	put.JointData <- function(etable, b, data, flag_migtype){
		colnamesAdd <- paste0("TransectBand_", flag_migtype, "_Veg1and2_",
							c("Abutting_rDensity", "rDensity"))
		if(flag_migtype==get("migtypes", envir = etr_vars)[1]){
			colnamesAdd <- c(colnamesAdd, paste0("TransectBand_", flag_migtype, "_",
							c("BSE_rDensity", "TF_rDensity", "BSEandTF_rDensity", "Veg1_rDensityToBSE", "Veg2_rDensityToTF")))
		}
		res <- matrix(NA, nrow=max(1, nrow(etable)), ncol=length(colnamesAdd), dimnames=list(NULL, colnamesAdd))

		res[b, 1] <- raster::cellStats(data$Veg1ABUTVeg2$grid, 'sum', na.rm = TRUE) / ncells
		Veg1_rDens <- etable[b, colnames(etable) == paste0("TransectBand_", flag_migtype, "_Veg1_rDensity")]
		Veg2_rDens <- etable[b, colnames(etable) == paste0("TransectBand_", flag_migtype, "_Veg2_rDensity")]
		res[b, 2] <- Veg1_rDens + Veg2_rDens
		if(flag_migtype==get("migtypes", envir = etr_vars)[1]){
			res[b, 3] <- sum(data$BSE$level3_density)
			res[b, 4] <- sum(data$TF$level3_density)
			res[b, 5] <- res[b, 3] + res[b, 4]
			res[b, 6] <- Veg1_rDens / res[b, 3]
			res[b, 7] <- Veg2_rDens / res[b, 4]
		}
		
		etable <- cbind(etable, res)
		return(etable)
	}
	etable <- put.JointData(etable=etable, b=b, data=data, flag_migtype=flag_migtype)

	return(etable)
}


tabulate_InterZoneLocation <- function(etable, b, data, colFlag, flag_migtype){
	iz_N <- length(data)
	listname <- names(data)
	for (iz in seq_len(iz_N)) {
		etable <- put.ListData1Level_TableData(etable=etable, b=b, data=data[[iz]], colFlag=colFlag, listname=listname[iz], flag_migtype=flag_migtype)
	}
	
	return(etable)
}

get.EcolBoundary_ListData2Levels_TableData <- function(etable, b, data, colFlag, flag_migtype){
	iz_N <- length(data)
	listname <- names(data)
	for (iz in seq_len(iz_N)) {
		etable <- put.ListData1Level_TableData(etable=etable, b=b, data=data[[iz]]$Veg1, colFlag=paste0(colFlag, "_Veg1"), listname=listname[iz], flag_migtype=flag_migtype)
		etable <- put.ListData1Level_TableData(etable=etable, b=b, data=data[[iz]]$Veg2, colFlag=paste0(colFlag, "_Veg2"), listname=listname[iz], flag_migtype=flag_migtype)
	}
	
	return(etable)
}

