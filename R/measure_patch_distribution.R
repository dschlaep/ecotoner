calc.EcolBoundary_PatchSizeDistribution <- function(eBIZ_data, gridGlobalVegClumps, GlobalVegClumpFreq, dist_AlongTransect, clumpAreaClasses_m2){
	lres <- list(Below_Area_m2=NA, Within_Area_m2=NA, Above_Area_m2=NA,
				Below_Clumps_N=NA, Within_Clumps_N=NA, Above_Clumps_N=NA,
				Below_Clumps_Area_Median_m2=NA, Within_Clumps_Area_Median_m2=NA, Above_Clumps_Area_Median_m2=NA,
				Below_Clumps_Area_SmallClass_N=NA, Below_Clumps_Area_MediumClass_N=NA, Below_Clumps_Area_LargeClass_N=NA, 
				Within_Clumps_Area_SmallClass_N=NA, Within_Clumps_Area_MediumClass_N=NA, Within_Clumps_Area_LargeClass_N=NA, 
				Above_Clumps_Area_SmallClass_N=NA, Above_Clumps_Area_MediumClass_N=NA, Above_Clumps_Area_LargeClass_N=NA)

	if(!is.na(eBIZ_data$start_m) & !is.na(eBIZ_data$end_m)){
		belowID <- 1
		startID <- findInterval(eBIZ_data$start_m, dist_AlongTransect, all.inside=TRUE)
		endID <- findInterval(eBIZ_data$end_m, dist_AlongTransect, all.inside=TRUE)
		aboveID <- length(dist_AlongTransect)
		belowIDs <- belowID:(startID-1)
		withinIDs <- startID:endID
		aboveIDs <- (endID+1):aboveID
	
		grid_zones <- raster(gridGlobalVegClumps)
		grid_zones[, belowIDs] <- 1
		grid_zones[, withinIDs] <- 2
		grid_zones[, aboveIDs] <- 3
		res_m <- raster::xres(gridGlobalVegClumps)
		
		temp.clumpsN <- zonal(x=gridGlobalVegClumps, z=grid_zones, fun=function(x, ...) length(unique(na.omit(x))))
		dimnames(temp.clumpsN) <- NULL
		temp.clumpSizes <- zonal(x=gridGlobalVegClumps, z=grid_zones, fun=function(x, ...) list(GlobalVegClumpFreq$count[GlobalVegClumpFreq$value %in% x] * res_m^2))
		temp.clumpSizesMedian <- sapply(temp.clumpSizes[, 2], FUN=function(x) median(unlist(x)))
		temp.clumpSizeClassesN <- sapply(temp.clumpSizes[, 2], FUN=function(x) c(sum((temp <- unlist(x)) <= clumpAreaClasses_m2[1]), sum(temp > clumpAreaClasses_m2[1] & temp <= clumpAreaClasses_m2[2]), sum(temp > clumpAreaClasses_m2[2])))
		if(nrow(temp.clumpsN) < 3){ #transfer values (guaranteeing good structure of clumpsN, clumpSizesMedian, and clumpSizeClassesN, in case there are no clumps in a zone)
			clumpsN <- matrix(0, nrow=3, ncol=2)
			clumpsN[, 1] <- 1:3 # ids of zones 1-3
			clumpSizesMedian <- rep(0, times=3)
			clumpSizeClassesN <- matrix(0, nrow=3, ncol=3)
			for(iz in 1:3) if(any(temp <- temp.clumpsN[, 1] == iz)){
				clumpsN[temp.clumpsN[temp, 1], 2] <- temp.clumpsN[temp, 2]
				clumpSizesMedian[temp.clumpsN[temp, 1]] <- temp.clumpSizesMedian[temp]
				clumpSizeClassesN[, temp.clumpsN[temp, 1]] <- temp.clumpSizeClassesN[, temp]
			}
		} else {
			clumpsN <- temp.clumpsN
			clumpSizesMedian <- temp.clumpSizesMedian
			clumpSizeClassesN <- temp.clumpSizeClassesN
		}

		lres <- list(Below_Area_m2=length(belowIDs) * (temp <- nrow(gridGlobalVegClumps) * res_m^2), Within_Area_m2=length(withinIDs) * temp, Above_Area_m2=length(aboveIDs) * temp,
					Below_Clumps_N=clumpsN[1, 2], Within_Clumps_N=clumpsN[2, 2], Above_Clumps_N=clumpsN[3, 2],
					Below_Clumps_Area_Median_m2=clumpSizesMedian[1], Within_Clumps_Area_Median_m2=clumpSizesMedian[2], Above_Clumps_Area_Median_m2=clumpSizesMedian[3],
					Below_Clumps_Area_SmallClass_N=clumpSizeClassesN[1, 1], Below_Clumps_Area_MediumClass_N=clumpSizeClassesN[2, 1], Below_Clumps_Area_LargeClass_N=clumpSizeClassesN[3, 1], 
					Within_Clumps_Area_SmallClass_N=clumpSizeClassesN[1, 2], Within_Clumps_Area_MediumClass_N=clumpSizeClassesN[2, 2], Within_Clumps_Area_LargeClass_N=clumpSizeClassesN[3, 2], 
					Above_Clumps_Area_SmallClass_N=clumpSizeClassesN[1, 3], Above_Clumps_Area_MediumClass_N=clumpSizeClassesN[2, 3], Above_Clumps_Area_LargeClass_N=clumpSizeClassesN[3, 3])
	}
	
	return(lres)
}


#' @export
InterZonePatchDistr <- function(i, b, im, ecotoner_settings, etband, template_etobs, flag_bfig, copy_FromMig1_TF, do_figures, ...) {
	#---5. Patch size distribution of vegetation on the band transect below, within, and above interaction zone
	eB_PatchSizeDist[[im]] <- vector(mode="list", length=2+length(stepsHullEdge(ecotoner_settings)))
	names(eB_PatchSizeDist[[im]]) <- c("Danz2012", "Eppinga2013", paste0("Gastner2009Step", stepsHullEdge(ecotoner_settings)))

#	eB_methods <- list(...)
stop()
	if (!copy_FromMig1_TF) {
		eB_PatchSizeDist[[im]]$Danz2012$Veg1 <- calc.EcolBoundary_PatchSizeDistribution(eBIZ_data=eB_InterZone[[im]]$Danz2012, gridGlobalVegClumps=etband$Veg[[im]]$Veg1$patches8, GlobalVegClumpFreq=etband$Veg$AllMigration$Veg1$patches8_freq, dist_AlongTransect=etband$Env$DistAlongXaxis_m, clumpAreaClasses_m2 = clumpAreaClasses_m2(ecotoner_settings))
		eB_PatchSizeDist[[im]]$Danz2012$Veg2 <- calc.EcolBoundary_PatchSizeDistribution(eBIZ_data=eB_InterZone[[im]]$Danz2012, gridGlobalVegClumps=etband$Veg[[im]]$Veg2$patches8, GlobalVegClumpFreq=etband$Veg$AllMigration$Veg2$patches8_freq, dist_AlongTransect=etband$Env$DistAlongXaxis_m, clumpAreaClasses_m2 = clumpAreaClasses_m2(ecotoner_settings))
		eB_PatchSizeDist[[im]]$Eppinga2013$Veg1 <- calc.EcolBoundary_PatchSizeDistribution(eBIZ_data=eB_InterZone[[im]]$Eppinga2013, gridGlobalVegClumps=etband$Veg[[im]]$Veg1$patches8, GlobalVegClumpFreq=etband$Veg$AllMigration$Veg1$patches8_freq, dist_AlongTransect=etband$Env$DistAlongXaxis_m, clumpAreaClasses_m2 = clumpAreaClasses_m2(ecotoner_settings))
		eB_PatchSizeDist[[im]]$Eppinga2013$Veg2 <- calc.EcolBoundary_PatchSizeDistribution(eBIZ_data=eB_InterZone[[im]]$Eppinga2013, gridGlobalVegClumps=etband$Veg[[im]]$Veg2$patches8, GlobalVegClumpFreq=etband$Veg$AllMigration$Veg2$patches8_freq, dist_AlongTransect=etband$Env$DistAlongXaxis_m, clumpAreaClasses_m2 = clumpAreaClasses_m2(ecotoner_settings))
		for (i_step in seq_along(stepsHullEdge(ecotoner_settings))) {
			eB_PatchSizeDist[[im]][[2+i_step]]$Veg1 <- calc.EcolBoundary_PatchSizeDistribution(eBIZ_data=eB_InterZone[[im]][[2+i_step]], gridGlobalVegClumps=etband$Veg[[im]]$Veg1$patches8, GlobalVegClumpFreq=etband$Veg$AllMigration$Veg1$patches8_freq, dist_AlongTransect=etband$Env$DistAlongXaxis_m, clumpAreaClasses_m2 = clumpAreaClasses_m2(ecotoner_settings))
			eB_PatchSizeDist[[im]][[2+i_step]]$Veg2 <- calc.EcolBoundary_PatchSizeDistribution(eBIZ_data=eB_InterZone[[im]][[2+i_step]], gridGlobalVegClumps=etband$Veg[[im]]$Veg2$patches8, GlobalVegClumpFreq=etband$Veg$AllMigration$Veg2$patches8_freq, dist_AlongTransect=etband$Env$DistAlongXaxis_m, clumpAreaClasses_m2 = clumpAreaClasses_m2(ecotoner_settings))
		}
		
		etable <- get.EcolBoundary_ListData2Levels_TableData(etable=etable, b=b, data=eB_PatchSizeDist[[im]], colFlag="PatchSizeDist", flag_migtype = migtype)
	} else {
		etable <- get.EcolBoundary_ListData2Levels_TableData(etable=etable, b=b, data=eB_PatchSizeDist[[1]], colFlag="PatchSizeDist", flag_migtype = migtype)
	}
	
	template_etobs$gETmeas[[b]] <- eB_PatchSizeDist
	template_etobs$etable[b] <- etable
	
	template_etobs
}

