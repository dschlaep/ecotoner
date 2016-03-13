calc_InterZoneLocation <- function(center, start, end, dist_AlongTransect, elev_AlongTransect, veg1_AlongTransect, veg2_AlongTransect, abut_AlongTransect) {
	center <- c(center)
	start <- c(start)
	end <- c(end)
	
	lres <- list(center_m = NA, start_m = NA, end_m = NA, width_m = NA,
				elev_atCenter_m = NA, elev_atStart_m = NA, elev_atEnd_m = NA, elev_max_m = NA, elev_min_m = NA, elev_totalGain_m = NA,
				slope_mpm = NA, slope_deg = NA, slopeRoughness_mpm = NA,
				veg1_fraction = NA, veg2_fraction = NA, abut_fraction = NA,
				veg1_below_densMean = NA, veg1_below_densSD = NA, veg1_within_densMean = NA, veg1_within_densSD = NA, veg1_above_densMean = NA, veg1_above_densSD = NA,
				veg2_below_densMean = NA, veg2_below_densSD = NA, veg2_within_densMean = NA, veg2_within_densSD = NA, veg2_above_densMean = NA, veg2_above_densSD = NA,
				abut_below_densMean = NA, abut_below_densSD = NA, abut_within_densMean = NA, abut_within_densSD = NA, abut_above_densMean = NA, abut_above_densSD = NA)

	
	if (!is.na(start) && !is.na(end)) {
		width <- abs(end - start)
		
		startID <- findInterval(start, dist_AlongTransect, all.inside = TRUE)
		endID <- findInterval(end, dist_AlongTransect, all.inside = TRUE)
		ids1 <- seq_len(startID)
		ids <- startID:endID
		ids_N <- length(ids)
		ids3v1 <- endID:length(veg1_AlongTransect)
		ids3v2 <- endID:length(veg2_AlongTransect)
		ids3a <- endID:length(abut_AlongTransect)
		
		elevC <- ifelse(is.finite(center), elev_AlongTransect[findInterval(center, dist_AlongTransect, all.inside = TRUE)], NA)
		elevS <- elev_AlongTransect[startID]
		elevE <- elev_AlongTransect[endID]
		elevMax <- max(elev_AlongTransect[ids])
		elevMin <- min(elev_AlongTransect[ids])
		slope <- abs(elevE - elevS) / width
		slopeD <- atan(slope) * 180/pi

		elevTotal <- sum(abs(elev_AlongTransect[ids][-1] - elev_AlongTransect[ids][-ids_N]))
		slopeRoughness <- (elevTotal - abs(elevE - elevS)) / width
	
		veg1_frac <- sum(veg1_AlongTransect[ids]) / ids_N
		veg2_frac <- sum(veg2_AlongTransect[ids]) / ids_N
		abut_frac <- sum(abut_AlongTransect[ids]) / ids_N
		
		sumy <- function(x) c(mean = mean(x), sd = sd(x))
		
		veg1_below_dens <- sumy(veg1_AlongTransect[ids1])
		veg1_within_dens <- sumy(veg1_AlongTransect[ids])
		veg1_above_dens <- sumy(veg1_AlongTransect[ids3v1])
		
		veg2_below_dens <- sumy(veg2_AlongTransect[ids1])
		veg2_within_dens <- sumy(veg2_AlongTransect[ids])
		veg2_above_dens <- sumy(veg2_AlongTransect[ids3v1])
		
		abut_below_dens <- sumy(abut_AlongTransect[ids1])
		abut_within_dens <- sumy(abut_AlongTransect[ids])
		abut_above_dens <- sumy(abut_AlongTransect[ids3v1])

		lres <- list(center_m = center, start_m = start, end_m = end, width_m = width,
					elev_atCenter_m = elevC, elev_atStart_m = elevS, elev_atEnd_m = elevE, elev_max_m = elevMax, elev_min_m = elevMin, elev_totalGain_m = elevTotal,
					slope_mpm = slope, slope_deg = slopeD, slopeRoughness_mpm = slopeRoughness,
					veg1_fraction = veg1_frac, veg2_fraction = veg2_frac, abut_fraction = abut_frac,
					veg1_below_densMean = veg1_below_dens["mean"], veg1_below_densSD = veg1_below_dens["sd"], veg1_within_densMean = veg1_within_dens["mean"], veg1_within_densSD = veg1_within_dens["sd"], veg1_above_densMean = veg1_above_dens["mean"], veg1_above_densSD = veg1_above_dens["sd"],
					veg2_below_densMean = veg2_below_dens["mean"], veg2_below_densSD = veg2_below_dens["sd"], veg2_within_densMean = veg2_within_dens["mean"], veg2_within_densSD = veg2_within_dens["sd"], veg2_above_densMean = veg2_above_dens["mean"], veg2_above_densSD = veg2_above_dens["sd"],
					abut_below_densMean = abut_below_dens["mean"], abut_below_densSD = abut_below_dens["sd"], abut_within_densMean = abut_within_dens["mean"], abut_within_densSD = abut_within_dens["sd"], abut_above_densMean = abut_above_dens["mean"], abut_above_densSD = abut_above_dens["sd"])
	}
	
	lres
}


#' @export
InterZoneLocation <- function(i, b, b_neigh, migtype, etband, etmeasure, flag_bfig, copy_FromMig1_TF, do_figures, dir_fig, ...) {
	#---4. Define location and contact/interaction zone of veg1 and veg2 for each ecological boundary

	dots <- list(...)
	if ("et_methods" %in% names(dots)) {
		et_methods <- dots["et_methods"]
	} else {
		stop("ecotoner::InterZoneLocation(): argument 'et_methods' is missing")
	}

	etmeasure$gETmeas[[migtype]]$etable[b, "Transect_ID"] <- i
	etmeasure$gETmeas[[migtype]]$etable[b, "Neighbor_Cells"] <- b_neigh
	
	step_names <- paste0("Step", stepsHullEdge(ecotoner_settings))

	eB_InterZone[[im]] <- vector(mode = "list", length = 2+length(stepsHullEdge(ecotoner_settings)))
	names(eB_InterZone[[im]]) <- c("Danz2012", "Eppinga2013", paste0("Gastner2009Step", stepsHullEdge(ecotoner_settings)))

	etmeasure$gETmeas[[migtype]] <- vector(mode = "list", length = length(step_lengths))
	step_names <- paste0("Step", step_lengths)
	names(etmeasure$gETmeas[[migtype]]) <- step_names

	if (!copy_FromMig1_TF) {
		ends <- c(NA, NA)
		if (!identical(class(eB_Danz[[im]]$Veg1VsDist_1D$sigmoidal$m), "try-error")) {
			ends[1] <- eB_Danz[[im]]$Veg1VsDist_1D$sigmoidal$center_x
		}
		if (!identical(class(eB_Danz[[im]]$Veg2VsDist_1D$sigmoidal$m), "try-error")) {
			ends[2] <- eB_Danz[[im]]$Veg2VsDist_1D$sigmoidal$center_x
		}
		eB_InterZone[[im]]$Danz2012 <- calc_InterZoneLocation(center = NA, start = ifelse(any(is.na(ends)), NA, min(ends)), end = ifelse(any(is.na(ends)), NA, max(ends)), dist_AlongTransect = etband$Env$DistAlongXaxis_m, elev_AlongTransect = etband$Env$elev$YMeans_ForEachX, veg1_AlongTransect = etband$Veg[[im]]$Veg1$density, veg2_AlongTransect = etband$Veg[[im]]$Veg2$density, abut_AlongTransect = etband$Veg[[im]]$Veg1ABUTVeg2$density)
		ends <- c((opt <- eB_Eppinga[[im]]$optim$pos_m) + eB_Eppinga[[im]]$adv_veg1$deltaFrontRunners_Mean_m, opt + eB_Eppinga[[im]]$adv_veg2$deltaFrontRunners_Mean_m)
		eB_InterZone[[im]]$Eppinga2013 <- calc_InterZoneLocation(center = opt, start = ifelse(any(is.na(ends)), NA, min(ends)), end = ifelse(any(is.na(ends)), NA, max(ends)), dist_AlongTransect = etband$Env$DistAlongXaxis_m, elev_AlongTransect = etband$Env$elev$YMeans_ForEachX, veg1_AlongTransect = etband$Veg[[im]]$Veg1$density, veg2_AlongTransect = etband$Veg[[im]]$Veg2$density, abut_AlongTransect = etband$Veg[[im]]$Veg1ABUTVeg2$density)
		for (i_step in seq_along(stepsHullEdge(ecotoner_settings))) {
			istep <- step_names[i_step]
			ends <- c(eB_Gastner[[im]][[istep]]$Veg1$stats$position_m, eB_Gastner[[im]][[istep]]$Veg2$stats$position_m)
			eB_InterZone[[im]][[2+i_step]] <- calc_InterZoneLocation(center = NA, start = ifelse(any(is.na(ends)), NA, min(ends)), end = ifelse(any(is.na(ends)), NA, max(ends)), dist_AlongTransect = etband$Env$DistAlongXaxis_m, elev_AlongTransect = etband$Env$elev$YMeans_ForEachX, veg1_AlongTransect = etband$Veg[[im]]$Veg1$density, veg2_AlongTransect = etband$Veg[[im]]$Veg2$density, abut_AlongTransect = etband$Veg[[im]]$Veg1ABUTVeg2$density)
		}
		
		if (do_figures) plot_interaction_zones(filename = file.path(dir_fig, paste0(flag_bfig, "InteractionZone_Comparison", migtype, ".pdf")), datLine = etline, eB_Env = etband$Env, eB_Veg = etband$Veg[[im]], datFit = eB_InterZone[[im]])
		etable <- tabulate_InterZoneLocation(etable = etable, b = b, data = eB_InterZone[[im]], colFlag = "InterZone", flag_migtype = migtype)
	} else {
		etable <- tabulate_InterZoneLocation(etable = etable, b = b, data = eB_InterZone[[1]], colFlag = "InterZone", flag_migtype = migtype)
	}
	
	template_etobs$gETmeas[[b]] <- eB_InterZone
	template_etobs$etable[b] <- etable
	
	template_etobs
}

