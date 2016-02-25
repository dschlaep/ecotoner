#' Definition: Ecological Boundary Zone
#'
#' @export
determine_ecotone_linear_extent <- function(dens_BSE, dens_Veg1, end1_toLeft, dens_TF, dens_Veg2, end2_toLeft, veg_density_low, veg_density_high, veg_density_extended_min) {
	cropUnion <- croptemp <- c(NA, NA)
	rmBSE <- rmTF <- NULL
	if (all(is.null(dens_BSE), is.null(dens_TF), !is.null(dens_Veg1), !is.null(dens_Veg2))) {
		cv1 <- calc_veg_decrease(dens_Veg1, end1_toLeft, veg_density_low, veg_density_high, veg_density_extended_min)
		cv2 <- calc_veg_decrease(dens_Veg2, end2_toLeft, veg_density_low, veg_density_high, veg_density_extended_min)
	}
	if (all(!is.null(dens_BSE), !is.null(dens_TF), is.null(dens_Veg1), is.null(dens_Veg2))) {
		cv1 <- calc_veg_decrease(dens_BSE, end1_toLeft, veg_density_low, veg_density_high, veg_density_extended_min)
		cv2 <- calc_veg_decrease(dens_TF, end2_toLeft, veg_density_low, veg_density_high, veg_density_extended_min)
	}
	if (all(!is.null(dens_BSE), !is.null(dens_TF), !is.null(dens_Veg1), !is.null(dens_Veg2))) {
		cv1 <- calc_veg_decrease(dens_Veg1, end1_toLeft, veg_density_low, veg_density_high, veg_density_extended_min)
		cv2 <- calc_veg_decrease(dens_Veg2, end2_toLeft, veg_density_low, veg_density_high, veg_density_extended_min)
		#no remaining BSE or TF vegetation should be more than veg_density_extended_min
		rmBSE <- calc_veg_decrease(dens_BSE-dens_Veg1, end1_toLeft, veg_density_low, veg_density_high, veg_density_extended_min)$rmain
		rmTF <- calc_veg_decrease(dens_TF-dens_Veg2, end2_toLeft, veg_density_low, veg_density_high, veg_density_extended_min)$rmain
	}
	
	
	if(all(!is.na(c(cv1$szone, cv1$ezone, cv2$szone, cv2$ezone)))) { #density pattern exists within type1 and type2
		#determine overall start:
		#	- leftmost position within minimum zone of veg2 such that
		#	- there is at least one point of the maximum zone of veg1 included (i.e., to the right) and that
		#	- all of veg1 between overall start and the rightmost point of the maximum zone of veg1 have sufficient density (i.e., > veg_density_extended_min)
		veg1max_rightmost <- tail(cv1$szone, n=1)
		tempj <- findInterval(veg1max_rightmost, temps <- c(0, cumsum(cv1$rmain$lengths)), all.inside=TRUE)
		if (cv1$rmain$values[tempj]) {
			veg1main_leftmost <- temps[tempj] + 1
			startpos <- intersect(cv2$szone, veg1main_leftmost:veg1max_rightmost)
			if (length(startpos) > 0) croptemp[1] <- startpos[1]
		}
		
		#determine overall end:
		#	- rightmost position within minimum zone of veg1 such that
		#	- overall start is on the left,
		#	- there is at least one point of the maximum zone of veg2 included (i.e., to the left) and that
		#	- all of veg2 between the lefmost point of the maximum zone of veg2 and overall end have sufficient density (i.e., > veg_density_extended_min)
		if (length(temp <- which(cv2$ezone > croptemp[1])) > 0) {
			veg2max_leftmost <- cv2$ezone[temp[1]]
			tempj <- findInterval(veg2max_leftmost, temps <- c(0, cumsum(cv2$rmain$lengths)), all.inside=TRUE)
			if (cv2$rmain$values[tempj]) {
				veg2main_rigthmost <- temps[tempj+1]
				endpos <- intersect(cv1$ezone, veg2max_leftmost:veg2main_rigthmost)
				if (length(endpos) > 0) croptemp[2] <- tail(endpos, n=1)
			}
		}
	}
	
	if (all(!is.na(croptemp)) && (any(is.null(rmBSE), is.null(rmTF)) || all(!is.null(rmBSE), !is.null(rmTF), !rmBSE$values, !rmTF$values))) {
		cropUnion <- croptemp
	}
	
	return(cropUnion)
}

calc_veg_decrease <- function(dens_Veg, end_toLeft, veg_density_low, veg_density_high, veg_density_extended_min, pk = 100) {

	if (requireNamespace("zoo", quietly = TRUE)) {
		ldens <- zoo::rollmean(dens_Veg, k = 3, fill = "extend", align = "center")
	} else {
		temp <- c(rep(dens_Veg[1], 2), dens_Veg, rep(tail(dens_Veg, 1), 2))
		ldens2 <- as.numeric(stats::filter(temp, filter = rep(1/5, 5), method = "convolution", sides = 2))
		ldens2 <- ldens2[-c(1:2, (length(ldens2) - 1):length(ldens2))]
	}
	
	ldensInt <- as.integer(round(ldens*pk)) #having matches of min and max accurate to (1/pk)-percent
	
	istart <- iend <- NA
	rhigh <- rle(dhigh <- ifelse(ldensInt >= veg_density_high * pk, TRUE, FALSE))
	rmain <- rle(ifelse(ldensInt >= veg_density_extended_min * pk, TRUE, FALSE))
	rlow <- rle(dlow <- ifelse(ldensInt <= veg_density_low * pk, TRUE, FALSE))
	#rhigh and rlow can be FALSE at the same locations, but they can never be TRUE at the same locations
	
	if (any(rhigh$values) && any(rlow$values)) {
		if (end_toLeft) {
			#first run of low values
			temps <- c(0, cumsum(rlow$lengths))	#indices of ends of low runs
			tempi <- which(rlow$values)[1]	#index of first low run in rle-object 
			istart <- (1 + temps[tempi]):temps[tempi + 1]
			#first high run after the last low run
			temps <- c(0, cumsum(rhigh$lengths))	#indices of ends of high runs
			tempj <- findInterval(max(istart), temps, all.inside=TRUE)	#index of rightmost low of last run among high runs; this will always be a high run of FALSE
			if (tempj < length(rhigh$values)) { #if not true, then low run is at the very right end of transect and no high run could be found
				tempi <- tempj + which(rhigh$values[(tempj + 1):length(rhigh$values)])[1]
				iend <- (1 + temps[tempi]):temps[tempi + 1]
			}
		} else {
			#last run of high values
			temps <- c(0, cumsum(rhigh$lengths))	#indices of ends of high runs
			tempi <- tail(which(rhigh$values), n=1)	#index of last high run in rle-object
			istart <- (1 + temps[tempi]):temps[tempi + 1]
			#last low run after the last high series
			temps <- c(0, cumsum(rlow$lengths))	#indices of ends of low runs
			tempj <- findInterval(max(istart), temps, all.inside=TRUE)	#index of rightmost high of last run among low runs; this will always be a low run of FALSE
			if (tempj < length(rlow$values)) { #if not true, then high run is at the very right end of transect and no low run could be found
				tempi <- tempj + tail(which(rlow$values[(tempj + 1):length(rlow$values)]), n=1)
				iend <- (1 + temps[tempi]):temps[tempi + 1]
			}
		}
	}
	
	list(szone = istart, ezone = iend, rmain = rmain)
}

