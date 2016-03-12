#------Eppinga, M.B., Pucko, C.A., Baudena, M., Beckage, B. & Molofsky, J. (2013) A new method to infer vegetation boundary movement from ‘snapshot’ data. Ecography, 36, 622-635.

#---Eppinga et al. 2013: 'Analytical analysis of vegetation boundary movement' (Fig. 1)
calc_Eppinga2013_optpos <- function(x, Veg1, Veg2){
	#in our cases, d(veg1) + d(veg2) <= 1
	#here, attempt 2: weigh each column by total column density to account for uneven distributed empty space
	mat1 <- raster::as.matrix(Veg1, maxpixels = raster::ncell(Veg1))
	mat2 <- raster::as.matrix(Veg2, maxpixels = raster::ncell(Veg2))
	matT <- raster::as.matrix(temp <- raster::overlay(Veg1, Veg2, fun = function(v1, v2) ifelse(!is.na(v1) | !is.na(v2), 1, NA)), maxpixels = raster::ncell(temp))

	sum1 <- sum(!is.na(mat1))
	sum2 <- sum(!is.na(mat2))
	cumsum.ColNotEmpty <- cumsum(apply(matT, 2, function(x) sum(!is.na(x))))
	sumT <- max(cumsum.ColNotEmpty)
	
	opt <- which(sum1 <= cumsum.ColNotEmpty)[1]
	
	list(pos_cell = opt, pos_m = max(x) / length(x) * opt,
		dr1 = sum1 / sumT, dr2 = sum2 / sumT,
		totd = sumT / raster::ncell(Veg1))
}

#---Eppinga et al. 2013: 'Analytical analysis of vegetation boundary movement' (Fig. 1)
calc_Eppinga2013_frontrunners <- function(veg, end_toLeft){
	mat <- raster::as.matrix(veg, maxpixels = raster::ncell(veg))
	fun <- if (end_toLeft) min else max
	xfront <- apply(mat, 1, function(v) ifelse(all(temp <- is.na(v)), NA, fun(which(!temp))))
	xfront[is.infinite(xfront)] <- NA
	xfront
}



#---Eppinga et al. 2013: 'Analytical analysis of vegetation boundary movement' and 'Statistical analyses'
calc_Eppinga2013_advancement <- function(x, veg, end_toLeft, optBoundary_m) {
	deltaF_m <- max(x) / length(x) * (calc_Eppinga2013_frontrunners(veg = veg, end_toLeft = end_toLeft)) - optBoundary_m
	TdeltaF <- transformation17(deltaF_m)
	bt_fmean <- backtransformation17(tm <- mean(TdeltaF, na.rm = TRUE))
	tsd <- sd(TdeltaF, na.rm = TRUE)
		
	list(deltaFrontRunners_T17 = TdeltaF, deltaFrontRunners_m = deltaF_m,
		deltaFrontRunners_Mean_m = bt_fmean,
		deltaFrontRunners_Mean_T17 = tm, deltaFrontRunners_SD_T17 = tsd)
}

#' @export
map_front_runners_Eppinga2013 <- function(filename, eB_Env, eB_Veg, datFit) {
	pdf(width=7, height=7, file=filename)
	par_old <- par(mfrow=c(1, 1), mar=c(0, 0.1, 1, 1), mgp=c(2, 0.5, 0), cex=cex <- 1.5)
	on.exit({par(par_old); dev.off()}, add = TRUE)

	ext1 <- raster::extent(eB_Env$elev$grid)
	xlim <- c(-1000, ext1@xmax)
	ylim <- c(-1000, ext1@ymax+1000)
	
	raster::image(eB_Env$elev$grid, col=gray(0:255/255), xlim=xlim, ylim=ylim, main="", xlab="", ylab="", asp=1, axes=FALSE)
	raster::image(eB_Veg$Veg1$grid, col=adjustcolor("red", alpha.f = 0.3), add=TRUE)
	raster::image(eB_Veg$Veg2$grid, col=adjustcolor("darkgreen", alpha.f = 0.3), add=TRUE)
	atx <- c((atx <- axTicks(1))[atx >= 0 & atx < ext1@xmax], ext1@xmax)
	axis(1, pos=0, at=atx)
	axis(2, pos=0, at=c(0, 2000, 4000, 6000))
	text(x=ext1@xmax/2, y=-strheight("0", units="user", cex=cex)*(0.5+2), labels="Transect length (m)")
	text(x=-strwidth("0", units="user", cex=cex)*(0.5+2.5), y=ext1@ymax/2, labels="Transect width (m)", srt=90)

	isnotna <- !is.na(datFit$adv_veg1$deltaFrontRunners_m)
	res_m <- raster::xres(eB_Env$elev$grid)
	points(x=((opt <- datFit$optim$pos_m) + datFit$adv_veg1$deltaFrontRunners_m)[isnotna], y=(ys <- (length(datFit$adv_veg1$deltaFrontRunners_m):1) * res_m)[isnotna], pch=".", col="magenta")
	isnotna <- !is.na(datFit$adv_veg2$deltaFrontRunners_m)
	points(x=(opt + datFit$adv_veg2$deltaFrontRunners_m)[isnotna], y=ys[isnotna], pch=".", col="green")
	segments(x0=opt, y0=ext1@ymin, x1=opt, y1=ext1@ymax, lwd=2, col="yellow")
	segments(x0=pos1 <- opt + datFit$adv_veg1$deltaFrontRunners_Mean_m, y0=ext1@ymin, x1=pos1, y1=ext1@ymax, lwd=2, col=adjustcolor("magenta", alpha.f = 0.7))
	segments(x0=pos2 <- opt + datFit$adv_veg2$deltaFrontRunners_Mean_m, y0=ext1@ymin, x1=pos2, y1=ext1@ymax, lwd=2, col=adjustcolor("green", alpha.f = 0.7))

	if (datFit$adv_stats$FrontsAdvBeyondOptBoundary) {
		p12 <- datFit$adv_stats$Front1Larger2_p < 0.05
		p21 <- datFit$adv_stats$Front2Larger1_p < 0.05
		if (p12) {
			arrows(x0=opt, y0=ext1@ymax+30, x1=pos1, y1=ext1@ymax+30, lwd=2, col="magenta")
			ptext <- paste0("adv(Veg1; BSE) > adv(Veg2; TF):\np < ", ifelse(datFit$adv_stats$Front1Larger2_p > 0.001, formatC(datFit$adv_stats$Front1Larger2_p, format="f", digits=3), "0.001"))
		}
		if (p21) {
			arrows(x0=opt, y0=ext1@ymax+30, x1=pos2, y1=ext1@ymax+30, lwd=2, col="green")
			ptext <- paste0("adv(Veg2; TF) > adv(Veg1; BSE):\np < ", ifelse(datFit$adv_stats$Front2Larger1_p > 0.001, formatC(datFit$adv_stats$Front2Larger1_p, format="f", digits=3), "0.001"))
		}
		if (!p12 & !p21) ptext <- paste0("adv(Veg1; BSE) = adv(Veg2; TF):\np > ", formatC(datFit$adv_stats$Front1Larger2_p, format="f", digits=3))
	} else {
		ptext <- paste0("adv() < optimal boundary")
	}
	text(x=ext1@xmax/2, y=ext1@ymax+strheight("0", units="user", cex=2/3*cex)*2, labels=ptext, cex=2/3)
	
	invisible()
}




#---Eppinga et al. 2013: 'Statistical analyses'
calc_Eppinga2013_stats <- function(deltaF1_T17, deltaF2_T17, seed = NULL) {
	#assumption that all frontrunners are y-row wise dispersed and do not originate from other sources
	advanced <- all(mean(deltaF1_T17, na.rm = TRUE) > 0, mean(deltaF2_T17, na.rm = TRUE) > 0)
	res <- list(FrontsAdvBeyondOptBoundary = advanced, FrontDiff_Mean_T17 = NA, FrontDiff_SE_T17 = NA,
			Front1Larger2_p = NA, Front2Larger1_p = NA,
			Front1Larger2_power = NA)
	
	if (advanced) {
		#Test if Veg1 advancement is larger than Veg2's: eq. 18
		if (!is.na(seed)) set.seed(seed)
		
		if (requireNamespace("simpleboot", quietly = TRUE) && requireNamespace("boot", quietly = TRUE)) {
			b <- simpleboot::two.boot(sample1 = deltaF1_T17, sample2 = deltaF2_T17, FUN = mean, R = R <- 10000, na.rm = TRUE)
		} else {
			warning("Package 'simpleboot' and/or 'boot' not installed: 'calc_Eppinga2013_stats' will not be estimated")
			return(res)
		}
		
		FrontDiff_T17.mean <- mean(b$t[,1], na.rm = TRUE) #added na.rm = TRUE, in case one vegetation type has no occurrence in some rows
		n <- sum(!is.na(deltaF1_T17 - deltaF2_T17))
		FrontDiff_T17.se <- sqrt(n / (n - 1)) * sd(b$t[,1], na.rm = TRUE) #standard error of the mean = f(n) * standard deviation of the bootstrap means
		FrontDiff_T17.sd <- FrontDiff_T17.se * sqrt(n)
		if (FrontDiff_T17.mean >= 0) { #Advance of veg1 is larger than of veg2; if both are advancing
			TAdv1LargerAdv2.p <- sum(Heaviside(abs(b$t) + abs(FrontDiff_T17.mean) - abs(b$t + FrontDiff_T17.mean))) / R
			TAdv2LargerAdv1.p <- 1
		} else {
			TAdv1LargerAdv2.p <- 1
			TAdv2LargerAdv1.p <- sum(Heaviside(abs(-b$t) + abs(-FrontDiff_T17.mean) - abs(-b$t - FrontDiff_T17.mean))) / R
		}
	
		#Power test: eq. 19
		tau <- 0.2 #effect size
		tcrit <- qt(0.95, df = n) #95% confidence
		TAdv1LargerAdv2.power <- 1 - (1/2 * (1 + erf((tcrit - tau * sqrt(n) / FrontDiff_T17.sd) / (sqrt(2) * FrontDiff_T17.sd))))

		res <- list(FrontsAdvBeyondOptBoundary = advanced, FrontDiff_Mean_T17 = FrontDiff_T17.mean, FrontDiff_SE_T17 = FrontDiff_T17.se,
				Front1Larger2_p = TAdv1LargerAdv2.p, Front2Larger1_p = TAdv2LargerAdv1.p,
				Front1Larger2_power = TAdv1LargerAdv2.power)
	}
		
	res
}

tabulate_Eppinga2013_advance <- function(etable, b, data, flag_migtype){
	colnamesAdd <- paste0(flag_migtype, "_Eppinga2013_",
						c("OptimalPosition_AlongXaxis_m",
						"Veg1_DeltaFront_Mean_m", "Veg1_DeltaFront_Mean_T17", "Veg1_DeltaFront_SD_T17", 
						"Veg2_DeltaFront_Mean_m", "Veg2_DeltaFront_Mean_T17", "Veg2_DeltaFront_SD_T17", 
						"FrontsAdvBeyondOptBoundary", "FrontDiff_Mean_T17", "FrontDiff_SE_T17", 
						"Front1Larger2_p", "Front2Larger1_p", "Front1Larger2_power"))
	res <- matrix(NA, nrow = max(1, nrow(etable)), ncol = length(colnamesAdd), dimnames = list(NULL, colnamesAdd))

	res[b, 1] <- data$optim$pos_m
	res[b, 2] <- data$adv_veg1$deltaFrontRunners_Mean_m
	res[b, 3] <- data$adv_veg1$deltaFrontRunners_Mean_T17
	res[b, 4] <- data$adv_veg1$deltaFrontRunners_SD_T17
	res[b, 5] <- data$adv_veg2$deltaFrontRunners_Mean_m
	res[b, 6] <- data$adv_veg2$deltaFrontRunners_Mean_T17
	res[b, 7] <- data$adv_veg2$deltaFrontRunners_SD_T17
	res[b, 8] <- data$adv_stats$FrontsAdvBeyondOptBoundary
	res[b, 9] <- data$adv_stats$FrontDiff_Mean_T17
	res[b, 10] <- data$adv_stats$FrontDiff_SE_T17
	res[b, 11] <- data$adv_stats$Front1Larger2_p
	res[b, 12] <- data$adv_stats$Front2Larger1_p
	res[b, 13] <- data$adv_stats$Front1Larger2_power
	
	cbind(etable, res)
}




#' @export
Eppinga2013Ecography <- function(i, b, migtype, ecotoner_settings, etband, etmeasure, flag_bfig, copy_FromMig1_TF, do_figures, ...) {
	#3b. Eppinga et al. 2013 Ecography: Location of boundary and front-runner distance
	#Objective: inference of vegetation boundary movement from one ‘snapshot’ (e.g. an aerial photograph or satellite image) in time
	#It is assumed that the current vegetation distribution is reflecting competitive interactions between communities over a longer time period (i.e. decades). Also, it is assumed that vegetation boundary movement is relatively slow as compared to fluctuations in environmental and meteorological conditions.
	#To meet these assumptions when applying the method to real ecosystems, we select snapshots of vegetation bound- aries that: 1) consist of two discrete communities that are characterized by plants with a relatively long lifespan (mostly trees and shrubs); 2) have been described in the literature, meaning that the direction of vegetation boundary move- ment is known; 3) are considered to be primarily driven by competitive interactions. 
	
	dots <- list(...)
	seed <- if ("seed" %in% names(dots)) dots["seed"] else NULL

	etmeasure$etable[b, "Transect_ID"] <- i
	etmeasure$etable[b, "Neighbor_Cells"] <- neighborhoods(ecotoner_settings)[b]
	
	if (!copy_FromMig1_TF) {
		etmeasure$gETmeas[[b]][[migtype]]$optim <- calc_Eppinga2013_optpos(x = etband$Env$DistAlongXaxis_m,
																Veg1 = etband$Veg[[migtype]]$Veg1$grid,
																Veg2 = etband$Veg[[migtype]]$Veg2$grid)
		etmeasure$gETmeas[[b]][[migtype]]$adv_veg1 <- calc_Eppinga2013_advancement(x = etband$Env$DistAlongXaxis_m,
																		veg = etband$Veg[[migtype]]$Veg1$grid,
																		end_toLeft = end_to_left(type_veg1(ecotoner_settings)),
																		optBoundary_m = etmeasure$gETmeas[[b]][[migtype]]$optim$pos_m)
		etmeasure$gETmeas[[b]][[migtype]]$adv_veg2 <- calc_Eppinga2013_advancement(x = etband$Env$DistAlongXaxis_m,
																		veg = etband$Veg[[migtype]]$Veg2$grid,
																		end_toLeft = end_to_left(type_veg2(ecotoner_settings)),
																		optBoundary_m = etmeasure$gETmeas[[b]][[migtype]]$optim$pos_m)
		etmeasure$gETmeas[[b]][[migtype]]$adv_stats <- calc_Eppinga2013_stats(deltaF1_T17 = etmeasure$gETmeas[[b]][[migtype]]$adv_veg1$deltaFrontRunners_T17,
																deltaF2_T17 = -etmeasure$gETmeas[[b]][[migtype]]$adv_veg2$deltaFrontRunners_T17,
																seed = seed)

		if(do_figures) map_front_runners_Eppinga2013(filename = file.path(dir_fig, paste0(flag_bfig, "Eppinga2013_FittedBoundaryAdvancement_", migtype, ".pdf")),
													eB_Env = etband$Env,
													eB_Veg = etband$Veg[[migtype]],
													datFit = etmeasure$gETmeas[[b]][[migtype]])				
	}
	
	etmeasure$etable <- tabulate_Eppinga2013_advance(etable = etmeasure$etable, b = b,
													data = etmeasure$gETmeas[[b]][[if (copy_FromMig1_TF) 1 else migtype]],
													flag_migtype = migtype)
		
	etmeasure
}
