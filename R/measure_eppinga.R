#------Eppinga, M.B., Pucko, C.A., Baudena, M., Beckage, B. & Molofsky, J. (2013) A new method to infer vegetation boundary movement from ‘snapshot’ data. Ecography, 36, 622-635.

#---Eppinga et al. 2013: 'Analytical analysis of vegetation boundary movement' (Fig. 1)
calc_Eppinga2013_optpos <- function(Veg1, Veg2){
	# in our cases, d(veg1) + d(veg2) <= 1
	# determine vegetation cover (cells that are not NA)
	ncells <- raster::ncell(Veg1)
	cov1 <- ncells - raster::cellStats(Veg1, "countNA")
	VegT <- raster::overlay(Veg1, Veg2, fun = function(v1, v2) ifelse(!is.na(v1) | !is.na(v2), 1, 0))
	covT <- raster::cellStats(VegT, "sum")
	
	# find location where veg1 would end if all of veg1 were arranged to the left and if transect width was reduced according to empty cells
	opt <- cov1 * raster::ncol(VegT)/ covT
	
	
	if (FALSE) {
		# plot illustrating how I adjusted the method for grids with 'empty' cells
		mat1 <- raster::as.matrix(Veg1, maxpixels = ncells)
		matT <- raster::as.matrix(VegT, maxpixels = ncells)
		
		add_transect_frame <- function(matT, res = 1) {
			axis(side = 1, pos = 0, at = c(axTicks(1), res * ncol(matT)))
			axis(side = 2, pos = 0, at = c(axTicks(2), res * nrow(matT)))
			rect(0, 0, res * ncol(matT), res* nrow(matT), lwd = 1)
			
			invisible(NULL)
		}	

		pdf(height = 4, width = 12, file = file.path(dir_fig, paste0(flag_bfig, "Eppinga2013_OptimalBoundaryLocation_", migtype, ".pdf")))
		op_old <- par(mfrow = c(1, 3), mar = c(2.5, 2.5, 1, 0.1), mgp = c(1.15, 0.15, 0), tcl = 0.5, cex = cex <- 1)
		
		# (a) map of transect band
		raster::image(Veg1, col = "red", main = "", xlab = "Transect length (m)", ylab = "Transect width (m)", asp = 1, axes = FALSE)
		raster::image(Veg2, col = "darkgreen", add = TRUE)
		add_transect_frame(matT, res = raster::xres(Veg1))
		mtext(side = 3, line = -1, text = "(a)", font = 2, adj = -0.075, cex = cex)
		
		# (b) vegetation separated from 'empty' cells
		ctveg <- apply(matT, 2, sum)
		cveg1 <- apply(mat1, 2, sum, na.rm = TRUE)
		
		barplot(ctveg, width = 1, space = 0, border = NA, col = "darkgreen", axes = FALSE, asp = 1, xlim = c(1, ncol(matT)), ylim = c(1, nrow(matT)), xlab = "Transect length (# column)", ylab = "Transect width (# row)")
		barplot(cveg1, width = 1, space = 0, border = NA, col = "red", axes = FALSE, add = TRUE)
		add_transect_frame(matT)
		mtext(side = 3, line = -1, text = "(b)", font = 2, adj = -0.075, cex = cex)

		# (c) vegetation optimally organized and 'empty' cells reduced transect width (spread fairly among vegetation types)
		y_empty <- covT / ncol(matT)
		
		plot(1, type = "n", axes = FALSE, asp = 1, xlim = c(1, ncol(matT)), ylim = c(1, nrow(matT)), xlab = "Transect length (# column)", ylab = "")
		rect(0, 0, ncol(matT), y_empty, border = NA, col = "darkgreen")
		rect(0, 0, opt, y_empty, border = NA, col = "red")
		add_transect_frame(matT)
		mtext(side = 3, line = -1, text = "(c)", font = 2, adj = -0.075, cex = cex)
		
		par(op_old)
		dev.off()
	}
	
	list(pos_cell = round(opt), pos_m = raster::xres(Veg1) * opt,
		dr1 = cov1 / covT, dr2 = 1 - cov1 / covT, totd = covT / ncells)
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
calc_Eppinga2013_advancement <- function(veg, end_toLeft, optBoundary_cell) {
	deltaF_m <- raster::xres(veg) * (calc_Eppinga2013_frontrunners(veg, end_toLeft) - optBoundary_cell)
	TdeltaF <- transformation17(deltaF_m)
	tm <- mean(TdeltaF, na.rm = TRUE)
	bt_fmean <- backtransformation17(tm)
		
	list(deltaFrontRunners_T17 = TdeltaF, deltaFrontRunners_m = deltaF_m,
		deltaFrontRunners_Mean_m = backtransformation17(tm),
		deltaFrontRunners_Mean_T17 = tm, deltaFrontRunners_SD_T17 = sd(TdeltaF, na.rm = TRUE))
}



#---Eppinga et al. 2013: 'Statistical analyses'
calc_Eppinga2013_stats <- function(deltaF1_T17, deltaF2_T17, seed = NULL) {
	#assumption that all frontrunners are y-row wise dispersed and do not originate from other sources
	advanced <- all(mean(deltaF1_T17, na.rm = TRUE) > 0, mean(deltaF2_T17, na.rm = TRUE) < 0)
	res <- list(FrontsAdvBeyondOptBoundary = advanced, FrontDiff_Mean_T17 = NA, FrontDiff_var_T17 = NA, advancement_p = NA, advancement_retro_power = NA)
	
	if (advanced) {
		if (requireNamespace("boot", quietly = TRUE)) {
			if (!is.na(seed)) set.seed(seed)
			
			# Test if vegetation types have same distance of front runners, i.e., test if difference is 0: eq. 18
			R <- 1e5L
			data <- cbind(deltaF1_T17, -deltaF2_T17)
			bmd <- boot::boot(data, boot_mean_of_diffs, R = R, stype = "i", sim = "ordinary", parallel = "no")
			
			#bmd_ci <- boot::boot.ci(bmd, conf = 0.95, type = "bca") # Test: is 0 contained in ci?
			
			# Calculate p-value (for H0: diff = 0); eq. 19
			# Direct interpretation of eq. 19: sum(Heaviside(abs(bmd$t) + abs(bmd$t0) - abs(bmd$t + bmd$t0))) / bmd$R
			ptol <- sqrt(.Machine$double.eps)
			ntol <- -sqrt(.Machine$double.neg.eps)
			advancement_p <- (if (bmd$t0 > ptol) sum(bmd$t <= ptol) else if (bmd$t0 < ntol) sum(bmd$t >= ntol) else bmd$R) / bmd$R

			# Retrospective power: eq. 20
			FrontDiff_T17_var <- as.numeric(var(bmd$t))
			tau <- 0.2 #effect size
			n <- sum(complete.cases(data))
			tcrit <- qt(0.95, df = n) #95% confidence
			advancement_retro_power <- 1 - (1/2 * (1 + erf((tcrit - tau * sqrt(n) / sqrt(FrontDiff_T17_var)) / (sqrt(2 * FrontDiff_T17_var)))))

			res <- list(FrontsAdvBeyondOptBoundary = advanced,
					FrontDiff_Mean_T17 = bmd$t0, FrontDiff_var_T17 = FrontDiff_T17_var,
					FrontDiff_Mean_m = backtransformation17(bmd$t0),
					boots_R = R, advancement_p = advancement_p, advancement_retro_power = advancement_retro_power)
		} else {
			warning("Package 'boot' not installed: 'calc_Eppinga2013_stats' will not be estimated")
		}
	}

	res
}


#---Output
#' @export
map_front_runners_Eppinga2013 <- function(filename, eB_Env, eB_Veg, datFit) {
#TODO(drs): make this function more general, e.g., axis(side = 2, ...) uses case-specific constants
	res_m <- raster::xres(eB_Env$elev$grid)

	pdf(width = 7, height = 7, file = filename)
	par_old <- par(mfrow = c(1, 1), mar = c(0, 0.1, 1, 1), mgp = c(2, 0.5, 0), cex = cex <- 1.5)
	on.exit({par(par_old); dev.off()}, add = TRUE)

	ext1 <- raster::extent(eB_Env$elev$grid)
	xlim <- c(-1000, ext1@xmax)
	ylim <- c(-1000, ext1@ymax + 1000)
	
	raster::image(eB_Env$elev$grid, col = gray(0:255/255), xlim = xlim, ylim = ylim, main = "", xlab = "", ylab = "", asp = 1, axes = FALSE)
	raster::image(eB_Veg$Veg1$grid, col = adjustcolor("red", alpha.f = 0.3), add = TRUE)
	raster::image(eB_Veg$Veg2$grid, col = adjustcolor("darkgreen", alpha.f = 0.3), add = TRUE)
	atx <- c((atx <- axTicks(1))[atx >= 0 & atx < ext1@xmax], ext1@xmax)
	axis(1, pos = 0, at = atx)
	axis(2, pos = 0, at = c(0, 2000, 4000, 6000))
	text(x = ext1@xmax/2, y = -strheight("0", units = "user", cex = cex) * (0.5 + 2), labels = "Transect length (m)")
	text(x = -strwidth("0", units = "user", cex = cex) * (0.5 + 2.5), y = ext1@ymax/2, labels = "Transect width (m)", srt = 90)

	# Optimal boundary
	segments(x0 = datFit$optim$pos_m, y0 = ext1@ymin, x1 = datFit$optim$pos_m, y1 = ext1@ymax, lwd = 2, col = "yellow")

	# Front runners
	ys <- res_m * (raster::nrow(eB_Env$elev$grid):1) - res_m / 2
	copt <- res_m * datFit$optim$pos_cell - res_m / 2
	isnotna <- !is.na(datFit$adv_veg1$deltaFrontRunners_m)
	points(x = (copt + datFit$adv_veg1$deltaFrontRunners_m)[isnotna], y = ys[isnotna], pch = 46, col = "magenta")
	isnotna <- !is.na(datFit$adv_veg2$deltaFrontRunners_m)
	points(x = (copt + datFit$adv_veg2$deltaFrontRunners_m)[isnotna], y = ys[isnotna], pch = 46, col = "green")

	if (datFit$adv_stats$FrontsAdvBeyondOptBoundary) {
		padv <- datFit$adv_stats$advancement_p < 0.05 && datFit$adv_stats$advancement_retro_power > 0.8
		p12 <- padv && datFit$adv_stats$FrontDiff_Mean_T17 > 0
		p21 <- padv && datFit$adv_stats$FrontDiff_Mean_T17 < 0
		p_str <- ifelse(datFit$adv_stats$advancement_p > 0.001, formatC(datFit$adv_stats$advancement_p, format = "f", digits = 3), "0.001")
		
		if (padv && abs(datFit$adv_stats$FrontDiff_Mean_m) > 2 * res_m) {
			# Mean front runner advancement
			pos1 <- copt + datFit$adv_veg1$deltaFrontRunners_Mean_m
			segments(x0 = pos1, y0 = ext1@ymin, x1 = pos1, y1 = ext1@ymax, lwd = 2, col = adjustcolor("magenta", alpha.f = 0.7))
			pos2 <- copt + datFit$adv_veg2$deltaFrontRunners_Mean_m
			segments(x0 = pos2, y0 = ext1@ymin, x1 = pos2, y1 = ext1@ymax, lwd = 2, col = adjustcolor("green", alpha.f = 0.7))
			arrows(x0 = copt, y0 = ext1@ymax + res_m,
				   x1 = if (p12) pos1 else if (p21) pos2 else 0, y1 = ext1@ymax + res_m,
				   lwd = 2, col = if (p12) "magenta" else if (p21) "green" else "white")
		}

		ptext <- paste0("adv(Veg1 = ", signif(datFit$adv_veg1$deltaFrontRunners_Mean_m, 2), " m) ",
						if (p12) ">" else if (p21) "<" else "=",
						" adv(Veg2 = ", signif(-datFit$adv_veg2$deltaFrontRunners_Mean_m, 2), " m):",
						"\np ", if (datFit$adv_stats$advancement_p < 0.05) "<" else ">=", " ", p_str,
						" based on ", datFit$adv_stats$boots_R, " bootstrap replicates,",
						"\nand with a retrospective power of ", signif(datFit$adv_stats$advancement_retro_power, 3))


	} else {
		ptext <- paste0("no advancement beyond optimal boundary location")
	}
	text(x = ext1@xmax / 2, y = ext1@ymax + 3 * strheight("0", units = "user", cex = 2 / 3 * cex), labels = ptext, cex = 2/3)
	
	invisible()
}


tabulate_Eppinga2013_advance <- function(etable, index, data){
	colnamesAdd <- paste0("Eppinga2013_",
						c("OptimalPosition_AlongXaxis_m",
						"Veg1_DeltaFront_Mean_m", "Veg1_DeltaFront_Mean_T17", "Veg1_DeltaFront_SD_T17", 
						"Veg2_DeltaFront_Mean_m", "Veg2_DeltaFront_Mean_T17", "Veg2_DeltaFront_SD_T17", 
						"FrontsAdvBeyondOptBoundary", "FrontDiff_Mean_T17", "FrontDiff_var_T17", "FrontDiff_Mean_m",
						"Bootstrap_reps", "advancement_p", "advancement_retrospective_power"))
	res <- vector("numeric", length(colnamesAdd))
	names(res) <- colnamesAdd
	
	res[1] <- data$optim$pos_m
	res[2] <- data$adv_veg1$deltaFrontRunners_Mean_m
	res[3] <- data$adv_veg1$deltaFrontRunners_Mean_T17
	res[4] <- data$adv_veg1$deltaFrontRunners_SD_T17
	res[5] <- data$adv_veg2$deltaFrontRunners_Mean_m
	res[6] <- data$adv_veg2$deltaFrontRunners_Mean_T17
	res[7] <- data$adv_veg2$deltaFrontRunners_SD_T17
	res[8] <- data$adv_stats$FrontsAdvBeyondOptBoundary
	res[9] <- data$adv_stats$FrontDiff_Mean_T17
	res[10] <- data$adv_stats$FrontDiff_var_T17
	res[11] <- data$adv_stats$FrontDiff_Mean_m
	res[12] <- data$adv_stats$boots_R
	res[13] <- data$adv_stats$advancement_p
	res[14] <- data$adv_stats$advancement_retro_power
	
	tabulate_merge_into_etable(etable, index, res)
}




#' @export
Eppinga2013Ecography <- function(i, b, migtype, ecotoner_settings, etband, etmeasure, copy_FromMig1_TF, do_figures, ...) {
	#3b. Eppinga et al. 2013 Ecography: Location of boundary and front-runner distance
	#Objective: inference of vegetation boundary movement from one ‘snapshot’ (e.g. an aerial photograph or satellite image) in time
	#It is assumed that the current vegetation distribution is reflecting competitive interactions between communities over a longer time period (i.e. decades). Also, it is assumed that vegetation boundary movement is relatively slow as compared to fluctuations in environmental and meteorological conditions.
	#To meet these assumptions when applying the method to real ecosystems, we select snapshots of vegetation bound- aries that: 1) consist of two discrete communities that are characterized by plants with a relatively long lifespan (mostly trees and shrubs); 2) have been described in the literature, meaning that the direction of vegetation boundary move- ment is known; 3) are considered to be primarily driven by competitive interactions. 
	
	dots <- list(...)
	seed <- if ("seed" %in% names(dots)) dots["seed"] else NULL
	if ("flag_bfig" %in% names(dots)) flag_bfig <- dots["flag_bfig"] else do_figures <- FALSE
	if ("dir_fig" %in% names(dots)) dir_fig <- dots["dir_fig"] else do_figures <- FALSE

	b_migtype <- (b - 1) * length(get("migtypes", envir = etr_vars)) + which(migtype == get("migtypes", envir = etr_vars))
	etmeasure$etable[b_migtype, "Transect_ID"] <- i
	etmeasure$etable[b_migtype, "Neighbor_Cells"] <- neighborhoods(ecotoner_settings)[b]
	etmeasure$etable[b_migtype, "Migration_Type"] <- migtype
	
	if (!copy_FromMig1_TF) {
		etmeasure$gETmeas[[b]][[migtype]]$optim <- calc_Eppinga2013_optpos(Veg1 = etband$Veg[[migtype]]$Veg1$grid,
																		   Veg2 = etband$Veg[[migtype]]$Veg2$grid)
		
		etmeasure$gETmeas[[b]][[migtype]]$adv_veg1 <- calc_Eppinga2013_advancement(veg = etband$Veg[[migtype]]$Veg1$grid,
																		end_toLeft = end_to_left(type_veg1(ecotoner_settings)),
																		optBoundary_cell = etmeasure$gETmeas[[b]][[migtype]]$optim$pos_cell)
		etmeasure$gETmeas[[b]][[migtype]]$adv_veg2 <- calc_Eppinga2013_advancement(veg = etband$Veg[[migtype]]$Veg2$grid,
																		end_toLeft = end_to_left(type_veg2(ecotoner_settings)),
																		optBoundary_cell = etmeasure$gETmeas[[b]][[migtype]]$optim$pos_cell)
		
		etmeasure$gETmeas[[b]][[migtype]]$adv_stats <- calc_Eppinga2013_stats(deltaF1_T17 = etmeasure$gETmeas[[b]][[migtype]]$adv_veg1$deltaFrontRunners_T17,
																		deltaF2_T17 = etmeasure$gETmeas[[b]][[migtype]]$adv_veg2$deltaFrontRunners_T17,
																		seed = seed)

		if(do_figures) map_front_runners_Eppinga2013(filename = file.path(dir_fig, paste0(flag_bfig, "Eppinga2013_FittedBoundaryAdvancement_", migtype, ".pdf")),
													eB_Env = etband$Env, eB_Veg = etband$Veg[[migtype]],
													datFit = etmeasure$gETmeas[[b]][[migtype]])				
	}
	
	etmeasure$etable <- tabulate_Eppinga2013_advance(etable = etmeasure$etable, index = b_migtype,
													data = etmeasure$gETmeas[[b]][[if (copy_FromMig1_TF) "AllMigration" else migtype]])
		
	etmeasure
}
