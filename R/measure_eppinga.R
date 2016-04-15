#------Eppinga, M.B., Pucko, C.A., Baudena, M., Beckage, B. & Molofsky, J. (2013) A new method to infer vegetation boundary movement from 'snapshot' data. Ecography, 36, 622-635.

version_Eppinga2013Ecography <- function() numeric_version("0.2.2")


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
	# front runner distance to optimal boundary in meters
	FR_dist_m <- raster::xres(veg) * (calc_Eppinga2013_frontrunners(veg, end_toLeft) - optBoundary_cell)
	FR_dist_T17 <- transformation17(FR_dist_m)
	FR_dist_mean_T17 <- mean(FR_dist_T17, na.rm = TRUE)
		
	list(FR_dist_T17 = FR_dist_T17, FR_dist_m = FR_dist_m,
		FR_dist_mean_m = backtransformation17(FR_dist_mean_T17),
		FR_dist_mean_T17 = FR_dist_mean_T17,
		FR_dist_sd_T17 = sd(FR_dist_T17, na.rm = TRUE))
}



#---Eppinga et al. 2013: 'Statistical analyses'
calc_Eppinga2013_stats <- function(FR_dist_T17_veg1, FR_dist_mean_T17_veg1, FR_dist_T17_veg2, FR_dist_mean_T17_veg2, seed = NULL) {
	#assumption that all frontrunners are y-row wise dispersed and do not originate from other sources

	res_names <- c("FrontsAdvBeyondOptBoundary", "FD_mean_T17", "FD_sd_T17", "FD_mean_m",
					"FD_boots_R",
					"FD_iidboot_mean", "FD_iidboot_bias", "FD_iidboot_se", "FD_iidboot_ci_type", "FD_iidboot_ci0_p", "FD_iidboot_freq_p",
					"FD_depboot_mean", "FD_depboot_bias", "FD_depboot_se", "FD_depboot_ci_type", "FD_depboot_ci0_p", "FD_depboot_freq_p", "FD_depboot_bstar", 
					"FD_WSRT_Z", "FD_WSRT_p", "FD_WSRT_midp",
					"FD_retro_power")
	res <- as.list(rep(NA, length(res_names)))
	names(res) <- res_names
	boot_ci_types <- c("normal", "basic", "student", "percent", "bca")	#see body(boot::boot.ci)
	
	res[["FrontsAdvBeyondOptBoundary"]] <- all(FR_dist_mean_T17_veg1 > 0, FR_dist_mean_T17_veg2 < 0)

	
	if (res[["FrontsAdvBeyondOptBoundary"]]) {
		if (!is.na(seed)) set.seed(seed)
		
		#---Test if vegetation types have advanced comparably
		#	i.e., test if difference between front runner distances (FD) of veg1 vs veg2 is 0: eq. 18
		data_T17 <- cbind(FR_dist_T17_veg1, -FR_dist_T17_veg2)
		delta_T17 <- apply(data_T17, 1,  function(x) x[1] - x[2])
		res[["FD_mean_T17"]] = mean(delta_T17, na.rm = TRUE)
		res[["FD_sd_T17"]] = sd(delta_T17, na.rm = TRUE)
		res[["FD_mean_m"]] = backtransformation17(res[["FD_mean_T17"]])
		res[["FD_boots_R"]] <- 1e5L # Eppinga et al. 2013: R = 1e5

		bmds <- list(iid = list("boot" = NULL, "ci" = NULL),
					 dep = list("boot" = NULL, "ci" = NULL))
		
		if (requireNamespace("boot", quietly = TRUE)) {
			#- Bootstrap approach assuming independent data (as used by Eppinga et al. 2013)
			bmds[["iid"]][["boot"]] <- boot::boot(data = data_T17,
										  statistic = indexed_mean_of_diffs,
										  R = res[["FD_boots_R"]],
										  sim = "ordinary", stype = "i",
										  parallel = "no")
			
			# adjusted bootstrap percentile (BCa) interval
			bmds[["iid"]][["ci"]] <- boot::boot.ci(bmds[["iid"]][["boot"]],
										conf = c(0.95, 0.99, 0.999), type = "bca")
		} else {
			warning("Package 'boot' not installed: 'calc_Eppinga2013_stats' will not estimate iid bootstrap")
		}
		
		if (requireNamespace("boot", quietly = TRUE) && requireNamespace("np", quietly = TRUE)) {
			# Stationary block bootstrap (Politis, D. N., and J. P. Romano. 1994. The Stationary Bootstrap. Journal of the American Statistical Association 89:1303-1313.)
			#	with optimal mean block length (Patton, A., D. N. Politis, and H. White. 2009. Correction to "Automatic Block-Length Selection for the Dependent Bootstrap" by D. Politis and H. White (vol 23, pg 53, 2004). Econometric Reviews 28:372-375.)

			if (anyNA(delta_T17)) {
				# Multiple imputations of time-series data
				# 'np::b.star' calls acf() and ccf() with default value for 'na.action', i.e., 'na.fail'
				# 'np::b.star' fails if anyNA(data); data contain NAs if a transect row has no cells of a vegetation type
				if (requireNamespace("Amelia", quietly = TRUE)) {
					im_T17 <- cbind(seq_len(nrow(data_T17)), data_T17)
					am_T17 <- Amelia::amelia(im_T17, m = 10, ts = 1, splinetime = 6,
											p2s = 0, parallel = "no")
				
					est_BstarSB <- mean(sapply(am_T17$imputations, function(x)
											np::b.star(apply(x, 1,  function(x) x[2] - x[3]), round = TRUE)[, "BstarSB"]))
				} else {
					warning("Package 'Amelia' not installed: 'calc_Eppinga2013_stats' cannot estimate optimal block-length for dependent bootstrap with missing data; poor approximation based on complete-cases used instead")
					est_BstarSB <- np::b.star(delta_T17[complete.cases(delta_T17)], round = TRUE)[, "BstarSB"]
				}
			} else {
				est_BstarSB <- np::b.star(delta_T17, round = TRUE)[, "BstarSB"]
			}
			
			res[["FD_depboot_bstar"]] <- min(nrow(data_T17), est_BstarSB)
			bmds[["dep"]][["boot"]] <- boot::tsboot(tseries = data_T17,
										 statistic = indexed_mean_of_diffs,
										 R = res[["FD_boots_R"]],
										 sim = "geom", l = res[["FD_depboot_bstar"]], endcorr = TRUE, n.sim = nrow(data_T17),
										 orig.t = TRUE, parallel = "no")
		
			# BCa and studentized CI don't apply for tsboot objects; use instead percentile method
			bmds[["dep"]][["ci"]] <- boot::boot.ci(bmds[["dep"]][["boot"]],
										conf = c(0.95, 0.99, 0.999), type = "perc")
		} else {
			warning("Package 'boot' and/or 'np' not installed: 'calc_Eppinga2013_stats' will not estimate dependent bootstrap")
		}

		# Extract bootstrap data
		ptol <- sqrt(.Machine$double.eps)
		ntol <- -sqrt(.Machine$double.neg.eps)
		for (ib in names(bmds)) if (!is.null(bmds[[ib]][["boot"]])) {
			res[[paste0("FD_", ib, "boot_mean")]] = mean(bmds[[ib]][["boot"]]$t, na.rm = TRUE)
			res[[paste0("FD_", ib, "boot_bias")]] <- res[[paste0("FD_", ib, "boot_mean")]] - res[["FD_mean_T17"]]
			res[[paste0("FD_", ib, "boot_se")]] <- as.numeric(sqrt(var(bmds[[ib]][["boot"]]$t, na.rm = TRUE)))

			# Test approach 1a: Is 0 contained in ci?
			i_item <- 4 # the i-th item in the list that is returned by boot::boot.ci
			res[[paste0("FD_", ib, "boot_ci_type")]] <- which(names(bmds[[ib]][["ci"]])[i_item] == boot_ci_types)
			conf <- bmds[[ib]][["ci"]][[i_item]]
			pid <- !(as.integer(apply(conf[, 4:5], 1, function(x) sum(x > ptol))) == 1)
			res[[paste0("FD_", ib, "boot_ci0_p")]] <- 1 - if (sum(pid) > 0) max(conf[pid, "conf"]) else 0 # this represents steps of 1, 0.05, 0.01, and 0.001 as upper bound of the p-value
			
			# Test approach 1b: Calculate p-value (for H0: diff = 0); eq. 19
			# Direct interpretation of eq. 19: sum(Heaviside(abs(bmdiid$t) + abs(bmdiid$t0) - abs(bmdiid$t + bmdiid$t0))) / bmdiid$R
			res[[paste0("FD_", ib, "boot_freq_p")]] <- (if (bmds[[ib]][["boot"]]$t0 > ptol) sum(bmds[[ib]][["boot"]]$t <= ptol) else if (bmds[[ib]][["boot"]]$t0 < ntol) sum(bmds[[ib]][["boot"]]$t >= ntol) else bmds[[ib]][["boot"]]$R) / bmds[[ib]][["boot"]]$R
		}

			
		if (requireNamespace("coin", quietly = TRUE)) {
			# Test approach 2: exact Wilcoxon signed rank test (with Pratt correction of zeros)
			wsrt <- coin::wilcoxsign_test(FR_dist_T17_veg1 ~ FR_dist_T17_veg2, distribution = "exact", alternative = "two.sided")
			res[["FD_WSRT_Z"]] <- coin::statistic(wsrt, type = "test")
			res[["FD_WSRT_p"]] <- coin::pvalue(wsrt)
			res[["FD_WSRT_midp"]] <- coin::midpvalue(wsrt)
		} else {
			warning("Package 'coin' not installed: 'calc_Eppinga2013_stats' will not calculate Wilcoxon signed rank test")
		}
			
		#---Retrospective power: Eppinga et al. 2013: eq. 20
		tau <- 0.2 #effect size
		n <- sum(complete.cases(data_T17))
		tcrit <- qt(0.95, df = n) #95% confidence
		res[["FD_retro_power"]] <- 1 - (1/2 * (1 + erf((tcrit - tau * sqrt(n) / res[["FD_sd_T17"]]) / (sqrt(2) * res[["FD_sd_T17"]]))))
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
	isnotna <- !is.na(datFit$adv_veg1$FR_dist_m)
	points(x = (copt + datFit$adv_veg1$FR_dist_m)[isnotna], y = ys[isnotna], pch = 46, cex = 2, col = "magenta")
	isnotna <- !is.na(datFit$adv_veg2$FR_dist_m)
	points(x = (copt + datFit$adv_veg2$FR_dist_m)[isnotna], y = ys[isnotna], pch = 46, cex = 2, col = "green")

	if (datFit$adv_stats$FrontsAdvBeyondOptBoundary) {
		if (!is.na(datFit$adv_stats$FD_depboot_freq_p) && !is.null(datFit$adv_stats$FD_depboot_freq_p)) {
			p_max2 <- max(datFit$adv_stats$FD_depboot_freq_p, na.rm = TRUE)
			p_max <- max(c(datFit$adv_stats$FD_depboot_ci0_p, p_max2), na.rm = TRUE) # 'FD_XXXboot_ci0_p' only comes in steps of 1, 0.05, 0.01, and 0.001 as upper bound of the p-value
			ttype <- "Stationary bootstrap p"
		} else if (!is.na(datFit$adv_stats$FD_iidboot_freq_p) && !is.null(datFit$adv_stats$FD_iidboot_freq_p)) {
			p_max2 <- max(datFit$adv_stats$FD_iidboot_freq_p, na.rm = TRUE)
			p_max <- max(c(datFit$adv_stats$FD_iidboot_ci0_p, p_max2), na.rm = TRUE) # 'FD_XXXboot_ci0_p' only comes in steps of 1, 0.05, 0.01, and 0.001 as upper bound of the p-value
			ttype <- "Ordinary iid bootstrap p"
		} else if (!is.na(datFit$adv_stats$FD_WSRT_p) && !is.null(datFit$adv_stats$FD_WSRT_p)) {
			p_max <- p_max2 <- max(datFit$adv_stats$FD_WSRT_p, na.rm = TRUE)
			ttype <- "Wilcoxon signed rank p"
		} else {
			p_max <- 1
			ttype <- "No test p"
		}
		p_value <- if (abs(p_max - 1) < sqrt(.Machine$double.eps)) p_max2 else p_max
		padv <- p_max < 0.05
		p12 <- padv && datFit$adv_stats$FD_mean_T17 > 0
		p21 <- padv && datFit$adv_stats$FD_mean_T17 < 0
		p_str <- ifelse(p_max > 0.001, formatC(p_value, format = "f", digits = 3), "0.001")
		
		if (padv && abs(datFit$adv_stats$FD_mean_m) > 2 * res_m) {
			# Mean front runner advancement
			pos1 <- copt + datFit$adv_veg1$FR_dist_mean_m
			segments(x0 = pos1, y0 = ext1@ymin, x1 = pos1, y1 = ext1@ymax, lwd = 2, col = adjustcolor("magenta", alpha.f = 0.7))
			pos2 <- copt + datFit$adv_veg2$FR_dist_mean_m
			segments(x0 = pos2, y0 = ext1@ymin, x1 = pos2, y1 = ext1@ymax, lwd = 2, col = adjustcolor("green", alpha.f = 0.7))
			arrows(x0 = copt, y0 = ext1@ymax + res_m,
				   x1 = if (p12) pos1 else if (p21) pos2 else 0, y1 = ext1@ymax + res_m,
				   lwd = 2, col = if (p12) "magenta" else if (p21) "green" else "white")
		}

		ptext <- paste0("adv(Veg1 = ", signif(datFit$adv_veg1$FR_dist_mean_m, 2), " m) ",
						if (p12) ">" else if (p21) "<" else "=",
						" adv(Veg2 = ", signif(-datFit$adv_veg2$FR_dist_mean_m, 2), " m):",
						"\n", ttype, if (p_max < 0.05) " < " else " >= ", p_str,
						" with a retrospective power of ", signif(datFit$adv_stats$FD_retro_power, 3))


	} else {
		ptext <- paste0("no advancement beyond optimal boundary location")
	}
	text(x = ext1@xmax / 2, y = ext1@ymax + 2 * strheight("0", units = "user", cex = 2 / 3 * cex), labels = ptext, cex = 2/3)
	
	invisible()
}


tabulate_Eppinga2013_advance <- function(etable, index, data) {
	temp1 <- unlist(data$adv_veg1[c("FR_dist_mean_m", "FR_dist_mean_T17", "FR_dist_sd_T17")])
	names(temp1) <- paste0("Veg1_", names(temp1))
	temp2 <- unlist(data$adv_veg2[c("FR_dist_mean_m", "FR_dist_mean_T17", "FR_dist_sd_T17")])
	names(temp2) <- paste0("Veg2_", names(temp2))
	
	tabulate_merge_into_etable(etable, index,
		data = c(unlist(data$optim), temp1, temp2, unlist(data$adv_stats)))
}




#' @export
Eppinga2013Ecography <- function(i, b, migtype, ecotoner_settings, etband, etmeasure, copy_FromMig1_TF, do_figures, ...) {
	#3b. Eppinga et al. 2013 Ecography: Location of boundary and front-runner distance
	#Objective: inference of vegetation boundary movement from one 'snapshot' (e.g. an aerial photograph or satellite image) in time
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
	
	etmeasure$gETmeas[[b]][[migtype]]$meta <- list()
	etmeasure$gETmeas[[b]][[migtype]]$meta$method <- "Eppinga2013Ecography"
	etmeasure$gETmeas[[b]][[migtype]]$meta$version <- version_Eppinga2013Ecography()
	etmeasure$gETmeas[[b]][[migtype]]$meta$copy_FromMig1_TF <- copy_FromMig1_TF

	
	if (!copy_FromMig1_TF) {
		etmeasure$gETmeas[[b]][[migtype]]$optim <- calc_Eppinga2013_optpos(Veg1 = etband$Veg[[migtype]]$Veg1$grid,
																		   Veg2 = etband$Veg[[migtype]]$Veg2$grid)
		
		etmeasure$gETmeas[[b]][[migtype]]$adv_veg1 <- calc_Eppinga2013_advancement(veg = etband$Veg[[migtype]]$Veg1$grid,
																		end_toLeft = end_to_left(type_veg1(ecotoner_settings)),
																		optBoundary_cell = etmeasure$gETmeas[[b]][[migtype]]$optim$pos_cell)
		etmeasure$gETmeas[[b]][[migtype]]$adv_veg2 <- calc_Eppinga2013_advancement(veg = etband$Veg[[migtype]]$Veg2$grid,
																		end_toLeft = end_to_left(type_veg2(ecotoner_settings)),
																		optBoundary_cell = etmeasure$gETmeas[[b]][[migtype]]$optim$pos_cell)
		
		etmeasure$gETmeas[[b]][[migtype]]$adv_stats <- calc_Eppinga2013_stats(FR_dist_T17_veg1 = etmeasure$gETmeas[[b]][[migtype]]$adv_veg1$FR_dist_T17,
																		FR_dist_mean_T17_veg1 = etmeasure$gETmeas[[b]][[migtype]]$adv_veg1$FR_dist_mean_T17,
																		FR_dist_T17_veg2 = etmeasure$gETmeas[[b]][[migtype]]$adv_veg2$FR_dist_T17,
																		FR_dist_mean_T17_veg2 = etmeasure$gETmeas[[b]][[migtype]]$adv_veg2$FR_dist_mean_T17,
																		seed = seed)

		if(do_figures) map_front_runners_Eppinga2013(filename = file.path(dir_fig, paste0(flag_bfig, "Eppinga2013_FittedBoundaryAdvancement_", migtype, ".pdf")),
													eB_Env = etband$Env, eB_Veg = etband$Veg[[migtype]],
													datFit = etmeasure$gETmeas[[b]][[migtype]])				
	}
	
	etmeasure$etable <- tabulate_Eppinga2013_advance(etable = etmeasure$etable, index = b_migtype,
													data = etmeasure$gETmeas[[b]][[if (copy_FromMig1_TF) "AllMigration" else migtype]])
		
	etmeasure
}
