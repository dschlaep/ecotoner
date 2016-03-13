
#------Danz, N.P., Frelich, L.E., Reich, P.B. & Niemi, G.J. (2012) Do vegetation boundaries display smooth or abrupt spatial transitions along environmental gradients? Evidence from the prairie–forest biome boundary of historic Minnesota, USA. Journal of Vegetation Science, n/a-n/a.



#---Danz et al. 2012: 'Spatial change in vegetation across the boundary', 'Spatial change in climate across the boundary', and 'Spatial vegetation/climate relationships across the boundary'
#' Estimates the abruptness of outcome y against x
#' 
#' @param x A numeric vector. Predictor along the ecotone transect from left-right. The length of x has to correspond to the ecotone transect. For instance, distance along the transect.
#' @param z A numeric vector. The outcome along the ecotone transect, e.g., vegetation density.
#' 
calc_Danz2012_abruptness_1D <- function(x, z, seed = NULL) {
	stopifnot(length(x) == length(z))
	
	if (!is.na(seed)) set.seed(seed)

	#Variation from published method: always use scaled and centered data
	xmin <- min(x, na.rm = TRUE)
	xmax <- max(x, na.rm = TRUE)
	x <- c(scale(x, center = xmin, scale = xmax - xmin))
	ymin <- min(z, na.rm = TRUE)
	ymax <- max(z, na.rm = TRUE)
	z <- c(scale(z, center = ymin, scale = ymax - ymin))
	
	#choice between linear vs sigmoidal only reflects how transect is defined and not necessarily how abrupt it is, i.e., whether 'shoulders' of vegetation density at low and high values are included in transect or not
	ms <- try(nls(z ~ sigmoidal(x, b, c), start = list(b = 1, c = 0.5)), silent = TRUE)
	reps <- 0
	while (inherits(ms, "try-error") && reps <= 50) {
		ms <- try(nls(z ~ sigmoidal(x, b, c), start = list(b = runif(1), c = runif(1))), silent = TRUE)
		reps <- reps + 1
	}
	
	ml <- try(lm(z ~ x), silent = TRUE)
	
	mbest <- "none"
	deltaAIC <- abb_ms <- abb_ml <- abbScaled_ms <- abbScaled_ml <- width_ms <- x_5percOfy_ms <- x_95percOfy_ms <- width_ml <- x_5percOfy_ml <- x_95percOfy_ml <- cenx_ms <- cenx_ml <- r2ms <- r2ml <- AICms <- AICml <- NA
	ceny <- ymin + 0.5 * (ymax - ymin)
	if (!inherits(ml, "try-error")) {
		abb_ml <- coef(ml)[2]
		abbScaled_ml <- abb_ml / (xmax - xmin)
		cenx_ml <- xmin + (0.5 - coef(ml)[1]) / coef(ml)[2] * (xmax - xmin)
		x_95percOfy_ml <- xmin + (0.95 - coef(ml)[1]) / coef(ml)[2] * (xmax - xmin)
		x_5percOfy_ml <- xmin + (0.05 - coef(ml)[1]) / coef(ml)[2] * (xmax - xmin)
		width_ml <- abs(x_95percOfy_ml - x_5percOfy_ml)
		AICml <- AIC(ml)
		r2ml <- summary(ml)$r.squared
	}
	if (!inherits(ms, "try-error")) {
		abb_ms <- ms$m$getPars()[1]
		abbScaled_ms <- abb_ms / (xmax - xmin)
		cenx_ms <- xmin + ms$m$getPars()[2] * (xmax - xmin)
		x_95percOfy_ms <- xmin + sigmoidal_inv(0.95, ms$m$getPars()[1], ms$m$getPars()[2]) * (xmax - xmin)
		x_5percOfy_ms <- xmin + sigmoidal_inv(0.05, ms$m$getPars()[1], ms$m$getPars()[2]) * (xmax - xmin)
		width_ms <- abs(x_95percOfy_ms - x_5percOfy_ms)
		AICms <- AIC(ms)
		r2ms <- 1 - sum(residuals(ms)^2) / sum((mean(z) - z)^2)	#residuals(ms) ==  predict(ms) - z
		if (!inherits(ml, "try-error")) {
			deltaAIC <- AICms - AICml
			mbest <- ifelse(deltaAIC < 0, "sigmoidal", "linear")
		}
	}
	lms <- list(m = ms, r2 = r2ms, aic = AICms, center_x = cenx_ms, abruptness_b = abb_ms, abruptness_bscaled = abbScaled_ms, width_m = width_ms, x_5percOfy = x_5percOfy_ms, x_95percOfy = x_95percOfy_ms)
	lml <- list(m = ml, r2 = r2ml, aic = AICml, center_x = cenx_ml, abruptness_b = abb_ml, abruptness_bscaled = abbScaled_ml, width_m = width_ml, x_5percOfy = x_5percOfy_ml, x_95percOfy = x_95percOfy_ml)
	lb <- if (!is.na(deltaAIC)) {
				if (deltaAIC < 0) lms else lml
			} else {
				list(m = NA, r2 = NA, aic = NA, center_x = NA, abruptness_b = NA, abruptness_bscaled = NA, width_m = NA, x_5percOfy = NA, x_95percOfy = NA)
			}
	
	list(sigmoidal = lms, linear = lml, best = lb, mbest_byAIC = mbest, deltaAIC = deltaAIC, center_y = ceny, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
}


tabulate_Danz2012_abruptness_1D <- function(etable, b, data, flag_migtype) {

	place_set1 <- function (etable, b, y, yname, vsname, flag_migtype) {
		colnamesAdd <- paste0(flag_migtype, "_Danz2012_", yname, "_", vsname, "_", 
						c("mbest_byAIC", "deltaAIC", "y_at50percOfy"))

		res <- data.frame(matrix(NA, nrow = max(1, nrow(etable)), ncol = length(colnamesAdd), dimnames = list(NULL, colnamesAdd)))
		res[b, 1] <- y["mbest_byAIC"]
		res[b, 2:3] <- unlist(y[c("deltaAIC", "center_y")])

		cbind(etable, res)
	}
	
	place_set2 <- function (etable, b, y, yname, vsname, mname, flag_migtype) {
		colnamesAdd <- paste0(flag_migtype, "_Danz2012_", yname, "_", vsname, "_", mname, "_",
						c("x_at05percOfy", "x_at50percOfy", "x_at95percOfy", "Width_m", "Abruptness_b", "Abruptness_bscaled", "AIC", "r2"))

		res <- matrix(NA, nrow = max(1, nrow(etable)), ncol = length(colnamesAdd), dimnames = list(NULL, colnamesAdd))
		res[b, seq_along(colnamesAdd)] <- unlist(y[c("x_5percOfy", "center_x", "x_95percOfy", "width_m", "abruptness_b", "abruptness_bscaled", "aic", "r2")])

		cbind(etable, res)
	}
	
	sets <- list(c("ElevVsDist_1D", "Elev", "vsDist"),
				c("Veg1VsDist_1D", "Veg1", "vsDist"), c("Veg1VsElev_1D", "Veg1", "vsElev"),
				c("Veg2VsDist_1D", "Veg2", "vsDist"), c("Veg2VsElev_1D", "Veg2", "vsElev"))
	mnames <- c("sigmoidal", "linear", "best")

	for (set in sets) {
		etable <- place_set1(etable, b, y = data[[set[1]]], yname = set[2], vsname = set[3], flag_migtype)
		for (mname in mnames) {
			etable <- place_set2(etable, b, y = data[[set[1]]][[mname]], yname = set[2], vsname = set[3], mname, flag_migtype)
		}
	}
	
	etable
}


#' @export
plot_Danz2012_abruptness_1D <- function(filename, eB_Env, eB_Veg, datFit) {

	pdf(width = 7, height = 8, file = filename)
	par_old <- par(mfrow = c(nr <- 2, 1), mar = c(2.5, 2.5, 0.5, 2.5), mgp = c(1.5, 0.5, 0), cex = cex <- 1.5)
	on.exit({par(par_old); dev.off()}, add = TRUE)


	#Panel (a): Danz et al. 2012 Fig. 3 and 4
	emin <- min(eB_Env$elev$YMeans_ForEachX)
	emax <- max(eB_Env$elev$YMeans_ForEachX)
	
	plot(eB_Env$DistAlongXaxis_m, eB_Env$elev$YMeans_ForEachX, type = "l", lwd = 2, xlab = "Transect length (m)", ylab = "Elevation (m)", axes = FALSE)
	at1 <- c(axTicks(1), max(eB_Env$DistAlongXaxis_m))
	if((temp <- tail(at1, n = 2))[1] > temp[2]) at1 <- c(at1[1:(length(at1)-2)], at1[length(at1)])
	axis(1, pos = emin, at = at1)
	axis(2, pos = 0, at = round(c(emin, axTicks(2), emax)))

	lines(eB_Env$DistAlongXaxis_m, transform_01_to_R(eB_Veg$Veg1$density, emin, emax), col = "red")
	lines(eB_Env$DistAlongXaxis_m, transform_01_to_R(eB_Veg$Veg2$density, emin, emax), col = "darkgreen")
	axis(4, pos = max(at1), at = at <- seq(emin, emax, length = 11), labels = round(transform_R_to_01(at, emin, emax), 1))
	mtext(text = "Vegetation density", side = 4, cex = cex, line = 1)
	mtext(text = "(a)", line = -1, cex = cex, adj = 0.05)
	legend("top", bty = "n", cex = 1/3*cex, legend = c("Elevation", "Big sagebrush ecosystem", "Temperate Forest", "Sigmoidal fit, center, & 90% range"), col = c("black", "red", "darkgreen", "black"), lty = c(1, 1, 1, 2), lwd = c(2, 1, 1, 2))

	#Veg vs Dist: BSE
	if (!inherits(datFit$Veg1VsDist_1D$sigmoidal$m, "try-error")) {
		ypred <- datFit$Veg1VsDist_1D$ymin + predict(datFit$Veg1VsDist_1D$sigmoidal$m) * (datFit$Veg1VsDist_1D$ymax - datFit$Veg1VsDist_1D$ymin)
		lines(eB_Env$DistAlongXaxis_m, transform_01_to_R(ypred, emin, emax), lty = 2, lwd = 2, col = "red")
		line_y <- par()$usr[3] + 0.02*(par()$usr[4]-par()$usr[3])
		arr_x0 <- datFit$Veg1VsDist_1D$sigmoidal$center_x
		arr_y0 <- transform_01_to_R(datFit$Veg1VsDist_1D$ymin + 0.5 * (datFit$Veg1VsDist_1D$ymax - datFit$Veg1VsDist_1D$ymin), emin, emax)
		arrows(x0 = arr_x0, y0 = arr_y0, x1 = arr_x0, y1 = line_y, lwd = 2, col = "red", length = 0.03*par()$pin[1]/nr)
		seg_x0 <- datFit$Veg1VsDist_1D$sigmoidal$x_5percOfy
		seg_y0 <- transform_01_to_R(datFit$Veg1VsDist_1D$ymin + 0.05 * (datFit$Veg1VsDist_1D$ymax - datFit$Veg1VsDist_1D$ymin), emin, emax)
		seg_x1 <- datFit$Veg1VsDist_1D$sigmoidal$x_95percOfy
		seg_y1 <- transform_01_to_R(datFit$Veg1VsDist_1D$ymin + 0.95 * (datFit$Veg1VsDist_1D$ymax - datFit$Veg1VsDist_1D$ymin), emin, emax)
		segments(x0 = seg_x0, y0 = line_y, x1 = seg_x1, y1 = line_y, lwd = 2, col = "red")
		segments(x0 = seg_x0, y0 = seg_y0, x1 = seg_x0, y1 = line_y, lwd = 2, col = "red")
		segments(x0 = seg_x1, y0 = seg_y1, x1 = seg_x1, y1 = line_y, lwd = 2, col = "red")
	}
	#Veg vs Dist: TF
	if (!inherits(datFit$Veg2VsDist_1D$sigmoidal$m, "try-error")) {
		ypred <- datFit$Veg2VsDist_1D$ymin + predict(datFit$Veg2VsDist_1D$sigmoidal$m) * (datFit$Veg2VsDist_1D$ymax - datFit$Veg2VsDist_1D$ymin)
		lines(eB_Env$DistAlongXaxis_m, transform_01_to_R(ypred, emin, emax), lty = 2, lwd = 2, col = "darkgreen")
		line_y <- par()$usr[3] + 0.03*(par()$usr[4]-par()$usr[3])
		arr_x0 <- datFit$Veg2VsDist_1D$sigmoidal$center_x
		arr_y0 <- transform_01_to_R(datFit$Veg2VsDist_1D$ymin + 0.5 * (datFit$Veg2VsDist_1D$ymax - datFit$Veg2VsDist_1D$ymin), emin, emax)
		arrows(x0 = arr_x0, y0 = arr_y0, x1 = arr_x0, y1 = line_y, lwd = 2, col = "darkgreen", length = 0.03*par()$pin[1]/nr)
		seg_x0 <- datFit$Veg2VsDist_1D$sigmoidal$x_5percOfy
		seg_y0 <- transform_01_to_R(datFit$Veg2VsDist_1D$ymin + 0.05 * (datFit$Veg2VsDist_1D$ymax - datFit$Veg2VsDist_1D$ymin), emin, emax)
		seg_x1 <- datFit$Veg2VsDist_1D$sigmoidal$x_95percOfy
		seg_y1 <- transform_01_to_R(datFit$Veg2VsDist_1D$ymin + 0.95 * (datFit$Veg2VsDist_1D$ymax - datFit$Veg2VsDist_1D$ymin), emin, emax)
		segments(x0 = seg_x0, y0 = line_y, x1 = seg_x1, y1 = line_y, lwd = 2, col = "darkgreen")
		segments(x0 = seg_x0, y0 = seg_y0, x1 = seg_x0, y1 = line_y, lwd = 2, col = "darkgreen")
		segments(x0 = seg_x1, y0 = seg_y1, x1 = seg_x1, y1 = line_y, lwd = 2, col = "darkgreen")
	}

	#Panel (b): Danz et al. 2012 Fig. 3 inset and Fig. 6
	#Veg vs Elev: BSE
	plot(sort(eB_Env$elev$YMeans_ForEachX), eB_Veg$Veg1$density[order(eB_Env$elev$YMeans_ForEachX)], xlim = c(min(eB_Env$elev$YMeans_ForEachX), max(eB_Env$elev$YMeans_ForEachX)), ylim = c(0, 1), type = "l", col = "red", xlab = "Elevation (m)", ylab = "Vegetation density", axes = FALSE)
	at2 <- c(emin, axTicks(1), emax)
	if ((temp <- tail(at2, n = 2))[1] > temp[2]) at2 <- c(at2[1:(length(at2)-2)], at2[length(at2)])
	if ((temp <- head(at2, n = 2))[1] > temp[2]) at2 <- c(at2[1], at2[3:length(at2)])
	axis(1, pos = 0, at = round(at2))
	axis(2, pos = emin)

	if (!inherits(datFit$Veg1VsElev_1D$sigmoidal$m, "try-error")) {
		ypred <- datFit$Veg1VsElev_1D$ymin + predict(datFit$Veg1VsElev_1D$sigmoidal$m) * (datFit$Veg1VsElev_1D$ymax - datFit$Veg1VsElev_1D$ymin)
		lines(eB_Env$elev$YMeans_ForEachX, ypred, col = "red", lty = 2, lwd = 2)
		line_y <- par()$usr[3] + 0.02*(par()$usr[4]-par()$usr[3])
		arr_x0 <- datFit$Veg1VsElev_1D$sigmoidal$center_x
		arr_y0 <- datFit$Veg1VsElev_1D$ymin + 0.5 * (datFit$Veg1VsElev_1D$ymax - datFit$Veg1VsElev_1D$ymin)
		arrows(x0 = arr_x0, y0 = arr_y0, x1 = arr_x0, y1 = line_y, lwd = 2, col = "red", length = 0.03*par()$pin[1]/nr)
		seg_x0 <- datFit$Veg1VsElev_1D$sigmoidal$x_5percOfy
		seg_y0 <- datFit$Veg1VsElev_1D$ymin + 0.05 * (datFit$Veg1VsElev_1D$ymax - datFit$Veg1VsElev_1D$ymin)
		seg_x1 <- datFit$Veg1VsElev_1D$sigmoidal$x_95percOfy
		seg_y1 <- datFit$Veg1VsElev_1D$ymin + 0.95 * (datFit$Veg1VsElev_1D$ymax - datFit$Veg1VsElev_1D$ymin)
		segments(x0 = seg_x0, y0 = line_y, x1 = seg_x1, y1 = line_y, lwd = 2, col = "red")
		segments(x0 = seg_x0, y0 = seg_y0, x1 = seg_x0, y1 = line_y, lwd = 2, col = "red")
		segments(x0 = seg_x1, y0 = seg_y1, x1 = seg_x1, y1 = line_y, lwd = 2, col = "red")
	}
	
	#Veg vs Elev: TF
	lines(sort(eB_Env$elev$YMeans_ForEachX), eB_Veg$Veg2$density[order(eB_Env$elev$YMeans_ForEachX)], col = "darkgreen")
	
	if (!inherits(datFit$Veg2VsElev_1D$sigmoidal$m, "try-error")) {
		ypred <- datFit$Veg2VsElev_1D$ymin + predict(datFit$Veg2VsElev_1D$sigmoidal$m) * (datFit$Veg2VsElev_1D$ymax - datFit$Veg2VsElev_1D$ymin)
		lines(eB_Env$elev$YMeans_ForEachX, ypred, col = "darkgreen", lty = 2, lwd = 2)
		line_y <- par()$usr[3] + 0.03*(par()$usr[4]-par()$usr[3])
		arr_x0 <- datFit$Veg2VsElev_1D$sigmoidal$center_x
		arr_y0 <- datFit$Veg2VsElev_1D$ymin + 0.5 * (datFit$Veg2VsElev_1D$ymax - datFit$Veg2VsElev_1D$ymin)
		arrows(x0 = arr_x0, y0 = arr_y0, x1 = arr_x0, y1 = line_y, lwd = 2, col = "darkgreen", length = 0.03*par()$pin[1]/nr)
		seg_x0 <- datFit$Veg2VsElev_1D$sigmoidal$x_5percOfy
		seg_y0 <- datFit$Veg2VsElev_1D$ymin + 0.05 * (datFit$Veg2VsElev_1D$ymax - datFit$Veg2VsElev_1D$ymin)
		seg_x1 <- datFit$Veg2VsElev_1D$sigmoidal$x_95percOfy
		seg_y1 <- datFit$Veg2VsElev_1D$ymin + 0.95 * (datFit$Veg2VsElev_1D$ymax - datFit$Veg2VsElev_1D$ymin)
		segments(x0 = seg_x0, y0 = line_y, x1 = seg_x1, y1 = line_y, lwd = 2, col = "darkgreen")
		segments(x0 = seg_x0, y0 = seg_y0, x1 = seg_x0, y1 = line_y, lwd = 2, col = "darkgreen")
		segments(x0 = seg_x1, y0 = seg_y1, x1 = seg_x1, y1 = line_y, lwd = 2, col = "darkgreen")
	}

	mtext(text = "(b)", line = -0.5, cex = cex, adj = 0.05)
		
	invisible()
}



#' @export
Danz2012JVegSci_1D <- function(i, b, migtype, ecotoner_settings, etband, etmeasure, copy_FromMig1_TF, do_figures, ...) {
	#---3. Ecological boundaries
	#3a. Danz et al. 2012 J.Veg.Sci.: Shape of vegetation boundary in relation to environmental conditions
	#Objective: of our boundary analysis was to evaluate whether the transition from prairie to forest across the boundary resulted from a smooth or abrupt climatic gradient, i.e. whether the transition followed pattern ‘(a)’ or pattern ‘(b)’ in Fig. 3. We used three analytical tactics to address this objective:
	#	(1) description of the spatial pattern of vegetation transition across the boundary,
	#	(2) evaluation of whether the climate gradient P – PET followed a steeper or shallower pattern of change across the boundary, and
	#	(3) direct modelling of the vegetation–climate relationship across the boundary.
	#Here, use elevation as surrogate driving environmental variable

	dots <- list(...)
	seed <- if ("seed" %in% names(dots)) dots["seed"] else NULL
	if ("flag_bfig" %in% names(dots)) flag_bfig <- dots["flag_bfig"] else do_figures <- FALSE
	if ("dir_fig" %in% names(dots)) dir_fig <- dots["dir_fig"] else do_figures <- FALSE

	etmeasure$etable[b, "Transect_ID"] <- i
	etmeasure$etable[b, "Neighbor_Cells"] <- neighborhoods(ecotoner_settings)[b]
	etmeasure$etable[b, "Migration_Type"] <- migtype
	
	if (!copy_FromMig1_TF) {
		etmeasure$gETmeas[[b]][[migtype]]$ElevVsDist_1D <- calc_Danz2012_abruptness_1D(x = etband$Env$DistAlongXaxis_m,
															z = etband$Env$elev$YMeans_ForEachX, seed = seed)
		etmeasure$gETmeas[[b]][[migtype]]$Veg1VsDist_1D <- calc_Danz2012_abruptness_1D(x = etband$Env$DistAlongXaxis_m,
															z = etband$Veg[[migtype]]$Veg1$density, seed = seed)
		etmeasure$gETmeas[[b]][[migtype]]$Veg2VsDist_1D <- calc_Danz2012_abruptness_1D(x = etband$Env$DistAlongXaxis_m,
															z = etband$Veg[[migtype]]$Veg2$density, seed = seed)
		etmeasure$gETmeas[[b]][[migtype]]$Veg1VsElev_1D <- calc_Danz2012_abruptness_1D(x = etband$Env$elev$YMeans_ForEachX,
															z = etband$Veg[[migtype]]$Veg1$density, seed = seed)
		etmeasure$gETmeas[[b]][[migtype]]$Veg2VsElev_1D <- calc_Danz2012_abruptness_1D(x = etband$Env$elev$YMeans_ForEachX,
															z = etband$Veg[[migtype]]$Veg2$density, seed = seed)

		if (do_figures) plot_Danz2012_abruptness_1D(filename = file.path(dir_fig, paste0(flag_bfig, "Danz2012_FittedBoundaryShapes_1D_", migtype, ".pdf")),
												eB_Env = etband$Env,
												eB_Veg = etband$Veg[[migtype]],
												datFit = etmeasure$gETmeas[[b]][[migtype]])				
	} 
	
	etmeasure$etable <- tabulate_Danz2012_abruptness_1D(etable = etmeasure$etable, b = b,
														data = etmeasure$gETmeas[[b]][[if (copy_FromMig1_TF) "AllMigration" else migtype]],
														flag_migtype = migtype)

	etmeasure	
}

