
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

calc_Danz2012_abruptness_2D <- function(x1, z1, x2, z2, seed = NULL) {
	stopifnot(raster::compareRaster(x2, z2), length(x1) == length(z1))
raster::origin(z2) <- raster::origin(x2)
	
	if (!is.na(seed)) set.seed(seed)
	
	# logistic regression on z ~ x with repeated 'rows' along the transect
	# idea by ld; implementation by drs
	x2p <- raster::rasterToPoints(x2)
	if (cellStats(z2, "countNA")) z2 <- raster::calc(z2, function(z) ifelse(is.na(z), 0, z))
	if (cellStats(x2, "countNA")) z2 <- raster::mask(z2, x2)
	z2p <- raster::rasterToPoints(z2)
stopifnot(x2p[, "y"] == z2p[, "y"])
	
	veg2 <- z2p[, "layer"]
	grad2 <- x2p[, "layer"]
	
	# model fits
	dats <- list('1d' = list(x = x1, y = z1), '2d' = list(x = x2p[, "layer"], y = z2p[, "layer"]))
	
	fits <- lapply(dats, function(dat) {
				linr <- try(lm(y ~ x, data = dat), silent = TRUE)
				# Zuur et al. 2009: binomial GLM: "We assume that Yi is binomial distributed with probability pi and ni = 1 independent trials"
				# Zuur et al. 2009: "The logit and probit link functions assume that you have approximately an equal number of zeros and ones"
				logr <- try(glm(y ~ x, data = dat, family = quasibinomial), silent = TRUE)
				# Zuur et al. 2009: "The clog–log may be an option if you have consid- erably more zeros than ones, or vice versa; the sigmoidal curve is asymmetrical."
				logr_cll <- try(glm(y ~ x, data = dat, family = quasibinomial(link = "cloglog")), silent = TRUE)
				sigm <- {i <- 1
							repeat {
								fit <- try(nls(y ~ sigmoidal(x, b, c), data = dat, start = list(b = runif(1, -0.01, 0.01), c = runif(1, -0.01, 0.01))), silent = TRUE)
								if (i > 50 || !inherits(fit, "try-error")) break
								i <- i + 1
							}
							fit}
			))
	
	lm_fit_1d <- lm(z1 ~ x1)
	logr_fit_2d <- glm(veg2 ~ grad2, family = quasibinomial)
	logr_cll_fit_2d <- glm(veg2 ~ grad2, family = quasibinomial(link = "cloglog"))
	logr_fit_1d <- glm(z1 ~ x1, family = quasibinomial)

	sig_fit_1d <- try(nls(z1 ~ sigmoidal(x1, b, c), start = list(b = 0.01, c = 0.005)), silent = TRUE)
	reps <- 0
	while (inherits(sig_fit_1d, "try-error") && reps <= 50) {
		sig_fit_1d <- try(nls(z1 ~ sigmoidal(x1, b, c), start = list(b = runif(1), c = runif(1))), silent = TRUE)
		reps <- reps + 1
	}
	sig_fit_2d <- try(nls(veg2 ~ sigmoidal(grad2, b, c), start = list(b = 0.01, c = 0.005)), silent = TRUE)
	reps <- 0
	while (inherits(sig_fit_2d, "try-error") && reps <= 50) {
		sig_fit_2d <- try(nls(veg2 ~ sigmoidal(grad2, b, c), start = list(b = runif(1), c = runif(1))), silent = TRUE)
		reps <- reps + 1
	}

xt <- seq(min(x2p[, "layer"]), max(x2p[, "layer"]), length.out = 100)
logr_pred_2d <- predict(logr_fit_2d, newdata = data.frame(grad2 = xt), type = "response", se.fit = TRUE)
#logr_cll_pred_2d <- predict(logr_cll_fit_2d, newdata = data.frame(grad2 = xt), type = "response", se.fit = TRUE)
logr_pred_1d <- predict(logr_fit_1d, newdata = data.frame(x1 = xt), type = "response", se.fit = TRUE)
sig_pred_1d <- predict(sig_fit_1d, newdata = data.frame(x1 = xt), type = "response")
sig_pred_2d <- predict(sig_fit_2d, newdata = data.frame(grad2 = xt), type = "response")

plot(xt, logr_pred_2d$fit, ylim = c(0, 1), type = "l")
polygon(x = c(xt, rev(xt)), y = with(logr_pred_2d, c(fit - 2 * se.fit, rev(fit + 2 * se.fit))), border = NA, col = adjustcolor("black", alpha.f = 0.3))

polygon(x = c(xt, rev(xt)), y = with(logr_pred_1d, c(fit - 2 * se.fit, rev(fit + 2 * se.fit))), border = NA, col = adjustcolor("gray", alpha.f = 0.3))
lines(xt, logr_pred_1d$fit, lty = 2, col = "black")

lines(xt, sig_pred_2d, lty = 1, col = "blue")
lines(xt, sig_pred_1d, lty = 2, col = "lightblue")

polygon(x = c(xt, rev(xt)), y = with(logr_cll_pred_2d, c(fit - 2 * se.fit, rev(fit + 2 * se.fit))), border = NA, col = adjustcolor("green", alpha.f = 0.3))
lines(xt, logr_cll_pred_2d$fit, lty = 1, col = "darkgreen")

rug(jitter(grad2[veg2 == 0]), side = 1, col = adjustcolor("lightblue", alpha.f = 0.1))
rug(jitter(grad2[veg2 == 1]), side = 3, col = adjustcolor("orange", alpha.f = 0.1))
temp <- density(grad2[veg2 == 0]); lines(temp$x, temp$y * 0.1 / max(temp$y), col = "lightblue")
temp <- density(grad2[veg2 == 1]); lines(temp$x, 1 - temp$y * 0.1 / max(temp$y), col = "orange")
lines(x1, z1, col = "red")
	
	# choice between linear vs sigmoidal only reflects how transect is defined and not necessarily how abrupt it is, i.e., whether 'shoulders' of vegetation density at low and high values are included in transect or not
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


tabulate_Danz2012_abruptness <- function(etable, b, data, flag_migtype) {

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
	
	sets <- list(c("ElevVsDist", "Elev", "vsDist"),
				c("Veg1VsDist", "Veg1", "vsDist"), c("Veg1VsElev", "Veg1", "vsElev"),
				c("Veg2VsDist", "Veg2", "vsDist"), c("Veg2VsElev", "Veg2", "vsElev"))
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
Danz2012JVegSci <- function(i, b, migtype, ecotoner_settings, etband, etmeasure, flag_bfig, copy_FromMig1_TF, do_figures, ...) {
	#---3. Ecological boundaries
	#3a. Danz et al. 2012 J.Veg.Sci.: Shape of vegetation boundary in relation to environmental conditions
	#Objective: of our boundary analysis was to evaluate whether the transition from prairie to forest across the boundary resulted from a smooth or abrupt climatic gradient, i.e. whether the transition followed pattern ‘(a)’ or pattern ‘(b)’ in Fig. 3. We used three analytical tactics to address this objective:
	#	(1) description of the spatial pattern of vegetation transition across the boundary,
	#	(2) evaluation of whether the climate gradient P – PET followed a steeper or shallower pattern of change across the boundary, and
	#	(3) direct modelling of the vegetation–climate relationship across the boundary.
	#Here, use elevation as surrogate driving environmental variable

	dots <- list(...)
	seed <- if ("seed" %in% names(dots)) dots["seed"] else NULL

	etmeasure$etable[b, "Transect_ID"] <- i
	etmeasure$etable[b, "Neighbor_Cells"] <- neighborhoods(ecotoner_settings)[b]
	
	if (!copy_FromMig1_TF) {
		etmeasure$gETmeas[[b]][[migtype]]$ElevVsDist <- calc_Danz2012_abruptness_1D(x = etband$Env$DistAlongXaxis_m,
															z = etband$Env$elev$YMeans_ForEachX, seed = seed)
		etmeasure$gETmeas[[b]][[migtype]]$Veg1VsDist <- calc_Danz2012_abruptness_1D(x = etband$Env$DistAlongXaxis_m,
															z = etband$Veg[[migtype]]$Veg1$density, seed = seed)
		etmeasure$gETmeas[[b]][[migtype]]$Veg2VsDist <- calc_Danz2012_abruptness_1D(x = etband$Env$DistAlongXaxis_m,
															z = etband$Veg[[migtype]]$Veg2$density, seed = seed)
		etmeasure$gETmeas[[b]][[migtype]]$Veg1VsElev <- calc_Danz2012_abruptness_1D(x = etband$Env$elev$YMeans_ForEachX,
															z = etband$Veg[[migtype]]$Veg1$density, seed = seed)
		etmeasure$gETmeas[[b]][[migtype]]$Veg2VsElev <- calc_Danz2012_abruptness_1D(x = etband$Env$elev$YMeans_ForEachX,
															z = etband$Veg[[migtype]]$Veg2$density, seed = seed)

if (FALSE) {
		grid_dist <- raster::raster(etband$Env$elev$grid)
		grid_dist[] <- rep(etband$Env$DistAlongXaxis_m, times = bandTransect_width_cellN(ecotoner_settings))
		
		etmeasure$gETmeas[[b]][[migtype]]$Veg1VsDist_2D <- calc_Danz2012_abruptness_2D(
															x1 = etband$Env$DistAlongXaxis_m, z1 = etband$Veg[[migtype]]$Veg1$density,
															x2 = grid_dist, z2 = etband$Veg[[migtype]]$Veg1$grid,
															seed = seed)
		etmeasure$gETmeas[[b]][[migtype]]$Veg2VsDist_2D <- calc_Danz2012_abruptness_2D(
															x1 = etband$Env$DistAlongXaxis_m, z1 = etband$Veg[[migtype]]$Veg2$density,
															x2 = grid_dist, z2 = etband$Veg[[migtype]]$Veg2$grid,
															seed = seed)
		etmeasure$gETmeas[[b]][[migtype]]$Veg1VsElev_2D <- calc_Danz2012_abruptness_2D(
															x1 = etband$Env$elev$YMeans_ForEachX, z1 = etband$Veg[[migtype]]$Veg1$density,
															x2 = etband$Env$elev$grid, z2 = etband$Veg[[migtype]]$Veg1$grid,
															seed = seed)
		etmeasure$gETmeas[[b]][[migtype]]$Veg2VsElev_2D <- calc_Danz2012_abruptness_2D(
															x1 = etband$Env$elev$YMeans_ForEachX, z1 = etband$Veg[[migtype]]$Veg2$density,
															x2 = etband$Env$elev$grid, z2 = etband$Veg[[migtype]]$Veg2$grid,
															seed = seed)
}


		if (do_figures) plot_Danz2012_abruptness(filename = file.path(dir_fig, paste0(flag_bfig, "Danz2012_FittedBoundaryShapes_", migtype, ".pdf")),
												eB_Env = etband$Env,
												eB_Veg = etband$Veg[[migtype]],
												datFit = etmeasure$gETmeas[[b]][[migtype]])				
	} 
	
	etmeasure$etable <- tabulate_Danz2012_abruptness(etable = etmeasure$etable, b = b,
														data = etmeasure$gETmeas[[b]][[if (copy_FromMig1_TF) "AllMigration" else migtype]],
														flag_migtype = migtype)

	etmeasure	
}
