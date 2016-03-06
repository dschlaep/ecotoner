
#------Danz, N.P., Frelich, L.E., Reich, P.B. & Niemi, G.J. (2012) Do vegetation boundaries display smooth or abrupt spatial transitions along environmental gradients? Evidence from the prairie–forest biome boundary of historic Minnesota, USA. Journal of Vegetation Science, n/a-n/a.



#---Danz et al. 2012: 'Spatial change in vegetation across the boundary', 'Spatial change in climate across the boundary', and 'Spatial vegetation/climate relationships across the boundary'
#' Estimates the abruptness of outcome y against x
#' 
#' @param x A numeric vector. Predictor along the ecotone transect from left-right. The length of x has to correspond to the ecotone transect. For instance, distance along the transect.
#' @param z A numeric vector. The outcome along the ecotone transect, e.g., vegetation density.
#' 
calc_Danz2012_abruptness_2D <- function(x1, z1, x2, z2, end_toLeft, xlab, fname, seed = NULL) {
	# logistic regression on z ~ x with repeated 'rows' along the transect
	# idea by ld; implementation by drs

	# Literature:
	#	Zuur, A. F., E. N. Ieno, N. J. Walker, A. A. Savellev, and G. M. Smith. 2009. Mixed effects models and extensions in ecology with R. Springer.
	#		- binomial GLM: "We assume that Yi is binomial distributed with probability pi and ni = 1 independent trials"
	#		- "The logit and probit link functions assume that you have approximately an equal number of zeros and ones"
	#		- "The clog–log may be an option if you have considerably more zeros than ones, or vice versa; the sigmoidal curve is asymmetrical."
	
	#---Test input
	if (!is.na(seed)) set.seed(seed)

	if (!(length(x1) == length(z1))) {
		stop("ecotoner::calc_Danz2012_abruptness_2D(): arguments 'x1' and 'z1' describe the same transect and thus must have the same length")
	}
	if (!(inherits(x2, "RasterLayer") && inherits(z2, "RasterLayer")) && !(inherits(x2, "matrix") && inherits(z2, "matrix"))) {
		stop("ecotoner::calc_Danz2012_abruptness_2D(): both argument(s) 'x2' and 'z2' must be either matrices or a rasters")
	}
	if (inherits(x2, "RasterLayer") && inherits(x2, "RasterLayer") && !raster::compareRaster(x2, z2)) {
		stop("ecotoner::calc_Danz2012_abruptness_2D(): argument(s) 'x2' and 'z2' describe the same transect and must be comparable rasters")
	} else
	if (inherits(x2, "matrix") && inherits(z2, "matrix") && !(identical(dim(x2), dim(z2)) && ncol(x2) >= 3)) {
		stop("ecotoner::calc_Danz2012_abruptness_2D(): argument(s) 'x2' and 'z2' describe the same transect and must be comparable matrices with at least 3 columns: x, y, value")
	}	
	
	#---Begin function calculations
	if (inherits(x2, "RasterLayer") && inherits(x2, "RasterLayer")) {
#TODO(drs): why are origins not exactly the same?
raster::origin(z2) <- raster::origin(x2)

		x2p <- raster::rasterToPoints(x2)
		if (raster::cellStats(z2, "countNA")) z2 <- raster::calc(z2, function(z) ifelse(is.na(z), 0, z))
		if (raster::cellStats(x2, "countNA")) z2 <- raster::mask(z2, x2)
		z2p <- raster::rasterToPoints(z2)
	
		veg2 <- z2p[, "layer"]
		grad2 <- x2p[, "layer"]
		rows2 <- x2p[, "y"]
	} else if (inherits(x2, "matrix") && inherits(z2, "matrix")) {
		isnotna <- !(is.na(z2[, 3]) | is.na(x2[, 3]))
		veg2 <- z2[isnotna, 3]
		grad2 <- x2[isnotna, 3]
		rows2 <- x2[isnotna, 2]
	}
	
	# model fits
	xt <- seq(min(grad2), max(grad2), length.out = 100)
	x1s <- scale(x1)
	grad2s <- scale(grad2)
	dats <- list('1D' = list(x = as.numeric(x1s), y = z1, r = rep(NA, length(x1)), w = rep(length(unique(rows2)), length(x1)),
							newdata = data.frame(x = scale(xt, center = attr(x1s, "scaled:center"), scale = attr(x1s, "scaled:scale")))),
				 '2D' = list(x = as.numeric(grad2s), y = veg2, r = factor(rows2),
				 			newdata = data.frame(x = scale(xt, center = attr(grad2s, "scaled:center"), scale = attr(grad2s, "scaled:scale"))))
				)
					
	
	fits <- lapply(dats, function(dat) {
				res <- list()
				# Model: P(Y_ij = 1) = link(b_0 + b_1 * X_ij) with X_ij = environmental gradient
				# Assumption: Y_ij are independent with i = along transect and j = 'rows' along width of transect band
#				res[["LM"]] <- m_glm(family = gaussian(link = "identity"), dat)
				res[["logistic GLM"]] <- m_glm(family = binomial(link = "logit"), dat)
				res[["logistic GLM with cloglog"]] <- m_glm(family = binomial(link = "cloglog"), dat)
				res[["sigmoidal"]] <- m_sig(dat)
				
				if (!anyNA(dat[["r"]])) {
					if (requireNamespace("lme4", quietly = TRUE)) {
						# Model: random intercept of rows representing repeated measures of one-cell wide transects
						#	P(Y_ij = 1|B_j) = link(b_0 + b_1 * X_ij + b_2 * B_j) with B_j = rows as blocks
						# Assumption: 
						res[["logistic GLMM"]] <- m_glmm(family = binomial(link = "logit"), dat)
					} else {
						warning("Package 'lme4' not installed: 'calc_Danz2012_abruptness_2D' will not estimate abruptness with GLMMs")
					}
				}
				
				res
			})
	
	plot_Danz2012_abruptness_2D(filename = fname, xt, x1, z1, grad2, veg2, fits, data[["newdata"]], end_toLeft, xlab)
	
	fits
}


#' @export
plot_Danz2012_abruptness_2D <- function(filename, xt, x1, z1, grad2, veg2, fits, newdata, end_toLeft, xlab) {
	m_N <- sapply(fits, length)
	m_labels <- paste(rep(names(fits), m_N), unlist(lapply(fits, names)))
	m_ltys <- unlist(lapply(m_N, seq_len))
	m_cols <- list(line = c(rep(c("blue", "black"), each = min(m_N)), rep("purple", abs(diff(m_N)))),
					range = c(rep(c("darkgray", "darkgreen"), each = min(m_N)), rep("orchid", abs(diff(m_N)))))
	
	pdf(width = 8, height = 6, file = filename)
	par_old <- par(mfrow = c(nr <- 1, 1), mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.5, 0.5, 0), cex = cex <- 1.5)
	on.exit({par(par_old); dev.off()}, add = TRUE)

	op <- par()
	plot(x1, z1, xlab = xlab, ylab = "Response", xlim = c(min(xt), max(xt)), ylim = c(0, 1), type = "l", col = "red")
	rug(jitter(grad2[veg2 == 0]), side = 1, col = adjustcolor("lightblue", alpha.f = 0.1), quiet = TRUE)
	rug(jitter(grad2[veg2 == 1]), side = 3, col = adjustcolor("orange", alpha.f = 0.1), quiet = TRUE)
	temp <- density(grad2[veg2 == 0]); lines(temp$x, temp$y * 0.1 / max(temp$y), col = "lightblue")
	temp <- density(grad2[veg2 == 1]); lines(temp$x, 1 - temp$y * 0.1 / max(temp$y), col = "orange")

	k <- 1
	temp <- lapply(fits, function(fit) lapply(fit, function(mfit) {
						if (mfit[["quals"]][["isConv"]] & !is.null(mfit[["preds"]])) {
							if (inherits(mfit[["preds"]], "list")) {
								lines(xt, mfit[["preds"]][["fit"]], lwd = 2, lty = m_ltys[k], col = m_cols[["line"]][k])
								polygon(x = c(xt, rev(xt)),
										y = c(mfit[["preds"]][["fit"]] - 2 * mfit[["preds"]][["se.fit"]], rev(mfit[["preds"]][["fit"]] + 2 * mfit[["preds"]][["se.fit"]])),
										border = NA, col = adjustcolor(m_cols[["range"]][k], alpha.f = 0.3))
							} else {
								lines(xt, mfit[["preds"]], lwd = 2, lty = m_ltys[k], col = m_cols[["line"]][k])
							}
						}
						k <<- k + 1
			}))
			
	legend(x = if (end_toLeft) "topleft" else "topright", inset = c(0, 0.025), bty = "n",
				legend = m_labels, cex = 0.75,
				lty = m_ltys, col = m_cols[["line"]],
				pch = 22, pt.cex = 1.5, pt.lwd = 0, pt.bg = adjustcolor(m_cols[["range"]], alpha.f = 0.3))
		
	invisible()
}



#' @export
Danz2012JVegSci_2D <- function(i, b, migtype, ecotoner_settings, etband, etmeasure, flag_bfig, copy_FromMig1_TF, do_figures, ...) {
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
		grid_dist <- raster::raster(etband$Env$elev$grid)
		grid_dist[] <- rep(etband$Env$DistAlongXaxis_m, times = bandTransect_width_cellN(ecotoner_settings))
		
		etmeasure$gETmeas[[b]][[migtype]]$Veg1VsDist_2D <- calc_Danz2012_abruptness_2D(
															x1 = etband$Env$DistAlongXaxis_m, z1 = etband$Veg[[migtype]]$Veg1$density,
															x2 = grid_dist, z2 = etband$Veg[[migtype]]$Veg1$grid,
															end_toLeft = end_to_left(type_veg1(ecotoner_settings)),
															xlab = "Distance along transect (m)", fname = file.path(dir_fig, paste0(flag_bfig, "Danz2012_FittedBoundaryShapes_", migtype, "_Veg1VsDist.pdf")),
															seed = seed)
		etmeasure$gETmeas[[b]][[migtype]]$Veg2VsDist_2D <- calc_Danz2012_abruptness_2D(
															x1 = etband$Env$DistAlongXaxis_m, z1 = etband$Veg[[migtype]]$Veg2$density,
															x2 = grid_dist, z2 = etband$Veg[[migtype]]$Veg2$grid,
															end_toLeft = end_to_left(type_veg2(ecotoner_settings)),
															xlab = "Distance along transect (m)", fname = file.path(dir_fig, paste0(flag_bfig, "Danz2012_FittedBoundaryShapes_", migtype, "_Veg2VsDist.pdf")),
															seed = seed)
		etmeasure$gETmeas[[b]][[migtype]]$Veg1VsElev_2D <- calc_Danz2012_abruptness_2D(
															x1 = etband$Env$elev$YMeans_ForEachX, z1 = etband$Veg[[migtype]]$Veg1$density,
															x2 = etband$Env$elev$grid, z2 = etband$Veg[[migtype]]$Veg1$grid,
															end_toLeft = end_to_left(type_veg1(ecotoner_settings)),
															xlab = "Elevation (m)", fname = file.path(dir_fig, paste0(flag_bfig, "Danz2012_FittedBoundaryShapes_", migtype, "_Veg1VsElev.pdf")),
															seed = seed)
		etmeasure$gETmeas[[b]][[migtype]]$Veg2VsElev_2D <- calc_Danz2012_abruptness_2D(
															x1 = etband$Env$elev$YMeans_ForEachX, z1 = etband$Veg[[migtype]]$Veg2$density,
															x2 = etband$Env$elev$grid, z2 = etband$Veg[[migtype]]$Veg2$grid,
															end_toLeft = end_to_left(type_veg2(ecotoner_settings)),
															xlab = "Elevation (m)", fname = file.path(dir_fig, paste0(flag_bfig, "Danz2012_FittedBoundaryShapes_", migtype, "_Veg2VsElev.pdf")),
															seed = seed)


		if (do_figures) plot_Danz2012_abruptness_1D(filename = file.path(dir_fig, paste0(flag_bfig, "Danz2012_FittedBoundaryShapes_", migtype, ".pdf")),
												eB_Env = etband$Env,
												eB_Veg = etband$Veg[[migtype]],
												datFit = etmeasure$gETmeas[[b]][[migtype]])				
	} 
	
	etmeasure$etable <- tabulate_Danz2012_abruptness_1D(etable = etmeasure$etable, b = b,
														data = etmeasure$gETmeas[[b]][[if (copy_FromMig1_TF) "AllMigration" else migtype]],
														flag_migtype = migtype)

	etmeasure	
}
