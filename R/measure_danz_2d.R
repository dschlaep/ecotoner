
#------Danz, N.P., Frelich, L.E., Reich, P.B. & Niemi, G.J. (2012) Do vegetation boundaries display smooth or abrupt spatial transitions along environmental gradients? Evidence from the prairie–forest biome boundary of historic Minnesota, USA. Journal of Vegetation Science, n/a-n/a.



#---Danz et al. 2012: 'Spatial change in vegetation across the boundary', 'Spatial change in climate across the boundary', and 'Spatial vegetation/climate relationships across the boundary'
#' Estimates the abruptness of outcome y against x
#' 
#' @param x A numeric vector. Predictor along the ecotone transect from left-right. The length of x has to correspond to the ecotone transect. For instance, distance along the transect.
#' @param z A numeric vector. The outcome along the ecotone transect, e.g., vegetation density.
#' 
calc_Danz2012_abruptness_2D <- function(x1, z1, x2, z2, seed = NULL) {
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
#TODO(drs): should I use data-splitting (e.g., cross-validation) to estimate model performance?
	dats <- list('1D' = list(x = as.numeric(x1s), y = z1, r = rep(NA, length(x1)), w = rep(length(unique(rows2)), length(x1)),
							newdata = data.frame(x = scale(xt, center = attr(x1s, "scaled:center"), scale = attr(x1s, "scaled:scale"))),
							is_binary = is_binary(z1),
							scaled_scale = attr(x1s, "scaled:scale")),
				 '2D' = list(x = as.numeric(grad2s), y = veg2, r = factor(rows2),
				 			newdata = data.frame(x = scale(xt, center = attr(grad2s, "scaled:center"), scale = attr(grad2s, "scaled:scale"))),
				 			is_binary = is_binary(veg2),
							scaled_scale = attr(grad2s, "scaled:scale"))
				)
	plot_data <- list('1D' = list(x = x1, y = z1), '2D' = list(x = grad2, y = veg2), newdata = xt)
	
	# GLM
	# 	Model: P(Y_ij = 1) = link(b_0 + b_1 * X_ij) with X_ij = environmental gradient
	# 	Assumption: Y_ij are independent with i = along transect and j = 'rows' along width of transect band
	# GLMM
	#	Model: random intercept of rows representing repeated measures of one-cell wide transects
	#			P(Y_ij = 1|B_j) = link(b_0 + b_1 * X_ij + b_2 * B_j) with B_j = rows as blocks
	
	dom <- list(#list(tag = "LM", cond = "TRUE", fun = "m_glm", family = "gaussian", link = "identity"),
				list(tag = "logistic GLM", cond = "TRUE", fun = "m_glm", family = "binomial", link = "logit"),
				list(tag = "logistic GLM with cloglog", cond = "TRUE", fun = "m_glm", family = "binomial", link = "cloglog"),
				list(tag = "sigmoidal NLS", cond = "TRUE", fun = "m_sig", family = NA, link = NA),
				list(tag = "logistic GLMM", cond = "!anyNA(dats[[k]][['r']])", fun = "m_glmm", family = "binomial", link = "logit")
				)
	
	preds <- fits <- vector("list", length = length(dats))
	names(preds) <- names(fits) <- names(dats)
	
	for (k in seq_along(dats)) {
		preds[[k]] <- fits[[k]] <- list()

		for (j in seq_along(dom)) if (eval(parse(text = dom[[j]][["cond"]]))) {
			f <- match.fun(dom[[j]][["fun"]])
			if ("family" %in% names(formals(f))) {
				temp <- f(family = match.fun(dom[[j]][["family"]])(link = dom[[j]][["link"]]), dats[[k]])
			} else {
				temp <- f(dats[[k]])
			}
			
			preds[[k]][[dom[[j]][["tag"]]]] <- modifyList(temp[["preds"]], list(isConv = temp[["quals"]][["isConv"]]))
			fits[[k]][[dom[[j]][["tag"]]]] <- temp[c("quals", "coef1", "perf")]
			
			# rescale coef1 to original scale
			temp <- fits[[k]][[dom[[j]][["tag"]]]][["coef1"]]
			if (length(temp) > 0 && !anyNA(temp)) fits[[k]][[dom[[j]][["tag"]]]][["coef1"]] <- temp / dats[[k]][["scaled_scale"]]
		}
	}
		
	list(fits = fits, plot_preds = preds, plot_data = plot_data)
}

	
add_Danz2012_abruptness_2D_panel <- function(preds, data, end_toLeft, xlab, panel = NULL, x_ann = TRUE, y_ann = TRUE, add_legend = TRUE) {
	m_N <- sapply(preds, length)
	m_labels <- paste(rep(names(preds), m_N), unlist(lapply(preds, names)))
	m_ltys <- unlist(lapply(m_N, seq_len))
	m_cols <- list(line = c(rep(c("blue", "black"), each = min(m_N)), rep("purple", abs(diff(m_N)))),
					range = c(rep(c("darkgray", "darkgreen"), each = min(m_N)), rep("orchid", abs(diff(m_N)))))

	# Start plot and add observed mean data
	plot(data[["1D"]][["x"]], data[["1D"]][["y"]],
			xlim = c(min(data[["newdata"]]), max(data[["newdata"]])), ylim = c(0, 1), 
			xlab = if (x_ann) xlab else "", ylab = if (y_ann) "Response" else "",
			type = "l", col = "red", axes = FALSE)
	axis(side = 1, labels = x_ann)
	axis(side = 2, labels = y_ann)
	
	# Add more observed data
	rug(jitter(data[["2D"]][["x"]][data[["2D"]][["y"]] == 0]), side = 1, col = adjustcolor("lightblue", alpha.f = 0.1), quiet = TRUE)
	rug(jitter(data[["2D"]][["x"]][data[["2D"]][["y"]] == 1]), side = 3, col = adjustcolor("orange", alpha.f = 0.1), quiet = TRUE)
	temp <- density(data[["2D"]][["x"]][data[["2D"]][["y"]] == 0]); lines(temp$x, temp$y * 0.1 / max(temp$y), col = "lightblue")
	temp <- density(data[["2D"]][["x"]][data[["2D"]][["y"]] == 1]); lines(temp$x, 1 - temp$y * 0.1 / max(temp$y), col = "orange")

	# Add model predictions
	k <- 1
	temp <- lapply(preds, function(fit) lapply(fit, function(mfit) {
						if (as.logical(mfit[["isConv"]])) {
							lines(data[["newdata"]], mfit[["fit"]], lwd = 2, lty = m_ltys[k], col = m_cols[["line"]][k])
							if (!is.null(mfit[["se.fit"]])) {
								polygon(x = c(data[["newdata"]], rev(data[["newdata"]])),
										y = c(mfit[["fit"]] - 2 * mfit[["se.fit"]], rev(mfit[["fit"]] + 2 * mfit[["se.fit"]])),
										border = NA, col = adjustcolor(m_cols[["range"]][k], alpha.f = 0.3))
							}
						}
						k <<- k + 1
			}))
	
	# Add annotations
	if (!is.null(panel)) mtext(side = 3, adj = 0.025, text = panel, font = 2)
	
	if (add_legend) {
		legend(x = if (end_toLeft) "topleft" else "topright", inset = c(0, 0.025), bty = "n",
				legend = m_labels, cex = par("cex") * 0.75,
				lty = m_ltys, col = m_cols[["line"]],
				pch = 22, pt.cex = 1.5, pt.lwd = 0, pt.bg = adjustcolor(m_cols[["range"]], alpha.f = 0.3))
	}
	
	invisible()
}
	
#' @export
plot_Danz2012_abruptness_2D <- function(filename, xlab, preds1, preds2, data1, data2, end_toLeft1, end_toLeft2) {
	pdf(height = 4.5 + 0.5, width = 2 * 5 + 0.5, file = filename)
	layout(mat = matrix(c(0, 1, 2, 0, 0, 0), nrow = 2, ncol = 3, byrow = TRUE),
			heights = c(4.5, 0.5), widths = c(0.5, 5, 5))
	par_old <- par(mar = c(0.5, 0.5, 1, 0.5), mgp = c(1.5, 0.5, 0), cex = cex <- 1.25, xpd = NA)
	on.exit({par(par_old); dev.off()}, add = TRUE)

	add_Danz2012_abruptness_2D_panel(preds = preds1, data = data1, end_toLeft = end_toLeft1, xlab = xlab, panel = "(a)")
	add_Danz2012_abruptness_2D_panel(preds = preds2, data = data2, end_toLeft = end_toLeft2, xlab = xlab, panel = "(b)",
				x_ann = TRUE, y_ann = FALSE,
				add_legend = !identical(lapply(preds1, names), lapply(preds2, names)))
		
	invisible()
}


tabulate_Danz2012_abruptness_2D <- function(etable, index, data) {
#TODO(drs): consider replacing this function with 'put.ListData1Level_TableData'
	temp <- unlist(data)
	cname <- paste("Danz2012", names(temp), sep = "_")
	cname_exist <- match(cname, colnames(etable), nomatch = 0)

	icol <- cname_exist > 0
	if (any(icol)) etable[index, cname_exist[icol]] <- temp[icol]

	icol <- cname_exist == 0
	if (any(icol)) {
		res <- as.data.frame(matrix(NA, nrow = max(1, nrow(etable)), ncol = sum(icol), dimnames = list(NULL, cname[icol])))
		res[index, ] <- temp[icol]
		etable <- cbind(etable, res)
	}
	
	etable
}



#' @export
Danz2012JVegSci_2D <- function(i, b, migtype, ecotoner_settings, etband, etmeasure, copy_FromMig1_TF, do_figures, ...) {
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
	
	b_migtype <- (b - 1) * length(get("migtypes", envir = etr_vars)) + which(migtype == get("migtypes", envir = etr_vars))
	etmeasure$etable[b_migtype, "Transect_ID"] <- i
	etmeasure$etable[b_migtype, "Neighbor_Cells"] <- neighborhoods(ecotoner_settings)[b]
	etmeasure$etable[b_migtype, "Migration_Type"] <- migtype
	
	if (!copy_FromMig1_TF) {
		# Vegetation versus elevation
		temp1 <- calc_Danz2012_abruptness_2D(x1 = etband$Env$elev$YMeans_ForEachX, z1 = etband$Veg[[migtype]]$Veg1$density,
												x2 = etband$Env$elev$grid, z2 = etband$Veg[[migtype]]$Veg1$grid, seed = seed)
		temp2 <- calc_Danz2012_abruptness_2D(x1 = etband$Env$elev$YMeans_ForEachX, z1 = etband$Veg[[migtype]]$Veg2$density,
												x2 = etband$Env$elev$grid, z2 = etband$Veg[[migtype]]$Veg2$grid, seed = seed)
		etmeasure$gETmeas[[b]][[migtype]]$Veg1VsElev_2D <- temp1$fits
		etmeasure$gETmeas[[b]][[migtype]]$Veg2VsElev_2D <- temp2$fits
		
		if (do_figures) plot_Danz2012_abruptness_2D(filename = file.path(dir_fig, paste0(flag_bfig, "Danz2012_FittedBoundaryShapes_2D_", migtype, "_VegsVsElev.pdf")),
													xlab = "Elevation (m)",
													preds1 = temp1$plot_preds, preds2 = temp2$plot_preds,
													data1 = temp1$plot_data, data2 = temp2$plot_data,
													end_toLeft1 = end_to_left(type_veg1(ecotoner_settings)), end_toLeft2 = end_to_left(type_veg2(ecotoner_settings)))
	} 
	
	etmeasure$etable <- tabulate_Danz2012_abruptness_2D(etable = etmeasure$etable, index = b_migtype,
														data = etmeasure$gETmeas[[b]][[if (copy_FromMig1_TF) "AllMigration" else migtype]])

	etmeasure	
}
