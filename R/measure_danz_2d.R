#TODO(drs): pass 'doms' etc via an option method to Danz2012; same approach to set options from user for other methods too; account for options in 'what_measure'
#TODO(drs): create S4 class for measure methods

version_Danz2012JVegSci_2D <- function() numeric_version("0.2.1")

#' Estimation of the abruptness of the change of outcome z along the \sQuote{rows} of 1- and 2-dimensional x
#'
#' Implementation of the methods by Danz et al. (2012). 
#'
#' \code{doms} contains the named items 
#'		\code{tag}: descriptive character string of the model; 
#'		\code{cond}: character string that evaluates to a logical value expressing under which conditions the model should be calculated; 
#'		\code{fun}: model function as character string; 
#'		\code{family}: error distribution function as character string (see \code{\link[stats]{family}}); 
#'		\code{link}: model link function as character string (see \code{\link[stats]{make.link}}); 
#'		\code{random}: a character string representing the random effects if any (see \code{\link[lme4]{glmer}} for calls to \code{m_glmm_lme4}, respectively \code{\link[MASS]{glmmPQL}} for calls to \code{m_glmm_PQL}); 
#'		\code{correlation}: \code{NULL} or a character string representing the residual correlation structure for \code{m_glmm_PQL} (see \code{\link[MASS]{glmmPQL}}); 
#'		\code{ytrans}: optional - a function applied to outcome z before model fitting; 
#'		\code{ytransinv}: optional - the backtransformation of \code{ytrans} applied after model fitting.
#'
#' @param doms A named list of named options as a list that defines models which are applied to the data. See details.
#' @param use_dims A logical vector of length 2. Value of first/second element indicates if 1/2-dimensional models are estimated.
#' @param x1d A numeric vector. The 1-dimensional predictor along the ecotone transect from left-right. The length of x has to correspond to the ecotone transect. For instance, distance along the transect.
#' @param z1d A numeric vector. The 1-dimensional outcome along the ecotone transect, e.g., vegetation density.
#' @param x2d A \code{\link{matrix}} or \linkS4class{RasterLayer}. The 2-dimensional predictor representing the ecotone transect band.
#' @param z2d A \code{\link{matrix}} or \linkS4class{RasterLayer}. The 2-dimensional outcome representing the ecotone transect band.
#' @param seed An integer value, \code{NULL}, or \code{NA}: \code{NA} will do nothing; the other two are passed to \code{\link{set.seed}}.
#' @param include_lm A logical value. If \code{TRUE}, then a linear model is added to \code{doms}.
#'
#' @references Danz, N. P., L. E. Frelich, P. B. Reich, and G. J. Niemi. 2012. Do vegetation boundaries display smooth or abrupt spatial transitions along environmental gradients? Evidence from the prairie-forest biome boundary of historic Minnesota, USA. Journal of Vegetation Science 24:1129-1140.
#' @seealso \code{\link{Danz2012JVegSci_2D}}, and the possible model functions \code{\link[ecotoner]{m_glm}}, \code{\link[ecotoner]{m_glmm_lme4}}, and \code{\link[ecotoner]{m_sig}}
#' @return A named list with three items \sQuote{fits}, \sQuote{plot_preds}, and \sQuote{plot_data}
#' @export
calc_Danz2012_abruptness_2D <- function(doms, use_dims, x1d, z1d, x2d, z2d, seed = NULL, include_lm = FALSE) {
	# logistic regression on z ~ x with repeated 'rows' along the transect
	# idea by ld; implementation by drs

	# Literature:
	#	Zuur, A. F., E. N. Ieno, N. J. Walker, A. A. Savellev, and G. M. Smith. 2009. Mixed effects models and extensions in ecology with R. Springer.
	#		- binomial GLM: "We assume that Yi is binomial distributed with probability pi and ni = 1 independent trials"
	#		- "The logit and probit link functions assume that you have approximately an equal number of zeros and ones"
	#		- "The clog-log may be an option if you have considerably more zeros than ones, or vice versa; the sigmoidal curve is asymmetrical."

	#---Test input
	if (!is.na(seed)) set.seed(seed)

	if (!(length(x1d) == length(z1d))) {
		stop("ecotoner::calc_Danz2012_abruptness_2D(): arguments 'x1d' and 'z1d' describe the same transect and thus must have the same length")
	}
	if (!(inherits(x2d, "RasterLayer") && inherits(z2d, "RasterLayer")) && !(inherits(x2d, "matrix") && inherits(z2d, "matrix"))) {
		stop("ecotoner::calc_Danz2012_abruptness_2D(): both argument(s) 'x2d' and 'z2d' must be either matrices or rasters")
	}
	if (inherits(x2d, "RasterLayer") && inherits(x2d, "RasterLayer") && !raster::compareRaster(x2d, z2d)) {
		stop("ecotoner::calc_Danz2012_abruptness_2D(): argument(s) 'x2d' and 'z2d' describe the same transect and must be comparable rasters")
	} else
	if (inherits(x2d, "matrix") && inherits(z2d, "matrix") && !(identical(dim(x2d), dim(z2d)) && ncol(x2d) >= 3)) {
		stop("ecotoner::calc_Danz2012_abruptness_2D(): argument(s) 'x2d' and 'z2d' describe the same transect and must be comparable matrices with at least 3 columns: x, y, value")
	}	
	
	#---Begin function calculations
	# Interpret NAs as absences of iveg
	if (inherits(z2d, "RasterLayer")) {
		z2d <- raster::calc(z2d, function(x) ifelse(is.na(x), 0, x))
	} else {
		z2d[is.na(z2d)] <- 0
	}
	dat2d <- transect_to_long(x2d, z2d)
	if (nrow(dat2d) == 0 || all(is.na(dat2d))) {
		stop("ecotoner::calc_Danz2012_abruptness_2D(): arguments 'x2d' and 'z2d' combined contain only NAs")
	}	

	# model fits
	xt <- seq(min(dat2d[, "x"], na.rm = TRUE), max(dat2d[, "x"], na.rm = TRUE), length.out = 100)
	x1d_scaled <- scale(x1d)
	x2d_scaled <- scale(dat2d[, "x"])

	dats <- list('1D' = list(x = as.numeric(x1d_scaled), y = z1d,
							r = rep(NA, length(x1d)), c = rep(NA, length(x1d)), 
							w = rep(length(unique(dat2d[, "rows"])), length(x1d)),
							newdata = data.frame(x = scale(xt, center = attr(x1d_scaled, "scaled:center"), scale = attr(x1d_scaled, "scaled:scale"))),
							is_binary = is_binary(z1d),
							scaled_scale = attr(x1d_scaled, "scaled:scale")),
				 '2D' = list(x = as.numeric(x2d_scaled), y = dat2d[, "y"],
				 			r = factor(dat2d[, "rows"]), c = factor(dat2d[, "cols"]),
				 			newdata = data.frame(x = scale(xt, center = attr(x2d_scaled, "scaled:center"), scale = attr(x2d_scaled, "scaled:scale"))),
				 			is_binary = is_binary(dat2d[, "y"]),
							scaled_scale = attr(x2d_scaled, "scaled:scale"))
				)
	if (!is.null(use_dims)) dats <- dats[use_dims]
	plot_data <- list('1D' = data.frame(x = x1d, y = z1d), '2D' = as.data.frame(dat2d[, c("x", "y")]), newdata = xt)
	
	# GLM
	# 	Model: P(Y_ij = 1) = link(b_0 + b_1 * X_ij) with X_ij = environmental gradient
	# 	Assumption: Y_ij are independent with i = along transect and j = 'rows' along width of transect band
	# GLMM
	#	Model: random intercept of rows representing repeated measures of one-cell wide transects
	#			P(Y_ij = 1|B_j) = link(b_0 + b_1 * X_ij + b_2 * B_j) with B_j = rows as blocks
	
	if (is.null(doms)) 
		doms <- list(lGLM = list(tag = "logistic GLM", cond = "TRUE", fun = "m_glm", family = "binomial", link = "logit"))
	if (include_lm) doms <- modifyList(list(lm = list(tag = "LM", cond = "TRUE", fun = "m_glm", family = "gaussian", link = "identity")), doms)
	
	grid <- if (any(sapply(doms, function(m) identical(m$fun, "m_glmm_RAC"))) && inherits(z2d, "RasterLayer")) z2d else NULL
	
	preds <- fits <- vector("list", length = length(dats))
	names(preds) <- names(fits) <- names(dats)
	
	for (k in seq_along(dats)) {
		preds[[k]] <- fits[[k]] <- list()

		for (j in seq_along(doms)) if (eval(parse(text = doms[[j]][["cond"]]))) {
			# call function
			fargs <- list(data. = dats[[k]])
			if (!is.null(doms[[j]][["family"]]))
				fargs <- c(fargs, list(family = match.fun(doms[[j]][["family"]])(link = doms[[j]][["link"]])))
			fargs <- c(fargs, list(random = doms[[j]][["random"]],
									correlation = doms[[j]][["correlation"]],
									ytrans = doms[[j]][["ytrans"]], ytransinv = doms[[j]][["ytransinv"]],
									grid = grid))
			res <- do.call(doms[[j]][["fun"]], args = fargs)
			if (inherits(res[["m"]], "try-error"))
				stop("ecotoner::calc_Danz2012_abruptness_2D(): ", res[["m"]])

			# copy result
			preds[[k]][[doms[[j]][["tag"]]]] <- modifyList(res[["preds"]], list(isConv = res[["quals"]][["isConv"]]))
			fits[[k]][[doms[[j]][["tag"]]]] <- res[c("quals", "coef1", "perf")]
			
			# rescale coef1 to original scale
			cres <- fits[[k]][[doms[[j]][["tag"]]]][["coef1"]]
			if (length(cres) > 0 && !anyNA(cres)) fits[[k]][[doms[[j]][["tag"]]]][["coef1"]] <- cres / dats[[k]][["scaled_scale"]]
		}
	}
		
	list(fits = fits, plot_preds = preds, plot_data = plot_data)
}

	
add_Danz2012_abruptness_2D_panel <- function(preds, data, end_toLeft, xlab, ylab = NULL, panel = NULL, x_ann = TRUE, y_ann = TRUE, add_legend = TRUE, cex) {
	m_N <- sapply(preds, length)
	m_labels <- paste(rep(names(preds), m_N), unlist(lapply(preds, names)))
	m_ltys <- unlist(lapply(m_N, seq_len))
	if (length(m_N) > 1) {
		ncol1 <- min(m_N)
		ncol2 <- abs(diff(m_N))
	} else {
		ncol1 <- floor(m_N / 3)
		ncol2 <- m_N - 3 * ncol1
	}
	m_cols <- list(line = c(rep(c("blue", "darkgreen", "black"), each = ncol1), rep("purple", ncol2)),
					range = c(rep(c("dodgerblue", "forestgreen", "darkgray"), each = ncol1), rep("orchid", ncol2)))

	# Start plot and add observed mean data
	plot(data[["1D"]][["x"]], data[["1D"]][["y"]],
			xlim = c(min(data[["newdata"]]), max(data[["newdata"]])), ylim = c(0, 1), 
			xlab = "", ylab = "", type = "l", col = "red", axes = FALSE)
	
	# Add more observed data
	index_y0 <- data[["2D"]][["y"]] == 0
	sum_y0 <- sum(index_y0)
	index_y1 <- !index_y0
	sum_y1 <- length(index_y0) - sum_y0
	if (sum_y0 > 0) rug(jitter(data[["2D"]][["x"]][index_y0]), side = 1, col = adjustcolor("lightblue", alpha.f = 0.1), quiet = TRUE)
	if (sum_y1 > 0) rug(jitter(data[["2D"]][["x"]][index_y1]), side = 3, col = adjustcolor("orange", alpha.f = 0.1), quiet = TRUE)
	if (sum_y0 > 2) {
		temp <- density(data[["2D"]][["x"]][index_y0])
		lines(temp$x, temp$y * 0.1 / max(temp$y), col = "lightblue")
	}
	if (sum_y1 > 2) {
		temp <- density(data[["2D"]][["x"]][index_y1])
		lines(temp$x, 1 - temp$y * 0.1 / max(temp$y), col = "orange")
	}

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
	if (missing(cex)) cex <- par("cex")
	axis(side = 1, labels = x_ann)
	axis(side = 2, labels = y_ann)
	if (!is.null(panel)) mtext(side = 3, adj = 0.025, text = panel, font = 2, xpd = NA)
	if (x_ann) mtext(side = 1, text = xlab, line = par("mgp")[1], cex = cex, xpd = NA)
	if (y_ann) mtext(side = 2, text = if (is.null(ylab)) "Response" else ylab, line = par("mgp")[1], cex = cex, xpd = NA)
	
	if (add_legend) {
		x_pos <- if (any(grepl("1D LM", m_labels))) {
					if (end_toLeft) "topleft" else "bottomleft"
				 } else {
					if (end_toLeft) "topleft" else "topright"
				 }
		legend(x = x_pos, inset = c(0, 0.025), bty = "n",
				legend = m_labels, cex = cex * 0.5,
				lty = m_ltys, col = m_cols[["line"]],
				pch = 22, pt.cex = 1.25, pt.lwd = 0, pt.bg = adjustcolor(m_cols[["range"]], alpha.f = 0.3))
	}
	
	invisible()
}
	
plot_Danz2012_abruptness_2D <- function(filename, xlab, preds1, preds2, data1, data2, end_toLeft1, end_toLeft2, responses_equal = TRUE) {
	pdf(height = 4.5 + 0.5, width = 2 * 5 + 0.5, file = filename)
	layout(mat = matrix(c(0, 1, 2, 0, 0, 0), nrow = 2, ncol = 3, byrow = TRUE),
			heights = c(4.5, 0.5), widths = c(0.5, 5, 5))
	par_old <- par(mar = c(0.5, 0.5, 1, 0.5), mgp = c(1.5, 0.5, 0), cex = 1.25, xpd = NA)
	on.exit({par(par_old); dev.off()}, add = TRUE)

	add_Danz2012_abruptness_2D_panel(preds = preds1, data = data1, end_toLeft = end_toLeft1, xlab = xlab, panel = "(a)",
				add_legend = responses_equal)
	add_Danz2012_abruptness_2D_panel(preds = preds2, data = data2, end_toLeft = end_toLeft2, xlab = xlab, panel = "(b)",
				y_ann = !responses_equal,
				add_legend = !identical(lapply(preds1, names), lapply(preds2, names)) || !responses_equal)
		
	invisible()
}


tabulate_Danz2012_abruptness_2D <- function(etable, index, data) {
	data <- unlist(data)
	names(data) <- paste("Danz2012", names(data), sep = "_")

	tabulate_merge_into_etable(etable, index, data)
}



#' Estimation of abruptness of the vegetation-transition along an environmental gradient of an ecotone transect
#'
#' Implementation of the methods by Danz et al. (2012). 
#'
#' @param i An integer value. The number of the ecotone transect; a value smaller than or equal to \code{transect_N(ecotoner_settings)}.
#' @param b An integer value. The index of the neighborhood window size; ; a value smaller than or equal to \code{neighbors_N(ecotoner_settings)}.
#' @param migtype A character string; one of \code{get("migtypes", envir = etr_vars)}.
#' @param ecotoner_settings An object of the class \linkS4class{EcotonerSettings}.
#' @param etband The object of the ecotone transect according to index \code{i} that was returned by the function \code{\link[ecotoner]{detect_ecotone_transects_from_searchpoint}}.
#' @param etmeasure The object of the ecotone transect according to index \code{i} that is created by the function \code{\link[ecotoner]{detect_ecotone_transects_from_searchpoint}}.
#' @param copy_FromMig1_TF A logical value. If \code{TRUE} and \code{migtype} is not \dQuote{AllMigration}, then output is copied from \dQuote{AllMigration} instead of re-calculated.
#' @param do_figures A logical value. If \code{TRUE}, then figures are generated and saved to disk.
#' @param ... Further arguments. Currently implemented are \code{seed} (to seed the random number generator), \code{flag_bfig} (character string to include in file names of figures), and \code{dir_fig} (character string of the path where figures are to be stored on disk).
#'
#' @return The object \code{etmeasure} with output of the call to \code{\link{calc_Danz2012_abruptness_2D}} based on \code{doms}, \code{b}, and \code{migtype} integrated.
#' @seealso \code{\link{calc_Danz2012_abruptness_2D}}
#' @references Danz, N. P., L. E. Frelich, P. B. Reich, and G. J. Niemi. 2012. Do vegetation boundaries display smooth or abrupt spatial transitions along environmental gradients? Evidence from the prairie-forest biome boundary of historic Minnesota, USA. Journal of Vegetation Science 24:1129-1140.
#' @export
Danz2012JVegSci_2D <- function(i, b, migtype, ecotoner_settings, etband, etmeasure, copy_FromMig1_TF, do_figures, ...) {
	dots <- list(...)
	seed <- if ("seed" %in% names(dots)) dots["seed"] else NULL
	if ("flag_bfig" %in% names(dots)) flag_bfig <- dots["flag_bfig"] else do_figures <- FALSE
	if ("dir_fig" %in% names(dots)) dir_fig <- dots["dir_fig"] else do_figures <- FALSE
	
	b_migtype <- (b - 1) * length(get("migtypes", envir = etr_vars)) + which(migtype == get("migtypes", envir = etr_vars))
	etmeasure$etable[b_migtype, "Transect_ID"] <- i
	etmeasure$etable[b_migtype, "Neighbor_Cells"] <- neighborhoods(ecotoner_settings)[b]
	etmeasure$etable[b_migtype, "Migration_Type"] <- migtype
	
	etmeasure$gETmeas[[b]][[migtype]]$meta <- list()
	etmeasure$gETmeas[[b]][[migtype]]$meta$method <- "Danz2012JVegSci_2D"
	etmeasure$gETmeas[[b]][[migtype]]$meta$version <- version_Danz2012JVegSci_2D()
	etmeasure$gETmeas[[b]][[migtype]]$meta$copy_FromMig1_TF <- copy_FromMig1_TF
	
	if (!copy_FromMig1_TF) {
		# Vegetation versus elevation
		doms <- list(	lm = list(tag = "LM", cond = "TRUE", fun = "m_glm", family = "gaussian", link = "identity")
						, lGLM = list(tag = "logistic GLM", cond = "TRUE", fun = "m_glm", family = "binomial", link = "logit")
						, bGLMc = list(tag = "binomial GLM with cloglog", cond = "TRUE", fun = "m_glm", family = "binomial", link = "cloglog")
						, bGLMci = list(tag = "binomial GLM with cloglog and 1 - y", cond = "TRUE", fun = "m_glm", family = "binomial", link = "cloglog", ytrans = function(x) 1 - x, ytransinv = function(x) 1 - x)
						, sNLS = list(tag = "sigmoidal NLS", cond = "TRUE", fun = "m_sig")
						, lGLMMr = list(tag = "logistic lme4-GLMM with random rows", cond = "!anyNA(dats[[k]][['r']])", fun = "m_glmm_lme4", family = "binomial", link = "logit", random = "(x|r)")
#						, lGLMMrc = list(tag = "logistic lme4-GLMM with random rows and columns", cond = "!anyNA(dats[[k]][['r']]) && !anyNA(dats[[k]][['c']])", fun = "m_glmm_lme4", family = "binomial", link = "logit", random = "(x|r) + (x|c)")
						, lGLMMPQLr = list(tag = "logistic PQL-GLMM with random rows", cond = "!anyNA(dats[[k]][['r']])", fun = "m_glmm_PQL", family = "binomial", link = "logit", random = "~ x|r")
#						, lGLMMPQLrSpher = list(tag = "logistic PQL-GLMM with random rows and spherical correlation in residuals", cond = "!anyNA(dats[[k]][['r']]) && !anyNA(dats[[k]][['c']])", fun = "m_glmm_PQL", family = "binomial", link = "logit", random = "~ x|r", correlation = "nlme::corSpher(form = ~ r + c, nugget = TRUE)")
						, lGLMMrRAC = list(tag = "logistic GLMM with random rows and RAC", cond = "!anyNA(dats[[k]][['r']])", fun = "m_glmm_RAC", family = "binomial", link = "logit", random = "(x|r)")
					)
		use_dims <- c('1D' = FALSE, '2D' = TRUE)

		temp1 <- calc_Danz2012_abruptness_2D(doms, use_dims,
					x1d = etband$Env$elev$YMeans_ForEachX, z1d = etband$Veg[[migtype]]$Veg1$density,
					x2d = etband$Env$elev$grid, z2d = etband$Veg[[migtype]]$Veg1$grid, seed = seed)
		temp2 <- calc_Danz2012_abruptness_2D(doms, use_dims,
					x1d = etband$Env$elev$YMeans_ForEachX, z1d = etband$Veg[[migtype]]$Veg2$density,
					x2d = etband$Env$elev$grid, z2d = etband$Veg[[migtype]]$Veg2$grid, seed = seed)
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
